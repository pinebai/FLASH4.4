!!****if* source/physics/sourceTerms/Heatexchange/HeatexchangeMain/Spitzer/Heatexchange
!!
!! NAME
!!
!!  Heatexchange
!!
!!
!! SYNOPSIS
!!
!!   call Heatexchange ( integer(IN) :: blockCount, 
!!                       integer(IN) :: blockList(blockCount), 
!!                       real(IN)    ::  dt  )    
!!
!! DESCRIPTION
!!
!!  Apply thermal heat exchange among temperature components to all
!!  blocks in the specified list. This implementation uses the Spitzer
!!  ion/electron equilibration time.
!!
!! ARGUMENTS
!!
!!   blockCount -- dimension of blockList
!!   blockList -- array of blocks where componenets should exchange
!!                heat
!!   dt  --       passed to the internal hx_burner module  
!!
!! PARAMETERS
!!
!!  useHeatexchange -- Boolean, True.  Turns on Heatexchange unit
!!
!! NOTES
!!
!!  You can define RECOVER_FROM_NEG_TION (see below) to force recovery
!!  from negative ion temperature anomalies such as occur with electron
!!  entropy advection in some situations.
!!***
subroutine Heatexchange ( blockCount, blockList, dt )
  
  use Grid_interface, ONLY: Grid_getBlkIndexLimits, Grid_getBlkPtr, &
       Grid_releaseBlkPtr, Grid_notifySolnDataUpdate
  use Eos_interface, ONLY: Eos_wrapped, Eos, Eos_getAbarZbar
  use Timers_interface, ONLY : Timers_start, Timers_stop  
  use Heatexchange_data,ONLY: hx_useHeatexchange, hx_ieTimeCoef, hx_navo, &
                              hx_logLevel
  use Driver_interface, ONLY: Driver_abortFlash
  use Logfile_interface, ONLY: Logfile_stampMessage

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Eos.h"
#include "Heatexchange.h"

  ! Arguments:
  real,    INTENT(in) :: dt
  integer, INTENT(in) :: blockCount
  integer, INTENT(in), DIMENSION(blockCount) :: blockList

  ! Local Variables:
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer :: i, j, k, n
  integer :: blockID, thisBlock
  real :: rho ! Mass density (g/cc)
  real :: tele, tele_new ! Electron temperature (K)
  real :: tion, tion_new ! Ion temperature (K)
  real :: cvele ! Electron specfic heat (ergs/g/K)
  real :: cvion ! Ion specfic heat (ergs/g/K)
  real :: eqtime ! Ion/Electron equilibration time (s)
  real :: erate ! Rate at which specific energy is transfered (ergs/g/s)
  real :: cvratio ! Ratio of specific heats
  real, pointer, dimension(:,:,:,:) :: solnData
  real :: A, B

  ! EOS data to compute cvele:
  logical, dimension(EOS_VARS+1:EOS_NUM) :: mask
  real, dimension(EOS_NUM) :: eos_arr
  real :: massfrac(NSPECIES)
  real :: abar, zbar

  real :: eele_old
  real :: eele_new
  real :: eion_old
  real :: eion_new
  real :: ratio

  integer, parameter :: maxWarningsDefault = 3
  integer, save      :: maxWarnings = maxWarningsDefault
  integer, save      :: warningsIssued = 0
  logical            :: issueWarning

  character(len=MAX_STRING_LENGTH) :: errmsg

!=========================================================================

  ! Check useHeatexchange flag
  if (.not. hx_useHeatexchange) return

  call Grid_notifySolnDataUpdate(&
       (/EION_VAR,EELE_VAR,TION_VAR,TELE_VAR,ENER_VAR,EINT_VAR,&
         PION_VAR,PELE_VAR,PRES_VAR/) )

  ! START TIMERS
  call Timers_start("heatXchg")

  ! BEGIN LOOP OVER BLOCKS PASSED IN
  do thisBlock = 1, blockCount
     blockID = blockList(thisBlock)

     !GET DIMENSION AND COORD POSITIONS
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
     
     !GET POINTER TO SOLUTION DATA
     call Grid_getBlkPtr(blockID,solnData)

#ifdef RECOVER_FROM_NEG_TION
#define SAFE_TION tion
#else
#define SAFE_TION max(0.0,tion)
#endif
     
     !LOOP OVER CURRENT BLOCK ZONES
     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

              rho = solnData(DENS_VAR, i, j, k)
              tion = solnData(TION_VAR,i,j,k)
              tele = solnData(TELE_VAR,i,j,k)

              ! Compute abar and zbar:
              call Eos_getAbarZbar(solnData(:,i,j,k),abar=abar,zbar=zbar)

              eos_arr(EOS_DENS) = solnData(DENS_VAR, i,j,k)
              eos_arr(EOS_TEMPION) = SAFE_TION
              eos_arr(EOS_TEMPELE) = tele
              eos_arr(EOS_TEMPRAD) = solnData(TRAD_VAR,i,j,k)
              eos_arr(EOS_TEMP)    = tele
              eos_arr(EOS_ABAR) = abar
              eos_arr(EOS_ZBAR) = zbar

              ! Get the electron specific heat:
              mask = .false.
              mask(EOS_CVELE) = .true.
              mask(EOS_CVION) = .true.
              mask(EOS_DET)   = .true.
              massfrac = solnData(SPECIES_BEGIN:SPECIES_BEGIN+NSPECIES-1,i,j,k)
              call Eos(MODE_DENS_TEMP_GATHER,1,eos_arr,massfrac,mask)
              cvele = max(0.0,eos_arr(EOS_CVELE))
              cvion = eos_arr(EOS_CVION)
              cvratio = cvele/cvion

#ifdef CVIO_VAR
              solnData(CVIO_VAR,i,j,k) = cvion
#endif

#ifdef CVEL_VAR
              solnData(CVEL_VAR,i,j,k) = cvele
#endif              

              if(cvion < 0.0 .or. cvele < 0.0) then
                 call Logfile_stampMessage('[Heatexchange] Negative specific heat detected:')

                 write(errmsg, '(a,1pe13.6)') '[Heatexchange] tion (K): ', tion 
                 call Logfile_stampMessage(errmsg)

                 write(errmsg, '(a,1pe13.6)') '[Heatexchange] tele (K): ', tele
                 call Logfile_stampMessage(errmsg)

                 write(errmsg, '(a,1pe13.6)') '[Heatexchange] dens (g/cc): ', solnData(DENS_VAR,i,j,k)
                 call Logfile_stampMessage(errmsg)

                 write(errmsg, '(a,1pe13.6)') '[Heatexchange] cvion (ergs/K/g): ', cvion
                 call Logfile_stampMessage(errmsg)

                 write(errmsg, '(a,1pe13.6)') '[Heatexchange] cvele (ergs/K/g): ', cvele
                 call Logfile_stampMessage(errmsg)

                 call Logfile_stampMessage("[Heatexchange] Species:")
                 do n = SPECIES_BEGIN, SPECIES_BEGIN+NSPECIES-1
                    write(errmsg, '(a,1pe15.6)') '[Heatexchange] ', solnData(n,i,j,k)
                    call Logfile_stampMessage(errmsg)
                 end do

                 call Driver_abortFlash('[Heatexchange] Negative specific heat. See log file')
              end if
                            
              call hx_ieEquilTime(&
                   zbar, &
                   abar, &
                   tele, SAFE_TION, &
                   rho*hx_navo/abar, &
                   eqtime)

              eqtime = eqtime * hx_ieTimeCoef

              A = (tion + cvratio * tele)/(1 + cvratio)
              B = (tele - tion)/(1 + cvratio)
#ifdef RECOVER_FROM_NEG_TION
              if (tion < 0.0) then
                 B = A / cvratio
              end if
#endif

              tele_new = A + B*exp(-(1+cvratio)*dt/eqtime)
              tion_new = A - cvratio*B*exp(-(1+cvratio)*dt/eqtime)
              
              ! Use the temperature change to compute new internal
              ! energies. This will ensure energy conservation.
              solnData(EION_VAR,i,j,k) = solnData(EION_VAR,i,j,k) + cvion*(tion_new-tion)
              solnData(EELE_VAR,i,j,k) = solnData(EELE_VAR,i,j,k) + cvele*(tele_new-tele)

              if ( solnData(EELE_VAR,i,j,k) < 0.0 .or. &
                   solnData(EION_VAR,i,j,k) < 0.0 ) then

                 issueWarning = .FALSE.
                 if (hx_logLevel .GE. HX_LOGLEVEL_INFO_FALLBACK) then
                    if (hx_logLevel .GE. HX_LOGLEVEL_INFO_ALL) maxWarnings = -1
                    if (maxWarnings < 0) then
                       issueWarning = .TRUE.
                    else
                       if (warningsIssued < maxWarnings) issueWarning = .TRUE.
                       if (warningsIssued == maxWarnings) then
                          warningsIssued = warningsIssued + 1
                          if (hx_logLevel > HX_LOGLEVEL_INFO_FALLBACK) &
                               call Logfile_stampMessage("[Heatexchange] Note: Further warnings suppressed", .true.)
                       end if
                    end if
                 end if

                 if ( issueWarning .AND. hx_logLevel > HX_LOGLEVEL_INFO_FALLBACK+1) then
                    call Logfile_stampMessage("[Heatexchange] Warning: Negative internal energy detected", .true.)
                    call Logfile_stampMessage("[Heatexchange] Correcting for this...", .true.)
                 end if

                 eele_old = solnData(EELE_VAR,i,j,k) - cvele*(tele_new-tele)
                 eion_old = solnData(EION_VAR,i,j,k) - cvion*(tion_new-tion)
                 eele_new = solnData(EELE_VAR,i,j,k)
                 eion_new = solnData(EION_VAR,i,j,k)
                 
                 if (issueWarning) then
                    write(errmsg, '(a,2(1pe13.6))') '[Heatexchange] OLD EION/TION: ', eion_old, tion
                    call Logfile_stampMessage(errmsg, .true.)

                    write(errmsg, '(a,2(1pe13.6))') '[Heatexchange] OLD EELE/TELE: ', eele_old, tele
                    call Logfile_stampMessage(errmsg, .true.)

                    write(errmsg, '(a,2(1pe13.6))') '[Heatexchange] BAD EION/TION: ', eion_new, tion_new
                    call Logfile_stampMessage(errmsg, .true.)

                    write(errmsg, '(a,2(1pe13.6))') '[Heatexchange] BAD EELE/TELE: ', eele_new, tele_new
                    call Logfile_stampMessage(errmsg, .true.)
                 end if

                 ! Negative internal energy. This can happen because
                 ! the specific heat can jump a lot (it may in fact be
                 ! discontinuous). In this case, the temperatures must
                 ! be used to update the internal energies...

                 eos_arr = 0.0
                 eos_arr(EOS_DENS) = solnData(DENS_VAR, i,j,k)
                 eos_arr(EOS_TEMPION) = tion_new
                 eos_arr(EOS_TEMPELE) = tele_new
                 eos_arr(EOS_TEMPRAD) = solnData(TRAD_VAR,i,j,k)
                 eos_arr(EOS_TEMP)    = tele

                 mask = .false.
                 mask(EOS_EINTELE) = .true.
                 mask(EOS_EINTION) = .true.
                 massfrac = solnData(SPECIES_BEGIN:SPECIES_BEGIN+NSPECIES-1,i,j,k)
                 call Eos(MODE_DENS_TEMP_GATHER,1,eos_arr,massfrac,mask)
                 eele_new = eos_arr(EOS_EINTELE)
                 eion_new = eos_arr(EOS_EINTION)

                 if (issueWarning) then
                    write(errmsg, '(a,2(1pe13.6))') '[Heatexchange] CORRECTED EION: ', eion_new
                    call Logfile_stampMessage(errmsg, .true.)

                    write(errmsg, '(a,2(1pe13.6))') '[Heatexchange] CORRECTED EELE: ', eele_new
                    call Logfile_stampMessage(errmsg, .true.)
                 end if

                 ! Now we have to scale the ion/electron internal
                 ! energy to ensure that we have energy
                 ! conservation...
                 
                 ratio = (eele_old + eion_old) / (eele_new + eion_new)
                 solnData(EION_VAR,i,j,k) = eion_new * ratio
                 solnData(EELE_VAR,i,j,k) = eele_new * ratio

                 if (issueWarning) then
                    write(errmsg, '(a,2(1pe13.6))') '[Heatexchange] FINAL EION: ', solnData(EION_VAR,i,j,k)
                    call Logfile_stampMessage(errmsg, .true.)

                    write(errmsg, '(a,2(1pe13.6))') '[Heatexchange] FINAL EELE: ', solnData(EELE_VAR,i,j,k)
                    call Logfile_stampMessage(errmsg, .true.)
                    call Logfile_stampMessage("", .true.)

                    warningsIssued = warningsIssued + 1
                 end if
              end if
           enddo
           
        enddo
     end do

     ! RELEASE MEMORY/POINTERS
     call Grid_releaseBlkPtr(blockID,solnData)
     
  end do

  do thisBlock = 1, blockCount
     blockID = blockList(thisBlock)

     !GET DIMENSION AND COORD POSITIONS
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

     ! MAKE HYDRO CONSISTENT WITH UPDATED INTERNAL ENERGY
     call Eos_wrapped(MODE_DENS_EI_GATHER,blkLimits,blockID)
     
  end do
  
  call Timers_stop("heatXchg")  
  return
  
end subroutine Heatexchange
