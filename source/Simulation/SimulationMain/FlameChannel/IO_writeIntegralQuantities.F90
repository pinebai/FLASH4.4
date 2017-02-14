!!****if* source/Simulation/SimulationMain/FlameChannel/IO_writeIntegralQuantities
!!
!! NAME
!!
!!  IO_writeIntegralQuantities
!!
!! SYNOPSIS
!!
!!  call IO_writeIntegralQuantities(integer(in) :: isfirst,
!!                                  real(in) :: simtime)
!!
!! DESCRIPTION
!!
!!   Compute the values of integral quantities (eg. total energy)
!!   and write them to an ASCII file.  If this is the initial step,
!!   create the file and write a header to it before writing the data.
!!
!!
!! ARGUMENTS
!!
!!   isFirst - if 1 then write header info plus data, otherwise just write data
!!
!!   simTime - simulation time
!!
!!
!!***

! See source/IO/IOMain/IO_writeIntegralQuantities.F90
! for more on API and original (example) subroutine
!
! This version has some additional metrics specific to the RT Flame
! in a channel (burning rate, surface area)
!
! Dean Townsley 2008


!!REORDER(4):solnData

!! So tired of trying to figure out which index of lsum,gsum is what
#define IO_DENS       1
#define IO_MOMX       2
#define IO_MOMY       3
#define IO_MOMZ       4
#define IO_ENER       5
#define IO_KE         6
#define IO_TURB_KE    7
#define IO_EINT       8
#define IO_MASS_ASH   9
#define IO_FLAM_DENS  10
#define IO_FLAM_TURB  11
#define IO_FLAM_VOL   12
#define IO_OUTF_MOMX  13
#define IO_OUTF_VOL   14
#define IO_FLAM_SURF1 15
#define IO_FLAM_SURF2 16
#define IO_FLAM_SURF3 17

subroutine IO_writeIntegralQuantities (isFirst, simTime)

  use IO_data, ONLY : io_restart, io_statsFileName, io_globalMe, io_globalComm
  use Grid_interface, ONLY : Grid_getListOfBlocks, &
    Grid_getBlkIndexLimits, Grid_getBlkPtr, Grid_getSingleCellVol, &
    Grid_releaseBlkPtr, Grid_fillGuardCells, Grid_getBlkBC
  use IO_interface, ONLY : IO_setScalar
  use Driver_data, ONLY : dr_dtOld
  use fl_effData, ONLY : fl_effDeltae, fl_eff_sumy_u, fl_eff_sumy_b, &
                         fl_eff_ye_u, fl_eff_ye_b
  use Eos_interface, ONLY : Eos
  use Eos_data, ONLY : eos_maxNewton, eos_tol
  use Simulation_data, ONLY : sim_rhoAmbient, sim_last_burned_mass, &
                              sim_ymin, sim_ymax, sim_zmin, sim_zmax, &
                              sim_inflowVx, sim_flamespeed, sim_variableInflow, &
                              sim_crossArea
  use ut_contourSurfaceInterface, ONLY: ut_contourSurfaceAreaBlock

  implicit none

#include "Flash_mpi.h"
#include "constants.h"
#include "Flash.h"
#include "Eos.h"
  
  real, intent(in) :: simTime

  integer, intent(in) :: isFirst

  integer :: lb, count
  
  integer :: funit = 99
  integer :: error
  
  character (len=MAX_STRING_LENGTH), save :: fname 
  
  integer, dimension(MAXBLOCKS) :: blockList

  integer, dimension(HIGH,MDIM) :: blkLimits, blkLimitsGC, faces

  integer, parameter ::  nGlobalSum = 17  ! Number of globally-summed quantities
  real, dimension(nGlobalSum) :: lsum, gsum ! Local and Global summed quantities

  integer :: i, j, k
  real :: dvol             !, del(MDIM)
  real, DIMENSION(:,:,:,:), POINTER :: solnData

  ! setting flam_max = 1.e-2 allows the density to change too
  ! much for an accurate estimate of the burning rate
  ! however, 1.e-3 is sometimes too small to actually sample
  ! any volume
  real, parameter :: flam_min = 1.e-6, flam_max = 1.e-2
  integer, parameter :: nlevels = 3
  real, dimension(nlevels) :: isolevels, blkAreas
  integer, dimension(MDIM) :: point
  real :: dt, brate, fspd, inflowVx, momx

  logical, dimension(NUNK_VARS) :: gcMask
  logical :: isOutflow

  real, dimension(EOS_NUM) :: eosData
  real :: p_target, rho_guess, eint_target, q_flam, flam, enthalpy
  real :: abar, zbar
  real :: err1, err2
  integer :: iter1, iter2
  logical, dimension(EOS_VARS+1:EOS_NUM) :: eosMask

  isolevels(1) = 0.1
  isolevels(2) = 0.5
  isolevels(3) = 0.9

  eosMask(:) = .false.
  eosMask(EOS_DPD) = .true.

#if NDIM > 1
  gcMask(:) = .false.
  gcMask(FLAM_MSCALAR) = .true.

  call Grid_fillGuardCells(CENTER, ALLDIR,&
       maskSize=NUNK_VARS, mask=gcMask,makeMaskConsistent=.true.)
#endif
  

  ! Sum quantities over all locally held leaf-node blocks.
  gsum  = 0.
  lsum = 0.
  
  call Grid_getListOfBlocks(LEAF, blockList, count)
  
  do lb = 1, count
     !get the index limits of the block
     call Grid_getBlkIndexLimits(blockList(lb), blkLimits, blkLimitsGC)

     ! get a pointer to the current block of data
     call Grid_getBlkPtr(blockList(lb), solnData)

     ! see if block is on an outflow boundary
     call Grid_getBlkBC(blockList(lb), faces)
     isOutFlow = (faces(LOW,IAXIS) == USER_DEFINED).or. &
                 (faces(LOW,IAXIS) == OUTFLOW)

     ! Sum contributions from the indicated blkLimits of cells.
     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
              
              point(IAXIS) = i
              point(JAXIS) = j
              point(KAXIS) = k

!! Get the cell volume for a single cell
              call Grid_getSingleCellVol(blockList(lb), EXTERIOR, point, dvol)
     
              ! mass   
#ifdef DENS_VAR
              lsum(IO_DENS) = lsum(IO_DENS) + solnData(DENS_VAR,i,j,k)*dvol 
#endif           


#ifdef DENS_VAR
#ifdef VELX_VAR      
              ! momentum
              momx = solnData(DENS_VAR,i,j,k) * solnData(VELX_VAR,i,j,k) * dvol
              lsum(IO_MOMX) = lsum(IO_MOMX) + momx 

#endif
#ifdef VELY_VAR      

              lsum(IO_MOMY) = lsum(IO_MOMY) + solnData(DENS_VAR,i,j,k) * & 
                   &                                solnData(VELY_VAR,i,j,k)*dvol
           
#endif
#ifdef VELZ_VAR      
              lsum(IO_MOMZ) = lsum(IO_MOMZ) + solnData(DENS_VAR,i,j,k) * & 
                   &                                solnData(VELZ_VAR,i,j,k)*dvol
#endif

              ! total energy
#ifdef ENER_VAR
              lsum(IO_ENER) = lsum(IO_ENER) + solnData(ENER_VAR,i,j,k) * & 
                   &                                solnData(DENS_VAR,i,j,k)*dvol
#endif
           
#ifdef VELX_VAR      
#ifdef VELY_VAR      
#ifdef VELZ_VAR      
              ! kinetic energy
              lsum(IO_KE) = lsum(IO_KE) + 0.5*solnData(DENS_VAR,i,j,k) * & 
                   &                             (solnData(VELX_VAR,i,j,k)**2+ & 
                   &                              solnData(VELY_VAR,i,j,k)**2+ & 
                   &                              solnData(VELZ_VAR,i,j,k)**2)*dvol           

#endif
#endif
#endif

#ifdef TURB_VAR
              ! turbulent kinetic energy
              lsum(IO_TURB_KE) = lsum(IO_TURB_KE) + 0.5*solnData(DENS_VAR,i,j,k) * &
                   &                  solnData(TURB_VAR,i,j,k)**2 * dvol
#endif

#ifdef EINT_VAR
              ! internal energy
              lsum(IO_EINT) = lsum(IO_EINT) + solnData(DENS_VAR,i,j,k) * & 
                   &                                solnData(EINT_VAR,i,j,k)*dvol
#endif

              flam = solnData(FLAM_MSCALAR,i,j,k)
              lsum(IO_MASS_ASH) = lsum(IO_MASS_ASH) + flam * &
                   &                      solnData(DENS_VAR,i,j,k)*dvol

              if (( flam > flam_min ) .and. ( flam < flam_max )) then

                 !! find average density just ahead of flame
                 !lsum(IO_FLAM_DENS) = lsum(IO_FLAM_DENS) + solnData(DENS_VAR,i,j,k)*dvol
                 !lsum(IO_FLAM_VOL) = lsum(IO_FLAM_VOL) + dvol

                 ! estimate the unburned density in low-Mach limit
                 ! eint + p/rho - q flam = const, p = const
                 p_target = solnData(PRES_VAR,i,j,k)
                 q_flam = fl_effDeltae * flam
                 abar = 1.0 / ( (1.0-flam)*fl_eff_sumy_u + flam*fl_eff_sumy_b )
                 zbar = abar * ( (1.0-flam)*fl_eff_ye_u + flam*fl_eff_ye_b )
                 enthalpy = solnData(EINT_VAR,i,j,k) - q_flam +  &
                    solnData(PRES_VAR,i,j,k) / solnData(DENS_VAR,i,j,k)
                 eint_target = solnData(EINT_VAR,i,j,k) - q_flam
                 rho_guess = p_target / (enthalpy - eint_target)
                 err1 = abs(eint_target - solnData(EINT_VAR,i,j,k))/eint_target
!                 err1 = abs(rho_guess - solnData(DENS_VAR,i,j,k))/rho_guess

                 iter1 = 0
                 do while (err1 > eos_tol)
                    iter1 = iter1 + 1

                    do iter2 = 1, eos_maxNewton

                       eosData(EOS_EINT) = eint_target
                       eosData(EOS_DENS) = rho_guess
                       eosData(EOS_ABAR) = abar
                       eosData(EOS_ZBAR) = zbar
                       call Eos(MODE_DENS_EI,1, eosData, mask=eosMask)

                       rho_guess = min(10.0*eosData(EOS_DENS), &
                                   max(0.1*eosData(EOS_DENS), &
                         eosData(EOS_DENS) + (p_target - eosData(EOS_PRES)) / &
                         eosData(EOS_DPD)))

                       err2 = abs(rho_guess - eosData(EOS_DENS)) / rho_guess

                       if (err2 <= eos_tol) exit
                    enddo

                    if (err2 > eos_tol) then
                       err1 = err2
                       exit
                    endif

                    eint_target = enthalpy - p_target / rho_guess
                    err1 = abs(eint_target - eosData(EOS_EINT)) / eint_target

!                    eosData(EOS_PRES) = p_target
!                    eosData(EOS_DENS) = rho_guess
!                    eosData(EOS_ABAR) = abar
!                    eosData(EOS_ZBAR) = zbar
!                    call Eos(MODE_DENS_PRES,1, eosData)
!
!                    rho_guess = p_target / ( enthalpy - eosData(EOS_EINT) )
!
!                    err1 = abs(rho_guess - eosData(EOS_DENS)) / rho_guess

                    if (iter1 >= eos_maxNewton) exit

                 enddo
                 
                 if (err1 > eos_tol) &
                    call Driver_abortFlash("[IO_writeIntegralQuantities] Unburned density estimate did not converge")

                 lsum(IO_FLAM_DENS) = lsum(IO_FLAM_DENS) + rho_guess*dvol
                 lsum(IO_FLAM_VOL) = lsum(IO_FLAM_VOL) + dvol
#ifdef TURB_VAR
                 lsum(IO_FLAM_TURB) = lsum(IO_FLAM_TURB) + solnData(TURB_VAR,i,j,k)*dvol
#endif

              endif

#ifdef VELX_VAR
              if (isOutflow .and. (i == blkLimits(LOW,IAXIS))) then
                 ! outflow conserves mass, so we can sum up the last cell on
                 ! grid since it is equal to the first guard cell along the
                 ! boundary
                 lsum(IO_OUTF_MOMX) = lsum(IO_OUTF_MOMX) + momx
                 lsum(IO_OUTF_VOL) = lsum(IO_OUTF_VOL) + dvol
              endif
#endif


#endif ! ifdef DENS_VAR

           enddo
        enddo
     enddo

#if NDIM > 1
     call ut_contourSurfaceAreaBlock(nlevels,isolevels,solnData(FLAM_MSCALAR,:,:,:), &
                                     blkLimits,blockList(lb),blkAreas)
     lsum(IO_FLAM_SURF1) = lsum(IO_FLAM_SURF1) + blkAreas(1)
     lsum(IO_FLAM_SURF2) = lsum(IO_FLAM_SURF2) + blkAreas(2)
     lsum(IO_FLAM_SURF3) = lsum(IO_FLAM_SURF3) + blkAreas(3)
#endif

     call Grid_releaseBlkPtr(blockList(lb), solnData)

  enddo
 
  
  ! Now the MASTER_PE sums the local contributions from all of
  ! the processors and writes the total to a file.
  
  call MPI_Reduce (lsum, gsum, nGlobalSum, FLASH_REAL, MPI_SUM, & 
       &                MASTER_PE, io_globalComm, error)

  dt = dr_dtOld
  if (io_globalMe == MASTER_PE) then ! gsum only valid on master node

     if (gsum(IO_OUTF_VOL) > 0.0) then
        gsum(IO_OUTF_MOMX) = gsum(IO_OUTF_MOMX)/gsum(IO_OUTF_VOL)
     else
        gsum(IO_OUTF_MOMX) = 0.0
     endif

     if (sim_last_burned_mass==-1.0) then
        brate = sim_flamespeed * sim_rhoAmbient * sim_crossArea
        ! for first time step, set to steady state solution
        ! fspd = brate / gsum(IO_FLAM_DENS) / gsum(IO_FLAM_SURF1)
     else ! defined average mass flow in to be positive 
        brate = (gsum(IO_MASS_ASH)-sim_last_burned_mass)/(2*dt) - &
                gsum(IO_OUTF_MOMX) * sim_crossArea 
!        brate = sim_flamespeed * sim_rhoAmbient * sim_crossArea - &
!                (gsum(9)-sim_last_burned_mass)/(2*dt)
     endif

     ! average density just ahead of flame
     if (gsum(IO_FLAM_VOL) > 0.0) then
        gsum(IO_FLAM_DENS) = gsum(IO_FLAM_DENS)/gsum(IO_FLAM_VOL)
        gsum(IO_FLAM_TURB) = gsum(IO_FLAM_TURB)/gsum(IO_FLAM_VOL)
     else ! let's just use what it should be
        gsum(IO_FLAM_DENS) = sim_rhoAmbient
     endif

     ! calculate flame front propagation speed: s = mb / rho / BurnArea
     fspd = brate / gsum(IO_FLAM_DENS) / gsum(IO_FLAM_SURF1)

     ! update inflow rate
     if (sim_variableInflow) then
        inflowVx = -brate / sim_rhoAmbient / sim_crossArea
        ! maximum change of 10% of current value (inflowVx < 0)
        inflowVx = min(0.9*sim_inflowVx, max(1.1*sim_inflowVx, inflowVx))
     else
        inflowVx = sim_inflowVx
     endif

     sim_last_burned_mass = gsum(IO_MASS_ASH)


     ! create the file from scratch if it is a not a restart simulation, 
     ! otherwise append to the end of the file
     if (isfirst == 0) then
        open (funit, file=trim(io_statsFileName), position='APPEND')
     else 
        if (.NOT. io_restart) then
           open (funit, file=trim(io_statsFileName)) 
           write (funit, 10)               &
                '#time                     ', &
                'mass                      ', &
                'x-momentum                ', &
                'y-momentum                ', & 
                'z-momentum                ', &
                'E_total                   ', &
                'E_kinetic                 ', &
                'E_turbulent               ', &
                'E_internal                ', &
                'Burned Mass               ', &
                'dens_burning_ave          ', &
                'turb_burning_ave          ', &
                'db_ave samplevol          ', &
                'inflow-momentum_ave       ', &
                'inflow samplevol          ', &
                'Burning rate              ', &
                'effective front speed     ', &
                'surface area flam=0.1     ', &
                'surface area flam=0.5     ', &
                'surface area flam=0.9     ', &
                'updated inflow rate       '

10         format (2x,50(a25, :, 1X))

        else
           open (funit, file=trim(io_statsFileName), position='APPEND')
           write (funit, 11) 
11         format('# simulation restarted')
        endif
     endif

     ! Write the global sums to the file.
     write (funit, 12) simtime, gsum(IO_DENS:IO_OUTF_VOL), brate, fspd, &
        gsum(IO_FLAM_SURF1:IO_FLAM_SURF3), inflowVx
12   format (1x, 50(es25.18, :, 1x))
 
     close (funit)          ! Close the file.
     
  else ! not MASTER_PE

     sim_last_burned_mass = 0.0

  endif

  if (sim_variableInflow) then

     call MPI_Bcast (inflowVx,1,FLASH_REAL,MASTER_PE,io_globalComm,error)
     sim_inflowVx = inflowVx

  endif

  call MPI_Barrier (io_globalComm, error)
  
  !=============================================================================
  
  return
end subroutine IO_writeIntegralQuantities



