!!****if* source/physics/sourceTerms/Heatexchange/HeatexchangeMain/ConstCoulomb/Heatexchange
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
!!  Apply thermal heat exchange among temperature components
!!  to all blocks in the specified list.
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
!!
!!***


subroutine Heatexchange ( blockCount, blockList, dt )
  
  use Grid_interface, ONLY  : Grid_getBlkIndexLimits, Grid_getBlkPtr, &
       Grid_releaseBlkPtr, Grid_notifySolnDataUpdate
  use Eos_interface, ONLY   : Eos_wrapped
  use Timers_interface, ONLY : Timers_start, Timers_stop
  
  use Heatexchange_data,ONLY: hx_useHeatexchange, hx_coulombLog, hx_c13, hx_c23, &
       hx_singleSpeciesA, hx_singleSpeciesZ, &
       hx_kBoltzmann, hx_Avogadro, hx_eMassInUAmu, hx_eleCharge

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Eos.h"

  !args
  integer, INTENT(in)                        :: blockCount
  integer, INTENT(in), DIMENSION(blockCount)  :: blockList
  real,    INTENT(in)                        :: dt

  ! locals
  integer                    :: i, j, k
  integer                    :: blockID, thisBlock
  real                       :: temp, dens, eint, pres
  real                       :: temp1,temp2,temp3, eint1,eint2,eint3
  real                       :: t12diff, t13diff, t23diff
  real                       :: bsEint1,bsEint2,bsEint3
  real                       :: Ye
  real                       :: flame
  real                       :: qdot, q1dot, q2dot, q3dot, qbar
  real                       :: dynamicZ, relA, ni, nePerDens, ionMassInUAmu, memi
  real                       :: numFactor, chargeEtcFactor, massFactor, temperFactor
  real                       :: nuj, lamx
  logical                    :: changedZone

  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  logical :: getGuardCells = .true.
  real, pointer, dimension(:,:,:,:)                    :: solnData
#ifndef EINT_VAR
  real :: energyKinetic
#endif




  ! Check useHeatexchange flag
  if (.not. hx_useHeatexchange) return

  call Grid_notifySolnDataUpdate(CENTER)

  ! START TIMERS
  call Timers_start("heatXchg")
  



  ! BEGIN LOOP OVER BLOCKS PASSED IN
  do thisBlock = 1, blockCount
     
     blockID = blockList(thisBlock)
     changedZone = .FALSE.
     
     !GET DIMENSION AND COORD POSITIONS
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
     
     !GET POINTER TO SOLUTION DATA
     call Grid_getBlkPtr(blockID,solnData)
     
     !LOOP OVER CURRENT BLOCK ZONES
     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

              dens    = solnData(DENS_VAR,i,j,k)
              temp1    = solnData(TION_VAR,i,j,k)
              temp2    = solnData(TELE_VAR,i,j,k)
              temp3    = solnData(TRAD_VAR,i,j,k)
              t12diff  = temp1 - temp2
              t13diff  = temp1 - temp3
              t23diff  = temp2 - temp3

#ifdef EINT_VAR
              eint    = solnData(EINT_VAR,i,j,k)
#else
              energyKinetic = solnData(VELX_VAR,i,j,k)**2
#if NDIM >= 2
              energyKinetic = energyKinetic + solnData(VELY_VAR,i,j,k)**2
#endif
#if NDIM > 2
              energyKinetic = energyKinetic + solnData(VELZ_VAR,i,j,k)**2
#endif
              eint    = solnData(ENER_VAR,i,j,k) - 0.5*energyKinetic
#endif

              eint1   = solnData(EION_VAR,i,j,k)
              eint2   = solnData(EELE_VAR,i,j,k)
              eint3   = solnData(ERAD_VAR,i,j,k)

              pres    = solnData(PRES_VAR,i,j,k)

              q1dot   = 0.e0
              q2dot   = 0.e0
              q3dot   = 0.e0

              bsEint1 = 1; bsEint2 = 1; bsEint3 = 1
              ! For the following see Imamura, Durisen, Lamb, Weast, ApJ 1987, (A5) and (A6) 

              dynamicZ = hx_singleSpeciesZ
!!$              Zp = dynamicZ + 1
!!$              Zinv = 1.0/dynamicZ; ZpInv = 1.0/Zp
              relA = hx_singleSpeciesA
              Ye = dynamicZ / relA
              ni = dens * hx_Avogadro / relA
              nePerDens = hx_Avogadro * Ye
              numFactor = 8*sqrt(2*PI)/3.0
              chargeEtcFactor = ni * dynamicZ**2 * hx_eleCharge**4 * hx_kBoltzmann**(-1.5)
              ionMassInUAmu = relA - dynamicZ * hx_eMassInUAmu
              memi = hx_eMassInUAmu / ionMassInUAmu
              massFactor = sqrt(hx_eMassInUAmu)/ionMassInUAmu * sqrt(hx_Avogadro) !DEV: not right for MKS units!

              temperFactor = (temp2 + memi * temp1)**(-1.5)

              nuj = numFactor * chargeEtcFactor * hx_coulombLog * massFactor * temperFactor
              lamx = 1.5 * nePerDens * hx_kBoltzmann * t12diff * nuj


              q3dot =   bsEint1 * hx_c13 * t13diff  +  bsEint2 * hx_c23 * t23diff
              q2dot =   bsEint1 * lamx              -  bsEint3 * hx_c23 * t23diff
              q1dot = -(bsEint2 * lamx              +  bsEint3 * hx_c13 * t13diff)
              qdot = q1dot + q2dot + q3dot

              changedZone = .TRUE.

                 ! UPDATE SOLUTION DATA

!!$              solnData(ENUC_VAR,i,j,k)     = qdot 
              solnData(EION_VAR,i,j,k)     = solnData(EION_VAR,i,j,k) + q1dot*dt
              solnData(EELE_VAR,i,j,k)     = solnData(EELE_VAR,i,j,k) + q2dot*dt
              solnData(ERAD_VAR,i,j,k)     = solnData(ERAD_VAR,i,j,k) + q3dot*dt

!!$              solnData(ENER_VAR,i,j,k)     = solnData(ENER_VAR,i,j,k) + qdot*dt
!!$              solnData(EINT_VAR,i,j,k)     = solnData(EINT_VAR,i,j,k) + qdot*dt
           enddo
           
        enddo
     enddo


     ! MAKE HYDRO CONSISTENT WITH UPDATED INTERNAL ENERGY
     if (changedZone) then
        call Eos_wrapped(MODE_DENS_EI_GATHER,blkLimits,blockID)
     end if
     
     ! RELEASE MEMORY/POINTERS
     call Grid_releaseBlkPtr(blockID,solnData)
  end do
  
  call Timers_stop("heatXchg")
  
  return
  
end subroutine Heatexchange
