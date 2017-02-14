!!****if* source/Simulation/SimulationMain/ShafranovShock/Heatexchange
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
!!  Apply thermal head exchange among temperature components
!!  to all blocks in specified list.
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
       Grid_releaseBlkPtr
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
  logical                    :: changedZone, changedBlock

  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC, eosLimits
  logical :: getGuardCells = .true.
  real, pointer, dimension(:,:,:,:)                    :: solnData
#ifndef EINT_VAR
  real :: energyKinetic
#endif

#define M_ELECTRON 9.1093897E-28
#define M_NUCLEON 1.6726231E-24

#define SHAFRANOV

  ! Check useHeatexchange flag
  if (.not. hx_useHeatexchange) return

  ! START TIMERS
  call Timers_start("heatXchg")
  



  ! BEGIN LOOP OVER BLOCKS PASSED IN
  do thisBlock = 1, blockCount
     
     blockID = blockList(thisBlock)
     changedBlock = .FALSE.
     
     !GET DIMENSION AND COORD POSITIONS
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
     eosLimits(LOW,IAXIS) = blkLimits(HIGH,IAXIS) + 1
     eosLimits(HIGH,IAXIS) = blkLimits(LOW,IAXIS) - 1
#if NDIM >= 2
     eosLimits(LOW,JAXIS) = blkLimits(HIGH,JAXIS) + 1
     eosLimits(HIGH,JAXIS) = blkLimits(LOW,JAXIS) - 1
#else
     eosLimits(:,JAXIS) = blkLimits(:,JAXIS)
#endif
#if NDIM > 2
     eosLimits(LOW,KAXIS) = blkLimits(HIGH,KAXIS) + 1
     eosLimits(HIGH,KAXIS) = blkLimits(LOW,KAXIS) - 1
#else
     eosLimits(:,KAXIS) = blkLimits(:,KAXIS)
#endif

     !GET POINTER TO SOLUTION DATA
     call Grid_getBlkPtr(blockID,solnData)
     
     !LOOP OVER CURRENT BLOCK ZONES
     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
              changedZone = .FALSE.

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
              lamx = (1.5 * nePerDens * hx_kBoltzmann * t12diff * nuj)

#ifdef SHAFRANOV
              ni   = dens*hx_Avogadro/(hx_singleSpeciesA)

              ! To reproduce exact Carlo's forumula
              ni   = dens/(hx_singleSpeciesZ*M_ELECTRON + hx_singleSpeciesA*M_NUCLEON)

              ! Electron-Ion relaxation time,  p. 421 of [Zel'Dovich, Y.B. and Raizer, Y.P.(2002)] eq 6.120
              nuj  =  (252.0*hx_singleSpeciesA*(temp2**1.5))/(ni*(hx_singleSpeciesZ**2)*hx_coulombLog)

              ! Heat exchange, p. 509 of [Zel'Dovich, Y.B. and Raizer, Y.P.(2002)] eq 7.31
              lamx = (1.5*hx_kBoltzmann*hx_singleSpeciesZ*ni*t12diff)/nuj  ! lamx is energy exchange per unit time per unit volume.

              lamx = lamx / dens ! We store specific energy in FLASH, energy exchange per unit time per unit mass.

#endif
              q3dot =   bsEint1 * hx_c13 * t13diff  +  bsEint2 * hx_c23 * t23diff
              q2dot =   bsEint1 * lamx              -  bsEint3 * hx_c23 * t23diff
              q1dot = -(bsEint2 * lamx              +  bsEint3 * hx_c13 * t13diff)
              qdot = q1dot + q2dot + q3dot

              if (abs(q1dot*dt) > abs(eint1*1.0e-19)) then
                 changedZone = .TRUE.
              else if (abs(q2dot*dt) > abs(eint2*1.0e-19)) then
                 changedZone = .TRUE.
              else if (abs(q3dot*dt) > abs(eint3*1.0e-19)) then
                 changedZone = .TRUE.
              end if


              if (changedZone) then
                 eosLimits(LOW,IAXIS)  = min(i, eosLimits(LOW,IAXIS))
                 eosLimits(HIGH,IAXIS) = max(i, eosLimits(HIGH,IAXIS))
                 eosLimits(LOW,JAXIS)  = min(j, eosLimits(LOW,JAXIS))
                 eosLimits(HIGH,JAXIS) = max(j, eosLimits(HIGH,JAXIS))
                 eosLimits(LOW,KAXIS)  = min(k, eosLimits(LOW,KAXIS))
                 eosLimits(HIGH,KAXIS) = max(k, eosLimits(HIGH,KAXIS))
                 changedBlock = .TRUE.
              end if

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
     if (changedBlock) then
        call Eos_wrapped(MODE_DENS_EI_GATHER,eosLimits,blockID)
     end if
     
     ! RELEASE MEMORY/POINTERS
     call Grid_releaseBlkPtr(blockID,solnData)
  end do
  
  call Timers_stop("heatXchg")
  
  return
  
end subroutine Heatexchange
