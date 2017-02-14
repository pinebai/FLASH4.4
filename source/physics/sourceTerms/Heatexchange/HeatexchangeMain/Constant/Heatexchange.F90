!!****if* source/physics/sourceTerms/Heatexchange/HeatexchangeMain/Constant/Heatexchange
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
!!   blockList -- array of blocks where components should exchange
!!                internal energy
!!   dt  --       current time step
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
  
  use Heatexchange_data,ONLY: hx_useHeatexchange, hx_c12, hx_c13, hx_c23, &
       hx_meshMe, hx_useRadTrans

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
  real                       :: temp, dens, eint
  real                       :: temp1,temp2,temp3, eint1,eint2,eint3
  real                       :: t12diff, t13diff, t23diff
  real                       :: bsEint1,bsEint2,bsEint3
  real                       :: sumyi, ye
  real                       :: flame
  real                       :: qdot, q1dot, q2dot, q3dot, qbar
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

              q1dot   = 0.e0
              q2dot   = 0.e0
              q3dot   = 0.e0

              bsEint1 = 1; bsEint2 = 1; bsEint3 = 1
              q3dot =   bsEint1 * hx_c13 * t13diff  +  bsEint2 * hx_c23 * t23diff
              q2dot =   bsEint1 * hx_c12 * t12diff  -  bsEint3 * hx_c23 * t23diff
              q1dot = -(bsEint2 * hx_c12 * t12diff  +  bsEint3 * hx_c13 * t13diff)
              qdot = q1dot + q2dot + q3dot

              changedZone = .TRUE.

                 ! UPDATE SOLUTION DATA

!!$              solnData(ENUC_VAR,i,j,k)     = qdot 
              solnData(EION_VAR,i,j,k)     = solnData(EION_VAR,i,j,k) + q1dot*dt
              solnData(EELE_VAR,i,j,k)     = solnData(EELE_VAR,i,j,k) + q2dot*dt

              if(.not. hx_useRadTrans) then
                 solnData(ERAD_VAR,i,j,k)     = solnData(ERAD_VAR,i,j,k) + q3dot*dt
              end if

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
