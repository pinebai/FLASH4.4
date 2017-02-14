!!****if* source/Grid/GridMain/Grid_renormMassScalars
!!
!! NAME
!!
!!  Grid_renormMassScalars
!!
!!
!! SYNOPSIS
!!
!!  Grid_renormMassScalars(integer(IN)  :: blkLimits(2,MDIM),
!!                         real,pointer :: solnData(:,:,:,:))
!!
!! DESCRIPTION
!!
!!  Renormalize the various Mass Scalar's in groups so they sum to 1.
!!
!! ARGUMENTS
!!
!!  blkLimits - the index limits for internal zones of the block to renormalize
!!  solnData -  Pointer to the block to be renormalized
!!
!!***

! solnData depends on the ordering on unk
!!REORDER(4): solnData


#ifdef DEBUG
#define DEBUG_MS
#endif

subroutine Grid_renormMassScalars(blkLimits,solnData)
  use Simulation_interface, ONLY : Simulation_getRenormGroup
  implicit none

#include "constants.h"
#include "Flash.h"

  integer, intent(in), dimension(2,MDIM)::blkLimits
  real,pointer :: solnData(:,:,:,:)

  integer :: i, j, k, n, gp
  
  !! +1 ensures compiler doesn't complain about 0'sized array
  !! MASS_SCALAR_GROUP numbering starts from 1 (0 means no normalization)
  real :: group_sum(NMASS_SCALAR_GROUPS+1)

  !! Common case, speed it up
  if (NMASS_SCALAR_GROUPS .eq. 0) return

  do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
     do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
               
           group_sum = 0.e0

           do n = MASS_SCALARS_BEGIN, MASS_SCALARS_END
              !! which group do I belong to?
              call Simulation_getRenormGroup(n,gp)
              if (gp .gt. 0) then
                 !! add to group tally
                 group_sum(gp) = group_sum(gp) + solnData(n,i,j,k)
              endif
           enddo

           !! Check group sum's to see if they are messed up
           !! No check done here. Should we check if sum is really small?
           do gp = 1, NMASS_SCALAR_GROUPS+1
              group_sum(gp) = 1.e0 / group_sum(gp) !! Invert sum
           enddo

           !! scale each group back
           do n = MASS_SCALARS_BEGIN, MASS_SCALARS_END
              !! which group do I belong to?
              call Simulation_getRenormGroup(n,gp)
              if (gp .gt. 0) then
                 solnData(n,i,j,k) = solnData(n,i,j,k)*group_sum(gp)
              endif
           enddo
 

        enddo !! end i loop
     enddo !! end j loop
  enddo !! end k loop
  
  return
end subroutine Grid_renormMassScalars

