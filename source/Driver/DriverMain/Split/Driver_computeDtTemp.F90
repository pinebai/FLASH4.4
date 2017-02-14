!!****if* source/Driver/DriverMain/Split/Driver_computeDtTemp
!!
!! NAME
!!
!!  Driver_computeDtTemp
!!
!!
!! SYNOPSIS
!!
!!  Driver_computeDtTemp() 
!!                       real(in) :: dt,
!!                       real(inout) :: dtTemperature,  
!!                       integer(inout) :: dtMinLoc(5), 
!!                       integer(in) :: blockID)
!!
!!
!!
!! DESCRIPTION
!!
!!  compute the timestep from limiting the temperature change to
!!  be smaller than temp_factor
!!
!!                                     T
!!             dt    =  temp_factor * --- * dt
!!               new                  dT      old
!!
!!  Note: during the hydro evolution, the temperature at the top of a sweep 
!!        is stored in unksm(itemp_old,...).  Since the timestep routine is
!!        called after refinements are done, and unksm does not get refined,
!!        the data in that array will not correspond to unk after refining.  
!!        To circumvent this problem, unksm is refilled with T/dT just before
!!        refining.  Any invalid blocks are flagged with a -1 (i.e. they are
!!        not leaf blocks).
!!
!! ARGUMENTS
!!
!!   myPE             --  current Processor 
!!   dt               --  current timestep
!!   dtTemperature    --  limiting timestep due to temperature
!!   dtMinLoc(5)      --  array to hold limiting zone info:  index[1-3] are zone
!!                        indices (i,j,k), index[4]=block ID, index[5]=processor 
!!                        number. The zone indices indicate the zone containing
!!                        the particle whose velocity set the timestep.
!!   blockID          --  current block
!! 
!! NOTES
!!  Since this constraint is computed after a timestep is taken, it does 
!!  not guarantee that we have kept the temperature change below temp_factor,
!!  it simply creates a prediction of a timestep that, based on the 
!!  previous one, will meet this criteria.
!!
!!***

subroutine Driver_computeDtTemp(myPE, dt, dtTemperature, dtMinLoc, blockID)

  use Driver_data, ONLY : dr_tempFactor
  use Grid_interface, ONLY : Grid_getBlkPtr, &
    Grid_getBlkIndexLimits, Grid_releaseBlkPtr

  implicit none

#include "Flash.h"
#include "constants.h"


  integer, intent(in) :: blockID
  integer, intent(inout) ::  dtMinLoc(5)  
  integer, intent(in) :: myPE
  real, intent(inout) :: dtTemperature
  real, intent(in)    :: dt
 
  real :: ldtTemperature
  

  integer :: i,j,k, temploc(5)

  real :: dtemp, dtemp_min

  real, DIMENSION(:,:,:,:), POINTER :: dataPtr

  integer, save :: itemp_old

  integer,dimension(2,MDIM)::blkLimits,blkLimitsGC



! compute the smallest T/dT
  dtemp_min = 1.e30

  call Grid_getBlkPtr(blockID, dataPtr, SCRATCH)

  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)


  do k = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
     do j = blkLimitsGC(LOW,JAXIS), blkLimitsGC(HIGH,JAXIS)
        do i = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH,IAXIS)

#ifdef OTMP_SCRATCH_GRID_VAR
  !         dtemp = dataPtr(OTMP_SCRATCH_GRID_VAR,i,j,k)
#endif

           ! since the grid may have refined since we computed dT/T, we flagged the
           ! bad data (non-leaf blocks) with -1.
           if (dtemp .LT. dtemp_min .AND. dtemp .GE. 0.e0) then
              dtemp_min = dtemp
              temploc(1) = i
              temploc(2) = j
              temploc(3) = k
              temploc(4) = blockID
              temploc(5) = MyPE
           endif
               
        enddo
     enddo
  enddo

! now compute the temperature timestep
  ldtTemperature = dr_tempFactor*dtemp_min*dt

  if (ldtTemperature .LT. dtTemperature) then
     dtTemperature= ldtTemperature
     dtMinLoc = temploc
  endif
  
  call Grid_releaseBlkPtr(blockID, dataPtr,SCRATCH)
  return
end subroutine Driver_computeDtTemp

               
      




