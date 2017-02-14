!!****if* source/IO/IOMain/hdf5/IO_writeRays
!!
!! NAME
!!
!!  IO_writeRays
!!
!! SYNOPSIS
!!
!!  call IO_writeRays(integer(in) :: numrays,
!!                    integer(in) :: raytags,
!!                    real(in) :: posbuffer,
!!                    real(in) :: powerbuffer,
!!                    integer(in) :: numpos)
!!
!! DESCRIPTION
!!
!! Write rays 
!! 
!! ARGUMENTS
!!
!!   numrays : number of rays  
!!
!!   raytags : tags of rays 
!!
!!   posbuffer : position buffer  
!!
!!   powerbuffer : power buffer 
!!
!!   numpos : 
!!
!!
!!
!!***

subroutine IO_writeRays(numRays, rayTags, posBuffer, powerBuffer, numPos)
  use IO_data, ONLY: io_wrotePlot,       &
                     io_oldPlotFileName, &
                     io_meshComm,        &
                     io_meshMe,          &
                     io_meshNumProcs,    &
                     io_rayFileID

  use Driver_interface, ONLY: Driver_abortFlash
  implicit none

#include "constants.h"
#include "Flash_mpi.h"

  integer, intent(in) :: numRays
  integer, intent(in) :: rayTags(:)
  real,    intent(in) :: posBuffer(:,:,:)
  real,    intent(in) :: powerBuffer(:,:)
  integer, intent(in) :: numPos(:)

  integer :: i, j, local_count, nints, ierr, global_count, start_pos, count
  integer, allocatable :: procPos(:)
  real, allocatable :: tags(:), xpos(:), ypos(:), zpos(:), power(:)

  allocate(procPos(io_meshNumProcs))

  if(.not. io_wrotePlot) then
     call Driver_abortFlash("[IO_writeRays] IO_writeRays should only be called after a plot")
  end if
  
  ! Before writing the ray data, each processor has to compute the
  ! total number of ray positions it has to write. This information
  ! must be shared with all other processors.
  local_count = 0
  do i = 1, numRays
     local_count = local_count + numPos(i)
  end do

  ! Now, share with everyone else...
  nints = 1
  call MPI_Allgather(local_count, nints, MPI_INTEGER, procPos,&
       nints, MPI_INTEGER, io_meshComm,ierr)

  ! Compute the total number of ray positions globally
  global_count = 0
  do i = 1, io_meshNumProcs
     global_count = global_count + procPos(i)
  end do  

  ! Compute the number of ray positions owned by processors whose rank
  ! is less than mine. This information is used to determine where in
  ! the RayPos dataset each processor writes:
  start_pos = 0
  if(io_meshMe > 0) then
     do i = 1, io_meshMe
        start_pos = start_pos + procPos(i)
     end do
  end if
  deallocate(procPos)

  ! Collect all ray data into 5 1D arrays:
  allocate(tags(local_count))
  allocate(xpos(local_count))
  allocate(ypos(local_count))
  allocate(zpos(local_count))
  allocate(power(local_count))

  count = 0
  do i = 1, numRays
     do j = 1, numPos(i)
        count = count + 1
        tags(count) = rayTags(i)
        xpos(count) = posBuffer(i,j,IAXIS)
        ypos(count) = posBuffer(i,j,JAXIS)
        zpos(count) = posBuffer(i,j,KAXIS)
        power(count) = powerBuffer(i,j)
     end do
  end do
   
  ! Write the data:
  call io_h5write_raydata(io_rayFileID, global_count, ierr, start_pos, &
       local_count, tags, xpos, ypos, zpos, power, io_meshMe)
  if(ierr < 0) then
     call Driver_abortFlash("[IO_writeRays] Error in io_h5write_raydata")
  end if

  deallocate(tags)
  deallocate(xpos)
  deallocate(ypos)
  deallocate(zpos)
  deallocate(power)

end subroutine IO_writeRays
