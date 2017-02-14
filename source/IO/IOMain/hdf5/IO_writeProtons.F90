!!****if* source/IO/IOMain/hdf5/IO_writeProtons
!!
!! NAME
!!
!!  IO_writeProtons
!!
!! SYNOPSIS
!!
!!  call IO_writeProtons (integer (in) :: numberOfProtons,
!!                        integer (in) :: tags,
!!                        real    (in) :: points,
!!                        integer (in) :: pointCount)
!!
!! DESCRIPTION
!!
!!  Writes a collection of IO protons to the HDF5 plot file into the 'ProtonData'
!!  dataset of that plot file.
!! 
!! ARGUMENTS
!!
!!   numberOfProtons : number of protons 
!!   tags            : the proton tags
!!   points          : the proton points  
!!   pointCount      : the number of points for each proton
!!
!!***

subroutine IO_writeProtons (numberOfProtons, tags, points, pointCount)

  use IO_data,          ONLY: io_wrotePlot,       &
                              io_meshComm,        &
                              io_meshMe,          &
                              io_meshNumProcs,    &
                              io_protonFileID

  use Driver_interface, ONLY: Driver_abortFlash

  implicit none

#include "constants.h"
#include "Flash_mpi.h"

  integer, intent (in) :: numberOfProtons
  integer, intent (in) :: tags       (:)
  real,    intent (in) :: points     (:,:,:)
  integer, intent (in) :: pointCount (:)

  integer :: count
  integer :: error
  integer :: globalPoints, localPoints, startPoint
  integer :: point
  integer :: proton

  real    :: tag

  integer, allocatable :: procPoints (:)
  real,    allocatable :: t (:), x (:), y (:), z (:)

  if (.not. io_wrotePlot) then
       call Driver_abortFlash ("[IO_writeProtons] IO_writeProtons should only be called after a plot")
  end if
  
  allocate (procPoints (1:io_meshNumProcs))
!
!
!     ...Before writing the proton data, each processor has to know the total
!        number of proton points it has to write. This information must be
!        shared with all other processors to compute the total number of proton
!        points globally. Also compute the number of proton points owned by
!        processors whose rank is less than mine. This information os used
!        to determine where in the 'ProtonData' dataset of the HDF5 plot file
!        each processor writes.
!
!
  localPoints = sum (pointCount (1:numberOfProtons))

  call MPI_Allgather (localPoints, &
                      1,           &
                      MPI_INTEGER, &
                      procPoints,  &
                      1,           &
                      MPI_INTEGER, &
                      io_meshComm, &
                      error        )

  globalPoints = sum (procPoints (1:io_meshNumProcs))

  if (io_meshMe > 0) then
      startPoint = sum (procPoints (1:io_meshMe))
  else
      startPoint = 0
  end if

  deallocate (procPoints)
!
!
!     ...Collect all proton data into 4 1D arrays and write to the HDF5 plot file.
!
!
  allocate (t (1:localPoints))
  allocate (x (1:localPoints))
  allocate (y (1:localPoints))
  allocate (z (1:localPoints))

  point = 0
  do proton = 1, numberOfProtons
     tag = real (tags (proton))
     count = pointCount (proton)
     t (point+1:point+count) = tag
     x (point+1:point+count) = points (1:count,proton,IAXIS)
     y (point+1:point+count) = points (1:count,proton,JAXIS)
     z (point+1:point+count) = points (1:count,proton,KAXIS)
     point = point + count
  end do
   
  call io_h5write_protondata (io_protonFileID, &
                              io_meshMe,       &
                              localPoints,     &
                              globalPoints,    &
                              startPoint,      &
                              error,           &
                              t, x, y, z       )

  if (error < 0) then
      call Driver_abortFlash ("[IO_writeProtons] Error in io_h5write_protonData")
  end if

  deallocate (t,x,y,z)

end subroutine IO_writeProtons
