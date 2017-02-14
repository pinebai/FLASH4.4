!!****if* source/Particles/ParticlesInitialization/pt_findTagOffset
!!
!! NAME
!!    pt_findTagOffset
!!
!! SYNOPSIS
!!
!!    call pt_findTagOffset(integer(in) :: newCount, integer(OUT) :: tagOffset)
!!
!! DESCRIPTION
!!    Add a unique tag each to particles generated during evolution.
!!    The algorithm first finds out the sum of number of particles 
!!    being added to all the processors to the left of MyPE. It then
!!    adds pt_startTagNumber (the largest tag number in use) to it.
!!    This number acts as the offset for the tag numbers being assigned
!!    to the newly added particles.
!!
!! NOTES
!!    This method of tag generation will work for up to 10^14 
!!    particles in a simulation.
!!
!!
!!
!!
!!***

subroutine pt_findTagOffset(newCount,tagOffset)

  use Particles_data, ONLY : pt_meshMe,pt_meshNumProcs,&
       pt_meshComm, pt_startTagNumber, pt_resetTag

  implicit none

  include "Flash_mpi.h"
  
  integer, intent(IN) :: newCount
  integer, intent(OUT) :: tagOffset
  
  integer ::  count, ierr, MyTagStart, maxTagOwner, maxTag,i
  
  integer, allocatable ::  localNumParticles(:)

  if(pt_resetTag) then  
     pt_startTagNumber=0
  end if

  allocate(localNumParticles(pt_meshNumProcs))

  count = 1
  call MPI_Allgather(newCount, count, MPI_INTEGER, localNumParticles,&
       count, MPI_INTEGER, pt_meshComm,ierr)

  MyTagStart=0
  do i = 1,pt_meshMe
     MyTagStart=MyTagStart+localNumParticles(i)
  end do
  tagOffset=MyTagStart+pt_startTagNumber

  deallocate(localNumParticles)

  ! Keep pt_startTagNumber properly updated, if it may be used in another call to pt_findTagOffset.
  if (.NOT. pt_resetTag) then
     maxTagOwner=pt_meshNumProcs-1
     if(pt_meshMe==maxTagOwner)pt_startTagNumber=tagOffset+newCount
     call MPI_Bcast(pt_startTagNumber ,count,MPI_INTEGER,&
          maxTagOwner, pt_meshComm, ierr)
  end if
  return
end subroutine pt_findTagOffset
