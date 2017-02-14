!!****if* source/Particles/ParticlesMain/active/Sink/Particles_sinkSumAttributes
!!
!! NAME
!!
!!  Particles_sinkSumAttributes
!!
!! SYNOPSIS
!!
!!  call Particles_sinkSumAttributes(  real(OUT) :: sums(:),
!!                                   integer(in) :: attribs(:),
!!                         OPTIONAL, integer(in) :: factor)
!!
!! DESCRIPTION
!!
!!  Compute global sums over all sink particles for the attributes in
!!  the list ATTRIBS, optionally multiplied by the value of another
!!  particle attribute FACTOR.
!!
!!
!! ARGUMENTS
!!
!!   sums - the sums are returned here
!!
!!   attribs - list of sink particle properties
!!
!!   factor  - a sink particle property to use as an optional factor
!!
!! NOTES
!!
!!   written by Klaus Weide, 2014
!!
!!***

subroutine Particles_sinkSumAttributes(sums, attribs, factor)

  use Particles_sinkData, ONLY : particles_global, localnpf
  use pt_sinkInterface, only: pt_sinkGatherGlobal
  use Particles_data, ONLY : pt_globalMe

  implicit none

#include "constants.h"

  real,intent(OUT)   :: sums(:)
  integer,intent(in) :: attribs(:)
  integer,intent(in),OPTIONAL :: factor


  integer            :: i, iout
  logical            :: doMult
  real               :: fact
  integer, save      :: MyPE, MasterPE
  logical, save      :: firstCall = .TRUE.

  integer :: gather_nprops
  integer, allocatable, dimension(:) :: gather_propinds

  doMult = .FALSE.

  gather_nprops = size(attribs)

  if (present(factor)) then
     if (factor > 0) then
        doMult = .TRUE.
        if (ALL(attribs(:).NE.factor)) gather_nprops = gather_nprops + 1
     end if
  end if

  allocate(gather_propinds(gather_nprops))
  gather_propinds(1:size(attribs)) = attribs
  if (doMult .AND. gather_nprops > size(attribs)) then
     gather_propinds(gather_nprops) = factor
  end if


  if (firstCall) then

     MyPE     = pt_globalMe
     MasterPE = MASTER_PE
     firstCall = .false.

  endif

  ! exchange particle information across CPUs (this needs to be called by all processors)
  call pt_sinkGatherGlobal(gather_propinds, gather_nprops)

  do iout = 1, size(sums)
     sums(iout) = 0.0
  end do

  ! only the master processor sums the data of all particles
  if (pt_globalMe .NE. MasterPE) return

  ! now write the actual data to the file
  do i = 1, localnpf

     if (doMult) then
        fact = particles_global(factor,i)
     else
        fact = 1.0
     end if

     do iout = 1, size(sums)
        sums(iout) = sums(iout) + fact * particles_global(gather_propinds(iout),i)
     end do

  enddo

  return

end subroutine Particles_sinkSumAttributes
