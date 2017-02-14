!!****f* source/IO/IO_writeRays
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
  implicit none

  integer, intent(in) :: numRays
  integer, intent(in) :: rayTags(:)
  real, intent(in) :: posBuffer(:,:,:)
  real, intent(in) :: powerBuffer(:,:)
  integer, intent(in) :: numPos(:)

end subroutine IO_writeRays
