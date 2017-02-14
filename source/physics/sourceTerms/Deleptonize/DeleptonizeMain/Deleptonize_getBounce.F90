!!****if* source/physics/sourceTerms/Deleptonize/DeleptonizeMain/Deleptonize_getBounce
!!
!! NAME
!!  
!!  Deleptonize_getBounce 
!!
!!
!! SYNOPSIS
!! 
!!  call Deleptonize_getBounce (logical(OUT) :: postBounce,
!!                              real(OUT) :: bounceTime)
!!  
!!  
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!  postBounce : logical stating whether bounce has been deteced
!!
!!***

subroutine Deleptonize_getBounce(postBounce,bounceTime)

  use Deleptonize_data, ONLY : delep_postBounce, delep_bounceTime

  implicit none
  
  logical, intent(OUT) :: postBounce
  real, intent(OUT) :: bounceTime

  postBounce = delep_postBounce
  bounceTime = delep_bounceTime

  return
end subroutine Deleptonize_getBounce
