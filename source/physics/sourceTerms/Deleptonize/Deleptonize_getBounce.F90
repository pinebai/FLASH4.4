!!****f* source/physics/sourceTerms/Deleptonize/Deleptonize_getBounce
!!
!! NAME
!!  
!!  Deleptonize_getBounce 
!!
!!
!! SYNOPSIS
!! 
!!  call Deleptonize_getBounce (logical(OUT) :: postBounce,
!!                              real(OUT) :: bouncTime)
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

  implicit none
  
  logical, intent(OUT) :: postBounce
  real, intent(OUT) :: bounceTime

  postBounce = .FALSE.
  bounceTime = 0.

  return
end subroutine Deleptonize_getBounce
