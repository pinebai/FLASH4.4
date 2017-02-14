!!****if* source/physics/sourceTerms/Flame/FlameMain/Flame_finalize
!!
!! NAME
!!
!!  Flame_finalize
!!
!! SYNOPSIS
!!
!!  call Flame_finalize()
!!
!! DESCRIPTION
!!
!! Dean Townsley 2008
!! Finalize
!!
!! ARGUMENTS
!!
!!   No arguments
!!
!!
!!
!!***

subroutine Flame_finalize()

  use Flame_data, ONLY: fl_useFlame

  implicit none

  if(fl_useFlame) then 
     call fl_fsFinalize
     call fl_effFinalize
  end if

  return

end subroutine Flame_finalize
