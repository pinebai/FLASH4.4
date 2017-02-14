!!****if* source/physics/sourceTerms/Flame/FlameSpeed/turbulent/fl_fsTFIFunc
!!
!! NAME
!!
!!  fl_fsTFIFunc
!!
!! SYNOPSIS
!!
!!  call fl_fsTFIFunc(real(in) :: x)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   x : input
!!
!!
!!
!!***

real function fl_fsTFIFunc(x)
   implicit none
   real, intent(in) :: x

   fl_fsTFIFunc = x

end function fl_fsTFIFunc
