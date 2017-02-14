!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/method/LowTemp/op_setPEcoeffsAij4
!!
!! NAME
!!
!!  op_setPEcoeffsAij4
!!
!! SYNOPSIS
!!
!!  call op_setPEcoeffsAij4 ()
!!
!! DESCRIPTION
!!
!!  This routine sets the A(i,j,4) expansion coefficients for determining the
!!  photoelectric cross sections according to the F.Biggs and R.Lighthill report:
!!
!!       Analytical Approximations for X-Ray Cross Sections II
!!       Frank Biggs and Ruth Lighthill
!!       Weapons Effects Research Department
!!       Sandia Laboratories, December 1971
!!
!!  The data is NOT the updated data from the 1988 update of the report!
!!  This routine can only be called after all the necessary arrays have
!!  been allocated. It sets the coefficients in three separate routines.
!!  Providing only one routine to set the entire coefficients leads to
!!  exaggerated compilation times!
!!
!! ARGUMENTS
!!
!!***
subroutine op_setPEcoeffsAij4 ()

  use op_interface,     ONLY : op_setPEcoeffsAij4atomsZ43,  &
                               op_setPEcoeffsAij4atomsZ73,  &
                               op_setPEcoeffsAij4atomsZ100

  implicit none
!
!
!   ...Set the coefficient array in three portions.
!
!
  call op_setPEcoeffsAij4atomsZ43  ()      ! sets coefficients for atoms:  1 =< Z =< 43
  call op_setPEcoeffsAij4atomsZ73  ()      ! sets coefficients for atoms: 44 =< Z =< 73
  call op_setPEcoeffsAij4atomsZ100 ()      ! sets coefficients for atoms: 74 =< Z =< 100
!
!
!   ...Ready! 
!
!
  return
end subroutine op_setPEcoeffsAij4
