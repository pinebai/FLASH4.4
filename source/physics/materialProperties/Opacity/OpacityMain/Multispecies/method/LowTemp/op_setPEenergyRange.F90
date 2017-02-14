!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/method/LowTemp/op_setPEenergyRange
!!
!! NAME
!!
!!  op_setPEenergyRange
!!
!! SYNOPSIS
!!
!!  call op_setPEenergyRange ()
!!
!! DESCRIPTION
!!
!!  This routine sets the energy range (boundaries) for determining the photoelectric
!!  cross sections according to the F.Biggs and R.Lighthill report:
!!
!!       Analytical Approximations for X-Ray Cross Sections II
!!       Frank Biggs and Ruth Lighthill
!!       Weapons Effects Research Department
!!       Sandia Laboratories, December 1971
!!
!!  The data is NOT the updated data from the 1988 update of the report!
!!  This routine can only be called after all the necessary arrays have
!!  been allocated. It sets the energy ranges in three separate routines.
!!  Providing only one routine to set the entire energy ranges leads to
!!  exaggerated compilation times!
!!
!! ARGUMENTS
!!
!!***
subroutine op_setPEenergyRange ()

  use op_interface,     ONLY : op_setPEenergyRangeAtomsZ43,  &
                               op_setPEenergyRangeAtomsZ73,  &
                               op_setPEenergyRangeAtomsZ100

  implicit none
!
!
!   ...Set the energy range array in three portions.
!
!
  call op_setPEenergyRangeAtomsZ43  ()      ! sets energy range for atoms:  1 =< Z =< 43
  call op_setPEenergyRangeAtomsZ73  ()      ! sets energy range for atoms: 44 =< Z =< 73
  call op_setPEenergyRangeAtomsZ100 ()      ! sets energy range for atoms: 74 =< Z =< 100
!
!
!   ...Ready! 
!
!
  return
end subroutine op_setPEenergyRange
