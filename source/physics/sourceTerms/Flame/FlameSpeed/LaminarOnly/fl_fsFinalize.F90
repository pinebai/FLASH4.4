!!****if* source/physics/sourceTerms/Flame/FlameSpeed/LaminarOnly/fl_fsFinalize
!!
!! NAME
!!
!!  fl_fsFinalize
!!
!! SYNOPSIS
!!
!!  call fl_fsFinalize()
!!
!! DESCRIPTION
!!
!! Dean Townsley 2008
!!
!!
!! ARGUMENTS
!!
!!   No arguments
!!
!!
!!
!!***


subroutine fl_fsFinalize()

   implicit none

   call fl_fsLaminarFinalize
   call fl_fsTFIFinalize

   return

end subroutine fl_fsFinalize
