!!****f* source/Grid/Grid_conserveField
!!
!! NAME
!!
!!  Grid_conserveField
!!
!! SYNOPSIS
!!
!!  Grid_conserveField ()
!!
!! ARGUMENTS
!!
!!
!! DESCRIPTION
!! 
!!  Corrects electric fields at refinement jump boundaries to make 
!!  sure electric fields at the common boundaries are consistent 
!!  at both refined and coarse meshes.
!!
!!  On the uniform grid this routine is not needed.
!!
!!***



subroutine Grid_conserveField ()


  implicit none

  return

end subroutine Grid_conserveField
