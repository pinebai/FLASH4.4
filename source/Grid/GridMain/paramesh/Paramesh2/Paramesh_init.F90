!!****if* source/Grid/GridMain/paramesh/Paramesh2/Paramesh_init
!!
!! NAME
!!  mesh_intialize
!! 
!! SYNOPSIS
!!  call Paramesh_init()
!!
!! DESCRIPTION
!!  Simple interface to mesh initialization routine.
!!  Each mesh package will have a different copy of this function.
!!
!! USED BY
!!  Grid_initMesh
!!***
subroutine Paramesh_init()
  use paramesh_interfaces , ONLY : amr_initialize
  implicit none
  call amr_initialize()

!!  call init_flash_physicaldata()
end subroutine Paramesh_init
