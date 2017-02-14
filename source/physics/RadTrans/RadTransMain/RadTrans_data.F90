!!****if* source/physics/RadTrans/RadTransMain/RadTrans_data
!!
!!  NAME 
!!
!!  RadTrans_data
!!
!!  SYNOPSIS
!!   use RadTrans_data
!!
!!  DESCRIPTION 
!!    Stores local data for the RadTrans unit
!!
!!***
module RadTrans_data
  use RadTrans_interfaceTypeDecl, ONLY: RadTrans_dbgContext_t
  implicit none
  
  integer, save :: rt_meshMe   ! Process rank in the mesh communicator
  integer, save :: rt_globalMe ! Global process rank
  integer, save :: rt_globalNumProcs
  logical, save :: rt_useRadTrans
  
  ! Physical constants:
  real, save :: rt_radconst ! Radiation constant
  real, save :: rt_speedlt ! Speed of light
  real, save :: rt_speedlt3 ! a third thereof
  real, save :: rt_boltz ! Boltzmann constant

  real, save :: rt_dtFactor ! Coefficient for RadTrans time step

  ! The number of replicated meshes active in this simulation
  integer, save :: rt_meshCopyCount

  ! The mesh number of this process
  integer, save :: rt_acrossMe

  ! The across communicator
  integer, save :: rt_acrossComm

  ! Global communicator for all processes
  integer, save :: rt_globalComm

  ! Structure that holds context information on the current operation,
  ! for debugging
  type(RadTrans_dbgContext_t),save,target :: rt_dbgContext

end module RadTrans_data
