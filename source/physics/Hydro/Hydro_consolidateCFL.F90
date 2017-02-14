!!****f* source/physics/Hydro/Hydro_consolidateCFL
!!
!! NAME
!!
!!  Hydro_consolidateCFL
!!
!! SYNOPSIS
!!
!!  call Hydro_consolidateCFL()
!!
!! DESCRIPTION
!!
!!  Consolidate the CFL factors of different MPI tasks by a global
!!  reduction operation.
!!
!!  Different MPI tasks may have computed different values for an updated Hydro
!!  CFL factor. This routine should be called before Hydro_computeDt is called
!!  to compute an updated Hydro time step.
!!
!! ARGUMENTS
!!
!!   No arguments
!!
!! SIDE EFFECTS
!!
!!  The CFL factor used by the Hydro unit is updated.
!!
!! NOTES
!!
!!  This routine must be called collectively by all MPI tasks in the global
!!  Hydro communicator.
!!
!!  For split Hydro implementations, there is only a stub implementation of
!!  this interface which does not do any computation or communication.
!!  That is okay since unsplit Hydro implementations are not expected to
!!  modify the CFL factor that is given as a runtime parameter.
!!***

subroutine Hydro_consolidateCFL()

  implicit none

end subroutine Hydro_consolidateCFL
