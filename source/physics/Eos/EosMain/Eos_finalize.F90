!!****if* source/physics/Eos/EosMain/Eos_finalize
!!
!! NAME
!!
!!  Eos_finalize
!!
!! SYNOPSIS
!!
!!  Eos_finalize()
!!
!! DESCRIPTION
!!
!!  Deallocates any memory that has been allocated in the Eos Unit
!!  and prepares the unit for shutdown
!!
!!
!!***


subroutine Eos_finalize()

  use eos_localInterface, ONLY: eos_tabFinalize
  implicit none

  call eos_tabFinalize

end subroutine Eos_finalize
