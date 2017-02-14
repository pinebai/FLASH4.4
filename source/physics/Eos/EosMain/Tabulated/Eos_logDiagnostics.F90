!!****if* source/physics/Eos/EosMain/Tabulated/Eos_logDiagnostics
!!
!! NAME
!!
!!  Eos_logDiagnostics
!!
!! SYNOPSIS
!!
!!  call Eos_logDiagnostics(logical, intent(IN)  :: force)
!!
!! DESCRIPTION
!! 
!! Stub
!!
!! ARGUMENTS
!!
!!   force : logical switch
!!
!!
!!
!!***

subroutine Eos_logDiagnostics(force)
  use eos_tabInterface,ONLY: eos_tabLogOutsideCounts, eos_tabZeroOutsideCounts
  implicit none
  logical, intent(IN) :: force

  call eos_tabLogOutsideCounts(force)
  call eos_tabZeroOutsideCounts()

end subroutine Eos_logDiagnostics
