!!****if* source/Grid/GridSolvers/Grid_setSolverDbgContextInfo
!!
!! NAME
!!
!!  Grid_setSolverDbgContextInfo
!!
!! SYNOPSIS
!!
!!  call Grid_setSolverDbgContextInfo(integer(in),OPTIONAL  :: component,
!!                                    integer(in),OPTIONAL  :: group)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   component : 
!!
!!   group : 
!!
!!
!!
!!***

subroutine Grid_setSolverDbgContextInfo(component,group)
  use gr_solversData, ONLY: gr_solversDbgContext
  implicit none
  integer,intent(in),OPTIONAL :: component, group

  if (.NOT. present(group) .AND. .NOT. present(component)) then
     ! clear the info
     gr_solversDbgContext%component = 0
     gr_solversDbgContext%group     = 0
  else if (present(group)) then
     gr_solversDbgContext%group     = group
     if (present(component)) then
        gr_solversDbgContext%component = component
     else
        gr_solversDbgContext%component = 3
     end if
  else
     gr_solversDbgContext%component = component
  end if

  gr_solversDbgContext%libErrCode = 0
  gr_solversDbgContext%flashErrCode = 0

end subroutine Grid_setSolverDbgContextInfo
