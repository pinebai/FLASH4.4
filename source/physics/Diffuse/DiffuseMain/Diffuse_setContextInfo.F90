!!****if* source/physics/Diffuse/DiffuseMain/Diffuse_setContextInfo
!!
!! NAME
!!
!!  Diffuse_setContextInfo
!!
!! SYNOPSIS
!!
!!  call Diffuse_setContextInfo(integer(in),OPTIONAL  :: group,
!!                              integer(in),OPTIONAL  :: component)
!!
!! DESCRIPTION
!!
!! Sets the context info
!!
!! ARGUMENTS
!!
!!   group : group number 
!!
!!   component : component number
!!
!!
!!
!!***

subroutine Diffuse_setContextInfo(group,component)
  use Diffuse_data, ONLY: diff_dbgContext
  implicit none
  integer,intent(in),OPTIONAL :: component, group

  if (.NOT. present(group) .AND. .NOT. present(component)) then
     ! clear the info
     diff_dbgContext%component = 0
     diff_dbgContext%group     = 0
  else if (present(group)) then
     diff_dbgContext%group     = group
     if (present(component)) then
        diff_dbgContext%component = component
     else
        diff_dbgContext%component = 3
     end if
  else
     diff_dbgContext%component = component
  end if

end subroutine Diffuse_setContextInfo
