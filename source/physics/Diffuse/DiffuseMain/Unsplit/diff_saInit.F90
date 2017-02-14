!!****if* source/physics/Diffuse/DiffuseMain/Unsplit/diff_saInit
!!
!! NAME
!!
!!  diff_saInit
!!
!!
!! SYNOPSIS
!!
!!  call diff_saInit()
!!
!! Description
!!
!!  Initializes local data for Unit Diffuse defined in Module diff_saData.
!!  All the variables here are initialized by calling the
!!  RuntimeParameters_get subroutine. These data variables are for
!!  Unit Scope ->  Diffuse.
!!
!! ARGUMENTS
!!
!!  none  
!!
!! PARAMETERS
!!
!!    diff_scaleFactThermSaTempDiff
!!        factor by which the solution is scaled.
!!    diff_scaleFactThermSaTime
!!        factor by which diffusion time step is scaled.
!!***

subroutine diff_saInit
  use RuntimeParameters_interface, ONLY: RuntimeParameters_get,RuntimeParameters_mapStrToInt
  use diff_saData, ONLY: diff_scaleFactThermSaTempDiff, diff_scaleFactThermSaTime, &
                         diff_eleDomainBC, diff_thetaImplct, diff_updEint
  use diff_saData, ONLY: diff_ionThetaImplct
  use diff_saData, ONLY: diff_ionDomainBC
  implicit none

#include "constants.h"

  character(len=MAX_STRING_LENGTH) :: xl_bcString,xr_bcString
  character(len=MAX_STRING_LENGTH) :: yl_bcString,yr_bcString
  character(len=MAX_STRING_LENGTH) :: zl_bcString,zr_bcString

  call RuntimeParameters_get('diff_scaleFactThermSaTempDiff',diff_scaleFactThermSaTempDiff)

  call RuntimeParameters_get('diff_scaleFactThermSaTime',diff_scaleFactThermSaTime)

  !get the boundary conditions stored as strings in the flash.par file

  ! Ion conduction
  call RuntimeParameters_get("diff_ionXlBoundaryType", xl_bcString)
  call RuntimeParameters_get("diff_ionXrBoundaryType", xr_bcString)
  call RuntimeParameters_get("diff_ionYlBoundaryType", yl_bcString)
  call RuntimeParameters_get("diff_ionYrBoundaryType", yr_bcString)
  call RuntimeParameters_get("diff_ionZlBoundaryType", zl_bcString)
  call RuntimeParameters_get("diff_ionZrBoundaryType", zr_bcString)

  !map the string boundary conditions to integer constants defined in constants.h
  call RuntimeParameters_mapStrToInt(xl_bcString,diff_ionDomainBC(1))
  call RuntimeParameters_mapStrToInt(xr_bcString,diff_ionDomainBC(2))
  call RuntimeParameters_mapStrToInt(yl_bcString,diff_ionDomainBC(3))
  call RuntimeParameters_mapStrToInt(yr_bcString,diff_ionDomainBC(4))
  call RuntimeParameters_mapStrToInt(zl_bcString,diff_ionDomainBC(5))
  call RuntimeParameters_mapStrToInt(zr_bcString,diff_ionDomainBC(6))

 
  ! Ele conduction
  call RuntimeParameters_get("diff_eleXlBoundaryType", xl_bcString)
  call RuntimeParameters_get("diff_eleXrBoundaryType", xr_bcString)
  call RuntimeParameters_get("diff_eleYlBoundaryType", yl_bcString)
  call RuntimeParameters_get("diff_eleYrBoundaryType", yr_bcString)
  call RuntimeParameters_get("diff_eleZlBoundaryType", zl_bcString)
  call RuntimeParameters_get("diff_eleZrBoundaryType", zr_bcString)

  !map the string boundary conditions to integer constants defined in constants.h
  call RuntimeParameters_mapStrToInt(xl_bcString,diff_eleDomainBC(1))
  call RuntimeParameters_mapStrToInt(xr_bcString,diff_eleDomainBC(2))
  call RuntimeParameters_mapStrToInt(yl_bcString,diff_eleDomainBC(3))
  call RuntimeParameters_mapStrToInt(yr_bcString,diff_eleDomainBC(4))
  call RuntimeParameters_mapStrToInt(zl_bcString,diff_eleDomainBC(5))
  call RuntimeParameters_mapStrToInt(zr_bcString,diff_eleDomainBC(6))

  call RuntimeParameters_get("diff_thetaImplct", diff_thetaImplct)
  call RuntimeParameters_get("diff_ionThetaImplct", diff_ionThetaImplct)
  
  call RuntimeParameters_get("diff_updEint",diff_updEint)

end subroutine diff_saInit
