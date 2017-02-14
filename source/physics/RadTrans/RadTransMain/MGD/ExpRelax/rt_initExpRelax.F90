!!****if* source/physics/RadTrans/RadTransMain/MGD/ExpRelax/rt_initExpRelax
!!
!!  NAME 
!!
!!  rt_initExpRelax
!!
!!  SYNOPSIS
!!
!!  call rt_initExpRelax()
!!
!!  DESCRIPTION 
!!    Initialize additional data for using the ExpRelax solver for radiative transfer
!!
!!***
subroutine rt_initExpRelax
  use rt_data, ONLY: rt_useMGD, rt_mgdthetaC, rt_mgdthetaD, rt_mgdthetaImplct, &
       rt_tightIonCoupling
  use rt_expData, ONLY: rt_expRelaxMaxIter
  use RadTrans_data, ONLY: rt_useRadTrans, rt_meshCopyCount, &
       rt_useRadTrans,                                       &
       rt_meshMe, rt_acrossMe
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  implicit none

#include "constants.h"

  
  call RuntimeParameters_get('rt_expRelaxMaxIter', rt_expRelaxMaxIter)
  if (rt_useMGD .AND. .NOT. rt_useRadTrans) then
     if (rt_meshME==MASTER_PE .AND. &
          rt_acrossMe==MASTER_PE) then
        print *,' Forcing rt_useMGD to FALSE because useRadTrans is FALSE!'
     end if
     rt_useMGD = .FALSE.
  end if

end subroutine rt_initExpRelax
