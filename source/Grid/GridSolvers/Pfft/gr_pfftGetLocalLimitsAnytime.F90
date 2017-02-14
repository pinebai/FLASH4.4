!!****if* source/Grid/GridSolvers/Pfft/gr_pfftGetLocalLimitsAnytime
!!
!! NAME 
!!
!!   gr_pfftGetLocalLimitsAnytime
!!
!! SYNOPSIS
!!
!!   gr_pfftGetLocalLimitsAnytime(integer(IN) :: axis1,
!!                          integer(IN) :: axis2,
!!                          integer(IN) :: meaxis,
!!                          integer(IN) :: currentGridShape(MDIM),
!!                          integer(IN) :: baseDatType,
!!                 OPTIONAL,integer(INOUT) :: currentLocalLimits(LOW:HIGH,MDIM))
!!
!! DESCRIPTION 
!!
!!  Find out the local index limits of the transformed
!!  data in a processor. This imformation is necessary to
!!  calculate the derivatives. Here the situation is complicated
!!  by the fact that after transformation, the order of physical
!!  dimension is <JAXIS,IAXIS> in 2D and <KAXIS,IAXIS,JAXIS> in 3D.
!!  So the ordering of dimensions in the array is different from the
!!  physical dimension order. The first two arguments passed to this 
!!  function provide the mapping.
!!
!!  Results are returned in the currentLocalLimits(:,axis1) slice of
!!  the INTENT(OUT) argument, currentLocalLimits.
!!
!! ARGUMENTS
!!
!!   axis1  - the first array dimension of the transformed data
!!   axis2  - the corresponding physical data dimension
!!   meaxis - the axis along which the information is sought
!!   currentGridShape - local shape per processor.
!!   baseDatType - self explanatory
!!   currentLocalLimits - results of this call are returned in
!!            currentLocalLimits(:,axis1), other elements of 
!!            currentLocalLimits are left unmodified, if the
!!            optional argument is present.
!!
!! SIDE EFFECTS
!!
!!  Results of this call are returned in pfft_localLimits(:,axis1)
!!  if the optional argument currentLocalLimits is not present.
!!  Other elements of pfft_localLimits are always left unmodified.
!!
!!***

!currentGridShape is the local shape per processor.

subroutine gr_pfftGetLocalLimitsAnytime(axis1,axis2,meaxis,currentGridShape,baseDatType,currentLocalLimits)

#include "constants.h"
#include "Pfft.h"

  use gr_pfftData, ONLY : pfft_me,pfft_globalLen, pfft_localLimits, &
       pfft_pclBaseDatType, pfft_outLen

  implicit none

  integer,intent(IN), target :: axis1
  integer,intent(IN) :: axis2
  integer,intent(IN), OPTIONAL, target ::  meaxis
  integer, dimension(MDIM), intent(IN), OPTIONAL, target :: currentGridShape
  integer, intent(IN), OPTIONAL, target :: baseDatType
  integer, intent(INOUT), OPTIONAL, target :: currentLocalLimits(LOW:HIGH,MDIM)

  integer :: len,s,e
  integer,pointer :: localLimPtr(:,:), curGridShapePtr(:)
  integer,pointer :: meaxisPtr,baseDatTypePtr

  if (present(meaxis)) then
     meaxisPtr => meaxis
  else
     meaxisPtr => axis1
  end if
  if (present(currentGridShape)) then
     curGridShapePtr => currentGridShape
  else
     curGridShapePtr => pfft_outLen
  end if
  if (present(baseDatType)) then
     baseDatTypePtr => baseDatType
  else
     baseDatTypePtr => pfft_pclBaseDatType(axis2)
  end if
  if (present(currentLocalLimits)) then
     localLimPtr => currentLocalLimits
  else
     localLimPtr => pfft_localLimits
  end if

  len=pfft_globalLen(axis2)
  if (baseDatTypePtr==PFFT_PCLDATA_COMPLEX .OR. &
       baseDatTypePtr==PFFT_PCLDATA_COMPLEX_STUFFED .OR. &
       baseDatTypePtr==PFFT_PCLDATA_COMPLEX_EXTENDED) then
     len=len/2
     if(baseDatTypePtr==PFFT_PCLDATA_COMPLEX_EXTENDED) then
        len=len+1
     end if
  end if
  
  s=pfft_me(meaxisPtr)*curGridShapePtr(axis1)

  localLimPtr(LOW,axis1) = min(s, len)
  e = s + curGridShapePtr(axis1) - 1
  if (e .GE. len) e = len - 1
  localLimPtr(HIGH,axis1) = e

end subroutine gr_pfftGetLocalLimitsAnytime
