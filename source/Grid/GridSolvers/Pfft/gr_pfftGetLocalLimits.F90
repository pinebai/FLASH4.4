!!****if* source/Grid/GridSolvers/Pfft/gr_pfftGetLocalLimits
!!
!! NAME 
!!
!!   gr_pfftGetLocalLimits
!!
!! SYNOPSIS
!!
!!   gr_pfftGetLocalLimits(integer(IN) :: axis1,
!!                         integer(IN) :: axis2)
!!
!! DESCRIPTION 
!!
!!  Find out the local index limits of the transformed
!!  data in a processor. This imformation is necessary to
!!  calculate the derivatives. Here the situation is complicated
!!  by the fact that after transformation, the order of physical
!!  dimension is <JAXIS,IAXIS> in 2D and <KAXIS,IAXIS,JAXIS> in 3D.
!!  So the ordering of dimensions in the array is different from the
!!  physical dimension order. The two arguments passed to this 
!!  function provide the mapping.
!!  
!! ARGUMENTS
!!
!!   axis1  - the first array dimension of the transformed data
!!   axis2  - the corresponding physical data dimension
!!
!! SIDE EFFECTS
!!
!!   Modifies pfft_localLimits(:,axis1).  This is how this subroutine
!!   communicates its results.
!!
!! NOTES
!!
!!   Indices in pfft_localLimit are 0-based.
!!
!!***

subroutine gr_pfftGetLocalLimits(axis1,axis2)

#include "constants.h"
#include "Pfft.h"

  use gr_pfftData, ONLY : pfft_me,pfft_globalLen,pfft_localLimits,pfft_outLen,&
       pfft_transformType, pfft_pclBaseDatType
  
  implicit none

  integer,intent(IN) :: axis1, axis2

  integer :: len,s,e


  len=pfft_globalLen(axis2)
  if(pfft_transformType(axis2)==PFFT_REAL2C .or. &
     pfft_transformType(axis2)==PFFT_REAL2C_STUFF .or. &
     pfft_transformType(axis2)==PFFT_REAL2C_EXTEND) then
     len=len/2
     if(pfft_transformType(axis2)==PFFT_REAL2C_EXTEND) then
        len=len+1
     end if
  end if
  
  s=pfft_me(axis1)*pfft_outLen(axis1)
  if(s>len) then
     pfft_localLimits(LOW:HIGH,axis1)= (/0,-1/)
  else
     pfft_localLimits(LOW,axis1)=s
     e=s+pfft_outLen(axis1)-1
     if(e>=len)e=len-1
     pfft_localLimits(HIGH,axis1)=e
  end if
!!$  print*,'pfft_globalLen(',axis2,')=',pfft_globalLen(axis2)
!!$  print*,'pfft_outLen(',axis1,')=',pfft_outLen(axis1)
!!$  print*,'pfft_LocalLimits(:,',axis1,')=',pfft_localLimits(:,axis1)
end subroutine gr_pfftGetLocalLimits
