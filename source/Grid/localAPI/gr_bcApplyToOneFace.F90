!!****if* source/Grid/localAPI/gr_bcApplyToOneFace
!!
!! NAME
!!  gr_bcApplyToOneFace
!!
!! SYNOPSIS
!!
!!  gr_bcApplyToOneFace(integer(IN) :: axis,
!!                      integer(IN) :: bcType,
!!                      integer(IN) :: gridDataStruct,
!!                      integer(IN) :: varCount,
!!                      integer(IN) :: regionType(MDIM)
!!                      integer(IN) :: blkLimits(LOW:HIGH,MDIM)
!!                      integer(IN) :: blkLimitsGC(LOW:HIGH,MDIM)
!!                      integer(IN) :: blockHandle,
!!                      integer(IN) :: idest)
!!  
!! DESCRIPTION 
!!
!!  
!! 
!! ARGUMENTS
!!  
!!    axis           - the direction for applying BC, one of IAXIS, JAXIS, or KAXIS
!!    bcType         - the type of boundary condition
!!    gridDataStruct - In PM3 and PM4 it can have values (CENTER,FACEX,FACEY,FACEZ,WORK),
!!                     in UG (CENTER,FACEX,FACEY,FACEZ), and in PM2 (CENTER,WORK).
!!    varCount       - the number of variable in the data structure specified in gridDataStruct
!!    regionType     - The part of the block that is involved in the boundary condition. This integer
!!                     array can have values (LEFT_EDGE, CENTER, RIGHT_EDGE, WHOLE_VEC and NO_VEC)
!!                     for each of the three dimensions of the physical data.
!!                     LEFT_EDGE implies guard cells along lower face. If this value is specified
!!                     for the axis along which the BC are being applies, that means that we 
!!                     are applying BC to the lowerface, if it is one of the other axes, then
!!                     the BC are being applied to one of the corners. Same is true of RIGHT_EDGE
!!                     except that implies the upperface. CENTER, WHOLE_VEC and NO_VEC values
!!                     are valid only for dimensions other than the one specified in "axis". 
!!                     CENTER implies only the interior cells, WHOLE_VEC implies all cells in 
!!                     the block and NO_VEC implies that the correspoding dimension is not
!!                     a part of the region. Normally this value is most likely to be used
!!                     along KAXIS in a 2D problems, and JAXIS and KAXIS in a 1D problem
!!    blkLimits      - the endpoints of the block without including the guardcells
!!    blkLimitsGC    - the endpoints of the block without including the guardcells
!!    blockHandle        - local block number
!!    idest         - this is useful in paramesh, when the boundary conditions
!!                   are being applied to the workspace array WORK
!!
!! NOTES
!!  A specific direction is required in axis - no ALLDIR at this point.
!!
!!***

subroutine gr_bcApplyToOneFace(axis,bcType,gridDataStruct,varCount,&
     regionType,blkLimits,blkLimitsGC,blockHandle,idest)

  implicit none
#include "constants.h"
  
  integer, intent(in) :: axis,bcType,gridDataStruct,varCount,blockHandle,idest
  integer,dimension(MDIM),intent(IN) :: regionType
  integer,dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimitsGC,blkLimits

end subroutine gr_bcApplyToOneFace
