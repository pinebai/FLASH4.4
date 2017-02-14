module ut_contourSurfaceInterface
interface
  subroutine ut_contourSurfaceAreaBlock(nlevels, isolevels, ctrData, blkLimits, blkIndex, areas)
     implicit none
#include "Flash.h"
#include "constants.h"

     integer,                     intent(IN)  :: nlevels
     real,dimension(nlevels),     intent(IN)  :: isolevels
     real,dimension(:,:,:),       intent(IN)  :: ctrData
     integer,dimension(HIGH,MDIM),intent(IN)  :: blkLimits
     integer,                     intent(IN)  :: blkIndex
     real,dimension(nlevels),     intent(OUT) :: areas

  end subroutine ut_contourSurfaceAreaBlock
end interface
end module ut_contourSurfaceInterface
