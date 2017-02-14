!!****h* source/physics/sourceTerms/Heatexchange/Heatexchange_interface
!!
!! NAME
!!   Heatexchange_interface
!!
!! SYNOPSIS
!!   use Heatexchange_interface
!!
!! DESCRIPTION
!! This is the header file for the Heatexchange unit that defines its
!! public interfaces.
!!***
Module Heatexchange_interface
#include "constants.h"

  interface
    subroutine Heatexchange_computeDt (blockID,  &
                                  blkLimits,blkLimitsGC,        &
                                  solnData,   &
                                  dt_heatXchg, dt_minloc )
     

      integer, intent(IN) :: blockID 
      integer, intent(IN),dimension(2,MDIM)::blkLimits,blkLimitsGC
      real,INTENT(INOUT)    :: dt_heatXchg
      integer,INTENT(INOUT)    :: dt_minloc(5)
      real, pointer, dimension(:,:,:,:) :: solnData

    end subroutine Heatexchange_computeDt
  end interface

  interface
    subroutine Heatexchange_finalize()
    end subroutine Heatexchange_finalize
  end interface

  interface
    subroutine Heatexchange_init(restart)
      
      logical, intent(IN) :: restart
    end subroutine Heatexchange_init
  end interface

  interface
    subroutine Heatexchange(blockCount, blockList, dt )
      implicit none
      integer, INTENT(in)                        :: blockCount
      integer, INTENT(in), DIMENSION(blockCount)  :: blockList
      real,    INTENT(in)                        :: dt
    end subroutine Heatexchange
  end interface


end Module Heatexchange_interface
