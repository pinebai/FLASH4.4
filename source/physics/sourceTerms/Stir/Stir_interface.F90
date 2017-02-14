!!****h* source/physics/sourceTerms/Stir/Stir_interface
!!
!! This is the header file for the Stir module that defines its
!! public interfaces.
!!***
Module Stir_interface

#include "constants.h"

  interface Stir
    subroutine Stir(blockCount,blockList,dt)
      integer, intent(IN) :: blockCount
      integer,dimension(blockCount), intent(IN) :: blockList
      real,intent(IN) :: dt
    end subroutine Stir
  end interface

  interface Stir_finalize
    subroutine Stir_finalize()
    end subroutine Stir_finalize
  end interface

  interface Stir_init
    subroutine Stir_init(restart)
      
      logical, intent(IN) :: restart
    end subroutine Stir_init
  end interface

  interface Stir_computeDt

     subroutine Stir_computeDt(blockID,                 & 
                           blkLimits,blkLimitsGC,        &
                           solnData,                     &
                           dt_stir, dt_minloc)

       !! arguments
       integer, intent(IN)   :: blockID
       integer, intent(IN),dimension(2,MDIM)::blkLimits,blkLimitsGC
       real, pointer           :: solnData(:,:,:,:) 
       real, intent(INOUT)     :: dt_stir
       integer, intent(INOUT)  :: dt_minloc(5)
       
     end subroutine Stir_computeDt

  end interface
end Module Stir_interface


