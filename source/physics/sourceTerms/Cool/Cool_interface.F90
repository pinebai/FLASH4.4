!!****h* source/physics/sourceTerms/Cool/Cool_interface
!!
!! This is the header file for the cool module that defines its
!! public interfaces.
!!***
Module Cool_interface
#include "constants.h"
#include "Flash.h"
  
  interface 
     subroutine Cool_computeDt (block_no, &
          blkLimits,blkLimitsGC,  &
          solnData,   &
          dt_check, dt_minloc )
       
       integer, intent(IN) :: block_no
       integer, intent(IN),dimension(2,MDIM)::blkLimits,blkLimitsGC
       real,INTENT(INOUT)    :: dt_check
       integer,INTENT(INOUT)    :: dt_minloc(5)
       real, pointer, dimension(:,:,:,:) :: solnData
     end subroutine Cool_computeDt
  end interface
  
  interface 
     subroutine Cool(blockCount,blockList,dt, time)
       integer, intent(IN) :: blockCount
       integer,dimension(blockCount), intent(IN) :: blockList
       real,intent(IN) :: dt, time
     end subroutine Cool
  end interface
  
  interface 
     subroutine Cool_finalize()
     end subroutine Cool_finalize
  end interface
  
  interface 
     subroutine Cool_init()
       
     end subroutine Cool_init
  end interface

  interface
     subroutine Cool_unitTest( fileUnit, perfect)
         integer, intent(IN) ::  fileUnit
         logical, intent(INOUT) :: perfect
     end subroutine Cool_unitTest
  end interface

end Module Cool_interface


