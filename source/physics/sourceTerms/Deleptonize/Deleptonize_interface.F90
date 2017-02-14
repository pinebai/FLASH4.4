!!****h* source/physics/sourceTerms/Deleptonize/Deleptonize_interface
!!
!! This is the header file for the heat module that defines its
!! public interfaces.
!!***
Module Deleptonize_interface
#include "constants.h"
#include "Flash.h"

  interface Deleptonize
     subroutine Deleptonize (blockCount,blockList,dt,time)
       integer,intent(IN) :: blockCount
       integer,dimension(blockCount),intent(IN)::blockList
       real,intent(IN) :: dt,time
     end subroutine Deleptonize
  end interface

  interface Deleptonize_init
     subroutine Deleptonize_init()
     end subroutine Deleptonize_init
  end interface


  interface Deleptonize_finalize
     subroutine Deleptonize_finalize ()
     end subroutine Deleptonize_finalize
  end interface

  interface Deleptonize_getBounce
     subroutine Deleptonize_getBounce(postBounce,bounceTime)
       logical, intent(OUT) :: postBounce
       real, intent(OUT) :: bounceTime
     end subroutine Deleptonize_getBounce
  end interface

end Module Deleptonize_interface
