!!****h* source/physics/sourceTerms/Ionize/Ionize_interface
!!
!!   Ionize's public interfaces
!!***

Module Ionize_interface
#include "constants.h"
#include "Flash.h"
#include "Ionize.h"

  interface
    subroutine Ionize(blockCount, blockList, dt, time)      
      implicit none
      integer, intent(IN) :: blockCount
      integer, dimension(blockCount), intent(IN) :: blockList
      real, intent(IN) :: dt, time
    end subroutine Ionize
  end interface

  interface
    subroutine Ionize_init()
      implicit none
      
    end subroutine Ionize_init
  end interface

  interface
    subroutine Ionize_finalize()    
      implicit none
    end subroutine Ionize_finalize
  end interface

  interface
     subroutine Ionize_equil(tx,nel,nion,delem)
       implicit none
        real, intent(INOUT) :: tx
        integer, intent(IN) :: nel
        integer, intent(OUT) :: nion
        real, dimension(ION_NIMAX), intent(OUT) :: delem
     end subroutine Ionize_equil
  end interface
end Module
