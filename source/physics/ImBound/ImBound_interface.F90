!****h* source/physics/ImBound/Imbound_interface
!!
!! NAME
!!
!!  Imbound_interface
!!
!! SYNOPSIS
!!
!!  use Imbound_interface
!!
!! DESCRIPTION
!!
!! This is the header file for the Immersed boundary (IB)
!! unit that defines its public interfaces.
!!
!!***

Module ImBound_interface

  implicit none

#include "Flash.h"


  interface !ImBound_init
     subroutine ImBound_init(restart)
       implicit none
       logical, INTENT(IN) :: restart
     end subroutine ImBound_init
  end interface

  interface !ImBound 
     subroutine ImBound(blockCount, blockList, dt, forcflag)
       implicit none
       !! ---- Argument List ----------------------------------
       integer, INTENT(IN) :: blockCount
       integer, INTENT(IN), dimension(MAXBLOCKS) :: blockList
       real, INTENT(IN) :: dt
       integer, INTENT(IN) :: forcflag
       !! -----------------------------------------------------
     end subroutine ImBound
  end interface 

 interface  !ImBound_finalize
   subroutine ImBound_finalize()
   implicit none
   end subroutine ImBound_finalize
 end interface 

  interface
     subroutine ImBound_setData()
     end subroutine ImBound_setData
  end interface

 
end Module ImBound_interface
