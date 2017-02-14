!!****if* source/Grid/GridSolvers/Pfft/MeshReconfiguration/PtToPt/gr_pfftNodeObject
!!
!! NAME
!!  gr_pfftNodeObject
!!
!! SYNOPSIS
!!
!!  use gr_pfftNodeObject
!!
!! DESCRIPTION
!!
!!   Derived data types and other data structures for mapping PM to  
!!   Pencil UG 
!!   
!!***
module gr_pfftNodeObject
#include "constants.h"
  implicit none
  include "Flash_mpi.h"

  !We will create an MPI derived datatype representing node_metadata.
  !The tagID is required because a single processor can send
  !more than 1 MPI message to a particular processor.
  type node_metadata
     integer :: flashProcID, flashBlockID, pfftProcID, tagID
     integer, dimension(1:MDIM) :: flashStartPos, flashEndPos, &
          pfftStartPos, pfftEndPos
  end type node_metadata


  !This is our node object.
  type node
     type(node_metadata) :: metadata
     type(node), pointer :: next, prev

     !Buffer containing grid data (Must be a pointer as in a type).
     real, pointer, dimension(:) :: buf

     !MPI status & request objects for this node.
     integer, dimension(MPI_STATUS_SIZE) :: status
     integer :: request

     !This is a debugging mode which was created when 
     !problems were encountered using MPI derived datatypes on BG/P.
#ifdef MANUALLY_PACK_MPI_MESSAGE
     integer, dimension(4+(4*MDIM)) :: packedData
#endif
  end type node

contains 

  !The user must provide: create_node, destroy_node and print_node.
  subroutine create_node(item)
    use Driver_interface, ONLY : Driver_abortFlash
    implicit none
    type(node), pointer :: item
    integer :: err    
    allocate(item, STAT=err)
    if (err /= 0) then
       call Driver_abortFlash &
            ("[create_node]: No more heap memory")
    end if
    nullify(item % next, item % prev)
  end subroutine create_node

  !We have a custom deallocation subroutine because our nodes 
  !contain lots of pieces of data which must be deallocated appropiately.
  subroutine destroy_node(item)
    implicit none
    type(node), pointer  :: item    
    deallocate(item % buf)
    deallocate(item)
    nullify(item)
  end subroutine destroy_node

  subroutine print_node(unitNumber, item)    
    implicit none
    integer, intent(IN) :: unitNumber
    type(node), pointer :: item
    write(unitNumber,*) "[print_node]: Not yet in place"
  end subroutine print_node
end module gr_pfftNodeObject
