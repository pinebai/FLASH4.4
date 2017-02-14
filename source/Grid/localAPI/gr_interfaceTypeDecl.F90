!!****ih* source/Grid/localAPI/gr_interfaceTypeDecl
!!
!! NAME
!!
!!  gr_interfaceTypeDecl
!!
!! SYNOPSIS
!!
!!  use gr_interfaceTypeDecl
!!
!! DESCRIPTION
!!
!!  Contains derived data type declarations used by gr_interface.  
!!  We need to place the derived data types in a module in order that 
!!  different subroutines have access to the same type declaration.
!!
!!***

module gr_interfaceTypeDecl

#include "constants.h"
#include "Flash.h"

  implicit none
  type BlockRegion_t
     integer :: numNegh  !No. neighbors at this guard cell region. 0 = external boundary.
     integer, dimension(BLKNO:TYPENO, 2**(NDIM-1)) :: details
  end type BlockRegion_t

  type AllBlockRegions_t
     type (BlockRegion_t), dimension(1:3, 1:1+2*K2D, 1:1+2*K3D) :: regionInfo
  end type AllBlockRegions_t

  type gr_solversDbgContext_t
     integer :: component
     integer :: group
     integer :: libErrCode
     integer :: flashErrCode    ! 1 for ERROR, 2 for INFO
     integer :: retriable       ! 0 for NO, 1 for YES
  end type gr_solversDbgContext_t

end module gr_interfaceTypeDecl
