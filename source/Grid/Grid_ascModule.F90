!!****f* source/Grid/Grid_ascModule
!!
!! NAME
!!  Grid_ascModule
!!
!! SYNOPSIS
!!
!!  use Grid_ascModule
!!
!! DESCRIPTION 
!!  
!!  Module with data and some accessor routines for allocatable scratches.
!!
!!   
!!***

Module Grid_ascModule

  implicit none


contains
  subroutine Grid_ascStart()
  end subroutine Grid_ascStart

  subroutine Grid_ascAllocMem(gds,var1,nvars, nGuardCtr,nGuardFaceN,nGuardFaceT, &
                                leafBlocks, arrayRank, highSize)

    integer, intent(IN) :: gds,var1,nvars
    integer, OPTIONAL, intent(IN) :: nGuardCtr,nGuardFaceN,nGuardFaceT
    integer, OPTIONAL, intent(IN), dimension(:) :: leafBlocks
    integer, OPTIONAL, intent(IN) :: arrayRank
    integer, OPTIONAL, intent(IN) :: highSize
  end subroutine Grid_ascAllocMem

  subroutine Grid_ascDeallocMem(gds, arrayRank)

    integer, OPTIONAL, intent(IN) :: gds
    integer, OPTIONAL, intent(IN) :: arrayRank

  end subroutine Grid_ascDeallocMem

end Module Grid_ascModule
