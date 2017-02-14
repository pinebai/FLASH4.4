!!****if* source/Grid/GridMain/gr_getCellVol
!!
!! NAME
!!  gr_getCellVol
!!
!! SYNOPSIS
!!
!!  gr_getCellVol(integer(IN) :: xb,
!!                integer(IN) :: xe,
!!                integer(IN) :: yb,
!!                integer(IN) :: ye,
!!                integer(IN) :: zb,
!!                integer(IN) :: ze,
!!                integer(IN) :: blockID,
!!                real(xb:xe,yb:ye,zb:ze) (OUT) :: dataBlock)
!!  
!! DESCRIPTION 
!!  
!!  This routine calculate the volume of 
!!  the specified cells (bounded by xb:xe,yb:ye,zb:ze) of a given block.
!!
!! ARGUMENTS 
!!
!!  xb        : starting cell along IAXIS
!!  xe        : last cell along IAXIS
!!  yb        : starting cell along JAXIS
!!  ye        : last cell along JAXIS
!!  zb        : starting cell along KAXIS
!!  ze        : last cell along KAXIS
!!  blockID   : my block number
!!  dataBlock : storage for returning calculated values
!!
!!
!!***

subroutine gr_getCellVol(xb,xe,yb,ye,zb,ze,blockID,dataBlock)

  use Grid_data, ONLY : gr_geometry
  use Grid_interface, ONLY : Grid_getDeltas, Grid_getSingleCellVol

  implicit none
#include "Flash.h"
#include "constants.h"

  integer,intent(IN) :: xb,xe,yb,ye,zb,ze,blockID
  real,dimension(xb:xe,yb:ye,zb:ze),intent(OUT)::dataBlock
  integer,dimension(MDIM) :: point
  real,dimension(MDIM) :: del
  integer :: i, j, k

#ifdef DEBUG_GRID
  print*,xb,xe,yb,ye,zb,ze
#endif
  if (gr_geometry==CARTESIAN) then
     call Grid_getDeltas(blockID,del)
     dataBlock=del(IAXIS)
     do i = 2,NDIM
        dataBlock = dataBlock*del(i)
     end do
#ifdef DEBUG_GRID
     print*,'in here it is ',maxval(dataBlock),del
#endif
  else
     do k=zb,ze
        point(KAXIS) = k
        do j=yb,ye
           point(JAXIS) = j
           do i=xb,xe
              point(IAXIS) = i
              call Grid_getSingleCellVol(blockID, EXTERIOR, point, &
                   dataBlock(i,j,k))
           end do
        end do
     end do
  end if
end subroutine gr_getCellVol
