!!****if* source/Grid/GridMain/paramesh/paramesh4/gr_extendedGetCellCoords
!!
!! NAME
!!  gr_extendedGetCellCoords
!!
!! SYNOPSIS
!!
!!  gr_extendedGetCellCoords(integer(IN) :: axis,
!!                      integer(IN):: blockID, 
!!                      integer(IN):: pe, 
!!                      integer(IN):: edge, 
!!                      logical(IN):: guardcell, 
!!                      real(OUT)  :: coordinates(size),
!!                      integer(IN):: size)
!!  
!!  
!! DESCRIPTION
!!
!!    This subroutine is an accessor function that gets the coordinates of
!!    the cells in a given block. This is a generalized variant of
!!    Grid_getCellCoords, for use within the Grid unit only:
!!    it can be called for blocks local to the executing processor or,
!!    at least in the Paramesh3 and 4 implementations of the Grid unit,
!!    for remote blocks for which cached information is currently locally
!!    available.
!!
!!    The block about which information is requested is identified by
!!    a pair (blockID, pe). All other arguments are as for Grid_getCellCoords
!!    and are passed to Grid_getCellCoords if it is called by this routine.
!!
!!    Coordinates are retrieved one axis at a time, 
!!    meaning you can get the i, j, _or_ k coordinates with one call.  
!!    If you want all the coordinates, all axes, you
!!    need to call gr_extendedGetCellCoords 3 times, one for each axis.
!!    The code carries coordinates at cell centers as well as faces.
!!    It is possible to get coordinates for CENTER, only LEFT_EDGE,
!!    only RIGHT_EDGE or for all FACES along a dimension.
!!
!!
!!
!!
!! ARGUMENTS
!!            
!!   axis - specifies the integer index coordinates of the cells being retrieved.
!!          axis can have one of three different values, IAXIS, JAXIS or KAXIS 
!!          (defined in constants.h as 1,2 and 3)
!!
!!   blockID - integer block number
!!
!!   pe      - processor where block (or a cached copy of block info) resides
!!
!!   edge - integer value with one of four values, 
!!          LEFT_EDGE, RIGHT_EDGE, CENTER or FACES
!!          The edge argument specifies what side of the zone to get, 
!!          the CENTER point, the LEFT_EDGE  or the RIGHT_EDGE of the zone.
!!          FACES gets the left and right face of each cell, but since 
!!          two adjacent cells have a common face, there are only N+1
!!          unique values if N is the number of cells.
!!
!!   guardcell - logical input. If true coordinates for guardcells are returned
!!          along with the interior cells, if false, only the interior coordinates 
!!          are returned.
!!
!!          
!!   coordinates - The array holding the data returning the coordinate values
!!                 coordinates must be at least as big as "size" (see below)
!!           
!!   size - integer specifying the size of the coordinates array.
!!          if edge = CENTER/LEFT_EDGE/RIGHT_EDGE then
!!                If guardcell true then size =  interior cells + 2*guardcells
!!                otherwise size = number of interior cells
!!          If edge=FACES 
!!                If guardcell true then size =  interior cells + 2*guardcells+1
!!                otherwise size = number of interior cells+1
!!
!!               
!!
!!
!!  NOTES
!!   variables that start with "gr_" are variables of Grid unit scope
!!   and are stored in the fortran module Grid_data. Variables are not
!!   starting with gr_ are local variables or arguments passed to the 
!!   routine.
!!
!!  SEE ALSO
!!
!!   Grid_getCellCoords
!!
!!***

#ifdef DEBUG
#define DEBUG_GRID
#endif

subroutine gr_extendedGetCellCoords(axis, blockID, pe, edge, guardcell, coordinates, size)

  use Grid_data, ONLY : gr_meshMe, &
       gr_delta, gr_imin, gr_jmin, gr_kmin
  use tree, ONLY: bnd_box, laddress, strt_buffer, last_buffer, lnblocks, lrefine, lrefine_max
  use Grid_interface, ONLY : Grid_getCellCoords
  use Driver_interface, ONLY : Driver_abortFlash

#include "constants.h"
#include "Flash.h"

  implicit none

  integer, intent(in) :: axis,blockID,pe, edge
  integer, intent(in) :: size
  logical, intent(in) :: guardcell
  real,intent(out), dimension(size) :: coordinates

  integer :: bOffset,eOffset,calcSize,numGuard

  integer :: accessBlk, accessPE ! what we use here to access coordinate information
  integer :: iblk, cornerID, stride, indLeft, indRight
  logical :: lfound
  real :: half_Delta, bnd_box_left, domain_left, indCenter
  real :: coordDeviation
  integer :: j

    ! The following logic copied from Tests/amr_1blk_bcset.F90 in Paramesh4 distribution. - KW
    ! This is heavily dependent on Paramesh internals.
  if (pe.EQ.gr_meshMe) then
     accessBlk = blockID
     accessPE  = gr_meshMe
  else
     lfound = .false.
     iblk = strt_buffer
     do while(.not.lfound.and.iblk.le.last_buffer)
        if(blockID.eq.laddress(1,iblk).and.pe.eq.laddress(2,iblk) ) then
           accessBlk = iblk
           accessPE  = gr_meshMe
           lfound = .true.
        endif
        iblk = iblk+1
     enddo
     if(.not.lfound) then
        call Driver_abortFlash('Paramesh error: gr_extendedGetCellCoords: '// & 
             &      ' remote block is not in list received on this pe')
     endif
  endif

  ! Do some error checking here
  

  if(((NDIM==1).and.(axis/=IAXIS)).or.((NDIM==2).and.(axis==KAXIS))) then
     bOffset = 0
     eOffset = 0
     numGuard = 0
  else
     if(guardcell) then
        bOffset = 0
        eOffset = 2*NGUARD
        numGuard = NGUARD
     else
        boffset = NGUARD
        eoffset = NGUARD
        numGuard = 0
     end if
  end if
  if(axis == IAXIS) then
     calcSize = NXB+2*numGuard
  else if(axis == JAXIS) then
     calcSize = NYB+2*K2D*numGuard
  else if(axis == KAXIS) then
     calcSize = NZB+2*K3D*numGuard
  end if

#ifdef DEBUG_GRID
  if (edge==FACES) calcSize = calcSize+1
  print*,' get coordinates', axis, blockID, edge, guardcell,size
  if((blockID<1).or.(blockID>MAXBLOCKS)) then
     call Driver_abortFlash("gr_extendedGetCellCoords :invalid blockID ")
  end if
  if(.not.((edge==LEFT_EDGE).or.(edge==RIGHT_EDGE).or.&
       &(edge==CENTER) .or. (edge==FACES))) then
     call Driver_abortFlash("Get Coords : invalid edge specification, must be CENTER, &
          LEFT_EDGE, RIGHT_EDGE, or FACES")
  end if

!!!  This can be refined further to make it geometry specific

  if(.not.((axis==IAXIS).or.(axis==JAXIS).or.(axis==KAXIS))) then
     call Driver_abortFlash("Get Coords : invalid axis, must be IAXIS, JAXIS or KAXIS ")
  end if
  
  if(size < calcSize)then
     call Driver_abortFlash("Get Coords : size of output array too small")
  end if
#endif

  if (accessBlk.LE.lnblocks) then !if it is a local block on this PE...
     call Grid_getCellCoords(axis, blockID, edge, guardcell, coordinates, size)
#ifndef DEBUG_GRID
     ! Return results of this call directly if NOT debugging.
     ! Otherwise, the response from this call is used for debugging.
     return
#endif
  endif

  ! The following logic copied from gr_updateData. This may not be the most
  ! efficient way of doing it, but copying should ensure that we get bit-by-bit
  ! identical results for cell center coordinates. - KW

  if (axis==IAXIS) then
     domain_left = gr_imin
  else if (axis==JAXIS) then
     domain_left = gr_jmin
  else
     domain_left = gr_kmin
  end if
  ! the variable delta stores dx, dy, and dz values for
  ! each level of refinement given an initial domain size
  ! set in Grid_init
  half_Delta = gr_delta(axis,lrefine_max)*0.5

  bnd_box_left = bnd_box(1,axis,accessBlk)
  stride = 2**(lrefine_max - lrefine(accessBlk))

  !cornerID is a unique integer value
  !assigned to each block
  cornerID = (bnd_box_left-domain_left+half_Delta)/&
       gr_delta(axis,lrefine_max) + 1

  ! calculate the left most index, use this to calculate coords
  indLeft = cornerID-stride*NGUARD -1
  do j=1,boffset+size
     indCenter = indLeft + stride*0.5
     indRight = indLeft + stride

     if (j > boffset) then
#ifdef DEBUG_GRID
        if (accessBlk.LE.lnblocks) then ! compare if possible for debugging, to go away later - KW

           if (edge==CENTER) then
              coordDeviation = (coordinates(j-boffset) - (domain_left + gr_delta(axis,lrefine_max)*indCenter))
           else if (edge==LEFT_EDGE .OR. edge==FACES) then
              coordDeviation = (coordinates(j-boffset) - (domain_left + gr_delta(axis,lrefine_max)*indLeft))
           else if (edge==RIGHT_EDGE) then
              coordDeviation = (coordinates(j-boffset) - (domain_left + gr_delta(axis,lrefine_max)*indRight))
           end if

           if (abs(coordDeviation) .GE. spacing(coordinates(j-boffset))) then
#undef REAL_FORMAT
#define REAL_FORMAT "(1PG33.26)"
              print*,'pe,block,accblk,axis,j',pe,blockID,accessBlk,axis,j,'coords differ!!!', &
                   & (coordDeviation .GT. spacing(coordinates(j-boffset)))
              print '('//REAL_FORMAT//',a6,3('//REAL_FORMAT//',1x))',coordinates(j-boffset),' .NE. &
                   &',(domain_left + gr_delta(axis,lrefine_max)*indCenter),&
                   coordDeviation,spacing(coordinates(j-boffset))
              if (abs(coordDeviation) .GT. spacing(coordinates(j-boffset))) &
                   call Driver_abortFlash &
                   &('[gr_extendedGetCellCoords] Two ways to calculate cell coordinates give significantly different answers!')
           end if
        endif
#endif

        if (edge==CENTER) then
           coordinates(j-boffset) = domain_left + gr_delta(axis,lrefine_max)*indCenter
        else if (edge==LEFT_EDGE .OR. edge==FACES) then
           coordinates(j-boffset) = domain_left + gr_delta(axis,lrefine_max)*indLeft
        else if (edge==RIGHT_EDGE) then
           coordinates(j-boffset) = domain_left + gr_delta(axis,lrefine_max)*indRight
        end if
     end if

     indLeft = indRight
  end do

  return
end subroutine gr_extendedGetCellCoords





