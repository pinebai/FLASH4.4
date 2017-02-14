!!****if* source/Grid/GridMain/Grid_getRowData
!!
!!
!! NAME
!!  Grid_getRowData
!!
!! SYNOPSIS
!!
!!  Grid_getRowData(integer(IN) :: blockID,
!!                  integer(IN) :: gridDataStruct,
!!                  integer(IN) :: structIndex,
!!                  integer(IN) :: beginCount, 
!!                  integer(IN) :: row,
!!                  integer(IN) :: startingPos(MDIM),
!!                  real(OUT)   :: datablock(dataSize),
!!                  integer(IN) :: dataSize)
!!  
!! DESCRIPTION 
!!  
!!  This routine allows the user to get an entire row or any contiguous
!!  part of a row of data depending on the arguments passed. The user is also
!!  allowed to specify if index counting should begin at the exterior edge
!!  of the block, (that is including guardcells)
!!  or the interior edge of the block. 
!!
!!  The data could be either a single variable in one of the grid data 
!!  structures such as cell centered/face centered/scratch, or it could
!!  be quantities derived from the size of the cell, such as cell volume/
!!  face area for all the cells in the specified portion of the block
!!
!!  
!! ARGUMENTS 
!!
!!  blockID : the local blockid
!!
!!  gridDataStruct : integer value specifying the type of data desired.
!!             Valid options are either one of the Grid data structures,
!!             or one of the derived quantities such as the cell volume
!!             or the cell area.
!!
!!             The options are defined in constants.h and they are :
!!                   CENTER cell centered variables
!!                   FACEX  face centered variable on faces along IAXIS
!!                   FACEY  face centered variable on faces along JAXIS
!!                   FACEZ  face centered variable on faces along IAXIS
!!                   WORK   single, cell centered variable, valid only
!!                          for paramesh
!!                   SCRATCH scratch space that can fit cell and face centered variables
!!                   SCRATCH_CTR scratch space for cell centered variables
!!                   SCRATCH_FACEX scratch space facex variables
!!                   SCRATCH_FACEY scratch space facey variables
!!                   SCRATCH_FACEZ scratch space facez variables
!!                   CELL_VOLUME volumes of specified cells 
!!                   CELL_FACEAREA face area for the specified cells on 
!!                                 specified face 
!!                   
!!  structIndex : integer value that provides index into the specified data
!!                structure. When gridDataStruct is one of the grid data structures
!!                structIndex translates to a specific variable in that 
!!                data structure, for example: DENS_VAR or PRES_VAR (define in
!!                Flash.h) if gridDataStruct = CENTER.
!!
!!                For gridDataStruct=CELL_VOLUME, it has no meaning
!!
!!                For gridDataStruct = CELL_FACEAREA, it can take one of 
!!                (ILO_FACE,IHI_FACE,JLO_FACE,JHI_FACE,KLO_FACE,KHI_FACE)
!!                defined in constants.h
!!  
!!
!!  beginCount : tells the routine where to start indice counting.  beginCount can
!!               be set to INTERIOR or EXTERIOR.  If INTERIOR is specified
!!               guardcell indices are not included and index 1 is the first interior cell. 
!!               If EXTERIOR is specified
!!               the first index, 1, is the left most guardcell.  See example
!!               below for more explanation.  (For most of the FLASH architecture code,
!!               we use EXTERIOR.  Some physics routines, however, find it helpful 
!!               only to work on the internal parts of the blocks (without
!!               guardcells) and wish to keep loop indicies  
!!               going from 1 to NXB without having to worry about finding 
!!               the correct offset for the number of guardcells.) 
!!               (INTERIOR and EXTERIOR are defined in constants.h)
!!
!!  row :      specifies the integer index coordinates of the cells
!!             The options are IAXIS, JAXIS or KAXIS defined in constants.h
!!
!!  startingPos(MDIM):
!!           specifies the starting position in each dimension of row
!!           row of data being fetched.
!!   
!!           startingPos(1) = i
!!           startingPos(2) = j
!!           startingPos(3) = k
!!
!!           For a 2d problem startingPos(3) is ignored
!!           For a 1d problem startingPos(3) and startingPos(2) are ignored
!!
!!
!!  datablock : a real array containing the data
!!              The dimensions for datablock are 
!!              datablock(datasize) Various compilers require the dimensions of
!!              datablock to be specified explicitly.  They are defined by 
!!              the next argument "dataSize".  
!!
!!
!!  dataSize : an integer specifying the dimensions for datablock
!!          
!!          dataSize holds the number of cells in the row to return
!!          
!!
!! EXAMPLE  
!!  
!!    Here is a 3d block example for getting an entire row of data, 
!!    including guardcells for each block on a local processor.  
!!    We will get the entire row in the 
!!    i direction at position j=5 and k=7.  beginCount is set to EXTERIOR
!!    meaning that the first cell of the block including guardcells is 
!!    index = 1.  If in this example, NGUARD = 4 then position 5 is the
!!    first interior cell in a dimension.
!! 
!!    (Can't really draw 3d for this example, but this
!!     picture is meant to help show where 
!!     counting begins when "beginCount" is set to EXTERIOR)
!!    
!!     j
!!    
!!     1         - - - - - - - -  
!!     2         - - - - - - - - 
!!     3         - - - - - - - - 
!!     4         - - - - - - - - 
!!     5 - - - -|-|-|-|-|-|-|-|-|- - - -
!!     6 - - - -|-|-|-|-|-|-|-|-|- - - -
!!     7 - - - -|-|-|-|-|-|-|-|-|- - - -
!!     8 - - - -|-|-|-|-|-|-|-|-|- - - -
!!     9 - - - -|-|-|-|-|-|-|-|-|- - - -
!!    10 - - - -|-|-|-|-|-|-|-|-|- - - -
!!    11 - - - -|-|-|-|-|-|-|-|-|- - - -
!!    12 - - - -|-|-|-|-|-|-|-|-|- - - -
!!    13         - - - - - - - - 
!!    14         - - - - - - - - 
!!    15         - - - - - - - - 
!!    16         - - - - - - - - 
!! i     1 2 3 4 5 6 7 8 9 10111213141516 
!! 
!!    #include "Flash.h"
!!    #include "constants.h"
!!
!!    ...
!!       
!!      integer ::    startingPos(MDIM)
!!      integer ::    dataSize
!!      integer ::    blockID
!!      real    ::    dataBlock(:)
!!       
!!          startingPos(1) = 1    
!!          startingPos(2) = 5
!!          startingPos(3) = 7
!!
!!          dataSize = blkLimitsGC(IAXIS) !This is equivalent to NXB + 2*NGUARD
!!                                        ! in FIXEDBLOCKSIZE mode
!!
!!          allocate(datablock(dataSize))
!!
!!          do blockID = 1, localNumBlocks
!!  
!!             call Grid_getRowData(blockID, CENTER, DENS_VAR, EXTERIOR, IAXIS, &
!!                               startingPos, dataBlock, dataSize)
!!  
!!          end do
!!
!!
!!  
!!    In this 2d block example we will get face area for x direction
!!    lower face of part of a row 
!!    for each block on a local processor.
!!    beginCount is set to INTERIOR, meaning that all the startPos indices
!!    will start where index 1 is the first interior cell of the block.
!!    In this example we will get just 3 cells in a row in the i direction
!!    beginning at position i= 4 and j position= 5
!!
!!    (hard to draw, but this is the idea, stars (*) are the cells to return
!!     notice the where indice counting starts when beginCount is set to INTERIOR)
!!            - - - - - - - - 
!!            - - - - - - - - 
!!            - - - - - - - - 
!!     j      - - - - - - - - 
!!     1 ----|-|-|-|-|-|-|-|-|----
!!     2 ----|-|-|-|-|-|-|-|-|----
!!     3 ----|-|-|-|-|-|-|-|-|----
!!     4 ----|-|-|-|-|-|-|-|-|----
!!     5 ----|-|-|-|*|*|*|-|-|----
!!     6 ----|-|-|-|-|-|-|-|-|----
!!     7 ----|-|-|-|-|-|-|-|-|----
!!     8 ----|-|-|-|-|-|-|-|-|----
!!            - - - - - - - - 
!!            - - - - - - - - 
!!            - - - - - - - - 
!!            - - - - - - - - 
!!         i  1-2-3-4 5-6-7-8 
!!
!!
!! 
!!    #include "Flash.h"
!!    #include "constants.h"
!!
!!    ...
!!       
!!      integer ::    startingPos(MDIM)
!!      integer ::    dataSize
!!      integer ::    blockID
!!      real    ::    dataBlock(:)
!!       
!!          startingPos(1) = 4    
!!          startingPos(2) = 5
!!          startingPos(3) = 1 !this is ignored since only 2d
!!
!!          dataSize = 3 !just getting 3 cells
!!
!!          allocate(datablock(dataSize))
!!
!!          do blockID = 1, localNumBlocks
!!  
!!             call Grid_getRowData(blockID, CELL_FACEAREA, ILO_FACE, INTERIOR, IAXIS, &
!!                               startingPos, dataBlock, dataSize)
!!  
!!          end do
!!
!!
!!
!!
!!
!!
!! 
!!
!!***


!!REORDER(5): unk, facevar[xyz]
!!REORDER(4): solnData

#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif

subroutine Grid_getRowData(blockID, gridDataStruct, structIndex, beginCount, &
     row, startingPos, datablock, dataSize)

  use Grid_data, ONLY : gr_iguard, gr_jguard, gr_kguard
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface, ONLY : Grid_getBlkPtr,Grid_releaseBlkPtr
  use gr_interface, ONLY : gr_getInteriorBlkPtr, gr_releaseInteriorBlkPtr


  implicit none

#include "constants.h"
#include "Flash.h"

  integer, intent(in) :: blockID, structIndex, beginCount, row, gridDataStruct
  integer, dimension(MDIM), intent(in) :: startingPos
  integer, intent(in) :: dataSize
  real, dimension(datasize),intent(out) :: datablock
  real,allocatable,dimension(:,:,:) :: cellvalues
  real, pointer, dimension(:,:,:,:) :: solnData
  integer :: xb,xe,yb,ye,zb,ze

  integer :: i, j ,k
  integer,dimension(MDIM)::begOffset,dataLen
  integer :: imax, jmax, kmax

  logical :: isget
  logical :: getIntPtr

#ifdef DEBUG_GRID
  isget = .true.
  call gr_checkDataType(blockID,gridDataStruct,imax,jmax,kmax,isget)

  !verify requested row is available given number of dims in problem
  if((row==KAXIS) .and. (NDIM < 3)) then
     print *, "Error: Grid_getRowData "
     call Driver_abortFlash("you have requested to return the KAXIS in a 1d or 2d problem")
  end if
  

 


  !verify requested row is available given number of dims in problem
  if((row==JAXIS) .and. (NDIM < 2)) then
     print *, "Error: Grid_getRowData"
     call Driver_abortFlash("you have requested to return the JAXIS in a 1d problem")
  end if
     

  !verify we have a valid blockid
  if((blockid<1).or.(blockid>MAXBLOCKS)) then
     print*,"Error: Grid_getRowData : invalid blockid "
     call Driver_abortFlash("Get_getRowData : invalid blockid ")
  end if
  
  !verify beginCount is set to a valid value
  if((beginCount /= INTERIOR) .and. (beginCount /= EXTERIOR)) then
     print *, "Error: Grid_getRowData: beginCount set to improper value"
     print *, "beginCount must = INTERIOR or EXTERIOR (defined in constants.h)"
     call Driver_abortFlash("beginCount must = INTERIOR or EXTERIOR (defined in constants.h)")
  end if
     

  !verify that there is enough space in datablock
  if ((row==IAXIS) .and. (dataSize > imax)) then
     print *, "Error: Grid_getRowData: dataSize too big"
     print *,"You are requesting more cells than block has in a dimension"
     call Driver_abortFlash("Grid_getRowData: dataSize too big")
  end if

  if ((row==JAXIS) .and. (dataSize > jmax)) then
     print *, "Error: Grid_getRowData: dataSize too big"
     print *,"You are requesting more cells than block has in a dimension"
     call Driver_abortFlash("Grid_getRowData: dataSize too big")
  end if

  if ((row==KAXIS) .and. (dataSize > kmax)) then
     print *, "Error: Grid_getRowData: dataSize too big"
     print *,"You are requesting more cells than block has in a dimension"
     call Driver_abortFlash("Grid_getRowData: dataSize too big")
  end if



  !verify that there is enough space in datablock
  if ((dataSize < 1)) then 
     
     print *, "Error: Grid_getRowData: dataSize too small"
     print *,"You are requesting more < 1 cell in a dimension of block, 1 is the min"
     call Driver_abortFlash("Grid_getRowData: dataSize too small")
  end if
  



  !verify that indicies aren't too big or too small for the block
  if(beginCount == EXTERIOR) then
    
     if (startingPos(1) > imax) then
        call Driver_abortFlash("Grid_getRowData startingPos(1) index larger than block")
     end if

     if ((NDIM > 1) .and. (startingPos(2) > jmax)) then
        call Driver_abortFlash("Grid_getRowData startingPos(2) index larger than block")
     end if
    
     if ((NDIM > 2) .and. (startingPos(3) > kmax)) then
        call Driver_abortFlash("Grid_getRowData startingPos(3) index larger than block")
     end if
    
     if (startingPos(1) < 1) then
        call Driver_abortFlash("Grid_getRowData startingPos(1) index smaller than 1")
     end if

     if ((NDIM > 1) .and. (startingPos(2) > jmax)) then
        call Driver_abortFlash("Grid_getRowData startingPos(2) index smaller than 1")
     end if
    
     if ((NDIM > 2) .and. (startingPos(3) > kmax)) then
        call Driver_abortFlash("Grid_getRowData startingPos(3) index smaller than 1")
     end if
        
  else !beginCount == INTERIOR

     if ((startingPos(1) + gr_iguard -1) > imax) then
        call Driver_abortFlash("Grid_getRowData startingPos(1) index larger than block")
     end if

     if ((NDIM > 1) .and. ((startingPos(2) + gr_jguard -1) > jmax)) then
        call Driver_abortFlash("Grid_getRowData startingPos(2) index larger than block")
     end if
    
     if ((NDIM > 2) .and. ((startingPos(3) + gr_kguard -1) > kmax)) then
        call Driver_abortFlash("Grid_getRowData startingPos(3) index larger than block")
     end if
    
     if (startingPos(1) < 1) then
        call Driver_abortFlash("Grid_getRowData startingPos(1) index smaller than 1")
     end if

     if ((NDIM > 1) .and. (startingPos(2) < 1)) then
        call Driver_abortFlash("Grid_getRowData startingPos(2) index smaller than 1")
     end if
    
     if ((NDIM > 2) .and. (startingPos(3) <  1)) then
        call Driver_abortFlash("Grid_getRowData startingPos(3) index smaller than 1")
     end if

  end if
  
  

  !more verification of indicies
  !check size and starting pos
  if(beginCount == EXTERIOR) then
     if(row == IAXIS) then
        if ((startingPos(IAXIS) + dataSize -1) > imax) then
           print *, "Error: Grid_getRowData"
           call Driver_abortFlash("Grid_getRowData indicies too large")
        end if
     end if

     if(row == JAXIS) then
        if ((startingPos(JAXIS) + dataSize -1) > jmax) then
           print *, "Error: Grid_getRowData"
           call Driver_abortFlash("Grid_getRowData indicies too large")
        end if
     end if

     if(row == KAXIS) then
        if ((startingPos(KAXIS) + dataSize -1) > kmax) then
           print *, "Error: Grid_getRowData"
           call Driver_abortFlash("Grid_getRowData indicies too large")
        end if
     end if
  else if(beginCount == INTERIOR) then
     if(row == IAXIS) then
        if ((startingPos(IAXIS) + dataSize + gr_iguard -1) > imax) then
           print *, "Error: Grid_getRowData"
           call Driver_abortFlash("Grid_getRowData indicies too large")
        end if
     end if

     if(row == JAXIS) then
        if ((startingPos(JAXIS) + dataSize + gr_jguard -1) > jmax) then
           print *, "Error: Grid_getRowData"
           call Driver_abortFlash("Grid_getRowData indicies too large")
        end if
     end if

     if(row == KAXIS) then
        if ((startingPos(KAXIS) + dataSize + gr_kguard -1) > kmax) then
           print *, "Error: Grid_getRowData"
           call Driver_abortFlash("Grid_getRowData indicies too large")
        end if
     end if
  end if
#endif
  dataLen=0
  dataLen(row)=dataSize
  call gr_getDataOffsets(blockID,gridDataStruct,startingPos,dataLen,beginCount,begOffset,getIntPtr)

  i=1
  j=1
  k=1
  

     !set the starting and ending position in x dir
  i = startingPos(IAXIS) + begOffset(IAXIS)
  if(NDIM > 1) j = startingPos(JAXIS) + begOffset(JAXIS)
  if(NDIM > 2) k = startingPos(KAXIS) + begOffset(KAXIS)
  if((gridDataStruct==CELL_VOLUME).or.(gridDataStruct==CELL_FACEAREA)) then
     xb=i;yb=j;zb=k;xe=xb;ye=yb;ze=zb
     if(row==IAXIS) xe=xb+datasize-1
     if(row==JAXIS) ye=yb+datasize-1
     if(row==KAXIS) ze=zb+datasize-1
     allocate(cellvalues(xb:xe,yb:ye,zb:ze))
     if(gridDataStruct==CELL_VOLUME) then
        call gr_getCellVol(xb,xe,yb,ye,zb,ze,blockID,cellvalues)
     else
        call gr_getCellFaceArea(xb,xe,yb,ye,zb,ze,&
             structIndex,blockID,cellvalues)
     end if
     if(row==IAXIS) datablock(:)=cellvalues(xb:xe,yb,zb)
     if(row==JAXIS) datablock(:)=cellvalues(xb,yb:ye,zb)
     if(row==KAXIS) datablock(:)=cellvalues(xb,yb,zb:ze)
     deallocate(cellvalues)
  elseif(getIntPtr) then
     call gr_getInteriorBlkPtr(blockID,solnData,gridDataStruct)
     if(row==IAXIS)datablock(:) = solnData(structIndex,i:i+datasize-1,j,k)
     if(row==JAXIS)datablock(:) = solnData(structIndex,i,j:j+datasize-1,k)
     if(row==KAXIS)datablock(:) = solnData(structIndex,i,j,k:k+datasize-1)
     call gr_releaseInteriorBlkPtr(blockID,solnData,gridDataStruct)
  else
     call Grid_getBlkPtr(blockID,solnData,gridDataStruct)
!!$     if(gridDataStruct==SCRATCH) then
!!$        if(row==IAXIS)datablock(:) = solnData(i:i+datasize-1,j,k,structIndex)
!!$        if(row==JAXIS)datablock(:) = solnData(i,j:j+datasize-1,k,structIndex)
!!$        if(row==KAXIS)datablock(:) = solnData(i,j,k:k+datasize-1,structIndex)
!!$     else
!!$     end if
     if(row==IAXIS)datablock(:) = solnData(structIndex,i:i+datasize-1,j,k)
     if(row==JAXIS)datablock(:) = solnData(structIndex,i,j:j+datasize-1,k)
     if(row==KAXIS)datablock(:) = solnData(structIndex,i,j,k:k+datasize-1)
     call Grid_releaseBlkPtr(blockID,solnData,gridDataStruct)

  end if
  return
end subroutine Grid_getRowData
