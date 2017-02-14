!!****f* source/Grid/Grid_getPlaneData
!!
!! NAME
!!  Grid_getPlaneData
!!
!! SYNOPSIS
!!
!!  Grid_getPlaneData(integer(IN) :: blockID,
!!                    integer(IN) :: gridDataStruct,
!!                    integer(IN) :: structIndex,
!!                    integer(IN) :: beginCount, 
!!                    integer(IN) :: plane,
!!                    integer(IN) :: startingPos(MDIM),
!!                    real(OUT)   :: datablock(dataSize(1),dataSize(2)),
!!                    integer(IN) :: dataSize(2))
!!  
!! DESCRIPTION 
!!  
!!  This routine allows the user to get an entire plane or any contiguous
!!  part of a plane of data depending on the arguments passed. The user is also
!!  allowed to specify if indice counting should begin at the exterior edge
!!  of the block, (that is including guardcells)
!!  or the interior edge of the block 
!!
!!  The data could be either a single variable in one of the grid data 
!!  structures such as cell centered/face centered/scratch, or it could
!!  be quantities derived from the size of the cell, such as cell volume/
!!  face area for all the cells in the specified portion of the block
!!
!!  For a 3d problem an XZ, XY, or YZ plane can be returned
!!  For a 2d problem the entire block can be returned
!!  For a 1d problem an error is returned
!!  
!! ARGUMENTS 
!!
!!  blockID : the local blockid
!!
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
!!                structure. When dataType is one of the grid data structures
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
!!  
!!
!!  beginCount : tells the routine where to start index counting.  beginCount can
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
!!  plane :    specifies the plane
!!             The options are XYPLANE, XZPLANE or YZPLANE defined in constants.h
!!             For a 2d problem this argument is ignored.  
!! 
!!  startingPos(MDIM):
!!           specifies the starting position in each dimension of 
!!           the plane of data being fetched.
!!   
!!           startingPos(1) = i
!!           startingPos(2) = j
!!           startingPos(3) = k
!!
!!           If a problem is 2 dimensions startingPos(3) is irrelevant and
!!           ignored.  If a problem is only 1 dimension this routine doesn't
!!           make any sense and an error is returned.
!!
!!
!!  datablock : a 2 dimensional array containing the data returned.
!!              The dimensions for datablock are 
!!              datablock(range1, range2) 
!!              Various compilers require the dimensions of
!!              datablock to be specified explicitly.  They are defined by 
!!              the next argument "dataSize".  
!!
!!
!!  dataSize : an integer array specifying the dimensions for datablock
!!          
!!          dataSize(1) holds the number of cells to return in 1st plane dim
!!
!!          dataSize(2) holds the number of cells to return in 2nd plane dim
!!                    
!!          Order in dataSize does matter, x comes before y and z.  y comes
!!          before z.  For example, to get an XZPLANE, get the x size dataSize(1) 
!!          and the z size in dataSize(2)                     
!!          
!!
!!
!! EXAMPLE  
!!
!!    EXAMPLE 1:  
!!    Here is a 3d block example getting an entire XZPlane, including guardcells
!!    For each block on a local processor, we will get the entire XZPlane
!!    where y = 5.  beginCount is set to EXTERIOR
!!    meaning that the first cell of the block including guardcells is 
!!    index = 1.  If in this example, the number of guardcells along
!!    all the dimensions is 4 then position 5 is the
!!    first interior cell in a dimension.
!! 
!!
!!    (Can't really draw 3d for this example, but this
!!     picture is meant to help show where 
!!     counting begins when "beginCount" is set to EXTERIOR)
!!    
!!     j
!!    
!!    16         - - - - - - - -  
!!    15         - - - - - - - - 
!!    14         - - - - - - - - 
!!    13         - - - - - - - - 
!!    12 - - - -|-|-|-|-|-|-|-|-|- - - -
!!    11 - - - -|-|-|-|-|-|-|-|-|- - - -
!!    10 - - - -|-|-|-|-|-|-|-|-|- - - -
!!     9 - - - -|-|-|-|-|-|-|-|-|- - - -
!!     8 - - - -|-|-|-|-|-|-|-|-|- - - -
!!     7 - - - -|-|-|-|-|-|-|-|-|- - - -
!!     6 - - - -|-|-|-|-|-|-|-|-|- - - -
!!     5 - - - -|-|-|-|-|-|-|-|-|- - - -
!!     4         - - - - - - - - 
!!     3         - - - - - - - - 
!!     2         - - - - - - - - 
!!     1         - - - - - - - - 
!! i     1 2 3 4 5 6 7 8 9 10111213141516 
!!
!!
!!    #include "Flash.h"
!!    #include "constants.h"
!!
!!    ...
!!       
!!      integer ::    startingPos(MDIM)
!!      integer ::    dataSize(2)
!!      integer ::    blockID
!!      real    ::    dataBlock(:,:)
!!       
!!          startingPos(1) = 1 !getting the entire blkLimitsGC of dim   
!!          startingPos(2) = 5
!!          startingPos(3) = 1 !getting the entire blkLimitsGC of dim
!!
!!          dataSize(1) = blkLimitsGC(HIGH,IAXIS) !This is equivalent to NXB + 2*NGUARD
!!                                                !in FIXEDBLOCKSIZE MODE
!!          dataSize(2) = blkLimitsGC(HIGH,JAXIX) !This is equivalent to NZB + 2*NGUARD
!!                                                !in FIXEDBLOCKSIZE MODE
!!
!!          allocate(datablock(dataSize(1), dataSize(2)))
!!
!!          do blockID = 1, localNumBlocks
!!  
!!             call Grid_getPlaneData(blockID, CENTER, DENS_VAR, XZPLANE, EXTERIOR, &
!!                               startingPos, dataBlock, dataSize)
!!  
!!          end do
!!
!!
!!    EXAMPLE 2:  
!!    In this 2d block example we will get cell volumes for a part of 
!!    a plane for each block on the
!!    local processor.
!!    beginCount is set to INTERIOR, meaning that all the startPos indices
!!    will start where index 1 is the first interior cell of the block.
!!    In this example we will get just 6 cells in the plane, 3 in the i
!!    direction, 2 in the j direction. i start position = 4, j start 
!!    position = 5
!!
!!    (Hard to draw, but this is the idea, stars (*) are the cells to return.
!!     Notice where index counting starts when beginCount is set to INTERIOR.)
!!            - - - - - - - - 
!!            - - - - - - - - 
!!            - - - - - - - - 
!!     j      - - - - - - - - 
!!     8 ----|-|-|-|-|-|-|-|-|----
!!     7 ----|-|-|-|-|-|-|-|-|----
!!     6 ----|-|-|-|*|*|*|-|-|----
!!     5 ----|-|-|-|*|*|*|-|-|----
!!     4 ----|-|-|-|-|-|-|-|-|----
!!     3 ----|-|-|-|-|-|-|-|-|----
!!     2 ----|-|-|-|-|-|-|-|-|----
!!     1 ----|-|-|-|-|-|-|-|-|----
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
!!      integer ::    dataSize(2)
!!      integer ::    blockID
!!      real    ::    dataBlock(:,:)
!!       
!!          startingPos(1) = 4    
!!          startingPos(2) = 5
!!          startingPos(3) = 1 !this doesn't really matter since only 2d
!!
!!          dataSize(1) = 3 !just getting 3 cells in the i dir
!!          dataSize(2) = 2 !just getting 2 cells in the i dir
!!
!!          allocate(datablock(dataSize(1), dataSize(2)))
!!
!!          do blockID = 1, localNumBlocks
!!  
!!             call Grid_getPlaneData(blockID, CELL_VOLUME, 0, XZPLANE, INTERIOR, &
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

subroutine Grid_getPlaneData(blockid, gridDataStruct, structIndex, beginCount, &
     plane, startingPos, datablock, dataSize)

  implicit none

  integer, intent(in) :: blockid, structIndex, beginCount, plane, gridDataStruct
  integer, dimension(MDIM), intent(in) :: startingPos
  integer, dimension(2), intent(in) :: dataSize
  real, dimension(datasize(1), dataSize(2)),intent(out) :: datablock

  datablock=0.0
  return
end subroutine Grid_getPlaneData
