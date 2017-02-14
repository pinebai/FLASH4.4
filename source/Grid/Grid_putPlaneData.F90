!!****f* source/Grid/Grid_putPlaneData
!!
!! NAME
!!  Grid_putPlaneData
!!
!! SYNOPSIS
!!
!!  Grid_putPlaneData(integer(IN) :: blockID,
!!                    integer(IN) :: gridDataStruct,
!!                    integer(IN) :: structIndex,
!!                    integer(IN) :: beginCount, 
!!                    integer(IN) :: plane,
!!                    integer(IN) :: startingPos(MDIM),
!!                    real(IN)   :: datablock(dataSize(1),dataSize(2)),
!!                    integer(IN) :: dataSize(2))
!!  
!! DESCRIPTION 
!!  
!!  Puts a plane of simulation data for a single variable into the 
!!  specified data structure
!!  
!!  This routine allows the user to put an entire plane or any contiguous
!!  part of a plane of data depending on the arguments passed. The user is also
!!  allowed to specify if index counting should begin at the exterior edge
!!  of the block (that is, including guardcells)
!!  or the interior edge of the block 
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
!!  gridDataStruct : integer value specifying data structure. 
!!                   The options are defined in constants.h and they are :
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
!!
!!  structIndex : integer value that specifies which variable to put into storage.
!!             for example: DENS_VAR, PRES_VAR as defined in Flash.h 
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
!!          before z.  For example, to put an XZPLANE, get the x size in
!!          dataSize(1) and the z size in dataSize(2)                     
!!          
!!
!!
!! EXAMPLE  
!!
!!    EXAMPLE 1:  
!!    Here is a 3d block example putting an entire XZPlane, including guardcells
!!    For each block on a local processor, we will put the entire XZPlane
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
!!          startingPos(IAXIS) = 1 !getting the entire blkLimitsGC of dim   
!!          startingPos(JAXIS) = 5
!!          startingPos(KAXIS) = 1 !getting the entire blkLimitsGC of dim
!!
!!          dataSize(1) = blkLimitsGC(HIGH,IAXIS) !This is equivalent to NXB + 2*NGUARD
!!                                                !in FIXEDBLOCKSIZE MODE
!!          dataSize(2) = blkLimitsGC(HIGH,KAXIS) !This is equivalent to NZB + 2*NGUARD
!!                                                !in FIXEDBLOCKSIZE MODE
!!
!!          allocate(datablock(dataSize(1), dataSize(2)))
!!
!!          do blockID = 1, localNumBlocks
!!  
!!             call Grid_putPlaneData(blockID, CENTER, DENS_VAR, XZPLANE, EXTERIOR, &
!!                               startingPos, dataBlock, dataSize)
!!  
!!          end do
!!
!!
!!    EXAMPLE 2:  
!!    In this 2d block example we will put part of a plane for each block on the
!!    local processor.
!!    beginCount is set to INTERIOR, meaning that all the startPos indices
!!    will start where index 1 is the first interior cell of the block.
!!    In this example we will put just 6 cells in the plane, 3 in the i
!!    direction, 2 in the j direction. i start position = 4, j start 
!!    position = 5
!!
!!    (hard to draw, but this is the idea, stars (*) are the cells to return
!!     notice the where indice counting starts when beginCount is set to INTERIOR)
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
!!          startingPos(IAXIS) = 4    
!!          startingPos(JAXIS) = 5
!!          startingPos(KAXIS) = 1 !this doesn't really matter since only 2d
!!
!!          dataSize(1) = 3 !just putting 3 cells in the i dir
!!          dataSize(2) = 2 !just putting 2 cells in the i dir
!!
!!          allocate(datablock(dataSize(1), dataSize(2)))
!!
!!          do blockID = 1, localNumBlocks
!!  
!!             call Grid_putPlaneData(blockID, CENTER, PRES_VAR, XZPLANE, INTERIOR, &
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

subroutine Grid_putPlaneData(blockid, gridDataStruct, structIndex, beginCount, &
     plane, startingPos, datablock, dataSize)

  implicit none

  integer, intent(in) :: blockid, structIndex, beginCount, plane, gridDataStruct
  integer, dimension(MDIM), intent(in) :: startingPos
  integer, dimension(2), intent(in) :: dataSize
  real, dimension(datasize(1), dataSize(2)),intent(in) :: datablock

  return
end subroutine Grid_putPlaneData
