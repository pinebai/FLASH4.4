!!****f* source/Grid/Grid_getPointData
!!
!! NAME
!!  Grid_getPointData
!!
!! SYNOPSIS
!!
!!  Grid_getPointData(integer(IN) :: blockID,
!!                    integer(IN) :: gridDataStruct,
!!                    integer(IN) :: structIndex,
!!                    integer(IN) :: beginCount, 
!!                    integer(IN) :: position(MDIM),
!!                    real(OUT)   :: datablock)
!!                    
!!  
!! DESCRIPTION 
!!  
!!  Gets a point of data from the specified block
!!  The data could be either a single variable in one of the grid data 
!!  structures such as cell centered/face centered/scratch, or it could
!!  be quantities derived from the size of the cell, such as cell volume/
!!  face area for all the cells in the specified portion of the block
!!  
!!  The user is
!!  allowed to specify if indice counting should begin at the exterior edge
!!  of the block, (that is including guardcells)
!!  or the interior edge of the block 
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
!!
!!  position(MDIM):
!!           specifies the point to return
!!   
!!           position(1) = i
!!           position(2) = j
!!           position(3) = k
!!
!!           If a problem is only 2d, position(3) is ignored.  For 1d problems
!!           position(2) and position(3) are ignored.
!!
!!
!!  datablock : a real value containing the data
!!
!!
!!
!!
!! EXAMPLE  
!!  
!!    Here is a 3d block example getting a point of data, 
!!    for each block on a local processor.  We will get the point at 
!!    position i=5, j=6 and k=7.  beginCount is set to EXTERIOR
!!    meaning that the first cell of the block including guardcells is 
!!    index = 1.  If in this example, NGUARD = 4 then position 5 is the
!!    first interior cell in a dimension.
!!
!! 
!!    (very hard to draw this, especially the j axis. 
!!     picture is just meant to help show where 
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
!!
!!    #include "Flash.h"
!!    #include "constants.h"
!!
!!    ...
!!       
!!      integer ::    position(MDIM)
!!      integer ::    blockID
!!      real    ::    dataBlock
!!
!!          position(1) = 5    
!!          position(2) = 6
!!          position(3) = 7
!!
!!          
!!
!!          do blockID = 1, localNumBlocks
!!  
!!             call Grid_getPointData(blockID, CENTER, DENS_VAR, EXTERIOR, &
!!                               position, dataBlock)
!!  
!!          end do
!!
!!
!!  
!!    In this 2d block example we will get the face area of upper face
!!    of a single cell for each block on a local processor.
!!    beginCount is set to INTERIOR, meaning that all the position indices
!!    will start where index 1 is the first interior cell of the block.
!!    In this example we will get a point where i=4, j=5
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
!!     5 ----|-|-|-|*|-|-|-|-|----
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
!!      integer ::    position(MDIM)
!!      integer ::    blockID
!!      real    ::    dataBlock
!!       
!!          position(1) = 4    
!!          position(2) = 5
!!          position(3) = 1 !ignored since only 2d problem
!!
!!
!!
!!          do blockID = 1, localNumBlocks
!!  
!!             call Grid_getPointData(blockID, CELL_FACEAREA, IHI_FACE, INTERIOR, &
!!                               position, dataBlock)
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

subroutine Grid_getPointData(blockid, gridDataStruct, structIndex, beginCount, &
     position, datablock)


  implicit none

  integer, intent(in) :: blockid, structIndex, beginCount, gridDataStruct
  integer, dimension(MDIM), intent(in) :: position
  real, intent(out) :: datablock

  datablock = 0.0
  return
end subroutine Grid_getPointData
