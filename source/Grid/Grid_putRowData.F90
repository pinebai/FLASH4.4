!!****f* source/Grid/Grid_putRowData
!!
!!
!! NAME
!!  Grid_putRowData
!!
!! SYNOPSIS
!!
!!  Grid_putRowData(integer(IN) :: blockID,
!!                  integer(IN) :: gridDataStruct,
!!                  integer(IN) :: structIndex,
!!                  integer(IN) :: beginCount, 
!!                  integer(IN) :: row,
!!                  integer(IN) :: startingPos(MDIM),
!!                  real(IN)   :: datablock(dataSize),
!!                  integer(IN) :: dataSize)
!!  
!! DESCRIPTION 
!!  
!!  Puts a row of simulation data for a single variable into the specified
!!  data structure.
!!  
!!  This routine allows the user to put an entire row or any contiguous
!!  part of a row of data depending on the arguments passed. The user is also
!!  allowed to specify if index counting should begin at the exterior edge
!!  of the block (that is including guardcells),
!!  or the interior edge of the block. 
!!
!!  
!! ARGUMENTS 
!!
!!  blockID : the local blockid
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
!!             for example, DENS_VAR, PRES_VAR as defined in Flash.h 
!!  
!!
!!  beginCount : tells the routine where to start index counting.  beginCount can
!!               be set to INTERIOR or EXTERIOR.  If INTERIOR is specified
!!               guardcell indices are not included and index 1 is the first interior cell. 
!!               If EXTERIOR is specified
!!               the first index, 1, is the leftmost guardcell.  See example
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
!!  dataSize : an integer specifying the size for datablock,
!!          
!!          dataSize specifies the number of cells in the row to set.
!!          
!!
!! EXAMPLE  
!!  
!!    Here is a 3d block example for putting an entire row of data, 
!!    including guardcells for each block on a local processor.  
!!    We will put the entire row in the 
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
!!          dataSize = blkLimitsGC(2,IAXIS) !This is equivalent to NXB + 2*NGUARD
!!                                          ! in FIXEDBLOCKSIZE mode
!!
!!          allocate(datablock(dataSize))
!!
!!          ! ...
!!          ! code to fill the elements of datablock ...
!!          ! ...
!!
!!          do blockID = 1, localNumBlocks
!!  
!!             call Grid_putRowData(blockID, CENTER, DENS_VAR, EXTERIOR, IAXIS, &
!!                               startingPos, dataBlock, dataSize)
!!  
!!          end do
!!
!!
!!  
!!    In this 2d block example we will put part of a row for the pressure variable
!!    for each block on a local processor.
!!    beginCount is set to INTERIOR, meaning that all the startPos indices
!!    will start where index 1 is the first interior cell of the block.
!!    In this example we will put just 3 cells in a row in the i direction
!!    beginning at position i= 4 and j position= 5
!!
!!    (Hard to draw, but this is the idea, stars (*) are the cells to modify.
!!     Notice where index counting starts when beginCount is set to INTERIOR.)
!!            - - - - - - - - 
!!            - - - - - - - - 
!!            - - - - - - - - 
!!     j      - - - - - - - - 
!!     8 ----|-|-|-|-|-|-|-|-|----
!!     7 ----|-|-|-|-|-|-|-|-|----
!!     6 ----|-|-|-|-|-|-|-|-|----
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
!!      integer ::    dataSize
!!      integer ::    blockID
!!      real    ::    dataBlock(:)
!!       
!!          startingPos(1) = 4    
!!          startingPos(2) = 5
!!          startingPos(3) = 1 !this is ignored since only 2d
!!
!!          dataSize = 3 !just putting 3 cells
!!
!!          allocate(datablock(dataSize))
!!
!!          ! ...
!!          ! code to fill the 3 elements of datablock ...
!!          ! ...
!!
!!          do blockID = 1, localNumBlocks
!!  
!!             call Grid_putRowData(blockID, CENTER, PRES_VAR, INTERIOR, IAXIS, &
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


subroutine Grid_putRowData(blockID, gridDataStruct, structIndex, beginCount, &
     row, startingPos, datablock, dataSize)


  implicit none

  integer, intent(in) :: blockID, structIndex, beginCount, row, gridDataStruct
  integer, dimension(MDIM), intent(in) :: startingPos
  integer, intent(in) :: dataSize
  real, dimension(datasize),intent(in) :: datablock


  return
end subroutine Grid_putRowData
