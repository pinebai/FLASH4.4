!!****if* source/Grid/GridMain/Grid_putPlaneData
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

!!REORDER(5): unk, facevar[xyz]
!!REORDER(4): solnData

#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif

subroutine Grid_putPlaneData(blockid, gridDataStruct, structIndex, beginCount, &
     plane, startingPos, datablock, dataSize)

  use Grid_data, ONLY : gr_iguard, gr_jguard, gr_kguard
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface, ONLY : Grid_getBlkPtr,Grid_releaseBlkPtr
  use gr_interface, ONLY : gr_getInteriorBlkPtr, gr_releaseInteriorBlkPtr

  implicit none

#include "constants.h"
#include "Flash.h"
#ifdef Grid_releaseBlkPtr
! disabling per-block drift logging because this routine gets called too much
! and clogs the drift logs with too much noise.  see: drift
#undef Grid_releaseBlkPtr
#endif


  integer, intent(in) :: blockid, structIndex, beginCount, plane, gridDataStruct
  integer, dimension(MDIM), intent(in) :: startingPos
  integer, dimension(2), intent(in) :: dataSize
  real, dimension(datasize(1), dataSize(2)),intent(in) :: datablock
  real, pointer, dimension(:,:,:,:) :: solnData


  integer :: i, var, xb, xe, yb, ye, zb, ze, x, y, z
  integer,dimension(MDIM)::begOffset,dataLen


  integer :: imax, jmax, kmax
  logical :: isget
  logical :: getIntPtr

#ifdef DEBUG_GRID
  isget=.true.
  call gr_checkDataType(blockID,gridDataStruct,imax,jmax,kmax,isget)
  
  !plane specific stuff
  if(NDIM == 1) then
     print *, "Error: Grid_putPlaneData"
     call Driver_abortFlash("Grid_putPlaneData.  Can not put plane data for 1d problem")
  end if

  if((plane == XZPLANE) .and. (NDIM < 3)) then
     print *, "Error: Grid_putPlaneData"
     call Driver_abortFlash("Grid_putPlaneData.  Can not put xzplane data for 2d problem")
  end if

  if((plane == YZPLANE) .and. (NDIM < 3)) then
     print *, "Error: Grid_putPlaneData"
     call Driver_abortFlash("Grid_putPlaneData.  Can not put yzplane data for 2d problem")
  end if



  !verify we have a valid blockid
  if((blockid<1).or.(blockid>MAXBLOCKS)) then
     print*,"Error: Grid_putPlaneData : invalid blockid "
     call Driver_abortFlash("Put_putPlaneData : invalid blockid ")
  end if
 
 
  
  !verify beginCount is set to a valid value
  if((beginCount /= INTERIOR) .and. (beginCount /= EXTERIOR)) then
     print *, "Error: Grid_putPlaneData: beginCount set to improper value"
     print *, "beginCount must = INTERIOR or EXTERIOR (defined in constants.h)"
     call Driver_abortFlash("beginCount must = INTERIOR or EXTERIOR (defined in constants.h)")
  end if

  !verify that dataSize isn't too big

  if (plane == XYPLANE .and. (dataSize(1) > imax .or. dataSize(2) > jmax)) then
     print *, "Error: Grid_putPlaneData: dataSize(1) or dataSize(2) too big"
     print *,"You are requesting more cells than block has in a dimension"
     call Driver_abortFlash("Grid_putPlaneData: dataSize(1) or dataSize(2) too big")
  end if

  if (plane==XZPLANE .and. &
     (dataSize(1) > imax .or. &
     dataSize(2) > kmax)) then
     print *, "Error: Grid_putPlaneData: dataSize(1) or dataSize(2) too big"
     print *,"You are requesting more cells than block has in a dimension"
     call Driver_abortFlash("Grid_putPlaneData: dataSize(1) or dataSize(2) too big")
  end if

  if ((plane==YZPLANE) .and. &
     ((dataSize(1) > jmax) .or. &
     (dataSize(2) > kmax))) then
     print *, "Error: Grid_putPlaneData: dataSize(1) or dataSize(2) too big"
     print *,"You are requesting more cells than block has in a dimension"
     call Driver_abortFlash("Grid_putPlaneData: dataSize(1) or dataSize(2) too big")
  end if




  !verify that there is enough space in datablock
  if ((dataSize(1)  < 1) .or. &
       (dataSize(2) < 1)) then
     
     print *, "Error: Grid_putPlaneData: dataSize(1) or (2) too small"
     print *,"You are requesting more < 1 cell in a dimension of block, 1 is the min"
     call Driver_abortFlash("Grid_putPlaneData: dataSize(1) or (2) too small")
  end if
  



  !verify that indicies aren't too big or too small for the block
  if(beginCount == EXTERIOR) then
    
     if (startingPos(1) > imax) then
        call Driver_abortFlash("Grid_putPlaneData startingPos(1) index larger than block")
     end if

     if ((NDIM > 1) .and. (startingPos(2) > jmax)) then
        call Driver_abortFlash("Grid_putPlaneData startingPos(2) index larger than block")
     end if
    
     if ((NDIM > 2) .and. (startingPos(3) > kmax)) then
        call Driver_abortFlash("Grid_putPlaneData startingPos(3) index larger than block")
     end if
    
     if (startingPos(1) < 1) then
        call Driver_abortFlash("Grid_putPlaneData startingPos(1) index smaller than 1")
     end if

     if ((NDIM > 1) .and. (startingPos(2) < 1)) then
        call Driver_abortFlash("Grid_putPlaneData startingPos(2) index smaller than 1")
     end if
    
     if ((NDIM > 2) .and. (startingPos(3) < 1)) then
        call Driver_abortFlash("Grid_putPlaneData startingPos(3) index smaller than 1")
     end if
        
  else !beginCount == INTERIOR

     if ((startingPos(1) + gr_iguard -1) > imax) then
        call Driver_abortFlash("Grid_putPlaneData startingPos(1) index larger than block")
     end if

     if ((NDIM > 1) .and. ((startingPos(2) + gr_jguard -1) > jmax)) then
        call Driver_abortFlash("Grid_putPlaneData startingPos(2) index larger than block")
     end if
    
     if ((NDIM > 2) .and. ((startingPos(3) + gr_kguard -1) > kmax)) then
        call Driver_abortFlash("Grid_putPlaneData startingPos(3) index larger than block")
     end if
    
     if (startingPos(1) < 1) then
        call Driver_abortFlash("Grid_putPlaneData startingPos(1) index smaller than 1")
     end if

     if ((NDIM > 1) .and. (startingPos(2) < 1)) then
        call Driver_abortFlash("Grid_putPlaneData startingPos(2) index smaller than 1")
     end if
    
     if ((NDIM > 2) .and. (startingPos(3) < 1)) then
        call Driver_abortFlash("Grid_putPlaneData startingPos(3) index smaller than 1")
     end if

  end if
  

  !more verification of indicies
  !check size and starting pos
  if(beginCount == EXTERIOR) then
     if(plane == XYPLANE) then
        if ((startingPos(IAXIS) + dataSize(1) -1) > imax) then
           print *, "Error: Grid_putPlaneData"
           call Driver_abortFlash("Grid_putPlaneData indicies too large")
        end if
        if ((startingPos(JAXIS) + dataSize(2) -1) > jmax) then
           print *, "Error: Grid_putPlaneData"
           call Driver_abortFlash("Grid_putPlaneData indicies too large")
        end if
     end if

     if(plane == XZPLANE) then
        if ((startingPos(IAXIS) + dataSize(1) -1) > imax) then
           print *, "Error: Grid_putPlaneData"
           call Driver_abortFlash("Grid_putPlaneData indicies too large")
        end if
        if ((startingPos(KAXIS) + dataSize(2) -1) > kmax) then
           print *, "Error: Grid_putPlaneData"
           call Driver_abortFlash("Grid_putPlaneData indicies too large")
        end if
     end if

     if(plane == YZPLANE) then
        if ((startingPos(JAXIS) + dataSize(1) -1) > jmax) then
           print *, "Error: Grid_putPlaneData"
           call Driver_abortFlash("Grid_putPlaneData indicies too large")
        end if
        if ((startingPos(KAXIS) + dataSize(2) -1) > kmax) then
           print *, "Error: Grid_putPlaneData"
           call Driver_abortFlash("Grid_putPlaneData indicies too large")
        end if
     end if

  !if INTERIOR counting, check same things   
  else   if(beginCount == INTERIOR) then
     if(plane == XYPLANE) then
        if ((startingPos(IAXIS) + dataSize(1) + gr_iguard -1) > imax) then
           print *, "Error: Grid_putPlaneData"
           call Driver_abortFlash("Grid_putPlaneData indicies too large")
        end if
        if ((startingPos(JAXIS) + dataSize(2) + gr_jguard -1) > jmax) then
           print *, "Error: Grid_putPlaneData"
           call Driver_abortFlash("Grid_putPlaneData indicies too large")
        end if
     end if

     if(plane == XZPLANE) then
        if ((startingPos(IAXIS) + dataSize(1) + gr_iguard -1) > imax) then
           print *, "Error: Grid_putPlaneData"
           call Driver_abortFlash("Grid_putPlaneData indicies too large")
        end if
        if ((startingPos(KAXIS) + dataSize(2) + gr_kguard -1) > kmax) then
           print *, "Error: Grid_putPlaneData"
           call Driver_abortFlash("Grid_putPlaneData indicies too large")
        end if
     end if

     if(plane == YZPLANE) then
        if ((startingPos(JAXIS) + dataSize(1) + gr_jguard -1) > jmax) then
           print *, "Error: Grid_putPlaneData"
           call Driver_abortFlash("Grid_putPlaneData indicies too large")
        end if
        if ((startingPos(KAXIS) + dataSize(2) + gr_kguard -1) > kmax) then
           print *, "Error: Grid_putPlaneData"
           call Driver_abortFlash("Grid_putPlaneData indicies too large")
        end if
     end if
  end if

#endif  

  call gr_getDataOffsets(blockID,gridDataStruct,startingPos,dataLen,beginCount,begOffset,getIntPtr)


  yb=1
  ye=1
  zb=1
  ze=1
  xb = startingPos(IAXIS) + begOffset(IAXIS)
  if(NDIM>1)yb = startingPos(JAXIS) + begOffset(JAXIS)
  if(NDIM>2)zb = startingPos(KAXIS) + begOffset(KAXIS)
  
  if (plane == XYPLANE) then
     xe = xb + dataSize(1) -1
     if(NDIM>1)ye = yb + dataSize(2) -1
     ze=zb
  elseif(plane == XZPLANE) then
     xe = xb + dataSize(1) -1
     if(NDIM>2)ze = zb + dataSize(2) -1
     ye = yb
  elseif(plane == YZPLANE) then
     xe = xb
     if(NDIM>1)ye = yb + dataSize(1) -1
     if(NDIM>2)ze = zb + dataSize(2) -1
  else
     call Driver_abortFlash("Grid_getPlaneData : invalid plane spec")
  end if

  if(getIntPtr) then
     call gr_getInteriorBlkPtr(blockID,solnData,gridDataStruct)
     if(plane==XYPLANE)solnData(structIndex,xb:xe,yb:ye,zb) = datablock(:,:)
     if(plane==XZPLANE)solnData(structIndex,xb:xe,yb,zb:ze) = datablock(:,:)
     if(plane==YZPLANE)solnData(structIndex,xb,yb:ye,zb:ze) = datablock(:,:)
     call gr_releaseInteriorBlkPtr(blockID,solnData,gridDataStruct)
  else
     call Grid_getBlkPtr(blockID,solnData,gridDataStruct)
!!$     if(gridDataStruct==SCRATCH) then
!!$        if(plane==XYPLANE)solnData(xb:xe,yb:ye,zb,structIndex) = datablock(:,:)
!!$        if(plane==XZPLANE)solnData(xb:xe,yb,zb:ze,structIndex) = datablock(:,:)
!!$        if(plane==YZPLANE)solnData(xb,yb:ye,zb:ze,structIndex) = datablock(:,:)
!!$     else
!!$     end if
     if(plane==XYPLANE)solnData(structIndex,xb:xe,yb:ye,zb) = datablock(:,:)
     if(plane==XZPLANE)solnData(structIndex,xb:xe,yb,zb:ze) = datablock(:,:)
     if(plane==YZPLANE)solnData(structIndex,xb,yb:ye,zb:ze) = datablock(:,:)
     call Grid_releaseBlkPtr(blockID,solnData,gridDataStruct)
  end if

  return
end subroutine Grid_putPlaneData
