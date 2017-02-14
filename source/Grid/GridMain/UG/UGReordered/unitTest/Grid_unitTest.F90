!!****if* source/Grid/GridMain/UG/UGReordered/unitTest/Grid_unitTest
!!
!! NAME
!!
!!  Grid_unitTest
!!
!! SYNOPSIS
!!
!!  Grid_unitTest(integer(in):: fileUnit,
!!                logical(inout)::perfect)
!!
!! DESCRIPTION
!!
!!  This routine tests the various functionality of the reordered UG (Uniform Grid) implementation.
!!  Specifically it tests guardcell filling with the routine
!!  Grid_fillGuardCells.  It also tests getting and putting data into
!!  the unk data structure with the routine
!!  Grid_get(Point/Row/Plane/Blk)Data and
!!  Grid_put(Point/Row/Plane/Blk)Data.  If all tests pass then a file
!!  usually named unitTest_000x is written with the line, "All
!!  tests conformed to expected results".  If the tests failed,
!!  various error messages are written to the file.
!!
!! ARGUMENTS
!!
!!  fileUnit : logical unit number for file in which to write error messages
!!  perfect : indicates if all tests passed or not
!!
!!
!!***

subroutine Grid_unitTest(fileUnit,perfect)

  use Grid_data, ONLY : gr_domainBC
  use Grid_interface, ONLY : Grid_fillGuardCells, &
    Grid_getBlkIndexLimits, Grid_getGlobalIndexLimits, &
    Grid_getBlkData, Grid_getBlkCornerID, Grid_dump, Grid_getPointData, Grid_putPointData, Grid_getRowData, &
    Grid_getPlaneData
                        

  implicit none

#include "Flash.h"
#include "constants.h"
  
  integer, intent(in)           :: fileUnit ! Output to file
  logical, intent(inout)        :: perfect  ! Flag to indicate errors

  logical :: bool = .true.
  integer :: n=1, iGuard = 4
  integer :: xSize,ySize,zSize,i,j,k
  integer :: xOffset,yOffset,zOffset
  integer :: yOffset0,zOffset0
  real    :: twopi,xPi,yPi,zPi,tempX,tempY,tempZ,temp
  real    :: thisError, maxError
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer, dimension(MDIM) :: globalRange, index, stride, axis = 1
  integer,dimension(NUNK_VARS)::vars

  integer :: blkID=1

  integer, dimension(MDIM) :: errorCoords
  
  integer :: startingPos(3), dataSizeBlk(3), dataSizeRow, dataSizePlane(2)
  real, allocatable :: datablockBlk(:,:,:)
  real    :: datablockPoint
  real, allocatable :: datablockRow(:)
  real, allocatable :: datablockPlane(:,:)


  call Grid_fillGuardCells(CENTER,ALLDIR)

  call Grid_getBlkIndexLimits(blkID,blkLimits,blkLimitsGC)

  xSize=blkLimitsGC(HIGH,IAXIS)
  ySize=blkLimitsGC(HIGH,JAXIS)
  zSize=blkLimitsGC(HIGH,KAXIS)

  allocate(dataBlockBlk(xSize, ySize, zSize))

  dataSizeBlk(1) = xSize
  dataSizeBlk(2) = ySize
  dataSizeBlk(3) = zSize

  axis = 1

  vars(1)=DENS_VAR

  call Grid_getGlobalIndexLimits(globalRange)
  print*,'the global blkLimits is ',globalRange
  twopi = 8.0*atan(1.0)
  xPi = twopi/globalRange(IAXIS)
  yPi = twopi/globalRange(JAXIS)
  zPi = twopi/globalRange(KAXIS)

  print *, "before Grid_getBlkData"

  call Grid_getBlkData(1, CENTER, DENS_VAR, EXTERIOR, axis, &
       dataBlockBlk, dataSizeBlk)

  print *, "after Grid_getBlkData"

  call Grid_getBlkCornerID(blkID,index,stride)

  !offset = |dist. from LH edge| - |width of guard cells            |   
  xOffset =  (index(IAXIS) - 1)  - (blkLimits(LOW,IAXIS) - blkLimitsGC(LOW,IAXIS))
  yOffset =  (index(JAXIS) - 1)  - (blkLimits(LOW,JAXIS) - blkLimitsGC(LOW,JAXIS))
  zOffset =  (index(KAXIS) - 1)  - (blkLimits(LOW,KAXIS) - blkLimitsGC(LOW,KAXIS))



  do k=blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
    ! offset = |dist. from LH edge| - |width of guard cells            | + iter.  
    zOffset0 = (index(KAXIS) - 1)   - (blkLimits(LOW,KAXIS) - blkLimitsGC(LOW,KAXIS)) + (k-1)
    zOffset  = zOffset0
    ! if we are on a boundary, shift the offset by the width of the
    ! domain in the appropriate dimension to represent the zone from
    ! which data for the real guard cell would be copied.

    ! check for LH boundary
    if (zOffset .lt. 0) then
      select case (gr_domainBC(LOW, KAXIS))
        case (PERIODIC)
          zOffset = zOffset + globalRange(KAXIS)
          print*,'The low K boundary condition is periodic.',k,zOffset0,zOffset
        case (REFLECTING)
           zOffset = abs(zOffset)-1
           print*,'The low K boundary condition is reflecting.',k,zOffset0,zOffset
        case (OUTFLOW)
           zOffset = 0
           print*,'The low K boundary condition is outflow.',k,zOffset0,zOffset
        case default
           print*,'The low K boundary condition is OTHER, which SHOULD NOT HAPPEN!'
      endselect
    ! check for RH boundary
    elseif (zOffset .ge. globalRange(KAXIS)) then
      select case (gr_domainBC(HIGH,KAXIS))
        case (PERIODIC)
          zOffset = zOffset - globalRange(KAXIS)
          print*,'The high K boundary condition is periodic.',k,zOffset0,zOffset
        case (REFLECTING)
           zOffset = globalRange(KAXIS) - (zOffset - globalRange(KAXIS) + 1)
           print*,'The high K boundary condition is reflecting.',k,zOffset0,zOffset
        case (OUTFLOW)
           zOffset = globalRange(KAXIS) - 1
          print*,'The high K boundary condition is outflow.',k,zOffset0,zOffset
        case default
           print*,'The high K boundary condition is OTHER, which SHOULD NOT HAPPEN!'
      endselect
    endif
    tempZ = cos(zOffset * zPi)

    do j=blkLimitsGC(LOW,JAXIS), blkLimitsGC(HIGH,JAXIS)
      ! offset = |dist. from LH edge| - |width of guard cells            | + iter.  
      yOffset0 = (index(JAXIS) - 1)   - (blkLimits(LOW,JAXIS) - blkLimitsGC(LOW,JAXIS)) + (j-1)
      yOffset  = yOffset0
      ! check for LH boundary
      if (yOffset .lt. 0) then
        select case (gr_domainBC(LOW, JAXIS))
          case (PERIODIC)
             yOffset = yOffset + globalRange(JAXIS)
             print*,'   The low J boundary condition is periodic.',j,yOffset0,yOffset
          case (REFLECTING)
             yOffset = abs(yOffset)-1
             print*,'   The low J boundary condition is reflecting.',j,yOffset0,yOffset
          case (OUTFLOW)
             yOffset = 0
             print*,'   The low J boundary condition is outflow.',j,yOffset0,yOffset
          case default
             print*,'** The low J boundary condition is OTHER, which SHOULD NOT HAPPEN!'
        endselect
      ! check for RH boundary
      elseif (yOffset .ge. globalRange(JAXIS)) then
        select case (gr_domainBC(HIGH,JAXIS))
          case (PERIODIC)
             yOffset = yOffset - globalRange(JAXIS)
             print*,'   The high J boundary condition is periodic.',j,yOffset0,yOffset
          case (REFLECTING)
             yOffset = globalRange(JAXIS) - (yOffset - globalRange(JAXIS) + 1)
             print*,'   The high J boundary condition is reflecting.',j,yOffset0,yOffset
          case (OUTFLOW)
             yOffset = globalRange(JAXIS) - 1
             print*,'   The high J boundary condition is outflow.',j,yOffset0,yOffset
          case default
             print*,'** The high J boundary condition is OTHER, which SHOULD NOT HAPPEN!'
        endselect
      endif
      tempY = cos(yOffset * yPi)

      do i=blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH,IAXIS)
        !offset = |dist. from LH edge| - |width of guard cells            | + iter.
        xOffset = (index(IAXIS) - 1)  -  (blkLimits(LOW,IAXIS) - blkLimitsGC(LOW,IAXIS)) + (i-1)
        ! check for LH boundary
        if (xOffset .lt. 0) then
          select case (gr_domainBC(LOW,IAXIS))
            case (PERIODIC)
              xOffset = xOffset + globalRange(IAXIS)
            case (REFLECTING)
              xOffset = abs(xOffset)-1
            case (OUTFLOW)
              xOffset = 0
            case default
               print*,'***** The low I boundary condition is OTHER, which SHOULD NOT HAPPEN!'
          endselect
        ! check for RH boundary
        elseif (xOffset .ge. globalRange(IAXIS)) then
          select case (gr_domainBC(HIGH, IAXIS))
            case (PERIODIC)
              xOffset = xOffset - globalRange(IAXIS)
            case (REFLECTING)
              xOffset = globalRange(IAXIS) - (xOffset - globalRange(IAXIS) + 1)
            case (OUTFLOW)
              xOffset = globalRange(IAXIS) - 1
            case default
               print*,'***** The high I boundary condition is OTHER, which SHOULD NOT HAPPEN!'
          endselect
        endif
        tempX = sin(xOffset * xPi)

        if (.not.isCorner(i,j,k) ) then
          temp = tempX*tempY*tempZ
          if(abs(dataBlockBlk(i,j,k)-temp) .ge. 0.00000000001) then
             perfect = .false.
             write (fileUnit,'("at (x,y,z)=(", I0,",", I0,",", I0, ")")') i,j,k
             write (fileUnit,'("the grid returned: ", F20.17)') dataBlockBlk(i,j,k)
             write (fileUnit,'("versus an expected value of: ", F20.17)') temp
             write (fileUnit,'("size of error: ", F20.17)') abs(dataBlockBlk(i,j,k)-temp)
             write (fileUnit,*)
          endif
        endif
      enddo
    enddo
  enddo
  call Grid_dump(vars,1,1,.true.)
  deallocate(dataBlockBlk)


  print *, "new stuff"

  !!New stuff katie added
  !!Test Grid_getBlkData for 1 variable, INTERIOR count, NOGUARDCELLS

  print *, "testing Grid_getBlkData for 1 variable INTERIOR"
  !print *, "calling Grid_dump"
  !call Grid_dump(vars, 1, 1, .false.)

  startingPos(1) = 1
  startingPos(2) = 1
  startingPos(3) = 1

  dataSizeBlk(1) = NXB
  dataSizeBlk(2) = NYB
  dataSizeBlk(3) = NZB

  
  allocate(datablockBlk(dataSizeBlk(1), dataSizeBlk(2), dataSizeBlk(3)))
    
  call Grid_getBlkData(1, CENTER, DENS_VAR, INTERIOR, startingPos, datablockBlk, dataSizeBlk)
  print *, "printing datablockBlk"
  do i=1, NZB
     print *, "face ", i
     print '(8F5.1)', datablockBlk(:,:,i)
  end do
  deallocate(datablockBlk)


  

  !Now test Grid_getBlkData, 1 var, EXTERIOR with guardcells
  print *, "testing Grid_getBlkData, 1 var, EXTERIOR"

  !print *, "calling Grid_dump"
  !call Grid_dump(vars, 1, 1, .true.)

  startingPos(1) = 1
  startingPos(2) = 1
  startingPos(3) = 1

  dataSizeBlk(1) = GRID_IHI_GC
  dataSizeBlk(2) = GRID_JHI_GC
  dataSizeBlk(3) = GRID_KHI_GC
  
  allocate(datablockBlk(dataSizeBlk(1), dataSizeBlk(2), dataSizeBlk(3)))
    
  call Grid_getBlkData(1, CENTER, DENS_VAR, EXTERIOR, startingPos, datablockBlk, dataSizeBlk)
  print *, "printing datablockBlk"

  do i=1, GRID_KHI_GC
     print *, "face = ", i
     print '(16F5.1)', datablockBlk(:,:,i)
  end do
  deallocate(datablockBlk)








  !Now test Grid_getBlkData, PRES_VAR, EXTERIOR
  print *, "testing Grid_getBlkData, ALLVAR, EXTERIOR"

 
  startingPos(1) = 1
  startingPos(2) = 1
  startingPos(3) = 1

  dataSizeBlk(1) = GRID_IHI_GC
  dataSizeBlk(2) = GRID_JHI_GC
  dataSizeBlk(3) = 1
  
  allocate(datablockBlk(dataSizeBlk(1), dataSizeBlk(2), dataSizeBlk(3)))
    
  print *, "calling Grid_getBlkData"
  call Grid_getBlkData(1, CENTER, PRES_VAR, EXTERIOR, startingPos, datablockBlk, dataSizeBlk)
  

  vars(1)=PRES_VAR
  print *, "calling Grid_dump"
  call Grid_dump(vars, 1, 1, .true.)
     
  print *, "printing datablockBlk"
  print '(16F7.3)', datablockBlk(:,:,:)

 
  deallocate(datablockBlk)











  !Now test Grid_getPointData, 1 var, EXTERIOR
  print *, "testing Grid_getPointData, 1 var, EXTERIOR"

  startingPos(1) = 6
  startingPos(2) = 6
  startingPos(3) = 6

    
  print *, "calling Grid_getPointData"
  call Grid_getPointData(1, CENTER, DENS_VAR, EXTERIOR, startingPos, datablockPoint)

  print *, datablockPoint









  !Now test Grid_getPointData, 1 var, INTERIOR
  print *, "testing Grid_getPointData, 1 var, INTERIOR"

  startingPos(1) = 2
  startingPos(2) = 1
  startingPos(3) = 1
    
  print *, "calling Grid_getPointData"
  call Grid_getPointData(1, CENTER, DENS_VAR, INTERIOR, startingPos, datablockPoint)

  print *, datablockPoint
 



  



  !Now test Grid_getPointData, PRES_VAR, INTERIOR
  print *, "testing Grid_getPointData, PRES_VAR, INTERIOR"

  startingPos(1) = 7
  startingPos(2) = 2
  startingPos(3) = 1

  print *, "calling Grid_getPointData"
  call Grid_getPointData(1, CENTER, PRES_VAR, INTERIOR, startingPos, datablockPoint)


  print '(11F7.3)', datablockPoint
 







  !Now test Grid_getPointData, 1 var INTERIOR
  print *, "testing Grid_getPointData, 1 var, INTERIOR"

  startingPos(1) = 7
  startingPos(2) = 12
  startingPos(3) = 1

    
  print *, "calling Grid_getPointData"
  call Grid_getPointData(1, CENTER, DENS_VAR, INTERIOR, startingPos, datablockPoint)

  print '(11F7.1)', datablockPoint
 



  !Now test Grid_putPointData, 1 var INTERIOR
  print *, "testing Grid_putPointData, 1 var, INTERIOR"

  startingPos(1) = 5
  startingPos(2) = 6
  startingPos(3) = 1

  datablockPoint = 9999.9
  
  print *, "calling Grid_putPointData"
  call Grid_putPointData(1, CENTER, DENS_VAR, INTERIOR, startingPos, datablockPoint)

  print '(11F7.1)', datablockPoint

  print *, "dumping grid"
  vars(1)=DENS_VAR
  call Grid_dump(vars, 1, 1, .false.)







  !Now test Grid_getRowData, 1 var, INTERIOR
  print *, "testing Grid_getRowData, 1 var, INTERIOR, I AXIS"

  startingPos(1) = 2
  startingPos(2) = 4
  startingPos(3) = 6

  dataSizeRow = 3


  allocate(datablockRow(dataSizeRow))
    
  print *, "calling Grid_getRowData"
  call Grid_getRowData(1, CENTER, DENS_VAR, INTERIOR, IAXIS, startingPos, datablockRow, dataSizeRow)


  print '(8F9.1)', datablockRow
 
  deallocate(datablockRow)




#if NDIM > 1
  !Now test Grid_getRowData, 1 var, INTERIOR JAXIS
  print *, "testing Grid_getRowData, 1 var, INTERIOR, JAXIS"

  startingPos(1) = 5
  startingPos(2) = 1
  startingPos(3) = 1

  dataSizeRow = NYB


  allocate(datablockRow(dataSizeRow))
    
  print *, "calling Grid_getRowData"
  call Grid_getRowData(1, CENTER, DENS_VAR, INTERIOR, JAXIS, startingPos, datablockRow, dataSizeRow)

  print '(11F8.1)', datablockRow
 
  deallocate(datablockRow)

#endif


  !Now test Grid_getRowData, 1 var, EXTERIOR IAXIS
  print *, "testing Grid_getRowData, 1 var, EXTERIOR, IAXIS"

  startingPos(1) = 1
  startingPos(2) = 1
  startingPos(3) = 1

  dataSizeRow = GRID_IHI_GC

  allocate(datablockRow(dataSizeRow))
    
  print *, "calling Grid_getRowData"
  call Grid_getRowData(1, CENTER, DENS_VAR, EXTERIOR, IAXIS, startingPos, datablockRow, dataSizeRow)

  print '(11F7.3)', datablockRow
 
  deallocate(datablockRow)





  !Now test Grid_getRowData, PRES_VAR, INTERIOR, IAXIS
  print *, "testing Grid_getRowData, INTERIOR, IAXIS"

  startingPos(1) = 1
  startingPos(2) = 1
  startingPos(3) = 1

  dataSizeRow = NXB

  allocate(datablockRow(dataSizeRow))
    
  print *, "calling Grid_getRowData"
  call Grid_getRowData(1, CENTER, PRES_VAR, INTERIOR, IAXIS, startingPos, datablockRow, dataSizeRow)

  print '(11F7.3)', datablockRow
  
  deallocate(datablockRow)





#if NDIM > 1
  !Now test Grid_getPlaneData, 1 var, INTERIOR, XYPlane
  print *, "testing Grid_getPlaneData, 1 var, INTERIOR, XYPlane"

  startingPos(1) = 1
  startingPos(2) = 1
  startingPos(3) = 1

  dataSizePlane(1) = NXB
  dataSizePlane(2) = NYB

  allocate(datablockPlane(dataSizePlane(1), dataSizePlane(2)))
    
  print *, "calling Grid_getPlaneData"
  call Grid_getPlaneData(1, CENTER, DENS_VAR, INTERIOR, XYPLANE, startingPos, datablockPlane, dataSizePlane)

  print '(8F9.1)', datablockPlane(:,:)

  deallocate(datablockPlane)






  !Now test Grid_getPlaneData, PRES_VAR, INTERIOR, XYPlane
  print *, "testing Grid_getPlaneData, all vars, INTERIOR, XYPlane"

  startingPos(1) = 1
  startingPos(2) = 1
  startingPos(3) = 1

  dataSizePlane(1) = NXB
  dataSizePlane(2) = NYB

  allocate(datablockPlane(dataSizePlane(1), dataSizePlane(2)))
    
  print *, "calling Grid_getPlaneData"
  call Grid_getPlaneData(1, CENTER, PRES_VAR, INTERIOR, XYPLANE, startingPos, datablockPlane, dataSizePlane)

  
  print '(8F7.3)', datablockPlane(:,:)
  
  deallocate(datablockPlane)





  !Now test Grid_getPlaneData, 1 var, EXTERIOR, XYPlane
  print *, "testing Grid_getPlaneData, 1 var, EXTERIOR, XYPlane"

  startingPos(1) = NGUARD + 1
  startingPos(2) = NGUARD + 1
  
  dataSizePlane(1) = NXB
  dataSizePlane(2) = NYB

  allocate(datablockPlane(dataSizePlane(1), dataSizePlane(2)))
    
  print *, "calling Grid_getPlaneData"
  call Grid_getPlaneData(1, CENTER, PRES_VAR, EXTERIOR, XYPLANE, startingPos, datablockPlane, dataSizePlane)
  
  print '(8F7.3)', datablockPlane(:,:)
  
  deallocate(datablockPlane)

#endif


  return

  contains

  ! figure out if we are on a corner (or edge) of the physical
  ! domain, where guard cell values don't matter. Notice we are
  ! never on a corner for a 1d simulation
  function isCorner(eye,jay,kay)

    logical :: isCorner
    integer, intent(in) :: eye,jay,kay
    logical :: eyeIsEdge, jayIsEdge, kayIsEdge

    eyeIsEdge=.false.
    jayIsEdge=.false.
    kayIsEdge=.false.

    if (eye .lt. blkLimits(LOW,IAXIS) .or. eye .gt. blkLimits(HIGH,IAXIS)) eyeIsEdge=.true.
    if (jay .lt. blkLimits(LOW,JAXIS) .or. jay .gt. blkLimits(HIGH,JAXIS)) jayIsEdge=.true.
    if (kay .lt. blkLimits(LOW,KAXIS) .or. kay .gt. blkLimits(HIGH,KAXIS)) kayIsEdge=.true.

    if ((eyeIsEdge .and. jayIsEdge) .or. &
        (eyeIsEdge .and. kayIsEdge) .or. &
        (jayIsEdge .and. kayIsEdge)) then
      isCorner = .true.
    else
      isCorner = .false.
    endif

    return

  end function isCorner
 
end subroutine Grid_unitTest
