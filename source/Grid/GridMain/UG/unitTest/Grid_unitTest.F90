!!****if* source/Grid/GridMain/UG/unitTest/Grid_unitTest
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
!!  This routine tests the internal operation of the UG (Uniform Grid) implementation.
!!  Specifically it tests guardcell filling with the routine
!!  Grid_fillGuardCells.  It also tests getting and putting data into
!!  the unk data structure with the routine
!!  Grid_get(Point/Row/Plane/Blk)Data and
!!  Grid_put(Point/Row/Plane/Blk)Data.  If all tests pass then a file
!!  usually named unitTest_000x is written with the line, "All
!!  tests conformed to expected results".  If the tests failed,
!!  various error messages may be written to the file.
!!
!! ARGUMENTS
!!
!!  fileUnit : logical unit number for file in which to write error messages
!!  perfect : indicates if all tests passed or not
!!
!! NOTES
!!
!!***

subroutine Grid_unitTest(fileUnit,perfect)

  use Grid_data, ONLY : gr_domainBC,gr_blkBC, gr_meshMe
  use Grid_interface, ONLY : Grid_fillGuardCells, &
    Grid_getBlkIndexLimits, Grid_getGlobalIndexLimits, &
    Grid_getBlkData, Grid_getBlkCornerID, Grid_dump, Grid_getCellCoords, Grid_getBlkPtr, Grid_releaseBlkPtr
  use gr_interfacesUT, ONLY : gr_testBoundary
  
  implicit none

#include "Flash.h"
#include "constants.h"
  
  integer, intent(in)           :: fileUnit ! Output to file
  logical, intent(inout)        :: perfect  ! Flag to indicate errors

  logical :: bool = .true.
  integer :: n=1, iGuard = 4
  integer :: xSize,ySize,zSize,i,j,k,ii,jj
  integer :: xOffset,yOffset,zOffset
  integer :: yOffset0,zOffset0
  real    :: xPi,yPi,zPi,tempX,tempY,tempZ,temp
  real    :: thisError, maxError, error
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC
  integer, dimension(MDIM) :: globalRange, index, stride, axis = 1
  integer,dimension(NUNK_VARS) :: vars

  integer, dimension(MDIM) :: errorCoords,strt,fin
  
  integer :: startingPos(3), dataSizeBlk(3), dataSizeRow, dataSizePlane(2)
  real, allocatable :: datablockBlk(:,:,:)
  real    :: datablockPoint
  real, allocatable :: datablockRow(:)
  real, allocatable :: datablockPlane(:,:)

  integer,dimension(MDIM) :: blkSize
  real, allocatable :: dataBlk(:,:,:)
  real, allocatable,dimension(:) :: iCoords, jCoords, kCoords
  real, pointer, dimension(:,:,:,:) :: solnData
  integer :: blockID = 1
  integer :: guard,bcType
  real :: twopi
  logical :: isFace
  twopi = 2.0*PI


  if (gr_meshMe < 100) then
     print*,'PE',gr_meshMe,"before Grid_fillGuardCells(CENTER_FACES,ALLDIR)"
  end if
  call Grid_fillGuardCells(CENTER_FACES,ALLDIR)
  if (gr_meshMe < 100) then
     print*,'PE',gr_meshMe,"after Grid_fillGuardCells(CENTER_FACES,ALLDIR)"
  end if

  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

  xSize=blkLimitsGC(HIGH,IAXIS)
  ySize=blkLimitsGC(HIGH,JAXIS)
  zSize=blkLimitsGC(HIGH,KAXIS)

  allocate(dataBlockBlk(xSize, ySize, zSize))

  dataSizeBlk(1) = xSize
  dataSizeBlk(2) = ySize
  dataSizeBlk(3) = zSize


  axis = 1

  call Grid_getGlobalIndexLimits(globalRange)
  print*,'the global blkLimits is ',globalRange

  xPi = twopi/globalRange(IAXIS)
  yPi = twopi/globalRange(JAXIS)
  zPi = twopi/globalRange(KAXIS)

  if (gr_meshMe < 100) then
     print*,'PE',gr_meshMe,"before Grid_getBlkData for DENS_VAR"
  end if

  call Grid_getBlkData(1, CENTER, DENS_VAR, EXTERIOR, axis, &
       dataBlockBlk, dataSizeBlk)

  if (gr_meshMe < 100) then
     print*,'PE',gr_meshMe,"after Grid_getBlkData for DENS_VAR"
  end if

  call Grid_getBlkCornerID(blockID,index,stride)

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
  vars(1)=DENS_VAR
  call Grid_dump(vars,1,1,.true.)
  deallocate(dataBlockBlk)

  if (gr_meshMe < 100) then
     print*,'PE',gr_meshMe,'Testing of guard cell contents for CENTER variable DENS_VAR all done.'
  end if
  
  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC,FACEX)
  blkSize(:)=blkLimitsGC(HIGH,:)-blkLimitsGC(LOW,:)+1
  allocate(dataBlk(blkSize(1),blkSize(2),blkSize(3)))
  allocate(iCoords(blkSize(IAXIS)))
  allocate(jCoords(blkSize(JAXIS)))
  allocate(kCoords(blkSize(KAXIS)))

  call Grid_getBlkPtr(blockID,solnData,FACEX)
  if (gr_meshMe < 100) then
     print*,'PE',gr_meshMe,'Testing of guard cell contents for first FACEX variable...'
     if (gr_meshMe==MASTER_PE) print*,'shape of dataBlk is         ',shape(dataBlk)
  end if
  call Grid_getCellCoords(IAXIS,blockID,FACES,.true.,iCoords,blkSize(IAXIS))
  call Grid_getCellCoords(JAXIS,blockID,CENTER,.true.,jCoords,blkSize(JAXIS))
  call Grid_getCellCoords(KAXIS,blockID,CENTER,.true.,kCoords,blkSize(KAXIS))
  do k=blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
     do j=blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
        do i=blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
           dataBlk(i,j,k)=cos(twopi*iCoords(i))*cos(twopi*jCoords(j))&
                         *cos(twopi*kCoords(k))
        end do
     end do
  end do

  do i = 1, NDIM
     do j = LOW,HIGH
        if((gr_blkBC(j,i)/=PERIODIC).and.(gr_blkBC(j,i)/=NOT_BOUNDARY)) then
           strt = blkLimitsGC(LOW,:)
           fin  = blkLimitsGC(HIGH,:)
           if(j==LOW)then
              strt(i)=blkLimitsGC(LOW,i)
              fin(i)=blkLimits(LOW,i)-1
           else
              strt(i)=blkLimits(HIGH,i)+1
              fin(i)=blkLimitsGC(HIGH,i)
           end if
           bcType=gr_blkBC(j,i)
           isFace= (i==IAXIS)
           guard = abs(blkLimits(j,i)-blkLimitsGC(j,i))
           call gr_testBoundary(dataBlk,strt,fin,guard,i,bcType,j,isFace)
        end if
     end do
  end do

  do i = 1,1 !! NFACE_VARS
     error=maxval(abs(dataBlk(:,:,:)-solnData(i,:,:,:)))
     if(error > 1.e-12)perfect=.false.
     print*,'and the error for face variable #',i,' is',error
     if(error > 1.e-12) then
        print*,'Location of maximum error at',maxloc(abs(dataBlk(:,:,:)-solnData(i,:,:,:)))
901     format(a,/,(16(17(1x,1PG8.1),/)))
        print 901,'solnData array is',               solnData(i,:,:,:)
        print 901,'dataBlk array is' ,dataBlk(:,:,:)
        print 901,'error array is'   ,dataBlk(:,:,:)-solnData(i,:,:,:)
     end if
  end do

  call Grid_releaseBlkPtr(blockID,solnData,FACEX)
  deallocate(iCoords)
  deallocate(jCoords)
  deallocate(kCoords)
  deallocate(dataBlk)

#if (NDIM>1)
  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC,FACEY)
  blkSize(:)=blkLimitsGC(HIGH,:)-blkLimitsGC(LOW,:)+1
  allocate(dataBlk(blkSize(1),blkSize(2),blkSize(3)))
  allocate(iCoords(blkSize(IAXIS)))
  allocate(jCoords(blkSize(JAXIS)))
  allocate(kCoords(blkSize(KAXIS)))
  
  call Grid_getBlkPtr(blockID,solnData,FACEY)
  if (gr_meshMe < 100) then
     print*,'PE',gr_meshMe,'Testing of guard cell contents for FACEY variables...'
     if (gr_meshMe==MASTER_PE) print*,'shape of dataBlk is         ',shape(dataBlk)
  end if
  call Grid_getCellCoords(IAXIS,blockID,CENTER,.true.,iCoords,blkSize(IAXIS))
  call Grid_getCellCoords(JAXIS,blockID,FACES,.true.,jCoords,blkSize(JAXIS))
  call Grid_getCellCoords(KAXIS,blockID,CENTER,.true.,kCoords,blkSize(KAXIS))
  do k=blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
     do j=blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
        do i=blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
           dataBlk(i,j,k)=cos(twopi*iCoords(i))*cos(twopi*jCoords(j))&
                         *cos(twopi*kCoords(k))
        end do
     end do
  end do

  do i = 1,NDIM
     do j = LOW,HIGH
        if((gr_blkBC(j,i)/=PERIODIC).and.(gr_blkBC(j,i)/=NOT_BOUNDARY)) then
           strt = blkLimitsGC(LOW,:)
           fin  = blkLimitsGC(HIGH,:)
           if(j==LOW)then
              strt(i)=blkLimitsGC(LOW,i)
              fin(i)=blkLimits(LOW,i)-1
           else
              strt(i)=blkLimits(HIGH,i)+1
              fin(i)=blkLimitsGC(HIGH,i)
           end if
           bcType=gr_blkBC(j,i)
           isFace= (i==JAXIS)
           guard = abs(blkLimits(j,i)-blkLimitsGC(j,i))
           call gr_testBoundary(dataBlk,strt,&
                fin,guard,i,bcType,j,isFace)
        end if
     end do
  end do


  do i = 1, NFACE_VARS
     error=maxval(abs(dataBlk(:,:,:)-solnData(i,:,:,:)))
     if(error > 1.e-12)perfect=.false.
     print*,'and the error for face variable #',i,' is',error
     if(error > 1.e-12) then
        print*,'Location of maximum error at',maxloc(abs(dataBlk(:,:,:)-solnData(i,:,:,:)))
     end if
  end do

  call Grid_releaseBlkPtr(blockID,solnData,FACEY)
  deallocate(iCoords)
  deallocate(jCoords)
  deallocate(kCoords)
  deallocate(dataBlk)
#endif
#if (NDIM >2)
  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC,FACEZ)
  blkSize(:)=blkLimitsGC(HIGH,:)-blkLimitsGC(LOW,:)+1
  allocate(dataBlk(blkSize(1),blkSize(2),blkSize(3)))
  allocate(iCoords(blkSize(IAXIS)))
  allocate(jCoords(blkSize(JAXIS)))
  allocate(kCoords(blkSize(KAXIS)))
  
  call Grid_getBlkPtr(blockID,solnData,FACEZ)
  if (gr_meshMe < 100) then
     print*,'PE',gr_meshMe,'Testing of guard cell contents for FACEZ variables...'
     if (gr_meshMe==MASTER_PE) print*,'shape of dataBlk is         ',shape(dataBlk)
  end if
  call Grid_getCellCoords(IAXIS,blockID,CENTER,.true.,iCoords,blkSize(IAXIS))
  call Grid_getCellCoords(JAXIS,blockID,CENTER,.true.,jCoords,blkSize(JAXIS))
  call Grid_getCellCoords(KAXIS,blockID,FACES,.true.,kCoords,blkSize(KAXIS))
  do k=blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
     do j=blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
        do i=blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
           dataBlk(i,j,k)=cos(twopi*iCoords(i))*cos(twopi*jCoords(j))&
                         *cos(twopi*kCoords(k))
        end do
     end do
  end do
  do i = 1,NDIM
     do j = LOW,HIGH
        if((gr_blkBC(j,i)/=PERIODIC).and.(gr_blkBC(j,i)/=NOT_BOUNDARY)) then
           strt = blkLimitsGC(LOW,:)
           fin  = blkLimitsGC(HIGH,:)
           if(j==LOW)then
              strt(i)=blkLimitsGC(LOW,i)
              fin(i)=blkLimits(LOW,i)-1
           else
              strt(i)=blkLimits(HIGH,i)+1
              fin(i)=blkLimitsGC(HIGH,i)
           end if
           bcType=gr_blkBC(j,i)
           isFace= (i==KAXIS)
           guard = abs(blkLimits(j,i)-blkLimitsGC(j,i))
           call gr_testBoundary(dataBlk,strt,&
                fin,guard,i,bcType,j,isFace)
        end if
     end do
  end do
!!  print*,dataBlk
  do i=1,NFACE_VARS
     error=maxval(abs(dataBlk(:,:,:)-solnData(i,:,:,:)))
     if(error > 1.e-12)perfect=.false.
     print*,'and the error for face variable #',i,' is',error
     if(error > 1.e-12) then
        print*,'Location of maximum error at',maxloc(abs(dataBlk(:,:,:)-solnData(i,:,:,:)))
     end if
  end do
  call Grid_releaseBlkPtr(blockID,solnData,FACEZ)
  deallocate(iCoords)
  deallocate(jCoords)
  deallocate(kCoords)
  deallocate(dataBlk)

#endif



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
