!!****if* source/Grid/unitTest/Grid_unitTest
!!
!! NAME
!!
!!  Grid_unitTest
!!
!! SYNOPSIS
!!
!!  call Grid_unitTest(integer(in):: fileUnit,
!!                     logical(inout)::perfect  )
!!
!! DESCRIPTION
!!
!!  This unit test exercises the data accessing functions of the Grid unit
!!  The routine has direct access to all the mesh data structures such as 
!!  "unk", "facevarx" etc. It uses the Grid_getBlk/Point/RowData functions 
!!  to fetch some or all of the block data, and then compares it with
!!  the corresponding section of the appropriate array.
!!
!! ARGUMENTS
!!  fileUnit - open f90 write unit
!!  perfect - returns a true if the test passed, false otherwise
!!
!!***

subroutine Grid_unitTest(fileUnit,perfect)


                        
  use physicaldata, ONLY : unk,facevarx,facevary,facevarz
  use Grid_interface, ONLY : Grid_getBlkData, &
                             Grid_getPointData, Grid_getRowData,&
                             Grid_getLocalNumBlks
  use Grid_data, ONLY: gr_ilo, gr_ihi, gr_jlo, gr_jhi, gr_klo, gr_khi, &
       gr_iloGC, gr_ihiGC, gr_jloGC, gr_jhiGC, gr_kloGC, gr_khiGC, &
       gr_iguard, gr_jguard, gr_kguard, gr_meshMe
  implicit none

#include "Flash.h"
#include "constants.h"

  
  integer, intent(in)           :: fileUnit ! Output to file
  logical, intent(inout)        :: perfect  ! Flag to indicate errors


  real, allocatable, dimension(:,:,:,:) :: dataBlock
  integer :: localNumBlks

  integer, dimension(MDIM) :: errorCoords
  integer,dimension(NUNK_VARS) :: vars

  integer,dimension(MDIM) :: grd,startingPos,dataSizeBlk
  integer,dimension(MDIM-1)::dataSizePlane
  integer :: dataSizeRow
  real, allocatable :: datablockBlk(:,:,:)
  real, allocatable :: datablockPlane(:,:)
  real, allocatable :: datablockRow(:)
  real :: datablockPoint

  integer :: i,j,ib,ie,jb,je,kb,ke
  real :: error

  error = 0.0
  grd(IAXIS)=gr_ihiGC-gr_ihi
  grd(JAXIS)=gr_jhiGC-gr_jhi
  grd(KAXIS)=gr_khiGC-gr_khi

  call Grid_getLocalNumBlks(localNumBlks)

  startingPos(1) = 1
  startingPos(2) = 1
  startingPos(3) = 1

  dataSizeBlk(1)=gr_ihi-gr_ilo+1-startingPos(1)+1
  dataSizeBlk(2)=gr_jhi-gr_jlo+1-startingPos(2)+1
  dataSizeBlk(3)=gr_khi-gr_klo+1-startingPos(3)+1
  ib=startingPos(1)+grd(1)
  ie=ib+dataSizeBlk(1)-1
  jb=startingPos(2)+grd(2)
  je=jb+dataSizeBlk(2)-1
  kb=startingPos(3)+grd(3)
  ke=kb+dataSizeBlk(3)-1


  allocate(datablockBlk(dataSizeBlk(1), dataSizeBlk(2), dataSizeBlk(3)))

  if (gr_meshMe.EQ.MASTER_PE) print *, "testing Grid_getBlkData for 1 variable INTERIOR"
  do i = 1,localNumBlks
     call Grid_getBlkData(i, CENTER, DENS_VAR, INTERIOR, startingPos, &
          datablockBlk, dataSizeBlk)
     error = error+maxval(abs(datablockBlk(:,:,:)-&
       unk(DENS_VAR,ib:ie,jb:je,kb:ke,i)))
  end do
  if (gr_meshMe.EQ.MASTER_PE) print*,'the cumulative error is ',error
  if (gr_meshMe.EQ.MASTER_PE) print*,'and the min-max val of density are',maxval(datablockBlk(:,:,:)),&
       minval(datablockBlk(:,:,:))
  
  deallocate(datablockBlk)
  
  !Now test Grid_getBlkData, 1 var, EXTERIOR with guardcells

  startingPos(1) = 1
  startingPos(2) = 1
  startingPos(3) = 9

  dataSizeBlk(1) = gr_ihiGC-startingPos(1)+1
  dataSizeBlk(2) = gr_jhiGC-startingPos(2)+1
  dataSizeBlk(3) = gr_khiGC-startingPos(3)+1

  ib=startingPos(1)
  ie=ib+dataSizeBlk(1)-1
  jb=startingPos(2)
  je=jb+dataSizeBlk(2)-1
  kb=startingPos(3)
  ke=kb+dataSizeBlk(3)-1
  
  allocate(datablockBlk(dataSizeBlk(1), dataSizeBlk(2), dataSizeBlk(3)))
  if (gr_meshMe.EQ.MASTER_PE) print *, "testing Grid_getBlkData, 1 var, EXTERIOR"
  do i = 1,localNumBlks
     call Grid_getBlkData(i, CENTER, DENS_VAR, EXTERIOR, startingPos, &
          datablockBlk, dataSizeBlk)
     error = error+maxval(abs(datablockBlk(:,:,:)-&
          unk(DENS_VAR,ib:ie,jb:je,kb:ke,i)))
  end do
  if (gr_meshMe.EQ.MASTER_PE) print*,'the cumulative error is ',error
  if (gr_meshMe.EQ.MASTER_PE) print*,'and the maximum val of density is ',maxval(abs(datablockBlk(:,:,:)))
  
  
  deallocate(datablockBlk)


  !Now test Grid_getBlkData, ALLVARS, EXTERIOR

 
  startingPos(1) = 1
  startingPos(2) = 1
  startingPos(3) = 1

  dataSizeBlk(1) = gr_ihiGC-startingPos(1)+1
  dataSizeBlk(2) = gr_jhiGC-startingPos(2)+1
  dataSizeBlk(3) = gr_khiGC-startingPos(3)+1

  ib=startingPos(1)
  ie=ib+dataSizeBlk(1)-1
  jb=startingPos(2)
  je=jb+dataSizeBlk(2)-1
  kb=startingPos(3)
  ke=kb+dataSizeBlk(3)-1
  
  if (gr_meshMe.EQ.MASTER_PE) print *, "testing Grid_getBlkData, ALLVAR, EXTERIOR"
  allocate(datablockBlk(dataSizeBlk(1), dataSizeBlk(2), dataSizeBlk(3)))
  
  do i = 1,localNumBlks
     do j = 1,NUNK_VARS
        call Grid_getBlkData(i,CENTER, j, EXTERIOR, startingPos, &
             datablockBlk, dataSizeBlk)
        error=error+(maxval(abs(&
             datablockBlk(:,:,:)-unk(j,ib:ie,jb:je,kb:ke,i))))
     end do
  end do
  if (gr_meshMe.EQ.MASTER_PE) print*,'the cumulative error is', error
  if (gr_meshMe.EQ.MASTER_PE) print*,'and the maximum value is ',maxval(abs(datablockBlk(:,:,:)))
     
  deallocate(datablockBlk)


  !Now test Grid_getPointData, 1 var, EXTERIOR
  startingPos(1) = 5
  startingPos(2) = 6
  startingPos(3) = 7

  if (gr_meshMe.EQ.MASTER_PE) print *, "testing Grid_getPointData, 1 var, EXTERIOR, from inside the blk"
  ib=startingPos(1)
  ie=ib
  jb=startingPos(2)
  je=jb
  kb=startingPos(3)
  ke=kb

  do i = 1,localNumblks
     call Grid_getPointData(i, CENTER, DENS_VAR, EXTERIOR, startingPos, &
          datablockPoint)

     error = error + abs(datablockpoint- &
       unk(DENS_VAR,ib,jb,kb,i))
  end do
  if (gr_meshMe.EQ.MASTER_PE) print*,'the cumulative error is ',error

  startingPos(1) = 7
  startingPos(2) = 12
  startingPos(3) = 13
  if (gr_meshMe.EQ.MASTER_PE) print *, "testing Grid_getPointData, 1 var, EXTERIOR, from guardcells"

  ib=startingPos(1)
  ie=ib
  jb=startingPos(2)
  je=jb
  kb=startingPos(3)
  ke=kb

  do i = 1,localNumblks
     call Grid_getPointData(i, CENTER, DENS_VAR, EXTERIOR, startingPos, &
          datablockPoint)

     error = error + abs(datablockpoint- &
       unk(DENS_VAR,ib,jb,kb,i))
  end do
  if (gr_meshMe.EQ.MASTER_PE) print*,'the cumulative error is ',error
  if (gr_meshMe.EQ.MASTER_PE) print *, datablockPoint
     


  !Now test Grid_getPointData, 1 var, INTERIOR
  if (gr_meshMe.EQ.MASTER_PE) print *, "testing Grid_getPointData, 1 var, INTERIOR"

  startingPos(1) = 7
  startingPos(2) = 3
  startingPos(3) = 4

  ib=startingPos(1)+grd(1)
  ie=ib
  jb=startingPos(2)+grd(2)
  je=jb
  kb=startingPos(3)+grd(3)
  ke=kb

  
  do i = 1,localNumblks
     call Grid_getPointData(i, CENTER, DENS_VAR, INTERIOR, startingPos, &
          datablockPoint)

     error = error + abs(datablockpoint- &
       unk(DENS_VAR,ib,jb,kb,i))
  end do
  if (gr_meshMe.EQ.MASTER_PE) print*,'the cumulative error is ',error





  !Now test Grid_getRowData, 1 var, INTERIOR
  if (gr_meshMe.EQ.MASTER_PE) print *, "testing Grid_getRowData, 1 var, INTERIOR, I AXIS"

  startingPos(1) = 3
  startingPos(2) = 1
  startingPos(3) = 4

  dataSizeRow = gr_ihi-gr_ilo+1-startingPos(1)+1

  ib=startingPos(1)+grd(1)
  ie=ib+dataSizeRow-1
  jb=startingPos(2)+grd(2)
  je=jb
  kb=startingPos(3)+grd(3)
  ke=kb




  allocate(datablockRow(dataSizeRow))
    
  do i = 1,localNumblks
     call Grid_getRowData(i, CENTER, DENS_VAR, INTERIOR, IAXIS, startingPos, &
          datablockRow, dataSizeRow)

     error = error+maxval(abs(datablockRow-unk(DENS_VAR,ib:ie,jb,kb,i)))
     
  end do
  if (gr_meshMe.EQ.MASTER_PE) print*,'error is ',error
 
  deallocate(datablockRow)

 


  !Now test Grid_getRowData, 1 var, INTERIOR JAXIS
  if (gr_meshMe.EQ.MASTER_PE) print *, "testing Grid_getRowData, 1 var, INTERIOR, JAXIS"

  startingPos(1) = 1
  startingPos(2) = 5
  startingPos(3) = 1

  dataSizeRow = gr_jhi-gr_jlo+1-startingPos(2)+1

  ib=startingPos(1)+grd(1)
  ie=ib
  jb=startingPos(2)+grd(2)
  je=jb+dataSizeRow-1
  kb=startingPos(3)+grd(3)
  ke=kb

  allocate(datablockRow(dataSizeRow))
    

  do i = 1,localNumblks
     call Grid_getRowData(i, CENTER, DENS_VAR, INTERIOR, JAXIS, startingPos, &
           datablockRow, dataSizeRow)
     
     error = error+maxval(abs(datablockRow-unk(DENS_VAR,ib,jb:je,kb,i)))
  end do
  if (gr_meshMe.EQ.MASTER_PE) print*,'error is ',error
  
  deallocate(datablockRow)
  




  startingPos(1) = 2
  startingPos(2) = 4
  startingPos(3) = 6

  dataSizeRow = gr_ihiGC-startingPos(1)+1

  ib=startingPos(1)
  ie=ib+dataSizeRow-1
  jb=startingPos(2)
  je=jb
  kb=startingPos(3)
  ke=kb

  allocate(datablockRow(dataSizeRow))
    
  if (gr_meshMe.EQ.MASTER_PE) print *, "testing Grid_getRowData, 1 var, EXTERIOR, IAXIS"
  do i = 1,localNumblks
     call Grid_getRowData(i, CENTER, DENS_VAR, EXTERIOR, IAXIS, &
          startingPos, datablockRow, dataSizeRow)
     error = error+maxval(abs(datablockRow-unk(DENS_VAR,ib:ie,jb,kb,i)))
  end do
  if (gr_meshMe.EQ.MASTER_PE) print*,'error is ',error
 
  deallocate(datablockRow)
  
  return
 
end subroutine Grid_unitTest
