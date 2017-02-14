!!****if* source/Grid/GridParticles/GridParticlesMapToMesh/Paramesh/gr_ptSameProcMap
!!
!! NAME
!!  gr_ptSameProcMap
!!
!! SYNOPSIS
!!
!!  gr_ptSameProcMap(integer,dimension(LOW:HIGH,MDIM), intent(IN) :: srcCoords, &
!!                   integer,dimension(LOW:HIGH,MDIM), intent(IN) :: destCoords, &
!!                   integer,dimension(BLKID:REFLEVELDIF), intent(IN) :: negh, &
!!                   integer, INTENT(in) :: varGrid)
!!
!! DESCRIPTION
!!
!!
!! This subroutine will copy the smeared values from the guard cells of
!! source block into the interior of the target block. The target block
!! is a neighbor existing on the same processor. 
!!
!! The input parameters specify the target block, and the coordinates 
!! of the source and destination region.  Data is taken from the source
!! region, operated on, and then copied directly to the destination region.
!! The operations include restriction (neighbor less refined), and 
!! prolongation (neighbor more refined).
!!
!!
!! ARGUMENTS
!!               srcCoords: The coordinates of the guard cell region in the source block.
!!               destCoords: The coordinates of the internal region in the destination block.
!!               negh:  An array containing information about the destination block.
!!               varGrid:   Index of gridded variable to receive interpolated
!!                              quantity
!!
!! 
!!***

!!REORDER(4): solnData

subroutine gr_ptSameProcMap (srcCoords, destCoords, negh, varGrid)

  use gr_ptData, ONLY : gr_ptBuf
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr
  use gr_ptInterface, ONLY : gr_ptProlongSmear
  use Driver_interface, ONLY : Driver_abortFlash
  implicit none

#include "constants.h"
#include "Flash.h"
#include "gr_ptMapToMesh.h"

  integer,dimension(LOW:HIGH,MDIM), intent(IN)  :: srcCoords, destCoords
  integer,dimension(BLKID:REFLEVELDIF), intent(IN):: negh
  integer,intent(IN) :: varGrid

  integer,dimension(MDIM) :: destStart,destEnd,srcStart,srcEnd
  real,dimension(:,:,:,:),pointer :: solnData
  integer, dimension(MDIM) :: upperLimit
  real, dimension(1:2,1:2,1:2) :: prolongedSection
  integer :: destBlkID, refLevelDif
  integer :: i,j,k, i1,j1,k1,i2,j2,k2, incSrc,incDest
  real :: val, srcToDestRefRatio


  destBlkID = negh(BLKID)
  refLevelDif = negh(REFLEVELDIF)
  
  !Copy the coordinates into smaller arrays which
  !are easier to work with.
  srcStart(:) = srcCoords(LOW,:)
  srcEnd(:) = srcCoords(HIGH,:)
  destStart(:) = destCoords(LOW,:)
  destEnd(:) = destCoords(HIGH,:)



  !write(*,"(a,i2,a,i2,a,i2,a,i2,a,i2,a,i2)") &
  !     " srcCoords: i=",srcStart(IAXIS),":",srcEnd(IAXIS),", j=",srcStart(JAXIS),":",&
  !     srcEnd(JAXIS),", k=",srcStart(KAXIS),":",srcEnd(KAXIS)
  !write(*,"(a,i2,a,i2,a,i2,a,i2,a,i2,a,i2)") &
  !     " destCoords: i=",destStart(IAXIS),":",destEnd(IAXIS),", j=",destStart(JAXIS),":",&
  !     destEnd(JAXIS),", k=",destStart(KAXIS),":",destEnd(KAXIS)

  
  call Grid_getBlkPtr(destBlkID, solnData, CENTER)

  srcToDestRefRatio = 2.0**(-refLevelDif)
  incSrc=int(max(1.0,srcToDestRefRatio))
  incDest=int(max(1.0,(1.0/srcToDestRefRatio)))


  !Determine the appropriate upper limit:
  if((refLevelDif == 0).or.(refLevelDif == 1)) then
     upperLimit(:) = srcEnd(1:MDIM)-srcStart(1:MDIM)+1
  else if(refLevelDif == -1) then
     upperLimit(:) = destEnd(1:MDIM)-destStart(1:MDIM)+1
  else
     call Driver_abortFlash("[gr_ptFromNegh]: Unrecognised refinement difference between source and dest blocks")
  end if


  do k = 1, upperLimit(KAXIS)

     k1 = srcStart(KAXIS) + ((k-1)*incSrc)
     k2 = destStart(KAXIS) + ((k-1)*incDest)
     
     do j = 1, upperLimit(JAXIS)

        j1 = srcStart(JAXIS) + ((j-1)*incSrc)
        j2 = destStart(JAXIS) + ((j-1)*incDest)

        do i = 1, upperLimit(IAXIS)

           i1 = srcStart(IAXIS) + ((i-1)*incSrc)
           i2 = destStart(IAXIS) + ((i-1)*incDest)
           
           
           if(refLevelDif == 0) then
              !Do a simple memory copy if the source and destination blocks
              !are at the same refinement.
              !write(*,*)"Neighbor at the same refinement."

              solnData(varGrid,i2,j2,k2)=solnData(varGrid,i2,j2,k2)+gr_ptBuf(i1,j1,k1)


           else if(refLevelDif == 1) then
              !Neighbor more refined.

              val = gr_ptBuf(i1,j1,k1)
              if(val /= 0.0) then

#ifdef DEBUG_GRIDMAPPARTICLES_VERBOSE
                 write(*,"(a,i2,a,i2,a,i2,a,i2,a,i2,a,i2,a,i2,a,i2,a,i2)") &
                      " To me: Prolong from i=",i1, ", j=", j1, ", k=", k1, &
                      " into iRange=",i2,":",i2+1,", jRange=",j2,":",j2+K2D,", kRange=",k2,":",k2+K3D
#endif

                 call gr_ptProlongSmear(i1,j1,k1,prolongedSection)

                 solnData(varGrid,i2:i2+1,j2:j2+K2D,k2:k2+K3D) = &
                      solnData(varGrid,i2:i2+1,j2:j2+K2D,k2:k2+K3D) + &
                      prolongedSection(1:2,1:1+K2D,1:1+K3D)
              end if


           else if(refLevelDif == -1) then
              !Neighbor less refined.

#ifdef DEBUG_GRIDMAPPARTICLES_VERBOSE             
              write(*,"(a,i2,a,i2,a,i2,a,i2,a,i2,a,i2,a,i2,a,i2,a,i2)") &
                   " To me: Restrict from iRange=",i1,":",i1+incSrc-1, &
                   ", jRange=",j1,":",j1+((incSrc-1)*K2D),", kRange=",k1,":",k1+((incSrc-1)*K3D), &
                   " into i=",i2,", j=",j2,", k=",k2
#endif              

              !Divide by 2.0**NDIM because the density from the smaller source cell 
              !is distributed into a larger cell space.
              solnData(varGrid,i2,j2,k2)=solnData(varGrid,i2,j2,k2) &
                   + (sum( gr_ptBuf(i1:i1+incSrc-1,j1:j1+((incSrc-1)*K2D), k1:k1+((incSrc-1)*K3D)) ) / &
                   2.0**NDIM)

           end if
           
        end do
     end do
  end do
  

  call Grid_releaseBlkPtr(destBlkID,solnData,CENTER)


  return

end subroutine gr_ptSameProcMap
