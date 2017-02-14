!!****if* source/physics/Hydro/HydroMain/split/MHD_8Wave/divbDiffuse/hy_8wv_divb
!!
!! NAME
!!
!!  hy_8wv_divb
!!
!!
!! SYNOPSIS
!!
!!  hy_8wv_divb(
!!           integer(IN) :: blockCount,
!!           integer(IN) :: blockList(blockCount),
!!           real(IN)    :: dt)
!!
!!
!! DESCRIPTION
!!
!!  Cleans divb using the parabolic diffusion method of
!!  Marder, J. Comput. Phys., 68, 48, 1987.
!!
!!  The routine also calculates current densities and divergence of velocity
!!  components if they are defined.
!!
!! ARGUMENTS
!!
!!  blockCount - the number of blocks in blockList
!!  blockList  - array holding local IDs of blocks on which to advance
!!  dt         - time step
!!
!!***

!!REORDER(4): U, scratchData, DB

subroutine hy_8wv_divb(blockCount,blockList,dt)

  use Hydro_data,     ONLY : hy_bref, hy_xref
  use Grid_interface, ONLY : Grid_fillGuardCells, Grid_getBlkIndexLimits, &
                             Grid_getBlkPtr, Grid_releaseBlkPtr, Grid_getDeltas

  implicit none

#include "Flash.h"
#include "constants.h"


  !!$ Argument list -------------------------------------
  integer, intent(IN) :: blockCount
  integer, intent(IN), dimension(blockCount) :: blockList
  real, intent(IN)    :: dt
  !!$ ---------------------------------------------------

  integer :: i,j,k,blockID,ierr,lb
  real :: dx,dy,dz,avisc
  real, pointer, dimension(:,:,:,:) :: U,scratchData
  real, dimension(MDIM) :: del
  real :: divb

#ifdef FIXEDBLOCKSIZE 
  integer, PARAMETER :: ibeg = GRID_ILO
  integer, PARAMETER :: jbeg = GRID_JLO
  integer, PARAMETER :: kbeg = GRID_KLO
  integer, PARAMETER :: iend = GRID_IHI
  integer, PARAMETER :: jend = GRID_JHI
  integer, PARAMETER :: kend = GRID_KHI
  integer, PARAMETER :: iSize = GRID_IHI_GC
  integer, PARAMETER :: jSize = GRID_JHI_GC
  integer, PARAMETER :: kSize = GRID_KHI_GC
  real, DIMENSION(IAXIS:KAXIS,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: DB
#else
  real,allocatable,dimension(:,:,:,:) :: DB
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits,blkLimitsGC
  integer :: iSize, jSize, kSize, ibeg, jbeg, kbeg, iend, jend, kend
#endif


  call Grid_fillGuardCells( CENTER, ALLDIR)

  ! Loop over leaf blocks
  do lb= 1,blockCount
     blockID = blockList(lb)

     call Grid_getDeltas(blockID, del)
     dx = del(IAXIS)
     dy = del(JAXIS)
     dz = del(KAXIS)

#ifndef FIXEDBLOCKSIZE
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
     iSize=blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
     jSize=blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
     kSize=blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1
     ibeg=blkLimits(LOW,IAXIS)
     jbeg=blkLimits(LOW,JAXIS)
     kbeg=blkLimits(LOW,KAXIS)
     iend=blkLimits(HIGH,IAXIS)
     jend=blkLimits(HIGH,JAXIS)
     kend=blkLimits(HIGH,KAXIS)
     allocate(DB(IAXIS:KAXIS,iSize,jSize,kSize),stat=ierr)
#endif

     DB = 0.

     ! Get block pointer
     call Grid_getBlkPtr(blockID,U,CENTER)
     call Grid_getBlkPtr(blockID,scratchData,SCRATCH)

     ! Compute viscous divB corrections
#if NDIM == 3
     do k = kbeg-1,kend+1
        do j = jbeg-1,jend+1
           do i = ibeg-1,iend+1
              
              avisc = 0.5/(1./(dx*dx)+1./(dy*dy)+1./(dz*dz))

              DB(IAXIS,i,j,k) = (avisc/dx)* &
                   ((U(MAGX_VAR,i+1,j,k)-2.*U(MAGX_VAR,i,j,k)+U(MAGX_VAR,i-1,j,k))/dx+ &
                    (U(MAGY_VAR,i+1,j+1, k )-U(MAGY_VAR,i+1,j-1, k )+ &
                     U(MAGY_VAR,i-1,j-1, k )-U(MAGY_VAR,i-1,j+1, k ))/(4.*dy)+ &
                    (U(MAGZ_VAR,i+1, j ,k+1)-U(MAGZ_VAR,i+1, j ,k-1)+ &
                     U(MAGZ_VAR,i-1, j ,k-1)-U(MAGZ_VAR,i-1, j ,k+1))/(4.*dz))

              DB(JAXIS,i,j,k) = (avisc/dy)* &
                   ((U(MAGY_VAR,i,j+1,k)-2.*U(MAGY_VAR,i,j,k)+U(MAGY_VAR,i,j-1,k))/dy+ &
                    (U(MAGX_VAR,i+1,j+1, k )-U(MAGX_VAR,i-1,j+1, k )+ &
                     U(MAGX_VAR,i-1,j-1, k )-U(MAGX_VAR,i+1,j-1, k ))/(4.*dx)+ &
                    (U(MAGZ_VAR, i ,j+1,k+1)-U(MAGZ_VAR, i ,j+1,k-1)+ &
                     U(MAGZ_VAR, i ,j-1,k-1)-U(MAGZ_VAR, i ,j-1,k+1))/(4.*dz))
              
              DB(KAXIS,i,j,k) = (avisc/dz)* &
                   ((U(MAGZ_VAR,i,j,k+1)-2.*U(MAGZ_VAR,i,j,k)+U(MAGZ_VAR,i,j,k-1))/dz+ &
                    (U(MAGX_VAR,i+1, j ,k+1)-U(MAGX_VAR,i-1, j ,k+1)+ &
                     U(MAGX_VAR,i-1, j ,k-1)-U(MAGX_VAR,i+1, j ,k-1))/(4.*dx)+ &
                    (U(MAGY_VAR, i ,j+1,k+1)-U(MAGY_VAR, i ,j-1,k+1)+ &
                     U(MAGY_VAR, i ,j-1,k-1)-U(MAGY_VAR, i ,j+1,k-1))/(4.*dy))
           end do
        end do
     end do
     
#else
     do j = jbeg-1,jend+1
        do i = ibeg-1,iend+1
           avisc = 0.5/(1./(dx*dx)+1./(dy*dy))

           DB(IAXIS,i,j,1) = (avisc/dx)* &
                ((U(MAGX_VAR,i+1,j,1)-2.*U(MAGX_VAR,i,j,1)+U(MAGX_VAR,i-1,j,1))/dx+ &
                 (U(MAGY_VAR,i+1,j+1,1)-U(MAGY_VAR,i+1,j-1,1)+ &
                  U(MAGY_VAR,i-1,j-1,1)-U(MAGY_VAR,i-1,j+1,1))/(4.*dy))
           
           DB(JAXIS,i,j,1) = (avisc/dy)* &
                ((U(MAGY_VAR,i,j+1,1)-2.*U(MAGY_VAR,i,j,1)+U(MAGY_VAR,i,j-1,1))/dy+ &
                 (U(MAGX_VAR,i+1,j+1,1)-U(MAGX_VAR,i-1,j+1,1)+ &
                  U(MAGX_VAR,i-1,j-1,1)-U(MAGX_VAR,i+1,j-1,1))/(4.*dx))
           
        end do
     end do
     
#endif
     
     ! Correct the field
     U(MAGX_VAR:MAGZ_VAR,:,:,:) = U(MAGX_VAR:MAGZ_VAR,:,:,:)+DB(IAXIS:KAXIS,:,:,:)

#ifdef MAGP_VAR
     U(MAGP_VAR,:,:,:) = 0.5*(U(MAGX_VAR,:,:,:)**2+U(MAGY_VAR,:,:,:)**2+U(MAGZ_VAR,:,:,:)**2)
#ifdef BETA_VAR
     U(BETA_VAR,:,:,:) = U(PRES_VAR,:,:,:)/U(MAGP_VAR,:,:,:)
#endif
#endif

#if NDIM == 3
     ! Compute and copy DivB into the database
     do k = kbeg,kend
        do j = jbeg,jend
           do i = ibeg,iend

              U(DIVB_VAR,i,j,k) = ((U(MAGX_VAR,i+1, j , k )-U(MAGX_VAR,i-1, j , k ))/dx+ &
                                   (U(MAGY_VAR, i ,j+1, k )-U(MAGY_VAR, i ,j-1, k ))/dy+ &
                                   (U(MAGZ_VAR, i , j ,k+1)-U(MAGZ_VAR, i , j ,k-1))/dz)/ &
                              ((abs(U(MAGX_VAR,i+1, j , k ))+abs(U(MAGX_VAR,i-1, j , k )))/dx+ &
                               (abs(U(MAGY_VAR, i ,j+1, k ))+abs(U(MAGY_VAR, i ,j-1, k )))/dy+ &
                               (abs(U(MAGZ_VAR, i , j ,k+1))+abs(U(MAGZ_VAR, i , j ,k-1)))/dz+ &
                             1.e-10*hy_bref/hy_xref)
#ifdef DIVV_VAR
              U(DIVV_VAR,i,j,k) =&
                   0.5*((U(VELX_VAR,i+1,j,  k  )-U(VELX_VAR,i-1,j,  k  ))/dx &
                       +(U(VELY_VAR,i,  j+1,k  )-U(VELY_VAR,i,  j-1,k  ))/dy &
                       +(U(VELZ_VAR,i,  j,  k+1)-U(VELZ_VAR,i,  j,  k-1))/dz)
#endif
#ifdef CURX_VAR
              U(CURX_VAR,i,j,k)=&
                   0.5*( (U(MAGZ_VAR,i,j+1,k)-U(MAGZ_VAR,i,j-1,k))/dy &
                        -(U(MAGY_VAR,i,j,k+1)-U(MAGY_VAR,i,j,k-1))/dz)
#endif
#ifdef CURY_VAR
              U(CURY_VAR,i,j,k)=&
                    0.5*(-(U(MAGZ_VAR,i+1,j,k)-U(MAGZ_VAR,i-1,j,k))/dx &
                         +(U(MAGX_VAR,i,j,k+1)-U(MAGX_VAR,i,j,k-1))/dz)
#endif
#ifdef CURZ_VAR
              U(CURZ_VAR,i,j,k)=&
                   0.5*( (U(MAGY_VAR,i+1,j,k)-U(MAGY_VAR,i-1,j,k))/dx &
                        -(U(MAGX_VAR,i,j+1,k)-U(MAGX_VAR,i,j-1,k))/dy)
#endif

           end do
        end do
     end do

#elif NDIM == 2
     do j = jbeg,jend
        do i = ibeg,iend

           U(DIVB_VAR,i,j,1) = ((U(MAGX_VAR,i+1, j ,1)-U(MAGX_VAR,i-1, j ,1))/dx+ &
                                (U(MAGY_VAR, i ,j+1,1)-U(MAGY_VAR, i ,j-1,1))/dy) / &
                           ((abs(U(MAGX_VAR,i+1, j ,1))+abs(U(MAGX_VAR,i-1, j ,1)))/dx+ &
                            (abs(U(MAGY_VAR, i ,j+1,1))+abs(U(MAGY_VAR, i ,j-1,1)))/dy+ &
                          1.e-10*hy_bref/hy_xref)

#ifdef DIVV_VAR
           U(DIVV_VAR,i,j,1) =&
                0.5*((U(VELX_VAR,i+1,j,  1  )-U(VELX_VAR,i-1,j,  1  ))/dx &
                    +(U(VELY_VAR,i,  j+1,1  )-U(VELY_VAR,i,  j-1,1  ))/dy) 
#endif
#ifdef CURX_VAR
              U(CURX_VAR,i,j,1)=&
                   0.5*(U(MAGZ_VAR,i,j+1,1)-U(MAGZ_VAR,i,j-1,1))/dy
#endif
#ifdef CURY_VAR
              U(CURY_VAR,i,j,1)=&
                   -0.5*(U(MAGZ_VAR,i+1,j,1)-U(MAGZ_VAR,i-1,j,1))/dx
#endif
#ifdef CURZ_VAR
              U(CURZ_VAR,i,j,1)=&
                   0.5*((U(MAGY_VAR,i+1,j,1)-U(MAGY_VAR,i-1,j,1))/dx &
                       -(U(MAGX_VAR,i,j+1,1)-U(MAGX_VAR,i,j-1,1))/dy)
#endif
#ifdef VECZ_VAR
                 !! advance Az to next time step
                 !! NOTE : This is for 2D only   
                 !! dA/dt = -E = VxB - magVisc*J       

                 U(VECZ_VAR,i,j,1) = U(VECZ_VAR,i,j,1)+&     
                      dt*( U(VELX_VAR,i,j,1)*U(MAGY_VAR,i,j,1)&
                          -U(VELY_VAR,i,j,1)*U(MAGX_VAR,i,j,1)&     
                     -0.5*scratchData(RESI_SCRATCH_GRID_VAR,i,j,1)*&
                     ( (U(MAGY_VAR,i+1,j,1)- U(MAGY_VAR,i-1,j,1))/dx&     
                      -(U(MAGX_VAR,i,j+1,1)- U(MAGX_VAR,i,j-1,1))/dy))    
#endif
        end do
     end do
#else ! NDIM == 1
     do i = ibeg,iend
        U(DIVB_VAR,i,1,1) = 0.
#ifdef CURX_VAR
        U(CURX_VAR,i,1,1) = 0.
#endif
#ifdef CURY_VAR
        U(CURY_VAR,i,1,1)=-0.5*(U(MAGZ_VAR,i+1,1,1)-U(MAGZ_VAR,i-1,1,1))/dx
#endif
#ifdef CURZ_VAR
        U(CURZ_VAR,i,1,1)= 0.5*(U(MAGY_VAR,i+1,1,1)-U(MAGY_VAR,i-1,1,1))/dx
#endif
     enddo
#endif
     
     ! Release pointer
     call Grid_releaseBlkPtr(blockID,U,CENTER)
     call Grid_releaseBlkPtr(blockID,scratchData,SCRATCH)

#ifndef FIXEDBLOCKSIZE     
     deallocate(DB)
#endif
  end do

end subroutine hy_8wv_divb
