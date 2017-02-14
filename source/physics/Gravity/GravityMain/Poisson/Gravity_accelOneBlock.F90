!!****if* source/physics/Gravity/GravityMain/Poisson/Gravity_accelOneBlock
!!
!! NAME
!!
!!  Gravity_accelOneBlock
!!
!!
!! SYNOPSIS
!!
!!  Gravity_accelOneBlock(integer, intent(in) :: blockID, 
!!                        integer, intent(in) :: ngcellcomp,
!!                        real(:,:,:,:)),intent(out) :: gvec, 
!!                        integer, intent(in),optional :: potentialIndex)
!!                      
!!                      
!!
!! DESCRIPTION
!!
!!  Compute components of the zone-averaged gravitational
!!  acceleration for this block.  Include ngcell layers outside
!!  block interior.
!!
!!  This routine computes the gravitational acceleration for
!!  zones in a given block. First-order
!!  finite-volume differencing is used everywhere.  It is assumed
!!  here that the requisite number of guard cells have peen appropriately
!!  filled for the variable containting the gravitational potential.
!!
!!  Dean Townsley 2008
!!  Contributed to Flash Center at the University of Chicago 2008
!!
!! ARGUMENTS
!!
!!  blockID            -  The local identifier of the block to work on
!!  gvec(:,:,:,:)   -  Array to receive gravitational acceleration
!!                        as as NDIM-dimensional vector.  It is assumed
!!                        the the space provided is the size of the block
!!                        plus all guard cells.  The first index is the vector
!!                        component and the latter are cell indices.
!!  ngcellcomp         -  Number of layers outside of block interior to
!!                        compute gravity
!!  potentialIndex     -  if specified,  Variable # to take as potential.
!!                        Default is GPOT_VAR for the potential stored in the
!!                        gpot slot of unk, which should correspond to the
!!                        potential at the current timestep.
!!
!!
!!***

subroutine Gravity_accelOneBlock ( blockID, ngcellcomp, gvec, potentialIndex)


  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, &
                             Grid_getDeltas, Grid_getBlkIndexLimits, &
                             Grid_getGeometry
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

#include "Flash.h"
#include "constants.h"

  integer, intent(in)                    :: blockID,  ngcellcomp
  real, dimension(:,:,:,:), intent(out)  :: gvec
  integer, intent(in),optional           :: potentialIndex

  real, pointer, dimension(:,:,:,:) :: solnData
  integer, dimension(2,MDIM)        :: blkLimits, blkLimitsGC
  real, dimension(MDIM)             :: delta, inv_2delta
  integer         :: sizeI, sizeJ, sizeK
  integer         :: potVar, geom, istat
  real, dimension(:,:,:), allocatable :: gpot(:,:,:)
  real, dimension(:), allocatable     :: inv_r, inv_sintheta
  integer :: i,j,k
  
  !==================================================

  ! If a variable index is explicitly specified, assume that as the potential
  ! otherwise use the default current potential GPOT_VAR  
  if(present(potentialIndex)) then
     potVar=potentialIndex
  else
     potVar=GPOT_VAR
  end if

  ! get block data and block info
  call Grid_getBlkPtr(blockID, solnData)
  call Grid_getDeltas(blockID, delta)
  inv_2delta(1:NDIM) = 0.5/delta(1:NDIM)
  call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC)
  call Grid_getGeometry(geom)

  sizeI = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
  sizeJ = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
  sizeK = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1

  ! pull out potential for better cache locality
  allocate( gpot(sizeI,sizeJ,sizeK), STAT=istat)
  if (istat/=0) call Driver_abortFlash("unable to allocate gpot in Gravity_accelOneBlock")
  gpot(:,:,:) = solnData(potVar, :,:,:)

  call Grid_releaseBlkPtr(blockID, solnData)

  ! now calculate gravity, depending on geometry
  select case (geom)
  case (CARTESIAN)

     do k = blkLimits(LOW,KAXIS)-K3D*ngcellcomp, blkLimits(HIGH,KAXIS)+K3D*ngcellcomp
        do j = blkLimits(LOW,JAXIS)-K2D*ngcellcomp, blkLimits(HIGH,JAXIS)+K2D*ngcellcomp
           do i = blkLimits(LOW,IAXIS)-ngcellcomp, blkLimits(HIGH,IAXIS)+ngcellcomp

              gvec(IAXIS,i,j,k) = (gpot(i+1,j,k)-gpot(i-1,j,k))*inv_2delta(IAXIS)
              if (NDIM>=2) gvec(JAXIS,i,j,k) = (gpot(i,j+1,k)-gpot(i,j-1,k))*inv_2delta(JAXIS)
              if (NDIM==3) gvec(KAXIS,i,j,k) = (gpot(i,j,k+1)-gpot(i,j,k-1))*inv_2delta(KAXIS)

           enddo
        enddo
     enddo

  case (CYLINDRICAL)

     ! need radius if theta axis is included
     if (NDIM==3) then
        allocate(inv_r(sizeI), STAT=istat)
        if (istat/=0) call Driver_abortFlash("unable to allocate inv_r in Gravity_accelOneBlock")
        inv_r(:) = 1.0/inv_r(:)
     endif
     do k = blkLimits(LOW,KAXIS)-K3D*ngcellcomp, blkLimits(HIGH,KAXIS)+K3D*ngcellcomp
        do j = blkLimits(LOW,JAXIS)-K2D*ngcellcomp, blkLimits(HIGH,JAXIS)+K2D*ngcellcomp
           do i = blkLimits(LOW,IAXIS)-ngcellcomp, blkLimits(HIGH,IAXIS)+ngcellcomp
           
              gvec(IAXIS,i,j,k) = (gpot(i+1,j,k)-gpot(i-1,j,k))*inv_2delta(IAXIS)
              if (NDIM>=2) gvec(JAXIS,i,j,k) = (gpot(i,j+1,k)-gpot(i,j-1,k))*inv_2delta(JAXIS)
              if (NDIM==3) gvec(KAXIS,i,j,k) = inv_r(i)*(gpot(i,j,k+1)-gpot(i,j,k-1))*inv_2delta(KAXIS)

           enddo
        enddo
     enddo
     if (NDIM==3) deallocate(inv_r)

  case (SPHERICAL)

     ! need radius if 2d or 3d
     if (NDIM>=2) then
        allocate(inv_r(sizeI), STAT=istat)
        if (istat/=0) call Driver_abortFlash("unable to allocate inv_r in Gravity_accelOneBlock")
        inv_r(:) = 1.0/inv_r(:)
     endif
     ! need sin_theta if 3d
     if (NDIM==3) then
        allocate(inv_sintheta(sizeJ), STAT=istat)
        if (istat/=0) call Driver_abortFlash("unable to allocate inv_sintheta in Gravity_accelOneBlock")
        inv_sintheta(:) = 1.0/inv_sintheta(:)
     endif
     do k = blkLimits(LOW,KAXIS)-K3D*ngcellcomp, blkLimits(HIGH,KAXIS)+K3D*ngcellcomp
        do j = blkLimits(LOW,JAXIS)-K2D*ngcellcomp, blkLimits(HIGH,JAXIS)+K2D*ngcellcomp
           do i = blkLimits(LOW,IAXIS)-ngcellcomp, blkLimits(HIGH,IAXIS)+ngcellcomp
           
              gvec(IAXIS,i,j,k) = (gpot(i+1,j,k)-gpot(i-1,j,k))*inv_2delta(IAXIS)
              if (NDIM>=2) gvec(JAXIS,i,j,k) = inv_r(i)*(gpot(i,j+1,k)-gpot(i,j-1,k))*inv_2delta(JAXIS)
              if (NDIM==3) gvec(KAXIS,i,j,k) = inv_r(i)*inv_sintheta(j)* &
                                                    (gpot(i,j,k+1)-gpot(i,j,k-1))*inv_2delta(KAXIS)

           enddo
        enddo
     enddo
     if (NDIM==3) deallocate(inv_sintheta)
     if (NDIM>=2) deallocate(inv_r)

  case default
     call Driver_abortFlash("unhandled geometry in Gravity_accelOneBlock")
  end select

  deallocate(gpot)
  
  return
   
end subroutine Gravity_accelOneBlock
