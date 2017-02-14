!!****if* source/Grid/GridSolvers/MultigridMC/poisson/poisson_mg_residual_leafs
!!
!! NAME
!!
!!  poisson_mg_residual_leafs
!!
!! SYNOPSIS
!!
!!  call poisson_mg_residual_leafs(:: irhs,
!!                                  :: ilhs,
!!                                  :: ires,
!!                                  :: flux_cons)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   irhs : 
!!
!!   ilhs : 
!!
!!   ires : 
!!
!!   flux_cons : 
!!
!! AUTOGENROBODOC
!!
!!
!!***

!*******************************************************************************

!  Routine:     mg_residual_leafs()

!  Description: Compute the residual of the equation to be solved using
!               multigrid, given a guess at the solution on a particular
!               level and the source term on leaf blocks.  This version
!               implements residuals for the Poisson equation.

!  Parameters:  
!               irhs        Right-hand side (source) of equation.
!               ilhs        Left-hand side (solution) of equation.
!               ires        Index of variable to receive the residual.

!               flux_cons   If nonzero, match fluxes at boundaries with
!                             finer blocks


subroutine poisson_mg_residual_leafs (irhs, ilhs, ires, flux_cons)

  !===============================================================================

  use gr_mgData

  use Grid_data, ONLY : gr_meshMe,gr_meshNumProcs

  use Grid_interface,    ONLY : Grid_getLocalNumBlks, &
                                Grid_getBlkPtr,       &
                                Grid_releaseBlkPtr,   &
                                Grid_getBlkPhysicalSize, &
                                Grid_getListOfBlocks, &
                                Grid_getBlkIndexLimits, &
                                Grid_setFluxHandling


  use physicaldata, only : flux_x,flux_y,flux_z,consv_fluxes,consv_flux_densities,&
                           diagonals
  use workspace, ONLY: work
  use tree, only : nodetype,lrefine
  use paramesh_dimensions, only : nguard_work
  use paramesh_interfaces, only :  amr_flux_conserve, amr_guardcell

  use Timers_interface, ONLY : Timers_start, Timers_stop

  implicit none
#include "constants.h"
#include "Multigrid.h"
#include "Flash.h"  
#include "mpif.h"

  integer :: level, irhs, ilhs, ires, flux_cons

  integer :: i, j, k, lb, ierr
  real    :: delx, dely, delz
  real    :: norm

  integer, parameter        :: MAXDIMS2 = MDIM
  integer, parameter        :: MFACES2  = 2*MAXDIMS2
  real, dimension(MAXDIMS2) :: size
  real, dimension(MFACES2)  :: neigh2, neigh_type

  real, allocatable, dimension(:,:,:,:) ::  mgfluxx, mgfluxy, mgfluxz

  real, pointer, save, dimension(:,:,:,:,:) :: flux_x_ptr, flux_y_ptr, flux_z_ptr

  real, pointer, save, dimension(:,:,:,:) :: unkt

  integer :: lnblocks2
  integer, save :: myPE2,num_PEs

  real dx,dy,dz,dxdz,dydz,dxdy

  integer, parameter :: nxb = NXB
  integer, parameter :: nyb = NYB
  integer, parameter :: nzb = NZB
  integer, parameter :: nguard = NGUARD
  integer, parameter :: ndim = NDIM

  logical, save :: first_call = .true.

  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC

  real, dimension(2,MDIM) :: boundBox

  integer blockcount,ii,jj,blockID

  real bsize(MDIM),coord(MDIM)

  integer blockList(MAXBLOCKS)

  !=============================================================================

!  call timer_start("residual")

!  if (leaf_only.ne.0) then
!     call timer_start("leaf_only")
!  else
!     call timer_start("not leaf_only")
!  end if

  if (first_call) then
     first_call = .false.
     myPE2 = gr_meshMe
     num_PEs = gr_meshNumProcs
     flux_x_ptr => flux_x
     flux_y_ptr => flux_y
     flux_z_ptr => flux_z
  end if


  call Grid_getLocalNumBlks(lnblocks2)

  allocate(mgfluxx(GRID_ILO_GC:GRID_IHI_GC,GRID_JLO_GC:GRID_JHI_GC,GRID_KLO_GC:GRID_KHI_GC,lnblocks2))
  allocate(mgfluxy(GRID_ILO_GC:GRID_IHI_GC,GRID_JLO_GC:GRID_JHI_GC,GRID_KLO_GC:GRID_KHI_GC,lnblocks2))
  allocate(mgfluxz(GRID_ILO_GC:GRID_IHI_GC,GRID_JLO_GC:GRID_JHI_GC,GRID_KLO_GC:GRID_KHI_GC,lnblocks2))

  mgfluxx = 0.
  mgfluxy = 0.
  mgfluxz = 0.

  flux_x_ptr(FLXMC_FLUX,:,:,:,:) = 0.
  flux_y_ptr(FLXMC_FLUX,:,:,:,:) = 0.
  flux_z_ptr(FLXMC_FLUX,:,:,:,:) = 0.

  !               Update boundary zones of the LHS array.
  !               Currently the mesh package doesn't supply a general boundary
  !               update routine, so we must copy the LHS into the "work" array,
  !               then update its boundaries.  This is handled by mg_bndry().

  ! Copy to work:
  call Grid_getListOfBlocks(LEAF,blockList,blockCount)
  do lb = 1, blockCount
     blockID = blockList(lb)
     ! Point to blocks center vars:
     call Grid_getBlkPtr(blockID,unkt,CENTER)
            
     do k = kli, kui
        do j = jli, jui
           do i = ili, iui
              work(i,j,k,blockID,1) = unkt(ilhs,i,j,k)
           enddo
        enddo
     enddo

     ! Release pointers:
     call Grid_releaseBlkPtr(blockID,unkt,CENTER)
  enddo

  ! Fill GuardCells:
  diagonals = .false.
  call amr_guardcell(myPE2,2,2,nlayersx=1,nlayersy=1,nlayersz=1,maxNodetype_gcWanted=1)
  diagonals = .true.


  !               Compute x, y, and z "fluxes."  Copy block boundary fluxes to
  !               the arrays used by the PARAMESH flux conservation routines.
  if (ndim == 1) then
     
     do lb = 1, blockCount
        
        blockID = blockList(lb)
        
        ! Get BlockSize:
        call Grid_getBlkPhysicalSize(blockID,size)
        delx = nxb/size(1)
        do i = ili-1, iui
           mgfluxx(i,jli:jui,kli:kui,blockID) = & 
                &          delx * (work(i+1,jli:jui,kli:kui,blockID,1) - & 
                &                  work(i,jli:jui,kli:kui,blockID,1))
        enddo
        
        flux_x_ptr(FLXMC_FLUX,1,:,:,blockID) = mgfluxx(ili-1,jli:jui,kli:kui,blockID)
        flux_x_ptr(FLXMC_FLUX,2,:,:,blockID) = mgfluxx(iui,jli:jui,kli:kui,blockID)
        
     enddo
     
  elseif (ndim == 2) then
     
     do lb = 1, blockCount
        
        blockID = blockList(lb)
        
        ! Get BlockSize:
        call Grid_getBlkPhysicalSize(blockID,size)
        
        delx = real(nxb)/size(1)
        dely = real(nyb)/size(2)
        do j = jli-1, jui
           do i = ili-1, iui
              mgfluxx(i,j,kli:kui,blockID) = & 
                   &            delx * (work(i+1,j,kli:kui,blockID,1) - & 
                   &                    work(i,j,kli:kui,blockID,1))
              mgfluxy(i,j,kli:kui,blockID) = & 
                   &            dely * (work(i,j+1,kli:kui,blockID,1) - & 
                   &                    work(i,j,kli:kui,blockID,1))
           enddo
        enddo
        
        
        flux_x_ptr(FLXMC_FLUX,1,:,:,blockID) = mgfluxx(ili-1,jli:jui,kli:kui,blockID)
        flux_x_ptr(FLXMC_FLUX,2,:,:,blockID) = mgfluxx(iui,jli:jui,kli:kui,blockID)
        flux_y_ptr(FLXMC_FLUX,:,1,:,blockID) = mgfluxy(ili:iui,jli-1,kli:kui,blockID)
        flux_y_ptr(FLXMC_FLUX,:,2,:,blockID) = mgfluxy(ili:iui,jui,kli:kui,blockID)
        
        
        
     enddo
     
  else
     
     select case(gr_mgDiffOpDiscretize)
        
     case(2)  ! 2nd order Central difference Solution
        
        call Timers_start("update_fluxes")
        
        do lb = 1, blockCount
           
           blockID = blockList(lb)
           
           
           ! Get BlockSize:
           call Grid_getBlkPhysicalSize(blockID,size)
           
           delx = real(nxb)/size(1)
           dely = real(nyb)/size(2)
           delz = real(nzb)/size(3)
           do k = kli-1, kui
              do j = jli-1, jui
                 do i = ili-1, iui
                    mgfluxx(i,j,k,blockID) = & 
                         &              delx * (work(i+1,j,k,blockID,1) - work(i,j,k,blockID,1))
                    mgfluxy(i,j,k,blockID) = & 
                         &              dely * (work(i,j+1,k,blockID,1) - work(i,j,k,blockID,1))
                    mgfluxz(i,j,k,blockID) = & 
                         &              delz * (work(i,j,k+1,blockID,1) - work(i,j,k,blockID,1))
                 enddo
              enddo
           enddo
           !endif
           flux_x_ptr(FLXMC_FLUX,1,:,:,blockID) = mgfluxx(ili-1,jli:jui,kli:kui,blockID) 
           flux_x_ptr(FLXMC_FLUX,2,:,:,blockID) = mgfluxx(iui,jli:jui,kli:kui,blockID)   
           flux_y_ptr(FLXMC_FLUX,:,1,:,blockID) = mgfluxy(ili:iui,jli-1,kli:kui,blockID) 
           flux_y_ptr(FLXMC_FLUX,:,2,:,blockID) = mgfluxy(ili:iui,jui,kli:kui,blockID)   
           flux_z_ptr(FLXMC_FLUX,:,:,1,blockID) = mgfluxz(ili:iui,jli:jui,kli-1,blockID) 
           flux_z_ptr(FLXMC_FLUX,:,:,2,blockID) = mgfluxz(ili:iui,jli:jui,kui,blockID)   
           
        enddo
        
        call Timers_stop("update_fluxes")
        
     case(4)    ! 4th order Central difference Solution
        do lb = 1, blockCount
           
           blockID = blockList(lb)
           
           ! Get BlockSize:
           call Grid_getBlkPhysicalSize(blockID,size)
           
           delx = nxb/size(1)
           dely = nyb/size(2)
           delz = nzb/size(3)
           do k = kli-2, kui+1
              do j = jli-2, jui+1
                 do i = ili-2, iui+1
                    mgfluxx(i,j,k,blockID) = delx * & 
                         &             ( 9./8.*(work(i+1,j,k,blockID,1) - work(i,j,k,blockID,1)) - &
                         &              1./24.*(work(i+2,j,k,blockID,1) - work(i-1,j,k,blockID,1)))
                    
                    mgfluxy(i,j,k,blockID) = dely * &
                         &             ( 9./8.*(work(i,j+1,k,blockID,1) - work(i,j,k,blockID,1)) - &
                         &              1./24.*(work(i,j+2,k,blockID,1) - work(i,j-1,k,blockID,1)))
                    
                    mgfluxz(i,j,k,blockID) = delz * &
                         &             ( 9./8.*(work(i,j,k+1,blockID,1) - work(i,j,k,blockID,1)) - &
                         &              1./24.*(work(i,j,k+2,blockID,1) - work(i,j,k-1,blockID,1)))
                 enddo
              enddo
           enddo
           flux_x_ptr(FLXMC_FLUX,1,:,:,blockID) = mgfluxx(ili-1,jli:jui,kli:kui,blockID) 
           flux_x_ptr(FLXMC_FLUX,2,:,:,blockID) = mgfluxx(iui,jli:jui,kli:kui,blockID)   
           flux_y_ptr(FLXMC_FLUX,:,1,:,blockID) = mgfluxy(ili:iui,jli-1,kli:kui,blockID) 
           flux_y_ptr(FLXMC_FLUX,:,2,:,blockID) = mgfluxy(ili:iui,jui,kli:kui,blockID)   
           flux_z_ptr(FLXMC_FLUX,:,:,1,blockID) = mgfluxz(ili:iui,jli:jui,kli-1,blockID) 
           flux_z_ptr(FLXMC_FLUX,:,:,2,blockID) = mgfluxz(ili:iui,jli:jui,kui,blockID)   
        enddo
        
     end select
     
     
     
  endif

  !-------------------------------------------------------------------------------
  !               Now use the PARAMESH flux conservation routine to enforce
  !               continuity of the first derivative of the solution by
  !               overwriting the coarse-grid boundary flux at coarse-fine
  !               interfaces with the fine-grid boundary flux.

  if (flux_cons /= 0) then

     !               Call the PARAMESH flux conservation routine.
     call Timers_start("residual_flux_conserve")
     call Grid_setFluxHandling("consv_flux_densities")
     call amr_flux_conserve(MyPE2, 0, 0)
     call Grid_setFluxHandling("consv_fluxes")
     call Timers_stop("residual_flux_conserve")
     
     if (ndim .eq. 1) then
        do lb = 1, blockCount
           blockID = blockList(lb)
           mgfluxx(ili-1,jli:jui,kli:kui,blockID) = flux_x_ptr(FLXMC_FLUX,1,:,:,blockID)
           mgfluxx(iui,jli:jui,kli:kui,blockID)   = flux_x_ptr(FLXMC_FLUX,2,:,:,blockID)
        enddo
     elseif (ndim .eq. 2) then
        do lb = 1, blockCount
           blockID = blockList(lb)
           mgfluxx(ili-1,jli:jui,kli:kui,blockID) = flux_x_ptr(FLXMC_FLUX,1,:,:,blockID)
           mgfluxx(iui,jli:jui,kli:kui,blockID)   = flux_x_ptr(FLXMC_FLUX,2,:,:,blockID)
           mgfluxy(ili:iui,jli-1,kli:kui,blockID) = flux_y_ptr(FLXMC_FLUX,:,1,:,blockID)
           mgfluxy(ili:iui,jui,kli:kui,blockID)   = flux_y_ptr(FLXMC_FLUX,:,2,:,blockID)
        enddo
     elseif (ndim .eq. 3) then
        do lb = 1, blockCount
           blockID = blockList(lb)
           mgfluxx(ili-1,jli:jui,kli:kui,blockID) = flux_x_ptr(FLXMC_FLUX,1,:,:,blockID) 
           mgfluxx(iui,jli:jui,kli:kui,blockID)   = flux_x_ptr(FLXMC_FLUX,2,:,:,blockID) 
           mgfluxy(ili:iui,jli-1,kli:kui,blockID) = flux_y_ptr(FLXMC_FLUX,:,1,:,blockID) 
           mgfluxy(ili:iui,jui,kli:kui,blockID)   = flux_y_ptr(FLXMC_FLUX,:,2,:,blockID) 
           mgfluxz(ili:iui,jli:jui,kli-1,blockID) = flux_z_ptr(FLXMC_FLUX,:,:,1,blockID) 
           mgfluxz(ili:iui,jli:jui,kui,blockID)   = flux_z_ptr(FLXMC_FLUX,:,:,2,blockID) 
        enddo
     endif
     
  endif

  !-------------------------------------------------------------------------------
  !               Now difference the fluxes and subtract from the RHS to
  !               obtain the residual.
  if (ndim == 1) then

     do lb = 1, blockCount
        blockID = blockList(lb)
        
        ! Get BlockSize:
        call Grid_getBlkPhysicalSize(blockID,size)
        
        ! Point to blocks center vars:
        call Grid_getBlkPtr(blockID,unkt,CENTER)      
        
        delx = real(nxb)/size(1)
        do i = ili, iui
           unkt(ires,i,jli:jui,kli:kui) = & 
                &          unkt(irhs,i,jli:jui,kli:kui) - & 
                &          delx * (mgfluxx(i,jli:jui,kli:kui,blockID) - & 
                &                  mgfluxx(i-1,jli:jui,kli:kui,blockID))
        enddo
        
        ! Point to blocks center vars:
        call Grid_releaseBlkPtr(blockID,unkt,CENTER)  
     enddo
     
  elseif (ndim == 2) then
     
     do lb = 1, blockCount
        blockID = blockList(lb)
        
        ! Get BlockSize:
        call Grid_getBlkPhysicalSize(blockID,size)
        
        ! Point to blocks center vars:
        call Grid_getBlkPtr(blockID,unkt,CENTER)      
        
        delx = real(nxb)/size(1)
        dely = real(nyb)/size(2)
        do j = jli, jui
           do i = ili, iui
              unkt(ires,i,j,kli:kui) = & 
                   &            unkt(irhs,i,j,kli:kui) - & 
                   &            delx * (mgfluxx(i,j,kli:kui,blockID) - & 
                   &                    mgfluxx(i-1,j,kli:kui,blockID)) - & 
                   &            dely * (mgfluxy(i,j,kli:kui,blockID) - & 
                   &                    mgfluxy(i,j-1,kli:kui,blockID))
           enddo
        enddo
        
        ! Point to blocks center vars:
        call Grid_releaseBlkPtr(blockID,unkt,CENTER) 
     enddo
     
  else
     
     select case(gr_mgDiffOpDiscretize)
        
     case(2)  ! 2nd order Central difference Solution
        
        call Timers_start("update_values")
        
        do lb = 1, blockCount
           blockID = blockList(lb)
           
           ! Get BlockSize:
           call Grid_getBlkPhysicalSize(blockID,size)
           
           ! Point to blocks center vars:
           call Grid_getBlkPtr(blockID,unkt,CENTER)  
           
           delx = real(nxb)/size(1)
           dely = real(nyb)/size(2)
           delz = real(nzb)/size(3)
           do k = kli, kui
              do j = jli, jui
                 do i = ili, iui
                    unkt(ires,i,j,k) = unkt(irhs,i,j,k) - & 
                         &              delx * (mgfluxx(i,j,k,blockID) - mgfluxx(i-1,j,k,blockID)) - & 
                         &              dely * (mgfluxy(i,j,k,blockID) - mgfluxy(i,j-1,k,blockID)) - & 
                         &              delz * (mgfluxz(i,j,k,blockID) - mgfluxz(i,j,k-1,blockID))
                 enddo
              enddo
           enddo
           
           ! Point to blocks center vars:
           call Grid_releaseBlkPtr(blockID,unkt,CENTER) 
        enddo
        
        call Timers_stop("update_values")
        
     case(4)    ! 4th order Central difference Solution
        do lb = 1, blockCount
           blockID = blockList(lb)
           
           ! Get BlockSize:
           call Grid_getBlkPhysicalSize(blockID,size)
           
           ! Point to blocks center vars:
           call Grid_getBlkPtr(blockID,unkt,CENTER)  
           
           delx = nxb/size(1)
           dely = nyb/size(2)
           delz = nzb/size(3)
           do k = kli, kui
              do j = jli, jui
                 do i = ili, iui
                    
                    unkt(ires,i,j,k) = unkt(irhs,i,j,k) - & 
                         &              delx * &
                         &             ( 9./8.*(mgfluxx(i,j,k,blockID) - mgfluxx(i-1,j,k,blockID))   - &
                         &              1./24.*(mgfluxx(i+1,j,k,blockID) - mgfluxx(i-2,j,k,blockID)))- &
                         &              dely * &
                         &             ( 9./8.*(mgfluxy(i,j,k,blockID) - mgfluxy(i,j-1,k,blockID))   - &
                         &              1./24.*(mgfluxy(i,j+1,k,blockID) - mgfluxy(i,j-2,k,blockID)))- &
                         &              delz * &
                         &             ( 9./8.*(mgfluxz(i,j,k,blockID) - mgfluxz(i,j,k-1,blockID))   - &
                         &              1./24.*(mgfluxz(i,j,k+1,blockID) - mgfluxz(i,j,k-2,blockID)))
                    
                    
                 enddo
              enddo
           enddo
           
           ! Point to blocks center vars:
           call Grid_releaseBlkPtr(blockID,unkt,CENTER) 
        enddo
        
        
     end select
     
  endif
  
  !===============================================================================

  deallocate(mgfluxx,mgfluxy,mgfluxz)

!  if (leaf_only.ne.0) then
!     call timer_stop("leaf_only")
!  else
!     call timer_stop("not leaf_only")
!  end if
!  call timer_stop("residual")
     
  return
end subroutine poisson_mg_residual_leafs



