!!****if* source/Grid/GridSolvers/MultigridMC/poisson/poisson_mg_residual
!!
!! NAME
!!
!!  poisson_mg_residual
!!
!! SYNOPSIS
!!
!!  call poisson_mg_residual(:: level,
!!                            :: irhs,
!!                            :: ilhs,
!!                            :: ires,
!!                            :: leaf_only,
!!                            :: flux_cons)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   level : 
!!
!!   irhs : 
!!
!!   ilhs : 
!!
!!   ires : 
!!
!!   leaf_only : 
!!
!!   flux_cons : 
!!
!! AUTOGENROBODOC
!!
!!
!!***

!*******************************************************************************

!  Routine:     mg_residual()

!  Description: Compute the residual of the equation to be solved using
!               multigrid, given a guess at the solution on a particular
!               level and the source term on that level.  This version
!               implements residuals for the Poisson equation.

!  Parameters:  level       Level to compute the residual on.
!               irhs        Right-hand side (source) of equation.
!               ilhs        Left-hand side (solution) of equation.
!               ires        Index of variable to receive the residual.
!               leaf_only   If nonzero, only compute residual on leaf nodes
!                             on the specified level.
!               flux_cons   If nonzero, match fluxes at boundaries with
!                             finer blocks


subroutine poisson_mg_residual (level, irhs, ilhs, ires, leaf_only, & 
                                flux_cons)

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


  use physicaldata, only : flux_x,flux_y,flux_z,consv_fluxes,consv_flux_densities
  use workspace, ONLY: work
  use tree, only : nodetype,lrefine
  use paramesh_dimensions, only : nguard_work
  use paramesh_interfaces, only :  amr_flux_conserve

  use Timers_interface, ONLY : Timers_start, Timers_stop

  implicit none
#include "constants.h"
#include "Multigrid.h"
#include "Flash.h"  
#include "mpif.h"

  integer :: level, irhs, ilhs, ires, leaf_only, flux_cons

  integer :: i, j, k, lb, ierr
  real    :: delx, dely, delz
  real    :: norm
  logical :: some_proc_needs_updating, my_proc_needs_updating
  logical, dimension(MAXBLOCKS) :: update_fluxes
  logical, dimension(MAXBLOCKS) :: update_values

  integer, parameter       :: MAXDIMS2 = 3
  integer, parameter       :: MFACES2 = 2*MAXDIMS2
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

  my_proc_needs_updating   = .false.
  some_proc_needs_updating = .false.

  do lb = 1, lnblocks2
     update_values(lb) = (lrefine(lb) == level)
     if (leaf_only /= 0) then
        update_values(lb) = update_values(lb) .and. (nodetype_save(lb) == 1)
     end if
     if (update_values(lb)) then
        my_proc_needs_updating = .true.
     end if
  end do

  call mpi_allreduce(my_proc_needs_updating, some_proc_needs_updating, &
       1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)

  if (some_proc_needs_updating) then

     do lb = 1, lnblocks2
        update_fluxes(lb) = (lrefine(lb) == level)
        if (flux_cons /= 0) then 
           update_fluxes(lb) = update_fluxes(lb) .or. (lrefine(lb) == level+1)
        end if
        if (leaf_only /= 0) then
           update_fluxes(lb) = update_fluxes(lb) .and. (nodetype_save(lb) == 1)
        end if
        if (update_fluxes(lb)) then
        end if
     end do
     

     !               Update boundary zones of the LHS array.
     !               Currently the mesh package doesn't supply a general boundary
     !               update routine, so we must copy the LHS into the "work" array,
     !               then update its boundaries.  This is handled by mg_bndry().

! ***
     call Timers_start("gr_mgBndry")
     call gr_mgBndry (level, ilhs, nguard_work, leaf_only, MG_COPY_UNK_TO_WORK, MG_STANDALONE) !MG_CONTINUE_SERIES
     call Timers_stop("gr_mgBndry")

     !               Compute x, y, and z "fluxes."  Copy block boundary fluxes to
     !               the arrays used by the PARAMESH flux conservation routines.
     if (ndim == 1) then

        do lb = 1, lnblocks2

           if (update_fluxes(lb)) then
              
              ! Get BlockSize:
              call Grid_getBlkPhysicalSize(lb,size)
              delx = nxb/size(1)
              do i = ili-1, iui
                 mgfluxx(i,jli:jui,kli:kui,lb) = & 
                      &          delx * (work(i+1,jli:jui,kli:kui,lb,1) - & 
                      &                  work(i,jli:jui,kli:kui,lb,1))
              enddo
           endif

           flux_x_ptr(FLXMC_FLUX,1,:,:,lb) = mgfluxx(ili-1,jli:jui,kli:kui,lb)
           flux_x_ptr(FLXMC_FLUX,2,:,:,lb) = mgfluxx(iui,jli:jui,kli:kui,lb)



        enddo

     elseif (ndim == 2) then

        do lb = 1, lnblocks2

           if (update_fluxes(lb)) then

              ! Get BlockSize:
              call Grid_getBlkPhysicalSize(lb,size)

              delx = real(nxb)/size(1)
              dely = real(nyb)/size(2)
              do j = jli-1, jui
                 do i = ili-1, iui
                    mgfluxx(i,j,kli:kui,lb) = & 
                         &            delx * (work(i+1,j,kli:kui,lb,1) - & 
                         &                    work(i,j,kli:kui,lb,1))
                    mgfluxy(i,j,kli:kui,lb) = & 
                         &            dely * (work(i,j+1,kli:kui,lb,1) - & 
                         &                    work(i,j,kli:kui,lb,1))
                 enddo
              enddo
           endif

           flux_x_ptr(FLXMC_FLUX,1,:,:,lb) = mgfluxx(ili-1,jli:jui,kli:kui,lb)
           flux_x_ptr(FLXMC_FLUX,2,:,:,lb) = mgfluxx(iui,jli:jui,kli:kui,lb)
           flux_y_ptr(FLXMC_FLUX,:,1,:,lb) = mgfluxy(ili:iui,jli-1,kli:kui,lb)
           flux_y_ptr(FLXMC_FLUX,:,2,:,lb) = mgfluxy(ili:iui,jui,kli:kui,lb)



        enddo

     else

        select case(gr_mgDiffOpDiscretize)

        case(2)  ! 2nd order Central difference Solution

        call Timers_start("update_fluxes")

        do lb = 1, lnblocks2 

           if (update_fluxes(lb)) then


              ! Get BlockSize:
              call Grid_getBlkPhysicalSize(lb,size)

              delx = (real(nxb)/size(1))
              dely = (real(nyb)/size(2))
              delz = (real(nzb)/size(3))
              do k = kli-1, kui
                 do j = jli-1, jui
                    do i = ili-1, iui
                       mgfluxx(i,j,k,lb) = & 
                            &   delx *     (work(i+1,j,k,lb,1) - work(i,j,k,lb,1)) !delx *
                       mgfluxy(i,j,k,lb) = & 
                            &   dely *     (work(i,j+1,k,lb,1) - work(i,j,k,lb,1)) !dely *
                       mgfluxz(i,j,k,lb) = & 
                            &   delz *     (work(i,j,k+1,lb,1) - work(i,j,k,lb,1)) !delz *
                    enddo
                 enddo
              enddo
           !endif
           flux_x_ptr(FLXMC_FLUX,1,:,:,lb) = mgfluxx(ili-1,jli:jui,kli:kui,lb) 
           flux_x_ptr(FLXMC_FLUX,2,:,:,lb) = mgfluxx(iui,jli:jui,kli:kui,lb)   
           flux_y_ptr(FLXMC_FLUX,:,1,:,lb) = mgfluxy(ili:iui,jli-1,kli:kui,lb) 
           flux_y_ptr(FLXMC_FLUX,:,2,:,lb) = mgfluxy(ili:iui,jui,kli:kui,lb)   
           flux_z_ptr(FLXMC_FLUX,:,:,1,lb) = mgfluxz(ili:iui,jli:jui,kli-1,lb) 
           flux_z_ptr(FLXMC_FLUX,:,:,2,lb) = mgfluxz(ili:iui,jli:jui,kui,lb)   
           endif

        enddo

        call Timers_stop("update_fluxes")

        case(4)    ! 4th order Central difference Solution
        do lb = 1, lnblocks2



           if (update_fluxes(lb)) then

              ! Get BlockSize:
              call Grid_getBlkPhysicalSize(lb,size)

              delx = nxb/size(1)
              dely = nyb/size(2)
              delz = nzb/size(3)
              do k = kli-2, kui+1
                 do j = jli-2, jui+1
                    do i = ili-2, iui+1
                       mgfluxx(i,j,k,lb) = delx * & 
                            &             ( 9./8.*(work(i+1,j,k,lb,1) - work(i,j,k,lb,1)) - &
                            &              1./24.*(work(i+2,j,k,lb,1) - work(i-1,j,k,lb,1)))
         
                       mgfluxy(i,j,k,lb) = dely * &
                            &             ( 9./8.*(work(i,j+1,k,lb,1) - work(i,j,k,lb,1)) - &
                            &              1./24.*(work(i,j+2,k,lb,1) - work(i,j-1,k,lb,1)))
 
                       mgfluxz(i,j,k,lb) = delz * &
                            &             ( 9./8.*(work(i,j,k+1,lb,1) - work(i,j,k,lb,1)) - &
                            &              1./24.*(work(i,j,k+2,lb,1) - work(i,j,k-1,lb,1)))
                    enddo
                 enddo
              enddo
           endif
           flux_x_ptr(FLXMC_FLUX,1,:,:,lb) = mgfluxx(ili-1,jli:jui,kli:kui,lb) 
           flux_x_ptr(FLXMC_FLUX,2,:,:,lb) = mgfluxx(iui,jli:jui,kli:kui,lb)   
           flux_y_ptr(FLXMC_FLUX,:,1,:,lb) = mgfluxy(ili:iui,jli-1,kli:kui,lb) 
           flux_y_ptr(FLXMC_FLUX,:,2,:,lb) = mgfluxy(ili:iui,jui,kli:kui,lb)   
           flux_z_ptr(FLXMC_FLUX,:,:,1,lb) = mgfluxz(ili:iui,jli:jui,kli-1,lb) 
           flux_z_ptr(FLXMC_FLUX,:,:,2,lb) = mgfluxz(ili:iui,jli:jui,kui,lb)   
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
        do i = 1, lnblocks2
           mgfluxx(ili-1,jli:jui,kli:kui,i) = flux_x_ptr(FLXMC_FLUX,1,:,:,i)
           mgfluxx(iui,jli:jui,kli:kui,i)   = flux_x_ptr(FLXMC_FLUX,2,:,:,i)
        enddo
        elseif (ndim .eq. 2) then
           do i = 1, lnblocks2 

              mgfluxx(ili-1,jli:jui,kli:kui,i) = flux_x_ptr(FLXMC_FLUX,1,:,:,i)
              mgfluxx(iui,jli:jui,kli:kui,i)   = flux_x_ptr(FLXMC_FLUX,2,:,:,i)
              mgfluxy(ili:iui,jli-1,kli:kui,i) = flux_y_ptr(FLXMC_FLUX,:,1,:,i)
              mgfluxy(ili:iui,jui,kli:kui,i)   = flux_y_ptr(FLXMC_FLUX,:,2,:,i)
                  
           enddo

        elseif (ndim .eq. 3) then
           do i = 1, lnblocks2
              if (update_values(i)) then
              mgfluxx(ili-1,jli:jui,kli:kui,i) = flux_x_ptr(FLXMC_FLUX,1,:,:,i) 
              mgfluxx(iui,jli:jui,kli:kui,i)   = flux_x_ptr(FLXMC_FLUX,2,:,:,i) 
              mgfluxy(ili:iui,jli-1,kli:kui,i) = flux_y_ptr(FLXMC_FLUX,:,1,:,i) 
              mgfluxy(ili:iui,jui,kli:kui,i)   = flux_y_ptr(FLXMC_FLUX,:,2,:,i) 
              mgfluxz(ili:iui,jli:jui,kli-1,i) = flux_z_ptr(FLXMC_FLUX,:,:,1,i) 
              mgfluxz(ili:iui,jli:jui,kui,i)   = flux_z_ptr(FLXMC_FLUX,:,:,2,i) 
              endif
           enddo
        endif

     endif

     !-------------------------------------------------------------------------------

     !               Now difference the fluxes and subtract from the RHS to
     !               obtain the residual.

     if (ndim == 1) then

        do lb = 1, lnblocks2
           if (update_values(lb)) then

              ! Get BlockSize:
              call Grid_getBlkPhysicalSize(lb,size)

              ! Point to blocks center vars:
              call Grid_getBlkPtr(lb,unkt,CENTER)      

              delx = real(nxb)/size(1)
              do i = ili, iui
                 unkt(ires,i,jli:jui,kli:kui) = & 
                      &          unkt(irhs,i,jli:jui,kli:kui) - & 
                      &          delx * (mgfluxx(i,jli:jui,kli:kui,lb) - & 
                      &                  mgfluxx(i-1,jli:jui,kli:kui,lb))
              enddo

              ! Point to blocks center vars:
              call Grid_releaseBlkPtr(lb,unkt,CENTER)  

           endif
        enddo

     elseif (ndim == 2) then

        do lb = 1, lnblocks2
           if (update_values(lb)) then

              ! Get BlockSize:
              call Grid_getBlkPhysicalSize(lb,size)

              ! Point to blocks center vars:
              call Grid_getBlkPtr(lb,unkt,CENTER)      

              delx = real(nxb)/size(1)
              dely = real(nyb)/size(2)
              do j = jli, jui
                 do i = ili, iui
                    unkt(ires,i,j,kli:kui) = & 
                         &            unkt(irhs,i,j,kli:kui) - & 
                         &            delx * (mgfluxx(i,j,kli:kui,lb) - & 
                         &                    mgfluxx(i-1,j,kli:kui,lb)) - & 
                         &            dely * (mgfluxy(i,j,kli:kui,lb) - & 
                         &                    mgfluxy(i,j-1,kli:kui,lb))
                 enddo
              enddo

              ! Point to blocks center vars:
              call Grid_releaseBlkPtr(lb,unkt,CENTER) 

           endif
        enddo

     else

        select case(gr_mgDiffOpDiscretize)

        case(2)  ! 2nd order Central difference Solution

        call Timers_start("update_values")

        do lb = 1, lnblocks2
           if (update_values(lb)) then


              ! Get BlockSize:
              call Grid_getBlkPhysicalSize(lb,size)

              ! Point to blocks center vars:
              call Grid_getBlkPtr(lb,unkt,CENTER)  

              delx = (real(nxb)/size(1))
              dely = (real(nyb)/size(2))
              delz = (real(nzb)/size(3))
              do k = kli, kui
                 do j = jli, jui
                    do i = ili, iui
                       unkt(ires,i,j,k) = unkt(irhs,i,j,k) - & 
                            &              delx * (mgfluxx(i,j,k,lb) - mgfluxx(i-1,j,k,lb)) - & 
                            &              dely * (mgfluxy(i,j,k,lb) - mgfluxy(i,j-1,k,lb)) - & 
                            &              delz * (mgfluxz(i,j,k,lb) - mgfluxz(i,j,k-1,lb))
                    enddo
                 enddo
              enddo

              ! Point to blocks center vars:
              call Grid_releaseBlkPtr(lb,unkt,CENTER) 

           endif
        enddo

        call Timers_stop("update_values")

        case(4)    ! 4th order Central difference Solution
        do lb = 1, lnblocks2
           if (update_values(lb)) then

              ! Get BlockSize:
              call Grid_getBlkPhysicalSize(lb,size)

              ! Point to blocks center vars:
              call Grid_getBlkPtr(lb,unkt,CENTER)  

              delx = nxb/size(1)
              dely = nyb/size(2)
              delz = nzb/size(3)
              do k = kli, kui
                 do j = jli, jui
                    do i = ili, iui

                       unkt(ires,i,j,k) = unkt(irhs,i,j,k) - & 
                            &              delx * &
                            &             ( 9./8.*(mgfluxx(i,j,k,lb) - mgfluxx(i-1,j,k,lb))   - &
                            &              1./24.*(mgfluxx(i+1,j,k,lb) - mgfluxx(i-2,j,k,lb)))- &
                            &              dely * &
                            &             ( 9./8.*(mgfluxy(i,j,k,lb) - mgfluxy(i,j-1,k,lb))   - &
                            &              1./24.*(mgfluxy(i,j+1,k,lb) - mgfluxy(i,j-2,k,lb)))- &
                            &              delz * &
                            &             ( 9./8.*(mgfluxz(i,j,k,lb) - mgfluxz(i,j,k-1,lb))   - &
                            &              1./24.*(mgfluxz(i,j,k+1,lb) - mgfluxz(i,j,k-2,lb)))


                    enddo
                 enddo
              enddo

              ! Point to blocks center vars:
              call Grid_releaseBlkPtr(lb,unkt,CENTER) 

           endif
        enddo
        

        end select




     endif

     !===============================================================================

  end if

  deallocate(mgfluxx,mgfluxy,mgfluxz)

!  if (leaf_only.ne.0) then
!     call timer_stop("leaf_only")
!  else
!     call timer_stop("not leaf_only")
!  end if
!  call timer_stop("residual")
     
  return
end subroutine poisson_mg_residual



