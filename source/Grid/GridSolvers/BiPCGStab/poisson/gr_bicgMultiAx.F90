!!****if* source/Grid/GridSolvers/BiPCGStab/poisson/gr_bicgMultiAx
!!
!! NAME
!!
!!  gr_bicgMultiAx
!!
!! SYNOPSIS
!!
!!  call gr_bicgMultiAx(integer(in) :: irhs,
!!                      integer(in) :: ilhs,
!!                      logical(in) :: gcellflg)
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
!!   gcellflg : 
!!
!! AUTOGENROBODOC
!!
!!
!!***

!*******************************************************************************

!  Routine:     gr_bicgMultiAx()

!  Description: Computes the product y=A*x


!  Parameters:  irhs        Right-hand side, index of x.
!               ilhs        Left-hand side, index of y.


subroutine gr_bicgMultiAx (irhs, ilhs, gcellflg)

  !===============================================================================

  use gr_bicgData, only : ili, iui, jli, jui, kli, kui

  use Grid_data, ONLY : gr_meshMe,gr_meshNumProcs

  use Grid_interface,    ONLY : Grid_getLocalNumBlks, &
                                Grid_getBlkPtr,       &
                                Grid_releaseBlkPtr,   &
                                Grid_getBlkPhysicalSize, &
                                Grid_getListOfBlocks, &
                                Grid_getBlkIndexLimits

  use gr_bicgInterface, ONLY : gr_bicgBndry  


  use physicaldata, only : flux_x,flux_y,flux_z
  use workspace, ONLY: work
  use tree, only : nodetype,lrefine
  use paramesh_dimensions, only : nguard_work
  use paramesh_interfaces, only :  amr_flux_conserve

  use Timers_interface, ONLY: Timers_start, Timers_stop

  implicit none
#include "constants.h"
#include "Flash.h"  
#include "mpif.h"

  integer, intent(in) :: irhs, ilhs
  logical, intent(in) :: gcellflg

  integer :: i, j, k, lb, ierr
  real    :: delx, dely, delz
!  real    :: norm
!  logical :: some_proc_needs_updating, my_proc_needs_updating
  logical, dimension(MAXBLOCKS) :: update_fluxes
!  logical, dimension(MAXBLOCKS) :: update_values

!  integer, parameter       :: MAXDIMS2 = 3
!  integer, parameter       :: MFACES2 = 2*MAXDIMS2
  real, dimension(MDIM) :: size
!  real, dimension(MFACES2)  :: neigh2, neigh_type



  real    :: mgfluxx(GRID_ILO_GC:GRID_IHI_GC,GRID_JLO_GC:GRID_JHI_GC,GRID_KLO_GC:GRID_KHI_GC, & 
       &                   MAXBLOCKS)
  real    :: mgfluxy(GRID_ILO_GC:GRID_IHI_GC,GRID_JLO_GC:GRID_JHI_GC,GRID_KLO_GC:GRID_KHI_GC, & 
       &                   MAXBLOCKS)
  real    :: mgfluxz(GRID_ILO_GC:GRID_IHI_GC,GRID_JLO_GC:GRID_JHI_GC,GRID_KLO_GC:GRID_KHI_GC, & 
       &                   MAXBLOCKS)

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


  if (first_call) then
     first_call = .false.
     myPE2 = gr_meshMe
     num_PEs = gr_meshNumProcs
     flux_x_ptr => flux_x
     flux_y_ptr => flux_y
     flux_z_ptr => flux_z
  end if
  call Grid_getLocalNumBlks(lnblocks2)


  flux_x_ptr(FLXBC_FLUX,:,:,:,:) = 0.
  flux_y_ptr(FLXBC_FLUX,:,:,:,:) = 0.
  flux_z_ptr(FLXBC_FLUX,:,:,:,:) = 0.

  update_fluxes(:) = .false.
  do lb = 1, lnblocks2     
    update_fluxes(lb) = (nodetype(lb) .eq. 1)
  enddo

! ***
  call Timers_start("gr_bicgBndry")
  call gr_bicgBndry (irhs,nguard_work,gcellflg)
  call Timers_stop("gr_bicgBndry")

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


           flux_x_ptr(FLXBC_FLUX,1,:,:,lb) = mgfluxx(ili-1,jli:jui,kli:kui,lb)
           flux_x_ptr(FLXBC_FLUX,2,:,:,lb) = mgfluxx(iui,jli:jui,kli:kui,lb)

        endif


     enddo

  elseif (ndim == 2) then

     do lb = 1, lnblocks2

        if (update_fluxes(lb)) then

           ! Get BlockSize:
           call Grid_getBlkPhysicalSize(lb,size)

           dx = size(1)/real(nxb)
           dy = size(2)/real(nyb)

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

           flux_x_ptr(FLXBC_FLUX,1,:,:,lb) = mgfluxx(ili-1,jli:jui,kli:kui,lb)*dy
           flux_x_ptr(FLXBC_FLUX,2,:,:,lb) = mgfluxx(iui,jli:jui,kli:kui,lb)*dy
           flux_y_ptr(FLXBC_FLUX,:,1,:,lb) = mgfluxy(ili:iui,jli-1,kli:kui,lb)*dx
           flux_y_ptr(FLXBC_FLUX,:,2,:,lb) = mgfluxy(ili:iui,jui,kli:kui,lb)*dx

        endif

     enddo
     
  else

     do lb = 1, lnblocks2


        if (update_fluxes(lb)) then

           ! Get BlockSize:
           call Grid_getBlkPhysicalSize(lb,size)

           dxdy =(size(1)*real(nxb)**(-1.)) * (size(2)*real(nyb)**(-1.))
           dydz =(size(2)*real(nyb)**(-1.)) * (size(3)*real(nzb)**(-1.))
           dxdz =(size(1)*real(nxb)**(-1.)) * (size(3)*real(nzb)**(-1.)) 

           delx = nxb/size(1)
           dely = nyb/size(2)
           delz = nzb/size(3)
           do k = kli-1, kui
              do j = jli-1, jui
                 do i = ili-1, iui
                       mgfluxx(i,j,k,lb) = & 
                            &              delx * (work(i+1,j,k,lb,1) - work(i,j,k,lb,1))
                       mgfluxy(i,j,k,lb) = & 
                            &              dely * (work(i,j+1,k,lb,1) - work(i,j,k,lb,1))
                       mgfluxz(i,j,k,lb) = & 
                            &              delz * (work(i,j,k+1,lb,1) - work(i,j,k,lb,1))
                 enddo
              enddo
           enddo
           flux_x_ptr(FLXBC_FLUX,1,:,:,lb) = mgfluxx(ili-1,jli:jui,kli:kui,lb)*dydz
           flux_x_ptr(FLXBC_FLUX,2,:,:,lb) = mgfluxx(iui,jli:jui,kli:kui,lb)*dydz
           flux_y_ptr(FLXBC_FLUX,:,1,:,lb) = mgfluxy(ili:iui,jli-1,kli:kui,lb)*dxdz
           flux_y_ptr(FLXBC_FLUX,:,2,:,lb) = mgfluxy(ili:iui,jui,kli:kui,lb)*dxdz
           flux_z_ptr(FLXBC_FLUX,:,:,1,lb) = mgfluxz(ili:iui,jli:jui,kli-1,lb)*dxdy
           flux_z_ptr(FLXBC_FLUX,:,:,2,lb) = mgfluxz(ili:iui,jli:jui,kui,lb)*dxdy

        endif

     enddo

  endif

  !-------------------------------------------------------------------------------
  !               Now use the PARAMESH flux conservation routine to enforce
  !               continuity of the first derivative of the solution by
  !               overwriting the coarse-grid boundary flux at coarse-fine
  !               interfaces with the fine-grid boundary flux.
  !               Call the PARAMESH flux conservation routine.
  call Timers_start("amr_flux_conserve")
  call amr_flux_conserve(MyPE2, 0, 0)
  call Timers_stop("amr_flux_conserve")

  if (ndim .eq. 1) then
     do i = 1, lnblocks2
        mgfluxx(ili-1,jli:jui,kli:kui,i) = flux_x_ptr(FLXBC_FLUX,1,:,:,i)
        mgfluxx(iui,jli:jui,kli:kui,i)   = flux_x_ptr(FLXBC_FLUX,2,:,:,i)
     enddo
  elseif (ndim .eq. 2) then
     do i = 1, lnblocks2

        ! Get BlockSize:
        call Grid_getBlkPhysicalSize(i,size)

        dx = real(nxb)/size(1)
        dy = real(nyb)/size(2)

        mgfluxx(ili-1,jli:jui,kli:kui,i) = flux_x_ptr(FLXBC_FLUX,1,:,:,i)*dy
        mgfluxx(iui,jli:jui,kli:kui,i)   = flux_x_ptr(FLXBC_FLUX,2,:,:,i)*dy
        mgfluxy(ili:iui,jli-1,kli:kui,i) = flux_y_ptr(FLXBC_FLUX,:,1,:,i)*dx
        mgfluxy(ili:iui,jui,kli:kui,i)   = flux_y_ptr(FLXBC_FLUX,:,2,:,i)*dx
                  
     enddo

  elseif (ndim .eq. 3) then
     do i = 1, lnblocks2

        ! Get BlockSize:
        call Grid_getBlkPhysicalSize(i,size)

        dxdy =(real(nxb)*size(1)**(-1.)) * (real(nyb)*size(2)**(-1.))
        dydz =(real(nyb)*size(2)**(-1.)) * (real(nzb)*size(3)**(-1.))
        dxdz =(real(nxb)*size(1)**(-1.)) * (real(nzb)*size(3)**(-1.))


        mgfluxx(ili-1,jli:jui,kli:kui,i) = flux_x_ptr(FLXBC_FLUX,1,:,:,i)*dydz
        mgfluxx(iui,jli:jui,kli:kui,i)   = flux_x_ptr(FLXBC_FLUX,2,:,:,i)*dydz
        mgfluxy(ili:iui,jli-1,kli:kui,i) = flux_y_ptr(FLXBC_FLUX,:,1,:,i)*dxdz
        mgfluxy(ili:iui,jui,kli:kui,i)   = flux_y_ptr(FLXBC_FLUX,:,2,:,i)*dxdz
        mgfluxz(ili:iui,jli:jui,kli-1,i) = flux_z_ptr(FLXBC_FLUX,:,:,1,i)*dxdy
        mgfluxz(ili:iui,jli:jui,kui,i)   = flux_z_ptr(FLXBC_FLUX,:,:,2,i)*dxdy
     enddo
  endif

  !-------------------------------------------------------------------------------
  !               Now difference the fluxes and subtract from the RHS to
  !               obtain the residual.

  if (ndim == 1) then

     do lb = 1, lnblocks2
        if (update_fluxes(lb)) then

           ! Get BlockSize:
           call Grid_getBlkPhysicalSize(lb,size)

           ! Point to blocks center vars:
           call Grid_getBlkPtr(lb,unkt,CENTER)      

           delx = nxb/size(1)
           do i = ili, iui
                 unkt(ilhs,i,jli:jui,kli:kui) = &  
                      &          delx * (mgfluxx(i,jli:jui,kli:kui,lb) - & 
                      &                  mgfluxx(i-1,jli:jui,kli:kui,lb))
           enddo
              
           ! Point to blocks center vars:
           call Grid_releaseBlkPtr(lb,unkt,CENTER)  

        endif
     enddo

  elseif (ndim == 2) then

     do lb = 1, lnblocks2
        if (update_fluxes(lb)) then

           ! Get BlockSize:
           call Grid_getBlkPhysicalSize(lb,size)

           ! Point to blocks center vars:
           call Grid_getBlkPtr(lb,unkt,CENTER)      

           delx = real(nxb)/size(1)
           dely = real(nyb)/size(2)
           do j = jli, jui
              do i = ili, iui
                    unkt(ilhs,i,j,kli:kui) = & 
                         &            delx * (mgfluxx(i,j,kli:kui,lb) - & 
                         &                    mgfluxx(i-1,j,kli:kui,lb)) + & 
                         &            dely * (mgfluxy(i,j,kli:kui,lb) - & 
                         &                    mgfluxy(i,j-1,kli:kui,lb))
              enddo
           enddo

           ! Point to blocks center vars:
           call Grid_releaseBlkPtr(lb,unkt,CENTER) 

        endif
     enddo

  else

     do lb = 1, lnblocks2
        if (update_fluxes(lb)) then


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
                    unkt(ilhs,i,j,k) = & 
                        &              delx * (mgfluxx(i,j,k,lb) - mgfluxx(i-1,j,k,lb)) + & 
                        &              dely * (mgfluxy(i,j,k,lb) - mgfluxy(i,j-1,k,lb)) + & 
                        &              delz * (mgfluxz(i,j,k,lb) - mgfluxz(i,j,k-1,lb))
                 enddo
              enddo
           enddo

           ! Point to blocks center vars:
           call Grid_releaseBlkPtr(lb,unkt,CENTER) 

        endif
     enddo

  endif

     !===============================================================================
     
  return
end subroutine gr_bicgMultiAx



