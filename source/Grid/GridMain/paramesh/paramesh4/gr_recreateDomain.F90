!!****if* source/Grid/GridMain/paramesh/paramesh4/gr_recreateDomain
!!
!! NAME
!!
!!  gr_recreateDomain
!!
!! 
!! SYNOPSIS
!!
!!  call gr_recreateDomain()
!!
!!
!! DESCRIPTION
!!
!!  Construct the top-level block structure.  Normally, when the Grid
!!  is created from scratch,   
!!  gr_createDomain() sets up an Nblockx * Nblocky * Nblockz array
!!  of top-level blocks.
!!  This is a pared-down variant of gr_createDomain appropriate when
!!  restarting from a checkpoint.  It is only really needed in simulations
!!  where Simulation_defineDomain is used to add obstacle blocks.
!!
!!
!! ARGUMENTS
!!
!!  NONE
!!
!! NOTES
!! 
!!  The user-defined Simulation_defineDomain should define the obstacles
!!  in the same ways as in previous runs up to the checkpoint from which
!!  the simulation is restarting. If Simulation_defineDomain is changed
!!  in mid-run, the resulting set of obstacle blocks may not be consistent
!!  any more with the properties, especially neighbors, of blocks read
!!  from the checkpoint file, and all bets are off.
!!
!! HISTORY
!!
!!  Started as a pared-down copy of gr_createDomain  2009-05-07  KW
!!***

#include "Flash.h"

#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif
   
!!
!!  NOTE: Support for removing floating point bias from the mesh discretization
!!  (controlled by unbiased_geometry switch) is not provided. This option is
!!  not recommended for general use but was found helpful in some cases.
!!  Originally suggested by Artur Gawryszczak, Copernicus Center,  Warsaw.


subroutine gr_recreateDomain
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_data, ONLY : gr_nBlockX, gr_nBlockY, gr_nBlockZ,&
                        gr_imin,gr_imax,gr_jmin,gr_jmax,gr_kmin,gr_kmax
  use Simulation_interface, ONLY : Simulation_defineDomain

  use tree, ONLY : boundary_box, boundary_index,mboundaries,nboundaries
  use gr_interface, ONLY : gr_packBCs

  implicit none
#include "constants.h"
#ifndef NBOUNDARIES
#define NBOUNDARIES 2*N_DIM
#endif

  real :: dx, dy, dz
  integer :: i, j, k
  


  integer, dimension(MDIM) :: nblks
  logical, allocatable,dimension(:,:,:) :: initialDomain
  integer, allocatable,dimension(:,:,:,:) :: boundaries
  integer :: bbox

!==============================================================================

!            Recreate the blocks, setting their coordinates and neighbor
!            information.


!!..initialize, calculate the number of initial blocks
  if((NDIM < 3).and.(gr_nBlockZ /=1 )) then
     print*,'Warning: setting Nblockz to 1 since not a 3d problem, you specified :',gr_nBlockZ
     gr_nBlockZ=1
  end if
  if((NDIM < 2).and.(gr_nBlockY /=1 )) then
     print*,'Warning : setting NblockY to 1 for 1D problem, you specified :',gr_nBlockY
     gr_nBlockY=1
  end if
  nblks(IAXIS)=gr_nBlockX
  nblks(JAXIS)=gr_nBlockY
  nblks(KAXIS)=gr_nBlockZ

  
  allocate(initialDomain(gr_nBlockX,gr_nBlocky,gr_nBlockZ))
  allocate(boundaries(MDIM*2,gr_nblockX,gr_nBlockY,gr_nBlockZ))

  call Simulation_defineDomain(initialDomain,boundaries,nblks)

  

  bbox = 2*NDIM+1
     
     
  dx = gr_imax - gr_imin
  dy = gr_jmax - gr_jmin
  dz = gr_kmax - gr_kmin

  !..set the deltas and loop ending points based on the number of dimensions

  dx   = dx / float(gr_nBlockX)

  if (NDIM > 1) dy = dy / float(gr_nBlockY)
  if (NDIM == 3)dz = dz / float(gr_nBlockZ)
  do k = 1,gr_nBlockZ
     do j = 1,gr_nBlockY
        do i = 1,gr_nBlockX
           if(.NOT. initialDomain(i,j,k)) then
              if (bbox .gt. NBOUNDARIES) then
                 call Driver_abortFlash('Too many boundary conditions - increase NBOUNDARIES!')
              end if
              if (bbox .gt. nboundaries) then
                 call Driver_abortFlash('Too many boundary conditions, found PARAMESH nboundaries < NBOUNDARIES!')
              end if
              if (bbox .gt. mboundaries) then
                 call Driver_abortFlash('Too many boundary conditions, found PARAMESH mboundaries < NBOUNDARIES!')
              end if
              boundary_box(1,1,bbox)=gr_imin+(i-1)*dx
              boundary_box(2,1,bbox)=gr_imin+i*dx
              boundary_index(bbox)=gr_packBCs( &
                   boundaries(1,i,j,k) , &
                   boundaries(2,i,j,k) , &
                   boundaries(3,i,j,k) , &
                   boundaries(4,i,j,k) , &
                   boundaries(5,i,j,k) , &
                   boundaries(6,i,j,k) )
              if(NDIM>1) then
                 boundary_box(1,2,bbox)=gr_jmin+(j-1)*dy
                 boundary_box(2,2,bbox)=gr_jmin+j*dy
              end if
              if(NDIM>2) then
                 boundary_box(1,3,bbox)=gr_kmin+(k-1)*dz
                 boundary_box(2,3,bbox)=gr_kmin+k*dz
              end if
              bbox = bbox + 1
           end if
        end do
     end do
  end do

  deallocate(initialDomain)
  deallocate(boundaries)
  
  return
end subroutine gr_recreateDomain
