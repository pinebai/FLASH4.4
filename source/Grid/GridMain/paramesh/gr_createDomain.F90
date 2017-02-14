!!****if* source/Grid/GridMain/paramesh/gr_createDomain
!!
!! NAME
!!
!!  gr_createDomain
!!
!! 
!! SYNOPSIS
!!
!!  gr_createDomain()
!!
!!
!! DESCRIPTION
!!
!!  Construct the top-level block structure.  This version of
!!  gr_createDomain() sets up an Nblockx * Nblocky * Nblockz array
!!  of top-level blocks.
!!
!! ARGUMENTS
!!
!!
!!***

#include "Flash.h"

#ifdef DEBUG_ALL
#define DEBUG_GRID
#ifndef FLASH_GRID_PARAMESH2
#define DEBUG_GRID_PARAMESH3
#endif
#endif

!!
!!  NOTE: Support for removing floating point bias from the mesh discretization
!!  (controlled by unbiased_geometry switch) is not provided. This option is
!!  not recommended for general use but was found helpful in some cases.
!!  Originally suggested by Artur Gawryszczak, Copernicus Center,  Warsaw.


subroutine gr_createDomain()
  use tree, ONLY : lnblocks,neigh, mfaces,lrefine_max, bnd_box
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_data, ONLY : gr_nBlockX, gr_nBlockY, gr_nBlockZ,gr_domainBC,&
                        gr_imin,gr_imax,gr_jmin,gr_jmax,gr_kmin,gr_kmax, &
                        gr_meshMe, gr_meshNumProcs
  use Simulation_interface, ONLY : Simulation_defineDomain

#ifndef FLASH_GRID_PARAMESH2
  use tree, ONLY : boundary_box, boundary_index,bsize,coord, mboundaries,nboundaries
  use gr_interface, ONLY : gr_createBlock, gr_packBCs
#endif

  implicit none
#include "constants.h"
#ifndef NBOUNDARIES
#define NBOUNDARIES 2*N_DIM
#endif


  integer :: i, j, k,n, i1, j1, k1,n1,m1,ii,s
  

  integer :: blockID, blockID1, mortIx ! 0-based
  integer :: myBlk ! 1-based
  
  integer :: gr_getIndex
  
  integer :: initBlocks,neighborBlock
  integer :: startblockID, endblockID ! 0-based
  integer, dimension(MDIM) :: pos,pos1,nblks
  logical, allocatable,dimension(:,:,:) :: initialDomain
  logical, allocatable :: mortonMask(:)
  integer, allocatable,dimension(:,:,:,:) :: boundaries
  real,    allocatable,dimension(:) :: coordI,coordJ,coordK
  integer :: bbox
  
!==============================================================================

!            Create the blocks, setting their coordinates and neighbor
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

  initBlocks = count(initialDomain)
  allocate(mortonMask(0:product(nblks)-1))
  
  ! compute mortonMask as initialDomain permuted to be in morton order
  do k=1, gr_nBlockZ
  do j=1, gr_nBlockY
  do i=1, gr_nBlockX
     pos(IAXIS) = i; pos(JAXIS) = j; pos(KAXIS) = k
     blockID = morton_ijk_to_index(nblks, pos-1)
     mortonMask(blockID) = initialDomain(i,j,k)
  end do
  end do
  end do

  startblockID = proc_to_block(gr_meshNumProcs, initBlocks, gr_meshMe)
  endblockID = proc_to_block(gr_meshNumProcs, initBlocks, gr_meshMe+1) - 1 ! dont worry if gr_meshMe+1 doesnt exist

  lnblocks = 0
  bbox = 2*NDIM+1

  if(startblockID <= endblockID) then
     allocate(coordI(0:gr_nBlockX))
     allocate(coordJ(0:gr_nBlockY))
     allocate(coordK(0:gr_nBlockZ))

     do i = 1,gr_nBlockX-1
        coordI(i) = ( gr_imin*(gr_nBlockX-i)  + gr_imax * i ) / real(gr_nBlockX)
     end do
     coordI(0) = gr_imin; coordI(gr_nBlockX) = gr_imax

     do j = 1,gr_nBlockY-1
        coordJ(j) = ( gr_jmin*(gr_nBlockY-j)  + gr_jmax * j ) / real(gr_nBlockY)
     end do
     coordJ(0) = gr_jmin; coordJ(gr_nBlockY) = gr_jmax

     do k = 1,gr_nBlockZ-1
        coordK(k) = ( gr_kmin*(gr_nBlockZ-k)  + gr_kmax * k ) / real(gr_nBlockZ)
     end do
     coordK(0) = gr_kmin; coordK(gr_nBlockZ) = gr_kmax
     
     ! Loop over all root blocks, act on only those assigned to this proc.
     do k=1, gr_nBlockZ
        do j=1, gr_nBlockY
           do i=1, gr_nBlockX
              if(initialDomain(i,j,k)) then
                 pos(IAXIS) = i; pos(JAXIS) = j; pos(KAXIS) = k
                 mortIx = morton_ijk_to_index(nblks, pos-1)
                 blockID = count(mortonMask(0:mortIx-1))
                 
                 if((startblockID <= blockID).and.(blockID <= endblockID)) then
                    myBlk = blockID - startblockID + 1
#ifdef DEBUG_GRID
                    print*,'I am creating block',gr_meshMe,myBlk,i,j,k
#endif
                    call gr_createBlock(coordI(i-1), coordI(i), &
                                        coordJ(j-1), coordJ(j), &
                                        coordK(k-1), coordK(k), &
                                        myBlk)
#ifndef FLASH_GRID_PARAMESH2
                    ! createBlock uses coordinate subtraction to compute widths
                    ! which is bad, lets just fix it here so the interface doesn't
                    ! have to change.
                    ! Actually, we aren't going through with this.  Leaving it
                    ! as commented out.
                    !bsize(1,myBlk) = (gr_imax-gr_imin)/real(nblks(1))
                    !bsize(2,myBlk) = (gr_jmax-gr_jmin)/real(nblks(2))
                    !bsize(3,myBlk) = (gr_kmax-gr_kmin)/real(nblks(3))
#endif

#ifdef DEBUG_GRID_PARAMESH3
                    print*,bsize(1:NDIM,myBlk),coord(1:NDIM,myBlk)
#endif
                    !! Now figure out the neighbors
                    !! If a block has index "1" along any dimension, then 
                    !! it is on the lowerface physical boundary 
                    !! along that dimension
                    !! and if its integer index is X/Y/Zblock 
                    !! then it is on the upperface
                    !! physical boundary of the block
                    do ii = 1,NDIM
                       do s=1, 2 ! LOW side, HIGH side
                          pos1 = pos
                          pos1(ii) = pos(ii) + merge(-1,1,s==1)
                          if(gr_domainBC(merge(LOW,HIGH,s==1),ii) == PERIODIC) then
                             pos1(ii) = 1 + modulo(pos1(ii)-1, nblks(ii))
                          end if
                          
                          n1 = 2*(ii-1) + s
                          if(pos1(ii) < 1 .or. nblks(ii) < pos1(ii)) then
                             neigh(:,n1,myBlk) = gr_domainBC(merge(LOW,HIGH,s==1),ii)
                          else
                             i1 = pos1(IAXIS); j1 = pos1(JAXIS); k1 = pos1(KAXIS)
                             if(initialDomain(i1,j1,k1)) then
                                mortIx = morton_ijk_to_index(nblks, pos1-1)
                                blockID1 = count(mortonMask(0:mortIx-1))
                                neigh(2,n1,myBlk) = block_to_proc(gr_meshNumProcs, initBlocks, blockID1)
                                neigh(1,n1,myBlk) = 1 + blockID1 - proc_to_block( &
                                  gr_meshNumProcs, initBlocks, neigh(2,n1,myBlk) &
                                )
                             else
                                neigh(:,n1,myBlk) = boundaries(n1-merge(-1,1,s==1),i1,j1,k1)
                             end if
                          end if
                       end do
                    end do
#ifdef DEBUG_GRID
                    print*,'for myblk',myblk,blockID
                    print*,neigh(1,:,myBlk)
#endif
                 end if
              else
#ifndef FLASH_GRID_PARAMESH2
                 if (bbox .gt. NBOUNDARIES) then
                    call Driver_abortFlash('Too many boundary conditions - increase NBOUNDARIES!')
                 end if
                 if (bbox .gt. nboundaries) then
                    call Driver_abortFlash('Too many boundary conditions, found PARAMESH nboundaries < NBOUNDARIES!')
                 end if
                 if (bbox .gt. mboundaries) then
                    call Driver_abortFlash('Too many boundary conditions, found PARAMESH mboundaries < NBOUNDARIES!')
                 end if
                 boundary_box(1,1,bbox)=coordI(i-1)
                 boundary_box(2,1,bbox)=coordI(i)
                 if(NDIM>1) then
                    boundary_box(1,2,bbox)=coordJ(j-1)
                    boundary_box(2,2,bbox)=coordJ(j)
                 end if
                 if(NDIM>2) then
                    boundary_box(1,3,bbox)=coordK(k-1)
                    boundary_box(2,3,bbox)=coordK(k)
                 end if
                 boundary_index(bbox)=gr_packBCs( &
                      boundaries(1,i,j,k) , &
                      boundaries(2,i,j,k) , &
                      boundaries(3,i,j,k) , &
                      boundaries(4,i,j,k) , &
                      boundaries(5,i,j,k) , &
                      boundaries(6,i,j,k) )
                 bbox = bbox + 1
#endif
              end if
           end do
        end do
     end do
     deallocate(coordI)
     deallocate(coordJ)
     deallocate(coordK)
  end if
  
  deallocate(initialDomain)
  deallocate(boundaries)
  deallocate(mortonMask)
  
  return
contains

  ! given global block id, compute owning processor
  ! np: global num procs
  ! nb: global num blocks
  ! b: 0-based global block id
  ! p: 0-based proc id
  function block_to_proc(np, nb, b) result(p)
    implicit none
    integer, intent(in) :: np, nb, b
    integer :: p, cut
    cut = mod(nb,np)*(1 + nb/np)
    if(b < cut) then
      p = b/(nb/np + 1)
    else
      p = (b-cut)/(nb/np) + mod(nb,np)
    end if
  end function
  
  ! given processor, compute first global block id it owns
  ! np: global num procs
  ! nb: global num blocks
  ! p: 0-based proc id, will work if called with np
  ! b0: 0-based global block id
  function proc_to_block(np, nb, p) result(b0)
    implicit none
    integer, intent(in) :: np, nb, p
    integer :: b0
    b0 = p*(nb/np) + min(p, mod(nb,np))
  end function
  
  ! rect: number of blocks along each dimension
  ! ijk: 0-based index coordinates into domain
  ! mort_ix: 0-based morton curve index
  function morton_ijk_to_index(rect, ijk) result(mort_ix)
    implicit none
    integer, intent(in) :: rect(MDIM), ijk(MDIM)
    integer :: mort_ix
    integer :: x(NDIM), box(NDIM)
    integer :: max_pow2, max_d, pow2, d
    x(:) = ijk(1:NDIM)
    box(:) = rect(1:NDIM)
    mort_ix = 0
    ! bisect box until it's just one element
    do while(.true.)
      ! find dim that can fit biggest pow2 strictly inside box
      max_pow2 = 0
      max_d = 0
      do d=1, NDIM
        pow2 = glb_pow2(box(d) - 1)
        if(pow2 >= max_pow2) then
          max_pow2 = pow2
          max_d = d
        end if
      end do
      ! found max_d, bisect with it
      if(max_pow2 == 0) then
        exit ! box is singular
      else if(x(max_d) < max_pow2) then
        box(max_d) = max_pow2
      else
        mort_ix = mort_ix + product((/(merge(max_pow2, box(I), I==max_d),I=1,NDIM)/))
        x(max_d) = x(max_d) - max_pow2
        box(max_d) = box(max_d) - max_pow2
      end if
    end do
  end function
  
  ! greatest lower bounding power of 2
  ! returns p2, the greatest power of 2 such that p2 <= x
  function glb_pow2(x) result(p2)
    implicit none
    integer, intent(in) :: x
    integer :: p2, i
    i = 1
    p2 = x
    do while(i < bit_size(x)) ! this should unroll
      p2 = ior(p2, ishft(p2,-i))
      i = ishft(i,1)
    end do
    p2 = p2 - ishft(p2,-1)
  end function
end subroutine gr_createDomain
