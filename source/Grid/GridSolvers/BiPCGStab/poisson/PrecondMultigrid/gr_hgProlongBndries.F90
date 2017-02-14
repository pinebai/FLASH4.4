!!****if* source/Grid/GridSolvers/BiPCGStab/poisson/PrecondMultigrid/gr_hgProlongBndries
!!
!! NAME
!!  gr_hgProlongBndries
!!
!! SYNOPSIS
!!  
!!  gr_hgProlongBndries(integer, intent(in) :: level
!!                      integer, intent(in) :: ifrom
!!                      integer, intent(in) :: ito
!!                      integer, intent(in) :: ichild)
!! DESCRIPTION
!!
!!  Prolong the boundary potentials present on a coarse level onto the
!!  face-guardcells of the children blocks at the next level.
!!  This occurs both on the faces and interiors of the parent blocks.
!!  For the faces, one interpolates the surface boundary values from the
!!  first layer of guardcells on the parent block.  For the second, one
!!  interpolates from the interior to get the approximated value at the
!!  bisection of the parent block.  These are put on the first layer of
!!  the child's boundary.  This scheme allows for seamless handling of
!!  external boundary values.  This routine is used in conjunction with a
!!  single-block solver to take an approximate solution on the coarse level
!!  and put it onto the fine levels.
!!  
!! ARGUMENTS
!!   
!!   ifrom  - the grid variable the data on the parent comes from
!!   ito    - the grid variable into which the face values are put
!!   level  - the level of the parents
!!   ichild - allows for only a single child block index to be filled
!!
!! NOTES
!!  
!!  Due to polynomial approximation's tendency to oscillate, the guardcell
!!  fills done here have two behaviors on the exterior blocks.  For
!!  the source term, the guardcells are extrapolated outwards.  For
!!  the residual terms, the guardcells are assumed 0 at the faces.
!! 
!! SEE ALSO
!!  
!!  gr_hgRestrict
!!
!!***

!!REORDER(5): unk

!==============================================================================

subroutine gr_hgProlongBndries(level, ifrom, ito, ichild)

#include "Flash.h"
#include "Multigrid.h"
#include "constants.h"

  use gr_hgData, ONLY: hg_myPE, nbbuf_prolong, nmax1, nmax2, n1off, n2off, n3off,  &
       & Pns, Px, Py, Pz, Pud, Pew, send_prolong_req, send_prolong_data, recv_prolong_data, &
       & gr_hgSaveNodetype, gr_hgSolnIndex
  use Driver_interface, ONLY : Driver_abortFlash
  use tree, ONLY : lnblocks,parent,child,neigh,nfaces,nchild,nodetype,lrefine
  use physicaldata, ONLY : unk
  use workspace, ONLY : work
  use Grid_data, ONLY : gr_meshComm
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use gr_hgInterface, ONLY: gr_hgBndry

  implicit none
  
#include "Flash_mpi.h"
  
  integer, intent(in)          :: ifrom, ito, level, ichild
  
  integer                      :: b, c, h, i, j, k, ierr1, ierr2, ierr3
  integer                      :: i1, i2, j1, j2, ii, jj, kk
  integer                      :: ierr, nsent
  integer                      :: blockID
  logical                      :: any_sent
  integer                      :: status(MPI_STATUS_SIZE)
  integer                      :: send_status(MPI_STATUS_SIZE,nchild*nbbuf_prolong)


!==============================================================================


  call Timers_start("gr_hgProlongBndries")

  ! Use (EXCHANGE_WORK, CONTINUE_SERIES) under the assumption that we're
  ! calling gr_hgProlongBndries immediately after a call to hg_solveLevel
  ! or hg_residual -- performance savings by not filling work() again
  

!!$  if (ifrom == gr_hgSolnIndex) then
!!$    call gr_hgBndry(level, ifrom, NGUARD, 0, MG_EXCHANGE_WORK, &
!!$                    MG_STANDALONE, .true.)  ! do extrapolation across the boundary
!!$  else
!!$    call gr_hgBndry(level, ifrom, NGUARD, 0, MG_EXCHANGE_WORK, &
!!$                    MG_STANDALONE, .false.)
!!$  endif

  nsent = 0
  h = 1
  do b=1,lnblocks
     ! gr_hgSaveNodetype is a copy of the original (pre-MG) node types.
     !  So this next line searches for blocks with children at this level
     if((lrefine(b)==level).and.(gr_hgSaveNodetype(b)>=PARENT_BLK)) then
        
        ! First, compute interpolated boundary values for each of this block's
        ! children.

        if (NDIM == 1) then
           
           do c = 1, nchild  ! over number of children from paramesh tree
              if ((ichild == 0) .or. (ichild == c)) then
                 
                 send_prolong_data(:,:,:,c,h) = 0.
                 
                 do ii = -2, 2
                    send_prolong_data(1,1,1,c,h) = &
                         send_prolong_data(1,1,1,c,h) + &
                         Pew(ii,1)*work(NGUARD+n1off(c)+1+ii,1,1,b,1)
                         ! Pew! is smelly prolongation east-west
                    send_prolong_data(1,1,2,c,h) = &
                         send_prolong_data(1,1,2,c,h) + &
                         Pew(ii,2)*work(NGUARD+n1off(c)+NXB/2+ii,1,1,b,1)
                 enddo
              
              endif
           enddo
        
        else if (NDIM == 2) then
        
           do c = 1, nchild
              if ((ichild == 0) .or. (ichild == c)) then
                 
                 send_prolong_data(:,:,:,c,h) = 0.
              
              ! N & S edges
              
                 do i = 1, NXB
                    do ii = -2, 2
                       do jj = -2, 2
                          send_prolong_data(i,1,3,c,h) = &
                               send_prolong_data(i,1,3,c,h) + &
                               Px(ii,i,c)*Pns(jj,1)*&
                               work(NGUARD+n1off(c)+1+(i-1)/2+ii,NGUARD+n2off(c)+1+jj,1,b,1)
                          send_prolong_data(i,1,4,c,h) = &
                               send_prolong_data(i,1,4,c,h) + &
                               Px(ii,i,c)*Pns(jj,2)*&
                               work(NGUARD+n1off(c)+1+(i-1)/2+ii,NGUARD+n2off(c)+NYB/2+jj,1,b,1)
                       enddo
                    enddo
                 enddo
                 
                 ! E & W edges
                 
                 do j = 1, NYB
                    do jj = -2, 2
                       do ii = -2, 2
                          send_prolong_data(j,1,1,c,h) = &
                               send_prolong_data(j,1,1,c,h) + &
                               Py(jj,j,c)*Pew(ii,1)*&
                               work(NGUARD+n1off(c)+1+ii,NGUARD+n2off(c)+1+(j-1)/2+jj,1,b,1)
                          send_prolong_data(j,1,2,c,h) = &
                               send_prolong_data(j,1,2,c,h) + &
                               Py(jj,j,c)*Pew(ii,2)*&
                               work(NGUARD+n1off(c)+NXB/2+ii,NGUARD+n2off(c)+1+(j-1)/2+jj,1,b,1)
                       enddo
                    enddo
                 enddo
                 
              endif
           enddo
           
        else ! NDIM == 3
           
           do c = 1, nchild
              if ((ichild == 0) .or. (ichild == c)) then
                 
                 send_prolong_data(:,:,:,c,h) = 0.
                 
                 ! N & S faces
                 
                 do k = 1, NZB
                    do i = 1, NXB
                       do kk = -2, 2
                          do ii = -2, 2
                             do jj = -2, 2
                                send_prolong_data(i,k,3,c,h) = send_prolong_data(i,k,3,c,h) + &
                                     Px(ii,i,c)*Pz(kk,k,c)*Pns(jj,1)* &
                                     work(NGUARD+n1off(c)+1+(i-1)/2+ii, &
                                     NGUARD+n2off(c)+1+jj, &
                                     NGUARD+n3off(c)+1+(k-1)/2+kk,b,1)
                                send_prolong_data(i,k,4,c,h) = send_prolong_data(i,k,4,c,h) + &
                                     Px(ii,i,c)*Pz(kk,k,c)*Pns(jj,2)* &
                                     work(NGUARD+n1off(c)+1+(i-1)/2+ii, &
                                     NGUARD+n2off(c)+NYB/2+jj, &
                                     NGUARD+n3off(c)+1+(k-1)/2+kk,b,1)
                             enddo
                          enddo
                       enddo
                    enddo
                 enddo
                 
                 ! E & W faces
                 
                 do k = 1, NZB
                    do j = 1, NYB
                       do kk = -2, 2
                          do jj = -2, 2
                             do ii = -2, 2
                                send_prolong_data(j,k,1,c,h) = send_prolong_data(j,k,1,c,h) + &
                                     Py(jj,j,c)*Pz(kk,k,c)*Pew(ii,1)* &
                                     work(NGUARD+n1off(c)+1+ii, &
                                     NGUARD+n2off(c)+1+(j-1)/2+jj, &
                                     NGUARD+n3off(c)+1+(k-1)/2+kk,b,1)
                                send_prolong_data(j,k,2,c,h) = send_prolong_data(j,k,2,c,h) + &
                                     Py(jj,j,c)*Pz(kk,k,c)*Pew(ii,2)* &
                                     work(NGUARD+n1off(c)+NXB/2+ii, &
                                     NGUARD+n2off(c)+1+(j-1)/2+jj, &
                                     NGUARD+n3off(c)+1+(k-1)/2+kk,b,1)
                             enddo
                          enddo
                       enddo
                    enddo
                 enddo
                 
                 ! U & D faces
                 
                 do j = 1, NYB
                    do i = 1, NXB
                       do jj = -2, 2
                          do ii = -2, 2
                             do kk = -2, 2
                                send_prolong_data(i,j,5,c,h) = send_prolong_data(i,j,5,c,h) + &
                                     Px(ii,i,c)*Py(jj,j,c)*Pud(kk,1)* &
                                     work(NGUARD+n1off(c)+1+(i-1)/2+ii, &
                                     NGUARD+n2off(c)+1+(j-1)/2+jj, &
                                     NGUARD+n3off(c)+1+kk,b,1)
                                send_prolong_data(i,j,6,c,h) = send_prolong_data(i,j,6,c,h) + &
                                     Px(ii,i,c)*Py(jj,j,c)*Pud(kk,2)* &
                                     work(NGUARD+n1off(c)+1+(i-1)/2+ii, &
                                     NGUARD+n2off(c)+1+(j-1)/2+jj, &
                                     NGUARD+n3off(c)+NZB/2+kk,b,1)
                             enddo
                          enddo
                       enddo
                    enddo
                 enddo
                 
              endif
           enddo
           
        endif
     
     ! Next, loop over children.  If a child is on this processor, set its
     ! boundary values directly.  If it is off-processor, send the data to
     ! the owning processor (non-blocking).
     
       any_sent = .false.
     
       do c = 1, nchild
          if ((ichild == 0) .or. (ichild == c)) then
               if (child(2,c,b) == hg_myPE) then  ! local child
                  blockID=child(1,c,b)
                  unk(ito,NGUARD,NGUARD*K2D+1:NGUARD*K2D+NYB,&
                     NGUARD*K3D+1:NGUARD*K3D+NZB,blockID) = &
                     send_prolong_data(1:NYB,1:NZB,1,c,h)
                unk(ito,NGUARD+NXB+1,NGUARD*K2D+1:NGUARD*K2D+NYB,&
                     NGUARD*K3D+1:NGUARD*K3D+NZB,blockID) = &
                     send_prolong_data(1:NYB,1:NZB,2,c,h)
                if (NDIM >= 2) then
                   unk(ito,NGUARD+1:NGUARD+NXB,((NGUARD-1)*K2D)+1,&
                        NGUARD*K3D+1:NGUARD*K3D+NZB,blockID) = &
                        send_prolong_data(1:NXB,1:NZB,3,c,h)
                   unk(ito,NGUARD+1:NGUARD+NXB,((NGUARD+NYB)*K2D)+1,&
                        NGUARD*K3D+1:NGUARD*K3D+NZB,blockID) = &
                        send_prolong_data(1:NXB,1:NZB,4,c,h)
                endif
                if (NDIM == 3) then
                   unk(ito,NGUARD+1:NGUARD+NXB,(NGUARD*K2D)+1:(NGUARD*K2D)+NYB,&
                        ((NGUARD-1)*K3D)+1,blockID) = send_prolong_data(1:NXB,1:NYB,5,c,h)
                   unk(ito,NGUARD+1:NGUARD+NXB,(NGUARD*K3D)+1:(NGUARD*K3D)+NYB,&
                        ((NGUARD+NZB)*K3D)+1,blockID) = send_prolong_data(1:NXB,1:NYB,6,c,h)
                endif
             else                            ! remote child
                any_sent = .true.
                nsent = nsent + 1
                call mpi_issend(send_prolong_data(1,1,1,c,h), nmax1*nmax2*nfaces, &
                     FLASH_REAL, child(2,c,b), child(1,c,b), &
                     gr_meshComm, send_prolong_req(nsent), ierr)
             endif
          endif
       enddo
      
       if (any_sent) h = h + 1
       if (h > nbbuf_prolong) &
            call Driver_abortFlash("Buffer space exceeded in gr_hgProlongBndries")
     end if
  enddo
  
  do b = 1, lnblocks
     if ((lrefine(b) == level+1) .and. (parent(2,b) /= hg_myPE)) then
        
        ! If parent is on another processor, receive the boundary data from
        ! the parent (blocking).
      
        
        call mpi_recv(recv_prolong_data(1,1,1), nmax1*nmax2*nfaces, &
             FLASH_REAL, parent(2,b), b, &
             gr_meshComm, status, ierr)
        unk(ito,NGUARD,NGUARD*K2D+1:NGUARD*K2D+NYB,&
             NGUARD*K3D+1:NGUARD*K3D+NZB,b) = &
             recv_prolong_data(1:NYB,1:NZB,1)
        unk(ito,NGUARD+NXB+1,NGUARD*K2D+1:NGUARD*K2D+NYB,&
             NGUARD*K3D+1:NGUARD*K3D+NZB,b) = &
             recv_prolong_data(1:NYB,1:NZB,2)
        if (NDIM >= 2) then
           unk(ito,NGUARD+1:NGUARD+NXB,(NGUARD-1)*K2D+1,&
                NGUARD*K3D+1:NGUARD*K3D+NZB,b) = &
                recv_prolong_data(1:NXB,1:NZB,3)
           unk(ito,NGUARD+1:NGUARD+NXB,(NGUARD+NYB)*K2D+1,&
                NGUARD*K3D+1:NGUARD*K3D+NZB,b) = &
                recv_prolong_data(1:NXB,1:NZB,4)
        endif
        if (NDIM == 3) then
           unk(ito,NGUARD+1:NGUARD+NXB,(NGUARD*K2D)+1:(NGUARD*K2D)+NYB,&
                ((NGUARD-1)*K3D)+1,b)           = recv_prolong_data(1:NXB,1:NYB,5)
           unk(ito,NGUARD+1:NGUARD+NXB,(NGUARD*K2D)+1:(NGUARD*K2D)+NYB,&
                ((NGUARD+NZB)*K3D)+1,b)     = recv_prolong_data(1:NXB,1:NYB,6)
        endif
        
        
     endif
  enddo
  
  call mpi_waitall(nsent, send_prolong_req, send_status, ierr)
  call Timers_stop("gr_hgProlongBndries")
  
  !===================================================================

  return
end subroutine gr_hgProlongBndries
