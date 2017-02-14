!!****if* source/Grid/GridSolvers/BiPCGStab/poisson/PrecondMultigrid/gr_hgSolveLevel
!!
!! NAME
!!  gr_hgSolveLevel
!!
!! SYNOPSIS
!!  call gr_hgSolveLevel(integer, intent(in) :: level,
!!                      integer, intent(in) :: gr_iSource,
!!                      integer, intent(in) :: gr_iSoln,
!!                      external            :: SolveBlock,
!!                      integer, intent(in) :: LeafFlag,
!!                      OPTIONAL,real(IN)   :: dt,
!!                      OPTIONAL,real(IN)   :: chi
!!                        OPTIONAL,real(IN)   :: theta)
!! 
!! DESCRIPTION
!!
!!  Block-per-block solve Poisson on all the blocks at level using boundary
!!  conditions on the exterior faces.  The second phase of this routine
!!  does two Jacobi sweeps on the first two interior cells for each block
!!  in order to dampen inter-block discontinuities.
!! 
!! ARGUMENTS
!!  
!!  level        - the level to solve on
!!  gr_iSource - the grid variable holding the source term
!!  gr_iSoln   - the grid variable in which the solution resides
!!  SolveBlock   - a function that can solve (or smooth) locally
!!  LeafFlag     - 0 => all blocks on the level solved
!!                 1 => only leaf blocks solved
!!                 2 => only parent blocks solved
!!  dt         - time step, to be passed down (maybe unused)
!!  chi        - a factor, to be passed down (maybe unused)
!!  theta      - Switch between Implicit/Explicit schemes (maybe unused)
!!
!! NOTES
!!
!!  The Jacobi sweeps SHOULD be replaced with something for better convergence.
!!
!!***

!!REORDER(5): unk

subroutine gr_hgSolveLevel(level, gr_iSource, gr_iSoln, SolveBlock, LeafFlag, dt, chi, theta)

!==================================================================

  use gr_hgData, ONLY: gr_hgBndTypes, &
       hg_ili, hg_iui, hg_jli, hg_jui, hg_kli, hg_kui
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface, ONLY : Grid_getDeltas
  use workspace, ONLY : work
  use tree, ONLY : lnblocks,lrefine,nodetype,bsize
  use physicaldata, ONLY : unk
  use Grid_data, ONLY : gr_meshComm  
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use gr_hgInterface, ONLY: gr_hgBndry

  implicit none
#include "Flash.h"
#include "constants.h"
#include "Multigrid.h"
#include "Flash_mpi.h"

  integer, intent(in) :: gr_iSource, gr_iSoln, level, LeafFlag
  real,intent(IN),OPTIONAL :: dt, chi,theta

  external               SolveBlock

  real, parameter              :: coeff2 = 2.0 !13./6.
  integer                      :: b, i, j, k, n
  integer                      :: nblockssolved
  real, dimension(NXB,NYB,NZB) :: soln
  real, dimension(MDIM)        :: block_size, zone_size, deltas
  real                         :: c, cx, cy, cz
  real                         :: avg, sum, lsum, vol, lvol, nbinv, bvol, cvol, bsum
  logical                      :: SolveThisBlock
  integer                      :: bnd_type, ierr
  
  integer, save                :: lrefine_min

  character                    :: level_timer_label*40

  integer       :: redblackpass, isweep, jsweep, ksweep
  logical, parameter :: relaxflag = .true.
  integer, parameter :: offin = 0
  integer, parameter :: offou = 0

  !====================================================================


  write(unit=level_timer_label, FMT='(A11,I2)') "solvelevel", level
  call Timers_start(level_timer_label)
 
  nblockssolved = 0

  do b = 1, lnblocks
     
     SolveThisBlock = (lrefine(b) == level)
     if (LeafFlag == MG_NODES_LEAF_ONLY) then
        SolveThisBlock = (SolveThisBlock .and. (nodetype(b) == LEAF))
     else if (LeafFlag == MG_NODES_PARENT_ONLY) then
        SolveThisBlock = (SolveThisBlock .and. (nodetype(b) == PARENT_BLK))
     endif

     if (SolveThisBlock) then
        nblockssolved = nblockssolved+1

        !Let's try a standard way of finding delta x rather than this convolution        
        call Grid_getDeltas(b,deltas) 

        do k = 1, NZB
           do j = 1, NYB
              do i = 1, NXB
                 soln(i,j,k) = unk(gr_iSource,i+NGUARD,j+K2D*NGUARD,k+K3D*NGUARD,b)
              enddo
           enddo
        enddo


        bnd_type = gr_hgBndTypes(2*NDIM-1) !low face highest dimension, should really check all directions...
        if (bnd_type == MG_BND_PERIODIC) then
           ! Only the coarsest mesh level is treated as periodic; the rest obtain
           ! boundary values by interpolation and are treated as Dirichlet.
           if (level == 1) then
              bnd_type = MG_BND_PERIODIC
           else
              bnd_type = MG_BND_DIRICHLET
           endif
        else if ((bnd_type == MG_BND_DIRICHLET) .or. &
                 (bnd_type == MG_BND_GIVENVAL)) then
           bnd_type = MG_BND_DIRICHLET
        else if ((bnd_type == MG_BND_NEUMANN)) then
           bnd_type = MG_BND_DIRICHLET
        else
           call Driver_abortFlash("gr_hgSolveLevel found an unrecognized bnd_type!")
        endif
        
        ! Adjust the source function to account for periodic boundary conditions.
        
        !    if (bnd_type == 0) then
        !
        !      avg = sum(soln) / (NXB*NYB*NZB)
        !      soln = soln - avg
        !
        !    endif
        
        ! Adjust the source function to account for given-value boundary
        ! conditions.

        if (bnd_type == MG_BND_DIRICHLET) then           

           soln(1,:,:)   = soln(1,:,:)   - coeff2*unk(gr_iSoln,NGUARD, &
                1+NGUARD*K2D:NYB+NGUARD*K2D, &
                1+NGUARD*K3D:NZB+NGUARD*K3D,b)/deltas(IAXIS)**2
           soln(NXB,:,:) = soln(NXB,:,:) - coeff2*unk(gr_iSoln,NGUARD+NXB+1, &
                1+NGUARD*K2D:NYB+NGUARD*K2D, &
                1+NGUARD*K3D:NZB+NGUARD*K3D,b)/deltas(IAXIS)**2
           
           if (NDIM >= 2) then
              soln(:,1,:)   = soln(:,1,:)   - coeff2*unk(gr_iSoln,NGUARD+1:NGUARD+NXB, &
                   (NGUARD-1)*K2D+1, &
                   1+NGUARD*K3D:NZB+NGUARD*K3D,b)/deltas(JAXIS)**2
              soln(:,NYB,:) = soln(:,NYB,:) - coeff2*unk(gr_iSoln,NGUARD+1:NGUARD+NXB, &
                   (NGUARD+NYB)*K2D+1, &
                   1+NGUARD*K3D:NZB+NGUARD*K3D,b)/deltas(JAXIS)**2
           endif


           if (NDIM == 3) then
              soln(:,:,1)   = soln(:,:,1)   - coeff2*unk(gr_iSoln,NGUARD+1:NGUARD+NXB, &
                   1+NGUARD*K2D:NYB+NGUARD*K2D, &
                   (NGUARD-1)*K3D+1,b )/deltas(KAXIS)**2
              soln(:,:,NZB) = soln(:,:,NZB) - coeff2*unk(gr_iSoln,NGUARD+1:NGUARD+NXB, &
                   1+NGUARD*K3D:NYB+NGUARD*K3D, &
                   (NGUARD+NZB)*K3D+1,b           )/deltas(KAXIS)**2
           endif
      
        endif 
        
        call Timers_start("fft")
        ! LBR Check SolveBlock = gr_hgPoissonSolveBlock
        call SolveBlock (soln, NXB, NYB, NZB, &
                         deltas(IAXIS), deltas(JAXIS), deltas(KAXIS), bnd_type, level)
        call Timers_stop("fft")        

        do k = 1, NZB
           do j = 1, NYB
              do i = 1, NXB
                 unk(gr_iSoln,i+NGUARD,j+K2D*NGUARD,k+K3D*NGUARD,b) = soln(i,j,k)
                 work(i+NGUARD,j+K2D*NGUARD,k+K3D*NGUARD,b,1) = soln(i,j,k)
              enddo
           enddo
        enddo
        
     endif
  enddo
  if (nblockssolved == 0) then
     call Timers_start("fft")   !trick to keep timers structure on different procs the same - KW
     call Timers_stop("fft")
  end if

!  call timer_stop("fft")
  
  !TEST - BOUNDARY RELAXATION
!  call timer_start("relaxation")

  
  do n = 1, 4
     
     if (n == 1) then
        call gr_hgBndry(level, gr_iSoln, 1, 0, MG_EXCHANGE_WORK, MG_BEGIN_SERIES, .false.)
     else
        call gr_hgBndry(level, gr_iSoln, 1, 0, MG_EXCHANGE_WORK, MG_CONTINUE_SERIES, .false.)
     endif
     
     do b = 1, lnblocks
        SolveThisBlock = (lrefine(b) == level)
        if (LeafFlag == 1) then
           SolveThisBlock = (SolveThisBlock .and. (nodetype(b) == 1))
        else if (LeafFlag == 2) then
           SolveThisBlock = (SolveThisBlock .and. (nodetype(b) /= 1))
        endif
        
        if (SolveThisBlock) then
           
           call Grid_getDeltas(b,deltas)
           
           if (NDIM == 1) then
              
              cx = 0.5
              c  = 0.5 * deltas(IAXIS)**2
              i = NGUARD+1
              do i = NGUARD+1, NGUARD+2
                 work(i,1,1,b,1) = cx*(work(i-1,1,1,b,1) + work(i+1,1,1,b,1)) - &
                      c*unk(gr_iSource,i,1,1,b)
                 unk(gr_iSoln,i,1,1,b) = work(i,1,1,b,1)
              enddo
              i = NGUARD+NXB
              do i = NGUARD+NXB-1, NGUARD+NXB
                 work(i,1,1,b,1) = cx*(work(i-1,1,1,b,1) + work(i+1,1,1,b,1)) - &
                      c*unk(gr_iSource,i,1,1,b)
                 unk(gr_iSoln,i,1,1,b) = work(i,1,1,b,1)
              enddo

           else if (NDIM == 2) then
              
              c  = 0.5 / (deltas(IAXIS)**2 + deltas(JAXIS)**2)
              cx = deltas(JAXIS)**2 * c
              cy = deltas(IAXIS)**2 * c
              c  = deltas(IAXIS)**2 * deltas(JAXIS)**2 * c
              
              do j = NGUARD+1, NGUARD+2
                 do i = NGUARD+1, NGUARD+NXB
                    work(i,j,1,b,1) = &
                         cx*(work(i-1,j,1,b,1) + work(i+1,j,1,b,1)) + &
                         cy*(work(i,j-1,1,b,1) + work(i,j+1,1,b,1)) - &
                         c*unk(gr_iSource,i,j,1,b)
                    unk(gr_iSoln,i,j,1,b) = work(i,j,1,b,1)
                 enddo
              enddo
              do j = NGUARD+NYB-1, NGUARD+NYB
                 do i = NGUARD+1, NGUARD+NXB
                    work(i,j,1,b,1) = &
                         cx*(work(i-1,j,1,b,1) + work(i+1,j,1,b,1)) + &
                         cy*(work(i,j-1,1,b,1) + work(i,j+1,1,b,1)) - &
                         c*unk(gr_iSource,i,j,1,b)
                    unk(gr_iSoln,i,j,1,b) = work(i,j,1,b,1)
                 enddo
              enddo
              do j = NGUARD+1, NGUARD+NYB
                 do i = NGUARD+1, NGUARD+2
                    work(i,j,1,b,1) = &
                         cx*(work(i-1,j,1,b,1) + work(i+1,j,1,b,1)) + &
                         cy*(work(i,j-1,1,b,1) + work(i,j+1,1,b,1)) - &
                         c*unk(gr_iSource,i,j,1,b)
                    unk(gr_iSoln,i,j,1,b) = work(i,j,1,b,1)
                 enddo
                 do i = NGUARD+NXB-1, NGUARD+NXB
                    work(i,j,1,b,1) = &
                         cx*(work(i-1,j,1,b,1) + work(i+1,j,1,b,1)) + &
                         cy*(work(i,j-1,1,b,1) + work(i,j+1,1,b,1)) - &
                         c*unk(gr_iSource,i,j,1,b)
                    unk(gr_iSoln,i,j,1,b) = work(i,j,1,b,1)
                 enddo
              enddo
              
           else ! NDIM == 3
              
              c  = 0.5 / (deltas(JAXIS)**2*deltas(KAXIS)**2 + &
                   deltas(IAXIS)**2*deltas(KAXIS)**2 + deltas(IAXIS)**2*deltas(JAXIS)**2)
              cx = deltas(JAXIS)**2*deltas(KAXIS)**2 * c
              cy = deltas(IAXIS)**2*deltas(KAXIS)**2 * c
              cz = deltas(IAXIS)**2*deltas(JAXIS)**2 * c
              c  = deltas(IAXIS)**2*deltas(JAXIS)**2*deltas(KAXIS)**2 * c

              

              if(relaxflag) then ! True, Jacobi

              do k = NGUARD+1, NGUARD+2
                 do j = NGUARD+1, NGUARD+NYB
                    do i = NGUARD+1, NGUARD+NXB
                       work(i,j,k,b,1) = &
                            cx*(work(i-1,j,k,b,1) + work(i+1,j,k,b,1)) + &
                            cy*(work(i,j-1,k,b,1) + work(i,j+1,k,b,1)) + &
                            cz*(work(i,j,k-1,b,1) + work(i,j,k+1,b,1)) - &
                            c*unk(gr_iSource,i,j,k,b)
                       unk(gr_iSoln,i,j,k,b) = work(i,j,k,b,1)
                    enddo
                 enddo
              enddo
              do k = NGUARD+NZB-1, NGUARD+NZB
                 do j = NGUARD+1, NGUARD+NYB
                    do i = NGUARD+1, NGUARD+NXB
                       work(i,j,k,b,1) = &
                            cx*(work(i-1,j,k,b,1) + work(i+1,j,k,b,1)) + &
                            cy*(work(i,j-1,k,b,1) + work(i,j+1,k,b,1)) + &
                            cz*(work(i,j,k-1,b,1) + work(i,j,k+1,b,1)) - &
                            c*unk(gr_iSource,i,j,k,b)
                       unk(gr_iSoln,i,j,k,b) = work(i,j,k,b,1)
                    enddo
                 enddo
              enddo

              else ! Red Black Gauss Seidel

              ksweep = 0
              do redblackpass = 1, 2
                 jsweep = ksweep
                 do k = NGUARD+1-offou, NGUARD+2+offin
                    isweep = jsweep                    
                    do j = NGUARD+1, NGUARD+NYB
                       do i = NGUARD+1+isweep, NGUARD+NXB, 2
                       work(i,j,k,b,1) = &
                            cx*(work(i-1,j,k,b,1) + work(i+1,j,k,b,1)) + &
                            cy*(work(i,j-1,k,b,1) + work(i,j+1,k,b,1)) + &
                            cz*(work(i,j,k-1,b,1) + work(i,j,k+1,b,1)) - &
                            c*unk(gr_iSource,i,j,k,b)
                       unk(gr_iSoln,i,j,k,b) = work(i,j,k,b,1)
                       enddo
                       isweep = 1 - isweep
                    enddo
                    jsweep = 1 - jsweep 
                 enddo
                 ksweep = 1 - ksweep        
              enddo

              ksweep = 0
              do redblackpass = 1, 2
                 jsweep = ksweep
                 do k = NGUARD+NZB-1-offin, NGUARD+NZB+offou
                    isweep = jsweep
                    do j = NGUARD+1, NGUARD+NYB
                       do i = NGUARD+1+isweep, NGUARD+NXB, 2
                          work(i,j,k,b,1) = &
                            cx*(work(i-1,j,k,b,1) + work(i+1,j,k,b,1)) + &
                            cy*(work(i,j-1,k,b,1) + work(i,j+1,k,b,1)) + &
                            cz*(work(i,j,k-1,b,1) + work(i,j,k+1,b,1)) - &
                            c*unk(gr_iSource,i,j,k,b)
                          unk(gr_iSoln,i,j,k,b) = work(i,j,k,b,1)
                       enddo
                       isweep = 1 - isweep
                    enddo
                    jsweep = 1 - jsweep
                 enddo
                 ksweep = 1 - ksweep
              enddo


              endif

           

              if (relaxflag) then ! True, Jacobi

              do k = NGUARD+1, NGUARD+NZB
                 do j = NGUARD+1, NGUARD+2
                    do i = NGUARD+1, NGUARD+NXB
                       work(i,j,k,b,1) = &
                            cx*(work(i-1,j,k,b,1) + work(i+1,j,k,b,1)) + &
                            cy*(work(i,j-1,k,b,1) + work(i,j+1,k,b,1)) + &
                            cz*(work(i,j,k-1,b,1) + work(i,j,k+1,b,1)) - &
                            c*unk(gr_iSource,i,j,k,b)
                       unk(gr_iSoln,i,j,k,b) = work(i,j,k,b,1)
                    enddo
                 enddo
                 do j = NGUARD+NYB-1, NGUARD+NYB
                    do i = NGUARD+1, NGUARD+NXB
                       work(i,j,k,b,1) = &
                            cx*(work(i-1,j,k,b,1) + work(i+1,j,k,b,1)) + &
                            cy*(work(i,j-1,k,b,1) + work(i,j+1,k,b,1)) + &
                            cz*(work(i,j,k-1,b,1) + work(i,j,k+1,b,1)) - &
                            c*unk(gr_iSource,i,j,k,b)
                       unk(gr_iSoln,i,j,k,b) = work(i,j,k,b,1)
                    enddo
                 enddo
                 do j = NGUARD+1, NGUARD+NYB
                    do i = NGUARD+1, NGUARD+2
                       work(i,j,k,b,1) = &
                            cx*(work(i-1,j,k,b,1) + work(i+1,j,k,b,1)) + &
                            cy*(work(i,j-1,k,b,1) + work(i,j+1,k,b,1)) + &
                            cz*(work(i,j,k-1,b,1) + work(i,j,k+1,b,1)) - &
                            c*unk(gr_iSource,i,j,k,b)
                       unk(gr_iSoln,i,j,k,b) = work(i,j,k,b,1)
                    enddo
                    do i = NGUARD+NXB-1, NGUARD+NXB
                       work(i,j,k,b,1) = &
                            cx*(work(i-1,j,k,b,1) + work(i+1,j,k,b,1)) + &
                            cy*(work(i,j-1,k,b,1) + work(i,j+1,k,b,1)) + &
                            cz*(work(i,j,k-1,b,1) + work(i,j,k+1,b,1)) - &
                            c*unk(gr_iSource,i,j,k,b)
                       unk(gr_iSoln,i,j,k,b) = work(i,j,k,b,1)
                    enddo
                 enddo
              enddo
              
              else ! Red Black Gauss Seidel


              jsweep = 0
              do redblackpass = 1, 2
                 ksweep = jsweep
                 do j = NGUARD+1-offou, NGUARD+2+offin
                    isweep = ksweep
                    do k = NGUARD+1, NGUARD+NZB
                       do i = NGUARD+1+isweep, NGUARD+NXB, 2
                          work(i,j,k,b,1) = &
                            cx*(work(i-1,j,k,b,1) + work(i+1,j,k,b,1)) + &
                            cy*(work(i,j-1,k,b,1) + work(i,j+1,k,b,1)) + &
                            cz*(work(i,j,k-1,b,1) + work(i,j,k+1,b,1)) - &
                            c*unk(gr_iSource,i,j,k,b)
                          unk(gr_iSoln,i,j,k,b) = work(i,j,k,b,1)
                       enddo
                       isweep = 1 - isweep
                    enddo
                    ksweep = 1 - ksweep
                 enddo
                 jsweep = 1 - jsweep
              enddo


              jsweep = 0
              do redblackpass = 1, 2
                 ksweep = jsweep
                 do j = NGUARD+NYB-1-offin, NGUARD+NYB+offou
                    isweep = ksweep
                    do k = NGUARD+1, NGUARD+NZB
                       do i = NGUARD+1+isweep, NGUARD+NXB, 2
                          work(i,j,k,b,1) = &
                            cx*(work(i-1,j,k,b,1) + work(i+1,j,k,b,1)) + &
                            cy*(work(i,j-1,k,b,1) + work(i,j+1,k,b,1)) + &
                            cz*(work(i,j,k-1,b,1) + work(i,j,k+1,b,1)) - &
                            c*unk(gr_iSource,i,j,k,b)
                          unk(gr_iSoln,i,j,k,b) = work(i,j,k,b,1)
                       enddo
                       isweep = 1 - isweep
                    enddo
                    ksweep = 1 - ksweep
                 enddo
                 jsweep = 1 - jsweep
              enddo




              isweep = 0
              do redblackpass = 1, 2
                 ksweep = isweep
                 do i = NGUARD+1-offou, NGUARD+2+offin
                    jsweep = ksweep
                    do k = NGUARD+1, NGUARD+NZB
                       do j = NGUARD+1+jsweep, NGUARD+NYB, 2
                          work(i,j,k,b,1) = &
                            cx*(work(i-1,j,k,b,1) + work(i+1,j,k,b,1)) + &
                            cy*(work(i,j-1,k,b,1) + work(i,j+1,k,b,1)) + &
                            cz*(work(i,j,k-1,b,1) + work(i,j,k+1,b,1)) - &
                            c*unk(gr_iSource,i,j,k,b)
                          unk(gr_iSoln,i,j,k,b) = work(i,j,k,b,1)
                       enddo
                       jsweep = 1 - jsweep
                    enddo
                    ksweep = 1 - ksweep
                 enddo
                 isweep = 1 - isweep
              enddo


              isweep = 0
              do redblackpass = 1, 2
                 ksweep = isweep
                 do i = NGUARD+NXB-1-offin, NGUARD+NXB+offou
                    jsweep = ksweep
                    do k = NGUARD+1, NGUARD+NZB
                       do j = NGUARD+1+jsweep, NGUARD+NYB, 2
                          work(i,j,k,b,1) = &
                            cx*(work(i-1,j,k,b,1) + work(i+1,j,k,b,1)) + &
                            cy*(work(i,j-1,k,b,1) + work(i,j+1,k,b,1)) + &
                            cz*(work(i,j,k-1,b,1) + work(i,j,k+1,b,1)) - &
                            c*unk(gr_iSource,i,j,k,b)
                          unk(gr_iSoln,i,j,k,b) = work(i,j,k,b,1)
                       enddo
                       jsweep = 1 - jsweep
                    enddo
                    ksweep = 1 - ksweep
                 enddo
                 isweep = 1 - isweep
              enddo

              endif


           endif         

        endif
     enddo
  enddo


!!$
!!$!  call timer_stop("relaxation")
!!$  !TEST - BOUNDARY RELAXATION
!!$!  call timer_stop("gr_hgSolveLevel")
!!$  
!!$  !TEST - NORMALIZATION
!!$  !low face highest dimension, should really check all directions...
!!$  if ( ALL(gr_hgBndTypes(1:2*NDIM) == MG_BND_PERIODIC .or. gr_hgBndTypes(1:2*NDIM) == MG_BND_NEUMANN) ) then
!!$     
!!$     lvol = 0.
!!$     lsum = 0.
!!$     nbinv = 1. / real(NXB)
!!$     if (NDIM >= 2) nbinv = nbinv / real(NYB)
!!$     if (NDIM == 3) nbinv = nbinv / real(NZB)
!!$     
!!$     do b = 1, lnblocks
!!$        if ((lrefine(b) == level) .or. &
!!$             ((nodetype(b) == 1) .and. (lrefine(b) < level))) then
!!$           block_size(1:MDIM) = bsize(1:MDIM,b)
!!$           bvol = block_size(1)
!!$           if (NDIM >= 2) bvol = bvol * block_size(2)
!!$           if (NDIM == 3) bvol = bvol * block_size(3)
!!$           cvol = bvol * nbinv
!!$           lvol = lvol + bvol
!!$           bsum = 0.
!!$           do k = hg_kli, hg_kui
!!$              do j = hg_jli, hg_jui
!!$                 do i = hg_ili, hg_iui
!!$                    bsum = bsum + work(i,j,k,b,1)
!!$                 enddo
!!$              enddo
!!$           enddo
!!$           lsum = lsum + bsum * cvol
!!$        endif
!!$     enddo
!!$     
!!$     call mpi_allreduce ( lsum, sum, 1, FLASH_REAL, &
!!$          MPI_SUM, gr_meshComm, ierr )
!!$     call mpi_allreduce ( lvol, vol, 1, FLASH_REAL, &
!!$          MPI_SUM, gr_meshComm, ierr )
!!$     
!!$     avg = sum / vol
!!$     
!!$     do b = 1, lnblocks
!!$        if (lrefine(b) == level) then
!!$           work(:,:,:,b,1) = work(:,:,:,b,1) - avg
!!$           unk(gr_iSoln,:,:,:,b) = unk(gr_iSoln,:,:,:,b) - avg
!!$        endif
!!$     enddo
!!$     
!!$  endif
!!$  !TEST - NORMALIZATION
  
  call Timers_stop(level_timer_label)

  return
end subroutine gr_hgSolveLevel
