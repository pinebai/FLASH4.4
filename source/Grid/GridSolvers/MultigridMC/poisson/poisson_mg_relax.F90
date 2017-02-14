!!****if* source/Grid/GridSolvers/MultigridMC/poisson/poisson_mg_relax
!!
!! NAME
!!
!!  poisson_mg_relax
!!
!! SYNOPSIS
!!
!!  call poisson_mg_relax(:: level,
!!                         :: irhs,
!!                         :: ilhs,
!!                         :: nsmooth)
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
!!   nsmooth : 
!!
!! AUTOGENROBODOC
!!
!!
!!***

!*******************************************************************************

!  Routine:     mg_relax()

!  Description: Perform a number of smoothing iterations for multigrid.
!               This version implements Gauss-Seidel with red-black ordering 
!               using the Poisson equation.

!  Parameters:  level       Level to relax on.
!               irhs        Right-hand side (source) of equation.
!               ilhs        Left-hand side (solution) of equation.  Receives
!                           the smoothed result.
!               nsmooth     Number of iterations to apply.


subroutine poisson_mg_relax (level, irhs, ilhs, nsmooth)

  !===============================================================================


  use gr_mgData, ONLY:  ili, iui, jli, jui, kli, kui, gr_mgDiffOpDiscretize

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  use Grid_interface, ONLY : Grid_getLocalNumBlks, Grid_getBlkPtr, &
    Grid_releaseBlkPtr, Grid_getBlkBoundBox, Grid_getBlkIndexLimits, &
    Grid_getListOfBlocks, Grid_getBlkPhysicalSize

  use Grid_data, ONLY : gr_meshMe

  use workspace, ONLY: work
  use tree, only : nodetype,lrefine
  use paramesh_dimensions, only : nguard_work

  use Timers_interface, ONLY : Timers_start, Timers_stop

  implicit none
#include "Multigrid.h"
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"

  integer       :: level, irhs, ilhs, nsmooth, idebug
  integer       :: firstpass = 1

  integer       :: iter, i, j, k, lb, ierr, Neff, lnblocks
  integer       :: redblackpass, isweep, jsweep, ksweep
  real          :: lerrnorm(2), errnorm(2), nxbinv, nybinv, nzbinv
  real          :: delx, dely, delz, error, prefac, critfac, pi
  real          :: idelx,idely,idelz

  real :: delydelz, delxdelz, delxdely, delxdelydelz

  logical       :: done

  character(len=256) :: str_buffer

  integer :: iterating_to_convergence_limit
  integer, parameter :: MAXDIMS = 3
  real, dimension(MAXDIMS) :: size

  logical, save :: first_call = .true.
  real, save    :: mgrid_smooth_tol
  integer, save :: Nblockx, Nblocky, Nblockz, myPE, masterPE

  real, pointer, dimension(:,:,:,:), save :: unk

  integer, parameter :: nxb = NXB
  integer, parameter :: nyb = NYB
  integer, parameter :: nzb = NZB
  integer, parameter :: ndim = NDIM


  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC

  real, dimension(2,MDIM) :: boundBox

  integer blockcount,ii,jj,blockID

  real bsize(MDIM),coord(MDIM)

  integer blockList(MAXBLOCKS)

  !               On the first call, check to make sure that the requested
  !               convergence criterion is reasonable.  The smoother relaxes
  !               long wavelengths (relative to the grid) most slowly.  The
  !               termination criterion stops the iteration when the iterate
  !               has stopped changing to within some tolerance.  This
  !               tolerance should not be larger than the convergence rate
  !               for the longest-wavelength mode.  The convergence rate
  !               used here is given by Briggs et al. (_A Multigrid Tutorial_)
  !               for Dirichlet boundary conditions.

  !               The way we do this comparison at the end may not be quite
  !               right... should probably compute the ratio of the new
  !               residual norm to the old residual norm.

!  call timer_start("relax")

  iterating_to_convergence_limit = nsmooth
  error    = 1.

  if (first_call) then
     call RuntimeParameters_get('mgrid_smooth_tol',    mgrid_smooth_tol)
     call RuntimeParameters_get('Nblockx',  Nblockx)
     call RuntimeParameters_get('Nblocky',  Nblocky)
     call RuntimeParameters_get('Nblockz',  Nblockz)
     myPE     = gr_meshMe 
     masterPE = MASTER_PE 
  end if

  if (first_call) then
     if (MyPE == MasterPE) then
        pi = PI
        i = 1      ! Will only try to converge via relaxation on level 1
        Neff = max( nxb*Nblockx*2**(i-1), & 
             &                  nyb*Nblocky*2**(i-1), & 
             &                  nzb*Nblockz*2**(i-1) )
        critfac = 2.*0.66667*sin(pi/(2.*Neff))**2
        if (mgrid_smooth_tol >= critfac) then
           write (str_buffer,*) 'mg_relax: termination', & 
                & ' tolerance of ', mgrid_smooth_tol, ' is larger', &
                & ' than the slowest', & 
                & ' convergence rate on level ', i, ' (', critfac, ').', &
                & ' Long-wavelength modes', & 
                & ' may not be able to converge on this level.'
           write(*,*) str_buffer
           !call stamp_logfile(str_buffer, 'warning')
        endif
     endif
     first_call = .false.
  endif

  !               Iteration loop.

  call Grid_getLocalNumBlks(lnblocks)

  iter    = 0
  nxbinv  = 1./nxb**2
  nybinv  = 1./nyb**2
  nzbinv  = 1./nzb**2
  done    = .false.

  
  !               Update boundary zones of the LHS array for this iteration.
  !               Currently the mesh package doesn't supply a general boundary
  !               update routine, so we must copy the LHS into the "work" array,
  !               then update its boundaries.  This is handled by mg_bndry().
  !               This relaxer assumes the LHS variable is defined on all blocks
  !               (not just leaf nodes), so we call mg_bndry with leaf_only == 0.
     
  ! copy unk into work and fill work's guardcells
! ***
  call Timers_start("gr_mgBndry")
  call gr_mgBndry (level, ilhs, 2, 0, MG_COPY_UNK_TO_WORK, MG_BEGIN_SERIES) 
  call Timers_stop("gr_mgBndry")

  do while ( .not. done )

     if (iter > iterating_to_convergence_limit) then 
        ! if we're iterating to convergence, 
        ! copy the old iterate for convergence testing.  
        do lb = 1, lnblocks
           if (lrefine(lb) == level) then

              ! Point to blocks center vars:
              call Grid_getBlkPtr(lb,unk,CENTER)

              do k = kli, kui
                 do j = jli, jui
                    do i = ili, iui
                       unk(ilhs,i,j,k) = work(i,j,k,lb,1)
                    enddo
                 enddo
              enddo

              ! Point to blocks center vars:
              call Grid_releaseBlkPtr(lb,unk,CENTER)


           endif
        enddo
     end if

     !               Perform the Gauss-Seidel iteration step. 
 
!     call timer_start("relax calc")

     if (ndim == 1) then
        do lb = 1, lnblocks
           if (lrefine(lb) == level) then
              call Grid_getBlkPhysicalSize(lb,size)
              delx = size(1)**2 * nxbinv

              ! Point to blocks center vars:
              call Grid_getBlkPtr(lb,unk,CENTER)

              isweep = 0
              do redblackpass = 1, 2
                 do i = ili+isweep, iui, 2  
                    work(i,jli:jui,kli:kui,lb,1) = & 
                         &  0.5 * ( work(i-1,jli:jui,kli:kui,lb,1) + & 
                         &          work(i+1,jli:jui,kli:kui,lb,1) - & 
                         &          delx*unk(irhs,i,jli:jui,kli:kui) )
                 end do
                 isweep = 1 - isweep 
              end do

              ! Point to blocks center vars:
              call Grid_releaseBlkPtr(lb,unk,CENTER)

           endif
        enddo

     elseif (ndim == 2) then

        do lb = 1, lnblocks
           if (lrefine(lb) == level) then
              call Grid_getBlkPhysicalSize(lb,size)

              delx = size(1)**2 * nxbinv
              dely = size(2)**2 * nybinv
              prefac = 0.5 / (delx+dely)

              ! Point to blocks center vars:
              call Grid_getBlkPtr(lb,unk,CENTER)
              jsweep = 0
              do redblackpass = 1, 2
                 isweep = jsweep
                 do j = jli, jui
                    do i = ili + isweep, iui, 2
                       work(i,j,kli:kui,lb,1) = & 
                            &              prefac * & 
                            &              ( dely*(work(i-1,j,kli:kui,lb,1) + & 
                            &                      work(i+1,j,kli:kui,lb,1)) + & 
                            &                delx*(work(i,j-1,kli:kui,lb,1) + & 
                            &                      work(i,j+1,kli:kui,lb,1)) - & 
                            &                delx*dely*unk(irhs,i,j,kli:kui) )
                    enddo
                    isweep = 1 - isweep
                 enddo
                 jsweep = 1 - jsweep
              enddo

              ! Point to blocks center vars:
              call Grid_releaseBlkPtr(lb,unk,CENTER)

           endif
        enddo




     else

        select case(gr_mgDiffOpDiscretize)

        case(2)    ! 2nd order Central difference Solution
        do lb = 1, lnblocks
           if (lrefine(lb) == level) then
              call Grid_getBlkPhysicalSize(lb,size)

              delx = size(1)**2 * nxbinv
              dely = size(2)**2 * nybinv
              delz = size(3)**2 * nzbinv

              delydelz = dely*delz
              delxdelz = delx*delz
              delxdely = delx*dely

              prefac = 0.5 / (delydelz+delxdelz+delxdely)

              delydelz = prefac*delydelz
              delxdelz = prefac*delxdelz
              delxdely = prefac*delxdely
              delxdelydelz = prefac*delx*dely*delz


              ! Point to blocks center vars:
              call Grid_getBlkPtr(lb,unk,CENTER)
              ksweep = 0
              do redblackpass = 1, 2
                 jsweep = ksweep
                 do k = kli, kui
                    isweep = jsweep
                    do j = jli, jui
                       do i = ili + isweep, iui, 2
                          work(i,j,k,lb,1) = & !&  prefac * & 
                               &  ( delydelz*(work(i-1,j,k,lb,1) + work(i+1,j,k,lb,1)) + & 
                               &    delxdelz*(work(i,j-1,k,lb,1) + work(i,j+1,k,lb,1)) + & 
                               &    delxdely*(work(i,j,k-1,lb,1) + work(i,j,k+1,lb,1)) - & 
                               &    delxdelydelz*unk(irhs,i,j,k) )
                       enddo
                       isweep = 1 - isweep
                    enddo
                    jsweep = 1 - jsweep
                 enddo
                 ksweep = 1 - ksweep
              enddo
              ! Point to blocks center vars:
              call Grid_releaseBlkPtr(lb,unk,CENTER)
           endif
        enddo

        case(4)    ! 4th order Central difference Solution
        do lb = 1, lnblocks
           if (lrefine(lb) == level) then
              call Grid_getBlkPhysicalSize(lb,size)

              idelx = (size(1)**2 * nxbinv)**(-1.) !1/dx**2
              idely = (size(2)**2 * nybinv)**(-1.) !1/dy**2
              idelz = (size(3)**2 * nzbinv)**(-1.) !1/dz**2

              prefac = 576./1460.*(idelx+idely+idelz)**(-1.)    

              ! Point to blocks center vars:
              call Grid_getBlkPtr(lb,unk,CENTER)
              ksweep = 0
              do redblackpass = 1, 2
                 jsweep = ksweep
                 do k = kli, kui
                    isweep = jsweep
                    do j = jli, jui
                       do i = ili + isweep, iui, 2

                          work(i,j,k,lb,1) = & 
                               &  prefac * & 
                               &  (0.001736111111111 * &
                               &  (idelx*(     (work(i-3,j,k,lb,1) + work(i+3,j,k,lb,1)) - &
                               &           54.*(work(i-2,j,k,lb,1) + work(i+2,j,k,lb,1)) + &
                               &          783.*(work(i-1,j,k,lb,1) + work(i+1,j,k,lb,1)))+ &
                               &   idely*(     (work(i,j-3,k,lb,1) + work(i,j+3,k,lb,1)) - &
                               &           54.*(work(i,j-2,k,lb,1) + work(i,j+2,k,lb,1)) + &
                               &          783.*(work(i,j-1,k,lb,1) + work(i,j+1,k,lb,1)))+ &
                               &   idelz*(     (work(i,j,k-3,lb,1) + work(i,j,k+3,lb,1)) - &
                               &           54.*(work(i,j,k-2,lb,1) + work(i,j,k+2,lb,1)) + &
                               &          783.*(work(i,j,k-1,lb,1) + work(i,j,k+1,lb,1)))) &
                               &   - unk(irhs,i,j,k) )

                       enddo
                       isweep = 1 - isweep
                    enddo
                    jsweep = 1 - jsweep
                 enddo
                 ksweep = 1 - ksweep
              enddo
              ! Point to blocks center vars:
              call Grid_releaseBlkPtr(lb,unk,CENTER)
           endif
        enddo

        end select


     endif

     !               Now compute the error norm (change in LHS from one iteration
     !               to the next).  This will be used to limit the number of
     !               iterations (see above).  Here we use the maximum norm rather
     !               than the L2 norm computed by mg_norm.

     !               However, don't start checking the convergence criterion
     !               until we've done ten iterations.  This saves work when doing
     !               the regular smoothing on fine levels, which isn't expected
     !               to produce convergence by itself.


     iter = iter + 1

     if (iter > iterating_to_convergence_limit) then

        lerrnorm(:) = 0.
      
        do lb = 1, lnblocks
           if (lrefine(lb) == level) then

              ! Point to blocks center vars:
              call Grid_getBlkPtr(lb,unk,CENTER)

              do k = kli, kui
                 do j = jli, jui
                    do i = ili, iui
                       lerrnorm(1) = & 
                            max( lerrnorm(1), & 
                            abs(unk(ilhs,i,j,k)-work(i,j,k,lb,1)) )
                       lerrnorm(2) = & 
                            max( lerrnorm(2), abs(unk(ilhs,i,j,k)) )
                    enddo
                 enddo
              enddo

              ! Point to blocks center vars:
              call Grid_releaseBlkPtr(lb,unk,CENTER)

           endif
        enddo

        call mpi_allreduce (lerrnorm, errnorm, 2, MPI_DOUBLE_PRECISION, & 
             MPI_MAX, MPI_COMM_WORLD, ierr)

        error = mgrid_smooth_tol*errnorm(2) - errnorm(1)

     endif

     !               End of loop.  Iterate until nsmooth iterations have been
     !               performed, or until the convergence criterion has been met.

     done = (iter == nsmooth) .or. ((iter > iterating_to_convergence_limit) .and. (error > 0.))
!     call timer_stop("relax calc")

     ! fill work's guardcells
     if (.not. done) then
     call Timers_start("gr_mgBndry")
     call gr_mgBndry (level, ilhs, nguard_work, 0, MG_EXCHANGE_WORK, MG_CONTINUE_SERIES) 
     call Timers_stop("gr_mgBndry")
     endif
  enddo

  do lb = 1, lnblocks
     if (lrefine(lb) == level) then

        ! Point to blocks center vars:
        call Grid_getBlkPtr(lb,unk,CENTER)

        do k = kli, kui
           do j = jli, jui
              do i = ili, iui
                 unk(ilhs,i,j,k) = work(i,j,k,lb,1)
              enddo
           enddo
        enddo

        ! Point to blocks center vars:
        call Grid_releaseBlkPtr(lb,unk,CENTER)


     endif
  enddo

  !===============================================================================

!  call timer_stop("relax")
  return
end subroutine poisson_mg_relax
