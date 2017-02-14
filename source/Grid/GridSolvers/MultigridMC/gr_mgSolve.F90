!!****if* source/Grid/GridSolvers/MultigridMC/gr_mgSolve
!!
!! NAME
!!
!!  gr_mgSolve
!!
!! SYNOPSIS
!!
!!  call gr_mgSolve(integer(in) :: isrc_dens,
!!                  integer(in) :: img_soln,
!!                  real(in) :: poisfact,
!!                  integer(in) :: img_src,
!!                  integer(in) :: img_res,
!!                  integer(in) :: img_corr,
!!                  integer(in) :: img_temp,
!!                  integer(in) :: img_temp2,
!!                  integer(in) :: bc_types,
!!                   :: mg_solve,
!!                   :: mg_residual,
!!                   :: mg_residual_leafs,
!!                   :: mg_relax)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   isrc_dens : 
!!
!!   img_soln : 
!!
!!   poisfact : 
!!
!!   img_src : 
!!
!!   img_res : 
!!
!!   img_corr : 
!!
!!   img_temp : 
!!
!!   img_temp2 : 
!!
!!   bc_types : 
!!
!!   mg_solve : 
!!
!!   mg_residual : 
!!
!!   mg_residual_leafs : 
!!
!!   mg_relax : 
!!
!! AUTOGENROBODOC
!!
!!
!!***

!*******************************************************************************

!  Routine:     gr_mgSolve()

!  Description: Main multigrid Poisson solver.  The solution is preconditioned
!               on the first call; afterward, the existing solution is taken
!               as the initial guess.

!               gr_mgSolve multiplies unk(isrc_dens) by poisfact and puts the
!               result into unk(img_src).  If the img_bc arg is 0 (dirichlet/
!               given-value), take boundary values for potential from boundaries
!               of unk(1st arg).  Then solve.  On output, solution goes into
!               unk(1st arg); unk(2nd arg) is not modified (unless 2nd
!               arg=img_src).

!  Parameters:  isrc_dens       Index for source array.  This is taken to be
!                               the density field; the source array to be used
!                               as the right-hand side of the Poisson equation
!                               is computed from this.
!               img_soln        Index for solution array.  The solution is
!                               written directly into this variable.
!               poisfact        Constant Poisson factor.
!               img_src, img_res, img_corr, img_temp
!                               Indices for work arrays to hold, respectively,
!                               the source (right-hand side), residual,
!                               correction, and temporary values.
!               img_bc          Boundary condition to apply to all boundaries.
!                                 0 = Isolated boundaries using James' algorithm
!                                 1 = Periodic boundaries
!                                 2 = Dirichlet boundaries
!                                 3 = Neumann boundaries


subroutine gr_mgSolve (isrc_dens, img_soln, poisfact, img_src, img_res, &
                       img_corr, img_temp, img_temp2, bc_types,mg_solve,&
                       mg_residual, mg_residual_leafs, mg_relax)

!===============================================================================

#include "Flash.h"

  use gr_mgData, ONLY: ile, iue, jle, jue, kle, kue, &
     mg_bnd_cond, mesh_lrefmin, mesh_lrefmax, interp_work,interp_mask_work_mg, &
     interp_mask_work_save,nodetype_save, &
     newchild_save, gr_mgDiffOpDiscretize

  use Grid_interface,    ONLY : Grid_getLocalNumBlks, &
                                Grid_getListOfBlocks, &
                                Grid_getBlkPtr,       &
                                Grid_releaseBlkPtr

  use Grid_data, ONLY : gr_meshMe

  use gr_mgInterface, ONLY : gr_mgInitSlv,gr_mgInitSrc, &
    gr_mgNorm, mg_restore_nodetypes, gr_mgCycle,     &
    gr_mgBndry

  use gr_interface, ONLY : gr_findMean

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  use workspace, ONLY: work, interp_mask_work
  use paramesh_dimensions, only : nguard_work,nvar_work

  use tree, only : maxblocks_tr

  use Timers_interface, ONLY : Timers_start, Timers_stop

  use paramesh_mpi_interfaces, ONLY : amr_morton_process

implicit none
#include "Multigrid.h"
#include "constants.h"

integer,intent(in) :: isrc_dens, img_soln, img_src, img_res, img_corr, img_temp, &
                 img_temp2
integer,intent(in) :: bc_types(:)
real,intent(in)    :: poisfact

external mg_solve, mg_residual, mg_residual_leafs, mg_relax

integer :: blockCount
integer :: blockList(MAXBLOCKS)

logical       :: done
integer       :: i, j, k, ii, jj, kk, lb, lnblocks, level
real          :: res_norm_ratio, res_norm_change, norm_old, norm_new, norm_rhs
real, save    :: mgrid_max_residual_norm, mgrid_max_iter_change
integer, save :: mgrid_max_vcycles, MyPE, MasterPE
logical, save :: mgrid_print_norm
logical, save :: first_call = .true.
logical, save :: first_step = .true.

integer :: blockID,img_bc

real, pointer, dimension(:,:,:,:) :: unk

real :: mean_soln

!===============================================================================

! Initialize
call Timers_start("Initialize_MultigridMC")

call Grid_getListOfBlocks(LEAF,blockList,blockCount)
if (first_call) then

  call RuntimeParameters_get('mgrid_max_residual_norm',    mgrid_max_residual_norm)
  call RuntimeParameters_get('mgrid_max_iter_change',    mgrid_max_iter_change)
  call RuntimeParameters_get('mgrid_max_vcycles',    mgrid_max_vcycles)
  call RuntimeParameters_get('mgrid_print_norm',    mgrid_print_norm)

  MyPE = gr_meshMe
  MasterPE = MASTER_PE

  allocate(interp_mask_work_mg(nvar_work))

  call RuntimeParameters_get("gr_mgDiffOpDiscretize",gr_mgDiffOpDiscretize)
  select case(gr_mgDiffOpDiscretize)
  case(2) ! 2nd Order Central differences
    interp_work = 2
  case(4) ! 4th Order Central differences
    interp_work = 4
    write(*,*) 'IN gr_mgSolve START gr_mgDiffOpDiscretize = 4'
  end select

  interp_mask_work_mg = interp_work
  allocate(interp_mask_work_save(nvar_work))

  allocate(nodetype_save(maxblocks_tr))
  allocate(newchild_save(maxblocks_tr))

  first_call = .false.
endif

interp_mask_work_save = interp_mask_work;
interp_mask_work = interp_mask_work_mg;

call gr_mgInitSlv(bc_types)
call gr_mgInitSrc (isrc_dens, poisfact, img_src, img_soln)

! Iterate
done = .false.
i    = 0
call gr_mgNorm (0, img_src, norm_rhs, 1)

call mg_residual_leafs (img_src, img_soln, img_res, 1)

!!$do level = mesh_lrefmax, mesh_lrefmin, -1
!!$   call Timers_start("poisson_mg_Residual_initial")
!!$   call mg_residual (level, img_src, img_soln, img_res, 1, 1)
!!$   call Timers_stop("poisson_mg_Residual_initial")
!!$enddo
!!$call mg_restore_nodetypes (mesh_lrefmax)

call gr_mgNorm (0, img_res, norm_old, 1)

call Timers_stop("Initialize_MultigridMC")

! if the solution is not good enough, then cycle
if (norm_old > mgrid_max_residual_norm * norm_rhs) then

   do while (.not. done)
      call Timers_start("gr_mgCycle")
      call gr_mgCycle (mesh_lrefmax, img_soln, img_src, img_res, img_corr, &
           img_temp, img_temp2, mg_solve, mg_residual, mg_relax)
      call Timers_stop("gr_mgCycle")

!!$      ! --------------------------------------------------------------
!!$      ! Substract Mean value of solution:
!!$      if (mg_bnd_cond == 0) then ! Substract solution mean.
!!$      call gr_findMean(img_soln,2,.false.,mean_soln)
!!$      call Grid_getListOfBlocks(LEAF,blockList,blockCount)
!!$      do lb = 1, blockCount
!!$      blockID = blockList(lb)
!!$      ! Point to blocks center vars:
!!$      call Grid_getBlkPtr(blockID,unk,CENTER)
!!$      do kk = kle, kue
!!$         do jj = jle, jue
!!$            do ii = ile, iue
!!$               unk(img_soln,ii,jj,kk) = unk(img_soln,ii,jj,kk) - mean_soln
!!$            enddo
!!$         enddo
!!$      enddo
!!$      ! Release pointers:
!!$      call Grid_releaseBlkPtr(blockID,unk,CENTER)
!!$      enddo
!!$      endif
!!$      ! --------------------------------------------------------------


!!$      do level = mesh_lrefmax, 1, -1
!!$         call Timers_start("poisson_mg_Residual")
!!$         call mg_residual (level, img_src, img_soln, img_res, 1, 1)
!!$         call Timers_stop("poisson_mg_Residual")
!!$      enddo
!!$      call mg_restore_nodetypes(mesh_lrefmax)

      call mg_residual_leafs (img_src, img_soln, img_res, 1)
      call gr_mgNorm (0, img_res, norm_new, 1)

      res_norm_ratio  = norm_new / norm_rhs
      res_norm_change = (norm_new - norm_old) / norm_old

      if ((i > 0) .and. (res_norm_change > 0. .and. &
      norm_new > mgrid_max_residual_norm) .and. (MyPE == MasterPE)) &
           print *, "gr_mgSolve:  WARNING:  V-cycles not converging"

      i = i + 1

      if ((mgrid_print_norm) .and. (MyPE == MasterPE)) then
               if (norm_rhs > 0.e0) then
                  write(*,'(a,i4,3(a,es9.2))')                 &
                    'cycle ', i,                               &
                    ' : res_norm_ratio = ', res_norm_ratio,    &
                    ' res norm = ',         norm_new,          &
                    ' new to old ratio = ', norm_new/norm_old
               else
                  write(*,'(a,i4,a,es9.2)')                    &
                     'cycle ', i,                              &
                     ' : res_norm_change = ', res_norm_change
               end if
      endif

      if ((i == mgrid_max_vcycles) .or. &
           (norm_new <= mgrid_max_residual_norm * norm_rhs) .or. &
           (abs(res_norm_change) < mgrid_max_iter_change) .and.  &
           (i .gt. 0)) &
           done = .true.

      norm_old = norm_new

   enddo

   ! Make sure the boundary zones are set properly
   do level = 1, mesh_lrefmax
      call Timers_start("gr_mgBndry")
      call gr_mgBndry (level, img_soln, nguard_work, 1, MG_COPY_UNK_TO_WORK, MG_CONTINUE_SERIES)
      call Timers_stop("gr_mgBndry")
   enddo

   ! Copy the updated guard cell data back into the database.
   call Grid_getListOfBlocks(LEAF,blockList,blockCount)

   do lb = 1, blockCount
      blockID = blockList(lb)
      ! Point to blocks center vars:
      call Grid_getBlkPtr(blockID,unk,CENTER)
      do k = kle, kue
         do j = jle, jue
            do i = ile, iue
               unk(img_soln,i,j,k) = work(i,j,k,blockID,1)
            enddo
         enddo
      enddo
      ! Release pointers:
      call Grid_releaseBlkPtr(blockID,unk,CENTER)
   enddo

   ! Substract mean value from solution:
   if (mg_bnd_cond == 0) then ! Substract solution mean.
     call gr_findMean(img_soln,2,.false.,mean_soln)
     do lb = 1, blockCount
      blockID = blockList(lb)
      ! Point to blocks center vars:
      call Grid_getBlkPtr(blockID,unk,CENTER)
      do k=kle,kue
         do j=jle,jue
            do i=ile,iue
               unk(img_soln,i,j,k) = unk(img_soln,i,j,k) - mean_soln
            enddo
         enddo
      enddo
      ! Release pointers:
      call Grid_releaseBlkPtr(blockID,unk,CENTER)
     enddo
   endif

end if

interp_mask_work = interp_mask_work_save

! Restore communication pattern to diagonals gcell filling:
call amr_morton_process()

!===============================================================================

return
end
