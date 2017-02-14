!!****if* source/Grid/GridSolvers/BiPCGStab/gr_bipcgstab
!!
!! NAME
!!
!!  gr_bipcgstab
!!
!! SYNOPSIS
!!
!!  call gr_bipcgstab(integer(in) :: isrc,
!!                    integer(in) :: isoln,
!!                    real    :: poisFact,
!!                    integer(in) :: bcro,
!!                    integer(in) :: bcri,
!!                    integer(in) :: bcvi,
!!                    integer(in) :: bcpi,
!!                    integer(in) :: bcsi,
!!                    integer(in) :: bczi,
!!                    integer(in) :: bcyi,
!!                    integer(in) :: bcTypes,
!!                    real(in), dimension(2,6)  :: bcValues)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   isrc : 
!!
!!   isoln : 
!!
!!   poisFact : 
!!
!!   bcro : 
!!
!!   bcri : 
!!
!!   bcvi : 
!!
!!   bcpi : 
!!
!!   bcsi : 
!!
!!   bczi : 
!!
!!   bcyi : 
!!
!!   bcTypes : 
!!
!!   bcValues : 
!!
!! AUTOGENROBODOC
!!
!!
!!***


subroutine gr_bipcgstab (isrc,isoln,poisFact, bcro, bcri, bcvi, bcpi, bcsi, bczi, bcyi, &
                      bcTypes,bcValues)


!===============================================================================

#include "Flash.h"

  use Grid_interface,    ONLY : Grid_getLocalNumBlks, &
                                Grid_getListOfBlocks, &
                                Grid_getBlkPtr,       &
                                Grid_releaseBlkPtr

  use Grid_data, ONLY : gr_meshMe

  use gr_bicgInterface, ONLY : gr_bicgInitSlv,gr_bicgInitSrc, &
    gr_bicgNorm, gr_bicgDotprod, gr_bicgMultiAx, gr_bicgApplyPrecond

  use gr_interface, ONLY : gr_findMean

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  
  use workspace, ONLY: work, interp_mask_work
  use paramesh_dimensions, only : nguard_work,nvar_work

  use gr_bicgData, ONLY : ili, iui, jli, jui, kli, kui, &
                          ile, iue, jle, jue, kle, kue, &
                          interp_work_bipcgstab,        &
                          interp_mask_work_save,        &
                          precond_init, gcellflg, bicg_bnd_cond
  use Timers_interface, ONLY: Timers_start, Timers_stop


implicit none
#include "constants.h"

integer, intent(in)       :: isrc, isoln, bcro, bcri, bcvi, bcpi, bcsi, bczi, bcyi 
integer, intent(in)       :: bcTypes(6)
real, intent(in), dimension(2,6) :: bcValues
real          :: poisFact

integer :: blockCount
integer :: blockList(MAXBLOCKS)

logical       :: done
integer       :: i, j, k, ii, jj, kk, lb, lnblocks, level
real          :: res_norm_ratio, res_norm_change, norm_old, norm_new, norm_rhs
real          :: numer, denom
real, save    :: bipcgs_max_residual_norm
integer, save :: bipcgs_max_cycles
integer, save :: MyPE, MasterPE
logical, save :: bipcgs_print_norm
logical, save :: firstCall = .true.

real :: rho, rho1,beta,alpha,wi

integer :: blockID, iter

real, pointer, dimension(:,:,:,:) :: unk

real :: mean_soln

integer TA(2),count_rate
real  :: ET 

!===============================================================================

ET = 0.

! Initialize:
call Grid_getListOfBlocks(LEAF,blockList,blockCount)
if (firstCall) then
 
  call RuntimeParameters_get('bipcgs_max_residual_norm',    bipcgs_max_residual_norm)
  call RuntimeParameters_get('bipcgs_max_cycles',    bipcgs_max_cycles)
  call RuntimeParameters_get('bipcgs_print_norm',    bipcgs_print_norm)

  MyPE = gr_meshMe 
  MasterPE = MASTER_PE

  allocate(interp_mask_work_save(nvar_work))

endif

interp_mask_work_save = interp_mask_work;
interp_mask_work = interp_work_bipcgstab;


! Get list of processors leaf blocks:
call Grid_getListOfBlocks(LEAF,blockList,blockCount)

!Initialize solver vars and source:
call gr_bicgInitSlv(bcTypes)
call gr_bicgInitSrc (isrc, isoln, poisFact)

if (firstCall) then
   if (MyPE .eq. MASTER_PE) write(*,*) 'gr_bipcgstab first call ...'
endif

! Get norm of the rhs:
call gr_bicgNorm (isrc, norm_rhs)
if (norm_rhs .eq. 0.0) norm_rhs = 1.

! Call Matmult for A*xo
call Timers_start("gr_bicgMultiAx")
gcellflg = .true.
call gr_bicgMultiAx(isoln, bcro, gcellflg)
call Timers_stop("gr_bicgMultiAx")

! ro = f - A*xo, vo=po=0
do lb = 1, blockCount
   blockID = blockList(lb)
   ! Point to blocks center vars:
   call Grid_getBlkPtr(blockID,unk,CENTER)
   do k=kli,kui
      do j=jli,jui
         do i=ili,iui
            unk(bcro,i,j,k) = unk(isrc,i,j,k) - unk(bcro,i,j,k)
                                                ! Contains A*xo
            !unk(bcvi,i,j,k) = 0.0
            !unk(bcpi,i,j,k) = 0.0
            ! Define bcri-1 in bcri
            unk(bcri,i,j,k) = unk(bcro,i,j,k)
         enddo                                  
      enddo
   enddo
   ! Release pointers:
   call Grid_releaseBlkPtr(blockID,unk,CENTER)
enddo

! Get norm of ro:
call gr_bicgNorm (bcri, norm_old)
norm_old = norm_old / norm_rhs

! ro hut is defined in bcro: rhoh=ro.

! Parameters:
rho1   = 1.
alpha  = 1.
wi     = 1.

! Initialize preconditioner set to true, to be used in first call to 
! gr_bicgApplyPrecond:
precond_init = .true.
! Start Main Loop:
do iter = 1, bipcgs_max_cycles

   ! rho = (roh,bcri-1)
   call Timers_start("gr_bicgDotProd")
   call gr_bicgDotprod(bcro,bcri, rho)
   call Timers_stop("gr_bicgDotProd")

   if (rho .eq. 0.0) then
      ! Write residual and exit:
      if ((bipcgs_print_norm) .and. (MyPE == MasterPE)) then
         write(*,'(a,i4,2(a,es9.2))')                    &
              'begin cyc', iter,                         &
              ' rho == 0.0, exit.'
      endif
      goto 919
   endif


   ! Pi: pi = ri-1 + beta*(pi-1 - wi-1*vi-1)
   if (iter .eq. 1) then
      do lb = 1, blockCount
         blockID = blockList(lb)
         ! Point to blocks center vars:
         call Grid_getBlkPtr(blockID,unk,CENTER)
         do k=kli,kui
            do j=jli,jui
               do i=ili,iui
                  unk(bcpi,i,j,k)=unk(bcri,i,j,k) 
               enddo
            enddo
         enddo
         ! Release pointers:
         call Grid_releaseBlkPtr(blockID,unk,CENTER)
      enddo
   else

      ! beta = (rhoi/rhoi-1)*(alpha/wi-1)
      beta = (rho/rho1)*(alpha/wi)
      do lb = 1, blockCount
         blockID = blockList(lb)
         ! Point to blocks center vars:
         call Grid_getBlkPtr(blockID,unk,CENTER)
         do k=kli,kui
            do j=jli,jui
               do i=ili,iui
                  unk(bcpi,i,j,k)=unk(bcri,i,j,k) + beta*(unk(bcpi,i,j,k) - wi*unk(bcvi,i,j,k)) 
               enddo
            enddo
         enddo
         ! Release pointers:
         call Grid_releaseBlkPtr(blockID,unk,CENTER)
      enddo
   endif


   if (gr_meshMe .eq. 0) CALL SYSTEM_CLOCK(TA(1),count_rate)

   ! Solve M*y =pi
   call gr_bicgApplyPrecond(bcyi, bcpi, bcTypes, bcValues)
   precond_init = .false. ! After first solvePrecond do not initialize preconditioner.

   if (gr_meshMe .eq. 0) then
      CALL SYSTEM_CLOCK(TA(2),count_rate)
      ET=ET+REAL(TA(2)-TA(1))/count_rate
   endif


   ! vi = A*y
   !gcellflg = .true. or .false. done on gr_bicgApplyPrecond
   call Timers_start("gr_bicgMultiAx")
   call gr_bicgMultiAx (bcyi, bcvi, gcellflg)
   call Timers_stop("gr_bicgMultiAx")   

   ! alpha = rho / (roh,vi)
   call Timers_start("gr_bicgDotProd")
   call gr_bicgDotprod(bcro,bcvi, denom)
   call Timers_stop("gr_bicgDotProd")
   alpha = rho / denom

   ! si = ri - alpha*vi
   do lb = 1, blockCount
      blockID = blockList(lb)
      ! Point to blocks center vars:
      call Grid_getBlkPtr(blockID,unk,CENTER)
      do k=kli,kui
         do j=jli,jui
            do i=ili,iui
               unk(bcsi,i,j,k)=unk(bcri,i,j,k) - alpha*unk(bcvi,i,j,k) 
            enddo
         enddo
      enddo
      ! Release pointers:
      call Grid_releaseBlkPtr(blockID,unk,CENTER)
   enddo  

   ! Early Convergence check on si
   call gr_bicgNorm (bcsi, norm_new)
   norm_new = norm_new / norm_rhs
   if ( norm_new .le. bipcgs_max_residual_norm) then  ! Check early convergence.
      ! Compute xi = xi + alpha*yi
      do lb = 1, blockCount
         blockID = blockList(lb)
         ! Point to blocks center vars:
         call Grid_getBlkPtr(blockID,unk,CENTER)
         do k=kli,kui
            do j=jli,jui
               do i=ili,iui
                  unk(isoln,i,j,k)=unk(isoln,i,j,k) + alpha*unk(bcyi,i,j,k) 
               enddo
            enddo
         enddo
         ! Release pointers:
         call Grid_releaseBlkPtr(blockID,unk,CENTER)
      enddo

      ! Write residual and exit:
      if ((bipcgs_print_norm) .and. (MyPE == MasterPE)) then
         write(*,'(a,i4,2(a,es9.2))')                    &
              'cyc1/2', iter,                            &
              ' res norm = ',         norm_new,          &
              ' new to old ratio = ', norm_new/norm_old
      endif

      goto 919

   endif

   if (gr_meshMe .eq. 0) CALL SYSTEM_CLOCK(TA(1),count_rate)

   ! Solve M*z =s
   call gr_bicgApplyPrecond(bczi, bcsi, bcTypes, bcValues)

   if (gr_meshMe .eq. 0) then
      CALL SYSTEM_CLOCK(TA(2),count_rate)
      ET=ET+REAL(TA(2)-TA(1))/count_rate
   endif


   ! Add to Xi: xi = xi-1 + alpha*yi    
   do lb = 1, blockCount
      blockID = blockList(lb)
      ! Point to blocks center vars:
      call Grid_getBlkPtr(blockID,unk,CENTER)
      do k=kli,kui
         do j=jli,jui
            do i=ili,iui
               unk(isoln,i,j,k)=unk(isoln,i,j,k) + alpha*unk(bcyi,i,j,k)
            enddo
         enddo
      enddo
      ! Release pointers:
      call Grid_releaseBlkPtr(blockID,unk,CENTER)
   enddo

   ! t = A*z  !!! From now to end of loop t is stored in bcyi !!!
   !gcellflg = .true. or .false. done on gr_bicgApplyPrecond
   call Timers_start("gr_bicgMultiAx")
   call gr_bicgMultiAx (bczi, bcyi, gcellflg)
                               !ti
   call Timers_stop("gr_bicgMultiAx")

   ! wi = (t,s)/(t,t)
   call Timers_start("gr_bicgDotProd")
   call gr_bicgDotprod(bcyi,bcsi, numer)
   call Timers_stop("gr_bicgDotProd")

   call Timers_start("gr_bicgDotProd")
   call gr_bicgDotprod(bcyi,bcyi, denom)
   call Timers_stop("gr_bicgDotProd")

   wi = numer / denom
   
   ! xi = xi + wi*zi, ri = si - wi*t
   do lb = 1, blockCount
      blockID = blockList(lb)
      ! Point to blocks center vars:
      call Grid_getBlkPtr(blockID,unk,CENTER)
      do k=kli,kui
         do j=jli,jui
            do i=ili,iui
               unk(isoln,i,j,k)=unk(isoln,i,j,k) + wi*unk(bczi,i,j,k)
               unk(bcri,i,j,k) =unk(bcsi,i,j,k)  - wi*unk(bcyi,i,j,k)
           enddo
         enddo
      enddo
      ! Release pointers:
      call Grid_releaseBlkPtr(blockID,unk,CENTER)
   enddo

   ! Check Convergence:
   call gr_bicgNorm (bcri, norm_new)
   norm_new = norm_new / norm_rhs

   ! Write residual:
   if ((bipcgs_print_norm) .and. (MyPE == MasterPE)) then
      write(*,'(a,i4,2(a,es9.2))')                    &
           'cycle ', iter,                            &
           ' res norm = ',         norm_new,          &
           ' new to old ratio = ', norm_new/norm_old
   endif

   if ((norm_new .le. bipcgs_max_residual_norm) .or.     &
       (iter .eq. bipcgs_max_cycles)) then  ! Check convergence.   
      ! Write residual and exit:
      if ((bipcgs_print_norm) .and. (MyPE == MasterPE)) then
         write(*,'(a,i4,2(a,es9.2))')                    &
              'cycle ', iter,                            &
              ' res norm = ',         norm_new,          &
              ' new to old ratio = ', norm_new/norm_old
      endif
      goto 919
   endif
    
   rho1     = rho
   norm_old = norm_new

enddo

919 continue

! Final Guardcell Fill
gcellflg = .true.
call gr_bicgBndry(isoln, nguard_work, gcellflg) ! Guardcells are filled in work array
! Copy solution from work to unk(isoln,:,:,:)
do lb = 1, blockCount
   blockID = blockList(lb)
   ! Point to blocks center vars:
   call Grid_getBlkPtr(blockID,unk,CENTER)
   do k=kle,kue
      do j=jle,jue
         do i=ile,iue
            unk(isoln,i,j,k)=work(i,j,k,blockID,1)  
         enddo
      enddo
   enddo
   ! Release pointers:
   call Grid_releaseBlkPtr(blockID,unk,CENTER)
enddo

! Substract mean value from solution:
if (bicg_bnd_cond == 0) then ! Substract solution mean.
  call gr_findMean(isoln,2,.false.,mean_soln)
  do lb = 1, blockCount
      blockID = blockList(lb)
      ! Point to blocks center vars:
      call Grid_getBlkPtr(blockID,unk,CENTER)
      do k=kle,kue
         do j=jle,jue
            do i=ile,iue 
               unk(isoln,i,j,k) = unk(isoln,i,j,k) - mean_soln
            enddo
         enddo
      enddo
      ! Release pointers:
      call Grid_releaseBlkPtr(blockID,unk,CENTER)
  enddo
endif

interp_mask_work = interp_mask_work_save


if (gr_meshMe .eq. 0) write(*,*) 'Elapsed Time in Preconditioner =',ET

if (firstCall) then
   if (MyPE .eq. MASTER_PE) write(*,*) 'gr_bipcgstab first call done'
   firstCall = .false.
endif


end subroutine gr_bipcgstab
