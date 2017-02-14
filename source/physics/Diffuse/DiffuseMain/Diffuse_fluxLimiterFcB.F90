!!****if* source/physics/Diffuse/DiffuseMain/Diffuse_fluxLimiterFcB
!!
!!  NAME 
!!
!!  Diffuse_fluxLimiterFcB
!!
!!  SYNOPSIS
!!
!! 
!!  call Diffuse_fluxLimiterFcB(integer(in) :: idcoef, 
!!                           integer(in) :: ifunc,
!!                           integer(in) :: ifl,
!!                           integer(in) :: mode,
!!                           integer(IN) :: blkcnt,
!!                           integer(IN) :: blklst(blkcnt),
!!                           logical(IN) :: returnCcB)
!!
!!  DESCRIPTION 
!!      This routine modifies the diffusion coefficient (idcoef) 
!!      and applies limiting.
!!
!! ARGUMENTS
!!
!!   idcoef  : index into solution vector, indicating a variable that holds
!!             the coefficient to which the limiter is applied.
!!             This variable is used for input and possibly also for output.
!!   ifunc   : index into solution vector giving the quantity whose flux
!!             is to be limited.
!!             This variable is used only for input.
!!   ifl     : index into solution vector giving the limiting values that
!!             are applied to the values of the variable given by
!!             idcoef so that the latter become less than or equal to the
!!             former.
!!             Overwritten with a cell-centered version of 3 times lambda
!!             (the flux limiter factor) if returnFlFactor is true.
!!   mode    : Flux limiter mode.
!!   blkcnt  : The number of blocks in the list
!!   blklst  : The list of blocks on which the solution must be updated.
!!   returnCcB : Indicates whether a modified (i.e., flux-limited)
!!               cell-centered version of the coefficient should be
!!               returned in the solution vector (in addition to the
!!               face-centered version which is always returned in
!!               the first variable of face-centered allocatable
!!               scratch buffers).
!!
!! SIDE EFFECTS
!!
!! Modifies the diffusion coefficient variable, indicated by the
!! argument idcoef, if returnCcB is true.
!!
!! Modifies the first two variables in face-centered allocatable
!! scratch buffers (gds type FACEX, FACEY, FACEZ), see
!! Grid_ascGetBlkPtr calls.
!!
!! If returnFlFactor=.TRUE., the flux limiter variable, indicated by
!! the argument ifl, is modified.
!! On output it will be overwritten with a value three times the flux
!! limiter factor lambda, which should be familiar from the literature
!! on flux limiters.  The value of a flux limiter factor lambda should
!! be in the range ( 0, 1/3 ].  Thus the value returned in variable
!! ifl of the solution vector will be in the range [ 0.0, 1.0 ].
!!
!! NOTES
!!
!!  If returnCcB is true, results are returned for interior cells
!!  and additionally one layer of guard cells.
!!  If returnCcB is false, results are only returned for the interior
!!  cells (including all their faces, for the face-centered results).
!!
!!  If returnCcB is true, valid input data in ifunc is required in two
!!  layers of guard cells.  If returnCcB is one, valid input data in
!!  ifunc is required in one layers of guard cells.  The caller is
!!  responsible for ensuring that these requirements are met, for
!!  example by calling Grid_fillGuardCells before calling this routine.
!!
!! SEE ALSO
!!
!!  Grid_ascGetBlkPtr
!!
!!***
subroutine Diffuse_fluxLimiterFcB(idcoef, ifunc, ifl, mode, blkcnt, blklst,returnCcB)
  use Grid_interface, ONLY: Grid_getBlkPtr, Grid_releaseBlkPtr, &
       Grid_ascGetBlkPtr,   Grid_ascReleaseBlkPtr, &
       Grid_getBlkIndexLimits, Grid_getDeltas
  use Diffuse_data, ONLY: diff_meshMe, diff_geometry
  use Driver_interface, ONLY: Driver_abortFlash
  implicit none

#include "Flash.h"
#include "constants.h"
#include "Eos.h"
  
  ! Arguments:
  integer, intent(in) :: idcoef
  integer, intent(in) :: ifunc
  integer, intent(in) :: ifl
  integer, intent(in) :: mode
  integer, intent(IN) :: blkcnt
  integer, intent(IN) :: blklst(blkcnt)
  logical, intent(IN) :: returnCcB

  integer :: lb, i, j, k
  real, pointer :: blkPtr(:,:,:,:), facBptr(:,:,:,:)
  integer :: blklim(2,MDIM), blklimgc(2,MDIM)
  integer :: dir, ioff, joff, koff, fgds
  real :: fl, fp1, fm1
  real :: delta(MDIM) ! The cell width
  real :: maggrad ! The magnitude of the gradient
  real, allocatable :: xcent(:), xfac(:)
  real, allocatable :: ycent(:), yfac(:)
  real, allocatable :: zcent(:)
  real :: gradi, gradj, gradk
  real :: R
  real :: dcoefFcB, fllmFcB, dcoefOldFcB, fllmOldFcB
  real, parameter :: hugeReal = 0.25*HUGE(R)
  real, save      :: sqrtHuge = 0.0
  logical, parameter :: returnFlFactor = .TRUE.
  integer, parameter :: iFactorB = 1,             iflFcB = 2
  integer :: kk1d, kk2d, kk3d

  kk1d = 0; kk2d = 0; kk3d = 0
  if (returnCcB) then
     kk1d = 1; kk2d = K2D; kk3d = K3D
  end if
  if (sqrtHuge==0.0) sqrtHuge = sqrt(hugeReal) !  Initialize this on first call.

  do lb = 1, blkcnt
     call Grid_getBlkIndexLimits(blklst(lb), blklim, blklimgc)
     call Grid_getBlkPtr(blklst(lb), blkPtr)
     call Grid_getDeltas(blklst(lb), delta)
     
     allocate(xcent(blklimgc(HIGH, IAXIS)))
     call Grid_getCellCoords(IAXIS, blklst(lb), CENTER, .true., &
          xcent, blklimgc(HIGH, IAXIS))

     if(diff_geometry == POLAR .AND. NDIM > 1) then
        allocate(xfac(blklim(LOW, IAXIS):blklim(HIGH, IAXIS)))
        call Grid_getCellCoords(IAXIS, blklst(lb), FACES, .false., &
             xfac, blklim(HIGH, IAXIS)-blklim(LOW, IAXIS)+1)
     end if

#if NDIM >= 2
     allocate(ycent(blklimgc(HIGH, JAXIS)))
     call Grid_getCellCoords(JAXIS, blklst(lb), CENTER, .true., &
          ycent, blklimgc(HIGH, JAXIS))
#endif
#if NDIM == 3
     allocate(zcent(blklimgc(HIGH, KAXIS)))          
     call Grid_getCellCoords(KAXIS, blklst(lb), CENTER, .true., &
          zcent, blklimgc(HIGH, KAXIS))
#endif

     do dir=IAXIS,NDIM
        ioff = 0; joff = 0; koff = 0
        if (dir==IAXIS) then
           ioff = 1; fgds = FACEX
        else if (dir==JAXIS) then
           joff = 1; fgds = FACEY
        else
           koff = 1; fgds = FACEZ
        end if
        call Grid_ascGetBlkPtr(blklst(lb),facBptr,fgds)
        do k = blklim(LOW,KAXIS)-kk3d, blklim(HIGH,KAXIS)+kk3d+koff
           do j = blklim(LOW,JAXIS)-kk2d, blklim(HIGH,JAXIS)+kk2d+joff
              do i = blklim(LOW,IAXIS)-kk1d, blklim(HIGH,IAXIS)+kk1d+ioff

!!$                 call diff_computMaggradLocal(maggrad,dir,i,j,k)
                 ! Compute the magnitude of the gradient:
                 if (dir==IAXIS) then
                    fm1 = blkPtr(ifunc, i-1, j, k)
                    fp1 = blkPtr(ifunc, i  , j, k)
                 else
                    fm1 = blkPtr(ifunc, i-1, j, k) + blkPtr(ifunc, i-1, j-joff, k-koff)
                    fp1 = blkPtr(ifunc, i+1, j, k) + blkPtr(ifunc, i+1, j-joff, k-koff)
                 endif
                    
                 gradi = ((fp1-fm1)/((4.0-3*ioff)*delta(IAXIS)))**2
                 maggrad = gradi

#if NDIM >= 2
                 if (dir==JAXIS) then
                    fm1 = blkPtr(ifunc, i, j-1, k)
                    fp1 = blkPtr(ifunc, i, j  , k)
                 else
                    fm1 = blkPtr(ifunc, i, j-1, k) + blkPtr(ifunc, i-ioff, j-1, k-koff)
                    fp1 = blkPtr(ifunc, i, j+1, k) + blkPtr(ifunc, i-ioff, j+1, k-koff)
                 end if
                 if(diff_geometry == POLAR) then
                    if (dir==JAXIS) then
                       gradj = ((fp1-fm1)/(delta(JAXIS) * xcent(i)) )**2
                    else if (dir==IAXIS) then
                       gradj = ((fp1-fm1)/(4.0*delta(JAXIS) * xfac(i) ))**2
                    else
                       gradj = ((fp1-fm1)/(4.0*delta(JAXIS) * xcent(i) ))**2
                    end if
                 else
                    gradj = ((fp1-fm1)/((4.0-3*joff)*delta(JAXIS)))**2
                 end if
                 maggrad = maggrad + gradj
#endif

#if NDIM == 3
                 if (dir==KAXIS) then
                    fm1 = blkPtr(ifunc, i, j, k-1)
                    fp1 = blkPtr(ifunc, i, j, k  )
                 else
                    fm1 = blkPtr(ifunc, i, j, k-1) + blkPtr(ifunc, i-ioff, j-joff, k-1)
                    fp1 = blkPtr(ifunc, i, j, k+1) + blkPtr(ifunc, i-ioff, j-joff, k+1)
                 end if
                 gradk = ((fp1-fm1)/((4.0-3*koff)*delta(KAXIS)))**2
                 maggrad = maggrad + gradk
#endif

                 maggrad = sqrt(maggrad)

                 dcoefOldFcB = 0.5 * ( blkPtr(idcoef,i,j,k) + blkPtr(idcoef,i-ioff,j-joff,k-koff) )
                 fllmOldFcB  = 0.5 * ( blkPtr(ifl   ,i,j,k) + blkPtr(ifl   ,i-ioff,j-joff,k-koff) )
                 dcoefFcB = dcoefOldFcB
                 fllmFcB  = fllmOldFcB
                 call diff_fluxLimiterLocal(maggrad, dcoefFcB, fllmFcB)
                 facBptr(iFactorB,i,j,k) = dcoefFcB
                 if (returnFlFactor) facBptr(iflFcB, i,j,k) = fllmFcB
              enddo
           enddo
        enddo
        call Grid_ascReleaseBlkPtr(blklst(lb),facBptr,fgds)
     end do

     deallocate(xcent)
     if (allocated(xfac))  deallocate(xfac)
     if (allocated(ycent)) deallocate(ycent)
     if (allocated(zcent)) deallocate(zcent)
     
     if (returnFlFactor .OR. returnCcB) then
        do dir=IAXIS,NDIM
           ioff = 0; joff = 0; koff = 0
           if (dir==IAXIS) then
              ioff = 1; fgds = FACEX
           else if (dir==JAXIS) then
              joff = 1; fgds = FACEY
           else
              koff = 1; fgds = FACEZ
           end if
           call Grid_ascGetBlkPtr(blklst(lb),facBptr,fgds)
           do k = blklim(LOW,KAXIS)-kk3d, blklim(HIGH,KAXIS)+kk3d
              do j = blklim(LOW,JAXIS)-kk2d, blklim(HIGH,JAXIS)+kk2d
                 do i = blklim(LOW,IAXIS)-kk1d, blklim(HIGH,IAXIS)+kk1d

                    if (returnFlFactor) then
                       if (dir==IAXIS) then
                          blkPtr (ifl,i,j,k) = 0.0
                       end if
                       blkPtr(ifl,i,j,k) = &
                            blkPtr(ifl,i,j,k) + facBptr(iflFcB,i,j,k) + facBptr(iflFcB,i+ioff,j+joff,k+koff)
                       if (dir==NDIM) then
                          blkPtr (ifl,i,j,k) = blkPtr(ifl,i,j,k) / (2*(NDIM))
                       end if
                    end if

                    if (returnCcB) then
#ifdef RETURNED_CCB_IS_AVG
                       if (dir==IAXIS) then
                          blkPtr (idcoef,i,j,k) = 0.0
                       end if
                       blkPtr(idcoef,i,j,k) = &
                            blkPtr(idcoef,i,j,k) + facBptr(iFactorB,i,j,k) + facBptr(iFactorB,i+ioff,j+joff,k+koff)
                       if (dir==NDIM) then
                          blkPtr (idcoef,i,j,k) = blkPtr(idcoef,i,j,k) / (2*(NDIM))
                       end if
#else
                       ! RETURNED_CCB_IS_MIN
                       blkPtr(idcoef,i,j,k) = &
                            min(blkPtr(idcoef,i,j,k), facBptr(iFactorB,i,j,k), facBptr(iFactorB,i+ioff,j+joff,k+koff))
#endif
                    end if
                 enddo
              enddo
           enddo
           call Grid_ascReleaseBlkPtr(blklst(lb),facBptr,fgds)
        end do
     end if

     call Grid_releaseBlkPtr(blklst(lb), blkPtr)
     
  end do

contains
  subroutine diff_fluxLimiterLocal(maggrad,dcoef,fllm)
    real, intent(in)    :: maggrad
    real, intent(inout) :: dcoef,fllm

    real :: dcoefOld
    real :: lambda3

    dcoefOld = dcoef
    if ( mode == FL_NONE .OR. dcoefOld == 0.0 ) then
       fllm = 1.0
       return
    end if

    ! Now compute the new diffusion coefficient:
    fl    = fllm

    select case(mode)
    case(FL_HARMONIC)
       dcoef = 1.0 / (1.0/dcoefOld + maggrad/(fl + 1.0D-100))
       lambda3 = dcoef / dcoefOld

    case(FL_MINMAX)
       dcoef = min(dcoefOld, fl/max(TINY(maggrad),TINY(maggrad)*fl,maggrad))
       lambda3 = dcoef / dcoefOld

    case(FL_LARSEN)
       dcoef = 1.0 / sqrt((1.0/dcoefOld)**2 + (maggrad/(fl + 1.0D-100))**2)
       lambda3 = dcoef / dcoefOld

    case(FL_LEVPOM)
       ! Using the flux limiter of Levermore & Pomraning (1981)

       if (maggrad == 0.0) then
          R = 0.0
       else if (maggrad .LE. sqrtHuge .AND. dcoefOld .LE. sqrtHuge .AND. &
            TINY(fl)*dcoefOld*maggrad > 1.0/HUGE(fl)) then
          R = dcoefOld * maggrad / max(TINY(fl)*dcoefOld*maggrad,fl) * 3.0
       else if (maggrad .LE. 1.25) then
          R = dcoefOld / max(TINY(fl)*dcoefOld,fl) * 3.0 * maggrad
       else
          R = dcoefOld / max(TINY(fl)*dcoefOld*maggrad,fl) * 3.0 * maggrad
       end if

       if (R > sqrtHuge) then
          lambda3 = 3.0/R
       else
          lambda3 = 3.0*(2.0 + R)/(6.0+3.0*R+R**2)
       endif

       dcoef = lambda3 * dcoefOld

    case DEFAULT
       call Driver_abortFlash("Invalid Flux limiter type")
    end select

!!$    dcoef = dcoef
    ! Overwrite ifl, normally FLLM_VAR, with 3 times the radiation
    ! flux limiter factor lambda,
    ! for diagnostics.
    fllm = lambda3
  end subroutine diff_fluxLimiterLocal
end subroutine Diffuse_fluxLimiterFcB
