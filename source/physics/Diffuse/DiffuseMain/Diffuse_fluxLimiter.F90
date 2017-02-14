!!****if* source/physics/Diffuse/DiffuseMain/Diffuse_fluxLimiter
!!
!!  NAME 
!!
!!  Diffuse_fluxLimiter
!!
!!  SYNOPSIS
!!
!! 
!!  call Diffuse_fluxLimiter(integer(in) :: idcoef, 
!!                           integer(in) :: ifunc,
!!                           integer(in) :: ifl,
!!                           integer(in) :: mode,
!!                           integer(IN) :: blkcnt,
!!                           integer(IN) :: blklst(blkcnt))
!!
!!  DESCRIPTION 
!!      This routine modifies the diffusion coefficient (idcoef) 
!!      and applies limiting.
!!
!! ARGUMENTS
!!
!!   idcoef  : index into solution vector, indicating a variable that holds
!!             the coefficient to which the limiter is applied.
!!             This variable is used for input and output.
!!   ifunc   : index into solution vector giving the quantity whose flux
!!             is to be limited.
!!             This variable is used only for input.
!!   ifl     : index into solution vector giving flux limiter variable.
!!             This variable is used for input and output.
!!   mode    : Flux limiter mode.
!!   blkcnt  : The number of blocks in the list
!!   blklst  : The list of blocks on which the solution must be updated.
!!
!! SIDE EFFECTS
!!
!! Modifies the diffusion coefficient variable, indicated by the
!! argument idcoef.
!!
!! Modifies the flux limiter variable, indicated by the argument ifl.
!! On output it will be overwritten with a value three times the flux
!! limiter factor lambda, which should be familiar from the literature
!! on flux limiters.  The value of a flux limiter factor lambda should
!! be in the range ( 0, 1/3 ].  Thus the value returned in variable
!! ifl of the solution vector will be in the range [ 0.0, 1.0 ].
!!
!! NOTES
!!
!!  Results for idcoef and ifl are returned for interior cells
!!  and additionally one layer of guard cells.
!!
!!  Valid input data in ifunc is required in one layers of guard
!!  cells.  The caller is responsible for ensuring that this
!!  requirement is met, for example by calling Grid_fillGuardCells
!!  before calling this routine.
!!***
subroutine Diffuse_fluxLimiter(idcoef, ifunc, ifl, mode, blkcnt, blklst)
  use Grid_interface, ONLY: Grid_getBlkPtr, Grid_releaseBlkPtr, &
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

  integer :: lb, i, j, k
  real, pointer :: blkPtr(:,:,:,:)
  integer :: blklim(2,MDIM), blklimgc(2,MDIM)
  real :: dcoef, dcoefOld
  real :: fl, fp1, fm1
  real :: delta(MDIM) ! The cell width
  real :: maggrad ! The magnitude of the gradient
  real, allocatable :: xcent(:)
  real, allocatable :: ycent(:)
  real, allocatable :: zcent(:)
  real :: gradi, gradj, gradk
  real :: R, lambda3
  real, parameter :: hugeReal = 0.25*HUGE(R)
  real, save      :: sqrtHuge = 0.0


  if (sqrtHuge==0.0) sqrtHuge = sqrt(hugeReal) !  Initialize this on first call.

  do lb = 1, blkcnt
     call Grid_getBlkPtr(blklst(lb), blkPtr)
     if(mode == FL_NONE) then
        blkPtr(ifl,:,:,:) = 1.0
        call Grid_releaseBlkPtr(blklst(lb), blkPtr)
        CYCLE
     end if
     call Grid_getBlkIndexLimits(blklst(lb), blklim, blklimgc)
     call Grid_getDeltas(blklst(lb), delta)
     
     allocate(xcent(blklimgc(HIGH, IAXIS)))
     call Grid_getCellCoords(IAXIS, blklst(lb), CENTER, .true., &
          xcent, blklimgc(HIGH, IAXIS))

     allocate(ycent(blklimgc(HIGH, JAXIS)))
     call Grid_getCellCoords(JAXIS, blklst(lb), CENTER, .true., &
          ycent, blklimgc(HIGH, JAXIS))

     allocate(zcent(blklimgc(HIGH, KAXIS)))          
     call Grid_getCellCoords(KAXIS, blklst(lb), CENTER, .true., &
          zcent, blklimgc(HIGH, KAXIS))

     do k = blklim(LOW,KAXIS)-K3D, blklim(HIGH,KAXIS)+K3D
        do j = blklim(LOW,JAXIS)-K2D, blklim(HIGH,JAXIS)+K2D
           do i = blklim(LOW,IAXIS)-1, blklim(HIGH,IAXIS)+1

              ! Compute the magnitude of the gradient:
              fm1 = blkPtr(ifunc, i-1, j, k)
              fp1 = blkPtr(ifunc, i+1, j, k)
              gradi = ((fp1-fm1)/(2.0*delta(IAXIS)))**2
              maggrad = gradi

#if NDIM >= 2
              fm1 = blkPtr(ifunc, i, j-1, k)
              fp1 = blkPtr(ifunc, i, j+1, k)
              if(diff_geometry == POLAR) then
                 gradj = ((fp1-fm1)/(2.0*delta(JAXIS)))**2 / xcent(i)
              else
                 gradj = ((fp1-fm1)/(2.0*delta(JAXIS)))**2
              end if
              maggrad = maggrad + gradj
#endif

#if NDIM == 3
              fm1 = blkPtr(ifunc, i, j, k-1)
              fp1 = blkPtr(ifunc, i, j, k+1)
              gradk = ((fp1-fm1)/(2.0*delta(KAXIS)))**2
              maggrad = maggrad + gradk
#endif

              maggrad = sqrt(maggrad)

              dcoefOld = blkPtr(idcoef, i, j, k)
              if (dcoefOld == 0.0) then
            !!     print*,'HELP! dcoefOld is 0.0!!! i=',i
                 blkPtr(ifl,i,j,k) = 1.0
                 CYCLE
              end if

              ! Now compute the new diffusion coefficient:
              fl    = blkPtr(ifl, i, j, k)

              select case(mode)
              case(FL_HARMONIC)
                 dcoef = 1.0 / (1.0/dcoefOld + maggrad/(fl + 1.0D-100))
                 lambda3 = dcoef / dcoefOld

              case(FL_MINMAX)
                 dcoef = min(dcoefOld, fl/(maggrad + 1.0d-100))
                 lambda3 = dcoef / dcoefOld

              case(FL_LARSEN)
                 dcoef = 1.0 / sqrt((1.0/dcoefOld)**2 + (maggrad/(fl + 1.0D-100))**2)
                 lambda3 = dcoef / dcoefOld

              case(FL_LEVPOM)
                 ! Using the flux limiter of Levermore & Pomraning (1981)

                 if (maggrad == 0.0) then
             !!       print*,'HELP! maggrad is 0.0!!! i=',i
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

              blkPtr(idcoef, i, j, k) = dcoef
              ! Overwrite ifl, normally FLLM_VAR, with 3 times the radiation
              ! flux limiter factor lambda,
              ! for diagnostics.
              blkPtr(ifl   , i, j, k) = lambda3
           enddo
        enddo
     enddo
     call Grid_releaseBlkPtr(blklst(lb), blkPtr)
     
     deallocate(xcent)
     deallocate(ycent)
     deallocate(zcent)
     
  end do

end subroutine Diffuse_fluxLimiter
