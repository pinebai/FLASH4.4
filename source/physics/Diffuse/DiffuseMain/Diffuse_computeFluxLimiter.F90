!!****if* source/physics/Diffuse/DiffuseMain/Diffuse_computeFluxLimiter
!!
!!  NAME 
!!
!!  Diffuse_computeFluxLimiter
!!
!!  SYNOPSIS
!!
!! 
!!  call Diffuse_computeFluxLimiter(integer(in) :: idcoef, 
!!                           integer(in) :: ifunc,
!!                           integer(in) :: ifl,
!!                           integer(in) :: iflOut,
!!                           integer(in) :: ieddi3,
!!                           integer(in) :: mode,
!!                           real(INOUT) :: solnData(:,lbUI:,lbUJ:,lbUK:),
!!                           integer(in),value :: lbIU,lbUJ,lbUK,
!!                           integer(IN) :: blockID,
!!                  OPTIONAL,integer(IN) :: gcLayers )
!!
!!  DESCRIPTION 
!!      This routine computes flux limiter factors and returns them
!!      in UNK variables. It does not apply the limiter by modifying
!!      diffusion coefficients or in any other way.
!!
!! ARGUMENTS
!!
!!   idcoef  : index into solution vector, indicating a variable that holds
!!             the coefficient to which the limiter is applied.
!!             This variable is used only for input.
!!   ifunc   : index into solution vector giving the quantity whose flux
!!             is to be limited. This gives either the quantity directly
!!             if iDenForIfunc is < 1, or must be multiplied by the variable
!!             given by iDenForIfunc.
!!             This variable is used only for input.
!!   iDenForIfunc : A density to multiply ifunc by, or -1.
!!             This variable is used only for input.
!!   ifl     : index into solution vector giving the limiting value that
!!             would be applied to the variable given by idcoef such that
!!             the latter is less than or equal to the former.
!!             This variable is used only for input.
!!   iflOut  : index into solution vector to put flux limiter factor.
!!             Note that this may be one of the variables given by the
!!             previous four dummy arguments that are used as inputs.
!!   mode     : Flux limiter mode.
!!   solnData : The block data on which to operate
!!   lbUI,lbUJ,lbUK: lower bounds of assumed-shape array solnData
!!   blockID  : Identifies the block on which to operate
!!   gcLayers : In how many layers of guard cells, in additionm to the
!!              interior cells, the flux limiter coefficients are desired.
!!              Note that one more layer is required in the variable given
!!              by ifl for input, for computing the gradient.
!!
!! SIDE EFFECTS
!!
!!  Modifies flux limiter variables (as requested by argument iflOut) in the solution
!!  storage UNK.  On output it will be overwritten with a value
!! three times the flux limiter factor lambda, which should be
!! familiar from the literature on flux limiters. The value
!! of lambda should be in the range [ 0, 1/3 ].
!!
!!***
subroutine Diffuse_computeFluxLimiter(idcoef, ifunc, iDenForIfunc, ifl, iflOut, ieddi3, mode, &
     solnData, lbUI,lbUJ,lbUK, &
     blockID, gcLayers)
  use Grid_interface, ONLY: Grid_getBlkIndexLimits, Grid_getDeltas, Grid_getCellCoords
  use Diffuse_data, ONLY: diff_meshMe, diff_geometry
  use Driver_interface, ONLY: Driver_abortFlash
  implicit none

#include "Flash.h"
#include "constants.h"
#include "Eos.h"
#include "FortranLangFeatures.fh"
  
  ! Arguments:
  integer, intent(in) :: idcoef
  integer, intent(in) :: ifunc
  integer, intent(in) :: iDenForIfunc
  integer, intent(in) :: ifl
  integer, intent(in) :: iflOut
  integer, intent(in) :: ieddi3
  integer, intent(in) :: mode
  integer, VALUE_INTENT(IN) :: lbUI,lbUJ,lbUK
  real,    intent(INOUT) :: solnData(:,lbUI:,lbUJ:,lbUK:)
  integer, intent(IN) :: blockID
  integer, intent(IN),OPTIONAL :: gcLayers

  integer :: lb, i, j, k
  integer :: blklim(2,MDIM), blklimgc(2,MDIM)
  integer :: ng
  real :: dcoef, dcoefOld
  real :: fl, fp1, fm1
  real :: delta(MDIM) ! The cell width
  real :: maggrad ! The magnitude of the gradient
  real, allocatable :: xcent(:)
  real, allocatable :: ycent(:)
  real, allocatable :: zcent(:)
  real :: gradi, gradj, gradk
  real :: R, lambda3, lambda
  real, parameter :: hugeReal = 0.25*HUGE(R)
  real, save      :: sqrtHuge = 0.0
  logical :: useIdenForFunc

  if(mode == FL_NONE) then
     solnData(iflOut,:,:,:) = 1.0
     if (ieddi3 > 0) solnData(ieddi3, :,:,:) = 1.0
     return
  end if

  if (sqrtHuge==0.0) sqrtHuge = sqrt(hugeReal) !  Initialize this on first call.

  if (iDenForIfunc > 0) then
     useIdenForFunc = .TRUE.
  else
     useIdenForFunc = .FALSE.
  end if

  if (present(gcLayers)) then
     ng = gcLayers
  else
     ng = 0
  end if

  call Grid_getBlkIndexLimits(blockID, blklim, blklimgc)
  call Grid_getDeltas(blockID, delta)

  allocate(xcent(blklimgc(HIGH, IAXIS)))
  call Grid_getCellCoords(IAXIS, blockID, CENTER, .true., &
       xcent, blklimgc(HIGH, IAXIS))

  allocate(ycent(blklimgc(HIGH, JAXIS)))
  call Grid_getCellCoords(JAXIS, blockID, CENTER, .true., &
       ycent, blklimgc(HIGH, JAXIS))

  allocate(zcent(blklimgc(HIGH, KAXIS)))          
  call Grid_getCellCoords(KAXIS, blockID, CENTER, .true., &
       zcent, blklimgc(HIGH, KAXIS))

  do k = blklim(LOW,KAXIS)-ng*K3D, blklim(HIGH,KAXIS)+ng*K3D
     do j = blklim(LOW,JAXIS)-ng*K2D, blklim(HIGH,JAXIS)+ng*K2D
        do i = blklim(LOW,IAXIS)-ng, blklim(HIGH,IAXIS)+ng

           ! Compute the magnitude of the gradient:
           fm1 = accessFuncData(ifunc, i-1, j, k)
           fp1 = accessFuncData(ifunc, i+1, j, k)
           gradi = ((fp1-fm1)/(2.0*delta(IAXIS)))**2
           maggrad = gradi

#if NDIM >= 2
           fm1 = accessFuncData(ifunc, i, j-1, k)
           fp1 = accessFuncData(ifunc, i, j+1, k)
           if(diff_geometry == POLAR) then
              gradj = ((fp1-fm1)/(2.0*delta(JAXIS)))**2 / xcent(i)
           else
              gradj = ((fp1-fm1)/(2.0*delta(JAXIS)))**2
           end if
           maggrad = maggrad + gradj
#endif

#if NDIM == 3
           fm1 = accessFuncData(ifunc, i, j, k-1)
           fp1 = accessFuncData(ifunc, i, j, k+1)
           gradk = ((fp1-fm1)/(2.0*delta(KAXIS)))**2
           maggrad = maggrad + gradk
#endif

           maggrad = sqrt(maggrad)

           dcoefOld = solnData(idcoef, i, j, k)
           if (dcoefOld == 0.0) then
              print*,'HELP! dcoefOld is 0.0!!! i=',i
              solnData(iflOut,i,j,k) = 1.0
              if (ieddi3 > 0) solnData(ieddi3, :,:,:) = 1.0 !???
              CYCLE
           end if

           ! Now compute the new diffusion coefficient:
           fl    = solnData(ifl, i, j, k)

           if (mode==FL_LEVPOM .OR. ieddi3 > 0) then ! need R in those cases...
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
           end if

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


              if (R > sqrtHuge) then
                 lambda3 = 3.0/R
              else
                 lambda3 = 3.0*(2.0 + R)/(6.0+3.0*R+R**2)
              endif

           case DEFAULT
              call Driver_abortFlash("Invalid Flux limiter type")
           end select

           if (ieddi3 > 0) then  ! Is this right only for FL_LEVPOM ?
              lambda = lambda3 / 3.0
!!$              solnData(ieddi, i, j, k) = 0.5 - (lambda * (0.5 + 0.5*lambda * R*R))
              solnData(ieddi3, i, j, k) = 1.5 - (lambda3 * (0.5 + 0.5*lambda * R*R))
!!$              print*,i,lambda3,lambda,R,solnData(ieddi3, i, j, k),ieddi3
           end if

           ! Overwrite iflOut, normally FLLM_VAR, with 3 times the radiation
           ! flux limiter factor lambda,
           ! for diagnostics.
           solnData(iflOut, i, j, k) = lambda3
        enddo
     enddo
  enddo

  deallocate(xcent)
  deallocate(ycent)
  deallocate(zcent)

contains
  pure real function accessFuncData(ifunc,i,j,k)
    integer,intent(IN) :: ifunc,i,j,k
    if (useIdenForFunc) then
       accessFuncData = solnData(ifunc,i,j,k) * solnData(iDenForIfunc,i,j,k)
    else
       accessFuncData = solnData(ifunc,i,j,k)
    end if
  end function accessFuncData
end subroutine Diffuse_computeFluxLimiter
