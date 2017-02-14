!!****if* source/Grid/GridMain/Grid_smoothVar
!!
!!  NAME 
!!
!!  Grid_smoothVar
!!
!!  SYNOPSIS
!!
!! 
!!  call Grid_smoothVar(integer(in) :: ivar, 
!!                           integer(in) :: ivarOut,
!!                           real(INOUT) :: solnData(:,lbUI:,lbUJ:,lbUK:),
!!                           integer(in),value :: lbIU,lbUJ,lbUK,
!!                  OPTIONAL,integer(in) :: smoothMethod,
!!                  OPTIONAL,integer(IN) :: gcLayers,
!!                  OPTIONAL,integer(IN) :: blockID ,
!!                  OPTIONAL,logical(IN) :: useMinSmoothVarVal,
!!                  OPTIONAL,real(IN)    :: minSmoothVarVal,
!!                  OPTIONAL,logical(IN) :: useMaxSmoothVarVal,
!!                  OPTIONAL,real(IN)    :: maxSmoothVarVal,
!!                  OPTIONAL,real(IN)    :: smoothCoeff )
!!
!!  DESCRIPTION 
!!      Smooths one variable in the data passsed in as array solnData.
!!
!! ARGUMENTS
!!
!!   ivar  : index into solution vector, indicating a variable to be smoothed.
!!   ivarOut: index into solution vector, indicating where the smoothed
!!            variable is to be stored.
!!            May be the same as ivar.
!!   solnData : The block data on which to operate
!!   lbUI,lbUJ,lbUK: Lower bounds of assumed-shape array solnData
!!   smoothMethod  : Method for smoothing.
!!   blockID  : Identifies the block on to whic hdata belongs; optional,
!!              for debugging.
!!   gcLayers : In how many layers of guard cells, in addition to the
!!              interior cells, the smoothed output is desired.
!!              Note that one more layer is required in the variable given
!!              by ivar for input to the smoothing.
!!
!! SIDE EFFECTS
!!
!!  Modifies a variable (as requested by argument ivarOut) in the solution
!!  storage UNK.
!!
!!***
subroutine Grid_smoothVar(ivar, ivarOut, &
     solnData, lbUI,lbUJ,lbUK, &
     blklim, smoothMethod, gcLayers, blockID,&
     useMinSmoothVarVal,&
     minSmoothVarVal,&
     useMaxSmoothVarVal,&
     maxSmoothVarVal,&
     smoothCoeff )

!!$  use Grid_interface, ONLY: Grid_getBlkIndexLimits, Grid_getDeltas, Grid_getCellCoords
  use Driver_interface, ONLY: Driver_abortFlash
  implicit none

#include "Flash.h"
#include "constants.h"
#include "FortranLangFeatures.fh"
  
  ! Arguments:
  integer, intent(in) :: ivar
  integer, intent(in) :: ivarOut
  integer, VALUE_INTENT(IN) :: lbUI,lbUJ,lbUK
  real,    intent(INOUT) :: solnData(:,lbUI:,lbUJ:,lbUK:)
  integer, intent(in)    :: blklim(2,MDIM)
  integer, intent(in),OPTIONAL :: smoothMethod
  integer, intent(IN),OPTIONAL :: gcLayers
  integer, intent(IN),OPTIONAL :: blockID
  logical, intent(IN),OPTIONAL :: useMinSmoothVarVal,useMaxSmoothVarVal
  real   , intent(IN),OPTIONAL :: minSmoothVarVal,maxSmoothVarVal
  real   , intent(IN),OPTIONAL :: smoothCoeff

  integer :: smoothMethodLoc
  logical :: useMinSmoothVarValLoc,useMaxSmoothVarValLoc
  real    :: minSmoothVarValLoc,maxSmoothVarValLoc
  real    :: smoothCoeffLoc
  integer :: lb, i, j, k
  integer :: ng
  real :: sval, corr
  real :: fl, fp1, fm1
  real :: delta(MDIM) ! The cell width
  integer :: outLim(2,MDIM)
  real :: maggrad ! The magnitude of the gradient
  real, allocatable :: xcent(:)
  real, allocatable :: ycent(:)
  real, allocatable :: zcent(:)
  real :: gradi, gradj, gradk
  real :: R, lambda3, lambda
  real,pointer :: s(:,:,:)

  smoothMethodLoc = SMOOTH_SOR
  if (present(smoothMethod)) smoothMethodLoc = smoothMethod
  useMinSmoothVarValLoc = .FALSE.
  if (present(useMinSmoothVarVal)) useMinSmoothVarValLoc = useMinSmoothVarVal
  minSmoothVarValLoc = 0.0
  if (present(minSmoothVarVal)) minSmoothVarValLoc = minSmoothVarVal
  useMaxSmoothVarValLoc = .FALSE.
  if (present(useMaxSmoothVarVal)) useMaxSmoothVarValLoc = useMaxSmoothVarVal
  maxSmoothVarValLoc = 1.0
  if (present(maxSmoothVarVal)) maxSmoothVarValLoc = maxSmoothVarVal
  smoothCoeffLoc = 1.0
  if (present(smoothCoeff)) then
     smoothCoeffLoc = smoothCoeff / real(2*NDIM)
  else
     smoothCoeffLoc = 1.0 / real(2*NDIM)
  endif

  if(smoothMethodLoc == SMOOTH_NONE) then
     if (ivarOut .NE. ivar) then
        solnData(ivarOut,:,:,:) = solnData(ivar,:,:,:)
     end if
     return
  end if

  if (present(gcLayers)) then
     ng = gcLayers
  else
     ng = 0
  end if

!!$  call Grid_getBlkIndexLimits(blockID, blklim, blklimgc)
!!$  call Grid_getDeltas(blockID, delta)
!!$
!!$  allocate(xcent(blklimgc(HIGH, IAXIS)))
!!$  call Grid_getCellCoords(IAXIS, blockID, CENTER, .true., &
!!$       xcent, blklimgc(HIGH, IAXIS))
!!$
!!$  allocate(ycent(blklimgc(HIGH, JAXIS)))
!!$  call Grid_getCellCoords(JAXIS, blockID, CENTER, .true., &
!!$       ycent, blklimgc(HIGH, JAXIS))
!!$
!!$  allocate(zcent(blklimgc(HIGH, KAXIS)))          
!!$  call Grid_getCellCoords(KAXIS, blockID, CENTER, .true., &
!!$       zcent, blklimgc(HIGH, KAXIS))

  outLim = blkLim
  outLim(LOW,IAXIS)  = outLim(LOW,IAXIS)  - ng
  outLim(HIGH,IAXIS) = outLim(HIGH,IAXIS) + ng
#if NDIM > 1
  outLim(LOW,JAXIS)  = outLim(LOW,JAXIS)  - ng
  outLim(HIGH,JAXIS) = outLim(HIGH,JAXIS) + ng
#endif
#if NDIM > 2
  outLim(LOW,KAXIS)  = outLim(LOW,KAXIS)  - ng
  outLim(HIGH,KAXIS) = outLim(HIGH,KAXIS) + ng
#endif

  if (ivarOut == ivar) then
     allocate(s(outLim(LOW,IAXIS):outLim(HIGH,IAXIS), &
                outLim(LOW,JAXIS):outLim(HIGH,JAXIS), &
                outLim(LOW,KAXIS):outLim(HIGH,KAXIS)  ))
  else

     call AssoPtr(s,solnData(ivarOut,outLim(LOW,IAXIS):outLim(HIGH,IAXIS), &
                                     outLim(LOW,JAXIS):outLim(HIGH,JAXIS), &
                                     outLim(LOW,KAXIS):outLim(HIGH,KAXIS)  ), &
                    outLim(LOW,IAXIS), &
                    outLim(LOW,JAXIS), &
                    outLim(LOW,KAXIS)  )
  end if

!!$  do k = blklim(LOW,KAXIS)-ng*K3D, blklim(HIGH,KAXIS)+ng*K3D
!!$     do j = blklim(LOW,JAXIS)-ng*K2D, blklim(HIGH,JAXIS)+ng*K2D
!!$        do i = blklim(LOW,IAXIS)-ng, blklim(HIGH,IAXIS)+ng
  do k = outLim(LOW,KAXIS), outLim(HIGH,KAXIS)
     do j = outLim(LOW,JAXIS), outLim(HIGH,JAXIS)
        do i = outLim(LOW,IAXIS), outLim(HIGH,IAXIS)

!!$           ! Compute the magnitude of the gradient:
!!$           fm1 = accessFuncData(ifunc, i-1, j, k)
!!$           fp1 = accessFuncData(ifunc, i+1, j, k)
!!$           gradi = ((fp1-fm1)/(2.0*delta(IAXIS)))**2
!!$           maggrad = gradi
!!$
!!$#if NDIM >= 2
!!$           fm1 = accessFuncData(ifunc, i, j-1, k)
!!$           fp1 = accessFuncData(ifunc, i, j+1, k)
!!$           if(gr_geometry == POLAR) then
!!$              gradj = ((fp1-fm1)/(2.0*delta(JAXIS)))**2 / xcent(i)
!!$           else
!!$              gradj = ((fp1-fm1)/(2.0*delta(JAXIS)))**2
!!$           end if
!!$           maggrad = maggrad + gradj
!!$#endif
!!$
!!$#if NDIM == 3
!!$           fm1 = accessFuncData(ifunc, i, j, k-1)
!!$           fp1 = accessFuncData(ifunc, i, j, k+1)
!!$           gradk = ((fp1-fm1)/(2.0*delta(KAXIS)))**2
!!$           maggrad = maggrad + gradk
!!$#endif
!!$
!!$           maggrad = sqrt(maggrad)

           if (useMaxSmoothVarValLoc) then
              
              if (maxval(solnData(ivar,i-1:i+1,j-K2D:j+K2D,k-K3D:k+K3D))>maxSmoothVarValLoc) then
                 sval = solnData(ivar,i, j, k)
                 s(i, j, k) = sval
                 CYCLE
              end if
           end if
           if (useMinSmoothVarValLoc) then
              if (minval(solnData(ivar,i-1:i+1,j-K2D:j+K2D,k-K3D:k+K3D))<minSmoothVarValLoc) then
                 sval = solnData(ivar,i, j, k)
                 s(i, j, k) = sval
                 CYCLE
              end if
           end if

           select case(smoothMethodLoc)
           case(SMOOTH_3POINT)
              sval =  solnData(ivar,i-1, j, k) + solnData(ivar,i, j, k) + solnData(ivar,i+1, j, k)
#if NDIM > 1
              sval = sval + solnData(ivar,i,j-1,k) + solnData(ivar,i,j,k) + solnData(ivar,i,j+1,k)
#endif
#if NDIM == 3
              sval = sval + solnData(ivar,i,j,k-1) + solnData(ivar,i,j,k) + solnData(ivar,i,j,k+1)
#endif
              sval = sval / (3.0*real(NDIM))
              
           case(SMOOTH_3CPOINT)
              sval =  sum(solnData(ivar,i-1:i+1,j-K2D:j+K2D,k-K3D:k+K3D))
              sval = sval / real(3**(NDIM))
              
           case(SMOOTH_SOR)
              corr = real(-2*NDIM)*solnData(ivar,i, j, k) + solnData(ivar,i-1, j, k) + solnData(ivar,i+1, j, k)
#if NDIM > 1
              corr = corr + solnData(ivar,i,j-1,k) + solnData(ivar,i,j+1,k)
#endif
#if NDIM == 3
              corr = corr + solnData(ivar,i,j,k-1) + solnData(ivar,i,j,k+1)
#endif
              sval =  solnData(ivar,i, j, k) + smoothCoeffLoc * corr

           case(SMOOTH_HARMONIC_SOR)
              corr = real(-2*NDIM)*accessFuncData(ivar,i, j, k) + accessFuncData(ivar,i-1, j, k) + accessFuncData(ivar,i+1, j, k)
#if NDIM > 1
              corr = corr + accessFuncData(ivar,i,j-1,k) + accessFuncData(ivar,i,j+1,k)
#endif
#if NDIM == 3
              corr = corr + accessFuncData(ivar,i,j,k-1) + accessFuncData(ivar,i,j,k+1)
#endif
              sval =  accessFuncData(ivar,i, j, k) + smoothCoeffLoc * corr
              sval = 1.0 / sval

           case DEFAULT
              call Driver_abortFlash("Invalid Smoothing method")
           end select

           s(i, j, k) = sval
        enddo
     enddo
  enddo

!!$  deallocate(xcent)
!!$  deallocate(ycent)
!!$  deallocate(zcent)

  if (ivarOut == ivar) then
     solnData(ivarOut,outLim(LOW,IAXIS):outLim(HIGH,IAXIS), &
                      outLim(LOW,JAXIS):outLim(HIGH,JAXIS), &
                      outLim(LOW,KAXIS):outLim(HIGH,KAXIS)  ) &
                  = s(outLim(LOW,IAXIS):outLim(HIGH,IAXIS), &
                      outLim(LOW,JAXIS):outLim(HIGH,JAXIS), &
                      outLim(LOW,KAXIS):outLim(HIGH,KAXIS)  )
     deallocate(s)
  end if

contains
  real function accessFuncData(ifunc,i,j,k)
    integer,intent(IN) :: ifunc,i,j,k
    if (smoothMethodLoc==SMOOTH_HARMONIC_SOR) then
       accessFuncData = 1.0/solnData(ifunc,i,j,k)
    else
       accessFuncData = solnData(ifunc,i,j,k)
    end if
  end function accessFuncData
  subroutine AssoPtr(p, d, lbI,lbJ,lbK)
    real,POINTER_INTENT_OUT :: p(:,:,:)
    integer,VALUE_INTENT(in) :: lbI,lbJ,lbK
    real,   intent(in),target :: d(lbI:,lbJ:,lbK:)
!    print*,'AssoPtr: LBOUND(d)',LBOUND(d)
    p => d
!    print*,'AssoPtr: LBOUND(p)',LBOUND(p)
    return
  end subroutine AssoPtr
end subroutine Grid_smoothVar
