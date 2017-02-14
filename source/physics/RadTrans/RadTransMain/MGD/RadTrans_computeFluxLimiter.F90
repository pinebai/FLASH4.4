!!****if* source/physics/RadTrans/RadTransMain/MGD/RadTrans_computeFluxLimiter
!!
!! NAME
!!
!!  RadTrans_computeFluxLimiter
!!
!! SYNOPSIS
!!
!!  call RadTrans_computeFluxLimiter(integer(in) :: ifl,
!!                                   integer(in) :: iflOut,
!!                                   integer(in) :: ieddi3,
!!                                   real(INOUT) :: solnData,
!!                                   integer(IN) :: blockID,
!!                                   integer(IN),OPTIONAL :: gcLayers)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   ifl : 
!!
!!   iflOut : 
!!
!!   ieddi3 : 
!!
!!   solnData : 
!!
!!   blockID : ID of block in current processor
!!
!!   gcLayers : 
!!
!! AUTOGENROBODOC
!!
!!
!!***

#include "Flash.h"

subroutine RadTrans_computeFluxLimiter(ifl, iflOut, ieddi3, solnData, blockID,gcLayers)
  use Grid_interface,    ONLY: Grid_getBlkIndexLimits
  use Diffuse_interface, ONLY: Diffuse_computeFluxLimiter
  use Opacity_interface, ONLY: Opacity
  use rt_data, ONLY: rt_useMGD, rt_mgdDomainBC, rt_mgdBounds, rt_mgdFlMode, &
       rt_mgdFlCoef, rt_mgdNumGroups, rt_mgdBcVals, rt_mgdthetaImplct, &
       rt_timeGroups, rt_groupBarrier, rt_gcMaskSize, rt_gcMask
  use RadTrans_data, ONLY: rt_boltz, rt_speedlt, rt_speedlt3, rt_radconst, rt_meshMe, &
       rt_meshCopyCount, rt_acrossMe, rt_globalComm , &
       rt_dbgContext

  implicit none

#include "constants.h"

  ! Arguments:
  integer, intent(in) :: ifl
  integer, intent(in) :: iflOut
  integer, intent(in) :: ieddi3
  real,    intent(INOUT) :: solnData(:,1:,1:,1:)
  integer, intent(IN) :: blockID
  integer, intent(IN),OPTIONAL :: gcLayers

  integer :: blkLimitsGC(LOW:HIGH,MDIM), blkLimits(LOW:HIGH,MDIM)
  integer :: lb, g, gloc, gvar, i, j, k
  integer :: nlayers, ng

  real    :: absorb_opac, emit_opac, trans_opac
  real    :: lambda, lambda3, pMat
  real    :: C  ! Speed of light [cm/s]

  C  = rt_speedlt

  if (present(gcLayers)) then
     nlayers = gcLayers
  else
     nlayers = 0
  end if
  ng = nlayers + 1

  !!2. CALL FLUX LIMITER

  do gloc = 1, NONREP_NLOCS(rt_acrossMe, rt_meshCopyCount, rt_mgdNumGroups)
     gvar = MGDR_NONREP_LOC2UNK(gloc)
     g = NONREP_LOC2GLOB(gloc, rt_acrossMe, rt_meshCopyCount)


     ! Set the various terms that go into the diffusion calculation
     ! and compute the emission term:

     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

     do k = blkLimits(LOW,KAXIS)-ng*K3D, blkLimits(HIGH,KAXIS)+ng*K3D
        do j = blkLimits(LOW,JAXIS)-ng*K2D, blkLimits(HIGH,JAXIS)+ng*K2D
           do i = blkLimits(LOW,IAXIS)-ng, blkLimits(HIGH,IAXIS)+ng

              ! Get the opacities:

              absorb_opac = 0.0
              emit_opac = 0.0
              trans_opac = 0.0

              call Opacity(solnData(:,i,j,k), g, absorb_opac, emit_opac, trans_opac)

              ! Set the diffusion coefficient before application of flux limiter:
              solnData(COND_VAR, i, j, k) = rt_speedlt3/max(TINY(trans_opac)*rt_speedlt3,trans_opac)

              ! Set the value of the flux limit:
              solnData(ifl,i,j,k) = rt_mgdFlCoef * C * solnData(gvar,i,j,k) * solnData(DENS_VAR,i,j,k)
              solnData(ifl,i,j,k) = max(0.0, solnData(ifl,i,j,k))
           enddo
        enddo
     enddo

     ! Compute the flux limiter: returns 3*lambda (flux limiter factor) in iflOut, which may be the same variable as ifl

     call Diffuse_computeFluxLimiter(COND_VAR, gvar,DENS_VAR,ifl, iflOut, ieddi3, &
                                     rt_mgdFlMode, solnData,1,1,1,blockID,nlayers)

  enddo

end subroutine RadTrans_computeFluxLimiter
