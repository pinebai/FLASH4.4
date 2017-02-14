!!****if* source/physics/Gravity/GravityMain/Poisson/BHTree/Gravity_bhNodeContrib
!!
!! NAME
!!
!!  Gravity_bhNodeContrib
!!
!!
!! SYNOPSIS
!!
!!   call Gravity_bhNodeContrib(
!!                          real(in)       :: node(:),
!!                          integer(in)    :: trLevel,
!!                          integer(in)    :: refLevel,
!!                          real(in)       :: dr(MDIM+2),
!!                          integer(in)    :: blockno,
!!                          integer(in)    :: point(MDIM),
!!                          integer(in)    :: blkLimits(2,MDIM),
!!                          real,pointer   :: solnData(:,:,:,:)
!!        )
!!
!! DESCRIPTION
!!
!!  Calculates contribution of the node to the gravitational
!!  potential and adds it to the solnData array (with index grv_defaultGpotVar).
!!
!! ARGUMENTS
!!
!!  node        : array of the node of the tree, whose
!!                contribution is added to the gravitational potential
!!  trLevel     : level of the node within the block-tree
!!  refLevel    : refinement level of the block with the node
!!  dr          : (1:MDIM) - position vector from the point-of-calculation to the node
!!                (MDIM+1) - square of the magnitude of the position vector
!!                (MDIM+2) - inverted magnitude of the position vector
!!  blockno     : number of block into which the point-of-calculation belongs
!!  point       : indeces of the point-of-calculation in the block
!!  blkLimits   : limits of indeces in the block
!!  solnData    : solution data from the grid
!!
!!
!!***

subroutine Gravity_bhNodeContrib(node, trLevel, refLevel, dr, blockno, point, blkLimits, solnData)
  use Gravity_data, ONLY : useGravity, grv_bhNewton, grv_bhUseEwald, &
    grv_defaultGpotVar, grv_bhIM, grv_bhIBM, grv_bhTreeLevels
  use Gravity_interface, ONLY : Gravity_bhEwaldAccV42, Gravity_bhEwaldPotV42, &
    Gravity_bhEwaldAcc, Gravity_bhEwaldPot
  implicit none
#include "constants.h"
#include "Flash.h"
#include "FortranLangFeatures.fh"
  real, dimension(:), intent(IN) :: node
  integer, intent(IN) :: trLevel, refLevel, blockno
  integer, dimension(MDIM), intent(IN) :: point
  real, dimension(MDIM+2), intent(IN) :: dr
  integer, dimension(2,MDIM), intent(IN)   :: blkLimits
  real, DIMENSION(:,:,:,:), POINTER_INTENT_IN :: solnData
  real :: ndMass
#ifdef GRAV_TREE_ACC
  real, dimension(MDIM) :: gacc
#endif

  if (.not. useGravity) return

  if (trLevel == grv_bhTreeLevels) then
    ndMass = node(grv_bhIBM)
  else
    ndMass = node(grv_bhIM)
  endif

  ! add the contribution to the potential
  if (grv_bhUseEwald) then

    !print*,'Calling Ewald field, dr = ',dr

#ifdef GRAV_TREE_EWALD_V42
    ! calculate gravitational acceleration using the original Ewald field
#ifdef GRAV_TREE_ACC
    gacc = grv_bhNewton * ndMass &
    &    * sign(Gravity_bhEwaldAccV42(abs(dr(IAXIS)), abs(dr(JAXIS)), abs(dr(KAXIS))),dr(1:MDIM))
    solnData(GACX_VAR, point(IAXIS), point(JAXIS), point(KAXIS)) &
    & = solnData(GACX_VAR, point(IAXIS), point(JAXIS), point(KAXIS)) &
    & + gacc(IAXIS)
    solnData(GACY_VAR, point(IAXIS), point(JAXIS), point(KAXIS)) &
    & = solnData(GACY_VAR, point(IAXIS), point(JAXIS), point(KAXIS)) &
    & + gacc(JAXIS)
    solnData(GACZ_VAR, point(IAXIS), point(JAXIS), point(KAXIS)) &
    & = solnData(GACZ_VAR, point(IAXIS), point(JAXIS), point(KAXIS)) &
    & + gacc(KAXIS)
#endif
    ! potential is always calculated
    solnData(grv_defaultGpotVar, point(IAXIS), point(JAXIS), point(KAXIS)) &
    & = solnData(grv_defaultGpotVar, point(IAXIS), point(JAXIS), point(KAXIS)) &
    & - grv_bhNewton*ndMass * Gravity_bhEwaldPotV42(abs(dr(IAXIS)), &
    &   abs(dr(JAXIS)), abs(dr(KAXIS)))
#else
    ! use Taylor expansion of the Ewald field
#ifdef GRAV_TREE_ACC
    gacc = grv_bhNewton * ndMass &
    &    * sign(Gravity_bhEwaldAcc(abs(dr(IAXIS)), abs(dr(JAXIS)), abs(dr(KAXIS)), dr(MDIM+2)),dr(1:MDIM))
    solnData(GACX_VAR, point(IAXIS), point(JAXIS), point(KAXIS)) &
    & = solnData(GACX_VAR, point(IAXIS), point(JAXIS), point(KAXIS)) &
    & + gacc(IAXIS)
    solnData(GACY_VAR, point(IAXIS), point(JAXIS), point(KAXIS)) &
    & = solnData(GACY_VAR, point(IAXIS), point(JAXIS), point(KAXIS)) &
    & + gacc(JAXIS)
    solnData(GACZ_VAR, point(IAXIS), point(JAXIS), point(KAXIS)) &
    & = solnData(GACZ_VAR, point(IAXIS), point(JAXIS), point(KAXIS)) &
    & + gacc(KAXIS)
#endif
    ! potential is always calculated because sinks need it
    ! calculate gravitational potential
    ! it is calculated even with GRAV_TREE_ACC is true, because
    ! sink particles need potential
    solnData(grv_defaultGpotVar, point(IAXIS), point(JAXIS), point(KAXIS)) &
    & = solnData(grv_defaultGpotVar, point(IAXIS), point(JAXIS), point(KAXIS)) &
    & - grv_bhNewton*ndMass * Gravity_bhEwaldPot(abs(dr(IAXIS)), &
    &   abs(dr(JAXIS)), abs(dr(KAXIS)), dr(MDIM+2))

#endif
  else
    ! calculate gravitational acceleration directly
#ifdef GRAV_TREE_ACC
    gacc = grv_bhNewton*ndMass*dr(MDIM+2)*dr(MDIM+2)*dr(MDIM+2) &
    &    * dr(1:MDIM)
    solnData(GACX_VAR, point(IAXIS), point(JAXIS), point(KAXIS)) &
    & = solnData(GACX_VAR, point(IAXIS), point(JAXIS), point(KAXIS)) &
    & + gacc(IAXIS)
    solnData(GACY_VAR, point(IAXIS), point(JAXIS), point(KAXIS)) &
    & = solnData(GACY_VAR, point(IAXIS), point(JAXIS), point(KAXIS)) &
    & + gacc(JAXIS)
    solnData(GACZ_VAR, point(IAXIS), point(JAXIS), point(KAXIS)) &
    & = solnData(GACZ_VAR, point(IAXIS), point(JAXIS), point(KAXIS)) &
    & + gacc(KAXIS)
    !print *, "NC: ", gacc, solnData(GACX_VAR:GACZ_VAR, point(IAXIS), point(JAXIS), point(KAXIS))
#endif
    ! potential is always calculated
    solnData(grv_defaultGpotVar, point(IAXIS), point(JAXIS), point(KAXIS)) &
    & = solnData(grv_defaultGpotVar, point(IAXIS), point(JAXIS), point(KAXIS)) &
    & - grv_bhNewton*ndMass*dr(MDIM+2)
  endif

  return
end subroutine Gravity_bhNodeContrib
