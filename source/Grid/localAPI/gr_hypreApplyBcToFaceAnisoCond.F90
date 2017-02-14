!!****if* source/Grid/localAPI/gr_hypreApplyBcToFaceAnisoCond
!!
!!  NAME
!!
!!  gr_hypreApplyBcToFaceAnisoCond
!!
!!  SYNOPSIS
!!
!!  call gr_hypreApplyBcToFaceAnisoCond (integer, intent(IN) :: blkLimitsGC(2,MDIM),
!!                              integer, intent(IN) :: blkLimits(2,MDIM),
!!                              integer, intent(IN) :: part,
!!                              integer, intent(IN) :: var,
!!                              integer, intent(IN) :: iFactorB,
!!                              integer, intent(IN) :: bcType,
!!                              integer, intent(IN) :: direction,
!!                              real,    intent(IN) :: bcValue(2),
!!                              real,    intent(IN) :: dt,
!!                              real,    intent(IN) :: theta,
!!                              real,    intent(IN) :: del,
!!                              integer, intent(IN) :: Lower(MDIM),
!!                              real,    intent(IN) :: scalefactor,
!!                              real,    intent(IN) :: faceArea(:,:,:),
!!                              real,    intent(IN) :: solnVec(:,:,:,:))
!!
!!  DESCRIPTION
!!      This routines modifies Matrix A (as in AX=B) or the right hand side B
!!      to account for the applied physical boundary conditions. The
!!      applied boundary conditions might or might not require changes
!!      to the RHS, particularly the way the diffusion operator has been
!!      implemented, OUTFLOW BC ends up as a no-op.
!!
!! ARGUMENTS
!!
!!   blkLimitsGC  : an array that holds the lower and upper indices of the
!!                  section of block with the guard cells.
!!   blkLimits    : an array that holds the lower and upper indices of the section
!!                  of block without the guard cells.
!!   part         : HYPRE part to which the block belongs to. Determined based
!!                  on its refinement level.
!!   var          : Variable associated with the HYPRE GRID (always set to zero).
!!   iFactorB     : Conductivity or Opacitiy.
!!   bcType       : GRID_PDE_BND_NEUMANN, VACUUM,
!!                  GRID_PDE_BND_DIRICHLET, GRID_PDE_BND_GIVENGRAD
!!   direction    : LOWER or UPPER face (i.e. ILO_FACE, IHI_FACE, JLO_FACE ...)
!!   bcValue      : used when GRID_PDE_BND_DIRICHLET/GRID_PDE_BND_GIVENGRAD is specified as bcType
!!                  GRID_PDE_BND_DIRICHLET
!!                    bcValue(1) -> Value of Temperature/Energy-Density specified on boundary.
!!                    bcValue(2) -> Value of Conductivity/Opacity specified on boundary.
!!                  GRID_PDE_BND_GIVENGRAD
!!                     bcValue(1) -> Value of gradient as specified on boundary.
!!                     bcValue(2) -> not used.
!!   dt           : Global time step.
!!   theta        : Varies scheme (0-> Explicit, 1-> backward euler, 0.5 -> Crank Nicholson
!!   del          : dx/dy/dz.
!!   Lower        : Lower limits of the domain (Global, computed based on
!!                  strides and block CornerID).
!!   scalefactor  : used to scale the RHS properly. Used with GIVENRAD, DIRICHLET only.
!!   faceArea     : Areas of the cells on the face (ILO_FACE, JLO_FACE ....).
!!                  At the same refinement level, cell values differ only for non-Cartesian
!!                  coordinates.
!!   solnVec      : slice of UNK corresponding to a particular block (blockID).
!!
!! SIDE EFFECTS
!!
!!
!! NOTES
!!
!!
!!***

!!REORDER(4): solnVec

subroutine gr_hypreApplyBcToFaceAnisoCond(blkLimits,blkLimitsGC,part,var,iFactorB,bcType,direction, &
     bcValue, dt, theta, del, Lower, scalefactor, faceArea, solnVec, blockID)
    
  implicit none

#include "Flash.h"  
#include "constants.h"
  
  integer, intent(IN) :: blkLimits (2,MDIM) 
  integer, intent(IN) :: blkLimitsGC (2,MDIM)
  integer, intent(IN) :: part
  integer, intent(IN) :: var
  integer, intent(IN) :: iFactorB
  integer, intent(IN) :: bcType
  integer, intent(IN) :: direction
  real,    intent(IN) :: bcValue(2)
  real,    intent(IN) :: dt
  real,    intent(IN) :: theta
  real,    intent(IN) :: del
  integer, intent(IN) :: Lower(MDIM), blockID 
  real,    intent(IN) :: scalefactor
  real,    intent(IN) :: faceArea(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
                                  blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
                                  blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS))
  real,    intent(IN) :: solnVec(NUNK_VARS, blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
       blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
       blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS))   
    

 
  return
  
end subroutine gr_hypreApplyBcToFaceAnisoCond

