!!****f* source/physics/sourceTerms/Heatexchange/Heatexchange_computeDt
!!
!! NAME
!!
!!  Heatexchange_computeDt
!!
!! SYNOPSIS
!!
!!  Heatexchange_computeDt(integer(IN) :: blockID,
!!                 )
!!                 integer(IN) :: blkLimits(2,MDIM)
!!                 integer(IN) :: blkLimitsGC(2,MDIM)
!!                 real,pointer::  solnData(:,:,:,:),   
!!                 real(OUT)   :: dt_heatXchg, 
!!                 real(OUT)   :: dt_minloc(5)) 
!!
!!
!! DESCRIPTION
!!
!!  compute a burning timestep limiter, by trying to force the energy
!!  generation from burning to be smaller than the internal energy
!!  in a zone.
!!
!!   The timestep limiter would be:
!!
!!                                      eintN
!!             dt     =  hx_dtFactor * ------
!!               burn                   qNdot
!!
!!  qNdot is energy/mass/s, so the time factor is already in there, and we
!!  are doing
!!
!!             
!!                                      eintN
!!             dt     =  hx_dtFactor * ------   * dt
!!               burn                   qNdot*dt
!!
!!  qNdot*dt is the amount of energy / mass deposited in a zone by burning. 
!!  eintN is the internal energy / mass in that zone.  If qNdot*dt is 2x
!!  eintN, then we want a timestep that is half the size.  
!!
!!  hx_dtFactor is a prefactor to scaling the timestep.  In general, we aim
!!  for qNdot*dt < hx_dtFactor * eintN.  For good coupling between the hydro
!!  and the burner, hx_dtFactor should be < 1.
!!
!!
!! ARGUMENTS
!!
!!  blockID       --  local block ID
!!  
!!  blkLimits     --  the indices for the interior endpoints of the block
!!  blkLimitsGC   --  the indices for endpoints including the guardcells
!!  solnData      --  the physical, solution data from grid
!!  dt_heatXchg       --  variable to hold timestep constraint
!!  dt_minloc(5)  --  array to hold limiting zone info:  zone indices
!!                    (i,j,k), block ID, PE number
!!
!!
!! PARAMETERS
!!
!!  hx_dtFactor    A parameter, such that qNdot*dt < hx_dtFactor * eint,
!!                  that is, the energy release from burning divided by
!!                  the internal energy in that zone is < hx_dtFactor.
!!
!! SEE ALSO
!!
!!  Driver_computeDt
!!
!!
!!***

subroutine Heatexchange_computeDt(blockID,   &
                           blkLimits,blkLimitsGC,        &
                           solnData,   &
                           dt_heatXchg, dt_minloc)

  implicit none

#include "constants.h"
#include "Flash.h"

  !! arguments
  integer, intent(IN)   :: blockID 
  integer, intent(IN),dimension(2,MDIM)::blkLimits,blkLimitsGC
  real, pointer           :: solnData(:,:,:,:) 
  real, intent(INOUT)     :: dt_heatXchg
  integer, intent(INOUT)  :: dt_minloc(5)

  return

end subroutine Heatexchange_computeDt
