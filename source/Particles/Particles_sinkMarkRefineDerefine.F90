!!****f* source/Particles/Particles_sinkMarkRefineDerefine
!!
!! NAME
!!
!!  Particles_sinkMarkRefineDerefine
!!
!! SYNOPSIS
!!
!!  call Particles_sinkMarkRefineDerefine()
!!
!! DESCRIPTION
!!
!!  This routine takes care of grid refinement based on Jeans length and sink particles.
!!  If the local density exceeds a given value that is computed based on Jeans analysis
!!  the block containing that cell is marked for refinement. Refinement and derefinement
!!  are triggered based on the number of cells per Jeans length, which the user must
!!  supply (jeans_ncells_ref and jeans_ncells_deref). Good values for these parameters are
!!  jeans_ncells_ref = 32 and jeans_ncells_deref = 64 (Federrath et al. 2011, ApJ 731, 62),
!!  but the user can choose any real number, where jeans_ncells_ref <= 2*jeans_ncells_deref.
!!  If sink particles are present, they must be at the highest level of AMR, so this routine
!!  also flags all cells within the sink particle accretion radius for refinement to the
!!  highest level.
!!
!! ARGUMENTS
!!
!! NOTES
!!
!!   written by Robi Banerjee, 2007-2008
!!   modified by Christoph Federrath, 2008-2015
!!   ported to FLASH3.3/4 by Chalence Safranek-Shrader, 2010-2012
!!   modified by Nathan Goldbaum, 2012
!!   refactored for FLASH4 by John Bachan, 2012
!!   cleaned and improved by Christoph Federrath, 2013-2015
!!
!!***

subroutine Particles_sinkMarkRefineDerefine()
  implicit none
end subroutine Particles_sinkMarkRefineDerefine
