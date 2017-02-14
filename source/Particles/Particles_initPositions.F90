!!****f* source/Particles/Particles_initPositions
!!
!! NAME
!!    Particles_initPositions
!!
!! SYNOPSIS
!!
!!    Particles_initPositions( logical, INTENT(inout) :: partPosInitialized,
!!                             logical, INTENT(out) :: updateRefine)
!!
!! DESCRIPTION
!!    Initialize particle locations.  There are the versions currently
!!    available. 
!!
!!      Lattice creates particles which are evenly distributed in space
!!      WithDensity creates particles distributed according to the global
!!         distribution of density.
!!      FromFile reads them in from a file
!!
!! ARGUMENTS
!!
!!
!!  partPosInitialized : boolean. Upon entry it is true if all particles have been
!!            successfully intialized, and there is no need for further
!!            action. Upon exit it has is true if it successfully initialized
!!            all the particles locally.
!!            The routine for configuring the initial AMR grid
!!            "gr_expandDomain" determines if all particles initialization
!!            was globally successful.
!!            When AMR is not being used, or when refinement based upon
!!            Particles is turned off, this routine is called directly 
!!            from Driver_initFlash with partPosInitialized set to .false. so that
!!            particles are initialized. In such instances, the value
!!            at output has no meaning.          
!! updateRefine : is set to true if the implementation wishes to
!!                retain the already initialized particles;
!!                is set to false if the implementation wants particles
!!                to be reinitialized instead as the grid refines.
!!
!! PARAMETERS
!!
!!    Each initialization version has parameters.  See the Config files
!!
!!***


subroutine Particles_initPositions (  partPosInitialized,updateRefine)


  implicit none


  logical, INTENT(INOUT) :: partPosInitialized
  logical, INTENT(OUT) :: updateRefine
  partPosInitialized=.true.

  updateRefine=.false.

end subroutine Particles_initPositions


