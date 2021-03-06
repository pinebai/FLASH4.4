!!****if* source/Driver/DriverMain/Driver_finalizeFlash
!!
!! NAME
!!  Driver_finalizeFlash
!!
!! SYNOPSIS
!!  Driver_finalizeFlash()
!!
!! DESCRIPTION
!!
!!  Calls all the unit finalize routines
!!  which may need
!!  memory deallocated etc before the run end.
!!  Order does matter.
!!
!!***


subroutine Driver_finalizeFlash()

  use Eos_interface, ONLY : Eos_finalize
  use Driver_interface, ONLY : Driver_finalizeSourceTerms
  use NSE_interface, ONLY: NSE_finalize
  use PlasmaState_interface, ONLY: PlasmaState_finalize
  use RuntimeParameters_interface, ONLY : RuntimeParameters_finalize
  use Multispecies_interface, ONLY : Multispecies_finalize
  use Particles_interface, ONLY : Particles_finalize
  use Grid_interface, ONLY : Grid_finalize
  use Hydro_interface, ONLY : Hydro_finalize
  use Driver_data, ONLY: dr_globalMe, dr_restart
  use Simulation_interface, ONLY : Simulation_finalize
  use ProtonImaging_interface, ONLY : ProtonImaging_finalize
  use ProtonEmission_interface, ONLY : ProtonEmission_finalize
  use IO_interface, ONLY : IO_finalize
  use Cosmology_interface, ONLY: Cosmology_finalize
  use Timers_interface, ONLY: Timers_finalize
  use Gravity_interface, ONLY: Gravity_finalize
  use IncompNS_interface, ONLY: IncompNS_finalize
  use Imbound_interface, ONLY: Imbound_finalize
  use SolidMechanics_interface, ONLY: SolidMechanics_finalize
implicit none
#include "mpif.h"

  integer :: ierr
 
  
!!$  call Profiler_finalize()
!!$  
  call RuntimeParameters_finalize()

  call Multispecies_finalize()

  call Driver_finalizeSourceTerms( dr_restart) ! some of these exist only as stubs

  call NSE_finalize()             ! NSE material property

  call PlasmaState_finalize()     ! Plasma state utility unit

  call Grid_finalize()            ! Grid package
 
  call Particles_finalize()       ! Particles
  
  call Hydro_finalize()           ! Hydrodynamics
  
  call Eos_finalize()             ! Equation of State

  call ProtonImaging_finalize()   ! Proton Imaging

  call ProtonEmission_finalize()  ! Proton Emission

  call Cosmology_finalize()       ! Cosmology

  call Gravity_finalize()         ! Gravity

  call IncompNS_finalize()        ! INS 

  call SolidMechanics_finalize()  ! Solid Mechanics

  call ImBound_finalize()         ! Im Boundaries

  call IO_finalize()

  call Simulation_finalize()

  call Timers_finalize()

  call MPI_Finalize(ierr)

  return
end subroutine Driver_finalizeFlash








