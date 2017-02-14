!!****if* source/Simulation/SimulationMain/unitTest/SinkMomTest/Simulation_init
!!
!! NAME
!!
!!  Simulation_init
!!
!! SYNOPSIS
!!
!!  call Simulation_init()
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   No arguments
!!
!! AUTOGENROBODOC
!!
!!
!!***


subroutine Simulation_init()

  use Simulation_data
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Particles_sinkData
  use pt_sinkInterface, only: pt_sinkCreateParticle, pt_sinkGatherGlobal
  use Driver_interface, ONLY : Driver_getMype

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Particles.h"

  integer :: myPE, pno, blockID
  real    :: pt
  logical :: restart

  logical, parameter :: Debug = .false.


  ! do nothing on restart
  call RuntimeParameters_get("restart", restart)
  if (restart) return


  if (Debug) print *, 'Simulation_init() entering.'

  call Driver_getMype(GLOBAL_COMM, myPE)
  sim_globalMe = myPE

  call RuntimeParameters_get('xmin',sim_xMin)
  call RuntimeParameters_get('ymin',sim_yMin)
  call RuntimeParameters_get('zmin',sim_zMin)
  call RuntimeParameters_get('xmax',sim_xMax)
  call RuntimeParameters_get('ymax',sim_yMax)
  call RuntimeParameters_get('zmax',sim_zMax)
  
  ! place initial sink particle

  if (sim_globalMe == MASTER_PE) then

    call RuntimeParameters_get("sim_sink_x", sim_sink_x)
    call RuntimeParameters_get("sim_sink_y", sim_sink_y)
    call RuntimeParameters_get("sim_sink_z", sim_sink_z)
    call RuntimeParameters_get("sim_sink_vx", sim_sink_vx)
    call RuntimeParameters_get("sim_sink_vy", sim_sink_vy)
    call RuntimeParameters_get("sim_sink_vz", sim_sink_vz)
    call RuntimeParameters_get("sim_sink_mass", sim_sink_mass)

    blockID = 1
    pt = 0.0

    pno = pt_sinkCreateParticle(sim_sink_x, sim_sink_y, sim_sink_z, pt, blockID, sim_globalMe)

    particles_local(VELX_PART_PROP, 1) = sim_sink_vx
    particles_local(VELY_PART_PROP, 1) = sim_sink_vy
    particles_local(VELZ_PART_PROP, 1) = sim_sink_vz
    particles_local(MASS_PART_PROP, 1) = sim_sink_mass

    write(*,'(A,4(1X,ES16.9),3I8)') "initial sink particle created (x, y, z, pt, blockID, MyPE, tag): ", &
      & sim_sink_x, sim_sink_y, sim_sink_z, pt, blockID, sim_globalMe, int(particles_local(TAG_PART_PROP,pno))

  endif

  call RuntimeParameters_get("sim_massTol", sim_massTol)
  call RuntimeParameters_get("sim_momXTol", sim_momXTol)
  call RuntimeParameters_get("sim_momYTol", sim_momYTol)
  call RuntimeParameters_get("sim_momZTol", sim_momZTol)


  call pt_sinkGatherGlobal()

  sim_testInitialized = .FALSE.

  if (Debug) print *, 'Simulation_init() done.'

end subroutine Simulation_init
