!!****if* source/Particles/ParticlesMain/active/Sink/pt_sinkGetSubCycleTimeStep
!!
!! NAME
!!
!!  pt_sinkGetSubCycleTimeStep
!!
!! SYNOPSIS
!!
!!  call pt_sinkGetSubCycleTimeStep(real(out) :: dt,
!!                                  real(in)  :: dt_global,
!!                                  real(in)  :: local_min_radius,
!!                                  real(in)  :: local_max_accel)
!!
!! DESCRIPTION
!!
!!  Computes sink particle time step for subcycling in case of close encounters and
!!  highly eccentric orbits based on the minimum distance and the maximum acceleration
!!  between sinks.
!!
!! ARGUMENTS
!!
!!   dt - time step to be returned by this routine
!!
!!   dt_global - the current global time step
!!
!!   local_min_radius - minimum distance between sinks
!!
!!   local_max_accel - maximum acceleration between sinks
!!
!! NOTES
!!
!!   written by Christoph Federrath, 2008-2012
!!   ported to FLASH3.3/4 by Chalence Safranek-Shrader, 2010-2012
!!   modified by Nathan Goldbaum, 2012
!!   refactored for FLASH4 by John Bachan, 2012
!!
!!***

subroutine pt_sinkGetSubCycleTimeStep(dt, dt_global, local_min_radius, local_max_accel)

    use Particles_sinkData
    use RuntimeParameters_interface, ONLY : RuntimeParameters_get
    use Driver_data, ONLY : dr_globalMe
    
    implicit none

#include "constants.h"
#include "Flash.h"
#include "Particles.h"
    include "Flash_mpi.h"
    
    real, intent(out)  :: dt
    real, intent(in)   :: dt_global, local_min_radius, local_max_accel
    real, save         :: dtmin, sink_subdt_factor
    integer            :: i, ierr
    real               :: local_dt, max_part_vel
    logical, save      :: first_call = .true.
    logical, parameter :: debug = .false.

    if (first_call)  then
       call RuntimeParameters_get("dtmin", dtmin)
       call RuntimeParameters_get("sink_subdt_factor", sink_subdt_factor)
       first_call = .false.
    end if

    local_dt = dt_global

    ! time step based on grav accel
    if (local_min_radius .gt. 0 .and. local_max_accel .gt. 0) then
       local_dt = min(local_dt, sink_subdt_factor*sqrt(local_min_radius/local_max_accel))
    end if

    ! get maximum particle velocity
    max_part_vel = 0.0

    if (sink_AdvanceSerialComputation) then
      do i = 1, localnpf
        max_part_vel = max(max_part_vel, abs(particles_global(VELX_PART_PROP, i)))
        max_part_vel = max(max_part_vel, abs(particles_global(VELY_PART_PROP, i)))
        max_part_vel = max(max_part_vel, abs(particles_global(VELZ_PART_PROP, i)))
      end do
    else
      do i = 1, localnp
        max_part_vel = max(max_part_vel, abs(particles_local(VELX_PART_PROP, i)))
        max_part_vel = max(max_part_vel, abs(particles_local(VELY_PART_PROP, i)))
        max_part_vel = max(max_part_vel, abs(particles_local(VELZ_PART_PROP, i)))
      end do
    endif

    if(debug) then
       print*, dr_globalMe, "get subcycle time: localnp=", localnp
       print*, dr_globalMe, "get subcycle time: max_part_vel=", max_part_vel
       print*, dr_globalMe, "get subcycle time: local_min_radius=", local_min_radius
    endif

    ! time step based on particle velocities
    if(local_min_radius .gt. 0 .and. max_part_vel .gt. 0) then
       local_dt = min(local_dt, sink_subdt_factor*local_min_radius/max_part_vel)
    end if

    ! no zero time step...
    if (local_dt .le. 0 ) then
       local_dt = dtmin
       print*, "WARNING... zero sink particle timestep!"
    end if

    ! avoid MPI_ALLREDUCE during subcycling if sink_AdvanceSerialComputation = .true.
    if (sink_AdvanceSerialComputation) then
      dt = local_dt ! is already the global min dt, in this case
    else
      ! get smallest timestep from all processors
      call MPI_ALLREDUCE(local_dt, dt, 1, FLASH_REAL, MPI_MIN, MPI_COMM_WORLD, ierr)
    endif

    if (debug) then
       if(dr_globalMe .eq. MASTER_PE) then
          print*, "get sub cycycle time - timestep=", local_dt
       endif
    end if

    return

end subroutine pt_sinkGetSubCycleTimeStep
