!!****if* source/Particles/ParticlesMain/active/Sink/Particles_sinkComputeDt
!!
!! NAME
!!
!!  Particles_sinkComputeDt
!!
!! SYNOPSIS
!!
!!  call Particles_sinkComputeDt(integer, INTENT(in) :: blockid,
!!                               real, INTENT(inout) :: dt_sink,
!!                            integer, INTENT(inout) :: dt_minloc)
!!
!! DESCRIPTION
!!
!!  Constrains the global timestep based on sink particle velocities and
!!  gravitational accelerations.
!!
!! ARGUMENTS
!!
!!   blockid - ID of block in current processor
!!
!!   dt_sink - dt constrained by sinks
!!
!!   dt_minloc - location of the cell that constrains the sink time step
!!
!! NOTES
!!
!!   written by Christoph Federrath, 2008-2015
!!   ported to FLASH3.3/4 by Chalence Safranek-Shrader, 2010-2012
!!   modified by Nathan Goldbaum, 2012
!!   refactored for FLASH4 by John Bachan, 2012
!!   added acceleration timestep constraint; Christoph Federrath, 2013
!!   added sinks off-domain support; Christoph Federrath, 2015
!!
!!***

subroutine Particles_sinkComputeDt(blockID,dt_sink,dt_minloc)

   ! Timestep correction for sink particles, idential to that
   ! of other particles. Sink particles can not travel more than
   ! a fraction sink_dt_factor across a cell in a single
   ! timestep.

   use Particles_sinkData, ONLY : localnp, particles_local, &
                  & sink_offDomainSupport, localnpf, particles_global
   use Grid_interface, ONLY : Grid_getBlkBoundBox, Grid_getDeltas
   use Driver_data, ONLY : dr_globalMe
   use RuntimeParameters_interface, ONLY : RuntimeParameters_get

   implicit none

#include "Flash.h"
#include "constants.h"

   integer, INTENT(in)    :: blockID
   real, INTENT(inout)    :: dt_sink
   integer, INTENT(inout) :: dt_minloc(5)
   logical, save          :: first_call = .true.
   real                   :: dtx, dty, dtz, dtnew, velxabs, velyabs, velzabs, accxabs, accyabs, acczabs
   real                   :: delta(MDIM), boundBox(2,MDIM)
   real, save             :: sink_dt_factor
   integer                :: i
   real, parameter        :: tiny_number = tiny(real(1.0))

   if (first_call) then
      call RuntimeParameters_get("sink_dt_factor", sink_dt_factor)
      first_call =.false.
   end if

   call Grid_getBlkBoundBox(blockID,boundBox) ! physical bounding box of the block
   call Grid_getDeltas(blockID,delta) ! cell size dx, dy, dz of the block

   if (sink_offDomainSupport) then

     do i = 1, localnpf

       if ( (particles_global(POSX_PART_PROP,i).ge.boundBox(1,1)-2.0*delta(1)) .and. &
            (particles_global(POSX_PART_PROP,i).le.boundBox(2,1)+2.0*delta(1)) .and. &
            (particles_global(POSY_PART_PROP,i).ge.boundBox(1,2)-2.0*delta(2)) .and. &
            (particles_global(POSY_PART_PROP,i).le.boundBox(2,2)+2.0*delta(2)) .and. &
            (particles_global(POSZ_PART_PROP,i).ge.boundBox(1,3)-2.0*delta(3)) .and. &
            (particles_global(POSZ_PART_PROP,i).le.boundBox(2,3)+2.0*delta(3))  ) then

         velxabs = abs(particles_global(VELX_PART_PROP,i))+tiny_number*delta(1)
         accxabs = abs(particles_global(ACCX_PART_PROP,i))+tiny_number*delta(1)
         dtx = min(delta(1)/velxabs,sqrt(delta(1)/accxabs))

         velyabs = abs(particles_global(VELY_PART_PROP,i))+tiny_number*delta(2)
         accyabs = abs(particles_global(ACCY_PART_PROP,i))+tiny_number*delta(2)
         dty = min(delta(2)/velyabs,sqrt(delta(2)/accyabs))

         velzabs = abs(particles_global(VELZ_PART_PROP,i))+tiny_number*delta(3)
         acczabs = abs(particles_global(ACCZ_PART_PROP,i))+tiny_number*delta(3)
         dtz = min(delta(3)/velzabs,sqrt(delta(3)/acczabs))

         dtnew = sink_dt_factor * min(dtx, dty, dtz)
         if (dtnew .lt. dt_sink) then
            dt_sink = dtnew
            ! info about where tstep restriction took place
            dt_minloc(1) = int((particles_global(POSX_PART_PROP,i)-boundBox(1,1))/delta(1)) + NGUARD+1
            dt_minloc(2) = int((particles_global(POSY_PART_PROP,i)-boundBox(1,2))/delta(2)) + NGUARD+1
            dt_minloc(3) = int((particles_global(POSZ_PART_PROP,i)-boundBox(1,3))/delta(3)) + NGUARD+1
            dt_minloc(4) = blockID
            dt_minloc(5) = dr_globalMe
         end if

       end if ! particle inside or near block ?

     end do ! loop over global particles

   else ! normal sinks (no off-domain support)

     do i = 1, localnp

       if (int(particles_local(BLK_PART_PROP,i)) .eq. blockID) then

         velxabs = abs(particles_local(VELX_PART_PROP,i))+tiny_number*delta(1)
         accxabs = abs(particles_local(ACCX_PART_PROP,i))+tiny_number*delta(1)
         dtx = min(delta(1)/velxabs,sqrt(delta(1)/accxabs))

         velyabs = abs(particles_local(VELY_PART_PROP,i))+tiny_number*delta(2)
         accyabs = abs(particles_local(ACCY_PART_PROP,i))+tiny_number*delta(2)
         dty = min(delta(2)/velyabs,sqrt(delta(2)/accyabs))

         velzabs = abs(particles_local(VELZ_PART_PROP,i))+tiny_number*delta(3)
         acczabs = abs(particles_local(ACCZ_PART_PROP,i))+tiny_number*delta(3)
         dtz = min(delta(3)/velzabs,sqrt(delta(3)/acczabs))

         dtnew = sink_dt_factor * min(dtx, dty, dtz)
         if (dtnew .lt. dt_sink) then
            dt_sink = dtnew
            ! info about where tstep restriction took place
            dt_minloc(1) = int((particles_local(POSX_PART_PROP,i)-boundBox(1,1))/delta(1)) + NGUARD+1
            dt_minloc(2) = int((particles_local(POSY_PART_PROP,i)-boundBox(1,2))/delta(2)) + NGUARD+1
            dt_minloc(3) = int((particles_local(POSZ_PART_PROP,i)-boundBox(1,3))/delta(3)) + NGUARD+1
            dt_minloc(4) = blockID
            dt_minloc(5) = dr_globalMe
         end if

       end if ! particle in block ?

     end do ! loop over local particles

   end if ! off-domain support ?

   return

end subroutine Particles_sinkComputeDt
