!!****if* source/Particles/ParticlesMain/active/Sink/pt_sinkFindList
!!
!! NAME
!!
!!  pt_sinkFindList
!!
!! SYNOPSIS
!!
!!  call pt_sinkFindList(real, intent(IN) :: x,
!!                       real, intent(IN) :: y,
!!                       real, intent(IN) :: z,
!!                       real, intent(IN) :: rad,
!!                       logical, intent(IN) :: create_part,
!!                       integer, dimension(maxsinks), intent(OUT) :: pindex_found,
!!                       integer, intent(OUT) :: np_found)
!!
!! DESCRIPTION
!!
!!  Searches the global and local sink particle lists to find any sink particle within
!!  a given search radius around any given position (x,y,z) and returns the list indexes
!!  of the found particels in pindex_found and the number found in np_found.
!!  If called with create_part = .false. then only particles_global(ipx, ipy, ipz, iptag)
!!  must be up-to-date, while calling it with create_part = .false. requires that a full
!!  pt_sinkGatherGlobal() was called before, such that the whole particles_global list is
!!  up-to-date.
!!
!! ARGUMENTS
!!
!!   x - x center of search sphere
!!
!!   y - y center of search sphere
!!
!!   z - z center of search sphere
!!
!!   rad - radius of search sphere
!!
!!   create_part - logical flag indicating whether to create a local dummy
!!                 particle or not, in case the found particle did not
!!                 reside on the local processor and was thus only found in
!!                 the global sink particle list.
!!
!!   pindex_found - list of found particle indexes to be returned by this routine
!!
!!   np_found - number of found particles
!!
!! NOTES
!!
!!   written by Robi Banerjee, 2007-2008
!!   modified by Christoph Federrath, 2008-2012
!!   ported to FLASH3.3/4 by Chalence Safranek-Shrader, 2010-2012
!!   modified by Martin Schroen, 2011
!!   modified by Nathan Goldbaum, 2012
!!   refactored for FLASH4 by John Bachan, 2012
!!
!!***

subroutine pt_sinkFindList(x, y, z, rad, create_part, pindex_found, np_found)

  use Particles_sinkData
  use Cosmology_interface, ONLY : Cosmology_getRedshift
  use Driver_interface, ONLY : Driver_abortFlash
  use pt_sinkInterface, ONLY: pt_sinkCorrectForPeriodicBCs
  use RuntimeParameters_interface, ONLY: RuntimeParameters_get

  implicit none

#include "Flash.h"
#include "Particles.h"
#include "constants.h"

  real, intent(IN)    :: x, y, z, rad
  logical, intent(IN) :: create_part

  integer, dimension(maxsinks), intent(OUT) :: pindex_found
  integer, intent(OUT)                      :: np_found
  integer                 :: found, found2, pno, lp
  real                    :: dist, redshift, onePlusRedshift, dx, dy, dz

  logical, save :: first_call = .true.

  character(len=80), save :: grav_boundary_type

  if (first_call) then
     call RuntimeParameters_get("grav_boundary_type", grav_boundary_type)
     first_call = .false.
  end if

  call Cosmology_getRedshift(redshift)
  onePlusRedshift = redshift + 1.0

  np_found = 0

  do pno = 1, localnpf

     found = 0

     dx = x - particles_global(ipx,pno)
     dy = y - particles_global(ipy,pno)
     dz = z - particles_global(ipz,pno)

     if (grav_boundary_type .eq. "periodic") call pt_sinkCorrectForPeriodicBCs(dx, dy, dz)

     dist = sqrt(dx**2 + dy**2 + dz**2)

     if (dist .le. rad) found = pno

     if (found .gt. 0) then

        np_found = np_found + 1

        found2 = 0
        do lp = 1, localnp
          if (int(particles_local(iptag,lp)) .eq. int(particles_global(iptag,pno))) found2 = lp
        enddo

        if (found2 .gt. 0) then
          pindex_found(np_found) = found2
        else if (create_part) then
          ! it was a particle in the global list, so create a dummy particle if desired
          localnp = localnp + 1
          if (localnp .gt. sink_maxSinks) &
            call Driver_abortFlash('sink_findList: Sink particle number exceeds sink_maxSinks. Increase.')
          particles_local(:,localnp) = particles_global(:,found)
          particles_local(ipblk,localnp) = NONEXISTENT
          pindex_found(np_found) = localnp
        else
          ! do this even if found <= localnp (means, particle is already created)
          pindex_found(np_found) = found ! index in global list
        endif

     endif ! found .gt. 0

  end do ! loop over all particles in global list

  return

end subroutine pt_sinkFindList
