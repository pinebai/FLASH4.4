!!****if* source/Particles/ParticlesMain/active/Sink/pt_sinkParticleMerging
!!
!! NAME
!!
!!  pt_sinkParticleMerging
!!
!! SYNOPSIS
!!
!!  call pt_sinkParticleMerging(real, intent(IN) :: dt)
!!
!! DESCRIPTION
!!
!!  Handles merging of sink particles, if the user switches it on.
!!  This subroutine is only marginally tested and might not work properly
!!  in all cases. Watch out for sink particle multiplication after merging.
!!  Eventually this needs to be done in a similar way as in SinkMergingAfterCreation(),
!!  which handles merging right after sink creation in symmetric setups. Sink particle
!!  merging is discouraged, but this routine provides the basic methods to do it, if the
!!  user really thinks they want to do it.
!!
!! ARGUMENTS
!!
!!   dt - the current simulation time step
!!
!! NOTES
!!
!!   written by Christoph Federrath, 2008-2012
!!   ported to FLASH3.3/4 by Chalence Safranek-Shrader, 2010-2012
!!   modified by Nathan Goldbaum, 2012
!!   refactored for FLASH4 by John Bachan, 2012
!!
!!***

subroutine pt_sinkParticleMerging(dt)

  use Particles_sinkData
  use pt_sinkInterface, only: pt_sinkGatherGlobal
  use Driver_data, ONLY : dr_globalMe
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get
  use Particles_interface, ONLY : Particles_sinkMoveParticles
  use Driver_interface, ONLY : Driver_abortFlash
  use Cosmology_interface, ONLY : Cosmology_getRedshift
  implicit none

#include "constants.h"
#include "Flash.h"
#include "GridParticles.h"
#include "Particles.h"
  include "Flash_mpi.h"
  
  real, intent(IN)    :: dt

  logical, save       :: first_call = .true.
  integer, save       :: MyPE, MasterPE
  real, save          :: accretion_radius, Newton

  integer             :: lp1, lp2, ptag_merged_to, ptag_merged_from, lpkeep, lpdelete, ierr
  real                :: rad, egrav, ekin, etot, dx, dy, dz
  real                :: accretion_radius_comoving, pm_del, px_old, py_old, pz_old, vx_del, vy_del, vz_del
  real                :: px1, px2, py1, py2, pz1, pz2, pm1, pm2, pvx1, pvx2, pvy1, pvy2, pvz1, pvz2
  real                :: pvcmx, pvcmy, pvcmz, vrad, new_mass, old_mass
  real                :: redshift, onePlusRedshift, onePlusRedshift2, onePlusRedshift3
  logical             :: sink_merging_happened, sink_merging_happened_on_masterPE

  character(len=80), save :: grav_boundary_type

  if (first_call) then

     MyPE = dr_globalMe
     MasterPE = MASTER_PE

     call RuntimeParameters_get("sink_accretion_radius", accretion_radius)
     call PhysicalConstants_get("Newton", Newton)
     call RuntimeParameters_get("grav_boundary_type", grav_boundary_type)

     first_call = .false.

  end if

  call Cosmology_getRedshift(redshift)
  onePlusRedshift = 1.0 + redshift
  onePlusRedshift2 = onePlusRedshift * onePlusRedshift
  onePlusRedshift3 = onePlusRedshift2 * onePlusRedshift

  accretion_radius_comoving = accretion_radius * onePlusRedshift

  call pt_sinkGatherGlobal()

  sink_merging_happened_on_masterPE = .false.
  sink_merging_happened = .false.

  ! only the master proc takes care of merging the particles
  ! (only work on the global particles list)
  if (MyPE .eq. MasterPE) then

     print*, "Starting to merge sink particles on MasterPE. Localnpf = ", localnpf
     lp1 = 1
     do while (lp1 .le. localnpf)

        lp2 = lp1 + 1
        do while (lp2 .le. localnpf)

           px1 = particles_global(ipx,lp1)
           py1 = particles_global(ipy,lp1)
           pz1 = particles_global(ipz,lp1)
           px2 = particles_global(ipx,lp2)
           py2 = particles_global(ipy,lp2)
           pz2 = particles_global(ipz,lp2)

           dx = px2 - px1
           dy = py2 - py1
           dz = pz2 - pz1

           if (grav_boundary_type .eq. "periodic") call pt_sinkCorrectForPeriodicBCs(dx, dy, dz)

           rad = sqrt(dx**2 + dy**2 + dz**2)

           if (rad .LT. accretion_radius_comoving) then

              pm1  = particles_global(ipm,lp1)
              pvx1 = particles_global(ipvx,lp1)
              pvy1 = particles_global(ipvy,lp1)
              pvz1 = particles_global(ipvz,lp1)
              pm2  = particles_global(ipm,lp2)
              pvx2 = particles_global(ipvx,lp2)
              pvy2 = particles_global(ipvy,lp2)
              pvz2 = particles_global(ipvz,lp2)

              if (rad .gt. 0) then
                 vrad = (dx*(pvx2-pvx1)+dy*(pvy2-pvy1)+dz*(pvz2-pvz1))/rad
              else
                 vrad = -1.E99 ! this will lead to definite merging if rad = 0
              endif

              if (vrad .lt. 0) then

                 ! The gravitational energy is calculated relative to the accretion radius, such that it is set to zero at r_accr.
                 ! The particles are not allowed to merge, if they have enough kinetic energy to reach a greater distance from each other
                 ! than their accretion radii. This means that they can only be merged, if they are bound within the accretion radii of each other.

                 if (rad .gt. 0 ) then
                    egrav = -Newton*pm1*pm2*(1.0/rad-1.0/accretion_radius_comoving)
                 else
                    egrav = -1.0e99
                 end if
                 egrav = egrav * onePlusRedshift3
                 new_mass = pm1 + pm2

                 pvcmx = (pm1*pvx1+pm2*pvx2)/new_mass
                 pvcmy = (pm1*pvy1+pm2*pvy2)/new_mass
                 pvcmz = (pm1*pvz1+pm2*pvz2)/new_mass
                 ekin  = 0.5*( pm1*((pvx1-pvcmx)**2+(pvy1-pvcmy)**2+(pvz1-pvcmz)**2) + &
                      pm2*((pvx2-pvcmx)**2+(pvy2-pvcmy)**2+(pvz2-pvcmz)**2)    )
                 ekin = ekin

                 etot = egrav + ekin

                 ! merge particles if total energy < 0
                 if (etot .lt. 0.0) then

                    if (pm1 .ge. pm2) then

                       ptag_merged_to   = int(particles_global(iptag,lp1))
                       ptag_merged_from = int(particles_global(iptag,lp2))
                       old_mass = pm1
                       lpkeep = lp1
                       lpdelete = lp2

                    else

                       ptag_merged_to   = int(particles_global(iptag,lp2))
                       ptag_merged_from = int(particles_global(iptag,lp1))
                       old_mass = pm2
                       lpkeep = lp2
                       lpdelete = lp1
                       dx = -dx
                       dy = -dy
                       dz = -dz

                    end if

                    ! mass and accretion rate
                    pm_del = particles_global(ipm,lpdelete)
                    particles_global(ipm,lpkeep) = new_mass
                    particles_global(ipmdot,lpkeep) = (particles_global(ipm,lpkeep) - & 
                         particles_global(iold_pmass,lpkeep))/dt

                    ! center of mass
                    px_old = particles_global(ipx,lpkeep)
                    py_old = particles_global(ipy,lpkeep)
                    pz_old = particles_global(ipz,lpkeep)
                    particles_global(ipx,lpkeep) = particles_global(ipx,lpkeep) + pm_del*dx/new_mass
                    particles_global(ipy,lpkeep) = particles_global(ipy,lpkeep) + pm_del*dy/new_mass
                    particles_global(ipz,lpkeep) = particles_global(ipz,lpkeep) + pm_del*dz/new_mass

                    ! center of mass velocity
                    particles_global(ipvx,lpkeep) = pvcmx
                    particles_global(ipvy,lpkeep) = pvcmy
                    particles_global(ipvz,lpkeep) = pvcmz

                    ! angular momentum
                    vx_del = particles_global(ipvx,lpdelete)
                    vy_del = particles_global(ipvy,lpdelete)
                    vz_del = particles_global(ipvz,lpdelete)
                    particles_global(iplx,lpkeep) = particles_global(iplx,lpkeep) + pm_del*(dy*vz_del-dz*vy_del) - &
                                    particles_global(ipm,lpkeep) * &
                                    ( (particles_global(ipy,lpkeep)-py_old)*particles_global(ipvz,lpkeep) - &
                                      (particles_global(ipz,lpkeep)-pz_old)*particles_global(ipvy,lpkeep)  )
                    particles_global(iply,lpkeep) = particles_global(iply,lpkeep) + pm_del*(dz*vx_del-dx*vz_del) - &
                                    particles_global(ipm,lpkeep) * &
                                    ( (particles_global(ipz,lpkeep)-pz_old)*particles_global(ipvx,lpkeep) - &
                                      (particles_global(ipx,lpkeep)-px_old)*particles_global(ipvz,lpkeep)  )
                    particles_global(iplz,lpkeep) = particles_global(iplz,lpkeep) + pm_del*(dx*vy_del-dy*vx_del) - &
                                    particles_global(ipm,lpkeep) * &
                                    ( (particles_global(ipx,lpkeep)-px_old)*particles_global(ipvy,lpkeep) - &
                                      (particles_global(ipy,lpkeep)-py_old)*particles_global(ipvx,lpkeep)  )

                    ! copy the last one in global list
                    particles_global(:,lpdelete) = particles_global(:,localnpf)
                    localnpf = localnpf - 1 ! the global list lost one particle

                    sink_merging_happened_on_masterPE = .true.

                    lp1 = 1
                    lp2 = 1

                 end if    ! bound in accretion accretion

              end if    ! vrad < 0

           end if    ! rad < accretion_radius

           lp2 = lp2 + 1
        end do

        lp1 = lp1 + 1
     end do

     print*, "done merging sink particles. merging happened=", sink_merging_happened_on_masterPE

  end if   ! myPE = MasterPE

  ! spread the word that merging happened
  call MPI_ALLREDUCE (sink_merging_happened_on_masterPE, sink_merging_happened, 1, FLASH_LOGICAL, FLASH_LOR, MPI_COMM_WORLD, ierr)
  
  ! all processors do this, if sink_merging_happened_on_masterPE = .true.
  if (sink_merging_happened) then

     ! delete all local particles on all processors
     do lp1 = 1, localnp
        particles_local(ipblk,lp1) = NONEXISTENT
        n_empty = n_empty + 1 ! delete
        localnp = localnp - 1 ! reduce local number of particles by one
     enddo

     call Particles_sinkSortParticles() ! I think we don't need that

     if (localnp .ne. 0) call Driver_abortFlash("pt_sinkParticleMerging: error while merging 3")

     ! make all particles in updated global list local for the master processor
     if (MyPE .eq. MasterPE) then

        if (localnpf .GT. sink_maxSinks) &
             call Driver_abortFlash('pt_sinkParticleMerging: Total number of sink particles exceeds sink_maxSinks.')

        do lp2 = 1, localnpf
           particles_local(:,lp2) = particles_global(:,lp2)
           n_empty = n_empty - 1 ! create
           localnp = localnp + 1 ! increase the local number of particles by one
        enddo

        if (localnp .ne. localnpf) call Driver_abortFlash("pt_sinkParticleMerging: error while merging 4")

     end if

     call Particles_sinkMoveParticles(.false.)
     call Particles_sinkSortParticles() ! I think we don't need that

     ! update particle's cpu info
     do lp1 = 1, localnp
        particles_local(ipcpu,lp1) = MyPE
     enddo

  end if

  return

end subroutine pt_sinkParticleMerging
