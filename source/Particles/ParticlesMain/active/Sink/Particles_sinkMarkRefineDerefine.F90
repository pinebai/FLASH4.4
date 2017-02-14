!!****if* source/Particles/ParticlesMain/active/Sink/Particles_sinkMarkRefineDerefine
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

  use RuntimeParameters_interface, only: RuntimeParameters_get

  implicit none
  
  logical, save :: first_call = .true.
  logical, save :: gr_refineOnSinkParticles
  logical, save :: gr_refineOnJeansLength

  if (first_call) then
     call RuntimeParameters_get("refineOnJeansLength", gr_refineOnJeansLength)
     call RuntimeParameters_get("refineOnSinkParticles", gr_refineOnSinkParticles)
     first_call = .false.
  end if
  
  !! Jeans Length:
  if (gr_refineOnJeansLength) call mark_blocks(3)
  !! Sink Particles:
  if (gr_refineOnSinkParticles) call mark_blocks(4)
  
  return
  
contains
  
  subroutine mark_blocks(input)
    use tree
    use paramesh_dimensions
    use physicaldata, ONLY : unk
    use Grid_data, ONLY : gr_maxRefine
    use Cosmology_interface, ONLY : Cosmology_getRedshift
    use Grid_interface, ONLY : Grid_getListOfBlocks, Grid_getBlkPhysicalSize, & 
         Grid_getCellCoords, Grid_getBlkIndexLimits, Grid_getBlkCenterCoords
    use RuntimeParameters_interface, ONLY : RuntimeParameters_get
    use PhysicalConstants_interface, ONLY : PhysicalConstants_get
    use Grid_interface, ONLY : Grid_getBlkPhysicalSize
    use Particles_sinkData, ONLY : localnpf, particles_global, maxsinks
    use pt_sinkInterface, only: pt_sinkGatherGlobal, pt_sinkCorrectForPeriodicBCs
    implicit none

    include "Flash_mpi.h"

#include "constants.h"
#include "Flash.h"
#include "Multispecies.h"

    ! Arguments

    integer, intent(IN) :: input

  ! Local data

    integer :: lb, kp, jp, ip, p
    real :: density, redshift

    real, dimension(:), allocatable :: xc, yc, zc
    integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
    integer :: size_x, size_y, size_z
    logical :: p_found

    real, dimension(MDIM) :: blockSize
    logical, save :: first_call = .true.
    real, save    :: accretion_radius
    real          :: accretion_radius_comoving, rad, distx, disty, distz
    real          :: dxhalf, dyhalf, dzhalf

    character(len=80), save :: grav_boundary_type

    real          :: jeans_min(MAXBLOCKS)
    real          :: cs2, maxd, jeans_numb
    real, save    :: jeans_ncells_ref, jeans_ncells_deref
    real, save    :: Newton
    real          :: comoving_density, oneplusz, oneplusz_neg2, oneplusz_cu

    integer, parameter :: gather_nprops = 3
    integer, dimension(gather_nprops), save :: gather_propinds = &
      (/ integer :: POSX_PART_PROP, POSY_PART_PROP, POSZ_PART_PROP /)

  !-------------------------------------------------------------------------------

    if (first_call) then

       call RuntimeParameters_get("sink_accretion_radius", accretion_radius)

       call RuntimeParameters_get("jeans_ncells_ref", jeans_ncells_ref)
       call RuntimeParameters_get("jeans_ncells_deref", jeans_ncells_deref)

       call PhysicalConstants_get("Newton", Newton)

       call RuntimeParameters_get("grav_boundary_type", grav_boundary_type)

       first_call = .false.

    end if


    call Cosmology_getRedshift(redshift)

    oneplusz = redshift + 1.0
    oneplusz_neg2 = oneplusz**(-2.0)
    oneplusz_cu = oneplusz**3.0

    accretion_radius_comoving = accretion_radius * oneplusz


    SELECT CASE(input)

    CASE(1)   ! OVERDESNTIY REFINEMENT

    ! not implemented in this version


    CASE(2)    ! UNDERDENSITY REFINEMENT

    ! not implemented in this version


    CASE(3)   ! JEANS LENGTH REFINEMENT

       do lb = 1, lnblocks

          jeans_min(lb) = 1.0e99

          if (nodetype(lb) .eq. 1) then

             call Grid_getBlkIndexLimits(lb, blkLimits, blkLimitsGC)

             call Grid_getBlkPhysicalSize(lb, blockSize)

             maxd = blockSize(1) / real(NXB)
             maxd = max(maxd, blockSize(2) / real(NYB))
             maxd = max(maxd, blockSize(3) / real(NZB))

             maxd = maxd / oneplusz

             do kp = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
                do jp = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
                   do ip = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

                      comoving_density = unk(DENS_VAR,ip,jp,kp,lb)

                      cs2 = unk(PRES_VAR,ip,jp,kp,lb) / comoving_density
                      cs2 = cs2 * oneplusz_neg2
                      density = comoving_density * oneplusz_cu

                      jeans_numb = sqrt(PI*cs2 / Newton / density) / maxd

                      if (jeans_numb .lt. jeans_min(lb)) then
                         jeans_min(lb) = jeans_numb
                      end if

                   enddo
                enddo
             enddo

          endif ! type of block

       enddo ! blocks

       ! label blocks for refinement/derefinement

       do lb = 1, lnblocks

          if (nodetype(lb) .eq. 1) then

              ! derefinement
              if ( (.not.refine(lb)) .and. (.not.stay(lb)) .and. & 
                   (jeans_min(lb) .gt. jeans_ncells_deref) ) then
                 derefine(lb) = .true.
              else
                 derefine(lb) = .false.
              end if

              ! refinement
              if ((jeans_min(lb) .lt. jeans_ncells_ref) .and. (lrefine(lb) .lt. gr_maxRefine)) then
                 derefine(lb) = .false.
                 refine(lb) = .true.
                 stay(lb) = .false.
              end if

              ! stay if within jeans_ncells_ref and jeans_ncells_deref
              if ((jeans_min(lb) .ge. jeans_ncells_ref) .and. (jeans_min(lb) .le. jeans_ncells_deref)) then
                 derefine(lb) = .false.
                 refine(lb) = .false.
                 stay(lb) = .true.
              end if

          end if       ! leaf blocks
       end do        ! blocks


    CASE(4)   ! SINK PARTICLE REFINEMENT

       ! update particles_global array
       call pt_sinkGatherGlobal(gather_propinds, gather_nprops)

       ! Any cell within accretion_radius of sink particle should be at the
       ! highest refinement level (its block, to be precise)

       do lb = 1, lnblocks

          if (nodetype(lb).eq.1) then

             call Grid_getBlkIndexLimits(lb, blkLimits, blkLimitsGC)
             size_x = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS) + 1
             size_y = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS) + 1
             size_z = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS) + 1

             allocate(xc(size_x))
             allocate(yc(size_y))
             allocate(zc(size_z))

             call Grid_getCellCoords(IAXIS, lb, CENTER, .true., xc, size_x)
             call Grid_getCellCoords(JAXIS, lb, CENTER, .true., yc, size_y)
             call Grid_getCellCoords(KAXIS, lb, CENTER, .true., zc, size_z)

             call Grid_getBlkPhysicalSize(lb, blockSize)

             dxhalf = blockSize(1)/real(NXB)/2.0
             dyhalf = blockSize(2)/real(NYB)/2.0
             dzhalf = blockSize(3)/real(NZB)/2.0

             ! see if particle is inside block
             do p = 1, localnpf
               if ( (particles_global(POSX_PART_PROP,p).ge.(xc(blkLimits(LOW,IAXIS))-dxhalf)) .and. &
                    (particles_global(POSX_PART_PROP,p).le.(xc(blkLimits(HIGH,IAXIS))+dxhalf)) .and. &
                    (particles_global(POSY_PART_PROP,p).ge.(yc(blkLimits(LOW,JAXIS))-dyhalf)) .and. &
                    (particles_global(POSY_PART_PROP,p).le.(yc(blkLimits(HIGH,JAXIS))+dyhalf)) .and. &
                    (particles_global(POSZ_PART_PROP,p).ge.(zc(blkLimits(LOW,KAXIS))-dzhalf)) .and. &
                    (particles_global(POSZ_PART_PROP,p).le.(zc(blkLimits(HIGH,KAXIS))+dzhalf)) ) then
                    if (lrefine(lb) .lt. gr_maxRefine) then
                        derefine(lb) = .false.
                        refine(lb) = .true.
                        stay(lb) = .false.
                    else
                        derefine(lb) = .false.
                        refine(lb) = .false.
                        stay(lb) = .true.
                    end if
               end if
             end do

             do kp = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
                do jp = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
                   do ip = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

                      ! cell within accretion radius?
                      p_found = .false.
                      do p = 1, localnpf
                         distx = xc(ip) - particles_global(POSX_PART_PROP,p)
                         disty = yc(jp) - particles_global(POSY_PART_PROP,p)
                         distz = zc(kp) - particles_global(POSZ_PART_PROP,p)
                         if (grav_boundary_type .eq. "periodic") call pt_sinkCorrectForPeriodicBCs(distx, disty, distz)
                         rad = sqrt(distx**2 + disty**2 + distz**2)
                         if (rad .le. accretion_radius_comoving) then
                            p_found = .true.
                         end if
                      end do

                      ! derefinement
                      if ((.not.p_found) .and. (.not.refine(lb)) .and. (.not.stay(lb))) then
                         derefine(lb) = .true.
                      else
                         derefine(lb) = .false.
                      end if

                      ! refinement if cell/block is within sink radius
                      if (p_found .and. (lrefine(lb) .lt. gr_maxRefine)) then
                         derefine(lb) = .false.
                         refine(lb) = .true.
                         stay(lb) = .false.
                      end if

                      ! stay at highest level if it already is
                      if (p_found .and. (lrefine(lb) .ge. gr_maxRefine)) then
                          derefine(lb) = .false.
                          refine(lb) = .false.
                          stay(lb) = .true.
                      end if

                   end do
                end do
             end do

             deallocate(xc)
             deallocate(yc)
             deallocate(zc)

          end if      ! nodetype

       end do         ! loop over blocks

    END SELECT

    return

  end subroutine mark_blocks

end subroutine Particles_sinkMarkRefineDerefine
