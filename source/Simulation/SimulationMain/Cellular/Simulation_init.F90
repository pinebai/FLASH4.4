!!****if* source/Simulation/SimulationMain/Cellular/Simulation_init
!!
!! NAME
!!
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!
!!  Simulation_init()
!!
!!
!! DESCRIPTION
!!
!!   Initialize general solution data for a planar detonation,
!!   perturbed with some noise.  This detonation should be unstable
!!   to becoming cellular at the detonation front.
!!
!! ARGUMENTS
!!
!!
!!
!! PARAMETERS
!!
!!    xhe4               mass fraction of he4
!!    xc12               mass fraction of c12
!!    xo16               mass fraction of o16
!!    rhoAmbient         density of the cold upstream material 
!!    tempAmbient        temperature of the cold upstream material
!!    velxAmbient        x-velocity of the cold upstream material
!!    rhoPerturb         density of the post shock material
!!    tempPerturb        temperature of the post shock material
!!    velxPerturb        x-velocity of the post shock material
!!    radiusPerturb      distance below which the perturbation is applied
!!    xCenterPerturb     origin of the of the perturbation
!!    yCenterPerturb     origin of the of the perturbation
!!    zCenterPerturb     origin of the of the perturbation
!!    usePseudo1d        .true. for a 1d initial configuration, with the ??
!!                          copied along the y and z directions
!!                       .false. for a spherical configuration
!!    noiseAmplitude     amplitude of the white noise added to the perturbation
!!    noiseDistance      distances above and below radiusPerturb get noise added
!!    xmax               boundary of domain
!!    xmin               boundary of domain
!!    ymax               boundary of domain
!!    ymin               boundary of domain
!!    zmax               boundary of domain
!!    zmin               boundary of domain
!!
!!    smallx             smallest allowed abundance
!!    smlrho             smallest allowed density
!!
!!  NOTES
!!    Species used in this simulation are HE4 (helium-4), C12 (carbon-12), O16 (oxygen-16)
!!
!! NOTES
!!   See paper: Timmes, FX; Zingale, M; Olson, K; Fryxell, B; The Astrophysical
!!               Journal, Nov. 10, 2000 : 543: 938-954
!!
!!***

subroutine Simulation_init()

  use Simulation_data, ONLY: sim_smallRho, sim_smallx, sim_radiusPerturb, sim_usePseudo1d, &
     sim_xhe4, sim_xc12, sim_xo16, &
     sim_rhoAmbient, sim_tempAmbient, sim_velxAmbient, &
     sim_rhoPerturb, sim_tempPerturb, sim_velxPerturb, &
     sim_noiseAmplitude, sim_noiseDistance, &
     sim_xCenterPerturb, sim_yCenterPerturb, sim_zCenterPerturb, &
     sim_xmin, sim_xmax, sim_ymin, sim_ymax, sim_zmin, sim_zmax , sim_meshMe
  use Driver_interface, ONLY : Driver_abortFlash, Driver_getMype
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Logfile_interface, ONLY : Logfile_stamp

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Multispecies.h"
#include "Eos.h"



! -----------------------------------------------------------------------------------

#ifndef FIXEDBLOCKSIZE
      call Driver_abortFlash("SORRY, Cellular not defined for non fixed block size")

! for non-fbs, you'll need to call Grid_getBlkIndexLimits to get the extent,
!  then allocate the arrays massFraction, x,y,z
  ! compute the maximum length of a vector in each coordinate direction
  ! (including guardcells)
!  q = max( blockExtent(HIGH,IAXIS), &
!           blockExtent(HIGH,JAXIS), &
!           blockExtent(HIGH,KAXIS) )
!  real, dimension(q,NSPECIES) ::  massFraction
!  real, dimension(q) :: xCoordsCell, yCoordsCell, zCoordsCell
  real             rvec(blockExtent(HIGH,IAXIS)*blockExtent(HIGH,JAXIS)*blockExtent(HIGH,KAXIS))
#endif  

     !-----------------------------------------------------------------------------
     ! grab the parameters relevant for this problem
     !-----------------------------------------------------------------------------
     call Driver_getMype(MESH_COMM, sim_meshMe)
     call RuntimeParameters_get('smallx', sim_smallx)
     call RuntimeParameters_get('smlrho', sim_smallRho)
     call RuntimeParameters_get('radiusPerturb', sim_radiusPerturb)
     call RuntimeParameters_get('xhe4', sim_xhe4)
     call RuntimeParameters_get('xc12', sim_xc12)
     call RuntimeParameters_get('xo16', sim_xo16)
     call RuntimeParameters_get('rhoAmbient', sim_rhoAmbient)
     call RuntimeParameters_get('rhoPerturb', sim_rhoPerturb)
     call RuntimeParameters_get('tempAmbient', sim_tempAmbient)
     call RuntimeParameters_get('tempPerturb', sim_tempPerturb)
     call RuntimeParameters_get('velxAmbient', sim_velxAmbient)
     call RuntimeParameters_get('velxPerturb', sim_velxPerturb)
     call RuntimeParameters_get('noiseAmplitude', sim_noiseAmplitude)
     call RuntimeParameters_get('noiseDistance', sim_noiseDistance)
     call RuntimeParameters_get('usePseudo1d', sim_usePseudo1d)
     call RuntimeParameters_get('xCenterPerturb', sim_xCenterPerturb)
     call RuntimeParameters_get('yCenterPerturb', sim_yCenterPerturb)
     call RuntimeParameters_get('zCenterPerturb', sim_zCenterPerturb)
     call RuntimeParameters_get('xmin', sim_xmin)
     call RuntimeParameters_get('xmax', sim_xmax)
     call RuntimeParameters_get('ymin', sim_ymin)
     call RuntimeParameters_get('ymax', sim_ymax)
     call RuntimeParameters_get('zmin', sim_zmin)
     call RuntimeParameters_get('zmax', sim_zmax)

!! Driver_abort stamps to logfile
     if (sim_meshMe == MASTER_PE) then
        !  dump the parameters
        if (sim_rhoAmbient .LT. 1.e4*sim_smallRho) then
           print *, 'Warning: ambient density is close to ', & 
                &              'cutoff density'
           print *, 'reset smlrho to be less than ',  & 
                &               1.e-4*sim_rhoAmbient
           call Driver_abortFlash('ERROR Simulation_init: smlrho is too high')
        endif

        if (sim_xCenterPerturb .GT. sim_xmax .OR. sim_xCenterPerturb .LT. sim_xmin) then
           print *, 'Error: xCenterPerturb must fall between xmin and xmax'
           call Driver_abortFlash('ERROR Simulation_init: xCenterPerturb must fall between xmin and xmax')
        endif

        if (sim_yCenterPerturb .GT. sim_ymax .OR. sim_yCenterPerturb .LT. sim_ymin) then
           print *, 'Error: yCenterPerturb must fall between ymin and ymax'
           call Driver_abortFlash('ERROR Simulation_init: yCenterPerturb must fall between ymin and ymax')
        endif

        if ((sim_zCenterPerturb .GT. sim_zmax .OR. sim_zCenterPerturb .LT. sim_zmin)  & 
             &              .AND. NDIM .EQ. 3) then
           print *, 'Error: zCenterPerturb must fall between zmin and zmax'
           call Driver_abortFlash('ERROR Simulation_init: zCenterPerturb must fall between zmin and zmax')
        endif

        call Logfile_stamp( "initializing for cellular detonation", &
             "[Simulation_init]")
        print *, 'flash: ', NDIM,  & 
             &           ' dimensional cellular detonation initialization'
        print *, ' '
        print *, 'helium mass fraction    = ', sim_xhe4
        print *, 'carbon mass fraction    = ', sim_xc12
        print *, 'oxygen mass fraction    = ', sim_xo16
        print *, ' '
        print *, 'upstream density        = ', sim_rhoAmbient
        print *, 'upstream temperature    = ', sim_tempAmbient
        print *, 'upstream velocity       = ', sim_velxAmbient
        print *, ' '
        print *, 'post shock density      = ', sim_rhoPerturb
        print *, 'post shock temperature  = ', sim_tempPerturb
        print *, 'post shock velocity     = ', sim_velxPerturb
        print *, ' '
        print *, 'post shock distance     = ', sim_radiusPerturb
        print *, 'x center                = ', sim_xCenterPerturb
        print *, 'y center                = ', sim_yCenterPerturb
        print *, 'z center                = ', sim_zCenterPerturb
        print *, ' '
        print *, 'noise amplitude         = ', sim_noiseAmplitude
        print *, 'noise distance          = ', sim_noiseDistance
        print *, ' '
     endif

  return

end subroutine Simulation_init


