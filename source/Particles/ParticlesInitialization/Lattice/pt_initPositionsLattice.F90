!!****if* source/Particles/ParticlesInitialization/Lattice/pt_initPositionsLattice
!!
!! NAME
!!    pt_initPositionsLattice
!!
!! SYNOPSIS
!!
!!    call pt_initPositionsLattice(integer(in)  :: blockID,
!!                                 logical(OUT) :: success)
!!
!! DESCRIPTION
!!    Initialize particle locations.  This version sets up particles
!!      which are evenly distributed in space along the axes.  Distribution
!!      in non-Cartesian coordinate systems is not regular.
!!
!! ARGUMENTS
!!
!!  blockID:        local block ID containing particles to create
!!  success:        returns .TRUE. if positions for all particles
!!                  that should be assigned to this block have been
!!                  successfully initialized.
!!
!! PARAMETERS
!!
!!    pt_numX:      number of particles along physical x-axis of domain
!!    pt_numY:      number of particles along physical y-axis of domain
!!    pt_numZ:      number of particles along physical z-axis of domain
!!    pt_initialXMin, pt_initialXMax:  physical domain to initialize with particles in X
!!    pt_initialYMin, pt_initialYMax:  physical domain to initialize with particles in Y
!!    pt_initialZMin, pt_initialZMax:  physical domain to initialize with particles in Z
!!    pt_initialRadius: Radius of region about center of the domain to
!!                      which particle initialization is confined.  Ignored
!!                      if negative.
!!
!!***


subroutine pt_initPositionsLattice (blockID,success)

  use Particles_data, ONLY:  pt_numLocal, particles, pt_maxPerProc, &
       pt_xmin, pt_ymin, pt_zmin,pt_xmax, pt_ymax, pt_zmax,&
       pt_posAttrib,pt_velNumAttrib, pt_velAttrib,pt_typeInfo, pt_meshMe
  use Driver_interface, ONLY : Driver_abortFlash
       
  use Particles_data, ONLY : pt_numX, pt_numY, pt_numZ,&
       pt_initialZMin, pt_initialZMax, &
       pt_initialXMin, pt_initialXMax, pt_initialYMin, pt_initialYMax, &
       pt_initialRadius
  use Particles_data, ONLY : pt_geometry

  use Grid_interface, ONLY : Grid_getBlkBoundBox, Grid_mapMeshToParticles
  use Logfile_interface, ONLY:  Logfile_stamp

  implicit none
#include "constants.h"
#include "Flash.h"
#include "Particles.h"

  integer, INTENT(in) :: blockID
  logical,intent(OUT) :: success

  integer       :: i, j, k
  integer       :: p
  integer       :: blockType,mapType
  logical       :: IsInBlock,IsInSphere
  real          :: xpos, ypos, zpos, rpos, bxl, byl, bzl, bxu, byu, bzu
  real          :: xvel, yvel, zvel

! NOTE dxParticle is particle spacing, not grid spacing
  real, dimension(MDIM) :: dxParticle = 0.0
  real, dimension(2,MDIM):: boundBox
  integer :: part_props=NPART_PROPS
!----------------------------------------------------------------------

  !       Initialization now done in Particles_init.
  
  !        Particle slot number
  
  p = pt_numLocal

  ! Get grid geometry for this block 

  call Grid_getBlkBoundBox(blockID,boundBox)
  bxl = boundBox(LOW,IAXIS)
  bxu = boundBox(HIGH,IAXIS)
  
  if (NDIM >= 2) then
     byl = boundBox(LOW,JAXIS)
     byu = boundBox(HIGH,JAXIS)
  endif
  if (NDIM == 3) then
     bzl = boundBox(LOW,KAXIS)
     bzu = boundBox(HIGH,KAXIS)
  endif
  
  ! determine particle spacing as dxParticle
  dxParticle = 0.0   ! initialize to zero for all dimensions
  dxParticle(IAXIS) = (pt_initialXMax - pt_initialXMin) / pt_numX
  if (NDIM > 1) dxParticle(JAXIS) = (pt_initialYMax - pt_initialYMin) / pt_numY
  if (NDIM == 3) dxParticle(KAXIS) = (pt_initialZMax - pt_initialZMin) / pt_numZ
  
  
  !! initialization in case of lower dimensionality
  zpos = 0.0
  ypos = 0.0
  xpos = 0.0
  
  loop_x:  do i = 1, pt_numX
     xpos = (i-0.5)*dxParticle(IAXIS) + pt_initialXMin
     IsInBlock = (xpos >= bxl) .and. (xpos < bxu)
     if (.not. IsInBlock) cycle loop_x   !! skip rest of statements if not in block
     
     loop_y:  do j = 1, pt_numY
        if (NDIM >= 2) then
           ypos = (j-0.5)*dxParticle(JAXIS) + pt_initialYMin
           IsInBlock = (ypos >= byl) .and. (ypos < byu)
           if (.not. IsInBlock) cycle loop_y
        endif
        
        loop_z:  do k = 1, pt_numZ
           if (NDIM == 3) then
              zpos = (k-0.5)*dxParticle(KAXIS) + pt_initialZMin  ! location of particle in z
              IsInBlock = (zpos >= bzl) .and. (zpos < bzu)
           endif

!! Restriction to sphere only implemented for CARTESIAN, CYLINDRICAL, or SPHERICAL
!! coordinates.
           
           if ( pt_geometry == CARTESIAN ) then

             rpos = sqrt(xpos**2+ypos**2+zpos**2)

           else if ( pt_geometry == CYLINDRICAL ) then

             rpos = sqrt(xpos**2+ypos**2)

           else if ( pt_geometry == SPHERICAL ) then

             rpos = xpos

           else
           
             rpos = -1.0

           endif
           
           IsInSphere = ( pt_initialRadius <= 0.0 ) .or. ( rpos <= pt_initialRadius )

           if (IsInBlock .and. IsInSphere) then
              p = p + 1
              !! Check space allocation
              if (p > pt_maxPerProc) then
                 print *,' '
                 print *,'Block', blockID,' on',pt_meshMe,' would get too many additional particles;'
                 print *,'                 proc',pt_meshMe,' already had', pt_numLocal,' particles.'
                 print *,'PARAMETER pt_maxPerProc is set to ',pt_maxPerProc
                 call Logfile_stamp(pt_meshMe,&
                      "[pt_initPositionsLattice] This proc would get too many additional particles")
                 success=.false.
                 pt_numLocal=0
                 return
              endif
              !! particle is defined, set up data structure
              
              
              particles(BLK_PART_PROP,p) = real(blockID)
              particles(PROC_PART_PROP,p) = real(pt_meshMe)
#ifdef MASS_PART_PROP
              particles(MASS_PART_PROP,p) = 1.
#endif
              particles(POSX_PART_PROP,p) = xpos
              particles(POSY_PART_PROP,p)  = ypos
              particles(POSZ_PART_PROP,p)  = zpos
              
              
           endif   !! end of IsInBlock .and. IsInSphere is true
           
        enddo loop_z
     enddo loop_y
  enddo loop_x

  !       Setting the particle database local number of particles
  pt_numLocal = p
  !       Now initialize velocity properties for the new particles
  mapType=pt_typeInfo(PART_MAPMETHOD,1)
  call Grid_mapMeshToParticles(particles,&
       part_props,BLK_PART_PROP, pt_numLocal,&
       pt_posAttrib,pt_velNumAttrib,pt_velAttrib,mapType)

  success=.true.


  return

!----------------------------------------------------------------------
  
end subroutine pt_initPositionsLattice


