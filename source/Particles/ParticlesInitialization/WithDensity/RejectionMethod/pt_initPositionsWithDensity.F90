!!****if* source/Particles/ParticlesInitialization/WithDensity/RejectionMethod/pt_initPositionsWithDensity
!!
!! NAME
!!
!!  pt_initPositionsWithDensity
!!
!! SYNOPSIS
!!
!!  call pt_initPositionsWithDensity(integer(IN)  :: blockID,
!!                                   logical(OUT) :: success)
!!
!! DESCRIPTION
!!
!!    Initialize particle locations.  This version sets up passive tracer
!!      particles that are distributed randomly according to the gas density     
!!
!! ARGUMENTS
!!
!!   blockID:  ID of block in current processor
!!
!!  success:        returns .TRUE. if positions for all particles
!!                  that should be assigned to this block have been
!!                  successfully initialized.
!!
!! PARAMETERS
!!
!!  pt_numParticlesWanted:   integer  Approximate number of tracer particles to use
!!                                throughout domain ??
!!
!!***
!===============================================================================

subroutine pt_initPositionsWithDensity (blockID,success)

  use Particles_data, ONLY:  pt_meshMe, pt_numLocal, pt_geometry, &
       pt_maxPerProc, particles, pt_velNumAttrib, pt_velAttrib,pt_typeInfo
  use Particles_data, ONLY : pt_numParticlesWanted, pt_totalMass, &
       pt_totalVolume, pt_averageDensity

  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, &
    Grid_getBlkIndexLimits, Grid_getSingleCellVol, Grid_getBlkBoundBox, &
    Grid_getBlkPhysicalSize, Grid_getDeltas
  use Driver_interface, ONLY:  Driver_abortFlash
  use Logfile_interface, ONLY:  Logfile_stamp
  use Particles_interface, ONLY : Particles_mapFromMesh
  use pt_interface, ONLY : pt_initLocal
  implicit none

#include "constants.h"
#include "Flash.h"
#include "Particles.h"

  integer, INTENT(in) :: blockID
  logical, INTENT(out) :: success

  real, save :: sumRatio = 0.0
  integer       :: newParticlesThisBlock
  real          :: maxDensity, blockMass, blockVolume, ratioVolume, ratioDensity
  real          :: dx, dy, dz, dvol, rho
  real          :: xpos, ypos, zpos, xvel, yvel, zvel, radius
  integer       :: p, b, bb, i, j, k, numParticlesThisBlock, tag
  logical       :: accept

  real, dimension(:,:,:,:), pointer :: solnData

  real    :: urp  ! random number generator

  integer :: sizeX, sizeY, sizeZ,partID, ind, myType
  integer :: blockList(MAXBLOCKS)
  integer, dimension(MDIM) :: point  ! indices for a given location
  real, dimension(MDIM) :: blockSize
  real, dimension(2,MDIM) :: boundBox
  integer, dimension(2,MDIM):: blkLimits, blkLimitsGC
  real, dimension(MDIM) :: pos, delta, partAttribVec
  integer, parameter :: part_props=NPART_PROPS
  !----------------------------------------------------------------

  ! Access the mesh data for this block.
  call Grid_getBlkPtr(blockID,solnData)
  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
  call Grid_getBlkBoundBox(blockID,boundBox)
  call Grid_getDeltas(blockID,delta)
  !----------------------------------------------------------------------

  ! Compute the amount of mass in this block and find the maximum density.
  myType=1
  blockMass = 0.
  blockVolume = 0.
  maxDensity = 0.

!!#define solnData(i,x,y,z) 1.0
  do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
     point(3) = k
     do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
        point(2) = j
        do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
           point(1) = i
           call Grid_getSingleCellVol(blockID,EXTERIOR,point,dvol)
           blockMass = blockMass + solnData(DENS_VAR,i,j,k)*dvol
           blockVolume  = blockVolume + dvol
           maxDensity = max(maxDensity, solnData(DENS_VAR,i,j,k))
        enddo
     enddo
  enddo
  

  ! Use the amount of mass in the block to determine the number of particles the
  ! block will get.  It is weighted by the 
!!  Previous version did not put any particles in small blocks due to truncation
!  numParticlesThisBlock = int(anint(pt_numParticlesWanted * (blockMass / pt_totalMass)))
  ratioDensity = (blockMass / blockVolume) / pt_averageDensity
  ratioVolume = blockVolume / pt_totalVolume
  ! New version calculates the exact same thing, but with much more work ;)
  numParticlesThisBlock = int(anint(pt_numParticlesWanted * ratioVolume * ratioDensity))
!  sumRatio = sumRatio + ratioVolume*ratioDensity
!  write(*,100)blockID,ratioDensity,ratioVolume,ratioVolume*ratioDensity,sumRatio,newParticlesThisBlock
!  write(*,200)"Old",0.0,0.0,blockMass / pt_totalMass,sumRatio,numParticlesThisBlock

 100 format(I5,4G12.4,I5)
 200 format(A5,4G12.4,I5)

  numParticlesThisBlock = int(anint(pt_numParticlesWanted * (blockMass / pt_totalMass)))

  !  Check that this requested number isn't going to blow data allocation
  success=.true.
  if (pt_numLocal + numParticlesThisBlock > pt_maxPerProc) then
     call Logfile_stamp(numParticlesThisBlock,"This block mass would generate too many additional particles:")
     call Grid_releaseBlkPtr(blockID,solnData) !cleanup before premature return
     success=.false.
     pt_numLocal=0
     return
  end if
 
  ! Initialize the particles for this block.

  call Grid_getBlkPhysicalSize(blockID,blockSize)  ! physical size of the block
  call Grid_getBlkBoundBox(blockID,boundBox)       ! physical corner of the block
  call Grid_getDeltas(blockID,delta)
  p = 0
  do while (p < numParticlesThisBlock)

     ! Initialize particle positions using the rejection method:  choose uniform
     ! random x, y, and z coordinates and density.  If the density lies below the
     ! mesh density for the zone in which the coordinates lie, keep the particle;
     ! otherwise, choose another set of random values.  In determining the random
     ! choice box dimension for the density, use the max density times 1.01 to make
     ! sure we don't have any problems with roundoff in the highest-density zones
     ! in the block.
     !  NOTE that some whole complicated (and LBR doesn't understand it) procedure
     !  has been called in Particles_initWithDensity to initialize the random seed.

     call random_number (harvest=urp)
     if (pt_geometry == CYLINDRICAL .and. NDIM == 2) then
        ! use surface density in 2d plane (annulus mass)
        rho = urp * boundBox(HIGH,IAXIS) * maxDensity * 1.01
     else
        rho  = urp * maxDensity * 1.01     ! generate a random density less than maximum over this block
     endif

     ! Calculate X position and index, needed in all geometries
     call random_number (harvest=urp)

     xpos = urp * blockSize(1)             
     i    = int(xpos/delta(IAXIS)) + blkLimits(LOW,IAXIS)
     xpos = xpos + boundBox(LOW,IAXIS)

     if (NDIM >= 2) then
        call random_number (harvest=urp)
        ypos = urp * blockSize(JAXIS)
        j    = int(ypos/delta(JAXIS)) + blkLimits(LOW,JAXIS)
        ypos = ypos + boundBox(LOW,JAXIS)
     else
        ypos = 0.
        j    = blkLimits(LOW,JAXIS)
     endif

     if (NDIM == 3) then
        call random_number (harvest=urp)
        zpos = urp * blockSize(3)
        k    = int(zpos/delta(KAXIS)) + blkLimits(LOW,KAXIS)
        zpos = zpos + boundBox(LOW,KAXIS)
     else
        zpos = 0.
        k    = blkLimits(LOW,KAXIS)
     endif

     pos(IAXIS)=xpos
     pos(JAXIS)=ypos
     pos(KAXIS)=zpos

     if (pt_geometry == CYLINDRICAL .and. NDIM == 2) then
        ! use surface density in 2d plane (annulus mass)
        radius = boundBox(LOW,IAXIS)+(i-blkLimits(LOW,IAXIS)+0.5)*delta(IAXIS)
        accept = (rho <= radius*solnData(DENS_VAR,i,j,k))
     else
        accept = (rho <= solnData(DENS_VAR,i,j,k))
     endif
     if (accept) then
        
        ! now set up 
        p    = p + 1
        partID=pt_numLocal+p
        
        particles(BLK_PART_PROP,partID) = real(blockID)
        particles(PROC_PART_PROP,partID) = real(pt_meshMe)
        particles(POSX_PART_PROP:POSZ_PART_PROP,partID) = pos(IAXIS:KAXIS)
        
        
        ! Now do velocities.  Originally this code (and FLASH2) had xvel=yvel=zvel=0
        ! Alan suggested that it should mimic what Lattice did, and initialize from the mesh

        call Particles_mapFromMesh(pt_typeInfo(PART_MAPMETHOD,myType),&
             pt_velNumAttrib,pt_velAttrib,pos,&
             boundbox,delta,solnData,partAttribVec)
        
        do ind=1,pt_velNumAttrib
           particles(pt_velAttrib(PART_DS_IND,ind),partID) = partAttribVec(ind)
        end do
 

     endif
     
  enddo
  
  !---------------------------------------------------------------------
  ! Cleanup

  ! Release the mesh data pointer for this block.
  call Grid_releaseBlkPtr(blockID,solnData)

  ! Set the particle database local number of particles.

  pt_numLocal = pt_numLocal + p

  return

end subroutine pt_initPositionsWithDensity


