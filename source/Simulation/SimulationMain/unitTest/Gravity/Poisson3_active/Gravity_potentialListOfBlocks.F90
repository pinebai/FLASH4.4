!!****if* source/Simulation/SimulationMain/unitTest/Gravity/Poisson3_active/Gravity_potentialListOfBlocks
!!
!!  NAME 
!!
!!     Gravity_potentialListOfBlocks
!!
!!  SYNOPSIS
!!
!!  call Gravity_potentialListOfBlocks(integer(IN) :: blockCount,
!!                                     integer(IN) :: blockList(blockCount),
!!                            optional,integer(IN) :: potentialIndex)
!!
!!  DESCRIPTION 
!!      This routine computes the gravitational potential for the gravity
!!      implementations (i.e., various Poisson implementations) which make
!!      use of it in computing the gravitational acceleration.
!!
!!      Supported boundary conditions are isolated (0) and
!!      periodic (1).  The same boundary conditions are applied
!!      in all directions.
!!
!! ARGUMENTS
!!
!!   blockCount   : The number of blocks in the list
!!   blockList(:) : The list of blocks on which to calculate potential
!!
!!
!!
!!***

!!REORDER(4): solnVec

subroutine Gravity_potentialListOfBlocks(blockCount,blockList, potentialIndex)

  use Gravity_data, ONLY : grav_poisfact, grav_temporal_extrp, grav_boundary, &
       grav_unjunkPden, &
       useGravity, updateGravity, grv_meshComm, grv_meshMe
  use Cosmology_interface, ONLY : Cosmology_getRedshift, &
       Cosmology_getOldRedshift
  use Driver_interface, ONLY : Driver_abortFlash
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Particles_interface, ONLY: Particles_updateGridVar
  use Grid_interface, ONLY : GRID_PDE_BND_PERIODIC, GRID_PDE_BND_NEUMANN, &
       GRID_PDE_BND_ISOLATED, GRID_PDE_BND_DIRICHLET, &
       Grid_getBlkPtr, Grid_releaseBlkPtr, &
       Grid_solvePoisson

  implicit none

#include "Flash.h"
#include "constants.h"
#include "Flash_mpi.h"

  integer,intent(IN) :: blockCount
  integer,dimension(blockCount),intent(IN) :: blockList
  integer, intent(IN), optional :: potentialIndex


  real, POINTER, DIMENSION(:,:,:,:) :: solnVec

  integer       :: ierr

  real          :: redshift, oldRedshift
  real          :: scaleFactor, oldScaleFactor
  real          :: invscale, rescale
  integer       :: lb
  integer       :: bcTypes(6)
  real          :: bcValues(2,6) = 0.
  logical       :: updateGrav
  integer       :: density


  !NOTE: The difference between this version of Gravity_potentialListOfBlocks.F90
  !and the default Gravity_potentialListOfBlocks.F90 are the additional code blocks marked 
  !with "ADDITION".
  !----------------------- ADDITION ---------------------------
  real, external :: getMassMesh, getMassParticles
  real :: massNotConserved, massInMeshBeforeParticleSmear, massInParticles, &
       massInMeshAfterParticleSmear
  !----------------------- ADDITION ---------------------------


  call Cosmology_getRedshift(redshift)
  call Cosmology_getOldRedshift(oldRedshift)

  scaleFactor = 1./(1.+redshift)
  oldScaleFactor = 1./(1.+oldRedshift)

  invscale = 1./scaleFactor**3

  ! Rescaling factor to try and keep initial guess at potential close to
  ! final solution (in cosmological simulations).  Source term in Poisson
  ! equation has 1/a(t)^3 in it; and in linear theory (Omega=1, matter dom.)
  ! the comoving density increases as a(t), so comoving peculiar potential
  ! (which is what we are calculating here) should vary as 1/a(t)^2.  For
  ! noncosmological simulations this has no effect, since oldscale = scale = 1.

  rescale = (oldScaleFactor/scaleFactor)**2

  !=========================================================================

  if(.not.useGravity) return

  !call Particles_updateGravity(updateGrav)  change to this
  !call ParticleUpdateGravity(updateGrav)

  if(.not.updateGravity) return

  call Timers_start("gravity Barrier")
  call MPI_Barrier (grv_meshComm, ierr)
  call Timers_stop("gravity Barrier")

  call Timers_start("gravity")

  bcTypes = grav_boundary
  where (bcTypes == PERIODIC)
     bcTypes = GRID_PDE_BND_PERIODIC
  elsewhere (bcTypes == ISOLATED)
     bcTypes = GRID_PDE_BND_ISOLATED
  elsewhere (bcTypes == DIRICHLET)
     bcTypes = GRID_PDE_BND_DIRICHLET
  elsewhere (bcTypes == OUTFLOW)
     bcTypes = GRID_PDE_BND_NEUMANN
  end where
  bcValues = 0.

  if (grav_temporal_extrp) then

     call Driver_abortFlash("shouldn't be here right now")
     !call extrp_initial_guess( igpot, igpol, igpot )

  else 

     do lb = 1, blockCount
        call Grid_getBlkPtr(blocklist(lb), solnVec)
#ifdef GPOL_VAR
        solnVec(GPOL_VAR,:,:,:) = solnVec(GPOT_VAR,:,:,:)
#endif
        solnVec(GPOT_VAR,:,:,:) = solnVec(GPOT_VAR,:,:,:) * rescale
        call Grid_releaseBlkPtr(blocklist(lb), solnVec)
     enddo

  endif


  !----------------------- ADDITION ---------------------------
  !Check the total mass in the DENS_VAR mesh element + the total 
  !mass in all the particles.  PDEN_VAR is not considered because
  !it is zeroed on entry to Grid_mapParticlesToMesh.
  massInMeshBeforeParticleSmear = getMassMesh(grv_meshComm, DENS_VAR)
  massInParticles = getMassParticles(grv_meshComm)

  if(grv_meshMe == 0) then
     print *, ""
     print *, "BEFORE PARTICLE SMEAR..."
     print *, "Total mass in the DENS_VAR mesh element =", massInMeshBeforeParticleSmear
     print *, "Total mass in the particles =", massInParticles
  end if
  !----------------------- ADDITION ---------------------------


  ! Poisson is solved with the total density of PDEN_VAR + DENS_VAR 
  density=DENS_VAR
  ! This only gets called if there are active particles.
#ifdef PDEN_VAR
#ifdef MASS_PART_PROP
  call Particles_updateGridVar(MASS_PART_PROP, PDEN_VAR)


  density = PDEN_VAR
#ifdef DENS_VAR           
  do lb = 1, blockCount
     call Grid_getBlkPtr(blocklist(lb), solnVec)
     solnVec(density,:,:,:) = solnVec(density,:,:,:) + &
          solnVec(DENS_VAR,:,:,:)
     call Grid_releaseBlkPtr(blocklist(lb), solnVec)
  enddo
#endif
#endif
#endif


  !----------------------- ADDITION ---------------------------
  !The particles have been mapped to the mesh and the density in the 
  !DENS_VAR mesh element has been added to the PDEN_VAR mesh element. 
  !Therefore the Poisson solve should be approximately the same as the 
  !Poisson solve in the original Poisson3 unitTest.
  massInMeshAfterParticleSmear = getMassMesh(grv_meshComm, density)
  
  massNotConserved = abs((massInMeshBeforeParticleSmear + massInParticles) - &
       massInMeshAfterParticleSmear)
  
  if(grv_meshMe == 0) then
     print *, "AFTER PARTICLE SMEAR..."
     print *, "Total mass in the density mesh element =", massInMeshAfterParticleSmear
     print *, "Mass not conserved =", massNotConserved
     print *, ""

     if (massNotConserved > 1.0E-10) &
          call Driver_abortFlash("Lost mass exceeds error tolerance... exiting")
  end if

  call MPI_Barrier (grv_meshComm, ierr)
  !----------------------- ADDITION ---------------------------


  invscale=grav_poisfact*invscale
  call Grid_solvePoisson (GPOT_VAR, density, bcTypes, bcValues, &
       invscale)

  ! Un-junk PDEN if it exists and if requested.

#ifdef PDEN_VAR
  if (grav_unjunkPden) then
     density = PDEN_VAR
#ifdef DENS_VAR           
     do lb = 1, blockCount
        call Grid_getBlkPtr(blocklist(lb), solnVec)
        solnVec(density,:,:,:) = solnVec(density,:,:,:) - solnVec(DENS_VAR,:,:,:)
        call Grid_releaseBlkPtr(blocklist(lb), solnVec)
     enddo
#endif
  end if
#endif

#ifdef USEBARS
  call MPI_Barrier (grv_meshComm, ierr)
#endif  
  call Timers_stop ("gravity")


  return
end subroutine Gravity_potentialListOfBlocks



!----------------------- ADDITION ---------------------------
real function getMassMesh(comm, meshDensityVar)

  use Particles_data, ONLY : pt_numLocal, pt_maxPerProc, particles
  use Grid_interface, ONLY : Grid_getPointData, Grid_getBlkIndexLimits, &
       Grid_getSingleCellVol, &
       Grid_getListOfBlocks

#include "Flash.h"
#include "constants.h"

  implicit none
#include "Flash_mpi.h"

  integer, intent(IN) :: comm, meshDensityVar
  integer :: i, j, k, iBlk, blockID, iPart
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits,blkLimitsGC
  integer, dimension(MDIM) :: point
  real :: cellDensity, cellVolume, meshMass, meshMassSum, particleMassSum, &
       globalMeshMassSum
  integer :: blkCount, ierr
  integer, dimension(MAXBLOCKS):: blkList

  meshMass = 0.0
  meshMassSum = 0.0
  particleMassSum = 0.0


  call Grid_getListOfBlocks(LEAF,blkList,blkCount)

  do iBlk = 1, blkCount
     blockID = blkList(iBlk)
     
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
     
     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        point(3) = k
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           point(2) = j
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
              point(1) = i
              
              !Extract the density that has been placed in each cell.
              call Grid_getPointData(blockID, CENTER, meshDensityVar, EXTERIOR, point, cellDensity)
              
              !Extract the volume of each cell
              call Grid_getSingleCellVol(blockID, EXTERIOR, point, cellVolume)
              
              meshMass = cellDensity * cellVolume
              meshMassSum = meshMassSum + meshMass
              
           end do
        end do
     end do
  end do

  call mpi_allreduce(meshMassSum, globalMeshMassSum, 1, FLASH_REAL, MPI_SUM, &
       comm, ierr)

  getMassMesh = globalMeshMassSum

end function getMassMesh



real function getMassParticles(comm)

  use Particles_data, ONLY : pt_numLocal, pt_maxPerProc, particles

#include "Flash.h"
#include "constants.h"

  implicit none
#include "Flash_mpi.h"

  integer, intent(IN) :: comm
  integer :: i, j, k, iBlk, blockID, iPart
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits,blkLimitsGC
  integer, dimension(MDIM) :: point
  real :: cellDens, cellPden, cellVolume, meshMass, meshMassSum, particleMassSum, &
       globalParticleMassSum
  integer :: blkCount, ierr
  integer, dimension(MAXBLOCKS):: blkList

  meshMass = 0.0
  meshMassSum = 0.0
  particleMassSum = 0.0

  do iPart = 1, pt_numLocal
     particleMassSum = particleMassSum + particles(MASS_PART_PROP,iPart)
  end do
  
  call mpi_allreduce(particleMassSum, globalParticleMassSum, 1, FLASH_REAL, MPI_SUM, &
       comm, ierr)

  getMassParticles = globalParticleMassSum

end function getMassParticles
!----------------------- ADDITION ---------------------------
