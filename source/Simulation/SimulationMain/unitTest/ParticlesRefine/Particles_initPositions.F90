!!****if* source/Simulation/SimulationMain/unitTest/ParticlesRefine/Particles_initPositions
!!
!! NAME
!!    Particles_initPositions
!!
!! SYNOPSIS
!!
!!    Particles_initPositions( logical, INTENT(out) :: success,
!!                             logical, INTENT(out) :: updateRefine)
!!
!! DESCRIPTION
!!    Initialize particle locations.  This version is designed for the 
!!    particle based poistest problem.  The particle locations are
!!    initialized according to the density given by a 
!!    Gaussian distribution function.
!!
!! ARGUMENTS
!!
!!    success : boolean indicating whether particles positions were 
!!              successfully initialized. This is not really relevant
!!              for this version of the routine.
!! updateRefine : is true if the routine wishes to retain the already
!!                initialized particles instead of reinitializing them
!!                as the grid refines.
!!
!!
!!***

!!REORDER(4): solnData

subroutine Particles_initPositions (success,updateRefine)

  use Particles_data, ONLY:  pt_numLocal, particles, pt_maxPerProc, &
       pt_posInitialized, pt_meshMe

  use Grid_interface, ONLY : Grid_getListOfBlocks, &
       Grid_getBlkPtr, Grid_releaseBlkPtr,Grid_getCellCoords,&
       Grid_getBlkIndexLimits, Grid_getBlkData

  use Simulation_data, ONLY : sim_ptMass, sim_densityThreshold, sim_print

  implicit none

#include "constants.h"
#include "Flash.h"

  logical, intent(INOUT) :: success
  logical, intent(OUT) :: updateRefine

  integer       :: i, j, k, lb,ii
  integer       :: structInd=1
  integer       :: nParticles
  integer       :: blockID

  integer       :: blkCount, ps, pe
  integer,dimension(MAXBLOCKS) :: blkList
  integer,dimension(MDIM) :: pos,dSize
  integer,dimension(LOW:HIGH, MDIM) :: blkLimits,blkLimitsGC
  real,dimension(:,:,:,:),pointer :: solnData
  real,allocatable,dimension(:) :: xCenter,yCenter,zCenter
  real,allocatable,dimension(:,:,:) :: cellVols
  real :: temp, err, mass, extraMass
  logical :: gcell=.true.

!----------------------------------------------------------------------




  updateRefine=.false.
  if(success) return

  pt_numLocal=0

  call Grid_getListOfBlocks(LEAF,blkList,blkCount)
  blockID=1
  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
  dSize(1:MDIM)=blkLimitsGC(HIGH,1:MDIM)-blkLimitsGC(LOW,1:MDIM)+1
  allocate(xCenter(dSize(IAXIS)))
  allocate(yCenter(dSize(JAXIS)))
  allocate(zCenter(dSize(JAXIS)))
  allocate(cellVols(dSize(IAXIS),dSize(JAXIS),dSize(KAXIS)))
  xCenter=0.0; yCenter=0.0; zCenter=0.0
  pos(1:MDIM)=blkLimitsGC(LOW,1:MDIM)
  do lb = 1, blkCount
     blockID=blkList(lb)
     call Grid_getBlkPtr(blockID, solnData, CENTER)
     call Grid_getCellCoords(IAXIS,blockID,CENTER,gcell,xCenter,dSize(IAXIS))
     if(NDIM>1)&
          call Grid_getCellCoords(JAXIS,blockID,CENTER,gcell,yCenter,dSize(JAXIS))
     if(NDIM>2)&
          call Grid_getCellCoords(KAXIS,blockID,CENTER,gcell,zCenter,dSize(KAXIS))

     call Grid_getBlkData(blockID, CELL_VOLUME, -1, EXTERIOR, pos, cellVols, dSize)
     do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
              mass=solnData(DENS_VAR,i,j,k)-sim_densityThreshold
              if(mass>0) then
                 mass =mass*cellVols(i,j,k)
                 nParticles=mass/sim_ptMass
                 if(nParticles==0)then
                    nParticles=1
                    extraMass=mass
                 else
                    extraMass=mass+(1-nParticles)*sim_ptMass
                 end if
                 if(pt_numLocal+nParticles>pt_maxPerProc) then
                    call Driver_abortFlash("too many particles")
                    call Grid_releaseBlkPtr(blockID, solnData)
                    if (pt_meshMe.EQ.MASTER_PE) &
                         print*,'[Particles_initPositions] no space to add more particles'
                    deallocate(xCenter)
                    deallocate(yCenter)
                    deallocate(zCenter)
                    deallocate(cellVols)
                    return
                 end if
                 
                 ps=pt_numLocal+1
                 pt_numLocal=pt_numLocal+nParticles
                 Particles(POSX_PART_PROP,ps:pt_numLocal)=xCenter(i)
                 Particles(POSY_PART_PROP,ps:pt_numLocal)=yCenter(j)
                 Particles(POSZ_PART_PROP,ps:pt_numLocal)=zCenter(k)
                 Particles(MASS_PART_PROP,ps:pt_numLocal)=sim_ptMass
                 Particles(BLK_PART_PROP,ps:pt_numLocal)=blockID
                 particles(PROC_PART_PROP,ps:pt_numLocal) = real(pt_meshMe)
                 Particles(MASS_PART_PROP,pt_numLocal)=extraMass
              end if
           end do
        end do
     end do
     call Grid_releaseBlkPtr(blockID, solnData)
  end do
  if (pt_meshMe.EQ.MASTER_PE) print*,'the local number of particles is ',pt_numLocal
  
  
  success = .true.
  pt_posInitialized = success
  
  deallocate(xCenter)
  deallocate(yCenter)
  deallocate(zCenter)
  deallocate(cellVols)
  return
  
end subroutine Particles_initPositions
