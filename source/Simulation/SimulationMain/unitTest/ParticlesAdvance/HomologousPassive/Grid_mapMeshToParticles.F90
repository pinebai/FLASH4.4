!!****if* source/Simulation/SimulationMain/unitTest/ParticlesAdvance/HomologousPassive/Grid_mapMeshToParticles
!!
!! NAME
!!
!!  Grid_mapMeshToParticles
!!
!! SYNOPSIS
!!
!!  Grid_mapMeshToParticles(real(inout) :: particles(part_props,numParticles),
!!                          integer(in) :: part_props,
!!                          integer(in) :: numParticles,
!!                          integer(in) :: posAttrib(MDIM),
!!                          integer(in) :: numAttrib,
!!                          integer(in) :: attrib(2,numAttrib),
!!                          integer(in) :: mapType,
!!               optional,  integer(in) :: gridDataStruct)
!!
!! DESCRIPTION
!!
!!  Routine to map a quantity defined on the mesh onto the particle positions.
!!
!!  This is a special overriding implementation for testing particle time
!!  advancement methods.  Its major difference from the regular implementation
!!  (in the GridSolvers subunit) is this: The regular implementation always
!!  calls Particles_mapFromMesh for mapping grid variable varGrid to particle
!!  property propPart. This implementation replaces this behavior by the
!!  following iff the requested grid variable is one of VEL{X,Y,Z}_VAR:
!!  It calls sim_ptAnaGetVelComponent rather than Particles_mapFromMesh,
!!  to get analytically computed velocity components rather than velocity
!!  components mapped (i.e., interpolated) from the grid fluid fields.
!!
!! ARGUMENTS
!!
!!     particles:     Data structure containing particles information
!!     part_props : number of particle attributes
!!     numParticles : the number of particles on my proc
!!     posAttrib           : particles data structure indices pointing
!!                           to particle positions
!!     numAttrib           : number of attributes that need to be mapped
!!     attrib              : list of attributes and their corresponding
!!                           mesh data structure indices
!!                           processor
!!     mapType            : method for mapping grid quantities to particles
!!     gridDataStruct: The Grid data structure that varGrid refers to; one of
!!                    CENTER, FACEX, FACEY, FACEZ, GRIDVAR.  
!!                    If this argument is not present, the default is CENTER, 
!!                    so that attrib indicates a list of variable in UNK.
!!
!!  
!!***
!*******************************************************************************

subroutine Grid_mapMeshToParticles (particles, part_props,part_blkID,&
                                    numParticles,posAttrib,&
                                    numAttrib, attrib,mapType,&
                                    gridDataStruct)


  use Driver_interface, ONLY : Driver_abortFlash, Driver_getSimTime
  use Grid_interface, ONLY : Grid_getBlkPtr,Grid_releaseBlkPtr,&
       Grid_getLocalNumBlks, Grid_getBlkBoundBox, Grid_getDeltas
  use Particles_interface, ONLY : Particles_mapFromMesh

  use sim_interface
  use Simulation_data, ONLY : sim_fakeMapMeshToParticles
  use Grid_data, ONLY: gr_meshMe
  use Particles_data, ONLY: pt_velNumAttrib, pt_velAttrib !shouldn't do this
  !in a Grid_* file, but this is under Simulation, so I hope it may be excused - KW

  implicit none

#include "Flash.h"
#include "constants.h"
#include "Particles.h"

  integer, INTENT(in) :: part_props, numParticles, numAttrib, part_blkID
  real, INTENT(inout),dimension(part_props,numParticles) :: particles
  integer,dimension(MDIM), intent(IN) :: posAttrib
  integer, dimension(2,numAttrib),INTENT(in) :: attrib
  integer,  INTENT(IN) :: mapType
  integer, optional, intent(IN) :: gridDataStruct

  integer                       :: varGrid, propPart
  integer                       :: xPart, yPart, zPart
  integer :: i,j,k,currentBlk,blkCount,inBlk, prevBlk

  real, pointer, dimension(:,:,:,:) :: solnVec
  real, dimension(LOW:HIGH,MDIM) :: bndBox
  real,dimension(MDIM) :: delta, pos
  real, dimension(numAttrib) :: partAttribVec
  logical :: doFake

  real :: t

  doFake = (sim_fakeMapMeshToParticles .AND. numAttrib .GE. max(1,pt_velNumAttrib) .AND. numattrib .LE. MDIM)
  if (doFake) &
       doFake = ALL(attrib(GRID_DS_IND,1:pt_velNumAttrib) == pt_velAttrib(GRID_DS_IND,1:pt_velNumAttrib))
  if (present(gridDataStruct)) then
     if (gridDataStruct .NE. CENTER) doFake = .FALSE.
  end if

  if (doFake) call Driver_getSimTime(t)

  if(numParticles>0 .AND. .NOT. doFake) then

     call Grid_getLocalNumBlks(blkCount)

     currentBlk=int(particles(BLK_PART_PROP,1))
     prevBlk=currentBlk
     call Grid_getBlkPtr(currentBlk,solnVec,gridDataStruct)
     call Grid_getBlkBoundBox(currentBlk,bndBox)
     call Grid_getDeltas(currentBlk,delta)
     do i = 1, numParticles
#ifdef DEBUG_GRIDPARTICLES
        if((particles(BLK_PART_PROP, i) < 0) .or. (particles(BLK_PART_PROP, i) > blkCount)) then
           call Driver_abortFlash("BLK_PART_PROP out of bounds")
        end if
#endif
        
        currentBlk=int(particles(BLK_PART_PROP,i))
        if(currentBlk /= prevBlk)then
           call Grid_releaseBlkPtr(prevBlk,solnVec,gridDataStruct)
           call Grid_getBlkPtr(currentBlk,solnVec,gridDataStruct)
           call Grid_getBlkBoundBox(currentBlk,bndBox)
           call Grid_getDeltas(currentBlk,delta)
        end if
        do j = 1,MDIM
           pos(j)=particles(posAttrib(j),i)
        end do
        
        call Particles_mapFromMesh (mapType,numAttrib, attrib,&
             pos, bndBox,delta,solnVec, partAttribVec)
        do j = 1,numAttrib
           particles(attrib(PART_DS_IND,j),i)=partAttribVec(j)
        end do
        prevBlk=currentBlk
     enddo
     call Grid_releaseBlkPtr(currentBlk,solnVec,gridDataStruct)


  else if (doFake .AND. (numParticles>0)) then
     
     do i = 1, numParticles
#ifdef DEBUG_GRIDPARTICLES
        if((particles(BLK_PART_PROP, i) < 0) .or. (particles(BLK_PART_PROP, i) > blkCount)) then
           call Driver_abortFlash("BLK_PART_PROP out of bounds")
        end if
#endif
        xPart = posAttrib(IAXIS)
        yPart = posAttrib(JAXIS)
        zPart = posAttrib(KAXIS)
        do j = 1,numAttrib
           varGrid  = attrib(GRID_DS_IND,j)
           propPart = attrib(PART_DS_IND,j)
           if (varGrid > 0 .AND. propPart > 0) then
!!$              print 999,'INFO    Grid_mapMeshToParticles analytic on proc',gr_meshMe,&
!!$                   j, numAttrib, varGrid, propPart,part_props,i
              call sim_ptAnaGetVelComponent(particles, numParticles, &
                   varGrid,xPart,yPart,zPart,t,propPart,i)
           else
999           format (A,I6,': attrib #',I2,'/',I2,' in this call has varGrid=',I2,' propPart=',I2,&
                   '/',I2,' (particle #',I2,')')
              print 999,'WARNING Grid_mapMeshToParticles analytic on proc',gr_meshMe,&
                   j, numAttrib, varGrid, propPart,part_props,i
        end if
        end do
     enddo

  end if


end subroutine Grid_mapMeshToParticles

