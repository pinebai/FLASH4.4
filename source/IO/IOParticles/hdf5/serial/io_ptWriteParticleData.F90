!!****if* source/IO/IOParticles/hdf5/serial/io_ptWriteParticleData
!!
!! NAME
!!
!! io_ptWriteParticleData
!!
!!
!! SYNOPSIS
!!
!! io_ptWriteParticleData(integer(in) :: fileID,
!!                     integer(in) :: globalNumParticles,
!!                     integer(in) :: localNumParticles,
!!                     integer(in) :: particleOffset,
!!                     character(in)(len=OUTPUT_PROP_LENGTH):: partAttributeLabels(NPART_PROPS)
!!                     logical(in) :: particlesToCheckpoint))
!!
!!
!!
!! DESCRIPTION
!!
!!    This routine writes out the particle data in a separate hdf5 file
!!    It calls  ioh5_write_particles
!!    This is being done to make particles easy to debug and since the particles
!!    are not associated with the mesh data it makes since to separate it from
!!    the checkpoint files and plotfiles
!!
!!
!! ARGUMENTS 
!!
!!   
!!   fileID - file handler to open file, in this case the fortran logical unit number
!!   
!!   globalNumParticles - the total number of particles in the simulation
!!   
!!   localNumParticles - number of particles on a given processor
!!   
!!   particleOffset - global particle offset, in this IO implementation
!!                    particleOffset is not used.  Each proc writes its own
!!                    particles to a unique file (either checkpoint or particle plotfile)
!!                    so no global offset is needed.  (Kept in the interface for
!!                    consistency with other IO implementations like pnetcdf and hdf5.)
!!
!!   partAttributeLabels - string label names for the particle properties.   
!!                         (like 'xpos', 'yvel', 'tag' etc)
!!
!!   particlesToCheckpoint - logical value indicating if particles should be written
!!                           to a checkpoint file or to a plotfile. (not used in this 
!!                           implementation)
!!
!!
!! NOTES
!!   To control particle output there are a few different parameters
!!   
!!
!!
!!***


subroutine io_ptWriteParticleData( fileID, globalNumParticles, &
     localNumParticles, particleOffset, partAttributeLabels, particlesToCheckpoint)
  
  
  use IO_data, ONLY : io_realParmNames, io_realParmValues, io_numRealParms, &
       io_intParmNames, io_intParmValues, io_numIntParms, &
       io_logParmNames, io_logParmValues, io_numLogParms, &
       io_strParmNames, io_strParmValues, io_numStrParms, &
       io_realScalarNames, io_realScalarValues, io_numRealScalars, &
       io_intScalarNames, io_intScalarValues, io_numIntScalars, &
       io_logScalarNames, io_logScalarValues, io_numLogScalars, &
       io_strScalarNames, io_strScalarValues, io_numStrScalars, &
       io_logToIntScalarValues, io_logToIntParmValues, &
       io_geometry, io_setupCall, io_buildDir, io_flashRelease, &
       io_fileCreationTime, io_buildDate, io_buildMachine, io_cflags, &
       io_fflags, io_setupTimeStamp, io_buildTimeStamp, io_outputSplitNum, &
       io_fileFormatVersion, io_globalMe, io_globalNumProcs

  use Grid_interface, ONLY : Grid_getLocalNumBlks, Grid_sortParticles

  use Grid_data, ONLY : gr_globalOffset, gr_globalNumBlocks

  use Particles_data, ONLY : particles, pt_maxPerProc

  implicit none

#include "constants.h"
#include "Flash_mpi.h"
#include "Flash.h"

  integer, intent(in) ::  fileID
  integer, intent(in) :: globalNumParticles, localNumParticles, particleOffset
  character (len=OUTPUT_PROP_LENGTH), intent(in) :: partAttributeLabels(NPART_PROPS)
  logical, intent(in) :: particlesToCheckpoint

 ! for logfile output
  character(len=MAX_STRING_LENGTH), dimension(2,2) :: strBuff
  integer :: ierr, blocksPerFile, localNumBlocks
  real, allocatable :: particlest(:,:)

  integer :: jproc, localnumparticlest, status(MPI_STATUS_SIZE)
  integer :: blockOffset, localNumBlockst, lb, beginIndex, endIndex
  integer,parameter :: particleTypes=1
  integer,dimension(MAXBLOCKS,particleTypes) :: particlesPerBlk, particlesPerBlkt
  integer :: partOffset
  integer :: l_numParticles


  l_numParticles=localNumParticles


  allocate(particlest(NPART_PROPS, pt_maxPerProc))

  call Grid_getLocalNumBlks(localNumBlocks)

!! This call returns particles sorted by block numbers and also the 
!! the number of particles that each block has.

  call Grid_sortParticles(particles,NPART_PROPS,l_numParticles,particleTypes,&
       pt_maxPerProc,particlesPerBlk, BLK_PART_PROP)

  if((io_globalMe == MASTER_PE) .and. (.not. particlesToCheckpoint)) then
    
    call io_h5write_header(io_globalMe, NPART_PROPS, fileID, io_geometry, &
         partAttributeLabels, io_setupCall, io_fileCreationTime, &
         io_flashRelease, io_buildDate, io_buildDir, io_buildMachine, &
         io_cflags, io_fflags, io_setupTimeStamp, io_buildTimeStamp, &
         io_fileFormatVersion, io_outputSplitNum)

    call io_prepareListsWrite()

         !! write the runtime parameters
     call io_h5write_runtime_parameters(io_globalMe, &
          fileID, &
          io_numRealParms, &
          io_realParmNames, &
          io_realParmValues, &
          io_numIntParms, &
          io_intParmNames, &
          io_intParmValues, &
          io_numStrParms, &
          io_strParmNames, &
          io_strParmValues, &
          io_numLogParms, &
          io_logParmNames, &
          io_logToIntParmValues, &
          io_outputSplitNum)
     
     !! write the scalars
     call io_h5write_scalars(io_globalMe, &
          fileID, &
          io_numRealScalars, &
          io_realScalarNames, &
          io_realScalarValues, &
          io_numIntScalars, &
          io_intScalarNames, &
          io_intScalarValues, &
          io_numStrScalars, &
          io_strScalarNames, &
          io_strScalarValues, &
          io_numLogScalars, &
          io_logScalarNames, &
          io_logToIntScalarValues, &
          io_outputSplitNum)



     call io_finalizeListsWrite()

  end if !(io_globalMe == MASTER_PE) .and. (.not. particlesToCheckpoint)

  !-----------------------------------------------------------------------------
  ! loop over all of the processors.  All the data is moved to processor 0 for
  ! storage using MPI sends and receives.
  !-----------------------------------------------------------------------------
  partOffset = 0  !in serial version, better to calculate than to communicate
  blockOffset = 0

  do jproc = 0,io_globalNumProcs-1
     
     if (io_globalMe == MASTER_PE) then
        
        ! fetch localNumblocks from other processors
        if (jproc /= 0) then
           
           call MPI_RECV (localNumParticlest,1,FLASH_INTEGER,jproc, & 
                1,MPI_COMM_WORLD,status,ierr)
           
           if (localNumParticlest > 0) then

              call MPI_RECV(particlest(1,1), NPART_PROPS*localNumParticlest, FLASH_REAL, &
                   jproc, 2, MPI_COMM_WORLD, status, ierr)

           end if !(if localNumParticlest > 0)

           
           call MPI_RECV (localNumBlockst,1,FLASH_INTEGER,jproc, & 
                3,MPI_COMM_WORLD,status,ierr)
           

           if (localNumBlockst > 0) then

              call MPI_RECV (particlesPerBlkt(1,1), localNumBlockst,FLASH_INTEGER,jproc, & 
                   4,MPI_COMM_WORLD,status,ierr)

           end if
           
        else
               
           localNumParticlest = localNumParticles
           particlest(:,1:localNumParticlest) = particles(:,1:localNumParticlest)
           localNumBlockst = localNumBlocks
           particlesPerBlkt(1:localNumBlockst,1) = particlesPerBlk(1:localNumBlockst,1)

        end if


        ! Write particles
        if (localNumParticlest > 0) then

           call io_h5write_particles(io_globalMe, &
                fileID, &
                globalNumParticles, &
                localNumParticlest, &
                NPART_PROPS, &
                partOffset, &
                particlest, &
                partAttributeLabels)      

        end if !if localNumParticlest > 0

        if(localNumBlockst > 0) then


           call io_h5write_localnp(io_globalMe, &
                fileID, &
                particlesPerBlkt, &
                localNumBlockst, &
                gr_globalNumBlocks, &
                blockOffset)

           
        endif

           
        !--------------------------------------------------------------------
        ! end local block loop
        !--------------------------------------------------------------------
        
        ! increment the global block number -- 
        !we just wrote localNumBlockst blocks from
        ! processor jproc to the output file
        partOffset = partOffset + localNumParticlest
        blockOffset = blockOffset + localNumBlockst

        
     else ! if (Io_GlobalMe == MASTER_PE)
        
        if (jproc == Io_GlobalMe) then
           
           call MPI_SEND(localNumParticles, 1, FLASH_INTEGER, 0, & 
                1, MPI_COMM_WORLD, ierr)

           if (localNumParticles > 0) then

              call MPI_SEND(particles(1,1), NPART_PROPS*localNumParticles, FLASH_REAL, &
                   0, 2, MPI_COMM_WORLD, ierr)

           endif !localNumParticles > 0
           
           call MPI_SEND(localNumBlocks, 1, FLASH_INTEGER, 0, & 
                3, MPI_COMM_WORLD, ierr)

           if (localNumBlocks > 0) then

               call MPI_SEND(particlesPerBlk(1,1), localNumBlocks, FLASH_INTEGER, 0, & 
                    4, MPI_COMM_WORLD, ierr)

           end if

        end if !jproc == Io_GlobalMe

     end if                 ! if Io_GlobalMe == MASTER_PE

     call MPI_BARRIER (MPI_COMM_WORLD, ierr)


     !------------------------------------------------------------------------
     ! end processor loop
     !------------------------------------------------------------------------
  end do
  deallocate(particlest)
  return

end subroutine io_ptWriteParticleData
