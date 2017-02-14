!!****if* source/IO/IOParticles/hdf5/parallel/io_ptWriteParticleData
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

#include "constants.h"
#include "Flash.h"
#include "io_flash.h"

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
       io_logToIntScalarValues, io_logToIntParmValues, io_outputSplitNum, &
       io_geometry, io_setupCall, io_buildDir, io_flashRelease, &
       io_fileCreationTime, io_buildDate, io_buildMachine, io_cflags, &
       io_fflags, io_setupTimeStamp, io_buildTimeStamp, io_comm, io_outputSplitNum, &
       io_splitNumBlks, io_splitParts, io_fileFormatVersion, io_globalMe
  use IOParticles_data, ONLY : io_writeParticleSubset, io_writeParticleAll

  use Grid_data, ONLY : gr_globalOffset, gr_globalNumBlocks
  use Grid_interface, ONLY : Grid_getLocalNumBlks, Grid_sortParticles
  use io_ptInterface, ONLY : io_ptCreateSubset
  use Particles_data, ONLY : particles,pt_maxPerProc
  use Logfile_interface, ONLY : Logfile_stamp
  use Driver_interface, ONLY : Driver_abortFlash
#ifdef USE_IO_C_INTERFACE
  use iso_c_binding, ONLY : c_loc
  use io_c_interface, ONLY : io_create_dataset, io_xfer_cont_slab
#else
#define c_loc(x) x
#endif

  implicit none
  include "Flash_mpi.h"

  integer, intent(in) ::  fileID
  integer, intent(in) :: globalNumParticles, localNumParticles, particleOffset
  character (len=OUTPUT_PROP_LENGTH), intent(in) :: partAttributeLabels(NPART_PROPS)
  logical, intent(in) :: particlesToCheckpoint
  integer,parameter :: particleTypes=1
  integer,target,dimension(MAXBLOCKS,particleTypes) :: particlesPerBlk


  ! for logfile output
  character(len=MAX_STRING_LENGTH), dimension(2,2) :: strBuff
  integer :: ierr, localNumBlocks, lb, beginIndex, endIndex

  !for split outpout
  integer :: partSplitOffset, localPartOffset, localOffset, splitOffset

  !-----------------------------------------------------------------------------
  !Variables used when writing a custom subset.
  !-----------------------------------------------------------------------------
  logical :: moreSubsetsRemain
  real, target, allocatable, dimension(:,:) :: particlest
  integer :: err, s, subsetGlobalNumParticles, subsetParticleOffset
  integer, parameter :: dims = 2, nothingVariable = 0
  integer, dimension(dims) :: subsetSize
  character(len=MAX_STRING_LENGTH) :: subsetName, subsetLabelName
  character(len=OUTPUT_PROP_LENGTH), target, dimension(NPART_PROPS) :: subsetLabels
  integer :: l_numParticles, typeMatchedXfer
  character(len=*), parameter :: localnp_str = "localnp"

  !The particle files are always double precision.
  typeMatchedXfer = 0
  
  l_numParticles=localNumParticles
  call Grid_getLocalNumBlks(localNumBlocks)

  !! This call returns particles sorted by block numbers and also the 
  !! the number of particles that each block has.

  call Grid_sortParticles(particles,NPART_PROPS,l_numParticles,particleTypes,&
       pt_maxPerProc, particlesPerBlk,BLK_PART_PROP)


  !!have to add in an offset for splitting particles among files.
  if(io_outputSplitNum > 1) then
     call MPI_ALLREDUCE(particleOffset, partSplitOffset,1, FLASH_INTEGER, &
          MPI_MIN, io_comm, ierr)
     call MPI_ALLREDUCE(localNumParticles, io_splitParts,1,FLASH_INTEGER, &
          MPI_SUM, io_comm, ierr)
     localPartOffset = particleOffset - partSplitOffset

     call MPI_ALLREDUCE(gr_globalOffset, splitOffset, 1, FLASH_INTEGER, &
          MPI_MIN, io_comm, ierr)
     call MPI_ALLREDUCE(localNumBlocks, io_splitNumBlks, 1, FLASH_INTEGER, &
          MPI_SUM, io_comm, ierr)
     localOffset = gr_globalOffset - splitOffset
  else
     localPartOffset = particleOffset
     localOffset = gr_globalOffset
     io_splitNumBlks = gr_globalNumBlocks
     io_splitParts = globalNumParticles
  endif

  !If we are outputting a particle plotfile and not a checkpoint, we need 
  !to include some overall metadata, especially scalars.  
  !specifically adding sim_info, runtime parameters and scalar values


  if( .not. particlesToCheckpoint) then

     call io_h5write_header(io_globalMe, NPART_PROPS, fileID, io_geometry, &
          partAttributeLabels, io_setupCall, io_fileCreationTime, &
          io_flashRelease, io_buildDate, io_buildDir, io_buildMachine, &
          io_cflags, io_fflags, io_setupTimeStamp, io_buildTimeStamp, &
          io_fileFormatVersion, io_outputSplitNum)

     call io_prepareListsWrite()

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

  end if ! .not. particlesToCheckpoint



  call io_create_dataset(&
       io_globalMe, &
       fileID, &
       IO_FILE_HDF5, &
       IO_FLASH_INT, &
       1, &
       (/io_splitNumBlks/), &
       localnp_str, &
       len_trim(localnp_str))

  !Treat the 2d particlePerBlk array as 1d.  We are only
  !interested in the first dimension.
  call io_xfer_cont_slab(io_globalMe, &
       fileID, &
       IO_FILE_HDF5, &
       IO_WRITE_XFER, &
       typeMatchedXfer, &
       localnp_str, &
       len_trim(localnp_str), &
       IO_FLASH_INT, &
       (/MAXBLOCKS/), &
       (/0/), &
       (/localNumBlocks/), &
       (/localOffset/), &
       (/localNumBlocks/), &
       1, &
       c_loc(particlesPerBlk(1,1)), err)
  if (err /= 0) then
     call Driver_abortFlash("Error writing localnp")
  end if


  !-----------------------------------------------------------------------------
  !Runtime parameter io_writeParticleSubset controls whether a subset of the
  !particles array is written to a particle file.  Default is .false..
  !Note: A checkpoint file NEVER contains a particle subset.
  !-----------------------------------------------------------------------------
  if (.not. particlesToCheckpoint) then
     customSubset: if (io_writeParticleSubset) then

        if (localNumParticles > 0) then
           !A subset is guaranteed to contain <= localNumParticles.
           allocate (particlest(NPART_PROPS, localNumParticles), STAT=err)
           if (err /= 0) then
              call Driver_abortFlash &
                   ("[io_ptWriteParticleData]: particlest alloc fail")
           end if
        end if


        s = 0; moreSubsetsRemain = .true.
        eachSubset: do while (moreSubsetsRemain .eqv. .true.)

           !There is the possibility to create an infinite loop.  Stamp a 
           !message every iteration to make it easy to detect user created
           !infinite loops.
           call Logfile_stamp( 'On user-defined subset '//char(48+s), &
                "[io_ptWriteParticleData]")

           !It is the users responsibility to fill particlest array.  The user
           !must also set three OUT arguments describing subset "s" and one OUT
           !argument indicating whether there are more subsets after subset "s".
           call io_ptCreateSubset ( s, NPART_PROPS, &
                localNumParticles, partAttributeLabels, particlest, &
                subsetSize, subsetName, subsetLabelName, subsetLabels, &
                moreSubsetsRemain)

           !Check user-returned subset size is sane.
           if (any (subsetSize(:) < 0)) then
              print *, "[io_ptWriteParticleData]: subset size:", subsetSize
              call Driver_abortFlash &
                   ("[io_ptWriteParticleData]: invalid subset size")
           end if

           !io_writeParticleAll means we want to write complete particles'
           !data to the particle file.  Check that the standard dataset names
           !do not conflict with the user selected dataset names.
           !We do not need need to worry about a conflict if we are
           !only writing a subset.
           if (io_writeParticleAll .and. &
                (subsetName == "tracer particles" .or. &
                subsetLabelName == "particle names")) then
              call Driver_abortFlash &
                   ("[io_ptWriteParticleData]: Dataset naming conflict")
           end if

           !Find global size of subset and io_globalMe's offset into the global subset.
           call io_getParticleOffset ( subsetSize(2), &
                subsetGlobalNumParticles, subsetParticleOffset)


           !The temporary vectors must be specified in C-order.
           if (subsetGlobalNumParticles > 0) then

              !Particle labels:
              !Add one to the number of dimensions so that the particle
              !labels appear in standard FLASH file layout.
              call io_create_dataset(&
                   io_globalMe, &
                   fileID, &
                   IO_FILE_HDF5, &
                   IO_FLASH_STRING, &
                   dims+1, &
                   (/subsetSize(1),1,OUTPUT_PROP_LENGTH/), &
                   subsetLabelName, &
                   len_trim(subsetLabelName))

              call io_xfer_cont_slab(&
                   io_globalMe, &
                   fileID, &
                   IO_FILE_HDF5, &
                   IO_WRITE_XFER_MASTER_PE, &
                   typeMatchedXfer, &
                   subsetLabelName, &
                   len_trim(subsetLabelName), &
                   IO_FLASH_STRING, &
                   (/subsetSize(1),1,OUTPUT_PROP_LENGTH/), &
                   (/0,0,0/), &
                   (/subsetSize(1),1,OUTPUT_PROP_LENGTH/), &
                   (/0,0,0/), &
                   (/subsetSize(1),1,OUTPUT_PROP_LENGTH/), &
                   dims+1, &
                   c_loc(subsetLabels(1)), err)
              if (err /= 0) then
                 call Driver_abortFlash("Error writing particle labels")
              end if


              !Particle data:
              call io_create_dataset(&
                   io_globalMe, &
                   fileID, &
                   IO_FILE_HDF5, &
                   IO_FLASH_DOUBLE, &
                   dims, &
                   (/subsetGlobalNumParticles,subsetSize(1)/), &
                   subsetName, &
                   len_trim(subsetName))
              call io_xfer_cont_slab(&
                   io_globalMe, &
                   fileID, &
                   IO_FILE_HDF5, &
                   IO_WRITE_XFER, &
                   typeMatchedXfer, &
                   subsetName, &
                   len_trim(subsetName), &
                   IO_FLASH_DOUBLE, &
                   (/localNumParticles,NPART_PROPS/), &
                   (/0,0/), &
                   (/subsetSize(2),subsetSize(1)/), &
                   (/subsetParticleOffset,0/), &
                   (/subsetSize(2),subsetSize(1)/), &
                   dims, &
                   c_loc(particlest(1,1)), err)
              if (err /= 0) then
                 call Driver_abortFlash("Error writing particle data")
              end if
           end if

           s = s + 1
        end do eachSubset


        if (localNumParticles > 0) then
           deallocate (particlest, STAT=err)
           if (err /= 0) then
              call Driver_abortFlash &
                   ("[io_ptWriteParticleData]: particlest dealloc fail")
           end if
        end if

     end if customSubset
  end if


  !-----------------------------------------------------------------------------
  !Runtime parameter io_writeParticleAll controls whether the entire particle
  !array is written to a particle file.  Default is .true..  Note: A checkpoint
  !file ALWAYS contains the entire particle array.
  !-----------------------------------------------------------------------------
  if (particlesToCheckpoint .or. &
       (.not.particlesToCheckpoint .and. io_writeParticleAll)) then

     call io_h5write_particles(io_globalMe, &
          fileID, &
          io_splitParts, &
          localNumParticles, &
          NPART_PROPS, &
          localPartOffset, &
          particles, &
          partAttributeLabels)
  end if

  return

end subroutine io_ptWriteParticleData
