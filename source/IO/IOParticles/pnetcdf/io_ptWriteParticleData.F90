!!****if* source/IO/IOParticles/pnetcdf/io_ptWriteParticleData
!!
!! NAME
!!
!! io_ptWriteParticleData
!!
!!
!! SYNOPSIS
!!
!! io_ptWriteParticleData(integer(in) :: fileID,
!!                      integer(in) :: globalNumParticles,
!!                      integer(in) :: localNumParticles,
!!                      integer(in) :: particleOffset,
!!                     character(in)(len=OUTPUT_PROP_LENGTH):: partAttributeLabels(NPART_PROPS)
!!                      logical(in) :: particlesToCheckpoint))
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
!!                           to a checkpoint file or to a plotfile. 
!!
!!
!! NOTES
!!   To control particle output there are a few different parameters
!!   
!!
!!***


subroutine io_ptWriteParticleData( fileID, globalNumParticles, &
     localNumParticles, particleOffset, partAttributeLabels, particlesToCheckpoint)

#include "constants.h"
#include "Flash.h"
#include "io_flash.h"

  use IO_data, ONLY : io_globalMe, io_meshMe, io_acrossMe
  use Grid_data, ONLY : gr_globalOffset, gr_globalNumBlocks
  use Grid_interface, ONLY : Grid_getLocalNumBlks, Grid_sortParticles

  use Particles_data, ONLY : particles,pt_maxPerProc
#ifdef USE_IO_C_INTERFACE
  use iso_c_binding, ONLY : c_loc
  use io_c_interface, ONLY : io_xfer_cont_slab
#else
#define c_loc(x) x
#endif

  implicit none
  include "Flash_mpi.h"

  integer, intent(in) :: fileID
  integer, intent(in) :: globalNumParticles, localNumParticles, particleOffset
  logical, intent(in) :: particlesToCheckpoint
  character (len=OUTPUT_PROP_LENGTH), intent(in) :: partAttributeLabels(NPART_PROPS)
  integer :: varid, partToChkInt

  integer :: ierr, localNumBlocks, lb, beginIndex, endIndex

  integer,parameter :: particleTypes=1
  integer,allocatable :: particlesPerBlk(:,:)
  integer :: l_numParticles, localnpElements

  integer, parameter :: libType = IO_FILE_PNETCDF
  integer, parameter :: xferType = IO_WRITE_XFER
  character(len=*), parameter :: particles_str = "particles", &
       localnp_str = "localnp"
  integer :: typeMatchedXfer, err


  if(globalNumParticles > 0) then

     !The particle files are always double precision.
     typeMatchedXfer = 0

     call Grid_getLocalNumBlks(localNumBlocks)
     allocate(particlesPerBlk(MAXBLOCKS,particleTypes))

     !! This call returns particles sorted by block numbers and also the 
     !! the number of particles that each block has.
     l_numParticles=localNumParticles     
     call Grid_sortParticles(particles,NPART_PROPS,l_numParticles,&
          particleTypes,pt_maxPerProc,&
          particlesPerBlk, BLK_PART_PROP)

     
     !change logical value to int value for c routines 
     if(particlesToCheckpoint) then
        partToChkInt = 1
     else
        partToChkInt = 0
     end if


     !varid is returned from io_ncmpi_write_part_dims.  We 
     !no longer use varid - instead we use ncmpi_inq_varid to convert
     !a variable name into a variable ID.
     call io_ncmpi_write_part_dims(fileID, &
          varid, &
          gr_globalNumBlocks, &
          NPART_PROPS, &
          globalNumParticles, &
          partToChkInt, &
          partAttributeLabels)        
     
#ifdef FIXEDBLOCKSIZE
     localnpElements = localNumBlocks
#else
     if (io_meshMe == MASTER_PE .and. io_acrossMe == 0) then
        particlesPerBlk(1,1) = globalNumParticles
        localnpElements = localNumBlocks
     else
        localnpElements = 0
     end if
#endif

     !Treat the 2d particlePerBlk array as 1d.  We are only
     !interested in the first dimension.
     call io_xfer_cont_slab(io_globalMe, &
          fileID, &
          libType, &
          xferType, &
          typeMatchedXfer, &
          localnp_str, &
          len_trim(localnp_str), &
          IO_FLASH_INT, &
          (/MAXBLOCKS/), &
          (/0/), &
          (/localnpElements/), &
          (/gr_globalOffset/), &
          (/localnpElements/), &
          1, &
          c_loc(particlesPerBlk(1,1)), err)
     if (err /= 0) then
        call Driver_abortFlash("Error writing localnp")
     end if


     call io_xfer_cont_slab(io_globalMe, &
          fileID, &
          libType, &
          xferType, &
          typeMatchedXfer, &
          particles_str, &
          len_trim(particles_str), &
          IO_FLASH_DOUBLE, &
          (/pt_maxPerProc,NPART_PROPS/), &
          (/0,0/), &
          (/localNumParticles,NPART_PROPS/), &
          (/particleOffset,0/), &
          (/localNumParticles,NPART_PROPS/), &
          2, &
          c_loc(particles(1,1)), err)
     if (err /= 0) then
        call Driver_abortFlash("Error writing particles")
     end if


     deallocate(particlesPerBlk)

  end if !if globalNumParticles > 0

end subroutine io_ptWriteParticleData
