!!****if* source/IO/IOParticles/direct/UG/io_writeParticleData
!!
!! NAME
!!
!! io_ptWRiteParticleData
!!
!!
!!
!! SYNOPSIS
!!
!! io_ptWRiteParticleData(integer(in)           :: myPE,
!!                      integer(in)           :: numProcs,
!!                      integer(in)           :: fileID,
!!                      integer(in)           :: globalNumParticles,
!!                      integer(in)           :: localNumParticles,
!!                      integer(in)           :: particleOffset,
!!                     character(in)(len=OUTPUT_PROP_LENGTH):: partAttributeLabels(NPART_PROPS)
!!                      logical(in)           :: particlesToCheckpoint)
!!
!!
!!
!! DESCRIPTION
!!
!!  This routine does a plain old fortran write to output particles
!!  The particles can be written to a checkpoint file or to their own
!!  particle plotfile.  (The argument 'fileID' initializes either a 
!!  checkpoint or particle plotfile in io_ptWRiteCheckpoint or io_ptWRiteParticles)
!!  Particles written in double precision to a checkpoint file and 
!!  single precision for a particle plot file.
!!
!!  Each processor writes its own data to a file!  Beware this can lead to many files!
!!
!!  Note the order in which values are written out as post processing will need
!!  to understand the order and size of values exactly
!!
!!
!!
!! ARGUMENTS 
!!
!!   myPE - current processor number
!!   
!!   numProcs - number of processors in the simulation
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


subroutine io_ptWRiteParticleData( fileID, globalNumParticles, &
     localNumParticles, particleOffset, partAttributeLabels, particlesToCheckpoint)


  use Particles_data, ONLY : particles
  use Driver_interface, ONLY : Driver_getSimTime, Driver_getNStep

  implicit none

#include "constants.h"
#include "Flash_mpi.h"
#include "Flash.h"

  integer, intent(in) ::  fileID
  integer, intent(in) :: globalNumParticles, localNumParticles, particleOffset
  logical, intent(in) :: particlesToCheckpoint
  character(len=OUTPUT_PROP_LENGTH),intent(in) :: partAttributeLabels(NPART_PROPS)

  real  :: localSimTime
  integer :: localnStep
  
 ! for logfile output
  character(len=MAX_STRING_LENGTH), dimension(2,2) :: strBuff

  integer, parameter :: single = SELECTED_REAL_KIND(p=6)
  real (kind=single), allocatable :: particlesSingle(:,:)



  !if the particles are being written to the checkpoint file
  if(particlesToCheckpoint) then
     write(fileID)localNumParticles
     
     if(localNumParticles > 0) then
        !double precision for checkpoint files
        write(fileID)particles(:,1:localNumParticles)
     end if

  else !particles being written to their own particle plotfile

     !allocate space for single precision particle data
     !only writing out xpos, ypos, zpos, xvel, yvel, zvel
     allocate(particlesSingle(7, localNumParticles))

     write(fileID)localNumParticles

     call Driver_getSimTime(localSimTime)
     call Driver_getNStep(localnStep)
     
     if(localNumParticles > 0) then
        write(fileID)localSimTime
        write(fileID)localnStep
        
        particlesSingle(1,1:localNumParticles) = &
             real(particles(POSX_PART_PROP,1:localNumParticles), kind = single)
        
        particlesSingle(2,1:localNumParticles) = &
             real(particles(POSY_PART_PROP,1:localNumParticles), kind = single)
        
        particlesSingle(3,1:localNumParticles) = &
             real(particles(POSZ_PART_PROP,1:localNumParticles), kind = single)
        
        
        particlesSingle(4,1:localNumParticles) = &
             real(particles(VELX_PART_PROP,1:localNumParticles), kind = single)
        
        particlesSingle(5,1:localNumParticles) = &
             real(particles(VELY_PART_PROP,1:localNumParticles), kind = single)
        
        particlesSingle(6,1:localNumParticles) = &
             real(particles(VELZ_PART_PROP,1:localNumParticles), kind = single)

        particlesSingle(7,1:localNumParticles) = &
             real(particles(TAG_PART_PROP,1:localNumParticles), kind = single)

   
        !single precision for particle plotfiles
        write(fileID)particlesSingle

     end if

     deallocate(particlesSingle)

  end if

  

  return

end subroutine io_ptWRiteParticleData
