!!****if* source/IO/IOParticles/pnetcdf/io_ptReadParticleData
!!
!! NAME
!!
!! io_ptReadParticleData
!!
!!
!! SYNOPSIS
!!
!! io_ptReadParticleData()
!!
!!
!!
!! DESCRIPTION
!!
!!    This routine reads out the particle data in a separate hdf5 file
!!    It calls  ioh5_read_particles
!!    This is being done to make particles easy to debug and since the particles
!!    are not associated with the mesh data it makes since to separate it from
!!    the checkpoint files and plotfiles
!!
!!
!! ARGUMENTS
!!
!!
!! NOTES
!!
!!
!!***

subroutine io_ptReadParticleData()

#include "constants.h"
#include "Flash.h"
#include "io_flash.h"

  use IO_data, ONLY : io_chkptFileID, &
       io_outputSplitNum, io_globalMe, io_meshNumProcs, io_meshMe
  use Driver_interface, ONLY : Driver_abortFlash
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Particles_interface, ONLY : Particles_putLocalNum
  use Grid_interface, ONLY : Grid_getLocalNumBlks

  use Particles_data, ONLY : particles, pt_maxPerProc, pt_posInitialized

  use Grid_data, ONLY : gr_globalOffset
  use IO_interface, ONLY : IO_getScalar  
#ifdef USE_IO_C_INTERFACE
  use iso_c_binding, ONLY : c_loc
  use io_c_interface, ONLY : io_xfer_cont_slab
#else
#define c_loc(x) x
#endif

  implicit none
  include "Flash_mpi.h"

  integer :: reLocalNumParticles, ierr, particleOffset
  integer, save :: globalNumParticles !must be saved for IBM compilers
  integer :: localNumBlocks
  integer :: lb, startIndex, endIndex

  integer,target,allocatable :: particlesPerBlk(:)

  logical :: useParticles
  integer, parameter :: libType = IO_FILE_PNETCDF
  integer, parameter :: xferType = IO_READ_XFER
  character(len=*), parameter :: particles_str = "particles", &
       localnp_str = "localnp"
  integer :: typeMatchedXfer, err,  particlesPerBlkSize

  !The particle files are always double precision.
  typeMatchedXfer = 0


  pt_posInitialized = .false. !So it is set even if we return early.

  !need to get previous runtime parameter storing max particles per proc
  call RuntimeParameters_get("pt_maxPerProc", pt_maxPerProc)


  call RuntimeParameters_get("useParticles", useParticles)

  if(.not. useParticles) then
     return
  end if

  !allocate particles data structure
  allocate (particles(NPART_PROPS,pt_maxPerProc), stat=ierr)
  if (ierr /= 0) then
     call Driver_abortFlash("io_ptReadParticleData:  could not allocate particle array")
  endif

  !particles must be initialized or the entire particles algorithm will fail
  particles = NONEXISTENT

  call IO_getScalar("globalNumParticles", globalNumParticles)

  call Grid_getLocalNumBlks(localNumBlocks)

  reLocalNumParticles = 0 !Needed for Particles_putLocalNum call.


  if(globalNumParticles > 0) then
  
     particlesPerBlkSize = max(1, localNumBlocks)
     allocate(particlesPerBlk(particlesPerBlkSize))
     particlesPerBlk = 0

#ifdef FIXEDBLOCKSIZE
     call io_xfer_cont_slab(io_globalMe, &
          io_chkptFileID, &
          libType, &
          xferType, &
          typeMatchedXfer, &
          localnp_str, &
          len_trim(localnp_str), &
          IO_FLASH_INT, &
          (/particlesPerBlkSize/), &
          (/0/), &
          (/localNumBlocks/), &
          (/gr_globalOffset/), &
          (/localNumBlocks/), &
          1, &
          c_loc(particlesPerBlk(1)), err)
     if (err /= 0) then
        call Driver_abortFlash("Error reading localnp")
     end if

     !now find the new reLocalNumParticles
     do lb=1, localNumBlocks
        reLocalNumParticles = particlesPerBlk(lb) + reLocalNumParticles
     end do
#else
     reLocalNumParticles=globalNumParticles/io_meshNumProcs
     if(io_meshMe==(io_meshNumProcs-1))&
          reLocalNumParticles=globalNumParticles-reLocalNumParticles*io_meshMe
#endif
     
     if (reLocalNumParticles > pt_maxPerProc) then
        call Driver_abortFlash &
             ('[io_ptReadParticleData] ERROR: too many particles on this proc; increase pt_maxPerProc')
     end if
     
     
     !now get the particle offset
     call io_getParticleOffset(reLocalNumParticles, globalNumParticles, particleOffset)


     call io_xfer_cont_slab(io_globalMe, &
          io_chkptFileID, &
          libType, &
          xferType, &
          typeMatchedXfer, &
          particles_str, &
          len_trim(particles_str), &
          IO_FLASH_DOUBLE, &
          (/pt_maxPerProc,NPART_PROPS/), &
          (/0,0/), &
          (/reLocalNumParticles,NPART_PROPS/), &
          (/particleOffset,0/), &
          (/reLocalNumParticles,NPART_PROPS/), &
          2, &
          c_loc(particles(1,1)), err)
     if (err /= 0) then
        call Driver_abortFlash("Error reading particles")
     end if
     

     !reset particles BLK_PART_PROP because it could have changed on restart
     startIndex = 1
     do lb=1, localNumBlocks
        
        if(particlesPerBlk(lb) > 0) then
           endIndex = startIndex + particlesPerBlk(lb)
           
           particles(BLK_PART_PROP,startIndex:endIndex-1) = lb
           startIndex = endIndex
           
        end if
     end do
     
     deallocate(particlesPerBlk)

     pt_posInitialized = .true.
  end if !if globalNumParticles > 0


  call Particles_putLocalNum(reLocalNumParticles)

  return

end subroutine io_ptReadParticleData
