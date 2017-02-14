!!****if* source/IO/IOParticles/hdf5/parallel/io_ptReadParticleData
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
!! This routine reads the particle data from a HDF5 checkpoint file.
!!
!!
!! ARGUMENTS
!!
!!
!!
!! NOTES
!!
!! This function contains two different strategies for reading particles
!! from file.  The first is to read all particles' properties in one go
!! and the second is to read a single particle property at a time.
!! The first method is faster, but can only be used when all properties 
!! are the same in the file and the FLASH binary.  The second method is
!! slower and we have observed terrible performance in large scale BG/P
!! runs when using it with independent I/O.  Collective I/O should improve
!! things for method 2, but it will always be slower than method 1.
!!
!! Previously we used io_h5read_particles and io_h5read_single_part_prop 
!! functions for method 1 and 2.  These functions deadlocked 
!! when using collective HDF5 optimization and reading 0 particles
!! from some processors.  We now use io_xfer_cont_slab which
!! is a generic function that works even when some processors read 
!! 0 particles.  It has an optional Fortran 2003 interface which
!! guarantees interoperability and it can be used by defining 
!! USE_IO_C_INTERFACE.
!!
!! There are still some limitations when using io_h5read_localnp
!! function.  It will fail when there are more blocks than the
!! size of the localnp dataset.  This is sized in the checkpoint
!! file to be the number of processors in the original run.  
!! There is the same situation when using NoFbs and restarting
!! on more processors.
!!
!!***

#include "constants.h"
#include "Flash.h"
#include "io_flash.h"

subroutine io_ptReadParticleData()

  use IO_data, ONLY : io_outputSplitNum, &
       io_chkptFileID, io_comm, io_splitNumBlks, io_splitParts, io_globalMe,&
       io_meshNumProcs, io_meshMe
  use Simulation_interface, ONLY : Simulation_mapIntToStr
  use Driver_interface, ONLY : Driver_abortFlash
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Particles_interface, ONLY : Particles_putLocalNum
  use Grid_interface, ONLY : Grid_getLocalNumBlks
  use Grid_data, ONLY : gr_globalNumBlocks, gr_globalOffset
  use Particles_data, ONLY : particles, pt_maxPerProc, pt_posInitialized
  use IO_interface, ONLY : IO_getScalar
#ifdef USE_IO_C_INTERFACE
  use iso_c_binding, ONLY : c_loc
  use io_c_interface, ONLY : io_xfer_cont_slab
#else
#define c_loc(x) x
#endif

  implicit none
  include "Flash_mpi.h"


  integer :: ierr, i, particleOffset
  integer, save :: globalNumParticles !done for IBM compilers
  integer :: localNumBlocks, blkid, lb, startIndex, endIndex

  integer :: reLocalNumParticles  !restart local number of particles, could be different from 
  !what was written out
  integer,target,allocatable :: particlesPerBlk(:)

  logical :: useParticles, outside

  real,dimension(2,MDIM) :: bndBox !for debugging only ...

  !for file splitting
  integer :: localOffset, splitOffset, partSplitOffset, localPartOffset 

  !for property-by-property read-in.
  integer :: fileNumPartProps, j
  character(len=24), allocatable :: filePropNames(:) 
!!$  integer :: fileToCurrentMap(NPART_PROPS)
  integer,allocatable :: fileToCurrentMap(:)
  integer :: propIndex
  character(len=24) :: propString
  logical :: allParticlePropsSame
  integer :: memOffset, fileOffset
  integer, parameter :: libType = IO_FILE_HDF5
  integer, parameter :: xferType = IO_READ_XFER
  integer, parameter :: memType = IO_FLASH_DOUBLE
  integer, parameter :: partArrayDims = 2
  character(len=*), parameter :: dsetName = "tracer particles", &
       localnp_str = "localnp"
  integer, parameter :: dsetLen = len_trim(dsetName)
  integer :: typeMatchedXfer, err, particlesPerBlkSize

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
  if (.NOT.allocated(particles)) &
       allocate (particles(NPART_PROPS,pt_maxPerProc), stat=ierr)
  if (ierr /= 0) &
     call Driver_abortFlash("io_ptReadParticleData:  could not allocate particle array")

  !particles must be initialized or the entire particles algorithm will fail
  particles = NONEXISTENT

  call IO_getScalar("globalNumParticles", globalNumParticles)
  if (io_outputSplitNum /= 1) then
     call IO_getScalar("splitNumParticles", io_splitParts)
     call IO_getScalar("splitNumBlocks", io_splitNumBlks)
  else
     io_splitParts = globalNumParticles
     io_splitNumBlks = gr_globalNumBlocks
  endif

  call Grid_getLocalNumBlks(localNumBlocks)
  
  call MPI_ALLREDUCE(gr_globalOffset, splitOffset, 1, FLASH_INTEGER, &
       MPI_MIN, io_comm, ierr);
  localOffset = gr_globalOffset - splitOffset

  reLocalNumParticles = 0


  if(globalNumParticles > 0 ) then

     particlesPerBlkSize = max(1, localNumBlocks)
     allocate(particlesPerBlk(particlesPerBlkSize))
     particlesPerBlk = 0

#ifdef FIXEDBLOCKSIZE
     !return an array particlePerBlk holding the number 
     !of particles on each blk on the local proc
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
          (/localOffset/), &
          (/localNumBlocks/), &
          1, &
          c_loc(particlesPerBlk(1)), err)
     if (err /= 0) then
        call Driver_abortFlash("Error reading localnp")
     end if

     !now find the newLocalNumParticles
     do lb=1, localNumBlocks
        reLocalNumParticles = particlesPerBlk(lb) + reLocalNumParticles
     end do
#else
     reLocalNumParticles=globalNumParticles/io_meshNumProcs
     if(io_meshMe==(io_meshNumProcs-1))&
          reLocalNumParticles=globalNumParticles-reLocalNumParticles*io_meshMe
     print*,'particles are ',globalNumParticles,reLocalNumParticles
#endif

     if (reLocalNumParticles > pt_maxPerProc) then
        call Driver_abortFlash &
             ('[io_ptReadParticleData] ERROR: too many particles on this proc; increase pt_maxPerProc')
     end if


     !now get the particle offset
     call io_getParticleOffset( reLocalNumParticles, io_splitParts, particleOffset)

     call MPI_ALLREDUCE(particleOffset, partSplitOffset, 1, FLASH_INTEGER, &
          MPI_MIN, io_comm, ierr)
     localPartOffset = particleOffset - partSplitOffset



     !DEV: Changes begin here.  --PR
     !grab the number of particles in the file.

     call io_h5read_num_props(io_chkptFileID, fileNumPartProps)

     allocate(filePropNames(fileNumPartProps))
     allocate(fileToCurrentMap(fileNumPartProps))
     fileToCurrentMap = NONEXISTENT



     call io_h5read_particle_names(io_chkptFileID, filePropNames, fileNumPartProps);


     !generate mapping
     do i = 1,NPART_PROPS

        !iterate over part
        call Simulation_mapIntToStr(i, propString, MAPBLOCK_PART)
        do j = 1, fileNumPartProps

           if(propString .eq. filePropNames(j)) then
              fileToCurrentMap(j) = i
              exit
           end if

        end do

     end do


     !Test whether the particle attributes are identical in the 
     !current simulation & in the checkpoint file.
     allParticlePropsSame = .false.
     if (NPART_PROPS == fileNumPartProps) then
        do i = 1, NPART_PROPS
           if (fileToCurrentMap(i) == i) then
              allParticlePropsSame = .true.
           else
              allParticlePropsSame = .false.
              exit
           end if
        end do
     end if


     !We got extremely poor performance on BG/P during the single property
     !at a time particle read.  This is because there is a strided read in 
     !io_h5read_single_part_prop and (at the time) we used independent I/O.
     !The collective I/O optimizations will improve performance, but it
     !will always be slower than the contiguous read with collective I/O
     !optimizations.
     if (allParticlePropsSame .eqv. .true.) then
        !This is the fast I/O read.
        if (io_globalMe == 0) then
           print *, "[io_ptReadParticleData]: Starting contiguous particle read:"
        end if

        !In HDF5 storage order:
        !Dimension 0 is the number of particles.
        !Dimension 1 is the number of properties.

        call io_xfer_cont_slab(io_globalMe, &
             io_chkptFileID, &
             libType, &
             xferType, &
             typeMatchedXfer, &
             dsetName, &
             dsetLen, &
             memType, &
             (/pt_maxPerProc,NPART_PROPS/), &
             (/0,0/), &
             (/reLocalNumParticles,NPART_PROPS/), &
             (/localPartOffset,0/), &
             (/reLocalNumParticles,NPART_PROPS/), &
             partArrayDims, &
             c_loc(particles(1,1)), err)
        if (err /= 0) then
           call Driver_abortFlash("Error reading particles")
        end if

        !The commented out function is no longer used because it deadlocks           
        !in collective HDF5 mode when a processor reads 0 particles.
        !call io_h5read_particles(io_chkptFileID, &
        !     particles, &
        !     reLocalNumParticles, &
        !     NPART_PROPS, &
        !     localPartOffset)
     else

        if (io_globalMe == 0) then
           print *, "[io_ptReadParticleData]: Starting property by property particle read:"
           print *, "[io_ptReadParticleData]: (io_comm == MPI_COMM_WORLD):", & 
                (io_comm == MPI_COMM_WORLD)
        end if

        do i = 1, fileNumPartProps

           if(fileToCurrentMap(i) .eq. NONEXISTENT ) then
              if (io_globalMe == 0) then
                  print *, "[io_ptReadParticleData]: Skip chkpoint attribute:", &
                       i, "...this is not in the simulation"
              end if
              !force iteration
              cycle
           end if
           
           if (io_globalMe == 0) then
              print *, "[io_ptReadParticleData]: Reading data from chkpoint attribute:", &
                   i, "and storing in simulation attribute:", fileToCurrentMap(i)
           end if


           !Zero based property offset in memory and file.
           memOffset = fileToCurrentMap(i) - 1
           fileOffset = i - 1

           call io_xfer_cont_slab( &
                io_globalMe, &
                io_chkptFileID, &
                libType, &
                xferType, &
                typeMatchedXfer, &
                dsetName, &
                dsetLen, &
                memType, &
                (/pt_maxPerProc,NPART_PROPS/), &
                (/0,memOffset/), &
                (/reLocalNumParticles,1/), &
                (/localPartOffset,fileOffset/), &
                (/reLocalNumParticles,1/), &
                partArrayDims, &
                c_loc(particles(1,1)), err)
           if (err /= 0) then
              call Driver_abortFlash("Error reading particles")
           end if

           !The commented out function is no longer used because it deadlocks           
           !in collective HDF5 mode when a processor reads 0 particles.
           !call io_h5read_single_part_prop(io_chkptFileID, &
           !     particles, &
           !     relocalnumparticles, &
           !     fileToCurrentMap(i), &
           !     i, &
           !     fileNumPartProps, &
           !     NPART_PROPS, &
           !     localPartOffset)

           call MPI_BARRIER(io_comm, ierr)
           !print *, particles(i,:)
        end do

     end if


     deallocate(fileToCurrentMap)

     if (io_globalMe == 0) then
        print *, "[io_ptReadParticleData]: Finished reading particles from file."
     end if




     !DEV: MAKE SURE THIS IS RIGHT! vvvv -PR
     !I think this can be replaced anyway, with the new particles movement. --PR
     !reset particles BLK_PART_PROP because it could have changed on restart.
     !I am also pretty sure this is going to have to change for particles
     startIndex = 1
     do lb=1, localNumBlocks

        if(particlesPerBlk(lb) > 0) then
           endIndex = startIndex + particlesPerBlk(lb)

           particles(BLK_PART_PROP,startIndex:endIndex-1) = lb
           startIndex = endIndex

        end if
     end do


     deallocate(filePropNames)

     deallocate(particlesPerBlk)

     pt_posInitialized = .true.
  end if !if globalNumParticles > 0

  call Particles_putLocalNum(reLocalNumParticles)

#ifdef DEBUG_GRIDPARTICLES
  !check to see if all the particles have a valid BLK_PART_PROP
  do i=1, reLocalNumParticles

     blkid = int(particles(BLK_PART_PROP,i))

     if((blkid < 0) .or. (blkid > localNumBlocks)) then
        call Driver_abortFlash("io_pteadParticleData, blkid out of bounds")
     end if

!!$     !do an expensive check here to see if all particles are within
!!$     call Grid_getBlkBoundbox(int(particles(BLK_PART_PROP,i)), bndBox)
!!$     call gr_particleOutsideBndBox(bndBox, particles(POSX_PART_PROP:POSZ_PART_PROP, i), outside)
!!$     if(outside) then
!!$        print *, "particle outside bndbox ", i
!!$        call Driver_abortFlash("Error: io_ptReadParticleData, particle outside bndBox")
!!$     end if

  end do
#endif



  return

end subroutine io_ptReadParticleData
