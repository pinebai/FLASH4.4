!!****if* source/IO/IOMain/hdf5/serial/UG/io_readData
!!
!! NAME
!!
!!  io_readData
!!
!!
!! SYNOPSIS
!!
!!  io_readData()
!!
!!
!! DESCRIPTION
!!
!!  This is the reading counterpart to io_writeData.  It reads an HDF5
!!  file and distributes it to the processors to restart a simulation.
!! 
!!  this is the reading counterpart to checkpoint_wr.  It eats the HDF
!!  output and distributes it to the processors.
!!
!!  Subroutine to read checkpoint file using AMR package.
!!  Currently reads are done serially by processor 0 and data sent to
!!  other processors.
!!
!! ARGUMENTS
!!
!!
!!
!!***

!!REORDER(5): unk, facevar[xyz], unkBuf, face[XYZ]Buf


subroutine io_readData()


  use IO_data, ONLY : io_globalMe, io_globalNumProcs,io_globalComm,&
        io_realParmNamesPrev, io_realParmValuesPrev, io_numRealParmsPrev, &
       io_intParmNamesPrev, io_intParmValuesPrev, io_numIntParmsPrev, &
       io_logParmNamesPrev, io_logParmValuesPrev, io_numLogParmsPrev, &
       io_strParmNamesPrev, io_strParmValuesPrev, io_numStrParmsPrev, &
       io_realScalarNames, io_realScalarValues, io_numRealScalars, &
       io_intScalarNames, io_intScalarValues, io_numIntScalars, &
       io_logScalarNames, io_logScalarValues, io_numLogScalars, &
       io_strScalarNames, io_strScalarValues, io_numStrScalars, &
       io_logToIntScalarValues, io_logToIntParmValuesPrev, io_unklabels, &
       io_ilo, io_ihi, io_jlo, io_jhi, io_klo, io_khi, &
       io_checkpointFileNumber, io_baseName, io_outputSplitNum, io_comm, &
       io_chkptFileID, io_faceXVarLabels, io_faceYVarLabels, io_faceZVarLabels
  use Driver_interface, ONLY : Driver_abortFlash
  use RuntimeParameters_interface, ONLY : RuntimeParameters_bcast
  use Logfile_interface, ONLY : Logfile_stamp
  use Grid_interface, ONLY : Grid_putLocalNumBlks
  use IO_interface, ONLY : IO_getScalar

  use physicalData, only : unk, facevarx, facevary, facevarz 

  implicit none

#include "Flash_mpi.h"
#include "constants.h"
#include "Flash.h"


      
! create a character variable to hold the string representation of the block
! number.  Note this is set to be 4 characters long (i.e. max = 9999).  
  character (len=4) :: fnumStr

  character (len=MAX_STRING_LENGTH) :: filename
  integer :: localNumBlocks

  integer :: blockID, globalNumBlocks, globalOffset 
  integer :: jproc, i, j
  
  integer :: ierr

  integer :: status(MPI_STATUS_SIZE)
      
  character(len=4) :: recordLabel

  real, allocatable :: unkBuf(:,:,:,:,:)
  real, allocatable :: faceXBuf(:,:,:,:,:)
  real, allocatable :: faceYBuf(:,:,:,:,:)
  real, allocatable :: faceZBuf(:,:,:,:,:)


  ! for logfile output
  character(len=MAX_STRING_LENGTH), save, allocatable, dimension(:,:) :: strBuf
  character(len=16)                                    :: numToStr


  ! particle stuff.
  integer, save :: MaxParticlesPerProc
  integer       :: ReNumParts, particle_offset
  integer       :: localnpt, n_part_msgs, pcount, pmsgcount
                            


  allocate(unkBuf(UNK_VARS_BEGIN:UNK_VARS_END, NXB, NYB, NZB, 1))
#if(NFACE_VARS > 0)
    allocate(faceXBuf(NFACE_VARS,NXB+1,NYB,NZB,1))
    if(NDIM .gt. 1) allocate(faceYBuf(NFACE_VARS,NXB,NYB+1,NZB,1))
    if(NDIM .gt. 2) allocate(faceZBuf(NFACE_VARS,NXB,NYB,NZB+1,1))
#endif
    

  call io_getOutputName(io_checkpointFileNumber, "hdf5", "_chk_", filename, .false.)  
  

  if (io_globalMe == MASTER_PE) then
     call io_h5open_file_for_read(io_chkptFileID, filename, io_comm, io_outputSplitNum)  

     write (*,*) '[io_readData] Opening ', trim(filename), &
          ' for restart'
     write(*,FMT="(1X,A11)", ADVANCE="NO") 'Progress:  '
     allocate (strBuf(2, 2))
     write (strBuf(1,1), "(A)") "type"
     write (strBuf(1,2), "(A)") "checkpoint"
     write (strBuf(2,1), "(A)") "name"
     write (strBuf(2,2), "(A)") trim(filename)
     call Logfile_stamp( strBuf, 2, 2, "[io_readData] file opened")
     deallocate(strBuf)


     call io_h5read_header(io_globalMe, io_chkptFileID, io_unklabels, io_outputSplitNum) 

     !allocate space for scalar and parameter lists
     call io_prepareListsRead()

     call io_h5read_runtime_parameters(io_chkptFileID, &
          io_numRealParmsPrev, &
          io_realParmNamesPrev, &
          io_realParmValuesPrev, &
          io_numIntParmsPrev, &
          io_intParmNamesPrev, &
          io_intParmValuesPrev, &
          io_numStrParmsPrev, &
          io_strParmNamesPrev, &
          io_strParmValuesPrev, &
          io_numLogParmsPrev, &
          io_logParmNamesPrev, &
          io_logToIntParmValuesPrev)
     
     
     call io_h5read_scalars(io_chkptFileID, &
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
          io_logToIntScalarValues)
     
     ! add the lists
     call io_finalizeListsRead()

  endif

  call MPI_BARRIER(io_globalComm,ierr)  

  !all procs call broadcast
  call RuntimeParameters_bcast()

  call io_bcastScalars()


  call IO_getScalar("globalNumBlocks", globalNumBlocks)

  if (globalNumBlocks /= io_globalNumProcs) then
     call Driver_abortFlash("globalNumBlocks does not equal io_globalNumProcs &
        &  UG requires 1 block per proc in fixedblocksize mode")
  end if

  call io_checkBlockShape(globalNumBlocks)

  globalOffset = io_globalMe
  localNumBlocks = 1

  call Grid_putLocalNumBlks(localNumBlocks)

  if (io_globalMe == MASTER_PE) then
     
     do jproc=0, io_globalNumProcs-1
        
        globalOffset = jproc

        do i = UNK_VARS_BEGIN,UNK_VARS_END
           recordLabel = io_unklabels(i)
           
           call io_h5read_unknowns(io_chkptFileID, & 
                NXB, & 
                NYB, & 
                NZB, & 
                unkBuf(i,:,:,:,:), &
                recordLabel,      &
                localNumBlocks, &
                globalNumBlocks, &
                globalOffset)
        enddo

    !!READ IN FACEVARS
#if(NFACE_VARS > 0)
          do i = 1,NFACE_VARS
       
        call io_h5read_unknowns(io_chkptFileID, & 
                 NXB+1, & 
                 NYB, & 
                 NZB, & 
                 faceXBuf(i,:,:,:,:), &
                 io_faceXVarLabels(i),      &
                 localNumBlocks, &
                 globalNumBlocks, &
                 globalOffset)
        
            if(NDIM .gt. 1) then

          call io_h5read_unknowns(io_chkptFileID, & 
                   NXB, & 
                   NYB+1, & 
                   NZB, & 
                   faceYBuf(i,:,:,:,:), &
                   io_faceYVarLabels(i), &
                   localNumBlocks, &
                   globalNumBlocks, &
                   globalOffset)
        
        end if

        if(NDIM .gt. 2) then
      
          call io_h5read_unknowns(io_chkptFileID, & 
                   NXB, & 
                   NYB, & 
                   NZB+1, & 
                   faceZBuf(i,:,:,:,:), &
                   io_faceZvarLabels(i), &
                   localNumBlocks, &
                   globalNumBlocks, &
                   globalOffset)
        
        end if
          end do
#endif
        
    if (jproc /= MASTER_PE) then
           
           do i = UNK_VARS_BEGIN,UNK_VARS_END
              call MPI_SEND (unkBuf(i,:,:,:,:), &
                   NXB*NYB*NZB, & 
                   FLASH_REAL, jproc, 7+i, io_globalComm, ierr)
           end do
           
      !HANDLE FACEVAR SEND
#if(NFACE_VARS > 0)
            do i = 1,NFACE_VARS
              call MPI_SEND(faceXBuf(i,:,:,:,:), &
               (NXB+1)*NYB*NZB, &
           FLASH_REAL, jproc, &
           7+i+NUNK_VARS, &
           io_globalComm, ierr)
              if(NDIM .gt. 1) then
                call MPI_SEND(faceYBuf(i,:,:,:,:), &
                 NXB*(NYB+1)*NZB, &
             FLASH_REAL, jproc, &
             7+i+NUNK_VARS+NFACE_VARS, &
             io_globalComm, ierr)
          end if
          if(NDIM .gt. 2) then
                call MPI_SEND(faceZBuf(i,:,:,:,:), &
                 NXB*NYB*(NZB+1), &
             FLASH_REAL, jproc, &
             7+i+NUNK_VARS+(NFACE_VARS*2), &
             io_globalComm, ierr) 
          end if
            end do
      
#endif 

        else

           unk(:,io_ilo:io_ihi, io_jlo:io_jhi, io_klo:io_khi,:) = & 
                unkBuf(:,1:NXB,1:NYB,1:NZB,:)
#if(NFACE_VARS > 0)
             facevarx(:,io_ilo:io_ihi+1, io_jlo:io_jhi, io_klo:io_khi,:) = & 
                faceXBuf(:,1:NXB+1,1:NYB,1:NZB,:)
             if(NDIM .gt. 1)&
           facevary(:,io_ilo:io_ihi, io_jlo:io_jhi+1, io_klo:io_khi,:) = & 
                 faceYBuf(:,1:NXB,1:NYB+1,1:NZB,:)
             if(NDIM .gt. 2) &
           facevarz(:,io_ilo:io_ihi, io_jlo:io_jhi, io_klo:io_khi+1,:) = & 
                  faceZBuf(:,1:NXB,1:NYB,1:NZB+1,:)
#endif
         
end if
    
     end do !end jproc loop
     
  else !io_globalMe /= MASTER_PE

     do i = UNK_VARS_BEGIN,UNK_VARS_END
        call MPI_RECV(unkBuf(i,:,:,:,:), &
             NXB*NYB*NZB, & 
             FLASH_REAL, MASTER_PE, 7+i, & 
             io_globalComm, status, ierr)
        
     end do

     unk(:,io_ilo:io_ihi, io_jlo:io_jhi, io_klo:io_khi,:) = & 
          unkBuf(:,1:NXB,1:NYB,1:NZB,:)
     
     !!RECEIVE FACEVARS
#if(NFACE_VARS > 0)
       do i=1,NFACE_VARS
          
      call MPI_RECV(faceXBuf(i,:,:,:,:), &
           (NXB+1)*NYB*NZB, &
           FLASH_REAL, MASTER_PE, &
           7+i+NUNK_VARS, &
           io_globalComm, status, ierr)
           
         if(NDIM .gt. 1) then
           
       call MPI_RECV(faceYBuf(i,:,:,:,:), &
            NXB*(NYB+1)*NZB, &
            FLASH_REAL, MASTER_PE, &
            7+i+NUNK_VARS+NFACE_VARS, &
            io_globalComm, status, ierr)
     end if

     if(NDIM .gt. 2) then

       call MPI_RECV(faceZBuf(i,:,:,:,:), &
            NXB*NYB*(NZB+1), &
            FLASH_REAL, MASTER_PE, &
            7+i+NUNK_VARS+(NFACE_VARS*2), &
            io_globalComm, status, ierr)
     end if
       
       end do !end i=1,NFACE_VARS
       

     facevarx(:,io_ilo:io_ihi+1, io_jlo:io_jhi, io_klo:io_khi,:) = & 
          faceXBuf(:,1:NXB+1,1:NYB,1:NZB,:)
     if(NDIM .gt. 1) &
       facevary(:,io_ilo:io_ihi, io_jlo:io_jhi+1, io_klo:io_khi,:) = & 
          faceYBuf(:,1:NXB,1:NYB+1,1:NZB,:)
     if(NDIM .gt. 2) &
       facevarz(:,io_ilo:io_ihi, io_jlo:io_jhi, io_klo:io_khi+1,:) = & 
          faceZBuf(:,1:NXB,1:NYB,1:NZB+1,:)
       
#endif
     
  end if
  
  deallocate(unkBuf)

#if(NFACE_VARS > 0)
    deallocate(faceXBuf)
    if(NDIM .gt. 1) deallocate(faceYBuf)
    if(NDIM .gt. 2) deallocate(faceZBuf)
#endif


  call MPI_BARRIER(io_globalComm,ierr)

  return
  
end subroutine io_readData


