!!****if* source/IO/IOMain/hdf5/serial/PM/io_readData
!!
!! NAME
!!
!!  io_readData
!!
!!
!! SYNOPSIS
!!
!!  call io_readData()
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
!!   none
!!
!!***


!!REORDER(5): unk, facevar[xyz], unkBuf, face[XYZ]Buf

#include "constants.h"
#include "Flash.h"
#include "io_flash.h"

subroutine io_readData()

  use io_typeInterface, ONLY : io_xfer_tree_data
 
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
       io_baseName, io_checkpointFileNumber, io_outputSplitNum, io_comm, io_chkptFileID, &
       io_faceXVarLabels, io_faceYVarLabels, io_faceZVarLabels, &
       tree_data_t
       
  use Driver_interface, ONLY : Driver_abortFlash
  use RuntimeParameters_interface, ONLY : RuntimeParameters_bcast
  use Logfile_interface, ONLY : Logfile_stamp
  use Grid_interface, ONLY : Grid_putLocalNumBlks, Grid_receiveInputData
  use IO_interface, ONLY : IO_getScalar
      
  
  use Grid_data, ONLY : gr_gid, gr_nToLeft, gr_globalOffset, &
       gr_globalNumBlocks

#ifdef FLASH_GRID_PARAMESH
  use tree, ONLY : maxblocks_tr, nfaces, nchild, nodetype, &
       lrefine, bnd_box, coord, bsize, neigh, parent, child, &
       lrefine_max
#ifdef FLASH_GRID_PARAMESH3OR4
  use tree, ONLY : MFLAGS, which_child, bflags
  use Grid_data, ONLY : gr_gsurr_blks, gr_is_gsurr_blks_initialized
#endif
#endif

  use physicaldata, ONLY : unk, facevarx, facevary, facevarz

!  use mpi
  implicit none

#include "Flash_mpi.h"

  type(tree_data_t) :: tree_data
      
  ! create a character variable to hold the string representation of the block
  ! number.  Note this is set to be 4 characters long (i.e. max = 9999).  
  character (len=4) :: fnumStr
  character (len=MAX_STRING_LENGTH) :: filename
  integer :: localNumBlocks, localNumBlockst

  integer :: blockID, ngid 
  integer :: jproc, i, j
  
  integer :: alnblocks !avg local num blocks 
  integer :: lnblockst !temp local num blocks
  integer :: ierr
  integer :: procBlocks
  integer :: globalOffsett

  ! temporary
  integer, target :: lrefinet(MAXBLOCKS), nodetypet(MAXBLOCKS)
  integer, target :: gidt(nfaces+1+nchild,MAXBLOCKS)
  real, target :: coordt(MDIM,MAXBLOCKS), bsizet(MDIM,MAXBLOCKS)
  real, target ::bnd_boxt(2,MDIM,MAXBLOCKS)

!  real :: unkBuf(NUNK_VARS, NXB, NYB, NZB, MAXBLOCKS)
  real,allocatable :: unkBuf(:,:,:,:,:)
  real,allocatable :: faceXBuf(:,:,:,:,:)
  real,allocatable :: faceYBuf(:,:,:,:,:)
  real,allocatable :: faceZBuf(:,:,:,:,:)

  integer :: xx, yy
  
  integer :: status(MPI_STATUS_SIZE)
      
  character(len=4) :: recordLabel

  ! for logfile output
  character(len=MAX_STRING_LENGTH), save, allocatable, dimension(:,:) :: strBuf
  character(len=16)                                    :: numToStr
  
#ifdef FLASH_GRID_PARAMESH3OR4
  !data structures specific to PM3f
  integer, target :: bflagst(mflags, MAXBLOCKS), which_childt(MAXBLOCKS), &
       gsurr_blkst(2,1+(K1D*2),1+(K2D*2),1+(K3D*2),MAXBLOCKS)
#endif
  integer, parameter :: presentDims = MDIM

#ifdef FLASH_GRID_PARAMESH4DEV_SURR_BLKS_OPTION
  logical, parameter :: do_gsurr_blks_read = .true.
#else
  logical, parameter :: do_gsurr_blks_read = .false.
#endif

  tree_data % bnd_box => bnd_boxt
  tree_data % coord => coordt
  tree_data % bsize => bsizet
  tree_data % gid => gidt
  tree_data % nodetype => nodetypet
  tree_data % lrefine => lrefinet
#ifdef FLASH_GRID_PARAMESH3OR4
  tree_data % bflags => bflagst
  tree_data % which_child => which_childt
#else
  nullify(tree_data % bflags)
  nullify(tree_data % which_child)
#endif
  nullify(tree_data % procnumber)
# ifdef FLASH_GRID_PARAMESH4DEV_SURR_BLKS_OPTION
  tree_data % gsurr_blks => gsurr_blkst
# else
  nullify(tree_data % gsurr_blks)
# endif

  call io_getOutputName(io_checkpointFileNumber, "hdf5", "_chk_", filename, .false.)


  if (io_globalMe == MASTER_PE) then
     write (*,*) '[io_readData] Opening ', trim(filename), &
          ' for restart'
     write(*,FMT="(1X,A11)", ADVANCE="NO") 'Progress'
     allocate (strBuf(2, 2))
     write (strBuf(1,1), "(A)") "type"
     write (strBuf(1,2), "(A)") "checkpoint"
     write (strBuf(2,1), "(A)") "name"
     write (strBuf(2,2), "(A)") trim(filename)
     call Logfile_stamp( strBuf, 2, 2, "[io_readData] file opened")
     deallocate(strBuf)
  end if


  if (io_globalMe == MASTER_PE) then
     call io_h5open_file_for_read(io_chkptFileID, filename, io_comm, io_outputSplitNum)  

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




  call IO_getScalar("globalNumBlocks", gr_globalNumBlocks)
  call io_checkBlockShape(gr_globalNumBlocks)
  
  ! compute the approximate number of blocks per processor
  alnblocks = int(gr_globalNumBlocks/io_globalNumProcs) + 1
  
  yy = (io_globalNumProcs*alnblocks) - gr_globalNumBlocks
  xx = io_globalNumProcs - yy
  
  ! check for error
  if ((alnblocks > MAXBLOCKS .AND. xx /= 0) .OR. &
       (alnblocks > MAXBLOCKS+1 .AND. xx == 0)) then
     
     print *
     print *, '********** ERROR in CHECKPOINT_FILENAME_RE ************'
     print *
     print *,' Number of blocks per processor exceeds maxblocks.'
     print *,' Suggest you reset maxblocks to a larger number or'
     print *,' run on a larger number of processors.'
     print *,' globalNumBlocks, io_globalNumProcs = ', gr_globalNumBlocks, io_globalNumProcs
     print *
     
     allocate (strBuf(5, 2))
     write (strBuf(1,1), "(A)") "type"
     write (strBuf(1,2), "(A)") "checkpoint"
     write (strBuf(2,1), "(A)") "name"
     write (strBuf(2,2), "(A)") trim(filename)
     
     write (strBuf(3,1), "(A)") "error"
     write (strBuf(3,2), "(A)") 'number of blocks per processor exceeds maxblocks'
     
     write (numToStr, "(I6)") gr_globalNumBlocks
     write (strBuf(4,1), "(A)") "globalNumBlocks"
     write (strBuf(4,2), "(A)") trim(adjustl(numToStr))
     
     write (numToStr, "(I6)") io_globalNumProcs
     write (strBuf(5,1), "(A)") "io_globalNumProcs"
     write (strBuf(5,2), "(A)") trim(adjustl(numToStr))
     
     call Logfile_stamp( strBuf, 5, 2, "[io_readData] error")
     if (allocated(strBuf)) deallocate(strBuf)
     
     call Driver_abortFlash("[io_readData] ERROR: &
          & number of blocks per processor exceeds maxblocks")

  end if
  

! loop over all the processor numbers and figure out how many blocks are
! stored to the left of the processor -- this is a little tricky
  
  gr_nToLeft(0) = 0
  
  do i = 0, io_globalNumProcs - 2
     if (i .LT. xx) then
        procBlocks = alnblocks
     else
        procBlocks = alnblocks - 1
     endif
     
     if (alnblocks .EQ. 0) then
        if (i .LT. gr_globalNumBlocks) then
           procBlocks = 1
        else
           procBlocks = 0
        end if
     end if
     
     ! we have the number of blocks on proc i, the number of blocks on i+1 is
     ! the number of blocks on i + the number of blocks left of i
     if (i .EQ. 0) then
        gr_nToLeft(i+1) = procBlocks
     else
        gr_nToLeft(i+1) = procBlocks + gr_nToLeft(i)
     endif
  enddo
  


  ! figure out how many blocks are on the current proc.
  if (io_globalMe < xx) then
     localNumBlocks = alnblocks
  else
     localNumBlocks = alnblocks - 1
  endif
  
  if (alnblocks .EQ. 0) then
     if (io_globalMe < gr_globalNumBlocks) then
        localNumBlocks = 1
     else
        localNumBlocks = 0
     end if
  end if
  
  ! compute the offset into the dataspace in the HDF5 file
  gr_globalOffset = gr_nToLeft(io_globalMe)

  call Grid_putLocalNumBlks(localNumBlocks)
  

  if (io_globalMe == MASTER_PE) then
     
     do jproc=0, io_globalNumProcs -1
        
        if (jproc /= MASTER_PE) then
           
           call MPI_RECV(localNumBlockst, 1, FLASH_INTEGER, & 
                jproc, 1, io_globalComm, status, ierr)
           
           call MPI_RECV(globalOffsett, 1, FLASH_INTEGER, & 
                jproc, 2, io_globalComm, status, ierr)

        else
           localNumBlockst = localNumBlocks
           globalOffsett = gr_globalOffset
        end if


        call io_xfer_tree_data(tree_data, io_chkptFileID, IO_FILE_HDF5, &
             IO_READ_XFER_MASTER_PE, &
             localNumBlockst, globalOffsett, presentDims)


        if (localNumBlockst > 0) then
           allocate(unkBuf(NUNK_VARS, NXB, NYB, NZB, localNumBlockst))
           do i = UNK_VARS_BEGIN,UNK_VARS_END
              recordLabel = io_unklabels(i)

              call io_h5read_unknowns(io_chkptFileID, & 
                                    NXB, & 
                                    NYB, & 
                                    NZB, & 
                                    unkBuf(i,:,:,:,:), &
                                    recordLabel,      &
                                    localNumBlockst, &
                                    gr_globalNumBlocks, &
                                    globalOffsett)
           enddo
#ifdef FLASH_GRID_PARAMESH3OR4
       !Face-centered variables
#if (NFACE_VARS > 0)
        allocate(faceXBuf(NFACE_VARS, NXB+1, NYB, NZB, localNumBlockst))
        if(NDIM .gt. 1) allocate(faceYBuf(NFACE_VARS,NXB,NYB+1,NZB,localNumBlockst))
        if(NDIM .gt. 2) allocate(faceZBuf(NFACE_VARS,NXB,NYB,NZB+1,localNumBlockst))
        
        do i = 1, NFACE_VARS
           
           call io_h5read_unknowns(io_chkptFileID,&
                NXB+1, &
                NYB, &
                NZB, &
                faceXBuf(i,:,:,:,:), &
                io_faceXVarLabels(i), &
                localNumBlockst,&
                gr_globalNumBlocks, &
                globalOffsett)
           if(NDIM .gt. 1) then
              
              call io_h5read_unknowns(io_chkptFileID,&
                   NXB, &
                   NYB+1, &
                   NZB, &
                   faceYBuf(i,:,:,:,:), &
                   io_faceYVarLabels(i), &
                   localNumBlockst,&
                   gr_globalNumBlocks, &
                   globalOffsett)
              
              if(NDIM .gt. 2) then
                 
                 call io_h5read_unknowns(io_chkptFileID,&
                      NXB, &
                      NYB, &
                      NZB+1, &
                      faceZBuf(i,:,:,:,:), &
                      io_faceZVarLabels(i), &
                      localNumBlockst,&
                      gr_globalNumBlocks, &
                      globalOffsett)
                 
              end if
           end if
        end do
#endif
#endif
      
           if (jproc /= MASTER_PE) then
              
              call MPI_SEND(coordt(1,1), MDIM*localNumBlockst, & 
                   FLASH_REAL, jproc, 2, io_globalComm, ierr)

              call MPI_SEND(bsizet(1,1), MDIM*localNumBlockst, & 
                   FLASH_REAL, jproc, 3, io_globalComm, ierr)

              call MPI_SEND(bnd_boxt(1,1,1), 2*MDIM*localNumBlockst, & 
                   FLASH_REAL, jproc, 4, io_globalComm, ierr)

              call MPI_SEND(lrefinet(1), localNumBlockst, & 
                   FLASH_INTEGER, jproc, 5, io_globalComm, ierr)

              call MPI_SEND(nodetypet(1), localNumBlockst, & 
                   FLASH_INTEGER, jproc, 6, io_globalComm, ierr)

              call MPI_SEND (gidt(1,1), localNumBlockst*(nfaces+1+nchild), & 
                   FLASH_INTEGER, jproc, 7, io_globalComm, ierr)

#ifdef FLASH_GRID_PARAMESH3OR4
              call MPI_SEND(which_childt(1), localNumBlockst, & 
                   FLASH_INTEGER, jproc, 8, io_globalComm, ierr)

              call MPI_SEND(bflagst(1,1), localNumBlockst*MFLAGS, & 
                   FLASH_INTEGER, jproc, 9, io_globalComm, ierr)

              if (do_gsurr_blks_read) then
                 call MPI_SEND(gsurr_blkst(1,1,1,1,1), &
                      2*(1+(K1D*2))*(1+(K2D*2))*(1+(K3D*2))*localNumBlockst, &
                      FLASH_INTEGER, jproc, 9999, io_globalComm, ierr)
              end if
#endif

              do i = UNK_VARS_BEGIN,UNK_VARS_END
                 call MPI_SEND (unkBuf(i,:,:,:,1:localNumBlockst), &
                      NXB*NYB*NZB*localNumBlockst, & 
                      FLASH_REAL, jproc, 9+i, io_globalComm, ierr)
              end do
#ifdef FLASH_GRID_PARAMESH3OR4

#if(NFACE_VARS > 0)
        do i = 1, NFACE_VARS
           call MPI_SEND (faceXBuf(i,:,:,:,1:localNumBlockst), &
                          (NXB+1)*NYB*NZB*localNumBlockst, &
                  FLASH_REAL, jproc, &
                  9+NUNK_VARS+i, io_globalComm, ierr)
        if(NDIM .gt. 1) then
           call MPI_SEND (faceYBuf(i,:,:,:,1:localNumBlockst), &
                          NXB*(NYB+1)*NZB*localNumBlockst, &
                  FLASH_REAL, jproc, &
                  9+NUNK_VARS+NFACE_VARS+i, &
                  io_globalComm, ierr)
          if(NDIM .gt. 2) then
             call MPI_SEND (faceZBuf(i,:,:,:,1:localNumBlockst), &
                            NXB*NYB*(NZB+1)*localNumBlockst, &
                    FLASH_REAL, jproc, &
                    9+NUNK_VARS+(NFACE_VARS*2)+i, &
                    io_globalComm, ierr)
                  end if
        end if
        end do
#endif
#endif

           else
              
              coord(:,1:localNumBlockst) = coordt(:,1:localNumBlockst)
              bsize(:,1:localNumBlockst) = bsizet(:,1:localNumBlockst)
              bnd_box(:,:,1:localNumBlockst) = bnd_boxt(:,:,1:localNumBlockst)
              lrefine(1:localNumBlockst) = lrefinet(1:localNumBlockst)
              nodetype(1:localNumBlockst) = nodetypet(1:localNumBlockst)
              gr_gid(:,1:localNumBlockst) = gidt(:,1:localNumBlockst)

#ifdef FLASH_GRID_PARAMESH3OR4
              which_child(1:localNumBlockst) = which_childt(1:localNumBlockst)
              bflags(:,1:localNumBlockst) = bflagst(:,1:localNumBlockst)
              if (do_gsurr_blks_read) then
                 gr_gsurr_blks(:,:,:,:,1:localNumBlockst) = &
                      gsurr_blkst(:,:,:,:,1:localNumBlockst)
              end if
#endif

              unk(:,io_ilo:io_ihi, io_jlo:io_jhi, io_klo:io_khi,1:localNumBlockst) & 
                   = unkBuf(:,1:NXB,1:NYB,1:NZB,1:localNumBlockst)
             
#if(NFACE_VARS > 0)
          facevarx(:,io_ilo:io_ihi+1, io_jlo:io_jhi, io_klo:io_khi,1:localNumBlockst) & 
                   = faceXBuf(:,1:NXB+1,1:NYB,1:NZB,1:localNumBlockst)
              if(NDIM .gt. 1) then
              facevary(:,io_ilo:io_ihi, io_jlo:io_jhi+1, io_klo:io_khi,1:localNumBlockst) & 
                   = faceYBuf(:,1:NXB,1:NYB+1,1:NZB,1:localNumBlockst)
              if(NDIM .gt. 2) then 
              facevarz(:,io_ilo:io_ihi, io_jlo:io_jhi, io_klo:io_khi+1,1:localNumBlockst) & 
                   = faceZBuf(:,1:NXB,1:NYB,1:NZB+1,1:localNumBlockst)

          endif
          endif
#endif

           end if !if(jproc /= MASTER_PE)
           
       deallocate(unkBuf)
#if(NFACE_VARS > 0) 
         deallocate(faceXBuf)
       if((NDIM .gt. 1) .AND. (NFACE_VARS .gt. 0)) deallocate(faceYBuf)
       if((NDIM .gt. 2) .AND. (NFACE_VARS .gt. 0)) deallocate(faceZBuf)
#endif
        end if !if(localNumBlockst > 0)
     end do !jproc= 0, io_globalNumProcs -1
     
  else !io_globalMe /= MASTER_PE
     
     call MPI_SEND(localNumBlocks, 1, FLASH_INTEGER, & 
          MASTER_PE, 1, io_globalComm, ierr)
     
     call MPI_SEND(gr_globalOffset, 1, FLASH_INTEGER, & 
          MASTER_PE, 2, io_globalComm, ierr)
     
     if (localNumBlocks > 0) then
        
        call MPI_RECV(coord(1,1), MDIM*localNumBlocks, &
             FLASH_REAL, MASTER_PE, 2, & 
             io_globalComm, status, ierr)
        
        call MPI_RECV(bsize(1,1), MDIM*localNumBlocks, & 
             FLASH_REAL, MASTER_PE, 3, & 
             io_globalComm, status, ierr)
        
        call MPI_RECV(bnd_box(1,1,1), 2*MDIM*localNumBlocks, & 
             FLASH_REAL, MASTER_PE, 4, & 
             io_globalComm, status, ierr)
        
        call MPI_RECV(lrefine(1), localNumBlocks, & 
             FLASH_INTEGER, MASTER_PE, 5, & 
             io_globalComm, status, ierr)
        
        call MPI_RECV(nodetype(1), localNumBlocks, & 
             FLASH_INTEGER, MASTER_PE, 6, & 
             io_globalComm, status, ierr)

        call MPI_RECV(gr_gid(1,1), localNumBlocks*(nfaces+1+nchild), & 
             FLASH_INTEGER, MASTER_PE, 7, & 
             io_globalComm, status, ierr)

#ifdef FLASH_GRID_PARAMESH3OR4
        call MPI_RECV(which_child(1), localNumBlocks, & 
             FLASH_INTEGER, MASTER_PE, 8, & 
             io_globalComm, status, ierr)
        
        call MPI_RECV(bflags(1,1), localNumBlocks*MFLAGS, & 
             FLASH_INTEGER, MASTER_PE, 9, & 
             io_globalComm, status, ierr)

        if (do_gsurr_blks_read) then
           call MPI_RECV(gr_gsurr_blks(1,1,1,1,1), &
                2*(1+(K1D*2))*(1+(K2D*2))*(1+(K3D*2))*localNumBlocks, &
                FLASH_INTEGER, MASTER_PE, 9999, & 
                io_globalComm, status, ierr)
        end if
#endif

        allocate(unkBuf(NUNK_VARS, NXB, NYB, NZB, localNumBlocks))
        do i = UNK_VARS_BEGIN,UNK_VARS_END
           call MPI_RECV(unkBuf(i,:,:,:,1:localNumBlocks), &
                NXB*NYB*NZB*localNumBlocks, & 
                FLASH_REAL, MASTER_PE, 9+i, & 
                io_globalComm, status, ierr)
           
        end do
        
        unk(:,io_ilo:io_ihi, io_jlo:io_jhi, io_klo:io_khi,1:localNumBlocks) = & 
             unkBuf(:,1:NXB,1:NYB,1:NZB,1:localNumBlocks)
        deallocate(unkBuf)

#ifdef FLASH_GRID_PARAMESH3OR4
    !Face-centered variables, receive end
#if(NFACE_VARS > 0)
      do i= 1, NFACE_VARS
      
        allocate(faceXBuf(NUNK_VARS, NXB+1, NYB, NZB, localNumBlocks))
           call MPI_RECV(faceXBuf(i,:,:,:,1:localNumBlocks), &
                (NXB+1)*NYB*NZB*localNumBlocks, & 
                FLASH_REAL, MASTER_PE, &
        9+i+NUNK_VARS, & 
                io_globalComm, status, ierr)
         facevarx(i,io_ilo:io_ihi+1, io_jlo:io_jhi, io_klo:io_khi,1:localNumBlocks) = & 
         faceXBuf(i,1:NXB+1,1:NYB,1:NZB,1:localNumBlocks)
     deallocate(faceXBuf)

       if(NDIM .gt. 1) then
       
        allocate(faceYBuf(NUNK_VARS, NXB, NYB+1, NZB, localNumBlocks))
           call MPI_RECV(faceYBuf(i,:,:,:,1:localNumBlocks), &
                NXB*(NYB+1)*NZB*localNumBlocks, & 
                FLASH_REAL, MASTER_PE, &
        9+i+NUNK_VARS+NFACE_VARS, & 
                io_globalComm, status, ierr)
    facevary(i,io_ilo:io_ihi, io_jlo:io_jhi+1, io_klo:io_khi,1:localNumBlocks) = & 
         faceYBuf(i,1:NXB,1:NYB+1,1:NZB,1:localNumBlocks)
     deallocate(faceYBuf)

     if(NDIM .gt. 2) then

           allocate(faceZBuf(NUNK_VARS, NXB, NYB, NZB+1, localNumBlocks))
           call MPI_RECV(faceZBuf(i,:,:,:,1:localNumBlocks), &
                NXB*NYB*(NZB+1)*localNumBlocks, & 
                FLASH_REAL, MASTER_PE, &
        9+i+NUNK_VARS+(NFACE_VARS*2), & 
                io_globalComm, status, ierr)
           facevarz(i,io_ilo:io_ihi, io_jlo:io_jhi, io_klo:io_khi+1,1:localNumBlocks) = & 
           faceZBuf(i,1:NXB,1:NYB,1:NZB+1,1:localNumBlocks)
       deallocate(faceZBuf)
     
     end if
     
       end if

     end do

#endif
#endif
    
        
     end if !localNumBlocks > 0
     
  end if !io_globalMe == MASTER_PE
  
  !! construct the tree from the gid data
  ! neighbor data


  alnblocks = int(gr_globalNumBlocks/io_globalNumProcs) + 1

  yy = (io_globalNumProcs*alnblocks) - gr_globalNumBlocks
  xx = io_globalNumProcs - yy
 

#ifdef FLASH_GRID_PARAMESH3OR4
  if (do_gsurr_blks_read) then
     call MPI_Bcast (gr_is_gsurr_blks_initialized, 1, FLASH_LOGICAL, MASTER_PE, & 
          io_globalComm, ierr)
  else
     gr_is_gsurr_blks_initialized = .false.
  end if
#endif
  !Extract data from gr_gid and gr_gsurr_blks arrays
  call Grid_receiveInputData(localNumBlocks, alnblocks, xx)


  call MPI_BARRIER(io_globalComm,ierr)

  return  
end subroutine io_readData
