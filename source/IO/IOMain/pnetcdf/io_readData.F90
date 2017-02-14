!!****if* source/IO/IOMain/pnetcdf/io_readData
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
!!
!! DESCRIPTION
!!
!!  This is the reading counterpart to io_writeData.  It reads parallel
!!  netcdf 
!!  file and distributes it to the processors to restart a simulation.
!!
!!
!! ARGUMENTS
!! 
!!
!!***

#include "constants.h"
#include "Flash.h"
#include "io_flash.h"

subroutine io_readData()

#ifdef USE_IO_C_INTERFACE
  use iso_c_binding, ONLY : c_loc, c_char
  use io_c_interface, ONLY : io_xfer_cont_slab, io_ncmpi_read_file_format
#endif

  use Grid_data, ONLY : gr_nToLeft, gr_globalOffset, &
       gr_globalNumBlocks, gr_gid
  use Driver_interface, ONLY : Driver_abortFlash
  use Logfile_interface, ONLY : Logfile_stamp
  use Grid_interface, ONLY : Grid_putLocalNumBlks, Grid_receiveInputData

  use IO_data, ONLY : io_globalMe, io_globalComm,&
        io_realParmNamesPrev, io_realParmValuesPrev, io_numRealParmsPrev, &
       io_intParmNamesPrev, io_intParmValuesPrev, io_numIntParmsPrev, &
       io_logParmNamesPrev, io_logParmValuesPrev, io_numLogParmsPrev, &
       io_strParmNamesPrev, io_strParmValuesPrev, io_numStrParmsPrev, &
       io_realScalarNames, io_realScalarValues, io_numRealScalars, &
       io_intScalarNames, io_intScalarValues, io_numIntScalars, &
       io_logScalarNames, io_logScalarValues, io_numLogScalars, &
       io_strScalarNames, io_strScalarValues, io_numStrScalars, &
       io_logToIntScalarValues, io_logToIntParmValuesPrev, &
       io_ilo, io_ihi, io_jlo, io_jhi, io_klo, io_khi, io_unklabels, &
       io_checkpointFileNumber, io_baseName, io_chkptFileID, &
       io_globalMe, io_meshMe, io_meshNumProcs, tree_data_t
  use IO_interface, ONLY : IO_getScalar

#ifdef FLASH_GRID_PARAMESH
  use tree, ONLY : nodetype, bnd_box, coord, bsize, lrefine, &
       neigh, child, parent, nfaces, nchild
#ifdef FLASH_GRID_PARAMESH3OR4
  use tree, ONLY : which_child, bflags
  use Grid_data, ONLY : gr_gsurr_blks
#endif
#endif

  use io_typeInterface, ONLY : io_xfer_mesh_data, io_xfer_tree_data
  use physicaldata, ONLY : unk, facevarx, facevary, facevarz

  implicit none
#ifndef USE_IO_C_INTERFACE
#define c_loc(x) x
  integer, parameter :: c_char = KIND('A')
#endif

#include "Flash_mpi.h"

  type(tree_data_t) :: tree_data
  integer, parameter :: libType = IO_FILE_PNETCDF, &
       fileType = CHECKPOINTFILE
  integer :: blockID, localNumBlocks, procBlocks, ngid
  integer :: i, lb, stat, ierr, varid, j, xx, yy, alnblocks, k

  integer, allocatable :: procnumber(:)

  integer :: alocalNumBlocks
  
  character (len=4) :: fnumString
  character (len=MAX_STRING_LENGTH) :: filename
  character(len=MAX_STRING_LENGTH), save, allocatable, dimension(:,:) :: strBuff

  integer, parameter :: xferType = IO_READ_XFER
  integer :: fileID, fileFmt, presentDims

#ifndef FLASH_GRID_PARAMESH
  !Declare arrays that exist only in Paramesh simulations.
  integer, parameter :: nfaces = 2*NDIM, nchild = 2**NDIM
  integer, target, dimension(1) :: nodetype = (/LEAF/), lrefine = (/1/)
  real, target, dimension(2,MDIM,1) :: bnd_box
  real, target, dimension(MDIM,1) :: coord, bsize
#endif

  fileID = io_chkptFileID


  !create the filename
  write (fnumString, '(i4.4)') io_checkpointFileNumber
  filename = trim(io_baseName) // 'ncmpi_chk_'// fnumString
  


  call io_ncmpi_open_file_for_read(io_chkptFileID, filename)

  if (io_globalMe == MASTER_PE) then
     allocate (strBuff(2,2))
     print *, 'file: ', trim(filename), ' opened for restart'
     write (strBuff(1,1), "(A)") "type"
     write (strBuff(1,2), "(A)") "checkpoint"
     write (strBuff(2,1), "(A)") "name"
     write (strBuff(2,2), "(A)") trim(filename)
     call Logfile_stamp( strBuff, 2, 2, "[io_readData] file opened")
     deallocate(strBuff)
  end if
  
  
  

  call io_prepareListsRead()
  
  call io_ncmpi_read_header(io_globalMe, &
       io_chkptFileID, &
       NUNK_VARS, &
       io_unklabels, &
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
       io_logToIntParmValuesPrev, &
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
  call io_ncmpi_read_file_format(io_globalMe, io_chkptFileID, fileFmt)
  
  call io_finalizeListsRead()
  
  call IO_getScalar("globalnumblocks", gr_globalNumBlocks)
  
  
  
  
  !---------------------------------------------------------------------------
  ! compute the number of blocks on each processor -- this will be used to
  ! get the offset into the file for the parallel read
  !---------------------------------------------------------------------------

#ifdef FLASH_GRID_PARAMESH
  ! compute the approximate number of blocks per processor
  alnblocks = int(gr_globalNumBlocks/io_meshNumProcs) + 1
  
  ! check for error -- if the number of blocks we want to put on each
  ! processor is greater than maxblocks, then abort
  if (alnblocks .GT. MAXBLOCKS) then
     
     print *
     print *, '********** ERROR in READ_DATA ************'
     print *
     print *,' Number of blocks per processor exceeds maxblocks.'
     print *,' Suggest you reset maxblocks to a larger number or'
     print *,' run on a larger number of processors.'
     print *,' globalNumBlocks, numProcs = ', gr_globalNumBlocks, io_meshNumProcs
     print *
     
     call Driver_abortFlash('[READ_DATA] ERROR: num blocks per proc exceeds maxblocks')
       
  end if
    
  
  ! figure out the excess blocks
  yy = (io_meshNumProcs*alnblocks) - gr_globalNumBlocks
  xx = io_meshNumProcs - yy
  
  ! loop over all the processor numbers and figure out how many blocks are
  ! stored to the left of the processor -- this is a little tricky
  
  gr_nToLeft(0) = 0
  do i = 0, io_meshNumProcs - 2
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
  if (io_meshMe < xx) then
     localNumBlocks = alnblocks
  else
     localNumBlocks = alnblocks - 1
  endif
  
  if (alnblocks .EQ. 0) then
     if (io_meshMe < gr_globalNumBlocks) then
        localNumBlocks = 1
     else
        localNumBlocks = 0
     end if
  end if
  
  ! compute the offset into the dataspace in the HDF5 file
  gr_globalOffset = gr_nToLeft(io_meshMe)

#else
  localNumBlocks = 1
#endif

  call Grid_putLocalNumBlks(localNumBlocks)
  call io_checkBlockShape(localNumBlocks)


#ifdef FLASH_GRID_PARAMESH
  tree_data % bnd_box => bnd_box
  tree_data % coord => coord
  tree_data % bsize => bsize
  tree_data % gid => gr_gid
  tree_data % nodetype => nodetype
  tree_data % lrefine => lrefine
# ifdef FLASH_GRID_PARAMESH3OR4
  tree_data % bflags => bflags
  tree_data % which_child => which_child
# else
  nullify(tree_data % bflags)
  nullify(tree_data % which_child)
# endif
  nullify(tree_data % procnumber)
# ifdef FLASH_GRID_PARAMESH4DEV_SURR_BLKS_OPTION
  tree_data % gsurr_blks => gr_gsurr_blks
# else
  nullify(tree_data % gsurr_blks)
# endif

  call io_ncmpi_read_present_dims(fileID, presentDims)

  call io_xfer_tree_data(tree_data, &
       io_chkptFileID, libType, xferType, &
       localNumBlocks, gr_globalOffset, presentDims)

  !Extract data from gr_gid and gr_gsurr_blks arrays
  call Grid_receiveInputData(localNumBlocks, alnblocks, xx)
#endif


  call io_xfer_mesh_data(fileID, fileFmt, fileType, &
       libType, xferType, localNumBlocks, gr_globalOffset)



  if (io_globalMe == MASTER_PE) then
     allocate (strBuff(2, 2))
     write (strBuff(1,1), "(A)") "type"
     write (strBuff(1,2), "(A)") "checkpoint"
     write (strBuff(2,1), "(A)") "name"
     write (strBuff(2,2), "(A)") trim(filename)
     call Logfile_stamp( strBuff, 2, 2, "file_rd_close")
     deallocate(strBuff)
  end if


  call mpi_barrier (io_globalComm, ierr)
  if (io_globalMe == MASTER_PE) &
    print *, 'io_readData:  finished reading input file.'

  return
end subroutine io_readData
