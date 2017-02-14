!!****if* source/IO/IOMain/hdf5/parallel/PM_argonne/io_readData
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
!!  All reading is done using parallel HDF5 -- no explicit data movement
!!  is used.
!!
!!
!!***

#ifdef DEBUG_ALL
#define DEBUG_IO
#endif

#include "constants.h"
#include "Flash.h"
#include "io_flash.h"

subroutine io_readData()

#ifdef USE_IO_C_INTERFACE
  use iso_c_binding, ONLY : c_loc, c_char
  use io_c_type_interface, ONLY : io_h5_read_file_format
#endif

  use Grid_data, ONLY : gr_globalNumBlocks, gr_nToLeft, &
       gr_globalOffset, gr_gid
  use Driver_interface, ONLY : Driver_abortFlash
  use Logfile_interface, ONLY : Logfile_stamp
  use Grid_interface, ONLY : Grid_putLocalNumBlks, Grid_getBlkIndexLimits, &
       Grid_receiveInputData

  use IO_data, ONLY : io_globalMe, io_globalNumProcs,io_globalComm,&
        io_baseName, io_checkpointFileNumber, &
       io_realParmNames, io_realParmValues, io_numRealParms, &
       io_intParmNames, io_intParmValues, io_numIntParms, &
       io_logParmNames, io_logParmValues, io_numLogParms, &
       io_strParmNames, io_strParmValues, io_numStrParms, &
       io_realScalarNames, io_realScalarValues, io_numRealScalars, &
       io_intScalarNames, io_intScalarValues, io_numIntScalars, &
       io_logScalarNames, io_logScalarValues, io_numLogScalars, &
       io_strScalarNames, io_strScalarValues, io_numStrScalars, &
       io_logToIntScalarValues, io_logToIntParmValues, io_unklabels, &
       io_ilo, io_ihi, io_jlo, io_jhi, io_klo, io_khi, &
       io_outputSplitNum, io_comm, io_chkptFileID, &
       io_faceXVarLabels, io_faceYVarLabels, io_faceZVarLabels,&
       io_realParmNamesPrev, io_realParmValuesPrev, io_numRealParmsPrev, &
       io_intParmNamesPrev, io_intParmValuesPrev, io_numIntParmsPrev, &
       io_logParmNamesPrev, io_logParmValuesPrev, io_numLogParmsPrev, &
       io_strParmNamesPrev, io_strParmValuesPrev, io_numStrParmsPrev, &
       io_logToIntParmValuesPrev, io_splitNumBlks, io_globalNumProcs, &
       io_meshMe, io_meshNumProcs, tree_data_t
  use IO_interface, ONLY : IO_getScalar
  use io_typeInterface, ONLY : io_xfer_mesh_data, io_xfer_tree_data
#ifdef FLASH_GRID_PARAMESH
  use tree, ONLY : nodetype, bnd_box, coord, bsize, lrefine, &
       neigh, child, parent, nchild, nfaces
#ifdef FLASH_GRID_PARAMESH3OR4
  use tree, ONLY : which_child, bflags
  use Grid_data, ONLY : gr_gsurr_blks
#endif
#endif

  implicit none
  type(tree_data_t) :: tree_data
#ifndef USE_IO_C_INTERFACE
#define c_loc(x) x
  integer, parameter :: c_char = KIND('A')
#endif

#include "Flash_mpi.h"


  integer :: localNumBlocks, ngid
  integer :: alocalNumBlocks

 
  character (len=4) :: fnumStr
  character (len=MAX_STRING_LENGTH) :: filename
  character(len=MAX_STRING_LENGTH), save, allocatable, dimension(:,:) :: strBuff


  integer :: blockID, procBlocks, ierr

  integer :: i, lb, j, xx, yy, alnblocks

  integer, allocatable :: procnumber(:) !dimension(localNumBlocks)
  
  integer :: realGlobalNumBlocks  
  integer :: blkLimits(2,MDIM), blkLimitsGC(2, MDIM)

  integer :: splitOffset, localOffset


  !presentDims needs to depend on FILE_FORMAT_VERSION.
  integer, parameter :: xferType = IO_READ_XFER
  integer :: fileFmt, presentDims
  integer, parameter :: libType = IO_FILE_HDF5, &
       fileType = CHECKPOINTFILE

#ifndef FLASH_GRID_PARAMESH
  !Declare arrays that exist only in Paramesh simulations.
  integer, parameter :: nfaces = 2*NDIM, nchild = 2**NDIM
  integer, target, dimension(1) :: nodetype = (/LEAF/), lrefine = (/1/)
  real, target, dimension(2,MDIM,1) :: bnd_box
  real, target, dimension(MDIM,1) :: coord, bsize
#endif



  call io_getOutputName(io_checkpointFileNumber, "hdf5", "_chk_", filename, .false.)


  call io_h5open_file_for_read(io_chkptFileID, filename, io_comm, io_outputSplitNum)

  if (io_globalMe == MASTER_PE) then
       allocate (strBuff(2,2))
       print *, 'file: ', trim(filename), ' opened for restart'
       write (strBuff(1,1), "(A)") "type"
       write (strBuff(1,2), "(A)") "checkpoint"
       write (strBuff(2,1), "(A)") "name"
       write (strBuff(2,2), "(A)") trim(filename)
       call Logfile_stamp( strBuff, 2, 2, "[io_readData]")
  end if
  if (allocated(strBuff)) deallocate(strBuff)


 !!read in the total number of blocks, time, and timestep
  call io_h5read_header(io_globalMe, &
                        io_chkptFileID, &
                        io_unklabels, &
                        io_outputSplitNum)



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


  call io_finalizeListsRead()
  call io_h5_read_file_format(io_globalMe, io_chkptFileID, fileFmt)

  call IO_getScalar("globalNumBlocks", gr_globalNumBlocks)
#ifndef IO_FLASH_NOFBS_UG
  call io_checkBlockShape(gr_globalNumBlocks)
#endif
  if (io_outputSplitNum /= 1) then
     call IO_getScalar("splitNumBlocks", io_splitNumBlks)
  else
     io_splitNumBlks = gr_globalNumBlocks
  endif


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
     
     call Driver_abortFlash('[io_readData] ERROR: num blocks per proc exceeds maxblocks')
     
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

  !find our offset into a potentially split file:
  if(io_outputSplitNum > 1) then
     call MPI_ALLREDUCE(gr_globalOffset, splitOffset,1,FLASH_INTEGER,&
                        MPI_MIN, io_comm, ierr)
     localOffset = gr_globalOffset - splitOffset
     
  else
     localOffset = gr_globalOffset
  end if


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

  call io_h5read_present_dims(io_chkptFileID, presentDims)

  call io_xfer_tree_data(tree_data, &
       io_chkptFileID, libType, xferType, &
       localNumBlocks, localOffset, presentDims)

  !Extract data from gr_gid and gr_gsurr_blks arrays
  call Grid_receiveInputData(localNumBlocks, alnblocks, xx)
#endif


  call io_xfer_mesh_data(io_chkptFileID, fileFmt, fileType, &
       libType, xferType, localNumBlocks, localOffset)


  if (io_globalMe == MASTER_PE) &
    print *, 'read_data:  read ', gr_globalNumBlocks, ' blocks.'

  if (io_globalMe == MASTER_PE) then
     allocate (strBuff(2, 2))
     write (strBuff(1,1), "(A)") "type"
     write (strBuff(1,2), "(A)") "checkpoint"
     write (strBuff(2,1), "(A)") "name"
     write (strBuff(2,2), "(A)") trim(filename)
     call Logfile_stamp( strBuff, 2, 2, "[io_readData] file_closed")
     if (allocated(strBuff)) deallocate(strBuff)
  end if
 
  call MPI_BARRIER (io_globalComm, ierr)
  if (io_globalMe == MASTER_PE) &
    print *, 'io_readData:  finished reading input file.'

  return
end subroutine io_readData
