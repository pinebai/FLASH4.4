!!****if* source/IO/IOMain/direct/UG/io_writeData
!!
!! NAME
!!
!!  io_writeData
!!
!!
!! SYNOPSIS
!!
!!  io_writeData(integer(in) :: ioLun)
!!               
!!
!!
!!
!! DESCRIPTION
!!
!!  This routine implements a bare bones, basic, non parallel fortran IO.
!!  Each processor writes its data to a unique file.  If you are running
!!  on 64,0000 procs, you will get 64,000 files for each call to 
!!  IO_writeCheckpoint.
!!
!!  In this method, no 'meta data', (like data/time of run, compiler flags etc)
!!  are written.  Only simulation data like (time, dt nstep etc) and the grid data
!!  boundBox, blockCenterCoords, neighbor data (gid),
!!  block physical size and unk data are written to the file.
!!
!!  The order in which the data is written is set up to match the corresponding
!!  io_readData routine.  If you add additional data to output, you will need   
!!  to make corresponding changes in io_readData.
!! 
!!  
!!
!! ARGUMENTS
!! 
!!  ioLun - logical unit number for fortran write file
!!
!! NOTES
!!
!!  Using this IO implementation will be a pain for anyone doing post
!!  processing.  We highly recommend that you use hdf5 or parallel netcdf
!!  if it is at all possible.
!!
!!  You will only be able to restart using the same number of processors!
!!
!!***




subroutine io_writeData(ioLun)  

  use IO_data, ONLY : io_globalMe, io_globalNumProcs,  io_realParmNames, io_realParmValues, io_numRealParms, &
       io_intParmNames, io_intParmValues, io_numIntParms, &
       io_logParmNames, io_logParmValues, io_numLogParms, &
       io_strParmNames, io_strParmValues, io_numStrParms, &
       io_realScalarNames, io_realScalarValues, io_numRealScalars, &
       io_intScalarNames, io_intScalarValues, io_numIntScalars, &
       io_logScalarNames, io_logScalarValues, io_numLogScalars, &
       io_strScalarNames, io_strScalarValues, io_numStrScalars, &
       io_logToIntScalarValues, io_logToIntParmValues, io_unklabels, &
       io_ilo, io_ihi, io_jlo, io_jhi, io_klo, io_khi, io_geometry, &
       io_setupCall, io_buildDir, io_flashRelease, &
       io_fileCreationTime, io_buildDate, io_buildMachine, io_cflags, io_fflags, &
       io_setupTimeStamp, io_buildTimeStamp, io_doublePrecision, &
       io_unklabels, io_ilo, io_ihi, io_jlo, io_jhi, io_klo, io_khi, &
       io_plotVarStr, io_nPlotVars
  use IO_interface, ONLY : IO_getScalar
  use Grid_interface, ONLY : Grid_getLocalNumBlks, &
    Grid_getBlkBoundBox, Grid_getBlkCenterCoords, &
    Grid_getBlkPhysicalSize

  
  use Grid_data, ONLY : gr_gid, scratch

  use physicalData, only : unk

  implicit none
  
  include "mpif.h"
#include "constants.h"
#include "Flash.h"

  integer, intent(in) :: ioLun

  integer :: localNumBlocks
  integer :: i,j,k,var, stat, globalNumBlocks, globalOffset
  logical :: isPlotVar

  !get these values from scalar lists
  integer :: localnstep
  real    :: localdt, localtime

  real :: boundBox(2, MDIM)
  real :: blockCenterCoords(MDIM) !!only single dim because only 1 block per proc
  real :: bsize(MDIM) !!only single dime because only 1 block per proc


  ! allocate storage to hold a single variable information
  ! this should only be a small memory overhead
  integer, parameter :: single = SELECTED_REAL_KIND(p=6)
  real (kind=single) :: unkt(1,NXB,NYB,NZB,1)


  ! localNumBlocks should be 1
  call Grid_getLocalNumBlks(localNumBlocks)


  !if we are writing a checkpoint file
  if(io_doublePrecision) then
     !! call the generic function prepareLists to allocate and 
     !! fill the runtime parameter lists and the scalar lists
     call io_prepareListsWrite()
     
     !first write the number of parameters and scalars
     !(may seem like a strange order, but necessary to know
     !number of parms and scalars first in order to restart)

     write(ioLun)io_numRealParms
     write(ioLun)io_numIntParms
     write(ioLun)io_numStrParms
     write(ioLun)io_numLogParms
     
     write(ioLun)io_numRealScalars
     write(ioLun)io_numIntScalars
     write(ioLun)io_numStrScalars
     write(ioLun)io_numLogScalars
     

     !! write the runtime parameters  
     write(ioLun)io_realParmNames
     write(ioLun)io_realParmValues
     write(ioLun)io_intParmNames
     write(ioLun)io_intParmValues
     write(ioLun)io_strParmNames
     write(ioLun)io_strParmValues
     write(ioLun)io_logParmNames
     write(ioLun)io_logToIntParmValues
     
     !! write the scalars  
     write(ioLun)io_realScalarNames
     write(ioLun)io_realScalarValues
     write(ioLun)io_intScalarNames
     write(ioLun)io_intScalarValues
     write(ioLun)io_strScalarNames
     write(ioLun)io_strScalarValues
     write(ioLun)io_logScalarNames
     write(ioLun)io_logToIntScalarValues
     
     
     call io_finalizeListsWrite()

  else !we are writing a plotfile
     
     call IO_getScalar("time", localtime)
     call IO_getScalar("dt", localdt)
     call IO_getScalar("nstep", localnstep)

     write(ioLun)localtime
     write(ioLun)localdt
     write(ioLun)localnstep

  end if


  
    
  !! write gid. gr_gid is filled in Grid_sendOutputData
  !! gid contains neighbor data
  write(ioLun)gr_gid


  !! Write the bndbox
  !! For the uniform grid there is always one block per proc so
  !! the first entry, block_id in Grid_getBlkBoundBox is always 1
  call Grid_getBlkBoundBox(1, boundBox)
  write(ioLun)boundBox
  


  !store the coordinates
  call Grid_getBlkCenterCoords(1, blockCenterCoords)
  write(ioLun)blockCenterCoords


  ! store the block size
  call Grid_getBlkPhysicalSize(1, bsize)
  write(ioLun)bsize


  !--------------------------------------------------------------------------
  ! store the unknowns -- 
  !--------------------------------------------------------------------------

  do var = UNK_VARS_BEGIN,UNK_VARS_END

     !if we are writing a checkpoint file
     if(io_doublePrecision) then
                
        write(ioLun)unk(var,io_ilo:io_ihi, io_jlo:io_jhi, io_klo:io_khi,1) 

     else !we are writing a plotfile

        !check to see if variable should be written to the plotfile
        call io_isPlotVar(var, isPlotVar, MAPBLOCK_UNK)

        if (isPlotVar) then

           unkt(1,:,:,:,1) = real( &
                (unk(var,io_ilo:io_ihi,io_jlo:io_jhi,io_klo:io_khi, 1)), kind = single) 

           
           write(ioLun)unkt               
           
        end if
     endif
     
  end do

  

  !write the scratch grid vars if the user defines any in flash.par
  !we can use the same routine as when writing the unknowns.
  do i = SCRATCH_GRID_VARS_BEGIN,SCRATCH_GRID_VARS_END
     
     call io_isPlotVar(i, isPlotVar, MAPBLOCK_SCRATCH)
     if(isPlotVar) then
        
        if(io_doublePrecision) then
           
           write(ioLun)scratch(i,io_ilo:io_ihi, io_jlo:io_jhi, io_klo:io_khi,1) 

        else
           
           unkt(1,1:NXB,1:NYB,1:NZB,1) = real(scratch(i,io_ilo:io_ihi, io_jlo:io_jhi, io_klo:io_khi,1), kind = single) 
           
           write(ioLun)unkt
           
        end if
     end if
  end do




end subroutine io_writeData
