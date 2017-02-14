!!****if* source/IO/IOMain/pnetcdf/explicitCopy/io_xfer_mesh_data
!!
!! NAME
!!  io_xfer_mesh_data
!!
!! SYNOPSIS
!!
!!  io_xfer_mesh_data(integer(in) :: fileID,
!!                    integer(in) :: fileFmt,
!!                    integer(in) :: fileType,
!!                    integer(in) :: libType,
!!                    integer(in) :: xferType,
!!                    integer(in) :: localNumBlocks,
!!                    integer(in) :: globalBlockOffset)
!!
!!
!! DESCRIPTION
!!
!! This subroutine calls a C function to transfer mesh data
!! between memory and file.  It works by setting a pointer to 
!! UNK, FACEX, FACEY, FACEZ and SCRATCH arrays, and then passing the
!! address of the first grid element to the C transfer function.
!! It can transfer data from memory to file and from file to memory.
!!
!!
!! ARGUMENTS
!!
!! fileID: The HDF5 or pnetcdf file identifier (used directly by the libraries)
!! fileFmt: The FLASH file layout.  The standard layout that tools
!!          understand is 9 (one mesh variable per dataset), but there is also
!!          support for an experimental file layout 10 (all mesh variables in
!!          the same dataset)
!! fileType: The FLASH file type (checkpoint file or plot file)
!! libType: The library we are using (HDF5 or pnetcdf)
!! xferType: The direction of data transfer: memory to file or file to memory
!! localNumBlocks: The number of blocks on myPE being trasferred to/from file
!! globalBlockOffset: The read/write block offset in file.
!!
!!***

!!REORDER(5): unk, facevar[xyz], unkBuf, face[XYZ]Buf

#include "constants.h"
#include "Flash.h"
#include "io_flash.h"

subroutine io_xfer_mesh_data(fileID, fileFmt, fileType, &
     libType, xferType, localNumBlocks, globalBlockOffset)

  use Grid_data, ONLY : gr_globalNumBlocks, scratch
  use Driver_interface, ONLY : Driver_abortFlash
  use IO_data, ONLY : io_globalMe, io_globalNumProcs, io_globalComm,&
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
       io_type_matched_xfer, io_doublePrecision, io_plotfileGridQuantityDP, &
       io_plotGridVarStr, io_faceXVarLabels, io_faceYVarLabels, io_faceZVarLabels, &
       io_unkLabelsGlobal
  use physicaldata, ONLY : unk, facevarx, facevary, facevarz
#ifdef USE_IO_C_INTERFACE
  use io_c_interface, ONLY : io_ncmpi_name_to_id
#endif

  implicit none
  include "Flash_mpi.h"

  integer, intent(IN) :: fileID, fileFmt, fileType, libType, &
       xferType, localNumBlocks, globalBlockOffset
  integer :: i, varid, j, k
  integer :: localOffset

  real,allocatable :: unkBuf(:,:,:,:,:)
  real,allocatable :: faceXBuf(:,:,:,:,:)
  real,allocatable :: faceYBuf(:,:,:,:,:)
  real,allocatable :: faceZBuf(:,:,:,:,:)
  real, allocatable :: globalVarMin(:), globalVarMax(:)
  real, allocatable :: globalFaceXMin(:), globalFaceXMax(:)
  real, allocatable :: globalFaceYMin(:), globalFaceYMax(:)
  real, allocatable :: globalFaceZMin(:), globalFaceZMax(:)

  integer, parameter :: single = SELECTED_REAL_KIND(p=6)
  real (kind=single), allocatable :: unkt(:,:,:,:,:)
  real (kind=single) :: spMax, spMin
  logical :: includeVar
  integer :: idx
  character(len=MAX_STRING_LENGTH) :: s


  !Set the local offset to 0 when a process will write no data.
  !This stops an (overly cautious) index out of bounds assertion
  !error in pnetcdf library (version 1.2.0pre1).
  localOffset = globalBlockOffset
  if (localOffset >= gr_globalNumBlocks) then
     localOffset = 0
  end if


  allocate(unkBuf(1, NXB, NYB, NZB, 1:MAXBLOCKS))


  if (xferType == IO_WRITE_XFER) then

     allocate(globalVarMin(NUNK_VARS))
     allocate(globalVarMax(NUNK_VARS))

     !get the max and minimum variables
     call io_getVarExtrema(NUNK_VARS, globalVarMin, globalVarMax, CENTER)
     allocate(unkt(1, NXB, NYB, NZB, MAXBLOCKS))


     idx = 0
     do i = UNK_VARS_BEGIN,UNK_VARS_END

        idx = idx + 1 !Index to extract current mesh variable string.

        if(io_doublePrecision .or. io_plotfileGridQuantityDP) then

           if(.not. io_doublePrecision) then
              call io_isPlotVar(i, includeVar,MAPBLOCK_UNK)
              if(.not. includeVar) cycle
           end if

           unkBuf(1,1:NXB,1:NYB,1:NZB,1:localNumBlocks) = &
                unk(i,io_ilo:io_ihi, io_jlo:io_jhi, io_klo:io_khi,1:localNumBlocks)


           s = io_unkLabelsGlobal(idx)
           call io_ncmpi_name_to_id(fileID, s, len_trim(s), varid)
           call io_ncmpi_write_unknowns(fileID, & 
                varid, &
                NXB, & 
                NYB, & 
                NZB, & 
                unkBuf, & 
                globalVarMin(i), &
                globalVarMax(i), &
                localNumBlocks, &
                gr_globalNumBlocks,  & 
                localOffset)
        else

           call io_isPlotVar(i, includeVar,MAPBLOCK_UNK)

           if (includeVar) then

              unkt(1,1:NXB,1:NYB,1:NZB,1:localNumBlocks) = &
                   real(unk(i,io_ilo:io_ihi, io_jlo:io_jhi, io_klo:io_khi,1:localNumBlocks), kind = single)


              spMin = real(globalVarMin(i), kind = single)
              spMax = real(globalVarMax(i), kind = single)

              s = io_unkLabelsGlobal(idx)
              call io_ncmpi_name_to_id(fileID, s, len_trim(s), varid)
              call io_ncmpi_write_unknowns_sp(fileID, & 
                   varid, &
                   NXB, & 
                   NYB, & 
                   NZB, & 
                   unkt, & 
                   spMin, &
                   spMax, &
                   localNumBlocks, &
                   gr_globalNumBlocks,  & 
                   localOffset)

           end if
        end if

     end do


     deallocate(unkBuf)
     deallocate(unkt)
     allocate(unkBuf(NXB,NYB,NZB,1,1))
     allocate(unkt(NXB,NYB,NZB,1,1))

     deallocate(globalVarMin)
     deallocate(globalVarMax)

     allocate(globalVarMin(NSCRATCH_GRID_VARS))
     allocate(globalVarMax(NSCRATCH_GRID_VARS))

     !get the max and minimum variables
     call io_getVarExtrema(NSCRATCH_GRID_VARS, globalVarMin, globalVarMax, SCRATCH)



     !write the scratch grid vars if the user defines any in flash.par
     !we can use the same routine as when writing the unknowns.
     do i = SCRATCH_GRID_VARS_BEGIN,SCRATCH_GRID_VARS_END

        call io_isPlotVar(i, includeVar, MAPBLOCK_SCRATCH)
        if(includeVar) then

           if(io_doublePrecision .or. io_plotfileGridQuantityDP) then

              unkBuf(1:NXB,1:NYB,1:NZB,1,1) = scratch(i,io_ilo:io_ihi, io_jlo:io_jhi, io_klo:io_khi,1) 


              call io_ncmpi_write_scratch(fileID, & 
                   NXB, & 
                   NYB, & 
                   NZB, & 
                   unkBuf, & 
                   globalVarMin(i), &
                   globalVarMax(i), &
                   io_plotGridVarStr(i), &
                   localNumBlocks, &
                   gr_globalNumBlocks,  & 
                   localOffset)

           else

              unkt(1:NXB,1:NYB,1:NZB,1,1) = real(scratch(i,io_ilo:io_ihi, io_jlo:io_jhi, io_klo:io_khi,1), kind = single) 

              spMax = real(globalVarMax(i), kind = single)
              spMin = real(globalVarMin(i), kind = single)


              call io_ncmpi_write_scratch_sp(fileID, & 
                   NXB, & 
                   NYB, & 
                   NZB, & 
                   unkt, & 
                   spMin, &
                   spMax, &
                   io_plotGridVarStr(i), &
                   localNumBlocks, &
                   gr_globalNumBlocks,  & 
                   localOffset)

           end if
        end if
     end do

     deallocate(unkBuf)
     deallocate(unkt)

     deallocate(globalVarMin)
     deallocate(globalVarMax)

     !start handling facevar output
#if(NFACE_VARS > 0)
     allocate(globalFaceXMin(NFACE_VARS))
     allocate(globalFaceXMax(NFACE_VARS))
     call io_getVarExtrema(NFACE_VARS, globalFaceXMin, globalFaceXMax, FACEX)

     if(NDIM .gt. 1) then
        allocate(globalFaceYMin(NFACE_VARS))
        allocate(globalFaceYMax(NFACE_VARS))
        call io_getVarExtrema(NFACE_VARS, globalFaceYMin, globalFaceYMax, FACEY)
     end if

     if(NDIM .gt. 2) then
        allocate(globalFaceZMin(NFACE_VARS))
        allocate(globalFaceZMax(NFACE_VARS))
        call io_getVarExtrema(NFACE_VARS, globalFaceZMin, globalFaceZMax, FACEZ)
     end if

     allocate(faceXBuf(1,NXB+1,NYB,NZB,1:localNumBlocks))
     if(NDIM .gt. 1) allocate(faceYBuf(1,NXB,NYB+1,NZB,1:localNumBlocks))
     if(NDIM .gt. 2) allocate(faceZBuf(1,NXB,NYB,NZB+1,1:localNumBlocks))


     do i = 1,NFACE_VARS
        if(io_doublePrecision) then !only for checkpointing
           faceXBuf(1,1:NXB+1,1:NYB,1:NZB,1:localNumBlocks) = &
                facevarx(i,io_ilo:io_ihi+1, io_jlo:io_jhi, io_klo:io_khi,1:localNumBlocks) 

           varid = 0 !Unused

           call io_ncmpi_write_facevar(fileID, & 
                varid, &
                IAXIS, &
                NXB+1, & 
                NYB, & 
                NZB, & 
                faceXBuf, & 
                globalFaceXMin(i), &
                globalFaceXMax(i), &
                io_faceXVarLabels(i), &
                localNumBlocks, &
                gr_globalNumBlocks,  & 
                localOffset)

           if (NDIM .gt. 1) then

              faceYBuf(1,1:NXB,1:NYB+1,1:NZB,1:localNumBlocks) = &
                   facevary(i,io_ilo:io_ihi, io_jlo:io_jhi+1, io_klo:io_khi,1:localNumBlocks) 

              varid = 0 !Unused

              call io_ncmpi_write_facevar(fileID, & 
                   varid, &
                   JAXIS, &
                   NXB, & 
                   NYB+1, & 
                   NZB, & 
                   faceYBuf, & 
                   globalFaceYMin(i), &
                   globalFaceYMax(i), &
                   io_faceYVarLabels(i), &
                   localNumBlocks, &
                   gr_globalNumBlocks,  & 
                   localOffset)

           end if !NDIM .gt. 1

           if (NDIM .gt. 2) then
              faceZBuf(1,1:NXB,1:NYB,1:NZB+1,1:localNumBlocks) = &
                   facevarz(i,io_ilo:io_ihi, io_jlo:io_jhi, io_klo:io_khi+1,1:localNumBlocks) 

              varid = 0 !Unused

              call io_ncmpi_write_facevar(fileID, & 
                   varid, &
                   KAXIS, &
                   NXB, & 
                   NYB, & 
                   NZB+1, & 
                   faceZBuf, & 
                   globalFaceZMin(i), &
                   globalFaceZMax(i), &
                   io_faceZVarLabels(i), &
                   localNumBlocks, &
                   gr_globalNumBlocks,  & 
                   localOffset)

           end if !NDIM .gt. 2

        end if !io_doublePrecision

     end do !i = 1,NFACE_VARS

     deallocate(faceXBuf)
     if(NDIM .gt. 1) deallocate(faceYBuf)
     if(NDIM .gt. 2) deallocate(faceZBuf)

     deallocate(globalFaceXMin)
     deallocate(globalFaceXMax)
     if(NDIM .gt. 1) then
        deallocate(globalFaceYMin)
        deallocate(globalFaceYMax)
     end if
     if(NDIM .gt. 2) then
        deallocate(globalFaceZMin)
        deallocate(globalFaceZMax)
     end if


#endif
  end if


  if (xferType == IO_READ_XFER) then

     idx = 0
     do i = UNK_VARS_BEGIN,UNK_VARS_END

        idx = idx + 1
        s = io_unkLabelsGlobal(idx)
        call io_ncmpi_name_to_id(fileID, s, len_trim(s), varid)
        call io_ncmpi_read_unknowns(io_chkptFileID, &
             varid, &
             NXB, &
             NYB, &
             NZB, &
             unkBuf, &
             localNumBlocks, &
             gr_globalNumBlocks,  &
             localOffset)



        unk(i,io_ilo:io_ihi, io_jlo:io_jhi, io_klo:io_khi,1:MAXBLOCKS) = & 
             unkBuf(1,1:NXB,1:NYB,1:NZB,1:MAXBLOCKS)


     enddo


     deallocate(unkBuf)

#if (NFACE_VARS > 0)
     allocate(faceXBuf(1, NXB+1, NYB, NZB, 1:MAXBLOCKS))
     if(NDIM > 1) allocate(faceYBuf(1, NXB, NYB+1, NZB, 1:MAXBLOCKS))
     if(NDIM > 2) allocate(faceZBuf(1, NXB, NYB, NZB+1, 1:MAXBLOCKS))


     do i = 1,NFACE_VARS
        s = io_faceXVarLabels(i)
        call io_ncmpi_name_to_id(fileID, s, len_trim(s), varid)
        call io_ncmpi_read_unknowns(io_chkptFileID, &
             varid, &
             NXB + 1, &
             NYB, &
             NZB, &
             faceXBuf, &
             localNumBlocks, &
             gr_globalNumBlocks, &
             localOffset)

        facevarx(i,io_ilo:io_ihi+1,io_jlo:io_jhi,io_klo:io_khi,1:MAXBLOCKS) = &
             faceXBuf(1,1:NXB+1,1:NYB,1:NZB,1:MAXBLOCKS)
        if(NDIM > 1) then
           s = io_faceYVarLabels(i)
           call io_ncmpi_name_to_id(fileID, s, len_trim(s), varid)
           call io_ncmpi_read_unknowns(io_chkptFileID, &
                varid, &
                NXB, &
                NYB + 1, &
                NZB, &
                faceYBuf, &
                localNumBlocks, &
                gr_globalNumBlocks, &
                localOffset)

           facevary(i,io_ilo:io_ihi,io_jlo:io_jhi+1,io_klo:io_khi,1:MAXBLOCKS) = &
                faceYBuf(1,1:NXB,1:NYB+1,1:NZB,1:MAXBLOCKS)
        end if
        if(NDIM > 2) then
           s = io_faceZVarLabels(i)
           call io_ncmpi_name_to_id(fileID, s, len_trim(s), varid)
           call io_ncmpi_read_unknowns(io_chkptFileID, &
                varid, &
                NXB, &
                NYB, &
                NZB + 1, &
                faceZBuf, &
                localNumBlocks, &
                gr_globalNumBlocks, &
                localOffset)

           facevarz(i,io_ilo:io_ihi,io_jlo:io_jhi,io_klo:io_khi+1,1:MAXBLOCKS) = &
                faceZBuf(1,1:NXB,1:NYB,1:NZB+1,1:MAXBLOCKS)
        end if

     end do

     deallocate(faceXBuf)
     if(NDIM > 1) deallocate(faceYBuf)
     if(NDIM > 2) deallocate(faceZBuf)
#endif

  end if

end subroutine io_xfer_mesh_data
