!!****if* source/IO/IOTypes/io_getZeroBasedVarInfo
!!
!! NAME
!!  io_getZeroBasedVarInfo
!!
!! SYNOPSIS
!!
!!  io_getZeroBasedVarInfo(integer(IN) :: fileType,
!!                         integer(IN) :: gridDataStruct,
!!                         integer(OUT) :: numGridVars,
!!                         integer(OUT) :: numOutputGridVars,
!!                         integer(OUT) :: gridVarOffsets(MAX_MESH_VAR),
!!      char(OUT)(len=MAX_STRING_LENGTH) :: gridVarLabels(MAX_MESH_VAR))
!!
!!
!! DESCRIPTION
!!
!! This subroutine returns information about the variables that
!! will be written to file for both checkpoint and plot files.
!!
!!
!! ARGUMENTS
!!
!!  fileType: The file type, either checkpoint file or plot file.
!!  gridDataStruct: The grid data structure, e.g. UNK, FACEX
!!  numGridVars: The total number of mesh variables in gridDataStruct.
!!  numOutputGridVars: The total number of mesh variables in gridDataStruct
!!                     that will be written to file (depends on whether the
!!                     file type is checkpoint file or plot file).
!!  gridVarOffsets: The zero-based index of each variable that will be
!!                  written to file.
!!  gridVarLabels: The name of each variable that will be written to file.
!!
!! NOTES
!!
!! Given a file type and grid data struture this routine provides global
!! information about the corresponding mesh variables.  Note the word global!
!! This routine knows nothing about mesh replication.  It is the responsibility
!! of other subroutines to determine the local memory offset for the
!! global variables that belong to this process.
!!
!!***

#include "constants.h"
#include "Flash.h"
#include "io_flash.h"

subroutine io_getZeroBasedVarInfo(fileType, gridDataStruct, numGridVars, &
     numOutputGridVars, gridVarOffsets, gridVarLabels)

  use Driver_interface, ONLY : Driver_abortFlash
  use IO_data, ONLY : io_unkLabels, io_faceXVarLabels, io_faceYVarLabels, &
       io_faceZVarLabels, io_scratchGridVarlabels, io_nPlotVars, &
       io_nPlotGridVars, io_plotVar, io_plotGridVar, io_unkLabelsGlobal, io_plotVarStr
  use ut_qsortInterface, ONLY : ut_qsort

  implicit none
  integer, intent(IN) :: fileType, gridDataStruct
  integer, intent(OUT) :: numGridVars, numOutputGridVars 
  integer, dimension(MAX_MESH_VAR), intent(OUT) :: gridVarOffsets
  character (len=MAX_STRING_LENGTH), dimension(MAX_MESH_VAR), intent(OUT) :: &
       gridVarLabels
  integer :: s, i


  !Redundant pieces of code are excluded using the preprocessor.  The logic does 
  !not need the preprocessor fragments, but the IBM compiler does when compiling 
  !with array-bounds checking.  Without the preprocessor fragments we get:
  !"(S) Zero-sized arrays must not be subscripted." compilation error message.

  numGridVars = 0
  numOutputGridVars = 0

  select case(fileType)
  case (CHECKPOINTFILE)


     select case(gridDataStruct)


     case (CENTER)
        numGridVars = ubound(io_unkLabelsGlobal,1)
        numOutputGridVars = numGridVars
#if NUNK_VARS > 0
        do s = 1, numOutputGridVars
           gridVarOffsets(s) = s - 1
           gridVarLabels(s) = trim(io_unkLabelsGlobal(s))
        end do
#endif


     case (FACEX)
        numGridVars = NFACE_VARS
        numOutputGridVars = numGridVars
#if NFACE_VARS > 0
        do s = 1, numOutputGridVars
           gridVarOffsets(s) = s - 1
           gridVarLabels(s) = trim(io_faceXVarLabels(s))
        end do
#endif


     case (FACEY)
        numGridVars = NFACE_VARS
        if (NDIM > 1) then
           numOutputGridVars = numGridVars
#if NFACE_VARS > 0
           do s = 1, numOutputGridVars
              gridVarOffsets(s) = s - 1
              gridVarLabels(s) = trim(io_faceYVarLabels(s))
           end do
#endif
        else
           numOutputGridVars = 0
        end if


     case (FACEZ)
        numGridVars = NFACE_VARS
        if (NDIM > 2) then
           numOutputGridVars = numGridVars
#if NFACE_VARS > 0
           do s = 1, numOutputGridVars
              gridVarOffsets(s) = s - 1
              gridVarLabels(s) = trim(io_faceZVarLabels(s))
           end do
#endif
        else
           numOutputGridVars = 0
        end if


     case (SCRATCH)
        numGridVars = NSCRATCH_GRID_VARS
        numOutputGridVars = io_nPlotGridVars
#if NSCRATCH_GRID_VARS > 0
        if (numOutputGridVars > 0) then
           gridVarOffsets = pack(io_plotGridVar-1, io_plotGridVar /= NONEXISTENT)
           call ut_qsort(gridVarOffsets, numOutputGridVars, ascOrderArg=.true.)
           do s = 1, numOutputGridVars
              gridVarLabels(s) = &
                   trim(io_scratchGridVarlabels(gridVarOffsets(s) + 1))
           end do
        end if
#endif


     case DEFAULT
        call Driver_abortFlash ("Checkpoint file data structure not recognised")
     end select


  case (PLOTFILE)

     !NOTE: The arrays io_plotVarStr, io_plotGridVarStr, io_plotFaceVarStr are 
     !of no use because they contain the plot file string names in the 
     !order given in flash.par.  In order to make a fast derived data type 
     !over memory I pick up the plot file variables in memory order.  As such, 
     !to obtain the variable names in the correct order I use the original 
     !arrays io_unkLabels, io_scratchGridVarlabels, io_faceXVarLabels, 
     !io_faceYVarLabels, io_faceZVarLabels.

     select case(gridDataStruct)


     case (CENTER)
        numGridVars = ubound(io_unkLabelsGlobal,1)
        numOutputGridVars = io_nPlotVars
#if NUNK_VARS > 0
        if (numOutputGridVars > 0) then
           i = 1
           do s = 1, numGridVars
              !We assume that the first io_nPlotVars strings in the list io_plotVarStr
              !are all valid names; IO_init should guarantee this. - KW
              !not all strings in the first io_nPlotVars elements of
              !io_plotVarStr are guaranteed to be valid strings.
              if (any(io_unklabelsGlobal(s) == io_plotVarStr(1:io_nPlotVars))) then
                 gridVarOffsets(i) = s - 1
                 i = i + 1
              end if        
           end do

           call ut_qsort(gridVarOffsets, numOutputGridVars, ascOrderArg=.true.)
           do s = 1, numOutputGridVars
              gridVarLabels(s) = &
                   trim(io_unkLabelsGlobal(gridVarOffsets(s) + 1))
           end do
        end if
#endif


     case (FACEX, FACEY, FACEZ)
        numGridVars = NFACE_VARS 
        numOutputGridVars = 0 !Face plot files not yet implemented.


     case (SCRATCH)
        numGridVars = NSCRATCH_GRID_VARS
        numOutputGridVars = io_nPlotGridVars
#if NSCRATCH_GRID_VARS > 0
        if (numOutputGridVars > 0) then
           gridVarOffsets = pack(io_plotGridVar-1, io_plotGridVar /= NONEXISTENT)
           call ut_qsort(gridVarOffsets, numOutputGridVars, ascOrderArg=.true.)
           do s = 1, numOutputGridVars
              gridVarLabels(s) = &
                   trim(io_scratchGridVarlabels(gridVarOffsets(s) + 1))
           end do
        end if
#endif
      
  
     case DEFAULT
        call Driver_abortFlash("Plot file data structure not recognised")
     end select


  case DEFAULT
     call Driver_abortFlash("File type not recognised")
  end select

end subroutine io_getZeroBasedVarInfo
