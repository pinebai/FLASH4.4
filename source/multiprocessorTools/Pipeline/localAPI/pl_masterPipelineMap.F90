!!****if* source/multiprocessorTools/Pipeline/localAPI/pl_masterPipelineMap
!!
!! NAME
!!  
!!  pl_masterPipelineMap
!!
!! SYNOPSIS
!! 
!!  call pl_masterPipelineMap ()
!!
!! DESCRIPTION
!!
!!  This routine sets up the pipeline map at the master rank. The pipeline map consists
!!  of a two-dimensional integer array, which contains info about how many channels (#)
!!  each rank has and also all the processor ID's of these channels:
!!
!!               Rank| # ||    channel proc ID's
!!               -----------------------------------
!!                 0 | 4 || 11 |  5 |  6 |  7 |
!!                 1 | 1 || 13 |
!!                 2 | 3 ||  1 |  8 |  5 |  4 |
!!                 3 | 5 || 17 | 15 | 16 |  6 |  4 |
!!                  .............................
!!
!!  For convenience, the pipeline map is declared such that the ranks are located in the
!!  columns and the channels in the rows of the pipeline map array.
!!
!! ARGUMENTS
!!
!!  none
!!
!! NOTES
!!
!!  After exiting this routine, the maximum number of channels per processor has been
!!  set on all processors and the pipeline map has been allocated and created on the master.
!!
!!***

subroutine pl_masterPipelineMap ()

  implicit none

  return
end subroutine pl_masterPipelineMap
