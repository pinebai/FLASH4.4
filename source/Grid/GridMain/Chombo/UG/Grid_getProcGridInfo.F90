!!****if* source/Grid/GridMain/Chombo/UG/Grid_getProcGridInfo
!!
!! NAME
!!  Grid_getProcGridInfo
!!
!! SYNOPSIS
!!
!!  Grid_getProcGridInfo(integer(OUT)  :: comm(MDIM),
!!               integer(OUT)  :: procGrid(MDIM),
!!               integer(OUT)  :: me(MDIM))
!!               
!!  
!! DESCRIPTION 
!!
!!  For Uniform Grid, and for the purpose of running the pfft routines
!!  the processors are logically configured as a cartesian grid, where
!!  each grid point represents a processor. Each dimension has its own
!!  communicator which is split from MPI_COMM_WORLD. This routine returns
!!  the communicators, the number of processors in each communicator,
!!  and the identity of the processor in each communicators.
!! 
!!
!! ARGUMENTS
!!
!!  comm        -  communicator along each dimension
!!  procGrid    - number of processors in each communicator
!!  me          - identity of MyPE in each of the communicators
!!               
!!
!!***

subroutine Grid_getProcGridInfo(comm,procGrid,me)

  use Grid_data, ONLY : gr_me, gr_procGrid, gr_comm

#include "constants.h"

implicit none
  integer, dimension(MDIM), intent(out)  :: comm, procGrid,me

  comm(1:MDIM)=gr_comm(1:MDIM)
  procGrid(1:MDIM)= gr_procGrid(1:MDIM)
  me(1:MDIM) = gr_me(1:MDIM)

  return
end subroutine Grid_getProcGridInfo
