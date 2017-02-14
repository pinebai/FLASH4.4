!!****if* source/physics/Hydro/HydroMain/unsplit/MHD_StaggeredMesh/hy_uhd_biermannSource
!!
!! NAME
!!
!!  hy_uhd_biermannSource
!!
!! SYNOPSIS
!!
!!  hy_uhd_biermannSource( integer (IN) :: blockCount,
!!                         integer (IN) :: blockList(blockCount),
!!                         real    (IN) :: dt )
!!
!! DESCRIPTION
!!
!! Implement Biermann Battery Term as a source to the magnetic field. This is a stub. We are not supporting this yet. 
!!
!! ARGUMENTS
!!
!!  blockCount -  the number of blocks in blockList
!!  blockList  -  array holding local IDs of blocks on which to advance
!!  dt         -  timestep
!!***

Subroutine hy_uhd_biermannSource ( blockCount, blockList, dt )
  use Grid_interface, ONLY: Grid_getBlkIndexLimits, &
                            Grid_getDeltas,         &
                            Grid_getBlkPtr,         &
                            Grid_releaseBlkPtr

  use Hydro_data, ONLY : hy_biermannCoef,  &
                         hy_useBiermann,   &
                         hy_useBiermann1T, &
                         hy_bier1TZ,       &
                         hy_bier1TA,       &
                         hy_avogadro,      &
                         hy_qele,          &
                         hy_speedOfLight

  use hy_uhd_slopeLimiters, ONLY : minmod

  use Eos_interface, ONLY : Eos_wrapped
  
  implicit none

#include "Flash.h"
#include "constants.h"
#include "UHD.h"

  ! Arguments:
  integer, intent(IN) :: blockCount
  integer, intent(IN) :: blockList(blockCount)
  real,    intent(IN) :: dt

  return

End Subroutine hy_uhd_biermannSource
