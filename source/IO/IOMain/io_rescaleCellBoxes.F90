!!****if* source/IO/IOMain/io_rescaleCellBoxes
!!
!! NAME
!!
!!  io_rescaleCellBoxes
!!
!!
!! SYNOPSIS
!!
!!  call io_rescaleCellBoxes() 
!!
!!
!! DESCRIPTION
!!
!!  Apply linear transformations (stretching / shifting) to the cell
!!  coordinates.
!!
!!  To be called immediately after io_readData has filled in various
!!  Paramesh metadata arrays with information from a checkpoint file
!!  upon restart.
!!
!!  This is a hack.
!!
!!  It is the user's responsibility to adjust the xmin,xmax,ymin,...,zmax
!!  Runtime parameters appropriately when coordinate rescaling is used.
!!
!! ARGUMENTS
!!
!!  none
!!
!! NOTES
!!
!!  This manipulates data structures owned by PARAMESH. Therefore it does not
!!  do anything if a different Grid implementation is used.
!!***

#ifdef DEBUG_ALL
#define DEBUG_IO
#endif

#include "constants.h"
#include "Flash.h"

subroutine io_rescaleCellBoxes()

  use Driver_interface, ONLY : Driver_abortFlash
  use Logfile_interface, ONLY : Logfile_stamp

  use IO_data, ONLY : io_globalMe

#ifdef FLASH_GRID_PARAMESH
  use tree, ONLY : lnblocks, bnd_box, coord, bsize
#endif

  implicit none

  integer :: lb

  !! These should of course all become runtime parameters! - KW
  real,parameter :: io_restartRescaleAi = 0.0
  real,parameter :: io_restartRescaleAj = 0.0
  real,parameter :: io_restartRescaleAk = 0.0
  real,parameter :: io_restartRescaleBi = 1.0
  real,parameter :: io_restartRescaleBj = 1.0
  real,parameter :: io_restartRescaleBk = 1.0


  if (io_restartRescaleAi == 0.0 .AND. &
      io_restartRescaleAj == 0.0 .AND. &
      io_restartRescaleAk == 0.0 .AND. &
      io_restartRescaleBi == 1.0 .AND. &
      io_restartRescaleBj == 1.0 .AND. &
      io_restartRescaleBk == 1.0) then
     return                     ! Nothing to do, RETURN!
  end if
      

#ifndef FLASH_GRID_PARAMESH
  if (io_globalMe == MASTER_PE) then
     print*,"NOTE: request to rescale cell boxes is being ignored!"
     call Logfile_stamp( "Cell boxes are NOT being rescaled!", "[io_rescaleCellBoxes]")
  end if

#else

  if (io_globalMe == MASTER_PE) then
     print*,"Cell boxes are being rescaled!"
     call Logfile_stamp( "Cell boxes are being rescaled!", "[io_rescaleCellBoxes]")
  end if


  do lb=1,lnblocks
     bnd_box(1,1,lb) = scaleI(bnd_box(1,1,lb))
     bnd_box(2,1,lb) = scaleI(bnd_box(2,1,lb))
     coord(1,lb) = scaleI(coord(1,lb))
     bsize(1,lb) = bsize(1,lb) * io_restartRescaleBi

#if N_DIM > 1
     bnd_box(1,2,lb) = scaleJ(bnd_box(1,2,lb))
     bnd_box(2,2,lb) = scaleJ(bnd_box(2,2,lb))
     coord(2,lb) = scaleJ(coord(2,lb))
     bsize(2,lb) = bsize(2,lb) * io_restartRescaleBj
#endif

#if N_DIM > 2
     bnd_box(1,3,lb) = scaleK(bnd_box(1,3,lb))
     bnd_box(2,3,lb) = scaleK(bnd_box(2,3,lb))
     coord(3,lb) = scaleK(coord(3,lb))
     bsize(3,lb) = bsize(3,lb) * io_restartRescaleBk
#endif
  end do

contains

  real function scaleI(x)
    real x
    scaleI = io_restartRescaleAi + io_restartRescaleBi * x
  end function scaleI
  real function scaleJ(y)
    real y
    scaleJ = io_restartRescaleAj + io_restartRescaleBj * y
  end function scaleJ
  real function scaleK(z)
    real z
    scaleK = io_restartRescaleAk + io_restartRescaleBk * z
  end function scaleK

#endif

end subroutine io_rescaleCellBoxes
