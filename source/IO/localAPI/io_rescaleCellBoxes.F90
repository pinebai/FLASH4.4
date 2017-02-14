!!****if* source/IO/localAPI/io_rescaleCellBoxes
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


subroutine io_rescaleCellBoxes()

  implicit none


end subroutine io_rescaleCellBoxes
