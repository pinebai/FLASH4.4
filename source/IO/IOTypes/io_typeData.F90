!!****if* source/IO/IOTypes/io_typeData
!!
!! NAME
!!  io_typeData
!!
!! SYNOPSIS
!!
!!  use io_typeData
!!
!! DESCRIPTION 
!!  
!!  Holds all the IO data that is needed by the IOtype sub unit
!!
!! ARGUMENTS
!!
!!  none    
!!
!!
!!***

#include "constants.h"
#include "Flash.h"

!Used only with derived datatype I/O (names may soon change)
!The names will eventually change to something like io_typeXXX.
module io_typeData
  implicit none
  logical, save :: io_packMeshPlotWriteHDF5
  logical, save :: io_packMeshChkWriteHDF5
  logical, save :: io_packMeshChkReadHDF5
  logical, save :: io_asyncMeshPlotWritePnet
  logical, save :: io_asyncMeshChkWritePnet
  logical, save :: io_asyncMeshChkReadPnet
  logical, save :: io_useLegacyLabels
  integer, parameter :: io_legacyLabelLength = 4
end module io_typeData
