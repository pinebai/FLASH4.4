!!****if* source/IO/IOTypes/io_do_xfer
!!
!! NAME
!!  io_do_xfer
!!
!! SYNOPSIS
!!
!!  io_do_xfer(integer(in) :: xferType,
!!             integer(in) :: gridStruct,
!!             character(len=*)(in) :: dataset,
!!             logical(out) :: doXfer)
!!
!! DESCRIPTION
!!
!! This subroutine is used to determine whether my process should
!! transfer a specified mesh variable to/from file.  It is required
!! because the mesh variables that exist locally are a subset of 
!! the global mesh variables when using mesh replicatoin.
!!
!! ARGUMENTS
!!
!! xferType: The direction of data transfer:
!!           IO_WRITE_XFER - write, IO_READ_XFER - read.
!! gridStruct: The mesh data structure, e.g. unk, facex, facey, ... e.t.c.
!! dataset: The mesh variable string name.
!! doXfer: Whether or not my process should transfer the specified mesh
!!         variable to/from file - depends on whether the mesh variable
!!         exists locally.
!!
!!***


#include "constants.h"
#include "Flash.h"
#include "io_flash.h"

subroutine io_do_xfer(xferType, gridStruct, dataset, doXfer)
  use Simulation_interface, ONLY : Simulation_mapStrToInt
  use IO_data, ONLY : io_acrossMe, io_unkNonRep
  implicit none
  integer, intent(IN) :: xferType, gridStruct
  character(len=*), intent(IN) :: dataset
  logical, intent(OUT) :: doXfer
  integer :: varIndex

  doXfer = .true.
  if (gridStruct == CENTER) then
     !Mesh replication can only happen for UNK.
     doXfer = .false.
     !Check if local view contains a NONEXISTENT mesh variable.
     call Simulation_mapStrToInt(trim(dataset), varIndex, MAPBLOCK_UNK)
     if (varIndex /= NONEXISTENT) then
        if (xferType == IO_WRITE_XFER) then
           ! only write mesh replicated data from mesh 0
           if (io_acrossMe .eq. 0 .or. io_unkNonRep(varIndex) > 0) then
              doXfer = .true.
           end if
        else if (xferType == IO_READ_XFER) then
           ! everyone reads mesh data
           doXfer = .true.
        end if
     end if
  end if
end subroutine io_do_xfer
