!!****if* source/IO/IOMain/io_isPlotVar
!!
!! NAME
!!  io_isPlotVar
!!
!! SYNOPSIS
!!
!!  io_isPlotVar(integer(in)  :: var,
!!               logical(out) :: val,
!!               integer(in)  :: map)
!!
!! DESCRIPTION
!!  Given a variable ID determine if it is to be output to a plotfile.
!! 
!!  This routine can handle both unk vars and scratch grid vars.
!!  
!!  Variables that are outputted to plotfile must be specified in the
!!  flash.par as (ex) plot_var_1 = 'dens' for unk vars and as
!!  plot_grid_var_1 = "vrtz" for grid vars.
!!
!! ARGUMENTS
!!
!!  var - index of variable in the unk data structure
!!  val - logical returned value, true if var is plotvar, false otherwise
!!  map - integer value defined in constants.h options are MAPBLOCK_UNK or
!!        MAPBLOCK_SCRATCH
!!
!!***


subroutine io_isPlotVar(var, val, map)

  use IO_data, ONLY : io_plotVar, io_nPlotVars, io_plotGridVar, io_maxPlotGridVars
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

#include "constants.h"
#include "Flash.h"

  integer, intent(in) :: var, map
  logical, intent(inout) :: val
  integer :: i

  if(map == MAPBLOCK_UNK) then

     do i=1, io_nPlotVars
        if(var == io_plotVar(i)) then
           val = .true.
           return
        end if
     end do
     
     val = .false.

  else if(map == MAPBLOCK_SCRATCH) then
     
     do i=1, io_maxPlotGridVars
        if(var == io_plotGridVar(i)) then
           val = .true.
           return
        end if
     end do
     
     val = .false.

  else

     call Driver_abortFlash("Error: io_isPlotVar, unknown map block")

  end if

end subroutine io_isPlotVar
