!!****f* source/Simulation/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer(IN) :: blockId)
!!
!!
!!
!! DESCRIPTION
!!  This routine applies initial conditions of a specific simulation
!!  to the specified block.
!!
!! 
!! ARGUMENTS
!!
!!  blockId -         the number of the block to update
!!
!! PARAMETERS
!!
!!  eosModeInit -     after this routine sets up initial conditions,
!!                    the grid package calls Eos to insure the
!!                    values are thermodynamically consistent.  This
!!                    parameter controls the mode of application of
!!                    the Eos.  Its default is "dens_ie", and it can
!!                    be one of "dens_ie", "dens_pres", "dens_temp".
!!                    Setting this value to dens_temp, for instance,
!!                    would make it possible to leave this routine
!!                    having just initialized density and temperature,
!!                    leaving Eos to calculate the rest.
!!
!! SEE ALSO
!!
!!  Eos_wrapped
!!***

subroutine Simulation_initBlock(blockId)

  
  implicit none
  
  integer, intent(in) :: blockId
  
  return
end subroutine Simulation_initBlock
