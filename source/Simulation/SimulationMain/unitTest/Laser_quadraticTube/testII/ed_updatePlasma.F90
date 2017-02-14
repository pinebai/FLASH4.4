!!****if* source/Simulation/SimulationMain/unitTest/Laser_quadraticTube/testII/ed_updatePlasma
!!
!! NAME
!!
!!  ed_updatePlasma
!!
!! SYNOPSIS
!!
!!  call ed_updatePlasma (integer, (in)           :: blockCount, 
!!                        integer, (in)           :: blockList (:),
!!                        real     (in), optional :: scaleFact)
!!
!! DESCRIPTION
!!
!!  Overriding (empty) routine from the original version to avoid calling Eos with
!!  nonsensical data.
!!
!! ARGUMENTS
!!
!!  blockCount : Number of blocks on current processor
!!  blockList  : All block ID numbers
!!  scaleFact  : optional scaling factor
!!
!! NOTES
!!
!!***

subroutine ed_updatePlasma (blockCount,blockList,scaleFact)

  implicit none
  
  integer,           intent (in) :: blockCount
  integer,           intent (in) :: blockList (1:blockCount)
  real   , optional, intent (in) :: scaleFact
!
!
!    ...Ready!
!
!
  return
end subroutine ed_updatePlasma
