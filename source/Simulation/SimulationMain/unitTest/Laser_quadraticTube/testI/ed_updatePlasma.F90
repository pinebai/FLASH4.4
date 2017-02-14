!!****if* source/Simulation/SimulationMain/unitTest/Laser_quadraticTube/testI/ed_updatePlasma
!!
!! NAME
!!
!!  ed_updatePlasma
!!
!! SYNOPSIS
!!
!!  call ed_updatePlasma (integer, intent (in) :: blockCount, 
!!                        integer, intent (in) :: blockList (:))
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
!!
!! NOTES
!!
!!***

subroutine ed_updatePlasma (blockCount,blockList)

  implicit none
  
  integer, intent (in) :: blockCount
  integer, intent (in) :: blockList (1:blockCount)
!
!
!    ...Ready!
!
!
  return
end subroutine ed_updatePlasma
