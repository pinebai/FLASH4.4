!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserComm/Async/ed_commFinalize
!!
!! NAME
!!  ed_commFinalize
!!
!! SYNOPSIS
!!  call ed_commFinalize()
!!
!! DESCRIPTION
!!  Finalizes the laser communication.
!!
!! ARGUMENTS
!!  No arguments
!!
!!***

subroutine ed_commFinalize()
  use ed_commData, ONLY : ed_commLog
  !$ use ed_commData, ONLY : ed_commLock
  use Logfile_interface, ONLY : Logfile_close
  use UTCounter_sharedCounter, ONLY : UTCounter_finalize
  implicit none

  call UTCounter_finalize()
  if (ed_commLog) then
     call Logfile_close(.true.)
  end if
  
  !$ call omp_destroy_lock(ed_commLock)
end subroutine ed_commFinalize
