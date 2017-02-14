!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserComm/ed_commProgressTransport
!!
!!  NAME     
!!   ed_commProgressTransport
!!
!!  SYNOPSIS
!!   call ed_commProgressTransport(logical, optional, intent(IN) :: optionalForceProgress)
!!
!!  DESCRIPTION 
!!   This is a stub.  It is only needed when using asynchronous laser
!!   communication.
!!
!!  ARGUMENTS
!!   optionalForceProgress: Forces communication progress.
!!
!!***

subroutine ed_commProgressTransport(optionalForceProgress)
  implicit none
  logical, optional, intent(IN) :: optionalForceProgress
end subroutine ed_commProgressTransport
