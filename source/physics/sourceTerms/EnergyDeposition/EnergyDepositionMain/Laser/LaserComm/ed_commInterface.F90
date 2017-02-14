!!****ih* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserComm/ed_commInterface
!!
!! NAME
!!
!!  ed_commInterface
!!
!! SYNOPSIS
!!
!!  use ed_commInterface
!!
!!***

module ed_commInterface
  implicit none

  interface
     subroutine ed_commIsTransportDone(isTransportDone)
       implicit none
       logical, intent(OUT) :: isTransportDone
     end subroutine ed_commIsTransportDone
  end interface

  interface
     subroutine ed_commProgressTransport(optionalForceProgress)
       implicit none
       logical, optional, intent(IN) :: optionalForceProgress
     end subroutine ed_commProgressTransport
  end interface

  interface
     subroutine ed_commProgressPipeline(doFlush)
       implicit none
       logical, intent(IN) :: doFlush
     end subroutine ed_commProgressPipeline
  end interface

  interface
     subroutine ed_commProgressCounter()
       implicit none
     end subroutine ed_commProgressCounter
  end interface

  interface
     subroutine ed_commHandleOffBlkRay(ray)
       implicit none
       integer, intent(IN) :: ray
     end subroutine ed_commHandleOffBlkRay
  end interface

  interface
     subroutine ed_commIncrementDeadRays(numDeadRays)
       implicit none
       integer, intent(IN) :: numDeadRays
     end subroutine ed_commIncrementDeadRays
  end interface

  interface
     subroutine ed_commCheckTermination(finalCount)
       implicit none
       integer, intent(IN) :: finalCount
     end subroutine ed_commCheckTermination
  end interface

  interface
     subroutine ed_commCheckRayExchange()
       implicit none
     end subroutine ed_commCheckRayExchange
  end interface

  interface
     subroutine ed_commCheckRayLocation(fromRay, toRay)
       implicit none
       integer, intent(IN) :: fromRay, toRay
     end subroutine ed_commCheckRayLocation
  end interface

  interface
     subroutine ed_commGetNewRays()
       implicit none
     end subroutine ed_commGetNewRays
  end interface

  interface
     subroutine ed_commSortParticles()
       implicit none
     end subroutine ed_commSortParticles
  end interface

  interface
     !This is a C function
     subroutine ed_comm_sort_particles(rays, pRayCount)
       implicit none
       real, dimension(*), intent(INOUT) :: rays
       integer, intent(INOUT) :: pRayCount
     end subroutine ed_comm_sort_particles
  end interface

  interface
     subroutine ed_commInit
       implicit none
     end subroutine ed_commInit
  end interface

  interface
     subroutine ed_commFinalize()
       implicit none
     end subroutine ed_commFinalize
  end interface

  interface
     subroutine ed_commInitComm()
       implicit none
     end subroutine ed_commInitComm
  end interface
end module ed_commInterface
