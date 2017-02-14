!!****ih* source/Grid/GridParticles/GridParticlesMove/Sieve/gr_ptSieveInterface
!!
!! NAME 
!!
!!   gr_ptSieveInterface
!!
!! SYNOPSIS
!!   
!!   use gr_ptSieveInterface
!!
!! DESCRIPTION
!! 
!! This is the header file for the GridParticles
!! subunit that defines its interfaces.
!!
!!***

#include "constants.h"
#include "Flash.h"

module gr_ptSieveInterface
  implicit none

  interface
     subroutine gr_ptResetProcPair(haveNonLocalData, mustCommunicate)

       logical, intent(IN)  :: haveNonLocalData
       logical, intent(OUT) :: mustCommunicate
  
     end subroutine gr_ptResetProcPair
  end interface
  
  interface

     subroutine gr_ptNextProcPair(timesInLoop, src, dest,  haveNonLocalData, mustCommunicate)
       integer, intent(IN)  :: timesInLoop
       integer, intent(OUT) :: src, dest
       logical, intent(IN)  :: haveNonLocalData
       logical, intent(OUT) :: mustCommunicate
     end subroutine gr_ptNextProcPair
  end interface
     
end module gr_ptSieveInterface
