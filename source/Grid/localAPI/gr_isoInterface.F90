!!****ih* source/Grid/localAPI/gr_isoInterface
!!
!! NAME
!!   gr_isoInterface
!!
!! SYNOPSIS
!!   use gr_isoInterface
!!
!! DESCRIPTION
!! This is the header file for the GridSolvers isolated boundary
!!  conditions subunit that defines its interfaces.
!!
!!***

module gr_isoInterface
  implicit none
#include "constants.h"
#include "Flash.h"

  interface
     subroutine gr_isoFindMassCenter (idensvar)
       integer, intent(IN)  :: idensvar
     end subroutine gr_isoFindMassCenter
  end interface

  interface
     subroutine gr_isoImageBdry (iiden, iipot, poisfact)
       integer,intent(IN)       :: iiden, iipot
       real,intent(IN)          :: poisfact
     end subroutine gr_isoImageBdry
  end interface

  interface
     subroutine gr_isoImageMass (isoln, iiden)
       integer,intent(IN)       :: isoln, iiden
     end subroutine gr_isoImageMass
  end interface

  interface
     subroutine gr_isoMpoleInit
     end subroutine gr_isoMpoleInit
  end interface

  interface
     subroutine gr_isoMpoleMoments (idensvar)
       integer, intent(IN)   :: idensvar
     end subroutine gr_isoMpoleMoments
  end interface

  interface
     subroutine gr_isoMpolePotential (ipotvar, poisfact)
       integer, intent(IN)   :: ipotvar
       real, intent(IN)      :: poisfact
     end subroutine gr_isoMpolePotential
  end interface

  interface
     subroutine gr_isoSumLocal (lsum, nsum, blockID, idensvar)
       integer,intent(IN) :: blockID, idensvar, nsum
       real,intent(OUT)    :: lsum(nsum)
     end subroutine gr_isoSumLocal
  end interface

  interface
     subroutine gr_isoZoneMoments (xprime, yprime, zprime, zonemass)
       real, intent(IN) :: xprime, yprime, zprime, zonemass
     end subroutine gr_isoZoneMoments
  end interface

  interface
     subroutine gr_isoZonepotential (xprime, yprime, zprime, potential)
       real, intent(IN)   :: xprime, yprime, zprime
       real, intent(OUT)  :: potential
     end subroutine gr_isoZonepotential
  end interface


end module gr_isoInterface






