!!****if* source/diagnostics/ProtonImaging/ProtonImagingMain/oldOneTimeStep/ProtonIO/pi_IOinit
!!
!! NAME
!!
!!  pi_IOinit
!!
!! SYNOPSIS
!!
!!  call pi_IOinit ()
!!
!! DESCRIPTION
!!
!!  This routine initializes variables for plotting proton data. As this needs
!!  the exact number of total protons that will be emitted, this routine is only
!!  successful, if all the beams have been set up. The number of emitted protons
!!  might change slightly during the beams setup, so we cannot use the requested
!!  proton emission number values from the flash.par.
!!
!!***

subroutine pi_IOinit ()
  
  use Driver_interface,    ONLY : Driver_abortFlash

  use ProtonImaging_data,  ONLY : pi_beamsAreSetup,            &
                                  pi_IOmaxBlockCrossingNumber, &
                                  pi_IOmaxPointsPerBlock,      &
                                  pi_IOmaxProtonCount,         &
                                  pi_IOnumberOfProtons2Plot,   &
                                  pi_IOprotonPointCount,       &
                                  pi_IOprotonPoints,           &
                                  pi_IOprotonTags,             &
                                  pi_IOprotonWriteModulo,      &
                                  pi_maxProtonCount,           &
                                  pi_totalProtons2Launch,      &
                                  pi_useIOprotonPlot

  implicit none

#include "constants.h"
#include "Flash.h"
!
!
!     ...Check, if proton plotting is wanted. If not, return at once.
!        Abort also, if the beams are not set up.
!
!
  if (.not. pi_useIOprotonPlot) return
  if (.not. pi_beamsAreSetup) then
       call Driver_abortFlash ("pi_IOinit: Proton beams are not set up!")
  end if
!
!
!     ...Plotting is wanted. Check, if number of IO protons is ok.
!
!
  if (pi_IOnumberOfProtons2Plot < 1) then
      call Driver_abortFlash ("pi_IOinit: pi_IOnumberOfProtons2Plot is < 1")
  end if
!
!
!     ...Calculate the maximum number of points per IO proton and per block.
!        This number is calculated from the supplied max block crossing number
!        (the maximum number of times an IO proton can move between opposite
!        walls of the cell) and the maximum number of cells in either of the
!        3D directions.
!
!
  pi_IOmaxPointsPerBlock = pi_IOmaxBlockCrossingNumber * max (NXB,NYB,NZB)

  if (pi_IOmaxPointsPerBlock < 1) then
      call Driver_abortFlash ("pi_IOinit: pi_IOmaxPointsPerBlock is < 1")
  end if
!
!
!     ...From the number of wanted IO protons, determine the IO proton write
!        modulo base, i.e. all proton with tags multiple of this base will be
!        plotted. Determine the maximum number of IO protons per processor needed,
!        using the maximum number of protons per processor.
!
!
  pi_IOprotonWriteModulo = pi_totalProtons2Launch / pi_IOnumberOfProtons2Plot
  pi_IOprotonWriteModulo = max (1,pi_IOprotonWriteModulo)
  pi_IOmaxProtonCount    = pi_maxProtonCount / pi_IOprotonWriteModulo + 1
  pi_IOmaxProtonCount    = min (pi_maxProtonCount, pi_IOmaxProtonCount) 
!
!
!     ...Allocate the needed IO arrays.
!
!
  allocate (pi_IOprotonTags       (                          1:pi_IOmaxProtonCount      ))
  allocate (pi_IOprotonPointCount (                          1:pi_IOmaxProtonCount      ))
  allocate (pi_IOprotonPoints     (1:pi_IOmaxPointsPerBlock, 1:pi_IOmaxProtonCount, MDIM))
!
!
!    ...Ready!
!
!
  return
end subroutine pi_IOinit
