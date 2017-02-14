!!****if* source/physics/materialProperties/Opacity/OpacityMain/unitTest/Opacity_unitTest
!!
!! NAME
!!
!!  Opacity_unitTest
!!
!! SYNOPSIS
!!
!!  Opacity_unitTest (integer (in)    :: fileUnit,
!!                    logical (inout) :: perfect)
!!
!! DESCRIPTION
!!
!!  This is a unitTest setup for testing the Opacity unit. Normally called
!!  from Driver_evolveFlash within a Simulation unitTest. See for example
!!  the directory 'source/Simulation/SimulationMain/unitTest/Opacity'.
!!
!!  Performs a test run on several of the opacity subunits: constant,
!!  tabulated, etc
!!
!! ARGUMENTS
!!
!!  fileUnit : number of file unit for diagnostic output
!!  perfect  : will be set .true., if the test is correct.
!!             This is always set .true. in this implementation.
!!
!! NOTES
!!
!!  In the Simulation unit, you should include in your Config file the line
!!    REQUIRES physics/materialProperties/Opacity/OpacityMain/unitTest
!!  in order to enable this test.
!!
!!
!!***
subroutine Opacity_unitTest (fileUnit, perfect)

  use Driver_interface,     ONLY : Driver_abortFlash
  use Opacity_interface,    ONLY : Opacity
  use Grid_interface,       ONLY : Grid_getListOfBlocks,    &
                                   Grid_getBlkPtr,          &
                                   Grid_releaseBlkPtr,      &
                                   Grid_getBlkIndexLimits
  use Opacity_data,         ONLY : op_speciesWeights
  use op_numericsData,      ONLY : op_Avogadro

  implicit none

# include "constants.h"
# include "Flash.h"

  integer, intent (in)    :: fileUnit
  logical, intent (inout) :: perfect

  character (len=9) :: opacityKind

  integer :: blkCount
  integer :: block
  integer :: blockID
  integer :: cell
  integer :: cellIndex
  integer :: firstCellIndex
  integer :: totalCells

  real    :: cellmassDensity
  real    :: cellmassFraction
  real    :: cellNumberDensity
  real    :: cellTemperature
  real    :: opacityAbsorption
  real    :: opacityEmission
  real    :: opacityTransport

  integer :: n
  real    :: Output,T

  integer, dimension (MAXBLOCKS) :: blkList

  integer, dimension (2,MDIM)    :: blkLimits
  integer, dimension (2,MDIM)    :: blkLimitsGC

  real, pointer, dimension (:,:,:,:) :: solnData
!
!
!   ...Temporarily test some stuff.
!
!  
!  write (*,*) ' Opacity kind ? '
!  read  (*,*) opacityKind
!
!  opacityKind = "Planck"
!  opacityKind = "Rosseland"
!
!  do n = 1000000,1000000
!     T = real (n)
!     call op_KleinGroupOpacity (opacityKind,1,T,510000.,510010.,Output)
!     call op_BiggsGroupOpacity (opacityKind,1,T,10.,11.,Output)
!     write (*,*) ' T = ',T,' Opacity= ',Output
!  end do
!
!  return
!
!
!   ...Get and loop over all leaf blocks.
!
!  
  call Grid_getListOfBlocks (LEAF, blkList, blkCount)

  do block = 1,blkCount

     blockID = blkList (block)
     call Grid_getBlkPtr         (blockID, solnData)
     call Grid_getBlkIndexLimits (blockID, blkLimits, blkLimitsGC)

     firstCellIndex = blkLimits (LOW,IAXIS)

     totalCells = blkLimits (HIGH,IAXIS) - blkLimits (LOW,IAXIS) + 1

     do cell = 1,totalCells

        cellIndex        = firstCellIndex + cell - 1
        cellmassDensity  = solnData (     DENS_VAR,cellIndex,1,1)
        cellmassFraction = solnData (SPECIES_BEGIN,cellIndex,1,1)
        cellTemperature  = solnData (     TELE_VAR,cellIndex,1,1)

        cellNumberDensity = cellmassFraction * cellmassDensity * op_Avogadro / op_speciesWeights (1)

        write (*,*) ' cell # ',cell
        write (*,*) ' ---------------------------------------- '
        write (*,*) '         mass density (     g/cm^3) in cell  = ',cellmassDensity
        write (*,*) '       number density (# ions/cm^3) in cell  = ',cellNumberDensity
        write (*,*) '        mass fraction               in cell  = ',cellmassFraction
        write (*,*) ' electron temperature (          K) in cell  = ',cellTemperature
!
!
!   ...Get the tabulated opacities.
!
!  
        call Opacity (solnData (1:NUNK_VARS,cellIndex,1,1), &
                      1,                          &
                      opacityAbsorption,          &
                      opacityEmission,            &
                      opacityTransport            )

        write (*,*) ' opacityAbsorption (tabulated) = ',opacityAbsorption
        write (*,*) ' opacityEmission   (tabulated) = ',opacityEmission
        write (*,*) ' opacityTransport  (tabulated) = ',opacityTransport

     end do

     call Grid_releaseBlkPtr (blockID, solnData)

  end do

  perfect = .true.
!
!
!    ...Ready!
!
!
  return
end subroutine Opacity_unitTest
