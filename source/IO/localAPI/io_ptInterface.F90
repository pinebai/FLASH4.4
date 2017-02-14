!!****ih* source/IO/localAPI/io_ptInterface
!!
!! NAME
!!  io_ptInterface
!!
!! SYNOPSIS
!!  use io_ptInterface
!!
!! DESCRIPTION
!!
!! This is the interface module for particle IO.  This handles whatever
!! correctoions need to be placed in checkpoint files to the next particle 
!! output time, as well as initialzation of particles.  These must only be
!! called from within IO.
!! 
!!***

#include "constants.h"
#include "Flash.h"

module io_ptInterface

  interface 
     subroutine io_ptCorrectNextPartTime(simTime, oldTime)
        real, intent(IN) :: simTime
        real, intent(out) :: oldTime
     end subroutine io_ptCorrectNextPartTime
  end interface

  interface
     subroutine io_ptResetNextFile(savedNext)
       real, intent(IN) :: savedNext
     end subroutine io_ptResetNextFile
  end interface

  interface
     subroutine io_ptInit()
       
       implicit none 
       
     end subroutine io_ptInit
  end interface

  interface
     subroutine io_ptSendOutputData()
       implicit none

     end subroutine io_ptSendOutputData
  end interface

  interface
     subroutine io_ptWriteParticleData( fileID, &
          globalNumParticles, localNumParticles, particleOffset, &
          partAttributeLabels, particlesToCheckpoint)
       implicit none
       integer, intent(IN) :: fileID, globalNumParticles, &
            localNumParticles, particleOffset
       character (len=OUTPUT_PROP_LENGTH), intent(IN) :: partAttributeLabels(NPART_PROPS)
       logical, intent(IN):: particlesToCheckpoint
     end subroutine io_ptWriteParticleData
  end interface

  interface
     subroutine io_ptReadParticleData()
       implicit none
     end subroutine io_ptReadParticleData
  end interface
       
  interface
     subroutine io_ptCreateSubset ( subsetIndex, numProperties, &
          numParticles, partLabels, particlest, &
          subsetSize, subsetName, subsetLabelName, subsetLabels, &
          moreSubsetsRemain)   
       implicit none
       integer, intent(IN) :: subsetIndex, numProperties, &
            numParticles
       character(len=OUTPUT_PROP_LENGTH), dimension(numProperties), &
            intent(IN) :: partLabels
       real, dimension(numProperties,numParticles), intent(INOUT) :: particlest
       integer, dimension(2), intent(OUT) :: subsetSize
       character(len=MAX_STRING_LENGTH), intent(OUT) :: subsetName, subsetLabelName
       character(len=OUTPUT_PROP_LENGTH), dimension(numProperties), &
            intent(OUT) :: subsetLabels
       logical, intent(OUT) :: moreSubsetsRemain
     end subroutine io_ptCreateSubset
  end interface

end module io_ptInterface
