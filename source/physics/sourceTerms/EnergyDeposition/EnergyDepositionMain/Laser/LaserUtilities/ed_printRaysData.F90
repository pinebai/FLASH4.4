!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserUtilities/ed_printRaysData
!!
!! NAME
!!
!!  ed_printRaysData
!!
!! SYNOPSIS
!!
!!  call ed_printRaysData (integer (in) :: processorID)
!!
!! DESCRIPTION
!!
!!  Utility routine, which prints detailed info regarding the generated data for all rays
!!  on the current processor to a text file. The information is written out to a file named
!!  <basenm>LaserRaysDataPrint<processorID>.txt, where <basenm> is the runtime parameter for
!!  output file names and <processorID> is the current processor.
!!
!! ARGUMENTS
!!
!!  processorID : processor identification number
!!
!!***

subroutine ed_printRaysData (processorID)

  use EnergyDeposition_data,  ONLY : ed_baseName, &
                                     ed_rayCount, &
                                     ed_rays

  implicit none
   
#include "Flash.h"
#include "EnergyDeposition.h"
#include "constants.h"

  integer, intent (in) :: processorID

  character (len = 4                ) :: charPID
  character (len = MAX_STRING_LENGTH) :: fileName

  integer :: blockID
  integer :: fileUnit
  integer :: ray
  integer :: rayTag
  integer :: ut_getFreeFileUnit

  real    :: power
  real    :: rayX, rayY, rayZ
  real    :: velX, velY, velZ
!
!
!   ...Open the (processor specific) printout file.
!
!
  write (charPID,'(I4.4)') processorID

  fileUnit = ut_getFreeFileUnit ()
  fileName = trim (ed_baseName) // "LaserRaysPrint" // charPID // ".txt"

  open (fileUnit, file = fileName)
!
!
!     ...Print out the rays position.
!
!
  write (fileUnit,'(/)')
  write (fileUnit,'(a)') "                         RAY POSITIONS"
  write (fileUnit,'(/)')
  write (fileUnit,'(a)') " Ray Nr               X                Y                Z"
  write (fileUnit,'(a)') " -----------------------------------------------------------------"

  do ray = 1, ed_rayCount

     rayX = ed_rays (RAY_POSX,ray)
     rayY = ed_rays (RAY_POSY,ray)
     rayZ = ed_rays (RAY_POSZ,ray)

     write (fileUnit,'(i6,3es20.12)') ray , rayX , rayY , rayZ

  end do
!
!
!     ...Print out the rays velocity.
!
!
  write (fileUnit,'(/)')
  write (fileUnit,'(a)') "                         RAY VELOCITIES"
  write (fileUnit,'(/)')
  write (fileUnit,'(a)') " Ray Nr               X                Y                Z"
  write (fileUnit,'(a)') " -----------------------------------------------------------------"

  do ray = 1, ed_rayCount

     velX = ed_rays (RAY_VELX,ray)
     velY = ed_rays (RAY_VELY,ray)
     velZ = ed_rays (RAY_VELZ,ray)

     write (fileUnit,'(i6,3es20.12)') ray , velX , velY , velZ

  end do
!
!
!     ...Print out the rays block ID number and globally unique tag.
!
!
  write (fileUnit,'(/)')
  write (fileUnit,'(a)') "  RAY BLOCK ID  &  GLOBAL TAG"
  write (fileUnit,'(/)')
  write (fileUnit,'(a)') " Ray Nr  Block ID  Global Tag"
  write (fileUnit,'(a)') " ----------------------------"

  do ray = 1, ed_rayCount

     blockID = int (ed_rays (RAY_BLCK,ray))
     rayTag  = int (ed_rays (RAY_TAGS,ray))

     write (fileUnit,'(i6,i6,i10)') ray , blockID, rayTag

  end do
!
!
!     ...Print out the rays power.
!
!
  write (fileUnit,'(/)')
  write (fileUnit,'(a)') "        RAY POWERS"
  write (fileUnit,'(/)')
  write (fileUnit,'(a)') " Ray Nr            Power"
  write (fileUnit,'(a)') " ----------------------------"

  do ray = 1, ed_rayCount

     power = ed_rays (RAY_POWR,ray)

     write (fileUnit,'(i6,es20.12)') ray , power

  end do
!
!
!   ...Close the printout file.
!
!
  close (fileUnit)
!
!
!    ...Ready!
!
!  
  return
end subroutine ed_printRaysData
