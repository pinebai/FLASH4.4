!!****if* source/diagnostics/ProtonImaging/ProtonImagingMain/oldOneTimeStep/ProtonUtilities/pi_printProtonsData
!!
!! NAME
!!
!!  pi_printProtonsData
!!
!! SYNOPSIS
!!
!!  call pi_printProtonsData (character (len=*), intent (in) :: fileLabel,
!!                            integer,           intent (in) :: processorID)
!!
!! DESCRIPTION
!!
!!  Utility routine, which prints detailed info regarding the generated data for all protons
!!  on the current processor to a text file. The information is written out to a file named
!!  <basenm><fileLabel><processorID>.txt, where <basenm> is the runtime parameter for output
!!  file names, <fileLabel> the chosen label of the file and <processorID> is the current
!!  processor. The use of different labels for <fileLabel> allows printing of the proton data
!!  at different stages of the simulation, without losing previous printed proton data.
!!
!! ARGUMENTS
!!
!!  fileLabel   : the label of the printout file
!!  processorID : processor identification number
!!
!!***

subroutine pi_printProtonsData (fileLabel, processorID)

  use ProtonImaging_data,  ONLY : pi_baseName,    &
                                  pi_protonCount, &
                                  pi_protons

  implicit none
   
#include "Flash.h"
#include "ProtonImaging.h"
#include "constants.h"

  character (len=*), intent (in) :: fileLabel
  integer,           intent (in) :: processorID

  character (len = 4                ) :: charPID
  character (len = MAX_STRING_LENGTH) :: fileName

  integer :: beam
  integer :: blockID
  integer :: detector
  integer :: fileUnit
  integer :: proton
  integer :: protonTag
  integer :: ut_getFreeFileUnit

  real    :: posX, posY, posZ
  real    :: velX, velY, velZ
!
!
!   ...Immediate return, if no protons present.
!
!
  if (pi_protonCount < 1) then
      return
  end if
!
!
!   ...Open the (label and processor specific) printout file.
!
!
  write (charPID,'(I4.4)') processorID  ! this transforms a number N into 000N

  fileUnit = ut_getFreeFileUnit ()
  fileName = trim (pi_baseName) // trim (fileLabel) // charPID // ".txt"

  open (fileUnit, file = fileName)
!
!
!     ...Print out the protons position.
!
!
  write (fileUnit,'(/)')
  write (fileUnit,'(a)') "                         PROTON POSITIONS"
  write (fileUnit,'(/)')
  write (fileUnit,'(a)') " Proton Nr               X                Y                Z"
  write (fileUnit,'(a)') " --------------------------------------------------------------------"

  do proton = 1, pi_protonCount

     posX = pi_protons (PROTON_POSX,proton)
     posY = pi_protons (PROTON_POSY,proton)
     posZ = pi_protons (PROTON_POSZ,proton)

     write (fileUnit,'(i6,3x,3es20.12)') proton , posX , posY , posZ

  end do
!
!
!     ...Print out the protons velocity.
!
!
  write (fileUnit,'(/)')
  write (fileUnit,'(a)') "                         PROTON VELOCITIES"
  write (fileUnit,'(/)')
  write (fileUnit,'(a)') " Proton Nr               X                Y                Z"
  write (fileUnit,'(a)') " --------------------------------------------------------------------"

  do proton = 1, pi_protonCount

     velX = pi_protons (PROTON_VELX,proton)
     velY = pi_protons (PROTON_VELY,proton)
     velZ = pi_protons (PROTON_VELZ,proton)

     write (fileUnit,'(i6,3x,3es20.12)') proton , velX , velY , velZ

  end do
!
!
!     ...Print out the protons block ID number and globally unique tag.
!
!
  write (fileUnit,'(/)')
  write (fileUnit,'(a)') "  PROTON BLOCK ID  &  GLOBAL TAG"
  write (fileUnit,'(/)')
  write (fileUnit,'(a)') " Proton Nr  Block ID  Global Tag"
  write (fileUnit,'(a)') " -------------------------------"

  do proton = 1, pi_protonCount

     blockID   = int (pi_protons (PROTON_BLCK,proton))
     protonTag = int (pi_protons (PROTON_TAGS,proton))

     write (fileUnit,'(i6,3x,i6,i10)') proton , blockID, protonTag

  end do
!
!
!     ...Print out the protons beam number from where they originated and
!        their target detector number.
!
!
  write (fileUnit,'(/)')
  write (fileUnit,'(a)') "  PROTON BEAM NUMBER  &  DETECTOR NUMBER"
  write (fileUnit,'(/)')
  write (fileUnit,'(a)') " Proton Nr   Beam #   Detector #"
  write (fileUnit,'(a)') " -------------------------------"

  do proton = 1, pi_protonCount

     beam     = int (pi_protons (PROTON_BEAM,proton))
     detector = int (pi_protons (PROTON_DETC,proton))

     write (fileUnit,'(i6,3x,i6,5x,i6)') proton , beam, detector

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
end subroutine pi_printProtonsData
