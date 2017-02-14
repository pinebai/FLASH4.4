!!****if* source/diagnostics/ProtonImaging/localAPI/pi_printProtonsData
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

  implicit none
   
  character (len=*), intent (in) :: fileLabel
  integer,           intent (in) :: processorID

  return
end subroutine pi_printProtonsData
