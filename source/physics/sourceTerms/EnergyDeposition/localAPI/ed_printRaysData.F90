!!****if* source/physics/sourceTerms/EnergyDeposition/localAPI/ed_printRaysData
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

  implicit none
   
  integer, intent (in) :: processorID

  return
end subroutine ed_printRaysData
