!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/method/Constant/op_writeConstants
!!
!! NAME
!!
!!  op_writeConstants
!!
!! SYNOPSIS
!!
!!  call op_writeConstants ()
!!
!! DESCRIPTION
!!
!!  Prints out the constant Opacities for each species. Only those constants
!!  are printed which were actually generated for each species.
!!
!! ARGUMENTS
!!
!!***
subroutine op_writeConstants ()

  use Opacity_data,    ONLY : op_totalSpecies,       &
                              op_absorptionKind,     &
                              op_emissionKind,       &
                              op_transportKind

  use op_constantData, ONLY : op_absorptionConstant, &
                              op_emissionConstant,   &
                              op_transportConstant

  implicit none

#include "Opacity.h"
#include "constants.h"

  logical :: needConstant
  integer :: fileUnit, ut_getFreeFileUnit
  integer :: species
!
!
!   ...Open the printout file. 
!
!
  fileUnit = ut_getFreeFileUnit ()

  open (unit = fileUnit, &
        file = "opacity_printout_constants.txt", &
        form = 'formatted')
!
!
!   ...Print out the title. 
!
!
  write (fileUnit,*)
  write (fileUnit,*) "   OPACITY CONSTANTS PRINTOUT (in units of 1/g)"
  write (fileUnit,*)
!
!
!   ...Loop over all species. 
!
!
  do species = 1,op_totalSpecies

     needConstant =      (op_absorptionKind (species) == OP_CONSTANT) &
                    .or. (op_emissionKind   (species) == OP_CONSTANT) &
                    .or. (op_transportKind  (species) == OP_CONSTANT)

     if (needConstant) then

         write (fileUnit,*)
         write (fileUnit,*) ' ------- Species # ',species,' -------'
         write (fileUnit,*)

         if (op_absorptionKind (species) == OP_CONSTANT) then
             write (fileUnit,'(1X,A22,ES14.6)') 'Absorption Constant = ', &
                                                 op_absorptionConstant (species)
         end if

         if (op_emissionKind (species) == OP_CONSTANT) then
             write (fileUnit,'(1X,A22,ES14.6)') 'Emission   Constant = ', &
                                                 op_emissionConstant (species)
         end if

         if (op_transportKind (species) == OP_CONSTANT) then
             write (fileUnit,'(1X,A22,ES14.6)') 'Transport  Constant = ', &
                                                 op_transportConstant (species)
         end if

     end if

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
end subroutine op_writeConstants
