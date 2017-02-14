!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/method/Constcm2g/op_writeConstcm2g
!!
!! NAME
!!
!!  op_writeConstcm2g
!!
!! SYNOPSIS
!!
!!  call op_writeConstcm2g ()
!!
!! DESCRIPTION
!!
!!  Prints out the constant (in units of cm^2/g) Opacities for each
!!  species. Only those constants are printed which were actually
!!  generated for each species.
!!
!! ARGUMENTS
!!
!!***
subroutine op_writeConstcm2g ()

  use Opacity_data,    ONLY : op_totalSpecies,       &
                              op_absorptionKind,     &
                              op_emissionKind,       &
                              op_transportKind

  use op_constcm2gData, ONLY : op_absorptionConstcm2g, &
                               op_emissionConstcm2g,   &
                               op_transportConstcm2g

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
        file = "opacity_printout_constcm2g.txt", &
        form = 'formatted')
!
!
!   ...Print out the title. 
!
!
  write (fileUnit,*)
  write (fileUnit,*) "   OPACITY CONSTANTS PRINTOUT (in units of cm^2/g)"
  write (fileUnit,*)
!
!
!   ...Loop over all species. 
!
!
  do species = 1,op_totalSpecies

     needConstant =      (op_absorptionKind (species) == OP_CONSTCM2G) &
                    .or. (op_emissionKind   (species) == OP_CONSTCM2G) &
                    .or. (op_transportKind  (species) == OP_CONSTCM2G)

     if (needConstant) then

         write (fileUnit,*)
         write (fileUnit,*) ' ------- Species # ',species,' -------'
         write (fileUnit,*)

         if (op_absorptionKind (species) == OP_CONSTCM2G) then
             write (fileUnit,'(1X,A22,ES14.6)') 'Absorption Constant = ', &
                                                 op_absorptionConstcm2g (species)
         end if

         if (op_emissionKind (species) == OP_CONSTCM2G) then
             write (fileUnit,'(1X,A22,ES14.6)') 'Emission   Constant = ', &
                                                 op_emissionConstcm2g (species)
         end if

         if (op_transportKind (species) == OP_CONSTCM2G) then
             write (fileUnit,'(1X,A22,ES14.6)') 'Transport  Constant = ', &
                                                 op_transportConstcm2g (species)
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
end subroutine op_writeConstcm2g
