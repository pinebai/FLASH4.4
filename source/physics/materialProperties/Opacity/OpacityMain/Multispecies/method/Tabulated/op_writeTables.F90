!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/method/Tabulated/op_writeTables
!!
!! NAME
!!
!!  op_writeTables
!!
!! SYNOPSIS
!!
!!  call op_writeTables ()
!!
!! DESCRIPTION
!!
!!  Prints out the tabulated Opacities for each species. Only those tables
!!  are printed which were actually generated for each species.
!!
!! ARGUMENTS
!!
!!***
subroutine op_writeTables ()

  use Opacity_data,     ONLY : op_totalSpecies,   &
                               op_absorptionKind, &
                               op_emissionKind,   &
                               op_transportKind

  use op_tabulatedData, ONLY : op_useLogTables,          &
                               op_speciesMaxTempPATable, &
                               op_speciesMaxTempPETable, &
                               op_speciesMaxTempROTable, &
                               op_speciesMinTempPATable, &
                               op_speciesMinTempPETable, &
                               op_speciesMinTempROTable

  use op_interface,     ONLY : op_writeSpeciesTables

  implicit none

#include "Opacity.h"
#include "constants.h"

  logical :: needTable
  logical :: needPATable
  logical :: needPETable
  logical :: needROTable
  integer :: fileUnit
  integer :: ut_getFreeFileUnit
  integer :: species
!
!
!    ...Open the printout file.
!
!
  fileUnit = ut_getFreeFileUnit ()

  open (unit = fileUnit, &
        file = "opacity_printout_tables.txt", &
        form = 'formatted')
!
!
!   ...Print out the title for the table temperature boundaries. 
!
!
  write (fileUnit,*)
  write (fileUnit,*) "   OPACITY TABLES TEMPERATURE BOUNDARIES (in units of Kelvin, -1 means no info)"
  write (fileUnit,*)
  write (fileUnit,*)

  write (fileUnit,*)
  write (fileUnit,*) "       SPECIES #    PLANCK ABSORPTION min T   PLANCK ABSORPTION max T"

  do species = 1,op_totalSpecies
     write (fileUnit,'(10X,I3,11X,ES14.6,12X,ES14.6)') species,op_speciesMinTempPATable (species), &
                                                               op_speciesMaxTempPATable (species)
  end do

  write (fileUnit,*)
  write (fileUnit,*) "       SPECIES #      PLANCK EMISSION min T     PLANCK EMISSION max T"

  do species = 1,op_totalSpecies
     write (fileUnit,'(10X,I3,11X,ES14.6,12X,ES14.6)') species,op_speciesMinTempPETable (species), &
                                                               op_speciesMaxTempPETable (species)
  end do

  write (fileUnit,*)
  write (fileUnit,*) "       SPECIES #            ROSSELAND min T           ROSSELAND max T"

  do species = 1,op_totalSpecies
     write (fileUnit,'(10X,I3,11X,ES14.6,12X,ES14.6)') species,op_speciesMinTempROTable (species), &
                                                               op_speciesMaxTempROTable (species)
  end do
!
!
!   ...Print out the title for the tables. 
!
!
  write (fileUnit,*)
  write (fileUnit,*)
  write (fileUnit,*)
  write (fileUnit,*)

  if (op_useLogTables) then

      write (fileUnit,*) "   OPACITY TABLES PRINTOUT (Opacities    in log10 (cm^2/g))"
      write (fileUnit,*) "                           (Temperatures in log10 (Kelvin)"
      write (fileUnit,*) "                           (Densities    in log10 (#ions/cm^3))"

  else

      write (fileUnit,*) "   OPACITY TABLES PRINTOUT (Opacities    in cm^2/g)"
      write (fileUnit,*) "                           (Temperatures in Kelvin)"
      write (fileUnit,*) "                           (Densities    in #ions/cm^3)"

  end if

  write (fileUnit,*)
!
!
!   ...Loop over all species. 
!
!
  do species = 1,op_totalSpecies

     needPATable  =      (op_absorptionKind (species) == OP_TABULAR_PA) &
                    .or. (op_emissionKind   (species) == OP_TABULAR_PA) &
                    .or. (op_transportKind  (species) == OP_TABULAR_PA)

     needPETable  =      (op_absorptionKind (species) == OP_TABULAR_PE) &
                    .or. (op_emissionKind   (species) == OP_TABULAR_PE) &
                    .or. (op_transportKind  (species) == OP_TABULAR_PE)

     needROTable  =      (op_absorptionKind (species) == OP_TABULAR_RO) &
                    .or. (op_emissionKind   (species) == OP_TABULAR_RO) &
                    .or. (op_transportKind  (species) == OP_TABULAR_RO)

     needTable    =       needPATable &
                    .or.  needPETable &
                    .or.  needROTable

     if (needTable) then

         write (fileUnit,*)
         write (fileUnit,*) ' ------- Species # ',species,' -------'
         write (fileUnit,*)

         call op_writeSpeciesTables (fileUnit,    &
                                     species,     &
                                     needPATable, &
                                     needPETable, &
                                     needROTable  )
     end if

  end do
!
!
!    ...Close the printout file.
!
!
  close (fileUnit)
!
!
!    ...Ready!
!
!
  return
end subroutine op_writeTables
