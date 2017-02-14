!!****if* source/physics/materialProperties/Opacity/OpacityMain/OPAL/op_writeTables
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
!!  Prints out the tabulated Opacities for each hydrogen abundance. Only those tables
!!  are printed which were actually generated for each hydrogen abundance.
!!
!! ARGUMENTS
!!
!!***
subroutine op_writeTables ()

  use Opacity_data,     ONLY : op_opalNumHydrogenAbundances,   &
                               op_absorptionKind, &
                               op_emissionKind,   &
                               op_transportKind

  use op_opalData, ONLY : OP_OPAL_LOWT,OP_OPAL_HIGHT,    &
                               op_useLogTables

  use op_interface,     ONLY : op_writeSpeciesTables

  implicit none

#include "Opacity.h"
#include "constants.h"

  logical :: needTable
  logical :: needLowTable
  logical :: needHighTable
  logical :: needROTable
  integer :: fileUnit
  integer :: ut_getFreeFileUnit
  integer :: m
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
!   ...Loop over all hydrogen abundance. 
!
!
  do m = 1,op_opalNumHydrogenAbundances

     needLowTable  =      (op_absorptionKind (m) == OP_OPAL_LOWT) &
                    .or. (op_emissionKind   (m) == OP_OPAL_LOWT) &
                    .or. (op_transportKind  (m) == OP_OPAL_LOWT)

     needHighTable  =      (op_absorptionKind (m) == OP_OPAL_HIGHT) &
                    .or. (op_emissionKind   (m) == OP_OPAL_HIGHT) &
                    .or. (op_transportKind  (m) == OP_OPAL_HIGHT)

     needROTable  =      (op_absorptionKind (m) == OP_TABULAR_RO) &
                    .or. (op_emissionKind   (m) == OP_TABULAR_RO) &
                    .or. (op_transportKind  (m) == OP_TABULAR_RO)

     needTable    =       needLowTable &
                    .or.  needHighTable &
                    .or.  needROTable

     if (needTable) then

         write (fileUnit,*)
         write (fileUnit,*) ' ------- H Abund # ',m,' -------'
         write (fileUnit,*)

         call op_writeSpeciesTables (fileUnit,    &
                                     m,     &
                                     needLowTable, &
                                     needHighTable, &
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
