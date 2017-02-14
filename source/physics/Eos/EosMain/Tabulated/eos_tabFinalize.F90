!!****if* source/physics/Eos/EosMain/Tabulated/eos_tabFinalize
!!
!! NAME
!!
!!  eos_tabFinalize
!!
!! SYNOPSIS
!!
!!  call eos_tabFinalize ()
!!
!! DESCRIPTION
!!
!!  Clean up the Opacity unit.
!!
!! ARGUMENTS
!!
!!***
#include "Flash.h"

subroutine eos_tabFinalize ()

  use eos_tabData,   ONLY : eos_tableKind,              &
                             eos_tableName,              &
                             eos_groupName,              &
                             eos_tabIonizationKind,         &
                             eos_tabIntEnergyKind,           &
                             eos_tabHeatCpKind,          &
           eos_allTab, &
           eos_tabAllDiag, &
                             EOS_TAB_NALLTAB
  implicit none

  integer :: species, i, j

  do species = 1,NSPECIES
     do i = 1,EOS_TAB_NALLTAB
        ! Currently, the only way these pointers get different from null is by being allocated. - KW
        if (associated(eos_allTab(species)%tg(i)%table)) then
           do j = LBOUND(eos_allTab(species)%tg(i)%table,1), &
                UBOUND(eos_allTab(species)%tg(i)%table,1)
              if (associated(eos_allTab(species)%tg(i)%table(j)%table)) deallocate(eos_allTab(species)%tg(i)%table(j)%table)
           end do
           deallocate(eos_allTab(species)%tg(i)%table)
        end if
        if (associated(eos_allTab(species)%tg(i)%mgTable)) deallocate(eos_allTab(species)%tg(i)%mgTable)
        if (associated(eos_allTab(species)%tg(i)%td%Temperatures)) deallocate(eos_allTab(species)%tg(i)%td%Temperatures)
        if (associated(eos_allTab(species)%tg(i)%td%Densities)) deallocate(eos_allTab(species)%tg(i)%td%Densities)
     end do
  end do
  deallocate(eos_allTab)
  deallocate(eos_tabAllDiag)

  deallocate (eos_tableKind)
  deallocate (eos_tableName)
  deallocate (eos_tabIonizationKind)
  deallocate (eos_tabIntEnergyKind)
  deallocate (eos_tabHeatCpKind)

end subroutine eos_tabFinalize
