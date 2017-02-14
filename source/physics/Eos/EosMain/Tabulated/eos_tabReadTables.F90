!!****if* source/physics/Eos/EosMain/Tabulated/eos_tabReadTables
!!
!! NAME
!!
!!  eos_tabReadTables
!!
!! SYNOPSIS
!!
!!  call eos_tabReadTables (character (in) :: tableKind (len=80),
!!                      character (in) :: tableName (len=80),
!!                      character (in) :: groupName (len=80),
!!                      logical   (in) :: wanted(EOS_TAB_NCOMP,EOS_TAB_NALLTAB),
!!                 type(eosT_tableGroupDescT)(INOUT) :: td(:), 
!!                        tbZF,tbEN,tbPR,tbHC,tbEntr)
!!
!! DESCRIPTION
!!
!!  Reads in the necessary data for processing and interpolating tabulated
!!  opacities from specific tables. The routine calls appropriate subroutines
!!  according to the kind of EOS tables specified. Currently only the
!!  IONMIX opacities can be read in.
!!
!! ARGUMENTS
!!
!!  tableKind   : the kind of tabulated Opacity file where data is going to be read
!!  tableName   : the name of tabulated Opacity file where data is going to be read
!!
!!***

#include "constants.h"

subroutine eos_tabReadTables (tableKind,   &
                              tableName,   &
                              groupName,   &
                              wanted,      &
                              td, &
                              tbZF,tbEN,tbPR,tbHC,tbEntr)

  use Driver_interface,   ONLY : Driver_abortFlash
  use eos_tabInterface, ONLY : eos_tabReadIonmixTables, &
                               eos_tabReadIonmix4Tables

  use eos_tabData,ONLY: EOS_TAB_NCOMP,EOS_TAB_NALLTAB, &
                        EOS_TAB_FOR_ELE, &
                        EOS_TAB_FOR_MAT, &
                        EOS_TABVT_EN, &
                        EOS_TABVT_ZF, &
                        EOS_TABVT_HC, &
                        EOS_TABVT_PR, &
                        EOS_TABVT_ENTR, &
                        eosT_tableGroupDescT, &
                        eosT_oneVarTablePT
                             

  implicit none

  logical,            intent (in) :: wanted(EOS_TAB_NCOMP,EOS_TAB_NALLTAB)
  type(eosT_tableGroupDescT),intent(inout) :: td(:)
  type(eosT_oneVarTablePT),pointer,dimension(:) :: tbZF,tbEN,tbPR,tbHC,tbEntr
  character (len=80), intent (in) :: tableKind
  character (len=80), intent (in) :: tableName
  character (len=80), intent (in) :: groupName

!
!   ...Call the appropriate routine.
!
!
  if (tableKind == "IONMIX") then

      call eos_tabReadIonmixTables (tableName,   &
                                wanted, &
                                td, &
                                tbZF,tbEN,tbHC)

  else if (tableKind=="IONMIX4" .OR. tableKind=="IONMIX6") then


      call eos_tabReadIonmix4Tables (tableName,   &
                                wanted, &
                                td, &
                                tbZF,tbEN,tbPR,tbHC,tbEntr)

  else if (tableKind == "OPACPLOT") then

      call eos_tabReadOpacplotTables (tableName,   &
                                groupName,   &
                                wanted, &
                                td, &
                                tbZF,tbEN,tbPR,tbHC,tbEntr)

  else
      call Driver_abortFlash ('[eos_tabReadTables] ERROR: EOS table kind not recognized')
  end if
!
!
!   ...Ready! 
!
!
  return
end subroutine eos_tabReadTables
