!!****if* source/physics/Eos/EosMain/Tabulated/eos_tabWriteSpeciesENTable
!!
!! NAME
!!
!!  eos_tabWriteSpeciesENTable
!!
!! SYNOPSIS
!!
!!  call eos_tabWriteSpeciesENTable (integer (in) :: fileUnit,
!!                               integer (in) :: species)
!!
!! DESCRIPTION
!!
!!  Prints out the tabulated internal energy data for the current species.
!!
!! ARGUMENTS
!!
!!  fileUnit    : unit # for the output file
!!  species     : the species index
!!
!!***

#include "Eos.h"


subroutine eos_tabWriteSpeciesENTable (fileUnit,species)

  use Driver_interface,  ONLY : Driver_abortFlash
  use Logfile_interface,  ONLY : Logfile_stamp
  use Eos_Data,           ONLY : eos_meshMe
  use eos_tabData,      ONLY : EOS_TAB_NCOMP,          &
                               EOS_TAB_FOR_MAT, &
                               EOS_TAB_FOR_ION, &
                               EOS_TAB_FOR_ELE, &
                               EOS_TABVT_EN, &
                                eos_allTab, &
                                eosT_tableGroupDescT, &
                                eosT_oneVarTablePT

  implicit none

  integer, intent (in) :: fileUnit
  integer, intent (in) :: species

  integer :: b,d,g,t
  integer :: blocks
  integer :: dstart,dend

  integer :: nstepsDensity
  integer :: nstepsTemperature

  type(eosT_tableGroupDescT),pointer :: td
  type(eosT_oneVarTablePT),pointer :: tb
!
!
!   ...Check whether storage for table data has been allocated. 
!
!


  if (.NOT.associated(eos_allTab(species)%tg(EOS_TABVT_EN)%table)) then
      call Driver_abortFlash ('[eos_tabWriteSpeciesENTable] ERROR: no EN data tables')
  end if
!
!
!   ...Get the current temperature and density grid.
!
  td => eos_allTab(species)%tg(EOS_TABVT_EN)%td
  nstepsTemperature = td%ntemp
  nstepsDensity     = td%ndens
!
!
!   ...Print out the EN table in columns of 10 for each group. 
!
!
  blocks = nstepsDensity / 10

  if (nstepsDensity*nstepsTemperature .LE. 0) then
#ifdef DEBUG_EOS
     call Logfile_stamp("Zero-size EN tata tables???")
#endif
     return !RETURN IMMEDIATELY
  end if


  do g = LBOUND(eos_allTab(species)%tg(EOS_TABVT_EN)%table,1), &
       UBOUND(eos_allTab(species)%tg(EOS_TABVT_EN)%table,1)


     tb => eos_allTab(species)%tg(EOS_TABVT_EN)%table(g)
     if (.NOT.associated(tb%table)) then
#ifdef DEBUG_EOS
        if (eos_meshMe==0) &
             print*,'eos_tabWriteSpeciesENTable: .NOT.associated(tb%table) for g=', g,', species=',species
#endif
        cycle
     end if

#ifdef DEBUG_EOS
     print*,'LBOUND(tb%table,1),LBOUND(tb%table,2):',LBOUND(tb%table,1),LBOUND(tb%table,2)
     print*,'UBOUND(tb%table,1),UBOUND(tb%table,2):',UBOUND(tb%table,1),UBOUND(tb%table,2)
#endif

     write (fileUnit,*)
     if (g==EOS_TAB_FOR_MAT) then
        write (fileUnit,*) '  --  All Matter --'
     else
        select case (g)
        case(EOS_TAB_FOR_ION)
           write (fileUnit,*) '   -- Ions --'
        case(EOS_TAB_FOR_ELE)
           write (fileUnit,*) '   -- Electrons --'
        end select
     end if
     write (fileUnit,*)

     dend = 0
     do b = 1,blocks
        dstart = dend + 1
        dend = dend + 10

        write (fileUnit,'(7X,A9,1P,10(G13.4,3X))') 'Density -',(td%Densities (d),d=dstart,dend)
        write (fileUnit,'(1X,A11)') 'Temperature'
        write (fileUnit,'(1X,A11)') '     |     '

        do t = 1,nstepsTemperature
           write (fileUnit,'(1X,1P,G13.7,2X,10(2X,G14.7))') td%Temperatures (t), &
                                                         (tb%table (t,d),d=dstart,dend)
        end do
     end do

     dstart = dend + 1

     if (dstart <= nstepsDensity) then

         write (fileUnit,'(7X,A9,1P,10(G13.4,3X))') 'Density -',(td%Densities (d),d=dstart,nstepsDensity)
         write (fileUnit,'(1X,A11)') 'Temperature'
         write (fileUnit,'(1X,A11)') '     |     '

         do t = 1,nstepsTemperature
            write (fileUnit,'(2X,1P,G8.2,6X,10(2X,G14.7))') td%Temperatures (t), &
                                                          (tb%table (t,d),d=dstart,nstepsDensity)
         end do

     end if

  end do
!
!
!   ...Ready! 
!
!
  return
end subroutine eos_tabWriteSpeciesENTable
