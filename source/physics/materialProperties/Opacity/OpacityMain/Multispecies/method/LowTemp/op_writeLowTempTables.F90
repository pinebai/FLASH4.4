!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/method/LowTemp/op_writeLowTempTables
!!
!! NAME
!!
!!  op_writeLowTempTables
!!
!! SYNOPSIS
!!
!!  call op_writeLowTempTables ()
!!
!! DESCRIPTION
!!
!!  Prints out the tabulated low temperature Planck and Rosseland Opacities for each species.
!!
!! ARGUMENTS
!!
!!***
subroutine op_writeLowTempTables ()

  use Opacity_data,   ONLY : op_nEnergyGroups,         &
                             op_totalSpecies,          &
                             op_energyGroupBoundaries

  use op_lowTempData, ONLY : op_maxNstepsLowTemp,      &
                             op_tableLowTemp,          &
                             op_PlanckLowTempTables,   &
                             op_RosselandLowTempTables

  implicit none

  integer :: b,g,s,t
  integer :: blocks
  integer :: fileUnit, ut_getFreeFileUnit
  integer :: group
  integer :: sFirst,sLast
!
!
!    ...Open the printout file.
!
!
  fileUnit = ut_getFreeFileUnit ()

  open (unit = fileUnit, &
        file = "opacity_printout_lowTemp_tables.txt", &
        form = 'formatted')
!
!
!   ...Print out the title for the tables. 
!
!
  write (fileUnit,*)
  write (fileUnit,*) "   OPACITY LOW TEMPERATURE TABLES PRINTOUT (Opacities    in cm^2/g)"
  write (fileUnit,*) "                                           (Temperatures in Kelvin)"
  write (fileUnit,*)
!
!
!   ...Set the block structure of the tables in species columns of 10. 
!
!
  blocks = op_totalSpecies / 10
!
!
!   ...Printout the low temperature Planck Opacities. 
!
!
  write (fileUnit,*)
  write (fileUnit,*) "   LOW TEMPERATURE PLANCK OPACITIES"
  write (fileUnit,*)

  do g = 1,op_nEnergyGroups

     write (fileUnit,*)
     write (fileUnit,*)
     write (fileUnit,*)
     write (fileUnit,*)
     write (fileUnit,*) ' Energy Group ',op_energyGroupBoundaries (g), &
                                    'to',op_energyGroupBoundaries (g+1),'(eV)'
     write (fileUnit,*)

     sLast = 0

     do b = 1,blocks
        sFirst = sLast + 1
        sLast = sLast + 10

        write (fileUnit,'(7X,A9,10(7X,I3,7X))') 'Species -',(s,s=sFirst,sLast)
        write (fileUnit,'(1X,A11)') 'Temperature'
        write (fileUnit,'(1X,A11)') '     |     '

        do t = 1,op_maxNstepsLowTemp
           write (fileUnit,'(1X,ES9.2,6X,10(2X,ES14.6))') op_tableLowTemp (t), &
                                                         (op_PlanckLowTempTables (t,s,g),s=sFirst,sLast)
        end do
     end do

     sFirst = sLast + 1
     sLast  = op_totalSpecies

     if (sFirst <= sLast) then

         write (fileUnit,'(7X,A9,10(7X,I3,7X))') 'Species -',(s,s=sFirst,sLast)
         write (fileUnit,'(1X,A11)') 'Temperature'
         write (fileUnit,'(1X,A11)') '     |     '

         do t = 1,op_maxNstepsLowTemp
            write (fileUnit,'(1X,ES9.2,6X,10(2X,ES14.6))') op_tableLowTemp (t) , &
                                                          (op_PlanckLowTempTables (t,s,g),s=sFirst,sLast)
         end do

     end if

  end do
!
!
!   ...Printout the low temperature Rosseland Opacities. 
!
!
  write (fileUnit,*)
  write (fileUnit,*)
  write (fileUnit,*)
  write (fileUnit,*)
  write (fileUnit,*) "   LOW TEMPERATURE ROSSELAND OPACITIES"
  write (fileUnit,*)

  do g = 1,op_nEnergyGroups

     write (fileUnit,*)
     write (fileUnit,*)
     write (fileUnit,*)
     write (fileUnit,*)
     write (fileUnit,*) ' Energy Group ',op_energyGroupBoundaries (g), &
                                    'to',op_energyGroupBoundaries (g+1),'(eV)'
     write (fileUnit,*)

     sLast = 0

     do b = 1,blocks
        sFirst = sLast + 1
        sLast = sLast + 10

        write (fileUnit,'(7X,A9,10(7X,I3,7X))') 'Species -',(s,s=sFirst,sLast)
        write (fileUnit,'(1X,A11)') 'Temperature'
        write (fileUnit,'(1X,A11)') '     |     '

        do t = 1,op_maxNstepsLowTemp
           write (fileUnit,'(1X,ES9.2,6X,10(2X,ES14.6))') op_tableLowTemp (t), &
                                                         (op_RosselandLowTempTables (t,s,g),s=sFirst,sLast)
        end do
     end do

     sFirst = sLast + 1
     sLast  = op_totalSpecies

     if (sFirst <= sLast) then

         write (fileUnit,'(7X,A9,10(7X,I3,7X))') 'Species -',(s,s=sFirst,sLast)
         write (fileUnit,'(1X,A11)') 'Temperature'
         write (fileUnit,'(1X,A11)') '     |     '

         do t = 1,op_maxNstepsLowTemp
            write (fileUnit,'(1X,ES9.2,6X,10(2X,ES14.6))') op_tableLowTemp (t) , &
                                                          (op_RosselandLowTempTables (t,s,g),s=sFirst,sLast)
         end do

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
end subroutine op_writeLowTempTables
