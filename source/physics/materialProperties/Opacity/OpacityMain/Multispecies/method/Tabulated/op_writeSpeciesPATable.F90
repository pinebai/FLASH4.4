!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/method/Tabulated/op_writeSpeciesPATable
!!
!! NAME
!!
!!  op_writeSpeciesPATable
!!
!! SYNOPSIS
!!
!!  call op_writeSpeciesPATable (integer (in) :: fileUnit,
!!                               integer (in) :: species)
!!
!! DESCRIPTION
!!
!!  Prints out the tabulated Planck Absorption Opacities for the current species.
!!
!! ARGUMENTS
!!
!!  fileUnit    : unit # for the output file
!!  species     : the species index
!!
!!***
subroutine op_writeSpeciesPATable (fileUnit,species)

  use Driver_interface,  ONLY : Driver_abortFlash

  use Opacity_data,      ONLY : op_nEnergyGroups,             &
                                op_energyGroupBoundaries

  use op_tabulatedData,  ONLY : op_nstepsDensityPA,           &
                                op_nstepsTemperaturePA,       &
                                op_tabulatedEnergyBoundaries, &
                                op_tableDensityPA,            &
                                op_tableTemperaturePA,        &
                                op_PlanckAbsorptionTables,    &
                                op_species2PATableIndex

  implicit none

  integer, intent (in) :: fileUnit
  integer, intent (in) :: species

  integer :: b,d,g,t
  integer :: blocks
  integer :: dstart,dend
  integer :: indexPA
  integer :: nstepsDensity
  integer :: nstepsTemperature
!
!
!   ...Get the PA table location index. 
!
!
  indexPA = op_species2PATableIndex (species)

  if (indexPA == 0) then
      call Driver_abortFlash ('[op_writeSpeciesPATable] ERROR: no index to PA tables')
  end if
!
!
!   ...Get the current temperature and density grid.
!
  nstepsTemperature = op_nstepsTemperaturePA (indexPA)
  nstepsDensity     = op_nstepsDensityPA     (indexPA)
!
!
!   ...Print out the PA opacity table in columns of 10 for each group. 
!
!
  blocks = nstepsDensity / 10

  do g = 1,op_nEnergyGroups

     write (fileUnit,*)
     write (fileUnit,*)
     write (fileUnit,*)
     write (fileUnit,*)
     write (fileUnit,*) ' Energy Group ',op_tabulatedEnergyBoundaries (g), &
                                    'to',op_tabulatedEnergyBoundaries (g+1),'(eV)'
     write (fileUnit,*)

     dend = 0
     do b = 1,blocks
        dstart = dend + 1
        dend = dend + 10

        write (fileUnit,'(7X,A9,10(2X,ES12.2,2X))') 'Density -',(op_tableDensityPA (d,indexPA),d=dstart,dend)
        write (fileUnit,'(1X,A11)') 'Temperature'
        write (fileUnit,'(1X,A11)') '     |     '

        do t = 1,nstepsTemperature
           write (fileUnit,'(1X,ES9.2,6X,10(2X,ES14.6))') op_tableTemperaturePA (t,indexPA), &
                                                         (op_PlanckAbsorptionTables (t,d,g,indexPA),d=dstart,dend)
        end do
     end do

     dstart = dend + 1

     if (dstart <= nstepsDensity) then

         write (fileUnit,'(7X,A9,10(2X,ES12.2,2X))') 'Density -',(op_tableDensityPA (d,indexPA),d=dstart,nstepsDensity)
         write (fileUnit,'(1X,A11)') 'Temperature'
         write (fileUnit,'(1X,A11)') '     |     '

         do t = 1,nstepsTemperature
            write (fileUnit,'(1X,ES9.2,6X,10(2X,ES14.6))') op_tableTemperaturePA (t,indexPA), &
                                                          (op_PlanckAbsorptionTables (t,d,g,indexPA),d=dstart,nstepsDensity)
         end do

     end if

  end do
!
!
!   ...Ready! 
!
!
  return
end subroutine op_writeSpeciesPATable
