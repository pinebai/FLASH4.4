!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/method/Tabulated/op_writeSpeciesROTable
!!
!! NAME
!!
!!  op_writeSpeciesROTable
!!
!! SYNOPSIS
!!
!!  call op_writeSpeciesROTable (integer (in) :: fileUnit,
!!                               integer (in) :: species)
!!
!! DESCRIPTION
!!
!!  Prints out the tabulated Rosseland Opacities for the current species.
!!
!! ARGUMENTS
!!
!!  fileUnit    : unit # for the output file
!!  species     : the species index
!!
!!***
subroutine op_writeSpeciesROTable (fileUnit,species)

  use Driver_interface,  ONLY : Driver_abortFlash

  use Opacity_data,      ONLY : op_nEnergyGroups,             &
                                op_energyGroupBoundaries

  use op_tabulatedData,  ONLY : op_nstepsDensityRO,           &
                                op_nstepsTemperatureRO,       &
                                op_tabulatedEnergyBoundaries, &
                                op_tableDensityRO,            &
                                op_tableTemperatureRO,        &
                                op_RosselandTables,           &
                                op_species2ROTableIndex

  implicit none

  integer, intent (in) :: fileUnit
  integer, intent (in) :: species

  integer :: b,d,g,t
  integer :: blocks
  integer :: dstart,dend
  integer :: indexRO
  integer :: nstepsDensity
  integer :: nstepsTemperature
!
!
!   ...Get the RO table location index. 
!
!
  indexRO = op_species2ROTableIndex (species)

  if (indexRO == 0) then
      call Driver_abortFlash ('[op_writeSpeciesROTable] ERROR: no index to RO tables')
  end if
!
!
!   ...Get the current temperature and density grid.
!
  nstepsTemperature = op_nstepsTemperatureRO (indexRO)
  nstepsDensity     = op_nstepsDensityRO     (indexRO)
!
!
!   ...Print out the RO opacity table in columns of 10 for each group. 
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

        write (fileUnit,'(7X,A9,10(2X,ES12.2,2X))') 'Density -',(op_tableDensityRO (d,indexRO),d=dstart,dend)
        write (fileUnit,'(1X,A11)') 'Temperature'
        write (fileUnit,'(1X,A11)') '     |     '

        do t = 1,nstepsTemperature
           write (fileUnit,'(1X,ES9.2,6X,10(2X,ES14.6))') op_tableTemperatureRO (t,indexRO), &
                                                         (op_RosselandTables (t,d,g,indexRO),d=dstart,dend)
        end do
     end do

     dstart = dend + 1

     if (dstart <= nstepsDensity) then

         write (fileUnit,'(7X,A9,10(2X,ES12.2,2X))') 'Density -',(op_tableDensityRO (d,indexRO),d=dstart,nstepsDensity)
         write (fileUnit,'(1X,A11)') 'Temperature'
         write (fileUnit,'(1X,A11)') '     |     '

         do t = 1,nstepsTemperature
            write (fileUnit,'(1X,ES9.2,6X,10(2X,ES14.6))') op_tableTemperatureRO (t,indexRO), &
                                                          (op_RosselandTables (t,d,g,indexRO),d=dstart,nstepsDensity)
         end do

     end if

  end do
!
!
!   ...Ready! 
!
!
  return
end subroutine op_writeSpeciesROTable
