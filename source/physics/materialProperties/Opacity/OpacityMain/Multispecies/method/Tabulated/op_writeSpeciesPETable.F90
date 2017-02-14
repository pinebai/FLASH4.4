!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/method/Tabulated/op_writeSpeciesPETable
!!
!! NAME
!!
!!  op_writeSpeciesPETable
!!
!! SYNOPSIS
!!
!!  call op_writeSpeciesPETable (integer (in) :: fileUnit,
!!                               integer (in) :: species)
!!
!! DESCRIPTION
!!
!!  Prints out the tabulated Planck Emission Opacities for the current species.
!!
!! ARGUMENTS
!!
!!  fileUnit    : unit # for the output file
!!  species     : the species index
!!
!!***
subroutine op_writeSpeciesPETable (fileUnit,species)

  use Driver_interface,  ONLY : Driver_abortFlash

  use Opacity_data,      ONLY : op_nEnergyGroups,             &
                                op_energyGroupBoundaries

  use op_tabulatedData,  ONLY : op_nstepsDensityPE,           &
                                op_nstepsTemperaturePE,       &
                                op_tabulatedEnergyBoundaries, &
                                op_tableDensityPE,            &
                                op_tableTemperaturePE,        &
                                op_PlanckEmissionTables,      &
                                op_species2PETableIndex

  implicit none

  integer, intent (in) :: fileUnit
  integer, intent (in) :: species

  integer :: b,d,g,t
  integer :: blocks
  integer :: dstart,dend
  integer :: indexPE
  integer :: nstepsDensity
  integer :: nstepsTemperature
!
!
!   ...Get the PE table location index. 
!
!
  indexPE = op_species2PETableIndex (species)

  if (indexPE == 0) then
      call Driver_abortFlash ('[op_writeSpeciesPETable] ERROR: no index to PE tables')
  end if
!
!
!   ...Get the current temperature and density grid.
!
  nstepsTemperature = op_nstepsTemperaturePE (indexPE)
  nstepsDensity     = op_nstepsDensityPE     (indexPE)
!
!
!   ...Print out the PE opacity table in columns of 10 for each group. 
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

        write (fileUnit,'(7X,A9,10(2X,ES12.2,2X))') 'Density -',(op_tableDensityPE (d,indexPE),d=dstart,dend)
        write (fileUnit,'(1X,A11)') 'Temperature'
        write (fileUnit,'(1X,A11)') '     |     '

        do t = 1,nstepsTemperature
           write (fileUnit,'(1X,ES9.2,6X,10(2X,ES14.6))') op_tableTemperaturePE (t,indexPE), &
                                                         (op_PlanckEmissionTables (t,d,g,indexPE),d=dstart,dend)
        end do
     end do

     dstart = dend + 1

     if (dstart <= nstepsDensity) then

         write (fileUnit,'(7X,A9,10(2X,ES12.2,2X))') 'Density -',(op_tableDensityPE (d,indexPE),d=dstart,nstepsDensity)
         write (fileUnit,'(1X,A11)') 'Temperature'
         write (fileUnit,'(1X,A11)') '     |     '

         do t = 1,nstepsTemperature
            write (fileUnit,'(1X,ES9.2,6X,10(2X,ES14.6))') op_tableTemperaturePE (t,indexPE), &
                                                          (op_PlanckEmissionTables (t,d,g,indexPE),d=dstart,nstepsDensity)
         end do

     end if

  end do
!
!
!   ...Ready! 
!
!
  return
end subroutine op_writeSpeciesPETable
