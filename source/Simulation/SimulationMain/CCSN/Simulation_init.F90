!!****if* source/Simulation/SimulationMain/CCSN/Simulation_init
!!
!! NAME
!!
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!
!!  Simulation_init()
!!
!! DESCRIPTION
!!
!!  Initialization routine for CCSN problem.  Reads in 
!!  runtime parameters.  Reads 1D model file and maps
!!  variable names to FLASH unk names.
!!
!! ARGUMENTS
!!
!!   
!!
!! PARAMETERS
!!
!! NOTES
!!  
!!  This problem is described in, e.g.,
!!  Couch, S.M. 2013, ApJ, 765, 29
!!  Couch, S.M. 2013, ApJ, 775, 35
!!  Couch, S.M. & O'Connor, E.P. 2013, arXiv:1310.5728
!!
!!***
subroutine Simulation_init()
  use Simulation_data
  use Simulation_interface, ONLY : Simulation_mapIntToStr
  use Driver_interface, ONLY : Driver_getMype
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Logfile_interface, ONLY : Logfile_stampMessage
  use Grid_interface, ONLY : Grid_getGeometry

  implicit none

#include "constants.h"
#include "Flash.h"

  integer :: i
  integer :: geom

  integer, parameter :: max_stored_vars = 30
  real var_temp(max_stored_vars)

  logical :: dr_restart

  character (len=256) :: current_line

  integer :: j, ipos, NUNK_VARS_stored
!  integer :: var_key(max_stored_vars)
  integer :: var_key (NUNK_VARS)
  character (len=4) :: var_labels(max_stored_vars)

  call Driver_getMype(MESH_COMM, sim_meshMe)

  call Grid_getGeometry(geom)

  call RuntimeParameters_get( 'model_file',   model_file)
  call RuntimeParameters_get( 'smlrho',   sim_smlrho)
  call RuntimeParameters_get( 'smallt',   sim_smallt)
  call RuntimeParameters_get( 'smallx',   sim_smallx)
  call RuntimeParameters_get( 'vel_mult', sim_velMult)
  call RuntimeParameters_get( 'nsub', nsub)
  call RuntimeParameters_get( 'restart', sim_restart)

  do i = UNK_VARS_BEGIN, UNK_VARS_END
     call Simulation_mapIntToStr(i, unklabels(i), MAPBLOCK_UNK)
     call makeLowercase(unklabels(i))
  enddo

  ! open the file and read in the header 
  open(unit=2,file=model_file,status='old')
  read (2,'(a80)') current_line

  if (sim_meshMe == MASTER_PE) print *, 'file opened'

  ! read in the number of variables line
  read (2,'(a80)') current_line
  ipos = index(current_line,'=') + 1
  read (current_line(ipos:),*) nvar_stored
  if (sim_meshMe == MASTER_PE) print *,"read nvar_stored", nvar_stored

  if (NUNK_VARS .NE. nvar_stored .AND. sim_meshMe == MASTER_PE) then
     print *, ' '
     print *, 'Warning: the number of variables stored in the'
     print *, 'input file is different than the number of'
     print *, 'variables in the current version of FLASH.'
     print *, ' '
     print *, 'The variables in the file that are also defined'
     print *, 'in FLASH will be read in.  Any missing variables'
     print *, 'will be initialized to zero'
     print *, ' '
  endif

  if (sim_meshMe == MASTER_PE) then
     print *, "Vaiables in file:"
  endif
  do i = 1, nvar_stored
     read (2,'(a4)') var_labels(i)
     if (sim_meshMe == MASTER_PE) &
          print *, var_labels(i)
     call makeLowercase(var_labels(i))
  enddo

  do j = 1, NUNK_VARS
     var_key(j) = NONEXISTENT

     do i = 1, nvar_stored

        if (unklabels(j) == var_labels(i)) then
           var_key(j) = i
        endif

     enddo

     if (var_key(j) == NONEXISTENT) then
        if(sim_meshMe == MASTER_PE) then
           print *, 'Warning, variable: ', unklabels(j), & 
                ' not found in the input file.'
           print *, 'initializing ', unklabels(j), ' to 0'
           print *, ' '
        endif
     endif

  enddo

  do i = 1, n1d_max
     read(2,*,end=11) xzn(i), (var_temp(j),j=1,nvar_stored)

     ! put these in order, so model1d_var always contains the same variables 
     ! in the same spots
     do j = 1, NUNK_VARS
        if (var_key(j) /= NONEXISTENT) then
           model_1d(i,j) = var_temp(var_key(j))
        else
           model_1d(i,j) = 0.0
        endif

     enddo

  enddo

11 close(unit=2)

  n1d_total = 0

  do while (xzn(n1d_total+1) .NE. 0.00)
     n1d_total = n1d_total + 1
  enddo

  if (sim_meshMe .EQ. MASTER_PE) then
     print *, 'file read completed'
     print *, n1d_total, 'points read in'
  endif

end subroutine Simulation_init
