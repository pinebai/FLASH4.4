!!****if* source/Simulation/SimulationMain/unitTest/Gravity/BHTree-layer/Simulation_init
!!
!! NAME
!!
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!
!!  call Simulation_init()
!!
!! ARGUMENTS
!!
!!    none
!!
!! DESCRIPTION
!!
!!  Initializes all the data specified in Simulation_data.
!!  It calls the RuntimeParameters_get routine for initialization.
!!
!!
!!***

subroutine Simulation_init()
  
  use Simulation_data
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get
  use Driver_interface, ONLY : Driver_abortFlash, Driver_getMype, Driver_getComm

  implicit none
#include "Flash.h"
#include "constants.h"
#include "Flash_mpi.h"


  integer :: i, istat
  real :: xmax, xmin, ymax, ymin, zmax, zmin

  call Driver_getMype(GLOBAL_COMM, sim_MyPE)
  call Driver_getComm(GLOBAL_COMM, sim_Comm)

  call RuntimeParameters_get( 'sim_prof_file' , sim_prof_file )
  call RuntimeParameters_get( 'sim_zMidplane', sim_zMidplane )
  call RuntimeParameters_get( 'sim_dir', sim_dir )

  call RuntimeParameters_get( 'smlrho', smlrho)
  call RuntimeParameters_get( 'smallp', smallp)
  call RuntimeParameters_get( 'smallX', smallX)
  call RuntimeParameters_get( 'gamma_1', sim_gamma_1)
  call RuntimeParameters_get( 'xmax', xmax)
  call RuntimeParameters_get( 'xmin', xmin)
  call RuntimeParameters_get( 'ymax', ymax)
  call RuntimeParameters_get( 'ymin', ymin)
  call RuntimeParameters_get( 'zmax', zmax)
  call RuntimeParameters_get( 'zmin', zmin)
  call RuntimeParameters_get( 'sim_solutionErrorTolerance1', sim_solutionErrorTolerance1)
  call RuntimeParameters_get( 'sim_solutionErrorTolerance2', sim_solutionErrorTolerance2)
  sim_Lx = xmax - xmin
  sim_Ly = ymax - ymin
  sim_Lz = zmax - zmin
  

  call PhysicalConstants_get( 'Boltzmann', sim_boltz)
  call PhysicalConstants_get( 'proton mass', sim_mH)
  call PhysicalConstants_get( 'pi', sim_pi)

  sim_vecLen = 1
  sim_mode = MODE_DENS_PRES

  ! read number of lines in the radial profile file and communicate it
  if (sim_myPE .eq. MASTER_PE) then
    open(unit = 55, file = sim_prof_file, status = 'old')
    read(55,*) sim_nProfile
  endif
  call MPI_Bcast(sim_nProfile, 1, MPI_INTEGER, MASTER_PE, MPI_COMM_WORLD, istat)

  ! read radial profile from file
  ! & construct the radial samples needed for the initialization.
  allocate(sim_zProf(sim_nProfile), stat=istat)
  if (istat .ne. 0) call Driver_abortFlash("Could not allocate sim_zProf")
  allocate(sim_rhoProf(sim_nProfile), stat=istat)
  if (istat .ne. 0) call Driver_abortFlash("Could not allocate sim_rhoProf")
  allocate(sim_pProf(sim_nProfile), stat=istat)
  if (istat .ne. 0) call Driver_abortFlash("Could not allocate sim_pProf")
  allocate(sim_vProf(sim_nProfile), stat=istat)
  if (istat .ne. 0) call Driver_abortFlash("Could not allocate sim_vProf")
  allocate(sim_xProf(sim_nProfile), stat=istat)
  if (istat .ne. 0) call Driver_abortFlash("Could not allocate sim_xProf")
  allocate(sim_PhiProf(sim_nProfile), stat=istat)
  if (istat .ne. 0) call Driver_abortFlash("Could not allocate sim_PhiProf")

  if (sim_myPE .eq. MASTER_PE) then
    do i = 1,sim_nProfile
      read(55,*) sim_zProf(i), sim_rhoProf(i), sim_pProf(i), sim_vProf(i) &
      & , sim_xProf(i), sim_PhiProf(i)
    enddo
    endfile(unit = 55)
  endif

  call MPI_Bcast(sim_zProf,   sim_nProfile, MPI_DOUBLE_PRECISION, MASTER_PE, MPI_COMM_WORLD, istat)
  call MPI_Bcast(sim_rhoProf, sim_nProfile, MPI_DOUBLE_PRECISION, MASTER_PE, MPI_COMM_WORLD, istat)
  call MPI_Bcast(sim_pProf,   sim_nProfile, MPI_DOUBLE_PRECISION, MASTER_PE, MPI_COMM_WORLD, istat)
  call MPI_Bcast(sim_vProf,   sim_nProfile, MPI_DOUBLE_PRECISION, MASTER_PE, MPI_COMM_WORLD, istat)
  call MPI_Bcast(sim_xProf,   sim_nProfile, MPI_DOUBLE_PRECISION, MASTER_PE, MPI_COMM_WORLD, istat)
  call MPI_Bcast(sim_PhiProf, sim_nProfile, MPI_DOUBLE_PRECISION, MASTER_PE, MPI_COMM_WORLD, istat)

end subroutine Simulation_init
