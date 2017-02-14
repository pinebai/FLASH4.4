!!****if* source/Simulation/SimulationMain/unitTest/Gravity/BHTree/Simulation_init
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
  use sim_interface, ONLY : sim_LPMNS

  implicit none
#include "Flash.h"
#include "constants.h"
#include "Flash_mpi.h"


  integer :: i, istat
  real :: xmax, xmin, ymax, ymin, zmax, zmin
  integer :: MAXSPHL
  parameter (MAXSPHL = 200)
  real :: PM(0:MAXSPHL), PD(0:MAXSPHL)


  call Driver_getMype(GLOBAL_COMM, sim_MyPE)
  call Driver_getComm(GLOBAL_COMM, sim_Comm)


  call RuntimeParameters_get( 'sim_radprof_file' , sim_radprof_file )
  call RuntimeParameters_get( 'sim_xCenter',       sim_xCenter )
  call RuntimeParameters_get( 'sim_yCenter',       sim_yCenter )
  call RuntimeParameters_get( 'sim_zCenter',       sim_zCenter )
  call RuntimeParameters_get( 'sim_vx',            sim_vx )
  call RuntimeParameters_get( 'sim_vy',            sim_vy )
  call RuntimeParameters_get( 'sim_vz',            sim_vz )
  call RuntimeParameters_get( 'sim_pertType',      sim_pertType )
  call RuntimeParameters_get( 'sim_pertamp',       sim_pertamp )
  call RuntimeParameters_get( 'sim_velamp',        sim_velamp )
  call RuntimeParameters_get( 'sim_spharm_l1',     sim_spharm_l1 )
  call RuntimeParameters_get( 'sim_spharm_m1',     sim_spharm_m1 )
  call RuntimeParameters_get( 'sim_nSubZones',     sim_nSubZones )
  call RuntimeParameters_get( 'sim_solutionErrorTolerance1', sim_solutionErrorTolerance1)
  call RuntimeParameters_get( 'sim_solutionErrorTolerance2', sim_solutionErrorTolerance2)

  call RuntimeParameters_get( 'smlrho', smlrho)
  call RuntimeParameters_get( 'smallp', smallp)
  call RuntimeParameters_get( 'smallX', smallX)
  call RuntimeParameters_get( 'xmax', xmax)
  call RuntimeParameters_get( 'xmin', xmin)
  call RuntimeParameters_get( 'ymax', ymax)
  call RuntimeParameters_get( 'ymin', ymin)
  call RuntimeParameters_get( 'zmax', zmax)
  call RuntimeParameters_get( 'zmin', zmin)
  sim_Lx = xmax - xmin
  sim_Ly = ymax - ymin
  sim_Lz = zmax - zmin
  
  call RuntimeParameters_get( 'jeans_ref', jeans_ref)
  call RuntimeParameters_get( 'jeans_deref', jeans_deref)

  call PhysicalConstants_get( 'Boltzmann', sim_boltz)
  call PhysicalConstants_get( 'proton mass', sim_mH)
  call PhysicalConstants_get( 'pi', sim_pi)

  !! DEV: sim_vecLen and sim_mode are currently unused. - KW
  sim_vecLen = 1
  sim_mode = MODE_DENS_PRES

  ! read number of lines in the radial profile file and communicate it
  if (sim_myPE .eq. MASTER_PE) then
    open(unit = 55, file = sim_radprof_file, status = 'old')
    read(55,*) sim_nProfile
  endif
  call MPI_Bcast(sim_nProfile, 1, MPI_INTEGER, MASTER_PE, sim_comm, istat)

  ! read radial profile from file
  ! & construct the radial samples needed for the initialization.
  allocate(sim_rProf(sim_nProfile), stat=istat)
  if (istat .ne. 0) call Driver_abortFlash("Could not allocate sim_rProf")
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
      read(55,*) sim_rProf(i), sim_rhoProf(i), sim_pProf(i), sim_vProf(i) &
      & , sim_xProf(i), sim_PhiProf(i)
    enddo
    close(55)
  endif

  call MPI_Bcast(sim_rProf,   sim_nProfile, MPI_DOUBLE_PRECISION, MASTER_PE, sim_comm, istat)
  call MPI_Bcast(sim_rhoProf, sim_nProfile, MPI_DOUBLE_PRECISION, MASTER_PE, sim_comm, istat)
  call MPI_Bcast(sim_pProf,   sim_nProfile, MPI_DOUBLE_PRECISION, MASTER_PE, sim_comm, istat)
  call MPI_Bcast(sim_vProf,   sim_nProfile, MPI_DOUBLE_PRECISION, MASTER_PE, sim_comm, istat)
  call MPI_Bcast(sim_xProf,   sim_nProfile, MPI_DOUBLE_PRECISION, MASTER_PE, sim_comm, istat)
  call MPI_Bcast(sim_PhiProf, sim_nProfile, MPI_DOUBLE_PRECISION, MASTER_PE, sim_comm, istat)


  ! normalize associated Legendre polynom
  plgndr_max = -1d99
  plgndr_min = 1d99
  if (sim_pertType .eq. 2) then
    do i = 0, 1000
      call sim_LPMNS(sim_spharm_m1, sim_spharm_l1, i/1000., PM, PD)
      if (PM(sim_spharm_l1) .gt. plgndr_max) plgndr_max = PM(sim_spharm_l1)
      if (PM(sim_spharm_l1) .lt. plgndr_min) plgndr_min = PM(sim_spharm_l1)
      if (i .eq. 0) plgndr_norm = PM(sim_spharm_l1) ! normalization taken at theta = 0
    enddo
    !print *, "plgndr_min,max,norm = ", plgndr_min, plgndr_max, plgndr_norm
  endif

  sim_absErrMax = 0.0
  sim_relErrMax = 0.0

end subroutine Simulation_init
