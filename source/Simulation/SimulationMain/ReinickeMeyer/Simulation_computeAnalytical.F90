!!****if* source/Simulation/SimulationMain/ReinickeMeyer/Simulation_computeAnalytical
!!
!! NAME
!!
!!  Simulation_computeAnalytical
!!
!!
!! SYNOPSIS
!!
!!  call Simulation_computeAnalytical(integer(IN) :: blockID, 
!!                                    real(IN)    :: tcurr)
!!
!! DESCRIPTION
!!
!!  Compute an analytical solution for the conducting, 1D blast
!!  wave. Solution described in:
!!
!!  Reinicke, P., Meyer-ter-Vehn, J., The point explosion with heat conduction, 
!!  Phys. Fluids A, 1807, 3, 1991
!!  
!! This simulation works by calling the rmtv subroutine, which takes
!! kind of weird units, so there is some conversion of units which has
!! to go on.
!! 
!! The rmtv subroutine takes the following parameters:
!!
!! rpos    = desired radial position for the solution in cm
!! aval_in = power in thermal conductivity chi0 * rho**a * T**b
!! bval_in = power in thermal conductivity chi0 * rho**a * T**b
!! chi0    = coefficient in thermal conductivity chi0 * rho**a * T**b
!! gamma   = ratio of specific heats
!! bigamma = Gruneisen coefficient (gamma-1)*ener = pres/den = G*temp
!! rf      = position of the heat front in cm
!! xif     = dimensionless position of the heat front
!! xis     = dimensionless position of the shock front
!! beta0   = eigenvalue of the problem
!! g0      = heat front scaling parameter
!!
!! output:
!! den  = density  g/cm**3
!! tev  = temperature ev
!! ener = specific internal energy erg/g
!! pres = presssure erg/cm**3
!! vel  = velocity cm/s
!! 
!! ARGUMENTS
!!
!!  blockID -        the number of the block to initialize
!!  tcurr   -        current time
!!
!! SIDE EFFECTS
!!
!!  The analytical solution is computed and stored in an appropriate slot
!!  (or slots) in the solution vector, UNK.
!!
!! NOTES
!! 
!!  By default this is just a stub that doesn't do anything.
!!
!!***

subroutine Simulation_computeAnalytical(blockId, tcurr)

  use Grid_interface, ONLY: Grid_getBlkIndexLimits
  use Grid_interface, ONLY: Grid_getCellCoords
  use Grid_interface, ONLY: Grid_putPointData

  use Simulation_data, ONLY : sim_gamma
  use Simulation_data, ONLY : sim_boltz
  use Simulation_data, ONLY : sim_avo
  use Simulation_data, ONLY : sim_singleSpeciesA
  use Simulation_data, ONLY : sim_densityExponent
  use Simulation_data, ONLY : sim_temperatureExponent
  use Simulation_data, ONLY : sim_K0
  use Simulation_data, ONLY : sim_tinitial
  use Simulation_data, ONLY : sim_rfInit
  use Simulation_data, ONLY : sim_smallt
  use Simulation_data, ONLY : sim_smlrho

  implicit none

#include "constants.h"
#include "Flash.h"

  integer, intent(in) :: blockId
  real,    intent(in) :: tcurr

  ! RMTV_1D SUBROUTINE INPUT:
  real :: aval_in
  real :: bval_in
  real :: chi0
  real :: gamma
  real :: bigamma
  real :: rf
  real :: xif_in
  real :: xis
  real :: beta0_in
  real :: g0

  ! RMTV_1D SUBROUTINE OUTPUT:
  real :: rmtv_rho   = 0.
  real :: rmtv_tev   = 0.
  real :: rmtv_eint  = 0.
  real :: rmtv_pres  = 0.
  real :: rmtv_vel   = 0.
  real :: rmtv_velx  = 0.
  real :: rmtv_vely  = 0.
  real :: rmtv_velz  = 0.
  real :: rmtv_ek    = 0.

  integer :: i
  integer :: j
  integer :: k
  integer :: sizeX
  integer :: sizeY
  integer :: sizeZ
  integer :: blkLimits(2,MDIM)
  integer :: blkLimitsGC(2,MDIM)
  integer :: axis(MDIM)

  real :: radius
  real :: zeta
  real :: alpha
  real :: xx 
  real :: yy 
  real :: zz

  real, allocatable :: xCoord(:)
  real, allocatable :: yCoord(:)
  real, allocatable :: zCoord(:)

  ! UNIT CONVERSIONS:

  ! Conversion from ergs to jerks:
  real, parameter :: ERG_TO_JERK = 1.0 / 1.0e+16

  ! Conversion from Kelvin to keV:
  real, parameter :: K_TO_KEV = 8.617332478e-05 / 1000.0

  ! Conversion from eV to Kelvin:
  real, parameter :: EV_TO_KELVIN = 11604.505

  ! Conversion from seconds to shakes:
  real, parameter :: S_TO_SHAKE = 1.0e+08
  

  ! *************************************************
  ! *                                               *
  ! *     SET INPUT PARAMETERS FOR CALL TO RMTV     *
  ! *                                               *
  ! *************************************************

  ! A, B, and gamma are specified in the runtime parameters file:
  aval_in = sim_densityExponent
  bval_in = sim_temperatureExponent
  gamma   = sim_gamma

  ! Technically, these can also be varied, however, I can't seem to
  ! make a reasonable simulation work with anything other than the
  ! values below. For this reason, I thought that xif_in, xis, and
  ! beta0_in should not be varied until they are better understood...
  xif_in   = 2.0
  xis      = 1.0
  beta0_in = 71975339.999999911

  ! chi0 has units of:
  ! 
  ! energy/time/length/temperature**(bval_in+1)/density**(aval_in)
  ! 
  ! We have to convert from FLASH units of:
  !
  ! ergs/s/cm/K**(bval_in+1) / (g/cc)**aval_in
  !
  ! to:
  !
  ! jerks / shakes / cm / keV**(bval_in+1) / (g/cc)**aval_in
  !
  chi0 = sim_K0
  chi0 = chi0 * ERG_TO_JERK / S_TO_SHAKE / (K_TO_KEV)**(bval_in+1)

  ! Bigamma has units of energy / mass / temperature. We have to
  ! convert from FLASH units of:
  !
  ! ergs / g / K
  !
  ! to:
  !
  ! jerks / g / keV:
  ! 
  bigamma = sim_avo * sim_boltz / sim_singleSpeciesA
  bigamma = bigamma * ERG_TO_JERK / K_TO_KEV

  ! Compute zeta, which is invarient throughout the simulation
  ! (Section IV from paper):
  alpha = (2*bval_in - 2*aval_in + 1) / (2*bval_in - 5*aval_in + 3)
  zeta = sim_rfInit / xif_in / (sim_tinitial*S_TO_SHAKE)**alpha
  rf = xif_in * zeta * (tcurr*S_TO_SHAKE)**alpha

  ! Compute g0 (Eqn 10 from paper):
  g0 = 2*chi0 * (alpha*zeta**(1/alpha))**(2*bval_in - 1) / &
       bigamma**(bval_in+1) / beta0_in
  g0 = g0**(1/(1-aval_in))

  ! *********************************
  ! *                               *
  ! *     SET ANALYTIC SOLUTION     *
  ! *                               *
  ! *********************************

  ! Get cell center coordinates:
  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)
  
  sizeX = blkLimitsGC(HIGH,IAXIS) - blkLimitsGC(LOW,IAXIS) + 1
  allocate(xCoord(sizeX)); xCoord = 0.0
  
  sizeY = blkLimitsGC(HIGH,JAXIS) - blkLimitsGC(LOW,JAXIS) + 1
  allocate(yCoord(sizeY)); yCoord = 0.0
  
  sizeZ = blkLimitsGC(HIGH,KAXIS) - blkLimitsGC(LOW,KAXIS) + 1
  allocate(zCoord(sizeZ)); zCoord = 0.0

  call Grid_getCellCoords(IAXIS, blockId, CENTER, .true., xCoord, sizeX)

  if (NDIM >= 2) then
     call Grid_getCellCoords(JAXIS, blockId, CENTER, .true., yCoord, sizeY)
  end if
  
  if (NDIM == 3) then
     call Grid_getCellCoords(KAXIS, blockId, CENTER, .true., zCoord, sizeZ)
  endif
  
  ! Loop over cell in block and set the initial condition in each cell
  ! using the output from the rmtv subroutine...
  do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
     do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)    
        do i = blkLimits(LOW,IAXIS), blkLimits(HIGH, IAXIS)
           
           xx = 0.0
           yy = 0.0
           zz = 0.0
           
           xx = xCoord(i)
           if (NDIM >= 2) yy = yCoord(j)
           if (NDIM == 3) zz = zCoord(k)

           radius = sqrt( xx**2 + yy**2 + zz**2)
           
           call rmtv_1d(radius,                                      &
                aval_in,bval_in,chi0,gamma,bigamma,                  &
                rf,xif_in,xis,beta0_in,g0,                           &
                rmtv_rho,rmtv_tev,rmtv_eint,rmtv_pres,rmtv_vel)
           
           rmtv_rho = max(rmtv_rho, sim_smlrho)
           rmtv_tev = max(rmtv_tev * EV_TO_KELVIN, sim_smallt)
           
           rmtv_velX = rmtv_vel*xx/radius
           rmtv_vely = rmtv_vel*yy/radius
           rmtv_velz = rmtv_vel*zz/radius
           
           axis(IAXIS)=i
           axis(JAXIS)=j
           axis(KAXIS)=k

           call Grid_putPointData(blockId, CENTER, ARHO_VAR, EXTERIOR, axis, rmtv_rho)
           call Grid_putPointData(blockId, CENTER, ATMP_VAR, EXTERIOR, axis, rmtv_tev)
           call Grid_putPointData(blockId, CENTER, AVLX_VAR, EXTERIOR, axis, rmtv_velx)
           call Grid_putPointData(blockId, CENTER, AVLY_VAR, EXTERIOR, axis, rmtv_vely)
           call Grid_putPointData(blockId, CENTER, AVLZ_VAR, EXTERIOR, axis, rmtv_velz)
           
        enddo
     enddo
  enddo  
  
  deallocate(xCoord)
  deallocate(yCoord)
  deallocate(zCoord)
  
  return  
  
end subroutine Simulation_computeAnalytical
