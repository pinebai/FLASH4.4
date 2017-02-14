!!****if* source/Simulation/SimulationMain/HeatexchangeIonEle/Simulation_init
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
!!
!! DESCRIPTION
!!
!!  Initializes all the parameters needed for a particular simulation
!!
!!
!! ARGUMENTS
!!
!!  
!!
!! PARAMETERS
!!
!!  sim_xctr           Temperature peak center X-coordinate
!!  sim_yctr           Temperature peak center Y-coordinate
!!  sim_zctr           Temperature peak center Z-coordinate
!!
!!   sim_maxTolCoeff0,sim_maxTolCoeff1,sim_maxTolCoeff2,sim_maxTolCoeff3
!!***

subroutine Simulation_init()
  
  use Simulation_data
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get
  use Eos_interface, ONLY : Eos
  use Logfile_interface, ONLY : Logfile_stamp
  use Driver_interface, ONLY : Driver_getMype
  implicit none
#include "Flash.h"
#include "constants.h"
#include "Eos.h"   

  

  real, dimension(:), allocatable :: eosData

  real,external :: alogam
  real,parameter ::        kev     = 8.617385e-5
  real :: n, ne ! sic
  real :: cl, KelvinPerEV
  real :: gasConstant,gammam1Inv,cv,zbar,abar
  real :: dynamicZ, relA, Ye, eMassInUAmu, ionMassInUAmu
  integer :: vecLen, case
  logical, dimension(EOS_VARS+1:EOS_NUM) :: eosMask


  call RuntimeParameters_get('orientation', sim_orientation)
  call RuntimeParameters_get('rho_init', sim_rhoInit)
  call RuntimeParameters_get('toffset', sim_toffset)
  call RuntimeParameters_get('sim_Q', sim_Q)
  call RuntimeParameters_get('sim_tempBackground', sim_tempBackground)
  call RuntimeParameters_get('sim_xctr',sim_xCenter)
  call RuntimeParameters_get('sim_yctr',sim_yCenter)
  call RuntimeParameters_get('sim_zctr',sim_zCenter)
  call RuntimeParameters_get('smallp', sim_smallP)
  call RuntimeParameters_get('smallx', sim_smallX) 
  
  call RuntimeParameters_get('smallt', sim_anaSmallT)
  call RuntimeParameters_get('sim_analytical_tolerance', sim_anaTol)
  call RuntimeParameters_get('sim_analytical_maxNewton', sim_anaMaxNewton)

  call RuntimeParameters_get('gamma', sim_gamma)
  call Driver_getMype(MESH_COMM, sim_meshMe)  

  call Logfile_stamp( "initializing Conduction problem",  &
       "[Simulation_init]")

  call RuntimeParameters_get("sim_ionTemp",sim_ionTemp)
  call RuntimeParameters_get("sim_eleTemp",sim_eleTemp)
  call RuntimeParameters_get("sim_radTemp",sim_radTemp)

  call RuntimeParameters_get("cond_TemperatureExponent", sim_condTemperatureExponent)
  call RuntimeParameters_get("initialCondTemperatureExponent", sim_initialCondTemperatureExponent)
  if (sim_initialCondTemperatureExponent==-999) then
     sim_initialCondTemperatureExponent = sim_condTemperatureExponent
  end if
  call RuntimeParameters_get("cond_K0", sim_alpha)
  call PhysicalConstants_get("ideal gas constant", gasConstant)
  gammam1Inv = 1.0/(sim_gamma-1.0)
  sim_alpha = sim_alpha / (sim_rhoInit * (gammam1Inv * gasConstant))
  n = sim_condTemperatureExponent

  vecLen = 1
  allocate(eosData(vecLen * EOS_NUM))
  ! Get thermodynamic quantities and load them on the grid using
  ! the EOS routine

!! #define CALL_EOS_HERE
#ifdef CALL_EOS_HERE

  eosMask = .FALSE.
  eosMask(EOS_CV) = .TRUE.
  eosMask(EOS_DET) = .TRUE.

  eosData(EOS_TEMPION) = sim_ionTemp
  eosData(EOS_TEMPELE) = sim_eleTemp
  eosData(EOS_TEMPRAD) = sim_radTemp
  eosData(EOS_DENS) = sim_rhoInit

#ifdef MODE_DENS_TEMP_ALL
  CALL Eos(MODE_DENS_TEMP_GATHER, vecLen, eosData, mask=eosMask)
#else
  CALL Eos(MODE_DENS_TEMP, vecLen, eosData, mask=eosMask)
#endif
  cv = eosData(EOS_CV)
  abar = eosData(EOS_ABAR)
  zbar = eosData(EOS_ZBAR)
  print*,"cv for all by Eos call is",cv
  sim_CvIon = cv / (zbar + 1)
  sim_CvEle = zbar * sim_CvIon
#else
  call RuntimeParameters_get("eos_singleSpeciesA", abar)
  call RuntimeParameters_get("eos_singleSpeciesZ", zbar)
  cv = 1.5 * (zbar+1) / abar * gasConstant
  sim_CvIon = 1.5 / abar * gasConstant
  sim_CvEle = zbar * sim_CvIon
#endif
  if (sim_meshMe == MASTER_PE) &
  print*,"Cv for ion,ele,all are",sim_CvIon,sim_CvEle,cv,sim_CvIon+sim_CvEle
  deallocate(eosData)

  call PhysicalConstants_get("Avogadro", sim_Avogadro)
  call PhysicalConstants_get("Boltzmann", sim_kBoltzmann)
  call PhysicalConstants_get("electron charge",sim_eleCharge)
  sim_dynamicZ = zbar
  sim_relA = abar
  Ye = sim_dynamicZ / sim_relA
  call PhysicalConstants_get("electron mass",sim_eMassInUAmu,unitMass="amu")
  ionMassInUAmu = sim_relA - sim_dynamicZ * sim_eMassInUAmu
  sim_memi = sim_eMassInUAmu / ionMassInUAmu
  
  call RuntimeParameters_get("sim_schemeOrder", sim_schemeOrder )
  call RuntimeParameters_get("sim_maxTolCoeff0", sim_maxTolCoeff0 )
  call RuntimeParameters_get("sim_maxTolCoeff1", sim_maxTolCoeff1 )
  call RuntimeParameters_get("sim_maxTolCoeff2", sim_maxTolCoeff2 )
  call RuntimeParameters_get("sim_maxTolCoeff3", sim_maxTolCoeff3 )


  if (sim_meshMe == MASTER_PE) then
55   format(1x,'Coulomb Logarithm by NRL Plasma Formulary estimate:',1PG23.16,/ &
     2x,'(p 34, case',I2,')')
     KelvinPerEV = 1/kev
     print*,KelvinPerEV,' Kelvin per eV,',kev,', eV per Kelvin'
     ne = sim_rhoInit * sim_Avogadro * Ye
     if (sim_eleTemp > sim_ionTemp*sim_memi .AND. sim_eleTemp < 10*zbar**2*KelvinPerEV) then
        cl = 23 - 0.5 * log ( ne * zbar**2 / sim_eleTemp**3 )
        print 55, cl, 1
     else if (sim_ionTemp*sim_memi < 10*zbar**2*KelvinPerEV .AND. sim_eleTemp > 10*zbar**2*KelvinPerEV) then
        cl = 24 - log ( sqrt(ne) / sim_eleTemp )
        print 55, cl, 2
     else if (sim_eleTemp < sim_ionTemp*sim_memi) then
        cl = 24 - log ( sqrt(ne) / sim_eleTemp )
        print 55, cl, 3
     else
        print*,' Coulomb Logarithm by NRL Plasma Formulary estimate: No case applies!'
     end if
  end if

end subroutine Simulation_init





