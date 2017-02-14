!!****if* source/Simulation/SimulationMain/Chemistry_Test/Simulation_init
!!
!! NAME
!!  
!!  Simulation_init
!!
!! SYNOPSIS
!!
!!  Simulation_init()
!!
!!
!! DESCRIPTION
!!
!!  Initializes all the parameters needed for the Marcus' cluster problem
!!  
!!
!! ARGUMENTS
!!
!!
!! PARAMETERS
!!
!!
!!***

subroutine Simulation_init()


!!***used modules from FLASH.***
  
  use Simulation_data
  use PhysicalConstants_interface, ONLY: PhysicalConstants_get
  use RuntimeParameters_interface, ONLY: RuntimeParameters_get,RuntimeParameters_getPrev
  use Logfile_interface, ONLY : Logfile_stamp

  implicit none
#include "Flash.h"
#include "constants.h"




  real :: nblockx,nblocky,nblockz
  real :: sim_lrefine_max
 
!! These variables check to see if the Geometry is right 
  integer :: meshGeom
  logical :: validGeom


 

     print *,'Getting Physical Constants'
     call PhysicalConstants_get("ideal gas constant",sim_gasConst)
     
     call RuntimeParameters_get("xmin", sim_xMin)
     call RuntimeParameters_get("ymin", sim_yMin)
 
     call RuntimeParameters_get("xmax", sim_xMax)
     call RuntimeParameters_get("ymax", sim_yMax)

     print *,'Getting RuntimeParameters'

     call RuntimeParameters_get(nblockx,"sim_nblockx")
     call RuntimeParameters_get(nblocky,"sim_nblocky")

     if(NDIM == 3) then
        call RuntimeParameters_get("zmin", sim_zMin)
        call RuntimeParameters_get("zmax", sim_zMax)
        call RuntimeParameters_get(nblockz, "sim_nblockz" )
     endif


     call RuntimeParameters_get("sim_c_temp", sim_c_temp)
     call RuntimeParameters_get("sim_c_den", sim_c_den)

     call RuntimeParameters_get("sim_xH", sim_xH)
     call RuntimeParameters_get("sim_xHP", sim_xHP)
     call RuntimeParameters_get("sim_xHM", sim_xHM)
     call RuntimeParameters_get("sim_xD", sim_xD)
     call RuntimeParameters_get("sim_xDM", sim_xDM)
     call RuntimeParameters_get("sim_xDP", sim_xDP)
     call RuntimeParameters_get("sim_xHE", sim_xHE)
     call RuntimeParameters_get("sim_xHEP", sim_xHEP)
     call RuntimeParameters_get("sim_xHEPP", sim_xHEPP)
     call RuntimeParameters_get("sim_xH2", sim_xH2)
     call RuntimeParameters_get("sim_xH2P", sim_xH2P)
     call RuntimeParameters_get("sim_xHD", sim_xHD)
     call RuntimeParameters_get("sim_xHDP", sim_xHDP)
     call RuntimeParameters_get("sim_xD2", sim_xD2)
     call RuntimeParameters_get("sim_xD2P", sim_xD2P)
     call RuntimeParameters_get("sim_xELEC", sim_xELEC)     

     call RuntimeParameters_get("sim_pchem_time", sim_pchem_time)
     call RuntimeParameters_get("sim_cool_time", sim_cool_time)
     call RuntimeParameters_get("sim_meta", sim_meta)

    return

end subroutine Simulation_init






