!!****if* source/physics/sourceTerms/Burn/BurnMain/nuclearBurn/Burn_init
!!
!! NAME
!!  
!!  Burn_init
!!
!!
!! SYNOPSIS
!! 
!!  call Burn_init()
!!
!!  
!! DESCRIPTION
!!
!!  Initializes various runtime parameters and the specific alpha-chain nuclear burning 
!!  network for the Burn unit.
!!
!! ARGUMENTS
!!
!!  
!!
!! PARAMETERS
!!
!!  useBurn -- Boolean, True.  Turns on burning module
!!  useBurnTable -- Boolean, False.  Controls the generation of reaction rates.
!!                TRUE interpolates from a stored table; FALSE generates them
!!                analytically.
!!  useShockBurn -- Boolean, FALSE.  Controls whether burning is allowed inside
!!                a regime experiencing shocks
!!  algebra -- Integer, 1, [1,2].  Controls choice of linear algebra package used
!!                for matrix solution.  1=Ma28 sparse package, 2=Gift hardwired package.
!!  odeStepper -- Integer, 1, [1,2].  Controls time integration routines.
!!                1=Bader-Deuflhard variable order, 2=Rosenbrock 4th order
!!  nuclearTempMin/Max -- Real, 1.1E+8/1.0E+12.  Minimum and maximum temperature
!!                ranges where burning can occur
!!  nuclearDensMin/Max -- Real, 1.0E-10/1.0E+14.  Minimum and maximum density range
!!                where burning can occur.
!!  nuclearNI56Max -- Real, 1.0.  Maximum mass fraction of nickel where burning
!!                can occur.
!!  enucDtFactor -- Real, 1.0E+30.  Timestep limiter.  See Burn_computeDt for details.            
!!
!!  SEE ALSO
!!    bn_initNetwork
!!
!!***


subroutine Burn_init()

  use Burn_data, ONLY: bn_algebra, bn_odeStepper, bn_useBurnTable, bn_useBurn, bn_meshMe, &
     &    bn_useShockBurn, bn_smallx, &
     &    bn_nuclearTempMin, bn_nuclearTempMax, bn_nuclearDensMin, bn_nuclearDensMax, &
     &    bn_nuclearNI56Max, bn_enucDtFactor, bn_meshMe
  use Driver_interface, ONLY : Driver_abortFlash, Driver_getMype
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  implicit none

#include "constants.h"
#include "Flash.h"


! This strange initialization is to keep interface consistent with other sourceTerms
  call Driver_getMype(MESH_COMM, bn_meshMe)

  call RuntimeParameters_get("useBurn", bn_useBurn )
  if (.not. bn_useBurn) then
     write(6,*)'WARNING:  You have included the Burn unit but have set '
     write(6,*)'   the runtime parameter useBurn to FALSE'
     write(6,*)'   No burning will occur but Burn_init will continue.'
  end if

  call RuntimeParameters_get("algebra", bn_algebra )
  call RuntimeParameters_get("odeStepper", bn_odeStepper )
  call RuntimeParameters_get("useBurnTable", bn_useBurnTable )

  if ((bn_algebra .lt. 1) .or. (bn_algebra .gt. 2)) then
     write(6,*) 
     write(6,*) 'only algebra=1 = ma28'
     write(6,*) 'and  algebra=2 = gift are valid'
     write(6,*) 'and you have specified algebra=',bn_algebra
     write(6,*) 'error in routine Burn'
     call Driver_abortFlash('ERROR in Burn, wrong algebra')
  end if
  if ((bn_odeStepper .lt. 1) .or. (bn_odeStepper .gt. 2)) then
     write(6,*) 
     write(6,*) 'only odeStepper=1 = bader-deuflhard'
     write(6,*) 'and  odeStepper=2 = rosenbrock integration are valid'    
     write(6,*) 'and you have specified odeStepper=',bn_odeStepper
     write(6,*) 'error in routine Burn'
     call Driver_abortFlash('ERROR in Burn, wrong integration type')
  end if

  call RuntimeParameters_get ('useShockBurn', bn_useShockBurn)


  !    eos_order = (/DENS_VAR, TEMP_VAR, PRES_VAR, ENER_VAR, EINT_VAR,         &
  !         VELX_VAR, VELY_VAR, VELZ_VAR, GAME_VAR, GAMC_VAR,         &
  !         SPECIES_BEGIN/)

  call RuntimeParameters_get( 'smallx', bn_smallx)
  call RuntimeParameters_get( 'nuclearTempMin', bn_nuclearTempMin)
  call RuntimeParameters_get( 'nuclearTempMax', bn_nuclearTempMax)
  call RuntimeParameters_get( 'nuclearDensMin', bn_nuclearDensMin)
  call RuntimeParameters_get( 'nuclearDensMax', bn_nuclearDensMax)
  call RuntimeParameters_get( 'nuclearNI56Max', bn_nuclearNI56Max)

  call RuntimeParameters_get('enucDtFactor', bn_enucDtFactor)


  !!  Now initialize the network things
  call bn_initNetwork()

end subroutine Burn_init
