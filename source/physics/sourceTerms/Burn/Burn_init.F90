!!****f* source/physics/sourceTerms/Burn/Burn_init
!!
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


  implicit none

  

  return

end subroutine Burn_init
