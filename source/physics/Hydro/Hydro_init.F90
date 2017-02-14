!!****f* source/physics/Hydro/Hydro_init
!!
!! NAME
!!
!!  Hydro_init
!!
!!
!! SYNOPSIS
!!
!!  call Hydro_init()
!!  
!!
!! DESCRIPTION
!! 
!!  Initialize unit scope variables which are typically the runtime parameters.
!!  This must be called once by Driver_initFlash.F90 first. Calling multiple
!!  times will not cause any harm but is unnecessary.
!!
!! ARGUMENTS
!!
!!
!! PARAMETERS
!!
!!   These are the runtime parameters used in the split PPM Hydro 
!!   implementation.
!!
!!   To see the default parameter values and all the runtime parameters
!!   specific to your simulation check the "setup_params" file in your
!!   object directory.
!!   You might have overwritten these values with the flash.par values
!!   for your specific run.  
!!
!!    useHydro [BOOLEAN]
!!        Should any Hydro calculations be performed?
!!    updateHydroFluxes [BOOLEAN]
!!        whether fluxes computed by Hydro should be used to update the solution
!!    cfl [REAL]
!!    geometry [STRING]
!!        Grid geometry
!!    eosMode [STRING]
!!        the default Eos mode, usually MODE_DENS_EI, 
!!        where density and energy are provided to 
!!        calculate pressure and temperature
!!    irenorm
!!    flux_correct [BOOLEAN]
!!    hybrid_riemann [BOOLEAN]
!!    smlrho [REAL]
!!        Cutoff value for density
!!    smallp [REAL]
!!        Cutoff value for pressure
!!    eintSwitch [REAL]  Defined in Eos Unit
!!    cvisc [REAL]
!!        Artificial viscosity constant
!!    dp_sh_md
!!    epsiln [REAL]
!!    nriem [INTEGER]
!!    omg1
!!        PPM dissipation parameter omega1
!!    omg2
!!        PPM dissipation parameter omega2
!!    rieman_tol [REAL]
!!    small [REAL]
!!       Generic small value that can be used as floor where needed
!!    smallu [REAL]
!!       Cutoff value for velocity
!!    smallx [REAL]
!!       Cutoff value for abundances
!!    vgrid
!!    ppm_modifystates
!!    leveque
!!    igodu [INTEGER]
!!       Enable Guodunov method instead of PPM if set to 1
!!    iplm [INTEGER]
!!       Enable piecewise linear method instead of PPM if set to 1
!!    use_steepening [BOOLEAN]
!!    use_cma_flattening [BOOLEAN]
!!    use_cma_advection [BOOLEAN]
!!
!!***

subroutine Hydro_init()

  implicit none

end subroutine Hydro_init

