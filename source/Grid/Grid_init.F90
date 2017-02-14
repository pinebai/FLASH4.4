!!****f* source/Grid/Grid_init
!!
!! NAME
!!  Grid_init
!!
!! SYNOPSIS
!!
!!  Grid_init()
!!
!! DESCRIPTION
!!  Initialize the runtime parameters needed by 
!!  the grid unit.
!!
!! ARGUMENTS
!!
!! PARAMETERS
!!  
!!   These are the runtime parameters used by the basic Grid unit.
!!   Your specific implementation (Paramesh, UG) WILL have more 
!!   runtime parameters particular to that unit. Since parameters
!!   inheritance is cumulative, all of the ones described here will
!!   be included by all the specific implementations.
!!
!!   To see the default parameter values and all the runtime parameters
!!   specific to your simulation check the "setup_params" file in your
!!   object directory.
!!   You might have overwritten these values with the flash.par values
!!   for your specific run.  
!!
!!    bndPriorityOne [INTEGER]
!!        indicates which axis, x,y,z has priority in guardcell filling at
!!        the physical boundaries.
!!        Defines priority of x direction boundary filling, can be 1,2 or 3
!!    bndPriorityThree [INTEGER]
!!        defines priority of z direction boundary filling, can be 1,2 or 3
!!    bndPriorityTwo [INTEGER]
!!        defines priority of y direction boundary filling, can be 1,2 or 3
!!    convertToConsvdForMeshCalls [BOOLEAN]
!!        indicates if vars are converted from primitive (i.e., velocity) to conservative 
!!        (i.e., momentum) before Grid functions that interpolate data are invoked.
!!        With PARAMESH3 or later, convertToConsvdInMeshInterp should be used instead.
!!    geometry [STRING]
!!        Grid geometry, one of "cartesian", "spherical", "cylindrical", "polar"
!!    unbiased_geometry [BOOLEAN]
!!        attempt to remove floating point bias from geometry discretization
!!        NOT YET IMPLEMENTED
!!    xl_boundary_type [STRING]
!!        lower (left) boundary condition in x dir
!!    xmax [REAL]
!!        physical domain upper bound in x dir
!!    xmin [REAL]
!!        physical domain lower bound in x dir
!!    xr_boundary_type [STRING]
!!        upper (right) boundary condition in x dir
!!    yl_boundary_type [STRING]
!!        lower boundary condition in y dir
!!    ymax [REAL]
!!        physical domain upper bound in y dir
!!    ymin [REAL]
!!        physical domain lower bound in y dir
!!    yr_boundary_type [STRING]
!!        upper boundary condition in y dir
!!    zl_boundary_type [STRING]
!!        lower boundary condition in z dir
!!    zmax [REAL]
!!        physical domain lower bound in x dir
!!    zmin [REAL]
!!        physical domain lower bound in z dir
!!    zr_boundary_type [STRING]
!!        upper boundary condition in z dir
!!    eosMode [STRING]
!!        the default Eos mode, usually MODE_DENS_EI, 
!!        where density and energy are provided to 
!!        calculate pressure and temperature
!!    smallx [REAL]  
!!        cutoff value for abundances, used in Paramesh2 "monotonic" interpolation
!!        for mesh prolongation
!!    interpol_order [INTEGER]  
!!        the order of interpolation, used in Paramesh2 "monotonic" interpolation
!!        for mesh prolongation
!!
!!    useParticles [BOOLEAN]
!!        Whether to initialize and advance particles
!!    pt_maxPerProc [INTEGER]
!!        Maximum number of particles per processor
!!
!!***

subroutine Grid_init()

  implicit none

end subroutine Grid_init

