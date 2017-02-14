!!****f* source/physics/materialProperties/NSE/NSE_finalAtDens
!!
!! NAME
!!
!!    NSE_finalAtDens
!!
!!
!! SYNOPSIS
!!
!!    call NSE_finalAtDens( real(OUT)   :: qbar_nse,
!!                          real(OUT)   :: sumyi_nse,
!!                          real(OUT)   :: approxtemp,
!!                          real(OUT)   :: edot,
!!                          real(OUT)   :: Yedot,
!!                          real(IN)    :: Ye,
!!                          real(IN)    :: dens,
!!                          real(IN)    :: emq)
!!
!!
!! DESCRIPTION
!!  
!! This function gives the NSE final state for a constant density burn
!! via table lookup.  qbar is returned in MeV/nucleon.
!! Only the compositional information (qbar, 1/Abar, edot, Yedot) should be
!! used for hydrodynamics purposes.
!! Accurate thermodynamic prperties of the NSE final state (rho, T) should be
!! obtained by solving
!!        eint - qbar_nse = emq .
!! This is important due to the limited accuracy of the table interpolation
!! (i.e., the interpolated values cannot satisfy any constraints on their own.)
!!
!! ARGUMENTS
!!
!!  outputs  (properties of the NSE state that we want)
!!    qbar_nse    --  Average binding energy per nucleon (averaged over nse)
!!    sumyi_nse   --  Sum over the abundances. (Y = X/A = X/(#of nucleons))
!!    approxtemp  --  approximate temperature of NSE final state (see note above)
!!    edot        --  Neutrino and anti-neutrino loss rates [ergs/g]
!!    Yedot       --  Change in electron fraction wrt time (neutronization rate)
!!
!!  inputs  (these uniquely define the NSE final state)
!!    Ye          --  Electron fraction
!!    dens        --  Density in/out
!!    emq         --  "e-qbar" = Internal energy per gram minus nuclear binding energy
!!                    per nucleon.  The latter is converted to units of ergs/g by assuming
!!                    that each nucleon has mass 1 a.m.u
!!                    This is the total energy in initial and final state -- this quantity
!!                    is constant during an isochoric burn because there is no work done
!!
!! NOTES
!!
!!
!!   Takes place after burning is complete.  NSE expands after burning
!!   and the equilibrium state is a function of the density and the
!!   temperature.  Been burned, flame is done, but then the hot stuff
!!   expands.  This changes the composition and releases energy.
!!
!! SEE ALSO
!!
!!   NSE_interface : top-level interface file for description of subroutine function
!!
!! HISTORY 
!!   Dean Townsley, Alan Calder 2006-8
!!   original interpolation kernel by Alan Calder (2006)
!!   thanks to Flash Code Group (esp. Lynn Reid) for bugfixes
!!
!!   Dean Townsley 2008  stub for when unit is not included
!!***

subroutine NSE_finalAtDens(qbar_nse,sumyi_nse,approxtemp,edot,Yedot, Ye, dens, emq)
  implicit none
  real, intent(IN) :: Ye, dens, emq
  real, intent(OUT) :: qbar_nse,sumyi_nse,approxtemp,edot,Yedot

  qbar_nse = 0.0
  sumyi_nse = 1.0
  approxtemp = 1.0
  edot = 0.0
  Yedot = 0.0
end subroutine NSE_finalAtDens
