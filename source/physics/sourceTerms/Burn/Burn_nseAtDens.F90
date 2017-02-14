!!****f* source/physics/sourceTerms/Burn/Burn_nseAtDens
!!
!! NAME
!!
!!    Burn_nseAtDens
!!
!!
!! SYNOPSIS
!!
!!    Burn_nseAtDens( real(OUT)   :: qbar_nse,
!!                          real(OUT)   :: sumyi_nse,
!!                          real(OUT)   :: approxtemp,
!!                          real(OUT)   :: edot,
!!                          real(OUT)   :: Yedot,
!!                          real(IN)    :: Ye,
!!                          real(INOUT) :: dens,
!!                          real(INOUT) :: emq)
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
!!    qbar_nse    --  Average binding energy per nucleon (averaged over nse)
!!    sumyi_nse   --  Sum over the abundances. (Y = X/A = X/(#of nucleons))
!!    approxtemp  --  ?? Ask Dean
!!    edot        --  Neutrino and anti-neutrino loss rates [ergs/g]
!!    Yedot       --  Change in electron fraction wrt time (neutronization rate)
!!    Ye          --  Electron fraction
!!    dens        --  Density in/out
!!    emq         --  E-Q = Internal energy minus binding energy (per volume?)
!!
!! NOTES
!!   Takes place after burning is complete.  NSE expands after burning
!!   and the equilibrium state is a function of the density and the
!!   temperature.  Been burned, flame is done, but then the hot stuff
!!   expands.  This changes the composition and releases energy.
!!
!!***


subroutine Burn_nseAtDens(qbar_nse,sumyi_nse,approxtemp,edot,Yedot, Ye, dens, emq)
  implicit none
  
  real, intent(IN) :: Ye, dens, emq
  real, intent(OUT) :: qbar_nse,sumyi_nse,approxtemp,edot,Yedot
  
! NOTE Ivo suggested values that are closer to sensible than just plain zero
!  This is representative of 100% Ni56 at a temperature of 4x10^9K
  qbar_nse = 8.64 ! MeV/nucleon
  sumyi_nse = 1./56.
  approxtemp = 4.0E9 !Kelvin
  edot = 0.0
  Yedot = 0.0

  return
end subroutine Burn_nseAtDens
