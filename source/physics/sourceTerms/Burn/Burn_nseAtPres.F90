!!****f* source/physics/sourceTerms/Burn/Burn_nseAtPres
!!
!! NAME
!!
!!
!!
!! SYNOPSIS
!!
!!    Burn_nseAtPres( real(OUT)   :: qbar_nse,
!!                          real(OUT)   :: sumyi_nse,
!!                          real(OUT)   :: approxtemp,
!!                          real(OUT)   :: edot,
!!                          real(OUT)   :: Yedot,
!!                          real(IN)    :: Ye,
!!                          real(INOUT) :: pres,
!!                          real(INOUT) :: hmq)
!!
!!
!! DESCRIPTION
!!
!! this function gives the NSE final state for a constant pressure burn
!! via table lookup.  qbar is returned in MeV/nucleon
!! only the compositional information (qbar, 1/Abar, edot, Yedot) should be
!! used for hydrodynamics purposes.
!! accurate thermodynamic prperties of the NSE final state (rho, T) should be
!! obtained by solving
!! eint + P/rho - qbar_nse = hmq
!! This is important due to the limited accuracy of the table interpolation
!! (i.e. the interpolated values cannot satisfy any constraints on their own.)
!!  
!! ARGUMENTS
!!  
!!    qbar_nse    --
!!    sumyi_nse   --
!!    approxtemp  --
!!    edot        --
!!    Yedot       --
!!    Ye          --
!!    pres        --
!!    hmq         --  H minus Q = Specific Enthalpy minus binding energy = eint + P/rho
!!
!!***


subroutine Burn_nseAtPres(qbar_nse,sumyi_nse,approxtemp,edot,Yedot, Ye, pres, hmq)
  implicit none
  
  real, intent(IN)    :: Ye, pres, hmq
  real, intent(OUT)   :: qbar_nse,sumyi_nse,approxtemp,edot,Yedot

! NOTE Ivo suggested values that are closer to sensible than just plain zero
!  This is representative of 100% Ni56 at a temperature of 4x10^9K
  qbar_nse = 8.64 ! MeV/nucleon
  sumyi_nse = 1./56.
  approxtemp = 4.0E9 !Kelvin
  edot = 0.0
  Yedot = 0.0
  
  return

end subroutine Burn_nseAtPres
