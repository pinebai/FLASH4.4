!!****f* source/physics/materialProperties/NSE/NSE_finalAtPres
!!
!! NAME
!!
!!   NSE_finalAtPres
!!
!! SYNOPSIS
!!
!!   call  NSE_finalAtPres( real(OUT)   :: qbar_nse,
!!                          real(OUT)   :: sumyi_nse,
!!                          real(OUT)   :: approxtemp,
!!                          real(OUT)   :: edot,
!!                          real(OUT)   :: Yedot,
!!                          real(IN)    :: Ye,
!!                          real(IN)    :: pres,
!!                          real(IN)    :: hmq)
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
!! NOTES
!!   see top level interfaces file for description of subroutine function
!! SEE ALSO
!!   NSE_interface
!! HISTORY
!!   Dean Townsley, Alan Calder 2006-8
!!   original interpolation kernel by Alan Calder 2006
!!
!!    Dean Townsley 2008  stub for when unit is not included
!!***

subroutine NSE_finalAtPres(qbar_nse,sumyi_nse,approxtemp,edot,Yedot, Ye, pres, hmq)
  implicit none
  real, intent(IN) :: Ye, pres, hmq
  real, intent(OUT) :: qbar_nse,sumyi_nse,approxtemp,edot,Yedot

  qbar_nse = 0.0
  sumyi_nse = 1.0
  approxtemp = 1.0
  edot = 0.0
  Yedot = 0.0
end subroutine
