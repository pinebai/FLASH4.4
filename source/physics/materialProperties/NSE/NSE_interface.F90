! Interface definitions for NSE properties unit
!
! Dean Townsley 2008
!

module NSE_interface

  interface
     subroutine NSE_init()
     end subroutine
  end interface

  interface
     subroutine NSE_finalize()
     end subroutine
  end interface

  interface
     subroutine NSE_finalAtDens(qbar_nse,sumyi_nse,approxtemp,edot,Yedot, Ye, dens, emq)
        implicit none
        real, intent(IN) :: Ye, dens, emq
        real, intent(OUT) :: qbar_nse,sumyi_nse,approxtemp,edot,Yedot
        ! This function gives the NSE final state for a constant density burn
        ! via table lookup.  qbar is returned in MeV/nucleon.
        ! Only the compositional information (qbar, 1/Abar, edot, Yedot) should be
        ! used for hydrodynamics purposes.  the approximate temperature should
        ! be used only as an approximation.
        ! Accurate thermodynamic prperties of the NSE final state (T) should be
        ! obtained by solving
        !        eint - qbar_nse = emq .
        ! This is important due to the limited accuracy of the table interpolation
        ! (i.e., the interpolated values cannot satisfy any constraints on their own.)
        !
        ! argument description
        !
        !  inputs:  (these uniquely define the NSE final state)
        !    Ye          --  Electron fraction (# of electrons = # of protons   per nucleon)
        !    dens        --  Density of unburned fuel *and* requested final state
        !    emq         --  "e-qbar" = Internal energy per gram minus nuclear binding energy
        !                    per nucleon.  The latter is converted to units of ergs/g by assuming
        !                    that each nucleon has mass 1 a.m.u
        !                    This is the total energy in initial and final state -- this quantity
        !                    is constant during an isochoric burn because there is no work done
        !
        !  outputs:  (properties of the NSE state that we want)
        !    qbar_nse    --  Average binding energy per nucleon in final nse state (MeV/nucleon)
        !    sumyi_nse   --  Ion number (# of ion particles in fluid per nucleon) in NSE
        !    approxtemp  --  approximate temperature of NSE final state (see note above)
        !    edot        --  Neutrino and anti-neutrino loss rates [ergs/g]
        !    Yedot       --  Change in electron fraction wrt time (neutronization rate)
     end subroutine NSE_finalAtDens
  end interface

  interface
     subroutine NSE_finalAtPres(qbar_nse,sumyi_nse,approxtemp,edot,Yedot, Ye, pres, hmq)
       implicit none
       real, intent(IN)    :: Ye, pres, hmq
       real, intent(OUT)   :: qbar_nse,sumyi_nse,approxtemp,edot,Yedot
       ! this function gives the NSE final state for a constant pressure burn
       ! via table lookup.  qbar is returned in MeV/nucleon
       ! only the compositional information (qbar, 1/Abar, edot, Yedot) should be
       ! used for hydrodynamics purposes.
       ! accurate thermodynamic prperties of the NSE final state (rho, T) should be
       ! obtained by solving
       ! eint + P/rho - qbar_nse = hmq
       ! This is important due to the limited accuracy of the table interpolation
       ! (i.e. the interpolated values cannot satisfy any constraints on their own.)
       !
       ! argument description
       !   as for NSE_finalAtPres except
       !     pres  -- Pressure in erg/cc
       !     hmq   --  h - qbar = Specific Enthalpy per gram ( = eint+P/rho )
       !                        minus nuclear binding energy per nucleon.
       !                        nuclear binding energy is converted to erg/g bu assuming
       !                        each nucleon has mass 1 a.m.u.
     end subroutine NSE_finalAtPres
  end interface

end module NSE_interface
