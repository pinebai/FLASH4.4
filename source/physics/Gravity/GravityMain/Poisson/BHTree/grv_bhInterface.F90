!!****ih* source/physics/Gravity/GravityMain/Poisson/BHTree/grv_bhInterface
!!
!! This is an interface module for the Gravity unit that defines some
!! private interfaces specific to the BHTree implementation.
!!***

Module grv_bhInterface

  interface grv_bhGenEwaldField
    subroutine grv_bhGenEwaldField()
    end subroutine grv_bhGenEwaldField
  end interface

  interface grv_bhGenEwaldFieldLevelV42
    subroutine grv_bhGenEwaldFieldLevelV42(nx, ny, nz, ewald_xmax, ewald_ymax, ewald_zmax, &
           & ewald_periodicity, field_pot, field_acc)
      implicit none
      integer, intent(IN) :: ewald_periodicity, nx, ny, nz
      real, intent(IN) ::  ewald_xmax, ewald_ymax, ewald_zmax
      real , intent(INOUT) :: field_pot(-1:,-1:,-1:)
      real , intent(INOUT) :: field_acc(1:,-1:,-1:,-1:)
    end subroutine grv_bhGenEwaldFieldLevelV42
  end interface

  interface grv_bhGenEwaldFieldLevel
    subroutine grv_bhGenEwaldFieldLevel(nx, ny, nz, field_ewald)
      implicit none
      integer, intent(IN) :: nx, ny, nz
      real, intent(OUT) :: field_ewald(0:,0:,0:,0:)
    end subroutine grv_bhGenEwaldFieldLevel
  end interface

  interface grv_elintF
    real function grv_elintF(theta, k)
      implicit none
      real, intent(in) :: theta, k
    end function grv_elintF
  end interface

  interface grv_IntSimpson
    real function grv_IntSimpson(fce,k,x,y,z,a,n)
      implicit none
      real, external :: fce
      integer, intent(in) :: k,n
      real, intent(in) :: x,y,z,a
    end function grv_IntSimpson
  end interface

  interface grv_coef1P2I
    real function grv_coef1P2I(x,k,y,z)
      implicit none
      integer, intent(in) :: k
      real, intent(in) :: x,y,z
    end function grv_coef1P2I
  end interface

  interface grv_coef1P2I_der
    real function grv_coef1P2I_der(x,k,y,z)
      implicit none
      integer, intent(in) :: k
      real, intent(in) :: x,y,z
    end function grv_coef1P2I_der
  end interface

  interface grv_accExtrnPotRow
    subroutine grv_accExtrnPotRow(pos, sweepDir, blockID, numCells, grav)
      implicit none
      integer, dimension(2), intent(in) :: pos
      integer, intent(in)               :: sweepDir, blockID,  numCells
      real, intent(inout)               :: grav(numCells)
    end subroutine grv_accExtrnPotRow
  end interface

  interface grv_readExtrnPotential
    subroutine grv_readExtrnPotential()
      implicit none
    end subroutine grv_readExtrnPotential
  end interface


  interface grv_derfc
    DOUBLE PRECISION FUNCTION DERFC(X)
      DOUBLE PRECISION X
    END FUNCTION DERFC
  end interface

  interface grv_derfcx
    DOUBLE PRECISION FUNCTION DERFCX(X)
      DOUBLE PRECISION X
    END FUNCTION DERFCX
  end interface

  interface grv_derf
    DOUBLE PRECISION FUNCTION DERF(X)
      DOUBLE PRECISION X
    END FUNCTION DERF
  end interface

  interface grv_besj0
    DOUBLE PRECISION FUNCTION BESJ0(X)
      DOUBLE PRECISION X
    END FUNCTION BESJ0
  end interface

  interface grv_besj1
    DOUBLE PRECISION FUNCTION BESJ1(X)
      DOUBLE PRECISION X
    END FUNCTION BESJ1
  end interface

  interface
    function grv_bhAccShort(ewald_alpha,xni,yni,zni)
      implicit none
      real, intent(IN) :: ewald_alpha,xni,yni,zni
      real :: grv_bhAccShort(0:12)
    end
  end interface

  interface
    function grv_bhAccLong1P2I(Linv, ewald_dzeta, ewald_eta, hi, x, y, z)
      implicit none
      real, intent(IN) :: Linv, ewald_dzeta, ewald_eta, x, y, z
      integer, intent(IN) :: hi
      real :: grv_bhAccLong1P2I(0:12)
    end
  end interface

  interface
    function grv_bhAccLong2P1I(Linv, ewald_dzeta, ratio_pinv, hi, hj, x, y, z)
      implicit none
      real, intent(IN) :: Linv, ewald_dzeta, ratio_pinv, x, y, z
      integer, intent(IN) :: hi, hj
      real :: grv_bhAccLong2P1I(0:12)
    end
  end interface

  interface
    function grv_bhAccLong3P(Linv, ewald_dzeta, ratio_pinv1, ratio_pinv2, hi, hj, hk, x, y, z)
      implicit none
      real, intent(IN) :: Linv, ewald_dzeta, ratio_pinv1, ratio_pinv2, x, y, z
      integer, intent(IN) :: hi, hj, hk
      real :: grv_bhAccLong3P(0:12)
    end
  end interface

end Module grv_bhInterface
