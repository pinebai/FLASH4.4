!!****ih* source/numericalTools/RungeKutta/localAPI/rk_interface
!!
!! NAME
!!
!!  rk_interface
!!
!! SYNOPSIS
!!
!!  use rk_interface
!!
!! DESCRIPTION
!!
!!  This is the header file for the Runge Kutta unit that defines its private interfaces.
!!
!!***

Module rk_interface

  interface
    integer function rk_orderRKmethod (method)
      character (len=*), intent (in) :: method
    end function rk_orderRKmethod
  end interface

  interface
    subroutine rk_setButcherTableauRepository () 
    end subroutine rk_setButcherTableauRepository
  end interface

  interface
    subroutine rk_stepEA (f, n, s,                      &
                          x, y,                         &
                          eMin, ePower, eFrac, eBase,   &
                          htry,                         &
                                          hused, hnext, &
                                          yout, eout    )
      interface
        function f (x,y)
          real, intent (in) :: x
          real, intent (in) :: y (:)
          real              :: f (1:size (y))
        end function f
      end interface

      integer, intent (in)  :: n
      integer, intent (in)  :: s
      real,    intent (in)  :: x
      real,    intent (in)  :: y     (:)
      real,    intent (in)  :: eMin
      real,    intent (in)  :: ePower
      real,    intent (in)  :: eFrac
      real,    intent (in)  :: eBase (:)
      real,    intent (in)  :: htry
      real,    intent (out) :: hused
      real,    intent (out) :: hnext
      real,    intent (out) :: yout  (:)
      real,    intent (out) :: eout  (:)
    end subroutine rk_stepEA
  end interface

  interface
    subroutine rk_stepEAC (f, n, nc, s,                  &
                           x, y, ymin, ymax,             &
                           eMin, ePower, eFrac, eBase,   &
                           htry,                         &
                                           hused, hnext, &
                                           yout, eout    )
      interface
        function f (x,y)
          real, intent (in) :: x
          real, intent (in) :: y (:)
          real              :: f (1:size (y))
        end function f
      end interface

      interface
        function ymin (nc,y)
          integer, intent (in) :: nc
          real,    intent (in) :: y (:)
          real                 :: ymin (1:nc)
        end function ymin
      end interface

      interface
        function ymax (nc,y)
          integer, intent (in) :: nc
          real,    intent (in) :: y (:)
          real                 :: ymax (1:nc)
        end function ymax
      end interface

      integer, intent (in)  :: n
      integer, intent (in)  :: nc
      integer, intent (in)  :: s
      real,    intent (in)  :: x
      real,    intent (in)  :: y     (:)
      real,    intent (in)  :: eMin
      real,    intent (in)  :: ePower
      real,    intent (in)  :: eFrac
      real,    intent (in)  :: eBase (:)
      real,    intent (in)  :: htry
      real,    intent (out) :: hused
      real,    intent (out) :: hnext
      real,    intent (out) :: yout  (:)
      real,    intent (out) :: eout  (:)
    end subroutine rk_stepEAC
  end interface

end Module rk_interface
