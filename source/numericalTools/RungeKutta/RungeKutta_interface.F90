!!****h* source/numericalTools/RungeKutta/RungeKutta_interface
!!
!! NAME
!!
!!  RungeKutta_interface
!!
!! SYNOPSIS
!!
!!  use RungeKutta_interface
!!
!! DESCRIPTION
!!
!!  This is the header file for the Runge Kutta unit that defines its public interfaces.
!!
!!***

Module RungeKutta_interface

  interface
    subroutine RungeKutta_finalize ()
    end subroutine RungeKutta_finalize
  end interface

  interface
    subroutine RungeKutta_init ()
    end subroutine RungeKutta_init
  end interface

  interface
    subroutine RungeKutta_step (method,                      &
                                f, x, y,                     &
                                eFrac, eBase,                &
                                htry,                        &
                                               hused, hnext, &
                                               yout, eout    )
      interface
        function f (x,y)
          real, intent (in) :: x
          real, intent (in) :: y (:)
          real              :: f (1:size (y))
        end function f
      end interface

      character (len=*), intent (in)  :: method
      real,              intent (in)  :: x
      real,              intent (in)  :: y     (:)
      real,              intent (in)  :: eFrac
      real,              intent (in)  :: eBase (:)
      real,              intent (in)  :: htry
      real,              intent (out) :: hused
      real,              intent (out) :: hnext
      real,              intent (out) :: yout  (:)
      real,              intent (out) :: eout  (:)
    end subroutine RungeKutta_step
  end interface

  interface
    subroutine RungeKutta_stepConfined (method,                      &
                                        f, nc,                       &
                                        x, y, ymin, ymax,            &
                                        eFrac, eBase,                &
                                        htry,                        &
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

      character (len=*), intent (in)  :: method
      integer,           intent (in)  :: nc
      real,              intent (in)  :: x
      real,              intent (in)  :: y     (:)
      real,              intent (in)  :: eFrac
      real,              intent (in)  :: eBase (:)
      real,              intent (in)  :: htry
      real,              intent (out) :: hused
      real,              intent (out) :: hnext
      real,              intent (out) :: yout  (:)
      real,              intent (out) :: eout  (:)
    end subroutine RungeKutta_stepConfined
  end interface

  interface
    real function RungeKutta_stepSizeEstimate (method,f,x,y,eFrac,eBase,hmax)

      interface
        function f (x,y)
          real, intent (in) :: x
          real, intent (in) :: y (:)
          real              :: f (1:size (y))
        end function f
      end interface

      character (len=*), intent (in)  :: method
      real,              intent (in)  :: x
      real,              intent (in)  :: y     (:)
      real,              intent (in)  :: eFrac
      real,              intent (in)  :: eBase (:)
      real, optional,    intent (in)  :: hmax
    end function RungeKutta_stepSizeEstimate
  end interface

end Module RungeKutta_interface
