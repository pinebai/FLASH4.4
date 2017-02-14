

module ins_statsInterface

  implicit none

  interface
    subroutine ins_statsVelpTimeAvg(n)
      implicit none
      integer, intent(in) :: n
    end subroutine
  end interface

  interface
    subroutine ins_statsRestressesTimeavg(n)
      implicit none
      integer, intent(in) :: n
    end subroutine
  end interface

  interface
    subroutine IncompNS_statsIOExport(expt_flag)
      implicit none
      logical, intent(in) :: expt_flag
    end subroutine
  end interface

end module
