!!****ih* source/Grid/localAPI/gr_bhInterface
!!
!! This is the header file for the Tree module that defines its
!! public interfaces.
!!***
Module gr_bhInterface

!! Basic components of the tree code: initialization, building of the tree,
!! communication and tree walk

  interface gr_bhBuildTree
    subroutine gr_bhBuildTree()
    implicit none
    end subroutine gr_bhBuildTree
  end interface

  interface gr_bhComBlkProperties
    subroutine gr_bhComBlkProperties()
    implicit none
    end subroutine gr_bhComBlkProperties
  end interface

  interface gr_bhDestroyTree
    subroutine gr_bhDestroyTree()
    end subroutine gr_bhDestroyTree
  end interface

  interface gr_bhExchangeTrees
    subroutine gr_bhExchangeTrees()
    end subroutine gr_bhExchangeTrees
  end interface

  interface gr_bhFinalize
    subroutine gr_bhFinalize()
    implicit none
    end subroutine gr_bhFinalize
  end interface

  interface gr_bhInit
    subroutine gr_bhInit()
    implicit none
    end subroutine gr_bhInit
  end interface

  interface gr_bhInitFieldVar
    subroutine gr_bhInitFieldVar(gpotVar)
    implicit none
    integer, intent(IN) :: gpotVar
    end subroutine gr_bhInitFieldVar
  end interface

  interface gr_bhTreeWalk
    subroutine gr_bhTreeWalk(iterate)
    implicit none
    logical,intent(OUT) :: iterate
    end subroutine gr_bhTreeWalk
  end interface

  interface gr_bhFinalizeIter
    subroutine gr_bhFinalizeIter()
    implicit none
    end subroutine gr_bhFinalizeIter
  end interface

  interface gr_bhGetTreeNodeSize
    real function gr_bhGetTreeNodeSize(level)
      integer, intent(IN) :: level
    end function gr_bhGetTreeNodeSize
  end interface

end Module gr_bhInterface


