!****h* source/physics
!!
!! NAME
!!
!!
!! SYNOPSIS
!!
!!
!! DESCRIPTION
!!
!! This is the header file for the 
!! module that defines its public interfaces.
!!
!!***


Module SolidMechanics_interface

  implicit none

  interface !SolidMechanics_init
     subroutine SolidMechanics_init(restart)
       implicit none
       logical, INTENT(IN) :: restart
     end subroutine SolidMechanics_init
  end interface

  interface !SolidMechanics 
     subroutine SolidMechanics(selector_flag,restart,convflag_all)
       implicit none
       integer, intent(in)              :: selector_flag
       logical, intent(in),optional     :: restart
       integer, optional, intent(inout) :: convflag_all
     end subroutine SolidMechanics
  end interface 

 interface  !SolidMechanics_finalize
    subroutine SolidMechanics_finalize()
      implicit none
    end subroutine SolidMechanics_finalize
 end interface 

 interface
    subroutine SolidMechanics_computeDt(dt_solid)
      implicit none
      real, INTENT(INOUT) :: dt_solid
    end subroutine SolidMechanics_computeDt
 end interface

 
end Module SolidMechanics_interface
