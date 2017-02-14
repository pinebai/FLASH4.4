!!****if* source/physics/RadTrans/localAPI/rt_interface
!!
!! NAME
!!   rt_interface
!!
!! SYNOPSIS
!!   use rt_interface : ONLY
!!
!!  DESCRIPTION 
!!    Interface for internal RadTrans subroutines
!!
!!***
module rt_interface
  interface
     subroutine rt_init
       implicit none
     end subroutine rt_init
  end interface
end module rt_interface
