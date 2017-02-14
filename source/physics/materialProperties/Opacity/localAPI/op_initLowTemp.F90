!!****if* source/physics/materialProperties/Opacity/localAPI/op_initLowTemp
!!
!! NAME
!!
!!  op_initLowTemp
!!
!! SYNOPSIS
!!
!!  call op_initLowTemp ()
!!
!! DESCRIPTION
!!
!!  Inititalizes the section of the low temperature opacity unit. Several arrays are
!!  allocated here and initialized with all the data. The number of atomic elements
!!  and their corresponding atomic numbers for the current opacity run must be known
!!  at this stage, otherwise the program will stop with a message.
!!
!!  The goal of this routine is to get and hold only the data corresponding to the
!!  atomic elements needed in memeory. The big arrays containing the data for all the
!!  elements will only be alive during the time spent in this routine.
!!
!! ARGUMENTS
!!
!!***
subroutine op_initLowTemp ()

  implicit none

  return
end subroutine op_initLowTemp
