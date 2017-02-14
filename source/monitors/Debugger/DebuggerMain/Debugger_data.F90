!!****if* source/monitors/Debugger/DebuggerMain/Debugger_data
!!
!! NAME
!!  Debugger_data
!!
!! SYNOPSIS
!!
!!  use Debugger_data
!!
!! DESCRIPTION
!!
!!  Holds the data needed by the Debugger unit
!!
!!***

module Debugger_data
  implicit none
  logical, save :: dbg_doHeapCheck
  integer, save :: dbg_globalMe
end module Debugger_data
