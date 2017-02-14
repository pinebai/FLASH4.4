!!****if* source/Grid/GridBoundaryConditions/Flash3HSE/gr_bcHseData
!!
!! NAME
!!  gr_bcHseData
!!
!! SYNOPSIS
!!
!!  use gr_bcHseData
!!
!!
!!***

Module gr_bcHseData
  integer,parameter :: HSE_FORWARD = 1
  integer,parameter :: HSE_BACKWARD = 2
  integer,parameter :: HSE_CONSTENTR = 3
  integer,parameter :: HSE_CONSTTEMP = 4
  integer,parameter :: HSE_SETTEMP = 5

  character(len=80), save :: gr_bcHseGravDirec
  integer, save :: gr_bcHseDirection

  real, save :: gr_bcHseGravConst

end Module gr_bcHseData
