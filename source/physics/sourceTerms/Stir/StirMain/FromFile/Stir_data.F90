!!****if* source/physics/sourceTerms/Stir/StirMain/FromFile/Stir_data
!!
!! NAME
!!  Stir_data
!!
!! SYNOPSIS
!!  Stir_data()
!!
!! DESCRIPTION
!!  Stores the local data for Source Term: StirFromFile
!!
!! AUTHOR
!!  Christoph Federrath, 2008
!!
!!***

Module Stir_data

  integer, parameter :: st_maxmodes = 1000

  integer, save :: st_meshMe, st_meshComm, st_globalMe

  ! number of modes
  integer, save :: st_nmodes

  real,save, dimension(3,st_maxmodes) :: st_mode, st_aka, st_akb
  real,save, dimension(6*st_maxmodes) :: st_OUphases
  real,save, dimension(  st_maxmodes) :: st_ampl

  character (len=80), save :: st_infilename
  real, save               :: st_lastTimeUpdatedAccel, st_dtUpdateAccel

  logical,save  :: st_useStir, st_computeDt
  real,save     :: st_decay
  real,save     :: st_energy
  real,save     :: st_stirmin, st_stirmax
  real,save     :: st_solweight
  real,save     :: st_solweightnorm
  integer,save  :: st_spectform

end Module Stir_data
