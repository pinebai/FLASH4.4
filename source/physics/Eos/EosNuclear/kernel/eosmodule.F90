 module eosmodule

   implicit none
   
   integer,save :: nrho,ntemp,nye

   integer,save :: warn_from !warn from given reflevel

   real,save :: energy_shift = 0.0e0
   
   real,save :: precision = 1.0e-10

!   real, parameter :: e_zeroPoint = 9.2927e18
   real, save :: e_zeroPoint
!   real, parameter :: e_zeroPoint = 9.2927e19

! min-max values:
   real,save :: eos_rhomin,eos_rhomax
   real,save :: eos_yemin,eos_yemax
   real,save :: eos_tempmin,eos_tempmax

   real,save :: eos_table_tmax

   real,save :: t_max_hack = 240.0e0

! basics
   integer, parameter :: nvars = 19
   real,save,allocatable :: alltables(:,:,:,:)
  ! index variable mapping:
  !  1 -> logpress
  !  2 -> logenergy
  !  3 -> entropy
  !  4 -> munu
  !  5 -> cs2
  !  6 -> dedT
  !  7 -> dpdrhoe
  !  8 -> dpderho
  !  9 -> muhat
  ! 10 -> mu_e
  ! 11 -> mu_p
  ! 12 -> mu_n
  ! 13 -> xa
  ! 14 -> xh
  ! 15 -> xn
  ! 16 -> xp
  ! 17 -> abar
  ! 18 -> zbar
  ! 19 -> gamma

   real,allocatable,save :: logrho(:)
   real,allocatable,save :: logtemp(:)
   real,allocatable,save :: ye(:)

! constants
   real,parameter :: mev_to_erg = 1.60217733e-6
   real,parameter :: amu_cgs = 1.66053873e-24
   real,parameter :: amu_mev = 931.49432e0
   real,parameter :: pi = 3.14159265358979e0
   real,parameter :: ggrav = 6.672e-8
   real,parameter :: temp_mev_to_kelvin = 1.16045221e10
   real,parameter :: clight = 2.99792458e10
   real,parameter :: kb_erg = 1.380658e-16
   real,parameter :: kb_mev = 8.61738568e-11   
   real, parameter  :: arad = 7.56557705e-15

  character(len=80),save :: eos_file

  real, save :: avo
  real, parameter :: amu = 1.6605402e-24, &
                     h   = 6.6260755e-27

  real, parameter :: sioncon = (2.0e0 * pi * amu * kb_erg)/(h*h)

 end module eosmodule
