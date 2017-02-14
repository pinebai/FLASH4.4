!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/op_setAtomWeights
!!
!! NAME
!!
!!  op_setAtomWeights
!!
!! SYNOPSIS
!!
!!  call op_setAtomWeights ()
!!
!! DESCRIPTION
!!
!!  This routine sets the atomic weights in g/mole.
!!
!! ARGUMENTS
!!
!!***
subroutine op_setAtomWeights ()

  use Opacity_data,      ONLY : op_atomWeight
  use Driver_interface,  ONLY : Driver_abortFlash

  implicit none

  logical :: goodShape
!
!
!   ...Check dimensions of the atomic element weight array.
!
!
  goodShape = size (op_atomWeight) == 100

  if (.not.goodShape) then
       call Driver_abortFlash ('[op_setAtomWeights] ERROR: array op_atomWeight has bad shape')
  end if
!
!
!   ...Set the atomic weights.
!
!
  op_atomWeight (  1) = 1.008
  op_atomWeight (  2) = 4.003
  op_atomWeight (  3) = 6.941
  op_atomWeight (  4) = 9.012
  op_atomWeight (  5) = 10.811
  op_atomWeight (  6) = 12.011
  op_atomWeight (  7) = 14.007
  op_atomWeight (  8) = 15.999
  op_atomWeight (  9) = 18.998
  op_atomWeight ( 10) = 20.18

  op_atomWeight ( 11) = 22.99
  op_atomWeight ( 12) = 24.305
  op_atomWeight ( 13) = 26.982
  op_atomWeight ( 14) = 28.086
  op_atomWeight ( 15) = 30.974
  op_atomWeight ( 16) = 32.065
  op_atomWeight ( 17) = 35.453
  op_atomWeight ( 18) = 39.948
  op_atomWeight ( 19) = 39.098
  op_atomWeight ( 20) = 40.078

  op_atomWeight ( 21) = 44.956
  op_atomWeight ( 22) = 47.867
  op_atomWeight ( 23) = 50.942
  op_atomWeight ( 24) = 51.996
  op_atomWeight ( 25) = 54.938
  op_atomWeight ( 26) = 55.845
  op_atomWeight ( 27) = 58.933
  op_atomWeight ( 28) = 58.693
  op_atomWeight ( 29) = 63.546
  op_atomWeight ( 30) = 65.409

  op_atomWeight ( 31) = 69.723
  op_atomWeight ( 32) = 72.64
  op_atomWeight ( 33) = 74.922
  op_atomWeight ( 34) = 78.96
  op_atomWeight ( 35) = 79.904
  op_atomWeight ( 36) = 83.8
  op_atomWeight ( 37) = 85.468
  op_atomWeight ( 38) = 87.62
  op_atomWeight ( 39) = 88.906
  op_atomWeight ( 40) = 91.224

  op_atomWeight ( 41) = 92.906
  op_atomWeight ( 42) = 95.94
  op_atomWeight ( 43) = 98
  op_atomWeight ( 44) = 101.07
  op_atomWeight ( 45) = 102.906
  op_atomWeight ( 46) = 106.42
  op_atomWeight ( 47) = 107.868
  op_atomWeight ( 48) = 112.411
  op_atomWeight ( 49) = 114.818
  op_atomWeight ( 50) = 118.71

  op_atomWeight ( 51) = 121.76
  op_atomWeight ( 52) = 127.6
  op_atomWeight ( 53) = 126.904
  op_atomWeight ( 54) = 131.293
  op_atomWeight ( 55) = 132.905
  op_atomWeight ( 56) = 137.327
  op_atomWeight ( 57) = 138.905
  op_atomWeight ( 58) = 140.116
  op_atomWeight ( 59) = 140.908
  op_atomWeight ( 60) = 144.242

  op_atomWeight ( 61) = 145
  op_atomWeight ( 62) = 150.36
  op_atomWeight ( 63) = 151.964
  op_atomWeight ( 64) = 157.25
  op_atomWeight ( 65) = 158.925
  op_atomWeight ( 66) = 162.5
  op_atomWeight ( 67) = 164.93
  op_atomWeight ( 68) = 167.259
  op_atomWeight ( 69) = 168.934
  op_atomWeight ( 70) = 173.04

  op_atomWeight ( 71) = 174.967
  op_atomWeight ( 72) = 178.49
  op_atomWeight ( 73) = 180.948
  op_atomWeight ( 74) = 183.84
  op_atomWeight ( 75) = 186.207
  op_atomWeight ( 76) = 190.23
  op_atomWeight ( 77) = 192.217
  op_atomWeight ( 78) = 195.084
  op_atomWeight ( 79) = 196.967
  op_atomWeight ( 80) = 200.59

  op_atomWeight ( 81) = 204.383
  op_atomWeight ( 82) = 207.2
  op_atomWeight ( 83) = 208.98
  op_atomWeight ( 84) = 209
  op_atomWeight ( 85) = 210
  op_atomWeight ( 86) = 222
  op_atomWeight ( 87) = 223
  op_atomWeight ( 88) = 226
  op_atomWeight ( 89) = 227
  op_atomWeight ( 90) = 232.038

  op_atomWeight ( 91) = 231.036
  op_atomWeight ( 92) = 238.029
  op_atomWeight ( 93) = 237
  op_atomWeight ( 94) = 244
  op_atomWeight ( 95) = 243
  op_atomWeight ( 96) = 247
  op_atomWeight ( 97) = 247
  op_atomWeight ( 98) = 251
  op_atomWeight ( 99) = 252
  op_atomWeight (100) = 257
!
!
!   ...Ready! 
!
!
  return
end subroutine op_setAtomWeights
