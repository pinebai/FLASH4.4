!!****ih* source/physics/Eos/unitTest/eos_testData
!!
!! NAME
!!
!!  eos_testData
!!
!! 
!! SYNOPSIS
!!
!! use eos_testData
!!
!! DESCRIPTION
!!
!!***

module eos_testData

  implicit none

#include "constants.h"

  character(len=MAX_STRING_LENGTH), save :: eos_testPresModeStr, eos_testEintModeStr, eos_testTempModeStr
  integer, save :: eos_testPresMode, eos_testEintMode, eos_testTempMode

end module eos_testData
