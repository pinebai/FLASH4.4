!!****if* source/physics/Eos/unitTest/eos_initTest
!!
!! NAME
!!
!!  eos_initTest
!!
!!
!! SYNOPSIS
!!
!!  call eos_initTest()
!!  
!!
!! DESCRIPTION
!! 
!!  Initialize Eos unitTest scope variables which are typically from runtime parameters.
!!  The part of the unitTest located in the Simulation unit tree should arrange for this
!!  initialization routine to be called.
!!
!! ARGUMENTS
!!
!!
!!  none
!!  
!!
!! PARAMETERS
!!
!!
!!    eos_testPresMode [STRING]
!!    eos_testEintMode [STRING]
!!    eos_testTempMode [STRING]
!!***

subroutine eos_initTest()

  use RuntimeParameters_interface, ONLY: RuntimeParameters_get, &
                                         RuntimeParameters_mapStrToInt
  use eos_testData, ONLY: eos_testPresModeStr, &
                          eos_testEintModeStr, &
                          eos_testTempModeStr, &
                          eos_testPresMode, &
                          eos_testEintMode, &
                          eos_testTempMode
                          

  implicit none

  call RuntimeParameters_get ("eos_testPresMode", eos_testPresModeStr)
  call RuntimeParameters_mapStrToInt(eos_testPresModeStr, eos_testPresMode)
  call RuntimeParameters_get ("eos_testEintMode", eos_testEintModeStr)
  call RuntimeParameters_mapStrToInt(eos_testEintModeStr, eos_testEintMode)
  call RuntimeParameters_get ("eos_testTempMode", eos_testTempModeStr)
  call RuntimeParameters_mapStrToInt(eos_testTempModeStr, eos_testTempMode)


end subroutine eos_initTest
