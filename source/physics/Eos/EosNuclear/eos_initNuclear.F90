!!****if* source/physics/Eos/EosNuclear/eos_initNuclear
!!
!! NAME
!!
!!  eos_initNuclear
!!
!! SYNOPSIS
!!  
!!  subroutine eos_initNuclear()
!!                 
!!
!! DESCRIPTION
!!
!!  Initialization for the Nuclear EOS appropriate for 
!!  core-collapse supernova simulations.  
!!
!! ARGUMENTS
!!
!! NOTES
!!      Parts of this unit are released under a different license than the
!!      usual FLASH license.  Specifically, some subroutines in the kernel 
!!      directory are released under the Creative Commons 
!!      attribution-noncommercial-share alike license.  Basically, if you use this
!!      unit in your work, the license requires that you cite the two articles 
!!      mentioned below.  More details may be found here:  
!!      stellarcollapse.org/equationofstate.
!!
!!      * O'Connor, E.P., & Ott, C.D. 2010, CQGra, 27, 114103
!!      * Couch, S.M. 2013, ApJ, 765, 29
!!
!!
!!***

#include "Eos.h"
#include "constants.h"
subroutine eos_initNuclear()

  use Eos_data, ONLY : eos_type
  use eosmodule
  use eos_nucData
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get
  use Driver_interface, ONLY : Driver_getComm, Driver_getMype
  use IO_interface, ONLY : IO_setScalar, IO_getPrevScalar
  use Logfile_interface, ONLY : Logfile_stamp

  implicit none
  character(len=100) :: message
  integer :: error

  eos_type=EOS_NUC

  call Driver_getComm(MESH_COMM,eos_meshComm)
  call Driver_getMype(MESH_COMM,eos_meshMe)

  call PhysicalConstants_get("Avogadro", avo)

  call RuntimeParameters_get('eos_file', eos_file)
  call RuntimeParameters_get('restart', eos_restart)
  if (eos_restart) then
     call IO_getPrevScalar("postBounce", eos_postBounce, error)
     call IO_getPrevScalar("bounceTime", eos_bounceTime, error)
     if (error /= 0) then
        if (eos_meshMe == MASTER_PE) then
           write(message,*) "postBounce not in checkpoint. Using runtime parameters."
           print *, message
           call Logfile_stamp(message, '[eos_initNuclear]')
        end if
        call RuntimeParameters_get("postBounce", eos_postBounce)
        call RuntimeParameters_get("bounceTime", eos_bounceTime)
     end if
  else
     call RuntimeParameters_get("postBounce", eos_postBounce)
     call RuntimeParameters_get("bounceTime", eos_bounceTime)
     call IO_setScalar("postBounce", eos_postBounce)
     call IO_setScalar("bounceTime", eos_bounceTime)
  end if
  if (eos_meshMe == MASTER_PE .AND. eos_postBounce) then
     write(message,*) "Starting sim post-bounce! Bounce time = ", eos_bounceTime
     print *, message
     call Logfile_stamp(message, '[eos_initNuclear]')
  end if

  call readtable(eos_file)
  e_zeroPoint = energy_shift


  return
end subroutine eos_initNuclear
