!!****ih* source/physics/Eos/EosMain/multiTemp/Helmholtz/eos_helmInterface
!!
!! NAME
!!     eos_helmInterface
!!
!! SYNOPSIS
!!     use eos_helmInterface
!!
!! DESCRIPTION
!!
!! This is an interface module for internal use of
!! the Helmholtz Eos implementation.
!!
!!***

module eos_helmInterface
  implicit none
#include "Eos.h"

  interface
     subroutine eos_helm(eos_jlo,eos_jhi,mask,componentMask)
       integer, intent(IN) :: eos_jlo, eos_jhi
       logical,optional, dimension(EOS_VARS+1:EOS_NUM),INTENT(in)::mask
       integer,optional, dimension(N_EOS_TEMP),INTENT(in)::componentMask
     end subroutine eos_helm
  end interface
end module eos_helmInterface
