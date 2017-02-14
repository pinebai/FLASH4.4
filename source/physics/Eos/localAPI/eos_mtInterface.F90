!!****ih* source/physics/Eos/localAPI/eos_mtInterface
!!
!! NAME
!!     eos_mtInterface
!!
!! SYNOPSIS
!!     use eos_mtInterface
!!
!! DESCRIPTION
!!
!! This is an interface module for internal use of
!! multiTemp Eos implementations.
!!
!!***

module eos_mtInterface
  implicit none
#include "Eos.h"

  interface
     subroutine eos_byTemp(eos_jlo,eos_jhi,mask,componentMask,eleFrac,ggProdEle)
       integer, intent(IN) :: eos_jlo, eos_jhi
       logical,optional, dimension(EOS_VARS+1:EOS_NUM),INTENT(in)::mask
       integer,optional,target, dimension(N_EOS_TEMP),INTENT(in)::componentMask
       real,optional,INTENT(in) :: eleFrac,ggProdEle
     end subroutine eos_byTemp
  end interface

  interface
     subroutine eos_byTempMG(eos_jlo,eos_jhi,gamM1Ion,mask,componentMask,eleFrac,ggProdEle)
       integer, intent(IN) :: eos_jlo, eos_jhi
       real, intent(IN) :: gamM1Ion(eos_jlo:eos_jhi)
       logical,optional, dimension(EOS_VARS+1:EOS_NUM),INTENT(in)::mask
       integer,optional,target, dimension(N_EOS_TEMP),INTENT(in)::componentMask
       real,optional,INTENT(in) :: eleFrac,ggProdEle
     end subroutine eos_byTempMG
  end interface

  interface
     subroutine eos_multiTypeByTemp(mode,vecLen,eosData,vecBegin,vecEnd,massFrac,mask,componentMask)
       implicit none
       integer, INTENT(in) :: mode, vecLen
       real,INTENT(inout), dimension(EOS_NUM*vecLen) :: eosData 
       integer, INTENT(in) :: vecBegin,vecEnd
       real, optional, INTENT(in),dimension(:)    :: massFrac
       logical,optional,target, dimension(EOS_VARS+1:EOS_NUM),INTENT(in)::mask
       integer,optional, dimension(N_EOS_TEMP),INTENT(in)::componentMask
     end subroutine eos_multiTypeByTemp
  end interface

  interface
     subroutine eos_initPhysData()
     end subroutine eos_initPhysData
  end interface
end module eos_mtInterface
