!!****ih* source/physics/Eos/localAPI/eos_localInterface
!!
!! NAME
!!     eos_localInterface
!!
!! SYNOPSIS
!!     use eos_localInterface
!!
!! DESCRIPTION
!!
!! This is an interface module for internal use of
!! the Eos unit.
!!
!!***

module eos_localInterface
  implicit none
#include "Eos.h"
#include "Flash.h"
  interface
     subroutine eos_idealGamma(mode, vecLen, eosData, vecBegin,vecEnd, massFrac, mask, diagFlag)
       implicit none
       integer, INTENT(in) :: mode, vecLen
       real,INTENT(inout), dimension(EOS_NUM*vecLen) :: eosData 
       integer,optional,INTENT(in) :: vecBegin,vecEnd
       logical, optional, INTENT(in),target,dimension(EOS_VARS+1:EOS_NUM) :: mask
       real, optional, INTENT(in),dimension(NSPECIES*vecLen)    :: massFrac
       integer, optional, INTENT(out)    :: diagFlag
     end subroutine Eos_idealGamma
  end interface
  interface
     subroutine eos_idealGamma3T(mode, vecLen, eosData, vecBegin,vecEnd, &
          eosType, material, massFrac, mask, componentMask)
       implicit none
       integer, INTENT(in) :: mode, vecLen
       real,INTENT(inout), dimension(EOS_NUM*vecLen) :: eosData 
       integer,optional,INTENT(in) :: vecBegin,vecEnd
       integer,optional,INTENT(in) :: eosType
       integer,optional,INTENT(in) :: material
       logical, optional, INTENT(in),target,dimension(EOS_VARS+1:EOS_NUM) :: mask
       real, optional, INTENT(in),dimension(NSPECIES*vecLen)    :: massFrac
       integer,optional,INTENT(in), dimension(N_EOS_TEMP)::componentMask
     end subroutine Eos_idealGamma3T
  end interface

  interface
     subroutine eos_mgamma(mode, vecLen, eosData, vecBegin,vecEnd, massFrac, mask)
       integer, INTENT(in) :: mode, vecLen
       real,INTENT(inout), dimension(EOS_NUM*vecLen) :: eosData 
       integer,optional,INTENT(in) :: vecBegin,vecEnd
       logical, optional, INTENT(in),dimension(EOS_VARS+1:EOS_NUM) :: mask
       real, optional, INTENT(in),dimension(NSPECIES*vecLen)    :: massFrac
     end subroutine Eos_mgamma
  end interface

  interface
     subroutine eos_mtemp(mode, vecLen, eosData, massFrac, mask)
       integer, INTENT(in) :: mode, vecLen
       real,INTENT(inout), dimension(EOS_NUM*vecLen) :: eosData 
       logical, optional, INTENT(in),target,dimension(EOS_VARS+1:EOS_NUM) :: mask
       real, optional, INTENT(in),dimension(NSPECIES*vecLen)    :: massFrac
     end subroutine eos_mtemp
  end interface

  interface
     subroutine eos_helmholtz(mode, vecLen, eosData, massFrac, mask)
       integer, INTENT(in) :: mode, vecLen
       real,INTENT(inout), dimension(EOS_NUM*vecLen) :: eosData 
       logical, optional, INTENT(in),dimension(EOS_VARS+1:EOS_NUM) :: mask
       real, optional,target, INTENT(in),dimension(NSPECIES*vecLen)    :: massFrac
     end subroutine Eos_helmholtz
  end interface

  interface
     subroutine eos_tabulated(mode, vecLen, eosData, massFrac, mask)
       integer, INTENT(in) :: mode, vecLen
       real,INTENT(inout), dimension(EOS_NUM*vecLen) :: eosData 
       logical, optional, INTENT(in),target,dimension(EOS_VARS+1:EOS_NUM) :: mask
       real, optional, INTENT(in),dimension(NSPECIES*vecLen)    :: massFrac
     end subroutine Eos_tabulated
  end interface
  interface
     subroutine eos_tabIonmix(mode, vecLen, eosData, vecBegin, vecEnd, eosType, subtype, material, mask)
       implicit none
       integer, INTENT(in) :: mode, vecLen
       real,INTENT(inout), dimension(EOS_NUM*vecLen) :: eosData 
       integer,optional,INTENT(in) :: vecBegin,vecEnd
       integer,optional,INTENT(in) :: eosType,subtype
       integer,optional,INTENT(in) :: material
       logical, optional, INTENT(in),target,dimension(EOS_VARS+1:EOS_NUM) :: mask
     end subroutine Eos_tabIonmix
  end interface

  interface
     subroutine eos_nuclear(mode, vecLen, eosData, massFrac, mask)
       integer, INTENT(in) :: mode, vecLen
       real,INTENT(inout), dimension(EOS_NUM*vecLen) :: eosData 
       logical, optional, INTENT(in),target,dimension(EOS_VARS+1:EOS_NUM) :: mask
       real, optional, INTENT(in),dimension(NSPECIES*vecLen)    :: massFrac
     end subroutine eos_nuclear
  end interface

  interface 
     subroutine eos_initGamma()
     end subroutine eos_initGamma
  end interface

  interface 
     subroutine eos_initMgamma()
     end subroutine eos_initMgamma
  end interface

  interface 
     subroutine eos_initMtemp()
     end subroutine eos_initMtemp
  end interface

  interface 
     subroutine eos_initHelmholtz()
     end subroutine eos_initHelmholtz
  end interface

  interface 
     subroutine eos_initTabulated()
     end subroutine eos_initTabulated
  end interface
  interface 
     subroutine eos_tabFinalize()
     end subroutine eos_tabFinalize
  end interface

  interface 
     subroutine eos_initNuclear()
     end subroutine eos_initNuclear
  end interface

  interface 
     subroutine eos_externalComputeAbarZbar(solnScalars, abarData, zbarData)
       implicit none
       real, dimension(:,:), intent(in) :: solnScalars
       real, dimension(:), intent(out)  :: abarData, zbarData
     end subroutine eos_externalComputeAbarZbar
  end interface

end module eos_localInterface
