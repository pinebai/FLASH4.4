Module Eos_interface

  use Eos_nucInterface, ONLY: Eos_nucOneZone, Eos_nucDetectBounce

  implicit none

# include "Eos.h"
# include "Flash.h"
# include "constants.h"
#include "FortranLangFeatures.fh"
  
  interface
     subroutine Eos_guardCells(eosMode, blockID,corners,layers,skipSrl)
       integer,intent(IN) :: eosMode,blockID
       logical,intent(IN) :: corners
       integer,dimension(MDIM),optional,intent(IN) :: layers
       logical,optional, intent(IN) :: skipSrl
     end subroutine Eos_guardCells
  end interface
  
  interface Eos_wrapped
     subroutine Eos_wrapped(mode,range,blockNum,gridDataStruct)
       integer, intent(in) :: mode
       integer, dimension(2,MDIM), intent(in) :: range
       integer,intent(IN) :: blockNum
       integer,optional,intent(IN) :: gridDataStruct
     end subroutine Eos_wrapped
     subroutine Eos_arrayWrapped(mode,range,solnData,gridDataStruct)
       integer, intent(in) :: mode
       integer, dimension(2,MDIM), intent(in) :: range
       real, POINTER_INTENT_IN :: solnData(:,:,:,:)
       integer,optional,intent(IN) :: gridDataStruct
     end subroutine Eos_arrayWrapped
  end interface
  
  interface
     subroutine Eos(mode, vecLen, eosData,  massFrac, mask, vecBegin,vecEnd,diagFlag)
       integer, INTENT(in) :: mode, vecLen
       real,INTENT(inout), dimension(EOS_NUM*vecLen) :: eosData 
       real, optional, INTENT(in),dimension(NSPECIES*vecLen)    :: massFrac
       logical, optional, INTENT(in),target,dimension(EOS_VARS+1:EOS_NUM) :: mask     
       integer,optional,INTENT(in) :: vecBegin,vecEnd
       integer, optional, INTENT(out)    :: diagFlag
     end subroutine Eos
  end interface
  
  interface Eos_init
     subroutine Eos_init()
     end subroutine Eos_init
  end interface
  
  interface Eos_finalize
     subroutine Eos_finalize()
     end subroutine Eos_finalize
  end interface

  interface
     subroutine Eos_getParameters(eintSwitch,inputsAreUnchanged,inputTempIsGuess,constantGammaC,&
          inputMassFracNeeded,smalle,smallE1,smallE2,smallE3)
       real,OPTIONAL,intent(OUT) :: eintSwitch
       logical,OPTIONAL,intent(OUT) :: inputsAreUnchanged
       logical,OPTIONAL,intent(OUT) :: inputTempIsGuess
       logical,OPTIONAL,intent(OUT) :: constantGammaC
       logical,OPTIONAL,intent(OUT) :: inputMassFracNeeded
       real,OPTIONAL,intent(OUT) :: smalle
       real,OPTIONAL,intent(OUT) :: smallE1,smallE2,smallE3
     end subroutine Eos_getParameters
  end interface

  interface Eos_unitTest
     subroutine Eos_unitTest(fileUnit, perfect)
       integer, intent(in) :: fileUnit
       logical, intent(out) :: perfect
     end subroutine Eos_unitTest
  end interface

  interface Eos_getAbarZbar
     subroutine Eos_getAbarZbar(solnVec,abar,zbar,sumY,Ye,massFrac)
       implicit none
       real, OPTIONAL,dimension(NUNK_VARS),intent(IN) :: solnVec
       real, OPTIONAL,                    intent(OUT) :: abar, zbar, Ye, sumY
       real, OPTIONAL,dimension(NSPECIES), intent(IN) :: massFrac
     end subroutine Eos_getAbarZbar
     subroutine Eos_getAbarZbarArraySection(ifirstVar,solnVec,abar,zbar,sumY,Ye,massFrac)
       implicit none
       integer,                            intent(IN) :: ifirstVar
       real, OPTIONAL,                     intent(IN) :: solnVec(ifirstVar:NUNK_VARS)
       real, OPTIONAL,                    intent(OUT) :: abar, zbar, Ye, sumY
       real, OPTIONAL,dimension(NSPECIES), intent(IN) :: massFrac
     end subroutine Eos_getAbarZbarArraySection
  end interface

  interface
     subroutine Eos_getData(axis,pos,vecLen,solnData,gridDataStruct,eosData, massFrac,eosMask)
       implicit none
       integer, intent(in) :: axis, vecLen,gridDataStruct
       integer, dimension(MDIM), intent(in) :: pos
       real, dimension(:),intent(OUT) :: eosData
       real, pointer,dimension(:,:,:,:) :: solnData
       real,dimension(:),optional,intent(OUT) :: massFrac
       logical, optional, INTENT(INOUT),dimension(EOS_VARS+1:) :: eosMask     
     end subroutine Eos_getData
  end interface

  interface Eos_getTempData
     subroutine Eos_getTempData(axis,pos,vecLen,solnData,gridDataStruct,eosData,mode)
       implicit none
       integer, intent(in) :: axis, vecLen, gridDataStruct, mode
       integer, dimension(MDIM), intent(in) :: pos
       real, dimension(:),intent(OUT) :: eosData
       real, pointer:: solnData(:,:,:,:)
     end subroutine Eos_getTempData
     subroutine Eos_getTempDataFromVec(solnVec,eosData,mode)
       implicit none
       integer, intent(in) :: mode
       real, dimension(:),intent(OUT) :: eosData
       real, dimension(NUNK_VARS),intent(IN) :: solnVec
     end subroutine Eos_getTempDataFromVec
  end interface
  interface
     subroutine Eos_getTempDataF(axis,pos,vecLen,solnData,gridDataStruct,eosData,mode)
       implicit none
       integer, intent(in) :: axis, vecLen, gridDataStruct, mode
       integer, dimension(MDIM), intent(in) :: pos
       real, dimension(EOS_NUM*vecLen),intent(OUT) :: eosData
       real, pointer:: solnData(:,:,:,:)
     end subroutine Eos_getTempDataF
  end interface

  interface
     subroutine Eos_putData(axis,pos,vecLen,solnData,gridDataStruct,eosData)
       integer, intent(in) :: axis, vecLen, gridDataStruct
       integer, dimension(MDIM), intent(in) :: pos
       real, dimension(:),intent(IN) :: eosData
       real, pointer,dimension(:,:,:,:) :: solnData
     end subroutine Eos_putData
  end interface

  interface
     subroutine Eos_logDiagnostics(force)
       implicit none
       logical, intent(IN) :: force
     end subroutine Eos_logDiagnostics
  end interface

end Module Eos_interface
  
