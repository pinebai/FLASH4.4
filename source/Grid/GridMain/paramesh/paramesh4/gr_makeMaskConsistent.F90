!!****if* source/Grid/GridMain/paramesh/paramesh4/gr_makeMaskConsistent
!!
!! NAME
!!
!!  gr_makeMaskConsistent
!!
!! SYNOPSIS
!!
!!  call gr_makeMaskConsistent(integer,intent(IN)  :: gridDataStruct,
!!                             integer,intent(IN)  :: eosMode,
!!                             logical,intent(OUT)  :: needEos)
!!
!! DESCRIPTION
!!
!!  For the solvers native to FLASH, this routine provides the service of 
!!  ensuring that none of essential variables are masked out. For example 
!!  if mask value is true for a variable that is determined
!!  through application of Eos, then it will be made sure that the inputs to 
!!  Eos also get their values filled. The output parameter "needEos" indicates 
!!  to the calling routine that the calculated guardcells may not be 
!!  thermodynamically consistent, and the calling routine should apply Eos. 
!!
!!  This subroutine takes account of various runtime parameters and 
!!  configuration symbols in order to determine the correct logic, 
!!  including names and existence of grid variables,
!!      number of mass fractions, the Eos mode of the caller,
!!      whether EOS_YE mode (formerly known as "EOS lite") is used,
!!      convertToConsvdForMeshCalls or convertToConsvdInMeshInterp,
!!      and gridinterpolation
!!
!!  The logic here is conservative in the following sense: 
!!  Even if an EOS call could have only small side effects on variables that are
!!  basically input (like ENER_VAR or EINT_VAR in MODE_DENS_EI), or if an EOS call 
!!  change a variable that is basically an input to it only in unusual 
!!  circumstances (as, e.g., when ENER_VAR or EINT_VAR
!!  are inconsistent with each other or other variables), then this routine will
!!  determine that the Eos call should be made.
!!  
!!
!! ARGUMENTS
!!
!!   gridDataStruct : indicates a variable that the caller needs, as an index into unk
!!
!!   eosMode : the eosMode which the calling routine is using.
!!
!!   needEos :  indicates whether Eos should be called after
!!                 guard cell filling.
!!
!! NOTES
!!
!!  Currently this routine is applicable to only cell centered data, since
!!  none of the FLASH solvers need consistency check for face centered data.
!! 
!!  This routine is MESSY MESSY MESSY and very hard to read.
!!  If you don't trust it, turn off runtime parameter enableMaskedGcFill
!!  and worry no more about missing some variables in guard cell filling;
!!  however, the application may then spend some more time in PARAMESH Grid calls
!!  (including MPI calls), exchange more data between processing entities than
!!  necessary, and require larger communication buffers.
!!
!! HISTORY
!!  gr_fillGcMask  started  Jan-2007  KW
!!  some tuning             Jun-2007  KW      
!!  renamed gr_makeMaskConsistent Aug-2007  KW
!!  allow for missing EINT_VAR    Nov-2007  KW
!!  split off less-strict version Dec-2007  KW
!!  changed to gr_makeMaskConsistent December-2007 AD
!!  changes to take 3T variables into account  December-2011 KW
!!  introduced call to Grid_guardCellMaskHook September-2013 KW
!!
!!***


subroutine gr_makeMaskConsistent(gridDataStruct,eosMode,needEos)

  use Eos_interface, ONLY : Eos_getParameters
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface,   ONLY : Grid_guardCellMaskHook
  use Grid_data, ONLY : gr_vartypes, gr_convertToConsvdForMeshCalls,gr_convertToConsvdInMeshInterp


  use physicaldata, ONLY : gcell_on_cc
  implicit none

#include "constants.h"
#include "Flash.h"

  integer,intent(IN) :: gridDataStruct
  integer,intent(IN) :: eosMode
  logical,intent(INOUT) :: needEos


  integer :: iref

  logical :: EosOutput, RelevantGammaC, EnergyVar
  logical,save :: inputsAreUnchangedByEos,inputTempIsGuessForEos,constantGammaC,&
       inputMassFracNeededByEos
  real,save :: smallE = 0.0     !! used sort-of like a 'first_time' flag here; 
                                !! a valid smallE should be .ne. 0.


  logical :: dens=.true.


  !! The next sections processes the defined variables in the current simulation
  !! to determine the situations that need to be examined to make the mask 
  !! consistent. The defined constants that provide indices into the grid
  !! data structures have to carefully handled otherwise the simulations that
  !! haven't included the corresponding variables will fail at compile time.
  !! Doing them separately in a section makes sure that the part of code that 
  !! is examining the conditions is not interspersed with a number of ifdefs,
  !! or even too many conditional, making it cumbersome to read. This is 
  !! especially true of the HAVE_UNK_.... It eliminates the need to check for
  !! the existence of each individual EOS related variable before using it in
  !! the code.


  !! Density, pressure, temperature and energy must be defined in order to
  !! use the equation of state. If any of them in undefined, Eos cannot be used

#define HAVE_UNK_VARS_FOR_EOS
#ifndef TEMP_VAR
#undef HAVE_UNK_VARS_FOR_EOS
#endif
#ifndef PRES_VAR
#undef HAVE_UNK_VARS_FOR_EOS
#endif
#ifndef ENER_VAR
#undef HAVE_UNK_VARS_FOR_EOS
#endif
#ifndef DENS_VAR
#undef HAVE_UNK_VARS_FOR_EOS
#define DENS_VAR 1
  dens=.false.
#endif

  !! If all the EOS variables are defined, that is the simulation is likely to
  !! be using EOS, certain other important EOS related variables must be
  !! handled. The undefined variables are given the index of a defined variable
  !! that will have the same value in the mask. For example if mask(VELX_VAR)
  !! is true, the mask(VELY_VAR) and mask(VELZ_VAR) must be true too. If they
  !! don't exist, they are made to point to VELX_VAR

#ifdef HAVE_UNK_VARS_FOR_EOS

#ifndef EINT_VAR
#define EINT_VAR ENER_VAR
#endif

#ifndef VELY_VAR
#define VELY_VAR VELX_VAR
#endif

#ifndef VELZ_VAR
#define VELZ_VAR VELX_VAR
#endif

#ifdef USE_EOS_YE
#ifndef YE_MSCALAR
  call Driver_abortFlash("gr_makeMaskConsistent compiled in USE_EOS_YE mode, but no YE_MSCALAR is defined")
#endif

#endif
     !! endif for HAVE_UNK...
#endif  


#ifdef TELE_VAR
#define HAVE_UNK_VARS_FOR_3T
#endif

  needEos=.false.

  if((gridDataStruct/=CENTER).and.(gridDataStruct/=CENTER_FACES)) return
  
  !! First find out if there are variables that are output from an Eos call.

#ifdef HAVE_UNK_VARS_FOR_EOS
  if (smallE .EQ. 0.0) then
     call Eos_getParameters(inputsAreUnchanged=inputsAreUnchangedByEos,&
          inputTempIsGuess=inputTempIsGuessForEos,&
          constantGammaC=constantGammaC,&
          inputMassFracNeeded=inputMassFracNeededByEos,&
          smalle=smallE)
  end if

#ifdef HAVE_UNK_VARS_FOR_3T
!!$  select case (eosMode)
!!$     case(MODE_DENS_TEMP)
!!$     case(MODE_DENS_PRES)
!!$     case(MODE_DENS_PRES)
!!$  end select
  if (gcell_on_cc(TEMP_VAR)) then
     needEos =.not.(eosMode==MODE_EOS_NOP   .OR. &
                    eosMode==MODE_DENS_TEMP .OR. &
                    eosMode==MODE_DENS_TEMP_ALL .OR. &
                    eosMode==MODE_DENS_TEMP_EQUI .OR. &
                    eosMode==MODE_DENS_TEMP_COMP .OR. &
                    eosMode==MODE_DENS_TEMP_ION .OR. &
                    eosMode==MODE_DENS_TEMP_ELE .OR. &
                    eosMode==MODE_DENS_TEMP_RAD .OR. &
                    eosMode==MODE_DENS_PRES_COMP .OR. &
                    eosMode==MODE_DENS_PRES_ION .OR. &
                    eosMode==MODE_DENS_PRES_ELE .OR. &
                    eosMode==MODE_DENS_PRES_RAD .OR. &
                    eosMode==MODE_DENS_EI_COMP .OR. &
                    eosMode==MODE_DENS_EI_ION .OR. &
                    eosMode==MODE_DENS_EI_ELE .OR. &
                    eosMode==MODE_DENS_EI_RAD) .OR. needEos
  end if
  if (gcell_on_cc(TION_VAR)) then
     needEos =.not.(eosMode==MODE_EOS_NOP   .OR. &
                    eosMode==MODE_DENS_TEMP .OR. &
                    eosMode==MODE_DENS_TEMP_ALL .OR. &
                    eosMode==MODE_DENS_TEMP_GATHER .OR. &
                    eosMode==MODE_DENS_TEMP_COMP .OR. &
                    eosMode==MODE_DENS_TEMP_ION .OR. &
                    eosMode==MODE_DENS_TEMP_ELE .OR. &
                    eosMode==MODE_DENS_TEMP_RAD .OR. &
                    eosMode==MODE_DENS_PRES .OR. &
                    eosMode==MODE_DENS_PRES_ELE .OR. &
                    eosMode==MODE_DENS_PRES_RAD .OR. &
                    eosMode==MODE_DENS_EI .OR. &
                    eosMode==MODE_DENS_EI_ELE .OR. &
                    eosMode==MODE_DENS_EI_RAD) .OR. needEos
  end if
  if (gcell_on_cc(TELE_VAR)) then
     needEos =.not.(eosMode==MODE_EOS_NOP   .OR. &
                    eosMode==MODE_DENS_TEMP .OR. &
                    eosMode==MODE_DENS_TEMP_ALL .OR. &
                    eosMode==MODE_DENS_TEMP_GATHER .OR. &
                    eosMode==MODE_DENS_TEMP_COMP .OR. &
                    eosMode==MODE_DENS_TEMP_ION .OR. &
                    eosMode==MODE_DENS_TEMP_ELE .OR. &
                    eosMode==MODE_DENS_TEMP_RAD .OR. &
                    eosMode==MODE_DENS_PRES .OR. &
                    eosMode==MODE_DENS_PRES_ION .OR. &
                    eosMode==MODE_DENS_PRES_RAD .OR. &
                    eosMode==MODE_DENS_EI .OR. &
                    eosMode==MODE_DENS_EI_ION .OR. &
                    eosMode==MODE_DENS_EI_RAD) .OR. needEos
  end if
  if (gcell_on_cc(TRAD_VAR)) then
     needEos =.not.(eosMode==MODE_EOS_NOP   .OR. &
                    eosMode==MODE_DENS_TEMP .OR. &
                    eosMode==MODE_DENS_TEMP_ALL .OR. &
                    eosMode==MODE_DENS_TEMP_GATHER .OR. &
                    eosMode==MODE_DENS_TEMP_COMP .OR. &
                    eosMode==MODE_DENS_TEMP_ION .OR. &
                    eosMode==MODE_DENS_TEMP_ELE .OR. &
                    eosMode==MODE_DENS_TEMP_RAD .OR. &
                    eosMode==MODE_DENS_PRES .OR. &
                    eosMode==MODE_DENS_PRES_ION .OR. &
                    eosMode==MODE_DENS_PRES_ELE .OR. &
                    eosMode==MODE_DENS_EI .OR. &
                    eosMode==MODE_DENS_EI_ION .OR. &
                    eosMode==MODE_DENS_EI_ELE) .OR. needEos
  end if
  if (gcell_on_cc(PRES_VAR)) then
     needEos =.not.(eosMode==MODE_EOS_NOP   .OR. &
                    eosMode==MODE_DENS_PRES .OR. &
                    eosMode==MODE_DENS_PRES_ALL .OR. &
                    eosMode==MODE_DENS_PRES_COMP .OR. &
                    eosMode==MODE_DENS_PRES_ION .OR. &
                    eosMode==MODE_DENS_PRES_ELE .OR. &
                    eosMode==MODE_DENS_PRES_RAD .OR. &
                    eosMode==MODE_DENS_TEMP_COMP .OR. &
                    eosMode==MODE_DENS_TEMP_ION .OR. &
                    eosMode==MODE_DENS_TEMP_ELE .OR. &
                    eosMode==MODE_DENS_TEMP_RAD .OR. &
                    eosMode==MODE_DENS_EI_COMP .OR. &
                    eosMode==MODE_DENS_EI_ION .OR. &
                    eosMode==MODE_DENS_EI_ELE .OR. &
                    eosMode==MODE_DENS_EI_RAD) .OR. needEos
  end if
  if (gcell_on_cc(PION_VAR)) then
     needEos =.not.(eosMode==MODE_EOS_NOP   .OR. &
                    eosMode==MODE_DENS_PRES .OR. &
                    eosMode==MODE_DENS_PRES_ALL .OR. &
                    eosMode==MODE_DENS_PRES_COMP .OR. &
                    eosMode==MODE_DENS_PRES_ION .OR. &
                    eosMode==MODE_DENS_PRES_ELE .OR. &
                    eosMode==MODE_DENS_PRES_RAD .OR. &
                    eosMode==MODE_DENS_TEMP .OR. &
                    eosMode==MODE_DENS_TEMP_ELE .OR. &
                    eosMode==MODE_DENS_TEMP_RAD .OR. &
                    eosMode==MODE_DENS_EI .OR. &
                    eosMode==MODE_DENS_EI_ELE .OR. &
                    eosMode==MODE_DENS_EI_RAD) .OR. needEos
  end if
  if (gcell_on_cc(PELE_VAR)) then
     needEos =.not.(eosMode==MODE_EOS_NOP   .OR. &
                    eosMode==MODE_DENS_PRES .OR. &
                    eosMode==MODE_DENS_PRES_ALL .OR. &
                    eosMode==MODE_DENS_PRES_COMP .OR. &
                    eosMode==MODE_DENS_PRES_ION .OR. &
                    eosMode==MODE_DENS_PRES_ELE .OR. &
                    eosMode==MODE_DENS_PRES_RAD .OR. &
                    eosMode==MODE_DENS_TEMP .OR. &
                    eosMode==MODE_DENS_TEMP_ION .OR. &
                    eosMode==MODE_DENS_TEMP_RAD .OR. &
                    eosMode==MODE_DENS_EI .OR. &
                    eosMode==MODE_DENS_EI_ION .OR. &
                    eosMode==MODE_DENS_EI_RAD) .OR. needEos
  end if
  if (gcell_on_cc(PRAD_VAR)) then
     needEos =.not.(eosMode==MODE_EOS_NOP   .OR. &
                    eosMode==MODE_DENS_PRES .OR. &
                    eosMode==MODE_DENS_PRES_ALL .OR. &
                    eosMode==MODE_DENS_PRES_COMP .OR. &
                    eosMode==MODE_DENS_PRES_ION .OR. &
                    eosMode==MODE_DENS_PRES_ELE .OR. &
                    eosMode==MODE_DENS_PRES_RAD .OR. &
                    eosMode==MODE_DENS_TEMP .OR. &
                    eosMode==MODE_DENS_TEMP_ION .OR. &
                    eosMode==MODE_DENS_TEMP_ELE .OR. &
                    eosMode==MODE_DENS_EI .OR. &
                    eosMode==MODE_DENS_EI_ION .OR. &
                    eosMode==MODE_DENS_EI_ELE) .OR. needEos
  end if
  if (gcell_on_cc(EION_VAR)) then
     needEos =.not.(eosMode==MODE_EOS_NOP .OR. &
                    eosMode==MODE_DENS_EI .OR. &
                    eosMode==MODE_DENS_EI_ALL .OR. &
                    eosMode==MODE_DENS_EI_GATHER .OR. &
                    eosMode==MODE_DENS_EI_COMP .OR. &
                    eosMode==MODE_DENS_EI_ION .OR. &
                    eosMode==MODE_DENS_EI_ELE .OR. &
                    eosMode==MODE_DENS_EI_RAD .OR. &
                    eosMode==MODE_DENS_TEMP .OR. &
                    eosMode==MODE_DENS_TEMP_ELE .OR. &
                    eosMode==MODE_DENS_TEMP_RAD .OR. &
                    eosMode==MODE_DENS_PRES .OR. &
                    eosMode==MODE_DENS_PRES_ELE .OR. &
                    eosMode==MODE_DENS_PRES_RAD) .OR. needEos
  end if
  if (gcell_on_cc(EELE_VAR)) then
     needEos =.not.(eosMode==MODE_EOS_NOP .OR. &
                    eosMode==MODE_DENS_EI .OR. &
                    eosMode==MODE_DENS_EI_ALL .OR. &
                    eosMode==MODE_DENS_EI_GATHER .OR. &
                    eosMode==MODE_DENS_EI_COMP .OR. &
                    eosMode==MODE_DENS_EI_ION .OR. &
                    eosMode==MODE_DENS_EI_ELE .OR. &
                    eosMode==MODE_DENS_EI_RAD .OR. &
                    eosMode==MODE_DENS_TEMP .OR. &
                    eosMode==MODE_DENS_TEMP_ION .OR. &
                    eosMode==MODE_DENS_TEMP_RAD .OR. &
                    eosMode==MODE_DENS_PRES .OR. &
                    eosMode==MODE_DENS_PRES_ION .OR. &
                    eosMode==MODE_DENS_PRES_RAD) .OR. needEos
  end if
  if (gcell_on_cc(ERAD_VAR)) then
     needEos =.not.(eosMode==MODE_EOS_NOP .OR. &
                    eosMode==MODE_DENS_EI .OR. &
                    eosMode==MODE_DENS_EI_ALL .OR. &
                    eosMode==MODE_DENS_EI_GATHER .OR. &
                    eosMode==MODE_DENS_EI_COMP .OR. &
                    eosMode==MODE_DENS_EI_ION .OR. &
                    eosMode==MODE_DENS_EI_ELE .OR. &
                    eosMode==MODE_DENS_EI_RAD .OR. &
                    eosMode==MODE_DENS_TEMP .OR. &
                    eosMode==MODE_DENS_TEMP_ION .OR. &
                    eosMode==MODE_DENS_TEMP_ELE .OR. &
                    eosMode==MODE_DENS_PRES .OR. &
                    eosMode==MODE_DENS_PRES_ION .OR. &
                    eosMode==MODE_DENS_PRES_ELE) .OR. needEos
  end if
#else  
  !! If one of the output variables in the current eosMode is set to true
  !! then Eos call is needed.
  if (eosMode .NE. MODE_EOS_NOP) then
     needEos = (.not.(eosMode==MODE_DENS_TEMP)).and.(gcell_on_cc(TEMP_VAR))
     needEos = (.not.(eosMode==MODE_DENS_PRES)).and.(gcell_on_cc(PRES_VAR)).or.needEos
  end if
#endif

  !! If any energy variable is set to true then eos call is needed DEV:: why ?
  needEos= gcell_on_cc(ENER_VAR).or.gcell_on_cc(EINT_VAR).or.needEos 
  
  !! If any of the two gammas are true then eos call is needed
  needEos= (gcell_on_cc(GAMC_VAR) .AND. .NOT. constantGammaC).or.&
       gcell_on_cc(GAME_VAR).or.needEos
  
  if(needEos) then
     
     !! if eos call is needed, density is always input, so must be true
     gcell_on_cc(DENS_VAR)=.true. 
     
     !! temp is true if either it is input to eos or is used as initial guess
     gcell_on_cc(TEMP_VAR)=gcell_on_cc(TEMP_VAR).or.inputTempIsGuessForEos.or.&
          (eosMode==MODE_DENS_TEMP)
     
#ifdef HAVE_UNK_VARS_FOR_3T
     gcell_on_cc(TEMP_VAR)=gcell_on_cc(TEMP_VAR).or. &
          eosMode==MODE_DENS_TEMP_EQUI .OR. &
          eosMode==MODE_DENS_TEMP_ALL
     gcell_on_cc(TION_VAR)=gcell_on_cc(TION_VAR).or. &
          eosMode==MODE_DENS_TEMP_GATHER .OR. &
          eosMode==MODE_DENS_TEMP_ALL .OR. &
          eosMode==MODE_DENS_TEMP_COMP .OR. &
          eosMode==MODE_DENS_TEMP_ION .OR. &
          (inputTempIsGuessForEos .AND. &
             (eosMode==MODE_DENS_PRES_ALL .OR. &
              eosMode==MODE_DENS_PRES_COMP .OR. &
              eosMode==MODE_DENS_PRES_ION .OR. &
              eosMode==MODE_DENS_EI_ALL .OR. &
              eosMode==MODE_DENS_EI_GATHER .OR. &
              eosMode==MODE_DENS_EI_MAT_GATHER .OR. &
              eosMode==MODE_DENS_EI_MAT_GATHER_PRADSCALE .OR. &
              eosMode==MODE_DENS_EI_RECAL_GATHER .OR. &
              eosMode==MODE_DENS_EI_COMP .OR. &
              eosMode==MODE_DENS_EI_MAT_EQUI .OR. &
              eosMode==MODE_DENS_EI_ION))
     gcell_on_cc(TELE_VAR)=gcell_on_cc(TELE_VAR).or. &
          eosMode==MODE_DENS_TEMP_GATHER .OR. &
          eosMode==MODE_DENS_TEMP_ALL .OR. &
          eosMode==MODE_DENS_TEMP_COMP .OR. &
          eosMode==MODE_DENS_TEMP_ELE .OR. &
          (inputTempIsGuessForEos .AND. &
             (eosMode==MODE_DENS_PRES_ALL .OR. &
              eosMode==MODE_DENS_PRES_COMP .OR. &
              eosMode==MODE_DENS_PRES_ELE .OR. &
              eosMode==MODE_DENS_EI_ALL .OR. &
              eosMode==MODE_DENS_EI_GATHER .OR. &
              eosMode==MODE_DENS_EI_MAT_GATHER .OR. &
              eosMode==MODE_DENS_EI_MAT_GATHER_PRADSCALE .OR. &
              eosMode==MODE_DENS_EI_RECAL_GATHER .OR. &
              eosMode==MODE_DENS_EI_COMP .OR. &
              eosMode==MODE_DENS_EI_MAT_EQUI .OR. &
              eosMode==MODE_DENS_EI_ELE))
     gcell_on_cc(TRAD_VAR)=gcell_on_cc(TRAD_VAR).or. &
          eosMode==MODE_DENS_TEMP_GATHER .OR. &
          eosMode==MODE_DENS_TEMP_ALL .OR. &
          eosMode==MODE_DENS_TEMP_COMP .OR. &
          eosMode==MODE_DENS_TEMP_RAD
#endif

     !! pressure is true if input to eos
     gcell_on_cc(PRES_VAR)=gcell_on_cc(PRES_VAR).or.(eosMode==MODE_DENS_PRES)
     
#ifdef HAVE_UNK_VARS_FOR_3T
     gcell_on_cc(PRES_VAR)=gcell_on_cc(PRES_VAR).or. &
          eosMode==MODE_DENS_PRES_ALL
     gcell_on_cc(PION_VAR)=gcell_on_cc(PION_VAR).or. &
          eosMode==MODE_DENS_PRES_ALL .OR. &
          eosMode==MODE_DENS_PRES_COMP .OR. &
          eosMode==MODE_DENS_PRES_ION
     gcell_on_cc(PELE_VAR)=gcell_on_cc(PELE_VAR).or. &
          eosMode==MODE_DENS_PRES_ALL .OR. &
          eosMode==MODE_DENS_PRES_COMP .OR. &
          eosMode==MODE_DENS_PRES_ELE
     gcell_on_cc(PRAD_VAR)=gcell_on_cc(PRAD_VAR).or. &
          eosMode==MODE_DENS_PRES_ALL .OR. &
          eosMode==MODE_DENS_PRES_COMP .OR. &
          eosMode==MODE_DENS_PRES_RAD
#endif

     !! Generic true for any kind of Eos
     gcell_on_cc(EINT_VAR) = .TRUE.
#ifdef HAVE_UNK_VARS_FOR_3T
     gcell_on_cc(EION_VAR)=gcell_on_cc(EION_VAR).or. &
          eosMode==MODE_DENS_EI_ALL .OR. &
          eosMode==MODE_DENS_EI_GATHER .OR. &
          eosMode==MODE_DENS_EI_MAT_GATHER .OR. &
          eosMode==MODE_DENS_EI_MAT_GATHER_PRADSCALE .OR. &
          eosMode==MODE_DENS_EI_RECAL_GATHER .OR. &
          eosMode==MODE_DENS_EI_COMP .OR. &
          eosMode==MODE_DENS_EI_MAT_EQUI .OR. &
          eosMode==MODE_DENS_EI_ION
     gcell_on_cc(EELE_VAR)=gcell_on_cc(EELE_VAR).or. &
          eosMode==MODE_DENS_EI_ALL .OR. &
          eosMode==MODE_DENS_EI_GATHER .OR. &
          eosMode==MODE_DENS_EI_MAT_GATHER .OR. &
          eosMode==MODE_DENS_EI_MAT_GATHER_PRADSCALE .OR. &
          eosMode==MODE_DENS_EI_RECAL_GATHER .OR. &
          eosMode==MODE_DENS_EI_COMP .OR. &
          eosMode==MODE_DENS_EI_MAT_EQUI .OR. &
          eosMode==MODE_DENS_EI_ELE
     gcell_on_cc(ERAD_VAR)=gcell_on_cc(ERAD_VAR).or. &
          eosMode==MODE_DENS_EI_ALL .OR. &
          eosMode==MODE_DENS_EI_GATHER .OR. &
          eosMode==MODE_DENS_EI_MAT_GATHER .OR. &
          eosMode==MODE_DENS_EI_MAT_GATHER_PRADSCALE .OR. &
          eosMode==MODE_DENS_EI_RECAL_GATHER .OR. &
          eosMode==MODE_DENS_EI_COMP .OR. &
          eosMode==MODE_DENS_EI_MAT_EQUI .OR. &
          eosMode==MODE_DENS_EI_RAD
#ifdef FLLM_VAR
     gcell_on_cc(FLLM_VAR)=gcell_on_cc(FLLM_VAR).or. &
          eosMode==MODE_DENS_EI_MAT_GATHER .OR. &
          eosMode==MODE_DENS_EI_MAT_GATHER_PRADSCALE
#endif
#endif
     gcell_on_cc(ENER_VAR) = .TRUE.
     gcell_on_cc(VELX_VAR) = .TRUE.
     gcell_on_cc(VELY_VAR) = .TRUE.
     gcell_on_cc(VELZ_VAR) = .TRUE.
     
     !! Relevant only when using helmholtz eos in either form
#ifdef USE_EOS_YE
     !! cal says abar=1/sumy
     !! cal says zbar=ye / sumy and he claims sumy are never zero
     gcell_on_cc(SUMY_MSCALAR) = .TRUE.
     gcell_on_cc(YE_MSCALAR) = .TRUE.
#else  
     !! if USE_EOS_YE not defined
#if defined(SUMY_MSCALAR) && defined(YE_MSCALAR)
     gcell_on_cc(SUMY_MSCALAR) = .TRUE.
     gcell_on_cc(YE_MSCALAR) = .TRUE.
#endif
     !! mass fractions are needed because they are needed to get the density right
     if (inputMassFracNeededByEos) gcell_on_cc(SPECIES_BEGIN:SPECIES_END) = .TRUE.
#endif


  end if  !! finish the handling of Eos needs in applying mask
  
#endif  

  !! Special stuff that may be needed by specific included code units is done here.
  !! The advantage is that that gr_makeMaskConsistent itself does not have
  !! to be updated for every EOS variant.
  call Grid_guardCellMaskHook(gcell_on_cc,needEos)



  !! now handle other situations that may arise when density is included 
  !! in the simulation
  if(dens) then
     if(.not.gcell_on_cc(DENS_VAR)) then
        if (gr_convertToConsvdForMeshCalls .OR. gr_convertToConsvdInMeshInterp) then
           do iref=1,NUNK_VARS
              !! For mass-specific variable u, interpolation may be done on the 
              !! product u*rho, so density is needed, too.
              gcell_on_cc(DENS_VAR)=gcell_on_cc(DENS_VAR).or.&
                   (gcell_on_cc(iref).and.(gr_vartypes(iref)==VARTYPE_PER_MASS))
           end do
        end if
     end if
  endif
  

#ifdef GRID_WITH_MONOTONIC

  if (SPECIES_BEGIN .LT. SPECIES_END) then
     if (gr_convertToConsvdForMeshCalls .OR. gr_convertToConsvdInMeshInterp) then
        if (dens) then
           if(gcell_on_cc(DENS_VAR)) gcell_on_cc(SPECIES_BEGIN:SPECIES_END) = .TRUE.
        end if
     else
        if( any(gcell_on_cc(SPECIES_BEGIN:SPECIES_END)))&
             gcell_on_cc(SPECIES_BEGIN:SPECIES_END) = .TRUE.
     end if
  end if

#endif
end subroutine gr_makeMaskConsistent
