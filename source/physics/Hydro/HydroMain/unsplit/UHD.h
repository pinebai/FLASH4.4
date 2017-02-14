!!***h* source/physics/Hydro/HydroMain/unsplit/UHD.h
!!
!! This is the internal header file for the Unsplit MHD/Hydro units.
!!
!!***

#include "Flash.h"
#include "constants.h"

#define DIR_X 1
#define DIR_Y 2
#define DIR_Z 3


#define MINMOD  1
#define MC      2
#define HYBRID  3
#define VANLEER 4
#define LIMITED 5

#define HARTEN      1
#define HARTENHYMAN 2

#define UPDATE_NONE     0
#define UPDATE_INTERIOR 1
#define UPDATE_BOUND    2
#define UPDATE_ALL      3
#define UPDATE_SPECMS_INTERIOR 5
#define UPDATE_ALL_SPECMS_BOUND 6

#define FWDCONVERT 1
#define BWDCONVERT 2

#define WENO5 1
#define WENOZ 2

#define HY_VAR1_SCRATCHCTR_VAR 1
#define HY_VAR2_SCRATCHCTR_VAR 2
! The following six are currently only used for 3T MHD. It is not necessary to make
! the definitions conditional, this is just done here to catch cases of improper use.
#if defined(FLASH_USM_MHD) && defined(FLASH_UHD_3T)
#define HY_XN01_SCRATCHCTR_VAR 3
#define HY_XN02_SCRATCHCTR_VAR 4
#define HY_XN03_SCRATCHCTR_VAR 5
#define HY_XN04_SCRATCHCTR_VAR 6
#define HY_XN05_SCRATCHCTR_VAR 7
#define HY_XN06_SCRATCHCTR_VAR 8
#endif

!!--------------------------------------------------------------!!
!! [1] FOR UNSPLIT STAGGERED MHD IMPLEMENTATION ----------------!!
!!--------------------------------------------------------------!!
#ifdef FLASH_USM_MHD

#define UNDEFINED_PROL 0
#define INJECTION_PROL 1
#define BALSARA_PROL   2

#define ROE   1
#define LLF   2
#define HLL   3
#define HLLC  4
#define HLLD  5
#define MARQ  6
#define MARM  7
#define HYBR  8

#ifdef CURX_VAR
#define HY_NEED_EXTRA_GCFILL
#endif

#ifdef CURY_VAR
#define HY_NEED_EXTRA_GCFILL
#endif

#ifdef CURZ_VAR
#define HY_NEED_EXTRA_GCFILL
#endif

#ifdef DIVV_VAR
#define HY_NEED_EXTRA_GCFILL
#endif

!! PRIMITIVE VARIABLES FOR MHD
#define HY_DENS 1
#define HY_VELX 2
#define HY_VELY 3
#define HY_VELZ 4
#define HY_PRES 5
#define HY_MAGX 6
#define HY_MAGY 7
#define HY_MAGZ 8
#define HY_GAMC 9
#define HY_GAME 10
#define HY_EINT 11

! MHD 3T 
#ifdef FLASH_UHD_3T
#define HY_EELE 12
#define HY_EION 13
#define HY_ERAD 14
#define HY_GRAV 15
#define HY_TEMP 16
#define HY_ABAR 17
#define HY_ZBAR 18
#else
#define HY_GRAV 12
#define HY_TEMP 13
#define HY_ABAR 14
#define HY_ZBAR 15
#endif


!! CONSERVATIVE VARIABLES FOR MHD 
#define HY_XMOM 2
#define HY_YMOM 3
#define HY_ZMOM 4
#define HY_ENER 5
!! NOTE THAT THE REST OF CONSERVATIVE VARIABLES 
!! ARE ALREADY DEFINED IN PRIMITIVE VARIABLES

!! FLUX VARIABLES FOR MHD
#define HY_DENS_FLUX 1
#define HY_XMOM_FLUX 2
#define HY_YMOM_FLUX 3
#define HY_ZMOM_FLUX 4
#define HY_ENER_FLUX 5
#define HY_MAGX_FLUX 6
#define HY_MAGY_FLUX 7
#define HY_MAGZ_FLUX 8
#define HY_P_FLUX 9
#define HY_EINT_FLUX 10
#define HY_VOLU_FLUX 11

! MHD 3T
#ifdef FLASH_UHD_3T
#define HY_EELE_FLUX 12
#define HY_EION_FLUX 13
#define HY_ERAD_FLUX 14
#endif

!! WAVE STRUCTURES FOR MHD
#define HY_FASTLEFT 1
#define HY_ALFNLEFT 2
#define HY_SLOWLEFT 3
#define HY_ENTROPY  4
#define HY_SLOWRGHT 5
#define HY_ALFNRGHT 6
#define HY_FASTRGHT 7

!! EXTRA PARAMETERS FOR PURE MHD
#define HY_VARINUM HY_MAGZ

#ifdef FLASH_UHD_3T
#define HY_SCRATCH_NUM 14
#else
#define HY_SCRATCH_NUM 11
#endif

!! NOTE: the following scratch indices are the definitions of the LOCAL
!!       scratch arrays defined in hy_memScratchData.F90, NOT from Config file.
!! The following definitions are only applied to MHD_StaggeredMesh.
!! They replace global scratch array definitions maintained by the Grid unit.
!! These definitions will be used whether FLASH_UHD_NEED_SCRATCHVARS is
!! defined or not.

#ifndef FLASH_UHD_3T
!without 3T

#define HY_P01_FACEXPTR_VAR 1
#define HY_P02_FACEXPTR_VAR 2
#define HY_P03_FACEXPTR_VAR 3
#define HY_P04_FACEXPTR_VAR 4
#define HY_P05_FACEXPTR_VAR 5
#define HY_P06_FACEXPTR_VAR 6
#define HY_P07_FACEXPTR_VAR 7
#define HY_P08_FACEXPTR_VAR 8
#define HY_P09_FACEXPTR_VAR 9
#define HY_P10_FACEXPTR_VAR 10
#define HY_P11_FACEXPTR_VAR 11


#define HY_N01_FACEXPTR_VAR 12
#define HY_N02_FACEXPTR_VAR 13
#define HY_N03_FACEXPTR_VAR 14
#define HY_N04_FACEXPTR_VAR 15
#define HY_N05_FACEXPTR_VAR 16
#define HY_N06_FACEXPTR_VAR 17
#define HY_N07_FACEXPTR_VAR 18
#define HY_N08_FACEXPTR_VAR 19
#define HY_N09_FACEXPTR_VAR 20
#define HY_N10_FACEXPTR_VAR 21
#define HY_N11_FACEXPTR_VAR 22


#if NDIM >= 2
#define HY_P01_FACEYPTR_VAR 1
#define HY_P02_FACEYPTR_VAR 2
#define HY_P03_FACEYPTR_VAR 3
#define HY_P04_FACEYPTR_VAR 4
#define HY_P05_FACEYPTR_VAR 5
#define HY_P06_FACEYPTR_VAR 6
#define HY_P07_FACEYPTR_VAR 7
#define HY_P08_FACEYPTR_VAR 8
#define HY_P09_FACEYPTR_VAR 9
#define HY_P10_FACEYPTR_VAR 10
#define HY_P11_FACEYPTR_VAR 11


#define HY_N01_FACEYPTR_VAR 12
#define HY_N02_FACEYPTR_VAR 13
#define HY_N03_FACEYPTR_VAR 14
#define HY_N04_FACEYPTR_VAR 15
#define HY_N05_FACEYPTR_VAR 16
#define HY_N06_FACEYPTR_VAR 17
#define HY_N07_FACEYPTR_VAR 18
#define HY_N08_FACEYPTR_VAR 19
#define HY_N09_FACEYPTR_VAR 20
#define HY_N10_FACEYPTR_VAR 21
#define HY_N11_FACEYPTR_VAR 22


#if NDIM == 3
#define HY_P01_FACEZPTR_VAR 1
#define HY_P02_FACEZPTR_VAR 2
#define HY_P03_FACEZPTR_VAR 3
#define HY_P04_FACEZPTR_VAR 4
#define HY_P05_FACEZPTR_VAR 5
#define HY_P06_FACEZPTR_VAR 6
#define HY_P07_FACEZPTR_VAR 7
#define HY_P08_FACEZPTR_VAR 8
#define HY_P09_FACEZPTR_VAR 9
#define HY_P10_FACEZPTR_VAR 10
#define HY_P11_FACEZPTR_VAR 11


#define HY_N01_FACEZPTR_VAR 12
#define HY_N02_FACEZPTR_VAR 13
#define HY_N03_FACEZPTR_VAR 14
#define HY_N04_FACEZPTR_VAR 15
#define HY_N05_FACEZPTR_VAR 16
#define HY_N06_FACEZPTR_VAR 17
#define HY_N07_FACEZPTR_VAR 18
#define HY_N08_FACEZPTR_VAR 19
#define HY_N09_FACEZPTR_VAR 20
#define HY_N10_FACEZPTR_VAR 21
#define HY_N11_FACEZPTR_VAR 22

#endif
!end #if NDIM == 3
#endif
!end #if NDIM >= 2

#else
! with 3T
#define HY_P01_FACEXPTR_VAR 1
#define HY_P02_FACEXPTR_VAR 2
#define HY_P03_FACEXPTR_VAR 3
#define HY_P04_FACEXPTR_VAR 4
#define HY_P05_FACEXPTR_VAR 5
#define HY_P06_FACEXPTR_VAR 6
#define HY_P07_FACEXPTR_VAR 7
#define HY_P08_FACEXPTR_VAR 8
#define HY_P09_FACEXPTR_VAR 9
#define HY_P10_FACEXPTR_VAR 10
#define HY_P11_FACEXPTR_VAR 11
#define HY_P12_FACEXPTR_VAR 12
#define HY_P13_FACEXPTR_VAR 13
#define HY_P14_FACEXPTR_VAR 14

#define HY_N01_FACEXPTR_VAR 15
#define HY_N02_FACEXPTR_VAR 16
#define HY_N03_FACEXPTR_VAR 17
#define HY_N04_FACEXPTR_VAR 18
#define HY_N05_FACEXPTR_VAR 19
#define HY_N06_FACEXPTR_VAR 20
#define HY_N07_FACEXPTR_VAR 21
#define HY_N08_FACEXPTR_VAR 22
#define HY_N09_FACEXPTR_VAR 23
#define HY_N10_FACEXPTR_VAR 24
#define HY_N11_FACEXPTR_VAR 25
#define HY_N12_FACEXPTR_VAR 26
#define HY_N13_FACEXPTR_VAR 27
#define HY_N14_FACEXPTR_VAR 28


#if NDIM >= 2
#define HY_P01_FACEYPTR_VAR 1
#define HY_P02_FACEYPTR_VAR 2
#define HY_P03_FACEYPTR_VAR 3
#define HY_P04_FACEYPTR_VAR 4
#define HY_P05_FACEYPTR_VAR 5
#define HY_P06_FACEYPTR_VAR 6
#define HY_P07_FACEYPTR_VAR 7
#define HY_P08_FACEYPTR_VAR 8
#define HY_P09_FACEYPTR_VAR 9
#define HY_P10_FACEYPTR_VAR 10
#define HY_P11_FACEYPTR_VAR 11
#define HY_P12_FACEYPTR_VAR 12
#define HY_P13_FACEYPTR_VAR 13
#define HY_P14_FACEYPTR_VAR 14

#define HY_N01_FACEYPTR_VAR 15
#define HY_N02_FACEYPTR_VAR 16
#define HY_N03_FACEYPTR_VAR 17
#define HY_N04_FACEYPTR_VAR 18
#define HY_N05_FACEYPTR_VAR 19
#define HY_N06_FACEYPTR_VAR 20
#define HY_N07_FACEYPTR_VAR 21
#define HY_N08_FACEYPTR_VAR 22
#define HY_N09_FACEYPTR_VAR 23
#define HY_N10_FACEYPTR_VAR 24
#define HY_N11_FACEYPTR_VAR 25
#define HY_N12_FACEYPTR_VAR 26
#define HY_N13_FACEYPTR_VAR 27
#define HY_N14_FACEYPTR_VAR 28


#if NDIM == 3
#define HY_P01_FACEZPTR_VAR 1
#define HY_P02_FACEZPTR_VAR 2
#define HY_P03_FACEZPTR_VAR 3
#define HY_P04_FACEZPTR_VAR 4
#define HY_P05_FACEZPTR_VAR 5
#define HY_P06_FACEZPTR_VAR 6
#define HY_P07_FACEZPTR_VAR 7
#define HY_P08_FACEZPTR_VAR 8
#define HY_P09_FACEZPTR_VAR 9
#define HY_P10_FACEZPTR_VAR 10
#define HY_P11_FACEZPTR_VAR 11
#define HY_P12_FACEZPTR_VAR 12
#define HY_P13_FACEZPTR_VAR 13
#define HY_P14_FACEZPTR_VAR 14

#define HY_N01_FACEZPTR_VAR 15
#define HY_N02_FACEZPTR_VAR 16
#define HY_N03_FACEZPTR_VAR 17
#define HY_N04_FACEZPTR_VAR 18
#define HY_N05_FACEZPTR_VAR 19
#define HY_N06_FACEZPTR_VAR 20
#define HY_N07_FACEZPTR_VAR 21
#define HY_N08_FACEZPTR_VAR 22
#define HY_N09_FACEZPTR_VAR 23
#define HY_N10_FACEZPTR_VAR 24
#define HY_N11_FACEZPTR_VAR 25
#define HY_N12_FACEZPTR_VAR 26
#define HY_N13_FACEZPTR_VAR 27
#define HY_N14_FACEZPTR_VAR 28

#endif
!end #if NDIM == 3
#endif
!end #if NDIM >= 2

#endif
!end if ifndef FLASH_UHD_3T

#define HY_NSCRATCH_VARS (2*HY_SCRATCH_NUM)


!!--------------------------------------------------------------!!
!! [2] FOR UNSPLIT HYDRO IMPLEMENTATION ------------------------!!
!!--------------------------------------------------------------!!
#elif defined(FLASH_UHD_HYDRO)
#define ROE   1
#define LLF   2
#define HLL   3
#define HLLC  4
#define MARQ  5
#define MARM  7
#define HYBR  8

!! PRIMITIVE VARIABLES FOR PURE HYDRO
#define HY_DENS 1
#define HY_VELX 2
#define HY_VELY 3
#define HY_VELZ 4
#define HY_PRES 5
#define HY_GAMC 6
#define HY_GAME 7
#define HY_EINT 8

#ifdef FLASH_UHD_3T
#define HY_EELE 9
#define HY_EION 10
#define HY_ERAD 11
#define HY_GRAV 12
#define HY_TEMP 13
#define HY_ABAR 14
#define HY_ZBAR 15
#else
#define HY_GRAV 9
#define HY_TEMP 10
#define HY_ABAR 11
#define HY_ZBAR 12
#endif

!! CONSERVATIVE VARIABLES FOR PURE HYDRO
#define HY_XMOM 2
#define HY_YMOM 3
#define HY_ZMOM 4
#define HY_ENER 5

!! FLUX VARIABLES FOR PURE HYDRO
#define HY_DENS_FLUX 1
#define HY_XMOM_FLUX 2
#define HY_YMOM_FLUX 3
#define HY_ZMOM_FLUX 4
#define HY_ENER_FLUX 5
#define HY_P_FLUX 6
#define HY_EINT_FLUX 7
#define HY_VOLU_FLUX 8
#ifdef FLASH_UHD_3T
#define HY_EELE_FLUX 9
#define HY_EION_FLUX 10
#define HY_ERAD_FLUX 11
#endif

!! WAVE STRUCTURES FOR PURE HYDRO
#define HY_FASTLEFT 1
#define HY_SLOWLEFT 2
#define HY_ENTROPY  3
#define HY_SLOWRGHT 4
#define HY_FASTRGHT 5

!! EXTRA PARAMETERS FOR PURE HYDRO
#define HY_VARINUM HY_PRES

#ifdef FLASH_UHD_3T
#define HY_SCRATCH_NUM 11
#else
#define HY_SCRATCH_NUM 8
#endif

!! NOTE: the following scratch indices are the definitions of the LOCAL
!!       scratch arrays defined in hy_memScratchData.F90, NOT from Config file.
!! The following definitions are only applied to Hydro_Unsplit.
!! They replace global scratch array definitions maintained by the Grid unit.
!! These definitions will be used whether FLASH_UHD_NEED_SCRATCHVARS is defined
!! or not.

#ifndef FLASH_UHD_3T
!without 3T

#define HY_P01_FACEXPTR_VAR 1
#define HY_P02_FACEXPTR_VAR 2
#define HY_P03_FACEXPTR_VAR 3
#define HY_P04_FACEXPTR_VAR 4
#define HY_P05_FACEXPTR_VAR 5
#define HY_P06_FACEXPTR_VAR 6
#define HY_P07_FACEXPTR_VAR 7
#define HY_P08_FACEXPTR_VAR 8

#define HY_N01_FACEXPTR_VAR 9
#define HY_N02_FACEXPTR_VAR 10
#define HY_N03_FACEXPTR_VAR 11
#define HY_N04_FACEXPTR_VAR 12
#define HY_N05_FACEXPTR_VAR 13
#define HY_N06_FACEXPTR_VAR 14
#define HY_N07_FACEXPTR_VAR 15
#define HY_N08_FACEXPTR_VAR 16


#if NDIM >= 2
#define HY_P01_FACEYPTR_VAR 1
#define HY_P02_FACEYPTR_VAR 2
#define HY_P03_FACEYPTR_VAR 3
#define HY_P04_FACEYPTR_VAR 4
#define HY_P05_FACEYPTR_VAR 5
#define HY_P06_FACEYPTR_VAR 6
#define HY_P07_FACEYPTR_VAR 7
#define HY_P08_FACEYPTR_VAR 8

#define HY_N01_FACEYPTR_VAR 9
#define HY_N02_FACEYPTR_VAR 10
#define HY_N03_FACEYPTR_VAR 11
#define HY_N04_FACEYPTR_VAR 12
#define HY_N05_FACEYPTR_VAR 13
#define HY_N06_FACEYPTR_VAR 14
#define HY_N07_FACEYPTR_VAR 15
#define HY_N08_FACEYPTR_VAR 16


#if NDIM == 3
#define HY_P01_FACEZPTR_VAR 1
#define HY_P02_FACEZPTR_VAR 2
#define HY_P03_FACEZPTR_VAR 3
#define HY_P04_FACEZPTR_VAR 4
#define HY_P05_FACEZPTR_VAR 5
#define HY_P06_FACEZPTR_VAR 6
#define HY_P07_FACEZPTR_VAR 7
#define HY_P08_FACEZPTR_VAR 8

#define HY_N01_FACEZPTR_VAR 9
#define HY_N02_FACEZPTR_VAR 10
#define HY_N03_FACEZPTR_VAR 11
#define HY_N04_FACEZPTR_VAR 12
#define HY_N05_FACEZPTR_VAR 13
#define HY_N06_FACEZPTR_VAR 14
#define HY_N07_FACEZPTR_VAR 15
#define HY_N08_FACEZPTR_VAR 16

#endif
!end #if NDIM == 3
#endif
!end #if NDIM >= 2

#else
! with 3T
#define HY_P01_FACEXPTR_VAR 1
#define HY_P02_FACEXPTR_VAR 2
#define HY_P03_FACEXPTR_VAR 3
#define HY_P04_FACEXPTR_VAR 4
#define HY_P05_FACEXPTR_VAR 5
#define HY_P06_FACEXPTR_VAR 6
#define HY_P07_FACEXPTR_VAR 7
#define HY_P08_FACEXPTR_VAR 8
#define HY_P09_FACEXPTR_VAR 9
#define HY_P10_FACEXPTR_VAR 10
#define HY_P11_FACEXPTR_VAR 11

#define HY_N01_FACEXPTR_VAR 12
#define HY_N02_FACEXPTR_VAR 13
#define HY_N03_FACEXPTR_VAR 14
#define HY_N04_FACEXPTR_VAR 15
#define HY_N05_FACEXPTR_VAR 16
#define HY_N06_FACEXPTR_VAR 17
#define HY_N07_FACEXPTR_VAR 18
#define HY_N08_FACEXPTR_VAR 19
#define HY_N09_FACEXPTR_VAR 20
#define HY_N10_FACEXPTR_VAR 21
#define HY_N11_FACEXPTR_VAR 22


#if NDIM >= 2
#define HY_P01_FACEYPTR_VAR 1
#define HY_P02_FACEYPTR_VAR 2
#define HY_P03_FACEYPTR_VAR 3
#define HY_P04_FACEYPTR_VAR 4
#define HY_P05_FACEYPTR_VAR 5
#define HY_P06_FACEYPTR_VAR 6
#define HY_P07_FACEYPTR_VAR 7
#define HY_P08_FACEYPTR_VAR 8
#define HY_P09_FACEYPTR_VAR 9
#define HY_P10_FACEYPTR_VAR 10
#define HY_P11_FACEYPTR_VAR 11

#define HY_N01_FACEYPTR_VAR 12
#define HY_N02_FACEYPTR_VAR 13
#define HY_N03_FACEYPTR_VAR 14
#define HY_N04_FACEYPTR_VAR 15
#define HY_N05_FACEYPTR_VAR 16
#define HY_N06_FACEYPTR_VAR 17
#define HY_N07_FACEYPTR_VAR 18
#define HY_N08_FACEYPTR_VAR 19
#define HY_N09_FACEYPTR_VAR 20
#define HY_N10_FACEYPTR_VAR 21
#define HY_N11_FACEYPTR_VAR 22


#if NDIM == 3
#define HY_P01_FACEZPTR_VAR 1
#define HY_P02_FACEZPTR_VAR 2
#define HY_P03_FACEZPTR_VAR 3
#define HY_P04_FACEZPTR_VAR 4
#define HY_P05_FACEZPTR_VAR 5
#define HY_P06_FACEZPTR_VAR 6
#define HY_P07_FACEZPTR_VAR 7
#define HY_P08_FACEZPTR_VAR 8
#define HY_P09_FACEZPTR_VAR 9
#define HY_P10_FACEZPTR_VAR 10
#define HY_P11_FACEZPTR_VAR 11

#define HY_N01_FACEZPTR_VAR 12
#define HY_N02_FACEZPTR_VAR 13
#define HY_N03_FACEZPTR_VAR 14
#define HY_N04_FACEZPTR_VAR 15
#define HY_N05_FACEZPTR_VAR 16
#define HY_N06_FACEZPTR_VAR 17
#define HY_N07_FACEZPTR_VAR 18
#define HY_N08_FACEZPTR_VAR 19
#define HY_N09_FACEZPTR_VAR 20
#define HY_N10_FACEZPTR_VAR 21
#define HY_N11_FACEZPTR_VAR 22

#endif
!end #if NDIM == 3
#endif
!end #if NDIM >= 2

#endif
!end if ifndef FLASH_UHD_3T


#define HY_NSCRATCH_VARS (2*HY_SCRATCH_NUM)

#endif /* ifdef FLASH_UHD_HYDRO */

!!--------------------------------------------------------------!!
!! END OF MHD AND HYDRO DEFINITIONS ----------------------------!!
!!--------------------------------------------------------------!!


#ifndef FLASH_UHD_3T
#ifndef GRAVITY
#define HY_END_VARS HY_EINT
#else
#define HY_END_VARS HY_GRAV
#endif
#define HY_END_FLUX HY_VOLU_FLUX
#else
#ifndef GRAVITY
#define HY_END_VARS HY_ERAD
#else
#define HY_END_VARS HY_GRAV
#endif
#define HY_END_FLUX HY_ERAD_FLUX
#endif

#define HY_NSPEC NSPECIES+NMASS_SCALARS
#define HY_SPEC_BEG HY_END_VARS+1
#define HY_SPEC_END HY_END_VARS+NSPECIES+NMASS_SCALARS


!! DEFINE TOTAL NUMBERS OF VARIABLES AND WAVES
#define HY_WAVENUM  HY_FASTRGHT
#define HY_VARINUM1 (HY_VARINUM+1)
#define HY_VARINUM2 HY_VARINUM+2
#define HY_VARINUM3 HY_VARINUM+3
#define HY_VARINUM4 HY_VARINUM+4
#define HY_VARINUM7 HY_VARINUM+7

#ifndef FLASH_UHD_3T
#define HY_VARINUMMAX HY_VARINUM4
#else
#define HY_VARINUMMAX HY_VARINUM7
#endif
