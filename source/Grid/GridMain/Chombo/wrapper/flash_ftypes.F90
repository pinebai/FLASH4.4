!!****if* source/Grid/GridMain/Chombo/wrapper/flash_ftypes
!!
!! NAME
!!
!!  flash_ftypes
!!
!! SYNOPSIS
!!
!!  use flash_ftypes
!!
!! DESCRIPTION 
!!  
!!  These are the non-primitive datatypes used by FLASH to 
!!  interoperate with C++ functions that interact with Chombo library.
!!  Each datatype corresponds to a datatype in flash_ctypes.h
!!  
!! NOTES 
!!
!!  The interoperable types need to be placed in a separate
!!  module to work around a bug in gfortran:
!!  http://gcc.gnu.org/bugzilla/show_bug.cgi?id=40920
!!  http://groups.google.com/group/comp.lang.fortran/browse_thread/thread/220286db98888bb4#
!!  http://groups.google.com/group/comp.lang.fortran/browse_frm/thread/ddc211f2987326b8
!!   
!!***

#include "constants.h"
#include "flash_bool.h"

module flash_ftypes
  use iso_c_binding, only : c_double, c_int, c_ptr, c_char, c_bool
  implicit none

  type, bind(c) :: flash_ug_info_t
     real(c_double), dimension(MDIM) :: lowDomain
     real(c_double), dimension(MDIM) :: highDomain
     integer(c_int), dimension(MDIM) :: procGrid
     integer(c_int), dimension(MDIM) :: baseDomainSize
     integer(c_int), dimension(MDIM) :: guardCells
     integer(c_int), dimension(MDIM) :: domainBC
     integer(c_int), dimension(MAX_GRID_DATA_STRUCT_TMP) :: meshTypes
     integer(c_int), dimension(MAX_GRID_DATA_STRUCT_TMP) :: meshNumVars
     integer(c_int) :: verbosity
  end type flash_ug_info_t

  type, bind(c) :: flash_amr_info_t
     real(c_double), dimension(MDIM) :: lowDomain
     real(c_double), dimension(MDIM) :: highDomain
     real(c_double) :: BRMeshRefineFillRatio
     integer(c_int) :: BRMeshRefineBufferSize
     integer(c_int) :: BRMeshRefineBlockFactor
     integer(c_int), dimension(MDIM) :: maxBlockSize
     integer(c_int), dimension(MDIM) :: baseDomainSize
     integer(c_int), dimension(MDIM) :: guardCells
     integer(c_int), dimension(MDIM) :: domainBC
     integer(c_int), dimension(MAX_GRID_DATA_STRUCT_TMP) :: meshTypes
     integer(c_int), dimension(MAX_GRID_DATA_STRUCT_TMP) :: meshNumVars
     integer(c_int) :: maxRefineLevel
     integer(c_int) :: verbosity
     integer(c_int) :: quadCFInterp
     integer(c_int) :: fluxCorrect
     integer(c_int) :: refRatio;
     integer(c_int) :: restrictBeforeGhostExchange;
     integer(c_int) :: scaleFineFluxes;
  end type flash_amr_info_t

  type, bind(c) :: box_info_t
     real(c_double), dimension(MDIM) :: deltas     !Grid_getDeltas
     real(c_double), dimension(MDIM) :: bsize      !Grid_getBlkPhysicalSize
     real(c_double), dimension(MDIM) :: coord      !Grid_getBlkCenterCoords
     real(c_double), dimension(MDIM) :: lowBndbox  !Grid_getBlkBoundBox
     real(c_double), dimension(MDIM) :: highBndbox !Grid_getBlkBoundBox
     integer(c_int), dimension(MDIM) :: lowLimits  !Grid_getBlkIndexLimits
     integer(c_int), dimension(MDIM) :: highLimits !Grid_getBlkIndexLimits
     integer(c_int), dimension(MDIM) :: guardcells !Grid_getBlkIndexLimits
     integer(c_int), dimension(MDIM) :: corner     !Grid_getBlkCornerID
     integer(c_int), dimension(MDIM) :: stride     !Grid_getBlkCornerID
     integer(c_int) :: lrefine                     !Grid_getBlkRefineLevel
     integer(c_int), dimension(MDIM) :: isNextToLowDomain
     integer(c_int), dimension(MDIM) :: isNextToHighDomain
  end type box_info_t
  
    type, bind(c) :: named_vals_t
        integer(c_int) :: real_count
        type(c_ptr) :: real_names ! character(kind=c_char,len=MAX_STRING_LENGTH), pointer, dimension(:)
        type(c_ptr) :: real_vals ! real(c_double), pointer, dimension(:)
        
        integer(c_int) :: int_count
        type(c_ptr) :: int_names ! character(kind=c_char,len=MAX_STRING_LENGTH), pointer, dimension(:)
        type(c_ptr) :: int_vals ! integer(c_int), pointer, dimension(:)
        
        integer(c_int) :: str_count
        type(c_ptr) :: str_names ! character(kind=c_char,len=MAX_STRING_LENGTH), pointer, dimension(:)
        type(c_ptr) :: str_vals ! character(kind=c_char,len=MAX_STRING_LENGTH), pointer, dimension(:)
        
        integer(c_int) :: log_count
        type(c_ptr) :: log_names ! character(kind=c_char,len=MAX_STRING_LENGTH), pointer, dimension(:)
        type(c_ptr) :: log_vals ! integer(c_int), pointer, dimension(:)
    end type
end module flash_ftypes
