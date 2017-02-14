##
## Lines starting with ## are comments inside template file
## All other lines including empty lines are non-comments
## 
## This file is a template for generating an F90 subroutine
## For syntax of this file see "Readme.template"
##
## VALID VARIABLE NAMES FOR THIS TEMPLATE
##
## len_eos_lists        -> length of list of Eos roles, see Eos_map.h
## eos_unk              -> the map for unk
## eos_face             -> the map for face
## eos_scratch          -> the map for scratch
## eos_scratch_ctr      -> the map for scratch_ctr
## eos_scratch_facexvar -> the map for scratch_facevarx
## eos_scratch_faceyvar -> the map for scratch_facevary
## eos_scratch_facezvar -> the map for scratch_facevarz
##
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! File created at setup time.  DO NOT EDIT !!!

#include "constants.h"
#include "Flash.h"
#include "Eos.h"
#include "Eos_map.h"

#define EOSIN 0
#define EOSOUT 1

integer function eos_variableMap(gridDataStruct, eosRole, direction)
   use Driver_interface, ONLY : Driver_abortFlash
   implicit none
   integer, intent(IN) :: gridDataStruct, eosRole
   integer, intent(IN) :: direction ! 0 for IN, 1 for OUT

   !The size of eos_unk should be the same as EOSMAP_NUM_ROLES preprocessor symbol.
   integer, save, dimension(1:EOSMAP_NUM_ROLES,0:1) :: eosmap_unk, eosmap_face, &
        eosmap_scratch, eosmap_scratch_ctr, eosmap_scratch_facexvar, &
        eosmap_scratch_faceyvar, eosmap_scratch_facezvar
   integer :: templateEosSize, eosIndex
   logical, save :: firstCall = .true.

   eosIndex = NONEXISTENT

   !Check input arguments are sensible.
   if (eosRole < 1 .or. eosRole > EOSMAP_NUM_ROLES) then
      call Driver_abortFlash("[eos_variableMap]: EOS index out of bounds")
   end if
   if ((direction /= EOSIN) .and. (direction /= EOSOUT)) then
      call Driver_abortFlash ("[eos_variableMap]: " //&
           "Direction is neither 0 (for IN) nor 1 (for OUT).")
   end if


   if (firstcall) then
      firstCall = .false.

      templateEosSize = %(len_eos_lists|1)s
      if (templateEosSize /= EOSMAP_NUM_ROLES) then
         call Driver_abortFlash &
              ("[eos_variableMap]: Eos_map.h inconsistent with setup script")
      end if

      !This is the map for unk:
      %(eos_unk!\n      )s

      !This is the map for face:
      %(eos_face!\n      )s

      !This is the map for scratch:
      %(eos_scratch!\n      )s

      !This is the map for scratch_ctr:
      %(eos_scratch_ctr!\n      )s

      !This is the map for scratch_facevarx:
      %(eos_scratch_facexvar!\n      )s

      !This is the map for scratch_facevary:
      %(eos_scratch_faceyvar!\n      )s

      !This is the map for scratch_facevarz:
      %(eos_scratch_facezvar!\n      )s
   end if


   !Return the appropriate element from the appropriate data structure.
   select case (gridDataStruct)
   case (CENTER)
      eosIndex = eosmap_unk(eosRole,direction)
   case (FACEX,FACEY,FACEZ)
      eosIndex = eosmap_face(eosRole,direction)
   case (SCRATCH)
      eosIndex = eosmap_scratch(eosRole,direction)
   case (SCRATCH_CTR)
      eosIndex = eosmap_scratch_ctr(eosRole,direction)
   case (SCRATCH_FACEX)
      eosIndex = eosmap_scratch_facexvar(eosRole,direction)
   case (SCRATCH_FACEY)
      eosIndex = eosmap_scratch_faceyvar(eosRole,direction)
   case (SCRATCH_FACEZ)
      eosIndex = eosmap_scratch_facezvar(eosRole,direction)
#ifdef FLASH_GRID_PARAMESH
   case (WORK)
      call Driver_abortFlash &
        ("[eos_variableMap]: work structure not implemented")
#endif
   case default
      call Driver_abortFlash &
         ("[eos_variableMap]: gridDataStruct not recognised")
   end select
   eos_variableMap = eosIndex

end function eos_variableMap
