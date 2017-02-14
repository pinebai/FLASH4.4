##
## Lines starting with ## are comments inside template file
## All other lines including empty lines are non-comments
## 
## This file is a template for Generating an F90 subroutine
## For syntax of this file see "Readme.template"
##
## VALID VARIABLE NAMES FOR THIS TEMPLATE
##
## values -> list of strings
## keys -> list of integers
## select_key_from_strlwr_and_map -> "select ... end select" F90 code
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! File created at setup time.  DO NOT EDIT !!!
!!

subroutine Simulation_mapStrToInt(str,key,map)
    use Grid_interface, only: Grid_parseNonRep
    use Driver_interface, only: Driver_getMype, Driver_getNumProcs
    use RuntimeParameters_interface, only: RuntimeParameters_get
    implicit none
#include "constants.h"
#include "Flash.h"
    character(len=*), intent(in) :: str
    integer, intent(out) :: key 
    integer, intent(in) :: map

    integer, parameter :: locunk1(0:NONREP_COUNT) = NONREP_LOCUNK1
    character(*), parameter :: rpcount_flat = NONREP_RPCOUNT_FLAT
    integer, parameter :: rpcount_start(NONREP_COUNT+1) = NONREP_RPCOUNT_START
    
    integer :: mesh, meshes
    character(len=MAX_STRING_LENGTH) :: strlwr
    integer :: nonrep, glob, nglob
    
    call Driver_getMype(MESH_ACROSS_COMM, mesh)
    call Driver_getNumProcs(MESH_ACROSS_COMM, meshes)
    key = NONEXISTENT
    strlwr = str
    call makeLowercase(strlwr)
    
    call Grid_parseNonRep(strlwr(1:len(str)), nonrep, glob)
    if(nonrep .gt. 0) then
        call RuntimeParameters_get(rpcount_flat(rpcount_start(nonrep):rpcount_start(nonrep+1)-1), nglob)
        if(glob .gt. nglob .or. mesh .ne. NONREP_MESHOFGLOB(glob, meshes)) return ! NONEXISTENT
        key = locunk1(nonrep)-1 + NONREP_GLOB2LOC(glob, mesh, meshes)
        return
    end if
    
    %(select_key_from_strlwr_and_map!\n    )s

    if(key .ne. NONEXISTENT) key = mod(key,MAPBLOCKSIZE)
end subroutine Simulation_mapStrToInt
