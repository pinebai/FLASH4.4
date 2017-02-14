!!****if* source/Grid/GridMain/paramesh/paramesh4/gr_commSetUp
!!
!! NAME
!!  gr_commSetUp
!!
!! SYNOPSIS
!!
!!  gr_commSetUp(integer, intent(IN) :: gridDataStruct)
!!
!!  
!! DESCRIPTION 
!!
!!  This subroutine is a wrapper around the Paramesh 4 procedure named 
!!  mpi_amr_comm_setup.  It will call the mpi_amr_comm_setup 
!!  subroutine in a way that is consistent with the Paramesh grid data 
!!  structure being used in the current simulation.  The grid data 
!!  structure type is described in the passed gridDataStruct argument.
!!
!!
!! ARGUMENTS   
!!
!!  gridDataStruct - integer constant, defined in "constants.h", 
!!                   indicating which grid data structure 
!!                   variable's guardcells to fill.
!!                   Paramesh has 5 data structures for grid variables, the first
!!                   four include all physical variables defined on the mesh. The 
!!                   fifth one includes a single variable.
!!
!!                   unk                cell centered, 
!!                   facex,facey,facez  face centered along i,j,k 
!!                                      direction respectively
!!                   work               cell centered, single variable.
!!                   
!!                   valid values of gridDataStruct are  
!!                   CENTER             unk only
!!                   WORK               work 
!!                   FACEX              facex
!!                   FACEY              facey
!!                   FACEZ              facez
!!                   FACES              facex,facey, and facez
!!                   CENTER_FACES     unk,facex,facey,facez
!!
!!***
subroutine gr_commSetUp(gridDataStruct)

#include "constants.h"

  use paramesh_mpi_interfaces, ONLY : mpi_amr_comm_setup
  use Grid_data, ONLY : gr_meshMe, gr_meshNumProcs
  implicit none
  integer, intent(IN) :: gridDataStruct
  logical :: lcc, lfc, lec, lnc, lguard, lprolong, lflux, ledge, &
             lrestrict, lfulltree
  integer :: tag_offset, iopt
  
  if((gridDataStruct==CENTER).or.(gridDataStruct==CENTER_FACES)) then
     lcc = .true.
  else
     lcc = .false.
  endif
  
  if((gridDataStruct==FACEX).or.(gridDataStruct==FACEY)&
       .or.(gridDataStruct==FACEZ).or.(gridDataStruct==CENTER_FACES)) then
     lfc = .true.
  else
     lfc = .false.
  end if
  
  lec       = .false.
  lnc       = .false.
  
  lguard    = .true.
  lprolong  = .false.
  lflux     = .false.
  ledge     = .false.
  lrestrict = .false.
  lfulltree = .false. ! It seems that FLASH never calls mpi_amr_comm_setup
                      ! (directly or indirectly) with fulltree true. - KW
  
  tag_offset = 100
  iopt = 1
  
  call mpi_amr_comm_setup(gr_meshMe,gr_meshNumProcs, &
       lguard,lprolong, &
       lflux,ledge,lrestrict, lfulltree, &
       iopt,lcc,lfc,lec,lnc,tag_offset)
  return
end subroutine gr_commSetUp
