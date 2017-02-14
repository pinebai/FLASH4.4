!!****if* source/Multispecies/MultispeciesMain/ms_mapMSIndex
!!
!! NAME
!!
!!  ms_mapMSIndex
!!
!! SYNOPSIS
!!
!!  ms_mapMSIndex(integer, intent(in)  :: name,
!!                integer, intent(out)  :: msindex)
!!
!! DESCRIPTION
!!
!!  Maps the integer value of the name of the species to the correct
!!  index in the multispecies array.
!!  
!!  This is done under the hood as so not to make the user learn 2
!!  different naming schemes
!!
!! ARGUMENTS
!!
!!    name - name of species define in Flash.h ie NI56_SPEC
!!    msindex - the corresponding index in the multispecies array
!!
!!***



subroutine ms_mapMSIndex(name, msindex)

  implicit none

#include "Flash.h"

  integer, intent(in)   :: name
  integer, intent(out)  :: msindex

  msindex = name - NPROP_VARS

end subroutine ms_mapMSIndex
