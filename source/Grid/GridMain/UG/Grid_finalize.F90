!!****if* source/Grid/GridMain/UG/Grid_finalize
!!
!! NAME
!!
!!  Grid_finalize
!!
!!
!! SYNOPSIS
!!
!!  Grid_finalize()
!!
!!
!! DESCRIPTION
!!
!!  Deallocates memory that has been allocated in the Grid Unit
!!
!!***


subroutine Grid_finalize()

  use gr_bcInterface, ONLY : gr_bcFinalize
  use gr_ptInterface, ONLY : gr_ptFinalize
  use Grid_data, ONLY : gr_iCoords,gr_jCoords,gr_kCoords, gr_gid
  use physicalData, ONLY : unk,facevarx,facevary,facevarz
  use Grid_data, ONLY : scratch,scratch_ctr,&
       &scratch_facevarx,scratch_facevary,scratch_facevarz
  use gr_sbInterface, ONLY : gr_sbFinalize

  implicit none

#include "Flash.h"

  deallocate(gr_gid)

#ifndef FIXEDBLOCKSIZE  
  deallocate(gr_iCoords)
  deallocate(gr_jCoords)
  deallocate(gr_kCoords)
  deallocate(unk)
  deallocate(scratch)
  deallocate(facevarx)
  deallocate(facevary)
  deallocate(facevarz)
  deallocate(scratch_ctr)
  deallocate(scratch_facevarx)
  deallocate(scratch_facevary)
  deallocate(scratch_facevarz)

#endif
  call gr_solversFinalize()
  call gr_ptFinalize()
  call gr_bcFinalize()
  call gr_sbFinalize()
end subroutine Grid_finalize
