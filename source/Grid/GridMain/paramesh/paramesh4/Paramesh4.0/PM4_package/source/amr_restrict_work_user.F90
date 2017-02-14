!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      subroutine amr_restrict_work_user()

!------------------------------------------------------------------------
!
! This is a stub routine for performing the interpolation operation
! during restriction of data stored in 'work'. It is meant to be replaced 
! by a user written routine if the provided default routine is insufficient
! for your application.
!
!------------------------------------------------------------------------

      return
      end subroutine amr_restrict_work_user
