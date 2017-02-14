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

!----------------------------------------------------------------------
! This is a stub routine to hold the place of a user written 
! interpolation routine to used during restriction of cell centered data 
! from fine to course meshes.
!----------------------------------------------------------------------

      subroutine amr_restrict_unk_user()

      return
      end subroutine amr_restrict_unk_user
