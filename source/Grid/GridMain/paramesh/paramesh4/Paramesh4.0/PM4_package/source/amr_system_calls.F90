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
!
! Some useful system calls which may or may not be available on some systems.
!
!
! $RCSfile: amr_system_calls.F90,v $
! $Revision: 1.2 $
! $Date: 2007/01/23 16:44:48 $
!
!--------------------------------------------------------------
!
! abort
! This is an SGI IRIX or UNICOS command. use amr_abort instead.

!
! flush
! use amr_flush instead

!--------------------------------------------------------------




      subroutine amr_abort()
      use paramesh_comm_data
      implicit none
      include 'mpif.h'
      integer :: ierrorcode,ierr

      call mpi_abort(amr_mpi_meshComm,ierrorcode,ierr)

      return
      end subroutine amr_abort


      subroutine amr_flush(iunit)
      implicit none
      integer,intent(in) :: iunit

      return
      end subroutine amr_flush

