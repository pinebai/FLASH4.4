!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* mpi_source/mpi_amr_global_domain_limits
!! NAME
!!
!!   mpi_amr_global_domain_limits
!!
!! SYNOPSIS
!!
!!   call mpi_amr_global_domain_limits ()
!!
!! ARGUMENTS
!!
!!   No arguments.
!!  
!! INCLUDES
!!
!!   paramesh_preprocessor.fh
!!   mpif.h
!!
!! USES
!!
!!   paramesh_dimensions
!!   physicaldata
!!   tree
!!   paramesh_comm_data
!!
!! CALLS
!!
!!   No other Paramesh routines are called.
!!
!! RETURNS
!!
!!   Upon return the coordinate ranges for the computational domain are found
!!   and are placed in the variables 'grid_xmin', 'grid_xmax', 'grid_ymin', 
!!   'grid_ymax', 'grid_zmin' and 'grid_zmax'.
!!
!! DESCRIPTION
!! 
!!   This routine computes the coordinate ranges for the computational domain
!!   as it currently exists.
!! 
!! AUTHORS
!!
!!   Written by Peter MacNeice.
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

         subroutine mpi_amr_global_domain_limits

         use paramesh_dimensions
         use physicaldata
         use tree
         Use paramesh_comm_data

         implicit none

         include 'mpif.h'

         real :: xmin,ymin,zmin,xmax,ymax,zmax
         real :: xmin1,ymin1,zmin1,xmax1,ymax1,zmax1
         integer :: ierr


! Find the coordinate ranges
         xmin1 = minval(bnd_box(1,1,1:lnblocks))
         ymin1 = minval(bnd_box(1,2,1:lnblocks))
         zmin1 = minval(bnd_box(1,3,1:lnblocks))
         xmax1 = maxval(bnd_box(2,1,1:lnblocks))
         ymax1 = maxval(bnd_box(2,2,1:lnblocks))
         zmax1 = maxval(bnd_box(2,3,1:lnblocks))
         xmin = min(1.e30, xmin1)
         ymin = min(1.e30, ymin1)
         zmin = min(1.e30, zmin1)
         xmax = max(-1.e30, xmax1)
         ymax = max(-1.e30, ymax1)
         zmax = max(-1.e30, zmax1)
      call MPI_ALLREDUCE (xmin,grid_xmin,1, & 
     &                    amr_mpi_real, & 
     &                    MPI_MIN,amr_mpi_meshComm,ierr)
      call MPI_ALLREDUCE (ymin,grid_ymin,1, & 
     &                    amr_mpi_real, & 
     &                    MPI_MIN,amr_mpi_meshComm,ierr)
      call MPI_ALLREDUCE (zmin,grid_zmin,1, & 
     &                    amr_mpi_real, & 
     &                    MPI_MIN,amr_mpi_meshComm,ierr)
      call MPI_ALLREDUCE (xmax,grid_xmax,1, & 
     &                    amr_mpi_real, & 
     &                    MPI_MAX,amr_mpi_meshComm,ierr)
      call MPI_ALLREDUCE (ymax,grid_ymax,1, & 
     &                    amr_mpi_real, & 
     &                    MPI_MAX,amr_mpi_meshComm,ierr)
      call MPI_ALLREDUCE (zmax,grid_zmax,1, & 
     &                    amr_mpi_real, & 
     &                    MPI_MAX,amr_mpi_meshComm,ierr)

      return
      end subroutine mpi_amr_global_domain_limits
