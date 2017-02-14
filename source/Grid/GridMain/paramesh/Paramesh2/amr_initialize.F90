      subroutine amr_initialize


! $RCSfile: amr_initialize.F90,v $
! $Revision: 1.1 $
! $Date: 2004/05/20 01:03:44 $


! This subroutine initializes the amr package. It performs any
! initialization required by the package which is application
! independent.

!
! NOTE : This routine MUST BE the first executed code in your application!!!!
!

use physicaldata
      use tree
      implicit none
      include 'mpif.h'



      integer ierr,mype,nprocs,maxprocs

!----------------------------------------------------------------

! Call the machine/software environment specific initialization routine.
! Different versions of comm_start are provided for use with machines
! which run or mpi. Make sure to compile with the appropriate version
! for your environment.

!      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,mype,ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
      MaxProcs = nprocs
      lnblocks = 0

! initialize tree data structure
        bsize(:,:) = -1.
        lrefine(:) = -1
        nodetype(:) = -1
        refine(:) = .FALSE.
        derefine(:) = .FALSE.
        parent(:,:) = -1
        child(:,:,:) = -1
        coord(:,:) = -1.
        bnd_box(:,:,:) = -1.
        neigh(:,:,:) = -1
        empty(:) = 0

! initialize solution arrays
        unk(:,:,:,:,:) = 0.
	unk2(:,:,:)    = 0.
	unk3(:,:,:)    = 0.
	unk4(:,:,:)    = 0.
        facevarx(:,:,:,:,:) = 0.
        facevary(:,:,:,:,:) = 0.
        facevarz(:,:,:,:,:) = 0.

! initialize other arrays
        flux_x(:,:,:,:,:) = -9999.9
        flux_y(:,:,:,:,:) = -9999.9
        flux_z(:,:,:,:,:) = -9999.9
!       .....


! Initialization required for prolongation routines
      call amr_prolong_fun_init

      return
      end



