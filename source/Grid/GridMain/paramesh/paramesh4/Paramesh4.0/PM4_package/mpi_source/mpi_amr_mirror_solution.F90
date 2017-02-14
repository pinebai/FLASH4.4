!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* mpi_source/amr_mirror_solution
!! NAME
!!
!!   amr_mirror_solution
!!
!! SYNOPSIS
!!
!!   call amr_mirror_solution(mype, nprocs, mirror_blks, ix, iy, iz)
!!
!!   call amr_mirror_solution(integer, integer, integer, integer, integer, integer)
!!
!! ARGUMENTS
!!   
!!   integer, intent(in) :: mype  
!!     The calling processor.
!!
!!   integer, intent(in) :: nprocs
!!     The number for processors used.
!!
!!   integer, intent(inout) :: mirror_blks
!!     The list of mirror blocks found by amr_mirror_blks
!!
!!   integer, intent(in) :: ix, iy, iz
!!     ???
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
!!   mpi_morton
!!   paramesh_mpi_interfaces
!!   paramesh_comm_data
!!
!! CALLS
!! 
!!   morton_number
!!    
!! RETURNS
!!
!!  Nothing returned.
!!
!! DESCRIPTION
!!
!!  This routine adjusts the solution in one half of a symmetric
!!  grid. The grid must be symmetric for reflection in a plane 
!!  which bisects one of the coordinate directions. The symmetry
!!  plane is defined by the integers ix,iy,iz. Two of these integers
!!  must be zero, and the other either 1 or 2. If 1, then the low
!!  coordinate half of the grid sends data to the high coordinate
!!  half - if 2, the data is sent in the opposite direction.
!!
!! AUTHORS
!!
!!   Peter MacNeice (December 2003).
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

!#define DEBUG

      subroutine amr_mirror_solution(mype,nprocs,mirror_blks,ix,iy,iz)

      use paramesh_dimensions
      use physicaldata
      use tree
      use mpi_morton
      Use paramesh_comm_data

      use paramesh_mpi_interfaces, only : morton_number

      implicit none
      include 'mpif.h'


      integer, intent(in)    ::  mype,nprocs,ix,iy,iz
      integer, intent(in)    ::  mirror_blks(:,:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! local variables

      real    :: xmin,xmax,ymin,ymax,zmin,zmax
      real    :: xmid,ymid,zmid
      real    :: xtest,ytest,ztest,direction

      real    :: ttunk(nvar,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd)
      real    :: ttfacevarx(nbndvar,il_bnd:iu_bnd+1,jl_bnd:ju_bnd, & 
     &                     kl_bnd:ku_bnd)
      real    :: ttfacevary(nbndvar,il_bnd:iu_bnd,jl_bnd:ju_bnd+k2d, & 
     &                     kl_bnd:ku_bnd)
      real    :: ttfacevarz(nbndvar,il_bnd:iu_bnd,jl_bnd:ju_bnd, & 
     &                     kl_bnd:ku_bnd+k3d)
      real    :: ttunk_e_x(nbndvare,il_bnd:iu_bnd,jl_bnd:ju_bnd+k2d, & 
     &                     kl_bnd:ku_bnd+k3d)
      real    :: ttunk_e_y(nbndvare,il_bnd:iu_bnd+1,jl_bnd:ju_bnd, & 
     &                     kl_bnd:ku_bnd+k3d)
      real    :: ttunk_e_z(nbndvare,il_bnd:iu_bnd+1,jl_bnd:ju_bnd+k2d, & 
     &                     kl_bnd:ku_bnd)
      real    :: ttunk_n(nbndvarc,il_bnd:iu_bnd+1,jl_bnd:ju_bnd+k2d, & 
     &                     kl_bnd:ku_bnd+k3d)
      real    ::  sign_cc_vars(nvar)

      integer :: lb,i,j,k
      integer :: ii,jj,kk
      integer :: iprocs,isize,isrc,idest,itag
      integer :: ierr,no_of_recvs
      integer :: allocation_status

      logical :: send_half(maxblocks)

      integer,dimension (:),  allocatable :: recvrequest
      integer,dimension (:,:),allocatable :: recvstatus


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Set the sign change required during the symmetry mapping
! if any of the cell-centered quantities are components of vectors.
      sign_cc_vars = 1.

#ifdef EXAMPLE
      if(ix.gt.0) then
        sign_cc_vars(2) = -1.
      elseif(iy.gt.0) then
        sign_cc_vars(3) = -1.
      elseif(iz.gt.0) then
        sign_cc_vars(4) = -1.
      endif
#endif /* EXAMPLE */

!
! get xmin,ymin,zmin,xmax,ymax,zmax from storage
      xmin = grid_xmin
      ymin = grid_ymin
      zmin = grid_zmin
      xmax = grid_xmax
      ymax = grid_ymax
      zmax = grid_zmax
      xmid = .5*(grid_xmin + grid_xmax)
      ymid = .5*(grid_ymin + grid_ymax)
      zmid = .5*(grid_zmin + grid_zmax)

!--------------------------------------------------
!
#ifdef DEBUG
      write(*,*) 'xmin,ymin,zmin,xmax,ymax,zmax ', & 
     & xmin,ymin,zmin,xmax,ymax,zmax
#endif /* DEBUG */

      no_of_recvs = 0
      do lb=1,lnblocks
      if(nodetype(lb).eq.1) then

! Is block lb in the half of the grid which should send data?
! set logical send_half
        if(ix.eq.0) then
          xtest = -1.
        elseif(ix.eq.1) then
          xtest = (coord(1,lb)-grid_xmin)*(coord(1,lb)-xmid)
        elseif(ix.eq.2) then
          xtest = (coord(1,lb)-grid_xmax)*(coord(1,lb)-xmid)
        endif
        if(ndim.ge.2) then
        if(iy.eq.0) then
          ytest = -1.
        elseif(iy.eq.1) then
          ytest = (coord(2,lb)-grid_ymin)*(coord(2,lb)-ymid)
        elseif(iy.eq.2) then
          ytest = (coord(2,lb)-grid_ymax)*(coord(2,lb)-ymid)
        endif
        else
          ytest = -1.
        endif
        if(ndim.eq.3) then
        if(iz.eq.0) then
          ztest = -1.
        elseif(iz.eq.1) then
          ztest = (coord(3,lb)-grid_zmin)*(coord(3,lb)-zmid)
        elseif(iz.eq.2) then
          ztest = (coord(3,lb)-grid_zmax)*(coord(3,lb)-zmid)
        endif
        else
          ztest = -1.
        endif

        send_half(lb) = .false.
        if(xtest.lt.0..and.ytest.lt.0..and.ztest.lt.0.) & 
     &                          send_half(lb) = .true.

        if(.not.send_half(lb)) no_of_recvs = no_of_recvs+1
      endif
      enddo


      if(allocated(recvrequest)) deallocate( recvrequest )
      allocate ( recvrequest(no_of_recvs) )

      if(allocated(recvstatus)) deallocate( recvstatus )
      allocate ( recvstatus(MPI_STATUS_SIZE,no_of_recvs) )

!--------------------------------------------------
!
! Set size of messages for unk data
      isize = len_block
!
! All leaf blocks in receiving half post a non-blocking receive

      no_of_recvs = 0
      do lb = 1,lnblocks
      if(nodetype(lb).eq.1.and..not.send_half(lb)) then

        isrc = mirror_blks(2,lb)
#ifdef PM_UNIQUE_MPI_TAGS
        itag = mype*maxblocks+lb
#else
        itag = lb
#endif
        no_of_recvs = no_of_recvs+1
        call Mpi_Irecv(unk(1,1,1,1,lb),isize, & 
     &                 amr_mpi_real, & 
     &                 isrc ,itag,amr_mpi_meshComm, & 
     &                 recvrequest(no_of_recvs),ierr)
#ifdef DEBUG
      write(*,*) 'pe ',mype,' receive of size ',isize, & 
     &          ' from ',isrc,' tag ',itag,' posted'
#endif /* DEBUG  */
      endif
      enddo
!
! All leaf blocks in sending half post a blocking send
      do lb = 1,lnblocks
      if(nodetype(lb).eq.1.and.send_half(lb)) then
        idest = mirror_blks(2,lb)
#ifdef PM_UNIQUE_MPI_TAGS
        itag = idest*maxblocks+mirror_blks(1,lb)
#else
        itag = mirror_blks(1,lb)
#endif
        call MPI_Ssend(unk(1,1,1,1,lb),isize, & 
     &                 amr_mpi_real, & 
     &                 idest,itag,amr_mpi_meshComm,ierr)
#ifdef DEBUG
      write(*,*) 'pe ',mype,' unk of size ',isize, & 
     &          ' sent to ',idest,' tag ',itag
#endif /* DEBUG  */
      endif
      enddo

      if(no_of_recvs.gt.0)  & 
     &     call MPI_Waitall(no_of_recvs,recvrequest,recvstatus,ierr)

!
!  end of unk messages
!--------------------------------------------------
      if(nfacevar.gt.0) then

!
! Set size of messages for facevarx data
      isize = nbndvar*(iu_bnd+1)*ju_bnd*ku_bnd
!
! All leaf blocks in receiving half post a non-blocking receive

      no_of_recvs = 0
      do lb = 1,lnblocks
      if(nodetype(lb).eq.1.and..not.send_half(lb)) then

        isrc = mirror_blks(2,lb)
#ifdef PM_UNIQUE_MPI_TAGS
        itag = mype*maxblocks+lb
#else
        itag = lb
#endif
        no_of_recvs = no_of_recvs+1
        call Mpi_Irecv(facevarx(1,1,1,1,lb),isize, & 
     &                 amr_mpi_real, & 
     &                 isrc ,itag,amr_mpi_meshComm, & 
     &                 recvrequest(no_of_recvs),ierr)
#ifdef DEBUG
      write(*,*) 'pe ',mype,' receive of size ',isize, & 
     &          ' from ',isrc,' tag ',itag,' posted'
#endif /* DEBUG  */
      endif
      enddo
!
! All leaf blocks in sending half post a blocking send
      do lb = 1,lnblocks
      if(nodetype(lb).eq.1.and.send_half(lb)) then
        idest = mirror_blks(2,lb)
#ifdef PM_UNIQUE_MPI_TAGS
        itag = idest*maxblocks+mirror_blks(1,lb)
#else
        itag = mirror_blks(1,lb)
#endif
        call MPI_Ssend(facevarx(1,1,1,1,lb),isize, & 
     &                 amr_mpi_real, & 
     &                 idest,itag,amr_mpi_meshComm,ierr)
#ifdef DEBUG
      write(*,*) 'pe ',mype,' facevarx of size ',isize, & 
     &          ' sent to ',idest,' tag ',itag
#endif /* DEBUG  */
      endif
      enddo

      if(no_of_recvs.gt.0) & 
     &     call MPI_Waitall(no_of_recvs,recvrequest,recvstatus,ierr)

      endif
!  end of facevarx messages
!--------------------------------------------------
      if(nfacevar.gt.0.and.ndim.ge.2) then

!
! Set size of messages for facevary data
      isize = nbndvar*iu_bnd*(ju_bnd+k2d)*ku_bnd
!
! All leaf blocks in receiving half post a non-blocking receive

      no_of_recvs = 0
      do lb = 1,lnblocks
      if(nodetype(lb).eq.1.and..not.send_half(lb)) then

        isrc = mirror_blks(2,lb)
#ifdef PM_UNIQUE_MPI_TAGS
        itag = mype*maxblocks+lb
#else
        itag = lb
#endif
        no_of_recvs = no_of_recvs+1
        call Mpi_Irecv(facevary(1,1,1,1,lb),isize, & 
     &                 amr_mpi_real, & 
     &                 isrc ,itag,amr_mpi_meshComm, & 
     &                 recvrequest(no_of_recvs),ierr)
#ifdef DEBUG
      write(*,*) 'pe ',mype,' receive of size ',isize, & 
     &          ' from ',isrc,' tag ',itag,' posted'
#endif /* DEBUG  */
      endif
      enddo
!
! All leaf blocks in sending half post a blocking send
      do lb = 1,lnblocks
      if(nodetype(lb).eq.1.and.send_half(lb)) then
        idest = mirror_blks(2,lb)
#ifdef PM_UNIQUE_MPI_TAGS
        itag = idest*maxblocks+mirror_blks(1,lb)
#else
        itag = mirror_blks(1,lb)
#endif
        call MPI_Ssend(facevary(1,1,1,1,lb),isize, & 
     &                 amr_mpi_real, & 
     &                 idest,itag,amr_mpi_meshComm,ierr)
#ifdef DEBUG
      write(*,*) 'pe ',mype,' facevary of size ',isize, & 
     &          ' sent to ',idest,' tag ',itag
#endif /* DEBUG  */
      endif
      enddo

      if(no_of_recvs.gt.0) & 
     &     call MPI_Waitall(no_of_recvs,recvrequest,recvstatus,ierr)

      endif
!  end of facevary messages
!--------------------------------------------------
      if(nfacevar.gt.0.and.ndim.eq.3) then

!
! Set size of messages for facevarz data
      isize = nbndvar*iu_bnd*ju_bnd*(ku_bnd+k3d)
!
! All leaf blocks in receiving half post a non-blocking receive

      no_of_recvs = 0
      do lb = 1,lnblocks
      if(nodetype(lb).eq.1.and..not.send_half(lb)) then

        isrc = mirror_blks(2,lb)
#ifdef PM_UNIQUE_MPI_TAGS
        itag = mype*maxblocks+lb
#else
        itag = lb
#endif
        no_of_recvs = no_of_recvs+1
        call Mpi_Irecv(facevarz(1,1,1,1,lb),isize, & 
     &                 amr_mpi_real, & 
     &                 isrc ,itag,amr_mpi_meshComm, & 
     &                 recvrequest(no_of_recvs),ierr)
#ifdef DEBUG
      write(*,*) 'pe ',mype,' receive of size ',isize, & 
     &          ' from ',isrc,' tag ',itag,' posted'
#endif /* DEBUG  */
      endif
      enddo
!
! All leaf blocks in sending half post a blocking send
      do lb = 1,lnblocks
      if(nodetype(lb).eq.1.and.send_half(lb)) then
        idest = mirror_blks(2,lb)
#ifdef PM_UNIQUE_MPI_TAGS
        itag = idest*maxblocks+mirror_blks(1,lb)
#else
        itag = mirror_blks(1,lb)
#endif
        call MPI_Ssend(facevarz(1,1,1,1,lb),isize, & 
     &                 amr_mpi_real, & 
     &                 idest,itag,amr_mpi_meshComm,ierr)
#ifdef DEBUG
      write(*,*) 'pe ',mype,' facevarz of size ',isize, & 
     &          ' sent to ',idest,' tag ',itag
#endif /* DEBUG  */
      endif
      enddo

      if(no_of_recvs.gt.0) & 
     &     call MPI_Waitall(no_of_recvs,recvrequest,recvstatus,ierr)

      endif
!  end of facevarz messages
!--------------------------------------------------
      if(nvaredge.gt.0.and.ndim.ge.2) then
!
! Set size of messages for unk_e_x data
      isize = nbndvare*iu_bnd*(ju_bnd+k2d)*(ku_bnd+k3d)
!
! All leaf blocks in receiving half post a non-blocking receive

      no_of_recvs = 0
      do lb = 1,lnblocks
      if(nodetype(lb).eq.1.and..not.send_half(lb)) then

        isrc = mirror_blks(2,lb)
#ifdef PM_UNIQUE_MPI_TAGS
        itag = mype*maxblocks+lb
#else
        itag = lb
#endif
        no_of_recvs = no_of_recvs+1
        call Mpi_Irecv(unk_e_x(1,1,1,1,lb),isize, & 
     &                 amr_mpi_real, & 
     &                 isrc ,itag,amr_mpi_meshComm, & 
     &                 recvrequest(no_of_recvs),ierr)
#ifdef DEBUG
      write(*,*) 'pe ',mype,' receive of size ',isize, & 
     &          ' from ',isrc,' tag ',itag,' posted'
#endif /* DEBUG  */
      endif
      enddo
!
! All leaf blocks in sending half post a blocking send
      do lb = 1,lnblocks
      if(nodetype(lb).eq.1.and.send_half(lb)) then
        idest = mirror_blks(2,lb)
#ifdef PM_UNIQUE_MPI_TAGS
        itag = idest*maxblocks+mirror_blks(1,lb)
#else
        itag = mirror_blks(1,lb)
#endif
        call MPI_Ssend(unk_e_x(1,1,1,1,lb),isize, & 
     &                 amr_mpi_real, & 
     &                 idest,itag,amr_mpi_meshComm,ierr)
#ifdef DEBUG
      write(*,*) 'pe ',mype,' unk_e_x of size ',isize, & 
     &          ' sent to ',idest,' tag ',itag
#endif /* DEBUG  */
      endif
      enddo

      if(no_of_recvs.gt.0) & 
     &     call MPI_Waitall(no_of_recvs,recvrequest,recvstatus,ierr)

      endif
!  end of unk_e_x messages
!--------------------------------------------------
      if(nvaredge.gt.0.and.ndim.ge.2) then
!
! Set size of messages for unk_e_y data
      isize = nbndvare*(iu_bnd+1)*ju_bnd*(ku_bnd+k3d)
!
! All leaf blocks in receiving half post a non-blocking receive

      no_of_recvs = 0
      do lb = 1,lnblocks
      if(nodetype(lb).eq.1.and..not.send_half(lb)) then

        isrc = mirror_blks(2,lb)
#ifdef PM_UNIQUE_MPI_TAGS
        itag = mype*maxblocks+lb
#else
        itag = lb
#endif
        no_of_recvs = no_of_recvs+1
        call Mpi_Irecv(unk_e_y(1,1,1,1,lb),isize, & 
     &                 amr_mpi_real, & 
     &                 isrc ,itag,amr_mpi_meshComm, & 
     &                 recvrequest(no_of_recvs),ierr)
#ifdef DEBUG
      write(*,*) 'pe ',mype,' receive of size ',isize, & 
     &          ' from ',isrc,' tag ',itag,' posted'
#endif /* DEBUG  */
      endif
      enddo
!
! All leaf blocks in sending half post a blocking send
      do lb = 1,lnblocks
      if(nodetype(lb).eq.1.and.send_half(lb)) then
        idest = mirror_blks(2,lb)
#ifdef PM_UNIQUE_MPI_TAGS
        itag = idest*maxblocks+mirror_blks(1,lb)
#else
        itag = mirror_blks(1,lb)
#endif
        call MPI_Ssend(unk_e_y(1,1,1,1,lb),isize, & 
     &                 amr_mpi_real, & 
     &                 idest,itag,amr_mpi_meshComm,ierr)
#ifdef DEBUG
      write(*,*) 'pe ',mype,' unk_e_y of size ',isize, & 
     &          ' sent to ',idest,' tag ',itag
#endif /* DEBUG  */
      endif
      enddo

      if(no_of_recvs.gt.0) & 
     &     call MPI_Waitall(no_of_recvs,recvrequest,recvstatus,ierr)

      endif
!  end of unk_e_y messages
!--------------------------------------------------
      if(nvaredge.gt.0.and.ndim.eq.3) then
!
! Set size of messages for unk_e_z data
      isize = nbndvare*(iu_bnd+1)*(ju_bnd+k2d)*ku_bnd
!
! All leaf blocks in receiving half post a non-blocking receive

      no_of_recvs = 0
      do lb = 1,lnblocks
      if(nodetype(lb).eq.1.and..not.send_half(lb)) then

        isrc = mirror_blks(2,lb)
#ifdef PM_UNIQUE_MPI_TAGS
        itag = mype*maxblocks+lb
#else
        itag = lb
#endif
        no_of_recvs = no_of_recvs+1
        call Mpi_Irecv(unk_e_z(1,1,1,1,lb),isize, & 
     &                 amr_mpi_real, & 
     &                 isrc ,itag,amr_mpi_meshComm, & 
     &                 recvrequest(no_of_recvs),ierr)
#ifdef DEBUG
      write(*,*) 'pe ',mype,' receive of size ',isize, & 
     &          ' from ',isrc,' tag ',itag,' posted'
#endif /* DEBUG  */
      endif
      enddo
!
! All leaf blocks in sending half post a blocking send
      do lb = 1,lnblocks
      if(nodetype(lb).eq.1.and.send_half(lb)) then
        idest = mirror_blks(2,lb)
#ifdef PM_UNIQUE_MPI_TAGS
        itag = idest*maxblocks+mirror_blks(1,lb)
#else
        itag = mirror_blks(1,lb)
#endif
        call MPI_Ssend(unk_e_z(1,1,1,1,lb),isize, & 
     &                 amr_mpi_real, & 
     &                 idest,itag,amr_mpi_meshComm,ierr)
#ifdef DEBUG
      write(*,*) 'pe ',mype,' unk_e_z of size ',isize, & 
     &          ' sent to ',idest,' tag ',itag
#endif /* DEBUG  */
      endif
      enddo

      if(no_of_recvs.gt.0) & 
     &     call MPI_Waitall(no_of_recvs,recvrequest,recvstatus,ierr)

      endif
!  end of unk_e_z messages
!--------------------------------------------------
!
      if(nvarcorn.gt.0) then

!
! Set size of messages for unk_n data
      isize = nbndvarc*len_blockn
!
! All leaf blocks in receiving half post a non-blocking receive

      no_of_recvs = 0
      do lb = 1,lnblocks
      if(nodetype(lb).eq.1.and..not.send_half(lb)) then

        isrc = mirror_blks(2,lb)
#ifdef PM_UNIQUE_MPI_TAGS
        itag = mype*maxblocks+lb
#else
        itag = lb
#endif
        no_of_recvs = no_of_recvs+1
        call Mpi_Irecv(unk_n(1,1,1,1,lb),isize, & 
     &                 amr_mpi_real, & 
     &                 isrc ,itag,amr_mpi_meshComm, & 
     &                 recvrequest(no_of_recvs),ierr)
#ifdef DEBUG
      write(*,*) 'pe ',mype,' receive of size ',isize, & 
     &          ' from ',isrc,' tag ',itag,' posted'
#endif /* DEBUG  */
      endif
      enddo
!
! All leaf blocks in sending half post a blocking send
      do lb = 1,lnblocks
      if(nodetype(lb).eq.1.and.send_half(lb)) then
        idest = mirror_blks(2,lb)
#ifdef PM_UNIQUE_MPI_TAGS
        itag = idest*maxblocks+mirror_blks(1,lb)
#else
        itag = mirror_blks(1,lb)
#endif
        call MPI_Ssend(unk_n(1,1,1,1,lb),isize, & 
     &                 amr_mpi_real, & 
     &                 idest,itag,amr_mpi_meshComm,ierr)
#ifdef DEBUG
      write(*,*) 'pe ',mype,' unk_n of size ',isize, & 
     &          ' sent to ',idest,' tag ',itag
#endif /* DEBUG  */
      endif
      enddo

      if(no_of_recvs.gt.0) & 
     &     call MPI_Waitall(no_of_recvs,recvrequest,recvstatus,ierr)

      endif
!  end of unk_n messages
!--------------------------------------------------


!
! Flip data in receiving blocks

      do lb = 1,lnblocks
      if(nodetype(lb).eq.1.and..not.send_half(lb)) then


! unk
        if(nvar.gt.0) then
        ttunk(:,:,:,:) = unk(:,:,:,:,lb)
        do k=kl_bnd,ku_bnd
          kk = ku_bnd - ( k - kl_bnd )
          if(iz.eq.0.or.ndim.lt.3) kk = k
        do j=jl_bnd,ju_bnd
          jj = ju_bnd - ( j - jl_bnd )
          if(iy.eq.0.or.ndim.lt.2) jj = j
        do i=il_bnd,iu_bnd
          ii = iu_bnd - ( i - il_bnd )
          if(ix.eq.0) ii = i
          unk(:nvar,i,j,k,lb) = ttunk(:nvar,ii,jj,kk) & 
     &                           *sign_cc_vars(:nvar)
        enddo
        enddo
        enddo
        endif

! facevar
        direction = 1.
        if(ix.gt.0) direction = -1.
        if(nfacevar.gt.0) then
        ttfacevarx(:,:,:,:) = facevarx(:,:,:,:,lb)
        do k=kl_bnd,ku_bnd
          kk = ku_bnd - ( k - kl_bnd )
          if(iz.eq.0.or.ndim.lt.3) kk = k
        do j=jl_bnd,ju_bnd
          jj = ju_bnd - ( j - jl_bnd )
          if(iy.eq.0.or.ndim.lt.2) jj = j
        do i=il_bnd,iu_bnd+1
          ii = iu_bnd + 1 - ( i - il_bnd )
          if(ix.eq.0) ii = i
          facevarx(:nfacevar,i,j,k,lb) = ttfacevarx(:nfacevar,ii,jj,kk) & 
     &                                      *direction
        enddo
        enddo
        enddo

        if(ndim.ge.2) then
        direction = 1.
        if(iy.gt.0) direction = -1.
        ttfacevary(:,:,:,:) = facevary(:,:,:,:,lb)
        do k=kl_bnd,ku_bnd
          kk = ku_bnd - ( k - kl_bnd )
          if(iz.eq.0.or.ndim.lt.3) kk = k
        do j=jl_bnd,ju_bnd+k2d
          jj = ju_bnd + k2d - ( j - jl_bnd )
          if(iy.eq.0.or.ndim.lt.2) jj = j
        do i=il_bnd,iu_bnd
          ii = iu_bnd - ( i - il_bnd )
          if(ix.eq.0) ii = i
          facevary(:nfacevar,i,j,k,lb) = ttfacevary(:nfacevar,ii,jj,kk) & 
     &                                      *direction
        enddo
        enddo
        enddo
        endif

        if(ndim.eq.3) then
        direction = 1.
        if(iz.gt.0) direction = -1.
        ttfacevarz(:,:,:,:) = facevarz(:,:,:,:,lb)
        do k=kl_bnd,ku_bnd+k3d
          kk = ku_bnd + k3d - ( k - kl_bnd )
          if(iz.eq.0.or.ndim.lt.3) kk = k
        do j=jl_bnd,ju_bnd
          jj = ju_bnd - ( j - jl_bnd )
          if(iy.eq.0.or.ndim.lt.2) jj = j
        do i=il_bnd,iu_bnd
          ii = iu_bnd - ( i - il_bnd )
          if(ix.eq.0) ii = i
          facevarz(:nfacevar,i,j,k,lb) = ttfacevarz(:nfacevar,ii,jj,kk) & 
     &                                      *direction
        enddo
        enddo
        enddo
        endif
        endif

! unk_e
        if(nvaredge.gt.0.and.ndim.ge.2) then
        direction = 1.
        if(ix.gt.0) direction = -1.
        ttunk_e_x(:,:,:,:) = unk_e_x(:,:,:,:,lb)
        do k=kl_bnd,ku_bnd+k3d
          kk = ku_bnd + k3d - ( k - kl_bnd )
          if(iz.eq.0.or.ndim.lt.3) kk = k
        do j=jl_bnd,ju_bnd+k2d
          jj = ju_bnd + k2d - ( j - jl_bnd )
          if(iy.eq.0.or.ndim.lt.2) jj = j
        do i=il_bnd,iu_bnd
          ii = iu_bnd - ( i - il_bnd )
          if(ix.eq.0) ii = i
          unk_e_x(:nvaredge,i,j,k,lb) = ttunk_e_x(:nvaredge,ii,jj,kk) & 
     &                                      *direction
        enddo
        enddo
        enddo

        ttunk_e_y(:,:,:,:) = unk_e_y(:,:,:,:,lb)
        direction = 1.
        if(iy.gt.0) direction = -1.
        do k=kl_bnd,ku_bnd+k3d
          kk = ku_bnd + k3d - ( k - kl_bnd )
          if(iz.eq.0.or.ndim.lt.3) kk = k
        do j=jl_bnd,ju_bnd
          jj = ju_bnd - ( j - jl_bnd )
          if(iy.eq.0.or.ndim.lt.2) jj = j
        do i=il_bnd,iu_bnd+1
          ii = iu_bnd + 1 - ( i - il_bnd )
          if(ix.eq.0) ii = i
          unk_e_y(:nvaredge,i,j,k,lb) = ttunk_e_y(:nvaredge,ii,jj,kk) & 
     &                                      *direction
        enddo
        enddo
        enddo

        if(ndim.eq.3) then
        direction = 1.
        if(iz.gt.0) direction = -1.
        ttunk_e_z(:,:,:,:) = unk_e_z(:,:,:,:,lb)
        do k=kl_bnd,ku_bnd
          kk = ku_bnd - ( k - kl_bnd )
          if(iz.eq.0.or.ndim.lt.3) kk = k
        do j=jl_bnd,ju_bnd+k2d
          jj = ju_bnd + k2d - ( j - jl_bnd )
          if(iy.eq.0.or.ndim.lt.2) jj = j
        do i=il_bnd,iu_bnd+1
          ii = iu_bnd + 1 - ( i - il_bnd )
          if(ix.eq.0) ii = i
          unk_e_z(:nvaredge,i,j,k,lb) = ttunk_e_z(:nvaredge,ii,jj,kk) & 
     &                                      *direction
        enddo
        enddo
        enddo
        endif
        endif

! unk_n
        if(nvarcorn.gt.0) then
        ttunk_n(:,:,:,:) = unk_n(:,:,:,:,lb)
        do k=kl_bnd,ku_bnd+k3d
          kk = ku_bnd + k3d - ( k - kl_bnd )
          if(iz.eq.0.or.ndim.lt.3) kk = k
        do j=jl_bnd,ju_bnd+k2d
          jj = ju_bnd + k2d - ( j - jl_bnd )
          if(iy.eq.0.or.ndim.lt.2) jj = j
        do i=il_bnd,iu_bnd+1
          ii = iu_bnd + 1 - ( i - il_bnd )
          if(ix.eq.0) ii = i
          unk_n(:nvarcorn,i,j,k,lb) = ttunk_n(:nvarcorn,ii,jj,kk)
        enddo
        enddo
        enddo
        endif

      endif
      enddo
!--------------------------------------------------

      deallocate( recvrequest )
      deallocate( recvstatus )

      return
      end subroutine amr_mirror_solution
