!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

!#define DEBUG

      subroutine amr_block_boundary_tecplot(loop)


! This routine constructs lines which trace the bounding boxes of each grid
! block.


! include file to define physical qualities of the model and mesh
      use paramesh_dimensions
      use physicaldata

! include file defining the tree
      use tree
      use io
      Use paramesh_comm_data

      use paramesh_interfaces

      implicit none

      include 'mpif.h'

! local amr variables
      integer :: nprocs,mype,loop
      character (len=6) :: filenumber
      integer :: ierror
      integer :: status1(MPI_STATUS_SIZE)
      logical :: lflag

!---------------------------------------------------------------
!
! no of arc segments for each block wire mesh
        integer,parameter :: nboxedges = 15                      ! spherical
        integer,parameter :: pts_per_edge = 25
        integer,parameter :: no_of_arcs = pts_per_edge*nboxedges

        integer,dimension (:),  allocatable :: glnblocks
        integer :: rnodetype(maxblocks)

        real :: arc_x(no_of_arcs,maxblocks)
        real :: arc_y(no_of_arcs,maxblocks)
        real :: arc_z(no_of_arcs,maxblocks)
        real :: arc_xt(no_of_arcs,maxblocks)
        real :: arc_yt(no_of_arcs,maxblocks)
        real :: arc_zt(no_of_arcs,maxblocks)

        integer :: lb,iproc,iout,k

!---------------------------------------------------------------

	Call MPI_COMM_RANK(amr_mpi_meshComm, mype, ierror)
	Call MPI_COMM_SIZE(amr_mpi_meshComm, nprocs, ierror)


        if(mype.eq.0) write(*,*) 'Started block bounding box trace'

! Initialize some message buffers
        arc_xt = 0.
        arc_yt = 0.
        arc_zt = 0.
!
!
! Initial tecplot file header.

        if(mype.eq.0) then

        write(filenumber,"(I6.6)") loop
        iout = 10

        open(unit=iout,status='unknown', & 
     &       file=trim(output_dir)//'grid.dat.'//filenumber)

        write(iout,190)
190     format('TITLE = "Grid block boundaries"')
        write(iout,191)
191     format('VARIABLES = "x", "y", "z"')
        write(iout,*) 'VARIABLES = "x", "y", "z"'

        endif

!--------
        do lb=1,lnblocks
          call draw_block_boundary(lb,pts_per_edge,arc_x,arc_y,arc_z)
        enddo


! broadcast number of blocks on each processor
        if(allocated(glnblocks)) deallocate( glnblocks )
        allocate ( glnblocks(0:nprocs-1) )
        call MPI_Allgather(lnblocks, 1,MPI_INTEGER, & 
     &                   glnblocks,1,MPI_INTEGER, & 
     &                   amr_mpi_meshComm,ierror)


        do iproc = 0,nprocs-1

!
! get remote nodetype array
          if(iproc.eq.0) then
            if(mype.eq.0) then
              rnodetype(1:lnblocks) = nodetype(1:lnblocks)
            endif
          else
            if(mype.eq.0) then
              call mpi_recv(rnodetype,maxblocks, & 
     &                      MPI_INTEGER, & 
     &                      iproc,2+iproc*10, & 
     &                      amr_mpi_meshComm,status1,ierror)
            elseif(mype.eq.iproc) then
              call mpi_send(nodetype,maxblocks, & 
     &                        MPI_INTEGER, & 
     &                        0,2+mype*10, & 
     &                        amr_mpi_meshComm,ierror)
            endif
          endif

         Call MPI_BARRIER(amr_mpi_meshComm, ierror)

          if(iproc.eq.0) then
            if(mype.eq.0) then
              arc_xt = arc_x
              arc_yt = arc_y
              arc_zt = arc_z
            endif
          else

! First exchange arc_x
            if(mype.eq.0) then
                call mpi_recv(arc_xt,no_of_arcs*maxblocks, & 
     &                        amr_mpi_real, & 
     &                        iproc,1+iproc*10, & 
     &                        amr_mpi_meshComm,status1,ierror)
            elseif(mype.eq.iproc) then
                call mpi_send(arc_x,no_of_arcs*maxblocks, & 
     &                        amr_mpi_real, & 
     &                        0,1+mype*10, & 
     &                        amr_mpi_meshComm,ierror)
            endif


! Next exchange arc_y
            if(mype.eq.0) then
                call mpi_recv(arc_yt,no_of_arcs*maxblocks, & 
     &                        amr_mpi_real, & 
     &                        iproc,2+iproc*10, & 
     &                        amr_mpi_meshComm,status1,ierror)
            elseif(mype.eq.iproc) then
                call mpi_send(arc_y,no_of_arcs*maxblocks, & 
     &                        amr_mpi_real, & 
     &                        0,2+mype*10, & 
     &                        amr_mpi_meshComm,ierror)
            endif


! Next exchange arc_z
            if(mype.eq.0) then
                call mpi_recv(arc_zt,no_of_arcs*maxblocks, & 
     &                        amr_mpi_real, & 
     &                        iproc,3+iproc*10, & 
     &                        amr_mpi_meshComm,status1,ierror)
            elseif(mype.eq.iproc) then
                call mpi_send(arc_z,no_of_arcs*maxblocks, & 
     &                        amr_mpi_real, & 
     &                        0,3+mype*10, & 
     &                        amr_mpi_meshComm,ierror)
            endif

          endif

          Call MPI_BARRIER(amr_mpi_meshComm, ierror)

          if( mype.eq.0) then
            do lb=1,glnblocks(iproc)
            if(rnodetype(lb).eq.1) then
! write data for this block into output file
              write(iout,197) no_of_arcs
197           format('ZONE I=',i4,' F=POINT')
              do k=1,no_of_arcs
                write(iout,100) arc_xt(k,lb),arc_yt(k,lb),arc_zt(k,lb)
              enddo
             endif
             enddo
           endif                                  ! end of mype=0 if test

         Call MPI_BARRIER(amr_mpi_meshComm, ierror)
        enddo                              ! end of loop over iproc

!--------------------------------------------------

        if(mype.eq.0) then
        close(unit=iout)
        endif                                ! mype if test


        call mpi_iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG, & 
     &                  amr_mpi_meshComm,lflag,status1,ierror)
        if(lflag) then
          write(*,*) 'Extra send(s) posted : ', & 
     &               'detected in amr_block_boundary_tecplot'
          call amr_abort()
        endif

        call MPI_BARRIER(amr_mpi_meshComm, ierror)


100     format(3(2x,1pe15.7))


        return
        end subroutine amr_block_boundary_tecplot




        subroutine draw_block_boundary( & 
     &                  lb,pts_per_edge,arc_x,arc_y,arc_z)

! include file to define physical qualities of the model and mesh
        use paramesh_dimensions
        use physicaldata

! include file defining the tree
        use tree

        use paramesh_interfaces



        implicit none
        integer,intent(in) :: lb,pts_per_edge
        real,intent(inout) :: arc_x(375,maxblocks)
        real,intent(inout) :: arc_y(375,maxblocks)
        real,intent(inout) :: arc_z(375,maxblocks)
!        real,intent(inout) :: arc_x(:,:)
!        real,intent(inout) :: arc_y(:,:)
!        real,intent(inout) :: arc_z(:,:)
        real :: x,y,z,r,theta,phi
        real :: r1,r2,t1,t2,p1,p2,del
        integer :: i,i0

        if(curvilinear) then

        if(spherical_pm) then

! spherical block
        i0=0

        r1 = bnd_box(1,1,lb)
        r2 = bnd_box(2,1,lb)
        t1 = bnd_box(1,2,lb)
        t2 = bnd_box(2,2,lb)
        p1 = bnd_box(1,3,lb)
        p2 = bnd_box(2,3,lb)

!------------
! low r face
!
! 1 r1-p1 edge - forward theta direction
        r = r1
        phi = p1
        do i= 1,pts_per_edge
          del = (t2-t1)/real(pts_per_edge-1)
          theta = t1 + real(i-1)*del
          call convert_rthetaphi_to_xyz(r,theta,phi,x,y,z)
          arc_x(i0+i,lb) = x
          arc_y(i0+i,lb) = y
          arc_z(i0+i,lb) = z
        enddo
        i0 = i0+pts_per_edge
! 2 r1-t2 edge - forward phi direction
        r = r1
        theta = t2
        do i= 1,pts_per_edge
          del = (p2-p1)/real(pts_per_edge-1)
          phi = p1 + real(i-1)*del
          call convert_rthetaphi_to_xyz(r,theta,phi,x,y,z)
          arc_x(i0+i,lb) = x
          arc_y(i0+i,lb) = y
          arc_z(i0+i,lb) = z
        enddo
        i0 = i0+pts_per_edge
! 3 r1-p2 edge - reverse theta direction
        r = r1
        phi = p2
        do i= 1,pts_per_edge
          del = -(t2-t1)/real(pts_per_edge-1)
          theta = t2 + real(i-1)*del
          call convert_rthetaphi_to_xyz(r,theta,phi,x,y,z)
          arc_x(i0+i,lb) = x
          arc_y(i0+i,lb) = y
          arc_z(i0+i,lb) = z
        enddo
        i0 = i0+pts_per_edge
! 4 r1-t1 edge - reverse phi
        r = r1
        theta = t1
        do i= 1,pts_per_edge
          del = -(p2-p1)/real(pts_per_edge-1)
          phi = p2 + real(i-1)*del
          call convert_rthetaphi_to_xyz(r,theta,phi,x,y,z)
          arc_x(i0+i,lb) = x
          arc_y(i0+i,lb) = y
          arc_z(i0+i,lb) = z
        enddo
        i0 = i0+pts_per_edge
!
! 5 t1-p1 edge - forward r
        theta = t1
        phi = p1
        do i= 1,pts_per_edge
          del = (r2-r1)/real(pts_per_edge-1)
          r = r1 + real(i-1)*del
          call convert_rthetaphi_to_xyz(r,theta,phi,x,y,z)
          arc_x(i0+i,lb) = x
          arc_y(i0+i,lb) = y
          arc_z(i0+i,lb) = z
        enddo
        i0 = i0+pts_per_edge
! 6 r2-p1 edge - forward theta
        r = r2
        phi = t1
        do i= 1,pts_per_edge
          del = (t2-t1)/real(pts_per_edge-1)
          theta = t1 + real(i-1)*del
          call convert_rthetaphi_to_xyz(r,theta,phi,x,y,z)
          arc_x(i0+i,lb) = x
          arc_y(i0+i,lb) = y
          arc_z(i0+i,lb) = z
        enddo
        i0 = i0+pts_per_edge
! 7 t2-p1 edge - reverse r
        theta = t2
        phi = p1
        do i= 1,pts_per_edge
          del = -(r2-r1)/real(pts_per_edge-1)
          r = r2 + real(i-1)*del
          call convert_rthetaphi_to_xyz(r,theta,phi,x,y,z)
          arc_x(i0+i,lb) = x
          arc_y(i0+i,lb) = y
          arc_z(i0+i,lb) = z
        enddo
        i0 = i0+pts_per_edge
! 8 r1-t2 edge - forward phi
        r = r1 
        theta = t2
        do i= 1,pts_per_edge
          del = (p2-p1)/real(pts_per_edge-1)
          phi = p1 + real(i-1)*del
          call convert_rthetaphi_to_xyz(r,theta,phi,x,y,z)
          arc_x(i0+i,lb) = x
          arc_y(i0+i,lb) = y
          arc_z(i0+i,lb) = z
        enddo
        i0 = i0+pts_per_edge
!
! 9 t2-p2 edge - forward r
        phi = p2
        theta = t2
        do i= 1,pts_per_edge
          del = (r2-r1)/real(pts_per_edge-1)
          r = r1 + real(i-1)*del
          call convert_rthetaphi_to_xyz(r,theta,phi,x,y,z)
          arc_x(i0+i,lb) = x
          arc_y(i0+i,lb) = y
          arc_z(i0+i,lb) = z
        enddo
        i0 = i0+pts_per_edge
! 10 r2-p2 edge - reverse theta
        r = r2
        phi = t2
        do i= 1,pts_per_edge
          del = -(t2-t1)/real(pts_per_edge-1)
          theta = t2 + real(i-1)*del
          call convert_rthetaphi_to_xyz(r,theta,phi,x,y,z)
          arc_x(i0+i,lb) = x
          arc_y(i0+i,lb) = y
          arc_z(i0+i,lb) = z
        enddo
        i0 = i0+pts_per_edge
! 11 t1-p2 edge - reverse r
        theta = t1
        phi = p2
        do i= 1,pts_per_edge
          del = -(r2-r1)/real(pts_per_edge-1)
          r = r2 + real(i-1)*del
          call convert_rthetaphi_to_xyz(r,theta,phi,x,y,z)
          arc_x(i0+i,lb) = x
          arc_y(i0+i,lb) = y
          arc_z(i0+i,lb) = z
        enddo
        i0 = i0+pts_per_edge
! 12 t1-p2 edge - forward r (ie reverse direction of last edge)
        theta = t1
        phi = p2
        do i= 1,pts_per_edge
          del = (r2-r1)/real(pts_per_edge-1)
          r = r1 + real(i-1)*del
          call convert_rthetaphi_to_xyz(r,theta,phi,x,y,z)
          arc_x(i0+i,lb) = x
          arc_y(i0+i,lb) = y
          arc_z(i0+i,lb) = z
        enddo
        i0 = i0+pts_per_edge
! 13 r2-t1 edge  - reverse phi
        r = r2
        theta = t1
        do i= 1,pts_per_edge
          del = -(p2-p1)/real(pts_per_edge-1)
          phi = p2 + real(i-1)*del
          call convert_rthetaphi_to_xyz(r,theta,phi,x,y,z)
          arc_x(i0+i,lb) = x
          arc_y(i0+i,lb) = y
          arc_z(i0+i,lb) = z
        enddo
        i0 = i0+pts_per_edge
!
! 14 r2-p1 edge - forward theta
        r = r2
        phi = p1
        do i= 1,pts_per_edge
          del = (t2-t1)/real(pts_per_edge-1)
          theta = t1 + real(i-1)*del
          call convert_rthetaphi_to_xyz(r,theta,phi,x,y,z)
          arc_x(i0+i,lb) = x
          arc_y(i0+i,lb) = y
          arc_z(i0+i,lb) = z
        enddo
        i0 = i0+pts_per_edge
! 15 r2-t2 edge - forward phi
        r = r2
        theta = t2
        do i= 1,pts_per_edge
          del = (p2-p1)/real(pts_per_edge-1)
          phi = p1 + real(i-1)*del
          call convert_rthetaphi_to_xyz(r,theta,phi,x,y,z)
          arc_x(i0+i,lb) = x
          arc_y(i0+i,lb) = y
          arc_z(i0+i,lb) = z
        enddo
        i0 = i0+pts_per_edge
!-----

        endif

        endif                      ! curvilinear

        return
        end subroutine draw_block_boundary
