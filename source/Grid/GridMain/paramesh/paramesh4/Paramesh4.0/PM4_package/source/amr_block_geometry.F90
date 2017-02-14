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

         subroutine amr_block_geometry(lb,pe)

!-----------------------------------------------------------
!
! This routine computes cell volumes, areas and edge lengths
! for various grid geometries, for the specified local block lb.

!
!
! Cartesian :
!     coord 1      x
!     coord 2      y
!     coord 3      z
!
! Cylindrical :  (NOTE: Coordinates are in this order to support a 2-d coordinate
!                       system which is axisymmetric about the z axis
!     coord 1      r 
!     coord 2      z
!     coord 3      theta
!
! Spherical :
!     coord 1      r
!     coord 2      theta
!     coord 3      phi           (azimuthal)
!
! Polar (2D) :
!     coord 1      r 
!     coord 2      theta
!     coord 3      z             (has only 1 grid cell in this direction)
!
!
! Written : Peter MacNeice      December 2001
! Cylindrical Axisymmetric added by Sergey Pancheshnyi and Kevin Olson (April 2006)

!-----------------------------------------------------------

      use paramesh_dimensions
      use physicaldata
      use tree
      use workspace
      use paramesh_comm_data
      implicit none

      Include 'mpif.h'

      integer, intent(in) :: lb,pe

!-----------------------------------------------------------

! Local arrays and variables

         real :: cell_face_coord1w(ilw1:iuw1+1)
         real :: cell_face_coord2w(jlw1:juw1+k2d)
         real :: cell_face_coord3w(klw1:kuw1+k3d)
      

         real :: del
         real :: cell_vol_1  ,cell_vol_2  ,cell_vol_3
         real :: cell_area1_1,cell_area1_2,cell_area1_3
         real :: cell_area2_1,cell_area2_2,cell_area2_3
         real :: cell_area3_1,cell_area3_2,cell_area3_3
         real :: cell_leng1_1,cell_leng1_2,cell_leng1_3
         real :: cell_leng2_1,cell_leng2_2,cell_leng2_3
         real :: cell_leng3_1,cell_leng3_2,cell_leng3_3
 
         real,save :: cbnd_box(2,3)

         integer :: ierr_trap
         integer :: mype,ierr
         integer :: lb0,pe0,iloc
         integer :: i, j, k
         logical :: lfound
         real :: dx, dy, dz
         real :: xleft,yleft,zleft
         real :: eps

!-----------------------

         eps = tiny(eps)

         Call MPI_COMM_RANK(amr_mpi_meshComm, mype, ierr)

         ierr_trap = 0
         if (cartesian_pm) then
            ierr_trap = ierr_trap + 1
         endif
         if (spherical_pm) then
            ierr_trap = ierr_trap + 1
         endif
         if (cylindrical_pm) then
            ierr_trap = ierr_trap + 1
         endif
         if (polar_pm) then
            ierr_trap = ierr_trap + 1
            if(ndim.ne.2) ierr_trap = ierr_trap + 1
         endif
         if(ierr_trap.gt.1) then
           write(*,*) 'Paramesh ERROR : amr_block_geometry. ', & 
     &           'Inconsistent choice of curvilinear coord.'
           call amr_abort()
         endif

!-----------------------
!
        lfound = .false.
        if(pe.eq.mype) then
          if(lb.le.lnblocks) then
            lfound = .true.
            lb0 = lb
            pe0 = mype
          elseif(lb.ge.strt_buffer.and.lb.le.last_buffer) then
            lfound = .true.
            lb0 = lb
            pe0 = mype
          endif
        else
          iloc = strt_buffer
          do while( (iloc.le.last_buffer) .and. & 
     &              (.not.lfound)  )
            if(laddress(1,iloc).eq.lb.and. & 
     &         laddress(2,iloc).eq.pe ) lfound = .true.
            if(.not.lfound) iloc = iloc + 1
          enddo
          if(lfound) then
            lb0 = iloc
            pe0 = mype
          endif
        endif

        if(.not.lfound) then
          write(*,*) 'amr_block_geometry ERROR : blk ', & 
     &        lb,pe,' not found on pe ',mype, & 
     &        ' strt_buffer:last_buffer ',strt_buffer,last_buffer, & 
     &        ' laddress ',laddress(:,strt_buffer:last_buffer)
          call amr_abort()
        endif
        cbnd_box(:,:) = bnd_box(:,:,lb0)

!-----------------------
! compute coords of cell interfaces


! for first coordinate direction
         del = (cbnd_box(2,1)-cbnd_box(1,1))/real(nxb)
         do i = il_bnd1,iu_bnd1+1
           cell_face_coord1(i) = cbnd_box(1,1) + del*real(i-1-nguard)
         enddo
         do i = ilw1,iuw1+1
           cell_face_coord1w(i) = cbnd_box(1,1)  & 
     &                            + del*real(i-1-nguard_work)
         enddo

         dx = (cbnd_box(2,1)-cbnd_box(1,1))/real(nxb)

! for second coordinate direction
         cell_face_coord2 = 0.
         cell_face_coord2w = 0.
         if(ndim.ge.2) then
         del = (cbnd_box(2,2)-cbnd_box(1,2))/real(nyb)
         yleft = coord(2,lb0) - bsize(2,lb0)/2.
         do j = jl_bnd1,ju_bnd1+1
           cell_face_coord2(j) = cbnd_box(1,2) + del*real(j-1-nguard)
         enddo
         do j = jlw1,juw1+1
           cell_face_coord2w(j) = cbnd_box(1,2)  & 
     &                            + del*real(j-1-nguard_work)
         enddo
         endif

         dy = (cbnd_box(2,2)-cbnd_box(1,2))/real(nyb)

! for third coordinate direction
         cell_face_coord3 = 0.
         cell_face_coord3w = 0.
         if(ndim.eq.3) then
         del = (cbnd_box(2,3)-cbnd_box(1,3))/real(nzb)
         do k = kl_bnd1,ku_bnd1+1
           cell_face_coord3(k) = cbnd_box(1,3) + del*real(k-1-nguard)
         enddo
         do k = klw1,kuw1+1
           cell_face_coord3w(k) = cbnd_box(1,3)  & 
     &                            + del*real(k-1-nguard_work)
         enddo
         endif

         dz = (cbnd_box(2,3)-cbnd_box(1,3))/real(nzb)

!-----------------------

! Apply any user specified coordinate transformation
         call user_coord_transfm(lb0,pe0)


!-----------------------
! compute cell volumes

! Note the style used here to compute cell_vol. We
! specify dependence of cell_vol on coord 1 in cell_vol_1,
! specify dependence of cell_vol on coord 2 in cell_vol_2,
! specify dependence of cell_vol on coord 3 in cell_vol_3.
! This style is used throughout this routine.


! first cell volumes for use with UNK data structure
         do k = kl_bnd1,ku_bnd1
         do j = jl_bnd1,ju_bnd1
         do i = il_bnd1,iu_bnd1

           cell_vol_1 = 1.
           cell_vol_2 = 1.
           cell_vol_3 = 1.

           if (cartesian_pm) then
           cell_vol_1 = dx 
           if(ndim.ge.2) & 
     &       cell_vol_2 = dy  
           if(ndim.eq.3) & 
     &       cell_vol_3 = dz
           if(ndim == 2 .AND. l2p5d.eq.1) cell_vol_3 =  & 
     &                    cbnd_box(2,3)-cbnd_box(1,3)
           endif  ! End If (cartesian_pm)

           if (spherical_pm) then
           cell_vol_1 = ( cell_face_coord1(i+1)**3 -  & 
                          cell_face_coord1(i)**3 ) /3.
           if(ndim.ge.2) & 
             cell_vol_2 = cos( cell_face_coord2(j)   ) - & 
                          cos( cell_face_coord2(j+1) )
           if(ndim.eq.3) & 
              cell_vol_3 = cell_face_coord3(k+k3d) -  & 
                           cell_face_coord3(k)
           if(ndim == 2 .AND. l2p5d.eq.1) cell_vol_3 =  & 
                          cbnd_box(2,3)-cbnd_box(1,3)
           endif  ! End If (spherical_pm)

           if (cylindrical_pm) then
           cell_vol_1 = ( cell_face_coord1(i+1)**2 -  & 
                          cell_face_coord1(i)**2 )*.5
           if(ndim.ge.2) & 
             cell_vol_2 =  cell_face_coord2(j+1) - & 
                           cell_face_coord2(j)
           if(ndim.eq.3) & 
             cell_vol_3 =  cell_face_coord3(k+k3d) - & 
                           cell_face_coord3(k)
           if(ndim == 2 .AND. l2p5d.eq.1) cell_vol_3 =  & 
                          cbnd_box(2,3)-cbnd_box(1,3)
           endif  ! End If (cylindrical_pm)

           if (polar_pm) then
           cell_vol_1 = ( cell_face_coord1(i+1)**2 -  & 
     &                    cell_face_coord1(i)**2 )*.5
           if(ndim.ge.2)  & 
     &          cell_vol_2 =  cell_face_coord2(j+1) - & 
     &                        cell_face_coord2(j)
           endif

           cell_vol(i,j,k) = max(abs(cell_vol_1 * cell_vol_2  & 
     &                                          * cell_vol_3), & 
     &                           eps) 

         enddo
         enddo
         enddo



! now cell volumes for use with WORK data structure
         do k = klw1,kuw1
         do j = jlw1,juw1
         do i = ilw1,iuw1

           cell_vol_1 = 1.
           cell_vol_2 = 1.
           cell_vol_3 = 1.

           if (cartesian_pm) then
           cell_vol_1 = dx
           if(ndim.ge.2) & 
             cell_vol_2 = dy
           if(ndim.eq.3) & 
             cell_vol_3 = dz
           if(ndim == 2 .AND. l2p5d.eq.1) cell_vol_3 =  & 
                          cbnd_box(2,3)-cbnd_box(1,3)
           endif  ! End If (cartesian_pm)

           if (spherical_pm) then
           cell_vol_1 = ( cell_face_coord1w(i+1)**3 - & 
                          cell_face_coord1w(i)**3 ) /3.
           if(ndim.ge.2) & 
             cell_vol_2 = cos( cell_face_coord2w(j)   ) - & 
                          cos( cell_face_coord2w(j+1) )
           if(ndim.eq.3) & 
             cell_vol_3 = cell_face_coord3w(k+k3d) - & 
                          cell_face_coord3w(k)
           if(ndim == 2 .AND. l2p5d.eq.1) cell_vol_3 =  & 
                          cbnd_box(2,3)-cbnd_box(1,3)
           endif  ! End If (spherical_pm)

           if (polar_pm) then
           cell_vol_1 = ( cell_face_coord1w(i+1)**2 - & 
     &                    cell_face_coord1w(i)**2 )*.5
           if(ndim.ge.2) & 
     &       cell_vol_2 =  cell_face_coord2w(j+1) - & 
     &                     cell_face_coord2w(j)
           endif

           if (cylindrical_pm) then
           ! INT(rdr) * INT(dz) * INT(d theta)
           cell_vol_1 = ( cell_face_coord1w(i+1)**2 - & 
                          cell_face_coord1w(i)**2 )*.5
           if(ndim.ge.2) & 
             cell_vol_2 =  cell_face_coord2w(j+1) - & 
                           cell_face_coord2w(j)
           if(ndim.eq.3) & 
             cell_vol_3 =  cell_face_coord3w(k+k3d) - & 
                           cell_face_coord3w(k)
           if(ndim == 2 .AND. l2p5d.eq.1) cell_vol_3 =  & 
                          cbnd_box(2,3)-cbnd_box(1,3)
           endif  ! End If (cylindrical_pm)

           cell_vol_w(i,j,k) = max(abs(cell_vol_1 * cell_vol_2 & 
     &                                            * cell_vol_3), & 
     &                             eps)

         enddo
         enddo
         enddo



!-----------------------
! Compute cell face areas


! compute cell area of faces perpendicular to first coord axis
         do k = kl_bnd1,ku_bnd1
         do j = jl_bnd1,ju_bnd1
         do i = il_bnd1,iu_bnd1+1

           cell_area1_1 = 1.
           cell_area1_2 = 1.
           cell_area1_3 = 1.

           if (cartesian_pm) then
           cell_area1_1 =  1.
           if(ndim.ge.2) & 
             cell_area1_2 =  dy
           if(ndim.eq.3) & 
             cell_area1_3 = dz
           if(ndim == 2 .AND. l2p5d.eq.1) cell_area1_3 =  & 
                          cbnd_box(2,3)-cbnd_box(1,3)
           endif  ! End If (cartesian_pm)
           
           if (spherical_pm) then
           cell_area1_1 = cell_face_coord1(i)**2
           if(ndim.ge.2) & 
             cell_area1_2 = cos( cell_face_coord2(j)   ) - & 
                            cos( cell_face_coord2(j+1) )
           if(ndim.eq.3) & 
             cell_area1_3 = cell_face_coord3(k+k3d) - & 
                            cell_face_coord3(k)
           if(ndim == 2 .AND. l2p5d.eq.1) cell_area1_3 =  & 
                          cbnd_box(2,3)-cbnd_box(1,3)
           endif  ! End If (spherical_pm)

           if (polar_pm) then
           cell_area1_1 =  cell_face_coord1(i) 
           if(ndim.ge.2) & 
     &       cell_area1_2 =  cell_face_coord2(j+1) - & 
     &                       cell_face_coord2(j)
           endif

           if (cylindrical_pm) then ! perp to r
           ! INT(dz) * r*INT(d theta)
           cell_area1_1 =  cell_face_coord1(i) 
           if(ndim.ge.2) & 
             cell_area1_2 =  cell_face_coord2(j+1) - & 
                             cell_face_coord2(j)
           If (ndim == 2 .and. l2p5d == 1) cell_area1_3 = &
                           cbnd_box(2,3)-cbnd_box(1,3)
           if (ndim.eq.3) then
                cell_area1_3 =  cell_face_coord3(k+k3d) - & 
                                cell_face_coord3(k)
           end if

           end if

           cell_area1(i,j,k) = max(abs(cell_area1_1 * cell_area1_2  & 
     &                                              * cell_area1_3), & 
     &                             eps)

         enddo
         enddo
         enddo



! compute cell area of faces perpendicular to second coord axis
         do k = kl_bnd1,ku_bnd1
         do j = jl_bnd1,ju_bnd1+k2d
         do i = il_bnd1,iu_bnd1

           cell_area2_1 = 1.
           cell_area2_2 = 1.
           cell_area2_3 = 1.

           if (cartesian_pm) then
           cell_area2_1 = dx
           if(ndim.eq.3) & 
             cell_area2_3 = dz
           if(ndim == 2 .AND. l2p5d.eq.1) cell_area2_3 =  & 
                          cbnd_box(2,3)-cbnd_box(1,3)
           endif  ! End If (cartesiam_pm)

           if (spherical_pm) then
           cell_area2_1 = (cell_face_coord1(i+1)-cell_face_coord1(i)) & 
     &            *(cell_face_coord1(i)+cell_face_coord1(i+1))*.5
           cell_area2_2 = sin( cell_face_coord2(j) )
           if(ndim.eq.3) & 
             cell_area2_3 = cell_face_coord3(k+k3d) - & 
                            cell_face_coord3(k)
           if(ndim == 2 .AND. l2p5d.eq.1) cell_area2_3 =  & 
                          cbnd_box(2,3)-cbnd_box(1,3)
           endif  ! End If (spherical_pm)

           if (polar_pm) then
           cell_area2_1 =  cell_face_coord1(i+1) - & 
     &                     cell_face_coord1(i)
           endif

           if (cylindrical_pm) then ! perp to z
           ! INT(rdr) * INT(d theta)
           cell_area2_1 = ( cell_face_coord1(i+1)**2 -  & 
                            cell_face_coord1(i)**2 )*.5
           if(ndim.eq.3) & 
             cell_area2_3 =  cell_face_coord3(k+k3d) - & 
                             cell_face_coord3(k)
           If (ndim == 2 .and. l2p5d == 1) cell_area2_3 =                              &
                           cbnd_box(2,3)-cbnd_box(1,3)
           endif  ! End If (cylindrical_pm)

           cell_area2(i,j,k) = max(abs(cell_area2_1 * cell_area2_2  & 
     &                                              * cell_area2_3), & 
     &                             eps) 

         enddo
         enddo
         enddo



! compute cell area of faces perpendicular to third coord axis
         do k = kl_bnd1,ku_bnd1+k3d
         do j = jl_bnd1,ju_bnd1
         do i = il_bnd1,iu_bnd1

           cell_area3_1 =  1.
           cell_area3_2 =  1.
           cell_area3_3 =  1.
           
           if (cartesian_pm) then
           cell_area3_1 = dx
!           cell_area3_1 =  cell_face_coord1(i+1) -
!     .                     cell_face_coord1(i)
           if(ndim.ge.2) & 
     &       cell_area3_2 = dy
!     .       cell_area3_2 =  cell_face_coord2(j+1) -
!     .                       cell_face_coord2(j)
           endif

           if (spherical_pm) then
           cell_area3_1 = (cell_face_coord1(i+1)-cell_face_coord1(i)) & 
     &            *(cell_face_coord1(i)+cell_face_coord1(i+1))*.5
           if(ndim.ge.2) & 
     &       cell_area3_2 = cell_face_coord2(j+1) - & 
     &                      cell_face_coord2(j)
           endif
           
           if (polar_pm) then
           cell_area3_1 = ( cell_face_coord1(i+1)**2 -  & 
     &                      cell_face_coord1(i)**2 )*.5
           if(ndim.ge.2) & 
     &       cell_area3_2 =  cell_face_coord2(j+1) - & 
     &                       cell_face_coord2(j)
           endif

           if (cylindrical_pm) then  ! perp to theta
           ! INT(dr) * INT(dz)
           cell_area3_1 =  cell_face_coord1(i+1) - & 
     &                     cell_face_coord1(i)
           if(ndim.ge.2) & 
     &       cell_area3_2 =  cell_face_coord2(j+k2d) - & 
     &                       cell_face_coord2(j)
           end if

           cell_area3(i,j,k) = max(abs(cell_area3_1 * cell_area3_2  & 
     &                                              * cell_area3_3), & 
     &                             eps)

         enddo
         enddo
         enddo

!-----------------------
! Compute cell edge lengths


! compute edge length in direction of first coord axis
         do k = kl_bnd1,ku_bnd1+k3d
         do j = jl_bnd1,ju_bnd1+k2d
         do i = il_bnd1,iu_bnd1

           cell_leng1_1 =  1.
           cell_leng1_2 =  1.
           cell_leng1_3 =  1.

           if (cartesian_pm) then
           cell_leng1_1 = dx
!           cell_leng1_1 =  cell_face_coord1(i+1) -
!     .                     cell_face_coord1(i)
           endif

           if (spherical_pm) then
           cell_leng1_1 =  cell_face_coord1(i+1) - & 
     &                     cell_face_coord1(i)
           endif

           if (polar_pm) then
           cell_leng1_1 =  cell_face_coord1(i+1) - & 
     &                     cell_face_coord1(i)
           endif

           if (cylindrical_pm) then
           cell_leng1_1 =  cell_face_coord1(i+1) - & 
     &                     cell_face_coord1(i)
           endif

           cell_leng1(i,j,k) = max(cell_leng1_1 * cell_leng1_2  & 
     &                                          * cell_leng1_3, & 
     &                             eps)

         enddo
         enddo
         enddo


! compute edge length in direction of second coord axis
         do k = kl_bnd1,ku_bnd1+k3d
         do j = jl_bnd1,ju_bnd1
         do i = il_bnd1,iu_bnd1+1

           cell_leng2_1 =  1.
           cell_leng2_2 =  1.
           cell_leng2_3 =  1.

           if (cartesian_pm) then
           if(ndim.ge.2) & 
     &       cell_leng2_2 = dy
!     .       cell_leng2_2 =  cell_face_coord2(j+k2d) -
!     .                       cell_face_coord2(j)
           endif

           if (spherical_pm) then
           cell_leng2_1 =  cell_face_coord1(i)
           if(ndim.ge.2) & 
     &       cell_leng2_2 =  cell_face_coord2(j+1) - & 
     &                       cell_face_coord2(j)
           endif

           if (polar_pm) then
           cell_leng2_1 =  cell_face_coord1(i)
           if(ndim.ge.2) & 
     &       cell_leng2_2 =  cell_face_coord2(j+1) - & 
     &                       cell_face_coord2(j)
           endif

           if (cylindrical_pm) then
           if(ndim.ge.2) & 
     &       cell_leng2_2 =  cell_face_coord2(j+k3d) - & 
     &                       cell_face_coord2(j)
           if(ndim == 2 .AND. l2p5d.eq.1) cell_leng2_2 =  & 
     &                    cbnd_box(2,2)-cbnd_box(1,2)
           endif

           cell_leng2(i,j,k) = max(cell_leng2_1 * cell_leng2_2 & 
     &                                          * cell_leng2_3, & 
     &                             eps)

         enddo
         enddo
         enddo


! compute edge length in direction of third coord axis
         do k = kl_bnd1,ku_bnd1
         do j = jl_bnd1,ju_bnd1+k2d
         do i = il_bnd1,iu_bnd1+1

           cell_leng3_1 =  1.
           cell_leng3_2 =  1.
           cell_leng3_3 =  1.

           if (cartesian_pm) then
           if(ndim.eq.3) & 
             cell_leng3_3 = dz
           if(ndim == 2 .AND. l2p5d.eq.1) cell_leng3_3 =  & 
                          cbnd_box(2,3)-cbnd_box(1,3)
           endif  ! End If (cartesian_pm)

           if (spherical_pm) then
           cell_leng3_1 =  cell_face_coord1(i)
           cell_leng3_2 =  sin( cell_face_coord2(j) )
           if(ndim.eq.3) & 
             cell_leng3_3 =  cell_face_coord3(k+k3d) - & 
                             cell_face_coord3(k)
           if(ndim == 2 .AND. l2p5d.eq.1) cell_leng3_3 =  & 
                          cbnd_box(2,3)-cbnd_box(1,3)
           endif  ! End If (spherical_pm)

           if (cylindrical_pm) then
           cell_leng3_1 =  cell_face_coord1(i)
           if(ndim.eq.3) & 
             cell_leng3_3 =  cell_face_coord3(k+k3d) - & 
                             cell_face_coord3(k)
           If (ndim == 2 .and. l2p5d == 1) cell_leng3_3 =                              &
                           cbnd_box(2,3)-cbnd_box(1,3)
           end if  ! End If (cylindrical_pm)

           cell_leng3(i,j,k) = max(cell_leng3_1 * cell_leng3_2  & 
     &                                          * cell_leng3_3, & 
     &                             eps)

         enddo
         enddo
         enddo


!-----------------------
         return
         end subroutine amr_block_geometry


