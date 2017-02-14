!!****if* source/Simulation/SimulationMain/unitTest/PFFT_PoissonFD/PresDp3P
!!
!! NAME
!!
!!  PresDp3P
!!
!!
!! SYNOPSIS
!!
!! PresDp3P( integer (in) :: nx,
!!           integer (in) :: nx,
!!           integer (in) :: ny,
!!           integer (in) :: nz,
!!           real (in)    :: dx,
!!           real (in)    :: dy,
!!           real (in)    :: dz,
!!           real (inout)    :: fsrc,
!!           real (inout)    :: usol)
!!
!! DESCRIPTION
!!
!!   Solve the 3D Poisson equation, with periodic boundary conditions
!!   and 3D ffts. Used for turbulence runs.
!!   Written: Marcos Vanella, July 2007.
!!
!! ARGUMENTS
!!
!!   nx - size along iaxis
!!   ny - size along jaxis
!!   nz - size along kaxis
!!
!!   dx - grid spacing along iaxis
!!   dy - grid spacing along jaxis
!!   dz - grid spacing along kaxis
!!  fsrc- input
!!  usol- output
!!
!!***

      subroutine PresDp3P(nx,ny,nz,dx,dy,dz,fsrc,usol)

      
      implicit none

      integer, intent(in) :: nx,ny,nz
      real, intent(in) ::  dx,dy,dz
      real,intent(inout) ::  fsrc(nx,ny,nz),usol(nx,ny,nz)


      ! Local variables

      real, parameter :: pi = 3.14159265358979

      integer l, n2mh, i , j , k
      real, save ::  dx2q,dy2q,dz2q
      real, save, allocatable, dimension(:) :: kx,ky,kz,&
                                    wsavex,wsavey,wsavez     

      complex, save, allocatable, dimension(:,:,:) :: var_c

      logical, save :: firstcall = .true.

!      integer TA(2),count_rate
!      real*8  ET     


      ! Fill used data:
      if (firstcall) then

         allocate(var_c(nx,ny,nz))
         allocate(kx(nx),ky(ny),kz(nz))
         kx = 0.; ky = 0.; kz = 0.;


         !...............................................................
         !     Modified wave numbers
         !...............................................................


         ! X wavenumbers:
         dx2q = 1.0/(dx*dx);
         n2mh = (nx-2)/2;

         do l=1 , n2mh
            kx(l+1) = cos(2.*pi*real(l)/real(nx-2));
            kx(l+1) = 2.*(1.-kx(l+1))*dx2q;
         enddo
         do l=n2mh+1 , (nx-2)-1;
            kx(l+1) = kx(nx-1-l);
         enddo 


         ! Y wavenumbers:
         dy2q = 1.0/(dy*dy);
         n2mh = (ny-2)/2;

         do l=1 , n2mh
            ky(l+1) = cos(2.*pi*real(l)/real(ny-2));
            ky(l+1) = 2.*(1.-ky(l+1))*dy2q;
         enddo
         do l=n2mh+1 , (ny-2)-1
            ky(l+1) = ky(ny-1-l);
         enddo 

         ! Z wavenumbers:
         dz2q = 1.0/(dz*dz);
         n2mh = (nz-2)/2;

         do l=1 , n2mh
            kz(l+1) = cos(2.*pi*real(l)/real(nz-2));
            kz(l+1) = 2.*(1.-kz(l+1))*dz2q;
         enddo
         do l=n2mh+1 , (nz-2)-1
            kz(l+1) = kz(nz-1-l);
         enddo 


         ! Allocate wsave
         allocate(wsavex(4*(nx-2)+15))
         allocate(wsavey(4*(ny-2)+15))
         allocate(wsavez(4*(nz-2)+15))


         ! Initialize FFTS in each dir:
         call CFFTI(NX-2,WSAVEX)
         call CFFTI(NY-2,WSAVEY)
         call CFFTI(NZ-2,WSAVEZ)

         
         firstcall = .false.

      endif


      ! Dump source into complex var_c
      var_c = cmplx(fsrc) 

      ! From real to wave space ...
      ! First fft in the x direction:

!      CALL SYSTEM_CLOCK(TA(1),count_rate)
      do k = 2 , nz-1
        do j = 2 , ny-1
           
           call CFFTF(NX-2,var_c(2:nx-1,j,k),WSAVEX)

        enddo
      enddo
      var_c = var_c * cmplx(nx-2)**(-1.)

!      CALL SYSTEM_CLOCK(TA(2),count_rate)
!      ET=REAL(TA(2)-TA(1))/count_rate
!      write(*,*) 'First fft time =',ET
      

      ! Second fft in the y direction:
!      CALL SYSTEM_CLOCK(TA(1),count_rate)

      do k = 2 , nz-1
         do i = 2 ,nx-1

           call CFFTF(NY-2,var_c(i,2:ny-1,k),WSAVEY)

        enddo
      enddo        
      var_c = var_c * cmplx(ny-2)**(-1.)

!      CALL SYSTEM_CLOCK(TA(2),count_rate)
!      ET=REAL(TA(2)-TA(1))/count_rate
!      write(*,*) 'Second fft time =',ET


      ! Third fft in the z direction:
!      CALL SYSTEM_CLOCK(TA(1),count_rate)

      do j = 2 , ny-1
         do i = 2 ,nx-1

           call CFFTF(NZ-2,var_c(i,j,2:nz-1),WSAVEZ)

        enddo
      enddo                       
      var_c = var_c * cmplx(nz-2)**(-1.)    

!      CALL SYSTEM_CLOCK(TA(2),count_rate)
!      ET=REAL(TA(2)-TA(1))/count_rate
!      write(*,*) 'Third fft time =',ET      

      ! solution
      do k = 2 , nz-1
         do j = 2 , ny-1
            do i = 2 ,nx-1

       var_c(i,j,k)= -var_c(i,j,k)*cmplx((kx(i-1)+ky(j-1)+kz(k-1))**(-1.))  

               if ( abs(aimag(var_c(i,j,k))) > 10.**(-10.)) then
                  write(*,*) 'Scalar I, J , K position =',i-1,j-1,k-1
                  write(*,*) 'Scalar Poisson Val=',var_c(i,j,k)
               endif           
           
            enddo
         enddo
      enddo


      var_c(2,2,2) = cmplx(0.);


      ! From wave space to real space...
      ! Inverse ffts:
      ! First ifft in the x direction:
      do k = 2 , nz-1
        do j = 2 , ny-1
           
           call CFFTB(NX-2,var_c(2:nx-1,j,k),WSAVEX)

        enddo
      enddo

      ! Second ifft in the y direction:
      do k = 2 , nz-1
         do i = 2 ,nx-1

           call CFFTB(NY-2,var_c(i,2:ny-1,k),WSAVEY)

        enddo
      enddo        

      ! Third ifft in the z direction:
      do j = 2 , ny-1
         do i = 2 ,nx-1

           call CFFTB(NZ-2,var_c(i,j,2:nz-1),WSAVEZ)

        enddo
      enddo                       

      
      ! Dump solution into Usol:
      usol = real(var_c)


      return

      end subroutine
