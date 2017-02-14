!!****if* source/Particles/ParticlesMain/active/Sink/pt_sinkPrepareEwald
!!
!! NAME
!!
!!  pt_sinkPrepareEwald
!!
!! SYNOPSIS
!!
!!  call pt_sinkPrepareEwald()
!!
!! DESCRIPTION
!!
!!  Generates a cube with gravitational acceleration corrections due to periodic boundary
!!  conditions. It implements the Ewald decomposition applied to sink particle accelerations.
!!  (Ewald 1921; see e.g., Hernquist et al. 1991, ApJS 75, 231 and Suto 1993, PTP 90, 1173).
!!  The Ewald correction field is generated here for a given geometry of the computational
!!  domain with side lengths Lx, Ly, Lz in case of periodic boundary conditions (PBCs).
!!  Only fully periodic BCs in all directions are supported. The lookup table is then
!!  interpolated to find force corrections in pt_sinkEwaldCorrection(). The Ewald field
!!  only has to be generated once and will be written to file for later use on FLASH restart.
!!
!! ARGUMENTS
!!
!! NOTES
!!
!!   written by Christoph Federrath, 2014
!!
!!***

subroutine pt_sinkPrepareEwald()

  use Particles_sinkData
  use Driver_data, ONLY : dr_globalMe, dr_globalNumProcs
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Logfile_interface, ONLY : Logfile_stamp

  implicit none

  include "Flash_mpi.h"

#include "Flash.h"
#include "Particles.h"
#include "constants.h"

  integer :: MyPE, NumPEs
  integer :: nx, ny, nz, es_nx, es_ny, es_nz, es_nrx, es_nry, es_nrz
  integer :: es_nfx, es_nfy, es_nfz, i, j, k, ni, nj, nk, hi, hj, hk
  integer :: index, chunk, ierr
  real :: xmin, xmax, ymin, ymax, zmin, zmax, Lx, Ly, Lz, es_radius2
  real :: ewald_xmax, ewald_ymax, ewald_zmax, Linv, Linv2, alpha, ewc
  real :: ratio_p1, ratio_p2, ratio_pinv1, ratio_pinv2
  real :: cr1, cr2, cr3, cf1, cf2, cf3, ewald_term
  real :: x, y, z, xni, yni, zni, rni, alpharni, hsq, kx, r

  integer :: funit
  integer :: ut_getFreeFileUnit
  logical :: file_exists
  character(len=MAX_STRING_LENGTH) :: strbuff, grav_boundary_type, openflag

  real, allocatable, dimension(:,:,:) :: loc_EwaldFieldX, &
                                         loc_EwaldFieldY, &
                                         loc_EwaldFieldZ

  intrinsic erfc

! ============================================

  MyPE = dr_globalMe
  NumPEs = dr_globalNumProcs

  call RuntimeParameters_get("grav_boundary_type", grav_boundary_type)

  if (grav_boundary_type .ne. "periodic") return

  call RuntimeParameters_get("sink_EwaldFieldNx", sink_EwaldNx)
  call RuntimeParameters_get("sink_EwaldFieldNy", sink_EwaldNy)
  call RuntimeParameters_get("sink_EwaldFieldNz", sink_EwaldNz)
  call RuntimeParameters_get("sink_EwaldSeriesN", sink_EwaldSeriesN)
  call RuntimeParameters_get("sink_EwaldFileName", sink_EwaldFileName)
  call RuntimeParameters_get("xmin", xmin)
  call RuntimeParameters_get("xmax", xmax)
  call RuntimeParameters_get("ymin", ymin)
  call RuntimeParameters_get("ymax", ymax)
  call RuntimeParameters_get("zmin", zmin)
  call RuntimeParameters_get("zmax", zmax)

  nx = sink_EwaldNx
  ny = sink_EwaldNy
  nz = sink_EwaldNz

  allocate(sink_EwaldFieldX(0:nx-1, 0:ny-1, 0:nz-1))
  allocate(sink_EwaldFieldY(0:nx-1, 0:ny-1, 0:nz-1))
  allocate(sink_EwaldFieldZ(0:nx-1, 0:ny-1, 0:nz-1))

  Lx = xmax - xmin
  Ly = ymax - ymin
  Lz = zmax - zmin

  es_nx = sink_EwaldSeriesN
  es_ny = sink_EwaldSeriesN
  es_nz = sink_EwaldSeriesN

  es_radius2 = real(sink_EwaldSeriesN)*real(sink_EwaldSeriesN)

  ewald_xmax = 0.5*Lx
  ewald_ymax = 0.5*Ly
  ewald_zmax = 0.5*Lz

  ! used in pt_sinkEwaldCorrection
  sink_EwaldDxI = real(nx-1) / ewald_xmax
  sink_EwaldDyI = real(ny-1) / ewald_ymax
  sink_EwaldDzI = real(nz-1) / ewald_zmax

  Linv = 1.0/Lx 
  ratio_p1 = Ly/Lx
  ratio_p2 = Lz/Lx
  ratio_pinv1 = 1.0/ratio_p1
  ratio_pinv2 = 1.0/ratio_p2
  es_nrx = sink_EwaldSeriesN
  es_nry = ceiling(sink_EwaldSeriesN/ratio_p1)
  es_nrz = ceiling(sink_EwaldSeriesN/ratio_p2)
  es_nfx = sink_EwaldSeriesN
  es_nfy = ceiling(sink_EwaldSeriesN*ratio_p1)
  es_nfz = ceiling(sink_EwaldSeriesN*ratio_p2)
  cr1 = 1.0
  cr2 = ratio_p1**2
  cr3 = ratio_p2**2
  cf1 = 1.0
  cf2 = 1.0/cr2
  cf3 = 1.0/cr3

  Linv2 = Linv/(PI*ratio_p1*ratio_p2)
  alpha = 2.0*Linv
  ewc = (PI*Linv/alpha)**2

  ! check if file with Ewald fields already exists
  if (MyPE == MASTER_PE) inquire(file=trim(sink_EwaldFileName), exist=file_exists)
  call MPI_Bcast(file_exists, 1, FLASH_LOGICAL, MASTER_PE, MPI_COMM_WORLD, ierr)

  if (file_exists) then ! read Ewald field from file

    if (MyPE == MASTER_PE) then

      write (strbuff, '("File with Ewald correction field exists. Reading, size = ", i4, " x ", i4, " x ", i4)') &
             nx, ny, nz
      call Logfile_stamp(strbuff, "[Sink particles]")
      write (*,*) 'pt_sinkPrepareEwald: '//strbuff

      funit = ut_getFreeFileUnit ()
      open(unit=funit, file=trim(sink_EwaldFileName), status='old')
      read(funit,*) ! header
      do k = 0, nz-1
        do j = 0, ny-1
          do i = 0, nx-1
            read(funit,'(6(e24.17, 1x))') x, y, z, &
                  sink_EwaldFieldX(i,j,k), sink_EwaldFieldY(i,j,k), sink_EwaldFieldZ(i,j,k)
          enddo
          read(funit,*) ! blank line
        enddo
        read(funit,*) ! blank line
      enddo
      close(unit=funit)
    endif
    call MPI_Bcast(sink_EwaldFieldX, nx*ny*nz, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(sink_EwaldFieldY, nx*ny*nz, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(sink_EwaldFieldZ, nx*ny*nz, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)

  else ! generate the Ewald field

    allocate(loc_EwaldFieldX(0:nx-1, 0:ny-1, 0:nz-1))
    allocate(loc_EwaldFieldY(0:nx-1, 0:ny-1, 0:nz-1))
    allocate(loc_EwaldFieldZ(0:nx-1, 0:ny-1, 0:nz-1))
    if (MyPE == MASTER_PE) then
      write (strbuff, '("Computing Ewald correction field (takes time...), size = ", i4, " x ", i4, " x ", i4)') &
             nx, ny, nz
      call Logfile_stamp(strbuff, "[Sink particles]")
      write (*,*) 'pt_sinkPrepareEwald: '//strbuff
    endif

    ! reset fields
    sink_EwaldFieldX(:,:,:) = 0.0
    sink_EwaldFieldY(:,:,:) = 0.0
    sink_EwaldFieldZ(:,:,:) = 0.0
    loc_EwaldFieldX(:,:,:) = 0.0
    loc_EwaldFieldY(:,:,:) = 0.0
    loc_EwaldFieldZ(:,:,:) = 0.0

    do k = 0, nz-1
      if (MyPE == MASTER_PE) write (*,*) 'pt_sinkPrepareEwald: MASTER_PE working on k = ', k+1, ' of ', nz, '...'
      do j = 0, ny-1
        do i = 0, nx-1

          ! calculate 1D index for 3D field
          index = k*ny*nx + j*nx + i
          chunk = nx*ny*nz/NumPEs + 1

          if ((index .ge. MyPE*chunk) .and. (index .lt. (MyPE+1)*chunk)) then

            x = i * ewald_xmax / real(nx-1)
            y = j * ewald_ymax / real(ny-1)
            z = k * ewald_zmax / real(nz-1)

            ! first term
            do nk = -es_nrz, es_nrz
              do nj = -es_nry, es_nry
                do ni = -es_nrx, es_nrx

                  ! terms with non-negligible contributions must lie inside ellipse
                  if ((cr1*ni**2 + cr2*nj**2 + cr3*nk**2) .le. es_radius2) then

                    xni = x + ni*Lx
                    yni = y + nj*Ly
                    zni = z + nk*Lz
                    rni = sqrt(xni*xni + yni*yni + zni*zni)

                    ! potential
                    ! ewald_term = erfc(alpha*rni) / (rni + 1d-99)

                    ! acceleration
                    alpharni = alpha*rni
                    ewald_term = ( erfc(alpharni) + 2.0*alpharni/sqrt(PI)*exp(-alpharni**2) ) / &
                                 (rni + 1e-99)**3

                    loc_EwaldFieldX(i,j,k) = loc_EwaldFieldX(i,j,k) + ewald_term * xni
                    loc_EwaldFieldY(i,j,k) = loc_EwaldFieldY(i,j,k) + ewald_term * yni
                    loc_EwaldFieldZ(i,j,k) = loc_EwaldFieldZ(i,j,k) + ewald_term * zni

                  endif

                enddo ! ni
              enddo ! nj
            enddo ! nk
            
            ! second term
            do hk = -es_nfz, es_nfz
              do hj = -es_nfy, es_nfy
                do hi = -es_nfx, es_nfx

                  if ( ((hi**2 + hj**2 + hk**2) .gt. 0) .and. &
                       ((cf1*hi**2 + cf2*hj**2 + cf3*hk**2) .le. es_radius2) ) then

                    hsq = hi**2 + (hj*ratio_pinv1)**2 + (hk*ratio_pinv2)**2
                    kx = 2.0*PI*Linv * (hi*x + hj*y*ratio_pinv1 + hk*z*ratio_pinv2)

                    ! potential
                    ! ewald_term = Linv2 * exp(-ewc*hsq) * cos(kx) / hsq

                    ! acceleration
                    ewald_term = 2.0*Linv2**2/hsq * exp(-ewc*hsq) * sin(kx) 

                    loc_EwaldFieldX(i,j,k) = loc_EwaldFieldX(i,j,k) + ewald_term * hi
                    loc_EwaldFieldY(i,j,k) = loc_EwaldFieldY(i,j,k) + ewald_term * hj
                    loc_EwaldFieldZ(i,j,k) = loc_EwaldFieldZ(i,j,k) + ewald_term * hk

                  endif ! not zero

                enddo ! hi
              enddo ! hj
            enddo ! hk

            ! Construct the force correction (see e.g., Hernquist et al. 1991).
            ! This is relatively smooth and better for interpolation in pt_sinkEwaldCorrection().
            r = sqrt(x**2 + y**2 + z**2)
            loc_EwaldFieldX(i,j,k) = -loc_EwaldFieldX(i,j,k) + x / (r+1e-99)**3
            loc_EwaldFieldY(i,j,k) = -loc_EwaldFieldY(i,j,k) + y / (r+1e-99)**3
            loc_EwaldFieldZ(i,j,k) = -loc_EwaldFieldZ(i,j,k) + z / (r+1e-99)**3

          endif ! MyPE chunck

        enddo ! i
      enddo ! j
    enddo ! k

    if (MyPE == MASTER_PE) write (*,*) 'pt_sinkPrepareEwald: Now reducing Ewald field...'

    ! distribute chunks of the Ewald field to all CPUs
    call MPI_AllReduce(loc_EwaldFieldX, sink_EwaldFieldX, nx*ny*nz, &
                       FLASH_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_AllReduce(loc_EwaldFieldY, sink_EwaldFieldY, nx*ny*nz, &
                       FLASH_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_AllReduce(loc_EwaldFieldZ, sink_EwaldFieldZ, nx*ny*nz, &
                       FLASH_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)

    ! write the Ewald fields to file
    if (MyPE == MASTER_PE) then

      write (strbuff, '("Ewald correction field calculated.")')
      call Logfile_stamp(strbuff, "[Sink particles]")
      write (*,*) 'pt_sinkPrepareEwald: '//strbuff

      if (file_exists) then
         openflag = 'replace'
      else
         openflag = 'new'
      endif
      funit = ut_getFreeFileUnit ()
      open(unit=funit, file=trim(sink_EwaldFileName), status=openflag)
      write(funit,'(6(a24, 1x))'), '|dx|', '|dy|', '|dz|', &
                                   'force correction x', 'force correction y', 'force correction z'
      do k = 0, nz-1
        z = k * ewald_zmax / real(nz-1)
        do j = 0, ny-1
          y = j * ewald_ymax / real(ny-1)
          do i = 0, nx-1
            x = i * ewald_xmax / real(nx-1)
            write(funit,'(6(e24.17, 1x))') x, y, z, &
                   sink_EwaldFieldX(i,j,k), sink_EwaldFieldY(i,j,k), sink_EwaldFieldZ(i,j,k) 
          enddo
          write(funit,*) ! blank line
        enddo
        write(funit,*) ! blank line
      enddo
      close(unit=funit)

    endif ! MASTER_PE

    deallocate(loc_EwaldFieldX)
    deallocate(loc_EwaldFieldY)
    deallocate(loc_EwaldFieldZ)

  endif ! generate the ewald field

  return

end subroutine pt_sinkPrepareEwald
