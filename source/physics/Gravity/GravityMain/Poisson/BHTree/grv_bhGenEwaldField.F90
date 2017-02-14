!!****if* source/physics/Gravity/GravityMain/Poisson/BHTree/grv_bhGenEwaldField
!!
!! NAME
!!
!!  grv_bhGenEwaldField
!!
!!
!! SYNOPSIS
!!
!!   call grv_bhGenEwaldField()
!!
!! DESCRIPTION
!!
!!  Generates the Ewald field. It is done in parallel, the whole field is
!!  distributed among all cpus. The field is eventually saved. 
!!  In case grv_bhEwaldAlwaysGenerate is FALSE and
!!  the file with the Ewald field exists, the Ewald field is read from it.
!!
!!  There are two approaches how the Ewald field can be generated. The default one
!!  utilizes subtracting 1/r term which assures only one level of the EF is 
!!  sufficient. In order to gain some speedup during interpolation, the 
!!  EF is generated with its partial derivatives and then it is reorganised.
!!  Generating the EF is done by grv_bhGenEwaldFieldLevel which is called by this subroutine.
!!  Since it is relatively new and serious problems cannot be completely ruled out, 
!!  we retained the possibility to switch off the older approach. It can
!!  be done by macro bhtreeEwaldV42. In this case the EF is generated 
!!  on a nested grid (with the highest density at small distances), one level 
!!  of the EF grid is calculated by grv_bhGenEwaldFieldLevelV42 which is called 
!!  by this subroutine. Note that the EF computed by grv_bhGenEwaldFieldLevel is
!!  completely different from the EF computed by grv_bhGenEwaldFieldLevelV42.
!!
!! ARGUMENTS
!!
!!   none
!!
!!***




subroutine grv_bhGenEwaldField()
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Logfile_interface, ONLY : Logfile_stamp
  use Gravity_data, ONLY : grv_bhEwaldFieldNxV42, grv_bhEwaldFieldNyV42, grv_bhEwaldFieldNzV42, &
    grv_bhLx, grv_bhLy, grv_bhLz, grv_bhEwaldLMaxV42, &
    grv_bhDxI, grv_bhDyIV42, grv_bhDzIV42, grv_bhMinEFSizeV42, &
    grv_meshMe, grv_meshComm, grv_bhTreeEwaldAccV42, grv_bhTreeEwaldPotV42, &
    grv_bhEwaldSeriesN, grv_bhEwaldNRefV42, grv_bhDirectionQuadV42, &
    grv_bhGenEwaldAccV42, grv_bhGenEwaldPotV42, &
    grv_bhLogfield, grv_bhLogfieldDer, grv_bhDLogI, grv_bhDLog, &
    grv_bhPotConst, grv_bhLayerPotConst, grv_bhLayerAccConst, grv_bhTreeEwald, &
    grv_bhEwald_periodicity, grv_margin, grv_L2inv
  use grv_bhInterface, ONLY : grv_bhGenEwaldFieldLevelV42, grv_bhGenEwaldFieldLevel
  use Grid_interface, ONLY : Grid_getMinCellSizes
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none
#include "Flash.h"
#include "Flash_mpi.h"
#include "constants.h"

  real, parameter :: pi = PI
  integer, parameter :: nlogbins=2000
  integer :: i, j, k, l, m, ierr, nref
  integer :: grv_bhEwaldNPer
  integer :: jj, kk, nlogbin_min, nlogbin_max
  integer :: n_ewaldGrid(0:2)
  real :: xlog, r2min, r2max
  real :: xmax, xmin, ymax, ymin, zmax, zmin, dum
  real :: x, y, z, ewald_xmax, ewald_ymax, ewald_zmax
  real, dimension(MDIM) :: mcs
  real, allocatable :: field_prep_accV42(:,:,:,:)
  real, allocatable :: field_prep_potV42(:,:,:)
  real, allocatable :: field_ewald_prep(:,:,:,:)
  real :: l_ewaldGrid(0:2)
  logical :: file_exists, file_acc_exists, file_pot_exists
  logical :: grv_bhLinearInterpolOnlyV42, grv_bhEwaldAlwaysGenerate
  real :: grv_bhDX, L1inv

  character(len=MAX_STRING_LENGTH) :: grv_bhEwaldFNameAccV42, grv_bhEwaldFNamePotV42
  character(len=MAX_STRING_LENGTH) :: strBuff, grav_boundary_type_x &
  & , grav_boundary_type_y, grav_boundary_type_z, grv_bhEwaldFName
  character(len=8) :: grv_char8

  call RuntimeParameters_get("grav_boundary_type_x", grav_boundary_type_x)
  call RuntimeParameters_get("grav_boundary_type_y", grav_boundary_type_y)
  call RuntimeParameters_get("grav_boundary_type_z", grav_boundary_type_z)

  ! return if periodic boundaries are not used in any direction
  if (    (grav_boundary_type_x .ne. "periodic") &
  & .and. (grav_boundary_type_y .ne. "periodic") &
  & .and. (grav_boundary_type_z .ne. "periodic")) return

  ! read runtime parameters
  call RuntimeParameters_get("xmin", xmin)
  call RuntimeParameters_get("xmax", xmax)
  call RuntimeParameters_get("ymin", ymin)
  call RuntimeParameters_get("ymax", ymax)
  call RuntimeParameters_get("zmin", zmin)
  call RuntimeParameters_get("zmax", zmax)
  call RuntimeParameters_get("grv_bhEwaldSeriesN", grv_bhEwaldSeriesN)

  call RuntimeParameters_get("grv_bhEwaldAlwaysGenerate", grv_bhEwaldAlwaysGenerate)

  ! periodicity and orientation of the problem
  grv_bhEwald_periodicity = 0
  if (grav_boundary_type_x == "periodic") grv_bhEwald_periodicity = grv_bhEwald_periodicity + 1
  if (grav_boundary_type_y == "periodic") grv_bhEwald_periodicity = grv_bhEwald_periodicity + 2
  if (grav_boundary_type_z == "periodic") grv_bhEwald_periodicity = grv_bhEwald_periodicity + 4

  ! dimensions of the computational domain (needed also in the new approach in
  ! routine grv_bhGenEwaldFieldLevel.F90)
  grv_bhLx = xmax - xmin
  grv_bhLy = ymax - ymin
  grv_bhLz = zmax - zmin

! old implementation of the Ewald method

#ifdef GRAV_TREE_EWALD_V42
! read runtime parameters
  call RuntimeParameters_get("grv_bhEwaldFieldNxV42", grv_bhEwaldFieldNxV42)
  call RuntimeParameters_get("grv_bhEwaldFieldNyV42", grv_bhEwaldFieldNyV42)
  call RuntimeParameters_get("grv_bhEwaldFieldNzV42", grv_bhEwaldFieldNzV42)
  call RuntimeParameters_get("grv_bhEwaldNRefV42", grv_bhEwaldNRefV42)
  call RuntimeParameters_get("grv_bhEwaldFNameAccV42", grv_bhEwaldFNameAccV42)
  call RuntimeParameters_get("grv_bhEwaldFNamePotV42", grv_bhEwaldFNamePotV42)
  call RuntimeParameters_get("grv_bhLinearInterpolOnlyV42", grv_bhLinearInterpolOnlyV42)

  ! choose between linear and quadratic (quadratic only in one direction) 
  ! interpolation method (for routine Gravity_bhEwald)
  if (grv_bhLinearInterpolOnlyV42) then
    grv_bhDirectionQuadV42 = 0
  else
    if (grv_bhEwald_periodicity == 6) then
      grv_bhDirectionQuadV42 = 1
    else if (grv_bhEwald_periodicity == 5) then
      grv_bhDirectionQuadV42 = 2
    else if (grv_bhEwald_periodicity == 3) then
      grv_bhDirectionQuadV42 = 3
    else
      grv_bhDirectionQuadV42 = 0
    endif
  endif

! Both fields are calculated
  grv_bhGenEwaldAccV42=.True.
  grv_bhGenEwaldPotV42=.True.

  ! dimensions of one grid-cell of the largest Ewald field,
  ! used in grv_bhEwald.F90
  grv_bhEwaldLMaxV42 = max(grv_bhLx,grv_bhLy,grv_bhLz)
  grv_bhDxI = (grv_bhEwaldFieldNxV42-1) / grv_bhEwaldLMaxV42
  grv_bhDyIV42 = (grv_bhEwaldFieldNyV42-1) / grv_bhEwaldLMaxV42
  grv_bhDzIV42 = (grv_bhEwaldFieldNzV42-1) / grv_bhEwaldLMaxV42

  ! number of refinement levels of the Ewald field
  if (grv_bhEwaldNRefV42 <= 0) then
    call Grid_getMinCellSizes(mcs)
    grv_bhEwaldNRefV42 = 1
    do
      !if (grv_meshMe == MASTER_PE)print *, "NRef: ", grv_bhEwaldNRefV42, mcs, &
      !& "|", 1/grv_bhDxI, 1/grv_bhDyIV42, 1/grv_bhDzIV42
      grv_bhMinEFSizeV42 = grv_bhEwaldLMaxV42/(2**(grv_bhEwaldNRefV42-1))
      if (    (grv_bhMinEFSizeV42 < 0.5*mcs(IAXIS)) &
      & .and. (grv_bhMinEFSizeV42 < 0.5*mcs(JAXIS)) &
      & .and. (grv_bhMinEFSizeV42 < 0.5*mcs(KAXIS))) exit
      grv_bhEwaldNRefV42 = grv_bhEwaldNRefV42 + 1
    enddo
  else
    grv_bhMinEFSizeV42 = grv_bhEwaldLMaxV42/(2**(grv_bhEwaldNRefV42-1))
  endif

  write (strBuff, '("Ewald refinement level determined: ",i3)') grv_bhEwaldNRefV42 
  call Logfile_stamp( strBuff, "[BHTree]")

  ! allocate array for the Ewald field
  if (grv_bhGenEwaldAccV42) then
    allocate(grv_bhTreeEwaldAccV42(0:grv_bhEwaldNRefV42-1, IAXIS:KAXIS,  -1:grv_bhEwaldFieldNxV42, &
    &                        -1:grv_bhEwaldFieldNyV42, -1:grv_bhEwaldFieldNzV42), stat=ierr)
    if (ierr /= 0) call Driver_abortFlash("could not allocate grv_bhTreeEwaldAccV42 in grv_bhGenEwaldField.F90")
  endif

  if (grv_bhGenEwaldPotV42) then
    allocate(grv_bhTreeEwaldPotV42(0:grv_bhEwaldNRefV42-1, -1:grv_bhEwaldFieldNxV42, &
    &                        -1:grv_bhEwaldFieldNyV42, -1:grv_bhEwaldFieldNzV42), stat=ierr)
    if (ierr /= 0) call Driver_abortFlash("could not allocate grv_bhTreeEwaldPotV42 in grv_bhGenEwaldField.F90")
  endif

  ! check the existence of the ewald_field files with correct grv_bhEwaldNRefV42
  if (grv_meshMe == MASTER_PE) then
    ! check necessity and existence of the ewald_field file for acceleration
    if (grv_bhGenEwaldAccV42) then
      inquire(file=grv_bhEwaldFNameAccV42, exist=file_acc_exists)
      if (file_acc_exists .and. (.not. grv_bhEwaldAlwaysGenerate)) then
        open(unit=53, file=grv_bhEwaldFNameAccV42, status='old')
        read(53,*) nref
        if (nref /= grv_bhEwaldNRefV42) then
          file_acc_exists = .false.
          close(unit=53)
        endif
      else
        file_acc_exists = .false.
      endif
    else
      file_acc_exists = .false.
    endif
    
    ! check necessity and existence of the ewald_field file for potential
    if (grv_bhGenEwaldPotV42) then
      inquire(file=grv_bhEwaldFNamePotV42, exist=file_pot_exists)
      if (file_pot_exists .and. (.not. grv_bhEwaldAlwaysGenerate)) then
        open(unit=53, file=grv_bhEwaldFNamePotV42, status='old')
        read(53,*) nref
        if (nref /= grv_bhEwaldNRefV42) then
          file_pot_exists = .false.
          close(unit=53)
        endif
      else
        file_pot_exists = .false.
      endif
    else
      file_pot_exists = .false.
    endif

    ! check existence of necessary files
    file_exists = .true.
    if (grv_bhGenEwaldAccV42) then
      if (.not. file_acc_exists) then
        file_exists = .false.
      endif
    endif

    if (grv_bhGenEwaldPotV42) then
      if (.not. file_pot_exists) then
        file_exists = .false.
      endif
    endif
  endif

  call MPI_Bcast(file_exists, 1, FLASH_LOGICAL, MASTER_PE, grv_meshComm, ierr)
  call MPI_Bcast(file_acc_exists, 1, FLASH_LOGICAL, MASTER_PE, grv_meshComm, ierr)
  call MPI_Bcast(file_pot_exists, 1, FLASH_LOGICAL, MASTER_PE, grv_meshComm, ierr)

  if (file_exists) then
    ! ewald_field file with correct grv_bhEwaldNRefV42 is present and there is 
    ! no user request to regenerate it => read the ewald field from the file
    if (file_acc_exists) then
      write (strBuff, '("Reading Ewald field from file: ",a20)') grv_bhEwaldFNameAccV42
      call Logfile_stamp( strBuff, "[BHTree]")
      if (grv_meshMe == MASTER_PE) then
        do m = IAXIS, KAXIS
          do l = 0, grv_bhEwaldNRefV42-1
            do k = -1,grv_bhEwaldFieldNzV42
              do j = -1,grv_bhEwaldFieldNyV42
                do i = -1,grv_bhEwaldFieldNxV42
                  read(53,*) x, y, z, grv_bhTreeEwaldAccV42(l,m,i,j,k), dum
                enddo
                read(53,*)
              enddo
              read(53,*)
            enddo
            read(53,*)
          enddo
          read(53,*)
        enddo
      endif
    endif

    if (file_pot_exists) then
      write (strBuff, '("Reading Ewald field from file: ",a20)') grv_bhEwaldFNamePotV42
      call Logfile_stamp( strBuff, "[BHTree]")
      if (grv_meshMe == MASTER_PE) then
        do l = 0, grv_bhEwaldNRefV42-1
          do k = -1,grv_bhEwaldFieldNzV42
            do j = -1,grv_bhEwaldFieldNyV42
              do i = -1,grv_bhEwaldFieldNxV42
                read(53,*) x, y, z, grv_bhTreeEwaldPotV42(l,i,j,k), dum
              enddo
              read(53,*)
            enddo
            read(53,*)
          enddo
          read(53,*)
        enddo
      close(unit=53)
      endif
    endif

    if (file_acc_exists) then
      call MPI_Bcast(grv_bhTreeEwaldAccV42, 3*grv_bhEwaldFieldNxV42*grv_bhEwaldFieldNyV42*grv_bhEwaldFieldNzV42 &
      &  , FLASH_REAL, MASTER_PE, grv_meshComm, ierr)
    endif

    if (file_pot_exists) then
      call MPI_Bcast(grv_bhTreeEwaldPotV42, grv_bhEwaldFieldNxV42*grv_bhEwaldFieldNyV42*grv_bhEwaldFieldNzV42 &
      &  , FLASH_REAL, MASTER_PE, grv_meshComm, ierr)
    endif
    write (strBuff, '("Ewald correction field read and communicated.")')
    call Logfile_stamp( strBuff, "[BHTree]")

  else

    ! Ewald field will be generated
    write (strBuff, '("Generating Ewald correction field")')
    call Logfile_stamp( strBuff, "[BHTree]")
    if (grv_bhGenEwaldAccV42) then
      allocate(field_prep_accV42(IAXIS:KAXIS, -1:grv_bhEwaldFieldNxV42, -1:grv_bhEwaldFieldNyV42, &
      & -1:grv_bhEwaldFieldNzV42), stat=ierr)
    endif
    if (grv_bhGenEwaldPotV42) then
      allocate(field_prep_potV42(-1:grv_bhEwaldFieldNxV42, -1:grv_bhEwaldFieldNyV42, &
      & -1:grv_bhEwaldFieldNzV42), stat=ierr)
    endif

    if (ierr /= 0) call &
    &  Driver_abortFlash("could not allocate grv_bhTreeEwald in grv_bhGenEwaldField.F90")
    do l = 0, grv_bhEwaldNRefV42-1
      ewald_xmax = grv_bhEwaldLMaxV42/(2**l)
      ewald_ymax = grv_bhEwaldLMaxV42/(2**l)
      ewald_zmax = grv_bhEwaldLMaxV42/(2**l)

      ! call subroutine which generate ewald field with appropriate boundary conditions
      call grv_bhGenEwaldFieldLevelV42(grv_bhEwaldFieldNxV42, grv_bhEwaldFieldNyV42, &
           & grv_bhEwaldFieldNzV42, ewald_xmax, ewald_ymax, ewald_zmax, &
           & grv_bhEwald_periodicity, field_prep_potV42, field_prep_accV42)

      if (grv_bhGenEwaldAccV42) then
        grv_bhTreeEwaldAccV42(l,:,:,:,:) = field_prep_accV42(:,:,:,:)
      endif
      if (grv_bhGenEwaldPotV42) then
        grv_bhTreeEwaldPotV42(l,:,:,:) = field_prep_potV42(:,:,:)
      endif

    enddo
    if (grv_bhGenEwaldAccV42) then
      deallocate(field_prep_accV42)
    endif
    if (grv_bhGenEwaldPotV42) then
      deallocate(field_prep_potV42)
    endif

    ! write the ewald field into the file
    if (grv_meshMe == MASTER_PE) then
      if (grv_bhGenEwaldAccV42) then
      open(unit=53, file=grv_bhEwaldFNameAccV42, status='replace')
      write(53,*) grv_bhEwaldNRefV42

      do m = IAXIS, KAXIS
        write(53,*) 'component ', m
        do l = 0, grv_bhEwaldNRefV42-1
          ewald_xmax = grv_bhEwaldLMaxV42/(2**l)
          ewald_ymax = grv_bhEwaldLMaxV42/(2**l)
          ewald_zmax = grv_bhEwaldLMaxV42/(2**l)
          write(53,*) 'refinement level ', l
          do k = -1,grv_bhEwaldFieldNzV42
            z = (k * ewald_zmax) / (grv_bhEwaldFieldNzV42-1)
            do j = -1,grv_bhEwaldFieldNyV42
              y = (j *  ewald_ymax) / (grv_bhEwaldFieldNyV42-1)
              do i = -1,grv_bhEwaldFieldNxV42
                x = (i * ewald_xmax) / (grv_bhEwaldFieldNxV42-1)
                write(53,'(5(e17.10, 2x))') x, y, z, grv_bhTreeEwaldAccV42(l,m,i,j,k) &
                , 1.0/(sqrt(x*x+y*y+z*z)+1d-99)
              enddo
              write(53,*)
            enddo
            write(53,*)
          enddo
          write(53,*)
        enddo
        write(53,*)
      enddo
      close(unit=53)
      write (strBuff, '("Ewald correction field written into file:", a20)') grv_bhEwaldFNameAccV42
      call Logfile_stamp( strBuff, "[BHTree]")
      endif ! grv_bhGenEwaldAccV42

      if (grv_bhGenEwaldPotV42) then
      open(unit=53, file=grv_bhEwaldFNamePotV42, status='replace')
      write(53,*) grv_bhEwaldNRefV42

      do l = 0, grv_bhEwaldNRefV42-1
        ewald_xmax = grv_bhEwaldLMaxV42/(2**l)
        ewald_ymax = grv_bhEwaldLMaxV42/(2**l)
        ewald_zmax = grv_bhEwaldLMaxV42/(2**l)
        write(53,*) 'refinement level ', l
        do k = -1,grv_bhEwaldFieldNzV42
          z = (k * ewald_zmax) / (grv_bhEwaldFieldNzV42-1)
          do j = -1,grv_bhEwaldFieldNyV42
            y = (j *  ewald_ymax) / (grv_bhEwaldFieldNyV42-1)
            do i = -1,grv_bhEwaldFieldNxV42
              x = (i * ewald_xmax) / (grv_bhEwaldFieldNxV42-1)
              write(53,'(5(e17.10, 2x))') x, y, z, grv_bhTreeEwaldPotV42(l,i,j,k) &
              , 1.0/(sqrt(x*x+y*y+z*z)+1d-99)
            enddo
            write(53,*)
          enddo
          write(53,*)
        enddo
        write(53,*)
      enddo
      close(unit=53)
      write (strBuff, '("Ewald correction field written into file:", a20)') grv_bhEwaldFNamePotV42
      call Logfile_stamp( strBuff, "[BHTree]")
      endif ! grv_bhGenEwaldPotV42
    endif ! MASTER_PE

  endif ! ewald field generated

#else

! new implementation of the Ewald method

! get variables needed for this approach
  call RuntimeParameters_get("grv_bhEwaldFName", grv_bhEwaldFName)
  call RuntimeParameters_get("grv_bhEwaldNPer", grv_bhEwaldNPer)

!  grv_margin = 4*grv_bhEwaldNPer-1
  grv_margin = 4*grv_bhEwaldNPer

! compute basic parameters of desired Ewald field
  select case (grv_bhEwald_periodicity)
    case (1)
      l_ewaldGrid(0) = 0.5D0*grv_bhLx
      l_ewaldGrid(1) = min(4*l_ewaldGrid(0), grv_bhLy)
      l_ewaldGrid(2) = min(4*l_ewaldGrid(0), grv_bhLz)

      n_ewaldGrid(0) = grv_bhEwaldNPer
      n_ewaldGrid(1) = min(grv_margin,floor(0.5D0 + grv_bhEwaldNPer*l_ewaldGrid(1)/l_ewaldGrid(0)))
      n_ewaldGrid(2) = min(grv_margin,floor(0.5D0 + grv_bhEwaldNPer*l_ewaldGrid(2)/l_ewaldGrid(0)))

      grv_bhDxI = 2*grv_bhEwaldNPer / grv_bhLx

      L1inv = 1.0D0/grv_bhLx
      grv_L2inv = 2.0D0/grv_bhLx

      r2min = (2*grv_bhLx)**2
      r2max = (max(grv_bhLx,grv_bhLy,grv_bhLz))**2

    case (2)
      l_ewaldGrid(1) = 0.5D0*grv_bhLy
      l_ewaldGrid(0) = min(4*l_ewaldGrid(1), grv_bhLx)
      l_ewaldGrid(2) = min(4*l_ewaldGrid(1), grv_bhLz)

      n_ewaldGrid(0) = min(grv_margin,floor(0.5D0 + grv_bhEwaldNPer*l_ewaldGrid(0)/l_ewaldGrid(1)))
      n_ewaldGrid(1) = grv_bhEwaldNPer
      n_ewaldGrid(2) = min(grv_margin,floor(0.5D0 + grv_bhEwaldNPer*l_ewaldGrid(2)/l_ewaldGrid(1)))

      grv_bhDxI = 2*grv_bhEwaldNPer / grv_bhLy

      L1inv = 1.0D0/grv_bhLy
      grv_L2inv = 2.0D0/grv_bhLy

      r2min = (2*grv_bhLy)**2
      r2max = (max(grv_bhLx,grv_bhLy,grv_bhLz))**2

    case (4)
      l_ewaldGrid(2) = 0.5D0*grv_bhLz
      l_ewaldGrid(0) = min(4*l_ewaldGrid(2), grv_bhLx)
      l_ewaldGrid(1) = min(4*l_ewaldGrid(2), grv_bhLy)

      n_ewaldGrid(0) = min(grv_margin,floor(0.5D0 + grv_bhEwaldNPer*l_ewaldGrid(0)/l_ewaldGrid(2)))
      n_ewaldGrid(1) = min(grv_margin,floor(0.5D0 + grv_bhEwaldNPer*l_ewaldGrid(1)/l_ewaldGrid(2)))
      n_ewaldGrid(2) = grv_bhEwaldNPer

      grv_bhDxI = 2*grv_bhEwaldNPer / grv_bhLz

      L1inv = 1.0D0/grv_bhLz
      grv_L2inv = 2.0D0/grv_bhLz

      r2min = (2*grv_bhLz)**2
      r2max = (max(grv_bhLx,grv_bhLy,grv_bhLz))**2

    case (3)
      n_ewaldGrid(0) = grv_bhEwaldNPer
      n_ewaldGrid(1) = floor(0.5D0+grv_bhEwaldNPer*(grv_bhLy/grv_bhLx))
      n_ewaldGrid(2) = min(grv_margin,floor(0.5D0 + 2*grv_bhEwaldNPer*grv_bhLz/grv_bhLx))

      grv_bhDxI = 2*grv_bhEwaldNPer / grv_bhLx

! acceleration in the infinity
      grv_bhLayerAccConst = 2*pi/(grv_bhLx*grv_bhLy)
      grv_bhLayerPotConst = -grv_bhLayerAccConst

    case (5)
      n_ewaldGrid(0) = floor(0.5D0+grv_bhEwaldNPer*(grv_bhLx/grv_bhLz))
      n_ewaldGrid(1) = min(grv_margin,floor(0.5D0 + 2*grv_bhEwaldNPer*grv_bhLy/grv_bhLz))
      n_ewaldGrid(2) = grv_bhEwaldNPer

      grv_bhDxI = 2*grv_bhEwaldNPer / grv_bhLz

! acceleration in the infinity
      grv_bhLayerAccConst = 2*pi/(grv_bhLx*grv_bhLz)
      grv_bhLayerPotConst = -grv_bhLayerAccConst

    case (6)
      n_ewaldGrid(0) = min(grv_margin,floor(0.5D0 + 2*grv_bhEwaldNPer*grv_bhLx/grv_bhLy))
      n_ewaldGrid(1) = grv_bhEwaldNPer
      n_ewaldGrid(2) = floor(0.5D0+grv_bhEwaldNPer*(grv_bhLz/grv_bhLy))

      grv_bhDxI = 2*grv_bhEwaldNPer / grv_bhLy

! acceleration in the infinity
      grv_bhLayerAccConst = 2*pi/(grv_bhLy*grv_bhLz)
      grv_bhLayerPotConst = -grv_bhLayerAccConst

    case (7)
      n_ewaldGrid(0) = grv_bhEwaldNPer
      n_ewaldGrid(1) = floor(0.5D0 + grv_bhEwaldNPer*(grv_bhLy/grv_bhLx))
      n_ewaldGrid(2) = floor(0.5D0 + grv_bhEwaldNPer*(grv_bhLz/grv_bhLx))

      grv_bhDxI = 2*grv_bhEwaldNPer / grv_bhLx

! reset grv_margin to avoid problem in Gravity_bhEwaldAcc or Gravity_bhEwaldPot routines in the case 
! of computational domain more elongated than by factor 4
      grv_margin = maxval(n_ewaldGrid)

  end select

! for cylindrical symmetry constructs arrays estimating logarithm
  if ((grv_bhEwald_periodicity == 1).or.(grv_bhEwald_periodicity == 2).or.(grv_bhEwald_periodicity == 4)) then

    grv_bhDLog = (r2max-r2min)/nlogbins
    grv_bhDLogI = 1.0D0/grv_bhDLog
    nlogbin_min = floor(grv_bhDLogI*r2min)
    nlogbin_max = ceiling(grv_bhDLogI*r2max)

    allocate(grv_bhLogfield(nlogbin_min:nlogbin_max), stat=ierr)
    allocate(grv_bhLogfieldDer(nlogbin_min:nlogbin_max), stat=ierr)

! fill the arrays in appropriate extent
    do i = 0, nlogbin_max
      if (i.lt.nlogbin_min) cycle
      xlog = i*grv_bhDLog 
      grv_bhLogfield(i) = -L1inv*log(xlog)
      grv_bhLogfieldDer(i) = -L1inv/xlog
    enddo

  endif

! allocate array for the field minus 1/r term
  allocate(grv_bhTreeEwald(0:12, 0:n_ewaldGrid(0), 0:n_ewaldGrid(1), 0:n_ewaldGrid(2)), stat=ierr)
  if (ierr /= 0) call Driver_abortFlash("could not allocate grv_bhTreeEwald in grv_bhGenEwaldField.F90")

! read the coefficients from file
  if ((.not.grv_bhEwaldAlwaysGenerate)) then
    if (grv_meshMe == MASTER_PE) then
      inquire(file=grv_bhEwaldFName, exist=file_exists)
      if (file_exists) then

      open(unit=53, file=grv_bhEwaldFName, status='old')
        do i=0, n_ewaldGrid(0)
          do j=0, n_ewaldGrid(1)
           do k=0, n_ewaldGrid(2)

             read(53,10) grv_char8, x, y, z

             read(53,11)   grv_bhTreeEwald(0,i,j,k),grv_bhTreeEwald(1,i,j,k), &
                           &  grv_bhTreeEwald(2,i,j,k),grv_bhTreeEwald(3,i,j,k)

             read(53,12)   grv_bhTreeEwald(4,i,j,k),grv_bhTreeEwald(5,i,j,k), &
                           &  grv_bhTreeEwald(6,i,j,k)

             read(53,13)   grv_bhTreeEwald(7,i,j,k),grv_bhTreeEwald(8,i,j,k),grv_bhTreeEwald(9,i,j,k), &
                           &  grv_bhTreeEwald(10,i,j,k),grv_bhTreeEwald(11,i,j,k), grv_bhTreeEwald(12,i,j,k)

             read(53,*)
           enddo
         enddo
       enddo

      endif
      close(unit=53)
    endif

    write (strBuff, '("Coefficients in the Taylor expansion read from file")')
    call Logfile_stamp( strBuff, "[BHTree]")

    call MPI_Bcast(file_exists, 1 &
    &  , FLASH_LOGICAL, MASTER_PE, grv_meshComm, ierr)

! communicate the arrays
  if (file_exists) then
      call MPI_Bcast(grv_bhTreeEwald, 13*(grv_bhEwaldNPer+1)*(grv_bhEwaldNPer+1)*(grv_margin+1) &
      &  , FLASH_REAL, MASTER_PE, grv_meshComm, ierr)
    endif

! compute the coefficients
  else

    allocate(field_ewald_prep(0:12, 0:n_ewaldGrid(0), 0:n_ewaldGrid(1), 0:n_ewaldGrid(2)), stat=ierr)

      call grv_bhGenEwaldFieldLevel(n_ewaldGrid(0), n_ewaldGrid(1), &
           & n_ewaldGrid(2), field_ewald_prep)

      grv_bhTreeEwald(:,:,:,:) = field_ewald_prep(:,:,:,:)

    deallocate(field_ewald_prep)

! write the coefficients into file
  if (grv_meshMe == MASTER_PE) then
    if (grv_bhEwaldAlwaysGenerate) then
      open(unit=53, file=grv_bhEwaldFName, status='replace')
!    else
!      open(unit=53, file='ewald_coefficients_check', status='replace')
    endif

  grv_bhDX = 1.0D0/grv_bhDxI

  do i=0, n_ewaldGrid(0)
    do j=0, n_ewaldGrid(1)
      do k=0, n_ewaldGrid(2)

        x = i * grv_bhDX
        y = j * grv_bhDX
        z = k * grv_bhDX

        write(53,10) 'Point at', x, y, z

        write(53,11)   grv_bhTreeEwald(0,i,j,k),grv_bhTreeEwald(1,i,j,k), &
                      &  grv_bhTreeEwald(2,i,j,k),grv_bhTreeEwald(3,i,j,k)

        write(53,12)   grv_bhTreeEwald(4,i,j,k),grv_bhTreeEwald(5,i,j,k), &
                      &  grv_bhTreeEwald(6,i,j,k)

        write(53,13)   grv_bhTreeEwald(7,i,j,k),grv_bhTreeEwald(8,i,j,k),grv_bhTreeEwald(9,i,j,k), &
                      &  grv_bhTreeEwald(10,i,j,k),grv_bhTreeEwald(11,i,j,k), grv_bhTreeEwald(12,i,j,k)

        write(53,*)
      enddo
    enddo
  enddo
    close(unit=53)

    write (strBuff, '("Coefficients in the Taylor expansion written into file")')
    call Logfile_stamp( strBuff, "[BHTree]")

  endif

10 FORMAT(A,3ES16.8)
11 FORMAT(4ES16.8)
12 FORMAT(3ES16.8)
13 FORMAT(6ES16.8)

  endif ! if ((.not.grv_bhEwaldAlwaysGenerate)) then

#endif

  return
end subroutine grv_bhGenEwaldField

