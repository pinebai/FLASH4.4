!!****if* source/Simulation/SimulationMain/StirFromFile/IO_writeIntegralQuantities
!!
!!  NAME
!!    IO_writeIntegralQuantities
!!
!!  SYNOPSIS
!!    call IO_writeIntegralQuantities(integer(in) :: isFirst,
!!                                    real(in)    :: simTime)
!!
!!  DESCRIPTION
!!   Compute the values of integral quantities (eg. total energy)
!!   and write them to an ASCII file.  If this is the initial step,
!!   create the file and write a header to it before writing the data.
!!
!!   Presently, this supports 1, 2, and 3-d Cartesian geometry and 2-d
!!   cylindrical geometry (r,z).  More geometries can be added by
!!   modifying the volume of each zone (dvol).
!!
!!   Users should modify this routine if they want to store any
!!   quantities other than default values in the flash.dat.  Make sure
!!   to modify the nGlobalSum parameter to match the number of
!!   quantities written.  Also make sure to modify the header to match
!!   the names of quantities with those calculated in the lsum and
!!   gsum arrays.
!!
!!  ARGUMENTS
!!   isFirst - if 1 then write header info plus data, otherwise just write data
!!   simTime - simulation time
!!
!!  AUTHOR: Christoph Federrath, 2008-2013
!!
!!***

subroutine IO_writeIntegralQuantities (isFirst, simTime)

  use IO_data       , ONLY : io_globalMe, io_restart, io_statsFileName, io_globalComm
  use Grid_interface, ONLY : Grid_computeUserVars, Grid_getListOfBlocks, Grid_getBlkIndexLimits, &
                             Grid_getDeltas, Grid_getBlkPtr, Grid_releaseBlkPtr

  implicit none

#include "mpif.h"
#include "constants.h"
#include "Flash.h"

  real, intent(in)    :: simTime
  integer, intent(in) :: isFirst

  integer :: lb, count

  integer :: funit = 99
  integer :: error, ioStat

  integer :: blockList(MAXBLOCKS)
  integer :: blkLimits(HIGH, MDIM), blkLimitsGC(HIGH, MDIM)

  integer, parameter :: nGlobalSum = 39     ! number of globally-summed quantities
  real               :: lsum(0:nGlobalSum)  ! local summed quantities
  real               :: gsum(0:nGlobalSum)  ! global summed quantities
  integer, parameter :: vortnum = 18, magnum = 30

  integer :: i, j, k
  real    :: dvol, del(MDIM)
  real, DIMENSION(:,:,:,:), POINTER :: solnData

  real :: rms_Mach, rms_Mach_mw, rms_Mach_netto, rms_Mach_netto_mw, rms_Forcing, mean_temp, mean_pres
  real :: lmin_dens = 1E99, lmax_dens, gmin_dens, gmax_dens, mean_dens, rms_dens
  real :: lmax_alfvenspeed, gmax_alfvenspeed, alfven_speed
  real :: sigma_dens, mean_ln_dens, rms_ln_dens, sigma_ln_dens
  logical, parameter :: Debug = .false.


  if (Debug .and. (io_globalMe == MASTER_PE)) print *, 'IO_writeIntegralQuantities:  entering ...'

  ! Make sure the vorticity is up-to-date
  call Grid_computeUserVars()

  ! Sum quantities over all locally held leaf-node blocks.
  gsum(:) = 0.0
  lsum(:) = 0.0

  call Grid_getListOfBlocks(LEAF, blockList, count)

  do lb = 1, count

     !get the index limits of the block
     call Grid_getBlkIndexLimits(blockList(lb), blkLimits, blkLimitsGC)

     ! get a pointer to the current block of data
     call Grid_getBlkPtr(blockList(lb), solnData)

     !get the deltas of this block and compute the cell volume
     call Grid_getDeltas(blockList(lb), del)

#if NDIM == 1
     dvol = del(1)
#endif
#if NDIM == 2
     dvol = del(1)*del(2)
#endif
#if NDIM == 3
     dvol = del(1)*del(2)*del(3)
#endif

     ! initialize initial local min and max density with
     ! solnData(DENS_VAR,blkLimits(LOW,IAXIS),blkLimits(LOW,JAXIS),blkLimits(LOW,KAXIS))
     ! for subsequent comparison
     lmin_dens = solnData(DENS_VAR,blkLimits(LOW,IAXIS),blkLimits(LOW,JAXIS),blkLimits(LOW,KAXIS))
     lmax_dens = solnData(DENS_VAR,blkLimits(LOW,IAXIS),blkLimits(LOW,JAXIS),blkLimits(LOW,KAXIS))
     gmin_dens = 0.0
     gmax_dens = 0.0
     lmax_alfvenspeed = 0.0
     gmax_alfvenspeed = 0.0

     ! Sum contributions from the indicated range of cells.
     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

              ! total volume
              lsum(0) = lsum(0) + dvol

#ifdef DENS_VAR
              ! mass
              lsum(1) = lsum(1) + solnData(DENS_VAR,i,j,k)*dvol
#endif
#ifdef DENS_VAR
              ! momentum
#ifdef VELX_VAR
              lsum(2) = lsum(2) + solnData(DENS_VAR,i,j,k)*solnData(VELX_VAR,i,j,k)*dvol
#endif
#ifdef VELY_VAR
              lsum(3) = lsum(3) + solnData(DENS_VAR,i,j,k)*solnData(VELY_VAR,i,j,k)*dvol
#endif
#ifdef VELZ_VAR
              lsum(4) = lsum(4) + solnData(DENS_VAR,i,j,k)*solnData(VELZ_VAR,i,j,k)*dvol
#endif
              ! total energy
#ifdef ENER_VAR
              lsum(5) = lsum(5) + solnData(ENER_VAR,i,j,k)*solnData(DENS_VAR,i,j,k)*dvol
#endif
              ! kinetic energy
#if NDIM == 1
              lsum(6) = lsum(6) + 0.5*solnData(DENS_VAR,i,j,k)*(solnData(VELX_VAR,i,j,k)**2)*dvol
#endif
#if NDIM == 2
              lsum(6) = lsum(6) + 0.5*solnData(DENS_VAR,i,j,k)* &
                (solnData(VELX_VAR,i,j,k)**2+solnData(VELY_VAR,i,j,k)**2)*dvol
#endif
#if NDIM == 3
              lsum(6) = lsum(6) + 0.5*solnData(DENS_VAR,i,j,k)* &
                (solnData(VELX_VAR,i,j,k)**2+solnData(VELY_VAR,i,j,k)**2+solnData(VELZ_VAR,i,j,k)**2)*dvol
#endif
#ifdef EINT_VAR
              ! internal energy
              lsum(7) = lsum(7) + solnData(DENS_VAR,i,j,k)*solnData(EINT_VAR,i,j,k)*dvol
#endif
              ! velocity squared
#if NDIM == 1
              lsum(8) = lsum(8) + (solnData(VELX_VAR,i,j,k)**2)*dvol
#endif
#if NDIM == 2
              lsum(8) = lsum(8) + (solnData(VELX_VAR,i,j,k)**2+solnData(VELY_VAR,i,j,k)**2)*dvol
#endif
#if NDIM == 3
              lsum(8) = lsum(8) + (solnData(VELX_VAR,i,j,k)**2+ &
                                   solnData(VELY_VAR,i,j,k)**2+ &
                                   solnData(VELZ_VAR,i,j,k)**2)*dvol
#endif
#ifdef PRES_VAR
#ifdef GAME_VAR
              ! sound speed squared
              if (.true.) then
                lsum(9) = lsum(9) + solnData(PRES_VAR,i,j,k)/solnData(DENS_VAR,i,j,k)*dvol
              else
                lsum(9) = lsum(9) + solnData(GAME_VAR,i,j,k)*solnData(PRES_VAR,i,j,k)/solnData(DENS_VAR,i,j,k)*dvol
              endif
#endif
#endif
              ! velocity squared mass weighted
#if NDIM == 1
              lsum(10) = lsum(10) + (solnData(VELX_VAR,i,j,k)**2)*solnData(DENS_VAR,i,j,k)*dvol
#endif
#if NDIM == 2
              lsum(10) = lsum(10) + (solnData(VELX_VAR,i,j,k)**2+ &
                                     solnData(VELY_VAR,i,j,k)**2)*solnData(DENS_VAR,i,j,k)*dvol
#endif
#if NDIM == 3
              lsum(10) = lsum(10) + (solnData(VELX_VAR,i,j,k)**2+ &
                                     solnData(VELY_VAR,i,j,k)**2+ &
                                     solnData(VELZ_VAR,i,j,k)**2)*solnData(DENS_VAR,i,j,k)*dvol
#endif
#ifdef PRES_VAR
#ifdef GAME_VAR
              ! sound speed squared mass weighted
              if (.true.) then
                lsum(11) = lsum(11) + solnData(PRES_VAR,i,j,k)*dvol
              else
                lsum(11) = lsum(11) + solnData(GAME_VAR,i,j,k)*solnData(PRES_VAR,i,j,k)*dvol
              endif
#endif
#endif
              ! forcing squared

#ifdef ACCX_VAR
#if NDIM == 1
              lsum(12) = lsum(12) + (solnData(ACCX_VAR,i,j,k)**2)*dvol
#endif
#ifdef ACCY_VAR
#if NDIM == 2
              lsum(12) = lsum(12) + (solnData(ACCX_VAR,i,j,k)**2+solnData(ACCY_VAR,i,j,k)**2)*dvol
#endif
#ifdef ACCZ_VAR
#if NDIM == 3
              lsum(12) = lsum(12) + (solnData(ACCX_VAR,i,j,k)**2+ &
                                     solnData(ACCY_VAR,i,j,k)**2+ &
                                     solnData(ACCZ_VAR,i,j,k)**2)*dvol
#endif
#endif
#endif
#endif

#ifdef TEMP_VAR
              ! <temperature>
              lsum(13) = lsum(13) + solnData(TEMP_VAR,i,j,k)*dvol
#endif
#ifdef PRES_VAR
              ! <pressure>
              lsum(14) = lsum(14) + solnData(PRES_VAR,i,j,k)*dvol
#endif
              ! min density
              lmin_dens = min(lmin_dens, solnData(DENS_VAR,i,j,k))

              ! max density
              lmax_dens = max(lmax_dens, solnData(DENS_VAR,i,j,k))

              ! rms density
              lsum(15) = lsum(15) + solnData(DENS_VAR,i,j,k)**2*dvol

              ! < ln density >
              lsum(16) = lsum(16) + alog(solnData(DENS_VAR,i,j,k))*dvol

              ! rms ( ln density )
              lsum(17) = lsum(17) + (alog(solnData(DENS_VAR,i,j,k)))**2*dvol

#endif

#ifdef DENS_VAR
              !! Here are the vorticity, divergence, curl, div  accel (& massweighted)

              !! vorticity and divergence
#ifdef MVRT_VAR
              lsum(vortnum  ) = lsum(vortnum  ) + abs(solnData(MVRT_VAR,i,j,k))*dvol
              lsum(vortnum+1) = lsum(vortnum+1) + abs(solnData(MVRT_VAR,i,j,k))**2*dvol
              lsum(vortnum+2) = lsum(vortnum+2) + abs(solnData(MVRT_VAR,i,j,k))*solnData(DENS_VAR,i,j,k)*dvol
#endif
#ifdef DVVL_VAR
              lsum(vortnum+3) = lsum(vortnum+3) + abs(solnData(DVVL_VAR,i,j,k))*dvol
              lsum(vortnum+4) = lsum(vortnum+4) + abs(solnData(DVVL_VAR,i,j,k))**2*dvol
              lsum(vortnum+5) = lsum(vortnum+5) + abs(solnData(DVVL_VAR,i,j,k))*solnData(DENS_VAR,i,j,k)*dvol
#endif
              !! curl and divergence of the acceleration field
#ifdef RTRF_VAR
              lsum(vortnum+6) = lsum(vortnum+6) + abs(solnData(RTRF_VAR,i,j,k))*dvol
              lsum(vortnum+7) = lsum(vortnum+7) + abs(solnData(RTRF_VAR,i,j,k))**2*dvol
              lsum(vortnum+8) = lsum(vortnum+8) + abs(solnData(RTRF_VAR,i,j,k))*solnData(DENS_VAR,i,j,k)*dvol
#endif
#ifdef DVRF_VAR
              lsum(vortnum+ 9) = lsum(vortnum+ 9) + abs(solnData(DVRF_VAR,i,j,k))*dvol
              lsum(vortnum+10) = lsum(vortnum+10) + abs(solnData(DVRF_VAR,i,j,k))**2*dvol
              lsum(vortnum+11) = lsum(vortnum+11) + abs(solnData(DVRF_VAR,i,j,k))*solnData(DENS_VAR,i,j,k)*dvol
#endif

              !! This is where all magnetic quantities start
#ifdef MAGX_VAR
              ! magnetic energy
              lsum(magnum  ) = lsum(magnum  ) + ( (solnData(MAGX_VAR,i,j,k)**2+ &
                                                   solnData(MAGY_VAR,i,j,k)**2+ &
                                                   solnData(MAGZ_VAR,i,j,k)**2) ) / (8.0*PI)*dvol
#endif
#ifdef MAGX_VAR
              ! <Bx> volume weighted
              lsum(magnum+1) = lsum(magnum+1) + solnData(MAGX_VAR,i,j,k)*dvol
              ! rms Bx
              lsum(magnum+2) = lsum(magnum+2) + solnData(MAGX_VAR,i,j,k)*solnData(MAGX_VAR,i,j,k)*dvol
#endif
#ifdef MAGY_VAR
              ! <By> volume weighted
              lsum(magnum+3) = lsum(magnum+3) + solnData(MAGY_VAR,i,j,k)*dvol
              ! rms By
              lsum(magnum+4) = lsum(magnum+4) + solnData(MAGY_VAR,i,j,k)*solnData(MAGY_VAR,i,j,k)*dvol
#endif
#ifdef MAGZ_VAR
              ! <Bz> volume weighted
              lsum(magnum+5) = lsum(magnum+5) + solnData(MAGZ_VAR,i,j,k)*dvol
              ! rms Bz
              lsum(magnum+6) = lsum(magnum+6) + solnData(MAGZ_VAR,i,j,k)*solnData(MAGZ_VAR,i,j,k)*dvol
#endif
#ifdef PRES_VAR
#ifdef MAGX_VAR
#ifdef MAGY_VAR
#ifdef MAGZ_VAR
! #ifdef BETA_VAR
              ! mean plasma beta < 2P / B^2 >
!               lsum(magnum+7) = lsum(magnum+7) + solnData(BETA_VAR,i,j,k)*dvol
! #else
! #ifdef MAGP_VAR
!               lsum(magnum+7) = lsum(magnum+7) + (solnData(PRES_VAR,i,j,k)/solnData(MAGP_VAR,i,j,k))*dvol
! #else
              lsum(magnum+7) = lsum(magnum+7) + (8.0*PI*solnData(PRES_VAR,i,j,k) / &
                (solnData(MAGX_VAR,i,j,k)**2+solnData(MAGY_VAR,i,j,k)**2+solnData(MAGZ_VAR,i,j,k)**2))*dvol
! #endif
! #endif
              alfven_speed = sqrt( (solnData(MAGX_VAR,i,j,k)**2+solnData(MAGY_VAR,i,j,k)**2+solnData(MAGZ_VAR,i,j,k)**2) / &
                (4.0*PI*solnData(DENS_VAR,i,j,k)) )
              lmax_alfvenspeed = max(lmax_alfvenspeed, alfven_speed)
#endif
#endif
#endif
#endif
#ifdef DIVB_VAR
              ! mean div B
              lsum(magnum+8) = lsum(magnum+8) + solnData(DIVB_VAR,i,j,k)*dvol
              ! rms div B
              lsum(magnum+9) = lsum(magnum+9) + solnData(DIVB_VAR,i,j,k)*solnData(DIVB_VAR,i,j,k)*dvol
#endif

! ifdef DENS_VAR
#endif

           enddo
        enddo
     enddo

     call Grid_releaseBlkPtr(blockList(lb), solnData)

  enddo

  if (Debug .and. (io_globalMe == MASTER_PE)) print *, 'IO_writeIntegralQuantities:  local sums finished ...'

  ! Now the MASTER_PE sums the local contributions from all of
  ! the processors and writes the total to a file.
  call MPI_Reduce (lsum, gsum, nGlobalSum+1, MPI_Double_Precision, MPI_Sum, & 
       &                MASTER_PE, MPI_Comm_World, error)

  ! reduce min density
  call MPI_Reduce (lmin_dens, gmin_dens, 1, MPI_Double_Precision, MPI_Min, & 
       &                MASTER_PE, MPI_Comm_World, error)

  ! reduce max density
  call MPI_Reduce (lmax_dens, gmax_dens, 1, MPI_Double_Precision, MPI_Max, & 
       &                MASTER_PE, MPI_Comm_World, error)

  ! reduce max Alfven speed
  call MPI_Reduce (lmax_alfvenspeed, gmax_alfvenspeed, 1, MPI_Double_Precision, MPI_Max, & 
       &                MASTER_PE, MPI_Comm_World, error)

  if (Debug .and. (io_globalMe == MASTER_PE)) print *, 'IO_writeIntegralQuantities:  mpi_reduce finished ...'

  if (io_globalMe == MASTER_PE) then

        if (Debug) print *, 'IO_writeIntegralQuantities:  Now calculating and writing integral quantities...'

        ! calculate rms Mach and rms Mach mass weighted
        rms_Mach          = sqrt(gsum(8)/gsum(9))
        rms_Mach_mw       = sqrt(gsum(10)/gsum(11))
        rms_Mach_netto    = sqrt((gsum(8)-(gsum(2)**2+gsum(3)**2+gsum(4)**2)/(gsum(1)**2))/gsum(9))
        rms_Mach_netto_mw = sqrt((gsum(10)-(gsum(2)**2+gsum(3)**2+gsum(4)**2)/gsum(1))/gsum(11))
        rms_Forcing       = sqrt(gsum(12))
        mean_temp         = gsum(13) / gsum(0)
        mean_pres         = gsum(14) / gsum(0)

        mean_dens         = gsum(1) / gsum(0)
        rms_dens          = sqrt ( gsum(15) / gsum(0) )
        sigma_dens = 0.
        if (rms_dens .GT. mean_dens) sigma_dens = sqrt ( rms_dens**2.0 - mean_dens**2.0 )
        mean_ln_dens      = gsum(16) / gsum(0) - alog(mean_dens)
        rms_ln_dens       = gsum(17) / gsum(0) - 2.0*alog(mean_dens)*gsum(16)/gsum(0) + alog(mean_dens)**2.0
        if (rms_ln_dens .GT. 0) then
            rms_ln_dens = sqrt(rms_ln_dens)
        else
            rms_ln_dens = 0.
        endif
        sigma_ln_dens = 0.
        if (rms_ln_dens .GT. mean_ln_dens) sigma_ln_dens = sqrt ( rms_ln_dens**2.0 - mean_ln_dens**2.0 )

        gsum(vortnum  )    = gsum(vortnum  ) / gsum(0) ! <|vorticity|> volume weighted
        gsum(vortnum+1)    = gsum(vortnum+1) / gsum(0) ! <|vorticity|^2> volume weighted
        gsum(vortnum+2)    = gsum(vortnum+2) / gsum(1) ! <|vorticity|> mass weighted
        gsum(vortnum+3)    = gsum(vortnum+3) / gsum(0) ! <|div vel|> volume weighted
        gsum(vortnum+4)    = gsum(vortnum+4) / gsum(0) ! <|div vel|^2> volume weighted
        gsum(vortnum+5)    = gsum(vortnum+5) / gsum(1) ! <|div vel|> mass weighted
        gsum(vortnum+6)    = gsum(vortnum+6) / gsum(0) ! <|rot accel|> volume weighted
        gsum(vortnum+7)    = gsum(vortnum+7) / gsum(0) ! <|rot accel|^2> volume weighted
        gsum(vortnum+8)    = gsum(vortnum+8) / gsum(1) ! <|rot accel|> mass weighted
        gsum(vortnum+9)    = gsum(vortnum+9) / gsum(0) ! <|div accel|> volume weighted
        gsum(vortnum+10)    = gsum(vortnum+10) / gsum(0) ! <|div accel|^2> volume weighted
        gsum(vortnum+11)    = gsum(vortnum+11) / gsum(1) ! <|div accel|> mass weighted

        gsum(magnum+1)    = gsum(magnum+1) / gsum(0) ! <Bx> volume weighted
        gsum(magnum+2)    = sqrt(gsum(magnum+2) / gsum(0)) ! rms Bx
        gsum(magnum+3)    = gsum(magnum+3) / gsum(0) ! <By> volume weighted
        gsum(magnum+4)    = sqrt(gsum(magnum+4) / gsum(0)) ! rms By
        gsum(magnum+5)    = gsum(magnum+5) / gsum(0) ! <Bz> volume weighted
        gsum(magnum+6)    = sqrt(gsum(magnum+6) / gsum(0)) ! rms Bz
        gsum(magnum+7)    = gsum(magnum+7) / gsum(0) ! plasma beta
        gsum(magnum+8)    = gsum(magnum+8) / gsum(0) ! mean div B
        gsum(magnum+9)    = sqrt(gsum(magnum+9)/ gsum(0)) ! rms div B


     ioStat = 0
     open(funit, file=trim(io_statsFileName), position='APPEND', status='OLD', iostat=ioStat)
     if (ioStat .NE. 0) then
        !print *, 'FILE FOUND'
        open(funit, file=trim(io_statsFileName), position='APPEND')
     endif

     if (isFirst .EQ. 1 .AND. (.NOT. io_restart .or. ioStat .NE. 0)) then

         write (funit, 10)               &
           '#00_time                      ', &
           '#01_mass                      ', &
           '#02_x-momentum                ', &
           '#03_y-momentum                ', & 
           '#04_z-momentum                ', &
           '#05_E_total                   ', &
           '#06_E_kinetic                 ', &
           '#07_E_internal                ', &
           '#08_rms_Mach                  ', &
           '#09_rms_Mach_mw               ', &
           '#10_rms_Mach_netto            ', &
           '#11_rms_Mach_netto_mw         ', &
           '#12_rms_Forcing               ', &
           '#13_mean_temperature          ', &
           '#14_mean_pressure             ', &
           '#15_min_density               ', &
           '#16_max_density               ', &
           '#17_abs_vorticity             ', &
           '#18_abs_vorticity_sqr         ', &
           '#19_abs_vorticity_mw          ', &
           '#20_abs_div_vel               ', &
           '#21_abs_div_vel_sqr           ', &
           '#22_abs_div_vel_mw            ', &
           '#23_abs_rot_accel_norm        ', &
           '#24_abs_rot_accel_norm_sqr    ', &
           '#25_abs_rot_accel_norm_mw     ', &
           '#26_abs_div_accel             ', &
           '#27_abs_div_accel_sqr         ', &
           '#28_abs_div_accel_mw          ', &
           '#29_E_magnetic                ', &
           '#30_mean_Bx                   ', &
           '#31_rms_Bx                    ', &
           '#32_mean_By                   ', &
           '#33_rms_By                    ', &
           '#34_mean_Bz                   ', &
           '#35_rms_Bz                    ', &
           '#36_plasma_beta               ', &
           '#37_mean_divB                 ', &
           '#38_rms_divB                  ', &
           '#39_max_alfven_speed          ', &
           '#40_mean_dens                 ', &
           '#41_rms_dens                  ', &
           '#42_sigma_dens                ', &
           '#43_mean_ln_dens              ', &
           '#44_rms_ln_dens               ', &
           '#45_sigma_ln_dens             '


10         format (2x,50(a25, :, 1X))

     else if (isFirst .EQ. 1) then
        write (funit, 11)
11      format('# simulation restarted')
     endif


     write (funit, 12) simtime, &          !! time
                       gsum(1), &          !! mass
                       gsum(2), &          !! x momentum
                       gsum(3), &          !! y momentum
                       gsum(4), &          !! z momentum
                       gsum(5), &          !! total energy
                       gsum(6), &          !! kinetic energy
                       gsum(7), &          !! internal energy
                       rms_Mach, &         !! rms Mach number
                       rms_Mach_mw, &      !! rms Mach number (mass weighted)
                       rms_Mach_netto, &   !! rms Mach number (bulk motion corrected)
                       rms_Mach_netto_mw, &!! rms Mach number (bulk motion corrected, mass weighted)
                       rms_Forcing, &      !! rms random forcing
                       mean_temp, &        !! <temperature>
                       mean_pres, &        !! <pressure>
                       gmin_dens, &        !! min density
                       gmax_dens, &        !! max density
                       gsum(vortnum  ), &  !! magnitude of vorticity (vol weighted)
                       gsum(vortnum+1), &  !! magnitude of vorticity^2 (vol weighted)
                       gsum(vortnum+2), &  !! magnitude of vorticity (mass weighted)
                       gsum(vortnum+3), &  !! |divergence of velocity| (vol weighted)
                       gsum(vortnum+4), &  !! |divergence of velocity|^2 (vol weighted)
                       gsum(vortnum+5), &  !! |divergence of velocity| (mass weighted)
                       gsum(vortnum+6), &  !! magnitude of curl of random force (vol weighted)
                       gsum(vortnum+7), &  !! magnitude of curl of random force^2 (vol weighted)
                       gsum(vortnum+8), &  !! magnitude of curl of random force (mass weighted)
                       gsum(vortnum+9), &  !! |divergence of random force| (vol weighted)
                       gsum(vortnum+10), & !! |divergence of random force|^2 (vol weighted)
                       gsum(vortnum+11), & !! |divergence of random force| (mass weighted)
                       gsum(magnum  ), &   !! magnetic energy
                       gsum(magnum+1), &   !! <Bx> vol weighted
                       gsum(magnum+2), &   !! rms Bx
                       gsum(magnum+3), &   !! <By> vol weighted
                       gsum(magnum+4), &   !! rms By
                       gsum(magnum+5), &   !! <Bz> vol weighted
                       gsum(magnum+6), &   !! rms Bz
                       gsum(magnum+7), &   !! plasma beta
                       gsum(magnum+8), &   !! mean div B
                       gsum(magnum+9), &   !! rms div B
                       gmax_alfvenspeed, & !! max Alfven speed
                       mean_dens, &        !! mean density
                       rms_dens, &             !! rms density
                       sigma_dens, &       !! sigma density
                       mean_ln_dens, &     !! mean ln(density)
                       rms_ln_dens, &      !! rms ln(density)
                       sigma_ln_dens       !! sigma ln(density)


12   format (1x, 50(es25.18, :, 1x))

     close (funit)          ! Close the file.

  endif

  call MPI_Barrier (MPI_Comm_World, error)

  !=============================================================================

  if (Debug .and. (io_globalMe == MASTER_PE)) print *, 'IO_writeIntegralQuantities:  exiting ...'
  return

end subroutine IO_writeIntegralQuantities
