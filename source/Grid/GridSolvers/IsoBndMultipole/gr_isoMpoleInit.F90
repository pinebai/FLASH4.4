!!****if* source/Grid/GridSolvers/IsoBndMultipole/gr_isoMpoleInit
!!
!! NAME
!!
!!  gr_isoMpoleInit
!!
!! 
!! SYNOPSIS
!!
!!  gr_isoMpoleInit()
!!
!!
!! DESCRIPTION
!!
!!  Initialize the isolated boundary multipole solver.  Read in any of the
!!  runtime parameters for this solver.  All solver common data
!!  is stored in the mpole_common module
!!
!! NOTES
!!   The regular multipole solver (Grid/GridSolvers/Multipole) can handle
!!   the 3D axisymmetric case with the use of a parameter mpole_3daxisymmetric.
!!  This unit CANNOT, nor does it check for incorrect usage.
!!
!!***

subroutine gr_isoMpoleInit

  use gr_isoMpoleData  ! ONLY not required in init routines
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use PhysicalConstants_interface, ONLY:  PhysicalConstants_get
  use Driver_interface, ONLY: Driver_abortFlash
  use Grid_data, ONLY : gr_geometry, gr_meshMe,&
                        gr_imin,gr_imax,gr_jmin,gr_jmax,gr_kmin,gr_kmax

#include "Flash.h"
#include "constants.h"  

#ifdef FLASH_GRID_PARAMESH
  use Grid_data, ONLY : gr_nblockX, gr_nblockY,gr_nblockZ
#else
  use Grid_Data,ONLY : gr_gIndexSize
#endif

  implicit none

  
  real    :: Nxmax, Nymax, Nzmax, dxmin, dymin, dzmin
  real    :: factrl
  integer :: istat, m, l
  integer, save :: lrefine_max  
  logical :: geo_2dcyl
  character(len=128) :: str_buffer
  !========================================================================
  
  !               Initialize database keys, constants, etc.
  
  twopi  = 2.*PI
  fourpi = 4.*PI
  fourpi_inv = 1.0/fourpi
  Nint_inv = 1.0/float(Nint)
  
  call RuntimeParameters_get( "mpole_lmax", mpole_lmax)
  
  call RuntimeParameters_get( "quadrant", mpole_quadrant)
  
  call RuntimeParameters_get("mpole_lmax", mpole_lmax)

  call RuntimeParameters_get("octant",   mpole_octant)

  if (mpole_quadrant) then
     cylfactor = 2.0*twopi
  else
     cylfactor = twopi
  endif




!               Check if we support the requested grid geometry.
!               Only bounded grid geometries are supported, and
!               support for the two which use angular coordinates
!               (3D cylindrical and 3D spherical) is deferred for now.

  geo_2dcyl=.false.
  if ((NDIM == 3) .and. (gr_geometry == CARTESIAN)) then
     
     mpole_geometry = CARTESIAN
     mpole_mmax = mpole_lmax
 
  else if ((NDIM == 2) .and. (gr_geometry == CYLINDRICAL)) then

     mpole_geometry = CYLINDRICAL
     geo_2dcyl=.true.
     mpole_mmax = 0

  else if ((NDIM == 1) .and. (gr_geometry == SPHERICAL)) then
     
     mpole_geometry = SPHERICAL
     mpole_lmax = 0
     mpole_mmax = 0
     
  else

     call Driver_abortFlash ('gr_isoMpoleInit:  FATAL:  unsupported geometry!')
     
  endif
  
  ! we are only allowed to do a quadrant (i.e. enforce reflection symmetry
  ! about y=0) if we are in 2-d cylindrical coords
  if ((mpole_quadrant) .AND. (.not.geo_2dcyl)) then
     call Driver_abortFlash('ERROR: quadrant only allowed in 2-d cylindrical geometry')
  endif
  
  
  !======================================================================
  
  !               Maximum number of zones across each dimension, if each
  !               dimension were to become fully refined.  Also compute
  !               minimum zone spacings at the maximum refinement level.
  !               Assume the domain is a Cartesian box.
  
#ifdef FLASH_GRID_PARAMESH

  call RuntimeParameters_get("lrefine_max", lrefine_max)
  Nxmax = gr_nblockX * NXB * 2.**(lrefine_max-1)
  
  dxmin = (gr_imax - gr_imin) / Nxmax
  
  if (NDIM >= 2) then
     Nymax = gr_nblockY * NYB * 2.**(lrefine_max-1)
     dymin = (gr_jmax - gr_jmin) / Nymax
  else
     Nymax = 0
     dymin = 1.
  endif

  if (NDIM == 3) then
     Nzmax = gr_nblockZ * NZB * 2.**(lrefine_max-1)
     dzmin = (gr_kmax - gr_kmin) / Nzmax
  else
     Nzmax = 0
     dzmin = 1.
  endif
#else
  Nxmax=gr_gIndexSize(IAXIS)
  Nymax=0
  Nzmax=0
  if(NDIM>1)Nymax=gr_gIndexSize(JAXIS)
  if(NDIM>2)Nzmax=gr_gIndexSize(KAXIS)
#endif
!               Inverse of sample spacing to use for moment arrays.
  
  dsinv = 2. / (dxmin*dymin*dzmin)**(1./NDIM)
  
  !               Number of radial samples to use in moment arrays.
  
  qmax = 2 * int( sqrt(Nxmax**2 + Nymax**2 + Nzmax**2) ) + 3
  
!  if (gr_meshMe == 1) then
!     write (str_buffer,*) 'gr_isoMpoleInit:  using ', qmax, ' radial samples', &
!          & ' moment array:  ', & 
!          & qmax*2*2*(mpole_lmax+1)*(mpole_mmax+1), ' items'
!     call write_logfile(str_buffer)
!  endif
  
  !               Allocate moment arrays and other data structures.
  
  allocate ( Moment(0:qmax,1:2,1:2,0:mpole_lmax,0:mpole_mmax), stat=istat)
  if (istat>0) call Driver_abortFlash ("gr_isoMpoleInit: Moment() allocate failed!")
  
  allocate( Momtmp(0:qmax), stat=istat)
  if (istat>0) call Driver_abortFlash ("gr_isoMpoleInit: Momtmp() allocate failed!")
  
  allocate( costable(0:mpole_mmax), stat=istat)
  if (istat>0) call Driver_abortFlash ("gr_isoMpoleInit: costable() allocate failed!")
  
  allocate( sintable(0:mpole_mmax), stat=istat)
  if (istat>0) call Driver_abortFlash ("gr_isoMpoleInit: sintable() allocate failed!")
  
  allocate( rpower(0:mpole_lmax), stat=istat)
  if (istat>0) call Driver_abortFlash ("gr_isoMpoleInit: rpower() allocate failed!")
  
  allocate( Legk1(0:mpole_lmax,0:mpole_mmax), stat=istat)
  if (istat>0) call Driver_abortFlash ("gr_isoMpoleInit: Legk1() allocate failed!")
  
  allocate( rprinv(0:mpole_lmax), stat=istat)
  if (istat>0) call Driver_abortFlash ("gr_isoMpoleInit: rprinv() allocate failed!")
  
  allocate( Legk2(0:mpole_lmax,0:mpole_mmax), stat=istat)
  if (istat>0) call Driver_abortFlash ("gr_isoMpoleInit: Legk2() allocate failed!")
  
  allocate( Leg_fact(0:mpole_lmax,0:mpole_mmax), stat=istat)
  if (istat>0) call Driver_abortFlash ("gr_isoMpoleInit: leg_fact() allocate failed!")
  
  !                       Coefficients for Legendre polynomials.
  
  do m = 0, mpole_mmax
     do l = m+2, mpole_lmax
        Legk1(l,m) = real(2*l - 1) / real(l - m)
        Legk2(l,m) = real(l + m - 1) / real(l - m)
     enddo
  enddo
  
  if ((mpole_lmax > 0) .and. (mpole_mmax > 0)) then
     
     do l = 1, mpole_lmax
        factrl = 2.
        do m = 1, l
           factrl = factrl / real((l+m) * (l-m+1))
           Leg_fact(l,m) = factrl
        enddo
     enddo
     
  endif
  
  !====================================================================
  
  return
end subroutine gr_isoMpoleInit
