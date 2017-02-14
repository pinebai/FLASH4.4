!!******************************************************
!! see the MAIN routine for user changes (last routine)!
!!******************************************************

Module Stir_data
!!
!! DESCRIPTION
!!
!!  Stores the turbulence stirring data in a file that is then read by the FLASH code
!!  to generate the patterns for driving of the turbulence as a time sequence
!!  The driving sequence follows an Ornstein-Uhlenbeck (OU) process.
!!
!!  For example applications see Federrath et al. 2008, ApJ 688, L79,
!!  Federrath et al. (2010, A&A 512, A81)
!!
!! PARAMETERS
!!
!!    ndim [INTEGER]
!!        number of spatial dimensions
!!    xmin, xmax, ymin, ymax, zmin, zmax [REAL]
!!        physical coordinates of the turbulent box
!!    st_decay [REAL]
!!        autocorrelation time for forcing
!!    st_energy [REAL]
!!        energy input/mode
!!    st_stirmin, st_stirmax [REAL]
!!        minimum, maximum stirring *wavenumber*
!!    st_solweight [REAL]
!!        solenoidal weight (solenoidal field 1.0, compressive field 0.0)
!!    st_seed [INTEGER]
!!        random number generator seed
!!    st_spectform [INTEGER]
!!        0: Band, 1:Parabola
!!
!! AUTHOR: Christoph Federrath, 2008 (adopted from FLASH)
!!
!!***

integer, parameter :: st_maxmodes = 1000

! OU variance corresponding to decay time and energy input rate
real, save :: st_OUvar

! number of modes
integer, save :: st_nmodes

real,save, dimension(3,st_maxmodes) :: st_mode, st_aka, st_akb
real,save, dimension(6*st_maxmodes) :: st_OUphases
real,save, dimension(  st_maxmodes) :: st_ampl

integer, save :: ndim ! number of spatial dimensions
real, save    :: xmin, xmax, ymin, ymax, zmin, zmax
real,save     :: st_decay
real,save     :: st_energy
real,save     :: st_stirmin, st_stirmax
real,save     :: st_solweight
real,save     :: st_solweightnorm
integer,save  :: st_seed
integer,save  :: st_spectform

contains
  

!!******************************************************
subroutine init_stir()
!!******************************************************

  implicit none

  integer            :: ikxmin, ikxmax, ikymin, ikymax, ikzmin, ikzmax
  integer            :: ikx, iky, ikz
  real               :: kx, ky, kz, k, kc, Lx, Ly, Lz
  real, parameter    :: twopi = 6.283185307
  ! the amplitude of the modes at kmin and kmax for a parabolic Fourier spectrum wrt 1.0 at the centre kc
  real, parameter    :: amin = 0.0
  logical, parameter :: Debug = .false.

  ! initialize some variables, allocate random seed
  st_OUvar = sqrt(st_energy/st_decay)
  kc       = 0.5*(st_stirmin+st_stirmax)
  ! this makes the rms force const irrespective of the solenoidal weight
  if (ndim.eq.3) st_solweightnorm = sqrt(3.0/3.0)*sqrt(3.0)*1.0/sqrt(1.0-2.0*st_solweight+3.0*st_solweight**2.0)
  if (ndim.eq.2) st_solweightnorm = sqrt(3.0/2.0)*sqrt(3.0)*1.0/sqrt(1.0-2.0*st_solweight+2.0*st_solweight**2.0)
  if (ndim.eq.1) st_solweightnorm = sqrt(3.0/1.0)*sqrt(3.0)*1.0/sqrt(1.0-2.0*st_solweight+1.0*st_solweight**2.0)

  ikxmin = 0
  ikymin = 0
  ikzmin = 0

  ikxmax = 20
  ikymax = 0
  ikzmax = 0
  if (ndim.gt.1) ikymax = 20
  if (ndim.gt.2) ikzmax = 20

  Lx = xmax - xmin
  Ly = ymax - ymin
  Lz = zmax - zmin

  st_nmodes = 0

  do ikx = ikxmin, ikxmax
      kx = twopi * ikx / Lx

      do iky = ikymin, ikymax
          ky = twopi * iky / Ly

          do ikz = ikzmin, ikzmax
              kz = twopi * ikz / Lz

              k = sqrt( kx*kx+ky*ky+kz*kz )

              if ((k .ge. st_stirmin) .and. (k .le. st_stirmax)) then

                 if ((st_nmodes + 2**(ndim-1)) .gt. st_maxmodes) then

                    print *, 'init_stir:  st_nmodes = ', st_nmodes, ' maxstirmodes = ',st_maxmodes
                    print *, 'Too many stirring modes'
                    exit

                 endif

                 st_nmodes = st_nmodes + 1

                 if (st_spectform == 0) st_ampl(st_nmodes) = 1.0
                 if (st_spectform == 1) st_ampl(st_nmodes) = 4.0*(amin-1.0)/((st_stirmax-st_stirmin)**2.0)*((k-kc)**2.0)+1.0
                 if (Debug) print *, "init_stir:  st_ampl(",st_nmodes,") = ", st_ampl(st_nmodes)

                 st_mode(1,st_nmodes) = kx
                 st_mode(2,st_nmodes) = ky
                 st_mode(3,st_nmodes) = kz

                 if (ndim.gt.1) then

                   st_nmodes = st_nmodes + 1

                   if (st_spectform == 0) st_ampl(st_nmodes) = 1.0
                   if (st_spectform == 1) st_ampl(st_nmodes) = 4.0*(amin-1.0)/((st_stirmax-st_stirmin)**2.0)*((k-kc)**2.0)+1.0
                   if (Debug) print *, "init_stir:  st_ampl(",st_nmodes,") = ", st_ampl(st_nmodes)

                   st_mode(1,st_nmodes) = kx
                   st_mode(2,st_nmodes) =-ky
                   st_mode(3,st_nmodes) = kz

                 endif

                 if (ndim.gt.2) then

                   st_nmodes = st_nmodes + 1

                   if (st_spectform == 0) st_ampl(st_nmodes) = 1.0
                   if (st_spectform == 1) st_ampl(st_nmodes) = 4.0*(amin-1.0)/((st_stirmax-st_stirmin)**2.0)*((k-kc)**2.0)+1.0
                   if (Debug) print *, "init_stir:  st_ampl(",st_nmodes,") = ", st_ampl(st_nmodes)

                   st_mode(1,st_nmodes) = kx
                   st_mode(2,st_nmodes) = ky
                   st_mode(3,st_nmodes) =-kz

                   st_nmodes = st_nmodes + 1

                   if (st_spectform == 0) st_ampl(st_nmodes) = 1.0
                   if (st_spectform == 1) st_ampl(st_nmodes) = 4.0*(amin-1.0)/((st_stirmax-st_stirmin)**2.0)*((k-kc)**2.0)+1.0
                   if (Debug) print *, "init_stir:  st_ampl(",st_nmodes,") = ", st_ampl(st_nmodes)

                   st_mode(1,st_nmodes) = kx
                   st_mode(2,st_nmodes) =-ky
                   st_mode(3,st_nmodes) =-kz

                 endif

              endif

          enddo
      enddo
  enddo

  write (*,'(A,I4,A)') 'initialized ',st_nmodes,' modes for stirring.'
  if (st_spectform == 0) write (*,'(A,I2,A)') ' spectral form        = ', st_spectform, ' (Band)'
  if (st_spectform == 1) write (*,'(A,I2,A)') ' spectral form        = ', st_spectform, ' (Parabola)'
  write (*,'(A,ES10.3)') ' solenoidal weight    = ', st_solweight
  write (*,'(A,ES10.3)') ' st_solweightnorm     = ', st_solweightnorm
  write (*,'(A,ES10.3)') ' stirring energy      = ', st_energy
  write (*,'(A,ES10.3)') ' autocorrelation time = ', st_decay
  write (*,'(A,ES10.3)') ' minimum wavenumber   = ', st_stirmin
  write (*,'(A,ES10.3)') ' maximum wavenumber   = ', st_stirmax
  write (*,'(A,I8)')    ' random seed          = ', st_seed

  return

end subroutine init_stir


!!******************************************************
subroutine st_ounoiseinit(vector, vectorlength, variance)
!!******************************************************
!! initialize pseudo random sequence for the OU process

  implicit none

  real, intent (INOUT)    :: vector (:)
  integer, intent (IN)    :: vectorlength
  real, intent (IN)       :: variance
  real                    :: grnval
  integer                 :: i

  do i = 1, vectorlength
     call st_grn (grnval)
     vector (i) = grnval * variance
  end do

  return

end subroutine st_ounoiseinit


!!******************************************************
subroutine st_ounoiseupdate(vector, vectorlength, variance, dt, ts)
!!******************************************************
!! Subroutine updates a vector of real values according to an algorithm
!!   that generates an Ornstein-Uhlenbeck sequence.
!!
!! The sequence x_n is a Markov process that takes the previous value,
!!   weights by an exponential damping factor with a given correlation
!!   time "ts", and drives by adding a Gaussian random variable with
!!   variance "variance", weighted by a second damping factor, also
!!   with correlation time "ts". For a timestep of dt, this sequence
!!   can be written as :
!!
!!     x_n+1 = f x_n + sigma * sqrt (1 - f**2) z_n
!!
!! where f = exp (-dt / ts), z_n is a Gaussian random variable drawn
!! from a Gaussian distribution with unit variance, and sigma is the
!! desired variance of the OU sequence. (See Bartosch, 2001).
!!
!! The resulting sequence should satisfy the properties of zero mean,
!!   and stationary (independent of portion of sequence) RMS equal to
!!   "variance". Its power spectrum in the time domain can vary from
!!   white noise to "brown" noise (P (f) = const. to 1 / f^2).
!!
!! References :
!!    Bartosch, 2001
!! http://octopus.th.physik.uni-frankfurt.de/~bartosch/publications/IntJMP01.pdf
!!   Finch, 2004
!! http://pauillac.inria.fr/algo/csolve/ou.pdf
!!         Uhlenbeck & Ornstein, 1930
!! http://prola.aps.org/abstract/PR/v36/i5/p823_1
!!
!! Eswaran & Pope 1988
!!
!! ARGUMENTS
!!
!!   vector :       vector to be updated
!!   vectorlength : length of vector to be updated
!!   variance :     variance of the distribution
!!   dt :           timestep
!!   ts :           autocorrelation time
!!
!!***

  implicit none

  real, intent (INOUT) :: vector (:)
  integer, intent (IN) :: vectorlength
  real, intent (IN)    :: variance, dt, ts
  real                 :: grnval, damping_factor
  integer              :: i

  damping_factor = exp(-dt/ts)

  do i = 1, vectorlength
     call st_grn (grnval)
     vector (i) = vector (i) * damping_factor + variance *   &
          sqrt(1.0 - damping_factor**2.0) * grnval
  end do

  return

end subroutine st_ounoiseupdate


!!******************************************************
subroutine st_grn(grnval)
!!******************************************************
!!
!! DESCRIPTION
!!
!!  Subroutine draws a number randomly from a Gaussian distribution
!!    with the standard uniform distribution function "random_number"
!!    using the Box-Muller transformation in polar coordinates. The
!!    resulting Gaussian has unit variance.
!!
!!***

  implicit none

  real, intent (OUT) :: grnval
  real               :: pi, r1, r2, g1

  pi = 4.0*atan(1.0)
  r1 = 0.0; r2 = 0.0;
  r1 = ran1s(st_seed)
  r2 = ran1s(st_seed)
  g1 = sqrt(2.0*log(1.0/r1))*cos(2.0*pi*r2)

  grnval = g1

  return

end subroutine st_grn

!!******************************************************
function ran1s(idum)
!!******************************************************
  integer idum,IA,IM,IQ,IR,NTAB,NDIV
  real ran1s,AM,EPS,RNMX
  parameter (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836, &
  NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
  integer k,iy
  if (idum.le.0) then
    idum=max(-idum,1)
  endif
  k=idum/IQ
  idum=IA*(idum-k*IQ)-IR*k
  if (idum.lt.0) idum=idum+IM
  iy=idum
  ran1s=min(AM*iy,RNMX)
  return

end function ran1s


!!******************************************************
subroutine st_calcPhases()
!!******************************************************
!!
!! DESCRIPTION
!!
!!     This routine updates the stirring phases from the OU phases.
!!     It copies them over and applies the projection operator.
!!
!!***

  implicit none

  real                 :: ka, kb, kk, diva, divb, curla, curlb
  integer              :: i,j
  logical, parameter   :: Debug = .false.

  do i = 1, st_nmodes
     ka = 0.0
     kb = 0.0
     kk = 0.0
     do j=1, ndim
        kk = kk + st_mode(j,i)*st_mode(j,i)
        ka = ka + st_mode(j,i)*st_OUphases(6*(i-1)+2*(j-1)+1+1)
        kb = kb + st_mode(j,i)*st_OUphases(6*(i-1)+2*(j-1)+0+1)
     enddo
     do j=1, ndim

        diva  = st_mode(j,i)*ka/kk
        divb  = st_mode(j,i)*kb/kk
        curla = (st_OUphases(6*(i-1)+2*(j-1) + 0 + 1) - divb)
        curlb = (st_OUphases(6*(i-1)+2*(j-1) + 1 + 1) - diva)

        st_aka(j,i) = st_solweight*curla+(1.0-st_solweight)*divb
        st_akb(j,i) = st_solweight*curlb+(1.0-st_solweight)*diva

! purely compressive
!         st_aka(j,i) = st_mode(j,i)*kb/kk
!         st_akb(j,i) = st_mode(j,i)*ka/kk

! purely solenoidal
!         st_aka(j,i) = bjiR - st_mode(j,i)*kb/kk
!         st_akb(j,i) = bjiI - st_mode(j,i)*ka/kk

        if (Debug) then
                print *, 'st_mode(dim=',j,',mode=',i,') = ', st_mode(j,i)
                print *, 'st_aka (dim=',j,',mode=',i,') = ', st_aka(j,i)
                print *, 'st_akb (dim=',j,',mode=',i,') = ', st_akb(j,i)
        endif

     enddo
  enddo

  return

end subroutine st_calcPhases


!!******************************************************
subroutine write_forcing_file(outfile, nsteps, step, time, end_time)
!!******************************************************
!!
!! DESCRIPTION
!!
!!     writes time dependent modes, phases and amplitudes to file
!!
!!***

  implicit none

  character (len=80), intent(in) :: outfile
  integer, intent(in)            :: nsteps, step
  real, intent(in)               :: time, end_time
  logical, save                  :: first_call = .true.
  integer                        :: iostat

  if (first_call) then
     open (unit=42, file=outfile, iostat=iostat, status='REPLACE', action='WRITE', &
           access='SEQUENTIAL', form='UNFORMATTED')
     ! header contains number of times and number of modes, end time, autocorrelation time, ...
     if (iostat.eq.0) then
        write (unit=42) nsteps, st_nmodes, end_time, st_decay, st_energy, st_solweight, &
                        st_solweightnorm, st_stirmin, st_stirmax, st_spectform
     else
        write (*,'(2A)') 'could not create file for write. filename: ', trim(outfile)
     endif
     close (unit=42)
     first_call = .false.
  endif

  open (unit=42, file=outfile, iostat=iostat, status='OLD', action='WRITE', &
        position='APPEND', access='SEQUENTIAL', form='UNFORMATTED')

  if (iostat.eq.0) then
     write (unit=42) step, time, &
                     st_mode    (:, 1:  st_nmodes), &
                     st_aka     (:, 1:  st_nmodes), &
                     st_akb     (:, 1:  st_nmodes), &
                     st_ampl    (   1:  st_nmodes), &
                     st_OUphases(   1:6*st_nmodes)
  else
     write (*,'(2A)') 'could not open file for write. filename: ', trim(outfile)
  endif

  close (unit=42)

  return

end subroutine write_forcing_file


end Module Stir_data


!!******************************************************
!!*********************** MAIN *************************
!!******************************************************
program generate_forcing_file
!!******************************************************

  use Stir_data

  implicit none

  real, save               :: end_time, time, dt
  integer,save             :: step, nsteps, m
  character (len=80), save :: outfilename
  logical                  :: print_stuff_to_shell = .true.

  ! =====================================================================
  ! =========================== user changes ============================
  ! =====================================================================
  ndim = 3;                ! 3-dimensional forcing
  
  xmin = -0.5; xmax = 0.5; ! spatial x-coordinates of the box (xmax-xmin=L_box)
  ymin = -0.5; ymax = 0.5; ! spatial y-coordinates of the box (ymax-ymin=L_box)
  zmin = -0.5; zmax = 0.5; ! spatial z-coordinates of the box (zmax-zmin=L_box)
                           ! Must be adjusted to physical problem

  st_spectform = 1                   ! 0 is band, 1 is paraboloid

  st_decay     = 0.5                 ! autocorrelation time, T=L_box/(2V)
                                     ! Must be adjusted to physical problem

  st_energy    = 2e-3                ! driving amplitude ~ sqrt(energy/decay)
                                     ! Note that energy ~ decay^-3 ~ Mach^3 ~ L_box^-1
                                     ! Must be adjusted to physical problem

  st_stirmin   = 6.283               ! <~ 2 pi / L_box ; k=1.0
  st_stirmax   = 18.95               ! >~ 6 pi / L_box ; k=3.0
                                     ! minimum and maximum wavenumbers stirred
                                     ! Must be adjusted to physical problem

  st_solweight = 1.0                 ! 1.0 solenoidal; 0.0 compressive; 0.5 natural mixture
                                     ! See Federrath et al. (2010, A&A 512, A81) for details
                                     ! Must be adjusted to physical problem

  st_seed      = 140281              ! random seed

  end_time     = 5.                  ! forcing table end time
                                     ! Must be adjusted to physical problem

  nsteps       = 100.                ! total number of time dumps to update the
                                     ! forcing pattern until end_time
                                     ! Must be adjusted to physical problem

  outfilename  = 'forcingfile_seed140281.dat'   ! output filename read by FLASH

  print_stuff_to_shell = .false.     ! print additional information
  ! =====================================================================

  call init_stir()
  call st_ounoiseinit(st_OUphases, 6*st_nmodes, st_OUvar)

  write(*,'(A)') '**********************************************************************************************************'

  dt = end_time / nsteps
  time = 0.0
  do step = 0, nsteps
     call st_ounoiseupdate(st_OUphases, 6*st_nmodes, st_OUvar, dt, st_decay)
     call st_calcPhases()
     call write_forcing_file(outfilename, nsteps, step, time, end_time)
     if (print_stuff_to_shell) then
        write(*,'(A,I6,A,ES10.3)') 'step = ', step, '  time = ', time
        do m = 1, st_nmodes
           write(*,'(A,I4,A,3(1X,ES9.2),A,3(1X,ES9.2),A,3(1X,ES9.2),A,1(1X,ES9.2))') &
                 'm=', m, '  st_mode:', st_mode(:,m), '  st_aka:', st_aka(:,m), '  st_akb:', st_akb(:,m), '  st_ampl:', st_ampl(m)
        enddo
        write(*,'(A)') '**********************************************************************************************************'
     endif
     time = time + dt
  enddo

  write(*,'(3A,I6,A,ES10.3,A)') 'generate_forcing_file:  outputfile "', trim(outfilename), &
        '" containing ', nsteps, ' times with end time ', end_time, ' written. Finished.'

end program generate_forcing_file

!!******************************************************
