!!****if* source/Simulation/SimulationMain/Jeans/JeansAnalytical
!!
!! NAME
!! 
!!  JeansAnalytical
!!
!!
!! SYNOPSIS
!!
!!  JeansAnalytical
!!
!!
!! DESCRIPTION
!!  This is a stand alone serial program which calculates the analytical 
!!  solution to Jeans problem.  The values printed to file are scaled 
!!  so that their values can be easily visualised on a single graph.
!!
!!***
Program JeansAnalytical

implicit none

integer :: i
real, parameter :: pi = 3.14159265358979323846, p0 = 1.5E7, delta = 1.0E-3, &
     G = 6.67428E-8, k = 10.984, gamma = 1.6666666666667, rho0 = 1.5E7
real :: time, c0, L, omega, T, U, W, U0
character (len=*), PARAMETER :: fileName="JeansAnalytical.dat"

L = 0.5 * sqrt ((pi * gamma * rho0) / (G * (p0**2)))
c0 = sqrt ((gamma * rho0) / p0)
omega = sqrt (((c0**2) * (k**2)) - (4.0*pi*G*p0)) 

!U0 is the internal energy at time = 0.0.
U0 = -0.125 * p0 * (c0**2) * (delta**2) * (L**2) * (1.0 - cos(2.0*omega*0.0))

open(unit=50, file=filename)

do i = 0, 1250, 1

   time = real(i) / 500.0

   !T is the kinetic energy.
   T = (p0 * (delta**2) * (omega**2) * (L**2) * (1.0 - cos(2.0*omega*time))) / &
        (8.0 * (k**2))

   !U is the internal energy.
   U = -0.125 * p0 * (c0**2) * (delta**2) * (L**2) * (1.0 - cos(2.0*omega*time))

   !W is the potential energy.
   W = -((pi * G * (p0**2) * (delta**2) * (L**2)) / (2.0 * (k**2))) * (1.0 + cos(2.0*omega*time))

   write(50,*)time, T, U-U0, W*10.0

end do

close(50)

End Program JeansAnalytical


!Format the simulation data using:
!sed '1,2d' jeans.dat | awk '{ ut=$8; u0=2.9452223072102301e+07; ur=(ut-u0); printf "%.16e, %.16e, %.16e, %.16e\n", $1,$7,ur,$9*10.0 }' > SimulationJeans.dat
!Here, 2.9452223072102301e+07 is the internal energy at time=0 (first row in jeans.dat data file).

!Plot analytical against simulation data using gnuplot.
!plot "JeansAnalytical.dat" using 1:2 title "Kinetic Energy (analytical) T(t)" with lines lt 3, \
!"JeansAnalytical.dat" using 1:3 title "Internal Energy (analytical) U(t)-U(0)" with lines lt 4, \
!"JeansAnalytical.dat" using 1:4 title "Potential Energy (analytical) W(t)*10" with lines lt 1, \
!"JeansSimulation.dat" using 1:2 title "Kinetic Energy (simulation) T(t)" with points lt 3 pt 1 ps 2.5, \
!"JeansSimulation.dat" using 1:3 title "Internal Energy (simulation) U(t)-U(0)" with points lt 4 pt 2 ps 2.5, \
!"JeansSimulation.dat" using 1:4 title "Potential Energy (simulation) W(t)*10" with points lt 1 pt 3 ps 2.5
