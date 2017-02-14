! This tests the measurement of isosurface area

! The core routines that do do the marching cubes and surface measurement
! are fairly stand-alone, and so can be tested simply outside of flash

! This program creates a box and fills it with values so that an isosurface
! creates a quadrant of a spherical shell.  This is then moved around to different
! parts of the box to test that each edge, corner, and face is clipped properly.




program UnitTest_contourSurface

use ut_contourSurfaceInterface

   implicit none
#include "constants.h"

   real, dimension(:,:,:), allocatable :: datacube
   integer, parameter :: boxsize = 32  ! this is the interior box size 
   integer, parameter :: nguard = 4
   integer :: fullsize, i,j,k
   real :: x,y,z, center(3), radius, correctarea, areas(2), isolevels(2), maxerror, error
   integer :: blkLimits(HIGH,MDIM)


   fullsize = boxsize+2*nguard ! put "guard cells" on each side

   blkLimits(LOW,IAXIS) = nguard+1
   blkLimits(HIGH,IAXIS) = nguard+boxsize
   blkLimits(LOW,JAXIS) = nguard+1
   blkLimits(HIGH,JAXIS) = nguard+boxsize
   blkLimits(LOW,KAXIS) = nguard+1
   blkLimits(HIGH,KAXIS) = nguard+boxsize

   radius = 0.5*boxsize
   correctarea = PI*radius**2

   allocate( datacube(fullsize,fullsize,fullsize) )


   print *, "Testing surface measurement with box size: ", boxsize
   print *, ""

   print *, "Spherical quarter centered on x axis, halfway along box edge"

   center(1)=0.5*boxsize
   center(2)=0.0
   center(3)=0.0

   call do_check

   print *, ""
   print *, "************************************************************"
   print *, "Spherical quarter centered on y axis, halfway along box edge"

   center(1)=0.0
   center(2)=0.5*boxsize
   center(3)=0.0

   call do_check

   print *, "************************************************************"
   print *, "Spherical quarter along upper y bottom edge, halfway along box edge"

   center(1)=0.5*boxsize
   center(2)=boxsize
   center(3)=0.0

   call do_check

   print *, ""
   print *, "************************************************************"
   print *, "Spherical quarter along upper x bottom edge, halfway along box edge"

   center(1)=boxsize
   center(2)=0.5*boxsize
   center(3)=0.0

   call do_check

   print *, ""
   print *, "************************************************************"
   print *, "Spherical quarter along top edge on lower y, halfway along box edge"

   center(1)=0.5*boxsize
   center(2)=0.0
   center(3)=boxsize

   call do_check

   print *, ""
   print *, "************************************************************"
   print *, "Spherical quarter along top edge on upper y, halfway along box edge"

   center(1)=0.5*boxsize
   center(2)=boxsize
   center(3)=boxsize

   call do_check

   print *, ""
   print *, "max error (%): ", maxerror
   ! error is in percent
   if ( maxerror < 1.0 ) then
      print *, "SUCCESS"
      stop 0
   else
      print *, "FAILURE!!"
      stop 1
   endif


   deallocate( datacube )





   contains

      subroutine do_check

         do k=1,fullsize
            z = k-nguard-0.5
            do j=1,fullsize
               y = j-nguard-0.5
               do i=1,fullsize
                  ! origin is defined as edge of interior cells
                  x = i-nguard-0.5
                  datacube(i,j,k) = sqrt( (x-center(1))**2 + (y-center(2))**2 + (z-center(3))**2 )
               enddo
            enddo
         enddo

         isolevels(1) = radius
         isolevels(2) = 0.5*radius

         print *, "  Cube filled, calculating surface area..."

         call ut_contourSurfaceAreaBlock( 2, isolevels, datacube, blkLimits, 0, areas)

         error = abs(areas(1)-correctarea)/correctarea*100.0
         maxerror = max(maxerror,error)
         print *, " sphere all the way to corner:"
         print *, "    actual area : ", correctarea
         print *, "     calculated : ", areas(1)
         print *, "       error(%) : ", error
 
         print *, ""

         error = abs(areas(2)-correctarea*0.25)*4/correctarea*100.0
         maxerror = max(maxerror,error)
         print *, " sphere half the way to corner:"
         print *, "    actual area : ", correctarea*0.25
         print *, "     calculated : ", areas(2)
         print *, "       error(%) : ", error


      end subroutine
   



end program
