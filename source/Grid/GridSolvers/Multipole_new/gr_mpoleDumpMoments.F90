!!****if* source/Grid/GridSolvers/Multipole_new/gr_mpoleDumpMoments
!!
!! NAME
!!
!!  gr_mpoleDumpMoments
!!
!! SYNOPSIS
!!
!!  gr_mpoleDumpMoments ()
!!
!! DESCRIPTION
!!
!!  Utility routine to output the regular moment array and the irregular moment
!!  array to a text file. The information is written out to a file named
!!  <basenm>MomentDump.txt, where <basenm> is the runtime parameter for output
!!  file names. The file is appended at each time for each iteration.
!!
!! NOTES
!!
!!  At the end of each iteration time, the line --- finished present iteration ---
!!  is inserted to make post-processing easier. The original distinctive
!!  phrase "Chakka Khan Chakka Khan I feel for you this is the end of this timestep"
!!  has NOT been retained. Sorry, but the author is open for discussion about
!!  including it again.
!!
!!***

subroutine gr_mpoleDumpMoments ()

  use Grid_data,                   ONLY : gr_meshMe
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  use gr_mpoleData,                ONLY : gr_mpoleMaxL,    &
                                          gr_mpoleMaxM,    &
                                          gr_mpoleMaxQ,    &
                                          gr_mpoleMomentR, &
                                          gr_mpoleMomentI, &
                                          gr_mpoleQRadii

  implicit none

#include "constants.h"
#include "gr_mpole.h"
   
  interface
     integer function ut_getFreeFileUnit()
     end function ut_getFreeFileUnit
  end interface

  logical, save :: firstCall = .true.

  integer :: fileUnit
  integer :: L,M,n
  integer :: posBlank
  integer :: Q

  real    :: radius

  character (len=MAX_STRING_LENGTH), save :: baseName
  character (len=MAX_STRING_LENGTH), save :: fileName
!
!
!   ...Do the printout only on the master processor.
!
!
  if (gr_meshMe == MASTER_PE) then
!
!
!   ...Open the printout file.
!
!
      fileUnit = ut_getFreeFileUnit ()

      if (firstCall) then
          call RuntimeParameters_get ("basenm",baseName)
          posBlank = index (baseName,' ')
          fileName = baseName (:posBlank-1) // 'MomentDump.txt'
          open  (fileUnit, file=fileName)
          firstCall = .false.
      else
          open (fileUnit, file=fileName, position='APPEND')
      end if
!
!
!     ...Printout the Moments in their order they were created and stored.
!
!
      do Q = 1,gr_mpoleMaxQ

         radius = gr_mpoleQRadii (Q)

         write (fileUnit,*)
         write (fileUnit,'(A5,I6,A30,E20.12)') &
                         ' Q = ', Q,' Radius from Center of Mass = ',radius
         write (fileUnit,*)
         write (fileUnit,'(A6,A3,A3,A20,A20)') &
                         ' PART ',' L ',' M ','   REGULAR MOMENTS  ','  IRREGULAR MOMENTS '
         write (fileUnit,*)

         n = 0

         do M = 0,gr_mpoleMaxM
         do L = M,gr_mpoleMaxL
            n = n + 1
            write (fileUnit,'(A6,I3,I3,E20.12,E20.12)') &
                            '  Cos ',L,M,gr_mpoleMomentR (n,Q),gr_mpoleMomentI (n,Q)
         end do
         end do

         do M = 1,gr_mpoleMaxM
         do L = M,gr_mpoleMaxL
            n = n + 1
            write (fileUnit,'(A6,I3,I3,E20.12,E20.12)') &
                            '  Sin ',L,M,gr_mpoleMomentR (n,Q),gr_mpoleMomentI (n,Q)
         end do
         end do

      end do

      write (fileUnit,'(A38)') '  --- finished present iteration ---  '
      close (fileUnit)

  end if
!
!
!       ...Ready!
!
!
  
  return
end subroutine gr_mpoleDumpMoments
