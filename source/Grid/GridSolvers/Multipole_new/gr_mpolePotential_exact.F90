!!****if* source/Grid/GridSolvers/Multipole_new/gr_mpolePotential_exact
!!
!! NAME
!!
!!  gr_mpolePotential_exact
!!
!! SYNOPSIS
!!
!!  gr_mpolePotential_exact (integer, intent(in) :: idensvar,
!!                           integer, intent(in) :: ipotvar,
!!                           real,    intent(in) :: poissonFactor )
!!
!! DESCRIPTION
!!
!!  Computes the potential field in the most exact way possible in
!!  FLASH and compares it to the one obtained by using Moments.
!!  This can only be done if the run is on 1 processor only.
!!  A simple O(N^2) algorithm over all cells is employed. This
!!  routine is for testing purposes only. The # of processors = 1
!!  is checked and the calculation is aborted if this condition is
!!  not met.
!!
!!  The exact FLASH potential at a cell is obtained excluding the
!!  internal cell potential contribution.
!!
!!                       !!! IMPORTANT !!!
!!
!!  This routine is intended for internal testing purposes by the author
!!  of the code. It should never be activated for an actual FLASH
!!  application. It might also be outdated and in need of a fixup, since
!!  the routine does not get used often.
!!
!! ARGUMENTS
!!
!!  idensvar       : index to variable containing the density
!!  ipotvar        : index to variable containing the potential
!!  poissonFactor  : the name says it all 
!!
!!***

!!REORDER(4): solnData

subroutine gr_mpolePotential_exact (idensvar,ipotvar,poissonFactor)

  use Grid_data,                   ONLY : gr_meshNumProcs
  use Driver_interface,            ONLY : Driver_abortFlash
  use Logfile_interface,           ONLY : Logfile_stamp
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  use Grid_interface,              ONLY : Grid_getBlkPtr,        &
                                          Grid_releaseBlkPtr,    &
                                          Grid_getBlkBoundBox,   &
                                          Grid_getDeltas,        &
                                          Grid_getBlkIndexLimits

  use gr_mpoleData,                ONLY : gr_mpoleSymmetryPlane2D, &
                                          gr_mpoleTwoPi,           &
                                          gr_mpoleFourPi,          &
                                          gr_mpoleFourPiInv,       &
                                          gr_mpoleGeometry,        &
                                          gr_mpoleBlockCount,      &
                                          gr_mpoleBlockList

  implicit none
  
#include "constants.h"
#include "gr_mpole.h"
#include "Flash.h"
  
  integer, intent (in) :: idensvar
  integer, intent (in) :: ipotvar
  real,    intent (in) :: poissonFactor

  real    :: bndBoxILow_A,bndBoxILow_B
  real    :: bndBoxJLow_A,bndBoxJLow_B
  real    :: bndBoxKLow_A,bndBoxKLow_B
  real    :: cellVolume
  real    :: cellMass
  real    :: DeltaI_A,DeltaI_B
  real    :: DeltaJ_A,DeltaJ_B
  real    :: DeltaK_A,DeltaK_B
  real    :: DeltaIHalf_A,DeltaIHalf_B
  real    :: DeltaJHalf_A,DeltaJHalf_B
  real    :: DeltaKHalf_A,DeltaKHalf_B
  real    :: gr_mpoleGravityConstant
  real    :: potential_exact, potential_mpole
  real    :: potdiff, potdiff_max
  real    :: r
  real    :: Rcyl_A,Rcyl_B
  real    :: Rsph_A,Rsph_B
  real    :: xA,xB,yA,yB,zA,zB

  real    :: delta_A (1:MDIM)
  real    :: delta_B (1:MDIM)

  real    :: bndBox_A      (LOW:HIGH,MDIM)
  real    :: bndBox_B      (LOW:HIGH,MDIM)
  integer :: blkLimits_A   (LOW:HIGH,MDIM)
  integer :: blkLimits_B   (LOW:HIGH,MDIM)
  integer :: blkLimitsGC_A (LOW:HIGH,MDIM)
  integer :: blkLimitsGC_B (LOW:HIGH,MDIM)

  real, pointer :: solnData (:,:,:,:)

  integer :: blockNr_A,blockNr_B
  integer :: blockID_A,blockID_B
  integer :: fileUnit
  integer :: iA,jA,kA
  integer :: iB,jB,kB
  integer :: imax_A, jmax_A, kmax_A
  integer :: imax_B, jmax_B, kmax_B
  integer :: imin_A, jmin_A, kmin_A
  integer :: imin_B, jmin_B, kmin_B
  integer :: posBlank

  logical :: firstCall = .true.
  logical :: same_cube
  logical :: same_ring
  logical :: same_layer

  character (len=MAX_STRING_LENGTH), save :: baseName
  character (len=MAX_STRING_LENGTH), save :: fileName
!  
!  ==========================================================================
!
!
!        1) if # of processors different from 1 -> abort.
!        2) if symmetry present -> abort.
!
!
  if (gr_meshNumProcs /= 1) then
      call Driver_abortFlash ('[gr_mpolePotential_exact] ERROR: # processors > 1')
  end if

  if (gr_mpoleSymmetryPlane2D) then
      call Driver_abortFlash ('[gr_mpolePotential_exact] ERROR: gr_mpoleSymmetryPlane2D = true')
  end if
!
!
!      This assignment of some random file unit number is not good
!      programming practice! Better to have a unit file number management
!      program in FLASH.
!
!
  fileUnit = 237     ! Let's hope the file number is not in use.

  if (firstCall) then
      call RuntimeParameters_get ("basenm",baseName)
      posBlank = index (baseName,' ')
      fileName = baseName (:posBlank-1) // 'comparePotentials.txt'
      open  (fileUnit, file=fileName)
      firstCall = .false.
  else
      open (fileUnit, file=fileName, position='APPEND')
  end if

  write (fileUnit,*)
  write (fileUnit,'(A7,A7,A7,A7,A20,A20)') &
                  ' BLOCK ','   I   ','   J   ','   K   ','  EXACT POTENTIALS  ','  MPOLE POTENTIALS  '
  write (fileUnit,*)
!
!
!        O(N^2) algorithm over all locally held leaf blocks.
!
!
  gr_mpoleGravityConstant = poissonFactor * gr_mpoleFourPiInv

  potdiff_max = ZERO

  do blockNr_A = 1,gr_mpoleBlockCount
     blockID_A = gr_mpoleBlockList (blockNr_A)

     call Grid_getBlkBoundBox     (blockID_A,  bndBox_A)
     call Grid_getDeltas          (blockID_A,  delta_A )
     call Grid_getBlkIndexLimits  (blockID_A,  blkLimits_A, blkLimitsGC_A)

     imin_A         = blkLimits_A (LOW, IAXIS)
     jmin_A         = blkLimits_A (LOW, JAXIS)
     kmin_A         = blkLimits_A (LOW, KAXIS)  
     imax_A         = blkLimits_A (HIGH,IAXIS)
     jmax_A         = blkLimits_A (HIGH,JAXIS)
     kmax_A         = blkLimits_A (HIGH,KAXIS)
     DeltaI_A      = delta_A (IAXIS)
     DeltaJ_A      = delta_A (JAXIS)
     DeltaK_A      = delta_A (KAXIS)
     DeltaIHalf_A = DeltaI_A * HALF
     DeltaJHalf_A = DeltaJ_A * HALF
     DeltaKHalf_A = DeltaK_A * HALF
     bndBoxILow_A = bndBox_A (LOW,IAXIS)
     bndBoxJLow_A = bndBox_A (LOW,JAXIS)
     bndBoxKLow_A = bndBox_A (LOW,KAXIS)

     select case (gr_mpoleGeometry)

     case (GRID_3DCARTESIAN)

       zA = bndBoxKLow_A + DeltaKHalf_A
       do kA = kmin_A,kmax_A
          yA = bndBoxJLow_A + DeltaJHalf_A
          do jA = jmin_A,jmax_A
             xA = bndBoxILow_A + DeltaIHalf_A
             do iA = imin_A,imax_A
!
!
!        Evaluate now the exact FLASH potential at cell cube (xA,yA,zA).
!
!
                potential_exact = ZERO 

                do blockNr_B = 1,gr_mpoleBlockCount
                   blockID_B = gr_mpoleBlockList (blockNr_B)

                   call Grid_getBlkBoundBox     (blockID_B,  bndBox_B)
                   call Grid_getDeltas          (blockID_B,  delta_B )
                   call Grid_getBlkPtr          (blockID_B,  solnData)
                   call Grid_getBlkIndexLimits  (blockID_B,  blkLimits_B, blkLimitsGC_B)

                   imin_B         = blkLimits_B (LOW, IAXIS)
                   jmin_B         = blkLimits_B (LOW, JAXIS)
                   kmin_B         = blkLimits_B (LOW, KAXIS)  
                   imax_B         = blkLimits_B (HIGH,IAXIS)
                   jmax_B         = blkLimits_B (HIGH,JAXIS)
                   kmax_B         = blkLimits_B (HIGH,KAXIS)
                   DeltaI_B      = delta_B (IAXIS)
                   DeltaJ_B      = delta_B (JAXIS)
                   DeltaK_B      = delta_B (KAXIS)
                   DeltaIHalf_B = DeltaI_B * HALF
                   DeltaJHalf_B = DeltaJ_B * HALF
                   DeltaKHalf_B = DeltaK_B * HALF
                   bndBoxILow_B = bndBox_B (LOW,IAXIS)
                   bndBoxJLow_B = bndBox_B (LOW,JAXIS)
                   bndBoxKLow_B = bndBox_B (LOW,KAXIS)

                   zB = bndBoxKLow_B + DeltaKHalf_B
                   do kB = kmin_B,kmax_B
                      yB = bndBoxJLow_B + DeltaJHalf_B
                      do jB = jmin_B,jmax_B
                         xB = bndBoxILow_B + DeltaIHalf_B
                         do iB = imin_B,imax_B

                            same_cube = xA.eq.xB .and. yA.eq.yB .and. zA.eq.zB

                            if (.not.same_cube) then
                                cellVolume = DeltaI_B * DeltaJ_B * DeltaK_B
                                cellMass = solnData (idensvar,iB,jB,kB) * cellVolume
                                r = sqrt ((xA-xB)*(xA-xB) + (yA-yB)*(yA-yB) + (zA-zB)*(zA-zB))
                                potential_exact = potential_exact + cellMass / r
                            end if

                            xB = xB + DeltaI_B
                         end do
                         yB = yB + DeltaJ_B
                      end do
                      zB = zB + DeltaK_B
                   end do

                   call Grid_releaseBlkPtr (blockID_B, solnData)

                end do

                potential_exact = - gr_mpoleGravityConstant * potential_exact
!
!
!        The entire exact FLASH potential at cell (iA,jA,kA) is ready.
!        Compare with the one evaluated using the multipoles.
!
!
                call Grid_getBlkPtr          (blockID_A, solnData)

                potential_mpole = solnData (ipotvar,iA,jA,kA)

                call Grid_releaseBlkPtr      (blockID_A, solnData)

                potdiff = abs (potential_exact - potential_mpole)
                if (potdiff > potdiff_max) then
!                    write (fileUnit,'(4(2X,I3,2X),E20.12,E20.12)') &
!                                      blockNr_A,iA,jA,kA,   potential_exact,potential_mpole
                    potdiff_max = potdiff
                end if

                write (fileUnit,'(4(2X,I3,2X),E20.12,E20.12)') &
                                   blockNr_A,iA,jA,kA,   potential_exact,potential_mpole

                xA = xA + DeltaI_A
             end do
             yA = yA + DeltaJ_A
          end do
          zA = zA + DeltaK_A
       end do

       write (fileUnit,*) ' maximum potential difference = ',potdiff_max

     case (GRID_2DCYLINDRICAL)

       zA = bndBoxJLow_A + DeltaJHalf_A
       do jA = jmin_A,jmax_A
          Rcyl_A = bndBoxILow_A + DeltaIHalf_A
          do iA = imin_A,imax_A
!
!
!        Evaluate now the potential at point (zA,RA). This 'point' represents,
!        in the 2D cylindrical case, actually a thin ring of crosssection:
!
!                                   -------
!                                  |       |
!                   O------------->|       | gr_mpoleDelta z
!                        Rcyl     |       |
!                                   -------
!                                 gr_mpoleDelta Rcyl
!
!        revolving in a circle around the coordinate origin 'O'.
!
!
             potential_exact = ZERO 

             do blockNr_B = 1,gr_mpoleBlockCount
                blockID_B = gr_mpoleBlockList (blockNr_B)

                call Grid_getBlkBoundBox     (blockID_B,  bndBox_B)
                call Grid_getDeltas          (blockID_B,  delta_B )
                call Grid_getBlkPtr          (blockID_B,  solnData)
                call Grid_getBlkIndexLimits  (blockID_B,  blkLimits_B, blkLimitsGC_B)

                imin_B         = blkLimits_B (LOW, IAXIS)
                jmin_B         = blkLimits_B (LOW, JAXIS)
                imax_B         = blkLimits_B (HIGH,IAXIS)
                jmax_B         = blkLimits_B (HIGH,JAXIS)
                DeltaI_B      = delta_B (IAXIS)
                DeltaJ_B      = delta_B (JAXIS)
                DeltaIHalf_B = DeltaI_B * HALF
                DeltaJHalf_B = DeltaJ_B * HALF
                bndBoxILow_B = bndBox_B (LOW,IAXIS)
                bndBoxJLow_B = bndBox_B (LOW,JAXIS)

                zB = bndBoxJLow_B + DeltaJHalf_B
                do jB = jmin_B,jmax_B
                   Rcyl_B = bndBoxILow_B + DeltaIHalf_B
                   do iB = imin_B,imax_B

                      same_ring = zA.eq.zB .and. Rcyl_A.eq.Rcyl_B

                      if (.not.same_ring) then
                          cellVolume = Rcyl_B * gr_mpoleTwoPi * DeltaI_B * DeltaJ_B
                          cellMass = solnData (idensvar,iB,jB,1) * cellVolume
                          r = sqrt ((zA-zB)*(zA-zB) + (Rcyl_A - Rcyl_B)*(Rcyl_A - Rcyl_B))
                          potential_exact = potential_exact + cellMass / r
                      end if

                      Rcyl_B = Rcyl_B + DeltaI_B
                   end do
                   zB = zB + DeltaJ_B
                end do

             end do

             potential_exact = - gr_mpoleGravityConstant * potential_exact
!
!
!        The entire exact FLASH potential at cell (iA,jA) is ready.
!        Compare with the one evaluated using the multipoles.
!
!
             potential_mpole = solnData (ipotvar,iA,jA,1)

             write (fileUnit,'(2(2X,I3,2X),7X,E20.12,E20.12)') &
                               iA,jA,   potential_exact,potential_mpole

             Rcyl_A = Rcyl_A + DeltaI_A
          end do
          zA = zA + DeltaJ_A
       end do

     case (GRID_1DSPHERICAL)

       Rsph_A = bndBoxILow_A + DeltaIHalf_A
       do iA = imin_A,imax_A
!
!
!        Evaluate now the potential at point (Rsph_A). This 'point' represents,
!        in the 1D spherical case, actually a thin spherical layer of thickness
!        (delta Rsph) with radius Rsph.
!
!
          potential_exact = ZERO 

          do blockNr_B = 1,gr_mpoleBlockCount
             blockID_B = gr_mpoleBlockList (blockNr_B)

             call Grid_getBlkBoundBox     (blockID_B,  bndBox_B)
             call Grid_getDeltas          (blockID_B,  delta_B )
             call Grid_getBlkPtr          (blockID_B,  solnData)
             call Grid_getBlkIndexLimits  (blockID_B,  blkLimits_B, blkLimitsGC_B)

             imin_B         = blkLimits_B (LOW, IAXIS)
             imax_B         = blkLimits_B (HIGH,IAXIS)
             DeltaI_B      = delta_B (IAXIS)
             DeltaIHalf_B = DeltaI_B * HALF
             bndBoxILow_B = bndBox_B (LOW,IAXIS)

             Rsph_B = bndBoxILow_B + DeltaIHalf_B
             do iB = imin_B,imax_B

                same_layer = Rsph_A .eq. Rsph_B

                if (.not.same_layer) then
                    cellVolume = Rsph_B * Rsph_B * gr_mpoleFourPi * DeltaI_B
                    cellMass = solnData (idensvar,iB,1,1) * cellVolume
                    potential_exact = potential_exact + cellMass / Rsph_B
                end if

                Rsph_B = Rsph_B + DeltaI_B
             end do

          end do

          potential_exact = - gr_mpoleGravityConstant * potential_exact
!
!
!        The entire exact FLASH potential at cell (iA) is ready.
!        Compare with the one evaluated using the multipoles.
!
!
          potential_mpole = solnData (ipotvar,iA,1,1)

          write (fileUnit,'((2X,I3,2X),2(7X),E20.12,E20.12)') &
                            iA,   potential_exact,potential_mpole

          Rsph_A = Rsph_A + DeltaI_A
       end do

     end select
!
!
!        Next leaf block.
!
!
  end do

  write (fileUnit,'(A38)') '  --- finished present iteration ---  '
  close (fileUnit)

  return
end subroutine gr_mpolePotential_exact
