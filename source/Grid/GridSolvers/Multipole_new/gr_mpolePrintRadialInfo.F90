!!****if* source/Grid/GridSolvers/Multipole_new/gr_mpolePrintRadialInfo
!!
!! NAME
!!
!!  gr_mpolePrintRadialInfo
!!
!! SYNOPSIS
!!
!!  gr_mpolePrintRadialInfo ()
!!
!! DESCRIPTION
!!
!!  Utility routine, which prints detailed info regarding the established radial
!!  bin grid to a text file. The information is written out to a file named
!!  <basenm>RadialInfoPrint.txt, where <basenm> is the runtime parameter for
!!  output file names. The file is appended at each time for each iteration.
!!
!!***

subroutine gr_mpolePrintRadialInfo ()

  use Grid_data,                   ONLY : gr_meshMe

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  use gr_mpoleData,                ONLY : gr_mpoleGeometry,            &
                                          gr_mpoleQDampingR,           &
                                          gr_mpoleQDampingI,           &
                                          gr_mpoleQRadii,              &
                                          gr_mpoleDr,                  &
                                          gr_mpoleDrInnerZone,         &
                                          gr_mpoleMaxR,                &
                                          gr_mpoleMaxQ,                &
                                          gr_mpoleMaxRadialZones,      &
                                          gr_mpoleMinRadialZone,       &
                                          gr_mpoleZoneQmax,            &
                                          gr_mpoleZoneType,            &
                                          gr_mpoleInnerZoneExists,     &
                                          gr_mpoleInnerZoneMaxR,       &
                                          gr_mpoleInnerZoneQmax,       &
                                          gr_mpoleInnerZoneDrRadii,    &
                                          gr_mpoleInnerZoneQlower,     &
                                          gr_mpoleInnerZoneQupper,     &
                                          gr_mpoleInnerZoneSize,       &
                                          gr_mpoleInnerZoneResolution, &
                                          gr_mpoleOuterZoneExists,     &
                                          gr_mpoleOuterZoneQshift

  use gr_mpoleData,                ONLY : gr_mpoleRad2Deg,             &
                                          gr_mpoleXcenter,             &
                                          gr_mpoleYcenter,             &
                                          gr_mpoleZcenter,             &
                                          gr_mpoleRcenter,             &
                                          gr_mpolePhiCenter,           &
                                          gr_mpoleThetaCenter,         &
                                          gr_mpoleDomainXmin,          &
                                          gr_mpoleDomainYmin,          &
                                          gr_mpoleDomainZmin,          &
                                          gr_mpoleDomainRmin,          &
                                          gr_mpoleDomainPhiMin,        &
                                          gr_mpoleDomainThetaMin,      &
                                          gr_mpoleDomainXmax,          &
                                          gr_mpoleDomainYmax,          &
                                          gr_mpoleDomainZmax,          &
                                          gr_mpoleDomainRmax,          &
                                          gr_mpoleDomainPhiMax,        &
                                          gr_mpoleDomainThetaMax

  implicit none

#include "constants.h"
#include "gr_mpole.h"
   
  interface
     integer function ut_getFreeFileUnit()
     end function ut_getFreeFileUnit
  end interface

  character*11 :: type

  logical, save :: firstCall = .true.
  logical :: firstLine

  integer :: fileUnit
  integer :: n
  integer :: posBlank
  integer :: Q,Qlocal,QlocalMax
  integer :: zone

  real    :: rglobal

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
          fileName = baseName (:posBlank-1) // 'RadialInfoPrint.txt'
          open (fileUnit, file=fileName)
          firstCall = .false.
      else
          open (fileUnit, file=fileName, position='APPEND')
      end if
!
!
!     ...Printout multipole center of expansion info.
!
!
      write (fileUnit,*)
      write (fileUnit,*) '     MULTIPOLE CENTER OF EXPANSION LOCATION'
      write (fileUnit,*)

      select case (gr_mpoleGeometry)

        case (GRID_3DCARTESIAN)

          write (fileUnit,'(A,ES20.12)') ' Center X     = ',gr_mpoleXcenter
          write (fileUnit,'(A,ES20.12)') ' Center Y     = ',gr_mpoleYcenter
          write (fileUnit,'(A,ES20.12)') ' Center Z     = ',gr_mpoleZcenter
          write (fileUnit,'(A,ES20.12)') ' Domain Min X = ',gr_mpoleDomainXmin
          write (fileUnit,'(A,ES20.12)') ' Domain Min Y = ',gr_mpoleDomainYmin
          write (fileUnit,'(A,ES20.12)') ' Domain Min Z = ',gr_mpoleDomainZmin
          write (fileUnit,'(A,ES20.12)') ' Domain Max X = ',gr_mpoleDomainXmax
          write (fileUnit,'(A,ES20.12)') ' Domain Max Y = ',gr_mpoleDomainYmax
          write (fileUnit,'(A,ES20.12)') ' Domain Max Z = ',gr_mpoleDomainZmax

        case (GRID_3DCYLINDRICAL)

          write (fileUnit,'(A,ES20.12)') ' Center R         = ',gr_mpoleRcenter
          write (fileUnit,'(A,ES20.12)') ' Center Z         = ',gr_mpoleZcenter
          write (fileUnit,'(A,ES20.12)') ' Center Phi (deg) = ',gr_mpolePhiCenter    * gr_mpoleRad2Deg
          write (fileUnit,'(A,ES20.12)') ' Domain Min R     = ',gr_mpoleDomainRmin
          write (fileUnit,'(A,ES20.12)') ' Domain Min Z     = ',gr_mpoleDomainZmin
          write (fileUnit,'(A,ES20.12)') ' Domain Min Phi   = ',gr_mpoleDomainPhiMin * gr_mpoleRad2Deg
          write (fileUnit,'(A,ES20.12)') ' Domain Max R     = ',gr_mpoleDomainRmax
          write (fileUnit,'(A,ES20.12)') ' Domain Max Z     = ',gr_mpoleDomainZmax
          write (fileUnit,'(A,ES20.12)') ' Domain Max Phi   = ',gr_mpoleDomainPhiMax * gr_mpoleRad2Deg

        case (GRID_2DCYLINDRICAL)

          write (fileUnit,'(A,ES20.12)') ' Center R     = ',gr_mpoleRcenter
          write (fileUnit,'(A,ES20.12)') ' Center Z     = ',gr_mpoleZCenter
          write (fileUnit,'(A,ES20.12)') ' Domain Min R = ',gr_mpoleDomainRmin
          write (fileUnit,'(A,ES20.12)') ' Domain Min Z = ',gr_mpoleDomainZmin
          write (fileUnit,'(A,ES20.12)') ' Domain Max R = ',gr_mpoleDomainRmax
          write (fileUnit,'(A,ES20.12)') ' Domain Max Z = ',gr_mpoleDomainZmax

        case (GRID_2DSPHERICAL)

          write (fileUnit,'(A,ES20.12)') ' Center R           = ',gr_mpoleRcenter
          write (fileUnit,'(A,ES20.12)') ' Center Theta (deg) = ',gr_mpoleThetaCenter    * gr_mpoleRad2Deg
          write (fileUnit,'(A,ES20.12)') ' Domain Min R       = ',gr_mpoleDomainRmin
          write (fileUnit,'(A,ES20.12)') ' Domain Min Theta   = ',gr_mpoleDomainThetaMin * gr_mpoleRad2Deg
          write (fileUnit,'(A,ES20.12)') ' Domain Max R       = ',gr_mpoleDomainRmax
          write (fileUnit,'(A,ES20.12)') ' Domain Max Theta   = ',gr_mpoleDomainThetaMax * gr_mpoleRad2Deg

        case (GRID_1DSPHERICAL)

          write (fileUnit,'(A,ES20.12)') ' Center R (domain origin) = ',ZERO
          write (fileUnit,'(A,ES20.12)') ' Domain Min R             = ',gr_mpoleDomainRmin
          write (fileUnit,'(A,ES20.12)') ' Domain Max R             = ',gr_mpoleDomainRmax

      end select
!
!
!     ...Printout the radial info.
!
!
      write (fileUnit,*)
      write (fileUnit,*) '               ZONE GRID INFO'
      write (fileUnit,*)
      write (fileUnit,'(A,ES20.12)') ' gr_mpoleDr              = ',gr_mpoleDr
      write (fileUnit,'(A,ES20.12)') ' gr_mpoleDrInnerZone     = ',gr_mpoleDrInnerZone
      write (fileUnit,'(A,ES20.12)') ' gr_mpoleMaxR            = ',gr_mpoleMaxR
      write (fileUnit,'(A,ES20.12)') ' gr_mpoleInnerZoneMaxR   = ',gr_mpoleInnerZoneMaxR
      write (fileUnit,'(A,I10)')     ' gr_mpoleMaxQ            = ',gr_mpoleMaxQ
      write (fileUnit,'(A,I10)')     ' gr_mpoleInnerZoneQmax   = ',gr_mpoleInnerZoneQmax
      write (fileUnit,'(A,I10)')     ' gr_mpoleInnerZoneSize   = ',gr_mpoleInnerZoneSize
      write (fileUnit,'(A,I10)')     ' gr_mpoleMinRadialZone   = ',gr_mpoleMinRadialZone
      write (fileUnit,'(A,I10)')     ' gr_mpoleMaxRadialZones  = ',gr_mpoleMaxRadialZones
      write (fileUnit,'(A,I10)')     ' gr_mpoleOuterZoneQshift = ',gr_mpoleOuterZoneQshift
!
!
!     ...Inner zone radial info (if present).
!
!
      if (gr_mpoleInnerZoneExists) then

          write (fileUnit,*)
          write (fileUnit,*)
          write (fileUnit,*)
          write (fileUnit,*) '               INNER ZONE QUENTCHING ARRAYS '
          write (fileUnit,*)
          write (fileUnit,'(A)') '   # Dr InnerZone         Q lower   Q upper    '
          write (fileUnit,'(A)') ' ----------------------------------------------'

          do n = 1,gr_mpoleInnerZoneSize
             write (fileUnit,'(8X,I4,                      16X,I4,                      6X,I4)') &
                                   n, gr_mpoleInnerZoneQlower (n), gr_mpoleInnerZoneQupper (n)
          end do

          write (fileUnit,*)
          write (fileUnit,*)
          write (fileUnit,*)
          write (fileUnit,*) '               INNER ZONE RADII'
          write (fileUnit,*)
          write (fileUnit,'(A)') '   Q #     R (Dr InnerZone)       R          DampingR       DampingI     '
          write (fileUnit,'(A)') ' ------------------------------------------------------------------------'

          do Q = 1,gr_mpoleInnerZoneQmax
             write (fileUnit,'(2X,I4,  6X,ES13.6,  2X,ES13.6,  2X,ES13.6,  2X,ES13.6)')      &
                               Q, gr_mpoleInnerZoneDrRadii (Q),                              &
                                  gr_mpoleQRadii           (Q),                              &
                                  gr_mpoleQDampingR        (Q),                              &
                                  gr_mpoleQDampingI        (Q)
          end do

      end if
!
!
!     ...Outer zone radial info.
!
!
      if (gr_mpoleOuterZoneExists) then

          write (fileUnit,*)
          write (fileUnit,*)
          write (fileUnit,*)
          write (fileUnit,*) '               OUTER (statistical) ZONE'
          write (fileUnit,*)
          write (fileUnit,'(A)') '   zone #     type           Q         R          DampingR        DampingI '
          write (fileUnit,'(A)') ' --------------------------------------------------------------------------'

          do zone = gr_mpoleMinRadialZone, gr_mpoleMaxRadialZones

             QlocalMax = gr_mpoleZoneQmax (zone) - gr_mpoleZoneQmax (zone-1)
             firstLine = .true.

             if (gr_mpoleZoneType (zone) == ZONE_EXPONENTIAL) then

                 type = 'exponential'

                 do Qlocal = 1,QlocalMax

                    Q       = gr_mpoleZoneQmax (zone-1) + Qlocal + gr_mpoleOuterZoneQshift
                    rglobal = gr_mpoleQRadii (Q)

                    if (rglobal > gr_mpoleInnerZoneMaxR) then

                        if (firstLine) then
                            firstLine = .false.
                            write (fileUnit,'(4X,I2,  4X,A11,  1X,I8,  2X,ES13.6,  2X,ES13.6,  2X,ES13.6)') &
                                              zone, type, Q, rglobal, gr_mpoleQDampingR (Q), gr_mpoleQDampingI (Q)
                        else
                            write (fileUnit,'(4X,2X,  4X,11X,  1X,I8,  2X,ES13.6,  2X,ES13.6,  2X,ES13.6 )') &
                                              Q, rglobal, gr_mpoleQDampingR (Q), gr_mpoleQDampingI (Q)
                        end if
                    end if
                 end do

             else if (gr_mpoleZoneType (zone) == ZONE_LOGARITHMIC) then

                 type = 'logarithmic'

                 do Qlocal = 1,QlocalMax

                    Q       = gr_mpoleZoneQmax (zone-1) + Qlocal + gr_mpoleOuterZoneQshift
                    rglobal = gr_mpoleQRadii (Q)

                    if (rglobal > gr_mpoleInnerZoneMaxR) then

                        if (firstLine) then
                            firstLine = .false.
                            write (fileUnit,'(4X,I2,  4X,A11,  1X,I8,  2X,ES13.6,  2X,ES13.6,  2X,ES13.6)') &
                                              zone, type, Q, rglobal, gr_mpoleQDampingR (Q), gr_mpoleQDampingI (Q)
                        else
                            write (fileUnit,'(4X,2X,  4X,11X,  1X,I8,  2X,ES13.6,  2X,ES13.6,  2X,ES13.6 )') &
                                              Q, rglobal, gr_mpoleQDampingR (Q), gr_mpoleQDampingI (Q)
                        end if
                    end if
                 end do

             end if

             write (fileUnit,'(A)') ' --------------------------------------------------------------------------'

          end do

      end if

      write (fileUnit,*)
      write (fileUnit,'(A)') '  ### finished present iteration ###  '
      close (fileUnit)

  end if
!
!
!    ...Ready!
!
!  
  return
end subroutine gr_mpolePrintRadialInfo


