!!****if* source/Grid/GridSolvers/Multipole_new/gr_mpolePot2Dcylindrical
!!
!! NAME
!!
!!  gr_mpolePot2Dcylindrical
!!
!! SYNOPSIS
!!
!!  gr_mpolePot2Dcylindrical  (integer, intent(in) :: ipotvar)
!!
!! DESCRIPTION
!!
!!  Computes the potential field for a 2D cylindrical geometry
!!  using the mass moments already calculated. On output the variable
!!  indexed by ipotvar contains the potential. The calculations are
!!  entirely local to each processor, since each processor has a local
!!  copy of the moments.
!!
!! ARGUMENTS
!!
!!  ipotvar  : index to variable containing the potential
!!
!!***

!!REORDER(4): solnData

subroutine gr_mpolePot2Dcylindrical (ipotvar)

  use Grid_interface,    ONLY : Grid_getBlkPtr,         &
                                Grid_releaseBlkPtr,     &
                                Grid_getBlkBoundBox,    &
                                Grid_getDeltas,         &
                                Grid_getBlkIndexLimits

  use gr_mpoleData,      ONLY : gr_mpoleGravityConstant,        &
                                gr_mpoleNumberInv,              &
                                gr_mpoleDrInv,                  &
                                gr_mpoleDrInnerZoneInv,         &
                                gr_mpoleMaxL,                   &
                                gr_mpoleMaxM,                   &
                                gr_mpoleMaxLM,                  &
                                gr_mpoleMaxQ,                   &
                                gr_mpoleMaxRadialZones,         &
                                gr_mpoleMinRadialZone,          &
                                gr_mpoleZoneRmax,               &
                                gr_mpoleZoneQmax,               &
                                gr_mpoleZoneType,               &
                                gr_mpoleZoneScalarInv,          &
                                gr_mpoleZoneLogNormInv,         &
                                gr_mpoleZoneExponentInv,        &
                                gr_mpoleInnerZoneMaxR,          &
                                gr_mpoleInnerZoneDrRadii,       &
                                gr_mpoleInnerZoneQlower,        &
                                gr_mpoleInnerZoneQupper,        &
                                gr_mpoleInnerZoneResolution,    &
                                gr_mpoleInnerZoneResolutionInv, &
                                gr_mpoleOuterZoneQshift

  use gr_mpoleData,      ONLY : gr_mpoleRcenter,                &
                                gr_mpoleZcenter,                &
                                gr_mpoleQDampingR,              &
                                gr_mpoleQDampingI,              &
                                gr_mpoleMomentR,                &
                                gr_mpoleMomentI,                &
                                gr_mpoleBlockCount,             &
                                gr_mpoleBlockList

  implicit none
  
#include "Flash.h"
#include "constants.h"
#include "gr_mpole.h"
  
  integer, intent (in) :: ipotvar

  logical :: i2, j2
  logical :: innerZonePotential

  integer :: blockNr, blockID
  integer :: DrUnit
  integer :: imax, imin, iC, iCmax, iF, iFmax, iB
  integer :: jmax, jmin, jC, jCmax, jF, jFmax
  integer :: L
  integer :: Q, Qlocal, Qlower, Qupper
  integer :: type
  integer :: zone

  integer :: blkLimits   (LOW:HIGH,1:MDIM)
  integer :: blkLimitsGC (LOW:HIGH,1:MDIM)

  real    :: bndBoxILow, bndBoxJLow
  real    :: dampI, Idamp
  real    :: dampR, Rdamp
  real    :: DeltaI, DeltaIHalf
  real    :: DeltaJ, DeltaJHalf
  real    :: f,g,h
  real    :: facePotential
  real    :: Ic0, Ic1, Ic2, IcL
  real    :: Qfloat, QfracI, QfracR
  real    :: r, rR, rI, rlocal, rinDrs, rsqR, rinvI, rsqinvI
  real    :: Rc0, Rc1, Rc2, RcL
  real    :: Rcyl
  real    :: Rdamping, Idamping
  real    :: RdampingQuotient, IdampingQuotient
  real    :: RdotI, IdotR
  real    :: sclInv, lgnInv, expInv
  real    :: z, zR, zI

  real    :: delta           (1:MDIM)
  real    :: bndBox (LOW:HIGH,1:MDIM)

  real, pointer   :: solnData (:,:,:,:)
!  
!
!     ...Sum quantities over all locally held leaf blocks.
!
!
!$omp do schedule (static)
  do blockNr = 1,gr_mpoleBlockCount

     blockID = gr_mpoleBlockList (blockNr)

     call Grid_getBlkBoundBox     (blockID, bndBox)
     call Grid_getDeltas          (blockID, delta)
     call Grid_getBlkPtr          (blockID, solnData)
     call Grid_getBlkIndexLimits  (blockID, blkLimits, blkLimitsGC)

     imin       = blkLimits (LOW, IAXIS)
     jmin       = blkLimits (LOW, JAXIS)
     imax       = blkLimits (HIGH,IAXIS)
     jmax       = blkLimits (HIGH,JAXIS)

     iCmax      = imax - imin + 1      ! # of local (inner) cells in i direction
     jCmax      = jmax - jmin + 1      ! # of local (inner) cells in j direction
     iFmax      = iCmax + iCmax        !  max local (inner) face index in i direction (1st face -> index 0)
     jFmax      = jCmax + jCmax        !  max local (inner) face index in j direction (1st face -> index 0)

     DeltaI     = delta (IAXIS)
     DeltaJ     = delta (JAXIS)
     DeltaIHalf = DeltaI * HALF
     DeltaJHalf = DeltaJ * HALF

     bndBoxILow = bndBox (LOW,IAXIS)
     bndBoxJLow = bndBox (LOW,JAXIS)

     solnData (ipotvar , imin:imax , jmin:jmax , 1) = ZERO
!
!
!         ...The 2D cylindrical case:
!
!
!                               ------
!                             /        \
!                            /     z    \
!                           |\     |    /|
!                           | \    |   / |
!                           |   ------   |
!                           |      |     |
!                           |      |     |         Rcyl --> stored in i-index (FLASH x)
!                           |      ----->|            z --> stored in j-index (FLASH y)
!                           |       Rcyl |
!                           |            |
!                           |            |
!                           |   ------   |
!                           | /        \ |
!                           |/          \|
!                            \          /
!                             \        /
!                               ------
!
!            Here the convention used in FLASH is to store the Rcyl values into the
!            i-index and the z values into the j-index. As for the 3D axisymmetric
!            case, only M=0 moments have been calculated and thus only M=0 regular
!            and irregular solid harmonics:
!
!                          0 =< L =< gr_mpoleMaxL
!                               M = 0 (cosine part)
!
!            are calculated and the appropriate scalar products with the moments
!            formed. For details on how z and Rcyl relate to the FLASH indices,
!            please see the moment evaluation routine.
!
!            Explanation of 'QfracR' and 'QfracI':
!            ------------------------------------
!
!            For outer zone (statistical) face potentials, we have to extract the
!            proper fractions of the Q-bin moments. The following picture, showing
!            the relevant radial bins, clarifies what is being done:
!
!
!
!                                            Q-th bin
!
!              ----- MR(Q-1) ---- | ---- MR(Q)----- MI(Q) -----| ------ MI(Q+1)----
!                                                 |
!                                          QfracR | QfracI
!                                                 |
!                                              Qfloat
!
!
!            Here 'Qfloat' is the value obtained by the cell's face radial position
!            and from 'Qfloat' we determined the bin number 'Q'. The values of
!            'QfracR' and 'QfracI' denote the fractional location of 'Qfloat'
!            within the Q-th bin. Hence we always have: QfracR + QfracI = 1.
!            When forming the dot products between the current Q-th bin solid
!            harmonics and the Moments, we take only the relevant fraction of
!            the Moments in Q. This sees the Moments in Q as a kind of statistical
!            average in terms of the radii comprised within the Q-th bin. As an
!            example, we take the dot product between the regular solid harmonics
!            in the Q-th bin R(Q) with all the necessary irregular Moments MI.
!            We have:
!
!
!              RdotI = R(Q) dot { MI(Q+1) + QfracI * [MI(Q) - MI(Q+1)] }
!
!                    = R(Q) dot { [QfracR + QfracI] * MI(Q+1) + QfracI * [MI(Q) - MI(Q+1)] }
!
!                    = R(Q) dot { QfracR * I(Q+1) + QfracI * I(Q) }
!
!
!            The same applies for the dot product between the irregular solid
!            harmonics in the Q-th bin I(Q) and the regular Moments MR. We get:
!
!
!              IdotR = I(Q) dot { QfracI * MR(Q-1) + QfracR * MR(Q) }
!
!
!            For the inner zone face potentials, we combine the regular face moments
!            with the irregular Q-bin moments and the irregular face moments with the
!            regular (Q-1)-bin moments. This can readily be achieved by setting
!            QfracR = 0 and QfracI = 1.
!
!
!
     z = bndBoxJLow - gr_mpoleZcenter                             ! initial z-axis location of j face

     do jF = 0, jFmax                                             ! loop over local j face indices

      jC = int (jF/2) + 1                                         ! local (inner) largest j cell index for j face
      j2 = mod (jF,2) == 0 .and. jC > 1 .and. jC < jCmax + 1      ! 2 cells j and j-1 share j face?
      jC = jmin - 1 + min (jC , jCmax)                            ! change to global (inner + guard) j cell index
      iB = mod (jF+1,2)                                           ! initial local i face index (0 or 1)

      Rcyl = bndBoxILow + iB * DeltaIHalf - gr_mpoleRcenter       ! initial Rcyl-axis location of i face

      do iF = iB, iFmax, 2                                        ! loop over local i face indices

       iC = int (iF/2) + 1                                        ! local (inner) largest i cell index for i face
       i2 = mod (iF,2) == 0 .and. iC > 1 .and. iC < iCmax + 1     ! 2 cells i and i-1 share i face?
       iC = imin - 1 + min (iC , iCmax)                           ! change to global (inner + guard) i cell index

       r = sqrt (z * z + Rcyl * Rcyl)                             ! radial position of face location
!
!
!        ...Find the radial bin.
!
!
       innerZonePotential = r <= gr_mpoleInnerZoneMaxR

       if (innerZonePotential) then

           rinDrs = r * gr_mpoleDrInnerZoneInv
           DrUnit = int (ceiling (rinDrs))
           Qlower = gr_mpoleInnerZoneQlower (DrUnit)
           Qupper = gr_mpoleInnerZoneQupper (DrUnit)
           QfracR = ZERO
           QfracI = ONE

           do Q = Qlower,Qupper
              if (rinDrs <= gr_mpoleInnerZoneDrRadii (Q)) exit
           end do

       else

           do zone = gr_mpoleMinRadialZone, gr_mpoleMaxRadialZones
              if (r - gr_mpoleZoneRmax (zone) <= ZERO) exit
           end do

           rlocal = r - gr_mpoleZoneRmax    (zone - 1)
           type   = gr_mpoleZoneType        (zone)
           sclInv = gr_mpoleZoneScalarInv   (zone)
           expInv = gr_mpoleZoneExponentInv (zone)

           if (type == ZONE_EXPONENTIAL) then
               Qfloat = (rlocal * sclInv * gr_mpoleDrInv) ** expInv
           else if (type == ZONE_LOGARITHMIC) then
               lgnInv = gr_mpoleZoneLogNormInv (zone)
               Qfloat = expInv * log (rlocal * sclInv * gr_mpoleDrInv * lgnInv + ONE)
           end if

           Qlocal = ceiling (Qfloat)
           QfracI = real (Qlocal) - Qfloat
           QfracR = ONE - QfracI
           Q      = gr_mpoleZoneQmax (zone - 1) + Qlocal + gr_mpoleOuterZoneQshift

       end if
!
!
!        ...Calculate the (damped) moments.
!
!
       Rdamping         = gr_mpoleQDampingR (Q)
       Idamping         = gr_mpoleQDampingI (Q)
       IdampingQuotient = gr_mpoleQDampingI (Q+1) / gr_mpoleQDampingI (Q)

       rI = r * Rdamping                  ! irregular solid harmonics will cancel Rdamping
       rinvI = ONE / rI

       RcL = Idamping                     ! start undo 1/(Idamping)^(L+1) in irregular Q moments
       IcL = rinvI * Rdamping             ! start undo (Rdamping)^L in regular Q moments

       Idamp = IdampingQuotient           ! damping for irregular moments for L = 0
       Rdamp = ONE                        ! damping for   regular moments for L = 0

       RdotI = RcL * (QfracI * gr_mpoleMomentI (1,Q) + QfracR * Idamp * gr_mpoleMomentI (1,Q+1))
       IdotR = IcL * (QfracR * gr_mpoleMomentR (1,Q) + QfracI         * gr_mpoleMomentR (1,Q-1))

       facePotential = RdotI + IdotR      ! initialize face potential variable
!
!
!        ...Proceed, if we go beyond monopoles.
!
!
       if (gr_mpoleMaxL > 0) then

        RdampingQuotient = gr_mpoleQDampingR (Q) / gr_mpoleQDampingR (Q-1)

        zR = z * Idamping                 !  regular  solid harmonics will cancel Idamping
        rR = r * Idamping                 !  regular  solid harmonics will cancel Idamping
        zI = z * Rdamping                 ! irregular solid harmonics will cancel Rdamping

        rsqR = rR * rR
        rsqinvI = rinvI * rinvI

        Rc1 = RcL * zR                    !  regular  solid harmonic for cos 10
        Ic1 = IcL * zI * rsqinvI          ! irregular solid harmonic for cos 10

        Idamp = Idamp * IdampingQuotient  ! store damping for irregular moments for L = 1
        Rdamp = Rdamp * RdampingQuotient  ! store damping for   regular moments for L = 1
        dampI = Idamp                     ! initialize for L accumulation damping
        dampR = Rdamp                     ! initialize for L accumulation damping

        RdotI = Rc1 * (QfracI * gr_mpoleMomentI (2,Q) + QfracR * dampI * gr_mpoleMomentI (2,Q+1))
        IdotR = Ic1 * (QfracR * gr_mpoleMomentR (2,Q) + QfracI * dampR * gr_mpoleMomentR (2,Q-1))

        facePotential = facePotential + RdotI + IdotR
!
!
!        ...Proceed, if we go beyond dipoles.
!
!
        if (gr_mpoleMaxL > 1) then

         Rc0 = RcL                                       !  regular  solid harmonic for cos 00
         Ic0 = IcL                                       ! irregular solid harmonic for cos 00

         do L = 2,gr_mpoleMaxL

            h = real (L + L - 1)
            g = real ((L - 1) * (L - 1))
            f = gr_mpoleNumberInv (L) * gr_mpoleNumberInv (L)

            Rc2 = (h * zR * Rc1 - rsqR * Rc0) * f        !  regular  solid harmonic for cos L0
            Ic2 = (h * zI * Ic1 -    g * Ic0) * rsqinvI  ! irregular solid harmonic for cos L0

            dampI = dampI * IdampingQuotient             ! next L damping
            dampR = dampR * RdampingQuotient             ! next L damping

            RdotI = Rc2 * (QfracI * gr_mpoleMomentI (L+1,Q) + QfracR * dampI * gr_mpoleMomentI (L+1,Q+1))
            IdotR = Ic2 * (QfracR * gr_mpoleMomentR (L+1,Q) + QfracI * dampR * gr_mpoleMomentR (L+1,Q-1))

            facePotential = facePotential + RdotI + IdotR

            Rc0 = Rc1                                    !
            Rc1 = Rc2                                    ! swap, to accumulate next L regular
            Ic0 = Ic1                                    ! and irregular solid harmonic
            Ic1 = Ic2                                    !

         end do

        end if    ! gr_mpoleMaxL > 1 condition
       end if     ! gr_mpoleMaxL > 0 condition

       facePotential = - gr_mpoleGravityConstant * facePotential
!
!
!        ...Add the current face potential to the relevant cell(s) of the potential block.
!           Note, that a face can only be shared by up to 2 cells and no more. The most
!           common case (shared on the i-index) is tested first.
!
!
       if (i2) then

           solnData (ipotvar,iC-1,jC,1) = solnData (ipotvar,iC-1,jC,1) + facePotential
           solnData (ipotvar,iC  ,jC,1) = solnData (ipotvar,iC  ,jC,1) + facePotential

       else if (j2) then

           solnData (ipotvar,iC,jC-1,1) = solnData (ipotvar,iC,jC-1,1) + facePotential
           solnData (ipotvar,iC,jC  ,1) = solnData (ipotvar,iC,jC  ,1) + facePotential

       else

           solnData (ipotvar,iC,jC,1) = solnData (ipotvar,iC,jC,1) + facePotential

       end if

       Rcyl = Rcyl + DeltaI      ! face position increment along radial axis
      end do
      z = z + DeltaJHalf         ! face position increment along z-axis
     end do
!
!
!    ...Form the potential average in each cell.
!
!
     solnData (ipotvar,imin:imax,jmin:jmax,1) = FOURTH * solnData (ipotvar,imin:imax,jmin:jmax,1)
!
!
!    ...Get ready for retrieving next LEAF block for the current processor.
!
!
     call Grid_releaseBlkPtr (blockID, solnData)

  end do
!$omp end do
!
!
!    ...Ready!
!
!
  return
end subroutine gr_mpolePot2Dcylindrical
