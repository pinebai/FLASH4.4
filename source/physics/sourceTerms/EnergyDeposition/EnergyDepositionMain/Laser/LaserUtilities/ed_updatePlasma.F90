!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserUtilities/ed_updatePlasma
!!
!! NAME
!!
!!  ed_updatePlasma
!!
!! SYNOPSIS
!!
!!  call ed_updatePlasma (integer (in)           :: blockCount, 
!!                        integer (in)           :: blockList (:),
!!                        real    (in), optional :: scaleFact)
!!
!! DESCRIPTION
!!
!!  Updates the physical variables of all cells on all blocks owned by the
!!  current processor. This is done by calling Eos.
!!
!! ARGUMENTS
!!
!!  blockCount : number of blocks on current processor
!!  blockList  : all block ID numbers
!!  scaleFact  : optional scaling factor
!!
!! NOTES
!!
!!***

subroutine ed_updatePlasma (blockCount,blockList,scaleFact)

  use Eos_interface,  ONLY : Eos_wrapped

  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
                             Grid_getBlkPtr,         &
                             Grid_releaseBlkPtr

  use EnergyDeposition_data,  ONLY : ed_depoVar, ed_depoVarIsPerMass

  implicit none
  
#include "constants.h"
#include "Flash.h"
#include "EnergyDeposition.h"

  integer,           intent (in) :: blockCount
  integer,           intent (in) :: blockList (1:blockCount)
  real   , optional, intent (in) :: scaleFact

  integer :: block
  integer :: blockID
  integer :: ib,ie,jb,je,kb,ke

  integer :: blkLimits   (LOW:HIGH,MDIM)
  integer :: blkLimitsGC (LOW:HIGH,MDIM)

  real, pointer :: solnData (:,:,:,:)
!
!
!    ...loop over all blocks in current processor.
!
!
  do block = 1, blockCount

     blockID = blockList (block)

     call Grid_getBlkPtr         (blockID, solnData, CENTER)
     call Grid_getBlkIndexLimits (blockID, blkLimits, blkLimitsGC, CENTER)

     ib = blkLimits (LOW ,IAXIS)
     ie = blkLimits (HIGH,IAXIS)
     jb = blkLimits (LOW ,JAXIS)
     je = blkLimits (HIGH,JAXIS)
     kb = blkLimits (LOW ,KAXIS)
     ke = blkLimits (HIGH,KAXIS)


     if (ed_depoVarIsPerMass) then
        if (present(scaleFact)) then
           solnData (EELE_VAR,ib:ie,jb:je,kb:ke) =   solnData (EELE_VAR,  ib:ie,jb:je,kb:ke) &
                                       + scaleFact * solnData (ed_depoVar,ib:ie,jb:je,kb:ke)
        else
           solnData (EELE_VAR,ib:ie,jb:je,kb:ke) =   solnData (EELE_VAR,  ib:ie,jb:je,kb:ke) &
                                                   + solnData (ed_depoVar,ib:ie,jb:je,kb:ke)
        end if
     else
        if (present(scaleFact)) then
           where ((          solnData(DENS_VAR,  ib:ie,jb:je,kb:ke).NE.0.0) .OR. &
                  (scaleFact*solnData(ed_depoVar,ib:ie,jb:je,kb:ke).NE.0.0))
              solnData (EELE_VAR,ib:ie,jb:je,kb:ke) =   solnData (EELE_VAR,  ib:ie,jb:je,kb:ke) &
                                       + (scaleFact *   solnData (ed_depoVar,ib:ie,jb:je,kb:ke) &
                                                      / solnData (DENS_VAR,  ib:ie,jb:je,kb:ke))
           end where
        else
           where ((          solnData(DENS_VAR,  ib:ie,jb:je,kb:ke).NE.0.0) .OR. &
                  (          solnData(ed_depoVar,ib:ie,jb:je,kb:ke).NE.0.0))
              solnData (EELE_VAR,ib:ie,jb:je,kb:ke) =   solnData (EELE_VAR,  ib:ie,jb:je,kb:ke) &
                                                   + (  solnData (ed_depoVar,ib:ie,jb:je,kb:ke) &
                                                      / solnData (DENS_VAR,  ib:ie,jb:je,kb:ke))
           end where
        end if
     end if

     call Grid_releaseBlkPtr (blockID, solnData, CENTER)

     call Eos_wrapped  (MODE_DENS_EI_GATHER, blkLimits, blockID)

  enddo
!
!
!    ...Ready!
!
!
  return
end subroutine ed_updatePlasma
