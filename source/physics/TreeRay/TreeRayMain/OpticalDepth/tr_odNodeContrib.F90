!!****if* source/physics/TreeRay/TreeRayMain/OpticalDepth/tr_odNodeContrib
!!
!! NAME
!!
!!  tr_odNodeContrib
!!
!!
!! SYNOPSIS
!!
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!
!! RESULT
!!
!!
!!***

subroutine tr_odNodeContrib(node, trLevel, refLevel, ii, jj, kk, ins, iph, ith, dr)
  use TreeRay_data, ONLY : tr_ilNI, tr_intersectList, tr_bhIM, tr_bhIBM, &
    tr_bhTreeLevels, tr_bhCdMaps, tr_nPo4pi, tr_bhTreeNodeSize
  use tr_odData, ONLY : tr_odIH2, tr_odICO, tr_odIBH2, tr_odIBCO, &
    tr_odICDTO, tr_odICDH2, tr_odICDCO, tr_odCDTOIndex
  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
  real, dimension(:), intent(IN) :: node
  integer, intent(IN) :: trLevel, refLevel, ins, iph, ith
  integer, intent(IN) :: ii, jj, kk
  real, dimension(MDIM+2), intent(IN) :: dr
  integer :: i, ipix, ndLevel
  real :: columnto, columnH2, columnCO, weight, pixarea_i
  real :: ndSize, ndVol_i

  ! inverted area of the healpix beam at distance to which the node projects
  ! nPix / (4*pi*dr^2)
  ndLevel = trLevel + refLevel
  pixarea_i = tr_nPo4pi * dr(MDIM+2)**2
  ndSize = tr_bhTreeNodeSize(ndLevel)
  ndVol_i = 1.0/ndSize**3
  !print *, "NC1 = ", tr_nPo4pi, dr, node

  ! node column densities
  if (trLevel == tr_bhTreeLevels) then
    !columnto = node(tr_bhIBM) * pixarea_i
    columnto = node(tr_bhIBM) * (node(tr_bhIBM)*ndVol_i)**(tr_odCDTOIndex-1) * pixarea_i
#if CHEMISTRYNETWORK == 4 || CHEMISTRYNETWORK == 5
    columnH2 = node(tr_odIBH2) * pixarea_i
#endif
#if CHEMISTRYNETWORK == 5
    columnCO = node(tr_odIBCO) * pixarea_i
#endif
  else
    !columnto = node(tr_bhIM) * pixarea_i
    columnto = node(tr_bhIM)*(node(tr_bhIBM)*ndVol_i)**(tr_odCDTOIndex-1) * pixarea_i
    !print *, "NC bnode = ", node, pixarea_i, columnto
#if CHEMISTRYNETWORK == 4 || CHEMISTRYNETWORK == 5
    columnH2 = node(tr_odIH2) * pixarea_i
#endif
#if CHEMISTRYNETWORK == 5
    columnCO = node(tr_odICO) * pixarea_i
#endif
  endif

  !print *, "NC2: ", columnto, columnH2, columnCO, pixarea_i

  do i = 1, tr_ilNI
    ! check if the node intersect with the ray
    if (tr_intersectList(i,ins,iph,ith) < 0) exit

    ! determine the ray index (ipix) and the weight of the intersection
    ipix = floor(tr_intersectList(i,ins,iph,ith))
    weight = (tr_intersectList(i,ins,iph,ith) - ipix)/0.999

    tr_bhCdMaps(tr_odICDTO,ipix,ii,jj,kk) &
    & = tr_bhCdMaps(tr_odICDTO,ipix,ii,jj,kk) + columnto*weight
    !print *, "NC3: ", i, tr_bhCdMaps(tr_odICDTO,ipix,ii,jj,kk), columnto, weight
#if CHEMISTRYNETWORK == 4 || CHEMISTRYNETWORK == 5
    tr_bhCdMaps(tr_odICDH2,ipix,ii,jj,kk) &
    & = tr_bhCdMaps(tr_odICDH2,ipix,ii,jj,kk) + columnH2*weight
#endif                                                                       
#if CHEMISTRYNETWORK == 5                                                    
    tr_bhCdMaps(tr_odICDCO,ipix,ii,jj,kk) &
    & = tr_bhCdMaps(tr_odICDCO,ipix,ii,jj,kk) + columnCO*weight
#endif
  enddo
  !print *, "NC4: ", tr_bhCdMaps(tr_odICDTO,ipix,ii,jj,kk), columnto, weight

  return
end subroutine tr_odNodeContrib
