!!****if* source/physics/sourceTerms/Flame/FlameSpeed/laminar/CONe/fl_fsLaminarInit
!!
!! NAME
!!
!!  fl_fsLaminarInit
!!
!! SYNOPSIS
!!
!!  call fl_fsLaminarInit()
!!
!! DESCRIPTION
!!
!!
!!  Aaron Jackson, Alan Calder 2008
!!
!!   Allocates and fills table for laminar flame speed look-up through X12, X22,
!!   LOG(den) from Chamulak, Brown, & Timmes.
!!  
!!  table convention i --> log(rho), j --> X22, k --> X12
!!
!! ARGUMENTS
!!
!!   No arguments
!!
!!
!!
!!***


subroutine fl_fsLaminarInit()

  use fl_fsConeData, ONLY : nc12it, nne22it, nldent,      &
                            ilt, iut, jlt, jut, klt, kut, &
                            lda, lna, lca, &
                            c12ib, ne22ib, ldenb, stab, dstab
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

  character (len=50),save :: tablename

  integer :: i, j, k, istat

  real :: rtemp1, rtemp2, rtemp3

  call RuntimeParameters_get("fl_fsCONeTableName",tablename)

  open (unit=21,file=tablename,status='unknown',iostat=istat)
  if (istat /= 0) call Driver_abortFlash("Unable to open Ne22 Flame Table in fl_fsConeInitTable")

  read(21,*) nc12it
  read(21,*) nne22it
  read(21,*) nldent
  read(21,*) 
!  write(6,*)'number of c12is' , nc12it
!  write(6,*)'number of ne22is', nne22it 
!  write(6,*)'number of ldens' , nldent

  klt = 1
  kut = nc12it
  lca = (klt+kut)/2

  jlt = 1
  jut = nne22it
  lna = (jlt+jut)/2

  ilt = 1
  iut = nldent
  lda = (ilt+iut)/2


  allocate(c12ib(klt:kut),stat=istat)
  if (istat /= 0) call Driver_abortFlash("Cannot allocate c12ib in fl_fsConeInitTable")
  allocate(ne22ib(jlt:jut),stat=istat)
  if (istat /= 0) call Driver_abortFlash("Cannot allocate ne22ib in fl_fsConeInitTable")
  allocate(ldenb(ilt:iut),stat=istat)
  if (istat /= 0) call Driver_abortFlash("Cannot allocate ldenb in fl_fsConeInitTable")
  allocate(stab(ilt:iut,jlt:jut,klt:kut),stat=istat)
  if (istat /= 0) call Driver_abortFlash("Cannot allocate stab in fl_fsConeInitTable")
  allocate(dstab(ilt:iut,jlt:jut,klt:kut),stat=istat)
  if (istat /= 0) call Driver_abortFlash("Cannot allocate dstab in fl_fsConeInitTable")

!! read the table 
!! assume X12 spans to 1
  do k = klt,kut-(jut-jlt+1)
     do j = jlt,jut
        do i = ilt,iut
           read(21,*) rtemp1, rtemp2, rtemp3, stab(i,j,k), dstab(i,j,k)
        enddo
     enddo
  enddo
  do k = kut-(jut-jlt+1)+1,kut
     do j = jlt, jut
        if (j == kut-k+1) then
           do i = ilt,iut
              read(21,*) rtemp1, rtemp2, rtemp3, stab(i,j,k), dstab(i,j,k)
           enddo
        else
            stab(:,j,k) = 0.0e0
           dstab(:,j,k) = 0.0e0
        endif
     enddo
  enddo

!! get the lden, ne22i, and c12i arrays. 
    
  rewind(21)
  read(21,*)
  read(21,*)
  read(21,*)
  read(21,*)
  do k = klt,kut-(jut-jlt+1)
     do j = jlt,jut
        do i = ilt,iut
           read(21,*) rtemp1, rtemp2, ldenb(i)
        enddo
        ne22ib(j) = rtemp2
     enddo
     c12ib(k) = rtemp1
  enddo
  do k = kut-(jut-jlt+1)+1,kut
     do j = jlt,jut
        do i = ilt,iut
           if (j == kut-k+1) then
!! assume densities have already been filled above.
              read(21,*) rtemp1, rtemp2  !, ldenb(i)
           elseif (k == kut-(jut-jlt+1)+1) then
               stab(i,j,k) =  stab(i,j,kut-j+1)  !! from condition above
              dstab(i,j,k) = dstab(i,j,kut-j+1)
           else
               stab(i,j,k) =  stab(i,j,k-1) !! allows for extrapolation to
              dstab(i,j,k) = dstab(i,j,k-1)  !! higher X12
           endif
        enddo
!! ne22 is ordered opposite the rest of the table due to o16=0
!       if (j == kut-k+1) ne22ib(j) = rtemp2
     enddo
     c12ib(k) = rtemp1
  enddo

  close(21)

  return
end subroutine
