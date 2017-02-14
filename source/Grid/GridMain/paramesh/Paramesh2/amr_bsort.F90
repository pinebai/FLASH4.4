
      subroutine amr_bi_sort (list,gid,npp)


! $RCSfile: amr_bsort.F90,v $
! $Revision: 1.1 $
! $Date: 2004/05/20 01:03:44 $


use physicaldata
      use tree
      implicit none
      include 'mpif.h'


      integer mpp
      parameter (mpp=maxblocks_tr)

      integer list(mpp),dlist(2*mpp)
      integer npp,nppt,nppt2
      integer gid(mpp)
      integer gidt(mpp),gidt2(mpp)
      integer ix(mpp),iproc(mpp)
      integer dix(2*mpp),diproc(2*mpp)
      integer npp_max,npp_max2,near_p_2
      integer i,j,k
!!!
      integer ierr
      integer reqs(mpp),reqr(mpp)
      integer stats(MPI_STATUS_SIZE,mpp)
      integer statr(MPI_STATUS_SIZE,mpp)
      integer nsend,nrecv
      integer kk(mpp)
!!!
      integer mype,nprocs

      save nppt,nppt2,gidt,gidt2

      call MPI_COMM_RANK (MPI_COMM_WORLD,mype,ierr)
      call MPI_COMM_SIZE (MPI_COMM_WORLD,nprocs,ierr)

! find max local list length among processors

      nppt = npp
      npp_max = npp
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      call MPI_ALLREDUCE (npp_max,npp_max2,1,MPI_INTEGER, & 
     &                    MPI_MAX,MPI_COMM_WORLD,ierr)
      npp_max = npp_max2

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

! compute nearest power of 2

      near_p_2 = 1
      do while (near_p_2.lt.npp_max)

         near_p_2 = 2*near_p_2

      end do

! load dummy sorting list

      do i = 1,npp
         dlist(i) = list(i)
      end do
      do i = npp+1,mpp
         dlist(i) = 1000000000
      end do

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

! sort dummy list

      call amr_bsort(dlist,dix,diproc,near_p_2)

! compute global id for each item in the sorted list

! valid items in sorted list send global ids back to their original
! locations 

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      do i = 1,near_p_2
        gidt(i) = i + near_p_2*mype
      end do

! SEND gidt back to the old location of the sorted item and put it in
! gidt2

      nrecv = 0
      do i = 1,near_p_2
         if (diproc(i).eq.mype) then
            gidt2(dix(i)) = gidt(i)
          else
            nrecv = nrecv + 1
         end if
      end do

      do i = 1,nrecv
         call MPI_IRECV(kk(i), & 
     &        1, & 
     &        MPI_INTEGER, & 
     &        MPI_ANY_SOURCE, & 
     &        MPI_ANY_TAG, & 
     &        MPI_COMM_WORLD, & 
     &        reqr(i), & 
     &        ierr)
      end do

      nsend = 0
      do i = 1,near_p_2
         if (diproc(i).ne.mype) then
            nsend = nsend + 1
            call MPI_SSEND(gidt(i), & 
     &           1, & 
     &           MPI_INTEGER, & 
     &           diproc(i),     &  ! PE TO SEND TO
     &           dix(i),        &  ! THIS IS THE TAG
     &           MPI_COMM_WORLD, & 
     &           ierr)
         end if
      end do

      if (nrecv.gt.0) then
         call MPI_WAITALL (nrecv, reqr, statr, ierr)
         do i = 1,nrecv
            gidt2(statr(MPI_TAG,i)) = kk(i)
         end do
      end if

      do i = 1,near_p_2
         list(i) = dlist(i)
         gid(i) = gidt2(i)
      end do

      call MPI_BARRIER (MPI_COMM_WORLD,ierr)

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine amr_bsort (list,ix,iproc,npp)

! Written K. Olson, 6/97

use physicaldata
      use tree
      implicit none
      include 'mpif.h'


      integer npp
      integer list(npp)
      integer n,nprocs
      integer n_fold,fold_len,next_fold(maxblocks_tr)
      integer ix(npp),iproc(npp)
      integer i,j,k,ierr
      integer mype

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      if (npp.gt.100000) print *,' ERROR in bsort: npp > 100000 '

      n_fold = 0
      fold_len = 1
      call MPI_COMM_SIZE (MPI_COMM_WORLD,nprocs,ierr)
      call MPI_COMM_RANK (MPI_COMM_WORLD,mype,ierr)

      do i = 1,npp
         ix(i) = i
         iproc(i) = mype
      end do

      call amr_int_heapsort (list,ix,iproc,npp)

      do while (fold_len.le.nprocs/2)

         call amr_fold_and_sort (list,ix,iproc,fold_len,npp)
         
         n_fold = n_fold + 1
         next_fold(n_fold) = fold_len

         do j = 1,n_fold-1

            call amr_fold_and_sort (list,ix,iproc,next_fold(j),npp)

            n_fold = n_fold + 1
            next_fold(n_fold) = next_fold(j)

            call MPI_BARRIER(MPI_COMM_WORLD,ierr)

         end do

         fold_len = 2*fold_len

         call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      end do

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine amr_fold_and_sort (list,ix,iproc,fold_len,npp)

      implicit none

      include 'mpif.h'

      integer nx
      parameter (nx = 4*16384)

      integer npp
      integer list(npp),ix(npp),iproc(npp)
      integer list_temp(npp),ix_temp(npp),iproc_temp(npp)
      integer fold_len,list_len,mype,list_no,left_most,pivot
      integer delta,fetch_from
      integer i,j,k
      integer list_com(nx),ix_com(nx),iproc_com(nx)
      integer list_com2(nx),ix_com2(nx),iproc_com2(nx)
      integer ierr
      integer stat(MPI_STATUS_SIZE)

      common /comm_fold/ list_com,ix_com,iproc_com,list_com2,ix_com2, & 
     &     iproc_com2

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      if (npp.gt.nx) print *,' ERROR IN fold_and_sort: npp > nx '

      call MPI_COMM_RANK(MPI_COMM_WORLD,mype,ierr)

      list_len = fold_len*2
      list_no = int(mype/list_len)
      left_most = list_len*list_no
      pivot = left_most + fold_len - 1
      
      if (mype.le.pivot) then
         delta = pivot - mype + 1
         fetch_from = pivot + delta
      else
         delta = mype - pivot - 1
         fetch_from = pivot - delta
      end if

      do i = 1,npp
         list_com(i) = list(i)
         ix_com(i) = ix(i)
         iproc_com(i) = iproc(i)
      end do

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

! For MPI mearly send data from local proc to fetch_from
! then issue a recieve to recieve incoming data which has been
! sent to local prom (fetch_from = send_to !!!)

      call MPI_SENDRECV  & 
     & (list_com(1),npp,MPI_INTEGER,fetch_from,1, & 
     &  list_com2(1),npp,MPI_INTEGER,fetch_from,1, & 
     &  MPI_COMM_WORLD,stat,ierr)

      call MPI_SENDRECV  & 
     & (ix_com(1),npp,MPI_INTEGER,fetch_from,2, & 
     &  ix_com2(1),npp,MPI_INTEGER,fetch_from,2, & 
     &  MPI_COMM_WORLD,stat,ierr)

      call MPI_SENDRECV  & 
     & (iproc_com(1),npp,MPI_INTEGER,fetch_from,3, & 
     &  iproc_com2(1),npp,MPI_INTEGER,fetch_from,3, & 
     &  MPI_COMM_WORLD,stat,ierr)

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      do i = 1,npp
         k = npp - i + 1
         list_temp(i) = list_com2(k)
         ix_temp(i) = ix_com2(k)
         iproc_temp(i) = iproc_com2(k)
      end do

      if (mype.le.pivot) then

         do i = 1,npp

            if (list_temp(i).lt.list(i)) then
               list(i) = list_temp(i)
               ix(i) = ix_temp(i)
               iproc(i) = iproc_temp(i)
            end if

         end do

      else

         do i = 1,npp

            if (list_temp(i).gt.list(i)) then
               list(i) = list_temp(i)
               ix(i) = ix_temp(i)
               iproc(i) = iproc_temp(i)
            end if

         end do

      end if

      call amr_int_heapsort (list,ix,iproc,npp)

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      return
      end

      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      subroutine amr_int_heapsort(ra,ix,iproc,n) 
!From "Numerical Recipes", index added
      
      implicit none
      
      integer n
      integer ra(n)
      integer ix(n),iproc(n)
      
! IMPLICIT DECLARATIONS
      
      integer rra
      integer i,j,k,L,ir,iix,iiproc
      
      if (n .lt. 2) return
      L = n/2+1
      ir = n
 10   continue
      if (L .gt. 1) then
         L = L-1
         rra = ra(L)
         iix = ix(L)
         iiproc = iproc(L)
      else
         rra = ra(ir)
         iix = ix(ir)
         iiproc = iproc(ir)
         ra(ir) = ra(1)
         ix(ir) = ix(1)
         iproc(ir) = iproc(1)
         ir = ir-1
         if (ir .eq. 1) then
            ra(1) = rra
            ix(1) = iix
            iproc(1) = iiproc
            return
         end if
      end if
      i = L
      j = L+L
 20   if (j .le. ir) then
         if (j .lt. ir) then
            if (ra(j) .lt. ra(j+1))   j = j+1
         end if
         if (rra .lt. ra(j)) then
            ra(i) = ra(j)
            ix(i) = ix(j)
            iproc(i) = iproc(j)
            i = j
            j = j+j
         else
            j = ir+1
         end if
         go to 20
      end if
      ra(i) = rra
      ix(i) = iix
      iproc(i) = iiproc
      go to 10
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine amr_int_heapsort2(ra,ix,n) 
!From "Numerical Recipes", index added
      
      implicit none
      
      integer n
      integer ra(n)
      integer ix(n)
      
! IMPLICIT DECLARATIONS
      
      integer rra
      integer i,j,k,L,ir,iix
      
      if (n .lt. 2) return
      L = n/2+1
      ir = n
 10   continue
      if (L .gt. 1) then
         L = L-1
         rra = ra(L)
         iix = ix(L)
      else
         rra = ra(ir)
         iix = ix(ir)
         ra(ir) = ra(1)
         ix(ir) = ix(1)
         ir = ir-1
         if (ir .eq. 1) then
            ra(1) = rra
            ix(1) = iix
            return
         end if
      end if
      i = L
      j = L+L
 20   if (j .le. ir) then
         if (j .lt. ir) then
            if (ra(j) .lt. ra(j+1))   j = j+1
         end if
         if (rra .lt. ra(j)) then
            ra(i) = ra(j)
            ix(i) = ix(j)
            i = j
            j = j+j
         else
            j = ir+1
         end if
         go to 20
      end if
      ra(i) = rra
      ix(i) = iix
      go to 10
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine amr_iq_sort(n,arr,indx,indx2)

! From Numerical Recipes, page 324

       use Driver_interface, ONLY : Driver_abortFlash
      implicit none

      integer n,indx(n),indx2(n),M,NSTACK
      integer arr(n)
      parameter (M=7,NSTACK=500)
      integer i,ir,temp,j,jstack,k,l,istack(NSTACK)
      integer a,b,c

      jstack = 0
      l = 1
      ir = n
 1    if (ir-l.lt.M) then
         do j = l+1,ir

            a = arr(j)
            b = indx(j)
            c = indx2(j)
            
            do i = j-1,1,-1
               if (arr(i).le.a) go to 2
               
               arr(i+1) = arr(i)
               indx(i+1) = indx(i)
               indx2(i+1) = indx2(i)
               
            end do
            i = 0
            
 2          arr(i+1) = a
            indx(i+1) = b
            indx2(i+1) = c
            
         end do

         if (jstack.eq.0) return
         ir = istack(jstack)
         l = istack(jstack-1)
         jstack = jstack-2
         
      else

         k = (l+ir)/2
         
         temp = arr(k)
         arr(k) = arr(l+1)
         arr(l+1) = temp
         temp = indx(k)
         indx(k) = indx(l+1)
         indx(l+1) = temp
         temp = indx2(k)
         indx2(k) = indx2(l+1)
         indx2(l+1) = temp
            
         if (arr(l+1).gt.arr(ir)) then
            
            temp = arr(l+1)
            arr(l+1) = arr(ir)
            arr(ir) = temp
            temp = indx(l+1)
            indx(l+1) = indx(ir)
            indx(ir) = temp
            temp = indx2(l+1)
            indx2(l+1) = indx2(ir)
            indx2(ir) = temp

         end if
         if (arr(l).gt.arr(ir)) then
            
            temp = arr(l)
            arr(l) = arr(ir)
            arr(ir) = temp
            temp = indx(l)
            indx(l) = indx(ir)
            indx(ir) = temp
            temp = indx2(l)
            indx2(l) = indx2(ir)
            indx2(ir) = temp
            
         end if
         if (arr(l+1).gt.arr(l)) then
            
            temp = arr(l+1)
            arr(l+1) = arr(l)
            arr(l) = temp
            temp = indx(l+1)
            indx(l+1) = indx(l)
            indx(l) = temp
            temp = indx2(l+1)
            indx2(l+1) = indx2(l)
            indx2(l) = temp
            
         end if
         i = l+1
         j = ir
         
         a = arr(l)
         b = indx(l)
         c = indx2(l)
         
 3       continue
         i = i + 1
         if (arr(i).lt.a) go to 3
 4       continue
         j = j-1
         if (arr(j).gt.a) go to 4
         if (j.lt.i) go to 5
         
         temp = arr(i)
         arr(i) = arr(j)
         arr(j) = temp
         temp = indx(i)
         indx(i) = indx(j)
         indx(j) = temp
         temp = indx2(i)
         indx2(i) = indx2(j)
         indx2(j) = temp
         
         go to 3
         
 5       arr(l) = arr(j)
         arr(j) = a
         indx(l) = indx(j)
         indx(j) = b
         indx2(l) = indx2(j)
         indx2(j) = c
         
         jstack = jstack+2

         if (jstack.gt.NSTACK) then
!! Removed pause as it is a) deprecated in F90 and b) stupid in a parallel environment
!!		pause 'NSTACK too small in iq_sort'
		call Driver_abortFlash('NSTACK too small in amr_iq_sort')
         end if
         if (ir-i+1.ge.j-1) then
            istack(jstack) = ir
            istack(jstack-1) = i
            ir = j-1
         else
            istack(jstack) = j-1
            istack(jstack-1) = l
            l = i
         end if
         
      end if
      
      go to 1

      end
