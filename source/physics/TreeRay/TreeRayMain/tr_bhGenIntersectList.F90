!!****if* source/physics/TreeRay/TreeRayMain/tr_bhGenIntersectList
!!
!! NAME
!!
!!  tr_bhGenIntersectList
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

subroutine tr_bhGenIntersectList()
  use TreeRay_data, ONLY : tr_nPix, tr_nSide, tr_ilNTheta, tr_ilNPhi, tr_ilNNS, &
    tr_ilNR, tr_ilNI, tr_ilFinePix, tr_intersectList, tr_meshMe, tr_comm, &
    tr_numProcs, tr_ilNSSampFac, tr_ilNSSampFacI
  use healpix, ONLY : pix2ang_nest
  use Logfile_interface, ONLY : Logfile_stamp
  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
  integer :: istat, il, iclosest, max_il, i1d, chunk, ierr
  integer :: ith, iph, ins, ipix, ifpix, ispix, ir
  integer :: node_count, ray_count(0:tr_nPix-1)
  real :: theta, phi, ns, nshalf, thfpix, phfpix
  real :: x, y, z, r, xnode, ynode, znode, dotprod, max_dotprod
  character(len=MAX_STRING_LENGTH) :: strBuff, tr_ilFName = "intersect-list"
  real :: loc_iList(tr_ilNI, tr_ilNNS, 0:tr_ilNPhi, 0:tr_ilNTheta)

  tr_ilNSSampFacI = 1.0/tr_ilNSSampFac

  allocate(tr_intersectList(tr_ilNI, tr_ilNNS, 0:tr_ilNPhi, 0:tr_ilNTheta), stat=istat)
  tr_intersectList = -1.0
  max_il = 0
  loc_iList = -1.0

  chunk = 1 + ((tr_ilNTheta+1)*(tr_ilNPhi+1) / tr_numProcs)
  if (tr_meshMe == MASTER_PE) then
    print *, "tr_ilNTheta, tr_ilNPhi, tr_numProcs = ", tr_ilNTheta, tr_ilNPhi, tr_numProcs
    print *, "chunk = ", chunk
  endif

  do ith = 0,tr_ilNTheta
    theta = ith*PI/tr_ilNTheta
    do iph = 0,tr_ilNPhi
      phi = iph*2*PI/tr_ilNPhi

      i1d = ith*(tr_ilNPhi+1) + iph
      !print *, "i1d = ", i1d, ith, iph
      ! check if this point should be calculated on this processor
      if ((i1d >= tr_meshMe*chunk) .and. (i1d < (tr_meshMe+1)*chunk)) then
        !print *, "accepted i1d at: ", i1d, tr_meshMe

        xnode = sin(theta)*cos(phi)
        ynode = sin(theta)*sin(phi)
        znode = cos(theta)
        do ins = 1,tr_ilNNS
          ! relative sizes of nodes will be sampled up to 1.5
          nshalf = 0.5*tr_ilNSSampFac*real(ins)/tr_ilNNS
          node_count = 0.0
          ray_count = 0.0
          ! check if the central point is out of the node
          !if (  (abs(xnode) < ns) .and. (abs(ynode) < ns) &
          !.and. (abs(znode) < ns)) cycle
       
          do ipix = 0,tr_nPix-1
            do ispix = 0,tr_ilFinePix**2-1 
              ! calculate number of the pixel in fine pixelation
              ifpix = tr_ilFinePix**2*ipix + ispix
              !print *, "calling pix2ang: ", tr_nSide, tr_ilFinePix, ifpix, thfpix, phfpix
              call pix2ang_nest(tr_nSide*tr_ilFinePix, ifpix, thfpix, phfpix)
              do ir = 0,tr_ilNR
                r = ir*3.0/tr_ilNR
                x = r*sin(thfpix)*cos(phfpix)
                y = r*sin(thfpix)*sin(phfpix)
                z = r*cos(thfpix)
                if (  (abs(x-xnode) < nshalf) .and. (abs(y-ynode) < nshalf) &
                .and. (abs(z-znode) < nshalf)) then
                  node_count = node_count + 1
                  ray_count(ipix) = ray_count(ipix) + 1
                endif
       
              enddo
       
            enddo ! loop over sub-rays
       
          enddo ! loop over all rays
          !if (i1d == 362) print *, "after all rays: ", ith, iph, ins, theta, phi, real(ray_count)/real(node_count)
       
          ! store list of intersections between the node and rays into the array
          il = 1
          do ipix = 0,tr_nPix-1
            if (ray_count(ipix) > 0) then
              if (il > max_il) max_il = il
              if (il > tr_ilNI) call Driver_abortFlash("TreeRay_bhGenIntersectList: Too many intersection. Increase tr_ilNI.")
              loc_iList(il, ins, iph, ith) = (ray_count(ipix)*0.999)/node_count + ipix
              il = il + 1
            endif
          enddo
          !if (i1d == 362) print *, "loc_ilist = ", loc_iList(:, ins, iph, ith)
       
          ! in case no intersection is found, just add the closest ray
          if (node_count == 0) then
            max_dotprod = 0.0
            do ipix = 0,tr_nPix-1
              call pix2ang_nest(tr_nSide, ipix, thfpix, phfpix)
              x = sin(thfpix)*cos(phfpix)
              y = sin(thfpix)*sin(phfpix)
              z = cos(thfpix)
              dotprod = x*xnode + y*ynode + z*znode
              if (dotprod > max_dotprod) then
                max_dotprod = dotprod
                iclosest = ipix
              endif
            enddo
            loc_iList(1,ins,iph,ith) = real(iclosest) + 0.999
          endif
          !if (i1d == 362) print *, "loc_ilist2 = ", ins, iph, ith, loc_iList(:, ins, iph, ith)
          !print *, "loc_ilist3 = ", ins, iph, ith, loc_iList(:, ins, iph, ith)

          !if (loc_iList(1,ins,iph,ith) < 0) then
          !  print *, 'No intersection: ', ins, ith, iph, node_count, loc_iList(1,ins,iph,ith)
          !  stop
          !endif
       
        enddo ! loop over nodes of different sizes
      endif ! distribution of work among processors
    enddo
  enddo

  !print *, "before allreduce: ", tr_ilNI*tr_ilNNS*(tr_ilNPhi+1)*(tr_ilNTheta+1), size(tr_intersectList), size(loc_iList)
  !print *, "Bef iList1675: ", loc_iList(1, 16, 5, 7), tr_intersectList(1, 16, 5, 7)
  call MPI_Barrier(tr_comm, ierr)
  call MPI_AllReduce(loc_iList,tr_intersectList, &
  tr_ilNI*tr_ilNNS*(tr_ilNPhi+1)*(tr_ilNTheta+1), FLASH_REAL, MPI_MAX, tr_comm, ierr)  
  !print *, "Aft iList1675: ", loc_iList(1, 16, 5, 7), tr_intersectList(1, 16, 5, 7)
  !do ins = 1,tr_ilNNS
  !  do ith = 0,tr_ilNTheta
  !    do iph = 0,tr_ilNPhi
  !      if (tr_intersectList(1,ins,iph,ith) < 0) then
  !        print *, 'No intersection2: ', ins, ith, iph, node_count, tr_intersectList(1,ins,iph,ith)
  !        stop
  !      endif
  !    enddo
  !  enddo
  !enddo

  ! write tr_intersectList into file
  if (tr_meshMe == MASTER_PE) then
 
    open(unit=53, file=tr_ilFName, status='replace')
    write(53,*) "# ", tr_ilNNS, tr_ilNTheta, tr_ilNPhi
    do ins = 1,tr_ilNNS
      ns = tr_ilNSSampFac*real(ins)/tr_ilNNS
      do ith = 0,tr_ilNTheta
        theta = ith*PI/tr_ilNTheta
        do iph = 0,tr_ilNPhi
          phi = iph*2*PI/tr_ilNPhi
          xnode = sin(theta)*cos(phi)
          ynode = sin(theta)*sin(phi)
          znode = cos(theta)
          write(53,'(6(e12.5, 2x))',advance='no') ns, theta, phi, xnode, ynode, znode
          do il = 1, tr_ilNI
            write(53,'(e15.8, 2x)',advance='no') tr_intersectList(il, ins, iph, ith)
          enddo
          write(53,*)
        enddo
        write(53,*)
      enddo
      write(53,*)
    enddo
    close(unit=53)
    write (strBuff, '("Intersection list written into file:", a20)') tr_ilFName
    call Logfile_stamp( strBuff, "[TreeRay]")
    write (strBuff, '("Maximum IL = ", i3)') max_il
    call Logfile_stamp( strBuff, "[TreeRay]")
  endif ! MASTER_PE

end subroutine tr_bhGenIntersectList
