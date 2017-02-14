!
!   subroutine b_<Type>_sendrcv(maxbatchlen,ns,sproc,stag,sval,
!                                            nr,rproc,rtag,rval)
!
!   These subroutines do a batched send of ns messages to procs
!    sproc, with tags stag and values svals, followed by batched
!    recieves of nr messages from rproc, into arrays rtag and rval.
!
!   <TYPE> == real
!
!   maxbatchlen: 
!            INTEGER (in) maximum number of messages to batch per sent msg
!   msglen:  INTEGER (in) number of <TYPE>s per message (not incl. tag)
!
!   ns:      INTEGER  (in)  number of messages to send
!   sproc:   INTEGER  (in)  processor to send message (i=1..ns) to 
!   stag:    <TYPE>   (in)  tag to send with value to processor sproc(i)
!   sval:    <TYPE>   (in)  value to send to processor sproc(i)
!
!   nr:      INTEGER  (in)  number of messages to recieve
!   rproc:   INTEGER  (in)  processor to recieve message (i=1..nr) from
!   rtag:    <TYPE>   (in)  tag to recieve with value from processor rproc(i)
!   rval:    <TYPE>   (out) value to recieve from processor rproc(i)
!

!#define DEBUG_BATCH
!#define VERBOSE_BATCH
!#define ERRROR_CHECK


       subroutine b_dbl_sendrcv(maxbatchlen,msglen, & 
     &                          ns,sproc,stag,sval, & 
     &                          nr,rproc,rtag,rval)
       

       IMPLICIT NONE

       integer :: ns, nr
       integer :: maxbatchlen, msglen
       integer :: sproc(ns), rproc(nr)
       integer :: stag(ns),  rtag(nr)
       real :: sval(ns*msglen),  rval(nr*msglen)

       integer :: nsproc, nrproc, nprocs, mype
       integer :: nmsgssent, nmsgsrecv, pair

       integer :: i,j,k,p,n,m
       integer :: istat,ierr

       integer :: tag,val,npairs,nbuffs

       real, allocatable, dimension(:,:)  ::  msgs, rmsgs
       integer, allocatable :: rcount(:)
       integer, allocatable :: oproc(:),omsgs(:),omsgsin(:),omsgsleft(:)
       integer, allocatable :: iproc(:),imsgs(:),imsgsleft(:),procno(:)
       integer, allocatable :: ibatches(:)
       integer, allocatable :: nmsg(:), neighnum(:), sneighnum(:)
       integer, allocatable :: rreqs(:)
       logical, allocatable :: rdone(:)

        include 'mpif.h'

       integer stat(MPI_STATUS_SIZE)



#ifdef VERBOSE_BATCH
       print *,mype,': PAIR SENDS/RECVS', ns,nr
#endif

!CCC
!    RECIEVE POSTING PHASE
!CCC

       call MPI_COMM_SIZE (MPI_COMM_WORLD,nprocs,ierr)
       call MPI_COMM_RANK (MPI_COMM_WORLD,mype,ierr)



!
!   find out who we have to get data from, and how much.
!      iproc(:)   -- processor #s we have to get messages from.
!      imsgs(:)   -- how many [tag val] pairs we have to get for each proc.
!      ibatches(:)-- how many batches the pairs get aggregated into
!      imsgsleft(:)- how many pairs have yet to be batched
!    

       allocate( nmsg(nprocs), neighnum(nprocs), stat = istat )
       if (istat .ne. 0)call Driver_abortFlash('batchsend: alloc nmsg')
       allocate( sneighnum(nprocs), stat = istat )
       if (istat .ne. 0)call Driver_abortFlash('batchsend: alloc sneighnum')
       allocate( rdone(nr), stat = istat)
       if (istat .ne. 0)call Driver_abortFlash('batchsend: alloc rdone')

       nmsg(:) = 0
       neighnum(:) = 0
       nrproc = 0
       do i=1,nr
          rdone(i) = .false.
          if (nmsg(rproc(i)+1) .eq. 0) then
             nrproc = nrproc + 1
             neighnum(rproc(i)+1) = nrproc
          endif
          nmsg(rproc(i)+1) = nmsg(rproc(i)+1) + 1
       enddo

       allocate(iproc(nrproc), imsgs(nrproc), stat=istat)
       if (istat .ne. 0) call Driver_abortFlash('batchsend:lloc iproc,imsgs')
       allocate(ibatches(nrproc), stat=istat)
       if (istat .ne. 0) call Driver_abortFlash('batchsend:lloc ibatches')
       allocate(imsgsleft(nrproc),stat=istat)
       if (istat .ne. 0) call Driver_abortFlash('batchsend: alloc imsgsleft')

       j = 1
       nmsg(:) = 0
       do i=1,nr
          if (nmsg(rproc(i)+1) .eq. 0) then
             iproc(j) = rproc(i)
             j= j + 1
          endif
          nmsg(rproc(i)+1) = nmsg(rproc(i)+1) + 1
       enddo

       do i=1,nrproc
         p = iproc(i)
         imsgs(i) = nmsg(p+1)
         imsgsleft(i) = nmsg(p+1)
       enddo 

#ifdef DEBUG_BATCH 
       print *,mype,' has ', nrproc, ' neighbors to recv from' 
       print *,mype,' procs: ',(iproc(i),i=1,nrproc)
       print *,mype,' nmsgs: ',(imsgs(i),i=1,nrproc)
#endif

       nmsgsrecv = 0
       do i=1,nrproc
          if (imsgs(i) .le. maxbatchlen) then
              ibatches(i) = 1
          else
              call intdivrndup(imsgs(i),maxbatchlen,m)
              ibatches(i) = m
          endif
          nmsgsrecv = nmsgsrecv + ibatches(i)
       enddo

!
!   allocate one message per incoming batch
!
       allocate( rmsgs((1+msglen)*maxbatchlen+1,nmsgsrecv),stat=istat)
       if (istat .ne. 0) call Driver_abortFlash('b_in: alloc rmsgs')
       allocate( rcount(nmsgsrecv), stat=istat )
       if (istat .ne. 0) call Driver_abortFlash('b_in: alloc rcount')
       allocate( rreqs(nmsgsrecv), stat=istat )
       if (istat .ne. 0) call Driver_abortFlash('b_in: alloc rreqs')
       allocate( procno(nmsgsrecv), stat=istat )
       if (istat .ne. 0) call Driver_abortFlash('b_in: alloc procno')

!
!   Post the appropriate number of recieves
!

       n = 1
       do i=1,nrproc
           do j=1,ibatches(i)
              procno(n) = iproc(i)
              rcount(n) = ibatches(i)*2+1
#ifdef DEBUG_BATCH
              print *, mype, 'posting a recieve from ',iproc(i)
#endif
              call MPI_IRECV(rmsgs(1,n), & 
     &                      (1+msglen)*maxbatchlen+1, & 
     &                      MPI_DOUBLE_PRECISION, & 
     &                      iproc(i), & 
     &                      iproc(i), & 
     &                      MPI_COMM_WORLD, & 
     &                      rreqs(n), & 
     &                      istat)
               n = n + 1
           enddo
       enddo

!CCC
!    SEND PHASE
!CCC
!
!   find out who we have to send data to, and how much.
!      oproc(:)   -- processor #s we have to send messages to.
!      omsgs(:)   -- how many [tag val] pairs we have to send
!      omsgsin(:) -- how many pairs have already been batched into the
!                    current (pending) message
!      omsgsleft(:)- how many pairs have yet to be batched
!    

       nmsg(:) = 0
       sneighnum(:) = 0
       nsproc = 0
       do i=1,ns
          if (nmsg(sproc(i)+1) .eq. 0) then
             nsproc = nsproc + 1
             sneighnum(sproc(i)+1) = nsproc
          endif
          nmsg(sproc(i)+1) = nmsg(sproc(i)+1) + 1
       enddo

       allocate(oproc(nsproc), omsgs(nsproc), stat=istat)
       if (istat .ne. 0) call Driver_abortFlash('b_in: alloc oproc,omsgs')
       allocate(omsgsin(nsproc),omsgsleft(nsproc),stat=istat)
       if (istat .ne. 0) call Driver_abortFlash('b_in: alloc omsgsin')

       j = 1
       nmsg(:) = 0
       do i=1,ns
          if (nmsg(sproc(i)+1) .eq. 0) then
             oproc(j) = sproc(i)
             j= j + 1
          endif
          nmsg(sproc(i)+1) = nmsg(sproc(i)+1) + 1
       enddo

       do i=1,nsproc
         p = oproc(i)
         omsgs(i) = nmsg(p+1)
         omsgsleft(i) = nmsg(p+1)
         omsgsin(i) = 0
       enddo 

#ifdef DEBUG_BATCH
       print *,mype,' has ', nsproc, ' neighbors to send to' 
       print *,mype,' procs: ',(oproc(i),i=1,nsproc)
       print *,mype,' nmsgs: ',(omsgs(i),i=1,nsproc)
#endif

       nbuffs = 0
       do i=1,nsproc
          if (omsgs(i).le. maxbatchlen) then
              nbuffs = nbuffs + 1
          else
              call intdivrndup(omsgs(i),maxbatchlen,m)
              nbuffs = nbuffs + m
          endif
       enddo


!
!   allocate one message per `neighboring' processor
!
       allocate( msgs((1+msglen)*maxbatchlen+1,nsproc),stat=istat)
       if (istat .ne. 0) call Driver_abortFlash('b_in: alloc msgs')

!
!   send the messages.
!
!   go through the list of messages to send, and construct the message:
!      msg(proc #i) = [num of msgs] [tag val] [tag val] ....
!
!   if we hit the pack length maximum (MAXBATCHLEN), send off the message
!    immediately and start a new one; otherwise, send off a message only
!    when all the data bound for processor p is aggregated.
!
       nmsgssent = 1
       do i=1,ns
          p = sneighnum(sproc(i)+1)

          msgs(omsgsin(p)*(msglen+1) + 2,p) = REAL(stag(i))
          do j=1,msglen 
            msgs(omsgsin(p)*(msglen+1)+2+j,p) = sval(msglen*(i-1)+j)
          enddo

          omsgsin(p) = omsgsin(p) + 1

!   send a message if we've hit the max, 

          if (omsgsin(p) .eq. maxbatchlen) then
            
             msgs(1,p) = maxbatchlen
             call MPI_SSEND(msgs(1,p), & 
     &                      maxbatchlen*(msglen+1)+1, & 
     &                      MPI_DOUBLE_PRECISION, & 
     &                      sproc(i), & 
     &                      mype, & 
     &                      MPI_COMM_WORLD, & 
     &                      ierr)
#ifdef DEBUG_BATCH
             print *,mype,' sent ',maxbatchlen,' to ',sproc(i)
             do j=1,omsgsin(p)
             print *,mype,'      [',msgs((j-1)*(msglen+3)+1,p),':' & 
     &              ,(msgs((j-1)*(msglen+1)+2+k,p),k=1,msglen),']'
             enddo
#endif
             if (ierr .ne. 0) call Driver_abortFlash('b_int: bsend 1')
             omsgsleft(p) = omsgsleft(p) - omsgsin(p)
             omsgsin(p) = 0

             nmsgssent = nmsgssent + 1

!    or if we're done for this neighbor.

          else if (omsgsin(p) .eq. omsgsleft(p)) then
             msgs(1,p) = omsgsleft(p)
             call MPI_SSEND(msgs(1,p),  & 
     &                      omsgsleft(p)*(msglen+1)+1, & 
     &                      MPI_DOUBLE_PRECISION, & 
     &                      sproc(i), & 
     &                      mype, & 
     &                      MPI_COMM_WORLD, & 
     &                      ierr)
             if (ierr .ne. 0) call Driver_abortFlash('b_int: bsend 2')
#ifdef DEBUG_BATCH
             print *,mype,' sent ',omsgsin(p),' to ',sproc(i)
             do j=1,omsgsin(p)
             print *,mype,'      [',msgs((j-1)*(msglen+1)+2,p),':' & 
     &              ,(msgs((j-1)*(msglen+1)+2+k,p),k=1,msglen),']'
             enddo
#endif
             omsgsleft(p) = 0
             omsgsin(p) = 0
             nmsgssent = nmsgssent + 1
          endif

#ifdef ERRROR_CHECK
          if (nmsgssent .gt. (nbuffs+1)) then
            call Driver_abortFlash('b_int_sendrcv: nmsgs gt num buffs')
          endif
#endif

       enddo
#ifdef ERRROR_CHECK
       if (nmsgssent .ne. (nbuffs+1)) then
         call Driver_abortFlash('b_int_sendrcv: nmsgs ne num buffs')
       endif
#endif
!CCC
!   DO RECIEVES PHASE
!CCC
       
!
!   While there are messags outstanding,
!
!       WAIT for one
!       parse it into source and tag/val pairs
!       match them against incoming source/tags
!       store the value along with the appropriate tag
!       next.
!   

       do m=1,nmsgsrecv 
           call MPI_WAITANY(nmsgsrecv,rreqs,n,stat,istat)

           p = procno(n)                  ! proc num
           i = neighnum(p+1)              ! index into neighbor list

           npairs = rmsgs(1,n)            ! get # of pairs from message.

#ifdef DEBUG_BATCH
            print *,mype, 'Msg ',n,': Got ',npairs,' from ',p
#endif
           imsgsleft(i) = imsgsleft(i) - npairs
#ifdef ERRROR_CHECK
           if (imsgsleft(i) .lt. 0) then 
               print *,mype, 'Too many pairs from ',p,'!'
           endif
#endif

           do pair=1,npairs
              tag = INT(rmsgs((pair-1)*(msglen+1) + 2 ,n))

              j = 1
              do while ( ((rproc(j).ne.p).or.(rtag(j).ne.tag) & 
     &                   .or.rdone(j)) .and. (j.le.nr) )
                 j = j + 1
              enddo

              if (j .gt. nr) then 
                print *,mype,': ',p,tag,'Tag not found!!'
              endif
              if (rtag(j) .ne. tag) then 
                print *,mype,tag,rtag(j),'Tag not found!!'
              endif
                 
              rdone(j) = .true.
              do k=1,msglen
               rval( msglen*(j-1)+k )=rmsgs((pair-1)*(msglen+1)+k+2,n)
              enddo

#ifdef DEBUG_BATCH
         print *,' [',rtag(j),(rval(msglen*(j-1)+k+2),k=1,msglen,n),']'
#endif
           enddo

       enddo

#ifdef ERRROR_CHECK
       do i=1,nrproc
           if (imsgsleft(i).gt.0) then 
               print *,mype,'Too few pairs from ',rproc(i),'!'
           endif
       enddo
#endif       

#ifdef VERBOSE_BATCH
       print *,mype,': ACTUAL SENDS/RECVS', nmsgssent-1, nmsgsrecv
#endif


!CCC
!   CLEANUP
!CCC

       deallocate(msgs)
       deallocate(rmsgs)
       deallocate(oproc)
       deallocate(iproc)
       deallocate(procno)
       deallocate(omsgs)
       deallocate(imsgs)
       deallocate(nmsg)
       deallocate(neighnum)
       deallocate(sneighnum)

       deallocate(imsgsleft)
       deallocate(rreqs)
       deallocate(ibatches)
       deallocate(rcount)
       deallocate(omsgsin)
       deallocate(omsgsleft)

       deallocate(rdone)

       return
       end

