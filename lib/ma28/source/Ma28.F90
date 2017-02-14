!!****if* source/physics/utilities/solvers/LinearAlgebra/Ma28/Ma28
!!
!! NAME
!!
!!  ma28
!!
!! SYNOPSIS
!!
!!    ma28ad(n,nz,a,licn,irn,lirn,icn,u,ikeep,iw,w,iflag)
!!  
!! DESCRIPTION 
!!
!! this file contains the harwell ma28 sparse matrix routines
!! 
!! easy to use front end routines: 
!! routine ma28ad does the symbolic and numeric lu decomp 
!! routine ma28bd does the numeric lu decomp of ma28ad
!! routine ma28cd solves the system of equations directly
!! routine ma28id solves the system iteratively
!! routine ma28dd does some pointer work
!! 
!! these are the hardball routines:
!! routine ma30ad does core symbolic and numeric lu decomp 
!! routine ma30bd does the numeric decomp the sparse pattern
!! routine ma30cd solves the linear system
!! routine ma30dd does garbage collection
!! 
!! support hardball routines:
!! routine ma28int1 does some common block initialization
!! roytine ma28int2 does some common block initialization
!! routine ma28int3 does some common block initialization
!! routine mc20ad sort a matrix
!! routine mc23ad does the block triangularization pointers
!! routine mc22ad reorders the off diagonal blocks based on the pivot info
!! routine mc21a front end of mc21b
!! routine mc21b pernutes the rows to get a zero free diagonal
!! routine mc13d front end for mc13e
!! routine mc13e permutes a lower triangular block
!! routine mc24ad gets a growth rate of fillin
!! 
!! initialization routines (was block data routines)
!! routine ma28int1 initializes the ma28 routine flags
!! routine ma28int2 initializes the ma28 routine flags
!! routine ma28int3 initializes the ma28 routine flags
!! 
!! never called:
!! routine mc20bd
!!***


      subroutine ma28ad(n,nz,a,licn,irn,lirn,icn,u,ikeep,iw,w,iflag)
      implicit none
      save 
!..
!..this subroutine performs the lu factorization of a.
!..
!..input:
!..n     order of matrix ... not altered by subroutine
!..nz    number of non-zeros in input matrix ... not altered by subroutine
!..a     is a real array  length licn.  holds non-zeros of matrix on entry
!..      and non-zeros of factors on exit.  reordered by mc20a/ad and
!..      mc23a/ad and altered by ma30a/ad
!..licn  integer  length of arrays a and icn ... not altered by subroutine
!..irn   integer array of length lirn.  holds row indices on input.
!..      used as workspace by ma30a/ad to hold column orientation of matrix
!..lirn  integer  length of array irn ... not altered by the subroutine
!..icn   integer array of length licn.  holds column indices on entry
!..      and column indices of decomposed matrix on exit. reordered by
!..      mc20a/ad and mc23a/ad and altered by ma30a/ad.
!..u     real variable  set by user to control bias towards numeric or
!..      sparsity pivoting.  u=1.0 gives partial pivoting while u=0. does
!..      not check multipliers at all.  values of u greater than one are
!..      treated as one while negative values are treated as zero.  not
!..      altered by subroutine.
!..ikeep integer array of length 5*n  used as workspace by ma28a/ad
!..      it is not required to be set on entry and, on exit, it contains 
!..      information about the decomposition. it should be preserved between 
!..      this call and subsequent calls to ma28b/bd or ma28c/cd.
!..      ikeep(i,1),i=1,n  holds the total length of the part of row i
!..      in the diagonal block.
!..      row ikeep(i,2),i=1,n  of the input matrix is the ith row in
!..      pivot order.
!..      column ikeep(i,3),i=1,n  of the input matrix is the ith column
!..      in pivot order.
!..      ikeep(i,4),i=1,n  holds the length of the part of row i in
!..      the l part of the l/u decomposition.
!..      ikeep(i,5),i=1,n  holds the length of the part of row i in the
!..      off-diagonal blocks.  if there is only one diagonal block,
!..      ikeep(1,5) will be set to -1.
!..iw    integer array of length 8*n.  if the option nsrch .le. n is
!..      used, then the length of array iw can be reduced to 7*n.
!..w     real array  length n.  used by mc24a/ad both as workspace and to
!..      return growth estimate in w(1).  the use of this array by ma28a/ad
!..      is thus optional depending on common block logical variable grow.
!..iflag integer variable  used as error flag by routine.  a positive
!..      or zero value on exit indicates success.  possible negative
!..      values are -1 through -14.
!..
!..declare
      integer          n,nz,licn,lirn,iflag,irn(lirn),icn(licn), & 
     &                 ikeep(n,5),iw(n,8),i,j1,j2,jj,j,length, & 
     &                 newpos,move,newj1,jay,knum,ii,i1,iend
      real             a(licn),u,w(n)
!..
!..common and private variables. common block ma28f/fd is used merely
!..to communicate with common block ma30f/fd  so that the user
!..need not declare this common block in his main program.
!..
!..the common block variables are:
!..lp,mp    default value 6 (line printer).  unit number
!..         for error messages and duplicate element warning resp.
!..nlp,mlp  unit number for messages from ma30a/ad and
!..         mc23a/ad resp.  set by ma28a/ad to value of lp.
!..lblock   logical  default value true.  if true mc23a/ad is used
!..         to first permute the matrix to block lower triangular form.
!..grow     logical  default value true.  if true then an estimate
!..         of the increase in size of matrix elements during l/u
!..         decomposition is given by mc24a/ad.
!..eps,rmin,resid  variables not referenced by ma28a/ad.
!..irncp,icncp  set to number of compresses on arrays irn and icn/a 
!..minirn,minicn  minimum length of arrays irn and icn/a; for success on 
!..               future runs.
!..irank  integer   estimated rank of matrix.
!..mirncp,micncp,mirank,mirn,micn communicate between ma30f/fd and ma28f/fd 
!..                               values of above named variables with 
!..                               somewhat similar names.
!..abort1,abort2  logical variables with default value true.  if false
!..               then decomposition will be performed even if the matrix is
!..               structurally or numerically singular respectively.
!..aborta,abortb  logical variables used to communicate values of
!                 abort1 and abort2 to ma30a/ad.
!..abort  logical  used to communicate value of abort1 to mc23a/ad.
!..abort3  logical variable not referenced by ma28a/ad.
!..idisp  integer array  length 2.  used to communicate information
!..       on decomposition between this call to ma28a/ad and subsequent
!..       calls to ma28b/bd and ma28c/cd.  on exit, idisp(1) and
!..       idisp(2) indicate position in arrays a and icn of the
!..       first and last elements in the l/u decomposition of the
!..       diagonal blocks, respectively.
!..numnz  integer  structural rank of matrix.
!..num    integer  number of diagonal blocks.
!..large  integer  size of largest diagonal block.
!..
!..
      logical          grow,lblock,abort,abort1,abort2,abort3,aborta, & 
     &                 abortb,lbig,lbig1
      integer          idisp(2),lp,mp,irncp,icncp,minirn,minicn, & 
     &                 irank,ndrop,maxit,noiter,nsrch,istart, & 
     &                 ndrop1,nsrch1,nlp,mirncp,micncp,mirank, & 
     &                 mirn,micn,mlp,numnz,num,large,lpiv(10), & 
     &                 lnpiv(10),mapiv,manpiv,iavpiv,ianpiv,kountl, & 
     &                 ifirst
      real             tol,themax,big,dxmax,errmax,dres,cgce, & 
     &                 tol1,big1,upriv,rmin,eps,resid,zero
!..
      common /ma28ed/  lp,mp,lblock,grow
      common /ma28fd/  eps,rmin,resid,irncp,icncp,minirn,minicn, & 
     &                 irank,abort1,abort2
      common /ma28gd/  idisp
      common /ma28hd/  tol,themax,big,dxmax,errmax,dres,cgce, & 
     &                 ndrop,maxit,noiter,nsrch,istart,lbig
      common /ma30id/  tol1,big1,ndrop1,nsrch1,lbig1
      common /ma30ed/  nlp,aborta,abortb,abort3
      common /ma30fd/  mirncp,micncp,mirank,mirn,micn
      common /mc23bd/  mlp,numnz,num,large,abort
      common /lpivot/  lpiv,lnpiv,mapiv,manpiv,iavpiv,ianpiv,kountl
      data zero        /0.0e0/, ifirst/0/
!..
!..format statements
99990 format(1x,'error return from ma28a/ad because')
99991 format(1x,'error return from ma30a/ad')
99992 format(1x,'error return from mc23a/ad')
99993 format(1x,'duplicate element in position',i8,' ',i8, & 
     &          'with value ',1pe22.14)
99994 format (1x,i6,'element with value',1pe22.14,'is out of range',/, & 
     &        1x,'with indices',i8,' ',i8)
99995 format(1x,'error return from ma28a/ad; indices out of range')
99996 format(1x,'lirn too small = ',i10)
99997 format(1x,'licn too small = ',i10)
99998 format(1x,'nz non positive = ',i10)
99999 format(1x,'n out of range = ',i10)
!..
!..
!..initialization and transfer of information between common blocks
      if (ifirst .eq. 0) then
       ifirst = 1
       call ma28int1
       call ma28int2
       call ma28int3
      end if
      iflag  = 0
      aborta = abort1
      abortb = abort2
      abort  = abort1
      mlp    = lp
      nlp    = lp
      tol1   = tol
      lbig1  = lbig
      nsrch1 = nsrch
!..
!..upriv private copy of u is used in case it is outside
      upriv = u
!..
!..simple data check on input variables and array dimensions.
      if (n .gt. 0) go to 10
      iflag = -8
      if (lp .ne. 0) write (lp,99999) n
      go to 210
10    if (nz .gt. 0) go to 20
      iflag = -9
      if (lp .ne. 0) write (lp,99998) nz
      go to 210
20    if (licn .ge. nz) go to 30
      iflag = -10
      if (lp .ne. 0) write (lp,99997) licn
      go to 210
30    if (lirn .ge. nz) go to 40
      iflag = -11
      if (lp .ne. 0) write (lp,99996) lirn
      go to 210
!..
!..data check to see if all indices lie between 1 and n.
40    do 50 i=1,nz
       if (irn(i) .gt. 0 .and. irn(i) .le. n .and. icn(i) .gt. 0 .and. & 
     &     icn(i) .le. n) go to 50
       if (iflag .eq. 0 .and. lp .ne. 0) write (lp,99995)
       iflag = -12
       if (lp .ne. 0) write (lp,99994) i,a(i),irn(i),icn(i)
50    continue
      if (iflag .lt. 0) go to 220
!..
!..sort matrix into row order.
      call mc20ad(n,nz,a,icn,iw,irn,0)
!..
!..part of ikeep is used here as a work-array.  ikeep(i,2) is the last row to 
!..have a non-zero in column i.  ikeep(i,3) is the off-set of column i from 
!..the start of the row.
      do 60 i=1,n
       ikeep(i,2) = 0
       ikeep(i,1) = 0
60    continue
!..
!..check for duplicate elements .. summing any such entries and printing a 
!..warning message on unit mp. move is equal to the number of duplicate 
!..elements found; largest element in the matrix is themax; j1 is position in 
!..arrays of first non-zero in row.
      move   = 0
      themax = zero
      j1     = iw(1,1)
      do 130 i=1,n
       iend = nz + 1
       if (i .ne. n) iend = iw(i+1,1)
       length = iend - j1
       if (length .eq. 0) go to 130
       j2 = iend - 1
       newj1 = j1 - move
       do 120 jj=j1,j2
        j = icn(jj)
        themax = max(themax,abs(a(jj)))
        if (ikeep(j,2) .eq. i) go to 110
!..
!..first time column has ocurred in current row.
        ikeep(j,2) = i
        ikeep(j,3) = jj - move - newj1
        if (move .eq. 0) go to 120
!..
!..shift necessary because of previous duplicate element.
        newpos = jj - move
        a(newpos) = a(jj)
        icn(newpos) = icn(jj)
        go to 120
!..
!..duplicate element.
110     move = move + 1
        length = length - 1
        jay = ikeep(j,3) + newj1
        if (mp .ne. 0) write (mp,99993) i,j,a(jj)
        a(jay) = a(jay) + a(jj)
        themax = max(themax,abs(a(jay)))
120    continue
       ikeep(i,1) = length
       j1 = iend
130    continue
!..
!..knum is actual number of non-zeros in matrix with any multiple entries 
!..counted only once
      knum = nz - move
      if (.not.lblock) go to 140
!..
!..perform block triangularisation.
      call mc23ad(n,icn,a,licn,ikeep,idisp,ikeep(1,2), & 
     &            ikeep(1,3),ikeep(1,5),iw(1,3),iw)
      if (idisp(1) .gt. 0) go to 170
      iflag = -7
      if (idisp(1) .eq. -1) iflag = -1
      if (lp .ne. 0) write (lp,99992)
      go to 210
!..
!..block triangularization not requested. move structure to end of data arrays 
!..in preparation for ma30a/ad; set lenoff(1) to -1 and set permutation arrays.
140   do 150 i=1,knum
       ii = knum - i + 1
       newpos = licn - i + 1
       icn(newpos) = icn(ii)
       a(newpos) = a(ii)
150   continue
      idisp(1) = 1
      idisp(2) = licn - knum + 1
      do 160 i=1,n
       ikeep(i,2) = i
       ikeep(i,3) = i
160   continue
      ikeep(1,5) = -1
170   if (lbig) big1 = themax
      if (nsrch .le. n) go to 180
!..
!..perform l/u decomosition on diagonal blocks.
      call ma30ad(n,icn,a,licn,ikeep,ikeep(1,4),idisp, & 
     &           ikeep(1,2),ikeep(1,3),irn,lirn,iw(1,2),iw(1,3),iw(1,4), & 
     &           iw(1,5),iw(1,6),iw(1,7),iw(1,8),iw,upriv,iflag)
      go to 190
!..
!..this call if used if nsrch has been set less than or equal n; in this case, 
!..two integer work arrays of length can be saved.
180    call ma30ad(n,icn,a,licn,ikeep,ikeep(1,4),idisp, & 
     &           ikeep(1,2),ikeep(1,3),irn,lirn,iw(1,2),iw(1,3),iw(1,4), & 
     &           iw(1,5),iw,iw,iw(1,6),iw,upriv,iflag)
!..
!..transfer common block information.
190   minirn = max0(mirn,nz)
      minicn = max0(micn,nz)
      irncp = mirncp
      icncp = micncp
      irank = mirank
      ndrop = ndrop1
      if (lbig) big = big1
      if (iflag .ge. 0) go to 200
      if (lp .ne. 0) write (lp,99991)
      go to 210
!..
!..reorder off-diagonal blocks according to pivot permutation.
200   i1 = idisp(1) - 1
      if (i1 .ne. 0) call mc22ad(n,icn,a,i1,ikeep(1,5),ikeep(1,2), & 
     &                         ikeep(1,3),iw,irn)
      i1 = idisp(1)
      iend = licn - i1 + 1
!..
!..optionally calculate element growth estimate.
      if (grow) call mc24ad(n,icn,a(i1),iend,ikeep,ikeep(1,4),w)
!..
!..increment growth estimate by original maximum element.
      if (grow) w(1) = w(1) + themax
      if (grow .and. n .gt. 1) w(2) = themax
!..
!..set flag if the only error is due to duplicate elements.
      if (iflag .ge. 0 .and. move .ne. 0) iflag = -14
      go to 220
210   if (lp .ne. 0) write (lp,99990)
220   return
      end
!..
!..
!..
!..
!..
      subroutine ma28bd(n,nz,a,licn,ivect,jvect,icn,ikeep,iw,w,iflag)
      implicit none
      save 
!..
!..this subroutine factorizes a matrix with the same pattern as that
!..previously factorized by ma28a/ad.
!..
!..input is :
!..n      order of matrix  not altered by subroutine.
!..nz     number of non-zeros in input matrix  not altered by subroutine.
!..a      array  length licn.  holds non-zeros of matrix on entry and 
!..       non-zeros of factors on exit.  reordered by ma28d/dd and altered by 
!..       subroutine ma30b/bd.
!..licn   integer  length of arrays a and icn.  not altered by subroutine.
!..ivect,jvect  integer arrays of length nz.  hold row and column
!..       indices of non-zeros respectively.  not altered by subroutine.
!..icn    integer array of length licn.  same array as output from ma28a/ad.  
!..       unchanged by ma28b/bd.
!..ikeep  integer array of length 5*n.  same array as output from
!..       ma28a/ad.  unchanged by ma28b/bd.
!..iw     integer array  length 5*n.  used as workspace by ma28d/dd and
!..       ma30b/bd.
!..w      array  length n.  used as workspace by ma28d/dd,ma30b/bd and 
!..       (optionally) mc24a/ad.
!..iflag  error flag with positive or zero value indicating success.
!..
!..
!..declare
      integer          n,nz,licn,iw(n,5),iflag,ikeep(n,5),ivect(nz), & 
     &                 jvect(nz),icn(licn),i1,iend,idup
      real             a(licn),w(n)
!..
!..private and common variables: unless otherwise stated common block 
!..variables are as in ma28a/ad. those variables referenced by ma28b/bd are 
!..mentioned below.
!..
!..lp,mp  used as in ma28a/ad as unit number for error and
!..       warning messages, respectively.
!..nlp    variable used to give value of lp to ma30e/ed.
!..eps    real/real  ma30b/bd will output a positive value
!..       for iflag if any modulus of the ratio of pivot element to the
!..       largest element in its row (u part only) is less than eps (unless
!..       eps is greater than 1.0 when no action takes place).
!..rmin   variable equal to the value of this minimum ratio in cases where 
!..       eps is less than or equal to 1.0.
!..meps,mrmin variables used by the subroutine to communicate between common 
!..        blocks ma28f/fd and ma30g/gd.
!..
!..declare
      logical          grow,lblock,aborta,abortb,abort1,abort2,abort3, & 
     &                 lbig,lbig1
      integer          idisp(2),mp,lp,irncp,icncp,minirn,minicn,irank, & 
     &                 ndrop,maxit,noiter,nsrch,istart,nlp,ndrop1,nsrch1
      real             eps,meps,rmin,mrmin,resid,tol,themax,big,dxmax, & 
     &                 errmax,dres,cgce,tol1,big1
      common /ma28ed/  mp,lp,lblock,grow
      common /ma28fd/  eps,rmin,resid,irncp,icncp,minirn,minicn,irank, & 
     &                 abort1,abort2
      common /ma28gd/  idisp
      common /ma28hd/  tol,themax,big,dxmax,errmax,dres,cgce,ndrop, & 
     &                 maxit,noiter,nsrch,istart,lbig
      common /ma30ed/  nlp,aborta,abortb,abort3
      common /ma30gd/  meps,mrmin
      common /ma30id/  tol1,big1,ndrop1,nsrch1,lbig1
!..
!..formats
99994 format(1x,'error return from ma28b/bd because')
99995 format(1x,'error return from ma30b/bd')
99996 format(1x,'licn too small = ',i10)
99997 format(1x,'nz non positive = ',i10)
99998 format(1x,'n out of range = ',i10)
99999 format(1x,'error return from ma28b/bd with iflag=',i4,/, & 
     &       1x,i7,' entries dropped from structure by ma28a/ad')
!..
!..
!..check to see if elements were dropped in previous ma28a/ad call.
      if (ndrop .eq. 0) go to 10
      iflag = -15
      write (6,99999) iflag,ndrop
      go to 70
10    iflag = 0
      meps  = eps
      nlp   = lp
!..
!..simple data check on variables.
      if (n .gt. 0) go to 20
      iflag = -11
      if (lp .ne. 0) write (lp,99998) n
      go to 60
20    if (nz .gt. 0) go to 30
      iflag = -10
      if (lp .ne. 0) write (lp,99997) nz
      go to 60
30    if (licn .ge. nz) go to 40
      iflag = -9
      if (lp .ne. 0) write (lp,99996) licn
      go to 60
!..
!..
40     call ma28dd(n,a,licn,ivect,jvect,nz,icn,ikeep,ikeep(1,4), & 
     &             ikeep(1,5),ikeep(1,2),ikeep(1,3),iw(1,3),iw, & 
     &             w(1),iflag)
!..
!..themax is largest element in matrix
      themax = w(1)
      if (lbig) big1 = themax
!..
!..idup equals one if there were duplicate elements, zero otherwise.
      idup = 0
      if (iflag .eq. (n+1)) idup = 1
      if (iflag .lt. 0) go to 60
!..
!..perform row-gauss elimination on the structure received from ma28d/dd
      call ma30bd(n,icn,a,licn,ikeep,ikeep(1,4),idisp, & 
     &            ikeep(1,2),ikeep(1,3),w,iw,iflag)
!..
!..transfer common block information.
      if (lbig) big1 = big
      rmin = mrmin
      if (iflag .ge. 0) go to 50
      iflag = -2
      if (lp .ne. 0) write (lp,99995)
      go to 60
!..
!..optionally calculate the growth parameter.
50    i1   = idisp(1)
      iend = licn - i1 + 1
      if (grow) call mc24ad(n,icn,a(i1),iend,ikeep,ikeep(1,4),w)
!..
!..increment estimate by largest element in input matrix.
      if (grow) w(1) = w(1) + themax
      if (grow .and. n .gt. 1) w(2) = themax
!..
!..set flag if the only error is due to duplicate elements.
      if (idup .eq. 1 .and. iflag .ge. 0) iflag = -14
      go to 70
60    if (lp .ne. 0) write (lp,99994)
70    return
      end
!..
!..
!..
!..
!..
      subroutine ma28cd(n,a,licn,icn,ikeep,rhs,w,mtype)
      implicit none
      save 
!..
!..uses the factors from ma28a/ad or ma28b/bd to solve a system of equations
!..
!..input:
!..n     order of matrix  not altered by subroutine.
!..a     array  length licn.  the same array as most recent call to ma28a/ad 
!..      or ma28b/bd.
!..licn  length of arrays a and icn.  not altered by subroutine.
!..icn   integer array of length licn.  same array as output from ma28a/ad.  
!..      unchanged by ma28c/cd.
!..ikeep integer array of length 5*n.  same array as output from ma28a/ad.  
!..      unchanged by ma28c/cd.
!..rhs   array  length n.  on entry, it holds the right hand side.  
!..      on exit, the solution vector.
!..w     array  length n. used as workspace by ma30c/cd.
!..mtype integer  used to tell ma30c/cd to solve the direct equation
!..      (mtype=1) or its transpose (mtype .ne. 1).
!..
!..resid  variable returns maximum residual of equations where pivot was zero.
!..mresid variable used by ma28c/cd to communicate with ma28f/fd and ma30h/hd.
!..idisp  integer array ; the same as that used by ma28a/ad. un changed.
!..
!..declare
      logical          abort1,abort2
      integer          n,licn,idisp(2),icn(licn),ikeep(n,5), & 
     &                 irncp,icncp,minirn,minicn,irank,mtype
      real             a(licn),rhs(n),w(n),resid,mresid,eps,rmin
      common /ma28fd/  eps,rmin,resid,irncp,icncp,minirn,minicn, & 
     &                 irank,abort1,abort2
      common /ma28gd/  idisp
      common /ma30hd/  mresid
!..
!..this call performs the solution of the set of equations.
      call ma30cd(n,icn,a,licn,ikeep,ikeep(1,4),ikeep(1,5), & 
     &            idisp,ikeep(1,2),ikeep(1,3),rhs,w,mtype)
!..
!..transfer common block information.
      resid = mresid
      return
      end
!..
!..
!..
!..
!..
      subroutine ma28id(n,nz,aorg,irnorg,icnorg,licn,a,icn, & 
     &                  ikeep,rhs,x,r,w,mtype,prec,iflag)
      implicit none
      save 
!..
!..this subroutine uses the factors from an earlier call to ma28a/ad
!..or ma28b/bd to solve the system of equations with iterative refinement.
!..
!..parameters are:
!..
!..n    order of the matrix. it is not altered by the subroutine.
!..nz   number of entries in the original matrix.  not altered by subroutine.
!..     for this entry the original matrix must have been saved in
!..     aorg,irnorg,icnorg where entry aorg(k) is in row irnorg(k) and
!..     column icnorg(k), k=1,...nz.  information about the factors of a
!..     is communicated to this subroutine via the parameters licn, a, icn
!..     and ikeep where:
!..aorg   array of length nz.  not altered by ma28i/id.
!..irnorg array of length nz.  not altered by ma28i/id.
!..icnorg array of length nz.  not altered by ma28i/id.
!..licn   equal to the length of arrays a and icn. not altered
!..a    array of length licn. it must be unchanged since the last call
!..     to ma28a/ad or ma28b/bd. it is not altered by the subroutine.
!..icn, ikeep are the arrays (of lengths licn and 5*n, respectively) of
!..     the same names as in the previous all to ma28a/ad. they should be
!..     unchanged since this earlier call. not altered.
!..
!..other parameters are as follows:
!..rhs array of length n. the user must set rhs(i) to contain the
!..    value of the i th component of the right hand side. not altered.
!..
!..x   array of length n. if an initial guess of the solution is
!..    given (istart equal to 1), then the user must set x(i) to contain
!..    the value of the i th component of the estimated solution.  on
!..    exit, x(i) contains the i th component of the solution vector.
!..r   array of length n. it need not be set on entry.  on exit, r(i)
!..    contains the i th component of an estimate of the error if maxit
!..    is greater than 0.
!..w is an array of length n. it is used as workspace by ma28i/id.
!..mtype must be set to determine whether ma28i/id will solve a*x=rhs
!..     (mtype equal to 1) or at*x=rhs (mtype ne 1, zero say). not altered.
!..prec should be set by the user to the relative accuracy required. the
!..     iterative refinement will terminate if the magnitude of the
!..     largest component of the estimated error relative to the largest
!..     component in the solution is less than prec. not altered.
!..iflag is a diagnostic flag which will be set to zero on successful
!..      exit from ma28i/id, otherwise it will have a non-zero value. the
!..      non-zero value iflag can have on exit from ma28i/id are ...
!..      -16    indicating that more than maxit iteartions are required.
!..      -17    indicating that more convergence was too slow.
!..
!..declare
      integer          n,nz,licn,mtype,iflag,icnorg(nz),irnorg(nz), & 
     &                 ikeep(n,5),icn(licn),i,iterat,nrow,ncol
      real             a(licn),aorg(nz),rhs(n),r(n),x(n),w(n),prec, & 
     &                 d,dd,conver,zero
!..
!..common block communication
      logical          lblock,grow,lbig
      integer          lp,mp,ndrop,maxit,noiter,nsrch,istart
      real             tol,themax,big,dxmax,errmax,dres,cgce
      common /ma28ed/  lp,mp,lblock,grow
      common /ma28hd/  tol,themax,big,dxmax,errmax,dres,cgce, & 
     &                 ndrop,maxit,noiter,nsrch,istart,lbig
      data             zero /0.0e0/
!..
!..
!..formats
99998 format(1x,'error return from ma28i with iflag = ', i3,/, & 
     &       1x,'convergence rate of',1pe9.2,'too slow',/, & 
     &       1x,'maximum acceptable rate set to ',1pe9.2)
99999 format(1x,'error return from ma28i/id with iflag = ',i3,/, & 
     &       1x,'more than',i5,'iterations required')
!..
!..
!..initialization of noiter,errmax and iflag.
      noiter = 0
      errmax = zero
      iflag   = 0
!..
!..jump if a starting vector has been supplied by the user.
      if (istart .eq. 1) go to 20
!..
!..make a copy of the right-hand side vector.
      do 10 i=1,n
       x(i) = rhs(i)
10    continue
!..
!..find the first solution.
      call ma28cd(n,a,licn,icn,ikeep,x,w,mtype)
!..
!..stop the computations if   maxit=0.
20    if (maxit .eq. 0) go to 160
!..
!..calculate the max-norm of the first solution.
      dd = 0.0e0
      do 30 i=1,n
       dd = max(dd,abs(x(i)))
30    continue
      dxmax = dd
!..
!..begin the iterative process.
      do 120 iterat=1,maxit
       d = dd
!..
!..calculate the residual vector.
       do 40 i=1,n
        r(i) = rhs(i)
40     continue
       if (mtype .eq. 1) go to 60
       do 50 i=1,nz
        nrow = irnorg(i)
        ncol = icnorg(i)
        r(ncol) = r(ncol) - aorg(i)*x(nrow)
50     continue
       go to 80
!..
!..mtype=1.
60     do 70 i=1,nz
        nrow = irnorg(i)
        ncol = icnorg(i)
        r(nrow) = r(nrow) - aorg(i)*x(ncol)
70     continue
80     dres = 0.0e0
!..
!..find the max-norm of the residual vector.
       do 90 i=1,n
        dres = max(dres,abs(r(i)))
90     continue
!..
!..stop the calculations if the max-norm of the residual vector is zero.
       if (dres .eq. 0.0) go to 150
!..
!..calculate the correction vector.
       noiter = noiter + 1
       call ma28cd(n,a,licn,icn,ikeep,r,w,mtype)
!..
!..find the max-norm of the correction vector.
       dd = 0.0e0
       do 100 i=1,n
        dd = max(dd,abs(r(i)))
100    continue
!..
!..check the convergence.
       if (dd .gt. d*cgce .and. iterat .ge. 2) go to 130
       if (dxmax*10.0e0 + dd .eq. dxmax*10.0e0) go to 140
!..
!..attempt to improve the solution.
       dxmax = 0.0e0
       do 110 i=1,n
        x(i) = x(i) + r(i)
        dxmax = max(dxmax,abs(x(i)))
110    continue
!..
!..check the stopping criterion; end of iteration loop
       if (dd .lt. prec*dxmax) go to 140
120   continue
!..
!..more than maxit iterations required.
      iflag = -16
      write (lp,99999) iflag,maxit
      go to 140
!..
!..convergence rate unacceptably slow.
130   iflag = -17
      conver = dd/d
      write (lp,99998) iflag,conver,cgce
!..
!..the iterative process is terminated.
140   errmax = dd
150   continue
160   return
      end
!..
!..
!..
!..
!..
      subroutine ma28dd(n,a,licn,ivect,jvect,nz,icn,lenr,lenrl, & 
     &                  lenoff,ip,iq,iw1,iw,w1,iflag)
      implicit none
      save 
!..
!..this subroutine need never be called by the user directly.
!..it sorts the user's matrix into the structure of the decomposed
!..form and checks for the presence of duplicate entries or
!..non-zeros lying outside the sparsity pattern of the decomposition
!..it also calculates the largest element in the input matrix.
!..
!..declare
      logical          lblock,grow,blockl
      integer          n,licn,nz,iw(n,2),idisp(2),icn(licn),ivect(nz), & 
     &                 jvect(nz),ip(n),iq(n),lenr(n),iw1(n,3),lenrl(n), & 
     &                 lenoff(n),iflag,lp,mp,i,ii,jj,inew,jnew,iblock, & 
     &                 iold,jold,j1,j2,idisp2,idummy,jdummy,midpt,jcomp
      real             a(licn),zero,w1,aa
!..
      common /ma28ed/  lp,mp,lblock,grow
      common /ma28gd/  idisp
      data             zero /0.0e0/
!..
!..formats
99997 format(1x,'element',i6,' ',i6,' was not in l/u pattern')
99998 format(1x,'non-zero',i7,' ',i6,' in zero off-diagonal block')
99999 format(1x,'element ',i6,' with value ',1pe22.14,/, & 
     &       1x,'has indices',i8,' ',i8,'out of range')
!..
!..iw1(i,3)  is set to the block in which row i lies and the inverse 
!..permutations to ip and iq are set in iw1(.,1) and iw1(.,2) resp.
!..pointers to beginning of the part of row i in diagonal and off-diagonal 
!..blocks are set in iw(i,2) and iw(i,1) resp.
      blockl  = lenoff(1) .ge. 0
      iblock  = 1
      iw(1,1) = 1
      iw(1,2) = idisp(1)
      do 10 i=1,n
       iw1(i,3) = iblock
       if (ip(i) .lt. 0) iblock = iblock + 1
       ii = iabs(ip(i)+0)
       iw1(ii,1) = i
       jj = iq(i)
       jj = iabs(jj)
       iw1(jj,2) = i
       if (i .eq. 1) go to 10
       if (blockl) iw(i,1) = iw(i-1,1) + lenoff(i-1)
       iw(i,2) = iw(i-1,2) + lenr(i-1)
10    continue
!..
!..place each non-zero in turn into its correct location in the a/icn array.
      idisp2 = idisp(2)
      do 170 i=1,nz
       if (i .gt. idisp2) go to 20
       if (icn(i) .lt. 0) go to 170
20     iold = ivect(i)
       jold = jvect(i)
       aa = a(i)
!..
!..dummy loop for following a chain of interchanges. executed nz times.
       do 140 idummy=1,nz
        if (iold .le. n .and. iold .gt. 0 .and. jold .le. n .and.  & 
     &      jold .gt. 0) go to 30
        if (lp .ne. 0) write (lp,99999) i, a(i), iold, jold
        iflag = -12
        go to 180
30      inew = iw1(iold,1)
        jnew = iw1(jold,2)
!..
!..are we in a valid block and is it diagonal or off-diagonal?
        if (iw1(inew,3)-iw1(jnew,3)) 40, 60, 50
40      iflag = -13
        if (lp .ne. 0) write (lp,99998) iold, jold
        go to 180
50      j1 = iw(inew,1)
        j2 = j1 + lenoff(inew) - 1
        go to 110
!..
!..element is in diagonal block.
60      j1 = iw(inew,2)
        if (inew .gt. jnew) go to 70
        j2 = j1 + lenr(inew) - 1
        j1 = j1 + lenrl(inew)
        go to 110
70      j2 = j1 + lenrl(inew)
!..
!..binary search of ordered list  .. element in l part of row.
        do 100 jdummy=1,n
         midpt = (j1+j2)/2
         jcomp = iabs(icn(midpt)+0)
         if (jnew-jcomp) 80, 130, 90
80       j2 = midpt
         go to 100
90       j1 = midpt
100     continue
        iflag = -13
        if (lp .ne. 0) write (lp,99997) iold, jold
        go to 180
!..
!..linear search ... element in l part of row or off-diagonal blocks.
110     do 120 midpt=j1,j2
         if (iabs(icn(midpt)+0) .eq. jnew) go to 130
120     continue
        iflag = -13
        if (lp .ne. 0) write (lp,99997) iold, jold
        go to 180
!..
!..equivalent element of icn is in position midpt.
130     if (icn(midpt) .lt. 0) go to 160
        if (midpt .gt. nz .or. midpt .le. i) go to 150
        w1 = a(midpt)
        a(midpt) = aa
        aa = w1
        iold = ivect(midpt)
        jold = jvect(midpt)
        icn(midpt) = -icn(midpt)
140    continue
!..
150    a(midpt) = aa
       icn(midpt) = -icn(midpt)
       go to 170
160    a(midpt) = a(midpt) + aa
!..
!..set flag for duplicate elements; end of big loop
       iflag = n + 1
170   continue
!..
!..reset icn array  and zero elements in l/u but not in a. get max a element
180   w1 = zero
      do 200 i=1,idisp2
       if (icn(i) .lt. 0) go to 190
       a(i) = zero
       go to 200
190    icn(i) = -icn(i)
       w1 = max(w1,abs(a(i)))
200   continue
      return
      end
!..
!..
!..
!..
!..
      subroutine ma30ad(nn,icn,a,licn,lenr,lenrl,idisp,ip,iq, & 
     &                  irn,lirn,lenc,ifirst,lastr,nextr,lastc, & 
     &                  nextc,iptr,ipc,u,iflag)
      implicit none
      save 
!..
!..
!..if the user requires a more convenient data interface then the ma28
!..package should be used.  the ma28 subroutines call the ma30 routines after 
!..checking the user's input data and optionally using mc23a/ad to permute the 
!..matrix to block triangular form.
!..
!..this package of subroutines (ma30a/ad, ma30b/bd, ma30c/cd and ma30d/dd) 
!..performs operations pertinent to the solution of a general sparse n by n 
!..system of linear equations (i.e. solve ax=b). structually singular matrices 
!..are permitted including those with row or columns consisting entirely of 
!..zeros (i.e. including rectangular matrices).  it is assumed that the 
!..non-zeros of the matrix a do not differ widely in size. if necessary a 
!..prior call of the scaling subroutine mc19a/ad may be made.
!..
!..a discussion of the design of these subroutines is given by duff and reid 
!..(acm trans math software 5 pp 18-35,1979 (css 48)) while fuller details of 
!..the implementation are given in duff (harwell report aere-r 8730,1977).  
!..the additional pivoting option in ma30a/ad and the use of drop tolerances 
!..(see common block ma30i/id) were added to the package after joint work with 
!..duff,reid,schaumburg,wasniewski and zlatev, harwell report css 135, 1983.
!..
!..ma30a/ad performs the lu decomposition of the diagonal blocks of the
!..permutation paq of a sparse matrix a, where input permutations p1 and q1 
!..are used to define the diagonal blocks.  there may be non-zeros in the 
!..off-diagonal blocks but they are unaffected by ma30a/ad. p and p1 differ 
!..only within blocks as do q and q1. the permutations p1 and q1 may be found 
!..by calling mc23a/ad or the matrix may be treated as a single block by 
!..using p1=q1=i. the matrix non-zeros should be held compactly by rows, 
!..although it should be noted that the user can supply the matrix by columns
!..to get the lu decomposition of a transpose.
!..
!..this description of the following parameters should also be consulted for 
!..further information on most of the parameters of ma30b/bd and ma30c/cd:
!..
!..n    is an integer variable which must be set by the user to the order
!..     of the matrix.  it is not altered by ma30a/ad.
!..
!..icn  is an integer array of length licn. positions idisp(2) to
!..     licn must be set by the user to contain the column indices of
!..     the non-zeros in the diagonal blocks of p1*a*q1. those belonging
!..     to a single row must be contiguous but the ordering of column
!..     indices with each row is unimportant. the non-zeros of row i
!..     precede those of row i+1,i=1,...,n-1 and no wasted space is
!..     allowed between the rows.  on output the column indices of the
!..     lu decomposition of paq are held in positions idisp(1) to
!..     idisp(2), the rows are in pivotal order, and the column indices
!..     of the l part of each row are in pivotal order and precede those
!..     of u. again there is no wasted space either within a row or
!..     between the rows. icn(1) to icn(idisp(1)-1), are neither
!..     required nor altered. if mc23a/ad been called, these will hold
!..     information about the off-diagonal blocks.
!..
!..a    is a real/real array of length licn whose entries
!..     idisp(2) to licn must be set by the user to the  values of the
!..     non-zero entries of the matrix in the order indicated by  icn.
!..     on output a will hold the lu factors of the matrix where again
!..     the position in the matrix is determined by the corresponding
!..     values in icn. a(1) to a(idisp(1)-1) are neither required nor altered.
!..
!..licn is an integer variable which must be set by the user to the
!..     length of arrays icn and a. it must be big enough for a and icn
!..     to hold all the non-zeros of l and u and leave some "elbow
!..     room".  it is possible to calculate a minimum value for licn by
!..     a preliminary run of ma30a/ad. the adequacy of the elbow room
!..     can be judged by the size of the common block variable icncp. it
!..     is not altered by ma30a/ad.
!..
!..lenr is an integer array of length n.  on input, lenr(i) should
!..     equal the number of non-zeros in row i, i=1,...,n of the
!..     diagonal blocks of p1*a*q1. on output, lenr(i) will equal the
!..     total number of non-zeros in row i of l and row i of u.
!..
!..lenrl is an integer array of length n. on output from ma30a/ad,
!..      lenrl(i) will hold the number of non-zeros in row i of l.
!..
!..idisp is an integer array of length 2. the user should set idisp(1)
!..      to be the first available position in a/icn for the lu
!..      decomposition while idisp(2) is set to the position in a/icn of
!..      the first non-zero in the diagonal blocks of p1*a*q1. on output,
!..      idisp(1) will be unaltered while idisp(2) will be set to the
!..      position in a/icn of the last non-zero of the lu decomposition.
!..
!..ip    is an integer array of length n which holds a permutation of
!..      the integers 1 to n.  on input to ma30a/ad, the absolute value of
!..      ip(i) must be set to the row of a which is row i of p1*a*q1. a
!..      negative value for ip(i) indicates that row i is at the end of a
!..      diagonal block.  on output from ma30a/ad, ip(i) indicates the row
!..      of a which is the i th row in paq. ip(i) will still be negative
!..      for the last row of each block (except the last).
!..
!..iq    is an integer array of length n which again holds a
!..      permutation of the integers 1 to n.  on input to ma30a/ad, iq(j)
!..      must be set to the column of a which is column j of p1*a*q1. on
!..      output from ma30a/ad, the absolute value of iq(j) indicates the
!..      column of a which is the j th in paq.  for rows, i say, in which
!..      structural or numerical singularity is detected iq(i) is negated.
!..
!..irn  is an integer array of length lirn used as workspace by ma30a/ad.
!..
!..lirn is an integer variable. it should be greater than the
!..     largest number of non-zeros in a diagonal block of p1*a*q1 but
!..     need not be as large as licn. it is the length of array irn and
!..     should be large enough to hold the active part of any block,
!..     plus some "elbow room", the  a posteriori  adequacy of which can
!..     be estimated by examining the size of common block variable irncp.
!..
!..lenc,ifirst,lastr,nextr,lastc,nextc 
!..     are all integer arrays of length n which are used as workspace by 
!..     ma30a/ad.  if nsrch is set to a value less than or equal to n, then 
!..     arrays lastc and nextc are not referenced by ma30a/ad and so can be 
!..     dummied in the call to ma30a/ad.
!..
!..iptr,ipc are integer arrays of length n; used as workspace by ma30a/ad.
!..
!..u    is a real/real variable which should be set by the
!..     user to a value between 0. and 1.0. if less than zero it is
!..     reset to zero and if its value is 1.0 or greater it is reset to
!..     0.9999 (0.999999999 in d version).  it determines the balance
!..     between pivoting for sparsity and for stability, values near
!..     zero emphasizing sparsity and values near one emphasizing
!..     stability. we recommend u=0.1 as a posible first trial value.
!..     the stability can be judged by a later call to mc24a/ad or by
!..     setting lbig to .true.
!..
!..iflag is an integer variable. it will have a non-negative value if
!..      ma30a/ad is successful. negative values indicate error
!..      conditions while positive values indicate that the matrix has
!..      been successfully decomposed but is singular. for each non-zero
!..      value, an appropriate message is output on unit lp.  possible
!..      non-zero values for iflag are
!.. -1   the matrix is structually singular with rank given by irank in
!..      common block ma30f/fd.
!.. +1   if, however, the user wants the lu decomposition of a
!..      structurally singular matrix and sets the common block variable
!..      abort1 to .false., then, in the event of singularity and a
!..      successful decomposition, iflag is returned with the value +1
!..      and no message is output.
!.. -2   the matrix is numerically singular (it may also be structually
!..      singular) with estimated rank given by irank in common block ma30f/fd.
!.. +2   the  user can choose to continue the decomposition even when a
!..      zero pivot is encountered by setting common block variable
!..      abort2 to .false.  if a singularity is encountered, iflag will
!..      then return with a value of +2, and no message is output if the
!..      decomposition has been completed successfully.
!.. -3   lirn has not been large enough to continue with the
!..      decomposition.  if the stage was zero then common block variable
!..      minirn gives the length sufficient to start the decomposition on
!..      this block.  for a successful decomposition on this block the user
!..      should make lirn slightly (say about n/2) greater than this value.
!.. -4   licn not large enough to continue with the decomposition.
!.. -5   the decomposition has been completed but some of the lu factors
!..      have been discarded to create enough room in a/icn to continue
!..      the decomposition. the variable minicn in common block ma30f/fd
!..      then gives the size that licn should be to enable the
!..      factorization to be successful.  if the user sets common block
!..      variable abort3 to .true., then the subroutine will exit
!..      immediately instead of destroying any factors and continuing.
!.. -6   both licn and lirn are too small. termination has been caused by
!..      lack of space in irn (see error iflag= -3), but already some of
!..      the lu factors in a/icn have been lost (see error iflag= -5).
!..      minicn gives the minimum amount of space required in a/icn for
!..      decomposition up to this point.
!..
!..
!..declare
      logical          abort1,abort2,abort3,lbig
      integer          nn,licn,lirn,iptr(nn),pivot,pivend,dispc, & 
     &                 oldpiv,oldend,pivrow,rowi,ipc(nn),idisp(2), & 
     &                 colupd,icn(licn),lenr(nn),lenrl(nn),ip(nn), & 
     &                 iq(nn),lenc(nn),irn(lirn),ifirst(nn),lastr(nn), & 
     &                 nextr(nn),lastc(nn),nextc(nn),lpiv(10),lnpiv(10), & 
     &                 msrch,nsrch,ndrop,kk,mapiv,manpiv,iavpiv,ianpiv, & 
     &                 kountl,minirn,minicn,morei,irank,irncp,icncp, & 
     &                 iflag,ibeg,iactiv,nzrow,num,nnm1,i,ilast,nblock, & 
     &                 istart,irows,n,ising,lp,itop,ii,j1,jj,j,indrow, & 
     &                 j2,ipos,nzcol,nzmin,nz,isw,isw1,jcost, & 
     &                 isrch,ll,ijfir,idummy,kcost,ijpos,ipiv,jpiv,i1, & 
     &                 i2,jpos,k,lc,nc,lr,nr,lenpiv,ijp1,nz2,l,lenpp, & 
     &                 nzpc,iii,idrop,iend,iop,jnew,ifill,jdiff,jnpos, & 
     &                 jmore,jend,jroom,jbeg,jzero,idispc,kdrop, & 
     &                 ifir,jval,jzer,jcount,jdummy,jold
      real             a(licn),u,au,umax,amax,zero,pivrat,pivr, & 
     &                 tol,big,anew,aanew,scale
!..
      common /ma30ed/  lp,abort1,abort2,abort3
      common /ma30fd/  irncp,icncp,irank,minirn,minicn
      common /ma30id/  tol,big,ndrop,nsrch,lbig
      common /lpivot/  lpiv,lnpiv,mapiv,manpiv,iavpiv,ianpiv,kountl
!..
      data             umax/0.999999999e0/
      data             zero /0.0e0/
!..
!..
!..formats
99992 format(1x,'to continue set lirn to at least',i8)
99993 format(1x,'at stage',i5,'in block',i5,'with first row',i5, & 
     &          'and last row',i5)
99994 format(1x,'error return from ma30a/ad lirn and licn too small')
99995 format(1x,'error return from ma30a/ad ; lirn not big enough')
99996 format(1x,'error return from ma30a/ad ; licn not big enough')
99997 format(1x,'lu decomposition destroyed to create more space')
99998 format(1x,'error return from ma30a/ad;', & 
     &          'matrix is numerically singular')
99999 format(1x,'error return from ma30a/ad;', & 
     &          'matrix is structurally singular')
!..
!..
!..
!..initialize
      msrch = nsrch
      ndrop = 0
      do 1272 kk=1,10
       lnpiv(kk) = 0
       lpiv(kk)  = 0
1272  continue
      mapiv  = 0
      manpiv = 0
      iavpiv = 0
      ianpiv = 0
      kountl = 0
      minirn = 0
      minicn = idisp(1) - 1
      morei  = 0
      irank  = nn
      irncp  = 0
      icncp  = 0
      iflag  = 0
      u      = min(u,umax)
      u      = max(u,zero)
!..
!..ibeg is the position of the next pivot row after elimination step using it.
!..iactiv is the position of the first entry in the active part of a/icn.
!..nzrow is current number of non-zeros in active and unprocessed part of row 
!..file icn.
      ibeg   = idisp(1)
      iactiv = idisp(2)
      nzrow  = licn - iactiv + 1
      minicn = nzrow + minicn
!..
!..count the number of diagonal blocks and set up pointers to the beginnings of 
!..the rows. num is the number of diagonal blocks.
      num = 1
      iptr(1) = iactiv
      if (nn .eq. 1) go to 20
      nnm1 = nn - 1
      do 10 i=1,nnm1
       if (ip(i) .lt. 0) num = num + 1
       iptr(i+1) = iptr(i) + lenr(i)
10    continue
!..
!..ilast is the last row in the previous block.
20    ilast = 0
!..
!..lu decomposition of block nblock starts
!..each pass on this loop performs lu decomp on one of the diagonal blocks.
      do 1000 nblock=1,num
       istart = ilast + 1
       do 30 irows=istart,nn
        if (ip(irows) .lt. 0) go to 40
30     continue
       irows = nn
40     ilast = irows
!..
!..n is the number of rows in the current block.
!..istart is the index of the first row in the current block.
!..ilast is the index of the last row in the current block.
!..iactiv is the position of the first entry in the block.
!..itop is the position of the last entry in the block.
       n = ilast - istart + 1
       if (n .ne. 1) go to 90
!..
!..code for dealing with 1x1 block.
       lenrl(ilast) = 0
       ising = istart
       if (lenr(ilast) .ne. 0) go to 50
!..
!..block is structurally singular.
       irank = irank - 1
       ising = -ising
       if (iflag .ne. 2 .and. iflag .ne. -5) iflag = 1
       if (.not.abort1) go to 80
       idisp(2) = iactiv
       iflag = -1
       if (lp .ne. 0) write (lp,99999)
       go to 1120
!..
50     scale = abs(a(iactiv))
       if (scale .eq. zero) go to 60
       if (lbig) big = max(big,scale)
       go to 70
60     ising = -ising
       irank = irank - 1
       iptr(ilast) = 0
       if (iflag .ne. -5) iflag = 2
       if (.not.abort2) go to 70
       idisp(2) = iactiv
       iflag    = -2
       if (lp .ne. 0) write (lp,99998)
       go to 1120
70     a(ibeg)       = a(iactiv)
       icn(ibeg)     = icn(iactiv)
       iactiv        = iactiv + 1
       iptr(istart)  = 0
       ibeg          = ibeg + 1
       nzrow         = nzrow - 1
80     lastr(istart) = istart
       ipc(istart)   = -ising
       go to 1000
!..
!..non-trivial block.
90     itop = licn
       if (ilast .ne. nn) itop = iptr(ilast+1) - 1
!..
!..set up column oriented storage.
       do 100 i=istart,ilast
        lenrl(i) = 0
        lenc(i)  = 0
100    continue
       if (itop-iactiv .lt. lirn) go to 110
       minirn = itop - iactiv + 1
       pivot  = istart - 1
       go to 1100
!..
!..calculate column counts.
110    do 120 ii=iactiv,itop
        i       = icn(ii)
        lenc(i) = lenc(i) + 1
120    continue
!..
!..set up column pointers so that ipc(j) points to position after end of 
!..column j in column file.
       ipc(ilast) = lirn + 1
       j1         = istart + 1
       do 130 jj=j1,ilast
        j      = ilast - jj + j1 - 1
        ipc(j) = ipc(j+1) - lenc(j+1)
130    continue
       do 150 indrow=istart,ilast
        j1 = iptr(indrow)
        j2 = j1 + lenr(indrow) - 1
        if (j1 .gt. j2) go to 150
        do 140 jj=j1,j2
         j         = icn(jj)
         ipos      = ipc(j) - 1
         irn(ipos) = indrow
         ipc(j)    = ipos
140     continue
150    continue
!..
!..dispc is the lowest indexed active location in the column file.
       dispc  = ipc(istart)
       nzcol  = lirn - dispc + 1
       minirn = max0(nzcol,minirn)
       nzmin  = 1
!..
!..initialize array ifirst.  ifirst(i) = +/- k indicates that row/col k has i 
!..non-zeros.  if ifirst(i) = 0,there is no row or column with i non zeros.
       do 160 i=1,n
        ifirst(i) = 0
160    continue
!..
!..compute ordering of row and column counts. first run through columns (from 
!..column n to column 1).
       do 180 jj=istart,ilast
        j  = ilast - jj + istart
        nz = lenc(j)
        if (nz .ne. 0) go to 170
        ipc(j) = 0
        go to 180
170     if (nsrch .le. nn) go to 180
        isw        = ifirst(nz)
        ifirst(nz) = -j
        lastc(j)   = 0
        nextc(j)   = -isw
        isw1       = iabs(isw)
        if (isw .ne. 0) lastc(isw1) = j
180    continue
!..
!..now run through rows (again from n to 1).
       do 210 ii=istart,ilast
        i  = ilast - ii + istart
        nz = lenr(i)
        if (nz .ne. 0) go to 190
        iptr(i)  = 0
        lastr(i) = 0
        go to 210
190     isw        = ifirst(nz)
        ifirst(nz) = i
        if (isw .gt. 0) go to 200
        nextr(i) = 0
        lastr(i) = isw
        go to 210
200     nextr(i)   = isw
        lastr(i)   = lastr(isw)
        lastr(isw) = i
210    continue
!..
!..
!..start of main elimination loop
!..
!..first find the pivot using markowitz criterion with stability control.
!..jcost is the markowitz cost of the best pivot so far,.. this pivot is in 
!..row ipiv and column jpiv.
!..
       do 980 pivot=istart,ilast
        nz2   = nzmin
        jcost = n*n
!..
!..examine rows/columns in order of ascending count.
        do 340 l=1,2
         pivrat = zero
         isrch  = 1
         ll     = l
!..
!..a pass with l equal to 2 is only performed in the case of singularity.
         do 330 nz=nz2,n
          if (jcost .le. (nz-1)**2) go to 420
          ijfir = ifirst(nz)
          if (ijfir) 230, 220, 240
220       if (ll .eq. 1) nzmin = nz + 1
          go to 330
230       ll    = 2
          ijfir = -ijfir
          go to 290
240       ll = 2
!..
!.. scan rows with nz non-zeros.
          do 270 idummy=1,n
           if (jcost .le. (nz-1)**2) go to 420
           if (isrch .gt. msrch) go to 420
           if (ijfir .eq. 0) go to 280
!..
!..row ijfir is now examined.
           i     = ijfir
           ijfir = nextr(i)
!..
!..first calculate multiplier threshold level.
           amax = zero
           j1   = iptr(i) + lenrl(i)
           j2   = iptr(i) + lenr(i) - 1
           do 250 jj=j1,j2
            amax = max(amax,abs(a(jj)))
250        continue
           au    = amax*u
           isrch = isrch + 1
!..
!..scan row for possible pivots
           do 260 jj=j1,j2
            if (abs(a(jj)) .le. au .and. l .eq. 1) go to 260
            j     = icn(jj)
            kcost = (nz-1)*(lenc(j)-1)
            if (kcost .gt. jcost) go to 260
            pivr = zero
            if (amax .ne. zero) pivr = abs(a(jj))/amax
            if (kcost .eq. jcost .and. (pivr .le. pivrat .or.  & 
     &          nsrch .gt. nn+1)) go to 260
!..
!..best pivot so far is found.
            jcost = kcost
            ijpos = jj
            ipiv  = i
            jpiv  = j
            if (msrch .gt. nn+1 .and. jcost .le. (nz-1)**2) go to 420
            pivrat = pivr
260        continue
270       continue
!..
!..columns with nz non-zeros now examined.
280       ijfir = ifirst(nz)
          ijfir = -lastr(ijfir)
290       if (jcost .le. nz*(nz-1)) go to 420
          if (msrch .le. nn) go to 330
          do 320 idummy=1,n
           if (ijfir .eq. 0) go to 330
           j     = ijfir
           ijfir = nextc(ijfir)
           i1    = ipc(j)
           i2    = i1 + nz - 1
!..
!..scan column j
           do 310 ii=i1,i2
            i     = irn(ii)
            kcost = (nz-1)*(lenr(i)-lenrl(i)-1)
            if (kcost .ge. jcost) go to 310
!..
!..pivot has best markowitz count so far ... now check its suitability on 
!..numeric grounds by examining the other non-zeros in its row.
            j1 = iptr(i) + lenrl(i)
            j2 = iptr(i) + lenr(i) - 1
!..
!..we need a stability check on singleton columns because of possible problems 
!..with underdetermined systems.
            amax = zero
            do 300 jj=j1,j2
             amax = max(amax,abs(a(jj)))
             if (icn(jj) .eq. j) jpos = jj
300         continue
            if (abs(a(jpos)) .le. amax*u .and. l .eq. 1) go to 310
            jcost = kcost
            ipiv  = i
            jpiv  = j
            ijpos = jpos
            if (amax .ne. zero) pivrat = abs(a(jpos))/amax
            if (jcost .le. nz*(nz-1)) go to 420
310        continue
320       continue
330      continue
!..
!..in the event of singularity; must make sure all rows and columns are tested.
!..matrix is numerically or structurally singular; it will be diagnosed later.
         msrch = n
         irank = irank - 1
340     continue
!..
!..assign rest of rows and columns to ordering array. matrix is singular.
        if (iflag .ne. 2 .and. iflag .ne. -5) iflag = 1
        irank = irank - ilast + pivot + 1
        if (.not.abort1) go to 350
        idisp(2) = iactiv
        iflag = -1
        if (lp .ne. 0) write (lp,99999)
        go to 1120
350     k = pivot - 1
        do 390 i=istart,ilast
         if (lastr(i) .ne. 0) go to 390
         k        = k + 1
         lastr(i) = k
         if (lenrl(i) .eq. 0) go to 380
         minicn = max0(minicn,nzrow+ibeg-1+morei+lenrl(i))
         if (iactiv-ibeg .ge. lenrl(i)) go to 360
         call ma30dd(a,icn,iptr(istart),n,iactiv,itop,.true.)
!..
!..check now to see if ma30d/dd has created enough available space.
         if (iactiv-ibeg .ge. lenrl(i)) go to 360
!..
!..create more space by destroying previously created lu factors.
         morei = morei + ibeg - idisp(1)
         ibeg = idisp(1)
         if (lp .ne. 0) write (lp,99997)
         iflag = -5
         if (abort3) go to 1090
360      j1 = iptr(i)
         j2 = j1 + lenrl(i) - 1
         iptr(i) = 0
         do 370 jj=j1,j2
          a(ibeg)   = a(jj)
          icn(ibeg) = icn(jj)
          icn(jj)   = 0
          ibeg      = ibeg + 1
370      continue
         nzrow = nzrow - lenrl(i)
380      if (k .eq. ilast) go to 400
390     continue
400     k = pivot - 1
        do 410 i=istart,ilast
         if (ipc(i) .ne. 0) go to 410
         k      = k + 1
         ipc(i) = k
         if (k .eq. ilast) go to 990
410     continue
!..
!..the pivot has now been found in position (ipiv,jpiv) in location ijpos in 
!..row file. update column and row ordering arrays to correspond with removal
!..of the active part of the matrix.
420     ising = pivot
        if (a(ijpos) .ne. zero) go to 430
!..
!..numerical singularity is recorded here.
        ising = -ising
        if (iflag .ne. -5) iflag = 2
        if (.not.abort2) go to 430
        idisp(2) = iactiv
        iflag = -2
        if (lp .ne. 0) write (lp,99998)
        go to 1120
430     oldpiv = iptr(ipiv) + lenrl(ipiv)
        oldend = iptr(ipiv) + lenr(ipiv) - 1
!..
!..changes to column ordering.
        if (nsrch .le. nn) go to 460
        colupd = nn + 1
        lenpp  = oldend-oldpiv+1
        if (lenpp .lt. 4) lpiv(1) = lpiv(1) + 1
        if (lenpp.ge.4 .and. lenpp.le.6) lpiv(2) = lpiv(2) + 1
        if (lenpp.ge.7 .and. lenpp.le.10) lpiv(3) = lpiv(3) + 1
        if (lenpp.ge.11 .and. lenpp.le.15) lpiv(4) = lpiv(4) + 1
        if (lenpp.ge.16 .and. lenpp.le.20) lpiv(5) = lpiv(5) + 1
        if (lenpp.ge.21 .and. lenpp.le.30) lpiv(6) = lpiv(6) + 1
        if (lenpp.ge.31 .and. lenpp.le.50) lpiv(7) = lpiv(7) + 1
        if (lenpp.ge.51 .and. lenpp.le.70) lpiv(8) = lpiv(8) + 1
        if (lenpp.ge.71 .and. lenpp.le.100) lpiv(9) = lpiv(9) + 1
        if (lenpp.ge.101) lpiv(10) = lpiv(10) + 1
        mapiv  = max0(mapiv,lenpp)
        iavpiv = iavpiv + lenpp
        do 450 jj=oldpiv,oldend
         j        = icn(jj)
         lc       = lastc(j)
         nc       = nextc(j)
         nextc(j) = -colupd
         if (jj .ne. ijpos) colupd = j
         if (nc .ne. 0) lastc(nc) = lc
         if (lc .eq. 0) go to 440
         nextc(lc) = nc
         go to 450
440      nz  = lenc(j)
         isw = ifirst(nz)
         if (isw .gt. 0) lastr(isw) = -nc
         if (isw .lt. 0) ifirst(nz) = -nc
450     continue
!..
!..changes to row ordering.
460     i1 = ipc(jpiv)
        i2 = i1 + lenc(jpiv) - 1
        do 480 ii=i1,i2
         i  = irn(ii)
         lr = lastr(i)
         nr = nextr(i)
         if (nr .ne. 0) lastr(nr) = lr
         if (lr .le. 0) go to 470
         nextr(lr) = nr
         go to 480
470      nz = lenr(i) - lenrl(i)
         if (nr .ne. 0) ifirst(nz) = nr
         if (nr .eq. 0) ifirst(nz) = lr
480     continue
!..
!..move pivot to position lenrl+1 in pivot row and move pivot row to the 
!..beginning of the available storage. the l part and the pivot in the old 
!..copy of the pivot row is nullified while, in the strictly upper triangular 
!..part, the column indices, j say, are overwritten by the corresponding
!..entry of iq (iq(j)) and iq(j) is set to the negative of the displacement of 
!..the column index from the pivot entry.
        if (oldpiv .eq. ijpos) go to 490
        au          = a(oldpiv)
        a(oldpiv)   = a(ijpos)
        a(ijpos)    = au
        icn(ijpos)  = icn(oldpiv)
        icn(oldpiv) = jpiv
!..
!..check if there is space available in a/icn to hold new copy of pivot row.
490     minicn = max0(minicn,nzrow+ibeg-1+morei+lenr(ipiv))
        if (iactiv-ibeg .ge. lenr(ipiv)) go to 500
        call ma30dd(a,icn,iptr(istart),n,iactiv,itop,.true.)
        oldpiv = iptr(ipiv) + lenrl(ipiv)
        oldend = iptr(ipiv) + lenr(ipiv) - 1
!..
!..check now to see if ma30d/dd has created enough available space.
        if (iactiv-ibeg .ge. lenr(ipiv)) go to 500
!..
!..create more space by destroying previously created lu factors.
        morei = morei + ibeg - idisp(1)
        ibeg = idisp(1)
        if (lp .ne. 0) write (lp,99997)
        iflag = -5
        if (abort3) go to 1090
        if (iactiv-ibeg .ge. lenr(ipiv)) go to 500
!..
!..there is still not enough room in a/icn.
        iflag = -4
        go to 1090
!..
!..copy pivot row and set up iq array.
500     ijpos = 0
        j1    = iptr(ipiv)
        do 530 jj=j1,oldend
         a(ibeg)   = a(jj)
         icn(ibeg) = icn(jj)
         if (ijpos .ne. 0) go to 510
         if (icn(jj) .eq. jpiv) ijpos = ibeg
         icn(jj) = 0
         go to 520
510      k       = ibeg - ijpos
         j       = icn(jj)
         icn(jj) = iq(j)
         iq(j)   = -k
520      ibeg    = ibeg + 1
530     continue
!..
        ijp1       = ijpos + 1
        pivend     = ibeg - 1
        lenpiv     = pivend - ijpos
        nzrow      = nzrow - lenrl(ipiv) - 1
        iptr(ipiv) = oldpiv + 1
        if (lenpiv .eq. 0) iptr(ipiv) = 0
!..
!..remove pivot row (including pivot) from column oriented file.
        do 560 jj=ijpos,pivend
         j       = icn(jj)
         i1      = ipc(j)
         lenc(j) = lenc(j) - 1
!..
!..i2 is last position in new column.
         i2 = ipc(j) + lenc(j) - 1
         if (i2 .lt. i1) go to 550
         do 540 ii=i1,i2
          if (irn(ii) .ne. ipiv) go to 540
          irn(ii) = irn(i2+1)
          go to 550
540      continue
550      irn(i2+1) = 0
560     continue
        nzcol = nzcol - lenpiv - 1
!..
!..go down the pivot column and for each row with a non-zero add the 
!..appropriate multiple of the pivot row to it. we loop on the number of 
!..non-zeros in the pivot column since ma30d/dd may change its actual position.
        nzpc = lenc(jpiv)
        if (nzpc .eq. 0) go to 900
        do 840 iii=1,nzpc
         ii = ipc(jpiv) + iii - 1
         i  = irn(ii)
!..
!..search row i for non-zero to be eliminated, calculate multiplier, and place 
!..it in position lenrl+1 in its row. idrop is the number of non-zero entries 
!..dropped from row i because these fall beneath tolerance level.
         idrop = 0
         j1    = iptr(i) + lenrl(i)
         iend  = iptr(i) + lenr(i) - 1
         do 570 jj=j1,iend
          if (icn(jj) .ne. jpiv) go to 570
!..
!..if pivot is zero, rest of column is and so multiplier is zero.
          au = zero
          if (a(ijpos) .ne. zero) au = -a(jj)/a(ijpos)
          if (lbig) big = max(big,abs(au))
          a(jj)    = a(j1)
          a(j1)    = au
          icn(jj)  = icn(j1)
          icn(j1)  = jpiv
          lenrl(i) = lenrl(i) + 1
          go to 580
570      continue
!..
!..jump if pivot row is a singleton.
580      if (lenpiv .eq. 0) go to 840
!..
!..now perform necessary operations on rest of non-pivot row i.
         rowi = j1 + 1
         iop  = 0
!..
!..jump if all the pivot row causes fill-in.
         if (rowi .gt. iend) go to 650
!..
!..perform operations on current non-zeros in row i. innermost loop.
         lenpp = iend-rowi+1
         if (lenpp .lt. 4) lnpiv(1) = lnpiv(1) + 1
         if (lenpp.ge.4 .and. lenpp.le.6) lnpiv(2) = lnpiv(2) + 1
         if (lenpp.ge.7 .and. lenpp.le.10) lnpiv(3) = lnpiv(3) + 1
         if (lenpp.ge.11 .and. lenpp.le.15) lnpiv(4) = lnpiv(4) + 1
         if (lenpp.ge.16 .and. lenpp.le.20) lnpiv(5) = lnpiv(5) + 1
         if (lenpp.ge.21 .and. lenpp.le.30) lnpiv(6) = lnpiv(6) + 1
         if (lenpp.ge.31 .and. lenpp.le.50) lnpiv(7) = lnpiv(7) + 1
         if (lenpp.ge.51 .and. lenpp.le.70) lnpiv(8) = lnpiv(8) + 1
         if (lenpp.ge.71 .and. lenpp.le.100) lnpiv(9) = lnpiv(9) + 1
         if (lenpp.ge.101) lnpiv(10) = lnpiv(10) + 1
         manpiv = max0(manpiv,lenpp)
         ianpiv = ianpiv + lenpp
         kountl = kountl + 1
         do 590 jj=rowi,iend
          j = icn(jj)
          if (iq(j) .gt. 0) go to 590
          iop    = iop + 1
          pivrow = ijpos - iq(j)
          a(jj)  = a(jj) + au*a(pivrow)
          if (lbig) big = max(abs(a(jj)),big)
          icn(pivrow) = -icn(pivrow)
          if (abs(a(jj)) .lt. tol) idrop = idrop + 1
590      continue
!..
!..jump if no non-zeros in non-pivot row have been removed because these are 
!..beneath the drop-tolerance  tol.
         if (idrop .eq. 0) go to 650
!..
!..run through non-pivot row compressing row so that only non-zeros greater 
!..than tol are stored. all non-zeros less than tol are also removed from the 
!..column structure.
         jnew = rowi
         do 630 jj=rowi,iend
          if (abs(a(jj)) .lt. tol) go to 600
          a(jnew)   = a(jj)
          icn(jnew) = icn(jj)
          jnew      = jnew + 1
          go to 630
!..
!..remove non-zero entry from column structure.
600       j = icn(jj)
          i1 = ipc(j)
          i2 = i1 + lenc(j) - 1
          do 610 ii=i1,i2
           if (irn(ii) .eq. i) go to 620
610       continue
620       irn(ii) = irn(i2)
          irn(i2) = 0
          lenc(j) = lenc(j) - 1
          if (nsrch .le. nn) go to 630
!..
!..remove column from column chain and place in update chain.
          if (nextc(j) .lt. 0) go to 630
!..
!..jump if column already in update chain.
          lc       = lastc(j)
          nc       = nextc(j)
          nextc(j) = -colupd
          colupd   = j
          if (nc .ne. 0) lastc(nc) = lc
          if (lc .eq. 0) go to 622
          nextc(lc) = nc
          go to 630
622       nz = lenc(j) + 1
          isw = ifirst(nz)
          if (isw .gt. 0) lastr(isw) = -nc
          if (isw .lt. 0) ifirst(nz) = -nc
630      continue
         do 640 jj=jnew,iend
          icn(jj) = 0
640      continue
!..
!..the value of idrop might be different from that calculated earlier because, 
!..we may have dropped some non-zeros which were not modified by the pivot row.
         idrop   = iend + 1 - jnew
         iend    = jnew - 1
         lenr(i) = lenr(i) - idrop
         nzrow   = nzrow - idrop
         nzcol   = nzcol - idrop
         ndrop   = ndrop + idrop
650      ifill   = lenpiv - iop
!..
!..jump is if there is no fill-in.
         if (ifill .eq. 0) go to 750
!..
!..now for the fill-in.
         minicn = max0(minicn,morei+ibeg-1+nzrow+ifill+lenr(i))
!..
!..see if there is room for fill-in. get maximum space for row i in situ.
         do 660 jdiff=1,ifill
          jnpos = iend + jdiff
          if (jnpos .gt. licn) go to 670
          if (icn(jnpos) .ne. 0) go to 670
660      continue
!..
!..there is room for all the fill-in after the end of the row so it can be 
!..left in situ. next available space for fill-in.
         iend = iend + 1
         go to 750
!..
!..jmore spaces for fill-in are required in front of row.
670      jmore = ifill - jdiff + 1
         i1    = iptr(i)
!..
!..look in front of the row to see if there is space for rest of the fill-in.
         do 680 jdiff=1,jmore
          jnpos = i1 - jdiff
          if (jnpos .lt. iactiv) go to 690
          if (icn(jnpos) .ne. 0) go to 700
680      continue
690      jnpos = i1 - jmore
         go to 710
!..
!..whole row must be moved to the beginning of available storage.
700      jnpos = iactiv - lenr(i) - ifill
!..
!..jump if there is space immediately available for the shifted row.
710      if (jnpos .ge. ibeg) go to 730
         call ma30dd(a,icn,iptr(istart),n,iactiv,itop,.true.)
         i1    = iptr(i)
         iend  = i1 + lenr(i) - 1
         jnpos = iactiv - lenr(i) - ifill
         if (jnpos .ge. ibeg) go to 730
!..
!..no space available; try to create some by trashing previous lu decomposition.
         morei = morei + ibeg - idisp(1) - lenpiv - 1
         if (lp .ne. 0) write (lp,99997)
         iflag = -5
         if (abort3) go to 1090
!..
!..keep record of current pivot row.
         ibeg      = idisp(1)
         icn(ibeg) = jpiv
         a(ibeg)   = a(ijpos)
         ijpos     = ibeg
         do 720 jj=ijp1,pivend
          ibeg      = ibeg + 1
          a(ibeg)   = a(jj)
          icn(ibeg) = icn(jj)
  720    continue
         ijp1   = ijpos + 1
         pivend = ibeg
         ibeg   = ibeg + 1
         if (jnpos .ge. ibeg) go to 730
!..
!..this still does not give enough room.
         iflag = -4
         go to 1090
730      iactiv = min0(iactiv,jnpos)
!..
!..move non-pivot row i.
         iptr(i) = jnpos
         do 740 jj=i1,iend
          a(jnpos)   = a(jj)
          icn(jnpos) = icn(jj)
          jnpos      = jnpos + 1
          icn(jj)    = 0
740      continue
!..
!..first new available space.
         iend  = jnpos
750      nzrow = nzrow + ifill
!..
!..innermost fill-in loop which also resets icn.
         idrop = 0
         do 830 jj=ijp1,pivend
          j = icn(jj)
          if (j .lt. 0) go to 820
          anew  = au*a(jj)
          aanew = abs(anew)
          if (aanew .ge. tol) go to 760
          idrop  = idrop + 1
          ndrop  = ndrop + 1
          nzrow  = nzrow - 1
          minicn = minicn - 1
          ifill  = ifill - 1
          go to 830
760       if (lbig) big = max(aanew,big)
          a(iend)   = anew
          icn(iend) = j
          iend      = iend + 1
!..
!..put new entry in column file.
          minirn = max0(minirn,nzcol+lenc(j)+1)
          jend   = ipc(j) + lenc(j)
          jroom  = nzpc - iii + 1 + lenc(j)
          if (jend .gt. lirn) go to 770
          if (irn(jend) .eq. 0) go to 810
770       if (jroom .lt. dispc) go to 780
!..
!..compress column file to obtain space for new copy of column.
          call ma30dd(a,irn,ipc(istart),n,dispc,lirn,.false.)
          if (jroom .lt. dispc) go to 780
          jroom = dispc - 1
          if (jroom .ge. lenc(j)+1) go to 780
!..
!..column file is not large enough.
          go to 1100
!..
!..copy column to beginning of file.
780       jbeg   = ipc(j)
          jend   = ipc(j) + lenc(j) - 1
          jzero  = dispc - 1
          dispc  = dispc - jroom
          idispc = dispc
          do 790 ii=jbeg,jend
           irn(idispc) = irn(ii)
           irn(ii) = 0
           idispc  = idispc + 1
790       continue
          ipc(j) = dispc
          jend   = idispc
          do 800 ii=jend,jzero
           irn(ii) = 0
800       continue
810       irn(jend) = i
          nzcol     = nzcol + 1
          lenc(j)   = lenc(j) + 1
!..
!..end of adjustment to column file.
          go to 830
!..
820       icn(jj) = -j
830      continue
         if (idrop .eq. 0) go to 834
         do 832 kdrop=1,idrop
          icn(iend) = 0
          iend = iend + 1
832      continue
834      lenr(i) = lenr(i) + ifill
!..
!..end of scan of pivot column.
840     continue
!..
!..
!..remove pivot column from column oriented storage; update row ordering arrays.
        i1 = ipc(jpiv)
        i2 = ipc(jpiv) + lenc(jpiv) - 1
        nzcol = nzcol - lenc(jpiv)
        do 890 ii=i1,i2
         i       = irn(ii)
         irn(ii) = 0
         nz      = lenr(i) - lenrl(i)
         if (nz .ne. 0) go to 850
         lastr(i) = 0
         go to 890
850      ifir       = ifirst(nz)
         ifirst(nz) = i
         if (ifir) 860, 880, 870
860      lastr(i) = ifir
         nextr(i) = 0
         go to 890
870      lastr(i)    = lastr(ifir)
         nextr(i)    = ifir
         lastr(ifir) = i
         go to 890
880      lastr(i) = 0
         nextr(i) = 0
         nzmin    = min0(nzmin,nz)
890     continue
!..
!..restore iq and nullify u part of old pivot row. record the column 
!..permutation in lastc(jpiv) and the row permutation in lastr(ipiv).
900     ipc(jpiv)   = -ising
        lastr(ipiv) = pivot
        if (lenpiv .eq. 0) go to 980
        nzrow      = nzrow - lenpiv
        jval       = ijp1
        jzer       = iptr(ipiv)
        iptr(ipiv) = 0
        do 910 jcount=1,lenpiv
         j         = icn(jval)
         iq(j)     = icn(jzer)
         icn(jzer) = 0
         jval      = jval + 1
         jzer      = jzer + 1
910     continue
!..
!..adjust column ordering arrays.
        if (nsrch .gt. nn) go to 920
        do 916 jj=ijp1,pivend
         j  = icn(jj)
         nz = lenc(j)
         if (nz .ne. 0) go to 914
         ipc(j) = 0
         go to 916
914      nzmin = min0(nzmin,nz)
916     continue
        go to 980
920     jj = colupd
        do 970 jdummy=1,nn
         j = jj
         if (j .eq. nn+1) go to 980
         jj = -nextc(j)
         nz = lenc(j)
         if (nz .ne. 0) go to 924
         ipc(j) = 0
         go to 970
924      ifir     = ifirst(nz)
         lastc(j) = 0
         if (ifir) 930, 940, 950
930      ifirst(nz)  = -j
         ifir        = -ifir
         lastc(ifir) = j
         nextc(j)    = ifir
         go to 970
940      ifirst(nz) = -j
         nextc(j)   = 0
         go to 960
950      lc          = -lastr(ifir)
         lastr(ifir) = -j
         nextc(j)    = lc
         if (lc .ne. 0) lastc(lc) = j
960      nzmin = min0(nzmin,nz)
970     continue
980    continue
!..
!..that was the end of main elimination loop
!..
!..
!..reset iactiv to point to the beginning of the next block.
990    if (ilast .ne. nn) iactiv = iptr(ilast+1)
1000  continue
!..
!..that was the end of deomposition of block   
!..
!..
!..record singularity (if any) in iq array.
      if (irank .eq. nn) go to 1020
      do 1010 i=1,nn
       if (ipc(i) .lt. 0) go to 1010
       ising     = ipc(i)
       iq(ising) = -iq(ising)
       ipc(i)    = -ising
1010  continue
!..
!..
!..run through lu decomposition changing column indices to that of new order 
!..and permuting lenr and lenrl arrays according to pivot permutations.
1020  istart = idisp(1)
      iend = ibeg - 1
      if (iend .lt. istart) go to 1040
      do 1030 jj=istart,iend
       jold    = icn(jj)
       icn(jj) = -ipc(jold)
1030  continue
1040  do 1050 ii=1,nn
       i        = lastr(ii)
       nextr(i) = lenr(ii)
       iptr(i)  = lenrl(ii)
1050  continue
      do 1060 i=1,nn
       lenrl(i) = iptr(i)
       lenr(i)  = nextr(i)
1060  continue
!..
!..update permutation arrays ip and iq.
      do 1070 ii=1,nn
       i        = lastr(ii)
       j        = -ipc(ii)
       nextr(i) = iabs(ip(ii)+0)
       iptr(j)  = iabs(iq(ii)+0)
1070  continue
      do 1080 i=1,nn
       if (ip(i) .lt. 0) nextr(i) = -nextr(i)
       ip(i) = nextr(i)
       if (iq(i) .lt. 0) iptr(i) = -iptr(i)
       iq(i) = iptr(i)
1080  continue
      ip(nn)   = iabs(ip(nn)+0)
      idisp(2) = iend
      go to 1120
!..
!..
!..error returns
1090  idisp(2) = iactiv
      if (lp .eq. 0) go to 1120
      write (lp,99996)
      go to 1110
1100  if (iflag .eq. -5) iflag = -6
      if (iflag .ne. -6) iflag = -3
      idisp(2) = iactiv
      if (lp .eq. 0) go to 1120
      if (iflag .eq. -3) write (lp,99995)
      if (iflag .eq. -6) write (lp,99994)
1110  pivot = pivot - istart + 1
      write (lp,99993) pivot,nblock,istart,ilast
      if (pivot .eq. 0) write (lp,99992) minirn
1120  return
      end
!..
!..
!..
!..
!..
      subroutine ma30bd(n,icn,a,licn,lenr,lenrl,idisp,ip,iq,w,iw,iflag)
      implicit none
      save 
!..
!..ma30b/bd performs the lu decomposition of the diagonal blocks of a new 
!..matrix paq of the same sparsity pattern, using information from a previous 
!..call to ma30a/ad. the entries of the input matrix  must already be in their 
!..final positions in the lu decomposition structure.  this routine executes 
!..about five times faster than ma30a/ad.
!..
!..parameters (see also ma30ad):
!..n   is an integer variable set to the order of the matrix.
!..
!..icn is an integer array of length licn. it should be unchanged
!..    since the last call to ma30a/ad. it is not altered by ma30b/bd.
!..
!..a   is a real/real array of length licn the user must set
!..    entries idisp(1) to idisp(2) to contain the entries in the
!..    diagonal blocks of the matrix paq whose column numbers are held
!..    in icn, using corresponding positions. note that some zeros may
!..    need to be held explicitly. on output entries idisp(1) to
!..    idisp(2) of array a contain the lu decomposition of the diagonal
!..    blocks of paq. entries a(1) to a(idisp(1)-1) are neither
!..    required nor altered by ma30b/bd.
!..
!..licn is an integer variable which must be set by the user to the
!..     length of arrays a and icn. it is not altered by ma30b/bd.
!..
!..lenr,lenrl are integer arrays of length n. they should be
!..    unchanged since the last call to ma30a/ad. not altered by ma30b/bd.
!..
!..idisp is an integer array of length 2. it should be unchanged since
!..      the last call to ma30a/ad. it is not altered by ma30b/bd.
!..
!..ip,iq are integer arrays of length n. they should be unchanged
!..      since the last call to ma30a/ad. not altered by ma30b/bd.
!..
!..w   is a array of length n which is used as workspace by ma30b/bd.
!..
!..iw  is an integer array of length n which is used as workspace by ma30b/bd.
!..
!..iflag  is an integer variable. on output from ma30b/bd, iflag has
!..       the value zero if the factorization was successful, has the
!..       value i if pivot i was very small and has the value -i if an
!..       unexpected singularity was detected at stage i of the decomposition.
!..
!..declare
      logical          abort1,abort2,abort3,stab,lbig
      integer          n,licn,iflag,iw(n),idisp(2),pivpos,icn(licn), & 
     &                 lenr(n),lenrl(n),ip(n),iq(n),ndrop,nsrch,ising, & 
     &                 i,istart,ifin,ilend,j,ipivj,jfin,jay,jayjay, & 
     &                 jj,lp
      real             a(licn),w(n),au,eps,rowmax,zero,one,rmin,tol,big
!..
      common /ma30ed/  lp,abort1,abort2,abort3
      common /ma30id/  tol,big,ndrop,nsrch,lbig
      common /ma30gd/  eps,rmin
      data             zero /0.0e0/, one /1.0e0/
!..
!..formats
99999 format(1x,'error return from ma30b/bd singularity in row',i8)
!..
!..initialize
      stab  = eps .le. one
      rmin  = eps
      ising = 0
      iflag = 0
      do 10 i=1,n
       w(i) = zero
10    continue
!..
!..set up pointers to the beginning of the rows.
      iw(1) = idisp(1)
      if (n .eq. 1) go to 25
      do 20 i=2,n
       iw(i) = iw(i-1) + lenr(i-1)
20    continue
!..
!..start  of main loop
!..at step i, row i of a is transformed to row i of l/u by adding appropriate 
!..multiples of rows 1 to i-1. using row-gauss elimination.
!..istart is beginning of row i of a and row i of l.
!..ifin is end of row i of a and row i of u.
!..ilend is end of row i of l.
25    do 160 i=1,n
       istart = iw(i)
       ifin   = istart + lenr(i) - 1
       ilend  = istart + lenrl(i) - 1
       if (istart .gt. ilend) go to 90
!..
!..load row i of a into vector w.
       do 30 jj=istart,ifin
        j    = icn(jj)
        w(j) = a(jj)
30     continue
!..
!..add multiples of appropriate rows of  i to i-1  to row i.
!..ipivj is position of pivot in row j. 
       do 70 jj=istart,ilend
        j     = icn(jj)
        ipivj = iw(j) + lenrl(j)
        au    = -w(j)/a(ipivj)
        if (lbig) big = max(abs(au),big)
        w(j) = au
!..
!..au * row j (u part) is added to row i.
        ipivj = ipivj + 1
        jfin  = iw(j) + lenr(j) - 1
        if (ipivj .gt. jfin) go to 70
!..
!..innermost loop.
        if (lbig) go to 50
        do 40 jayjay=ipivj,jfin
         jay    = icn(jayjay)
         w(jay) = w(jay) + au*a(jayjay)
40      continue
        go to 70
50      do 60 jayjay=ipivj,jfin
         jay    = icn(jayjay)
         w(jay) = w(jay) + au*a(jayjay)
         big    = max(abs(w(jay)),big)
60      continue
70     continue
!..
!..reload w back into a (now l/u)
       do 80 jj=istart,ifin
        j     = icn(jj)
        a(jj) = w(j)
        w(j)  = zero
80     continue
!..
!..now perform the stability checks.
90     pivpos = ilend + 1
       if (iq(i) .gt. 0) go to 140
!..
!..matrix had singularity at this point in ma30a/ad.
!..is it the first such pivot in current block ?
       if (ising .eq. 0) ising = i
!..
!..does current matrix have a singularity in the same place ?
       if (pivpos .gt. ifin) go to 100
       if (a(pivpos) .ne. zero) go to 170
!..
!..it does .. so set ising if it is not the end of the current block
!..check to see that appropriate part of l/u is zero or null.
100    if (istart .gt. ifin) go to 120
       do 110 jj=istart,ifin
        if (icn(jj) .lt. ising) go to 110
        if (a(jj) .ne. zero) go to 170
110    continue
120    if (pivpos .le. ifin) a(pivpos) = one
       if (ip(i) .gt. 0 .and. i .ne. n) go to 160
!..
!..end of current block ... reset zero pivots and ising.
       do 130 j=ising,i
        if ((lenr(j)-lenrl(j)) .eq. 0) go to 130
        jj    = iw(j) + lenrl(j)
        a(jj) = zero
130    continue
       ising = 0
       go to 160
!..
!..matrix had non-zero pivot in ma30a/ad at this stage.
140    if (pivpos .gt. ifin) go to 170
       if (a(pivpos) .eq. zero) go to 170
       if (.not.stab) go to 160
       rowmax = zero
       do 150 jj=pivpos,ifin
        rowmax = max(rowmax,abs(a(jj)))
150    continue
       if (abs(a(pivpos))/rowmax .ge. rmin) go to 160
       iflag = i
       rmin  = abs(a(pivpos))/rowmax
160   continue
      go to 180
!..
!..error return
170   if (lp .ne. 0) write (lp,99999) i
      iflag = -i
180   return
      end
!..
!..
!..
!..
!..
      subroutine ma30cd(n,icn,a,licn,lenr,lenrl,lenoff,idisp,ip, & 
     &                  iq,x,w,mtype)
      implicit none
      save 
!..
!..
!..ma30c/cd uses the factors produced by ma30a/ad or ma30b/bd to solve
!..ax=b or a transpose x=b when the matrix p1*a*q1 (paq) is block lower 
!..triangular (including the case of only one diagonal block).
!..
!..parameters: 
!..n  is an integer variable set to the order of the matrix. it is not
!..   altered by the subroutine.
!..
!..icn is an integer array of length licn. entries idisp(1) to
!..    idisp(2) should be unchanged since the last call to ma30a/ad. if
!..    the matrix has more than one diagonal block, then column indices
!..    corresponding to non-zeros in sub-diagonal blocks of paq must
!..    appear in positions 1 to idisp(1)-1. for the same row those
!..    entries must be contiguous, with those in row i preceding those
!..    in row i+1 (i=1,...,n-1) and no wasted space between rows.
!..    entries may be in any order within each row. not altered by ma30c/cd.
!..
!..a  is a real/real array of length licn.  entries
!..   idisp(1) to idisp(2) should be unchanged since the last call to
!..   ma30a/ad or ma30b/bd.  if the matrix has more than one diagonal
!..   block, then the values of the non-zeros in sub-diagonal blocks
!..   must be in positions 1 to idisp(1)-1 in the order given by icn.
!..   it is not altered by ma30c/cd.
!..
!..licn  is an integer variable set to the size of arrays icn and a.
!..      it is not altered by ma30c/cd.
!..
!..lenr,lenrl are integer arrays of length n which should be
!..     unchanged since the last call to ma30a/ad. not altered by ma30c/cd.
!..
!..lenoff  is an integer array of length n. if the matrix paq (or
!..        p1*a*q1) has more than one diagonal block, then lenoff(i),
!..        i=1,...,n should be set to the number of non-zeros in row i of
!..        the matrix paq which are in sub-diagonal blocks.  if there is
!..        only one diagonal block then lenoff(1) may be set to -1, in
!..        which case the other entries of lenoff are never accessed. it is
!..        not altered by ma30c/cd.
!..
!..idisp  is an integer array of length 2 which should be unchanged
!..       since the last call to ma30a/ad. it is not altered by ma30c/cd.
!..
!c..ip,iq are integer arrays of length n which should be unchanged
!..       since the last call to ma30a/ad. they are not altered by ma30c/cd.
!..
!..x   is a real/real array of length n. it must be set by
!..    the user to the values of the right hand side vector b for the
!..    equations being solved.  on exit from ma30c/cd it will be equal
!..    to the solution x required.
!..
!..w  is a real/real array of length n which is used as
!..   workspace by ma30c/cd.
!..
! mtype is an integer variable which must be set by the user. if
!..     mtype=1, then the solution to the system ax=b is returned; any
!..     other value for mtype will return the solution to the system a
!..     transpose x=b. it is not altered by ma30c/cd.
!..
!..declare
      logical          neg,nobloc
      integer          n,licn,idisp(2),icn(licn),lenr(n),lenrl(n), & 
     &                 lenoff(n),ip(n),iq(n),mtype,ii,i,lt,ifirst, & 
     &                 iblock,ltend,jj,j,iend,j1,ib,iii,j2,jpiv, & 
     &                 jpivp1,ilast,iblend,numblk,k,j3,iback,lj2,lj1
      real             a(licn),x(n),w(n),wii,wi,resid,zero
      common /ma30hd/  resid
      data zero       /0.0e0/
!..
!..final value of resid is the max residual for inconsistent set of equations.
      resid = zero
!..
!..nobloc is .true. if subroutine block has been used previously and is .false. 
!..otherwise.  the value .false. means that lenoff will not be subsequently 
!..accessed.
      nobloc = lenoff(1) .lt. 0
      if (mtype .ne. 1) go to 140
!..
!..now solve   a * x = b. neg is used to indicate when the last row in a block 
!..has been reached.  it is then set to true whereafter backsubstitution is
!..performed on the block.
      neg = .false.
!..
!..ip(n) is negated so that the last row of the last block can be recognised.  
!..it is reset to its positive value on exit.
      ip(n) = -ip(n)
!..
!..preorder vector ... w(i) = x(ip(i))
      do 10 ii=1,n
       i     = ip(ii)
       i     = iabs(i)
       w(ii) = x(i)
10    continue
!..
!..lt is the position of first non-zero in current row of off-diagonal blocks.
!..ifirst holds the index of the first row in the current block.
!..iblock holds the position of the first non-zero in the current row
!..of the lu decomposition of the diagonal blocks.
      lt     = 1
      ifirst = 1
      iblock = idisp(1)
!..
!..if i is not the last row of a block, then a pass through this loop adds the 
!..inner product of row i of the off-diagonal blocks and w to w and performs 
!..forward elimination using row i of the lu decomposition.   if i is the last 
!..row of a block then, after performing these aforementioned operations, 
!..backsubstitution is performed using the rows of the block.
      do 120 i=1,n
       wi = w(i)
       if (nobloc) go to 30
       if (lenoff(i) .eq. 0) go to 30
!..
!..operations using lower triangular blocks. 
!..ltend is the end of row i in the off-diagonal blocks.
       ltend = lt + lenoff(i) - 1
       do 20 jj=lt,ltend
        j  = icn(jj)
        wi = wi - a(jj)*w(j)
20     continue
!..
!..lt is set the beginning of the next off-diagonal row.
!..set neg to .true. if we are on the last row of the block.
       lt = ltend + 1
30     if (ip(i) .lt. 0) neg = .true.
       if (lenrl(i) .eq. 0) go to 50
!..
!..forward elimination phase.
!..iend is the end of the l part of row i in the lu decomposition.
       iend = iblock + lenrl(i) - 1
       do 40 jj=iblock,iend
        j  = icn(jj)
        wi = wi + a(jj)*w(j)
40     continue
!..
!..iblock is adjusted to point to the start of the next row.
50     iblock = iblock + lenr(i)
       w(i)   = wi
       if (.not.neg) go to 120
!..
!..back substitution phase.
!..j1 is position in a/icn after end of block beginning in row ifirst
!..and ending in row i.
       j1 = iblock
!..
!..are there any singularities in this block?  if not, continue
       ib = i
       if (iq(i) .gt. 0) go to 70
       do 60 iii=ifirst,i
        ib = i - iii + ifirst
        if (iq(ib) .gt. 0) go to 70
        j1    = j1 - lenr(ib)
        resid = max(resid,abs(w(ib)))
        w(ib) = zero
60     continue
!..
!..entire block is singular.
       go to 110
!..
!..
!..each pass through this loop performs the back-substitution
!..operations for a single row, starting at the end of the block and
!..working through it in reverse order.
!..j2 is end of row ii. j1 is beginning of row ii. jpiv is the position of the 
!..pivot in row ii. jump out if row ii of u has no non-zeros.
70     do 100 iii=ifirst,ib
        ii     = ib - iii + ifirst
        j2     = j1 - 1
        j1     = j1 - lenr(ii)
        jpiv   = j1 + lenrl(ii)
        jpivp1 = jpiv + 1
        if (j2 .lt. jpivp1) go to 90
        wii = w(ii)
        do 80 jj=jpivp1,j2
         j   = icn(jj)
         wii = wii - a(jj)*w(j)
80      continue
        w(ii) = wii
90      w(ii) = w(ii)/a(jpiv)
100    continue
110    ifirst = i + 1
       neg    = .false.
120   continue
!..
!..reorder solution vector ... x(i) = w(iqinverse(i))
      do 130 ii=1,n
       i    = iq(ii)
       i    = iabs(i)
       x(i) = w(ii)
130   continue
      ip(n) = -ip(n)
      go to 320
!..
!..
!..now solve  atranspose * x = b. preorder vector ... w(i)=x(iq(i))
140   do 150 ii=1,n
       i     = iq(ii)
       i     = iabs(i)
       w(ii) = x(i)
150   continue
!..
!..lj1 points to the beginning the current row in the off-diagonal blocks.
!..iblock is initialized to point to beginning of block after the last one
!..ilast is the last row in the current block.
!..iblend points to the position after the last non-zero in the current block.
      lj1    = idisp(1)
      iblock = idisp(2) + 1
      ilast  = n
      iblend = iblock
!..
!..each pass through this loop operates with one diagonal block and
!..the off-diagonal part of the matrix corresponding to the rows
!..of this block.  the blocks are taken in reverse order and the
!..number of times the loop is entered is min(n,no. blocks+1).
      do 290 numblk=1,n
       if (ilast .eq. 0) go to 300
       iblock = iblock - lenr(ilast)
!..
!..this loop finds the index of the first row in the current block. it is 
!..first and iblock is set to the position of the beginning of this first row.
       do 160 k=1,n
        ii = ilast - k
        if (ii .eq. 0) go to 170
        if (ip(ii) .lt. 0) go to 170
        iblock = iblock - lenr(ii)
160    continue
170    ifirst = ii + 1
!..
!..j1 points to the position of the beginning of row i (lt part) or pivot
       j1 = iblock
!..
!..forward elimination. each pass through this loop performs the operations 
!..for one row of the block.  if the corresponding entry of w is zero then the
!..operations can be avoided.
       do 210 i=ifirst,ilast
        if (w(i) .eq. zero) go to 200
!..
!.. jump if row i singular.
        if (iq(i) .lt. 0) go to 220
!..
!..j2 first points to the pivot in row i and then is made to point to the
!..first non-zero in the u transpose part of the row.
!..j3 points to the end of row i.
        j2 = j1 + lenrl(i)
        wi = w(i)/a(j2)
        if (lenr(i)-lenrl(i) .eq. 1) go to 190
        j2 = j2 + 1
        j3 = j1 + lenr(i) - 1
        do 180 jj=j2,j3
         j    = icn(jj)
         w(j) = w(j) - a(jj)*wi
180     continue
190     w(i) = wi
200     j1   = j1 + lenr(i)
210    continue
       go to 240
!..
!..deals with rest of block which is singular.
220    do 230 ii=i,ilast
        resid = max(resid,abs(w(ii)))
        w(ii) = zero
230    continue
!..
!..back substitution. this loop does the back substitution on the rows of the 
!..block in the reverse order doing it simultaneously on the l transpose part
!..of the diagonal blocks and the off-diagonal blocks.
!..j1 points to the beginning of row i.
!..j2 points to the end of the l transpose part of row i.
240    j1 = iblend
       do 280 iback=ifirst,ilast
        i  = ilast - iback + ifirst
        j1 = j1 - lenr(i)
        if (lenrl(i) .eq. 0) go to 260
        j2 = j1 + lenrl(i) - 1
        do 250 jj=j1,j2
         j    = icn(jj)
         w(j) = w(j) + a(jj)*w(i)
250     continue
260     if (nobloc) go to 280
!..
!..operations using lower triangular blocks.
!..lj2 points to the end of row i of the off-diagonal blocks.
!..lj1 points to the beginning of row i of the off-diagonal blocks.
        if (lenoff(i) .eq. 0) go to 280
        lj2 = lj1 - 1
        lj1 = lj1 - lenoff(i)
        do 270 jj=lj1,lj2
         j = icn(jj)
         w(j) = w(j) - a(jj)*w(i)
270     continue
280    continue
       iblend = j1
       ilast = ifirst - 1
290   continue
!..
!..reorder solution vector ... x(i)=w(ipinverse(i))
300   do 310 ii=1,n
       i    = ip(ii)
       i    = iabs(i)
       x(i) = w(ii)
310   continue
320   return
      end
!..
!..
!..
!..
!..
!..
      subroutine ma30dd(a,icn,iptr,n,iactiv,itop,reals)
      implicit none
      save 
!..
!..this subroutine performs garbage collection operations on the arrays a, 
!..icn and irn. iactiv is the first position in arrays a/icn from which the 
!..compress starts.  on exit, iactiv equals the position of the first entry
!..in the compressed part of a/icn
!..
      logical          reals
      integer          n,itop,iptr(n),icn(itop), & 
     &                 irncp,icncp,irank,minirn,minicn,j,k,kn,kl, & 
     &                 jpos,iactiv
      real             a(itop)
      common /ma30fd/  irncp,icncp,irank,minirn,minicn
!..
!..
      if (reals) icncp = icncp + 1
      if (.not.reals) irncp = irncp + 1
!..
!..set the first non-zero entry in each row to the negative of the
!..row/col number and hold this row/col index in the row/col
!..pointer.  this is so that the beginning of each row/col can
!..be recognized in the subsequent scan.
      do 10 j=1,n
       k = iptr(j)
       if (k .lt. iactiv) go to 10
       iptr(j) = icn(k)
       icn(k) = -j
10    continue
      kn = itop + 1
      kl = itop - iactiv + 1
!..
!..go through arrays in reverse order compressing to the back so
!..that there are no zeros held in positions iactiv to itop in icn.
!..reset first entry of each row/col and pointer array iptr.
      do 30 k=1,kl
       jpos = itop - k + 1
       if (icn(jpos) .eq. 0) go to 30
       kn = kn - 1
       if (reals) a(kn) = a(jpos)
       if (icn(jpos) .ge. 0) go to 20
!..
!..first non-zero of row/col has been located
       j         = -icn(jpos)
       icn(jpos) = iptr(j)
       iptr(j)   = kn
20     icn(kn)   = icn(jpos)
30    continue
      iactiv = kn
      return
      end
!..
!..
!..
!..
!..
      subroutine ma28int1
      implicit none
      save 
!..
!..  
!..lp,mp are used by the subroutine as the unit numbers for its warning
!..      and diagnostic messages. default value for both is 6 (for line
!..      printer output). the user can either reset them to a different
!..      stream number or suppress the output by setting them to zero.
!..      while lp directs the output of error diagnostics from the
!..      principal subroutines and internally called subroutines, mp
!..      controls only the output of a message which warns the user that he
!..      has input two or more non-zeros a(i), . . ,a(k) with the same row
!..      and column indices.  the action taken in this case is to proceed
!..      using a numerical value of a(i)+...+a(k). in the absence of other
!..      errors, iflag will equal -14 on exit.
!..lblock is a logical variable which controls an option of first
!..       preordering the matrix to block lower triangular form (using
!..       harwell subroutine mc23a). the preordering is performed if lblock
!..       is equal to its default value of .true. if lblock is set to
!..       .false. , the option is not invoked and the space allocated to
!..       ikeep can be reduced to 4*n+1.
!..grow is a logical variable. if it is left at its default value of
!..     .true. , then on return from ma28a/ad or ma28b/bd, w(1) will give
!..     an estimate (an upper bound) of the increase in size of elements
!..     encountered during the decomposition. if the matrix is well
!..     scaled, then a high value for w(1), relative to the largest entry
!..     in the input matrix, indicates that the lu decomposition may be
!..     inaccurate and the user should be wary of his results and perhaps
!..     increase u for subsequent runs.  we would like to emphasise that
!..     this value only relates to the accuracy of our lu decomposition
!..     and gives no indication as to the singularity of the matrix or the
!..     accuracy of the solution.  this upper bound can be a significant
!..     overestimate particularly if the matrix is badly scaled. if an
!..     accurate value for the growth is required, lbig (q.v.) should be
!..     set to .true.
!..eps,rmin are real variables. if, on entry to ma28b/bd, eps is less
!..     than one, then rmin will give the smallest ratio of the pivot to
!..     the largest element in the corresponding row of the upper
!..     triangular factor thus monitoring the stability of successive
!..     factorizations. if rmin becomes very large and w(1) from
!..     ma28b/bd is also very large, it may be advisable to perform a
!..     new decomposition using ma28a/ad.
!..resid is a real variable which on exit from ma28c/cd gives the value
!..      of the maximum residual over all the equations unsatisfied because
!..      of dependency (zero pivots).
!..irncp,icncp are integer variables which monitor the adequacy of "elbow
!..     room" in irn and a/icn respectively. if either is quite large (say
!..     greater than n/10), it will probably pay to increase the size of
!..     the corresponding array for subsequent runs. if either is very low
!..     or zero then one can perhaps save storage by reducing the size of
!..     the corresponding array.
!..minirn,minicn are integer variables which, in the event of a
!..     successful return (iflag ge 0 or iflag=-14) give the minimum size
!..     of irn and a/icn respectively which would enable a successful run
!..     on an identical matrix. on an exit with iflag equal to -5, minicn
!..     gives the minimum value of icn for success on subsequent runs on
!..     an identical matrix. in the event of failure with iflag= -6, -4,
!..     -3, -2, or -1, then minicn and minirn give the minimum value of
!..     licn and lirn respectively which would be required for a
!..     successful decomposition up to the point at which the failure occurred.
!..irank is an integer variable which gives an upper bound on the rank of
!..      the matrix.
!..abort1 is a logical variable with default value .true.  if abort1 is
!..       set to .false.  then ma28a/ad will decompose structurally singular
!..       matrices (including rectangular ones).
!..abort2 is a logical variable with default value .true.  if abort2 is
!..       set to .false. then ma28a/ad will decompose numerically singular
!..       matrices.
!..idisp is an integer array of length 2. on output from ma28a/ad, the
!..      indices of the diagonal blocks of the factors lie in positions
!..      idisp(1) to idisp(2) of a/icn. this array must be preserved
!..      between a call to ma28a/ad and subsequent calls to ma28b/bd,
!..      ma28c/cd or ma28i/id.
!..tol is a real variable.  if it is set to a positive value, then any
!..    non-zero whose modulus is less than tol will be dropped from the
!..    factorization.  the factorization will then require less storage
!..    but will be inaccurate.  after a run of ma28a/ad with tol positive
!..    it is not possible to use ma28b/bd and the user is recommended to
!..    use ma28i/id to obtain the solution.  the default value for tol is 0.0.
!..themax is a real variable.  on exit from ma28a/ad, it will hold the
!..       largest entry of the original matrix.
!..big is a real variable. if lbig has been set to .true., big will hold
!..    the largest entry encountered during the factorization by ma28a/ad
!..     or ma28b/bd.
!..dxmax is a real variable. on exit from ma28i/id, dxmax will be set to
!..      the largest component of the solution.
!..errmax is a real variable.  on exit from ma28i/id, if maxit is
!..       positive, errmax will be set to the largest component in the
!..       estimate of the error.
!..dres is a real variable.  on exit from ma28i/id, if maxit is positive,
!..     dres will be set to the largest component of the residual.
!..cgce is a real variable. it is used by ma28i/id to check the
!..     convergence rate.  if the ratio of successive corrections is
!..     not less than cgce then we terminate since the convergence
!..     rate is adjudged too slow.
!..ndrop is an integer variable. if tol has been set positive, on exit
!..     from ma28a/ad, ndrop will hold the number of entries dropped from
!..     the data structure.
!..maxit is an integer variable. it is the maximum number of iterations
!..     performed by ma28i/id. it has a default value of 16.
!..noiter is an integer variable. it is set by ma28i/id to the number of
!..     iterative refinement iterations actually used.
!..nsrch is an integer variable. if nsrch is set to a value less than n,
!..     then a different pivot option will be employed by ma28a/ad.  this
!..     may result in different fill-in and execution time for ma28a/ad.
!..     if nsrch is less than or equal to n, the workspace array iw can be
!..     reduced in length.  the default value for nsrch is 32768.
!..istart is an integer variable. if istart is set to a value other than
!..     zero, then the user must supply an estimate of the solution to
!..     ma28i/id.  the default value for istart is zero.
!..lbig is a logical variable. if lbig is set to .true., the value of the
!..    largest element encountered in the factorization by ma28a/ad or
!..    ma28b/bd is returned in big.  setting lbig to .true.  will
!..    increase the time for ma28a/ad marginally and that for ma28b/bd
!..    by about 20%.  the default value for lbig is .false.
!..
!..declare
      logical          lblock,grow,abort1,abort2,lbig
      integer          lp,mp,irncp,icncp,minirn,minicn,irank, & 
     &                 ndrop,maxit,noiter,nsrch,istart
      real             eps,rmin,resid,tol,themax,big,dxmax, & 
     &                 errmax,dres,cgce
      common /ma28ed/  lp,mp,lblock,grow
      common /ma28fd/  eps,rmin,resid,irncp,icncp,minirn,minicn, & 
     &                 irank,abort1,abort2
      common /ma28hd/  tol,themax,big,dxmax,errmax,dres,cgce, & 
     &                 ndrop,maxit,noiter,nsrch,istart,lbig
!..
      eps    = 1.0e-4
      tol    = 0.0e0
      cgce   = 0.5e0
      maxit  = 16
      lp     = 6
      mp     = 6
      nsrch  = 32768
      istart = 0
      lblock = .true.
      grow   = .false.
      lbig   = .false.
      abort1 = .true.
      abort2 = .true.
      return
      end
!..
!..
!..
!..
!..
      subroutine ma28int2
      implicit none
      save 
!..
!..
!..common block ma30e/ed holds control parameters
!..    common /ma30ed/ lp, abort1, abort2, abort3
!..the integer lp is the unit number to which the error messages are
!..sent. lp has a default value of 6.  this default value can be
!..reset by the user, if desired.  a value of 0 suppresses all
!..messages.
!..the logical variables abort1,abort2,abort3 are used to control the
!..conditions under which the subroutine will terminate.
!..if abort1 is .true. then the subroutine will exit  immediately on
!..detecting structural singularity.
!..if abort2 is .true. then the subroutine will exit immediately on
!..detecting numerical singularity.
!..if abort3 is .true. then the subroutine will exit immediately when
!..the available space in a/icn is filled up by the previously decomposed, 
!..active, and undecomposed parts of the matrix.
!..
!..the default values for abort1,abort2,abort3 are set to .true.,.true.
!..and .false. respectively.
!..
!..
!..the variables in the common block ma30f/fd are used to provide the
!..user with information on the decomposition.
!..common /ma30fd/ irncp, icncp, irank, minirn, minicn
!..
!..irncp and icncp are integer variables used to monitor the adequacy
!..of the allocated space in arrays irn and a/icn respectively, by
!..taking account of the number of data management compresses
!..required on these arrays. if irncp or icncp is fairly large (say
!..greater than n/10), it may be advantageous to increase the size
!..of the corresponding array(s).  irncp and icncp are initialized
!..to zero on entry to ma30a/ad and are incremented each time the
!..compressing routine ma30d/dd is entered.
!..
!..icncp is the number of compresses on a/icn.
!..irncp is the number of compresses on irn.
!..
!..irank is an integer variable which gives an estimate (actually an
!..upper bound) of the rank of the matrix. on an exit with iflag
!..equal to 0, this will be equal to n.
!..
! minirn is an integer variable which, after a successful call to
!..ma30a/ad, indicates the minimum length to which irn can be
!..reduced while still permitting a successful decomposition of the
!..same matrix. if, however, the user were to decrease the length
!..of irn to that size, the number of compresses (irncp) may be
!..very high and quite costly. if lirn is not large enough to begin
!..the decomposition on a diagonal block, minirn will be equal to
!..the value required to continue the decomposition and iflag will
!..be set to -3 or -6. a value of lirn slightly greater than this
!..(say about n/2) will usually provide enough space to complete
!..the decomposition on that block. in the event of any other
!..failure minirn gives the minimum size of irn required for a
!..successful decomposition up to that point.
!..
!..minicn is an integer variable which after a successful call to
!..ma30a/ad, indicates the minimum size of licn required to enable
!..a successful decomposition. in the event of failure with iflag=
!..-5, minicn will, if abort3 is left set to .false., indicate the
!..minimum length that would be sufficient to prevent this error in
!..a subsequent run on an identical matrix. again the user may
!..prefer to use a value of icn slightly greater than minicn for
!..subsequent runs to avoid too many conpresses (icncp). in the
!..event of failure with iflag equal to any negative value except
!..-4, minicn will give the minimum length to which licn could be
!..reduced to enable a successful decomposition to the point at
!..which failure occurred.  notice that, on a successful entry
!..idisp(2) gives the amount of space in a/icn required for the
!..decomposition while minicn will usually be slightly greater
!..because of the need for "elbow room".  if the user is very
!..unsure how large to make licn, the variable minicn can be used
!..to provide that information. a preliminary run should be
!..performed with abort3 left set to .false. and licn about 3/2
!..times as big as the number of non-zeros in the original matrix.
!..unless the initial problem is very sparse (when the run will be
!..successful) or fills in extremely badly (giving an error return
!..with iflag equal to -4), an error return with iflag equal to -5
!..should result and minicn will give the amount of space required
!..for a successful decomposition.
!..
!..
!..common block ma30g/gd is used by the ma30b/bd entry only.
!..   common /ma30gd/ eps, rmin
! eps is a real/real variable. it is used to test for
!..small pivots. its default value is 1.0e-4 (1.0e-4 in d version).
!..if the user sets eps to any value greater than 1.0, then no
!..check is made on the size of the pivots. although the absence of
!..such a check would fail to warn the user of bad instability, its
!..absence will enable ma30b/bd to run slightly faster. an  a
!..posteriori  check on the stability of the factorization can be
!..obtained from mc24a/ad.
!..
!..rmin is a real/real variable which gives the user some
!..information about the stability of the decomposition.  at each
!..stage of the lu decomposition the magnitude of the pivot apiv
!..is compared with the largest off-diagonal entry currently in its
!..row (row of u), rowmax say. if the ratio min (apiv/rowmax)
!..where the minimum is taken over all the rows, is less than eps
!..then rmin is set to this minimum value and iflag is returned
!..with the value +i where i is the row in which this minimum
!..occurs.  if the user sets eps greater than one, then this test
!..is not performed. in this case, and when there are no small
!..pivots rmin will be set equal to eps.
!..
!..
!..common block ma30h/hd is used by ma30c/cd only.
!..   common /ma30hd/ resid
!..resid is a real/real variable. in the case of singular
!..or rectangular matrices its final value will be equal to the
!..maximum residual for the unsatisfied equations; otherwise its
!..value will be set to zero.
!..
!..
!..common  block ma30i/id controls the use of drop tolerances, the
!..modified pivot option and the the calculation of the largest
!..entry in the factorization process. this common block was added
!..to the ma30 package in february, 1983.
!..   common /ma30id/ tol, big, ndrop, nsrch, lbig
!..
!..tol is a real/real variable.  if it is set to a positive
!..value, then ma30a/ad will drop from the factors any non-zero
!..whose modulus is less than tol.  the factorization will then
!..require less storage but will be inaccurate.  after a run of
!..ma30a/ad where entries have been dropped, ma30b/bd  should not
!..be called.  the default value for tol is 0.0.
!..
!..big is a real/real variable.  if lbig has been set to
!...true., big will be set to the largest entry encountered during
!..the factorization.
!..ndrop is an integer variable. if tol has been set positive, on exit
!..from ma30a/ad, ndrop will hold the number of entries dropped
!..from the data structure.
!..
!..nsrch is an integer variable. if nsrch is set to a value less than
!..or equal to n, then a different pivot option will be employed by
!..ma30a/ad.  this may result in different fill-in and execution
!..time for ma30a/ad. if nsrch is less than or equal to n, the
!..workspace arrays lastc and nextc are not referenced by ma30a/ad.
!..the default value for nsrch is 32768.
!..lbig is a logical variable. if lbig is set to .true., the value of
!..the largest entry encountered in the factorization by ma30a/ad
!..is returned in big.  setting lbig to .true.  will marginally
!..increase the factorization time for ma30a/ad and will increase
!..that for ma30b/bd by about 20%.  the default value for lbig is
!...false.
!..
!..declare
      logical          abort1,abort2,abort3,lbig
      integer          lp,ndrop,nsrch  
      real             eps,rmin,tol,big
!..
      common /ma30ed/ lp,abort1,abort2,abort3
      common /ma30gd/ eps,rmin
      common /ma30id/ tol,big,ndrop,nsrch,lbig
!..
      eps    = 1.0e-4
      tol    = 0.0e0
      big    = 0.0e0
      lp     = 6
      nsrch  = 32768
      lbig   = .false.
      abort1 = .true.
      abort2 = .true.
      abort3 = .false.
      return
      end
!..
!..
!..
!..
!..
      subroutine ma28int3
      implicit none
      logical         abort
      integer         lp,numnz,num,large
      common /mc23bd/ lp,numnz,num,large,abort
      lp       = 6
      abort    = .false.
      return
      end
!..
!..
!..
!..
!..
      subroutine mc20ad(nc,maxa,a,inum,jptr,jnum,jdisp)
      implicit none
      save 
!..
!..sorts a matrix into row order
!..
      integer          nc,maxa,inum(maxa),jnum(maxa),jptr(nc),jdisp, & 
     &                 null,j,k,kr,ice,jce,ja,jb,i,loc,icep,jcep
      real             a(maxa),ace,acep
!..
!..go
      null = -jdisp
      do 60 j=1,nc
       jptr(j) = 0
60    continue
!..
!..count the number of elements in each column.
      do 120 k=1,maxa
       j       = jnum(k) + jdisp
       jptr(j) = jptr(j) + 1
120   continue
!..
!..set the jptr array
      k = 1
      do 150 j=1,nc
       kr      = k + jptr(j)
       jptr(j) = k
       k       = kr
150   continue
!..
!..reorder the elements into column order; an in-place sort of order maxa.
!.. jce is the current entry.
      do 230 i=1,maxa
       jce = jnum(i) + jdisp
       if (jce .eq. 0) go to 230
       ace = a(i)
       ice = inum(i)
!..
!..clear the location vacated.
       jnum(i) = null
!..
!..chain from current entry to store items.
       do 200 j=1,maxa
        loc        = jptr(jce)
        jptr(jce)  = jptr(jce) + 1
        acep       = a(loc)
        icep       = inum(loc)
        jcep       = jnum(loc)
        a(loc)     = ace
        inum(loc)  = ice
        jnum(loc)  = null
        if (jcep .eq. null) go to 230
        ace = acep
        ice = icep
        jce = jcep + jdisp
200    continue
230   continue
!..
!..reset jptr vector.
      ja = 1
      do 250 j=1,nc
       jb      = jptr(j)
       jptr(j) = ja
       ja      = jb
250   continue
      return
      end
!..
!..
!..
!..
!..
      subroutine mc23ad(n,icn,a,licn,lenr,idisp,ip,iq,lenoff,iw,iw1)
      implicit none
      save 
!..
!..performs the block triangularization
!..
!..declare
      logical          abort
      integer          n,licn,idisp(2),iw1(n,2),icn(licn),lenr(n), & 
     &                 ip(n),iq(n),lenoff(n),iw(n,5),lp,numnz,num, & 
     &                 large,i,ii,ibeg,iend,i1,i2,k,iblock,jnpos, & 
     &                 ilend,inew,irowe,irowb,leni,nz,j,jj,iold,jold, & 
     &                 jnew
      real             a(licn)
      common /mc23bd/  lp,numnz,num,large,abort
!..
!..formats
180   format(1x,'matrix is structurally singular, rank = ',i6)
200   format(1x,'licn not big enough increase by ',i6)
220   format(1x,'error return from mc23ad because')
!..
!..set pointers iw(*,1) to beginning of the rows and set lenoff equal to lenr.
      iw1(1,1)  = 1
      lenoff(1) = lenr(1)
      if (n .eq. 1) go to 20
      do 10 i=2,n
       lenoff(i) = lenr(i)
       iw1(i,1)  = iw1(i-1,1) + lenr(i-1)
10    continue
!..
!..idisp(1) points to the first position in a/icn after the off-diagonal blocks 
!..and untreated rows.
20    idisp(1) = iw1(n,1) + lenr(n)
!..
!..find row permutation ip to make diagonal zero-free.
      call mc21a(n,icn,licn,iw1,lenr,ip,numnz,iw)
!..
!..possible error return for structurally singular matrices.
      if (numnz .ne. n  .and.  abort) go to 170
!..
!..iw1(*,2) and lenr are permutations of iw1(*,1) and lenr/lenoff suitable for 
!..entry to mc13d since matrix with these row pointer and length arrays has 
!..maximum number of non-zeros on the diagonal.
      do 30 ii=1,n
       i         = ip(ii)
       iw1(ii,2) = iw1(i,1)
       lenr(ii)  = lenoff(i)
30    continue
!..
!..find symmetric permutation iq to block lower triangular form.
      call mc13d(n,icn,licn,iw1(1,2),lenr,iq,iw(1,4),num,iw)
      if (num .ne. 1) go to 60
!..
!..action taken if matrix is irreducible. whole matrix is just moved to the 
!..end of the storage.
      do 40 i=1,n
       lenr(i) = lenoff(i)
       ip(i)   = i
       iq(i)   = i
40    continue
      lenoff(1) = -1
!..
!..idisp(1) is the first position after the last element in the off-diagonal 
!..blocks and untreated rows.
      nz       = idisp(1)-1
      idisp(1) = 1
!..
!..idisp(2) is position in a/icn of the first element in the diagonal blocks.
      idisp(2) = licn - nz + 1
      large    = n
      if (nz .eq. licn) go to 230
      do 50 k=1,nz
       j       = nz - k + 1
       jj      = licn - k + 1
       a(jj)   = a(j)
       icn(jj) = icn(j)
50    continue
      go to 230
!..
!..data structure reordered. form composite row permutation:ip(i) = ip(iq(i)).
60    do 70 ii=1,n
       i        = iq(ii)
       iw(ii,1) = ip(i)
70    continue
      do 80 i=1,n
       ip(i) = iw(i,1)
80    continue
!..
!..run through blocks in reverse order separating diagonal blocks which are 
!..moved to the end of the storage.  elements in off-diagonal blocks are left 
!..in place unless a compress is necessary.
!..ibeg indicates the lowest value of j for which icn(j) has been
!..     set to zero when element in position j was moved to the
!..     diagonal block part of storage.
!..iend is position of first element of those treated rows which are in 
!..     diagonal blocks.
!..large is the dimension of the largest block encountered so far.
!..num is the number of diagonal blocks.
!..i1 is first row (in permuted form) of block iblock.
!..i2 is last row (in permuted form) of block iblock.
      ibeg  = licn + 1
      iend  = licn + 1
      large = 0
      do 150 k=1,num
       iblock = num - k + 1
       i1 = iw(iblock,4)
       i2 = n
       if (k .ne. 1) i2 = iw(iblock+1,4) - 1
       large = max0(large,i2-i1+1)
!..
!..go through the rows of block iblock in the reverse order.
       do 140 ii=i1,i2
        inew = i2 - ii + i1
!..
!..we now deal with row inew in permuted form (row iold in original matrix).
        iold = ip(inew)
!..
!..if there is space to move up diagonal block portion of row go to 110
        if (iend-idisp(1) .ge. lenoff(iold)) go to 110
!..
!..in-line compress.; moves separated off-diagonal elements and untreated rows 
!..to front of storage.
        jnpos = ibeg
        ilend = idisp(1)-1
        if (ilend .lt. ibeg) go to 190
        do 90 j=ibeg,ilend
         if (icn(j) .eq. 0) go to 90
         icn(jnpos) = icn(j)
         a(jnpos)   = a(j)
         jnpos      = jnpos + 1
90      continue
        idisp(1) = jnpos
        if (iend-jnpos .lt. lenoff(iold)) go to 190
        ibeg = licn + 1
!..
!..reset pointers to the beginning of the rows.
        do 100 i=2,n
         iw1(i,1) = iw1(i-1,1) + lenoff(i-1)
100     continue
!..
!..row iold is now split into diag. and off-diag. parts.
110     irowb = iw1(iold,1)
        leni  = 0
        irowe = irowb+lenoff(iold)-1
!..
!..backward scan of whole of row iold (in original matrix).
        if (irowe .lt. irowb) go to 130
        do 120 jj=irowb,irowe
         j    = irowe - jj + irowb
         jold = icn(j)
!..
!..iw(.,2) holds the inverse permutation to iq.; it was set to this in mc13d.
         jnew = iw(jold,2)
!..
!..if (jnew.lt.i1) then element is in off-diagonal block and so is left in situ.
         if (jnew .lt. i1) go to 120
!..
!..element is in diagonal block and is moved to the end of the storage.
         iend      = iend-1
         a(iend)   = a(j)
         icn(iend) = jnew
         ibeg      = min0(ibeg,j)
         icn(j)    = 0
         leni      = leni + 1
120     continue
        lenoff(iold) = lenoff(iold) - leni
130     lenr(inew)   = leni
140    continue
       ip(i2)        = -ip(i2)
150   continue
!..
!..resets ip(n) to positive value.
!..idisp(2) is position of first element in diagonal blocks.
      ip(n)    = -ip(n)
      idisp(2) = iend
!..
!..this compress used to move all off-diagonal elements to the front of storage.
      if (ibeg .gt. licn) go to 230
      jnpos = ibeg
      ilend = idisp(1) - 1
      do 160 j=ibeg,ilend
       if (icn(j) .eq. 0) go to 160
       icn(jnpos) = icn(j)
       a(jnpos)   = a(j)
       jnpos      = jnpos + 1
160   continue
!..
!..idisp(1) is first position after last element of off-diagonal blocks.
      idisp(1) = jnpos
      go to 230
!..
!..error return
170   if (lp .ne. 0) write(lp,180) numnz
      idisp(1) = -1
      go to 210
190   if (lp .ne. 0) write(lp,200) n
      idisp(1) = -2
210   if (lp .ne. 0) write(lp,220)
230   return
      end
!..
!..
!..
!..
!..
      subroutine mc22ad(n,icn,a,nz,lenrow,ip,iq,iw,iw1)
      implicit none
      save 
!..
!..reorders the off diagonal blocks based on the pivot information
!..
!..declare
      integer          n,nz,iw(n,2),icn(nz),lenrow(n),ip(n),iq(n), & 
     &                 iw1(nz),i,jj,iold,j2,length,j,ipos,jval, & 
     &                 ichain,newpos,jnum
      real             a(nz),aval
!..
!..go
      if (nz .le. 0) go to 1000
      if (n  .le. 0) go to 1000
!..
!..set start of row i in iw(i,1) and lenrow(i) in iw(i,2)
      iw(1,1) = 1
      iw(1,2) = lenrow(1)
      do 10 i=2,n
       iw(i,1) = iw(i-1,1) + lenrow(i-1)
       iw(i,2) = lenrow(i)
10    continue
!..
!..permute lenrow according to ip.  set off-sets for new position of row iold 
!..in iw(iold,1) and put old row indices in iw1 in positions corresponding to 
!..the new position of this row in a/icn.
      jj = 1
      do 20 i=1,n
       iold      = ip(i)
       iold      = iabs(iold)
       length    = iw(iold,2)
       lenrow(i) = length
       if (length .eq. 0) go to 20
       iw(iold,1) = iw(iold,1) - jj
       j2 = jj + length - 1
       do 15 j=jj,j2
        iw1(j) = iold
15     continue
       jj = j2 + 1
20    continue
!..
!..set inverse permutation to iq in iw(.,2).
      do 30 i=1,n
       iold       = iq(i)
       iold       = iabs(iold)
       iw(iold,2) = i
30    continue
!..
!..permute a and icn in place, changing to new column numbers.
!..main loop; each pass through this loop places a closed chain of column 
!..indices in their new (and final) positions ... this is recorded by
!..setting the iw1 entry to zero so that any which are subsequently
!..encountered during this major scan can be bypassed.
      do 200 i=1,nz
       iold = iw1(i)
       if (iold .eq. 0) go to 200
       ipos = i
       jval = icn(i)
!..
!..if row iold is in same positions after permutation go to 150.
       if (iw(iold,1) .eq. 0) go to 150
       aval = a(i)
!..
!..chain loop; each pass through this loop places one (permuted) column index
!..in its final position  .. viz. ipos.
!..newpos is the original position in a/icn of the element to be placed
!..in position ipos.  it is also the position of the next element in the chain.
       do 100 ichain=1,nz
        newpos = ipos + iw(iold,1)
        if (newpos .eq. i) go to 130
        a(ipos)   = a(newpos)
        jnum      = icn(newpos)
        icn(ipos) = iw(jnum,2)
        ipos      = newpos
        iold      = iw1(ipos)
        iw1(ipos) = 0
100    continue
130    a(ipos)   = aval
150    icn(ipos) = iw(jval,2)
200   continue
1000  return
      end
!..
!..
!..
!..
!..
      subroutine mc21a(n,icn,licn,ip,lenr,iperm,numnz,iw)
      implicit none
      save 
      integer n,licn,ip(n),icn(licn),lenr(n),iperm(n),iw(n,4),numnz
      call mc21b(n,icn,licn,ip,lenr,iperm,numnz,iw(1,1),iw(1,2),iw(1,3), & 
     &           iw(1,4))
      return
      end
!..
!..
!..
!..
!..
      subroutine mc21b(n,icn,licn,ip,lenr,iperm,numnz,pr,arp,cv,out)
      implicit none
      save 
!..
!..does a row permutation to make the diagonal zero free
!..
!..pr(i) is the previous row to i in the depth first search.
!..     it is used as a work array in the sorting algorithm.
!..     elements (iperm(i),i) i=1, ... n  are non-zero at the end of the
!..     algorithm unless n assignments have not been made.  in which case
!..(iperm(i),i) will be zero for n-numnz entries.
!..cv(i)  is the most recent row extension at which column i was visited.
!..arp(i) is one less than the number of non-zeros in row i
!..       which have not been scanned when looking for a cheap assignment.
!..out(i) is one less than the number of non-zeros in row i
!..       which have not been scanned during one pass through the main loop.
!..
!..declare
      integer n,licn,ip(n),icn(licn),lenr(n),iperm(n),pr(n),cv(n), & 
     &        arp(n),out(n),i,jord,j,in1,in2,k,ii,ioutk,j1,kk,numnz
!..
!..initialization of arrays.
      do 10 i=1,n
       arp(i)   = lenr(i)-1
       cv(i)    = 0
       iperm(i) = 0
10    continue
      numnz=0
!..
!..main loop. each pass round this loop either results in a new assignment
!..or gives a row with no assignment.
      do 130 jord=1,n
       j     = jord
       pr(j) = -1
       do 100 k=1,jord
!..
!..look for a cheap assignment
        in1 = arp(j)
        if (in1 .lt. 0) go to 60
        in2 = ip(j) + lenr(j)-1
        in1 = in2 - in1
        do 50 ii=in1,in2
         i = icn(ii)
         if (iperm(i) .eq. 0) go to 110
50      continue
!..
!..no cheap assignment in row.
!..begin looking for assignment chain starting with row j.
        arp(j) = -1
60      out(j) = lenr(j)-1
!..
!..c inner loop.  extends chain by one or backtracks.
        do 90 kk=1,jord
         in1 = out(j)
         if (in1 .lt. 0) go to 80
         in2 = ip(j)+lenr(j)-1
         in1 = in2 - in1
!..
!..forward scan.
         do 70 ii=in1,in2
          i = icn(ii)
          if (cv(i) .eq. jord) go to 70
!..
!..column i has not yet been accessed during this pass.
          j1      = j
          j       = iperm(i)
          cv(i)   = jord
          pr(j)   = j1
          out(j1) = in2-ii-1
          go to 100
70       continue
!..
!..backtracking step.
80       j = pr(j)
         if (j .eq. -1) go to 130
90      continue
100    continue
!..
!..new assignment is made.
110    iperm(i) = j
       arp(j)   = in2 - ii - 1
       numnz    = numnz + 1
       do 120 k=1,jord
        j = pr(j)
        if (j .eq. -1) go to 130
        ii       = ip(j) + lenr(j) - out(j) - 2
        i        = icn(ii)
        iperm(i) = j
120    continue
130   continue
!..
!..if matrix is structurally singular, we now complete the permutation iperm.
      if (numnz .eq. n) return
      do 140 i=1,n
       arp(i) = 0
140   continue
      k = 0
      do 160 i=1,n
       if (iperm(i) .ne. 0) go to 150
       k      = k + 1
       out(k) = i
       go to 160
150    j      = iperm(i)
       arp(j) = i
160   continue
      k = 0
      do 170 i=1,n
       if (arp(i) .ne. 0) go to 170
       k            = k+1
       ioutk        = out(k)
       iperm(ioutk) = i
170   continue
      return
      end
!..
!..
!..
!..
!..
      subroutine mc13d(n,icn,licn,ip,lenr,ior,ib,num,iw)
      implicit none
      save 
      integer n,licn,ip(n),icn(licn),lenr(n),ior(n),ib(n),iw(n,3),num
      call mc13e(n,icn,licn,ip,lenr,ior,ib,num,iw(1,1),iw(1,2),iw(1,3))
      return
      end
!..
!..
!..
!..
!..
      subroutine mc13e(n,icn,licn,ip,lenr,arp,ib,num,lowl,numb,prev)
      implicit none
      save 
!..
!.. arp(i) is one less than the number of unsearched edges leaving
!..        node i.  at the end of the algorithm it is set to a
!..        permutation which puts the matrix in block lower
!..        triangular form.
!..ib(i)   is the position in the ordering of the start of the ith
!..        block.  ib(n+1-i) holds the node number of the ith node
!..        on the stack.
!..lowl(i) is the smallest stack position of any node to which a path
!..        from node i has been found.  it is set to n+1 when node i
!..        is removed from the stack.
!..numb(i) is the position of node i in the stack if it is on
!..        it, is the permuted order of node i for those nodes
!..        whose final position has been found and is otherwise zero.
!..prev(i) is the node at the end of the path when node i was
!..        placed on the stack.
!..
!..declare
      integer n,licn,stp,dummy,ip(n),icn(licn),lenr(n),arp(n),ib(n), & 
     &        lowl(n),numb(n),prev(n),icnt,num,nnm1,j,iv,ist,i1,i2, & 
     &        ii,iw,ist1,lcnt,i,isn,k
!..
!..
!..icnt is number of nodes whose positions in final ordering have been found.
!..num is the number of blocks that have been found.
      icnt = 0
      num  = 0
      nnm1 = n + n-1
!..
!..initialization of arrays.
      do 20 j=1,n
       numb(j) = 0
       arp(j)  = lenr(j)-1
20    continue
!..
!..look for a starting node
!..ist is the number of nodes on the stack ... it is the stack pointer.
      do 120 isn=1,n
       if (numb(isn) .ne. 0) go to 120
       iv  = isn
       ist = 1
!..
!..put node iv at beginning of stack.
       lowl(iv) = 1
       numb(iv) = 1
       ib(n)    = iv
!..
!..the body of this loop puts a new node on the stack or backtracks.
       do 110 dummy=1,nnm1
        i1 = arp(iv)
!..
!..have all edges leaving node iv been searched.
        if (i1 .lt. 0) go to 60
        i2 = ip(iv) + lenr(iv) - 1
        i1 = i2 - i1
!..
!..look at edges leaving node iv until one enters a new node or all edges are 
!..exhausted.
        do 50 ii=i1,i2
         iw = icn(ii)
         if (numb(iw) .eq. 0) go to 100
         lowl(iv) = min0(lowl(iv),lowl(iw))
50      continue
!..
!..there are no more edges leaving node iv.
        arp(iv) = -1
!..
!..is node iv the root of a block.
60      if (lowl(iv) .lt. numb(iv)) go to 90
!..
!..order nodes in a block.
        num  = num + 1
        ist1 = n + 1 - ist
        lcnt = icnt + 1
!..
!..peel block off the top of the stack starting at the top and working down to 
!..the root of the block.
        do 70 stp=ist1,n
         iw       = ib(stp)
         lowl(iw) = n + 1
         icnt     = icnt + 1
         numb(iw) = icnt
         if (iw .eq. iv) go to 80
70      continue
80      ist     = n - stp
        ib(num) = lcnt
!..
!..are there any nodes left on the stack.
        if (ist .ne. 0) go to 90
!..
!..have all the nodes been ordered.
        if (icnt .lt. n) go to 120
        go to 130
!..
!..backtrack to previous node on path.
90      iw = iv
        iv = prev(iv)
!..
!..update value of lowl(iv) if necessary.
        lowl(iv) = min0(lowl(iv),lowl(iw))
        go to 110
!..
!..put new node on the stack.
100     arp(iv)  = i2 - ii - 1
        prev(iw) = iv
        iv       = iw
        ist      = ist+1
        lowl(iv) = ist
        numb(iv) = ist
        k        = n+1-ist
        ib(k)    = iv
110    continue
120   continue
!..
!..put permutation in the required form.
130   do 140 i=1,n
       ii      = numb(i)
       arp(ii) = i
140   continue
      return
      end
!..
!..
!..
!..
!..
      subroutine mc24ad(n,icn,a,licn,lenr,lenrl,w)
      implicit none
      save 
!..
!..computes the gwoth rate of fill in
!..
      integer          n,licn,icn(licn),lenr(n),lenrl(n),i,j0,j2,j1,jj,j
      real             a(licn),w(n),amaxl,wrowl,amaxu,zero
      data             zero/0.0e0/
!..
!..initialize
      amaxl = zero
      do 10 i=1,n
       w(i) = zero
10    continue
      j0 = 1
      do 100 i=1,n
       if (lenr(i) .eq. 0) go to 100
       j2=j0+lenr(i)-1
       if (lenrl(i) .eq. 0) go to 50
!..
!..calculation of 1-norm of l.
       j1 = j0 + lenrl(i) - 1
       wrowl=zero
       do 30 jj=j0,j1
        wrowl = wrowl + abs(a(jj))
30     continue
!..
!..amaxl is the maximum norm of columns of l so far found.
       amaxl = max(amaxl,wrowl)
       j0    = j1 + 1
!..
!..calculation of norms of columns of u (max-norms).
50     j0 = j0 + 1
       if (j0 .gt. j2) go to 90
       do 80 jj=j0,j2
        j    = icn(jj)
        w(j) = max(abs(a(jj)),w(j))
80     continue
90     j0 = j2 + 1
100   continue
!..
!..amaxu is set to maximum max-norm of columns of u.
      amaxu = zero
      do 200 i=1,n
       amaxu = max(amaxu,w(i))
200   continue
!..
!..grofac is max u max-norm times max l 1-norm.
      w(1) = amaxl*amaxu
      return
      end
!..
!..
!..
!..
!..
      subroutine mc20bd(nc,maxa,a,inum,jptr)
      implicit none
      save 
!..
!..never called
!..
!..
!..declare
      integer          nc,maxa,inum(maxa),jptr(nc),kmax,jj,j,klo,kor, & 
     &                 kdummy,ice,k,ik
      real             a(maxa),ace
!..
!..go
      kmax=maxa
      do 35 jj=1,nc
       j   = nc + 1 - jj
       klo = jptr(j)+1
       if (klo .gt. kmax) go to 30
       kor=kmax
!..
!..items kor, kor+1, .... ,kmax are in order
       do 25 kdummy=klo,kmax
        ace = a(kor-1)
        ice = inum(kor-1)
        do 10 k=kor,kmax
         ik = inum(k)
         if (iabs(ice) .le. iabs(ik)) go to 20
         inum(k-1) = ik
         a(k-1)    = a(k)
10      continue
        k         = kmax+1
20      inum(k-1) = ice
        a(k-1)    = ace
        kor       = kor-1
25     continue
30     kmax = klo - 2
35    continue
      return
      end
