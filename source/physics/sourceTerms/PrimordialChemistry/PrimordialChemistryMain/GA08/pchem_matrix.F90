!!***if* source/physics/sourceTerms/PrimordialChemistry/PrimordialChemistryMain/GA08/pchem_matrix
!!
!! NAME
!!  pchem_matrix
!!
!!
!! SYNOPSIS
!!
!!  Call the needed subroutines
!!
!!
!! DESCRIPTION
!!  This is going to hold the two new fuctions I need to 
!!  decompose and solve the dense matrix that I have
!!
!!  These two functions are from Numerical Recipies.
!!
!!  The first is matrix_ludcmp. 
!!
!!***

  subroutine ludcmp(a,n,np,indx,d)
    INTEGER n,np,indx(n),NMAX
    REAL d,a(np,np),TINY
    PARAMETER (NMAX=500, TINY=1.0e-20)

    INTEGER i,imax,j,k
    REAL aamax,dum,sum,vv(NMAX)
    d=1.
    do i=1,n
       aamax = 0.
       do j=1,n
          if (abs(a(i,j)) .gt. aamax) aamax=abs(a(i,j))
       enddo
       if (aamax .eq. 0.) stop 'singular matrix in ludcmp'
       vv(i) = 1./aamax
    enddo 

    do j=1,n
         do i=1,j-1
            sum = a(i,j)
            do k=1,i-1
               sum = sum-a(i,k)*a(k,j)
            enddo
            a(i,j) = sum
         enddo
         aamax = 0.
         do i=j,n
             sum=a(i,j)
             do k=1,j-1
                sum=sum-a(i,k)*a(k,j)
             enddo
             a(i,j) = sum
             dum = vv(i)*abs(sum)
             if (dum .ge. aamax) then
                imax = i
                aamax = dum
             endif
         enddo
         if (j .ne. imax) then
             do k=1,n
                  dum = a(imax,k)
                  a(imax,k) = a(j,k)
                  a(j,k) = dum
             enddo
             d=-d
             vv(imax)=vv(j)
         endif
         indx(j)=imax
         if(a(j,j) .eq. 0.) a(j,j) = TINY
         if (j .ne. n) then
            dum = 1./a(j,j)
            do i=j+1,n
                 a(i,j) = a(i,j)*dum
            enddo
         endif
     enddo
     return 
     END






  subroutine lubksb(a,n,np,indx,b)
  
  INTEGER n,np,indx(n)
  REAL a(np,np),b(n)
  INTEGER i,ii,j,ll
  REAL sum

  ii=0

  do i=1,n
       ll = indx(i)
       sum = b(ll)
       b(ll) = b(i)
       if (ii .ne. 0) then
          do j = ii,i-1
               sum = sum-a(i,j)*b(j)
          enddo
       else if (sum .ne. 0.) then
           ii=i
       endif
       b(i) = sum
  enddo
 

  do i=n,1,-1
       sum = b(i)
          do j=i+1,n
             sum = sum - a(i,j)*b(j)
          enddo
       b(i) = sum/a(i,i)
  enddo
  return
  END    


  subroutine leqs(a,b,n,np)
  implicit none
  save


  integer  n,np,n1,i,j,k,l,imax,jj
  real     a(np,np),b(np),r,c
 
!  print *, 'inside leqs'  

  n1=n-1
  do i=1,n
     r = abs(a(i,1))
     do j=2,n
        c = abs(a(i,j))
        if (r .lt. c) r=c
     enddo

     do j=1,n
        a(i,j) = a(i,j)/r
     enddo
     b(i) = b(i)/r
  enddo


  do j=1,n1
     l = j+1
     do i=l,n
        r = -a(i,j)
        if (r .eq. 0.0) goto 50
        r = r/a(j,j)
        do k=l,n
           a(i,k) = a(i,k) + r*a(j,k)
        enddo
        b(i) = b(i) + r*b(j)
50      continue
     enddo
  enddo


  b(n) = b(n)/a(n,n)
  do l=1,n1
     i = n-l
     r = 0.0e0
     imax = i+1
     do j = imax, n
        jj = i+n+1-j
         r = r+a(i,jj)*b(jj)
     enddo
     b(i) = (b(i)-r)/a(i,i)
  enddo
  return
  end
  