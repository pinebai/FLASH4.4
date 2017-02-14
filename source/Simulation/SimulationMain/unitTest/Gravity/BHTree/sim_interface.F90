Module sim_interface


  interface
    subroutine sim_computeError()
    end subroutine sim_computeError
  end interface

  interface
     subroutine sim_spharm(l, m, theta, phi, Ylm, dYdt, dYdp)
       implicit none
       integer l, m
       real theta, phi
       real Ylm, dYdt, dYdp
       real PM(0:l), PD(0:l), fact
     end subroutine sim_spharm
  end interface

  interface
     SUBROUTINE sim_LPMNS(M,N,X,PM,PD)
!
!       ========================================================
!       Purpose: Compute associated Legendre functions Pmn(x)
!                and Pmn'(x) for a given order
!       Input :  x --- Argument of Pmn(x)
!                m --- Order of Pmn(x),  m = 0,1,2,...,n
!                n --- Degree of Pmn(x), n = 0,1,2,...,N
!       Output:  PM(n) --- Pmn(x)
!                PD(n) --- Pmn'(x)
!       Copyright: Jianming Jin
!       http://jin.ece.uiuc.edu/routines/routines.html
!       ========================================================
!
       IMPLICIT none
       integer,intent(IN) :: M, N
       real,intent(IN)    :: X
       real,intent(OUT)   :: PM(0:N),PD(0:N)
     END SUBROUTINE sim_LPMNS
  end interface

end Module sim_interface

