!!****if* source/Simulation/SimulationMain/unitTest/Gravity/BHTree/sim_spharm
!!
!! NAME
!!
!!  sim_spharm
!!
!! SYNOPSIS
!!
!!  call sim_spharm(:: l,
!!                   :: m,
!!                   :: theta,
!!                   :: phi,
!!                   :: ylm,
!!                   :: dydt,
!!                   :: dydp)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   l : 
!!
!!   m : 
!!
!!   theta : 
!!
!!   phi : 
!!
!!   ylm : 
!!
!!   dydt : 
!!
!!   dydp : 
!!
!! AUTOGENROBODOC
!!
!!
!!***

subroutine sim_spharm(l, m, theta, phi, Ylm, dYdt, dYdp)
  use sim_interface, ONLY: sim_LPMNS
  implicit none
  integer l, m
  real theta, phi
  real Ylm, dYdt, dYdp
  real PM(0:l), PD(0:l), fact

  fact = 1.0
  
  ! no normalization neccesary - normalized to 1 elsewhere anyway
  !do i = l-m+1, l+m
  !  fact = fact*i
  !enddo
  !fact = sqrt((2*l+1)/(4*pi) / fact)

  call sim_LPMNS(m,l,cos(theta),PM,PD)
  
  Ylm  = fact * PM(l) * cos(m*phi)
  dYdt = fact * PD(l) * (-sin(theta)) * cos(m*phi)
  dYdp = fact * PM(l) * (-m*sin(m*phi))

end subroutine sim_spharm




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

  integer :: K
  real    :: X0, PM0, PMK, PM1, PM2
  DO K=0,N
     PM(K)=0.0D0
     PD(K)=0.0D0
  ENDDO
  IF (ABS(X).EQ.1.0D0) THEN
     DO K=0,N
        IF (M.EQ.0) THEN
           PM(K)=1.0D0
           PD(K)=0.5D0*K*(K+1.0)
           IF (X.LT.0.0) THEN
              PM(K)=(-1)**K*PM(K)
              PD(K)=(-1)**(K+1)*PD(K)
           ENDIF
        ELSE IF (M.EQ.1) THEN
           PD(K)=huge(0.0)
        ELSE IF (M.EQ.2) THEN
           PD(K)=-0.25D0*(K+2.0)*(K+1.0)*K*(K-1.0)
           IF (X.LT.0.0) PD(K)=(-1)**(K+1)*PD(K)
        ENDIF
     ENDDO
     RETURN
  ENDIF
  X0=ABS(1.0D0-X*X)
  PM0=1.0D0
  PMK=PM0
  DO K=1,M
     PMK=(2.0D0*K-1.0D0)*SQRT(X0)*PM0
     PM0=PMK
  ENDDO
  PM1=(2.0D0*M+1.0D0)*X*PM0
  PM(M)=PMK
  PM(M+1)=PM1
  DO K=M+2,N
     PM2=((2.0D0*K-1.0D0)*X*PM1-(K+M-1.0D0)*PMK)/(K-M)
     PM(K)=PM2
     PMK=PM1
     PM1=PM2
  ENDDO
  PD(0)=((1.0D0-M)*PM(1)-X*PM(0))/(X*X-1.0)  
  DO K=1,N
     PD(K)=(K*X*PM(K)-(K+M)*PM(K-1))/(X*X-1.0D0)
  ENDDO
  RETURN
END SUBROUTINE sim_LPMNS
