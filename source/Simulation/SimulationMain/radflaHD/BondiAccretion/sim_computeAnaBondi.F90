!!****if* source/Simulation/SimulationMain/radflaHD/BondiAccretion/sim_computeAnaBondi
!!
!! NAME
!!
!!  sim_computeAnaBondi
!!
!! SYNOPSIS
!!
!!  call sim_computeAnaBondi(real,intent(IN)  :: r,
!!                           real,intent(OUT)  :: velr,
!!                           real,intent(OUT)  :: rho)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   r : 
!!
!!   velr : 
!!
!!   rho : 
!!
!! AUTOGENROBODOC
!!
!!
!!***

subroutine sim_computeAnaBondiScaled(x,u,alpha)
  implicit none

  real,intent(IN) :: x
  real,intent(OUT) :: u,alpha

  real,parameter :: small_ug = 1.e-20
  real,parameter :: tol = 1.e-8
  real :: ug,gx,ff,f,u_final
  integer :: j
  integer :: max_iter=1000

  gx = (1./x)+2.*log(x)-(3./2.)+log(4.)

  ug = 1./(x*2.0)

  do j =1,max_iter

     f = ug**2/2.-log(ug)-gx
     ff = ug - 1.0/ug

     if (ff.eq.0.0) then
        u_final = 1.0
        alpha = (exp(3./2.)/4.)/(x**2)
        exit
     end if

     ug = ug - f/ff

     if(ug <= 0.)then
        ug = small_ug
     endif

!!     print*,'TEST 1',j,ug,f,ff,gx

     if((abs(f/ff)/ug).le.tol)then

        u_final = ug
        alpha = (exp(3./2.)/4.)/(u_final*x**2)

        exit

     endif

  enddo

  if(j > max_iter) then
     print*,'The NR methods did not converge for x =',x
     stop 'This is BAAAD!'
  end if

  u = u_final
  

end subroutine sim_computeAnaBondiScaled

subroutine sim_computeAnaBondi(r,velR,rho)

  use sim_interface, ONLY: sim_computeAnaBondiScaled
  use Simulation_data, ONLY: sim_rho_vac, sim_soundSpeedInf,&
       sim_bondiRadius

  implicit none

  real,intent(IN) :: r
  real,intent(OUT) :: velR,rho

  real :: u,alpha

  call sim_computeAnaBondiScaled(r/sim_bondiRadius,u,alpha)
  velR = u * (-sim_soundSpeedInf)
  rho = alpha * sim_rho_vac

end subroutine sim_computeAnaBondi
