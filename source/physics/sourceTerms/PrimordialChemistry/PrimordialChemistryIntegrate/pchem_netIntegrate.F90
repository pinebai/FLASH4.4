!!****if* source/physics/sourceTerms/PrimordialChemistry/PrimordialChemistryIntegrate/pchem_netIntegrate
!!
!! NAME
!!
!!  pchem_netIntegrate
!!
!! SYNOPSIS
!!   subroutine pchem_netIntegrate( real(IN)    :: start,
!!                         real(IN)    :: stptry,
!!                         real(IN)    :: stpmin,
!!                         real(IN)    :: stopp,
!!                         real(INOUT)    :: bc(:),
!!                         real(IN)     :: eps,
!!                         real(IN)     :: dxsav,
!!                         integer(IN)  :: kmax, 
!!                         real(OUT)    :: xrk(:),
!!                         real(OUT)    :: yrk(:,:),   ! called with elem
!!                         integer(IN)  :: xphys,
!!                         integer(IN)  :: yphys,
!!                         integer(IN)  :: xlogi,
!!                         integer(IN)  :: ylogi,   
!!                         integer(OUT) :: nok,
!!                         integer(OUT) :: nbad,
!!                         integer(OUT) :: kount,
!!                         real(IN)     :: odescal,
!!                         integer(IN)  :: iprint,
!!                         procedure(IN):: derivs,
!!                         procedure(IN):: jakob,
!!                         procedure(IN):: bjakob,
!!                         procedure(IN):: steper)   
!!   
!! DESCRIPTION
!!
!!   ode integrator for stiff odes with an analytic and sparse jacobian.
!!   
!! ARGUMENTS
!!
!!   start   - real(IN)     beginning integration point
!!   stptry  - real(IN)     suggested first step size
!!   stpmin  - real(IN)     minimum allowable step size
!!   stopp   - real(IN)     ending integration point
!!   bc(:)   - real(INOUT)     initial conditions, array of physical dimension yphys
!!   eps     - real(IN)     desired fraction error during the integration
!!   dxsav   - real(IN)     incremental vale of indep variable at which to store the solution
!!                          if zero, solution is stored at every intermediate point
!!                          if not zero, solution is done and saved at every dxsav point
!!   kmax    - integer(IN)  maximum number of solution points to store, kkmax < xphys
!!   xrk(:)  - real(OUT)    the independent variable solution 
!!                          array of physical dimension xphys, logical dimension xlogi
!!   yrk(:,:)- real(OUT)    the dependent variable solution 
!!                          array of physical dimension (yphys,xphys) with 
!!                          logical  dimension (ylogi,xlogi)
!!   xphys   - integer(IN)  physical size of array xrk,yrk
!!   yphys   - integer(IN)  physical size of array yrk
!!   xlogi   - integer(IN)  logical size of array xrk,yrk
!!   ylogi   - integer(IN)  logical size of array yrk
!!   nok     - integer(OUT) number of succesful steps taken
!!   nbad    - integer(OUT) number of bad steps taken, bad but retried and then succesful
!!   kount   - integer(OUT) total number of steps stored in arrays xrk and yrk
!!   odescal - real(IN)     error scaling factor 
!!   iprint  - integer(IN)  determines if the solution is printed as it evolves
!!   derivs  - procedure(IN) name of the routine that contains the odes
!!   jakob   - procedure(IN) name of the routine that contains the jacobian of the odes
!!   bjakob  - procedure(IN) name of the routine that sets the pointers of the sparse jacobian
!!   steper  - procedure(IN) name of the routine that will take a single step
!!   
!! NOTES
!!  this file uses 8 routines to drive the integration of 
!!  nuclear reaction networks with 2 choices of the time stepper
!!  and 2 choices for the linear pchem_algebra.
!!  routine pchem_netIntegrate drives the integration of the odes
!!  routine pchem_baderMa28 drives a bader-deuflhard step with ma28 pchem_algebra
!!  routine pchem_baderStepMa28 is a bader-deuflhard stepper with ma28 pchem_algebra
!!  routine pchem_baderGift drives a bader-deuflhard step with gift pchem_algebra
!!  routine pchem_baderStepGift is a bader-deuflhard stepper with ma28 pchem_algebra
!!  routine pchem_pzExtr does extrapolations for any of the Step_* routines
!!  routine pchem_rosenMa28 drives a rosenbrock step with ma28 pchem_algebra
!!  routine pchem_rosenGift drives a rosenbrock step with gift pchem_algebra
!!***   




!!---------------------------------------------------------------------------------

subroutine pchem_netIntegrate(start,stptry,stpmin,stopp,bc,  & 
     &                  eps,dxsav,kmax,   & 
     &                  xrk,yrk,xphys,yphys,xlogi,ylogi,  & 
     &                  nok,nbad,kount,odescal,iprint,  & 
     &                  derivs,jakob,bjakob,steper,jcounts)
   
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none


  ! For some mysterious reason, can't use the following line
  ! use bnNetwork_interface, ONLY: derivs, jakob, bjakob
  ! use pchem_integrateInterface, ONLY: steper.
  ! But you CAN directly put the interfaces into the file.  Go figure.
  ! Also note that the "use" lines have to be ABOVE implicit none, but
  !  the direct interfaces have to be BELOW it.
  interface   ! = pchem_network
     subroutine derivs(tt,y,dydt)   !! == pchem_network
       implicit none
       real, intent(IN) :: tt
       real, intent(INOUT), dimension(*)  :: y
       real, intent(OUT), dimension(*) :: dydt
     end subroutine derivs

     subroutine jakob(tt,y,dfdy,nzo,nDummy) ! = pchem_networkSparseJakob or pchem_networkDenseJakob
       implicit none                        ! See notes in bnNetwork_interface about this...
       integer, intent(IN) :: nzo, nDummy
       real, intent(IN)    :: tt
       real, intent(INOUT) :: y(*)
       real, intent(OUT)   :: dfdy(nzo,nDummy)
     end subroutine jakob

     subroutine bjakob(iloc,jloc,nzo,np)
       implicit none
       integer, intent(IN)  ::   iloc(*),jloc(*),np
       integer, intent(OUT) ::   nzo
     end subroutine bjakob

     subroutine steper(y,dydx,nv,x,htry,eps,yscal,hdid,hnext, & 
          &                       derivs,jakob,bjakob,jcounts)
       implicit none
       external               derivs,jakob,bjakob
       integer, intent(IN) :: nv
       real, intent(INOUT) :: y(nv)
       real, intent(IN)    :: dydx(nv), yscal(nv), htry, eps
       real, intent(OUT)   :: hdid, hnext
       real, intent(INOUT) :: x, jcounts
     end subroutine steper
  end interface

  !! arguments
  !! Note, can't give an INTENT statement with external functions
  !!  And if you don't use the interface above, you need the following line
  !  external                    derivs,jakob,bjakob, steper 

  integer, intent(IN)  :: xphys,yphys,xlogi,ylogi
  integer, intent(IN)  :: kmax, iprint
  real, intent(IN)     :: odescal, dxsav, eps
  real, intent(IN)     :: start, stptry, stpmin, stopp
  real, intent(INOUT), dimension(yphys) :: bc

  integer, intent(OUT) :: nok, nbad, kount
  real, intent(OUT), dimension(xphys)       :: xrk
  real, intent(OUT), dimension(yphys,xphys) :: yrk
  real, intent(INOUT) :: jcounts

  !! local declarations

  integer, parameter :: nmax=30   !! this should be 13 for Aprox13
  integer, parameter :: stpmax=10000
  integer            :: i,j,nstp 
  real, parameter    :: zero=0.0e0
  real, parameter    :: one=1.0e0
  real, parameter    :: tiny=1.0e-15
  real, save         :: yscal(nmax),y(nmax),dydx(nmax),    & 
       &                 x,xsav,h,hdid,hnext

  save

  !!   here are the format statements for printouts as we integrate 
100 format(1x,i4,1pe10.2) 
101 format(1x,1p12e10.2) 
102 format(1x,5(a,' ',1pe10.2))
103 format(1x,5(a,' ',i6))

  !------------------------------------------------------------------------------

  !!   initialize    
  if (ylogi .gt. yphys) then
     write(*,*) 'ylogi > yphys in routine pchem_netIntegrate' 
     call Driver_abortFlash('ERROR in pchem_netIntegrate: ylogi > yphys')
  end if

  if (yphys .gt. nmax)  then
     write(*,*) 'yphys > nmax in routine pchem_netIntegrate' 
     call Driver_abortFlash('ERROR in pchem_netIntegrate: yphys > nmax')
  end if

  x     = start    
  h     = sign(stptry,stopp-start)  
  nok   = 0  
  nbad  = 0 
  kount = 0    


  !!   store the first step  
  do i=1,ylogi 
     y(i) = bc(i)
!     print*, 'INSIDE BN_NI: y(',i,')= ', y(i)
      
  enddo
  xsav = x - 2.0e0 * dxsav 

  !!   take at most stpmax steps 
  do nstp=1,stpmax 


     !!   positive definite mass fractions 
     do i=1,ylogi 
        y(i) = max(y(i),1.0e-30) 
!        print *, 'Y(',i,')= ', y(i)
     enddo
     !! dummy procedure name, expands to pchem_network
     call derivs(x,y,dydx) 


     !!   scaling vector used to monitor accuracy 
     do i=1,ylogi 
        yscal(i) = max(odescal,abs(y(i))) 
     enddo

     !!   store intermediate results    
     if (kmax .gt. 0) then 
        if ( (abs(dxsav) - abs(x-xsav)) .le. tiny) then  
           if ( kount .lt. (kmax-1) ) then   
              kount         = kount+1   
              xrk(kount)    = x    
              do i=1,ylogi  
                 yrk(i,kount) = y(i) 
              enddo
              if (iprint .eq. 1) then 
                 write(*,100) kount,xrk(kount) 
                 write(*,101) (yrk(j,kount), j=1,ylogi) 
              end if
              xsav=x  
           end if
        end if
     end if


     !!   if the step can overshoot the stop point or the dxsav increment then cut it 
     if ((x+h-stopp)*(x+h-start).gt.zero) h = stopp - x   
     if (dxsav.ne.zero .and. h.gt.(xsav-x+dxsav)) h = xsav + dxsav-x 

     !!   do an integration step 
     !! dummy procedure names in this dummy procedure call
     !! for Aprox13, expands to steper = pchem_bader/rosen|Gift/Ma28 
     !! Passes on the names of 
     !! derivs = pchem_aprox13
     !! jakob  = pchem_saprox13 (for sparse ma28 solver) or pchem_daprox13 (for dense gift solver)
     !! bjakob = pchem_baprox13 (really used only for ma28 solver)
     !!   
!     do i=1,ylogi
!        print *, 'Before steper: y(',i,')= ', y(i)
!     enddo


     call steper(y,dydx,ylogi,x,h,eps,yscal,hdid,hnext, & 
          &             derivs,jakob,bjakob,jcounts)    
     if (hdid.eq.h) then 
        nok = nok+1    
     else  
        nbad = nbad+1  
     end if


     !!   this is the normal exit point, save the final step    
     if (nstp .eq. stpmax .or. (x-stopp)*(stopp-start) .ge. zero) then 
        do i=1,ylogi   
           bc(i) = y(i)  
        enddo
                

        if (kmax.ne.0) then    
           kount         = kount+1   
           xrk(kount)    = x    
           do i=1,ylogi  
              yrk(i,kount) = y(i)  
           end do
           if (iprint .eq. 1) then 
              write(*,100) kount,xrk(kount) 
              write(*,101) (yrk(j,kount), j=1,ylogi) 
           end if
        end if
        return   
     end if

     !!   set the step size for the next iteration; stay above stpmin 
     h = hnext 
     if (abs(hnext).lt.stpmin) then
        write(*,*) ' '
        write(*,102) 'hnext=',hnext,' stpmin=',stpmin
        write(*,*) 'hnext < stpmin in pchem_netIntegrate' 
        write(*,*)' '
        call Driver_abortFlash('ERROR in pchem_netIntegrate: hnext < stpmin')
     end if

     !!   back for another iteration or death 
  enddo

  ! crash here if no normal return occurred above
  write(*,*) 
  write(*,103)'stpmax=',stpmax
  write(*,*)'more than stpmax steps required in pchem_netIntegrate'  
  write(*,*) 
  call Driver_abortFlash('ERROR in pchem_netIntegrate: too many steps required')

end subroutine pchem_netIntegrate
