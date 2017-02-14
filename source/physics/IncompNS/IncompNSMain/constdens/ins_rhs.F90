
      SUBROUTINE ins_rhs2d(uni,vni,ru1,ix1,ix2,jy1,jy2,dx,dy,ru,rv)

  !***************************************************************
  ! This subroutine computes the discretization of the RHS of the 
  ! Helmholtz equation on a staggered uniform grid.
  !
  ! Input:  uni,vni     = velocity at timestep n
  !         ru1         = molecular viscosity
  !         ix1,ix2     = starting and ending x indices
  !         jy1,jy2     = starting and ending y indices
  !         dx,dy       = grid spacing in x and y directions
  !
  ! Output: ru,rv    = u and v momentum for Helmholtz RHS
  !**************************************************************

      implicit none
      INTEGER, INTENT(IN):: ix1, ix2, jy1, jy2
      REAL, INTENT(IN):: ru1, dx, dy
      REAL, DIMENSION(:,:,:), INTENT(IN):: uni, vni
      REAL, DIMENSION(:,:,:), INTENT(OUT):: ru, rv

      INTEGER:: i, j
      REAL:: dx1, dy1
      ! x-component variables
      REAL:: uxplus, uxminus, vxplus, vxminus, wxplus, wxminus, &
             uyplus, uyminus
      REAL:: dudxp, dudxm, dudyp, dudym, dvdxp, dvdxm
      REAL:: tvjp, tvjm
      REAL:: txxp, txxm, tyyp, tyym
      REAL:: txyp, txym
      ! new y-component variables
      REAL:: vyplus, vyminus
      REAL:: dvdyp, dvdym
      REAL:: tvip, tvim
      INTEGER, parameter :: kz1 = 1

      ! grid spacings
      dx1 = 1.0/dx
      dy1 = 1.0/dy

      !++++++++++  U-COMPONENT  ++++++++++
       do j = jy1,jy2
          do i = ix1,ix2+1
             ! get velocities at 1/2 locations
             uxplus = (uni(i+1,j,kz1) + uni(i,j,kz1))*0.5
             uxminus = (uni(i,j,kz1) + uni(i-1,j,kz1))*0.5

             vxplus = (vni(i,j+1,kz1) + vni(i-1,j+1,kz1))*0.5
             vxminus = (vni(i,j,kz1) + vni(i-1,j,kz1))*0.5

             uyplus = (uni(i,j+1,kz1) + uni(i,j,kz1))*0.5
             uyminus = (uni(i,j,kz1) + uni(i,j-1,kz1))*0.5

             ! get derivatives at 1/2 locations
             dudxp = (uni(i+1,j,kz1) - uni(i,j,kz1))*dx1
             dudxm = (uni(i,j,kz1) - uni(i-1,j,kz1))*dx1
             dudyp = (uni(i,j+1,kz1) - uni(i,j,kz1))*dy1
             dudym = (uni(i,j,kz1) - uni(i,j-1,kz1))*dy1
             dvdxp = (vni(i,j+1,kz1) - vni(i-1,j+1,kz1))*dx1
             dvdxm = (vni(i,j,kz1) - vni(i-1,j,kz1))*dx1

             ! flux of normal total stresses
             txxp = ru1*dudxp
             txxm = ru1*dudxm
             tyyp = ru1*dudyp
             tyym = ru1*dudym

             ! calculate RHS for u-momentum
             ru(i,j,kz1) =                                              &
                          - (uxplus*uxplus - uxminus*uxminus)*dx1       &! advection term
                          - (vxplus*uyplus - vxminus*uyminus)*dy1       &                          
                          + (txxp - txxm)*dx1                           &! diffusion - normal terms 
                          + (tyyp - tyym)*dy1
          enddo
       enddo

    !++++++++++  V-COMPONENT  ++++++++++

       do j = jy1,jy2+1
          do i = ix1,ix2
             ! get velocities at 1/2 locations
             vxplus = (vni(i+1,j,kz1) + vni(i,j,kz1))*0.5
             vxminus = (vni(i,j,kz1) + vni(i-1,j,kz1))*0.5

             vyplus = (vni(i,j+1,kz1) + vni(i,j,kz1))*0.5
             vyminus = (vni(i,j,kz1) + vni(i,j-1,kz1))*0.5

             uyplus = (uni(i+1,j,kz1) + uni(i+1,j-1,kz1))*0.5
             uyminus = (uni(i,j,kz1) + uni(i,j-1,kz1))*0.5

             ! get derivatives at 1/2 locations
             dvdxp = (vni(i+1,j,kz1) - vni(i,j,kz1))*dx1
             dvdxm = (vni(i,j,kz1) - vni(i-1,j,kz1))*dx1
             dvdyp = (vni(i,j+1,kz1) - vni(i,j,kz1))*dy1
             dvdym = (vni(i,j,kz1) - vni(i,j-1,kz1))*dy1
             dudyp = (uni(i+1,j,kz1) - uni(i+1,j-1,kz1))*dy1
             dudym = (uni(i,j,kz1) - uni(i,j-1,kz1))*dy1

             ! flux of normal total stresses
             txxp = ru1*dvdxp
             txxm = ru1*dvdxm
             tyyp = ru1*dvdyp
             tyym = ru1*dvdym

             ! calculate RHS for v-momentum
             rv(i,j,kz1) =                                              &
                          - (uyplus*vxplus - uyminus*vxminus)*dx1       &! advection term
                          - (vyplus*vyplus - vyminus*vyminus)*dy1       &
                          + (txxp - txxm)*dx1                           &! diffusion - normal terms
                          + (tyyp - tyym)*dy1
          enddo
       enddo

       END SUBROUTINE ins_rhs2d




      SUBROUTINE ins_rhs3d(uni,vni,wni,tv,ru1,      &
                           ix1,ix2,jy1,jy2,kz1,kz2, &
                           dx,dy,dz,ru,rv,rw)

  !*****************************************************************
  ! This subroutine computes the centered discretization of the RHS 
  ! of the momentum equation (advection + viscous terms) on a 
  ! staggered uniform grid based on the Paramesh grid structure.
  !
  ! Input:  uni,vni,wni = velocity at timestep n
  !         tv          = eddy viscosity
  !         ru1         = molecular viscosity
  !         ix1,ix2     = starting and ending x indices
  !         jy1,jy2     = starting and ending y indices
  !         kz1,kz2     = starting and ending z indices
  !         dx,dy,dz    = grid spacing in x, y, and z directions
  !
  ! Output: ru,rv,rw    = RHS of u, v, and w momentum equations
  !
  ! E. Balaras   July 1999
  ! P. Rabenold  August 2006
  !**************************************************************

      implicit none
      INTEGER, INTENT(IN):: ix1, ix2, jy1, jy2, kz1, kz2
      REAL, INTENT(IN):: ru1, dx, dy, dz
      REAL, DIMENSION(:,:,:), INTENT(IN):: uni, vni, wni, tv
      REAL, DIMENSION(:,:,:), INTENT(OUT):: ru, rv, rw

!!$      REAL, DIMENSION(nx+1,ny,nz), INTENT(IN):: uni
!!$      REAL, DIMENSION(nx,ny+1,nz), INTENT(IN):: vni
!!$      REAL, DIMENSION(nx,ny,nz+1), INTENT(IN):: wni
!!$      REAL, DIMENSION(nx,ny,nz)  , INTENT(IN):: tv
!!$
!!$      REAL, DIMENSION(nx+1,ny,nz), INTENT(OUT):: ru
!!$      REAL, DIMENSION(nx,ny+1,nz), INTENT(OUT):: rv
!!$      REAL, DIMENSION(nx,ny,nz+1), INTENT(OUT):: rw

      INTEGER:: i, j, k
      REAL:: dx1, dy1, dz1
      ! x-component variables
      REAL:: uxplus, uxminus, vxplus, vxminus, wxplus, wxminus, &
             uyplus, uyminus, uzplus, uzminus
      REAL:: dudxp, dudxm, dudyp, dudym, dudzp, dudzm, dvdxp, dvdxm, &
             dwdxp, dwdxm
      REAL:: tvjp, tvjm, tvkp, tvkm
      REAL:: txxp, txxm, tyyp, tyym, tzzp, tzzm
      REAL:: txyp, txym, txzp, txzm
      ! additional y-component variables
      REAL:: vyplus, vyminus, vzplus, vzminus, wyplus, wyminus
      REAL:: dvdyp, dvdym, dvdzp, dvdzm, dwdyp, dwdym
      REAL:: tvip, tvim
      REAL:: tyzp, tyzm
      ! additional z-component variables
      REAL:: wzplus, wzminus
      REAL:: dwdzp, dwdzm

      ! grid spacings
      dx1 = 1.0/dx
      dy1 = 1.0/dy
      dz1 = 1.0/dz


      !++++++++++  U-COMPONENT  ++++++++++
      do k = kz1,kz2
         do j = jy1,jy2
            do i = ix1,ix2+1

               ! get velocities at 1/2 locations
               uxplus  = (uni(i+1,j  ,k  ) + uni(i  ,j  ,k  ))*0.5
               uxminus = (uni(i  ,j  ,k  ) + uni(i-1,j  ,k  ))*0.5

               uyplus  = (uni(i  ,j+1,k  ) + uni(i  ,j  ,k  ))*0.5
               uyminus = (uni(i  ,j  ,k  ) + uni(i  ,j-1,k  ))*0.5

               uzplus  = (uni(i  ,j  ,k+1) + uni(i  ,j  ,k  ))*0.5
               uzminus = (uni(i  ,j  ,k  ) + uni(i  ,j  ,k-1))*0.5

               vxplus  = (vni(i  ,j+1,k  ) + vni(i-1,j+1,k  ))*0.5
               vxminus = (vni(i  ,j  ,k  ) + vni(i-1,j  ,k  ))*0.5

               wxplus  = (wni(i  ,j  ,k+1) + wni(i-1,j  ,k+1))*0.5
               wxminus = (wni(i  ,j  ,k  ) + wni(i-1,j  ,k  ))*0.5

               ! get derivatives at 1/2 locations
               dudxp = (uni(i+1,j  ,k  ) - uni(i  ,j  ,k  ))*dx1
               dudxm = (uni(i  ,j  ,k  ) - uni(i-1,j  ,k  ))*dx1
               dudyp = (uni(i  ,j+1,k  ) - uni(i  ,j  ,k  ))*dy1
               dudym = (uni(i  ,j  ,k  ) - uni(i  ,j-1,k  ))*dy1
               dudzp = (uni(i  ,j  ,k+1) - uni(i  ,j  ,k  ))*dz1
               dudzm = (uni(i  ,j  ,k  ) - uni(i  ,j  ,k-1))*dz1
               dvdxp = (vni(i  ,j+1,k  ) - vni(i-1,j+1,k  ))*dx1
               dvdxm = (vni(i  ,j  ,k  ) - vni(i-1,j  ,k  ))*dx1
               dwdxp = (wni(i  ,j  ,k+1) - wni(i-1,j  ,k+1))*dx1
               dwdxm = (wni(i  ,j  ,k  ) - wni(i-1,j  ,k  ))*dx1

               ! get nu_t

               ! ****** requires DIAGONALS for corner ghost cells ******

               tvjp = 0.25*(tv(i-1,j  ,k  ) + tv(i  ,j  ,k  ) + &
                            tv(i  ,j+1,k  ) + tv(i-1,j+1,k  ))
               tvjm = 0.25*(tv(i-1,j-1,k  ) + tv(i  ,j-1,k  ) + &
                            tv(i  ,j  ,k  ) + tv(i-1,j  ,k  ))
               tvkp = 0.25*(tv(i-1,j  ,k  ) + tv(i  ,j  ,k  ) + &
                            tv(i  ,j  ,k+1) + tv(i-1,j  ,k+1))
               tvkm = 0.25*(tv(i-1,j  ,k-1) + tv(i  ,j,  k-1) + &
                            tv(i  ,j  ,k  ) + tv(i-1,j  ,k  ))

               ! flux of normal total stresses
               txxp = (ru1 + 2.0*tv(i,j,k))*dudxp
               txxm = (ru1 + 2.0*tv(i-1,j,k))*dudxm
               tyyp = (ru1 + tvjp)*dudyp
               tyym = (ru1 + tvjm)*dudym
               tzzp = (ru1 + tvkp)*dudzp
               tzzm = (ru1 + tvkm)*dudzm

               ! flux of cross SGS stresses
               txyp = tvjp*dvdxp
               txym = tvjm*dvdxm
               txzp = tvkp*dwdxp
               txzm = tvkm*dwdxm

               ! calculate RHS for u-momentum
               ru(i,j,k) =                                          &                              
                           - (uxplus*uxplus - uxminus*uxminus)*dx1  &! advection term
                           - (vxplus*uyplus - vxminus*uyminus)*dy1  &
                           - (wxplus*uzplus - wxminus*uzminus)*dz1  &             
                           + (txxp - txxm)*dx1                      &! diffusion - normal terms
                           + (tyyp - tyym)*dy1                      &
                           + (tzzp - tzzm)*dz1                      &
                           + (txyp - txym)*dy1                      &! diffusion - cross terms
                           + (txzp - txzm)*dz1   
            enddo
         enddo
      enddo

      !++++++++++  V-COMPONENT  ++++++++++

      do k = kz1,kz2
         do j = jy1,jy2+1
            do i = ix1,ix2

               ! get velocities at 1/2 locations
               vxplus  = (vni(i+1,j  ,k  ) + vni(i  ,j  ,k  ))*0.5
               vxminus = (vni(i  ,j  ,k  ) + vni(i-1,j  ,k  ))*0.5

               vyplus  = (vni(i  ,j+1,k  ) + vni(i  ,j  ,k  ))*0.5
               vyminus = (vni(i  ,j  ,k  ) + vni(i  ,j-1,k  ))*0.5

               vzplus  = (vni(i  ,j  ,k+1) + vni(i  ,j  ,k  ))*0.5
               vzminus = (vni(i  ,j  ,k  ) + vni(i  ,j  ,k-1))*0.5

               uyplus  = (uni(i+1,j  ,k  ) + uni(i+1,j-1,k  ))*0.5
               uyminus = (uni(i  ,j  ,k  ) + uni(i  ,j-1,k  ))*0.5

               wyplus  = (wni(i  ,j  ,k+1) + wni(i  ,j-1,k+1))*0.5
               wyminus = (wni(i  ,j  ,k  ) + wni(i  ,j-1,k  ))*0.5

               ! get derivatives at 1/2 locations
               dvdxp = (vni(i+1,j  ,k  ) - vni(i  ,j  ,k  ))*dx1
               dvdxm = (vni(i  ,j  ,k  ) - vni(i-1,j  ,k  ))*dx1
               dvdyp = (vni(i  ,j+1,k  ) - vni(i  ,j  ,k  ))*dy1
               dvdym = (vni(i  ,j  ,k  ) - vni(i  ,j-1,k  ))*dy1
               dvdzp = (vni(i  ,j  ,k+1) - vni(i  ,j  ,k  ))*dz1
               dvdzm = (vni(i  ,j  ,k  ) - vni(i  ,j  ,k-1))*dz1
               dudyp = (uni(i+1,j  ,k  ) - uni(i+1,j-1,k  ))*dy1
               dudym = (uni(i  ,j  ,k  ) - uni(i  ,j-1,k  ))*dy1
               dwdyp = (wni(i  ,j  ,k+1) - wni(i  ,j-1,k+1))*dy1
               dwdym = (wni(i  ,j  ,k  ) - wni(i  ,j-1,k  ))*dy1

               ! get nu_t

               ! ****** requires DIAGONALS for corner ghost cells ******

               tvip = 0.25*(tv(i  ,j-1,k  ) + tv(i+1,j-1,k  ) + &
                            tv(i+1,j  ,k  ) + tv(i  ,j  ,k  ))
               tvim = 0.25*(tv(i-1,j-1,k  ) + tv(i  ,j-1,k  ) + &
                            tv(i  ,j  ,k  ) + tv(i-1,j  ,k  ))
               tvkp = 0.25*(tv(i  ,j-1,k  ) + tv(i  ,j-1,k+1) + &
                            tv(i  ,j  ,k+1) + tv(i  ,j  ,k  ))
               tvkm = 0.25*(tv(i  ,j-1,k-1) + tv(i  ,j-1,k  ) + &
                            tv(i  ,j  ,k  ) + tv(i  ,j  ,k-1))

               ! flux of normal total stresses
               txxp = (ru1 + tvip)*dvdxp
               txxm = (ru1 + tvim)*dvdxm
               tyyp = (ru1 + 2.0*tv(i,j,k))*dvdyp
               tyym = (ru1 + 2.0*tv(i,j-1,k))*dvdym
               tzzp = (ru1 + tvkp)*dvdzp
               tzzm = (ru1 + tvkm)*dvdzm

               ! flux of cross SGS stresses
               txyp = tvip*dudyp
               txym = tvim*dudym
               tyzp = tvkp*dwdyp
               tyzm = tvkm*dwdym

               ! calculate RHS for v-momentum
               rv(i,j,k) =                                           &
                           - (uyplus*vxplus - uyminus*vxminus)*dx1   &! advection term
                           - (vyplus*vyplus - vyminus*vyminus)*dy1   &
                           - (wyplus*vzplus - wyminus*vzminus)*dz1   &
                           + (txxp - txxm)*dx1                       &! diffusion - normal terms
                           + (tyyp - tyym)*dy1                       &
                           + (tzzp - tzzm)*dz1                       &
                           + (txyp - txym)*dx1                       &! diffusion - cross terms
                           + (tyzp - tyzm)*dz1
            enddo
         enddo
      enddo

      !++++++++++  W-COMPONENT  ++++++++++
      
      do k = kz1,kz2+1
         do j = jy1,jy2
            do i = ix1,ix2

               ! get velocities at 1/2 locations
               wxplus  = (wni(i+1,j  ,k  ) + wni(i  ,j  ,k  ))*0.5
               wxminus = (wni(i  ,j  ,k  ) + wni(i-1,j  ,k  ))*0.5
               
               wyplus  = (wni(i  ,j+1,k  ) + wni(i  ,j  ,k  ))*0.5
               wyminus = (wni(i  ,j  ,k  ) + wni(i  ,j-1,k  ))*0.5
               
               wzplus  = (wni(i  ,j  ,k+1) + wni(i  ,j  ,k  ))*0.5
               wzminus = (wni(i  ,j  ,k  ) + wni(i  ,j  ,k-1))*0.5

               uzplus  = (uni(i+1,j  ,k  ) + uni(i+1,j  ,k-1))*0.5
               uzminus = (uni(i  ,j  ,k  ) + uni(i  ,j  ,k-1))*0.5

               vzplus  = (vni(i  ,j+1,k  ) + vni(i  ,j+1,k-1))*0.5
               vzminus = (vni(i  ,j  ,k  ) + vni(i  ,j  ,k-1))*0.5

               ! get derivatives at 1/2 locations
               dwdxp = (wni(i+1,j  ,k  ) - wni(i  ,j  ,k  ))*dx1
               dwdxm = (wni(i  ,j  ,k  ) - wni(i-1,j  ,k  ))*dx1
               dwdyp = (wni(i  ,j+1,k  ) - wni(i  ,j  ,k  ))*dy1
               dwdym = (wni(i  ,j  ,k  ) - wni(i  ,j-1,k  ))*dy1
               dwdzp = (wni(i  ,j  ,k+1) - wni(i  ,j  ,k  ))*dz1
               dwdzm = (wni(i  ,j  ,k  ) - wni(i  ,j  ,k-1))*dz1
               dudzp = (uni(i+1,j  ,k  ) - uni(i+1,j  ,k-1))*dz1
               dudzm = (uni(i  ,j  ,k  ) - uni(i  ,j  ,k-1))*dz1
               dvdzp = (vni(i  ,j+1,k  ) - vni(i  ,j+1,k-1))*dz1
               dvdzm = (vni(i  ,j  ,k  ) - vni(i  ,j  ,k-1))*dz1

               ! get nu_t

               ! ****** requires DIAGONALS for corner ghost cells ******

               tvip = 0.25*(tv(i  ,j  ,k-1) + tv(i+1,j  ,k-1) + &
                            tv(i+1,j  ,k  ) + tv(i  ,j  ,k  ))
               tvim = 0.25*(tv(i-1,j  ,k-1) + tv(i  ,j  ,k-1) + &
                            tv(i  ,j  ,k  ) + tv(i-1,j  ,k  ))
               tvjp = 0.25*(tv(i  ,j  ,k-1) + tv(i  ,j  ,k  ) + &
                            tv(i  ,j+1,k  ) + tv(i  ,j+1,k-1))
               tvjm = 0.25*(tv(i  ,j-1,k-1) + tv(i  ,j-1,k  ) + & 
                            tv(i  ,j  ,k  ) + tv(i  ,j  ,k-1))

               ! flux of normal total stresses
               txxp = (ru1 + tvip)*dwdxp
               txxm = (ru1 + tvim)*dwdxm
               tyyp = (ru1 + tvjp)*dwdyp
               tyym = (ru1 + tvjm)*dwdym
               tzzp = (ru1 + 2.0*tv(i,j,k))*dwdzp
               tzzm = (ru1 + 2.0*tv(i,j,k-1))*dwdzm

               ! flux of cross SGS stresses
               txzp = tvip*dudzp
               txzm = tvim*dudzm
               tyzp = tvjp*dvdzp
               tyzm = tvjm*dvdzm

               ! calculate RHS for w-momentum
               rw(i,j,k) =                                          &
                          - (uzplus*wxplus - uzminus*wxminus)*dx1   &! advection term
                          - (vzplus*wyplus - vzminus*wyminus)*dy1   &
                          - (wzplus*wzplus - wzminus*wzminus)*dz1   &
                          + (txxp - txxm)*dx1                       &! diffusion - normal terms
                          + (tyyp - tyym)*dy1                       &
                          + (tzzp - tzzm)*dz1                       &
                          + (txzp - txzm)*dx1                       &! diffusion - cross terms
                          + (tyzp - tyzm)*dy1                       
            enddo
         enddo
      enddo

      END SUBROUTINE ins_rhs3d

