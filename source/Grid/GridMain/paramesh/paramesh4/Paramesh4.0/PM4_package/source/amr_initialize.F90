!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"


      subroutine amr_initialize 


!----------------------------------------------------------------


! This subroutine initializes the amr package. It performs any
! initialization required by the package which is application
! independent.

!
! NOTE : This routine MUST BE the first executed code in your application!!!!
!

      use paramesh_dimensions
      use physicaldata
      use tree
      use workspace
      use mpi_morton
      use timings
      use prolong_arrays
      use timings
      use paramesh_comm_data

      use paramesh_interfaces, only : amr_1blk_guardcell_reset, & 
     &                                amr_prolong_fun_init, & 
     &                                amr_bcset_init

      implicit none

      include 'mpif.h'

      integer :: nfield, nprocs, mype, maxprocs
      integer :: i
      real    :: test1,test2
      real    :: real_test_in(2),real_test_out(2)
      integer :: ierr
      double precision :: time1

!----------------------------------------------------------------

! START UP MPI

      call comm_start(maxprocs,nprocs,mype)
      Call MPI_BARRIER(amr_mpi_meshComm, ierr)

!-----Set the type for passing reals using MPI
#ifdef REAL8
      amr_mpi_real = MPI_DOUBLE_PRECISION
#else
      amr_mpi_real = MPI_REAL
#endif

!--
! test mpi real communitcation
      real_test_in(1) = 100. + real(mype)
      real_test_in(2) = 1000. + 2.*real(mype)

! reduce
      call mpi_real_allreduce(real_test_in, real_test_out, 2, & 
     &     amr_mpi_real, & 
     &     MPI_MAX, amr_mpi_meshComm, ierr)
      test1 = 100. + real(nprocs-1)
      test2 = 1000. + 2.*real(nprocs-1)
      if(real_test_out(1).ne.test1.or.real_test_out(2).ne.test2) then
        if(mype.eq.0) then
        write(*,*) 'PARAMESH ERROR: A test of mpi communication of ', & 
     &             'REAL type failed. Possible cause is using -r8 ', & 
     &             'without defining the preprocessor variable REAL8' & 
     &      ,' real_test_in ',real_test_in,' real_test_out ', & 
     &        real_test_out,' test1  test2 ', test1, test2
        endif
        Call MPI_BARRIER(amr_mpi_meshComm,ierr)
        call amr_abort()
      endif
!--

! allocate storage for paramesh arrays
! this allows for runtime definition of variables such as ,e.g., maxblocks

      call amr_set_runtime_parameters()

      if (timing_mpi) then
         time1 = mpi_wtime()
      end if

#ifdef LIBRARY

      clean_divb = .FALSE.

      if (no_permanent_guardcells) then
         npgs = 0
      else
         npgs = 1
      end if

      if (var_dt) then
         advance_all_levels = .true.
      end if

      if (curvilinear .and. .not.cartesian_pm) then
         consv_fluxes = .true.
         consv_flux_densities = .false.
         edge_value_integ = .true.
         edge_value = .false.
      elseif (.not.curvilinear) then
         curvilinear_conserve = .false.
         cartesian_pm = .false.
         cylindrical_pm = .false.
         spherical_pm = .false.
         polar_pm = .false.
      end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! To be put in paramesh runtime control file !!!!

      if (ndim == 2) then
         nzb = 1
      endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      nbedges = ndim*2**(ndim-1)
      k3d=(ndim-1)/2
      k2d=ndim/2
      k1d=1
      red_f = 0.25
      if (consv_fluxes) then
         if (ndim == 3) then
            red_f = 1.0
         elseif (ndim == 2) then
            red_f = 0.5
         elseif (ndim == 1) then
            red_f = 0.25
         endif
      end if
      nchild=2**ndim
      nfaces=2*ndim
      if (nboundaries < 2*ndim) then
         nboundaries=2*ndim
      endif

      if (nvar_work < 2) nvar_work = 2
      nbndvar=max(1,nfacevar)
      nbndvare=max(1,nvaredge)
      nbndvarc=max(1,nvarcorn)

      nfluxes=max(1,nfluxvar)
      nbndmax=max(nbndvar,nfluxes)

      nedgevar=max(nedgevar1,nvaredge)
      nedges=max(1,nedgevar)

      maxdim=max(nxb,nyb,nzb)
      gc_off_x=mod(nxb,2)
      gc_off_y=mod(nyb,2)
      gc_off_z=mod(nzb,2)
      il_bnd=1
      jl_bnd=1
      kl_bnd=1
      iu_bnd=nxb+2*nguard*npgs
      ju_bnd=nyb+2*nguard*npgs*k2d
      ku_bnd=nzb+2*nguard*npgs*k3d
      il_bndi=nguard*npgs+1
      iu_bndi=nguard*npgs+nxb
      jl_bndi=nguard*npgs*k2d+1
      ju_bndi=nguard*npgs*k2d+nyb
      kl_bndi=nguard*npgs*k3d+1
      ku_bndi=nguard*npgs*k3d+nzb
      len_block=iu_bnd*ju_bnd*ku_bnd*nvar
      len_blockfx=(iu_bnd+1)*ju_bnd*ku_bnd
      len_blockfy=iu_bnd*(ju_bnd+k2d)*ku_bnd
      len_blockfz=iu_bnd*ju_bnd*(ku_bnd+k3d)
      len_blockex=iu_bnd*(ju_bnd+k2d)*(ku_bnd+k3d)
      len_blockey=(iu_bnd+1)*ju_bnd*(ku_bnd+k3d)
      len_blockez=(iu_bnd+1)*(ju_bnd+k2d)*ku_bnd
      len_blockn=(iu_bnd+1)*(ju_bnd+k2d)*(ku_bnd+k3d)
      len_blockfxf=2*ju_bnd*ku_bnd
      len_blockfyf=iu_bnd*2*ku_bnd
      len_blockfzf=iu_bnd*ju_bnd*2
      il_bnd1=1
      jl_bnd1=1
      kl_bnd1=1
      iu_bnd1=nxb+2*nguard
      ju_bnd1=nyb+2*nguard*k2d
      ku_bnd1=nzb+2*nguard*k3d
      len_block1 = iu_bnd1*ju_bnd1*ku_bnd1*nvar
      len_blockfx1 = (iu_bnd1+1)*ju_bnd1*ku_bnd1 
      len_blockfy1 = iu_bnd1*(ju_bnd1+k2d)*ku_bnd1
      len_blockfz1 = iu_bnd1*ju_bnd1*(ku_bnd1+k3d)
      len_blockex1 = iu_bnd1*(ju_bnd1+k2d)*(ku_bnd1+k3d)
      len_blockey1 = (iu_bnd1+1)*ju_bnd1*(ku_bnd1+k3d)
      len_blockez1 = (iu_bnd1+1)*(ju_bnd1+1)*ku_bnd1 
      len_blockn1 = (iu_bnd1+1)*(ju_bnd1+1)*(ku_bnd1+k3d)
      ilw=1
      jlw=1
      klw=1
      ngw2=2*nguard_work
      iuw=nxb+ngw2*npgs
      juw=nyb+ngw2*npgs*k2d
      kuw=nzb+ngw2*npgs*k3d
      len_wblock=iuw*juw*kuw
      ilw1=1
      jlw1=1
      klw1=1
      iuw1=nxb+ngw2
      juw1=nyb+ngw2*k2d
      kuw1=nzb+ngw2*k3d
      len_wblock1=iuw1*juw1*kuw1
      if (ndim == 1) then
         nmax_lays = nxb/2
      end if
      if (ndim == 2) then
         nmax_lays = min(nxb/2,nyb/2)
      end if
      if (ndim == 3) then
         nmax_lays = min(nxb/2,nyb/2,nzb/2)
      end if

      mxblks_buf = max(1,maxblocks/100)
      maxblocks_alloc = maxblocks * 10

      maxblocksf = 1+(maxblocks-1)*min(1,nfacevar)
      maxblocksue = 1+(maxblocks-1)*min(1,nvaredge)
      maxblocksn = 1+(maxblocks-1)*min(1,nvarcorn)
      maxblocks_gt=(maxblocks-1)*(1-npgs)+1
      maxblocksf_gt=(maxblocksf-1)*(1-npgs)+1
      maxblocksue_gt=(maxblocksue-1)*(1-npgs)+1
      maxblocksn_gt=(maxblocksn-1)*(1-npgs)+1
      maxblocksfl=1+(maxblocks-1)*min(1,nfluxvar)
      maxblockse=1+(maxblocks-1)*min(1,nedgevar)

      allocate( & 
     & unk(nvar, & 
     &     il_bnd:iu_bnd, & 
     &     jl_bnd:ju_bnd, & 
     &     kl_bnd:ku_bnd, & 
     &     maxblocks))
      allocate(interp_mask_unk(nvar))
      allocate(interp_mask_unk_res(nvar))

      allocate(gcell_on_cc_pointer(nvar))
      allocate(gcell_on_cc(nvar))
      allocate(int_gcell_on_cc(nvar))

      allocate(gcell_on_fc_pointer(3,nbndvar))
      allocate(gcell_on_fc(3,nbndvar))
      allocate(int_gcell_on_fc(3,nbndvar))

      allocate(gcell_on_ec_pointer(3,nbndvare))
      allocate(gcell_on_ec(3,nbndvare))
      allocate(int_gcell_on_ec(3,nbndvare))

      allocate(gcell_on_nc_pointer(nbndvarc))
      allocate(gcell_on_nc(nbndvarc))
      allocate(int_gcell_on_nc(nbndvarc))

      allocate(checkp_on_cc(nvar))
      allocate(checkp_on_fc(3,nbndvar))
      allocate(checkp_on_ec(3,nbndvare))
      allocate(checkp_on_nc(nbndvarc))

      allocate( & 
     & facevarx(nbndvar, & 
     &          il_bnd:iu_bnd+1, & 
     &          jl_bnd:ju_bnd, & 
     &          kl_bnd:ku_bnd, & 
     &          maxblocksf))
      allocate( & 
     & facevary(nbndvar, & 
     &          il_bnd:iu_bnd, & 
     &          jl_bnd:ju_bnd+k2d, & 
     &          kl_bnd:ku_bnd, & 
     &          maxblocksf))
      allocate( & 
     & facevarz(nbndvar, & 
     &          il_bnd:iu_bnd, & 
     &          jl_bnd:ju_bnd, & 
     &          kl_bnd:ku_bnd+k3d, & 
     &          maxblocksf))
      allocate(interp_mask_facex(nbndvar))
      allocate(interp_mask_facey(nbndvar))
      allocate(interp_mask_facez(nbndvar))
      allocate(interp_mask_facex_res(nbndvar))
      allocate(interp_mask_facey_res(nbndvar))
      allocate(interp_mask_facez_res(nbndvar))
      allocate( & 
     & unk_e_x(nbndvare, & 
     &         il_bnd:iu_bnd, & 
     &         jl_bnd:ju_bnd+k2d, & 
     &         kl_bnd:ku_bnd+k3d, & 
     &         maxblocksue))
      allocate( & 
     & unk_e_y(nbndvare, & 
     &         il_bnd:iu_bnd+1, & 
     &         jl_bnd:ju_bnd, & 
     &         kl_bnd:ku_bnd+k3d, & 
     &         maxblocksue))
      allocate( & 
     & unk_e_z(nbndvare, & 
     &         il_bnd:iu_bnd+1, & 
     &         jl_bnd:ju_bnd+k2d, & 
     &         kl_bnd:ku_bnd, & 
     &         maxblocksue))
      allocate(interp_mask_ec(nbndvare))
      allocate(interp_mask_ec_res(nbndvare))
      allocate(unk_n(nbndvarc, & 
     &         il_bnd:iu_bnd+1,jl_bnd:ju_bnd+k2d, & 
     &         kl_bnd:ku_bnd+k3d, & 
     &         maxblocksn))
      allocate(interp_mask_nc(nbndvarc))
      allocate(interp_mask_nc_res(nbndvarc))
      allocate(unk1(nvar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1, & 
     &              kl_bnd1:ku_bnd1, npblks))
      allocate(facevarx1(nbndvar,il_bnd1:iu_bnd1+1, & 
     &                   jl_bnd1:ju_bnd1, & 
     &                   kl_bnd1:ku_bnd1,npblks))
      allocate(facevary1(nbndvar,il_bnd1:iu_bnd1, & 
     &                   jl_bnd1:ju_bnd1+k2d, & 
     &                   kl_bnd1:ku_bnd1,npblks))
      allocate(facevarz1(nbndvar,il_bnd1:iu_bnd1, & 
     &                   jl_bnd1:ju_bnd1, & 
     &                   kl_bnd1:ku_bnd1+k3d,npblks))
      allocate(unk_e_x1(nbndvare, & 
     &                  il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1+k2d, & 
     &                  kl_bnd1:ku_bnd1+k3d, & 
     &                  npblks))
      allocate(unk_e_y1(nbndvare, & 
     &                  il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1, & 
     &                  kl_bnd1:ku_bnd1+k3d, & 
     &                  npblks))
      allocate(unk_e_z1(nbndvare, & 
     &                  il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d, & 
     &                  kl_bnd1:ku_bnd1, & 
     &                  npblks))
      allocate(unk_n1(nbndvarc, & 
     &                il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d, & 
     &                kl_bnd1:ku_bnd1+k3d, & 
     &                npblks))
      allocate(unk1_fl(nvar, & 
     &                 il_bnd1:iu_bnd1+nxb+2*nguard, & 
     &                 jl_bnd1:ju_bnd1+(nyb+2*nguard)*k2d, & 
     &                 kl_bnd1:ku_bnd1+(nzb+2*nguard)*k3d ))
      allocate(facevarx1_fl(nbndvar, & 
     &                      il_bnd1:iu_bnd1+nxb+2*nguard+1, & 
     &                      jl_bnd1:ju_bnd1+(nyb+2*nguard)*k2d, & 
     &                      kl_bnd1:ku_bnd1+(nzb+2*nguard)*k3d ))
      allocate(facevary1_fl(nbndvar, & 
     &                      il_bnd1:iu_bnd1+nxb+2*nguard, & 
     &                      jl_bnd1:ju_bnd1+(nyb+2*nguard+1)*k2d, & 
     &                      kl_bnd1:ku_bnd1+(nzb+2*nguard)*k3d ))
      allocate(facevarz1_fl(nbndvar, & 
     &                      il_bnd1:iu_bnd1+nxb+2*nguard, & 
     &                      jl_bnd1:ju_bnd1+(nyb+2*nguard)*k2d, & 
     &                      kl_bnd1:ku_bnd1+(nzb+2*nguard+1)*k3d ))
      allocate(unk_e_x1_fl(nbndvare, & 
     &                     il_bnd1:iu_bnd1+nxb+2*nguard, & 
     &                     jl_bnd1:ju_bnd1+(nyb+2*nguard+1)*k2d, & 
     &                     kl_bnd1:ku_bnd1+(nzb+2*nguard+1)*k3d ))
      allocate(unk_e_y1_fl(nbndvare, & 
     &                     il_bnd1:iu_bnd1+nxb+2*nguard+1, & 
     &                     jl_bnd1:ju_bnd1+(nyb+2*nguard)*k2d, & 
     &                     kl_bnd1:ku_bnd1+(nzb+2*nguard+1)*k3d ))
      allocate(unk_e_z1_fl(nbndvare, & 
     &                     il_bnd1:iu_bnd1+nxb+2*nguard+1, & 
     &                     jl_bnd1:ju_bnd1+(nyb+2*nguard+1)*k2d, & 
     &                     kl_bnd1:ku_bnd1+(nzb+2*nguard)*k3d ))
      allocate(unk_n1_fl(nbndvarc, & 
     &                   il_bnd1:iu_bnd1+nxb+2*nguard+1, & 
     &                   jl_bnd1:ju_bnd1+(nyb+2*nguard+1)*k2d, & 
     &                   kl_bnd1:ku_bnd1+(nzb+2*nguard+1)*k3d ))
      allocate(time_loc(maxblocks_alloc))
      allocate(ldtcomplete(maxblocks_alloc))


      allocate( & 
     & flux_x(nfluxes,1:2, & 
     &        jl_bndi:ju_bndi,kl_bndi:ku_bndi,maxblocksfl))
      allocate( & 
     & flux_y(nfluxes,il_bndi:iu_bndi, & 
     &        1:2,kl_bndi:ku_bndi,maxblocksfl))
      allocate( & 
     & flux_z(nfluxes,il_bndi:iu_bndi, & 
     &        jl_bndi:ju_bndi,1:2,maxblocksfl))
      allocate( & 
     & tflux_x(nfluxes,1:2, & 
     &        jl_bndi:ju_bndi,kl_bndi:ku_bndi,maxblocksfl))
      allocate( & 
     & tflux_y(nfluxes,il_bndi:iu_bndi, & 
     &         1:2,kl_bndi:ku_bndi,maxblocksfl))
      allocate( & 
     & tflux_z(nfluxes,il_bndi:iu_bndi, & 
     &        jl_bndi:ju_bndi,1:2,maxblocksfl))

      allocate( & 
     & bedge_facex_y(nedges,1:2,jl_bnd:ju_bnd+1, & 
     &               kl_bnd:ku_bnd+1,maxblockse))
      allocate( & 
     & bedge_facex_z(nedges,1:2,jl_bnd:ju_bnd+1, & 
     &               kl_bnd:ku_bnd+1,maxblockse))
      allocate( & 
     & bedge_facey_x(nedges,il_bnd:iu_bnd+1,1:2, & 
     &               kl_bnd:ku_bnd+1,maxblockse))
      allocate( & 
     & bedge_facey_z(nedges,il_bnd:iu_bnd+1,1:2, & 
     &               kl_bnd:ku_bnd+1,maxblockse))
      allocate( & 
     & bedge_facez_x(nedges,il_bnd:iu_bnd+1, & 
     &               jl_bnd:ju_bnd+1,1:2,maxblockse))
      allocate( & 
     & bedge_facez_y(nedges,il_bnd:iu_bnd+1, & 
     &               jl_bnd:ju_bnd+1, & 
     &               1:2,maxblockse))
      allocate( & 
     & recvarx1e(nedges,1:2,jl_bnd:ju_bnd+1, & 
     &           kl_bnd:ku_bnd+1))
      allocate( & 
     & recvary1e(nedges,il_bnd:iu_bnd+1,1:2, & 
     &           kl_bnd:ku_bnd+1))
      allocate( & 
     & recvarz1e(nedges,il_bnd:iu_bnd+1,jl_bnd:ju_bnd+1, & 
     &           1:2))
      allocate( & 
     & recvarx2e(nedges,1:2,jl_bnd:ju_bnd+1, & 
     &           kl_bnd:ku_bnd+1))
      allocate( & 
     & recvary2e(nedges,il_bnd:iu_bnd+1,1:2, & 
     &           kl_bnd:ku_bnd+1))
      allocate( & 
     & recvarz2e(nedges,il_bnd:iu_bnd+1,jl_bnd:ju_bnd+1, & 
     &           1:2))
      allocate( & 
     & tbedge_facex_y(nedges,1:2,jl_bnd:ju_bnd+1, & 
     &                kl_bnd:ku_bnd+1,maxblockse))
      allocate( & 
     & tbedge_facex_z(nedges,1:2,jl_bnd:ju_bnd+1, & 
     &                kl_bnd:ku_bnd+1,maxblockse))
      allocate( & 
     & tbedge_facey_x(nedges,il_bnd:iu_bnd+1,1:2, & 
     &                kl_bnd:ku_bnd+1,maxblockse))
      allocate( & 
     & tbedge_facey_z(nedges,il_bnd:iu_bnd+1,1:2, & 
     &                kl_bnd:ku_bnd+1,maxblockse))
      allocate( & 
     & tbedge_facez_x(nedges,il_bnd:iu_bnd+1, & 
     &                jl_bnd:ju_bnd+1,1:2,maxblockse))
      allocate( & 
     & tbedge_facez_y(nedges,il_bnd:iu_bnd+1, & 
     &                jl_bnd:ju_bnd+1, & 
     &                1:2,maxblockse))


      allocate(recvarxf(nfluxes,1:2,jl_bndi:ju_bndi,kl_bndi:ku_bndi))
      allocate(recvaryf(nfluxes,il_bndi:iu_bndi,1:2,kl_bndi:ku_bndi))
      allocate(recvarzf(nfluxes,il_bndi:iu_bndi,jl_bndi:ju_bndi,1:2))
      allocate(bndtempx1(nfluxes,1:2,jl_bndi:ju_bndi,kl_bndi:ku_bndi))
      allocate(bndtempy1(nfluxes,il_bndi:iu_bndi,1:2,kl_bndi:ku_bndi))
      allocate(bndtempz1(nfluxes,il_bndi:iu_bndi,jl_bndi:ju_bndi,1:2))

      len_block_bndx=2*(ju_bndi-jl_bndi+1)*(ku_bndi-kl_bndi+1)
      len_block_bndy=2*(iu_bndi-il_bndi+1)*(ku_bndi-kl_bndi+1)
      len_block_bndz=2*(iu_bndi-il_bndi+1)*(ju_bndi-jl_bndi+1)
      len_block_ex=2*(ju_bnd+k2d)*(ku_bnd+k3d)
      len_block_ey=2*(iu_bnd+1  )*(ku_bnd+k3d)
      len_block_ez=2*(iu_bnd+1  )*(ju_bnd+k2d)

! tree data

      maxblocks_tr=10*maxblocks

      allocate(neigh(2,mfaces,maxblocks_tr))
      allocate(child(2,mchild,maxblocks_tr))
      allocate(which_child(maxblocks_tr))
      allocate(parent(2,maxblocks_tr))
      allocate(lrefine(maxblocks_tr))
      allocate(nodetype(maxblocks_tr))
      allocate(empty(maxblocks_tr))
      allocate(bflags(mflags,maxblocks_tr))
      allocate(newchild(maxblocks_tr))
      allocate(derefine(maxblocks_tr))
      allocate(refine(maxblocks_tr))
      allocate(stay(maxblocks_tr))
      allocate(work_block(maxblocks_tr))
      allocate(coord(mdim,maxblocks_tr))
      allocate(bsize(mdim,maxblocks_tr))
      allocate(bnd_box(2,mdim,maxblocks_tr))
      allocate(level_cell_sizes(mdim,maxlevels))
      allocate(laddress(1:2,1:maxblocks_alloc))
      allocate(surr_blks(3,3,1+2*k2d,1+2*k3d,maxblocks_alloc))
#ifdef SAVE_MORTS
      allocate(surr_morts(6,3,1+2*k2d,1+2*k3d,maxblocks_alloc))
#endif
      allocate(boundary_box(2,mdim,mboundaries))
      allocate(boundary_index(mboundaries))

! workspace data

      allocate(work(ilw:iuw,jlw:juw,klw:kuw,maxblocks, & 
     &               nvar_work))
      allocate(interp_mask_work(nvar_work))
      allocate(interp_mask_work_res(nvar_work))
      allocate(recvw(ilw:iuw,jlw:juw,klw:kuw))
      allocate(sendw(ilw:iuw,jlw:juw,klw:kuw))
      allocate(tempw(ilw:iuw,jlw:juw,klw:kuw))
      allocate(work1(ilw1:iuw1,jlw1:juw1,klw1:kuw1,npblks))
      allocate(work1_fl(ilw1:iuw1+nxb+2*nguard_work, & 
     &                  jlw1:juw1+(nyb+2*nguard_work)*k2d, & 
     &                  klw1:kuw1+(nzb+2*nguard_work)*k3d) )
      allocate(recvw1(ilw1:iuw1,jlw1:juw1,klw1:kuw1,npblks))
      allocate(tempw1(ilw1:iuw1,jlw1:juw1,klw1:kuw1))

! morton data

      allocate(mortonbnd(6,1:3,1:maxblocks))
      allocate(laddress_guard(1:2,1:maxblocks_alloc))
      allocate(laddress_prol(1:2,1:maxblocks_alloc))
      allocate(laddress_flux(1:2,1:maxblocks_alloc))
      allocate(laddress_restrict(1:2,1:maxblocks_alloc))

! prolong_arrays data

      allocate(prol_dx(il_bnd1:iu_bnd1))
      allocate(prol_dy(jl_bnd1:ju_bnd1))
      allocate(prol_dz(kl_bnd1:ku_bnd1))
      allocate(prol_indexx(2,il_bnd1:iu_bnd1,2))
      allocate(prol_indexy(2,jl_bnd1:ju_bnd1,2))
      allocate(prol_indexz(2,kl_bnd1:ku_bnd1,2))
      allocate(prol_f_dx(il_bnd1:iu_bnd1+1))
      allocate(prol_f_dy(jl_bnd1:ju_bnd1+k2d))
      allocate(prol_f_dz(kl_bnd1:ku_bnd1+k3d))
      allocate(prol_f_indexx(2,il_bnd1:iu_bnd1+1,2))
      allocate(prol_f_indexy(2,jl_bnd1:ju_bnd1+k2d,2))
      allocate(prol_f_indexz(2,kl_bnd1:ku_bnd1+k3d,2))
      allocate(prolw_dx(ilw1:iuw1))
      allocate(prolw_dy(jlw1:juw1))
      allocate(prolw_dz(klw1:kuw1))
      allocate(prolw_indexx(2,ilw1:iuw1,2))
      allocate(prolw_indexy(2,jlw1:juw1,2))
      allocate(prolw_indexz(2,klw1:kuw1,2))

! an array for timings
      allocate(timer_amr_1blk_to_perm(0:1+nvar_work))

#endif /* LIBRARY */


      allocate( & 
     &  gt_unk(nvar, & 
     &         il_bnd:iu_bnd,jl_bnd:ju_bnd, & 
     &         kl_bnd:ku_bnd,maxblocks_gt))
      if (no_permanent_guardcells) then
      allocate( & 
     & gt_facevarx(nbndvar,il_bnd:iu_bnd+1,jl_bnd:ju_bnd, & 
     &             kl_bnd:ku_bnd,maxblocksf_gt))
      allocate( & 
     & gt_facevary(nbndvar,il_bnd:iu_bnd,jl_bnd:ju_bnd+k2d, & 
     &             kl_bnd:ku_bnd,maxblocksf_gt))
      allocate( & 
     & gt_facevarz(nbndvar,il_bnd:iu_bnd,jl_bnd:ju_bnd, & 
     &              kl_bnd:ku_bnd+k3d,maxblocksf_gt))
      else ! no_permanent_guardcells
      allocate( & 
     & gt_facevarx(nbndvar,1:2,jl_bnd:ju_bnd, & 
     &             kl_bnd:ku_bnd,maxblocksf))
      allocate( & 
     & gt_facevary(nbndvar,il_bnd:iu_bnd,1:1+k2d, & 
     &             kl_bnd:ku_bnd,maxblocksf))
      allocate( & 
     & gt_facevarz(nbndvar,il_bnd:iu_bnd,jl_bnd:ju_bnd, & 
     &             1:1+k3d,maxblocksf))
      endif
      allocate( & 
     & gt_unk_e_x(nbndvare,il_bnd:iu_bnd,jl_bnd:ju_bnd+k2d, & 
     &             kl_bnd:ku_bnd+k3d,maxblocksue_gt))
      allocate( & 
     & gt_unk_e_y(nbndvare,il_bnd:iu_bnd+1,jl_bnd:ju_bnd, & 
     &            kl_bnd:ku_bnd+k3d,maxblocksue_gt))
      allocate( & 
     & gt_unk_e_z(nbndvare,il_bnd:iu_bnd+1, & 
     &            jl_bnd:ju_bnd+k2d, & 
     &            kl_bnd:ku_bnd,maxblocksue_gt))
      allocate( & 
     & gt_unk_n(nbndvarc,il_bnd:iu_bnd+1,jl_bnd:ju_bnd+k2d, & 
     &          kl_bnd:ku_bnd+k3d,maxblocksn_gt))

      if (var_dt .or. pred_corr) then
      allocate( & 
     & t_unk(nvar, & 
     &       il_bnd:iu_bnd, & 
     &       jl_bnd:ju_bnd, & 
     &       kl_bnd:ku_bnd, & 
     &       maxblocks))
      allocate( & 
     & tfacevarx(nbndvar, & 
     &          il_bnd:iu_bnd+1, & 
     &          jl_bnd:ju_bnd, & 
     &          kl_bnd:ku_bnd, & 
     &          maxblocksf))
      allocate( & 
     & tfacevary(nbndvar, & 
     &          il_bnd:iu_bnd, & 
     &          jl_bnd:ju_bnd+k2d, & 
     &          kl_bnd:ku_bnd, & 
     &          maxblocksf))
      allocate( & 
     & tfacevarz(nbndvar, & 
     &          il_bnd:iu_bnd, & 
     &          jl_bnd:ju_bnd, & 
     &          kl_bnd:ku_bnd+k3d, & 
     &          maxblocksf))
      allocate( & 
     & t_unk_e_x(nbndvare, & 
     &          il_bnd:iu_bnd, & 
     &          jl_bnd:ju_bnd+k2d, & 
     &          kl_bnd:ku_bnd+k3d, & 
     &          maxblocksue))
      allocate( & 
     & t_unk_e_y(nbndvare, & 
     &           il_bnd:iu_bnd+1, & 
     &           jl_bnd:ju_bnd, & 
     &           kl_bnd:ku_bnd+k3d, & 
     &           maxblocksue))
      allocate( & 
     & t_unk_e_z(nbndvare, & 
     &          il_bnd:iu_bnd+1, & 
     &          jl_bnd:ju_bnd+k2d, & 
     &          kl_bnd:ku_bnd, & 
     &          maxblocksue))
      allocate( & 
     & t_unk_n(nbndvarc, & 
     &         il_bnd:iu_bnd+1,jl_bnd:ju_bnd+k2d, & 
     &         kl_bnd:ku_bnd+k3d, & 
     &         maxblocksn))
      endif ! var_dt

      if (var_dt) then
      allocate( & 
     & ttflux_x(nfluxes,1:2,jl_bnd:ju_bnd, & 
     &          kl_bnd:ku_bnd,maxblocksfl))
      allocate( & 
     & ttflux_y(nfluxes,il_bnd:iu_bnd, & 
     &          1:2,kl_bnd:ku_bnd,maxblocksfl))
      allocate( & 
     & ttflux_z(nfluxes,il_bnd:iu_bnd, & 
     &          jl_bnd:ju_bnd,1:2,maxblocksfl))
      endif

      if (var_dt) then
      allocate( & 
     & ttbedge_facex_y(nedges,1:2,jl_bnd:ju_bnd+1, & 
     &                 kl_bnd:ku_bnd+1,maxblockse))
      allocate( & 
     & ttbedge_facex_z(nedges,1:2,jl_bnd:ju_bnd+1, & 
     &                 kl_bnd:ku_bnd+1,maxblockse))
      allocate( & 
     & ttbedge_facey_x(nedges,il_bnd:iu_bnd+1,1:2, & 
     &                 kl_bnd:ku_bnd+1,maxblockse))
      allocate( & 
     & ttbedge_facey_z(nedges,il_bnd:iu_bnd+1,1:2, & 
     &                 kl_bnd:ku_bnd+1,maxblockse))
      allocate( & 
     & ttbedge_facez_x(nedges,il_bnd:iu_bnd+1, & 
     &                jl_bnd:ju_bnd+1,1:2,maxblockse))
      allocate( & 
     & ttbedge_facez_y(nedges,il_bnd:iu_bnd+1, & 
     &                jl_bnd:ju_bnd+1, & 
     &                1:2,maxblockse))
      end if


      if (curvilinear) then
      allocate( & 
     & cell_vol(il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1,kl_bnd1:ku_bnd1))
      allocate(cell_area1(il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1, & 
     &                    kl_bnd1:ku_bnd1))
      allocate(cell_area2(il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1+k2d, & 
     &                    kl_bnd1:ku_bnd1))
      allocate(cell_area3(il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1, & 
     &                    kl_bnd1:ku_bnd1+k3d))
      allocate(cell_leng1(il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1+k2d, & 
     &                    kl_bnd1:ku_bnd1+k3d))
      allocate(cell_leng2(il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1, & 
     &                    kl_bnd1:ku_bnd1+k3d))
      allocate(cell_leng3(il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d, & 
     &                    kl_bnd1:ku_bnd1))
      allocate(cell_face_coord1(il_bnd1:iu_bnd1+1))
      allocate(cell_face_coord2(jl_bnd1:ju_bnd1+k2d))
      allocate(cell_face_coord3(kl_bnd1:ku_bnd1+k3d))
      allocate(cell_vol_w(ilw1:iuw1,jlw1:juw1,klw1:kuw1))
      end if

      allocate(ladd_strt(0:nprocs-1))
      allocate(ladd_end(0:nprocs-1))

! initialize tree data structure
        bsize(:,:) = -1.
        lrefine(:) = -1
        nodetype(:) = -1
        stay(:) = .TRUE.
        refine(:) = .FALSE.
        derefine(:) = .FALSE.
        parent(:,:) = -1
        child(:,:,:) = -1
        which_child(:) = -1
        coord(:,:) = -1.
        bnd_box(:,:,:) = -1.
        neigh(:,:,:) = -1
        empty(:) = 0
        bflags(:,:) = -1
        work_block(:) = 0.
        surr_blks(:,:,:,:,:) = -1
#ifdef SAVE_MORTS
        surr_morts(:,:,:,:,:) = -1
#endif

! initialize solution arrays
        unk(:,:,:,:,:) = 0.
        facevarx(:,:,:,:,:) = 0.
        facevary(:,:,:,:,:) = 0.
        facevarz(:,:,:,:,:) = 0.
        unk_e_x(:,:,:,:,:) = 0.
        unk_e_y(:,:,:,:,:) = 0.
        unk_e_z(:,:,:,:,:) = 0.
        unk_n(:,:,:,:,:) = 0.

! initialize boundary location arrays for mpi use.
        boundary_box(:,:,:nboundaries) = 0.
        boundary_index(:nboundaries) = -1


! Initialization required for prolongation routines
      call amr_prolong_fun_init

! Set default values for gcell logical control arrays
      do i = 1,nvar
        gcell_on_cc_pointer(i) = i
      enddo
      gcell_on_cc(:)     = .true.
      int_gcell_on_cc(:) = .true.

      do i = 1,nfacevar
        gcell_on_fc_pointer(:,i) = i
      enddo
      gcell_on_fc(:,:)     = .true.
      int_gcell_on_fc(:,:) = .true.

      do i = 1,nedgevar
        gcell_on_ec_pointer(:,i) = i
      enddo
      gcell_on_ec(:,:)     = .true.
      int_gcell_on_ec(:,:) = .true.

      do i = 1,nvarcorn
        gcell_on_nc_pointer(i) = i
      enddo
      gcell_on_nc(:)     = .true.
      int_gcell_on_nc(:) = .true.

      checkp_on_cc = .true.
      checkp_on_fc = .true.
      checkp_on_ec = .true.
      checkp_on_nc = .true.

! Set default values for interp_masks

      interp_mask_unk(:) = 1
      interp_mask_work(:) = 1
      interp_mask_facex(:) = 1
      interp_mask_facey(:) = 1
      interp_mask_facez(:) = 1
      interp_mask_ec(:) = 1
      interp_mask_nc(:) = 1

      interp_mask_unk_res(:) = 1
      interp_mask_work_res(:) = 1
      interp_mask_facex_res(:) = 1
      interp_mask_facey_res(:) = 1
      interp_mask_facez_res(:) = 1
      interp_mask_ec_res(:) = 1
      interp_mask_nc_res(:) = 1

! Initialize index array defining the variables which
! constitute any divergence free fields
      allocate(i_divf_fc_vars(3,nfield_divf))
      do nfield = 1,nfield_divf
        i_divf_fc_vars(:,nfield) = nfield
      enddo

! Initialization required for boundary condition routine
      call amr_bcset_init


      call amr_1blk_guardcell_reset

! Mark amr_gsurrounding_blks uncalled. This will be set to +1 if
! and when amr_gsurrounding_blks is called.
      gsurrblks_set = -1

! Initialize grid change marker flag
! Initial value = 1 reflects a changed grid.
      grid_changed = 1

! Initialize flag to detect if amr_checkpoint_re or amr_refine_derefine
! have been called. grid_analysed_mpi = +1 means tha one of them has, so
! the mpi version will have control info for any communication dependent
! routines.
      grid_analysed_mpi = -1

!
! Initialize mpi communication pattern id
      mpi_pattern_id = 0

! set state flags used by mpi communications
      lrestrict_in_progress = .false.
      lprolong_in_progress  = .false.
      lguard_in_progress    = .false.

! set state flags used by mpi block boundary info list
! -1 is unset. When set, it will have value 100.
      bc_block_neighs_status = -1

!
! Initialize some arrays used in controling mpi communications
      commatrix_recv(:) = 0
      commatrix_send(:) = 0

! initialize the instance counter
      instance = 0

! initialize no permanent_guardcells variables

      Call MPI_BARRIER(amr_mpi_meshComm,ierr)

      if (timing_mpi) then
         timer_amr_initialize =  timer_amr_initialize & 
     &                         + mpi_wtime() - time1
      endif

      return
      end subroutine amr_initialize
