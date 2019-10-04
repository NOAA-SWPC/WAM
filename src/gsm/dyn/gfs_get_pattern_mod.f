module gfs_get_pattern_mod
 use gfs_dyn_vert_def
 use gfs_dyn_resol_def
 use gfs_dyn_layout1
 use gfs_dyn_gg_def
 use gfs_dyn_date_def
 use namelist_dynamics_def
 use gfs_dyn_coordinate_def
 use gfs_dyn_stoch_data
 use gfs_dyn_patterngenerator
 use gfs_dyn_mpi_def, only: mc_comp,mpi_sum,mpi_real4,mpi_complex,mpi_r_io_r,mpi_real8
 use gfs_dyn_machine, only: r_kind => kind_io4, kind_evod
 use gfs_dyn_physcons, only: rerth=>con_rerth
 use gfs_dyn_tracer_const
 implicit none
 private

 public  get_pattern_sppt,get_pattern_shum,get_pattern_skeb,get_pattern_vc
 public dump_patterns,restore_patterns
 contains

subroutine get_pattern_sppt(&
           ls_node,ls_nodes,max_ls_nodes,&
           lats_nodes_a,global_lats_a,lonsperlar,&
           plnev_a,plnod_a,&
           sppt3d,dt)
! generate random pattern for SPPT.
! output array sppt3d contains pattern for latitudes on this task.
 implicit none

 integer, intent(in) :: ls_node(ls_dim,3),ls_nodes(ls_dim,nodes),&
   max_ls_nodes(nodes),lats_nodes_a(nodes),global_lats_a(latg),lonsperlar(latg)
 real(kind=kind_evod),intent(in) :: &
    plnev_a(len_trie_ls,latg2),plnod_a(len_trio_ls,latg2)

 real(kind=kind_phys),intent(out) :: sppt3d(lonf,lats_node_a_max,levs)
 integer i,j,k,l,lat,ierr,n

 real(kind_phys) dt
 real(kind_evod), dimension(lonf,lats_node_a):: sppt2d,wrk2d
 integer :: num2d
! logical lprint

 !real(r_kind), allocatable, dimension(:,:,:) :: workg,workg_out
 !real (kind=kind_io8)   glolal(lonf,lats_node_a)
 !real (kind=r_kind)   wrkga(lonf*latg)
 !integer kmsk0(lonf,lats_node_a)
 !kmsk0 = 0

 if (.not. allocated (spec_sppt_e)) then
    call init_stochdata(dt,ls_node)
 endif
 sppt2d = 0.
 do n=1,nsppt
 call patterngenerator_advance(spec_sppt_e(:,:,n),spec_sppt_o(:,:,n),rpattern_sppt(n))
! print *,'min/max spec_sppt_e',minval(spec_sppt_e),maxval(spec_sppt_e)
 call scalarspect_to_grid(&
           spec_sppt_e(:,:,n),spec_sppt_o(:,:,n),wrk2d,&
           ls_node,ls_nodes,max_ls_nodes,&
           lats_nodes_a,global_lats_a,lonsperlar,&
           plnev_a,plnod_a,1)
 sppt2d = sppt2d + wrk2d
 enddo
 !print *,'min/max sppt2d',minval(sppt2d),maxval(sppt2d)
 if (sppt_logit) sppt2d = (2./(1.+exp(sppt2d)))-1.
 !print *,'min/max sppt2d',minval(sppt2d),maxval(sppt2d)
 do k = 1, levs
    sppt3d(:,1:lats_node_a,k)= sppt2d*vfact_sppt(k)+ 1.0
 enddo
! write out data
 !allocate(workg(lonf,latg,levs))
 !allocate(workg_out(lonf,latg,levs))
 !workg = 0.
 !do k=1,levs
 !  CALL uninterpred(2,kmsk0,glolal,sppt3d(:,:,k),&
 !                   global_lats_a,lonsperlar)
 !  do j=1,lats_node_a
 !     lat=global_lats_a(ipt_lats_node_a-1+j)
 !     do i=1,lonf
 !        workg(i,lat,k) = glolal(i,j)
 !     enddo
 !  enddo
 !enddo
 !call mpi_reduce(workg,workg_out,lonf*latg*levs,&
 !                mpi_real4,mpi_sum,me_l_0,mc_comp,ierr)
 !if (me .eq. me_l_0) then
 !   open(77,form='unformatted',access='direct',recl=lonf*latg*levs)
 !   write(77,rec=1) workg_out
 !   close(77)
 !endif
 !deallocate(workg,workg_out)
 !call mpi_barrier(mc_comp,ierr)
 !call mpi_quit(9999)

end subroutine get_pattern_sppt

subroutine get_pattern_shum(&
           ls_node,ls_nodes,max_ls_nodes,&
           lats_nodes_a,global_lats_a,lonsperlar,&
           plnev_a,plnod_a,&
           shum3d,dt)

! generate random pattern for SHUM.
! output array shum3d contains pattern for latitudes on this task.

 implicit none

 integer, intent(in) :: ls_node(ls_dim,3),ls_nodes(ls_dim,nodes),&
   max_ls_nodes(nodes),lats_nodes_a(nodes),global_lats_a(latg),lonsperlar(latg)
 real(kind=kind_evod),intent(in) :: &
    plnev_a(len_trie_ls,latg2),plnod_a(len_trio_ls,latg2)

 real(kind=kind_phys),intent(out) :: shum3d(lonf,lats_node_a_max,levs)
    
 integer i,j,k,l,lat,ierr,n

 real(kind_phys) dt
 real(kind_evod), dimension(lonf,lats_node_a):: shum2d
 integer :: num2d

 !real(r_kind), allocatable, dimension(:,:,:) :: workg,workg_out
 !real (kind=kind_io8)   glolal(lonf,lats_node_a)
 !real (kind=r_kind)   wrkga(lonf*latg)
 !integer kmsk0(lonf,lats_node_a)
 !kmsk0 = 0

 if (.not. allocated (spec_shum_e)) then
    call init_stochdata(dt,ls_node)
 endif
 shum3d(:,:,:)=0.0
 do n=1,nshum
 call patterngenerator_advance(spec_shum_e(:,:,n),spec_shum_o(:,:,n),rpattern_shum(n))
 call scalarspect_to_grid(&
           spec_shum_e(:,:,n),spec_shum_o(:,:,n),shum2d,&
           ls_node,ls_nodes,max_ls_nodes,&
           lats_nodes_a,global_lats_a,lonsperlar,&
           plnev_a,plnod_a,1)
 do k = 1, levs
    shum3d(:,1:lats_node_a,k)=shum3d(:,1:lats_node_a,k)+shum2d*vfact_shum(k)
 enddo
 enddo
! write out data
!allocate(workg(lonf,latg,levs))
!allocate(workg_out(lonf,latg,levs))
!workg = 0.; workg_out = 0.
!do k=1,levs
!  CALL uninterpred(2,kmsk0,glolal,shum3d(:,:,k),&
!                   global_lats_a,lonsperlar)
!  do j=1,lats_node_a
!     lat=global_lats_a(ipt_lats_node_a-1+j)
!     do i=1,lonf
!        workg(i,lat,k) = glolal(i,j)
!     enddo
!  enddo
!enddo
!call mpi_reduce(workg,workg_out,lonf*latg*levs,&
!                mpi_real4,mpi_sum,me_l_0,mc_comp,ierr)
!if (me .eq. me_l_0) then
!   print *,'min/max shum out',minval(workg_out),maxval(workg_out)
!   open(77,form='unformatted',access='direct',recl=lonf*latg*levs)
!   write(77,rec=1) workg_out
!   close(77)
!endif
!deallocate(workg,workg_out)
!call mpi_barrier(mc_comp,ierr)
!call mpi_quit(9999)

end subroutine get_pattern_shum

subroutine get_pattern_skeb(vrtspec_e,&
                            vrtspec_o,&
                            divspec_e,&
                            divspec_o,&
           ls_node,ls_nodes,max_ls_nodes,&
           lats_nodes_a,global_lats_a,lonsperlar,&
           epsedn,epsodn,snnp1ev,snnp1od,&
           plnev_a,plnod_a,plnew_a,plnow_a,&
                            skeb3d_u,&
                            skeb3d_v,dt)

! generate random patterns for SKEB (stochastic kinetic energy backscatter).
! output arrays skeb3d_u,skeb3d_v contains u and v patterns for latitudes on this task.

 implicit none

 real(kind=kind_evod), intent(in), dimension(len_trie_ls,2,levs) :: &
  vrtspec_e,divspec_e
 real(kind=kind_evod), intent(in), dimension(len_trio_ls,2,levs) :: &
  vrtspec_o,divspec_o
 integer, intent(in) :: ls_node(ls_dim,3),ls_nodes(ls_dim,nodes),&
   max_ls_nodes(nodes),lats_nodes_a(nodes),global_lats_a(latg),lonsperlar(latg)
 real(kind=kind_evod),intent(in) ::  epsedn(len_trie_ls),&
  epsodn(len_trio_ls),snnp1ev(len_trie_ls),snnp1od(len_trio_ls),&
  plnev_a(len_trie_ls,latg2),plnod_a(len_trio_ls,latg2),&
  plnew_a(len_trie_ls,latg2),plnow_a(len_trio_ls,latg2)
 real(kind=kind_phys),intent(out) :: &
   skeb3d_u(lonf,lats_node_a_max,levs),&
   skeb3d_v(lonf,lats_node_a_max,levs)
 real(kind_phys), intent(in) :: dt
 ! locals
 real(kind=kind_phys),dimension(lonf,lats_node_a,levs) ::&
  udiffg,vdiffg,dissrate,ug,vg
 real(kind=kind_evod), dimension(len_trie_ls,2,levs) :: &
  vrtdiffspec_e,divdiffspec_e
 real(kind=kind_evod), dimension(len_trio_ls,2,levs) :: &
  vrtdiffspec_o,divdiffspec_o
 real(kind=kind_evod), dimension(len_trie_ls) :: &
  smoothfact_e,kenorm_e,wavenumsq_e
 real(kind=kind_evod), dimension(len_trio_ls) :: &
  smoothfact_o,kenorm_o,wavenumsq_o
 complex(r_kind), dimension((jcap+1)*(jcap+2)/2) :: workspec
 integer i,j,k,l,n,nn,locl,indev,indod,ierr,jbasev,jbasod,indlsod,indlsev,lat
 real(kind=kind_evod) :: globalvar,globalvar0
 include 'function_indlsod'
 include 'function_indlsev'

! real(r_kind) t0,t1,t2

 real(r_kind) rnn1,rnn0,epstiny,rnnmax
 logical lprint

 !real(r_kind), allocatable, dimension(:,:,:) :: workg,workg_out
 !integer kmsk0(lonf,lats_node_a)
 !real (kind=kind_io8)   glolal(lonf,lats_node_a)
 !real (kind=r_kind)   wrkga(lonf*latg)
 !kmsk0 = 0

 lprint = .false.

 epstiny = tiny(rnn1)

 rnn1 = skeb_diss_smooth ! use namelist parameter to define smoothing scale.
 rnn0 = rnn1*(rnn1+1.)
 smoothfact_e=1.; smoothfact_o=1. ! used to smooth dissipation estimate.
 kenorm_e=0.; kenorm_o=0. ! used to convert forcing pattern to wind field.
 do locl=1,ls_max_node
     l = ls_node(locl,1)
     jbasev = ls_node(locl,2)
     indev = indlsev(l,l)
     jbasod = ls_node(locl,3)
     indod = indlsod(l+1,l)
     do n=l,jcap,2
        rnn1 = n*(n+1.)
        smoothfact_e(indev) = exp(-(rnn1/rnn0))
        kenorm_e(indev) = sqrt(rnn1)/rerth
        indev = indev + 1
     enddo
     do n=l+1,jcap,2
        rnn1 = n*(n+1.)
        smoothfact_o(indod) = exp(-(rnn1/rnn0))
        kenorm_o(indod) = sqrt(rnn1)/rerth
        indod = indod + 1
     enddo
  enddo
  ! set the even and odd (n-l) terms of the top row to zero
  do locl=1,ls_max_node
     l = ls_node(locl,1)
     jbasev = ls_node(locl,2)
     jbasod = ls_node(locl,3)
     if (mod(l,2) .EQ. mod(jcap+1,2)) then
        smoothfact_e(indlsev(jcap+1,l)) = 0.
        kenorm_e(indlsev(jcap+1,l)) = 0.
     endif
     if (mod(l,2) .NE. mod(jcap+1,2)) then
        smoothfact_o(indlsod(jcap+1,l)) = 0.
        kenorm_o(indlsod(jcap+1,l)) = 0.
     endif
  enddo
  wavenumsq_e = ((kenorm_e*rerth)**2) ! n*(n+1)
  wavenumsq_o = ((kenorm_o*rerth)**2) 
! streamfunction norm (default is KE norm)
! perturbations are smaller scale
! kenorm_e = rerth*kenorm_e**2; kenorm_o = rerth*kenorm_o**2
! vorticity norm 
! perturbations are larger scale
! kenorm_e = 1./rerth; kenorm_o = 1./rerth

! compute vorticity gradient.
! numerical diffusion assumed proportional to magnitude of vorticity
! gradient, as in ECMWF specral stochastic backscatter implementation.

! fill divdiffspec with laplacian of vorticity
 do k=1,levs
   do n=1,2
     divdiffspec_e(:,n,k) = -wavenumsq_e*vrtspec_e(:,n,k)
     divdiffspec_o(:,n,k) = -wavenumsq_o*vrtspec_o(:,n,k)
   enddo
 enddo
 vrtdiffspec_e=0; vrtdiffspec_o=0.
! on return udiffg, vdiffg are vorticity gradient.
! (laplacian of vorticity takes place of divergence when computing
!  u and v, vorticity is zero.  In other words, vorticity
!  mimics velocity potential, divergent wind is gradient of vel. potential)
 call vrtdivspect_to_uvgrid(&
           divdiffspec_e,divdiffspec_o,vrtdiffspec_e,vrtdiffspec_o,&
           udiffg,vdiffg,&
           ls_node,ls_nodes,max_ls_nodes,&
           lats_nodes_a,global_lats_a,lonsperlar,&
           epsedn,epsodn,snnp1ev,snnp1od,plnev_a,plnod_a,levs)

! vorticity gradient divided by earth radius 
! gradient modulus squared, skeb coefficient O(1000)
! (ECMWF uses this)
 udiffg = (udiffg**2+vdiffg**2)/rerth**2
! RMS gradient skeb coefficient is O(10)-O(100)
 !udiffg = sqrt(udiffg**2+vdiffg**2)/rerth
 !if (me .eq. me_l_0) print *,'min/max dissrate',&
 ! minval(dissrate),maxval(dissrate)

! smooth the dissipation estimate.
! back to spectral space.
 call scalargrid_to_spect(&
            divdiffspec_e,divdiffspec_o,udiffg,&
            ls_node,ls_nodes,max_ls_nodes,&
            lats_nodes_a,global_lats_a,lonsperlar,&
            plnew_a,plnow_a,levs)
! smooth in spectral space.
 do k=1,levs
   do n=1,2
     divdiffspec_e(:,n,k) = smoothfact_e*divdiffspec_e(:,n,k)
     divdiffspec_o(:,n,k) = smoothfact_o*divdiffspec_o(:,n,k)
   enddo
 enddo
! back to grid
 call scalarspect_to_grid(&
           divdiffspec_e,divdiffspec_o,dissrate,&
           ls_node,ls_nodes,max_ls_nodes,&
           lats_nodes_a,global_lats_a,lonsperlar,&
           plnev_a,plnod_a,levs)

! t2 = mpi_wtime()
! if (me .eq. me_l_0) print *,'time to compute dissrate = ',t2-t1
! t1 = mpi_wtime()

 if (.not. allocated(spec_skeb_e)) then
   call init_stochdata(dt,ls_node)
 endif
 ! generate random streamfunction forcing patterns.
 skeb3d_u=0; skeb3d_v=0.
 do n=1,nskeb
 do k=1,levs
   call patterngenerator_advance(spec_skeb_e(:,:,k,n),spec_skeb_o(:,:,k,n),&
                                 rpattern_skeb(n))
 enddo
 ! don't modify spectral arrays used to evolve pattern
 vrtdiffspec_e = spec_skeb_e(:,:,:,n)
 vrtdiffspec_o = spec_skeb_o(:,:,:,n)

! apply successive applications of 1-2-1 filter in vertical to introduce vertical correlations.
 if (skeb_vfilt(n) > 0) then

!   if (me .eq. me_l_0 .and. lprint) print *,'applying 1-2-1 filter',skeb_vfilt,'times'

   do nn=1,skeb_vfilt(n)
      do k=2,levs-1
         divdiffspec_e(:,:,k) = vrtdiffspec_e(:,:,k+1)+ &
                        2.*vrtdiffspec_e(:,:,k)+&
                        vrtdiffspec_e(:,:,k-1)
         divdiffspec_o(:,:,k) = vrtdiffspec_o(:,:,k+1)+ &
                        2.*vrtdiffspec_o(:,:,k)+&
                        vrtdiffspec_o(:,:,k-1)
      enddo
      divdiffspec_e(:,:,1) =  &
      (1.+1./3.)*vrtdiffspec_e(:,:,2)+&
      2.*(1.+1./3.)*vrtdiffspec_e(:,:,1)
      divdiffspec_e(:,:,levs) = &
      (1.+1./3.)*vrtdiffspec_e(:,:,levs-1)+&
      2.*(1.+1./3.)*vrtdiffspec_e(:,:,levs)
      divdiffspec_o(:,:,1) =  &
      (1.+1./3.)*vrtdiffspec_o(:,:,2)+&
      2.*(1.+1./3.)*vrtdiffspec_o(:,:,1)
      divdiffspec_o(:,:,levs) = &
      (1.+1./3.)*vrtdiffspec_o(:,:,levs-1)+&
      2.*(1.+1./3.)*vrtdiffspec_o(:,:,levs)
      vrtdiffspec_e = 0.25*divdiffspec_e
      vrtdiffspec_o = 0.25*divdiffspec_o
   enddo
   ! inflate variance of random patterns back to pre-vertical filtered values
   do k=1,levs
      ! compute variance at each level before smoothing
      call computevarspec_eo(spec_skeb_e(:,:,k,n),spec_skeb_o(:,:,k,n),globalvar0,ls_node)
!      if ( me .EQ. me_l_0) print*,'gvar before smoothing',n,k,sqrt(globalvar0)
      ! compute variance at each level after smoothing
      call computevarspec_eo(vrtdiffspec_e(:,:,k),vrtdiffspec_o(:,:,k),globalvar,ls_node)
!      if ( me .EQ. me_l_0) print*,'gvar after smoothing',n,k,sqrt(globalvar)
      ! normalize back to original variance
      vrtdiffspec_e(:,:,k)=vrtdiffspec_e(:,:,k)*sqrt(globalvar0)/sqrt(globalvar)
      vrtdiffspec_o(:,:,k)=vrtdiffspec_o(:,:,k)*sqrt(globalvar0)/sqrt(globalvar)
      !call
      !computevarspec_eo(vrtdiffspec_e(:,:,k),vrtdiffspec_o(:,:,k),globalvar,ls_node)
      !! check to see if new variance make sense
      !if ( me .EQ. me_l_0) print*,'gvar after adjustment',n,k,sqrt(globalvar)
    enddo

 end if

! t2 = mpi_wtime()
! if (me .eq. me_l_0) print *,'time for vertical smoothing',me,t2-t1
!
! t1 = mpi_wtime()

! ke norm (convert streamfunction forcing to vorticity forcing)
 divdiffspec_e = 0; divdiffspec_o = 0.
 do k=1,levs
   do nn=1,2
     vrtdiffspec_e(:,nn,k) = kenorm_e*vrtdiffspec_e(:,nn,k)*vfact_skeb(k)
     vrtdiffspec_o(:,nn,k) = kenorm_o*vrtdiffspec_o(:,nn,k)*vfact_skeb(k)
   enddo
 enddo
 ! modulate u and v forcing by smoothed dissipation rate.
 call vrtdivspect_to_uvgrid(&
           divdiffspec_e,divdiffspec_o,vrtdiffspec_e,vrtdiffspec_o,&
           udiffg,vdiffg,&
           ls_node,ls_nodes,max_ls_nodes,&
           lats_nodes_a,global_lats_a,lonsperlar,&
           epsedn,epsodn,snnp1ev,snnp1od,plnev_a,plnod_a,levs)
 skeb3d_u(:,1:lats_node_a,:) = skeb3d_u(:,1:lats_node_a,:) + udiffg*dissrate
 skeb3d_v(:,1:lats_node_a,:) = skeb3d_v(:,1:lats_node_a,:) + vdiffg*dissrate
 enddo

!t2 = mpi_wtime()
!if (me .eq. me_l_0 .and. lprint) then
!  print *,'min/max dissrate',minval(dissrate),maxval(dissrate)
!  print *,'min/max stoch u forcing',minval(udiffg),maxval(udiffg)
!  print *,'min/max total stoch forcing',&
!      minval(skeb3d_u),maxval(skeb3d_u)

!  print *,'time to compute tendencies',me,t2-t1
!  print *,'total time',me,t2-t0
!endif

! write out dissipation estimate if random patterns are being written out.
! (for diagnostic purposes)
!allocate(workg(lonf,latg,levs))
!allocate(workg_out(lonf,latg,levs))
!if (me .eq. me_l_0) then
!   open(77,file='skeb.bin',form='unformatted')
!endif
!workg = 0.
!workg_out = 0.
!do k=1,levs
!  CALL uninterpred(2,kmsk0,glolal,skeb3d_v(:,:,k),&
!                   global_lats_a,lonsperlar)
!  do j=1,lats_node_a
!     lat=global_lats_a(ipt_lats_node_a-1+j)
!     do i=1,lonf
!        workg(i,lat,k) = glolal(i,j)
!     enddo
!  enddo
!enddo
!call mpi_reduce(workg,workg_out,lonf*latg*levs,&
!                mpi_real4,mpi_sum,me_l_0,mc_comp,ierr)
!if (me .eq. me_l_0) then
!!   print *,'min/max dissrate',minval(dissrate),maxval(dissrate)
!!   print *,'min/max stoch u forcing',minval(udiffg),maxval(udiffg)
!!   print *,'min/max total stoch forcing',&
!!      minval(skeb3d_u),maxval(skeb3d_u)
!!
!   do k=1,levs
!      write(77) workg_out(:,:,k)
!   enddo 
!endif
!
!if (me .eq. me_l_0) then
!   close(77)
!endif
!deallocate(workg,workg_out)
!call mpi_barrier(mc_comp,ierr)
!call mpi_quit(9999)

end subroutine get_pattern_skeb

subroutine get_pattern_vc(vrtspec_e,&
                            vrtspec_o,&
                            divspec_e,&
                            divspec_o,&
           ls_node,ls_nodes,max_ls_nodes,&
           lats_nodes_a,global_lats_a,lonsperlar,&
           epsedn,epsodn,snnp1ev,snnp1od,&
           plnev_a,plnod_a,plnew_a,plnow_a,&
                            vc3d_u,&
                            vc3d_v,dt)

! generate perturbations for vorticity confinment 
! output arrays vc3d_u,vc3d_v contains u and v patterns for latitudes on this task.

 implicit none

 real(kind=kind_evod), intent(in), dimension(len_trie_ls,2,levs) :: &
  vrtspec_e,divspec_e
 real(kind=kind_evod), intent(in), dimension(len_trio_ls,2,levs) :: &
  vrtspec_o,divspec_o
 integer, intent(in) :: ls_node(ls_dim,3),ls_nodes(ls_dim,nodes),&
   max_ls_nodes(nodes),lats_nodes_a(nodes),global_lats_a(latg),lonsperlar(latg)
 real(kind=kind_evod),intent(in) ::  epsedn(len_trie_ls),&
  epsodn(len_trio_ls),snnp1ev(len_trie_ls),snnp1od(len_trio_ls),&
  plnev_a(len_trie_ls,latg2),plnod_a(len_trio_ls,latg2),&
  plnew_a(len_trie_ls,latg2),plnow_a(len_trio_ls,latg2)
 real(kind=kind_phys),intent(out) :: &
   vc3d_u(lonf,lats_node_a_max,levs),&
   vc3d_v(lonf,lats_node_a_max,levs)
 real(kind_phys), intent(in) :: dt
 ! locals
 real(kind=kind_phys),dimension(lonf,lats_node_a,levs) ::&
  vrtgradx,vrtgrady,vrtg,vrtgradmod
 real(kind=kind_evod), dimension(len_trie_ls,2,levs) :: &
  vrtdiffspec_e,divdiffspec_e
 real(kind=kind_evod), dimension(len_trio_ls,2,levs) :: &
  vrtdiffspec_o,divdiffspec_o
 real(kind=kind_evod), dimension(len_trie_ls) :: &
  wavenumsq_e
 real(kind=kind_evod), dimension(len_trio_ls) :: &
  wavenumsq_o
 real(kind_evod), dimension(lonf,lats_node_a):: vcfact,wrk2d
 integer i,j,k,l,n,locl,indev,indod,ierr,jbasev,jbasod,indlsod,indlsev,lat
 include 'function_indlsod'
 include 'function_indlsev'

 real(r_kind) rnn1
 real(r_kind) si(levs+1),sl
 if (.not. allocated(vfact_vc)) allocate(vfact_vc(levs))
   do k=1,levs+1
      si(levs+2-k)= ak5(k)/101.3+bk5(k) ! si are now sigmas
   enddo
   do k=1,levs
     sl = 0.5*(si(k)+si(k+1))
     if (sl .lt. vc_sigtop1 .and. sl .gt. vc_sigtop2) then
       vfact_vc(k) = (sl-vc_sigtop2)/(vc_sigtop1-vc_sigtop2)
     else if (sl .lt. vc_sigtop2) then
       vfact_vc(k) = 0.0
     else
       vfact_vc(k) = 1.0
     endif
  enddo

 wavenumsq_e=0.; wavenumsq_o=0. ! used to convert forcing pattern to wind field.
 do locl=1,ls_max_node
     l = ls_node(locl,1)
     jbasev = ls_node(locl,2)
     indev = indlsev(l,l)
     jbasod = ls_node(locl,3)
     indod = indlsod(l+1,l)
     do n=l,jcap,2
        rnn1 = n*(n+1.)
        wavenumsq_e(indev) = rnn1
        indev = indev + 1
     enddo
     do n=l+1,jcap,2
        rnn1 = n*(n+1.)
        wavenumsq_o(indod) = rnn1
        indod = indod + 1
     enddo
  enddo
  ! set the even and odd (n-l) terms of the top row to zero
  do locl=1,ls_max_node
     l = ls_node(locl,1)
     jbasev = ls_node(locl,2)
     jbasod = ls_node(locl,3)
     if (mod(l,2) .EQ. mod(jcap+1,2)) then
        wavenumsq_e(indlsev(jcap+1,l)) = 0.
     endif
     if (mod(l,2) .NE. mod(jcap+1,2)) then
       wavenumsq_o(indlsod(jcap+1,l)) = 0.
     endif
  enddo

! compute vorticity gradient.

! fill divdiffspec with laplacian of vorticity
 do k=1,levs
   do n=1,2
     divdiffspec_e(:,n,k) = -wavenumsq_e*vrtspec_e(:,n,k)
     divdiffspec_o(:,n,k) = -wavenumsq_o*vrtspec_o(:,n,k)
   enddo
 enddo
 vrtdiffspec_e=0; vrtdiffspec_o=0.
! on return udiffg, vdiffg are vorticity gradient.
! (laplacian of vorticity takes place of divergence when computing
!  u and v, vorticity is zero.  In other words, vorticity
!  mimics velocity potential, divergent wind is gradient of vel. potential)
 call vrtdivspect_to_uvgrid(&
           divdiffspec_e,divdiffspec_o,vrtdiffspec_e,vrtdiffspec_o,&
           vrtgradx,vrtgrady,&
           ls_node,ls_nodes,max_ls_nodes,&
           lats_nodes_a,global_lats_a,lonsperlar,&
           epsedn,epsodn,snnp1ev,snnp1od,plnev_a,plnod_a,levs)
 call scalarspect_to_grid(&
           vrtspec_e,vrtspec_o,vrtg,&
           ls_node,ls_nodes,max_ls_nodes,&
           lats_nodes_a,global_lats_a,lonsperlar,&
           plnev_a,plnod_a,levs)
 vrtgradmod = sqrt(vrtgradx**2+vrtgrady**2) + 1.e-12
 if (.not. allocated (spec_vc_e) .and. vcamp(1) > 0.) then
    call init_stochdata(dt,ls_node) ! init for stochastic component
 endif
 if (nvc > 0) then
    ! stochastic component
    vcfact = 0
    do n=1,nvc
    call patterngenerator_advance(spec_vc_e(:,:,n),spec_vc_o(:,:,n),rpattern_vc(n))
    call scalarspect_to_grid(&
              spec_vc_e(:,:,n),spec_vc_o(:,:,n),wrk2d,&
              ls_node,ls_nodes,max_ls_nodes,&
              lats_nodes_a,global_lats_a,lonsperlar,&
              plnev_a,plnod_a,1)
    vcfact = vcfact + wrk2d
    enddo
    if (vc_logit) then
       vcfact = (2./(1.+exp(vcfact)))-1.
       vcfact = vc*(vcfact+1.)
    else
       vcfact = vc + vcfact  ! deterministic + stochastic
    endif
 else
    vcfact = vc ! purely deterministic
 endif
 do k=1,levs
    vc3d_u(:,1:lats_node_a,k) = vfact_vc(k)*vcfact*dt*vrtgrady(:,:,k)*abs(vrtg(:,:,k))/vrtgradmod(:,:,k)
    vc3d_v(:,1:lats_node_a,k) = -vfact_vc(k)*vcfact*dt*vrtgradx(:,:,k)*abs(vrtg(:,:,k))/vrtgradmod(:,:,k)
 enddo

end subroutine get_pattern_vc

! the routines below are spectral transform routines used internally.

subroutine gatherspec(spharmspec_out,spharmspec_e,spharmspec_o,ls_node)
   implicit none
   integer j,n,l,nm,ierr,ls_node(ls_dim,3),nn
   real(kind_evod), intent(in) :: spharmspec_e(len_trie_ls,2)
   real(kind_evod), intent(in) :: spharmspec_o(len_trio_ls,2)
   complex(r_kind), intent(out) :: spharmspec_out((jcap+1)*(jcap+2)/2)
   complex(r_kind) :: spharmspec((jcap+1)*(jcap+2)/2)
   integer :: idxspec(0:jcap,0:jcap)
   integer indlsod,indlsev,jbasev,jbasod
   include "function_indlsod"
   include "function_indlsev"
   nm = 0
   idxspec = 0
   do l=0,jcap
      do n=l,jcap
         nm = nm + 1
         idxspec(l,n) = nm
      enddo
   enddo
   spharmspec = 0.; spharmspec_out = 0.
   do j = 1, ls_max_node   
      l=ls_node(j,1) ! zonal wavenumber
      jbasev=ls_node(j,2)
      jbasod=ls_node(j,3)
      nn = indlsev(l,l)
      do n=l,jcap,2
         nm = idxspec(l,n)
         spharmspec(nm) = cmplx(spharmspec_e(nn,1),spharmspec_e(nn,2))
         nn = nn + 1
      enddo
      nn=indlsod(l+1,l)
      do n=l+1,jcap,2
         nm = idxspec(l,n)
         spharmspec(nm) = cmplx(spharmspec_o(nn,1),spharmspec_o(nn,2))
         nn = nn + 1
      enddo
   enddo
   ! mpi_reduce
   call mpi_reduce(spharmspec,spharmspec_out,(jcap+1)*(jcap+2)/2,&
                   mpi_complex,mpi_sum,me_l_0,mc_comp,ierr)
end subroutine gatherspec

subroutine computevarspec_eo(spharmspec_e,spharmspec_o,varout,ls_node)
   implicit none
   integer j,n,l,nm,ierr,ls_node(ls_dim,3),nn
   real(kind_evod), intent(in) :: spharmspec_e(len_trie_ls,2)
   real(kind_evod), intent(in) :: spharmspec_o(len_trio_ls,2)
   real(kind_evod), intent(out) :: varout
   real(kind_evod) :: varspec
   integer :: idxspec(0:jcap,0:jcap)
   integer indlsod,indlsev,jbasev,jbasod
   include "function_indlsod"
   include "function_indlsev"
   varspec = 0.; varout = 0.
   do j = 1, ls_max_node
      l=ls_node(j,1) ! zonal wavenumber
      jbasev=ls_node(j,2)
      jbasod=ls_node(j,3)
      nn = indlsev(l,l)
      do n=l,jcap,2
         nm = idxspec(l,n)
         if ( l .EQ. 0) varspec = varspec+0.5*(spharmspec_e(nn,1)**2+spharmspec_e(nn,2)**2)
         if ( l .NE. 0) varspec = varspec+spharmspec_e(nn,1)**2+spharmspec_e(nn,2)**2
         nn = nn + 1
      enddo
      nn=indlsod(l+1,l)
      do n=l+1,jcap,2
         nm = idxspec(l,n)
         if ( l .EQ. 0) varspec = varspec+0.5*(spharmspec_o(nn,1)**2+spharmspec_o(nn,2)**2)
         if ( l .NE. 0) varspec = varspec+spharmspec_o(nn,1)**2+spharmspec_o(nn,2)**2
         nn = nn + 1
      enddo
   enddo
   ! mpi_reduce
   call mpi_allreduce(varspec,varout,1,mpi_real8,mpi_sum,mc_comp,ierr)
end subroutine computevarspec_eo

subroutine subsetspec(spharmspec_in,spharmspec_e,spharmspec_o,ls_node)
   implicit none
   integer j,n,l,nm,ierr,ls_node(ls_dim,3),nn
   real(kind_evod), intent(out) :: spharmspec_e(len_trie_ls,2)
   real(kind_evod), intent(out) :: spharmspec_o(len_trio_ls,2)
   complex(r_kind), intent(in) :: spharmspec_in((jcap+1)*(jcap+2)/2)
   integer :: idxspec(0:jcap,0:jcap)
   integer indlsod,indlsev,jbasev,jbasod
   include "function_indlsod"
   include "function_indlsev"
   nm = 0
   idxspec = 0
   do l=0,jcap
      do n=l,jcap
         nm = nm + 1
         idxspec(l,n) = nm
      enddo
   enddo
   spharmspec_e = 0.; spharmspec_o=0.
   do j = 1, ls_max_node   
      l=ls_node(j,1) ! zonal wavenumber
      jbasev=ls_node(j,2)
      jbasod=ls_node(j,3)
      nn = indlsev(l,l)
      do n=l,jcap,2
         nm = idxspec(l,n)
         spharmspec_e(nn,1) = real(spharmspec_in(nm))
         spharmspec_e(nn,2) = imag(spharmspec_in(nm))
         nn = nn + 1
      enddo
      nn=indlsod(l+1,l)
      do n=l+1,jcap,2
         nm = idxspec(l,n)
         spharmspec_o(nn,1) = real(spharmspec_in(nm))
         spharmspec_o(nn,2) = imag(spharmspec_in(nm))
         nn = nn + 1
      enddo
   enddo
end subroutine subsetspec

subroutine scalarspect_to_grid(&
           trie_ls,trio_ls,datag,&
           ls_node,ls_nodes,max_ls_nodes,&
           lats_nodes_a,global_lats_a,lonsperlar,&
           plnev_a,plnod_a,nlevs)


      implicit none
      real(kind=kind_evod), intent(in) :: trie_ls(len_trie_ls,2,nlevs)
      real(kind=kind_evod), intent(in) :: trio_ls(len_trio_ls,2,nlevs)
      real(kind=kind_phys),  intent(out) :: datag(lonf,lats_node_a,nlevs)
      integer, intent(in) :: ls_node(ls_dim,3),ls_nodes(ls_dim,nodes),&
        nlevs,max_ls_nodes(nodes),lats_nodes_a(nodes),global_lats_a(latg),lonsperlar(latg)
      real(kind=kind_evod),intent(in) :: plnev_a(len_trie_ls,latg2),plnod_a(len_trio_ls,latg2)
! local vars
      real(kind=kind_evod) for_gr_a_1(lon_dim_a,nlevs,lats_dim_a)
      real(kind=kind_evod) for_gr_a_2(lonf,nlevs,lats_dim_a)
      integer              i,j,k
      integer              l,lan,lat
      integer              lons_lat

      call sumfln_slg_gg(trie_ls,&
                  trio_ls,&
                  lat1s_a,&
                  plnev_a,plnod_a,&
                  nlevs,ls_node,latg2,&
                  lats_dim_a,nlevs,for_gr_a_1,&
                  ls_nodes,max_ls_nodes,&
                  lats_nodes_a,global_lats_a,&
                  lats_node_a,ipt_lats_node_a,lon_dim_a,&
                  lonsperlar,lon_dim_a,latg,0)

      do lan=1,lats_node_a
         lat = global_lats_a(ipt_lats_node_a-1+lan)
         lons_lat = lonsperlar(lat)
         CALL FOUR_TO_GRID(for_gr_a_1(1,1,lan),for_gr_a_2(1,1,lan),&
                           lon_dim_a,lonf,lons_lat,nlevs)
      enddo  

      datag = 0.
      do lan=1,lats_node_a
        lat      = global_lats_a(ipt_lats_node_a-1+lan)
        lons_lat = lonsperlar(lat)
        do k=1,nlevs
          do i=1,lons_lat
            datag(i,lan,k) = for_gr_a_2(i,k,lan)
          enddo
        enddo
      enddo

      return
      end subroutine scalarspect_to_grid

 subroutine scalargrid_to_spect(&
           trie_ls,trio_ls,datag,&
           ls_node,ls_nodes,max_ls_nodes,&
           lats_nodes_a,global_lats_a,lonsperlar,&
           plnew_a,plnow_a,nlevs)

      implicit none
      real(kind=kind_evod), intent(out) :: trie_ls(len_trie_ls,2,nlevs)
      real(kind=kind_evod), intent(out) :: trio_ls(len_trio_ls,2,nlevs)
      real(kind=kind_phys),  intent(in) :: datag(lonf,lats_node_a,nlevs)
      integer, intent(in) :: ls_node(ls_dim,3),ls_nodes(ls_dim,nodes),&
        nlevs,max_ls_nodes(nodes),lats_nodes_a(nodes),global_lats_a(latg),lonsperlar(latg)
      real(kind=kind_evod),intent(in) :: plnew_a(len_trie_ls,latg2),plnow_a(len_trio_ls,latg2)
! local vars
      real(kind=kind_evod) for_gr_a_1(lon_dim_a,nlevs,lats_dim_a)
      real(kind=kind_evod) for_gr_a_2(lonf,nlevs,lats_dim_a)
      integer              i,j,k
      integer              l,lan,lat
      integer              lons_lat

      trie_ls = 0.; trio_ls = 0.

      do lan=1,lats_node_a
        lat = global_lats_a(ipt_lats_node_a-1+lan)
        lons_lat = lonsperlar(lat)
        do k=1,nlevs
          do i=1,lons_lat
            for_gr_a_2(i,k,lan) = datag(i,lan,k)
          enddo
        enddo
      enddo

      do lan=1,lats_node_a
         lat = global_lats_a(ipt_lats_node_a-1+lan)
         lons_lat = lonsperlar(lat)
         call grid_to_four(for_gr_a_2(1,1,lan),for_gr_a_1(1,1,lan),&
                           lonf,lon_dim_a,lons_lat,nlevs)
      enddo

      call four2fln_gg(lats_dim_a,nlevs,nlevs,for_gr_a_1,&
                    ls_nodes,max_ls_nodes,&
                    lats_nodes_a,global_lats_a,lon_dim_a,&
                    lats_node_a,ipt_lats_node_a,&
                    lat1s_a,lon_dim_a,latg,latg2,&
                    trie_ls(1,1,1), trio_ls(1,1,1),&
                    plnew_a, plnow_a,&
                    ls_node,0,&
                    nlevs,nlevs)

      return
      end subroutine scalargrid_to_spect

      subroutine vrtdivspect_to_uvgrid(&
           trie_di,trio_di,trie_ze,trio_ze,&
           uug,vvg,&
           ls_node,ls_nodes,max_ls_nodes,&
           lats_nodes_a,global_lats_a,lonsperlar,&
           epsedn,epsodn,snnp1ev,snnp1od,plnev_a,plnod_a,nlevs)

      implicit none
      real(kind=kind_evod), intent(in) :: trie_di(len_trie_ls,2,nlevs)
      real(kind=kind_evod), intent(in) :: trio_di(len_trio_ls,2,nlevs)
      real(kind=kind_evod), intent(in) :: trie_ze(len_trie_ls,2,nlevs)
      real(kind=kind_evod), intent(in) :: trio_ze(len_trio_ls,2,nlevs)
      real(kind=kind_phys),  intent(out) :: uug(lonf,lats_node_a,nlevs)
      real(kind=kind_phys),  intent(out) :: vvg(lonf,lats_node_a,nlevs)
      integer, intent(in) :: ls_node(ls_dim,3),ls_nodes(ls_dim,nodes),&
        nlevs,max_ls_nodes(nodes),lats_nodes_a(nodes),global_lats_a(latg),lonsperlar(latg)
      real(kind=kind_evod),intent(in) ::  epsedn(len_trie_ls),&
       epsodn(len_trio_ls),snnp1ev(len_trie_ls),snnp1od(len_trio_ls),&
       plnev_a(len_trie_ls,latg2),plnod_a(len_trio_ls,latg2)
! local vars
      real(kind=kind_evod) trie_ls(len_trie_ls,2,2*nlevs)
      real(kind=kind_evod) trio_ls(len_trio_ls,2,2*nlevs)
      real(kind=kind_evod) for_gr_a_1(lon_dim_a,2*nlevs,lats_dim_a)
      real(kind=kind_evod) for_gr_a_2(lonf,2*nlevs,lats_dim_a)
      integer              i,j,k
      integer              l,lan,lat
      integer              lons_lat
      real (kind=kind_evod) tx1

      do k=1,nlevs
        call dezouv(trie_di(1,1,k),       trio_ze(1,1,k),&
                    trie_ls(1,1,k), trio_ls(1,1,nlevs+k),&
                    epsedn,epsodn,snnp1ev,snnp1od,ls_node)
        call dozeuv(trio_di(1,1,k),       trie_ze(1,1,k),&
                    trio_ls(1,1,k), trie_ls(1,1,nlevs+k),&
                    epsedn,epsodn,snnp1ev,snnp1od,ls_node)
      enddo

      call sumfln_slg_gg(trie_ls,&
                  trio_ls,&
                  lat1s_a,&
                  plnev_a,plnod_a,&
                  2*nlevs,ls_node,latg2,&
                  lats_dim_a,2*nlevs,for_gr_a_1,&
                  ls_nodes,max_ls_nodes,&
                  lats_nodes_a,global_lats_a,&
                  lats_node_a,ipt_lats_node_a,lon_dim_a,&
                  lonsperlar,lon_dim_a,latg,0)

      do lan=1,lats_node_a
         lat = global_lats_a(ipt_lats_node_a-1+lan)
         lons_lat = lonsperlar(lat)
         CALL FOUR_TO_GRID(for_gr_a_1(1,1,lan),for_gr_a_2(1,1,lan),&
                           lon_dim_a,lonf,lons_lat,2*nlevs)
      enddo  

      uug = 0.; vvg = 0.
      do lan=1,lats_node_a
        lat      = global_lats_a(ipt_lats_node_a-1+lan)
        lons_lat = lonsperlar(lat)
        tx1      = 1. / coslat_a(lat)
        do k=1,nlevs
          do i=1,lons_lat
            uug(i,lan,k) = for_gr_a_2(i,k,lan) * tx1
            vvg(i,lan,k) = for_gr_a_2(i,nlevs+k,lan) * tx1
          enddo
        enddo
      enddo

      return
      end subroutine vrtdivspect_to_uvgrid
subroutine dump_patterns(sfile)
    implicit none
    character*9 :: sfile
    integer :: stochlun,k,n
    stochlun=99
    if (me .EQ. me_l_0) then
       if (nsppt > 0 .OR. nvc > 0 .OR. nshum > 0 .OR. nskeb > 0) then
          OPEN(stochlun,file=sfile//trim(ens_nam),form='unformatted')
          print*,'open ',sfile,' for output'
       endif
    endif
    if (nsppt > 0) then
       do n=1,nsppt 
       call write_pattern(rpattern_sppt(n),spec_sppt_e(:,:,n),spec_sppt_o(:,:,n),stochlun)
       enddo
    endif
    if (nvc > 0) then
       do n=1,nvc 
       call write_pattern(rpattern_vc(n),spec_vc_e(:,:,n),spec_vc_o(:,:,n),stochlun)
       enddo
    endif
    if (nshum > 0) then
       do n=1,nshum
       call write_pattern(rpattern_shum(n),spec_shum_e(:,:,n),spec_shum_o(:,:,n),stochlun)
       enddo
    endif
    if (nskeb > 0) then
       do n=1,nskeb
       do k=1,levs
          call write_pattern(rpattern_skeb(n),spec_skeb_e(:,:,k,n),spec_skeb_o(:,:,k,n),stochlun)
       enddo
       enddo
    endif
    close(stochlun)
 end subroutine dump_patterns
 subroutine restore_patterns(sfile)
    implicit none
    character*9 :: sfile
    integer :: stochlun,k,n
    stochlun=99
    if (me .EQ. me_l_0) then
       if (nsppt > 0 .OR. nvc > 0 .OR. nshum > 0 .OR. nskeb > 0) then
          OPEN(stochlun,file=sfile//trim(ens_nam),form='unformatted',status='old')
       endif
    endif
    if (nsppt > 0) then
       do n=1,nsppt
       call read_pattern(rpattern_sppt(n),spec_sppt_e(:,:,n),spec_sppt_o(:,:,n),stochlun)
       enddo
    endif
    if (nvc > 0) then
       do n=1,nvc
       call read_pattern(rpattern_vc(n),spec_vc_e(:,:,n),spec_vc_o(:,:,n),stochlun)
       enddo
    endif
    if (nshum > 0) then
       do n=1,nshum
       call read_pattern(rpattern_shum(n),spec_shum_e(:,:,n),spec_shum_o(:,:,n),stochlun)
       enddo
    endif
    if (nskeb > 0) then
       do n=1,nskeb
       do k=1,levs
          call read_pattern(rpattern_skeb(n),spec_skeb_e(:,:,k,n),spec_skeb_o(:,:,k,n),stochlun)
       enddo
       enddo
    endif
    close(stochlun)
 end subroutine restore_patterns
subroutine read_pattern(rpattern,pattern2d_e,pattern2d_o,lunptn)
   implicit none
   real(kind_evod), intent(inout) :: pattern2d_e(len_trie_ls,2)
   real(kind_evod), intent(inout) :: pattern2d_o(len_trio_ls,2)
   type(random_pattern), intent(inout) :: rpattern
   integer, intent(in) :: lunptn
   real(kind_evod),allocatable  :: pattern2d(:)
   integer nm,nn,ierr

   allocate(pattern2d(2*ndimspec))

   ! read only on root process, and send to all tasks
   if (me .EQ. me_l_0) then
      read(lunptn) pattern2d
      print*,'reading in random pattern (min/max/size)',minval(pattern2d),maxval(pattern2d),size(pattern2d)
   endif
   call mpi_bcast(pattern2d,2*ndimspec,mpi_r_io_r,me_l_0,mc_comp,ierr)
   ! subset
   do nn=1,len_trie_ls
      nm = rpattern%idx_e(nn)
      if (nm == 0) cycle
      pattern2d_e(nn,1) = pattern2d(nm)
      pattern2d_e(nn,2) = pattern2d(ndimspec+nm)
   enddo
   do nn=1,len_trio_ls
      nm = rpattern%idx_o(nn)
      if (nm == 0) cycle
      pattern2d_o(nn,1) = pattern2d(nm)
      pattern2d_o(nn,2) = pattern2d(ndimspec+nm)
   enddo
   !print*,'after scatter...',me,maxval(pattern2d_e),maxval(pattern2d_o) &
   ! ,minval(pattern2d_e),minval(pattern2d_o)
   deallocate(pattern2d)
 end subroutine read_pattern

 subroutine write_pattern(rpattern,pattern2d_e,pattern2d_o,lunptn)
   implicit none
   real(kind_evod), intent(in) :: pattern2d_e(len_trie_ls,2)
   real(kind_evod), intent(in) :: pattern2d_o(len_trio_ls,2)
   type(random_pattern), intent(inout) :: rpattern
   integer, intent(in) :: lunptn
   real(kind_evod), allocatable  :: pattern2d(:),pattern2d_out(:)
   integer nm,nn,ierr

   allocate(pattern2d(2*ndimspec),pattern2d_out(2*ndimspec))
   pattern2d=0.0
   pattern2d_out=0.0
   ! fill in apprpriate pieces of array
   !print*,'before collection...',me,maxval(pattern2d_e),maxval(pattern2d_o) &
   ! ,minval(pattern2d_e),minval(pattern2d_o)
   do nn=1,len_trie_ls
      nm = rpattern%idx_e(nn)
      if (nm == 0) cycle
      pattern2d(nm)          = pattern2d_e(nn,1)
      pattern2d(ndimspec+nm) = pattern2d_e(nn,2)
   enddo
   do nn=1,len_trio_ls
      nm = rpattern%idx_o(nn)
      if (nm == 0) cycle
      pattern2d(nm)          = pattern2d_o(nn,1)
      pattern2d(ndimspec+nm) = pattern2d_o(nn,2)
   enddo
   call mpi_reduce(pattern2d,pattern2d_out,2*ndimspec,&
                   mpi_r_io_r,mpi_sum,me_l_0,mc_comp,ierr)
  !  write only on root process
   if (me .EQ. me_l_0) then
      print*,'writing out random pattern (min/max/size)',&
      minval(pattern2d_out),maxval(pattern2d_out),size(pattern2d_out)
      !print*,'max/min pattern=',maxval(pattern2d_out),minval(pattern2d_out)
      write(lunptn) pattern2d_out
   endif
   deallocate(pattern2d,pattern2d_out)
 end subroutine write_pattern
end module gfs_get_pattern_mod
