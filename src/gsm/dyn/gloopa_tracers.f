       subroutine gloopa_hyb_slg (deltim,trie_ls,trio_ls,gzie_ln,gzio_ln
     &,                           ls_node,ls_nodes,max_ls_nodes
     &,                           lats_nodes_a,global_lats_a
     &,                           lonsperlat,ldfi,grid_gr_dfi
     &,                           epse,epso,epsedn,epsodn
     &,                           snnp1ev,snnp1od
     &,                           plnev_a,plnod_a,pddev_a
     &,                           pddod_a,plnew_a,plnow_a
     &,                           kdt,batah,lsout,end_step)
!    &                            global_times_a,kdt,batah,lsout,end_step)

      use gfs_dyn_machine     , only : kind_evod
      use gfs_dyn_resol_def   , only : jcap,jcap1,latg,latg2,ntrac,
     &                                 levh,levp1,levs,lnt2,lonf,lota,
     &                                 p_di,p_dlam,p_dphi,p_gz,p_q,
     &                                 p_rq,p_rt,p_te,p_uln,p_vln,
     &                                 p_w,p_x,p_y,p_ze,lotls,adiabatic
      use gfs_dyn_layout1     , only : ipt_lats_node_a,lat1s_a,
     &                                 lats_dim_a,lats_node_a,
     &                                 len_trie_ls,len_trio_ls,
     &                                 lon_dim_a,ls_dim,
     &                                 ls_max_node,me,nodes
      use gfs_dyn_gg_def      , only : rcs2_a,lats_nodes_h,global_lats_h
      use gfs_dyn_date_def    , only : spdmax, fhour
!my12272012
      use namelist_dynamics_def,only : gg_tracers,hybrid,ref_temp,ngptc,
     &                                 settls_dep3ds,settls_dep3dg,
     &                                 ref_pres,nemsio_in,nemsio_out,
     &                                 ldfi_spect
      use gfs_dyn_mpi_def     , only : kind_mpi,
     &                                 mpi_r_mpi,mpi_max,mc_comp
!    &                                 mpi_r_mpi,mpi_real8,mpi_max
      use gfs_dyn_physcons    , only : rerth => con_rerth, rd => con_rd,
     &                                 grav => con_g, omega => con_omega
      use layout_lag          , only : lat1s_h,lats_dim_h,lon_dim_h
      use layout_grid_tracers , only : rgt_a,rgt_h,xhalo,yhalo
      use gfs_dyn_coordinate_def, only : ak5,bk5
      use      pmgrid         , only : platd,plond
      use      slgshr         , only : dlam,lam,nlonex,ra,rdlam,rdlam6
      use gfs_dyn_dfi_mod
!     use gfs_dyn_io_header,    only: gz_grid
!
      implicit none

      type(gfs_dfi_grid_gr),intent(inout) :: grid_gr_dfi
      logical,intent(in) :: ldfi

      real(kind=kind_evod),intent(inout) ::
     &                           trie_ls(len_trie_ls,2,lotls),
     &                           trio_ls(len_trio_ls,2,lotls)
      real(kind=kind_evod), dimension(len_trie_ls) :: epse, epsedn,
     &                                                snnp1ev
      real(kind=kind_evod), dimension(len_trio_ls) :: epso, epsodn,
     &                                                snnp1od
      real(kind=kind_evod), dimension(len_trie_ls,latg2) ::
     &                        plnev_a, pddev_a, plnew_a
      real(kind=kind_evod), dimension(len_trio_ls,latg2) ::
     &                        plnod_a, pddod_a, plnow_a
      real(kind=kind_evod) :: gzie_ln(len_trie_ls,2),
     &                        gzio_ln(len_trio_ls,2),
!    &                        global_times_a(latg,nodes),
     &                        batah,deltim
      integer, dimension(nodes) :: max_ls_nodes, lats_nodes_a
      integer, dimension(latg)  :: global_lats_a, lonsperlat
      integer                   :: ls_node(ls_dim,3),
     &                             ls_nodes(ls_dim,nodes)
!    &,                            ndexev(len_trie_ls)
!    &,                            ndexod(len_trio_ls)
!
!Moo  real(kind=kind_evod), parameter ::  gama=6.5/1000.,
!Moo &                                    expon=rd*gama/grav,
!Moo &                                    psref=101.325,  tzero=288.
      real(kind=kind_evod), parameter ::  gama=0.0/1000.,
     &                                    expon=rd*gama/grav,
     &                                    zero=0.0d0

!    &                                    psref=101.325,  zero=0.0d0
!    &                                    psref=101.325,  tzero=300.

      logical lsout
      integer iter_slg,ifirst,i_1,i_2,kdt,
     &        kdt_save,indlsev,jbasev,ini_slg,ini_dp,
     &        njeff,lon,iblk,ngptcd,ngptcs,xlon,ylan,

     &        i,ierr,j,k,kap,kar,kat,kau,kav,kdrlam,
     &        kdrphi,kdtlam,kdtphi,kdulam,kduphi,kdvlam,
     &        kdvphi,ksd,ksplam,kspphi,ksq,ksrc,ksrg,kst,
     &        ksu,ksv,ksz,l,lan,lat,lmax,locl,lmax2p1,
     &        lons_lat,n,node,nvcn,jtem, ktem, kk, nn,
     &        s_lat_in_pe,n_lat_in_pe,
     &        ipt_ls,n_ggp,s_ggp,iter_max,ntr,kl1,kt1,i2

!     integer n_rgt_a

!     integer ,allocatable, dimension(:) :: lats_nodes_h,global_lats_h
!
      real(kind=kind_evod) :: r_rt          ! used in ln(ps) 'linarization'
     &,                       tref, tzero,rcs2_loc,psref
     &,                       consv, xvcn, reall
     &,                       lat_phi,degrad,hfdt, dummy

!
      real(kind=kind_evod),dimension(levs) :: ak,bk,pkref,tkref,
     &                                        spdmax_node
      real(kind=kind_mpi), dimension(levs) :: spdmax_nodesm
     &,                                       spdmax_nodem

      real(kind=kind_evod) ,allocatable, dimension(:) :: tb_k,etaint,
     &           etamid, omega_v, cosf1, rcos1

      real(kind=kind_evod) ,allocatable, dimension(:,:) ::
     &           dpdt_a,       sd_m,      bak_grid_a,
     &           phige,        phigo,
!    &                         grad_gzlam,
     &           gz_grid,      grad_gzlam,
     &           grad_gzphi,   rlcs2,      spdlat

      real(kind=kind_evod) ,allocatable, dimension(:,:,:) ::
     &            lammp_a,    phimp_a, sigmp_a,
     &            for_coef_a, dyn_coef_a,
     &            bak_coef_a, for_grid_a, dyn_grid_a,
     &            unlm,       unla,
     &            vnlm,       vnla,
     &            tnlm,       tnla,     z_m,
     &            ud3_a,      vd3_a,    td3_a,
     &            dpdt_d3_a,  ug_h,
     &            ug_m,       vg_h,
     &            vg_m,       ww_h,
     &            ww_m,       ud3_h1,   ud3_h2,
     &            vd3_h1,     vd3_h2,   td3_h1,     td3_h2,
     &            dpdt_d3_h1, dpdt_d3_h2
!

      save  gz_grid, phige, phigo, r_rt, tb_k, omega_v
     &,     grad_gzlam, grad_gzphi
     &,     lammp_a, phimp_a, sigmp_a, z_m,  sd_m, unlm, vnlm, tnlm
     &,     etamid,  etaint,  rcos1,   ug_m, vg_m, ww_m
     &,     s_lat_in_pe, n_lat_in_pe,  iter_max, rlcs2
     &,     ifirst, kdt_save, ini_slg, ini_dp

      logical go_forward,end_step,lprnt
      logical, save :: first
      data first/.true./

      data ifirst /1/,kdt_save/0/
!
      include 'function_indlsev'

!

!mjr  lon_dim_h  = lonf + 1 + 2*xhalo  !  even or odd
!

      lprnt  = .false.
      ngptcs = ngptc
      ngptcd = ngptc
      hfdt   = 0.5*deltim

      if (me == 0 .and. kdt <10 ) write(0,*)'ngptcd = ',ngptcd,
     &                                     ' ngptcs = ',ngptcs

      ksz     =                   1
      ksd     =                   1
      kst     =   levs          + 1
      ksq     = 2*levs          + 1
      ksplam  = 2*levs          + 2
      kspphi  = 2*levs          + 3

      ksu     = 2*levs          + 4
      ksv     = 3*levs          + 4
      ksrc    = 4*levs          + 4

      ksrg    = 2*levs          + 4

      kdtphi  =                   1
      kdrphi  =   levs          + 1
      kdtlam  =   levs          + 1
      kdrlam  = 2*levs +   levh + 1
      kdulam  = 2*levs + 2*levh + 1
      kdvlam  = 3*levs + 2*levh + 1
      kduphi  = 4*levs + 2*levh + 1
      kdvphi  = 5*levs + 2*levh + 1
      kau     =                   1
      kav     =   levs          + 1
      kat     = 2*levs          + 1
      kap     = 3*levs          + 1
      kar     = 3*levs          + 2

      allocate ( for_coef_a(lon_dim_a,(4*levs+levh+3),lats_dim_a) )
      allocate ( for_grid_a(lonf     ,(2*levs+levh+3),lats_dim_a) )

      if (ifirst == 1) then                      ! first ifirst
        for_grid_a = 0.0
        allocate (tb_k(levs))
        tref   = ref_temp
        tzero  = ref_temp
        psref  = ref_pres
        r_rt   = 1./(rd*tref)
!!!!    r_rt   = 0.0   ! for no ln(ps) linearization
!--------------------------------------------------
        degrad = 45./atan(1.)

        do k=1,levs
          ak(k)    = 0.5*(ak5(k)+ak5(k+1))
          bk(k)    = 0.5*(bk5(k)+bk5(k+1))
          pkref(k) = ak(k) + bk(k)*psref
          tkref(k) = tzero*(pkref(k)/psref)**expon

          tb_k(k)  = bk(k)*psref*expon*tkref(k)/pkref(k) !order is top2bot
          tb_k(k)  = - tb_k(k)/(rd*tzero)                !check tzero or tref
!!!!      tb_k(k)  = 0.0 !  "for no temp. perturbation"
        enddo
        if(me == 0)then
          do k=1,levs
            print 111, k, pkref(k),tkref(k),tb_k(k),ak(k),bk(k)
          enddo
          print 109
          print*,' deltim in gloopa_slg=',deltim,' batah=',batah
          print 109
        endif

109     format(1h ,'------------------------------------')
111     format('in gloopa k=',i3,2x,' pkref=',e10.3,2x,'tkref=',e10.3,
     & ' tb_k=',e12.4,2x,'ak=',e10.3,2x,'bk=',e10.3)

        call set_pmgrid_h(lonf,latg,levs)

        allocate(dlam(platd),rdlam(platd),rdlam6(platd))
        allocate(lam(plond,platd),nlonex(platd))

        ra = 1./rerth

        if (.not. allocated(lats_nodes_h))
     &                    allocate ( lats_nodes_h(nodes) )
        if (.not. allocated(lat1s_h))
     &                    allocate (lat1s_h(0:jcap))
        if (.not. allocated(global_lats_h))
     &                    allocate (global_lats_h(latg+2*yhalo*nodes))

        call getcon_lag(lats_nodes_a, global_lats_a,
     &                  lats_nodes_h, global_lats_h,
     &                  lonsperlat,   xhalo, yhalo)
!
      if (.not.allocated(rgt_h))
     &         allocate (rgt_h(lon_dim_h,levs,lats_dim_h,ntrac))

!     write(0,*)' rgt_h dims=',size(rgt_h,dim=1),size(rgt_h,dim=2),
!    &size(rgt_h,dim=3),size(rgt_h,dim=4)
!

        s_lat_in_pe = latg + 1
     &              - global_lats_a(ipt_lats_node_a-1+lats_node_a)
        s_lat_in_pe = s_lat_in_pe - yhalo + 2

!       write(0,*)'me1z=',me,' s_lat_in_pe=',s_lat_in_pe
!    &,' lats_dim_h=',lats_dim_h

        n_lat_in_pe = s_lat_in_pe + lats_dim_h - 1

!       write(0,*)'me1z=',me,' n_lat_in_pe=',n_lat_in_pe
!    &,' lats_dim_h=',lats_dim_h

!---------------------------------------------------
        allocate (  etamid(levs), etaint(levp1) )
        do k=1,levp1
          etaint(k) = ak5(k)/psref + bk5(k)
          if (me == 0) print*,'etaint_def  etaint(k) k=',etaint(k),k
        enddo
        do k=1,levs
          etamid(k) = (etaint(k) + etaint(k+1) ) * 0.5
          if (me == 0) print*,'etamid_def  etamid(k) k=',etamid(k),k
        enddo
        if (me == 0) then
          do k=2,levs             ! note k=1 not used
            print*,'detamd_def k=',k, etamid(k)-etamid(k-1)
          enddo
        endif
!---------------------------------------------------
        allocate ( omega_v(lats_dim_a) )
        allocate ( rlcs2(jcap1,lats_node_a) )
        allocate ( cosf1(lats_dim_a),rcos1(lats_dim_a) )
        iter_max = 0

        do lan=1,lats_node_a                        ! begin  lan loop 0
          lat      = global_lats_a(ipt_lats_node_a-1+lan)
          lons_lat = lonsperlat(lat)
          iter_max = max ( iter_max , (lons_lat+ngptcd-1)/ngptcd )
          ipt_ls   = min(lat,latg-lat+1)
!--------------------------------------------
          lmax     = min(jcap,lons_lat/2)
          do i=1,lmax+1
            if ( ipt_ls >= lat1s_a(i-1) ) then
               reall        = i - 1
               rlcs2(i,lan) = reall*rcs2_a(ipt_ls)*ra
            else
               rlcs2(i,lan) = zero
            endif
          enddo
!--------------------------------------------
          cosf1(lan)   = 1. / sqrt( rcs2_a(ipt_ls) )
          lat_phi      = degrad*acos(cosf1(lan))

          rcos1(lan)   = sqrt( rcs2_a(min(lat,latg-lat+1)) )
          omega_v(lan) = 2.*omega*rerth*cosf1(lan)
!         consv        = 2.*omega*rerth
          lat_phi      = degrad*acos(1./rcos1(lan))
!         if (me == 0) print*,' consv=2. omega rerth ',consv
!sela     omega_v(lan) = 0.

        enddo                                       ! end of lan loop 0
!***************************************************
        allocate (    phige(len_trie_ls,2) )
        allocate (    phigo(len_trio_ls,2) )
!       allocate (    gz_grid(lonf,lats_dim_a) )
        allocate (    gz_grid(lon_dim_a,lats_dim_a) )
        allocate ( grad_gzlam(lon_dim_a,lats_dim_a) )
        allocate ( grad_gzphi(lon_dim_a,lats_dim_a) )

!$omp parallel do private(i,j)
        do j=1,lats_dim_a
          do i=1,lon_dim_a
            gz_grid(i,j) = 0.0
          enddo
        enddo

        call get_topo_grid_grad (cosf1,
     &                           grad_gzlam,grad_gzphi,gzie_ln,gzio_ln,
     &                           gz_grid, phige,phigo,
     &                           for_coef_a,for_grid_a,
     &                           trie_ls,trio_ls,
     &                           ls_node,ls_nodes,max_ls_nodes,
     &                           lats_nodes_a,global_lats_a,
     &                           lonsperlat,epse,epso,
     &                           plnev_a,plnod_a,
     &                           0,r_rt) ! 1 for write ; 0 for no write
!***************************************************
        deallocate (cosf1)
        allocate (  ug_m(lon_dim_h,levs  ,lats_dim_h) )
        allocate (  vg_m(lon_dim_h,levs  ,lats_dim_h) )
        allocate (  ww_m(lon_dim_h,levs  ,lats_dim_h) )
        allocate ( lammp_a(lonf,levs,lats_dim_a) )
        allocate ( phimp_a(lonf,levs,lats_dim_a) )
        allocate ( sigmp_a(lonf,levs,lats_dim_a) )
        allocate (    unlm(lonf,levs,lats_dim_a) )
        allocate (    vnlm(lonf,levs,lats_dim_a) )
        allocate (    tnlm(lonf,levs,lats_dim_a) )
        allocate (     z_m(lonf,levs,lats_dim_a) )
        allocate (    sd_m(lon_dim_a,lats_dim_a) )
        if (.not. allocated(rgt_a))
     &          allocate (rgt_a(lonf,levs,lats_dim_a,ntrac))
!!
        ifirst = 0
      endif !  first ifirst

!!----------------------------------------------------------------
!for halo testing  ud3_h =0. ; vd3_h =0. ; td3_h =0. ; dpdt_d3_h =0.
!!----------------------------------------------------------------
   

!        if(kdt < 3)write(0,*)' me =',me,'  iter_max=',iter_max

      call delnpe(trie_ls(1,1,p_q   ), trio_ls(1,1,p_dphi),
     &            trie_ls(1,1,p_dlam), epse, epso,ls_node)
      call delnpo(trio_ls(1,1,p_q   ), trie_ls(1,1,p_dphi),
     &            trio_ls(1,1,p_dlam), epse, epso,ls_node)

!$omp parallel do shared(trie_ls,trio_ls)
!$omp+shared(epsedn,epsodn,snnp1ev,snnp1od,ls_node)
!$omp+private(k)
      do k=1,levs
         call dezouv(trie_ls(1,1,p_di +k-1), trio_ls(1,1,p_ze +k-1),
     &               trie_ls(1,1,p_uln+k-1), trio_ls(1,1,p_vln+k-1),
     &               epsedn,epsodn,snnp1ev,snnp1od,ls_node)
         call dozeuv(trio_ls(1,1,p_di +k-1), trie_ls(1,1,p_ze +k-1),
     &               trio_ls(1,1,p_uln+k-1), trie_ls(1,1,p_vln+k-1),
     &               epsedn,epsodn,snnp1ev,snnp1od,ls_node)
      enddo

!     call synctime()

      call sumfln_slg_gg(trie_ls(1,1,p_di),
     &                   trio_ls(1,1,p_di),
     &                   lat1s_a,
     &                   plnev_a,plnod_a,
     &                   4*levs+3,ls_node,latg2,
     &                   lats_dim_a,4*levs+levh+3,
     &                   for_coef_a,
     &                   ls_nodes,max_ls_nodes,
     &                   lats_nodes_a,global_lats_a,
     &                   lats_node_a,ipt_lats_node_a,lon_dim_a,
     &                   lonsperlat,lon_dim_a,latg,0)

      if (.not. gg_tracers .or. (kdt == 1 .and. .not. nemsio_in)) then
         call sumfln_slg_gg(trie_ls(1,1,p_rq),
     &                      trio_ls(1,1,p_rq),
     &                      lat1s_a,
     &                      plnev_a,plnod_a,
     &                      levh,ls_node,latg2,
     &                      lats_dim_a,4*levs+levh+3,
     &                      for_coef_a,
     &                      ls_nodes,max_ls_nodes,
     &                      lats_nodes_a,global_lats_a,
     &                      lats_node_a,ipt_lats_node_a,lon_dim_a,
     &                      lonsperlat,lon_dim_a,latg,4*levs+3)
      endif

      allocate ( dyn_coef_a(lon_dim_a,2*levs,lats_dim_a) )
      call sumder2_slg(trie_ls(1,1,p_te),
     &                 trio_ls(1,1,p_te),
     &                 lat1s_a,
     &                 pddev_a,pddod_a,
     &                 levs,ls_node,latg2,
     &                 lats_dim_a,2*levs,
     &                 dyn_coef_a,
     &                 ls_nodes,max_ls_nodes,
     &                 lats_nodes_a,global_lats_a,
     &                 lats_node_a,ipt_lats_node_a,lon_dim_a,
     &                 lonsperlat,lon_dim_a,latg,0)

      do k=1,levs
        spdmax_node(k) = zero
      enddo

      allocate ( dyn_grid_a(lonf, 2*levs, lats_dim_a) )
      allocate ( ug_h(lon_dim_h,levs, lats_dim_h) )
      allocate ( vg_h(lon_dim_h,levs, lats_dim_h) )
      allocate ( ww_h(lon_dim_h,levs, lats_dim_h) )

      ug_h = 0.0 ; vg_h = 0.0

      do lan=1,lats_node_a   !begin lan loop 1 (four and grid)
!       timer1   = rtc()
        lat      = global_lats_a(ipt_lats_node_a-1+lan)
        lons_lat = lonsperlat(lat)
        lmax     = min(jcap,lons_lat/2)
        lmax2p1  = lmax + lmax + 3
        ylan     = lats_node_a + 1 - lan
!$omp parallel do private(k,i,i2,kl1,kt1)
        do k=1,levs
          kl1 = kdtlam-1+k
          kt1 = kst   -1+k
          do i=1,lmax+1
            i2  = i + i
            dyn_coef_a(i2-1,kl1,lan) = - for_coef_a(i2 ,kt1,lan)
     &                               *   rlcs2(i,lan)
            dyn_coef_a(i2  ,kl1,lan) =   for_coef_a(i2-1,kt1,lan)
     &                               *   rlcs2(i,lan)
          enddo
          if (lmax2p1 <= lons_lat+2) then
            do i = lmax2p1, lons_lat+2
              dyn_coef_a(i,kl1,lan) = zero
!             dyn_coef_a(i,kdtphi-1+k,lan) = zero
            enddo
          endif
        end do
!!!!  synthesize  grids in South to north order
!     -----------------------------------------
        call four_to_grid(for_coef_a(1,ksd,lan),
     &                    for_grid_a(1,ksd,ylan),
     &                    lon_dim_a,lon_dim_a-2,lons_lat,2*levs+3)
        call four_to_grid(for_coef_a(1,ksu,lan),
     &                    ug_h(1+xhalo,1,ylan+yhalo),
     &                    lon_dim_a, lon_dim_h, lons_lat, levs)
        call four_to_grid(for_coef_a(1,ksv,lan),
     &                    vg_h(1+xhalo,1,ylan+yhalo),
     &                    lon_dim_a, lon_dim_h, lons_lat, levs)
        if ( .not. gg_tracers
     &       .or.  (kdt == 1 .and. .not. nemsio_in)) then
           call four_to_grid(for_coef_a(1,ksrc,lan),
     &                       for_grid_a(1,ksrg,ylan),
     &                       lon_dim_a,lon_dim_a-2,lons_lat,levh)
        endif
        call four_to_grid(dyn_coef_a(1,1,lan),
     &                    dyn_grid_a(1,1,ylan),
     &                    lon_dim_a,lon_dim_a-2,lons_lat,2*levs)
!       timer2=rtc()
!       global_times_a(lat,me+1)=timer2-timer1

      enddo !sela finished lan loop 1
!
!     if ((lsout .and. kdt > 0) .or. ldfi) THEN
      if (gg_tracers) then
        if (kdt == 1 .and. .not. nemsio_in ) then
!$omp parallel do private(i,k,kk,n,lan,ylan,lat,lons_lat)
          do lan=1,lats_node_a
            ylan     = yhalo + lan
            lat      = global_lats_a(ipt_lats_node_a+lats_node_a-lan)
            lons_lat = lonsperlat(lat)
            do n=1,ntrac
              do k=1,levs
                kk = ksrg + n*levs - k
                do i=1,lons_lat
                  rgt_h(xhalo+i,k,ylan,n) = max(0.,for_grid_a(i,kk,lan))
                enddo
              enddo
            enddo
          enddo
        elseif (kdt > 1 .and. adiabatic) then
!$omp parallel do private(i,k,kk,n,lan,ylan,lat,lons_lat)
          do lan=1,lats_node_a
            ylan     = yhalo + lan
            lat      = global_lats_a(ipt_lats_node_a+lats_node_a-lan)
            lons_lat = lonsperlat(lat)
            do n=1,ntrac
              do k=1,levs
                kk = levs+1-k
                do i=1,lons_lat
                  rgt_h(xhalo+i,k,ylan,n) = rgt_a(i,kk,lan,n)
                enddo
              enddo
            enddo
          enddo
        endif
            
!$omp parallel do private(i,k,kk,n,lan,ylan,lat,lons_lat)
        do lan=1,lats_node_a
          ylan     = yhalo + lan
          lat      = global_lats_a(ipt_lats_node_a+lats_node_a-lan)
          lons_lat = lonsperlat(lat)
          do n=1,ntrac
            do k=1,levs
              kk = ksrg + n*levs - k
              do i=1,lons_lat
                for_grid_a(i,kk,lan) = rgt_h(xhalo+i,k,ylan,n)
              enddo
            enddo
          enddo
        enddo
      endif              ! gg_tracers if
!     endif
!
!                                       digital filter component
!     write(0,*)' in gloopa me=',me,' kdt=',kdt,' ldfi=',ldfi,
!    &' before gfs_dficoll_dynamics_slg'

      IF (ldfi .and. .not.ldfi_spect) then
        call gfs_dficoll_dynamics_slg(for_grid_a,ug_h,vg_h,
     &                                kst, ksq, ksrg,
     &                                ak5, bk5,
     &                                global_lats_a,lonsperlat,
     &                                gz_grid,grid_gr_dfi)
      endif
!

      if (nemsio_out .and. (lsout .and. kdt > 0) ) then
        CALL wrtout_dynamics_slg_nemsio(fhour, for_grid_a,ug_h,vg_h,
     &                                  global_lats_a,lonsperlat,
     &                                  kst, ksq, ksrg, gz_grid)
      endif

      deallocate ( for_coef_a )
      deallocate ( dyn_coef_a )

      if(end_step) return

      allocate ( spdlat(levs,iter_max ) )

      if (kdt == 1) then ! calculate max meridional velocity

        do lan=1,lats_node_a   !sela begin lan loop 2
!         timer1   = rtc()
          lat      = global_lats_a(ipt_lats_node_a-1+lan)
          lons_lat = lonsperlat(lat)

!sela  ngptcd=(lons_lat+num_threads-1)/num_threads ! opt ngptc

!$omp parallel do schedule(dynamic,1) private(lon,njeff,iblk)
          do lon=1,lons_lat,ngptcd
            njeff = min(ngptcd,lons_lat-lon+1)
            iblk  = (lon-1)/ngptcd+ 1
            call gfidi_speed(lon_dim_a, lon_dim_h, njeff, lat,
     &                       vg_h(xhalo+lon,1,lats_node_a+1-lan+yhalo),
     &                       rcs2_a(min(lat,latg-lat+1)),spdlat(1,iblk))
          enddo   ! omp lon
          iblk = 1
          do lon=1,lons_lat,ngptcd
            do k=1,levs
              spdmax_node(k) = max(spdmax_node(k),spdlat(k,iblk))
            enddo
            iblk = iblk + 1
          enddo

        enddo ! finished  lan loop 2

        spdmax_nodem = spdmax_node

        call mpi_allreduce(spdmax_nodem,spdmax_nodesm,levs,mpi_r_mpi,
     &                     mpi_max,mc_comp,ierr)
        do k=1,levs
          spdmax(k) = sqrt(spdmax_nodesm(k))
        enddo

        if ( me == 0 ) then
          print 101,(spdmax(k),k=1,levs)
101       format(' v_max(01:10)=',10f5.0,:/' v_max(11:20)=',10f5.0,
     x         :/' v_max(21:30)=',10f5.0,:/' v_max(31:40)=',10f5.0,
     x         :/' v_max(41:50)=',10f5.0,:/' v_max(51:60)=',10f5.0,
     x         :/' v_max(61:70)=',10f5.0,:/' v_max(71:80)=',10f5.0,
     x         :/' v_max(81:90)=',10f5.0,:/' v_max(91:00)=',10f5.0)
        endif
      endif ! finished meridional velocity

      do k=1,levs
       spdmax_node(k) = zero
      enddo

      go_forward = .false.
      iter_slg   = 2
!
!
      if (kdt == 1 .or. first) then       
          iter_slg   = 4         !if(kdt == 1 )
          go_forward = .true.

!       if (kdt == 1) then
!         iter_slg   = 4         !if(kdt == 1 )
!         go_forward = .true.
!       endif

        ini_slg      = 1
        ini_dp       = 1

        first = .false.
      endif

      if (kdt > 1 .and. kdt < kdt_save) then
        iter_slg   = 4   !if(kdt > 1 .and. kdt < kdt_save)
        go_forward = .true.
        ini_slg    = 0
        ini_dp     = 1
      endif 

      if ( ini_slg == 1 .or. ini_dp == 1) then
        i_1 = 1          ! this call of slgscan must be made from 
        i_2 = 2          ! an unthreaded region  for initializing.
!   Note 3-D arrays do not need to be declared before initialization
!   (ini_slg .or. ini_dp = 1) -- dummy argument passed
!sk10022012
        call slgscan_h(i_1, i_2, s_lat_in_pe, n_lat_in_pe,
     &                 j, lat, lon_dim_h,
     &                 deltim, iter_slg, etamid, lats_node_a ,
     &                 ug_h, vg_h, ww_h, 
     &                 ug_m,vg_m,ww_m,
     &                 lammp_a,phimp_a,sigmp_a,me,
     &                 dummy, dummy,
     &                 dummy, dummy, dummy,
     &                 dummy, dummy, dummy,
     &                 dummy, dummy, dummy,
     &                 dummy, dummy, dummy,
     &                 global_lats_a,lonsperlat,ini_slg,ini_dp,.false.)
      endif

      ini_slg  = 0
      ini_dp   = 0

      kdt_save = kdt

      allocate ( unla(lonf,levs,lats_dim_a) )
      allocate ( vnla(lonf,levs,lats_dim_a) )
      allocate ( tnla(lonf,levs,lats_dim_a) )
      allocate ( dpdt_a(lon_dim_a,lats_dim_a) )

      allocate (     ud3_h1(lon_dim_h,levs  ,lats_dim_h) )
      allocate (     ud3_h2(lon_dim_h,levs  ,lats_dim_h) )
      allocate (     vd3_h1(lon_dim_h,levs  ,lats_dim_h) )
      allocate (     vd3_h2(lon_dim_h,levs  ,lats_dim_h) )
      allocate (     td3_h1(lon_dim_h,levs  ,lats_dim_h) )
      allocate (     td3_h2(lon_dim_h,levs  ,lats_dim_h) )
      allocate ( dpdt_d3_h1(lon_dim_h,levs,lats_dim_h) )
      allocate ( dpdt_d3_h2(lon_dim_h,levs,lats_dim_h) )

      if(hybrid) then !--------------  hybrid  ----------------------

        do lan=1,lats_node_a   !sela begin lan loop 3
!         timer1   = rtc()
          lat      = global_lats_a(ipt_lats_node_a+lats_node_a-lan)
          lons_lat = lonsperlat(lat)
          ylan     = yhalo + lan
          rcs2_loc =  rcs2_a(min(lat,latg-lat+1))

!$omp parallel do private(n,k,kk,i)
          do n=1,ntrac
            do k=1,levs
              kk = ksrg + n*levs - k
              do i=1,lons_lat
                rgt_h(xhalo+i,k,ylan,n) = for_grid_a(i,kk,lan)
              enddo
            enddo
          enddo


!$omp parallel do schedule(dynamic,1) private(lon)
!$omp+private(njeff,iblk,xlon,nvcn,xvcn)

          do lon=1,lons_lat,ngptcd
            njeff = min(ngptcd,lons_lat-lon+1)
            iblk  = (lon-1)/ngptcd + 1
            xlon  = xhalo + lon

            call gfidi_hyb_resonan
     &             (lon_dim_a-2, lon_dim_h, njeff, lat,
     &        for_grid_a( lon,ksd   , lan ),
     &        for_grid_a( lon,kst   , lan ),
     &              ug_h(xlon,1     ,ylan),
     &              vg_h(xlon,1     ,ylan),
     &              ug_m(xlon,1     ,ylan),
     &              vg_m(xlon,1     ,ylan),
     &              ww_h(xlon,1     ,ylan),
     &              ww_m(xlon,1     ,ylan),
!    &        for_grid_a( lon,ksrg  , lan),
     &        for_grid_a( lon,kspphi, lan),
     &        for_grid_a( lon,ksplam, lan),
     &        for_grid_a( lon,ksq   , lan),
     &            rcs2_loc,
     &                   etaint,
     &            spdlat(1,iblk),
     &            deltim,nvcn,xvcn,
     &        dyn_grid_a( lon,kdtphi, lan),
     &        dyn_grid_a( lon,kdtlam, lan),
     &             td3_h1(xlon,1     ,ylan),
     &             td3_h2(xlon,1     ,ylan),
     &             rgt_h (xlon,1     ,ylan,1),
     &             ud3_h1(xlon,1     ,ylan),
     &             ud3_h2(xlon,1     ,ylan),
     &             vd3_h1(xlon,1     ,ylan),
     &             vd3_h2(xlon,1     ,ylan),
     &               sd_m( lon       , lan),
     &                z_m( lon,1     , lan),
     &               unlm( lon,1     , lan),
     &               vnlm( lon,1     , lan),
     &               tnlm( lon,1     , lan),
     &               unla( lon,1     , lan),
     &               vnla( lon,1     , lan),
     &               tnla( lon,1     , lan),
     &             dpdt_a( lon       , lan),
     &         dpdt_d3_h1(xlon,1     ,ylan),
     &         dpdt_d3_h2(xlon,1     ,ylan),
     &         grad_gzlam( lon       , lan),
     &         grad_gzphi( lon       , lan),
     &            gz_grid( lon       , lan),
     &             go_forward,
     &             batah,gg_tracers,tb_k,kdt,r_rt,omega_v(lan))
          enddo   ! omp lon
!
          iblk = 1
          do lon=1,lons_lat,ngptcd
             do k=1,levs
               spdmax_node(k) = max(spdmax_node(k),spdlat(k,iblk))
             enddo
            iblk = iblk + 1
          enddo

!         timer2 = rtc()
!         global_times_a(lat,me+1)=global_times_a(lat,me+1)+timer2-timer1

        enddo !sela finshed lan loop 3

      else  !-------------  sigma branch ---------------------------

!      It appears that the sigma option is not coded - Moorthi
!      -------------------------------------------------------
!!$omp parallel do schedule(dynamic,1) private(lon,njeff,iblk)
!         do lon=1,lons_lat,ngptcd
!           njeff = min(ngptcd,lons_lat-lon+1)
!           iblk  = (lon-1)/ngptcd+ 1
!         enddo   !lon
!         write(0,*)' Sigma branch option not coded for SL Dynamics'

      endif ! -----------------------  hybrid if ------------------

      deallocate ( for_grid_a, dyn_grid_a, spdlat )


      call set_halos(ug_h,       lats_nodes_h, global_lats_h,
     &               lonsperlat, lon_dim_h,    xhalo, yhalo,
     &               levs, levs)
      call set_halos(vg_h,       lats_nodes_h, global_lats_h,
     &               lonsperlat, lon_dim_h,    xhalo, yhalo,
     &               levs, levs)
      call set_halos(ww_h,       lats_nodes_h, global_lats_h,
     &               lonsperlat, lon_dim_h,    xhalo, yhalo,
     &               levs, levs)
 
      call set_halos(ud3_h1,     lats_nodes_h, global_lats_h,
     &               lonsperlat, lon_dim_h,    xhalo, yhalo,
     &               levs, levs)
      call set_halos(ud3_h2,     lats_nodes_h, global_lats_h,
     &               lonsperlat, lon_dim_h,    xhalo, yhalo,
     &               levs, levs)
      call set_halos(vd3_h1,     lats_nodes_h, global_lats_h,
     &               lonsperlat, lon_dim_h,    xhalo, yhalo,
     &               levs, levs)
      call set_halos(vd3_h2,     lats_nodes_h, global_lats_h,
     &               lonsperlat, lon_dim_h,    xhalo, yhalo,
     &               levs, levs)
      call set_halos(td3_h1,     lats_nodes_h, global_lats_h,
     &               lonsperlat, lon_dim_h,    xhalo, yhalo,
     &               levs, levs)
      call set_halos(td3_h2,     lats_nodes_h, global_lats_h,
     &               lonsperlat, lon_dim_h,    xhalo, yhalo,
     &               levs, levs)

!     if (me == 1) then
!     write(1000+me,*)' SPHUM BEF set_halos for me=',me,' rgt_h=',
!    & rgt_h(xhalo+1,1:64,yhalo+lan,1)*1000
!     write(1000+me,*)' OZONE BEF set_halos for me=',me,' rgt_h=',
!    & rgt_h(xhalo+1,1:64,yhalo+lan,2)*1000000
!    &,' lon_dim_h=',lon_dim_h,' lats_dim=',lats_dim_a
!     endif

      do ntr=1,ntrac
        call set_halos(rgt_h(1,1,1,ntr), lats_nodes_h, global_lats_h,
     &                 lonsperlat,       lon_dim_h,    xhalo,yhalo,
     &                 levs, levs)
      enddo

      call set_halos(dpdt_d3_h1, lats_nodes_h, global_lats_h,
     &               lonsperlat, lon_dim_h,    xhalo, yhalo,
     &               levs, levs)
      call set_halos(dpdt_d3_h2, lats_nodes_h, global_lats_h,
     &               lonsperlat, lon_dim_h,    xhalo, yhalo,
     &               levs, levs)

200   format(' n.hemi me=',i4,2x,' n_ggp=',i4,'lons_lat=',i4)
201   format(' s.hemi me=',i4,2x,' s_ggp=',i4,'lons_lat=',i4)

      if (me == 0) then
        n_ggp    = lats_nodes_h(1) - yhalo
        lons_lat = lonsperlat(latg)

        if (kdt == 1) print 200,me,n_ggp,lons_lat

        do ntr=1,ntrac
          call sltini_nh_scalar(n_ggp, lon_dim_h, rgt_h(1,1,1,ntr),
     &                          lons_lat)
        enddo
        call sltini_n_uvw(n_ggp, lon_dim_h, ug_h, vg_h, lons_lat)
        call sltini_nh_scalar(n_ggp, lon_dim_h, ww_h, lons_lat)

        call sltini_n_uvw(n_ggp, lon_dim_h, ud3_h1, vd3_h1, lons_lat)
        call sltini_n_uvw(n_ggp, lon_dim_h, ud3_h2, vd3_h2, lons_lat)

        call sltini_nh_scalar(n_ggp, lon_dim_h, td3_h1(1,1,1), lons_lat)
        call sltini_nh_scalar(n_ggp, lon_dim_h, td3_h2(1,1,1), lons_lat)
        call sltini_nh_scalar(n_ggp, lon_dim_h, dpdt_d3_h1(1,1,1),
     &                        lons_lat)
        call sltini_nh_scalar(n_ggp, lon_dim_h, dpdt_d3_h2(1,1,1),
     &                        lons_lat)
      endif
      if (me == nodes-1) then
        s_ggp    = yhalo + 1
        lons_lat = lonsperlat(1)
        if (kdt  == 1) print 201,me,s_ggp,lons_lat

        do ntr=1,ntrac
          call sltini_sh_scalar(s_ggp, lon_dim_h, rgt_h(1,1,1,ntr),
     &                          lons_lat)
        enddo
        call sltini_s_uvw(s_ggp, lon_dim_h, ug_h, vg_h, lons_lat)
        call sltini_sh_scalar(s_ggp, lon_dim_h, ww_h, lons_lat)

        call sltini_s_uvw(s_ggp, lon_dim_h, ud3_h1, vd3_h1, lons_lat)
        call sltini_s_uvw(s_ggp, lon_dim_h, ud3_h2, vd3_h2, lons_lat)

        call sltini_sh_scalar(s_ggp, lon_dim_h, td3_h1(1,1,1), lons_lat)
        call sltini_sh_scalar(s_ggp, lon_dim_h, td3_h2(1,1,1), lons_lat)
        call sltini_sh_scalar(s_ggp, lon_dim_h, dpdt_d3_h1(1,1,1),
     &                        lons_lat)
        call sltini_sh_scalar(s_ggp, lon_dim_h, dpdt_d3_h2(1,1,1),
     &                        lons_lat)
      endif


!sela if(kdt.eq.1) then
!sela do j=1,lats_dim_h   !  lan loop scan
!sela
!sela   k = j
!sela   print*, 'lat,k ud3_h = ',lat,k
!sela   print*,'m ud3_h = ',k,minval(ud3_h(:,:,k)),maxval(ud3_h(:,:,k))
!sela   print*,'m vd3_h = ',k,minval(vd3_h(:,:,k)),maxval(vd3_h(:,:,k))
!sela   print*,'m td3_h = ',k,minval(td3_h(:,:,k)),maxval(td3_h(:,:,k))
!sela   print*,'m dpdt_d3_h = ',k,minval(dpdt_d3_h(:,:,k)),
!sela.          maxval(dpdt_d3_h(:,:,k))
!sela   print*,'m ug_h = ',k,minval(ug_h(:,:,k)),maxval(ug_h(:,:,k))
!sela   print*,'m vg_h = ',k,minval(vg_h(:,:,k)),maxval(vg_h(:,:,k))
!sela enddo
!sela stop
!sela endif
      
      if (me == 0 .and. kdt < 10) write(0,*)'iter_slg = ',iter_slg 

      rgt_a = -9999.0

      allocate (     ud3_a(lonf,levs,lats_dim_a) )
      allocate (     vd3_a(lonf,levs,lats_dim_a) )
      allocate (     td3_a(lonf,levs,lats_dim_a) )
      allocate ( dpdt_d3_a(lonf,levs,lats_dim_a) )

      do j=1,lats_node_a   !  lan loop scan
        lat      = latg+1-global_lats_a(ipt_lats_node_a+lats_node_a-j)
        lons_lat = lonsperlat(lat)

!$omp parallel do schedule(dynamic,1) private(lon,njeff,i_1,i_2)
        do lon=1,lons_lat,ngptcs
          njeff = min(ngptcs,lons_lat-lon+1)
          i_1   = lon
          i_2   = lon+njeff-1

          call slgscan_h(i_1, i_2, s_lat_in_pe, n_lat_in_pe,
     &                   j, lat, lon_dim_h,
     &                   deltim, iter_slg, etamid, lats_node_a,
     &                   ug_h, vg_h, ww_h,
     &                   ug_m, vg_m, ww_m,
     &                   lammp_a, phimp_a, sigmp_a, me,
     &                   rgt_h(1,1,1,1),    rgt_a(1,1,1,1),
     &                   ud3_h1(1,1,1),     ud3_h2(1,1,1), ud3_a(1,1,1),
     &                   vd3_h1(1,1,1),     vd3_h2(1,1,1), vd3_a(1,1,1),
     &                   td3_h1(1,1,1),     td3_h2(1,1,1), td3_a(1,1,1),
     &                   dpdt_d3_h1(1,1,1), dpdt_d3_h2(1,1,1),
     &                   dpdt_d3_a(1,1,1),  global_lats_a, lonsperlat,
     &                   ini_slg, ini_dp,   lprnt)
        enddo   ! lon

      enddo   !  finished lan loop scan

      deallocate ( ug_h,   vg_h,   ww_h,   ud3_h1, ud3_h2, vd3_h1,
     &             vd3_h2, td3_h1, td3_h2, dpdt_d3_h1, dpdt_d3_h2 )

      allocate ( bak_coef_a(lon_dim_a,lota,lats_dim_a)
     &,          bak_grid_a(lonf ,lota) )

      do lan=1,lats_node_a                      !sela begin lan loop 5
         jtem     = lats_node_a+1-lan
         lat      = latg + 1 - global_lats_a(ipt_lats_node_a-1+lan)
         lons_lat = lonsperlat(lat)

!$omp parallel do private(k,i,ktem)
         do k=1,levs
           ktem = levs+1-k
           do i=1,lons_lat
!sela       bak_grid_a(i,kau-1+k)= -omega_v(lan)+

             bak_grid_a(i,kau-1+k) = (ud3_a(i,k,jtem)
     &                             +  hfdt*unla(i,ktem,jtem))*rcos1(lan)
             bak_grid_a(i,kav-1+k) = (vd3_a(i,k,jtem)
     &                             +  hfdt*vnla(i,ktem,jtem))*rcos1(lan)
             bak_grid_a(i,kat-1+k) = td3_a(i,k,jtem)
     &                             + hfdt*tnla(i,ktem,jtem)

           enddo
         enddo
!$omp parallel do private(i)
         do i=1,lons_lat
           bak_grid_a(i,kap) = dpdt_a(i,jtem)
         enddo
         do k=1,levs
!$omp parallel do private(i)
           do i=1,lons_lat
             bak_grid_a(i,kap) = bak_grid_a(i,kap)
     &                         + dpdt_d3_a(i,levp1-k,jtem)
           enddo
         enddo
 
         call grid_to_four(bak_grid_a(1,kau), bak_coef_a(1,kau,lan),
     &                     lon_dim_a-2, lon_dim_a,lons_lat, 3*levs+1)

         if (.not. gg_tracers .or. lsout) then
           do ntr=1,ntrac
             nn = (ntr-1)*levs
!$omp parallel do private(k,i,kk)
             do k=1,levs
               kk = kar-1+nn+k
               do i=1,lons_lat
                 bak_grid_a(i,kk) = rgt_a(i,k,jtem,ntr)
               enddo
             enddo
           enddo
           call grid_to_four(bak_grid_a(1,kar), bak_coef_a(1,kar,lan),
     &                       lon_dim_a-2, lon_dim_a, lons_lat,levh)
         endif

      enddo                            ! finished lan loop 5

      deallocate ( bak_grid_a, unla, vnla, tnla
     &,            ud3_a, vd3_a, td3_a , dpdt_a, dpdt_d3_a )

!     call synctime()

      call four2fln_gg(lats_dim_a,lota,3*levs+1,bak_coef_a,
     &              ls_nodes,max_ls_nodes,
     &              lats_nodes_a,global_lats_a,lon_dim_a,
     &              lats_node_a,ipt_lats_node_a,
     &              lat1s_a,lon_dim_a,latg,latg2,
     &              trie_ls(1,1,p_w), trio_ls(1,1,p_w),
     &              plnew_a, plnow_a,
     &              ls_node,0,2*levs+1,3*levs+1)

      if(.not.gg_tracers .or. lsout) then
         call four2fln_gg(lats_dim_a,lota,levh,bak_coef_a,
     &              ls_nodes,max_ls_nodes,
     &              lats_nodes_a,global_lats_a,lon_dim_a,
     &              lats_node_a,ipt_lats_node_a,
     &              lat1s_a,lon_dim_a,latg,latg2,
     &              trie_ls(1,1,p_rt), trio_ls(1,1,p_rt),
     &              plnew_a, plnow_a,
     &              ls_node,3*levs+1,1,levh)
      endif
      deallocate ( bak_coef_a )

!!$omp parallel do shared(trie_ls,trio_ls,epse,epso,ls_node)
!!$omp+private(k)

!$omp parallel do private(k)
      do k=1,levs
         call uveodz(trie_ls(1,1,p_w  +k-1), trio_ls(1,1,p_x  +k-1),
     &               trie_ls(1,1,p_uln+k-1), trio_ls(1,1,p_vln+k-1),
     &               epse,epso,ls_node)
         call uvoedz(trio_ls(1,1,p_w  +k-1), trie_ls(1,1,p_x  +k-1),
     &               trio_ls(1,1,p_uln+k-1), trie_ls(1,1,p_vln+k-1),
     &               epse,epso,ls_node)
      enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!sela add tb to Y to complete temp. tendency. WATCH SIGN!!!!!!!!

!$omp parallel do private(k,i,ktem,kk)
      do k=1,levs
        ktem = p_y+k-1
        kk   = levp1-k
        do i=1,len_trie_ls
          trie_ls(i,1,ktem) = trie_ls(i,1,ktem) + tb_k(kk)*phige(i,1)
          trie_ls(i,2,ktem) = trie_ls(i,2,ktem) + tb_k(kk)*phige(i,2)
        enddo
        do i=1,len_trio_ls
          trio_ls(i,1,ktem) = trio_ls(i,1,ktem) + tb_k(kk)*phigo(i,1)
          trio_ls(i,2,ktem) = trio_ls(i,2,ktem) + tb_k(kk)*phigo(i,2)
        enddo
      enddo
!sela add tb to Y to complete temp. tendency. WATCH SIGN!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (r_rt == 0.0) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!$omp parallel do private(k,i,jtem,ktem)
        do k=1,levs
          ktem = p_x  +k-1
          jtem = p_uln+k-1
          do i=1,len_trie_ls
            trie_ls(i,1,ktem) = trie_ls(i,1,jtem)
            trie_ls(i,2,ktem) = trie_ls(i,2,jtem)
          enddo
          do i=1,len_trio_ls
            trio_ls(i,1,ktem) = trio_ls(i,1,jtem)
            trio_ls(i,2,ktem) = trio_ls(i,2,jtem)
          enddo
        enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      else
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!$omp parallel do private(k,i,jtem,ktem)
        do k=1,levs
          ktem = p_x  +k-1
          jtem = p_uln+k-1
          do i=1,len_trie_ls
            trie_ls(i,1,ktem) = trie_ls(i,1,jtem)
     &                        - batah*hfdt * trie_ls(i,1,P_gz)
            trie_ls(i,2,ktem) = trie_ls(i,2,jtem)
     &                        - batah*hfdt * trie_ls(i,2,P_gz)
          enddo
          do i=1,len_trio_ls
            trio_ls(i,1,ktem) = trio_ls(i,1,jtem)
     &                        - batah*hfdt * trio_ls(i,1,P_gz)
            trio_ls(i,2,ktem) = trio_ls(i,2,jtem)
     &                        - batah*hfdt * trio_ls(i,2,P_gz)
          enddo
        enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!$omp parallel do private(k,i,jtem,ktem)
      do k=1,levs
        ktem = p_w  +k-1
        jtem = p_vln+k-1
        do i=1,len_trie_ls
          trie_ls(i,1,ktem) = trie_ls(i,1,jtem)
          trie_ls(i,2,ktem) = trie_ls(i,2,jtem)
        enddo
        do i=1,len_trio_ls
          trio_ls(i,1,ktem) = trio_ls(i,1,jtem)
          trio_ls(i,2,ktem) = trio_ls(i,2,jtem)
        enddo
      enddo
      do locl=1,ls_max_node
              l=ls_node(locl,1)
         jbasev=ls_node(locl,2)
         if ( l == 0 ) then
            n = 0
            do k=1,levs
               trie_ls(indlsev(n,l),1,p_w+k-1) = zero
               trie_ls(indlsev(n,l),2,p_w+k-1) = zero
            enddo
         endif
      end do
      spdmax_nodem = spdmax_node
      call mpi_allreduce(spdmax_nodem,spdmax_nodesm,levs,mpi_r_mpi,
     &                   mpi_max,mc_comp,ierr)

!$omp parallel do private(k)
      do k=1,levs
        spdmax(k) = sqrt(spdmax_nodesm(k))
      enddo

!     n_rgt_a never used - so comment out
!     if (lsout) then
!       do k=1,levs
!         n_rgt_a = 0

!         do lan=1,lats_node_a   !sela begin neg. moist loop
!           lat      = global_lats_a(ipt_lats_node_a-1+lan)
!           lons_lat = lonsperlat(lat)
!           do lon=1,lons_lat
!             if(rgt_a(lon,k,lan,1) < 0.)then
!               n_rgt_a = n_rgt_a + 1
!             endif
!           enddo   ! omp lon
!         enddo     ! fin  lat loop


!       enddo       ! fin  levs loop
!     endif         !lsout


!     call mpi_bcast(spdmax,levs,mpi_real8,0,mc_comp,ierr)

      return
      end
