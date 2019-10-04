!
! ==========================================================================
!
      subroutine gfs_dficoll_dynamics(grid_gr,grid_gr_dfi)
!       
!*** Oct 2009, Jun Wang collect n+1 time step data for digital filter
!
      use gfs_dyn_layout1, only : lats_node_a_max
      use gfs_dyn_resol_def
      use gfs_dyn_dfi_mod, only : gfs_dfi_grid_gr
!
      implicit none
!
      REAL(KIND=KIND_EVOD)  :: grid_gr(lonf,lats_node_a_max,lotgr)
      type(gfs_dfi_grid_gr) :: grid_gr_dfi
!
!local vars
      integer :: kstr,kend,krq
!1: gz
      if(grid_gr_dfi%z_imp==1) then
        kstr=g_gz
        kend=g_gz
        grid_gr_dfi%hs(:,:,1)=grid_gr(:,:,kstr)
!      print*,'in dficoll,hs=',maxval(grid_gr_dfi%hs),
!     &  minval(grid_gr_dfi%hs) 
      endif
!2: ps
      if(grid_gr_dfi%ps_imp==1) then
        kstr=g_zq
        kend=g_zq
        grid_gr_dfi%ps(:,:,1)=grid_gr(:,:,kstr)
!      print*,'in dficoll,ps=',maxval(grid_gr_dfi%ps),
!     &   minval(grid_gr_dfi%ps)   
      endif
!3: t
      if(grid_gr_dfi%temp_imp==1) then
        kstr=g_t
        kend=kstr + levs - 1
        grid_gr_dfi%t(:,:,:)=grid_gr(:,:,kstr:kend)
!      print*,'in dficoll,t=',maxval(grid_gr_dfi%t),minval(grid_gr_dfi%t)
      endif
!4: u
      if(grid_gr_dfi%u_imp==1) then
        kstr=g_u
        kend=kstr + levs - 1
        grid_gr_dfi%u=grid_gr(:,:,kstr:kend)
!      print*,'in dficoll,u=',maxval(grid_gr_dfi%u),minval(grid_gr_dfi%u)
      endif
!5: v
      if(grid_gr_dfi%v_imp==1) then
        kstr=g_v
        kend=kstr + levs - 1
        grid_gr_dfi%v=grid_gr(:,:,kstr:kend)
!      print*,'in dficoll,v=',maxval(grid_gr_dfi%v),minval(grid_gr_dfi%v)
      endif
!6: tracer
!       print *,'in dficoll0,tracer=',size(grid_gr_dfi%tracer,3),
!     &   'ntrac=',ntrac,'levs=',levs,'g_rt=',g_rt
      if(grid_gr_dfi%tracer_imp==1) then
        kstr=g_rt
        kend=kstr + levs*ntrac - 1
        grid_gr_dfi%tracer=grid_gr(:,:,kstr:kend)
!        print*,'in dficoll,q=',maxval(grid_gr_dfi%tracer),
!     &   minval(grid_gr_dfi%tracer)
      endif
!       print *,'in dficoll'
!7: p
      if(grid_gr_dfi%p_imp==1) then
        kstr=g_p
        kend=kstr + levs - 1
        grid_gr_dfi%p=grid_gr(:,:,kstr:kend)
!      print*,'in dficoll,p=',maxval(grid_gr_dfi%p),minval(grid_gr_dfi%p)
      endif
!8: dp
      if(grid_gr_dfi%dp_imp==1) then
        kstr=g_dp
        kend=kstr + levs - 1
       grid_gr_dfi%dp=grid_gr(:,:,kstr:kend)
!      print*,'in dficoll,p=',maxval(grid_gr_dfi%p),minval(grid_gr_dfi%p)
      endif
!9: dpdt
      if(grid_gr_dfi%dpdt_imp==1) then
        kstr=g_dpdt
        kend=kstr + levs - 1
        grid_gr_dfi%dpdt=grid_gr(:,:,kstr:kend)
!      print*,'in dficoll,p=',maxval(grid_gr_dfi%dpdt),
!     &    minval(grid_gr_dfi%dpdt)
      endif

!       print *,'end of dfi coll'
!
      end subroutine gfs_dficoll_dynamics
