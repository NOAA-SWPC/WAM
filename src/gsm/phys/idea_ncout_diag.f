          Module idea_ncout_diag
!!
! fixed 3D fields
!         use layout1                            !lats_node_r, lats_node_r_max
          use resol_def, only : lonr, latr, levs !lonf, latg
!
! accumulation id performed at lonr, lats_node_r  with grid_fld%ps( lonr, lats_node_rmax)
!
          integer                              :: wamd_daysteps
          real,allocatable, dimension(:, :)    :: d_wps, d_wphis
!10-3d
          real,allocatable, dimension(:, :, :) :: d_wU, d_wV,  d_wT
          real,allocatable, dimension(:, :, :) :: d_dp, d_p3d, d_z3d
          real,allocatable, dimension(:, :, :) :: wAxz, wAyz, wEps, wKed
!
          real, allocatable, dimension(:)   ::
     &   ws24,  wc24, ws12, wc12, ws8, wc8
!
          real,allocatable, dimension(:, :, :)    ::  w3d
          real,allocatable, dimension(:, :)       ::  w2d
!
!
          real                              ::  Rfd   ! for daily mean values
          real                              ::  Rfd2  !  2*rfd

!
!================================================================================
          contains
!
          subroutine idea_init_diagav(deltim)
          use resol_def, only : lonr, latr, levs
          use layout1,only : n1 => lats_node_r, nmax=>lats_node_r_max
          implicit none   
          real, intent(in)  :: deltim

          integer ::  NT
          real    :: fr24, fr12, fr8, pi2day
          integer :: k
          NT = Nint(86400./deltim)              !480 points
          rfd = deltim/86400.
          rfd2 = 2.*rfd
          pi2day = 8.*atan(1.0)*rfd            ! 2*pi/Tday*dt
          fr24 = 1.*pi2day
          fr12 = 2.*pi2day
          fr8  = 3.*pi2day
!
          allocate( ws24(nt), wc24(nt), ws12(nt), wc12(nt))
          allocate( ws8(nt), wc8(nt) )
!
          do k=1, nt
            ws24(k)= sin(fr24*K)*rfd2 
            wc24(k)= cos(fr24*K)*rfd2 
            ws12(k)= sin(fr12*K)*rfd2 
            wc12(k)= cos(fr12*K)*rfd2 
             ws8(k)= sin(fr8*K)*rfd2 
             wc8(k)= cos(fr8*K)*rfd2 
          enddo
!
!  

           allocate(d_wps(lonr, n1), d_wphis(lonr,n1), w2d(lonr,n1) )
!         
           allocate(wAxz(lonr, n1, levs), wAyz(lonr,n1,levs))
           allocate(w3d(lonr, n1, levs),  wEps(lonr,n1,levs))
           allocate(wKed(lonr, n1, levs), d_wu(lonr,n1,levs))
           allocate(d_wv(lonr, n1, levs), d_wT(lonr,n1,levs))
!
           allocate(d_dP(lonr, n1, levs), d_P3d(lonr,n1,levs))   ! pressures
           allocate(d_Z3d(lonr, n1, levs))                       ! geo-height

!          real, dimension(lonr, latr, levs) ::  wQ, wO3, wOP, wO2
!          real, dimension(lonr, latr, levs) ::  u24s, u24c, u12s, u12c, u08s, u08c
!           allocate(u24s(lonr, n1, levs), u24c(lonr,n1,levs))
!           allocate(v24s(lonr, n1, levs), v24c(lonr,n1,levs))
!           allocate(t24s(lonr, n1, levs), t24c(lonr,n1,levs))

!           allocate(u12s(lonr, n1, levs), u12c(lonr,n1,levs))
!           allocate(v12s(lonr, n1, levs), v12c(lonr,n1,levs))
!           allocate(t12s(lonr, n1, levs), t12c(lonr,n1,levs))


!
          end subroutine idea_init_diagav

          end MODULE idea_ncout_diag
!
       subroutine idea_make_diagav
     &      (grid_fld, rih24, ind_lats, ind_lons)
!
!23456
!
          use layout1, only : lats_node_r, lats_node_r_max,  me
          use resol_def, only : lonr, latr, levs                  !lonf, latg
!          use  wam_diag3d_mod only : G3D_WAMD
          use  idea_ncout_diag
          use gfs_physics_gridgr_mod, only : Grid_Var_Data         
          use wam_pass_diag_types,    only : Gis_wam
          implicit none
!input
!          TYPE(G3D_WAMD)          :: gis_wam     !
          TYPE(Grid_Var_Data)     :: grid_fld    ! regular state
!          integer                 :: wam_daysteps
          integer, dimension(lats_node_r) :: ind_lats, ind_lons
          real                    :: rih24
!locals
          integer :: n1, n2, ntk      
!
          integer ::  i,j,k         
          n1 = lonr
          n2 = lats_node_r
          ntk = wamD_daysteps
!
!       print *, ' VAY-NCOUT ', n1, n2, ntk  
!       print *, ' VAY-NCOUT ', maxval(grid_fld%ps), minval(grid_fld%ps)
!
!       if (rih24.ne.0.0) then
!
!$omp parallel do private(i,j)
         do j=1,n2
            do i=1,lonr                   !ind_lons(j)
              w2d(i,j) = grid_fld%ps(i, j)
            enddo
         enddo      
!        
! 
       if (wamD_daysteps .gt. 0 ) then

!$omp parallel do private(i,j)
      do j=1,n2
        do i=1,lonr
        d_wps(i,j)=d_wps(i,j) +grid_fld%ps(i, j)*rfd
        d_wphis(i,j)=d_wphis(i,j) +grid_fld%z(i, j)*rfd
        enddo
      enddo


!$omp parallel do private(i,j,k)
      do k=1, levs
!
      do j=1,n2
        do i=1,lonr                                     !ind_lons(j)
        d_wu(i,j,k)=d_wu(i,j,k) +grid_fld%U(i, j,k)*rfd
        d_wv(i,j,k)=d_wv(i,j,k) +grid_fld%V(i, j,k)*rfd
        d_wt(i,j,k)=d_wt(i,j,k) +grid_fld%T(i, j,k)*rfd
!
         wAxz(i,j,k)=wAxz(i,j,k)+gis_wam%daxz(i, j,k)*rfd
         wAyz(i,j,k)=wAyz(i,j,k)+gis_wam%dayz(i, j,k)*rfd
         wEps(i,j,k)=weps(i,j,k)+gis_wam%deps(i, j,k)*rfd
         wKed(i,j,k)=wked(i,j,k)+gis_wam%dked(i, j,k)*rfd
!
!  wAxz, wAyz, wEps, wKed
!        wQ(i,j,k)=wQ(i,j,k)+grid_fld%tracers(1)%flds(i,j,k)*rfd
!        wO3(i,j,k)=wO3(i,j,k)+grid_fld%tracers(2)%flds(i,j,k)*rfd
!        wOP(i,j,k)=wOP(i,j,k)+grid_fld%tracers(4)%flds(i,j,k)*rfd
!        wO2(i,j,k)=wO2(i,j,k)+grid_fld%tracers(5)%flds(i,j,k)*rfd
!
! p-dp-dpdt
!
        d_dp(i,j,k)=d_dp(i,j,k)+grid_fld%dp(i, j,k)*rfd
        d_p3d(i,j,k)=d_p3d(i,j,k)+grid_fld%p(i, j,k)*rfd
        d_z3d(i,j,k)=d_z3d(i,j,k)+gis_wam%zgkm(i, j,k)*rfd	
        enddo
      enddo
      ENDDO    !k=1,levs
!               wam_daysteps = wam_daysteps+1
!         if (me==0) print *, ' VAY-sumtides'      
         ELSE
!$omp parallel do private(i,j)
      do j=1,n2
        do i=1,lonr
        d_wps(i,j)= grid_fld%ps(i, j)*rfd
        d_wphis(i,j)=grid_fld%z(i, j)*rfd
        enddo
      enddo
!$omp parallel do private(i,j,k)
      do k=1, levs
!
      do j=1,n2
        do i=1,lonr

        d_wu(i,j,k)=grid_fld%U(i, j,k)*rfd
        d_wv(i,j,k)=grid_fld%V(i, j,k)*rfd
        d_wt(i,j,k)=grid_fld%T(i, j,k)*rfd
         wAxz(i,j,k)=gis_wam%daxz(i, j,k)*rfd
         wAyz(i,j,k)=gis_wam%dayz(i, j,k)*rfd
         wEps(i,j,k)=gis_wam%deps(i, j,k)*rfd
         wKed(i,j,k)=gis_wam%dked(i, j,k)*rfd
!        wQ(i,j,k)=grid_fld%tracers(1)%flds(i,j,k)*rfd
!        wO3(i,j,k)=grid_fld%tracers(2)%flds(i,j,k)*rfd
!        wOP(i,j,k)=grid_fld%tracers(4)%flds(i,j,k)*rfd
!        wO2(i,j,k)=grid_fld%tracers(5)%flds(i,j,k)*rfd

        d_dp(i,j,k) = grid_fld%dp(i, j,k)*rfd
        d_p3d(i,j,k)= grid_fld%p(i, j,k) *rfd
        d_z3d(i,j,k)= gis_wam%zgkm(i, j,k)*rfd	
        enddo
      enddo
      
!
!       w2d(1:lonr, 1:n2) =grid_fld%U(1:lonr, 1:n2,k)
!       call Get_tidal_dec(w2d, u24s(:,:,k), u24c(:,:,k), 
!     & u12s(:,:,k), u12c(:,:,k), u08s(:,:,k), u08c(:,:,k), n1, n2, ntk)
!       w2d(1:lonr, 1:n2) =grid_fld%V(1:lonr, 1:n2,k)
!       call Get_tidal_dec(w2d, V24s(:,:,k), V24c(:,:,k), 
!     & V12s(:,:,k), V12c(:,:,k), V08s(:,:,k), V08c(:,:,k), n1, n2, ntk)
!       w2d(1:lonr, 1:n2) =grid_fld%T(1:lonr, 1:n2,k)
!      call Get_tidal_dec(w2d, T24s(:,:,k), T24c(:,:,k), 
!     & t12s(:,:,k), t12c(:,:,k), t08s(:,:,k), t08c(:,:,k), n1, n2, ntk)
!         
       ENDDO 


         if (me==0)  print *, ' VAY-init-DIAG'      
!        


      ENDIF

      END subroutine idea_make_diagav
!
