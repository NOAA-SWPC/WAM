          Module idea_ncout_phys
!
! fixed 3D fields
!         use layout1                            !lats_node_r, lats_node_r_max
          use resol_def, only : lonr, latr, levs !lonf, latg
!
! accumulation id performed at lonr, lats_node_r  with grid_fld%ps( lonr, lats_node_rmax)
!
          integer                              :: wam_daysteps
          real,allocatable, dimension(:, :)    ::  wps, wphis
          real,allocatable, dimension(:, :, :) :: wU, wV, wT
          real,allocatable, dimension(:, :, :) :: wQ, wO3, wOP, wO2
          real,allocatable, dimension(:, :, :) :: Wp3d, Wz3d, Wdp
          real,allocatable, dimension(:, :, :) :: 
     & u24s, u24c, u12s, u12c, u08s, u08c
          real,allocatable, dimension(:, :, :) ::
     &  v24s, v24c, v12s, v12c, v08s, v08c
          real,allocatable, dimension(:, :,:) ::
     &  t24s, t24c, t12s, t12c, t08s, t08c
          real,allocatable, dimension(:, :, :) ::
     &  oz24s, oz24c, oz12s, oz12c, oz08s, oz08c
          real,allocatable, dimension(:, :, :) ::
     &  op24s, op24c, op12s, op12c, op08s, op08c
          real,allocatable, dimension(:, :, :) ::
     &  om24s, om24c, om12s, om12c, om08s, om08c
!
          real,allocatable, dimension(:, :, :)    ::  w3d
          real,allocatable, dimension(:, :)       ::  w2d
!
!
          real, allocatable, dimension(:)   ::
     &   ws24,  wc24, ws12, wc12, ws8, wc8
          real                              ::  Rfd   ! for daily mean values
          real                              ::  Rfd2  !  2*rfd

!===============================================================================
! Total # of fields: 7-daily
!                  : 6-tidal x 6 =36
!
! 43- [3D ararys] + 2- [2D arrays]
!          gfs_physics_gridgr_mod.f:grid_fld%ps(dim1,dim2),  grid_fld%t(dim1,dim2,dim3),
!
!        [lonr, lats_node_r_max, levs]=(dim1,dim2,dim3)
!
!         grid_fld%tracers(n)%flds(dim1,dim2,dim3)
!
!================================================================================
          contains
!
          subroutine idea_init_daily(deltim)
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
!          deallocate( ws24, wc24, ws12, wc12, ws8, wc8 )
!          real, dimension(lonr, latr) ::  wps, wphis

           allocate(wps(lonr, n1), wphis(lonr,n1), w2d(lonr,n1) )
!          real, dimension(lonr, latr, levs) ::  wU, wV, wT
           allocate(wu(lonr, n1, levs), wv(lonr,n1,levs))
           allocate(w3d(lonr, n1, levs),  wT(lonr,n1,levs))
           allocate(wQ(lonr, n1, levs), wO3(lonr,n1,levs))
           allocate(wOP(lonr, n1, levs), wO2(lonr,n1,levs))
           allocate(wdP(lonr, n1, levs), wP3d(lonr,n1,levs))   ! pressures
           allocate(wZ3d(lonr, n1, levs))                      ! geo-height

!          real, dimension(lonr, latr, levs) ::  wQ, wO3, wOP, wO2
!          real, dimension(lonr, latr, levs) ::  u24s, u24c, u12s, u12c, u08s, u08c
           allocate(u24s(lonr, n1, levs), u24c(lonr,n1,levs))
           allocate(v24s(lonr, n1, levs), v24c(lonr,n1,levs))
           allocate(t24s(lonr, n1, levs), t24c(lonr,n1,levs))

           allocate(u12s(lonr, n1, levs), u12c(lonr,n1,levs))
           allocate(v12s(lonr, n1, levs), v12c(lonr,n1,levs))
           allocate(t12s(lonr, n1, levs), t12c(lonr,n1,levs))

           allocate(u08s(lonr, n1, levs), u08c(lonr,n1,levs))
           allocate(v08s(lonr, n1, levs), v08c(lonr,n1,levs))
           allocate(t08s(lonr, n1, levs), t08c(lonr,n1,levs))

!
          end subroutine idea_init_daily

          end Module idea_ncout_phys
!
       subroutine idea_make_daily(grid_fld, rih24, ind_lats, ind_lons)
!23456
          use layout1, only : lats_node_r, lats_node_r_max,  me
          use resol_def, only : lonr, latr, levs             !lonf, latg
          use gfs_physics_gridgr_mod, only : Grid_Var_Data
          use  idea_ncout_phys
          
          implicit none
!input
          TYPE(Grid_Var_Data)     :: grid_fld     ! dims: lonr, lats_node_r_max, levs, gis_phy%ntrac in gfs_physics_initialize_mod.f
!          integer                 :: wam_daysteps
          integer, dimension(lats_node_r) :: ind_lats, ind_lons
          real                    :: rih24
!locals
          integer :: n1, n2, ntk      
!
          integer ::  i,j,k         
          n1 = lonr
          n2 = lats_node_r
          ntk = wam_daysteps +1
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
!       print *,  'VAY-grid_fld%ps ', me, maxval(w2d), minval(w2d) 
       if (wam_daysteps .gt. 0 ) then
!$omp parallel do private(i,j)
      do j=1,n2
        do i=1,lonr
!              w2d(i,j) = grid_fld%ps(i, j)
        wps(i,j)=wps(i,j) +grid_fld%ps(i, j)*rfd
        wphis(i,j)=wphis(i,j) +grid_fld%z(i, j)*rfd
        enddo
      enddo


!$omp parallel do private(i,j,k)
      do k=1, levs
!
      do j=1,n2
        do i=1,lonr                                     !ind_lons(j)
        wu(i,j,k)=wu(i,j,k) +grid_fld%U(i, j,k)*rfd
        wv(i,j,k)=wv(i,j,k) +grid_fld%V(i, j,k)*rfd
        wt(i,j,k)=wt(i,j,k) +grid_fld%T(i, j,k)*rfd
        wQ(i,j,k)=wQ(i,j,k)+grid_fld%tracers(1)%flds(i,j,k)*rfd
        wO3(i,j,k)=wO3(i,j,k)+grid_fld%tracers(2)%flds(i,j,k)*rfd
        wOP(i,j,k)=wOP(i,j,k)+grid_fld%tracers(4)%flds(i,j,k)*rfd
        wO2(i,j,k)=wO2(i,j,k)+grid_fld%tracers(5)%flds(i,j,k)*rfd
!
! p-dp-dpdt
!
        wdp(i,j,k)=wdp(i,j,k)+grid_fld%dp(i, j,k)*rfd
        wp3d(i,j,k)=wp3d(i,j,k)+grid_fld%p(i, j,k)*rfd
        enddo
      enddo
!U
! 2./Ntime = 2*rdf
       w2d(1:lonr, 1:n2) =grid_fld%U(1:lonr, 1:n2,k)
            call Get_tidal_decsum(w2d, u24s(:,:,k), u24c(:,:,k), 
     & u12s(:,:,k), u12c(:,:,k), u08s(:,:,k), u08c(:,:,k), n1, n2, ntk)
!V
       w2d(1:lonr, 1:n2) =grid_fld%V(1:lonr, 1:n2,k)
!
            call Get_tidal_decsum(w2d, v24s(:,:,k), v24c(:,:,k), 
     & v12s(:,:,k), v12c(:,:,k), v08s(:,:,k), v08c(:,:,k), n1, n2, ntk)

       w2d(1:lonr, 1:n2)=grid_fld%T(1:lonr, 1:n2,k)

            call Get_tidal_decsum(w2d, T24s(:,:,k), T24c(:,:,k), 
     & t12s(:,:,k), t12c(:,:,k), t08s(:,:,k), t08c(:,:,k), n1, n2, ntk)
           enddo

!               wam_daysteps = wam_daysteps+1
         if (me==0) print *, ' VAY-sumtides'      
         ELSE
!$omp parallel do private(i,j)
      do j=1,n2
        do i=1,lonr
        wps(i,j)= grid_fld%ps(i, j)*rfd
        wphis(i,j)=grid_fld%z(i, j)*rfd
        enddo
      enddo
!$omp parallel do private(i,j,k)
      do k=1, levs
!
      do j=1,n2
        do i=1,lonr
        wu(i,j,k)=grid_fld%U(i, j,k)*rfd
        wv(i,j,k)=grid_fld%V(i, j,k)*rfd
        wt(i,j,k)=grid_fld%T(i, j,k)*rfd
        wQ(i,j,k)=grid_fld%tracers(1)%flds(i,j,k)*rfd
        wO3(i,j,k)=grid_fld%tracers(2)%flds(i,j,k)*rfd
        wOP(i,j,k)=grid_fld%tracers(4)%flds(i,j,k)*rfd
        wO2(i,j,k)=grid_fld%tracers(5)%flds(i,j,k)*rfd
        wdp(i,j,k)=grid_fld%dp(i, j,k)*rfd
        wp3d(i,j,k)=grid_fld%p(i, j,k)*rfd
        enddo
      enddo
!
       w2d(1:lonr, 1:n2) =grid_fld%U(1:lonr, 1:n2,k)
       call Get_tidal_dec(w2d, u24s(:,:,k), u24c(:,:,k), 
     & u12s(:,:,k), u12c(:,:,k), u08s(:,:,k), u08c(:,:,k), n1, n2, ntk)
       w2d(1:lonr, 1:n2) =grid_fld%V(1:lonr, 1:n2,k)
       call Get_tidal_dec(w2d, V24s(:,:,k), V24c(:,:,k), 
     & V12s(:,:,k), V12c(:,:,k), V08s(:,:,k), V08c(:,:,k), n1, n2, ntk)
       w2d(1:lonr, 1:n2) =grid_fld%T(1:lonr, 1:n2,k)
      call Get_tidal_dec(w2d, T24s(:,:,k), T24c(:,:,k), 
     & t12s(:,:,k), t12c(:,:,k), t08s(:,:,k), t08c(:,:,k), n1, n2, ntk)
!         
           enddo
! 
!         if (me==0)  print *, ' VAY-initides'      
!            wam_daysteps = 1


          ENDIF

      END subroutine idea_make_daily
!
      Subroutine Get_tidal_dec(T, T24s, T24c, t12s, t12c,
     &    t08s, t08c, n1, n2, ntk)
      use  idea_ncout_phys, only : ws24, wc24, ws12, wc12, ws8, wc8
          implicit none
          integer :: n1, n2, ntk
       real, dimension(n1, n2)::T, T24s, T24c, t12s, t12c, t08s, t08c
          T24s = ws24(ntK) * T
          T24c = wc24(ntK) * T

          T12s = ws12(ntK) * T
          T12c = wc12(ntK) * T

          T08s = ws8(ntK) * T
          T08c = wc8(ntK) * T
!
          END Subroutine Get_tidal_dec
!
      Subroutine Get_tidal_decsum(T, T24s, T24c, t12s, t12c, 
     &           t08s, t08c, n1, n2, ntk)
      use  idea_ncout_phys, only : ws24, wc24, ws12, wc12, ws8, wc8
      implicit none
      integer :: n1, n2, ntk
      real, dimension(n1, n2) ::  T,T24s, T24c, t12s, t12c, t08s, t08c
          T24s = ws24(ntK) * T + T24s 
          T24c = wc24(ntK) * T + T24c 

          T12s = ws12(ntK) * T + T12s
          T12c = wc12(ntK) * T + T12c

          T08s = ws8(ntK) * T + T08s
          T08c = wc8(ntK) * T + T08c
!
          END Subroutine Get_tidal_decsum
!
