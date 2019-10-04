      subroutine wrtout_dynamics_slg_nemsio(fhour,
     &                                      for_grid_a,ug_h,vg_h,
     &                                      global_lats_a,lonsperlat,
     &                                      kst, ksq, ksrg, gz_grid)

! Dec 2014 Weiyu Yang - modified for NEMSIO in the Semi-Lagrangian GSM model.
!----------------------------------------------------------------------------
!!
!! write out only grid values for nemsio
!!
      use gfs_dyn_resol_def,      ONLY : lonf, levs, levh, latg
      use gfs_dyn_layout1,        ONLY : lats_node_a, lats_dim_a,
     &                                   nodes_comp, me, ipt_lats_node_a
     &                                 , lon_dim_a
      use gfs_dyn_coordinate_def, ONLY : ak5, bk5
      use namelist_dynamics_def,  ONLY : ens_nam
      use gfs_dyn_mpi_def
      use gfs_dyn_physcons,       ONLY : fv => con_fvirt, grav => con_g
      use gfs_dyn_MACHINE,        ONLY : kind_evod, kind_grid
      use layout_lag          ,   ONLY : lats_dim_h, lon_dim_h
      use layout_grid_tracers ,   ONLY : xhalo, yhalo
      use gfs_dyn_gg_def,         ONLY : coslat_a

      implicit none
!
      real(kind=kind_evod), parameter :: gravi=1.0/grav
     &,                                  one=1.0, cb2pa=1000.0
     &,                                  qmin=1.e-10
      real(kind=kind_evod) fhour

      real(kind=kind_evod) :: for_grid_a(lonf,2*levs+levh+3,lats_dim_a),
     &                        ug_h(lon_dim_h,levs,lats_dim_h),
     &                        vg_h(lon_dim_h,levs,lats_dim_h),
     &                        gz_grid(lon_dim_a,lats_dim_a) 

      integer                 kst, ksq, ksrg
     &,                       lat, lan, lanh, lan1, km, ks
     &,                       lons_lat,kk,nfill,jlonf,
     &                        ierr,i,j,k,l,n,ih,ndig,kh,ioproc
      logical lfnhr
!     real(kind=8) t1,t2,t3,t4,t5,ta,tb,tc,td,te,tf,rtc,tx,ty
!
      real(kind=kind_evod), allocatable :: tfac(:,:)
      real(kind=kind_evod)  tx1

      character CFHOUR*40,CFORM*40,filename*255
!!
      INTEGER, dimension(latg) :: GLOBAL_lats_a,   lonsperlat
!
      real(kind=kind_grid), allocatable :: zsg(:,:), psg(:,:)
     &,                       dpg(:,:,:), ttg(:,:,:), uug(:,:,:)
     &,                       vvg(:,:,:), rqg(:,:,:)

      call mpi_barrier(mc_comp,ierr)

      ioproc = nodes_comp - 1

!sela set lfnhr to false for writing one step output etc.
      lfnhr = .true.    ! no output
      lfnhr = 3600*abs(fhour-nint(fhour)) <= 1
      IF(LFNHR) THEN
        KH   = NINT(FHOUR)
        NDIG = MAX(LOG10(KH+0.5)+1.,2.)
        WRITE(CFORM,'("(I",I1,".",I1,")")') NDIG,NDIG
        WRITE(CFHOUR,CFORM) KH
      ELSE
        KS   = NINT(FHOUR*3600)
        KH   = KS/3600
        KM   = (KS-KH*3600)/60
        KS   = KS-KH*3600-KM*60
        NDIG = MAX(LOG10(KH+0.5)+1.,2.)
        WRITE(CFORM,
     &      '("(I",I1,".",I1,",A1,I2.2,A1,I2.2)")') NDIG,NDIG
        WRITE(CFHOUR,CFORM) KH,':',KM,':',KS
      ENDIF
      if( nfill(ens_nam) == 0 ) then
        CFHOUR = CFHOUR(1:nfill(CFHOUR))
      else
        CFHOUR = CFHOUR(1:nfill(CFHOUR)) // ens_nam(1:nfill(ens_nam))
      endif
      if (me == ioproc)
     &write(0,*)' in wrtout_dynamics_slg cfhour=',cfhour,' ens_nam=',
     &  ens_nam,'fhour=',fhour,'lfnhr=',lfnhr
!
      call MPI_BARRIER(mpi_comm_all,ierr)
!
!***                    BUILD STATE ON EACH NODE                ********
! COMP tasks only assemble upper air state first then sfc state,
! then (only if liope)  flux state.
!
      if (.not. allocated(zsg)) allocate(zsg(lonf,lats_node_a))
      if (.not. allocated(psg)) allocate(psg(lonf,lats_node_a))
      if (.not. allocated(rqg)) allocate(rqg(lonf,lats_node_a,levh))
      if (.not. allocated(dpg)) allocate(dpg(lonf,lats_node_a,levs))
      if (.not. allocated(uug)) allocate(uug(lonf,lats_node_a,levs))
      if (.not. allocated(vvg)) allocate(vvg(lonf,lats_node_a,levs))
      if (.not. allocated(ttg)) allocate(ttg(lonf,lats_node_a,levs))

      if(mc_comp /= MPI_COMM_NULL) then
        if (.not. allocated(tfac)) allocate (tfac(lonf,levs))
!$omp parallel do private(lan,lan1,lanh,lat,lons_lat,tx1,i,k,kk,ih)
!$omp+private(tfac)
        do lan=1,lats_node_a
          lan1     = lats_node_a + 1 - lan
          lanh     = lan1 + yhalo
          lat      = global_lats_a(ipt_lats_node_a-1+lan)
          lons_lat = lonsperlat(lat)
          tx1      = one / coslat_a(lat)

          DO k = 1, levh
            kk = ksrg + k - 1
            DO i = 1, lons_lat
              rqg(i,lan,k) = for_grid_a(i,kk,lan1)
            ENDDO
          ENDDO
          DO i = 1, lons_lat
            psg(i, lan) = exp(for_grid_a(i,ksq,lan1))
            zsg(i, lan) = gz_grid(i, lan1) * gravi
          ENDDO
          do k=1,levs
            kk = levs - k + 1
            do i=1,lons_lat
              dpg(i,lan,k) = ak5(kk+1)-ak5(kk)
     &                     + (bk5(kk+1)-bk5(kk)) * psg(i,lan)
            enddo
          enddo
          do i=1,lons_lat
            psg(i,lan) = cb2pa*psg(i,lan)
          enddo
          do k=1,levs
            do i=1,lons_lat
              tfac(i,k)    = one / (one + fv*max(rqg(i,lan,k),qmin))
              dpg(i,lan,k) = cb2pa*dpg(i,lan,k)
            enddo
          enddo

          DO k = 1, levs
            kk = kst + k - 1
            DO i = 1, lons_lat
              ih = i + xhalo
              uug(i,lan,k) = ug_h(ih,k,lanh) * tx1
              vvg(i,lan,k) = vg_h(ih,k,lanh) * tx1
              ttg(i,lan,k) = for_grid_a(i,kk,lan1) * tfac(i,k)
            ENDDO
          ENDDO                 

        enddo
!
!  done with state build
!  NOW STATE IS ASSEMBLED ON EACH NODE.  GET EVERYTHING OFF THE COMPUTE
!  NODES (currently done with a send to the I/O task_
!  send state to I/O task.  All tasks
!
        call grid_collect (zsg,psg,uug,vvg,ttg,rqg,dpg,
     &                     global_lats_a,lonsperlat)
!
      endif
!
      end subroutine wrtout_dynamics_slg_nemsio
