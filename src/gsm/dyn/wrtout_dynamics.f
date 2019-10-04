      INTEGER FUNCTION nfill(C)
      implicit none
      integer j
      CHARACTER*(*) C
      NFILL = LEN(C)
      DO J=1,NFILL
        IF(C(J:J) == ' ') THEN
          NFILL = J-1
          RETURN
        ENDIF
      ENDDO
      RETURN
      END


      module gfs_dyn_write_state
!
! new module to supply domain information to the GFS output routines
!called by wrtout.
!
! May 2009 Jun Wang, modified to use write grid component
! Feb 2011 Henry Juang, modified to have options for mass_dp and ndsl advection
! Oct 2012 Jun Wang, add sigio output option
! Aug 2013 Henry Juang, add sigio output with ndsl
! Sep 2014 S Moorthi - some cleanup and optimization
! May 2017 Weiyu Yang - Modified for the WAM restart.
!
      use gfs_dyn_machine
      use gfs_dyn_resol_def
      implicit none
!
      real(kind=kind_io4), allocatable,target :: buff_mult_pieceg(:,:,:)
      real(kind=kind_io4), allocatable :: buff_mult_piecesg(:)
!
      real(kind=kind_io4), allocatable :: buff_mult_piece(:,:,:),
     &                                    buff_mult_pieces(:,:,:,:)
!     real(kind=kind_io4), allocatable :: buff_mult_piecef(:,:,:),
!    &                                    buff_mult_piecesf(:,:,:,:)
!     real(kind=kind_io4), allocatable :: buff_mult_piecea(:,:,:),
!    &                                    buff_mult_piecesa(:,:,:,:)
!     integer, allocatable :: ivar_global(:),ivar_global_a(:,:)

      integer, allocatable :: ivarg_global(:),ivarg_global_a(:,:)
!
      integer ngrid ,ngrida,ngridg
!     save ngrid,ngrida,buff_mult_piece,buff_mult_pieces,ivar_global
!    &,    ngridg,buff_mult_pieceg,buff_mult_piecesg,ivarg_global

      save ngrid,ngrida,buff_mult_piece,buff_mult_pieces,
     &     ngridg,buff_mult_pieceg,buff_mult_piecesg,ivarg_global
      end module gfs_dyn_write_state

      subroutine wrtout_dynamics(phour,fhour,zhour,idate,
     &                  TRIE_LS,TRIO_LS,grid_gr,
     &                  sl,si,
     &                  ls_node,ls_nodes,max_ls_nodes,
     &                  lats_nodes_a,global_lats_a,lonsperlat,
     &                  snnp1ev,snnp1od,
     &                  epsedn,epsodn,plnev_a,plnod_a,
     &                  epse  ,epso  ,plnew_a,plnow_a,
     &                  pdryini,sigf)
!!
!! write out only grid values for nemsio
!!
      use gfs_dyn_resol_def
      use gfs_dyn_layout1
      use gfs_dyn_coordinate_def
      use namelist_dynamics_def
      use gfs_dyn_mpi_def
      use gfs_dyn_gg_def
      use gfs_dyn_tracer_const
      use gfs_dyn_physcons, cp => con_cp 
     &                    , rd => con_rd, fv => con_fvirt
     &                    , rkappa => con_rocp
      use do_dynamics_mod
      implicit none
!
      real(kind=kind_evod) phour,fhour,zhour
!
      integer              idate(4),km,iostat,no3d,ks
      logical lfnhr
      real lat, lan
      real(kind=8) t1,t2,t3,t4,t5,ta,tb,tc,td,te,tf,rtc,tx,ty
      real timesum
!
      real(kind=kind_evod) sl(levs), si(levp1)
!
      integer              ls_node(ls_dim,3)
      integer              ls_nodes(ls_dim,nodes)
      integer, dimension(nodes) ::  max_ls_nodes, lats_nodes_a

      real(kind=kind_evod), allocatable :: tfac(:,:), sumq(:,:)
     &,                                    tki(:,:)
      real(kind=kind_evod), parameter :: one=1.0, cb2pa=1000.0
     &,                                  qmin=1.e-10
      real(kind=kind_evod)  tx1, tkrt0
      integer               lons_lat,nn,kk,nnl, nosig,nfill,jlonf,
     &                      ierr,i,j,k,l,lenrec,locl,n,node,
     &                      thermodyn_id_out,sfcpress_id_out
      character*5 sigf
      data timesum/0./
cc
      REAL(KIND=KIND_EVOD) TRIE_LS(LEN_TRIE_LS,2,lotls)
     &,                    TRIO_LS(LEN_TRIO_LS,2,lotls)
      REAL(KIND=KIND_grid) grid_gr(lonf*lats_node_a_max,lotgr)
!!
      character CFHOUR*40,CFORM*40,filename*255
      integer jdate(4),nzsig,ndigyr,ndig,kh,ioproc
!  ---  locals for iau:
      integer :: jdate_iau(4),idat(8),jdat(8)
      real (kind=kind_evod) :: zhour_iau,fhour_iau
      real :: rinc(5)
!!
      REAL (KIND=KIND_grid) pdryini
      INTEGER, dimension(latg) :: GLOBAL_lats_a,   lonsperlat
!
      real(kind=kind_evod), dimension(len_trie_ls) ::  epsedn, epse
     &,                                                snnp1ev
      real(kind=kind_evod), dimension(len_trio_ls) ::  epsodn, epso
     &,                                                snnp1od
!!
      real(kind=kind_evod), dimension(len_trie_ls,latg2) ::
     &                        plnev_a, plnew_a
      real(kind=kind_evod), dimension(len_trio_ls,latg2) ::
     &                        plnod_a, plnow_a
!!
      real(kind=kind_grid), allocatable :: zsg(:,:), psg(:,:)
     &,                       dpg(:,:,:), ttg(:,:,:), uug(:,:,:)
     &,                       vvg(:,:,:), rqg(:,:,:)
!!
!     real(kind=kind_mpi),allocatable :: trieo_ls_nodes_buf(:,:,:,:,:)
!     real(kind=kind_mpi),allocatable :: trieo_ls_node(:,:,:)
!     save trieo_ls_nodes_buf,trieo_ls_node
!     real(kind=8) tba,tbb,tbc,tbd
      integer iret
! for ndsl
      real(kind=kind_evod),allocatable:: trie_ls_rqt(:,:,:)
     &,                                  trio_ls_rqt(:,:,:)
     &,                                  rqt_gr_a_1(:,:)
     &,                                  rqt_gr_a_2(:,:)
!
!     t3=rtc()
!jw   call mpi_barrier(mpi_comm_all,ierr)

      call mpi_barrier(mc_comp,ierr)

!     ioproc = nodes_comp - 1
      ioproc = 0

!     t4=rtc()
!     tba=t4-t3
!jw      if(nodes_comp .lt. 1 .or. nodes_comp .gt. nodes) then
!jw        print *, '  NODES_COMP UNDEFINED, CANNOT DO I.O '
!jw        call mpi_finalize()
!jw         stop 333
!jw      endif
!
!     if(.not. allocated ( trieo_ls_node)) then
!       allocate ( trieo_ls_node  ( len_trie_ls_max+len_trio_ls_max,
!    &                            2, 3*levs+1*levh+1 ) )
!     endif

!     t3=rtc()
!jw      call shapeset (ls_nodes,max_ls_nodes,pdryini)
!jw      call MPI_BARRIER(mpi_comm_all,ierr)
!jw      call MPI_BARRIER(mpi_comp,ierr)
!     t4=rtc()
!     tbb=t4-t3
       
!     if (.not. allocated (trieo_ls_nodes_buf) )then
!       allocate( trieo_ls_nodes_buf ( len_trie_ls_max+len_trio_ls_max,
!    x                               2, 3*levs+1*levh+1, nodes,1 ) )
!     endif

!     t1=rtc()

!!
      JDATE  = IDATE
      ndigyr = 4
      IF (NDIGYR == 2) THEN
        JDATE(4) = MOD(IDATE(4)-1,100) + 1
      ENDIF

      if ( iau .and. fhour >= 6.0 ) then
        idat = 0
        idat(1) = idate(4)
        idat(2) = idate(2)
        idat(3) = idate(3)
        idat(5) = idate(1)
        rinc = 0.
        rinc(2) = 6.
        call w3movdat(rinc,idat,jdat)
        jdate_iau(4) = jdat(1)
        jdate_iau(2) = jdat(2)
        jdate_iau(3) = jdat(3)
        jdate_iau(1) = jdat(5)
        fhour_iau = fhour - 6.0
        zhour_iau = zhour - 6.0
        if ( zhour_iau < 0.0 ) zhour_iau = 0.0
      else
        jdate_iau = jdate
        fhour_iau = fhour
        zhour_iau = zhour
      endif

!sela set lfnhr to false for writing one step output etc.
      lfnhr = .true.    ! no output
!     lfnhr=3600*abs(fhour-nint(fhour)).le.1.or.phour.eq.0
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
     &print *,' in wrtout_dynamics cfhour=',cfhour,' ens_nam=',
     &  ens_nam,'fhour=',fhour,'lfnhr=',lfnhr
!
      nosig = 61
!!
!     t3=rtc()
      call MPI_BARRIER(mpi_comm_all,ierr)
!     t4=rtc()
!
C*** BUILD STATE ON EACH NODE ********
c build state on each node.   COMP tasks only
c assemble upair state first then sfc state,
c then (only if liope)  flux state.
!
      if(nemsio_out) then
!       t3=rtc()

        if (.not. allocated(zsg)) allocate(zsg(lonf,lats_node_a))
        if (.not. allocated(psg)) allocate(psg(lonf,lats_node_a))
        if (.not. allocated(rqg)) allocate(rqg(lonf,lats_node_a,levh))
        if (.not. allocated(dpg)) allocate(dpg(lonf,lats_node_a,levs))
        if (.not. allocated(uug)) allocate(uug(lonf,lats_node_a,levs))
        if (.not. allocated(vvg)) allocate(vvg(lonf,lats_node_a,levs))
        if (.not. allocated(ttg)) allocate(ttg(lonf,lats_node_a,levs))

        if(mc_comp .ne. MPI_COMM_NULL) then

          do lan=1,lats_node_a
            jlonf = (lan-1)*lonf
            zsg(1:lonf,lan) = grid_gr(jlonf+1:jlonf+lonf,g_gz)
          enddo
          do k=1,levh
            do lan=1,lats_node_a
              jlonf = (lan-1)*lonf
              rqg(1:lonf,lan,k)=
     &        grid_gr(jlonf+1:jlonf+lonf,g_rq-1+k)
            enddo
          enddo

          do lan=1,lats_node_a
            lat      = global_lats_a(ipt_lats_node_a-1+lan)
            lons_lat = lonsperlat(lat)
            tx1      = one / coslat_a(lat)
            jlonf    = (lan-1)*lonf

            if (gen_coord_hybrid) then
              psg(1:lons_lat,lan) = grid_gr(jlonf+1:jlonf+lons_lat,g_q)
            else
              psg(1:lons_lat,lan) = 
     &        exp(grid_gr(jlonf+1:jlonf+lons_lat,g_q))
            endif

            if (gen_coord_hybrid) then        ! for general sigma-theta-p hybrid
              if(mass_dp) then
                do k=1,levs
                  do i=1,lons_lat
                    dpg(i,lan,k) = grid_gr(i+jlonf,g_dp-1+k)
                  enddo
                enddo
              else
                if (.not. allocated(tki)) allocate (tki(lonf,levs+1))
                tki(:,1)       = 0.0
                tki(:,levs+1)  = 0.0
                do k=2,levs
                  do i=1,lons_lat
                    tkrt0 = ( grid_gr(i+jlonf,g_tt-1+k-1)
     &                       +grid_gr(i+jlonf,g_tt-1+k) )
     &                        /(thref(k-1)+thref(k))
                    tki (i,k) = ck5(k)*tkrt0**rkappa
                  enddo
                enddo
                do k=1,levs
                  do i=1,lons_lat
                    dpg(i,lan,k) = ak5(k)-ak5(k+1)+(bk5(k)-bk5(k+1))
     &                       * psg(i,lan) + tki(i,k) - tki(i,k+1)
                  enddo
                enddo
              endif
            else if( hybrid ) then            ! for sigma-p hybrid (ECWMF)
              do k=1,levs
                kk = levs - k + 1
                do i=1,lons_lat
                  dpg(i,lan,k) = ak5(kk+1)-ak5(kk)
     &                     + (bk5(kk+1)-bk5(kk)) * psg(i,lan)
                enddo
              enddo
            else                              ! For sigma coordinate
              do k=1,levs
                do i=1,lons_lat
                  dpg(i,lan,k) = (si(k) - si(k+1)) * psg(i,lan)
                enddo
              enddo
            endif
            if (.not. allocated(tfac)) allocate (tfac(lonf,levs))
            if (thermodyn_id == 3) then
              if (.not. allocated(sumq)) allocate (sumq(lonf,levs))
              do k=1,levs
                do i=1,lons_lat
                  tfac(i,k) = 0.0
                  sumq(i,k) = 0.0
                enddo
              enddo
              do nn=1,ntrac
                nnl = (nn-1)*levs
                if (cpi(nn) .ne. 0.0) then
                  do k=1,levs
                    do i=1,lons_lat
                      sumq(i,k) = sumq(i,k) + rqg(i,lan,nnl+k)
                      tfac(i,k) = tfac(i,k) + cpi(nn)*rqg(i,lan,nnl+k)
                    enddo
                  enddo
                endif
              enddo
              do k=1,levs
                do i=1,lons_lat
                  tfac(i,k) = (one-sumq(i,k))*cpi(0) + tfac(i,k)
                enddo
              enddo
            else
              do k=1,levs
                do i=1,lons_lat
                  tfac(i,k) = one + fv*max(rqg(i,lan,k),qmin)
                enddo
              enddo
            endif
            do k=1,levs
              do i=1,lons_lat
                uug(i,lan,k) = grid_gr(i+jlonf,g_uu-1+k) * tx1
                vvg(i,lan,k) = grid_gr(i+jlonf,g_vv-1+k) * tx1
                ttg(i,lan,k) = grid_gr(i+jlonf,g_tt-1+k) / tfac(i,k)
              enddo
            enddo
            do k=1,levs
              do i=1,lons_lat
                dpg(i,lan,k) = cb2pa*dpg(i,lan,k)
              enddo
            enddo
            do i=1,lons_lat
              psg(i,lan) = cb2pa*psg(i,lan)
            enddo

          enddo

        endif                 ! comp node
!
!  done with state build
!  NOW STATE IS ASSEMBLED ON EACH NODE.  GET EVERYTHING OFF THE COMPUTE
!  NODES (currently done with a send to the I/O task_
!  send state to I/O task.  All tasks
!
        call grid_collect (zsg,psg,uug,vvg,ttg,rqg,dpg,
     &                         global_lats_a,lonsperlat)
!
      endif
!
! 
      if(sigio_out) then       !add sigio out
!
!*** for enthalpy and ps
! keep enthalpy and ps variables before write
!
!         if(run_enthalpy) then
!          do k=1,levs
!            kk = P_TE + k - 1
!            trie_te(:,:,k) = trie_ls(:,:,kk)
!            trio_te(:,:,k) = trio_ls(:,:,kk)
!          enddo
!          trie_q (:,:) = trie_ls(:,:,P_Q)
!          trio_q (:,:) = trio_ls(:,:,P_Q)
!
!          direction=-1          ! from (enthalpy,ps) to (virttemp,lnps)
!          call spect_tv_enthalpy_ps
!!!   &       (direction,run_enthalpy,
!     &       (direction,
!     X        TRIE_LS(1,1,P_Q ), TRIO_LS(1,1,P_Q ),
!     X        TRIE_LS(1,1,P_TE), TRIO_LS(1,1,P_TE),
!     X        TRIE_LS(1,1,P_RQ), TRIO_LS(1,1,P_RQ),
!     &        ls_node,ls_nodes,max_ls_nodes,
!     &        lats_nodes_r,global_lats_r,lonsperlar,
!     &        plnev_r,plnod_r,plnew_r,plnow_r)
!
!        endif           ! (run enthalpy
!!
!       thermodyn_id_out = 1
!       if( gen_coord_hybrid ) then
!         sfcpress_id_out  = 2
!       else
!         sfcpress_id_out  = 1
!       endif

        thermodyn_id_out = thermodyn_id
        sfcpress_id_out  = sfcpress_id

!
! n time step spectral file
!
!        filename = sigf//trim(CFHOUR)
         filename = sigf//CFHOUR
         if (me == ioproc) print *,'in twrites_hst,filename=',
     &                            trim(filename)


         if( .not. ndslfv ) then

           if(iau) then
             call twrites_hst(filename,ioproc,fhour_iau,jdate_iau,
     &              ls_nodes,max_ls_nodes,trie_ls,trio_ls,
     &              trie_ls,trio_ls,
     &              thermodyn_id_out,sfcpress_id_out,pdryini)
           else
             call twrites_hst(filename,ioproc,fhour,idate,
     &              ls_nodes,max_ls_nodes,trie_ls,trio_ls,
     &              trie_ls,trio_ls,
     &              thermodyn_id_out,sfcpress_id_out,pdryini)
           end if

         else

           if ( .not. allocated (trie_ls_rqt) )then
             allocate( trie_ls_rqt(len_trie_ls,2,levh) )
           endif
           if ( .not. allocated (trio_ls_rqt) )then
             allocate( trio_ls_rqt(len_trio_ls,2,levh) )
           endif
           if ( .not. allocated (rqt_gr_a_1) )then
             allocate( rqt_gr_a_1(lonfx*levh,lats_dim_ext) )
           endif
           if ( .not. allocated (rqt_gr_a_2) )then
             allocate( rqt_gr_a_2(lonfx*levh,lats_dim_ext) )
           endif
           call do_dynamics_gridc2rqt(grid_gr,rqt_gr_a_2,
     &                               global_lats_a,lonsperlat)
           call grid_to_spect_rqt(rqt_gr_a_1,rqt_gr_a_2,
     &                         trie_ls_rqt,trio_ls_rqt,levh,
     &                         ls_node,ls_nodes,max_ls_nodes,
     &                 lats_nodes_a,global_lats_a,lonsperlat,
     &                   epse,epso,plnew_a,plnow_a)
           call twrites_hst(filename,ioproc,fhour,idate,
     &            ls_nodes,max_ls_nodes,trie_ls,trio_ls,
     &            trie_ls_rqt,trio_ls_rqt,
     &            thermodyn_id_out,sfcpress_id_out,pdryini)

         endif

         if (me == ioproc) print *,'finished end of sigio output for ',
     &     trim(filename)
!
!        if (runenthalpy) then
!! te
!          do k=1,levs
!            kk = P_TE + k - 1
!            trie_ls(:,:,kk) = trie_te(:,:,k)
!            trio_ls(:,:,kk) = trio_te(:,:,k)
!          enddo
!! ps
!          trie_ls(:,:,P_Q) = trie_q (:,:)
!          trio_ls(:,:,P_Q) = trio_q (:,:)
!!
!        endif      
!
      endif                                        !end sigio_out

!
      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE wrt_restart_dynamics(TRIE_LS,TRIO_LS,grid_gr,
     &        SI,fhour,idate,igen,pdryini,
     x        ls_node,ls_nodes,max_ls_nodes,
     &        global_lats_a,lonsperlat,lats_nodes_a,
     &        epse,epso,plnew_a,plnow_a,
     &        ens_nam,kdt,nfcstdate7)
!
      use gfs_dyn_resol_def
      use gfs_dyn_layout1
      use gfs_dyn_mpi_def
      use namelist_dynamics_def, only : ndslfv, lsidea
      use do_dynamics_mod
!
      implicit none
!
      real(kind=kind_evod) fhour, pdryini
      character (len=*)  :: ens_nam
      character (255)    :: filename
!
      integer              idate(4), igen
      INTEGER              LS_NODE (LS_DIM*3)
      integer              ls_nodes(ls_dim,nodes)
      integer              max_ls_nodes(nodes)
      integer step,kdt,nfcstdate7(7)
 
      real(kind=kind_evod) si(levp1)
!c
      REAL(KIND=KIND_EVOD) TRIE_LS(LEN_TRIE_LS,2,lotls)
      REAL(KIND=KIND_EVOD) TRIO_LS(LEN_TRIO_LS,2,lotls)
      REAL(KIND=KIND_grid) grid_gr(lonf,lats_node_a_max,lotgr)

      real(kind=kind_evod)  epse  (len_trie_ls)
      real(kind=kind_evod)  epso  (len_trio_ls)
      real(kind=kind_evod)   plnew_a(len_trie_ls,latg2)
      real(kind=kind_evod)   plnow_a(len_trio_ls,latg2)

! for ndsl
      real(kind=kind_evod),allocatable:: trie_ls_rqt(:,:,:)
     &,                                  trio_ls_rqt(:,:,:)
     &,                                  rqt_gr_a_1(:,:)
     &,                                  rqt_gr_a_2(:,:)
!!
      integer, dimension(latg) :: global_lats_a, lonsperlat
      INTEGER                     lats_nodes_a(nodes)
!
!-- local variables
      integer IOPROC, IPRINT, needoro, iret, nfill
!!
      IPRINT = 0
      IOPROC = nodes - 1
      if( ndslfv ) then
        if ( .not. allocated (trie_ls_rqt) )then
           allocate( trie_ls_rqt(len_trie_ls,2,levh) )
        endif
        if ( .not. allocated (trio_ls_rqt) )then
           allocate( trio_ls_rqt(len_trio_ls,2,levh) )
        endif
         if ( .not. allocated (rqt_gr_a_1) )then
           allocate( rqt_gr_a_1(lonfx*levh,lats_dim_ext) )
        endif
         if ( .not. allocated (rqt_gr_a_2) )then
           allocate( rqt_gr_a_2(lonfx*levh,lats_dim_ext) )
        endif
      endif
!
      if (me == ioproc) print *,'in restart,lonsperlat=',lonsperlat
! n-1 time step spectral file
!
        step = -1
        filename = 'SIGR1'

        if( .not. ndslfv ) then
            CALL TWRITES_rst(filename,ioproc,FHOUR,idate,
     &                  SI,LS_NODES,MAX_LS_NODES,step,
     &                  trie_ls,trio_ls,trie_ls,trio_ls)

        else

            call do_dynamics_gridm2rqt(grid_gr,rqt_gr_a_2,
     &                                 global_lats_a,lonsperlat)
            call grid_to_spect_rqt(rqt_gr_a_1,rqt_gr_a_2,
     &                           trie_ls_rqt,trio_ls_rqt,levh,
     &                           ls_node,ls_nodes,max_ls_nodes,
     &                   lats_nodes_a,global_lats_a,lonsperlat,
     &                     epse,epso,plnew_a,plnow_a)
            CALL TWRITES_rst(filename,ioproc,FHOUR,idate,
     &                  SI,LS_NODES,MAX_LS_NODES,step,
     &                  trie_ls,trio_ls,trie_ls_rqt,trio_ls_rqt)

        endif

        if (me == ioproc) print *,'1 end of twritero_rst,',
     &                             trim(filename)
!
! n time step spectral file
!
        step = 0
        filename = 'SIGR2'
      
        if( .not. ndslfv ) then
            CALL TWRITES_rst(filename,ioproc,FHOUR,idate,
     &                  SI,LS_NODES,MAX_LS_NODES,step,
     &                  trie_ls,trio_ls,trie_ls,trio_ls)

        else

            call do_dynamics_gridc2rqt(grid_gr,rqt_gr_a_2,
     &                                 global_lats_a,lonsperlat)
            call grid_to_spect_rqt(rqt_gr_a_1,rqt_gr_a_2,
     &                           trie_ls_rqt,trio_ls_rqt,levh,
     &                           ls_node,ls_nodes,max_ls_nodes,
     &                   lats_nodes_a,global_lats_a,lonsperlat,
     &                     epse,epso,plnew_a,plnow_a)
            CALL TWRITES_rst(filename,ioproc,FHOUR,idate,
     &                  SI,LS_NODES,MAX_LS_NODES,step,
     &                  trie_ls,trio_ls,trie_ls_rqt,trio_ls_rqt)

        endif

        if (me == ioproc) print *,'2 end of twritero_rst for ',
     &                             trim(filename)
      IF(lsidea) THEN

        CALL TWRITES_rst_idea('fort.1051',ioproc,FHOUR,idate,
     &              SI,LS_NODES,MAX_LS_NODES,trie_ls,trio_ls)
      END IF

! n-1 time step grid file
!
       filename = 'GRDR1'

       CALL TWRITEG_rst(filename,ioproc,FHOUR,idate,
     X                SI,pdryini,global_lats_a,lonsperlat,lats_nodes_a,
     &                grid_gr(1,1,g_qm),
     &                grid_gr(1,1,g_dpm),grid_gr(1,1,g_ttm),
     &                grid_gr(1,1,g_uum),grid_gr(1,1,g_vvm),
     &                grid_gr(1,1,g_rm),grid_gr(1,1,g_gz),
     &    kdt,nfcstdate7 )
       if (me == ioproc) print *,'1 end twriteg_rst,',trim(filename)
!
! n time step grid file
!

      filename = 'GRDR2'
      CALL TWRITEG_rst(filename,ioproc,FHOUR,idate,
     X                SI,pdryini,global_lats_a,lonsperlat,lats_nodes_a,
     &                grid_gr(1,1,g_q),
     &                grid_gr(1,1,g_dp),grid_gr(1,1,g_tt),
     &                grid_gr(1,1,g_uu),grid_gr(1,1,g_vv),
     &                grid_gr(1,1,g_rq),grid_gr(1,1,g_gz),
     &    kdt,nfcstdate7 )
      if (me == ioproc) print *,'2 end twriteg_rst,',trim(filename)
      call mpi_barrier(mpi_comm_all,iret)
!
      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE wrtlog_dynamics(phour,fhour,idate)
      use gfs_dyn_resol_def
      use namelist_dynamics_def
      implicit none

      integer idate(4),ndigyr,nolog
      integer ks,kh,km,ndig,nfill
      character CFHOUR*40,CFORM*40
      logical lfnhr
      real phour,fhour
c
c     CREATE CFHOUR

csela set lfnhr to false for writing one step output etc.
      lfnhr=.true.    ! no output
ccmr  lfnhr=.false.   !    output
!      lfnhr=3600*abs(fhour-nint(fhour)).le.1.or.phour.eq.0
      lfnhr=3600*abs(fhour-nint(fhour)).le.1
      IF(LFNHR) THEN
        KH=NINT(FHOUR)
        NDIG=MAX(LOG10(KH+0.5)+1.,2.)
        WRITE(CFORM,'("(I",I1,".",I1,")")') NDIG,NDIG
        WRITE(CFHOUR,CFORM) KH
        WRITE(CFORM,'("(I",I1,".",I1,")")') NDIG,NDIG
        WRITE(CFHOUR,CFORM) KH
      ELSE
        KS=NINT(FHOUR*3600)
        KH=KS/3600
        KM=(KS-KH*3600)/60
        KS=KS-KH*3600-KM*60
        NDIG=MAX(LOG10(KH+0.5)+1.,2.)
        WRITE(CFORM,
     &      '("(I",I1,".",I1,",A1,I2.2,A1,I2.2)")') NDIG,NDIG
        WRITE(CFHOUR,CFORM) KH,':',KM,':',KS
      ENDIF
      if( nfill(ens_nam) == 0 ) then
      CFHOUR = CFHOUR(1:nfill(CFHOUR))
      else
      CFHOUR = CFHOUR(1:nfill(CFHOUR)) // ens_nam(1:nfill(ens_nam))
      endif
!      print *,' in wrtlog_dynamics cfhour=',cfhour,' ens_nam=',ens_nam

      nolog=99
      OPEN(NOlog,FILE='LOG.F'//CFHOUR,FORM='FORMATTED')
      write(nolog,100)fhour,idate
100   format(' completed mrf fhour=',f10.3,2x,4(i4,2x))
      CLOSE(NOlog)

      RETURN
      END


      subroutine  shapeset (ls_nodes,max_ls_nodes,pdryini)
!
      use gfs_dyn_resol_def
      use gfs_dyn_layout1
      use namelist_dynamics_def
      use gfs_dyn_mpi_def
      implicit none
!
      integer              ls_nodes(ls_dim,nodes), max_ls_nodes(nodes)
      integer              ierr,j,k,l,lenrec,locl,n,node
     &,                    indjoff, indev, indod
!
      real(kind=kind_evod) gencode,order,ppid,realform
      real(kind=kind_evod) subcen,tracers,trun,vcid,vmid,vtid
!
      real(kind=kind_evod) dummy(201-levp1-levs)
      real(kind=kind_evod) ensemble(2),dummy2(18)
!
      real(kind=kind_io4)   tmps(4+nodes+jcap1*nodes)
      real(kind=kind_io4)   tmpr(3+nodes+jcap1*(nodes-1))
      REAL (KIND=KIND_grid) pdryini
!
      INTEGER              GLOBAL_lats_a(latg)
      INTEGER                 lonsperlat(latg)
!
      integer  il,ilen,i,msgtag,ls_diml,nodesl,ioproc, itmpr
                                                                                                        
!  Now define shape of the coefficients array
!  as a function of node. This will define how
!  to assemble the few wavenumbers on each node
!  into a full coefficient array.
!
       IOPROC = nodes
       IF (LIOPE) then
!199    format(' GWVX MAX_LS_NODES ',i20)
        if (me == 0.or. me == ioproc) then
        tmps    = 0.
        tmps(1) = PDRYINI
        tmps(2:nodes_comp+1) = max_ls_nodes(1:nodes_comp)
        tmps(nodes_comp+2)   = ls_dim
        tmps(nodes_comp+3)   = len_trie_ls_max
        tmps(nodes_comp+4)   = len_trio_ls_max
        il = nodes_comp+4

        do i=1,nodes_comp
          do j=1,ls_dim
             il = il+1
             tmps(il) = ls_nodes(j,i)
          enddo
        enddo
        ilen   = 4 + nodes_comp + jcap1*nodes_comp
        msgtag = 2345
        if(me == 0) then
            CALL mpi_send(tmps,ilen,MPI_R_IO,ioproc,
     &                msgtag,MPI_COMM_ALL,info)
           endif
        endif
!
        if (me == ioproc) then
         ilen   = 4+nodes_comp+jcap1*(nodes_comp)
         msgtag = 2345
             CALL mpi_recv(tmpr,ilen,MPI_R_IO,0,
     &                     msgtag,MPI_COMM_ALL,stat,info)

          itmpr = 3 + nodes + jcap1*(nodes-1)
          tmps(1:itmpr) = tmpr(1:itmpr)
          ls_nodes = 0
          pdryini  = tmps(1)
          max_ls_nodes(1:nodes_comp) = int(tmps(2:nodes_comp+1))
          ls_diml = int(tmps(nodes_comp+2))
          len_trie_ls_max = int(tmps(nodes_comp+3))
          len_trio_ls_max = int(tmps(nodes_comp+4))
          il = nodes_comp + 3 + 1
                                                                                                        
          do i=1,nodes_comp
            do j=1,ls_diml
              il=il+1
              ls_nodes(j,i)=int(tmps(il))
            enddo
          enddo
        endif
      ENDIF

      return
      end
