      subroutine get_w_z(grid_gr,
     &                 trie_ls,trio_ls,
     &                 LS_NODE,LS_NODES,MAX_LS_NODES,
     &                 LATS_NODES_A,GLOBAL_LATS_A,LONSPERLAT,
     &                 EPSE,EPSO,EPSEDN,EPSODN,
     &                 PLNEV_A,PLNOD_A,PLNEW_A,PLNOW_A,
     &                 PDDEV_A,PDDOD_A,SNNP1EV,SNNP1OD,
     &                 kdt, deltim, restart_step)
!!
! Program History Log:
! Mar 2015    Henry Juang	use existed variables to get w hydrostatically
! Feb 2016    Weiyu Yang        add u, v, t variables for the WAM-IPE coupling.
!                               The outputs for coupling are: uug --> u,
!                                                             vvg --> v,
!                                                             ttg --> t,
!                                                             wwg --> w,
!                                                             zzg --> z.
!                                                    rqg --> q, oz, clw, O1, O2
!                                                             n2g --> N2.
! Mar 2016    Sajal Kar         zzg & wwg for g(z) when var_g=true
! May 2016    Weiyu Yang        Convert the O, O2 and N2 into total number of density.
!----------------------------------------------
!!
      use gfs_dyn_resol_def
      use gfs_dyn_layout1
      use gfs_dyn_gg_def
      use gfs_dyn_vert_def
      use gfs_dyn_coordinate_def
      use gfs_dyn_date_def
      use namelist_dynamics_def
      use gfs_dyn_mpi_def
      use gfs_dyn_dfi_mod
      use gfs_dyn_physcons, only: p0 => con_p0, fv => con_fvirt
     &,                           re => con_rerth, g0 => con_g
     &,                           cpd => con_cp, rkap => con_rocp
     &,                           con_amw, con_rgas, con_amo2, con_avgd
     &,                           con_amo3, con_boltz
      use do_dynamics_mod
      use gfs_dyn_tracer_const, only: cpi
      use gci,                  only: grid_collect_ipe
      use get_variables_for_WAM_IPE_coupling
!!     
      IMPLICIT NONE

      INTEGER :: LONSPERLAT(LATG)
!!     
      REAL(KIND=KIND_EVOD) TRIE_LS(LEN_TRIE_LS,2,LOTls)
      REAL(KIND=KIND_EVOD) TRIO_LS(LEN_TRIO_LS,2,LOTls)
      REAL(KIND=KIND_GRID) GRID_GR(lonf*lats_node_a_max,lotgr)

      integer          ls_node(ls_dim,3), kdt
!
      INTEGER          LS_NODES(LS_DIM,NODES)
      INTEGER          MAX_LS_NODES   (NODES)
      INTEGER          LATS_NODES_A   (NODES)
      INTEGER          GLOBAL_LATS_A(LATG)
!
      REAL(KIND=KIND_EVOD)      EPSE(LEN_TRIE_LS)
      REAL(KIND=KIND_EVOD)      EPSO(LEN_TRIO_LS)
      REAL(KIND=KIND_EVOD)    EPSEDN(LEN_TRIE_LS)
      REAL(KIND=KIND_EVOD)    EPSODN(LEN_TRIO_LS)
      REAL(KIND=KIND_EVOD)   SNNP1EV(LEN_TRIE_LS)
      REAL(KIND=KIND_EVOD)   SNNP1OD(LEN_TRIO_LS)

      REAL(KIND=KIND_EVOD)   PLNEV_A(LEN_TRIE_LS,LATG2)
      REAL(KIND=KIND_EVOD)   PLNOD_A(LEN_TRIO_LS,LATG2)
      REAL(KIND=KIND_EVOD)   PLNEW_A(LEN_TRIE_LS,LATG2)
      REAL(KIND=KIND_EVOD)   PLNOW_A(LEN_TRIO_LS,LATG2)
      REAL(KIND=KIND_EVOD)   PDDEV_A(LEN_TRIE_LS,LATG2)
      REAL(KIND=KIND_EVOD)   PDDOD_A(LEN_TRIO_LS,LATG2)

      REAL(KIND=KIND_EVOD) SYN_GR_A_1(LONFX*LOTS,LATS_DIM_EXT)
      REAL(KIND=KIND_EVOD) DYN_GR_A_1(LONFX*LOTD,LATS_DIM_EXT)
      REAL(KIND=KIND_EVOD) ANL_GR_A_1(LONFX*LOTA,LATS_DIM_EXT)

      REAL(KIND=KIND_EVOD) SYN_GR_A_2(LONFX*LOTS,LATS_DIM_EXT)
      REAL(KIND=KIND_EVOD) DYN_GR_A_2(LONFX*LOTD,LATS_DIM_EXT)
      REAL(KIND=KIND_EVOD) ANL_GR_A_2(LONFX*LOTA,LATS_DIM_EXT)

      REAL(KIND=KIND_EVOD) SYN_GR_S_Z(LONFX     ,LATS_DIM_EXT)
      REAL(KIND=KIND_EVOD) ANL_GR_A_W(LONFX*levs,LATS_DIM_EXT)
      REAL(KIND=KIND_EVOD) ANL_GR_A_Z(LONFX*levs,LATS_DIM_EXT)
      REAL(KIND=KIND_EVOD) TEREF(LEVP1),CK5P(LEVP1)

      REAL(KIND=KIND_GRID) tfac(lonf, levs), sumq(lonf, levs)
      REAL(KIND=KIND_GRID) tx1
      REAL(KIND=KIND_GRID), PARAMETER :: qmin=1.e-10
!!     
      REAL (KIND=KIND_grid), parameter :: cons0=0.0d0, cons1=1.0d0,
     &                                    cons2=2.0d0, cons0p5 = 0.5d0,
     &                                    con_amn2=28.013
!sk
      REAL(KIND=KIND_GRID) zs,phi,grav, rkapi, rmdo1, rmdo2, rmdn2
      REAL(KIND=KIND_EVOD) pptt, ppttbz, mmm, con_amo1, avgdr, mmmavgd

      REAL(KIND=KIND_GRID) phis(lonf,lats_node_a)
      real, parameter:: g0re = g0*re, g0re2 = g0*re*re
      real           :: deltim
      logical :: restart_step

      INTEGER               I,J,K,L,LOCL,N

      integer               lan,lat,nt, levs3, levs4
      integer               lon_dim,lons_lat,node
      integer               njeff,lon,jlonf,ilan
     &,                     ngptcd, nn, nnl

      IF(.NOT. ALLOCATED(wwg)) THEN 
         ALLOCATE(wwg(lonf,lats_node_a,levs)) 
         wwg = 0.0 
      endif 
      IF(.NOT. ALLOCATED(zzg)) then  
         ALLOCATE(zzg(lonf,lats_node_a,levs)) 
         zzg = 0.0 
      endif 
      IF(.NOT. ALLOCATED(uug)) then  
         ALLOCATE(uug(lonf,lats_node_a,levs)) 
         uug = 0.0 
      endif 
      IF(.NOT. ALLOCATED(vvg)) then 
         ALLOCATE(vvg(lonf,lats_node_a,levs)) 
         vvg=0.0 
      endif 
      IF(.NOT. ALLOCATED(ttg)) then 
         ALLOCATE(ttg(lonf,lats_node_a,levs)) 
         ttg = 0.0 
      endif 
      IF(.NOT. ALLOCATED(rqg)) then  
         ALLOCATE(rqg(lonf,lats_node_a,levh)) 
         rqg = 0.0 
      endif 
      IF(.NOT. ALLOCATED(n2g)) then 
         ALLOCATE(n2g(lonf,lats_node_a,levs)) 
         n2g = 0.0 
      endif
      IF(.NOT. ALLOCATED(den)) then
         ALLOCATE(den(lonf,lats_node_a,levs))
         den = 0.0
      endif
      IF(.NOT. ALLOCATED(gmol)) then
         ALLOCATE(gmol(lonf,lats_node_a,levs))
         gmol = 0.0
      endif
      IF(.NOT. ALLOCATED(ppg)) then
         ALLOCATE(ppg(lonf,lats_node_a,levs))
         ppg = 0.0
      endif
      IF(.NOT. ALLOCATED(ps)) then
         ALLOCATE(ps(lonf,lats_node_a))
         ps = 0.0
      endif
!----------------------------------------------------------
      if (me < num_pes_fcst) then
        levs3 = levs * 3
        levs4 = levs * 4
        ngptcd = ngptc
        rkapi = 1.0 / rkap
        if( thermodyn_id.eq.3 ) then
          do k=1,levp1
            thref(k)=300.*cpd
            teref(k)=255.*cpd
          enddo
        else
          do k=1,levp1
            thref(k)=300.
            teref(k)=255.
          enddo
       endif
       ck5p(levp1)=ck5(levp1)
       do k=1,levp1
         ck5p(k)=ck5(k)*(teref(k)/thref(k))**rkapi
       enddo
!----------------------------------------------------------
! transform spectral to grid to syn
          call spect_to_grid(trie_ls,trio_ls, 
     &                       syn_gr_a_1,syn_gr_a_2,
     &                       ls_node,ls_nodes,max_ls_nodes,
     &                       lats_nodes_a,global_lats_a,lonsperlat,
     &                       epse,epso,epsedn,epsodn,
     &                       snnp1ev,snnp1od,plnev_a,plnod_a)

! ---------------------------------------------------------------------
! transform spectral deriative to dyn
        call spect_to_gridxy(trie_ls,trio_ls,
     &                       syn_gr_a_1,syn_gr_a_2,
     &                       dyn_gr_a_1,dyn_gr_a_2,
     &                       ls_node,ls_nodes,max_ls_nodes,
     &                       lats_nodes_a,global_lats_a,lonsperlat,
     &                       pddev_a,pddod_a)
! ------------------------------------------------------------------
! move zs in grid_gr to syn
        do lan=1,lats_node_a
          lat = global_lats_a(ipt_lats_node_a-1+lan)
          lons_lat = lonsperlat(lat)

          tx1      = 1.0 / coslat_a(lat)
          jlonf = (lan-1)*lonf

          do lon=1,lons_lat,ngptc
            njeff = min(ngptc,lons_lat-lon+1)
            if (height_dependent_g) then
!sk compute surface geopotential (phis) as in subroutine phi2z
              do i=lon,lon+njeff-1
                ilan=i+jlonf
                syn_gr_s_z(i,lan)=grid_gr(ilan,g_gz)
                zs=syn_gr_s_z(i,lan)
! Since later when calculate zzg, the anl_gr_a_z is phi/g+g0*zs, to
! avoid add zs twice, need to pre-remove the term g0*zs.  Weiyu.
!-------------------------------------------------------------------
!                phis(i,lan)=g0re*zs/(re+zs)
                phis(i,lan)=g0*zs*(re / (re+zs)-1)
              enddo
            else !g(z)=g0
              do i=lon,lon+njeff-1
                ilan=i+jlonf
                syn_gr_s_z(i,lan)=grid_gr(ilan,g_gz)
              enddo
            endif
          enddo

          do i=1,lons_lat
            ps(i,lan) = grid_gr(i+jlonf,g_q)
          enddo

          do i=1,lons_lat
            pptt = ak5(1) + bk5(1) * ps(i, lan) + ck5p(1)
            do k=1,levs
              ppg(i, lan, k) = pptt
              pptt         = ak5(k+1)+bk5(k+1)*ps(i, lan)+ck5p(k+1)
              ppg(i, lan, k) = (ppg(i, lan, k) + pptt) * 500.0  ! ~Pa.
            end do
          end do
 
          do k=1,levh
            do i=1,lons_lat
              rqg(i,lan,k) = max(grid_gr(i+jlonf,g_rq-1+k), qmin)
            enddo
          enddo

          do k=1,levs
            do i=1,lons_lat
              n2g(i, lan, k) = 1.0
            enddo
          enddo

! rqg containing:
! ntrac=5, nn=1, ==> q.
!          nn=2, ==> oz.
!          nn=3, ==> clw.
!          nn=4, ==> o1.
!          nn=5, ==> o2.
!------------------------
          do nn=1,ntrac
            if (cpi(nn) .ne. 0.0) then
              nnl = (nn-1)*levs
              do k=1,levs
                do i=1,lons_lat
                  n2g(i, lan, k) = n2g(i, lan, k) - rqg(i,lan,nnl+k)
                enddo
              enddo
            endif
          enddo

         if (thermodyn_id == 3) then
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
               tfac(i,k) = (1.0-sumq(i,k))*cpi(0) + tfac(i,k)
             enddo
           enddo
         else
           do k=1,levs
             do i=1,lons_lat
               tfac(i,k) = 1.0 + fv*max(rqg(i,lan,k),qmin)
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

! Convert O, O2 and N2 to the total number of density.
!-----------------------------------------------------
         con_amo1 = con_amo2 * 0.5        ! Mi for O1.
         avgdr    = con_avgd / con_rgas   ! average number / R.

         do k=1,levs
           do i=1,lons_lat

! Calculate the mean molecular mass.
!-----------------------------------
             mmm = rqg(i, lan, k) / con_amw                ! for q weight.
             mmm = mmm + rqg(i, lan, k+levs ) / con_amo3   ! add ozone weight.
             mmm = mmm + rqg(i, lan, k+levs3) / con_amo1   ! add     O weight.
             mmm = mmm + rqg(i, lan, k+levs4) / con_amo2   ! add    O2 weight.
             mmm = mmm + n2g(i, lan, k)       / con_amn2     ! add    N2 weight.
             mmm = 1.0 / mmm                               ! final mmm.
!             IF(i == 4. and. lan == 1 .and.me == 0) THEN
!               print*, 'in get_w_z, i, lan, k, mmm=',i, lan, k, mmm
!             END IF

             rmdo1 = mmm * avgdr         ! Md * average number / R.
             rmdo2 = rmdo1 / con_amo2    ! Md * average number / R / Mi_o2.
             rmdn2 = rmdo1 / con_amn2    ! Md * average number / R / Mi_n2.
             rmdo1 = rmdo1 / con_amo1    ! Md * average number / R / Mi_o1.
             
             pptt   = ppg(i, lan, k) / ttg(i, lan, k)  ! P / T.
             ppttbz = pptt / con_boltz
             rqg(i, lan, k+levs3) = rqg(i, lan, k+levs3) * rmdo1 * pptt
             rqg(i, lan, k+levs4) = rqg(i, lan, k+levs4) * rmdo2 * pptt
             n2g(i, lan, k)       = n2g(i, lan, k)       * rmdn2 * pptt

             den (i, lan, k) = mmm * 1.0e-3 / con_avgd * ppttbz ! unit ~ kg/m^3.
             gmol(i, lan, k) = mmm ! unit ~ g/mol (amu).
           end do
         end do
       enddo
!
! =============================
! -------------------------------------------------------------------
       do lan=1,lats_node_a  
!
           lat      = global_lats_a(ipt_lats_node_a-1+lan)
           lon_dim  = lon_dims_a(lan)
           lons_lat = lonsperlat(lat)
! --------------------------------------------------------------------
!$omp parallel do schedule(dynamic,1) private(lon)
!$omp+private(njeff)
           do lon=1,lons_lat,ngptcd
!!
              njeff = min(ngptcd,lons_lat-lon+1)
!
              if( gen_coord_hybrid ) then                        
                 if( thermodyn_id == 3 ) then                    
                   call gfidi_hyb_gc_h_w_z(lon_dim, njeff, lat,
     &               syn_gr_a_2(lon+(ksd   -1)*lon_dim,lan),
     &               syn_gr_a_2(lon+(kst   -1)*lon_dim,lan),
     &               syn_gr_a_2(lon+(ksu   -1)*lon_dim,lan),
     &               syn_gr_a_2(lon+(ksv   -1)*lon_dim,lan),
     &               syn_gr_a_2(lon+(ksr   -1)*lon_dim,lan),
     &               syn_gr_a_2(lon+(kspphi-1)*lon_dim,lan),
     &               syn_gr_a_2(lon+(ksplam-1)*lon_dim,lan),
     &               syn_gr_a_2(lon+(ksq   -1)*lon_dim,lan),
     &               syn_gr_a_2(lon+(kzsphi-1)*lon_dim,lan),
     &               syn_gr_a_2(lon+(kzslam-1)*lon_dim,lan),
     &               syn_gr_s_z(lon                   ,lan),
     &               rcs2_a(min(lat,latg-lat+1)),           
     &               dyn_gr_a_2(lon+(kdtphi-1)*lon_dim,lan),
     &               dyn_gr_a_2(lon+(kdtlam-1)*lon_dim,lan),
     &               dyn_gr_a_2(lon+(kdrphi-1)*lon_dim,lan),
     &               dyn_gr_a_2(lon+(kdrlam-1)*lon_dim,lan),
     &               anl_gr_a_w(lon,lan),
     &               anl_gr_a_z(lon,lan),me)

                 else
                   print *,' get_w_z error: not enthalpy '
                 endif
              else
                 print *,' get_w_z error: not gen_coord_hybrid '
              endif 
!
           enddo   !lon
! ---------------------------------------------------------------
       enddo   ! end of lan
!
! ===================================================
! move w in anl to grid_gr
       do lan=1,lats_node_a
         lat = global_lats_a(ipt_lats_node_a-1+lan)
         lon_dim = lon_dims_a(lan)
         lons_lat = lonsperlat(lat)
!$omp parallel do schedule(dynamic,1) private(lon)
!$omp+private(i,k,ilan,njeff)
         do lon=1,lons_lat,ngptc
           njeff = min(ngptc,lons_lat-lon+1)
           if (height_dependent_g) then
!sk compute zzg and grav as in subroutine phi2z, then compute wwg
             do k=1,levs
               do i=lon,lon+njeff-1
                 phi=g0*anl_gr_a_z(i+(k-1)*lon_dim,lan)
                 zzg(i,lan,k)=re*(phis(i,lan)+phi)
     &                        /(g0re-(phis(i,lan)+phi))
                 grav=g0re2/((re+zzg(i,lan,k))*(re+zzg(i,lan,k)))
                 wwg(i,lan,k)=(g0/grav)*anl_gr_a_w(i+(k-1)*lon_dim,lan)
               enddo
             enddo
           else   !g(z)=g0
             do k=1,levs
               do i=lon,lon+njeff-1
                 wwg(i,lan,k)=anl_gr_a_w(i+(k-1)*lon_dim,lan)
                 zzg(i,lan,k)=anl_gr_a_z(i+(k-1)*lon_dim,lan)
               enddo
             enddo
           endif
         enddo
       enddo

! -------------------------------------------------------------------
      endif ! only for fcst nodes

! Check outputs.
!---------------
!      print*, 'In get_w_z, uug = ', uug(2, 2, 4), uug(4, 1,  19)
!      print*, 'In get_w_z, vvg = ', vvg(2, 2, 4), vvg(4, 1,  19)
!      print*, 'In get_w_z, ttg = ', ttg(2, 2, 4), ttg(4, 1,  19)
!      print*, 'In get_w_z, wwg = ', wwg(2, 2, 4), wwg(4, 1,  19)
! At level 1.
!------------
!      print*, 'In get_w_z, zzg = ', zzg(2, 2, 1), zzg(4, 1,  1)
!      print*, 'In get_w_z, ppg = ', ppg(2, 2, 1), ppg(4, 1,  1)
! At level 99 and 149.
!---------------------
!      print*, 'at level=99'
!      print*, 'In get_w_z, zzg = ', zzg(2, 2, 99), zzg(4, 1,  99)
!      print*, 'In get_w_z, ppg = ', ppg(2, 2, 99), ppg(4, 1,  99)
!      print*, 'In get_w_z, o1g = ', rqg(2, 2, 549), rqg(4, 1,  549)
!      print*, 'In get_w_z, o2g = ', rqg(2, 2, 699), rqg(4, 1,  699)
!      print*, 'In get_w_z, n2g = ', n2g(2, 2, 99), n2g(4, 1,  99)
!      print*, 'at level=149'
!      print*, 'In get_w_z, zzg = ', zzg(2, 2, 149), zzg(4, 1,  149)
!      print*, 'In get_w_z, ppg = ', ppg(2, 2, 149), ppg(4, 1,  149)
!      print*, 'In get_w_z, o1g = ', rqg(2, 2, 599), rqg(4, 1,  599)
!      print*, 'In get_w_z, o2g = ', rqg(2, 2, 749), rqg(4, 1,  749)
!      print*, 'In get_w_z, n2g = ', n2g(2, 2, 149), n2g(4, 1,  149)
!
! The following is only to output for making figures and comparisons. WY.
!------------------------------------------------------------------------

! The following is to output the interface restart file for WAM-IPE coupling
! restart run.
!---------------------------------------------------------------------------
! Joe Schoonover ( Jan. 31, 2018 ) : Needed to change kdt to
! kdt+1 to write the neutrals at a multiple of FHRES, instead
! of a FHRES+dt
      IF(wam_ipe_cpl_rst_output .AND. kdt /= 0 .AND.
     &  MOD(NINT(deltim) * (kdt+1), NINT(FHRES) * 3600) == 0) THEN
        PRINT*,'Write out the WAM-IPE rst file needed for IPE, kdt=',
     &       kdt
        CALL grid_collect_ipe(wwg,zzg,uug,vvg,
     &                        ttg,rqg,n2g,global_lats_a,lonsperlat,
     &                        lats_nodes_a,kdt,deltim,restart_step)
      END IF

! The following is to the NetCDF diagnostic files.
!-------------------------------------------------
      IF(NC_output .AND.
     &  MOD(NINT(deltim) * kdt, DELOUT_NC) == 0) THEN
        CALL grid_collect_ipe(wwg,zzg,uug,vvg,
     &                        ttg,rqg,n2g,global_lats_a,lonsperlat,
     &                        lats_nodes_a,kdt,deltim,restart_step,
     &                        den,gmol)
      END IF

      END subroutine get_w_z
