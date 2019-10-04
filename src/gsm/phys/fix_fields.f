
      SUBROUTINE fix_fields(
     &                  LONSPERLAR,GLOBAL_LATS_R,XLON,XLAT,sfc_fld,
     &                  nst_fld,HPRIME,JINDX1,JINDX2,DDY,OZPLIN,CREAD,
     &                  cread_grd,cread_nst,nblck,phy_f3d,phy_f2d)
!!     
!! Code Revision
!! jan 26 2010 Jun Wang, added phy_f3d,phy_f2d read in from restart file
!! Nov    2010 S. Moorthi - nst model related changes
!! Mar    2013 Jun Wang  restart for idea
!! Aug    2015 Xu Li     change nst_fcst and nst_spinup to be nstf_name
!!                       introduce the depth mean SST, 
!!                       remove nst_reset_nonwater
!!                       add nemsio for nst file

!!
      use machine , only : kind_rad, kind_phys
      use funcphys                         
      use module_progtm             
      use resol_def
      use namelist_physics_def
      use layout1
      use gg_def
      use ozne_def
      use module_nst_water_prop, only: get_dtzm_2d
      use gfs_physics_sfc_flx_mod
      use gfs_physics_nst_var_mod
      use idea_composition,             only: pr_idea,gg,amgms,prsilvl
      IMPLICIT NONE
!!     
      TYPE(Sfc_Var_Data)  :: sfc_fld
      TYPE(Nst_Var_Data)  :: nst_fld
      CHARACTER (len=*)   :: CREAD
      CHARACTER (len=*)   :: CREAD_grd
      CHARACTER (len=*)   :: cread_nst
      INTEGER JINDX1(LATS_NODE_R),JINDX2(LATS_NODE_R)
      REAL (KIND=KIND_RAD) DDY(LATS_NODE_R)
      REAL (KIND=KIND_RAD) HPRIME(NMTVR,LONR,LATS_NODE_R)

      INTEGER IOZONDP
      REAL (kind=kind_rad) OZPLIN(LATSOZP,LEVOZP,pl_coeff,timeoz)
     &,                    XLON(LONR,LATS_NODE_R)
     &,                    XLAT(LONR,LATS_NODE_R)
     &,                    dtzm(lonr,lats_node_r)
       
      integer, dimension(latr) :: global_lats_r, lonsperlar
!
      integer              nblck
      real (kind=kind_rad) phy_f3d(NGPTC,LEVS,ntot3d,NBLCK,LATS_NODE_R)
     &,                    phy_f2d(lonr,lats_node_r,ntot2d)
      REAL (KIND=KIND_RAD) plyr(levs)
!
      real (kind=kind_phys) gaul(lats_node_r), pi
      real (kind=kind_phys) :: zsea1,zsea2

      real, PARAMETER:: RLAPSE=0.65E-2
      real dt_warm
      integer needoro, i, j, lat, k
      INTEGER NREAD, NREAD_NST
!!     
      call gfuncphys
      if (lsm == 0) then ! For OSU LSM
         CALL GRDDF
         CALL GRDKT
      endif
!!     
      IOZONDP = 0
      if (ntoz > 0) IOZONDP = 1
      NREAD   = 14
!     CREAD   = 'fort.14'
      sfc_fld%ORO = 0.
      NEEDORO     = 0

      if (fhini == fhrot) then
        if (me == 0) print *,' call read_sfc CREAD=',cread
     &,   'nemsio_in=',nemsio_in
        if(nemsio_in) then
          CALL read_sfc_nemsio(sfc_fld,NEEDORO,NREAD,
     &                         CREAD,GLOBAL_LATS_R,LONSPERLAR)
!!read in vcoord
          call read_vcoord(NREAD,CREAD_grd)
        else
!sfcio input
          call read_sfc(sfc_fld,NEEDORO,NREAD,
     &                  CREAD,GLOBAL_LATS_R,LONSPERLAR)
        endif
!nst   input
        if (nstf_name(1) > 0) then
          nst_fld%slmsk = sfc_fld%slmsk
          if ( nstf_name(2) > 0 ) then
          if (me == 0) print *,'call set_nst when nst_spinup = ',
     &                           nstf_name(2)
            CALL set_nst(sfc_fld%tsea,nst_fld)
          else
            NREAD_NST   = 15
!
!         Add nemsio for nst file here in the future
!
            if(nemsio_in) then
              if (me == 0) print *,'call read_nst_nemsio'
              CALL read_nst_nemsio(nst_fld,NREAD_NST,CREAD_NST,
     &                      GLOBAL_LATS_R,LONSPERLAR)
            else
              if (me == 0) print *,'call read_nst'
              CALL read_nst(nst_fld,NREAD_NST,CREAD_NST,
     &                      GLOBAL_LATS_R,LONSPERLAR)
            endif

!
!           Get SST (Tf + dTw - dTc) if NSSTm coupled
!
            if (  nstf_name(1) > 1 ) then   ! if NSSTM coupled on
              zsea1=0.001*real(nstf_name(4))
              zsea2=0.001*real(nstf_name(5))
              call get_dtzm_2d(nst_fld%xt,nst_fld%xz,nst_fld%dt_cool,
     &                         nst_fld%z_c,sfc_fld%slmsk,
     &                         zsea1,zsea2,lonr,lats_node_r,dtzm)
              do j = 1, lats_node_r
                do i = 1, lonr
                  if (sfc_fld%slmsk(i,j) == 0 ) then
                    sfc_fld%tsea(i,j) = nst_fld%tref(i,j) + dtzm(i,j)
     &              - (sfc_fld%oro(i,j)-sfc_fld%oro_uf(i,j)) * rlapse
                  endif
                enddo
              enddo
!
!    When AM and NST is not coupled, tsea (in surface file) ==> tref
!
            elseif (nstf_name(1) == 1) then
              nst_fld%tref = sfc_fld%tsea
            endif
          endif                        ! if ( nstf_name(2) > 0 ) then
        endif                          ! if ( nstf_name(1) > 0 ) then
      else
!
        if(lsidea) then
          if(.not.allocated(pr_idea)) then
            allocate(pr_idea(levs))
            allocate(gg(levs))
            allocate(prsilvl(levs+1))
          endif
        endif
!
        if (me .eq. 0) print *,' call read_sfc_r CREAD=',cread
        CALL read_sfc_r(cread,sfc_fld,phy_f2d,phy_f3d,
     &                  NGPTC,NBLCK,GLOBAL_LATS_R,LONSPERLAR,
     &                  NEEDORO,lsidea,pr_idea,gg,prsilvl,amgms)

        if ( nstf_name(1) > 0 ) then
          nst_fld%slmsk = sfc_fld%slmsk
          NREAD_NST   = 15
          if (me == 0) print *,' call read_nst_r CREAD=',cread_nst
          CALL read_nst_r(nst_fld,NREAD_NST,CREAD_NST,
     &                    GLOBAL_LATS_R,LONSPERLAR)
        endif
      endif
      NEEDORO=1
      CALL read_mtn_hprim_oz(sfc_fld%SLMSK,HPRIME,NEEDORO,sfc_fld%ORO,
     &                       sfc_fld%oro_uf,IOZONDP,OZPLIN,
     &                       GLOBAL_LATS_R,LONSPERLAR)
!      
!   Set up some interpolation coefficients for ozone forcing
!
      if (ntoz > 0) then
        pi = acos(-1.0)
        do j=1, lats_node_r
          lat = global_lats_r(ipt_lats_node_r-1+J)
          if (lat <= latr2) then
            gaul(j) = 90.0 - colrad_r(lat)*180.0/PI
          else
            gaul(j) = -(90.0 - colrad_r(lat)*180.0/PI)
          endif
!cselaif(me.eq.0) print*,'gau(j,1) gau(j,2)',gaul(j,1),gaul(j,2)
        enddo
        CALL SETINDXOZ(LATS_NODE_R,LATS_NODE_R,GAUL,
     &                 JINDX1,JINDX2,DDY)
      endif
!      
      CALL LONLAT_PARA(GLOBAL_LATS_R,XLON,XLAT,LONSPERLAR)
!!     
      RETURN
      END
