      SUBROUTINE dfi_fixwr(iflag, sfc_fld, nst_fld)

!
!***********************************************************************
!     PURPOSE:
!      save or retrieve fixed fields in digifilt
!
!     REVISION HISOTRY:
!     2011-12-05  J.Wang    Adopted from fixwr in /nwprod/sorc/global_fcst.fd
!     2014-09-18  S.Moorthi cleaned up and made nst model an option
!
!***********************************************************************
!
      use resol_def
      use namelist_physics_def, only : nstf_name
      use layout1
      use gfs_physics_sfc_flx_mod, only : Sfc_Var_Data
      use gfs_physics_nst_var_mod, ONLY : Nst_Var_Data

      implicit none

      TYPE(Sfc_Var_Data)        :: sfc_fld
      TYPE(Nst_Var_Data)        :: nst_fld

      integer,intent(in)        :: iflag

      real , allocatable :: SMC1(:,:,:), STC1(:,:,:),
     &                      HICE1(:,:),  FICE1(:,:),   TISFC1(:,:),
     &                      TSEA1(:,:),  weasd1(:,:),  TG31(:,:),
     &                      ZORL1(:,:),  CV1(:,:),     CVB1(:,:),
     &                      CVT1(:,:),   ALVSF1(:,:),  ALVWF1(:,:),
     &                      ALNSF1(:,:), ALNWF1(:,:),  SLMSK1(:,:),
     &                      VFRAC1(:,:), CANOPY1(:,:), F10M1(:,:),
     &                      VTYPE1(:,:), STYPE1(:,:),  FACSF1(:,:),
     &                      FACWF1(:,:), UUSTAR1(:,:), FFMM1(:,:),
     &                      FFHH1(:,:),  TPRCP1(:,:),  SRFLAG1(:,:),
     +                      SLC1(:,:,:), SNWDPH1(:,:), SLOPE1(:,:),
     +                      SHDMIN1(:,:),SHDMAX1(:,:), SNOALB1(:,:),
     &                      SNCOVR1(:,:),

! li added for NST components
     &                      xt1(:,:),xs1(:,:),xu1(:,:),xv1(:,:),xz1(:,:)
     &,                     zm1(:,:),xtts1(:,:),xzts1(:,:),dt_cool1(:,:)
     &,                     z_c1(:,:),c_01(:,:),c_d1(:,:),w_01(:,:)
     &,                     w_d1(:,:),d_conv1(:,:),ifd1(:,:),Tref1(:,:)
     &,                     Qrain1(:,:)

      logical first
      data    first/.true./
      save    first,SMC1,STC1,TSEA1,weasd1,TG31,ZORL1,CV1,CVB1,CVT1
      save    HICE1,FICE1,TISFC1                       ! FOR SEA-ICE - XW Nov04
!                                                      ! FOR NST   - XL Dec0r97
      save    xt1,xs1,xu1,xv1,xz1,zm1,xtts1,xzts1,
     &        dt_cool1,z_c1,c_01,c_d1,w_01,w_d1,
     &        d_conv1,ifd1,Tref1,Qrain1
      save    ALVSF1,ALVWF1,ALNSF1,ALNWF1,SLMSK1,VFRAC1,CANOPY1,F10M1
     &,       VTYPE1,STYPE1,FACSF1,FACWF1,UUSTAR1,FFMM1,FFHH1
!lu [+2L]: save (tprcp1,srflag1),(slc1,snwdph1,slope1,shdmin1,shdmax1,snoalb1)
     &,       TPRCP1,SRFLAG1
     &,       SLC1,SNWDPH1,SLOPE1,SHDMIN1,SHDMAX1,SNOALB1,SNCOVR1

      integer i,j,k

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!     print *,' enter fixwr '                                   ! hmhj
      if (first) then
        allocate (SMC1(lonr,lsoil,lats_node_r))
        allocate (STC1(lonr,lsoil,lats_node_r))
        allocate (HICE1(lonr,lats_node_r))              ! FOR SEA-ICE - XW Nov04
        allocate (FICE1(lonr,lats_node_r))              ! FOR SEA-ICE - XW Nov04
        allocate (TISFC1(lonr,lats_node_r))             ! FOR SEA-ICE - XW Nov04
                                                        ! FOR NST     - XL Dec09
        allocate (TSEA1(lonr,lats_node_r))
        allocate (weasd1(lonr,lats_node_r))
        allocate (TG31(lonr,lats_node_r))
        allocate (ZORL1(lonr,lats_node_r))
        allocate (CV1(lonr,lats_node_r))
        allocate (CVB1(lonr,lats_node_r))
        allocate (CVT1(lonr,lats_node_r))
        allocate (ALVSF1(lonr,lats_node_r))
        allocate (ALVWF1(lonr,lats_node_r))
        allocate (ALNSF1(lonr,lats_node_r))
        allocate (ALNWF1(lonr,lats_node_r))
        allocate (SLMSK1(lonr,lats_node_r))
        allocate (VFRAC1(lonr,lats_node_r))
        allocate (CANOPY1(lonr,lats_node_r))
        allocate (F10M1(lonr,lats_node_r))
        allocate (VTYPE1(lonr,lats_node_r))
        allocate (STYPE1(lonr,lats_node_r))
        allocate (FACSF1(lonr,lats_node_r))
        allocate (FACWF1(lonr,lats_node_r))
        allocate (UUSTAR1(lonr,lats_node_r))
        allocate (FFMM1(lonr,lats_node_r))
        allocate (FFHH1(lonr,lats_node_r))
Clu [+8L]: allocate (tprcp,srflag),(slc,snwdph,slope,shdmin,shdmax,snoalb)
        allocate (TPRCP1(lonr,lats_node_r))
        allocate (SRFLAG1(lonr,lats_node_r))
        allocate (SLC1(lonr,lsoil,lats_node_r))
        allocate (SNWDPH1(lonr,lats_node_r))
        allocate (SLOPE1(lonr,lats_node_r))
        allocate (SHDMIN1(lonr,lats_node_r))
        allocate (SHDMAX1(lonr,lats_node_r))
        allocate (SNOALB1(lonr,lats_node_r))
        allocate (SNCOVR1(lonr,lats_node_r))

        if (nstf_name(1) > 0) then
          allocate (xt1(lonr,lats_node_r))
          allocate (xs1(lonr,lats_node_r))
          allocate (xu1(lonr,lats_node_r))
          allocate (xv1(lonr,lats_node_r))
          allocate (xz1(lonr,lats_node_r))
          allocate (zm1(lonr,lats_node_r))
          allocate (xtts1(lonr,lats_node_r))
          allocate (xzts1(lonr,lats_node_r))

          allocate (dt_cool1(lonr,lats_node_r))
          allocate (z_c1(lonr,lats_node_r))
          allocate (c_01(lonr,lats_node_r))
          allocate (c_d1(lonr,lats_node_r))
          allocate (w_01(lonr,lats_node_r))
          allocate (w_d1(lonr,lats_node_r))
          allocate (d_conv1(lonr,lats_node_r))
          allocate (ifd1(lonr,lats_node_r))
          allocate (Tref1(lonr,lats_node_r))
          allocate (Qrain1(lonr,lats_node_r))
        endif
        first = .false.
      endif
!
      if(iflag == 1) then
!      print *,'fixwr,iflag=1,asve 3hr phys'
        do k=1,lsoil
          do j=1,lats_node_r
            do i=1,lonr
              smc1(i,k,j) = sfc_fld%smc(k,i,j)
              stc1(i,k,j) = sfc_fld%stc(k,i,j)
              slc1(i,k,j) = sfc_fld%slc(k,i,j)        !! Clu [+1L]: slc -> slc1
            enddo
          enddo
        enddo
        do j=1,lats_node_r
          do i=1,lonr
            hice1(i,j)    = sfc_fld%hice(i,j)           ! FOR SEA-ICE - XW Nov04
            fice1(i,j)    = sfc_fld%fice(i,j)           ! FOR SEA-ICE - XW Nov04
            tisfc1(i,j)   = sfc_fld%tisfc(i,j)          ! FOR SEA-ICE - XW Nov04
            tsea1(i,j)    = sfc_fld%tsea(i,j)
            weasd1(i,j)   = sfc_fld%weasd(i,j)
            tg31(i,j)     = sfc_fld%tg3(i,j)
            zorl1(i,j)    = sfc_fld%zorl(i,j)
            cv1(i,j)      = sfc_fld%cv(i,j)
            cvb1(i,j)     = sfc_fld%cvb(i,j)
            cvt1(i,j)     = sfc_fld%cvt(i,j)
            alvsf1(i,j)   = sfc_fld%alvsf(i,j)
            alvwf1(i,j)   = sfc_fld%alvwf(i,j)
            alnsf1(i,j)   = sfc_fld%alnsf(i,j)
            alnwf1(i,j)   = sfc_fld%alnwf(i,j)
            slmsk1(i,j)   = sfc_fld%slmsk(i,j)
            vfrac1(i,j)   = sfc_fld%vfrac(i,j)
            canopy1(i,j)  = sfc_fld%canopy(i,j)
            f10m1(i,j)    = sfc_fld%f10m(i,j)
            vtype1(i,j)   = sfc_fld%vtype(i,j)
            stype1(i,j)   = sfc_fld%stype(i,j)
            facsf1(i,j)   = sfc_fld%facsf(i,j)
            facwf1(i,j)   = sfc_fld%facwf(i,j)
            uustar1(i,j)  = sfc_fld%uustar(i,j)
            ffmm1(i,j)    = sfc_fld%ffmm(i,j)
            ffhh1(i,j)    = sfc_fld%ffhh(i,j)
Clu [+7L]: add (tprcp,srflag),(snwdph,slope,shdmin,shdmax,snoalb)
            tprcp1(i,j)   = sfc_fld%tprcp(i,j)
            srflag1(i,j)  = sfc_fld%srflag(i,j)
            snwdph1(i,j)  = sfc_fld%snwdph(i,j)
            slope1(i,j)   = sfc_fld%slope(i,j)
            shdmin1(i,j)  = sfc_fld%shdmin(i,j)
            shdmax1(i,j)  = sfc_fld%shdmax(i,j)
            snoalb1(i,j)  = sfc_fld%snoalb(i,j)
            sncovr1(i,j)  = sfc_fld%sncovr(i,j)
          enddo
        enddo
        if (nstf_name(1) > 0 ) then                     ! For NST
          do j=1,lats_node_r
            do i=1,lonr
              xs1(i,j)      = nst_fld%xs(i,j)
              xu1(i,j)      = nst_fld%xu(i,j)
              xv1(i,j)      = nst_fld%xv(i,j)
              xz1(i,j)      = nst_fld%xz(i,j)
              zm1(i,j)      = nst_fld%zm(i,j)
              xtts1(i,j)    = nst_fld%xtts(i,j)
              xzts1(i,j)    = nst_fld%xzts(i,j)

              dt_cool1(i,j) = nst_fld%dt_cool(i,j)
              z_c1(i,j)     = nst_fld%z_c(i,j)
              c_01(i,j)     = nst_fld%c_0(i,j)
              c_d1(i,j)     = nst_fld%c_d(i,j)
              w_01(i,j)     = nst_fld%w_0(i,j)
              w_d1(i,j)     = nst_fld%w_d(i,j)
              d_conv1(i,j)  = nst_fld%d_conv(i,j)
              Tref1(i,j)    = nst_fld%Tref(i,j)
              Qrain1(i,j)   = nst_fld%Qrain(i,j)
            enddo
          enddo
        endif
      elseif(iflag == 2) then
!      print *,'fixwr,iflag=2,set 3hr phys'
        do k=1,lsoil
          do j=1,lats_node_r
            do i=1,lonr
              sfc_fld%smc(k,i,j) = smc1(i,k,j)
              sfc_fld%stc(k,i,j) = stc1(i,k,j)
              sfc_fld%slc(k,i,j) = slc1(i,k,j)         !! Clu [+1L]: slc1 -> slc
            enddo
          enddo
        enddo
        do j=1,lats_node_r
          do i=1,lonr
            sfc_fld%hice(i,j)    = hice1(i,j)          ! FOR SEA-ICE - XW Nov04
            sfc_fld%fice(i,j)    = fice1(i,j)          ! FOR SEA-ICE - XW Nov04
            sfc_fld%tisfc(i,j)   = tisfc1(i,j)         ! FOR SEA-ICE - XW Nov04

            sfc_fld%tsea(i,j)    = tsea1(i,j)
            sfc_fld%weasd(i,j)   = weasd1(i,j)
            sfc_fld%tg3(i,j)     = tg31(i,j)
            sfc_fld%zorl(i,j)    = zorl1(i,j)
            sfc_fld%cv(i,j)      = cv1(i,j)
            sfc_fld%cvb(i,j)     = cvb1(i,j)
            sfc_fld%cvt(i,j)     = cvt1(i,j)
            sfc_fld%alvsf(i,j)   = alvsf1(i,j)
            sfc_fld%alvwf(i,j)   = alvwf1(i,j)
            sfc_fld%alnsf(i,j)   = alnsf1(i,j)
            sfc_fld%alnwf(i,j)   = alnwf1(i,j)
            sfc_fld%slmsk(i,j)   = slmsk1(i,j)
            sfc_fld%vfrac(i,j)   = vfrac1(i,j)
            sfc_fld%canopy(i,j)  = canopy1(i,j)
            sfc_fld%f10m(i,j)    = f10m1(i,j)
            sfc_fld%vtype(i,j)   = vtype1(i,j)
            sfc_fld%stype(i,j)   = stype1(i,j)
            sfc_fld%facsf(i,j)   = facsf1(i,j)
            sfc_fld%facwf(i,j)   = facwf1(i,j)
            sfc_fld%uustar(i,j)  = uustar1(i,j)
            sfc_fld%ffmm(i,j)    = ffmm1(i,j)
            sfc_fld%ffhh(i,j)    = ffhh1(i,j)
Clu [+7L]: add (tprcp,srflag),(snwdph,slope,shdmin,shdmax,snoalb)
            sfc_fld%tprcp(i,j)   = tprcp1(i,j)
            sfc_fld%srflag(i,j)  = srflag1(i,j)
            sfc_fld%snwdph(i,j)  = snwdph1(i,j)
            sfc_fld%slope(i,j)   = slope1(i,j)
            sfc_fld%shdmin(i,j)  = shdmin1(i,j)
            sfc_fld%shdmax(i,j)  = shdmax1(i,j)
            sfc_fld%snoalb(i,j)  = snoalb1(i,j)
            sfc_fld%sncovr(i,j)  = sncovr1(i,j)
          enddo
        enddo
        if (nstf_name(1) > 0 ) then                     ! For NST
          do j=1,lats_node_r
            do i=1,lonr
              nst_fld%xs(i,j)      = xs1(i,j)
              nst_fld%xu(i,j)      = xu1(i,j)
              nst_fld%xv(i,j)      = xv1(i,j)
              nst_fld%xz(i,j)      = xz1(i,j)
              nst_fld%zm(i,j)      = zm1(i,j)
              nst_fld%xtts(i,j)    = xtts1(i,j)
              nst_fld%xzts(i,j)    = xzts1(i,j)

              nst_fld%dt_cool(i,j) = dt_cool1(i,j)
              nst_fld%z_c(i,j)     = z_c1(i,j)
              nst_fld%c_0(i,j)     = c_01(i,j)
              nst_fld%c_d(i,j)     = c_d1(i,j)
              nst_fld%w_0(i,j)     = w_01(i,j)
              nst_fld%w_d(i,j)     = w_d1(i,j)
              nst_fld%d_conv(i,j)  = d_conv1(i,j)
              nst_fld%Tref(i,j)    = Tref1(i,j)
              nst_fld%Qrain(i,j)   = Qrain1(i,j)
            enddo
          enddo
        endif
      endif
!     print *,' leave fixwr '
      return
      end subroutine dfi_fixwr
