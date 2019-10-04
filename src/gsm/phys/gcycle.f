      SUBROUTINE GCYCLE(ME,LATS_NODE_R,LONSPERLAR,GLOBAL_LATS_R,
     &                  IPT_LATS_NODE_R,IDATE,FHOUR,FHCYC,
     &                  XLON ,XLAT,sfc_fld,ialb)
!
      USE MACHINE
      USE PHYSCONS, PI => con_PI
      use resol_def
      use gfs_physics_sfc_flx_mod
      use namelist_physics_def, only : use_ufo, nst_anl
      implicit none
!
      TYPE(Sfc_Var_Data)        :: sfc_fld
!
      INTEGER LONSPERLAR(LATR),LONS_LAT,IALB
      INTEGER GLOBAL_LATS_R(LATR),IPT_LATS_NODE_R
      INTEGER ME,NVAR,IDATE(4),LATS_NODE_R,LAT
      real fhour, fhcyc
      real XLON  (LONR,LATS_NODE_R),XLAT  (LONR,LATS_NODE_R)

!
!     Local variables
!     ---------------
      integer ilat,len,il,iq,ilon,i,j,lon,l
!
      real  RLA(LONR*LATS_NODE_R),           RLO(LONR*LATS_NODE_R),
     &      SLMASK(LONR*LATS_NODE_R),        OROG(LONR*LATS_NODE_R),
     &                                       OROG_UF(LONR*LATS_NODE_R),
     &      TSFFCS(LONR*LATS_NODE_R),        SNOFCS(LONR*LATS_NODE_R),
     &      ZORFCS(LONR*LATS_NODE_R),        ALBFCS(LONR*LATS_NODE_R,4),
     &      TG3FCS(LONR*LATS_NODE_R),        CNPFCS(LONR*LATS_NODE_R),
     &      SMCFCS(LONR*LATS_NODE_R,LSOIL),
     &      STCFCS(LONR*LATS_NODE_R,LSOIL),
     &      SLIFCS(LONR*LATS_NODE_R),        AISFCS(LONR*LATS_NODE_R),
     &      F10MFCS(LONR*LATS_NODE_R),       VEGFCS(LONR*LATS_NODE_R),
     &      VETFCS(LONR*LATS_NODE_R),        SOTFCS(LONR*LATS_NODE_R),
     &      ALFFCS(LONR*LATS_NODE_R,2),      CVFCS(LONR*LATS_NODE_R),
     &      CVBFCS(LONR*LATS_NODE_R),        CVTFCS(LONR*LATS_NODE_R),
!
     &      SMCFC1(LONR*LATS_NODE_R*LSOIL),
     &      STCFC1(LONR*LATS_NODE_R*LSOIL),
     &      ALBFC1(LONR*LATS_NODE_R*4),       ALFFC1(LONR*LATS_NODE_R*2)
!CluX add swdfcs, sihfcs, sicfcs
     +,     SWDFCS(LONR*LATS_NODE_R)
     +,     SIHFCS(LONR*LATS_NODE_R),SICFCS(LONR*LATS_NODE_R)
     &,     SITFCS(LONR*LATS_NODE_R)
!CluX add vmnfcs, vmxfcs, slpfcs, absfcs, slcfc1, slcfcs
     +,     VMNFCS(LONR*LATS_NODE_R),VMXFCS(LONR*LATS_NODE_R)
     +,     SLPFCS(LONR*LATS_NODE_R),ABSFCS(LONR*LATS_NODE_R)
     +,     SLCFC1(LONR*LATS_NODE_R*LSOIL)
     +,     SLCFCS(LONR*LATS_NODE_R,LSOIL)


      real  sig1t, pifac
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!     if (me .eq. 0) print *,' nlats=',nlats,' lonsinpe='
!    *,lonsinpe(0,1)

      sig1t = 0.0
!
      pifac = 180.0 / pi
      len = 0
      DO ilat=1,LATS_NODE_R  !-----BEGIN LATITUDE LOOP------------------
        LAT = GLOBAL_LATS_R(IPT_LATS_NODE_R-1+ILAT)
        LONS_LAT = LONSPERLAR(LAT)
        DO ILON=1,LONS_LAT   !-----BEGIN LONGITUDE LOOP-----------------

!     print *,' calling gcycle for ilat',ilat,' me=',me,' nlats='
!    *,nlats,' lonsinpe=',lonsinpe(:,ilat)
!     if (ilat .eq. nlats) stop
!
          len = len + 1
          RLA(len)      = XLAT(ilon,ilat) * pifac
          RLO(len)      = XLON(ilon,ilat) * pifac
          OROG(len)     = sfc_fld%ORO(ilon,ilat)
          OROG_UF(len)  = sfc_fld%ORO_UF(ilon,ilat)
          TSFFCS(len)   = sfc_fld%TSEA(ilon,ilat)
          SNOFCS(len)   = sfc_fld%weasd(ilon,ilat)
          ZORFCS(len)   = sfc_fld%ZORL(ilon,ilat)
          ALBFCS(len,1) = sfc_fld%ALVSF(ilon,ilat)
          ALBFCS(len,2) = sfc_fld%ALVWF(ilon,ilat)
          ALBFCS(len,3) = sfc_fld%ALNSF(ilon,ilat)
          ALBFCS(len,4) = sfc_fld%ALNWF(ilon,ilat)
          TG3FCS(len)   = sfc_fld%TG3(ilon,ilat)
          CNPFCS(len)   = sfc_fld%CANOPY(ilon,ilat)
          SMCFCS(len,:) = sfc_fld%SMC(:,ilon,ilat)
          STCFCS(len,:) = sfc_fld%STC(:,ilon,ilat)
          SLIFCS(len)   = sfc_fld%SLMSK(ilon,ilat)
          F10MFCS(len)  = sfc_fld%F10M(ilon,ilat)
          VEGFCS(len)   = sfc_fld%VFRAC(ilon,ilat)
          VETFCS(len)   = sfc_fld%VTYPE(ilon,ilat)
          SOTFCS(len)   = sfc_fld%STYPE(ilon,ilat)
          ALFFCS(len,1) = sfc_fld%FACSF(ilon,ilat)
          ALFFCS(len,2) = sfc_fld%FACWF(ilon,ilat)
          CVFCS(len)    = sfc_fld%CV(ilon,ilat)
          CVBFCS(len)   = sfc_fld%CVB(ilon,ilat)
          CVTFCS(len)   = sfc_fld%CVT(ilon,ilat)
!CluX add swdfcs, sihfcs, sicfcs
          SWDFCS(len)   = sfc_fld%SNWDPH(ilon,ilat)
          SIHFCS(len)   = sfc_fld%HICE(ilon,ilat)
          SICFCS(len)   = sfc_fld%FICE(ilon,ilat)
          SITFCS(len)   = sfc_fld%TISFC(ilon,ilat)
!CluX add slcfcs, vmnfcs, vmxfcs, slpfcs, absfcs
          SLCFCS(len,:) = sfc_fld%SLC(:,ilon,ilat)
          VMNFCS(len)   = sfc_fld%SHDMIN(ilon,ilat)
          VMXFCS(len)   = sfc_fld%SHDMAX(ilon,ilat)
          SLPFCS(len)   = sfc_fld%SLOPE(ilon,ilat)
          ABSFCS(len)   = sfc_fld%SNOALB(ilon,ilat)

!
          IF (SLIFCS(len) .LT. 0.1 .OR. SLIFCS(len) .GT. 1.5) THEN
             SLMASK(len) = 0
          ELSE
             SLMASK(len) = 1
          ENDIF

          IF (SLIFCS(len) .EQ. 2) THEN
            AISFCS(len) = 1.
          ELSE
            AISFCS(len) = 0.
          ENDIF

!     if (me .eq. 0)
!    &   print *,' len=',len,' rla=',rla(len),' rlo=',rlo(len)
        ENDDO                 !-----END LONGITUDE LOOP------------------------------

      ENDDO                   !-----END LATITUDE LOOP-------------------------------
!
! check
!     print *,' total points = ',len
!
      do l=1,lsoil
        il = (l-1)*len
        do i=1,len
          SMCFC1(il+i) = SMCFCS(i,l)
          STCFC1(il+i) = STCFCS(i,l)
          SLCFC1(il+i) = SLCFCS(i,l)
        enddo
      enddo
      do l=1,4
        il = (l-1)*len
        do i=1,len
          ALBFC1(il+i) = ALBFCS(i,l)
        enddo
      enddo
      do l=1,2
        il = (l-1)*len
        do i=1,len
          ALFFC1(il+i) = ALFFCS(i,l)
        enddo
      enddo
! check
!     call mymaxmin(slifcs,len,len,1,'slifcs')
!     call mymaxmin(slmask,len,len,1,'slmsk')
!
      CALL SFCCYCLE(101,LEN,LSOIL,SIG1T,fhcyc
     &,             idate(4), idate(2), idate(3), idate(1), fhour
     &,             RLA, RLO, SLMASK, OROG, OROG_UF, USE_UFO, nst_anl
     &,             SIHFCS,   SICFCS, SITFCS
     &,             SWDFCS,   SLCFC1
     &,             VMNFCS,   VMXFCS, SLPFCS, ABSFCS
     &,             TSFFCS,   SNOFCS, ZORFCS, ALBFC1, TG3FCS
     &,             CNPFCS,   SMCFC1, STCFC1, SLIFCS, AISFCS, F10MFCS
     &,             VEGFCS,   VETFCS, SOTFCS, ALFFC1
     &,             CVFCS,    CVBFCS, CVTFCS, me, nlunit, ialb)
!
      do l=1,lsoil
        il = (l-1)*len
        do i=1,len
          SMCFCS(i,l) = SMCFC1(il+i)
          STCFCS(i,l) = STCFC1(il+i)
          SLCFCS(i,l) = SLCFC1(il+i)
        enddo
      enddo
      do l=1,4
        il = (l-1)*len
        do i=1,len
          ALBFCS(i,l) = ALBFC1(il+i)
        enddo
      enddo
      do l=1,2
        il = (l-1)*len
        do i=1,len
          ALFFCS(i,l) = ALFFC1(il+i)
        enddo
      enddo
!
      il = 0
      DO ILAT=1,LATS_NODE_R  !-----BEGIN LATITUDE LOOP------------------
        LAT = GLOBAL_LATS_R(IPT_LATS_NODE_R-1+ILAT)
        LONS_LAT = LONSPERLAR(LAT)
!
        DO ILON=1,LONS_LAT   !-----BEGIN LONGITUDE LOOP-----------------
          il = il + 1
          sfc_fld%TSEA(ilon,ilat)   = TSFFCS(il)
          sfc_fld%weasd(ilon,ilat)  = SNOFCS(il)
          sfc_fld%ZORL(ilon,ilat)   = ZORFCS(il)
          sfc_fld%ALVSF(ilon,ilat)  = ALBFCS(il,1)
          sfc_fld%ALVWF(ilon,ilat)  = ALBFCS(il,2)
          sfc_fld%ALNSF(ilon,ilat)  = ALBFCS(il,3)
          sfc_fld%ALNWF(ilon,ilat)  = ALBFCS(il,4)
          sfc_fld%TG3(ilon,ilat)    = TG3FCS(il)
          sfc_fld%CANOPY(ilon,ilat) = CNPFCS(il)
          sfc_fld%SMC(:,ilon,ilat)  = SMCFCS(il,:)
          sfc_fld%STC(:,ilon,ilat)  = STCFCS(il,:)
          sfc_fld%SLMSK(ilon,ilat)  = SLIFCS(il)
          sfc_fld%F10M(ilon,ilat)   = F10MFCS(il)
          sfc_fld%VFRAC(ilon,ilat)  = VEGFCS(il)
          sfc_fld%VTYPE(ilon,ilat)  = VETFCS(il)
          sfc_fld%STYPE(ilon,ilat)  = SOTFCS(il)
          sfc_fld%FACSF(ilon,ilat)  = ALFFCS(il,1)
          sfc_fld%FACWF(ilon,ilat)  = ALFFCS(il,2)
          sfc_fld%CV(ilon,ilat)     = CVFCS(il)
          sfc_fld%CVB(ilon,ilat)    = CVBFCS(il)
          sfc_fld%CVT(ilon,ilat)    = CVTFCS(il)
          sfc_fld%SNWDPH(ilon,ilat) = SWDFCS(il)
          sfc_fld%HICE(ilon,ilat)   = SIHFCS(il)
          sfc_fld%FICE(ilon,ilat)   = SICFCS(il)
          sfc_fld%TISFC(ilon,ilat)  = SITFCS(il)
          sfc_fld%SLC(:,ilon,ilat)  = SLCFCS(il,:)
          sfc_fld%SHDMIN(ilon,ilat) = VMNFCS(il)
          sfc_fld%SHDMAX(ilon,ilat) = VMXFCS(il)
          sfc_fld%SLOPE(ilon,ilat)  = SLPFCS(il)
          sfc_fld%SNOALB(ilon,ilat) = ABSFCS(il)
!
        ENDDO     !-----END LONGITUDE LOOP------------------------------
!
      ENDDO       !-----END LATITUDE LOOP-------------------------------
!
!     if (me .eq. 0) print*,'executed gcycle during hour=',fhour
      
      RETURN
      END

