      SUBROUTINE read_mtn_hprim_oz(SLMSK,HPRIME,NEEDORO,ORO,oro_uf,
     &           iozondp,ozplin,global_lats_r,lonsperlar)
!
!***********************************************************************
!
      use resol_def, ONLY: latr, lonr, nmtvr, jcap
      use layout1,   ONLY: me, nodes, lats_node_r, nodes_comp
      use ozne_def,  ONLY: latsozp, levozp, timeoz, pl_coeff
      USE machine,   ONLY: kind_io8, kind_io4
      use namelist_physics_def, only: use_ufo, griboro
      implicit none

!
      integer iozondp
      integer, dimension(latr) :: global_lats_r, lonsperlar
      real (kind=kind_io8), dimension(lonr,lats_node_r) ::
     &                      slmsk, oro, oro_uf
!    &,                     buffo, buff2
      real (kind=kind_io8)  HPRIME(NMTVR,lonr,lats_node_r)
 
      real (kind=kind_io8) ozplin(latsozp,levozp,pl_coeff,timeoz)
 
      real (kind=kind_io8), allocatable :: buffo(:,:), buff2(:,:)
      real (kind=kind_io4), allocatable :: buff1(:,:), buffm(:,:,:)
      integer kmsk0(lonr,lats_node_r)
      integer i,j,k,nmtn,needoro,iop,iret
!
      kmsk0 = 0
!     iop = nodes_comp - 1
      iop = 0

      allocate (buffo(lonr,lats_node_r), buff2(lonr,lats_node_r))
      allocate (buff1(lonr,latr), buffm(lonr,latr,nmtvr))
!
!     Read HPRIME from file MTNVAR
!     ****************************
      nmtn = 24
      IF (me == iop) THEN   
        READ(nmtn) buffm
!!      do k=1,nmtvr
!!        write(200) buffm(:,:,k)
!!      enddo
      ENDIF
      DO k=1,nmtvr
        call split2d_phys(buffm(1,1,k),buffo,global_lats_r,iop)
        CALL interpred_phys(1,kmsk0,buffo,buff2,global_lats_r,
     &                      lonsperlar)
        do j=1,lats_node_r
          do i=1,lonr
            HPRIME(k,i,j) = buff2(i,j)
          enddo
        enddo
      ENDDO
 
!my jordan's mb
!sela  print *, ' (*j*)  nmtvr= ',nmtvr, 'reading hprime'
!my      DO j=1,lats_node_r
!my      DO i=1,lonr
!my      DO k=1,NMTVR
!my        IF(SLMSK(i,j).NE.1.) HPRIME(k,i,j) = 0.
!my      ENDDO
!my      ENDDO
!my      ENDDO
 

 
      if (iozondp == 1) CALL readoz_disprd(ozplin)
!
!     reading the grib orography and scattering the data
!
      if(needoro.eq.1) then

        IF( me == iop) then
          if (griboro) then
            CALL ORORD(101,lonr,latr,buff1,'orography')
            print *,'read grb orography'
          else
            open(101, file='orography', form='unformatted'
     &,                      status='old', iostat=iret)
            read(101) buff1
            print *,'read binary orography'
            close(101)
          endif
        endif
        call split2d_phys(buff1,buffo,global_lats_r,iop)
        CALL interpred_phys(1,kmsk0,buffo,oro,global_lats_r,lonsperlar)
      endif
!                                  read unfiltered orography
      if (use_ufo) then
        IF( me == iop) then
          if (griboro) then
            CALL ORORD(101,lonr,latr,buff1,'orography_uf')
            print *,'read grb orography_uf'
          else
            open(101, file='orography_uf', form='unformatted'
     &,                      status='old', iostat=iret)
            read(101) buff1
            print *,'read binary orography_uf'
            close(101)
          endif
        endif
        call split2d_phys(buff1,buffo,global_lats_r,iop)
        CALL interpred_phys(1,kmsk0,buffo,oro_uf,
     &                      global_lats_r,lonsperlar)
      else
        oro_uf = 0.0
      endif
      deallocate(buff1, buffm, buffo, buff2)

      RETURN
      END

      SUBROUTINE read_sfc(sfc_fld,NEEDORO,nread,
     &                    cfile,global_lats_r,lonsperlar)
!
!***********************************************************************
!
! !REVISION HISTORY:
!
!  2012/10/20  Jun Wang, Modified nwprod read_sfc.f(sfcio) for nems gfs
!
      use sfcio_module, ONLY: sfcio_head, sfcio_data, sfcio_realfill,
     &                        sfcio_srohdc, sfcio_axdata
      use resol_def,               ONLY: latr, latr2, lonr, lsoil
      use layout1,                 ONLY: me, nodes, lats_node_r
     &,                                  lats_node_r_max, nodes_comp
      use gfs_physics_sfc_flx_mod, ONLY: Sfc_Var_Data
      use namelist_soilveg ,       only: salp_data, snupx
      use physcons,                only: tgice => con_tice
      use machine,                 only: kind_io4, kind_io8

      implicit none
!
      TYPE(Sfc_Var_Data)        :: sfc_fld
      integer, dimension(latr)  :: global_lats_r, lonsperlar
      integer jump, needoro

!     real(kind=kind_io4) buff1(lonr,latr)
!     real(kind=kind_io8) buffo(lonr,lats_node_r)
!     real(kind=kind_io8) buff3(lonr,lats_node_r)

      real(kind=kind_io4), allocatable :: buff1(:,:)
      real(kind=kind_io8), allocatable :: buffo(:,:), buff3(:,:)
      integer nread,i,j,k,ij,idate(4),lonsfc,latsfc,lplsfc(latr2)
      character*(*) cfile
      integer kmsk(lonr,latr)
      CHARACTER*8 labfix(4)
      real t1,t2,timef,rsnow
      type(sfcio_head) head
      type(sfcio_data) data
      integer iret, vegtyp, iop
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!     t1=timef()

!     iop = nodes_comp - 1
      iop = 0
      allocate (buff1(lonr,latr), buffo(lonr,lats_node_r)
     &,         buff3(lonr,lats_node_r))
      IF (me == iop) THEN
        print *,'in read_sfcio nread=',nread,' cfile=',cfile
        call sfcio_srohdc(nread,cfile,head,data,iret)

        PRINT 99,nread,head%fhour,head%idate,
     &         head%lonb,head%latb,head%lsoil,head%ivs,iret,lats_node_r
99      FORMAT(1H ,'in fixio nread=',i3,2x,'HOUR=',f8.2,3x,'IDATE=',
     &  4(1X,I4),4x,'lonsfc,latsfc,lsoil,ivssfc,iret=',6i8)

        if(iret.ne.0) goto 5000
!       if(head%ivs.ne.200412.and.head%ivs.ne.200501) goto 5000
        if(head%lonb.ne.lonr) goto 5000
        if(head%latb.ne.latr) goto 5000
        if(head%lsoil.ne.lsoil) goto 5000

      ENDIF

      kmsk = 0

      if(me == iop) buff1 = data%tsea
      call split2d_phys(buff1, buffo,global_lats_r,iop)
      CALL interpred_phys(1,kmsk,buffo,sfc_fld%TSEA,
     &               global_lats_r,lonsperlar)

      DO K=1, LSOIL
        if(me == iop) buff1 = data%smc(:,:,k)
        call split2d_phys(buff1, buffo,global_lats_r,iop)
        CALL interpred_phys(1,kmsk,buffo,buff3,global_lats_r,lonsperlar)
        do j=1,lats_node_r
          do i=1,lonr
            sfc_fld%SMC(k,i,j) = buff3(i,j)
          enddo
        enddo
      ENDDO

      if(me == iop) buff1=data%weasd
      call split2d_phys(buff1, buffo,global_lats_r,iop)
      CALL interpred_phys(1,kmsk,buffo,sfc_fld%weasd,
     &               global_lats_r,lonsperlar)

      DO K = 1, LSOIL
        if(me == iop) buff1=data%stc(:,:,k)
        call split2d_phys(buff1, buffo,global_lats_r,iop)
        CALL interpred_phys(1,kmsk,buffo,buff3,global_lats_r,lonsperlar)
        do j=1,lats_node_r
          do i=1,lonr
            sfc_fld%STC(k,i,j) = buff3(i,j)
          enddo
        enddo
      ENDDO

      if(me == iop) buff1=data%tg3
      call split2d_phys(buff1, buffo,global_lats_r,iop)
      CALL interpred_phys(1,kmsk,buffo,sfc_fld%TG3,
     &               global_lats_r,lonsperlar)

      if(me == iop) buff1=data%zorl
      call split2d_phys(buff1, buffo,global_lats_r,iop)
      CALL interpred_phys(1,kmsk,buffo,sfc_fld%ZORL,
     &               global_lats_r,lonsperlar)

      sfc_fld%cv  = 0
      sfc_fld%cvb = 0
      sfc_fld%cvt = 0

      if(me == iop) buff1=data%alvsf
      call split2d_phys(buff1, buffo,global_lats_r,iop)
      CALL interpred_phys(1,kmsk,buffo,sfc_fld%ALVSF,
     &               global_lats_r,lonsperlar)
      if(me == iop) buff1=data%alvwf
      call split2d_phys(buff1, buffo,global_lats_r,iop)
      CALL interpred_phys(1,kmsk,buffo,sfc_fld%ALVWF,
     &               global_lats_r,lonsperlar)
      if(me == iop) buff1=data%alnsf
      call split2d_phys(buff1, buffo,global_lats_r,iop)
      CALL interpred_phys(1,kmsk,buffo,sfc_fld%ALNSF,
     &               global_lats_r,lonsperlar)
      if(me == iop) buff1=data%alnwf
      call split2d_phys(buff1, buffo,global_lats_r,iop)
      CALL interpred_phys(1,kmsk,buffo,sfc_fld%ALNWF,
     &               global_lats_r,lonsperlar)

!     The mask cannot be interpolated
      if(me == iop) buff1=data%slmsk
      call split2d_phys(buff1, buffo,global_lats_r,iop)
      CALL interpred_phys(1,kmsk,buffo,sfc_fld%SLMSK,
     &               global_lats_r,lonsperlar)

      if(me == iop) buff1=data%vfrac
      call split2d_phys(buff1, buffo,global_lats_r,iop)
      CALL interpred_phys(1,kmsk,buffo,sfc_fld%VFRAC,
     &               global_lats_r,lonsperlar)

      if(me == iop) buff1=data%canopy
      call split2d_phys(buff1, buffo,global_lats_r,iop)
      CALL interpred_phys(1,kmsk,buffo,sfc_fld%CANOPY,
     &               global_lats_r,lonsperlar)

      if(me == iop) buff1=data%f10m
      call split2d_phys(buff1, buffo,global_lats_r,iop)
      CALL interpred_phys(1,kmsk,buffo,sfc_fld%F10M,
     &               global_lats_r,lonsperlar)

      if(me == iop) buff1=data%vtype
      call split2d_phys(buff1, buffo,global_lats_r,iop)
      CALL interpred_phys(1,kmsk,buffo,sfc_fld%VTYPE,
     &               global_lats_r,lonsperlar)

      if(me == iop) buff1=data%stype
      call split2d_phys(buff1, buffo,global_lats_r,iop)
      CALL interpred_phys(1,kmsk,buffo,sfc_fld%STYPE,
     &               global_lats_r,lonsperlar)

      if(me == iop) buff1=data%facsf
      call split2d_phys(buff1, buffo,global_lats_r,iop)
      CALL interpred_phys(1,kmsk,buffo,sfc_fld%FACSF,
     &               global_lats_r,lonsperlar)
      if(me == iop) buff1=data%facwf
      call split2d_phys(buff1, buffo,global_lats_r,iop)
      CALL interpred_phys(1,kmsk,buffo,sfc_fld%FACWF,
     &               global_lats_r,lonsperlar)

!szunyogh 06/16/99
        if(me == iop) buff1=data%uustar
         call split2d_phys(buff1, buffo,global_lats_r,iop)
         CALL interpred_phys(1,kmsk,buffo,sfc_fld%UUSTAR,
     &               global_lats_r,lonsperlar)

        if(me == iop) buff1=data%ffmm
         call split2d_phys(buff1, buffo,global_lats_r,iop)
         CALL interpred_phys(1,kmsk,buffo,sfc_fld%FFMM,
     &                  global_lats_r,lonsperlar)

        if(me == iop) buff1=data%ffhh
         call split2d_phys(buff1, buffo,global_lats_r,iop)
         CALL interpred_phys(1,kmsk,buffo,sfc_fld%FFHH,
     &                  global_lats_r,lonsperlar)

!c-- XW: FOR SEA-ICE Nov04
!c-- XW: FOR SEA-ICE JAN15 change "sfc_fld%FICE(i,j)>=0.15"
!    Sea-ice (hice/fice) was added to the surface files.
         if(me == iop) buff1=data%hice
         call split2d_phys(buff1, buffo,global_lats_r,iop)
         CALL interpred_phys(1,kmsk,buffo,sfc_fld%HICE,
     &                  global_lats_r,lonsperlar)

         if(me == iop) buff1=data%fice
         call split2d_phys(buff1, buffo,global_lats_r,iop)
         CALL interpred_phys(1,kmsk,buffo,sfc_fld%FICE,
     &                  global_lats_r,lonsperlar)

         if(me == iop) buff1=data%tisfc
         call split2d_phys(buff1, buffo,global_lats_r,iop)
         CALL interpred_phys(1,kmsk,buffo,sfc_fld%TISFC,
     &                  global_lats_r,lonsperlar)
         if (lats_node_r > 0 )  then
           if (sfc_fld%tisfc(1,1) < 0.0)  then
             DO j=1,lats_node_r
!$omp parallel do private(i)
               DO i=1,LONR
                 sfc_fld%TISFC(i,j) = sfc_fld%TSEA(i,j)
                 IF(sfc_fld%SLMSK(i,j) >=  2. .AND.
!    &             sfc_fld%FICE(i,j)  >= 0.5) THEN
     &             sfc_fld%FICE(i,j)  >= 0.15) THEN
                   sfc_fld%TISFC(i,j) = (sfc_fld%TSEA(i,j)
     &            -tgice*(1.-sfc_fld%FICE(i,j))) / sfc_fld%FICE(i,j)
                   sfc_fld%TISFC(i,j)=MIN(sfc_fld%TISFC(i,j),tgice)
                 ENDIF
               ENDDO
             ENDDO
           endif
         endif

!-- XW: END SEA-ICE

!lu   11/10/2004
!*     surface files for GFS/Noah contain 8 additional records:
!*     tprcp, srflag, snwdph, slc, shdmin, shdmax, slope, snoalb

         if(me == iop) buff1=data%tprcp
         call split2d_phys(buff1, buffo,global_lats_r,iop)
         CALL interpred_phys(1,kmsk,buffo,sfc_fld%TPRCP,
     &                  global_lats_r,lonsperlar)

!* srflag
         if(me == iop) buff1=data%srflag
         call split2d_phys(buff1, buffo,global_lats_r,iop)
         CALL interpred_phys(1,kmsk,buffo,sfc_fld%SRFLAG,
     &                  global_lats_r,lonsperlar)

!* snwdph
         if(me == iop) buff1=data%snwdph
         call split2d_phys(buff1, buffo,global_lats_r,iop)
         CALL interpred_phys(1,kmsk,buffo,sfc_fld%SNWDPH,
     &                  global_lats_r,lonsperlar)

!* slc
         DO K=1, LSOIL
           if(me == iop) buff1=data%slc(:,:,k)
           call split2d_phys(buff1, buffo,global_lats_r,iop)
           CALL interpred_phys(1,kmsk,buffo,buff3,
     &                  global_lats_r,lonsperlar)
           do j=1,lats_node_r
             do i=1,lonr
               sfc_fld%SLC(k,i,j) = buff3(i,j)
             enddo
           enddo
         ENDDO

!* shdmin
         if(me == iop) buff1=data%shdmin
         call split2d_phys(buff1, buffo,global_lats_r,iop)
         CALL interpred_phys(1,kmsk,buffo,sfc_fld%SHDMIN,
     &                  global_lats_r,lonsperlar)

!* shdmax
         if(me == iop) buff1=data%shdmax
         call split2d_phys(buff1, buffo,global_lats_r,iop)
         CALL interpred_phys(1,kmsk,buffo,sfc_fld%SHDMAX,
     &                  global_lats_r,lonsperlar)

!* slope
         if(me == iop) buff1=data%slope
         call split2d_phys(buff1, buffo,global_lats_r,iop)
         CALL interpred_phys(1,kmsk,buffo,sfc_fld%SLOPE,
     &                  global_lats_r,lonsperlar)

!* snoalb
         if(me == iop) buff1=data%snoalb
         call split2d_phys(buff1, buffo,global_lats_r,iop)
         CALL interpred_phys(1,kmsk,buffo,sfc_fld%SNOALB,
     &                  global_lats_r,lonsperlar)
!     print *,' snoalb=',sfc_fld%snoalb(1,:)
!lu [+67L]: the addition of 8 Noah records ends here .........................

       if(needoro.eq.1) then
         if(me == iop) then
           buff1=data%orog
           needoro=1
           if(all(data%orog.ne.sfcio_realfill)) needoro=0
           print *,'read sfc orography'
         endif
         call split2d_phys(buff1, buffo,global_lats_r,iop)
         CALL interpred_phys(1,kmsk,buffo,sfc_fld%ORO,
     &                  global_lats_r,lonsperlar)
         call skip(needoro)
       endif
!
!Wei initialize snow fraction(weasd is in mm)
      DO j=1,lats_node_r
        DO i=1,LONR
          sfc_fld%SNCOVR(i,j) = 0.0
          if (sfc_fld%slmsk(i,j) > 0.001 ) then
            vegtyp = sfc_fld%VTYPE(i,j)
            if( vegtyp == 0 ) vegtyp = 7      !hmhj
!        if (vegtyp < 0 .or. vegtyp >13)
!    &   write(0,*)' vegtyp=',vegtyp, sfc_fld%VTYPE(i,j),' i=',' j=',j
!    &,' me=',me
            RSNOW  = 0.001*sfc_fld%weasd(i,j)/SNUPX(vegtyp)
            IF (0.001*sfc_fld%weasd(i,j) < SNUPX(vegtyp)) THEN
              sfc_fld%SNCOVR(i,j) = 1.0 - ( EXP(-SALP_DATA*RSNOW)
     &                                    - RSNOW*EXP(-SALP_DATA))
            ELSE
              sfc_fld%SNCOVR(i,j) = 1.0
            ENDIF
!           if (i == 1)
!    &       print*,SNUPX(vegtyp),SALP_DATA,sfc_fld%SNCOVR(i,j),
!    &       '************debug',sfc_fld%weasd(i,j),vegtyp,' j=',j
!    &,      ' snoalb1=',sfc_fld%snoalb(i,j)
!
          endif
        ENDDO
       ENDDO
!

      IF (me == iop) then
         call sfcio_axdata(data,iret)
!        t2=timef()
!        print *,'FIXIO TIME ',t2-t1,t1,t2
      endif
      deallocate (buff1, buffo, buff3)
!
      RETURN
 5000 PRINT *, ' ERROR IN INPUT IN FIXIO'
      STOP
      END

      SUBROUTINE read_vcoord(nread,cfile)
!
!***********************************************************************
!
! !REVISION HISTORY:
!
! 2016/03/01 Jun Wang read in vertical coordinate info from grid point 
!                     model lavel file
!
      USE machine,        only : kind_io4
      use coordinate_def, only : ak5,bk5,ck5
      use nemsio_module
!
      integer       nread
      character*(*) cfile
!
      integer idvc,k,dimz,iret
      logical gen_coord_hybrid, hybrid
      type(nemsio_gfile) gfile_in
      real(kind=kind_io4),allocatable :: vcoord(:,:,:)
!
!!read grid point model level file
      print *,' nread=',nread,' cfile=',cfile
      call nemsio_init()
!
      call nemsio_open(gfile_in,trim(cfile),'read',iret=iret)
!
      call nemsio_getfilehead(gfile_in,iret=iret,
     &     dimz=dimz,idvc=idvc)
!
      gen_coord_hybrid=.false.
      call nemsio_getheadvar(gfile_in,'GEN_COOR', 
     &     gen_coord_hybrid, iret=iret)
      hybrid=.false.
      call nemsio_getheadvar(gfile_in,'HYBRID',hybrid,iret=iret)  
!      print *,' in phys, read vcoord,hybrid=',hybrid,'dimz=',
!     &  dimz,'idvc=',idvc
      if (iret /= 0) then
        if (idvc==2) hybrid=.true.
      endif
!
      allocate(vcoord(dimz+1,3,2))
      call nemsio_getfilehead(gfile_in,iret=iret,vcoord=vcoord)
!
      call nemsio_close(gfile_in,iret)
! 
      if (gen_coord_hybrid) then
!   ak bk ck in file have the same order as model
        do k=1,dimz+1
          ak5(k) = vcoord(k,1,1)/1000.
          bk5(k) = vcoord(k,2,1)
          ck5(k) = vcoord(k,3,1)/1000.
        enddo
      else if (hybrid .and. idvc .eq. 2) then
        ck5=0.
        do k=1,dimz+1
          ak5(k) = vcoord(dimz+2-k,1,1)/1000.
          bk5(k) = vcoord(dimz+2-k,2,1)
        enddo
      endif
!
      deallocate(vcoord)

      end subroutine read_vcoord
!

      SUBROUTINE read_sfc_nemsio(sfc_fld,NEEDORO,nread,
     &                    cfile,global_lats_r,lonsperlar)
!
!***********************************************************************
!
! !REVISION HISTORY:
!
!  2011/02/26  Sarah Lu, Modified to read cold-start (from chgres) and warm
!              start (from replay) ICs
!
!      use sfcio_module, ONLY: sfcio_head, sfcio_data, sfcio_realfill,
!     &                        sfcio_srohdc, sfcio_axdata
      use resol_def,    ONLY: latr, latr2, lonr, lsoil
      use layout1,      ONLY: me, nodes, lats_node_r, lats_node_r_max
     &,                       nodes_comp
      use gfs_physics_sfc_flx_mod, ONLY: Sfc_Var_Data
      use namelist_soilveg ,       only: salp_data, snupx
      use physcons,     only : tgice => con_tice
      USE machine,      ONLY: kind_io4, kind_io8
      use nemsio_module
      implicit none
!
      TYPE(Sfc_Var_Data)        :: sfc_fld
      integer              global_lats_r(latr)
      integer                 lonsperlar(latr)

      integer jump
      integer needoro

      real(kind=kind_io4) buff1(lonr*latr),buff2(lonr,latr,LSOIL)
      real(kind=kind_io8) buffo(lonr,lats_node_r_max)
      real(kind=kind_io8) buff3(lonr,lats_node_r_max)
      integer nread,i,j,k,ij,idate7(7),lonsfc,latsfc,lplsfc(latr2)
      character*(*) cfile
      integer kmsk(lonr,latr),kmskcv(lonr,latr)
      CHARACTER*8 labfix(4)
      real t1,t2,rsnow
      real(4) fhour
      integer nfhour,nfminute,nfsecondn,nfsecondd
      type(nemsio_gfile) gfile_in
      integer iret, vegtyp,lonb4,latb4,nsoil4,ivs4
      integer size1, size2, size3, iop
      real(8) timef,stime,etime
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
      t1=timef()
!     iop = nodes_comp - 1
      iop = 0

      if(me == iop) print *,' nread=',nread,' cfile=',cfile
      call nemsio_init()
!
      stime=timef()
      call nemsio_open(gfile_in,trim(cfile),'read',iret=iret)
       etime=timef()
!        print *,'in read_sfc,nemsio open time=',etime-stime
!
      IF (me == iop) THEN

        call nemsio_getheadvar(gfile_in,'dimx',lonb4,iret=iret)
        call nemsio_getheadvar(gfile_in,'dimy',latb4,iret=iret)
        call nemsio_getheadvar(gfile_in,'nsoil',nsoil4,iret=iret)
        call nemsio_getheadvar(gfile_in,'ivs',ivs4,iret=iret)
        if(iret /= 0) then
          call nemsio_getheadvar(gfile_in,'ivssfc',ivs4,iret=iret)
        endif
        call nemsio_getfilehead(gfile_in,idate=idate7,nfhour=nfhour,    
     &       nfminute=nfminute,nfsecondn=nfsecondn,nfsecondd=nfsecondd,
     &       iret=iret)   
        if(iret/=0) print *,' after sfcio_srohdc,iret=',iret

        FHOUR     = real(nfhour,8)+real(nfminute,8)/60.+        
     &              real(nfsecondn,8)/(real(nfsecondd,8)*3600.)

!        PRINT 99,nread,head%fhour,head%idate,
!     &           head%lonb,head%latb,head%lsoil,head%ivs,iret
        write(*,99) nread,fhour,idate7(1:4),
     &           lonb4,latb4,nsoil4,ivs4,iret

 99     FORMAT(1H ,'in read_sfc_nemsio, nread=',i3,2x,'HOUR=',f8.2,3x,
     &   'IDATE=',4(1X,I4),4x,'lonsfc,latsfc,lsoil,ivssfc,iret=',5i8)

        if(iret.ne.0) goto 6000
!        if(lonb4.ne.lonr) goto 6000
!        if(latb4.ne.latr) goto 6000
        if(nsoil4.ne.lsoil) goto 6000

      ENDIF

      kmsk = 0
!
       stime=timef()
      if(me == iop) call nemsio_readrecv(gfile_in,'tmp','sfc',1,buff1,
     &                                   iret=iret)
      call split2d_phys(buff1, buffo,global_lats_r,iop)
      CALL interpred_phys(1,kmsk,buffo,sfc_fld%TSEA,
     &    global_lats_r,lonsperlar)

      DO K=1, LSOIL

        if(me == iop) call nemsio_readrecv(gfile_in,'smc','soil layer',
     &                                      k,buff1,iret=iret)
        call split2d_phys(buff1, buffo,global_lats_r,iop)
        CALL interpred_phys(1,kmsk,buffo,buff3,global_lats_r,lonsperlar)
        sfc_fld%SMC(k,:,:) = buff3(:,:)
      ENDDO

      if(me == iop) call nemsio_readrecv(gfile_in,'weasd','sfc',1,buff1,
     &                                   iret=iret)
      call split2d_phys(buff1, buffo,global_lats_r,iop)
      CALL interpred_phys(1,kmsk,buffo,sfc_fld%weasd,
     &               global_lats_r,lonsperlar)

      DO K = 1, LSOIL
        if(me == iop) call nemsio_readrecv(gfile_in,'stc','soil layer',
     &                                      k,buff1,iret=iret)
        call split2d_phys(buff1, buffo,global_lats_r,iop)
        CALL interpred_phys(1,kmsk,buffo,buff3,global_lats_r,lonsperlar)
        sfc_fld%STC(k,:,:) = buff3(:,:)
      ENDDO

      if(me == iop) call nemsio_readrecv(gfile_in,'tg3','sfc',1,buff1,
     &                                   iret=iret)
      call split2d_phys(buff1, buffo,global_lats_r,iop)
      CALL interpred_phys(1,kmsk,buffo,sfc_fld%TG3,
     &    global_lats_r,lonsperlar)

      if(me == iop) call nemsio_readrecv(gfile_in,'sfcr','sfc',1,buff1,
     &                                   iret=iret)
      call split2d_phys(buff1, buffo,global_lats_r,iop)
      CALL interpred_phys(1,kmsk,buffo,sfc_fld%ZORL,
     &    global_lats_r,lonsperlar)

      sfc_fld%cv  = 0
      sfc_fld%cvb = 0
      sfc_fld%cvt = 0

      if(me == iop) call nemsio_readrecv(gfile_in,'alvsf','sfc',1,buff1,
     &                                 iret=iret)
!      if(me==0) buff1=data%alvsf
      call split2d_phys(buff1, buffo,global_lats_r,iop)
      CALL interpred_phys(1,kmsk,buffo,sfc_fld%ALVSF,
     &               global_lats_r,lonsperlar)
      if(me == iop) call nemsio_readrecv(gfile_in,'alvwf','sfc',1,buff1,
     &                                   iret=iret)
!      if(me==0) buff1=data%alvwf
      call split2d_phys(buff1, buffo,global_lats_r,iop)
      CALL interpred_phys(1,kmsk,buffo,sfc_fld%ALVWF,
     &               global_lats_r,lonsperlar)
      if(me == iop) call nemsio_readrecv(gfile_in,'alnsf','sfc',1,buff1,
     &                                   iret=iret)
!      if(me==0) buff1=data%alnsf
      call split2d_phys(buff1, buffo,global_lats_r,iop)
      CALL interpred_phys(1,kmsk,buffo,sfc_fld%ALNSF,
     &               global_lats_r,lonsperlar)
      if(me == iop) call nemsio_readrecv(gfile_in,'alnwf','sfc',1,buff1,
     &                                   iret=iret)
!      if(me==0) buff1=data%alnwf
      call split2d_phys(buff1, buffo,global_lats_r,iop)
      CALL interpred_phys(1,kmsk,buffo,sfc_fld%ALNWF,
     &               global_lats_r,lonsperlar)

!     The mask cannot be interpolated
      if(me == iop) call nemsio_readrecv(gfile_in,'land','sfc',1,buff1,
     &                                   iret=iret)
!      if(me==0) buff1=data%slmsk
      call split2d_phys(buff1, buffo,global_lats_r,iop)
      CALL interpred_phys(1,kmsk,buffo,sfc_fld%SLMSK,
     &               global_lats_r,lonsperlar)

      if(me == iop) call nemsio_readrecv(gfile_in,'veg','sfc',1,buff1,
     &                                   iret=iret)
!      if(me==0) buff1=data%vfrac
      call split2d_phys(buff1, buffo,global_lats_r,iop)
      CALL interpred_phys(1,kmsk,buffo,sfc_fld%VFRAC,
     &               global_lats_r,lonsperlar)

      if(me == iop) call nemsio_readrecv(gfile_in,'cnwat','sfc',1,buff1,
     &                                   iret=iret)
!      if(me==0) buff1=data%canopy
      call split2d_phys(buff1, buffo,global_lats_r,iop)
      CALL interpred_phys(1,kmsk,buffo,sfc_fld%CANOPY,
     &               global_lats_r,lonsperlar)

      if(me == iop) call nemsio_readrecv(gfile_in,'f10m','sfc',1,buff1,
     &                                   iret=iret)
!      if(me==0) buff1=data%f10m
      call split2d_phys(buff1, buffo,global_lats_r,iop)
      CALL interpred_phys(1,kmsk,buffo,sfc_fld%F10M,
     &    global_lats_r,lonsperlar)

      if(me == iop) call nemsio_readrecv(gfile_in,'vtype','sfc',1,buff1,
     &                                   iret=iret)
!      if(me==0) buff1=data%vtype
      call split2d_phys(buff1, buffo,global_lats_r,iop)
      CALL interpred_phys(1,kmsk,buffo,sfc_fld%VTYPE,
     &               global_lats_r,lonsperlar)

      if(me == iop) call nemsio_readrecv(gfile_in,'sotyp','sfc',1,buff1,
     &                                   iret=iret)
!      if(me==0) buff1=data%stype
      call split2d_phys(buff1, buffo,global_lats_r,iop)
      CALL interpred_phys(1,kmsk,buffo,sfc_fld%STYPE,
     &               global_lats_r,lonsperlar)

      if(me == iop) call nemsio_readrecv(gfile_in,'facsf','sfc',1,buff1,
     &                                   iret=iret)
!      if(me==0) buff1=data%facsf
      call split2d_phys(buff1, buffo,global_lats_r,iop)
      CALL interpred_phys(1,kmsk,buffo,sfc_fld%FACSF,
     &               global_lats_r,lonsperlar)
      if(me == iop) call nemsio_readrecv(gfile_in,'facwf','sfc',1,buff1,
     &                                   iret=iret)
!      if(me==0) buff1=data%facwf
      call split2d_phys(buff1, buffo,global_lats_r,iop)
      CALL interpred_phys(1,kmsk,buffo,sfc_fld%FACWF,
     &               global_lats_r,lonsperlar)

!
      if(me == iop) call nemsio_readrecv(gfile_in,'fricv','sfc',1,buff1,
     &                                   iret=iret)
!        if(me==0) buff1=data%uustar
         call split2d_phys(buff1, buffo,global_lats_r,iop)
         CALL interpred_phys(1,kmsk,buffo,sfc_fld%UUSTAR,
     &               global_lats_r,lonsperlar)

      if(me == iop) call nemsio_readrecv(gfile_in,'ffmm','sfc',1,buff1,
     &                                   iret=iret)
!        if(me==0) buff1=data%ffmm
         call split2d_phys(buff1, buffo,global_lats_r,iop)
         CALL interpred_phys(1,kmsk,buffo,sfc_fld%FFMM,
     &                  global_lats_r,lonsperlar)

      if(me == iop) call nemsio_readrecv(gfile_in,'ffhh','sfc',1,buff1,
     &                                   iret=iret)
!        if(me==0) buff1=data%ffhh
         call split2d_phys(buff1, buffo,global_lats_r,iop)
         CALL interpred_phys(1,kmsk,buffo,sfc_fld%FFHH,
     &                  global_lats_r,lonsperlar)

!    Sea-ice (hice/fice) was added to the surface files.

      if(me == iop) call nemsio_readrecv(gfile_in,'icetk','sfc',1,buff1,
     &                                   iret=iret)
!         if(me==0) buff1=data%hice
         call split2d_phys(buff1, buffo,global_lats_r,iop)
         CALL interpred_phys(1,kmsk,buffo,sfc_fld%HICE,
     &                  global_lats_r,lonsperlar)

      if(me == iop) call nemsio_readrecv(gfile_in,'icec','sfc',1,buff1,
     &                                   iret=iret)
!         if(me==0) buff1=data%fice
         call split2d_phys(buff1, buffo,global_lats_r,iop)
         CALL interpred_phys(1,kmsk,buffo,sfc_fld%FICE,
     &                  global_lats_r,lonsperlar)

      if(me == iop) call nemsio_readrecv(gfile_in,'tisfc','sfc',1,buff1,
     &                                   iret=iret)
!         if(me==0) buff1=data%tisfc
         call split2d_phys(buff1, buffo,global_lats_r,iop)
         CALL interpred_phys(1,kmsk,buffo,sfc_fld%TISFC,
     &                  global_lats_r,lonsperlar)
         if (lats_node_r > 0 )  then
           if (sfc_fld%tisfc(1,1) < 0.0)  then
             DO j=1,lats_node_r
               DO i=1,LONR
                  sfc_fld%TISFC(i,j)= sfc_fld%TSEA(i,j)
                  IF(sfc_fld%SLMSK(i,j) >=  2. .AND.
!    &               sfc_fld%FICE(i,j)  >= 0.5) THEN
     &               sfc_fld%FICE(i,j)  >= 0.15) THEN
                     sfc_fld%TISFC(i,j) = (sfc_fld%TSEA(i,j)
     &              -tgice*(1.-sfc_fld%FICE(i,j))) / sfc_fld%FICE(i,j)
                   sfc_fld%TISFC(i,j)=MIN(sfc_fld%TISFC(i,j),tgice)
                   ENDIF
               ENDDO
             ENDDO
           endif
         endif


!*     surface files for GFS/Noah contain 8 additional records:
!*     tprcp, srflag, snwdph, slc, shdmin, shdmax, slope, snoalb

      if(me == iop) call nemsio_readrecv(gfile_in,'tprcp','sfc',1,buff1,
     &                                   iret=iret)
!         if(me==0) buff1=data%tprcp
         call split2d_phys(buff1, buffo,global_lats_r,iop)
         CALL interpred_phys(1,kmsk,buffo,sfc_fld%TPRCP,
     &                  global_lats_r,lonsperlar)

!* srflag
      if(me == iop) call nemsio_readrecv(gfile_in,'crain','sfc',1,buff1,
     &                                   iret=iret)
!         if(me==0) buff1=data%srflag
         call split2d_phys(buff1, buffo,global_lats_r,iop)
         CALL interpred_phys(1,kmsk,buffo,sfc_fld%SRFLAG,
     &                  global_lats_r,lonsperlar)

!* snwdph
      if(me == iop) call nemsio_readrecv(gfile_in,'snod','sfc',1,buff1,
     &                                   iret=iret)
!         if(me==0) buff1=data%snwdph
         call split2d_phys(buff1, buffo,global_lats_r,iop)
         CALL interpred_phys(1,kmsk,buffo,sfc_fld%SNWDPH,
     &                  global_lats_r,lonsperlar)

!* slc
         DO K=1, LSOIL
!         if(me==0) buff1=data%slc(:,:,k)
      if(me == iop) call nemsio_readrecv(gfile_in,'slc','soil layer',k,
     &                                   buff1,iret=iret)
         call split2d_phys(buff1, buffo,global_lats_r,iop)
         CALL interpred_phys(1,kmsk,buffo,buff3,
     &    global_lats_r,lonsperlar)
         sfc_fld%SLC(k,:,:) = buff3(:,:)
         ENDDO

!* shdmin
      if(me == iop) call nemsio_readrecv(gfile_in,'shdmin','sfc',1,buff1
     &,                                   iret=iret)
!         if(me==0) buff1=data%shdmin
         call split2d_phys(buff1, buffo,global_lats_r,iop)
         CALL interpred_phys(1,kmsk,buffo,sfc_fld%SHDMIN,
     &                  global_lats_r,lonsperlar)

!* shdmax
      if(me == iop) call nemsio_readrecv(gfile_in,'shdmax','sfc',1,buff1
     &,                                  iret=iret)
!         if(me==0) buff1=data%shdmax
         call split2d_phys(buff1, buffo,global_lats_r,iop)
         CALL interpred_phys(1,kmsk,buffo,sfc_fld%SHDMAX,
     &                  global_lats_r,lonsperlar)

!* slope
      if(me == iop) call nemsio_readrecv(gfile_in,'sltyp','sfc',1,buff1,
     &                                   iret=iret)
!         if(me==0) buff1=data%slope
         call split2d_phys(buff1, buffo,global_lats_r,iop)
         CALL interpred_phys(1,kmsk,buffo,sfc_fld%SLOPE,
     &                  global_lats_r,lonsperlar)

!* snoalb
      if(me == iop) call nemsio_readrecv(gfile_in,'salbd','sfc',1,buff1,
     &                                   iret=iret)
!         if(me==0) buff1=data%snoalb
         call split2d_phys(buff1, buffo,global_lats_r,iop)
         CALL interpred_phys(1,kmsk,buffo,sfc_fld%SNOALB,
     &                  global_lats_r,lonsperlar)
!     print *,' snoalb=',sfc_fld%snoalb(1,:)
!lu [+67L]: the addition of 8 Noah records ends here .........................

       if(needoro == 1) then
         if (me == iop) then
           call nemsio_readrecv(gfile_in,'orog','sfc',1,buff1,iret=iret)
!           buff1=data%orog
           needoro = 1
           if(all(buff1.ne.-9999.)) needoro=0
!           print *,'read sfc orography'
         endif
         call split2d_phys(buff1, buffo,global_lats_r,iop)
         CALL interpred_phys(1,kmsk,buffo,sfc_fld%ORO,
     &                  global_lats_r,lonsperlar)
         call skip(needoro)
       endif
       etime=timef()
!      print *,'in read_sfc,nemsio read time=',etime-stime
!
!Wei initialize snow fraction(weasd is in mm)
      DO j=1,lats_node_r
        DO i=1,LONR
          sfc_fld%SNCOVR(i,j) = 0.0
          if (sfc_fld%slmsk(i,j) > 0.001 .AND. 
     &        ABS(sfc_fld%VTYPE(i,j)) >= 0.5 ) then
            vegtyp = sfc_fld%VTYPE(i,j)
            RSNOW  = 0.001*sfc_fld%weasd(i,j)/SNUPX(vegtyp)
            IF (0.001*sfc_fld%weasd(i,j) < SNUPX(vegtyp)) THEN
              sfc_fld%SNCOVR(i,j) = 1.0 - ( EXP(-SALP_DATA*RSNOW)
     &                                    - RSNOW*EXP(-SALP_DATA))
            ELSE
              sfc_fld%SNCOVR(i,j) = 1.0
            ENDIF
!           if (i == 1)
!    &       print*,SNUPX(vegtyp),SALP_DATA,sfc_fld%SNCOVR(i,j),
!    &       '************debug',sfc_fld%weasd(i,j),vegtyp,' j=',j
!    &,      ' snoalb1=',sfc_fld%snoalb(i,j)
!
          endif
        ENDDO
       ENDDO
!
       IF (me == iop) then
!         call sfcio_axdata(data,iret)
         t2 = timef()
!         print *,'FIXIO TIME ',t2-t1,t1,t2
       endif
!
      call nemsio_close(gfile_in,iret=iret)
!
      call nemsio_finalize()
!
      RETURN
 6000 PRINT *, ' error in input in routine read_sfc_nemsio'
      STOP
      END
!
      SUBROUTINE read_nst(nst_fld, nread, cfile,
     &                   global_lats_r, lonsperlar)
!
!***********************************************************************
!
      use namelist_physics_def
      USE machine,        ONLY: kind_ior, kind_io8, kind_rad
      use nstio_module
      use resol_def
      use layout1
      use mpi_def
      use gfs_physics_nst_var_mod
      implicit none
!
      TYPE(Nst_Var_Data)       :: nst_fld
      integer              global_lats_r(latr)
      integer                 lonsperlar(latr)

!     real (kind=kind_io8) slmsk(lonr,lats_node_r),

      real(kind=kind_io4) buff1(lonr,latr)
      real(kind=kind_io8) buffo(lonr,lats_node_r_max)
      integer nread,i,j,k,ij,idate(4),lonnst,latnst,lplnst(latr2)
      character*(*) cfile
      integer kmsk(lonr,latr)
      CHARACTER*8 labfix(4)
      real t1,t2,timef
      type(nstio_head) head
      type(nstio_data) data
      integer iret,iop
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!     iop = nodes_comp - 1
      iop = 0
      t1=timef()

      print *,'read nst filem nread=',nread,'cfile=',cfile
      IF (me == 0) then
        call nstio_srohdc(nread,cfile,head,data,iret)

        PRINT 99,nread,head%fhour,head%idate,
     &     head%lonb,head%latb,head%lsea,head%ivn,iret,lats_node_r
99      FORMAT(1H ,'in fixio nread=',i3,2x,'HOUR=',f8.2,3x,'IDATE=',
     &  4(1X,I4),4x,'lonnst,latnst,lsea,ivsnst,iret=',6i8)

        if(iret.ne.0) goto 5000
        if(head%lonb.ne.lonr) goto 5000
        if(head%latb.ne.latr) goto 5000
        if(head%lsea.ne.lsea) goto 5000

      ENDIF

      kmsk=0
!
!     Assign ocnf(lonr,lats_node_r,nf_ocn)
!
      IF (me == 0)  buff1=data%xt
      call split2d_phys(buff1, buffo,global_lats_r,iop)
      CALL interpred_phys(1,kmsk,buffo,nst_fld%xt,
     &               global_lats_r,lonsperlar)

      IF (me == 0) buff1=data%xs
      call split2d_phys(buff1, buffo,global_lats_r,iop)
      CALL interpred_phys(1,kmsk,buffo,nst_fld%xs,
     &               global_lats_r,lonsperlar)

      IF (me == 0) buff1=data%xu
      call split2d_phys(buff1, buffo,global_lats_r,iop)
      CALL interpred_phys(1,kmsk,buffo,nst_fld%xu,
     &               global_lats_r,lonsperlar)

      IF (me == 0) buff1=data%xv
      call split2d_phys(buff1, buffo,global_lats_r,iop)
      CALL interpred_phys(1,kmsk,buffo,nst_fld%xv,
     &               global_lats_r,lonsperlar)

      IF (me == 0) buff1=data%xz
      call split2d_phys(buff1, buffo,global_lats_r,iop)
      CALL interpred_phys(1,kmsk,buffo,nst_fld%xz,
     &               global_lats_r,lonsperlar)

      IF (me == 0) buff1=data%zm
      call split2d_phys(buff1, buffo,global_lats_r,iop)
      CALL interpred_phys(1,kmsk,buffo,nst_fld%zm,
     &               global_lats_r,lonsperlar)

      IF (me == 0) buff1=data%xtts
      call split2d_phys(buff1, buffo,global_lats_r,iop)
      CALL interpred_phys(1,kmsk,buffo,nst_fld%xtts,
     &               global_lats_r,lonsperlar)

      IF (me == 0) buff1=data%xzts
      call split2d_phys(buff1, buffo,global_lats_r)
      CALL interpred_phys(1,kmsk,buffo,nst_fld%xzts,
     &               global_lats_r,lonsperlar)

      IF (me == 0) buff1=data%dt_cool
      call split2d_phys(buff1, buffo,global_lats_r,iop)
      CALL interpred_phys(1,kmsk,buffo,nst_fld%dt_cool,
     &               global_lats_r,lonsperlar)

      IF (me == 0) buff1=data%z_c
      call split2d_phys(buff1, buffo,global_lats_r,iop)
      CALL interpred_phys(1,kmsk,buffo,nst_fld%z_c,
     &               global_lats_r,lonsperlar)

      IF (me == 0) buff1=data%c_0
      call split2d_phys(buff1, buffo,global_lats_r,iop)
      CALL interpred_phys(1,kmsk,buffo,nst_fld%c_0,
     &               global_lats_r,lonsperlar)
      IF (me == 0) buff1=data%c_d
      call split2d_phys(buff1, buffo,global_lats_r,iop)
      CALL interpred_phys(1,kmsk,buffo,nst_fld%c_d,
     &               global_lats_r,lonsperlar)

      IF (me == 0) buff1=data%w_0
      call split2d_phys(buff1, buffo,global_lats_r,iop)
      CALL interpred_phys(1,kmsk,buffo,nst_fld%w_0,
     &               global_lats_r,lonsperlar)

      IF (me == 0) buff1=data%w_d
      call split2d_phys(buff1, buffo,global_lats_r,iop)
      CALL interpred_phys(1,kmsk,buffo,nst_fld%w_d,
     &               global_lats_r,lonsperlar)

      IF (me == 0) buff1=data%d_conv
      call split2d_phys(buff1, buffo,global_lats_r,iop)
      CALL interpred_phys(1,kmsk,buffo,nst_fld%d_conv,
     &               global_lats_r,lonsperlar)

      IF (me == 0) buff1=data%ifd
      call split2d_phys(buff1, buffo,global_lats_r,iop)
      CALL interpred_phys(1,kmsk,buffo,nst_fld%ifd,
     &               global_lats_r,lonsperlar)

      IF (me == 0) buff1=data%tref
      call split2d_phys(buff1, buffo,global_lats_r,iop)
      CALL interpred_phys(1,kmsk,buffo,nst_fld%tref,
     &               global_lats_r,lonsperlar)

      IF (me == 0) buff1=data%Qrain
      call split2d_phys(buff1, buffo,global_lats_r,iop)
      CALL interpred_phys(1,kmsk,buffo,nst_fld%Qrain,
     &               global_lats_r,lonsperlar)

!     IF (icolor.eq.2.and.me.eq.nodes-1) then
      IF (me == 0) then
         call nstio_axdata(data,iret)
         t2=timef()
         print *,'FIXIO for NST TIME ',t2-t1,t1,t2
      endif
!
      RETURN
 5000 PRINT *, ' ERROR IN INPUT IN read_nst'
      STOP
      END
!
      SUBROUTINE read_nst_nemsio(nst_fld,nread,cfile,
     &                    global_lats_r,lonsperlar)
!
!***********************************************************************
!
! !REVISION HISTORY:
!
!  2015/08/26  Xu Li, created to read nemsio NSST files
!
      use namelist_physics_def, ONLY: lsea
      use resol_def,    ONLY: latr, latr2, lonr
      use layout1,      ONLY: me, nodes, lats_node_r, lats_node_r_max
     &,                       nodes_comp
      use gfs_physics_nst_var_mod, ONLY: Nst_Var_Data
      USE machine,      ONLY: kind_io4, kind_io8
      use nemsio_module
      implicit none
!
      TYPE(Nst_Var_Data)        :: nst_fld
      integer              global_lats_r(latr)
      integer                 lonsperlar(latr)

      integer jump
      integer needoro

      real(kind=kind_io4) buff1(lonr*latr)
      real(kind=kind_io8) buffo(lonr,lats_node_r_max)
      integer nread,i,j,k,ij,idate7(7),lonsfc,latsfc,lplsfc(latr2)
      character*(*) cfile
      integer kmsk(lonr,latr)
      CHARACTER*8 labfix(4)
      real t1,t2
      real(4) fhour
      integer nfhour,nfminute,nfsecondn,nfsecondd
      type(nemsio_gfile) gfile_in
      integer iret,lonb4,latb4,nsea4,ivn4
      integer size1, size2, size3, iop
      real(8) timef,stime,etime
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
      t1=timef()
!     iop = nodes_comp - 1
      iop = 0

      if(me == iop) print *,' nread=',nread,' cfile=',cfile
      call nemsio_init()
!
      stime=timef()
      call nemsio_open(gfile_in,trim(cfile),'read',iret=iret)
      etime=timef()
         print *,'in read_nst_nemsio open time=',etime-stime
!
      IF (me == iop) THEN

        call nemsio_getheadvar(gfile_in,'dimx',lonb4,iret=iret)
        call nemsio_getheadvar(gfile_in,'dimy',latb4,iret=iret)
        call nemsio_getheadvar(gfile_in,'ivsnst',ivn4,iret=iret)

        call nemsio_getfilehead(gfile_in,idate=idate7,nfhour=nfhour,
     &       nfminute=nfminute,nfsecondn=nfsecondn,nfsecondd=nfsecondd,
     &       iret=iret)
        if(iret/=0) print *,' after nstio_srohdc,iret=',iret

        FHOUR     = real(nfhour,8)+real(nfminute,8)/60.+
     &              real(nfsecondn,8)/(real(nfsecondd,8)*3600.)

        PRINT 99,nread,fhour,idate7(1:4),
     &           lonb4,latb4,nsea4,ivn4,iret

99      FORMAT(1H ,'in read_nst_nemsio, nread=',i3,2x,'HOUR=',f8.2,3x,
     &   'IDATE=',4(1X,I4),4x,'lonsfc,latsfc,lsea,ivsnst,iret=',5i8)

        if(iret.ne.0) goto 7000
        if(lonb4.ne.lonr) goto 7000
        if(latb4.ne.latr) goto 7000
        if(nsea4.ne.lsea) goto 7000

      ENDIF

      kmsk = 0
!
      stime=timef()

      if(me == iop) call nemsio_readrecv(gfile_in,'slmsk','sfc',1,buff1,
     &                                   iret=iret)
      call split2d_phys(buff1, buffo,global_lats_r,iop)
      CALL interpred_phys(1,kmsk,buffo,nst_fld%slmsk,
     &    global_lats_r,lonsperlar)

      if(me == iop) call nemsio_readrecv(gfile_in,'xt','sfc',1,buff1,
     &                                   iret=iret)
      call split2d_phys(buff1, buffo,global_lats_r,iop)
      CALL interpred_phys(1,kmsk,buffo,nst_fld%xt,
     &               global_lats_r,lonsperlar)

      if(me == iop) call nemsio_readrecv(gfile_in,'xs','sfc',1,buff1,
     &                                   iret=iret)
      call split2d_phys(buff1, buffo,global_lats_r,iop)
      CALL interpred_phys(1,kmsk,buffo,nst_fld%xs,
     &               global_lats_r,lonsperlar)

      if(me == iop) call nemsio_readrecv(gfile_in,'xu','sfc',1,buff1,
     &                                   iret=iret)
      call split2d_phys(buff1, buffo,global_lats_r,iop)
      CALL interpred_phys(1,kmsk,buffo,nst_fld%xu,
     &               global_lats_r,lonsperlar)

      if(me == iop) call nemsio_readrecv(gfile_in,'xv','sfc',1,buff1,
     &                                   iret=iret)
      call split2d_phys(buff1, buffo,global_lats_r,iop)
      CALL interpred_phys(1,kmsk,buffo,nst_fld%xv,
     &               global_lats_r,lonsperlar)     
      if(me == iop) call nemsio_readrecv(gfile_in,'xz','sfc',1,buff1,

     &                                   iret=iret)
      call split2d_phys(buff1, buffo,global_lats_r,iop)
      CALL interpred_phys(1,kmsk,buffo,nst_fld%xz,
     &               global_lats_r,lonsperlar)

      if(me == iop) call nemsio_readrecv(gfile_in,'zm','sfc',1,buff1,
     &                                   iret=iret)
      call split2d_phys(buff1, buffo,global_lats_r,iop)
      CALL interpred_phys(1,kmsk,buffo,nst_fld%zm,
     &               global_lats_r,lonsperlar)

      if(me == iop) call nemsio_readrecv(gfile_in,'xtts','sfc',1,buff1,
     &                                   iret=iret)
      call split2d_phys(buff1, buffo,global_lats_r,iop)
      CALL interpred_phys(1,kmsk,buffo,nst_fld%xtts,
     &               global_lats_r,lonsperlar)

      if(me == iop) call nemsio_readrecv(gfile_in,'xzts','sfc',1,buff1,
     &                                   iret=iret)
      call split2d_phys(buff1, buffo,global_lats_r,iop)
      CALL interpred_phys(1,kmsk,buffo,nst_fld%xzts,
     &               global_lats_r,lonsperlar)

      if(me == iop) call nemsio_readrecv(gfile_in,'dtcool','sfc',1,
     &                                   buff1,iret=iret)
      call split2d_phys(buff1, buffo,global_lats_r,iop)
      CALL interpred_phys(1,kmsk,buffo,nst_fld%dt_cool,
     &               global_lats_r,lonsperlar)

      if(me == iop) call nemsio_readrecv(gfile_in,'zc','sfc',1,buff1,
     &                                   iret=iret)
      call split2d_phys(buff1, buffo,global_lats_r,iop)
      CALL interpred_phys(1,kmsk,buffo,nst_fld%z_c,
     &               global_lats_r,lonsperlar)

      if(me == iop) call nemsio_readrecv(gfile_in,'c0','sfc',1,buff1,
     &                                   iret=iret)
      call split2d_phys(buff1, buffo,global_lats_r,iop)
      CALL interpred_phys(1,kmsk,buffo,nst_fld%c_0,
     &               global_lats_r,lonsperlar)     
      if(me == iop) call nemsio_readrecv(gfile_in,'cd','sfc',1,buff1,

     &                                   iret=iret)
      call split2d_phys(buff1, buffo,global_lats_r,iop)
      CALL interpred_phys(1,kmsk,buffo,nst_fld%c_d,
     &               global_lats_r,lonsperlar)

      if(me == iop) call nemsio_readrecv(gfile_in,'w0','sfc',1,buff1,
     &                                   iret=iret)
      call split2d_phys(buff1, buffo,global_lats_r,iop)
      CALL interpred_phys(1,kmsk,buffo,nst_fld%w_0,
     &               global_lats_r,lonsperlar)

      if(me == iop) call nemsio_readrecv(gfile_in,'wd','sfc',1,buff1,
     &                                   iret=iret)
      call split2d_phys(buff1, buffo,global_lats_r,iop)
      CALL interpred_phys(1,kmsk,buffo,nst_fld%w_d,
     &               global_lats_r,lonsperlar)

      if(me == iop) call nemsio_readrecv(gfile_in,'dconv','sfc',1,buff1,
     &                                   iret=iret)
      call split2d_phys(buff1, buffo,global_lats_r,iop)
      CALL interpred_phys(1,kmsk,buffo,nst_fld%d_conv,
     &               global_lats_r,lonsperlar)

      if(me == iop) call nemsio_readrecv(gfile_in,'ifd','sfc',1,buff1,
     &                                   iret=iret)
      call split2d_phys(buff1, buffo,global_lats_r,iop)
      CALL interpred_phys(1,kmsk,buffo,nst_fld%ifd,
     &               global_lats_r,lonsperlar)

      if(me == iop) call nemsio_readrecv(gfile_in,'tref','sfc',1,buff1,
     &                                   iret=iret)
      call split2d_phys(buff1, buffo,global_lats_r,iop)
      CALL interpred_phys(1,kmsk,buffo,nst_fld%tref,
     &               global_lats_r,lonsperlar)     
      if(me == iop) call nemsio_readrecv(gfile_in,'qrain','sfc',1,buff1,
     &                                   iret=iret)
      call split2d_phys(buff1, buffo,global_lats_r,iop)
      CALL interpred_phys(1,kmsk,buffo,nst_fld%qrain,
     &               global_lats_r,lonsperlar)


       etime=timef()
!      print *,'in read_nst_nemsio read time=',etime-stime
!
       IF (me == iop) then
!         call nstio_axdata(data,iret)
          t2 = timef()
!         print *,'read_nst_nemsio TIME ',t2-t1,t1,t2
       endif
!
      call nemsio_close(gfile_in,iret=iret)
!
      call nemsio_finalize()
!
      RETURN
 7000 PRINT *, ' error in input in routine read_nst_nemsio'
      END           ! read_nst_nemsio

      SUBROUTINE set_nst(tsea, nst_fld)
c
c***********************************************************************
c
      use namelist_physics_def
      USE machine,     ONLY: kind_io8
      use resol_def
      use layout1
      use gfs_physics_nst_var_mod
      use module_nst_parameters, only: z_w_max
      use mpi_def
      implicit none
c
      TYPE(Nst_Var_Data)       :: nst_fld
      real (kind=kind_io8) tsea(lonr,lats_node_r)

      integer i,j,k
      real t1,t2,timef

c
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c
      t1=timef()

!      print *,'in set_nst start'
      nst_fld%xt      = 0.0
      nst_fld%xs      = 0.0
      nst_fld%xu      = 0.0
      nst_fld%xv      = 0.0
      nst_fld%xz      = z_w_max
      nst_fld%zm      = 0.0
      nst_fld%xtts    = 0.0
      nst_fld%xzts    = 0.0
      nst_fld%dt_cool = 0.0
      nst_fld%z_c     = 0.0
      nst_fld%c_0     = 0.0
      nst_fld%c_d     = 0.0
      nst_fld%w_0     = 0.0
      nst_fld%w_d     = 0.0
      nst_fld%d_conv  = 0.0
      nst_fld%ifd     = 0.0
      nst_fld%Tref(:,1:lats_node_r)= tsea(:,1:lats_node_r)
      nst_fld%Qrain   = 0.0
!
      t2=timef()
!      print *,'FIXIO for set_nst TIME ',t2-t1,t1,t2
!
      RETURN
      END
!
!***********************************************************************
!
      SUBROUTINE nst_reset_nonwater(tsea,nst_fld)
c
c***********************************************************************
c
      use resol_def
      USE machine,     ONLY: kind_io8
      use layout1
      use gfs_physics_nst_var_mod
      use module_nst_parameters, only: z_w_max
      use mpi_def
      implicit none
c
      TYPE(Nst_Var_Data)       :: nst_fld
      real (kind=kind_io8) tsea(lonr,lats_node_r)

      integer i,j
      real t1,t2,timef
c
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c
      t1=timef()

      do j = 1, lats_node_r
        do i = 1, lonr
          if ( nst_fld%slmsk(i,j) /= 0.0 ) then
            nst_fld%xt(i,j)      = 0.0
            nst_fld%xs(i,j)      = 0.0
            nst_fld%xu(i,j)      = 0.0
            nst_fld%xv(i,j)      = 0.0
            nst_fld%xz(i,j)      = z_w_max
            nst_fld%zm(i,j)      = 0.0
            nst_fld%xtts(i,j)    = 0.0
            nst_fld%xzts(i,j)    = 0.0
            nst_fld%dt_cool(i,j) = 0.0
            nst_fld%z_c(i,j)     = 0.0
            nst_fld%c_0(i,j)     = 0.0
            nst_fld%c_d(i,j)     = 0.0
            nst_fld%w_0(i,j)     = 0.0
            nst_fld%w_d(i,j)     = 0.0
            nst_fld%d_conv(i,j)  = 0.0
            nst_fld%ifd(i,j)     = 0.0
            nst_fld%Tref(i,j)    = tsea(i,j)
            nst_fld%Qrain(i,j)   = 0.0
          endif
        enddo
      enddo

            t2=timef()
!            print *,'FIXIO for nst_reset_nonwater TIME ',t2-t1,t1,t2
!
      RETURN
      END
!
!***********************************************************************
!
      subroutine interpred_phys(iord,kmsk,f,fi,global_lats_r,lonsperlar)
!!
      use resol_def,   ONLY: latr, lonr
      use layout1,     ONLY: ipt_lats_node_r, lats_node_r
      USE machine,     ONLY: kind_io8
      implicit none
!!
      integer              global_lats_r(latr)
      integer,intent(in):: iord
      integer,intent(in):: kmsk(lonr,lats_node_r)
      integer,intent(in):: lonsperlar(latr)
      real(kind=kind_io8),intent(in):: f(lonr,lats_node_r)
      real(kind=kind_io8),intent(out):: fi(lonr,lats_node_r)
      integer j,lons,lat
!!
      do j=1,lats_node_r
        lat=global_lats_r(ipt_lats_node_r-1+j)
        lons=lonsperlar(lat)
        if(lons /= lonr) then
          call intlon_phys(iord,1,1,lonr,lons,
     &                  kmsk(1,j),f(1,j),fi(1,j))
cjfe      fi(lons+1:lonr,j)=-9999.e9
          fi(lons+1:lonr,j)=0.
        else
          fi(:,j)=f(:,j)
        endif
      enddo
      end subroutine
c
c***********************************************************************
c
      subroutine intlon_phys(iord,imon,imsk,m1,m2,k1,f1,f2)
      use machine, ONLY: kind_io8
      implicit none
      integer,intent(in):: iord,imon,imsk,m1,m2
      integer,intent(in):: k1(m1)
      real (kind=kind_io8),intent(in):: f1(m1)
      real (kind=kind_io8),intent(out):: f2(m2)
      integer i2,in,il,ir
      real (kind=kind_io8) r,x1
      r=real(m1)/real(m2)
      do i2=1,m2
         x1=(i2-1)*r
         il=int(x1)+1
         ir=mod(il,m1)+1
          if(iord.eq.2.and.(imsk.eq.0.or.k1(il).eq.k1(ir))) then
            f2(i2)=f1(il)*(il-x1)+f1(ir)*(x1-il+1)
          else
            in=mod(nint(x1),m1)+1
            f2(i2)=f1(in)
          endif
      enddo
      end subroutine
c
c**********************************************************************
c
      SUBROUTINE readoz_disprd(ozplin)
 
      use ozne_def, ONLY: latsozp, levozp, timeoz, pl_coeff, kozpl
      USE machine,  ONLY: kind_phys, kind_io4
      implicit none
!!
      integer n,k,kk,i
      real (kind=kind_phys) ozplin(latsozp,levozp,pl_coeff,timeoz)
      real(kind=kind_io4) tempin(latsozp)
!
      DO I=1,timeoz
        do n=1,pl_coeff
          DO k=1,levozp
            READ(kozpl) tempin
            ozplin(:,k,n,i) = tempin(:)
          ENDDO
        enddo
      ENDDO
 
      RETURN
      END
!
!***********************************************************************
!
      SUBROUTINE ORORD(LUGB,IORO,JORO,ORO,FNOROG)
!
      use layout1, ONLY: me
      USE machine, ONLY: kind_io4, kind_io8
      implicit none
!!
      integer lugb, ioro, joro, kpdoro, ior, jor, i,k
      CHARACTER*(*) FNOROG
!
      real (kind=kind_io4) oro(ioro,joro)
      real (kind=kind_io8) orog(ioro,joro), blnm, bltm
      logical gausm
!
!     FNOROG = 'orography'
      kpdoro = 8
      IOR    = IORO
      JOR    = JORO
      CALL FIXRDG(LUGB,IOR,JOR,FNOROG,
     &            KPDORO,OROG,GAUSM,BLNM,BLTM,me)
!
      if (ior .ne. ioro .or. jor .ne. joro) then
         print *,' orography file not o.k. run aborted'
         call abort
      endif
      ORO = OROG
!
      RETURN
      END
!
!***********************************************************************
!
      subroutine split2d_phys(x,xl,global_lats_r,iop)
!
!***********************************************************************
!
      use resol_def,     ONLY: latr, lonr
      use layout1,       ONLY: me, nodes, lats_node_r, ipt_lats_node_r
      use mpi_def,       ONLY: info, mpi_r_io, mpi_comm_all
      USE machine,       ONLY: kind_io4, kind_io8
      implicit none
!!
      real(kind=kind_io4) x(lonr,latr)
      real (kind=kind_io8) xl(lonr,lats_node_r)
      real(kind=kind_io4) tmp(lonr,latr)
      integer global_lats_r(latr),iop
      integer nprocf,nodesr
!     integer maxfld,nprocf,nodesr
!     integer proc,j,lat,msgtag,nproc,i,msgtag1,buff,startlat,ierr
      integer proc,j,lat,nproc,i,buff,startlat,ierr
!     integer ifld/0/
!     save ifld
      real t1,t2,t3,t4,timef,ta,tb
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!!
      XL = 0.
!     ifld=ifld+1
      IF (me == iop) THEN
!
!         Sending the data
!         ----------------
!-- do not need to send data, all processores read the data
        tmp = 0.
        do j=1,latr
          do i=1,lonr
            tmp(i,j) = X(i,j)
          enddo
        enddo
      ENDIF

      call mpi_bcast(tmp,lonr*latr,MPI_R_IO,iop,MPI_COMM_ALL,info)
      call mpi_barrier(mpi_comm_all,info)

!-- get subdomain of data
      do j=1,lats_node_r
        lat = global_lats_r(ipt_lats_node_r-1+j)
        do i=1,lonr
          xl(i,j) = tmp(i,lat)
        enddo
      enddo

      return
      end
!
!***********************************************************************
!
      subroutine split2d_phys_r8(x,xl,global_lats_r,iop)
!
!***********************************************************************
!
      use resol_def,     ONLY: latr, lonr
      use layout1,       ONLY: me, nodes, lats_node_r, ipt_lats_node_r
      use mpi_def,       ONLY: info, mpi_r_io_r, mpi_comm_all
      USE machine,       ONLY: kind_io8
      implicit none
!!
      integer iop
      real(kind=kind_io8) x(lonr,latr)
      real (kind=kind_io8) xl(lonr,lats_node_r)
      real(kind=kind_io8) tmp(lonr,latr)
      integer global_lats_r(latr)
      integer nprocf,nodesr
!     integer maxfld,nprocf,nodesr
!     integer proc,j,lat,msgtag,nproc,i,msgtag1,buff,startlat,ierr
      integer proc,j,lat,nproc,i,buff,startlat,ierr
!     integer ifld/0/
!     save ifld
      real t1,t2,t3,t4,timef,ta,tb
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!!
      XL = 0.
!     ifld=ifld+1
      IF (me == iop) THEN
!
!         Sending the data
!         ----------------
!-- do not need to send data, all processores read the data
        tmp = 0.
        do j=1,latr
          do i=1,lonr
            tmp(i,j) = X(i,j)
          enddo
        enddo
      ENDIF

      call mpi_bcast(tmp,lonr*latr,MPI_R_IO_R,iop,MPI_COMM_ALL,info)
      print *,"in split2d_phys_r8 info=", info, "lonr, latr=", lonr, latr
      
      call mpi_barrier(mpi_comm_all,info)

!-- get subdomain of data
      do j=1,lats_node_r
        lat = global_lats_r(ipt_lats_node_r-1+j)
        do i=1,lonr
          xl(i,j) = tmp(i,lat)
        enddo
      enddo

      return
      end
!
!***********************************************************************
!
      SUBROUTINE skip(jump)
 
!*************************************************************************
 
      use resol_def
      use layout1
      use mpi_def
      implicit none
 
      integer jump,ipe
 
      ipe=0
 
      CALL MPI_BCAST(jump,1,MPI_INTEGER,ipe,MPI_COMM_ALL,info)
 
      RETURN
      END
!
!***********************************************************************
!
      subroutine excht(lats_nodes_r,global_lats_r,trcx,trcy)
!
!***********************************************************************
!
      use machine,   only: kind_io8
      use resol_def
      use layout1
      use mpi_def
      implicit none

      integer n,i,j,ierr,ilat,lat,node,nsend
      integer              global_lats_r(latr)
      integer              lats_nodes_r(nodes)
      real(kind=kind_io8)  trcx(lats_node_r,ntrac,2)
      real(kind=kind_io8)  trcy(latr,ntrac,2)
      real(kind=kind_io8) tmps(ntrac*2,lats_node_r_max,nodes)
      real(kind=kind_io8) tmpr(ntrac*2,lats_node_r_max,nodes)
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
      if (nodes /= 1) then
        do node=1,nodes
          do i=1,lats_node_r
            lat = global_lats_r(ipt_lats_node_r-1+i)
            do n=1,ntrac
              tmps(n,      i,node) = trcx(i,n,1)
              tmps(n+ntrac,i,node) = trcx(i,n,2)
            enddo
          enddo
        enddo
!
        nsend = (ntrac+ntrac)*lats_node_r_max

        call mpi_alltoall(tmps,nsend,mpi_r_def,
     &                    tmpr,nsend,mpi_r_def,
     &                    mc_comp,ierr)
!!
        ilat=1
        do node=1,nodes
          do i=1,lats_nodes_r(node)
             lat = global_lats_r(ilat)
             do n=1,ntrac
               trcy(lat,n,1) = tmpr(n,      i,node)
               trcy(lat,n,2) = tmpr(n+ntrac,i,node)
             enddo
             ilat = ilat + 1
          enddo
        enddo
!!
      else
        trcy = trcx
      endif
!!
      return
      end
!
!***********************************************************************
!
c
c***********************************************************************
c
c     SUBROUTINE EXCHA(lats_nodes_r,global_lats_r,X1,X2,Y1,Y2)
c
c***********************************************************************
c
c     use resol_def,  ONLY: latr
c     use layout1,    ONLY: nodes, lats_node_r_max, lats_node_r,
c    &                      ipt_lats_node_r
c     use mpi_def,    ONLY: mc_comp, mpi_r_def
c     USE machine,    ONLY: kind_io8
c     implicit none
c
c     integer n,i,j,ierr,ilat,lat,node,nsend
c     integer              global_lats_r(latr)
c     integer              lats_nodes_r(nodes)
c     real(kind=kind_io8) X1(lats_node_r),X2(lats_node_r)
c     real(kind=kind_io8) Y1(latr),Y2(latr)
cjfe  real(kind=kind_mpi) tmps(2,lats_node_r_max,nodes)
cjfe  real(kind=kind_mpi) tmpr(2,lats_node_r_max,nodes)
c     real(kind=kind_io8) tmps(2,lats_node_r_max,nodes)
c     real(kind=kind_io8) tmpr(2,lats_node_r_max,nodes)
c
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c
c     if (nodes.ne.1) then
c       do node=1,nodes
c         do i=1,lats_node_r
c          lat=global_lats_r(ipt_lats_node_r-1+i)
c          tmps(1,i,node)=X1(I)
c          tmps(2,i,node)=X2(I)
c         enddo
c       enddo
c!
c       nsend=2*lats_node_r_max
cjfe    call mpi_alltoall(tmps,nsend,MPI_R_MPI,
cjfe x                     tmpr,nsend,MPI_R_MPI,
cjfe x                     MC_COMP,ierr)
c       call mpi_alltoall(tmps,nsend,MPI_R_DEF,
c    x                     tmpr,nsend,MPI_R_DEF,
c    x                     MC_COMP,ierr)
c!
c       ilat=1
c       do node=1,nodes
c         do i=1,lats_nodes_r(node)
c            lat=global_lats_r(ilat)
c            Y1(lat)=tmpr(1,i,node)
c            Y2(lat)=tmpr(2,i,node)
c            ilat=ilat+1
c         enddo
c       enddo
c!
c     ELSE
c       Y1=X1
c       Y2=X2
c     ENDIF
c!
c     RETURN
c     END
c
c***********************************************************************
c
c     SUBROUTINE SUMLAT(n,X,nodes)
c
c***********************************************************************
c
c     use mpi_def,   ONLY: MC_COMP, MPI_R_DEF, info, mpi_sum
c     USE machine,   ONLY: kind_io8, kind_io4
c     implicit none
c
c     integer n,i,j,np,mr,nodes
c     real(kind=kind_io8) X(n),Y(N)
c     real(kind=kind_io4) Z(n)
c
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c
c     if (nodes.ne.1) then
c       DO i=1,n
c         Y(i)=X(i)
c       ENDDO
c       CALL mpi_allreduce(Y,X,n,MPI_R_DEF,MPI_SUM,
c    &                    MC_COMP   ,info)
c     endif
c       DO i=1,n
c         Z(i)=X(i)
c       ENDDO
c       DO i=1,n
c         X(i)=Z(i)
c       ENDDO
c!
c     RETURN
c     END
c
c***********************************************************************
c
      subroutine unsplit2d_phys(ioproc,x,xl,global_lats_r)
c
c***********************************************************************
c
      use resol_def,   ONLY: latr, lonr
      use layout1,     ONLY: me, lats_node_r, lats_node_r_max,
     &                       ipt_lats_node_r, nodes
      use mpi_def,     ONLY: info, mpi_comm_all, liope, mpi_r_io,
     &                       stat
      USE machine,     ONLY: kind_io4, kind_io8
      implicit none
!!
      real(kind=kind_io4) x(lonr,latr)
      real (kind=kind_io8) xl(lonr,lats_node_r)
      real(kind=kind_io4) tmp(lonr,latr+2)
      integer global_lats_r(latr),ipt_lats_node_rl,nodesr
      integer lats_nodes_rl
      integer maxfld,ioproc,nproct
      integer proc,j,lat,msgtag,nproc,i,buff,startlat,ierr
!     integer proc,j,lat,msgtag,nproc,i,msgtag1,buff,startlat,ierr
      integer ifldu/0/
      save ifldu
      integer illen,ncc
      data ncc/0/
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!!
!     X      = 0.
      maxfld = 50
      ifldu  = ifldu + 1
!!
      IF (me /= ioproc) THEN    ! all fcst node need to send data
!
!         Sending the data
!         ----------------
         tmp = 0.
         tmp(lonr,latr+1) = ipt_lats_node_r
         tmp(lonr,latr+2) = lats_node_r
!$omp parallel do private(i,j)
         do j=1,lats_node_r
            do i=1,lonr
              tmp(i,j) = XL(i,j)
            enddo
         enddo
         if (.NOT. LIOPE) then
           nodesr = nodes
         else
           nodesr = nodes+1
         endif
!        msgtag = 1000 + (me+1)*nodesr*maxfld + ifldu
         msgtag = 1000 + (me+1)*nodesr        + ifldu
         msgtag = me + 1
!     write(0,*)' calling mpi_send in unsplit2d_phys for'
!    &,' pe=',me,' msgtag=',msgtag
         call MPI_SEND(tmp(lonr,latr+1),1,MPI_R_IO,ioproc,
     &                 msgtag,MPI_COMM_ALL,info)
          call MPI_SEND(tmp(lonr,latr+2),1,MPI_R_IO,ioproc,
     &                  msgtag,MPI_COMM_ALL,info)
         illen = tmp(lonr,latr+2)
! send the local grid domain
         CALL mpi_send(tmp(1,1),illen*lonr,MPI_R_IO,ioproc,
     &                  msgtag,MPI_COMM_ALL,info)

      ELSE            !     for pes ioproc
        x = 0.
        if (.NOT.LIOPE) then
          nproct = nodes
!$omp parallel do private(i,j,lat)
          do j=1,lats_node_r
             lat = global_lats_r(ipt_lats_node_r-1+j)
             do i=1,lonr
                x(i,lat) = XL(i,j)
             enddo
          enddo
        else
          nproct = nodes - 1
        endif
        DO proc=1,nproct
          if (proc.ne.ioproc+1) then
!           msgtag = 1000 + proc*nodes*maxfld + ifldu
!           msgtag = 1000 + proc*nodesr       + ifldu
            msgtag = proc
            CALL mpi_recv(tmp(lonr,latr+1),1,MPI_R_IO,proc-1,
     &                    msgtag,MPI_COMM_ALL,stat,info)
            CALL mpi_recv(tmp(lonr,latr+2),1,MPI_R_IO,proc-1,
     &                    msgtag,MPI_COMM_ALL,stat,info)
            illen = tmp(lonr,latr+2)
            CALL mpi_recv(tmp(1,1),illen*lonr ,MPI_R_IO,proc-1,
     &                     msgtag,MPI_COMM_ALL,stat,info)
            if (.NOT.LIOPE) then
              ipt_lats_node_rl = tmp(lonr,latr+1)
              lats_nodes_rl    = tmp(lonr,latr+2)
            else
              ipt_lats_node_rl = tmp(lonr,lats_node_r_max+1)
              lats_nodes_rl    = tmp(lonr,lats_node_r_max+2)
            endif
!$omp parallel do private(i,j,lat)
            do j=1,lats_nodes_rl
              lat = global_lats_r(ipt_lats_node_rl-1+j)
              do i=1,lonr
                x(i,lat) = tmp(i,j)
              enddo
            enddo
          endif   !(proc.ne.ioproc+1)
        enddo

      ENDIF
      ncc = ncc + 1

!!
      return
      end
!
!***********************************************************************
!
      subroutine uninterpred(iord,kmsk,f,fi,global_lats_r,lonsperlar)
!!
      use resol_def,   ONLY: latr, lonr
      use layout1,     ONLY: lats_node_r, ipt_lats_node_r
      USE machine,     ONLY: kind_io8
      implicit none
!!
      integer              global_lats_r(latr)
      integer,intent(in):: iord
      integer,intent(in):: kmsk(lonr,lats_node_r)
      integer,intent(in):: lonsperlar(latr)
      real(kind=kind_io8),intent(out):: f(lonr,lats_node_r)
      real(kind=kind_io8),intent(in) :: fi(lonr,lats_node_r)
      integer j,lons,lat
!!
      do j=1,lats_node_r
        lat  = global_lats_r(ipt_lats_node_r-1+j)
        lons = lonsperlar(lat)
        if(lons .ne. lonr) then
          call intlon_phys(iord,1,1,lons,lonr,
     &                     kmsk(1,j),fi(1,j),f(1,j))
        else
          f(:,j) = fi(:,j)
        endif
      enddo
      end subroutine



      subroutine uninterprez(iord,kmsk,f,fi,global_lats_r,lonsperlar,  
     &                       buff_mult_piecea)
!!
      use resol_def,   ONLY: latr, lonr
      use layout1,     ONLY: lats_node_r, ipt_lats_node_r
      USE machine,     ONLY: kind_io4,kind_io8
      implicit none
!!
      integer,intent(in), dimension(latr) :: global_lats_r,lonsperlar
      integer,intent(in):: iord
      integer,intent(in):: kmsk(lonr,lats_node_r)
      real(kind=kind_io8),intent(out):: f(lonr,lats_node_r)
      real(kind=kind_io8),intent(in) :: fi(lonr,lats_node_r)
!
      real(kind=kind_io4),intent(inout) :: buff_mult_piecea
     &                                     (1:lonr,1:lats_node_r)
      integer i,j,lons,lat
!!
!$omp parallel do private(j,lat,lons)
      do j=1,lats_node_r
        lat  = global_lats_r(ipt_lats_node_r-1+j)
        lons = lonsperlar(lat)
        if(lons .ne. lonr) then
          call intlon_phys(iord,1,1,lons,lonr,
     &                     kmsk(1,j),fi(1,j),f(1,j))
        else
          f(:,j) = fi(:,j)
        endif
      enddo
!$omp parallel do private(i,j)
      do j=1,lats_node_r
        do i=1,lonr
          buff_mult_piecea(i,j) = f (i,j)
        end do
      end do
      end subroutine


      subroutine unsplit2z(ngridx,ngridt,x,global_lats_r)
!
!***********************************************************************
!
      use resol_def,   ONLY: lonr,latr
      use mod_state,   ONLY: ivar_global_a, buff_mult_pieces
      use layout1,     ONLY: me, nodes_comp
      use mpi_def,     ONLY: liope
      USE machine,     ONLY: kind_io4
      implicit none
!!
      real(kind=kind_io4) x(lonr,latr)
      integer             global_lats_r(latr),ngridx,ngridt

!     Locals
      integer i,j,proc,lat,lats_nodes_rl,ipt_lats_node_rl,illen,nd1,nd2
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!!
      X = 0.
!
      nd1 = 0
      DO proc=1,nodes_comp
        ipt_lats_node_rl = ivar_global_a(1,proc)
        lats_nodes_rl    = ivar_global_a(2,proc)
        nd2 = nd1 + lonr*lats_nodes_rl*(ngridx-1)
        do j=1,lats_nodes_rl
          lat = global_lats_r(ipt_lats_node_rl-1+j)
          do i=1,lonr
            x(i,lat) = buff_mult_pieces(nd2+i+(j-1)*lonr)
          enddo
        enddo
        nd1 = nd1 + lonr*lats_nodes_rl*ngridt
      enddo
!!
      return
      end
 
!
!***********************************************************************
!
      subroutine unsplit2d_phys_r(ioproc,x,xl,global_lats_r)
!
!***********************************************************************
!
      use resol_def,   ONLY: latr, lonr
      use layout1,     ONLY: me, lats_node_r, lats_node_r_max, 
     &                       ipt_lats_node_r, nodes
      use mpi_def,     ONLY: liope, info, stat, mpi_comm_all, 
     &                       mpi_r_io_r
      USE machine,     ONLY: kind_ior, kind_io8
      implicit none
!!
      real(kind=kind_ior) x(lonr,latr)
      real (kind=kind_io8) xl(lonr,lats_node_r)
      real(kind=kind_ior) tmp(lonr,latr+2)
      integer global_lats_r(latr),ipt_lats_node_rl,nodesr
      integer lats_nodes_rl
      integer maxfld,ioproc,nproct
      integer proc,j,lat,msgtag,nproc,i,buff,startlat,ierr
!     integer proc,j,lat,msgtag,nproc,i,msgtag1,buff,startlat,ierr
      integer ifldu/0/
      save ifldu
      integer illen,ncc
      data ncc/0/
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!!
!     X = 0.               ! commented by moorthi on 20051117
      maxfld = 50
      ifldu  = ifldu + 1
!!
      IF (me.ne.ioproc) THEN
!
!         Sending the data
!         ----------------
         tmp = 0.
         tmp(lonr,latr+1) = ipt_lats_node_r
         tmp(lonr,latr+2) = lats_node_r
!$omp parallel do private(i,j)
         do j=1,lats_node_r
           do i=1,lonr
             tmp(i,j) = XL(i,j)
           enddo
         enddo
         if (.NOT.LIOPE) then
           nodesr = nodes
         else
           nodesr = nodes + 1
         endif
!        msgtag = 1000 + (me+1)*nodesr*maxfld + ifldu
!        msgtag = 1000 + (me+1)*nodesr        + ifldu
         msgtag = me+1
!        write(0,*)'sending data to ioproc=',ioproc, 'mes=',msgtag
!    &,' me=',me
!        write(0,*)'sending data to ioproc=',ioproc, mes=',msgtag,
!     &    'nodesr=',nodesr,'maxfld=',maxfld,'ifldu=',ifldu,
!     &    'liope=',liope
         call MPI_SEND(tmp(lonr,latr+1),1,MPI_R_IO_R,ioproc,
     &                 msgtag,MPI_COMM_ALL,info)
         call MPI_SEND(tmp(lonr,latr+2),1,MPI_R_IO_R,ioproc,
     &                 msgtag,MPI_COMM_ALL,info)
         illen = tmp(lonr,latr+2)
         
! send the local grid domain
         CALL mpi_send(tmp(1,1),illen*lonr,MPI_R_IO_R,ioproc,
     &                 msgtag,MPI_COMM_ALL,info)

      ELSE        !     for pes ioproc
        x = 0.0               ! added by Moorthi on 2005111700
        if (.NOT.LIOPE) then
          nproct = nodes
!$omp parallel do private(i,j,lat)
          do j=1,lats_node_r
             lat = global_lats_r(ipt_lats_node_r-1+j)
             do i=1,lonr
                x(i,lat) = XL(i,j)
             enddo
          enddo
        else
          nproct = nodes - 1
        endif
        DO proc=1,nproct
          if (proc.ne.ioproc+1) then
!           msgtag = 1000 + proc*nodes*maxfld + ifldu
!           msgtag = 1000 + proc*nodesr       + ifldu
            msgtag = proc
!          print *,'receive data fm pe=',proc-1,'mes=',msgtag,
!    &     'nodes=',nodes,'maxfld=',maxfld,'ifldu=',ifldu
            CALL mpi_recv(tmp(lonr,latr+1),1,MPI_R_IO_R,proc-1,
     &                    msgtag,MPI_COMM_ALL,stat,info)
            CALL mpi_recv(tmp(lonr,latr+2),1,MPI_R_IO_R,proc-1,
     &                    msgtag,MPI_COMM_ALL,stat,info)
            illen = tmp(lonr,latr+2)
            CALL mpi_recv(tmp(1,1),illen*lonr ,MPI_R_IO_R,proc-1,
     &                   msgtag,MPI_COMM_ALL,stat,info)

            if (.NOT.LIOPE) then
              ipt_lats_node_rl = tmp(lonr,latr+1)
              lats_nodes_rl    = tmp(lonr,latr+2)
            else
              ipt_lats_node_rl = tmp(lonr,lats_node_r_max+1)
              lats_nodes_rl    = tmp(lonr,lats_node_r_max+2)
            endif

!!$omp parallel do private(i,j,lat)
            do j=1,lats_nodes_rl
              lat = global_lats_r(ipt_lats_node_rl-1+j)
!     write(0,*)' in unsplit2d_phys_r j=',j,' lat=',lat
              do i=1,lonr
                x(i,lat) = tmp(i,j)
              enddo
            enddo
          endif   !(proc.ne.ioproc+1)
        enddo
!!
      ENDIF
      ncc = ncc + 1
 
!!
      return
      end
c
c***********************************************************************
c
      subroutine split2d_phys_r(x,xl,global_lats_r)
c
c***********************************************************************
c
      use resol_def,      ONLY: latr, lonr
      use layout1,        ONLY: me, lats_node_r, ipt_lats_node_r, nodes
      use mpi_def,        ONLY: liope, mpi_comm_all, info,mpi_r_io_r
      USE machine,        ONLY: kind_ior, kind_io8
      implicit none
!!
      real(kind=kind_ior) x(lonr,latr)
      real (kind=kind_io8) xl(lonr,lats_node_r)
      real(kind=kind_ior) tmp(lonr,latr)
      integer global_lats_r(latr)
      integer nprocf,nodesr
!     integer maxfld,nprocf,nodesr
      integer proc,j,lat,nproc,i,buff,startlat,ierr
!     integer proc,j,lat,msgtag,nproc,i,msgtag1,buff,startlat,ierr
!     integer ifld/0/
!     save ifld
      real t1,t2,t3,t4,timef,ta,tb
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!!
      XL=0.
!     maxfld=50
!     ifld=ifld+1
!!
      IF (me == 0) THEN
        ta = timef()
        t3 = ta
!       DO proc=1,nodes-1
        do proc=1,1
!
!         Sending the data
!         ----------------
          tmp=0.
          do j=1,latr
            do i=1,lonr
              tmp(i,j)=X(i,j)
            enddo
          enddo
!Moor    msgtag=1000+proc*nodes*maxfld+ifld
          t1 = timef()
!sela    print *,' GWVX BROADCASTING FROM ',nodes-1
          call mpi_bcast
     &        (tmp,lonr*latr,MPI_R_IO_R,nodes-1,MPI_COMM_ALL,info)
          call mpi_comm_rank(MPI_COMM_ALL,i,info)
c         CALL mpi_send(tmp,lonr*latr,MPI_R_IO_R,proc-1,msgtag,
c     &                  MPI_COMM_ALL,info)
          t2 = timef()
!sela    print 102,t2-t1
 
 102      format(' SEND TIME ',f10.5)
        enddo
        t4 = timef()
      ELSE
        if (.NOT. LIOPE) then
          nodesr = nodes
        else
          nodesr = nodes+1
        endif
!Moor   msgtag=1000+(me+1)*nodesr*maxfld+ifld
!sela   print *,' GWVX BROADCASTREC  FROM ',nodesr-1
        call mpi_bcast
     &       (tmp,lonr*latr,MPI_R_IO_R,nodesr-1,MPI_COMM_ALL,info)
        call mpi_comm_rank(MPI_COMM_ALL,i,info)
!sela   print *,'GWVX IPT ',ipt
c       CALL mpi_recv(tmp,lonr*latr,MPI_R_IO_R,nodesr-1,
c    &                msgtag,MPI_COMM_ALL,stat,info)
        do j=1,lats_node_r
          lat=global_lats_r(ipt_lats_node_r-1+j)
          do i=1,lonr
             xl(i,j)=tmp(i,lat)
          enddo
        enddo
!!
      ENDIF
!!
!!     for pes nodes-1
      if (.NOT.LIOPE) then
        if (me.eq.nodes-1) then
          do j=1,lats_node_r
             lat=global_lats_r(ipt_lats_node_r-1+j)
             do i=1,lonr
                xl(i,j)=X(i,lat)
             enddo
          enddo
        endif
      endif
!!
      tb=timef()
      call mpi_comm_rank(MPI_COMM_ALL,i,info)
 
!sela  if(icolor.eq.2.and.me.eq.nodes-1)print 103,tb-ta,t4-t3
 103  format(' GLOBAL AND SEND TIMES  split2d_phys',2f10.5)
      return
      end

!
c***********************************************************************
c
      subroutine split2d_rst(x,xl,fieldsize,global_lats_r,lonsperlar)
c
c***********************************************************************
c
      use resol_def,      ONLY: latr, lonr
      use layout1,        ONLY: me, lats_node_r, ipt_lats_node_r, nodes
      use mpi_def,        ONLY: liope, mpi_comm_all, info,mpi_r_io_r
      USE machine,        ONLY: kind_ior, kind_io8
      implicit none
!!
!!
      integer,intent(in) :: fieldsize,global_lats_r(latr),
     &                      lonsperlar(latr)
      real(kind=kind_ior), intent(in)    :: x(fieldsize)
      real (kind=kind_io8),intent(inout) :: xl(lonr,lats_node_r)
      integer j,lat,i,lon
!      real t1,t2,t3,t4,timef,ta,tb
c
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c
!!
!--- get subdomain of data
       do j=1,lats_node_r
         lat = global_lats_r(ipt_lats_node_r-1+j)
         if(lat /= 1) then
           lon = sum(lonsperlar(1:lat-1))
         else
           lon = 0
         endif
!
         do i=1,lonsperlar(lat)
           xl(i,j) = X(lon+i)
         enddo
       enddo
!!
!sela  if(icolor.eq.2.and.me.eq.nodes-1)print 103,tb-ta,t4-t3
! 103  format(' GLOBAL AND SEND TIMES  split2d_phys',2f10.5)

      return
      end subroutine split2d_rst


!***********************************************************************
!
      SUBROUTINE read_sfc_r(cfile,sfc_fld,phy_f2d,phy_f3d,
     &           NGPTC,NBLCK,global_lats_r,lonsperlar,NEEDORO,
     &           lsidea,pr_idea,gg,prsilvl,amgms)
!
!***********************************************************************
!
! !REVISION HISTORY:
!
!  2011/9/26   Jun Wang add cv/cvt/cvb
!  2013/3/06   Jun Wang add idea variables to restart file

      use resol_def,      ONLY: latr, lonr, latr2, lsoil,levs
     &,                         ntot3d, ntot2d
      use layout1,        ONLY: me, nodes, lats_node_r,ipt_lats_node_r
      USE machine,        ONLY: kind_ior, kind_io8, kind_rad

      use gfs_physics_sfc_flx_mod, ONLY: Sfc_Var_Data
      use namelist_soilveg ,       only: salp_data, snupx
      use physcons,                only : tgice => con_tice
      use nemsio_module
!
      implicit none
!
      character(*),intent(in) :: cfile
      TYPE(Sfc_Var_Data),intent(inout) :: sfc_fld
      integer,intent(in)               :: global_lats_r(latr)
     &,                                   lonsperlar(latr)
     &,                                   NGPTC,NBLCK
      real(kind=kind_rad),intent(inout) ::
     &    phy_f2d(lonr,lats_node_r,ntot2d),
     &    phy_f3d(NGPTC,LEVS,ntot3d,NBLCK,lats_node_r)
      integer,intent(inout) :: needoro
      logical,intent(in)    :: lsidea
      real(kind=kind_rad),intent(inout) :: pr_idea(levs),gg(levs),
     &                                     amgms(levs), prsilvl(levs+1)
!
      integer jump

      real(kind=kind_io8) buff3(lonr,lats_node_r)
!
      real(kind=kind_ior),allocatable :: buff1(:)
!
      integer i,j,k,im,jm,idate(4),lplsfc(latr2)
      real t1,t2,timef,rsnow
!---
      type(nemsio_gfile) :: gfile
      integer iret, vegtyp,fieldsize,iblk,il,lons_lat,njeff,l,lat,lon
      character*2 nump2d,nump3d
      character(255) varname
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
      t1=timef()
!
      call nemsio_init()
!
      call nemsio_open(gfile,trim(cfile),'read',iret=iret)
!      print *,'after nemsio_open, iret=',iret
      if(iret/=0) then
        PRINT *, ' ERROR in input routine read_sfc_r'
        return
      endif
!
      call nemsio_getfilehead(gfile,dimx=im,dimy=jm,iret=iret)
      fieldsize=im*jm
      allocate(buff1(fieldsize))
!
      if(lsidea) then
        call nemsio_getheadvar(gfile,"pr_idea",pr_idea,iret=iret)
        if(iret/=0) print *,'idea, restart, get plyr wrong!'
        print *,'in read rst sfc, pr_idea=',pr_idea(1:levs)
        call nemsio_getheadvar(gfile,"gravity",gg,iret=iret)
        if(iret/=0) print *,'idea, restart, get gravity wrong!'
        call nemsio_getheadvar(gfile,"amgms",amgms,iret=iret)
        if(iret/=0) print *,'idea, restart, get amgms wrong!'
        call nemsio_getheadvar(gfile,"prsilvl",prsilvl,iret=iret)
        if(iret/=0) print *,'idea, restart, get prsilvl wrong!'
      endif

!
!-- tsea
      call nemsio_readrecv(gfile,'tmp','sfc',1,buff1,iret=iret)
      call split2d_rst(buff1,sfc_fld%TSEA,fieldsize,global_lats_r,
     &  lonsperlar)
!-- smc
      DO K=1, LSOIL
        call nemsio_readrecv(gfile,'smc','soil layer',k,buff1,iret=iret)
        call split2d_rst(buff1, sfc_fld%smc(k,:,:),fieldsize,
     &    global_lats_r,lonsperlar)
!        print *,'read inrst,smc=',sfc_fld%smc(k,1:5,1:5)
      ENDDO

!-- weasd
      call nemsio_readrecv(gfile,'weasd','sfc',1,buff1,iret=iret)
      call split2d_rst(buff1,sfc_fld%weasd,fieldsize,global_lats_r,
     &  lonsperlar)
!--stc
      DO K = 1, LSOIL
        call nemsio_readrecv(gfile,'stc','soil layer',k,buff1,iret=iret)
        call split2d_rst(buff1, sfc_fld%stc(k,:,:),fieldsize,
     &    global_lats_r,lonsperlar)
!        print *,'read inrst,stc=',sfc_fld%stc(k,1:5,1:5)
      ENDDO

!--tg3
      call nemsio_readrecv(gfile,'tg3','sfc',1,buff1,iret=iret)
      call split2d_rst(buff1,sfc_fld%tg3,fieldsize,global_lats_r,
     &  lonsperlar)
!        print *,'read inrst,tg3=',sfc_fld%tg3(1:3,1:3)
!--zorl
      call nemsio_readrecv(gfile,'sfcr','sfc',1,buff1,iret=iret)
      call split2d_rst(buff1,sfc_fld%zorl,fieldsize,global_lats_r,
     &  lonsperlar)
!
!--cv
      call nemsio_readrecv(gfile,'tcdc','convect-cld laye',1,
     &  buff1,iret=iret)
      call split2d_rst(buff1,sfc_fld%cv,fieldsize,global_lats_r,
     &  lonsperlar)
      if(me==0) print *,'read inrst from rst,cwafter cv=',
     &  maxval(sfc_fld%cv),minval(sfc_fld%cv)
!--cvb
      call nemsio_readrecv(gfile,'pres','convect-cld bot',1,
     &  buff1,iret=iret)
      call split2d_rst(buff1,sfc_fld%cvb,fieldsize,global_lats_r,
     &  lonsperlar)
!--cvt
      call nemsio_readrecv(gfile,'pres','convect-cld top',1,
     &  buff1,iret=iret)
      call split2d_rst(buff1,sfc_fld%cvt,fieldsize,global_lats_r,
     &  lonsperlar)
!        print *,'read inrst,cwafter cvt'
!-- alvsf
      call nemsio_readrecv(gfile,'alvsf','sfc',1,buff1,iret=iret)
      call split2d_rst(buff1,sfc_fld%alvsf,fieldsize,global_lats_r,
     &  lonsperlar)
!-- alvwf
      call nemsio_readrecv(gfile,'alvwf','sfc',1,buff1,iret=iret)
      call split2d_rst(buff1,sfc_fld%alvwf,fieldsize,global_lats_r,
     &  lonsperlar)
!-- alnsf
      call nemsio_readrecv(gfile,'alnsf','sfc',1,buff1,iret=iret)
      call split2d_rst(buff1,sfc_fld%alnsf,fieldsize,global_lats_r,
     &  lonsperlar)
!--alnwf
      call nemsio_readrecv(gfile,'alnwf','sfc',1,buff1,iret=iret)
      call split2d_rst(buff1,sfc_fld%alnwf,fieldsize,global_lats_r,
     &  lonsperlar)
!-- slmsk
      call nemsio_readrecv(gfile,'land','sfc',1,buff1,iret=iret)
      call split2d_rst(buff1,sfc_fld%slmsk,fieldsize,global_lats_r,
     &  lonsperlar)

!-- vfrac
      call nemsio_readrecv(gfile,'veg','sfc',1,buff1,iret=iret)
      call split2d_rst(buff1,sfc_fld%vfrac,fieldsize,global_lats_r,
     &  lonsperlar)
!-- canopy
      call nemsio_readrecv(gfile,'cnwat','sfc',1,buff1,iret=iret)
      call split2d_rst(buff1,sfc_fld%canopy,fieldsize,global_lats_r,
     &  lonsperlar)
!-- f10m
      call nemsio_readrecv(gfile,'f10m','10 m above gnd',1,buff1,
     &   iret=iret)
      call split2d_rst(buff1,sfc_fld%f10m,fieldsize,global_lats_r,
     &  lonsperlar)
!--vtype
      call nemsio_readrecv(gfile,'vtype','sfc',1,buff1,iret=iret)
      call split2d_rst(buff1,sfc_fld%vtype,fieldsize,global_lats_r,
     &  lonsperlar)
!-- stype
      call nemsio_readrecv(gfile,'sotyp','sfc',1,buff1,iret=iret)
      call split2d_rst(buff1,sfc_fld%stype,fieldsize,global_lats_r,
     &  lonsperlar)
!-- facsf
      call nemsio_readrecv(gfile,'facsf','sfc',1,buff1,iret=iret)
      call split2d_rst(buff1,sfc_fld%facsf,fieldsize,global_lats_r,
     &  lonsperlar)
!-- facwf
      call nemsio_readrecv(gfile,'facwf','sfc',1,buff1,iret=iret)
      call split2d_rst(buff1,sfc_fld%facwf,fieldsize,global_lats_r,
     &  lonsperlar)
!-- uustar (fricv)
      call nemsio_readrecv(gfile,'fricv','sfc',1,buff1,iret=iret)
      call split2d_rst(buff1,sfc_fld%uustar,fieldsize,global_lats_r,
     &  lonsperlar)
!-- ffhh
      call nemsio_readrecv(gfile,'ffhh','sfc',1,buff1,iret=iret)
      call split2d_rst(buff1,sfc_fld%ffhh,fieldsize,global_lats_r,
     &  lonsperlar)
!-- ffmm
      call nemsio_readrecv(gfile,'ffmm','sfc',1,buff1,iret=iret)
      call split2d_rst(buff1,sfc_fld%ffmm,fieldsize,global_lats_r,
     &  lonsperlar)
!-- hice
      call nemsio_readrecv(gfile,'icetk','sfc',1,buff1,iret=iret)
      call split2d_rst(buff1,sfc_fld%hice,fieldsize,global_lats_r,
     &  lonsperlar)
!-- fice
      call nemsio_readrecv(gfile,'icec','sfc',1,buff1,iret=iret)
      call split2d_rst(buff1,sfc_fld%fice,fieldsize,global_lats_r,
     &  lonsperlar)
!-- tisfc
      call nemsio_readrecv(gfile,'tisfc','sfc',1,buff1,iret=iret)
      call split2d_rst(buff1,sfc_fld%tisfc,fieldsize,global_lats_r,
     &  lonsperlar)
!        print *,'read inrst,tisfc=',sfc_fld%tisfc(1:3,1:3)
      if (lats_node_r > 0 )  then
        if (sfc_fld%tisfc(1,1) < 0.0) then
          DO j=1,lats_node_r
            DO i=1,LONR
               sfc_fld%TISFC(i,j) = sfc_fld%TSEA(i,j)
               IF(sfc_fld%SLMSK(i,j) >=  2. .AND.
!    &            sfc_fld%FICE(i,j)  >= 0.5) THEN
     &            sfc_fld%FICE(i,j)  >= 0.15) THEN
                  sfc_fld%TISFC(i,j) = (sfc_fld%TSEA(i,j)
     &           -tgice*(1.-sfc_fld%FICE(i,j))) / sfc_fld%FICE(i,j)
                  sfc_fld%TISFC(i,j) = MIN(sfc_fld%TISFC(i,j),tgice)
               ENDIF
            ENDDO
          ENDDO
        endif
      endif
!-- tprcp
      call nemsio_readrecv(gfile,'tprcp','sfc',1,buff1,iret=iret)
      call split2d_rst(buff1,sfc_fld%tprcp,fieldsize,global_lats_r,
     &  lonsperlar)
!-- srflag (crain)
      call nemsio_readrecv(gfile,'crain','sfc',1,buff1,iret=iret)
      call split2d_rst(buff1,sfc_fld%srflag,fieldsize,global_lats_r,
     &  lonsperlar)
!-- snwdph
      call nemsio_readrecv(gfile,'snod','sfc',1,buff1,iret=iret)
      call split2d_rst(buff1,sfc_fld%SNWDPH,fieldsize,global_lats_r,
     &  lonsperlar)
!-- slc
      DO K=1, LSOIL
        call nemsio_readrecv(gfile,'slc','soil layer',k,buff1,iret=iret)
        call split2d_rst(buff1,sfc_fld%slc(k,:,:),fieldsize,
     &    global_latS_r,lonsperlar)
!        print *,'read inrst,slc=',sfc_fld%slc(k,1:3,1:3)
      ENDDO
!-- shdmin
      call nemsio_readrecv(gfile,'shdmin','sfc',1,buff1,iret=iret)
      call split2d_rst(buff1,sfc_fld%shdmin,fieldsize,global_lats_r,
     &  lonsperlar)
!-- shdmax
      call nemsio_readrecv(gfile,'shdmax','sfc',1,buff1,iret=iret)
      call split2d_rst(buff1,sfc_fld%shdmax,fieldsize,global_lats_r,
     &  lonsperlar)
!-- slope (sltyp)
      call nemsio_readrecv(gfile,'sltyp','sfc',1,buff1,iret=iret)
      call split2d_rst(buff1,sfc_fld%slope,fieldsize,global_lats_r,
     &  lonsperlar)
!-- salbd
      call nemsio_readrecv(gfile,'salbd','sfc',1,buff1,iret=iret)
      call split2d_rst(buff1,sfc_fld%SNOALB,fieldsize,global_lats_r,
     &  lonsperlar)
!        print *,'read inrst,snoalb=',sfc_fld%snoalb(1:3,1:3)
!-- orog
      if(needoro.eq.1) then
        call nemsio_readrecv(gfile,'orog','sfc',1,buff1,iret=iret)
        needoro=1
        if(any(buff1.eq.-9999.)) needoro=0
!        print *,'read sfc orography'
        call split2d_rst(buff1,sfc_fld%oro,fieldsize,global_lats_r,
     &  lonsperlar)
        call skip(needoro)
      endif
!        print *,'read inrst,after orog'
!jw read sncovr from rstart file
!-- read in snow cover from restart file
      sfc_fld%SNCOVR = 0.0
      call nemsio_readrecv(gfile,'sncovr','sfc',1,buff1,iret=iret)
      if(iret==0)
     &call split2d_rst(buff1,sfc_fld%sncovr,fieldsize,global_lats_r,
     &  lonsperlar)
!        print *,'read inrst,snoalb=',sfc_fld%sncovr(38,3),
!     &    sfc_fld%weasd(38,3)
!
!-- ntot2d
      DO K=1, ntot2d
        write(nump2d,'(I2.2)')k
        varname='phyf2d_'//nump2d
        call nemsio_readrecv(gfile,trim(varname),'sfc',1,buff1,
     &    iret=iret)
!        print *,'read inrst,',trim(varname),'iret=',iret
        call split2d_rst(buff1,phy_f2d(:,:,k),fieldsize,global_lats_r,
     &    lonsperlar)
      ENDDO
!
!-- ntot3d == num_p3d+npdf3d+nshoc_3d+ncnvcld3d
      do k=1, ntot3d
        write(nump3d,'(I2.2)')k
        varname='phyf3d_'//nump3d
        do l=1, levs
          call nemsio_readrecv(gfile,trim(varname),'mid layer',L,
     &      buff1,iret=iret)
!        print *,'read inrst,phy_p3d,',trim(varname),'iret=',iret
          call split2d_rst(buff1,buff3,fieldsize,global_lats_r,
     &    lonsperlar)
!
          do j=1,lats_node_r
            lat = global_lats_r(ipt_lats_node_r-1+j)
            lons_lat = lonsperlar(lat)
            iblk = 0
            il   = 1
            do lon=1,lons_lat,NGPTC
              NJEFF=MIN(NGPTC,lons_lat-lon+1)
              iblk = iblk + 1
              do i=1,NJEFF
                phy_f3d(i,l,k,iblk,j) = buff3(il,j)
                il = il + 1
              enddo
            enddo
          enddo
!
        ENDDO
      ENDDO

      call nemsio_close(gfile)
      call nemsio_finalize()
!
      t2=timef()
!      print *,'FIXIO TIME ',t2-t1,t1,t2
!
      RETURN

      STOP
      END
!
      SUBROUTINE read_nst_r(nst_fld, nread, cfile,
     &                     global_lats_r, lonsperlar)
!
!***********************************************************************
!
      use namelist_physics_def
      USE machine,        ONLY: kind_ior, kind_io8, kind_rad
      use resol_def
      use layout1
      use mpi_def
      use gfs_physics_nst_var_mod
      use nemsio_module
      implicit none
!
      TYPE(Nst_Var_Data)       :: nst_fld
      integer              global_lats_r(latr)
      integer                 lonsperlar(latr)

!     real (kind=kind_io8) slmsk(lonr,lats_node_r),

      real(kind=kind_ior),allocatable :: buff1(:)
      real(kind=kind_io8) buffo(lonr,lats_node_r)
      integer nread,i,j,k,ij,idate(4),lonnst,latnst,lplnst(latr2)
      character*(*) cfile
      integer kmsk(lonr,latr)
      CHARACTER*8 labfix(4)
      real t1,t2,timef
!---
      type(nemsio_gfile) :: gfile
      integer iret, fieldsize, im, jm
      character(255) varname

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
      t1=timef()
!
      call nemsio_init()
!
      call nemsio_open(gfile,trim(cfile),'read',iret=iret)
!      print *,'after nemsio_open, iret=',iret
      if(iret /= 0) then
        PRINT *, ' ERROR in input routine read_sfc_r'
        return
      endif
!
      call nemsio_getfilehead(gfile,dimx=im,dimy=jm,iret=iret)
      fieldsize = im*jm
      allocate(buff1(fieldsize))
!

!-- xt
      call nemsio_readrecv(gfile,'xt','sfc',1,buff1,iret=iret)
      call split2d_rst(buff1,nst_fld%xt,fieldsize,global_lats_r,
     &  lonsperlar)

!-- xs
      call nemsio_readrecv(gfile,'xs','sfc',1,buff1,iret=iret)
      call split2d_rst(buff1,nst_fld%xs,fieldsize,global_lats_r,
     &  lonsperlar)

!-- xu
      call nemsio_readrecv(gfile,'xu','sfc',1,buff1,iret=iret)
      call split2d_rst(buff1,nst_fld%xu,fieldsize,global_lats_r,
     &  lonsperlar)

!-- xv
      call nemsio_readrecv(gfile,'xv','sfc',1,buff1,iret=iret)
      call split2d_rst(buff1,nst_fld%xv,fieldsize,global_lats_r,
     &  lonsperlar)

!-- xz
      call nemsio_readrecv(gfile,'xz','sfc',1,buff1,iret=iret)
      call split2d_rst(buff1,nst_fld%xz,fieldsize,global_lats_r,
     &  lonsperlar)

!-- zm
      call nemsio_readrecv(gfile,'zm','sfc',1,buff1,iret=iret)
      call split2d_rst(buff1,nst_fld%zm,fieldsize,global_lats_r,
     &  lonsperlar)

!-- xtts
!      call nemsio_readrecv(gfile,'xtts','sfc',1,buff1,iret=iret)
!     call split2d_rst(buff1,nst_fld%xtts,fieldsize,global_lats_r,
!    &  lonsperlar)

!-- xzts
!     call nemsio_readrecv(gfile,'xzts','sfc',1,buff1,iret=iret)
!     call split2d_rst(buff1,nst_fld%xzts,fieldsize,global_lats_r,
!    &  lonsperlar)

!-- dt_cool
      call nemsio_readrecv(gfile,'dtcool','sfc',1,buff1,iret=iret)
      call split2d_rst(buff1,nst_fld%dt_cool,fieldsize,global_lats_r,
     &  lonsperlar)

!-- z_c
      call nemsio_readrecv(gfile,'zc','sfc',1,buff1,iret=iret)
      call split2d_rst(buff1,nst_fld%z_c,fieldsize,global_lats_r,
     &  lonsperlar)

!-- c_0
      call nemsio_readrecv(gfile,'c0','sfc',1,buff1,iret=iret)
      call split2d_rst(buff1,nst_fld%c_0,fieldsize,global_lats_r,
     &  lonsperlar)

!-- c_d
      call nemsio_readrecv(gfile,'cd','sfc',1,buff1,iret=iret)
      call split2d_rst(buff1,nst_fld%c_d,fieldsize,global_lats_r,
     &  lonsperlar)

!-- xt
      call nemsio_readrecv(gfile,'w0','sfc',1,buff1,iret=iret)
      call split2d_rst(buff1,nst_fld%w_0,fieldsize,global_lats_r,
     &  lonsperlar)

!-- w_d
      call nemsio_readrecv(gfile,'wd','sfc',1,buff1,iret=iret)
      call split2d_rst(buff1,nst_fld%w_d,fieldsize,global_lats_r,
     &  lonsperlar)

!-- d_conv
      call nemsio_readrecv(gfile,'dconv','sfc',1,buff1,iret=iret)
      call split2d_rst(buff1,nst_fld%xt,fieldsize,global_lats_r,
     &  lonsperlar)

!-- ifd
      call nemsio_readrecv(gfile,'ifd','sfc',1,buff1,iret=iret)
      call split2d_rst(buff1,nst_fld%ifd,fieldsize,global_lats_r,
     &  lonsperlar)

!-- tref
      call nemsio_readrecv(gfile,'tref','sfc',1,buff1,iret=iret)
      call split2d_rst(buff1,nst_fld%tref,fieldsize,global_lats_r,
     &  lonsperlar)

!-- Qrain
      call nemsio_readrecv(gfile,'Qrain','sfc',1,buff1,iret=iret)
      call split2d_rst(buff1,nst_fld%Qrain,fieldsize,global_lats_r,
     &  lonsperlar)
!       print *,'in read_nst_r,qrain=',nst_fld%Qrain(1:3,1:3)

!
      call nemsio_close(gfile)
      call nemsio_finalize()
!
      t2=timef()
!      print *,'end of read_nst_r time ',t2-t1,t1,t2
!
      RETURN

      STOP
      END
!
!***********************************************************************
!
