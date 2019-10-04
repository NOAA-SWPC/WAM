!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
    subroutine post_run_gfs(wrt_int_state,mype,mpicomp,lead_write,      &
               mygridtype,mymaptype,mynsoil,nbdl,mynfhr,mynfmin)
!
!  revision history:
!     Jun 2010    J. Wang      Initial code
!     Jan 2012    J. Wang      Add aerosol fields
!     Feb 2012    J. Wang      setvar_aerfile is reset after post
!     Jan 2013    Sarah Lu     EL_MYJ changed to EL_PBL to be consistent
!                              with nceppost upgrade
!     28May2013   Sarah Lu     Specify iostatusD3D
!     07Nov2014   S. Moorthi   Threading, optimization, bug-fix
!     09Oct2015   S. Moorthi - adding imp_physics argument to MICROINIT
!
!-----------------------------------------------------------------------
!*** run post on quilt
!-----------------------------------------------------------------------
!
      use MODULE_WRITE_INTERNAL_STATE_GFS
      use CTLBLK_mod, only : komax,ifhr,ifmin,MODELNAME,datapd,fld_info, &
                             npset,grib,gocart_on,imp_physics
      use grib2_module, only : gribit2,num_pset,nrecout,first_grbtbl
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      type(WRITE_INTERNAL_STATE_GFS),intent(in) :: wrt_int_state
      integer,intent(in)                        :: mype
      integer,intent(in)                        :: mpicomp
      integer,intent(in)                        :: lead_write
      character(1),intent(in)                   :: mygridtype
      integer,intent(in)                        :: mymaptype
      integer,intent(in)                        :: mynsoil
      integer,intent(in)                        :: nbdl
      integer,intent(in)                        :: mynfhr
      integer,intent(in)                        :: mynfmin
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      integer,PARAMETER :: NFD=11,NBND=6
!
      integer N,NWTPG,IEOF,LCNTRL,ierr,jts,jte
      integer,allocatable  :: jstagrp(:),jendgrp(:)
      integer,save         :: kpo,kth,kpv
      real,dimension(komax),save :: po, th, pv 
      logical,save   :: LOG_POSTALCT=.false.
      logical,save   :: setvar_gridfile=.false.,setvar_flxfile=.false.
      logical,save   :: setvar_aerfile=.false.
      logical        :: Log_runpost
      character(255) :: post_fname*255

      integer,save   :: iostatusD3D=-1
!
!      print *,'in post_run start'
!-----------------------------------------------------------------------
!*** set up dimensions
!-----------------------------------------------------------------------
!
      MODELNAME = "GFS"
      NWTPG     = wrt_int_state%WRITE_TASKS_PER_GROUP
      JTS       = wrt_int_state%JSTART_WRITE(MYPE-lead_write+1)      !<-- Starting J of this write task's subsection
      JTE       = wrt_int_state%JEND_WRITE  (MYPE-lead_write+1)      !<-- Ending J of this write task's subsection

!      print *,'in post_run,jts=',jts,'jte=',jte,'nwtpg=',nwtpg, &
!        'log_postalct=',log_postalct,'jts=',jts,'jte=',jte,'nbdl=',nbdl
!
!-----------------------------------------------------------------------
!*** set up fields ro run post
!-----------------------------------------------------------------------
!
      IF(.not.LOG_POSTALCT) THEN
!
        ALLOCATE(JSTAGRP(NWTPG),JENDGRP(NWTPG))
!
        DO N=0,NWTPG-1
          JSTAGRP(N+1) = wrt_int_state%JSTART_WRITE(N+1)
          JENDGRP(N+1) = wrt_int_state%JEND_WRITE  (N+1)
        ENDDO
!      print *,'in post_run,jstagrp=',jstagrp,'jendgrp=',jendgrp

!-----------------------------------------------------------------------
!*** read namelist for pv,th,po
!-----------------------------------------------------------------------
!
         call read_postnmlt(kpo,kth,kpv,po,th,pv,wrt_int_state%nlunit, &
                            wrt_int_state%post_namelist)
!     write(0,*)'in post_run,aft nmlst po gocart_on=',gocart_on
!
!-----------------------------------------------------------------------
!*** allocate post variables
!-----------------------------------------------------------------------
!
        call post_alctvars(wrt_int_state%im,wrt_int_state%jm,        &
               wrt_int_state%lm,MYPE,wrt_int_state%WRITE_TASKS_PER_GROUP,   &
               mpicomp,mygridtype,mymaptype,wrt_int_state%post_gribversion,&
               MYNSOIL,         &
               LEAD_WRITE,JTS,JTE,JSTAGRP,JENDGRP)
!      print *,'in post_run,aft post_alctvars'
!
!-----------------------------------------------------------------------
!*** read namelist for pv,th,po
!-----------------------------------------------------------------------
!
!        call read_postnmlt(kpo,kth,kpv,po,th,pv)
!     write(0,*)'in post_run,aft nmlst po'
!
!-----------------------------------------------------------------------
!*** read namelist for pv,th,po
!-----------------------------------------------------------------------
!
        LOG_POSTALCT = .true.
        first_grbtbl = .true.
!
      ENDIF
!       
!-----------------------------------------------------------------------
!*** fill post variables with values from forecast results
!-----------------------------------------------------------------------
!
      ifhr  = mynfhr
      ifmin = mynfmin
      call set_postvars_gfs(wrt_int_state,mpicomp,nbdl,JTS,JTE)
      if(trim(wrt_int_state%FILENAME_BASE(NBDL))=='SIG.F') setvar_gridfile=.true.
      if(trim(wrt_int_state%FILENAME_BASE(NBDL))=='FLX.F') setvar_flxfile=.true.
      if(trim(wrt_int_state%FILENAME_BASE(NBDL))=='AER.F') setvar_aerfile=.true.

!       print *,'af set_postvars,setvar_gridfile=',setvar_gridfile,  &
!        'setvar_flxfile=',setvar_flxfile,'ifmin=',ifmin,'ifhr=',ifhr, &
!       'setvar_aerfile=',setvar_aerfile,'gocart_aer2post=',wrt_int_state%gocart_aer2post
!
      if(.not.wrt_int_state%gocart_aer2post) then
        Log_runpost = setvar_gridfile.and.setvar_flxfile
      else
        Log_runpost = setvar_gridfile.and.setvar_flxfile.and.setvar_aerfile
      endif
!
      if(Log_runpost) then
!  
        if(imp_physics==5 .or. imp_physics==85 .or. imp_physics==95) then
          call MICROINIT(imp_physics)
        endif
!
        if(grib=="grib2" .and. first_grbtbl) then
          call READ_xml()
        endif
!
        IEOF  = 0
        npset = 0
        do while( IEOF == 0)
!
          if(grib == "grib1") then
            CALL READCNTRL(kth,IEOF)
           else if(grib == "grib2") then
            npset = npset + 1
            CALL SET_OUTFLDS(kth,kpv,pv)
            if(allocated(datapd))deallocate(datapd)
            allocate(datapd(wrt_int_state%im(1),jte-jts+1,nrecout+10))
            datapd = 0.
            call get_postfilename(post_fname)
          endif
!
!         if ( IEOF.eq.0) CALL PROCESS(KTH,KPV,TH(1:KTH),PV(1:KPV))
          if ( IEOF == 0) CALL PROCESS(KTH,KPV,TH(1:KTH),PV(1:KPV),iostatusD3D)
!
          if(grib == "grib2") then
            call mpi_barrier(mpicomp,ierr)
            call gribit2(post_fname)
            if(allocated(datapd))deallocate(datapd)
            if(allocated(fld_info))deallocate(fld_info)
            if(npset >= num_pset) exit
          endif
!
        enddo
!
        if(grib == "grib1") then
          LCNTRL = 14
          rewind(LCNTRL)
        endif
!
        setvar_gridfile = .false.
        setvar_flxfile  = .false.
        setvar_aerfile  = .false.
!
      endif
!          write(0,*)'after readcntrl and process'

    end subroutine post_run_gfs
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
    subroutine set_postvars_gfs(wrt_int_state,mpicomp,nbdl,jts,jte)
!
!  revision history:
!     Jun 2010    J. Wang      Initial code
!     Jan 2012    Wang/Lu      set aerosol post variables from aer file
!     Feb 2012    Wang         fix grid information for grib2
!
!-----------------------------------------------------------------------
!*** set up post fields from nmint_state
!-----------------------------------------------------------------------
!
      use vrbls4d
      use vrbls3d
      use vrbls2d
      use soil
      use masks
      use ctlblk_mod
      use params_mod
      use gridspec_mod
      use lookup_mod
      use physcons
      use rqstfld_mod
!
      use MODULE_WRITE_INTERNAL_STATE_GFS
!
!-----------------------------------------------------------------------
!
      implicit none
!
      include 'mpif.h'
!
!-----------------------------------------------------------------------
!
      type(WRITE_INTERNAL_STATE_GFS),intent(in) :: wrt_int_state
      integer,intent(in)                        :: mpicomp,nbdl,jts,jte
!
!-----------------------------------------------------------------------
!
      integer I,ii,J,jj,L,LL,LEV,K,N,N1,N2,NPOSN_1,NPOSN_2, LENGTH
      integer NPOS_START,NPOS_END,indx_num,nfield,ierr,iret,igdout
      integer INDEXP,INDX_2D,INDX_2D2,INDX_2D3,iostatusD3D,itr
      integer RINC(5)
      character(NAME_MAXSTR) :: NAME,NAMEWL
      CHARACTER(6)           :: FMT='(I2.2)'
      CHARACTER(6)           :: MODEL_LEVEL
      CHARACTER(3)           :: CK
      CHARACTER(16)          :: layer,VarAerName,VarAerLevtyp
      REAL :: FACT,TLMH,RADI,TMP,ES,TV,RHOAIR,tem
      REAL,dimension(:,:,:),allocatable :: buf3d,fi
      REAL,dimension(:,:),  allocatable :: dummy,dummy2,buf
      REAL,dimension(:),    allocatable  :: slat,qstl
      real,external::FPVSNEW
!     real,   allocatable :: d2d(:,:), u2d(:,:), v2d(:,:), omga2d(:,:)
      real,   allocatable :: p2d(:,:), t2d(:,:), q2d(:,:),  qs2d(:,:),  &
                             cw2d(:,:), cfr2d(:,:)
!     real*8, allocatable :: pm2d(:,:), pi2d(:,:)
!
!-----------------------------------------------------------------------
!***  INTEGER SCALAR/1D HISTORY VARIABLES
!-----------------------------------------------------------------------
!
      maptype     = 4      ! default gaussian grid
      gridtype    = 'A'
      imp_physics = 99 !set GFS mp physics to 99 for Zhao scheme
      iostatusD3D = -1

!      print*,'MP_PHYSICS= ',imp_physics,'nbdl=',nbdl
!
! nems gfs has zhour defined 
      tprec   = 1.0*ifhr-wrt_int_state%ZHOUR
      tclod   = tprec
      trdlw   = tprec
      trdsw   = tprec
      tsrfc   = tprec
      tmaxmin = tprec
      td3d    = tprec
!

!      write(6,*) 'maptype and gridtype is ', maptype,gridtype
!
!-----------------------------------------------------------------------
! some fields are only in grid file, not in flx
!-----------------------------------------------------------------------
!
      allocate(buf(im,jsta:jend))

      inifields_if: if(trim(wrt_int_state%FILENAME_BASE(NBDL))=='SIG.F') then
!
        allocate(buf3d(im,jsta:jend,lm),fi(im,jm,2))
!
!set gaussian grid on write lead task:
        allocate(dummy(im,jm),dummy2(im,jm))
        if(me == 0) then
          allocate(slat(jm))
          call splat(4,jm,slat)
          radi = 180.0d0/(4.d0*atan(1.d0))
          tmp = 360.0 / im
!$omp parallel do private(i,j,tem)
          do  j=1,jm
            tem = asin(slat(j)) * radi
            do  i=1,im
              dummy(i,j)  = tem
              dummy2(i,j) = tmp*(i-1)
            enddo
          enddo
          deallocate(slat)
        endif
!scatter to other write tasks
        call mpi_scatterv(dummy(1,1),icnt,idsp,mpi_real                   &
                         ,gdlat(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,ierr)
        call mpi_scatterv(dummy2(1,1),icnt,idsp,mpi_real                   &
                         ,gdlon(1,jsta),icnt(me),mpi_real,0,MPI_COMM_COMP,ierr)
!       print *,'in ste var,gdlat=',maxval(gdlat(:,jsta:jend)),minval(gdlat(:,jsta:jend))
!       print *,'in ste var,gdlon=',maxval(gdlon(:,jsta:jend)),minval(gdlon(:,jsta:jend))
!
        CALL EXCH(gdlat(1,JSTA_2L))

!$omp parallel do private(i,j)
      do j = jsta, jend_m
        do i = 1, im-1
          DX( i,j) = ERAD*COS(GDLAT(I,J)*DTR)*(GDLON(I+1,J)-GDLON(I,J))*DTR
          DY(i,j)  = ERAD*(GDLAT(I,J)-GDLAT(I,J+1))*DTR  ! like A*DPH
        end do
      end do

!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=1,im
          F(I,J) = 1.454441e-4*sin(gdlat(i,j)*DTR)   ! 2*omeg*sin(phi)
        end do
      end do
!
! GFS does not output PD
!     pd=spval
!     pdtop=spval
!     pt=0.

!     pdtop = spval
!     pt    = 0.
      pt    = 10000.          ! this is for 100 hPa added by Moorthi
      pd    = spval           ! GFS does not output PD

      do j=jsta,jend
        do i=1,im
          pint(i,j,1) = PT
        end do
      end do
!jw GFS set omga
      omga = spval
!
! cloud water and ice mixing ratio  for zhao scheme
! need to look up old eta post to derive cloud water/ice from cwm
! Zhao scheme does not produce suspended rain and snow
      qqr = 0.
      qqs = 0.
      qqw = 0.
      qqi = 0.
      qqg = spval

      ttnd = spval
! GFS does not have surface specific humidity
      QS = SPVAL

! GFS does not have inst sensible heat flux
      twbs = SPVAL

! GFS does not have inst latent heat flux
      qwbs = SPVAL

!  GFS does not have time step and physics time step, make up ones since they
! are not really used anyway
      NPHS = 2.
      DT   = 80.
      DTQ2 = DT * NPHS  !MEB need to get physics DT
!jw      TSPH = 3600./DT   !MEB need to get DT
! All GFS time-averaged quantities are in 6 hour bucket
!      TPREC = 6.0
! GFS does not have accumulated total, gridscale, and convective precip, will use inst precip to
! derive in SURFCE.f

! GFS does not have similated precip
      lspa = spval

! GFS does not have convective cloud efficiency
      CLDEFI = SPVAL

! GFS does not have 10 m theta
      TH10 = SPVAL

! GFS does not have 10 m humidity
      Q10 = SPVAL

! GFS does not have TKE because it uses MRF scheme
      Q2 = SPVAL

! GFS does not have surface exchange coeff

! GFS does not have snow free albedo
      ALBASE = SPVAL

! GFS probably does not use zenith angle
      Czen   = spval
      CZMEAN = SPVAL

! GFS does not have inst surface outgoing longwave
      radot  = spval
!
! GFS does not have inst cloud fraction for high, middle, and low cloud
      cfr    = spval
      cfrach = spval
      cfracl = spval
      cfracm = spval
!
      if (gocart_on) then
        dust   = spval
      endif
!
      ICING_GFIP = spval
!
      sr = 0.
!
! GFS does not have inst ground heat flux
      grnflx = spval
! GFS doesn not yet output soil layer thickness, assign SLDPTH to be the same as nam

      SLDPTH(1) = 0.10
      SLDPTH(2) = 0.3
      SLDPTH(3) = 0.6
      SLDPTH(4) = 1.0
! GFS does not output time averaged convective and strat cloud fraction, set acfrcv to spval, n
! cfrcv to 1
      acfrcv = spval
      ncfrcv = 1.0
! GFS does not output time averaged cloud fraction, set acfrst to spval, ncfrst to 1
      acfrst = spval
      ncfrst = 1.0

! GFS does not have storm runoff
      ssroff = spval

! GFS does not have UNDERGROUND RUNOFF
      bgroff = spval
! GFS incoming sfc longwave has been averaged over 6 hr bucket, set ARDLW to 1
      ardlw  = 1.0
!     trdlw  = 6.0

! GFS does not have inst incoming sfc longwave
      rlwin = spval

! GFS does not have inst model top outgoing longwave
      rlwtoa = spval
! GFS does not have inst incoming sfc shortwave
      rswin  = spval

! GFS does not have inst incoming clear sky sfc shortwave
      rswinc = spval

! GFS does not have inst outgoing sfc shortwave
      rswout = spval

! GFS incoming sfc longwave has been averaged, set ARDLW to 1
      ardsw = 1.0
!     trdsw = 6.0
! GFS surface flux has been averaged, set  ASRFC to 1
      asrfc = 1.0
! GFS does not have snow phase change heat flux
      snopcx = spval
! GFS does not use total momentum flux
      sfcuvx = spval
! GFS does not have temperature tendency due to long wave radiation
      rlwtt  = spval

! GFS does not have temperature tendency due to long wave radiation
      rswtt  = spval

! GFS does not have temperature tendency due to latent heating from convection
      tcucn  = spval
      tcucns = spval

! set avrain to 1
      avrain = 1.0
      avcnvc = 1.0
      theat  = 6.0 ! just in case GFS decides to output T tendency

! GFS does not have temperature tendency due to latent heating from grid scale
      train  = spval

! GFS does not have soil moisture availability
      smstav = spval

! GFS does not have total soil moisture
      smstot = spval
! GFS does not have accumulated surface evaporation
      sfcevp = spval

! GFS does not have surface exchange coeefficient
      sfcexc = spval

! GFS does not have averaged accumulated snow
      acsnow = spval

! GFS does not have snow melt
      acsnom = spval

! GFS does not have sst????
      sst = spval
! GFS does not have mixing length
!     el_myj = spval
      el_pbl = spval

! GFS does not output exchange coefficient
      exch_h = spval


! GFS does not output humidity at roughness length
      qz0 = spval

! GFS does not output u at roughness length
      uz0 = spval

! GFS does not output humidity at roughness length
      vz0        = spval
      MAXRHSHLTR = SPVAL
      MINRHSHLTR = SPVAL
! GFS does not have deep convective cloud top and bottom fields
      htop       = spval
      HTOPD      = SPVAL
      HBOTD      = SPVAL
      HTOPS      = SPVAL
      HBOTS      = SPVAL
      CUPPT      = SPVAL
!
      F_ice      = spval
      F_rain     = spval
      F_RimeF    = spval
      dbz        = spval
!
      IDAT(5) = 0
      N2      = 0
      DO N=1,wrt_int_state%KOUNT_I1D(NBDL)                                    !<-- Loop through all scalar/1D integer data
!
        NPOSN_1 = (N-1)*NAME_MAXSTR+1
        NPOSN_2 = N*NAME_MAXSTR
        NAME    = wrt_int_state%NAMES_I1D_STRING(NBDL)(NPOSN_1:NPOSN_2)               !<-- The variable's name
        LENGTH  = wrt_int_state%LENGTH_DATA_I1D(N,NBDL)                            !<-- The variable's length in words
!        print *,'in set_var,I1d,name=',trim(name),'length=',length
!
        IF(LENGTH == 1)THEN
          N2 = N2 + 1
        ELSE
!
          DO N1=1,LENGTH
            N2 = N2 + 1
            IF(TRIM(NAME) == 'IDAT') THEN
              if(N1==1) IDAT(N1) = wrt_int_state%ALL_DATA_I1D(N2,NBDL)
              IF(N1==2) IDAT(N1) = wrt_int_state%ALL_DATA_I1D(N2,NBDL)
              IF(N1==3) IDAT(N1) = wrt_int_state%ALL_DATA_I1D(N2,NBDL)
              IF(N1==4) IDAT(N1) = wrt_int_state%ALL_DATA_I1D(N2,NBDL)
            ENDIF
          ENDDO
!
       ENDIF
!
     ENDDO
     SDAT(1) = IDAT(2)
     SDAT(2) = IDAT(3)
     SDAT(3) = IDAT(4)
     IHRST   = IDAT(1)
     IMIN    = IDAT(5)
!
!     print *,'in set_var,nbdl=',nbdl,'KOUNT_I1D=',wrt_int_state%KOUNT_I1D(NBDL),'idat=',idat
!
     endif inifields_if
!
!-----------------------------------------------------------------------------
! get post fields 
!-----------------------------------------------------------------------------
!
      N=0
      DO nfield=1,wrt_int_state%KOUNT_R2D(NBDL)
!
       NPOS_START=(NFIELD-1)*NAME_MAXSTR+1
       NPOS_END=NFIELD*NAME_MAXSTR
!       print *,'NPOS_START=',NPOS_START,'NPOS_END=',NPOS_END
       NAMEWL=wrt_int_state%NAMES_R2D_STRING(NBDL)(NPOS_START:NPOS_END)
       INDX_2D=index(NAMEWL,"_2D")
!       print *,'set_var,nfield=',nfield,'NAMEWL=',trim(NAMEWL),' INDX_2D=',INDX_2D,&
!        'NBDL=',NBDL,'KOUNT_R2D=',wrt_int_state%KOUNT_R2D(NBDL),'NPOS_START=',NPOS_START, &
!        'NPOS_END=',NPOS_END

!for 3D fields
       if(INDX_2D>0) then
         INDX_2D2=INDEX(NAMEWL,"_")
         if(INDX_2D2>0) NAME=NAMEWL(1:INDX_2D2-1)
         MODEL_LEVEL=NAMEWL(INDX_2D-2:INDX_2D-1)
         lev=(ICHAR(MODEL_LEVEL(1:1))-48)*10+ICHAR(MODEL_LEVEL(2:2))-48
         ll=lm-lev+1

!       print *,'set_var,3D vars, NAME=',trim(NAME),'nfield=',nfield,    &
!         'model_level=',MODEL_LEVEL,'lev=',lev

! model level T
         if(trim(name)=='tmp') then
!$omp parallel do private(j)
           do j=jsta,jend
             t(:,j,ll)=wrt_int_state%WRITE_SUBSET_R(:,j,nfield)
           enddo
!           print *,'in set_var,nbdl=',nbdl,'tmp=',maxval(t(:,jsta:jend,ll)), &
!               minval(t(:,jsta:jend,ll)),'ll=',ll
         endif
! model level q
          if(trim(name)=='spfh') then
!$omp parallel do private(j)
            do j=jsta,jend
              q(:,j,ll)=wrt_int_state%WRITE_SUBSET_R(:,j,nfield)
            enddo
!           print *,'in set_var,nbdl=',nbdl,'q=',maxval(q(:,jsta:jend,ll)), &
!               minval(q(:,jsta:jend,ll)),'ll=',ll
          endif
! model level u
          if(trim(name)=='ugrd') then
!$omp parallel do private(j)
            do j=jsta,jend
              uh(:,j,ll)=wrt_int_state%WRITE_SUBSET_R(:,j,nfield)
            enddo
          endif
! model level v
          if(trim(name)=='vgrd') then
!$omp parallel do private(j)
            do j=jsta,jend
              vh(:,j,ll)=wrt_int_state%WRITE_SUBSET_R(:,j,nfield)
            enddo
          endif
! model level pressure
          if(trim(name)=='pres') then
!$omp parallel do private(j)
            do j=jsta,jend
              pmid(:,j,ll)=wrt_int_state%WRITE_SUBSET_R(:,j,nfield)
            enddo
          endif
! GFS is on A grid and does not need PMIDV

! dp
          if(trim(name)=='dpres') then
!$omp parallel do private(j)
            do j=jsta,jend
              buf3d(:,j,ll)=wrt_int_state%WRITE_SUBSET_R(:,j,nfield)
            enddo
          endif
! ozone mixing ratio
          if(trim(name)=='o3mr') then
!$omp parallel do private(j)
            do j=jsta,jend
              o3(:,j,ll)=wrt_int_state%WRITE_SUBSET_R(:,j,nfield)
            enddo
          endif
!cloud water
          if(trim(name)=='clwmr') then
!$omp parallel do private(i,j)
            do j=jsta,jend
              do i=1,im
                cwm(i,j,ll) = wrt_int_state%WRITE_SUBSET_R(i,j,nfield)
                if (cwm(i,j,ll) < 1.0e-8) cwm(i,j,ll) = 0.0
                if(t(i,j,ll) < (TFRZ-15.) )then ! dividing cloud water from ice
                  qqi(i,j,ll) = cwm(i,j,ll)
                else
                  qqw(i,j,ll) = cwm(i,j,ll)
                end if
              enddo
            enddo
!         if (j.eq.jm/2 .and. mod(i,50).eq.0)
!     +   print*,'sample ',trim(VarName), ' after scatter= '
!     +   ,i,j,ll,cwm(i,j,ll)
          endif
! GFS does output omeg now
          if(trim(name)=='vvel') then
!$omp parallel do private(j)
            do j=jsta,jend
              omga(:,j,ll) = wrt_int_state%WRITE_SUBSET_R(:,j,nfield)
            enddo
!            print *,'in set_postvar,ll=',ll,'omga=',maxval(omga(1:im,jsta:jend,ll)), &
!              minval(omga(1:im,jsta:jend,ll))
          endif
          if (gocart_on) then
! GFS output dust in nemsio
          if(trim(name)=='du001') then
!$omp parallel do private(j)
            do j=jsta,jend
              dust(:,j,ll,1) = wrt_int_state%WRITE_SUBSET_R(:,j,nfield)
            enddo
!             print *,'for du001,setvar_post,dust1=',maxval(dust(:,jsta:jend,ll,1)), &
!             minval(dust(:,jsta:jend,ll,1)),'ll=',ll,'lev=',lev
!            print *,'ll=',ll,'dust(2,102,ll,1)=',dust(2,102,ll,1),'value=', &
!             maxval(dust(:,jsta:jend,ll,1)),minval(dust(:,jsta:jend,ll,1))
          endif
          if(trim(name)=='du002') then
!$omp parallel do private(j)
            do j=jsta,jend
              dust(:,j,ll,2) = wrt_int_state%WRITE_SUBSET_R(:,j,nfield)
            enddo
          endif
          if(trim(name)=='du003') then
!$omp parallel do private(j)
            do j=jsta,jend
              dust(:,j,ll,3) = wrt_int_state%WRITE_SUBSET_R(:,j,nfield)
            enddo
          endif
          if(trim(name)=='du004') then
!$omp parallel do private(j)
            do j=jsta,jend
              dust(:,j,ll,4) = wrt_int_state%WRITE_SUBSET_R(:,j,nfield)
            enddo
          endif
          if(trim(name)=='du005') then
!$omp parallel do private(j)
            do j=jsta,jend
              dust(:,j,ll,5) = wrt_int_state%WRITE_SUBSET_R(:,j,nfield)
            enddo
          endif
          endif    ! if gocart_in

!2D fields
      else
        INDX_2D2=INDEX(NAMEWL,"_")
        if(INDX_2D2>0) then
          NAME = trim(NAMEWL(1:INDX_2D2-1))
        else
          NAME = TRIM(NAMEWL)
        endif
        INDX_2D3 = INDEX(NAMEWL,"_",BACK=.true.)
        if(INDX_2D3>0) then
          layer = NAMEWL(INDX_2D3+1:)
        else
          layer = 'sfc'
        endiF
        itr=0
        if(INDEX(NAMEWL,"_ave") > 0) then
          itr = 3
        elseif(INDEX(NAMEWL,"_acc") > 0) then
          itr = 4
        elseif(INDEX(NAMEWL,"_win") > 0) then
          itr = 2
        endif

       print *,'set_var,2D vars, NAME=',trim(NAME),'layer=',trim(layer),'nfield=',nfield,'itr=',itr
!
       if(trim(wrt_int_state%FILENAME_BASE(NBDL))=='SIG.F' ) then
! Surface pressure  using nemsio
         if(trim(name)=='pres'.and.trim(layer)=='sfc') then
!$omp parallel do private(j)
           do j=jsta,jend
             pint(:,j,lp1)=wrt_int_state%WRITE_SUBSET_R(:,j,nfield)
           enddo
!           print *,'in set_postvar,psfc=',maxval(pint(1:im,jsta:jend,lp1)), &
!           minval(pint(1:im,jsta:jend,lp1))
!$omp parallel do private(i,j)
           do j=jsta,jend
             do i=1,im
!jw real4      ALPINT(I,J,LP1) = ALOG(PINT(I,J,LP1))
               ALPINT(I,J,LP1) = LOG(PINT(I,J,LP1))
!              if(i==ii .and. j==jj)print*,'sample PSFC',i,j,pint(i,j,lp1)
             end do
           end do
         endif
! Terrain height * G
         if(trim(name)=='hgt'.and.trim(layer)=='sfc') then
!$omp parallel do private(i,j)
           do j=jsta,jend
             do i=1,im
               fis(i,j) = wrt_int_state%WRITE_SUBSET_R(i,j,nfield)
               if (fis(i,j) < spval) then
                 fis(i,j) = fis(i,j) * con_g
               else
                 fis(i,j) = spval
               endif
             enddo
           enddo
!           print *,'in set_var,fis=',maxval(fis(:,jsta:jend)),minval(fis(:,jsta:jend))
         endif
       endif
!land
       if(trim(name)=='land') then
!$omp parallel do private(i,j)
         do j=jsta,jend
           do i=1,im
             sm(i,j) = wrt_int_state%WRITE_SUBSET_R(i,j,nfield)
             if (sm(i,J) < spval) sm(i,j) = 1.0 - sm(i,j)    ! convert to sea mask
             if (sm(i,J) > spval) sm(i,j) = spval
           enddo
         enddo
       endif
!icec
       if(trim(name)=='icec') then
!$omp parallel do private(i,j)
         do j=jsta,jend
           do i=1,im
             sice(i,j) = wrt_int_state%WRITE_SUBSET_R(i,j,nfield)
             if (sm(i,j) < spval .and. sm(i,j) == 0.0) sice(i,j) = 0.0 !specify sea ice=0 at land
           enddo
         enddo
       endif
! PBL height
       if(trim(name)=='hpbl') then
!$omp parallel do private(j)
         do j=jsta,jend
           pblh(:,j) = wrt_int_state%WRITE_SUBSET_R(:,j,nfield)
         enddo
       endif
! frictional velocity using nemsio
       if(trim(name)=='fricv') then
!$omp parallel do private(j)
         do j=jsta,jend
           ustar(:,j) = wrt_int_state%WRITE_SUBSET_R(:,j,nfield)
         enddo
       endif
! roughness length
       if(trim(name)=='sfcr') then
!$omp parallel do private(i,j)
         do j=jsta,jend
           do i=1,im
             z0(i,j) = wrt_int_state%WRITE_SUBSET_R(i,j,nfield)
             if (z0(i,j) < 1.0e-5) z0(i,j) = 0.0
           enddo
         enddo
       endif
! surface potential T
       if(trim(name)=='tmp'.and.trim(layer)=='sfc') then
!$omp parallel do private(i,j)
         do j=jsta,jend
           do i=1,im
             ths(i,j) = wrt_int_state%WRITE_SUBSET_R(i,j,nfield)
             if (ths(i,j) < spval) then
               ths(i,j) = ths(i,j) * (p1000/pint(i,j,lp1))**capa ! convert to THS
             else
               ths(i,j) = spval
             endif
             thz0(i,j) = ths(i,j)  ! GFS does not have THZ0, use THS to substitute
           enddo
         enddo
!         print *,'in set_var,tsfc=',maxval(ths(:,jsta:jend)),minval(ths(:,jsta:jend))
! GFS does not have THZ0, use THS to substitute
       endif

! convective precip in m per physics time step
       if(trim(name)=='cprat'.and.trim(layer)=='sfc') then
!$omp parallel do private(i,j)
         do j=jsta,jend
           do i=1,im
             avgcprate(i,j) = wrt_int_state%WRITE_SUBSET_R(i,j,nfield)
             if (avgcprate(i,j) < spval) then
               avgcprate(i,j) = avgcprate(i,j) * dtq2/1000. ! convert to m
             else
               avgcprate(i,j) = spval
             endif
             cprate(i,j) = avgcprate(i,j)
           enddo
         enddo
       endif
! precip rate in m per physics time step
       if(trim(name)=='prate'.and.trim(layer)=='sfc') then
!$omp parallel do private(i,j)
         do j=jsta,jend
           do i=1,im
             avgprec(i,j) = wrt_int_state%WRITE_SUBSET_R(i,j,nfield)
             if (avgprec(i,j) < spval) then
               avgprec(i,j) = avgprec(i,j) * dtq2/1000. ! convert to m
             else
               avgprec(i,j) = spval
             endif
             prec(i,j) = avgprec(i,j)
           enddo
         enddo
       endif
! inst snow water eqivalent
       if(trim(name)=='weasd'.and.trim(layer)=='sfc') then
!$omp parallel do private(j)
         do j=jsta,jend
           sno(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,nfield)
         enddo
!         print *,'in set_var,weasd=',maxval(sno(:,jsta:jend)),minval(sno(:,jsta:jend))
       endif
! snow depth in mm
       if(trim(name)=='snod'.and.trim(layer)=='sfc') then
!$omp parallel do private(i,j)
         do j=jsta,jend
           do i=1,im
             si(i,j) = wrt_int_state%WRITE_SUBSET_R(i,j,nfield)
             if (si(i,j) < spval) then
               si(i,j) = si(i,j) * 1000.0 ! convert to mm
             else
               si(i,j) = spval
             endif
           enddo
         enddo
       endif
! 2m T
       if(trim(name)=='tmp'.and.trim(layer)=='2 m above gnd') then
!$omp parallel do private(j)
         do j=jsta,jend
           tshltr(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,nfield)
         enddo
! GFS does not have 2m pres, estimate it, also convert t to theta
!$omp parallel do private(i,j)
         Do j=jsta,jend
           Do i=1,im
            PSHLTR(I,J) = pint(I,J,lm+1)*EXP(-0.068283/tshltr(i,j))
            tshltr(i,j) =  tshltr(i,j)*(p1000/PSHLTR(I,J))**CAPA ! convert to theta
!          if (j.eq.jm/2 .and. mod(i,50).eq.0)
!     +   print*,'sample 2m T and P after scatter= '
!     +   ,i,j,tshltr(i,j),pshltr(i,j)
          end do
        end do
!         print *,'in set_var,t2m=',maxval(tshltr(:,jsta:jend)),minval(tshltr(:,jsta:jend))
       endif
! 2m specific humidity
       if(trim(name)=='spfh'.and.trim(layer)=='2 m above gnd') then
!$omp parallel do private(j)
         do j=jsta,jend
           qshltr(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,nfield)
         enddo
       endif
! mid day avg albedo in fraction
       if(trim(name)=='albdo'.and.trim(layer)=='sfc') then
!$omp parallel do private(i,j)
         do j=jsta,jend
           do i=1,im
             avgalbedo(i,j) = wrt_int_state%WRITE_SUBSET_R(i,j,nfield)
             if (avgalbedo(i,j) < spval) then
               avgalbedo(i,j) = avgalbedo(i,j) * 0.01     ! convert to fraction
             else
                avgalbedo(i,j) = spval
             endif
           enddo
         enddo
       endif
! time averaged column cloud fraction
       if(trim(name)=='tcdc'.and.trim(layer)=='atmos col') then
!$omp parallel do private(i,j)
         do j=jsta,jend
           do i=1,im
             avgtcdc(i,j) = wrt_int_state%WRITE_SUBSET_R(i,j,nfield)
             if (avgtcdc(i,j) < spval) then
               avgtcdc(i,j) = avgtcdc(i,j) * 0.01   ! convert to fraction
             else
               avgtcdc(i,j) = spval
             endif
           enddo
         enddo
       endif
! maximum snow albedo in fraction
       if(trim(name)=='mxsalb'.and.trim(layer)=='sfc') then
!$omp parallel do private(i,j)
         do j=jsta,jend
           do i=1,im
             mxsnal(i,j) = wrt_int_state%WRITE_SUBSET_R(i,j,nfield)
             if (mxsnal(i,j) < spval) then
               mxsnal(i,j) = mxsnal(i,j) *0.01  ! convert to fraction
             else
               mxsnal(i,j) = 0.0
             endif
           enddo
         enddo
       endif
! ave high cloud fraction using nemsio
       if(trim(name)=='tcdc'.and.trim(layer)=='high cld lay') then
!$omp parallel do private(i,j)
         do j=jsta,jend
           do i=1,im
             avgcfrach(i,j) = wrt_int_state%WRITE_SUBSET_R(i,j,nfield)
             if (avgcfrach(i,j) < spval) then
               avgcfrach(i,j) = avgcfrach(i,j) * 0.01    ! convert to fraction
             else
               avgcfrach(i,j) = spval
             endif
           enddo
         enddo
       endif
! ave low cloud fraction using nemsio
       if(trim(name)=='tcdc'.and.trim(layer)=='low cld lay') then
!$omp parallel do private(i,j)
         do j=jsta,jend
           do i=1,im
             avgcfracl(i,j) = wrt_int_state%WRITE_SUBSET_R(i,j,nfield)
             if (avgcfracl(i,j) < spval) then
               avgcfracl(i,j) = avgcfracl(i,j) * 0.01   ! convert to fraction
             else
               avgcfracl(i,j) = spval
             endif
           enddo
         enddo
       endif
! ave middle cloud fraction using nemsio
       if(trim(name)=='tcdc'.and.trim(layer)=='mid cld lay') then
!$omp parallel do private(i,j)
         do j=jsta,jend
           do i=1,im
             avgcfracm(i,j) = wrt_int_state%WRITE_SUBSET_R(i,j,nfield)
             if (avgcfracm(i,j) < spval) then
               avgcfracm(i,j) = avgcfracm(i,j) * 0.01     ! convert to fraction
             else
               avgcfracm(i,j) = spval
             endif
           enddo
         enddo
       endif
! inst convective cloud fraction using nemsio
       if(trim(name)=='tcdc'.and.trim(layer)=='convect-cld laye') then
!$omp parallel do private(i,j)
         do j=jsta,jend
           do i=1,im
             cnvcfr(i,j) = wrt_int_state%WRITE_SUBSET_R(i,j,nfield)
             if (cnvcfr(i,j) < spval) then
               cnvcfr(i,j) = cnvcfr(i,j) * 0.01  ! convert to fraction
             else
               cnvcfr(i,j) = spval
             endif
           enddo
         enddo
       endif
! slope type using nemsio
       if(trim(name)=='sltyp'.and.trim(layer)=='sfc') then
!$omp parallel do private(i,j)
         do j=jsta,jend
           do i=1,im
             buf(i,j) = wrt_int_state%WRITE_SUBSET_R(i,j,nfield)
             if (buf(i,j) < spval) then
               islope(i,j) = nint(buf(i,j))
             else
               islope(i,j) = 99999
!              islope(i,j) = spval
             endif
           enddo
         enddo
!         print *,'islope=',maxval(islope(:,jsta:jend)),minval(islope(:,jsta:jend))
       endif
! plant canopy sfc wtr in m using nemsio
       if(trim(name)=='cnwat'.and.trim(layer)=='sfc') then
!$omp parallel do private(i,j)
         do j=jsta,jend
           do i=1,im
             cmc(i,j) = wrt_int_state%WRITE_SUBSET_R(i,j,nfield)
             if (cmc(i,j) < spval) then
               cmc(i,j) = cmc(i,j) * 0.001 ! convert from kg*m^2 to m
             else
               cmc(i,j) = spval
             endif
           enddo
         enddo
       endif
! vegetation fraction in fraction. using nemsio
       if(trim(name)=='veg'.and.trim(layer)=='sfc') then
!$omp parallel do private(i,j)
         do j=jsta,jend
           do i=1,im
             vegfrc(i,j) = wrt_int_state%WRITE_SUBSET_R(i,j,nfield)
             if (vegfrc(i,j) < spval) then
               vegfrc(i,j) = vegfrc(i,j) *0.01   ! convert to fraction
             else
               vegfrc(i,j) = 0.0
             endif
           enddo
         enddo
       endif
! liquid volumetric soil mpisture in fraction using nemsio
       if(trim(name)=='soill'.and.trim(layer)=='0-10 cm down') then
!$omp parallel do private(i,j)
         do j=jsta,jend
           do i=1,im
             sh2o(i,j,1) = wrt_int_state%WRITE_SUBSET_R(i,j,nfield)
             if (sh2o(i,j,1) >= spval) sh2o(i,j,1) = 0.0
           enddo
         enddo
       endif
       if(trim(name)=='soill'.and.trim(layer)=='10-40 cm down') then
!$omp parallel do private(i,j)
         do j=jsta,jend
           do i=1,im
             sh2o(i,j,2) = wrt_int_state%WRITE_SUBSET_R(i,j,nfield)
             if (sh2o(i,j,2) >= spval) sh2o(i,j,2) = 0.0
           enddo
         enddo
       endif
       if(trim(name)=='soill'.and.trim(layer)=='40-100 cm down') then
!$omp parallel do private(i,j)
         do j=jsta,jend
           do i=1,im
             sh2o(i,j,3) = wrt_int_state%WRITE_SUBSET_R(i,j,nfield)
             if (sh2o(i,j,3) >= spval) sh2o(i,j,3) = 0.0
           enddo
         enddo
       endif
       if(trim(name)=='soill'.and.trim(layer)=='100-200 cm down') then
!$omp parallel do private(i,j)
         do j=jsta,jend
           do i=1,im
             sh2o(i,j,4) = wrt_int_state%WRITE_SUBSET_R(i,j,nfield)
             if (sh2o(i,j,4) >= spval) sh2o(i,j,4) = 0.0
           enddo
         enddo
       endif
! volumetric soil moisture using nemsio
       if(trim(name)=='soilw'.and.trim(layer)=='0-10 cm down') then
!$omp parallel do private(i,j)
         do j=jsta,jend
           do i=1,im
             smc(i,j,1) = wrt_int_state%WRITE_SUBSET_R(i,j,nfield)
             if (smc(i,j,1) >= spval) smc(i,j,1) = 0.0
           enddo
         enddo
       endif
       if(trim(name)=='soilw'.and.trim(layer)=='10-40 cm down') then
!$omp parallel do private(i,j)
         do j=jsta,jend
           do i=1,im
             smc(i,j,2) = wrt_int_state%WRITE_SUBSET_R(i,j,nfield)
             if (smc(i,j,2) >= spval) smc(i,j,2) = 0.0
           enddo
         enddo
       endif
       if(trim(name)=='soilw'.and.trim(layer)=='40-100 cm down') then
!$omp parallel do private(i,j)
         do j=jsta,jend
           do i=1,im
             smc(i,j,3) = wrt_int_state%WRITE_SUBSET_R(i,j,nfield)
             if (smc(i,j,3) >= spval) smc(i,j,3) = 0.0
           enddo
         enddo
       endif
       if(trim(name)=='soilw'.and.trim(layer)=='100-200 cm down') then
!$omp parallel do private(i,j)
         do j=jsta,jend
           do i=1,im
             smc(i,j,4) = wrt_int_state%WRITE_SUBSET_R(i,j,nfield)
             if (smc(i,j,4) >= spval) smc(i,j,4) = 0.0
           enddo
         enddo
       endif
! soil temperature using nemsio
       if(trim(name)=='tmp'.and.trim(layer)=='0-10 cm down') then
!$omp parallel do private(i,j)
         do j=jsta,jend
           do i=1,im
             stc(i,j,1) = wrt_int_state%WRITE_SUBSET_R(i,j,nfield)
             if (stc(i,j,1) > spval) stc(i,j,1) = 0.0
           enddo
         enddo
       endif
       if(trim(name)=='tmp'.and.trim(layer)=='10-40 cm down') then
!$omp parallel do private(i,j)
         do j=jsta,jend
           do i=1,im
             stc(i,j,2) = wrt_int_state%WRITE_SUBSET_R(i,j,nfield)
             if (stc(i,j,2) > spval) stc(i,j,2) = 0.0
           enddo
         enddo
       endif
       if(trim(name)=='tmp'.and.trim(layer)=='40-100 cm down') then
!$omp parallel do private(i,j)
         do j=jsta,jend
           do i=1,im
             stc(i,j,3) = wrt_int_state%WRITE_SUBSET_R(i,j,nfield)
             if (stc(i,j,3) > spval) stc(i,j,3) = 0.0
           enddo
         enddo
       endif
       if(trim(name)=='tmp'.and.trim(layer)=='100-200 cm down') then
!$omp parallel do private(i,j)
         do j=jsta,jend
           do i=1,im
             stc(i,j,4) = wrt_int_state%WRITE_SUBSET_R(i,j,nfield)
             if (stc(i,j,4) > spval) stc(i,j,4) = 0.0
           enddo
         enddo
       endif
! time averaged incoming sfc longwave using nemsio
       if(trim(name)=='dlwrf'.and.trim(layer)=='sfc') then
!$omp parallel do private(j)
         do j=jsta,jend
           alwin(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,nfield)
         enddo
       endif
! time averaged outgoing sfc longwave using gfsio
       if(trim(name)=='ulwrf'.and.trim(layer)=='sfc'.and.itr==3) then
!$omp parallel do private(i,j)
         do j=jsta,jend
           do i=1,im
             alwout(i,j) = wrt_int_state%WRITE_SUBSET_R(i,j,nfield)
             if (alwout(i,j) < spval) then
               alwout(i,j) = - alwout(i,j)   ! CLDRAD puts a minus sign before gribbing
             else
               alwout(i,j) = spval
             endif
           enddo
         enddo
       endif
! time averaged outgoing model top longwave using gfsio
       if(trim(name)=='ulwrf'.and.trim(layer)=='nom. top') then
!$omp parallel do private(j)
         do j=jsta,jend
           alwtoa(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,nfield)
         enddo
       endif
! time averaged incoming sfc shortwave using gfsio
       if(trim(name)=='dswrf'.and.trim(layer)=='sfc') then
!$omp parallel do private(j)
         do j=jsta,jend
           aswin(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,nfield)
         enddo
       endif
! time averaged incoming sfc uv-b using getgb
       if(trim(name)=='duvb'.and.trim(layer)=='sfc') then
!$omp parallel do private(j)
         do j=jsta,jend
           auvbin(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,nfield)
         enddo
       endif
! time averaged incoming sfc clear sky uv-b using getgb
       if(trim(name)=='cduvb'.and.trim(layer)=='sfc') then
!$omp parallel do private(j)
         do j=jsta,jend
           auvbinc(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,nfield)
         enddo
       endif
! time averaged outgoing sfc shortwave using gfsio
       if(trim(name)=='uswrf'.and.trim(layer)=='sfc'.and.itr==3) then
!$omp parallel do private(i,j)
         do j=jsta,jend
           do i=1,im
             aswout(i,j) = wrt_int_state%WRITE_SUBSET_R(i,j,nfield)
             if (aswout(i,j) < spval) then
               aswout(i,j) = - aswout(i,j) ! CLDRAD puts a minus sign before gribbing
             else
               aswout(i,j) = spval
             endif
           enddo
         enddo
       endif
! time averaged model top outgoing shortwave
       if(trim(name)=='uswrf'.and.trim(layer)=='nom. top') then
!$omp parallel do private(j)
         do j=jsta,jend
           aswtoa(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,nfield)
         enddo
       endif
! time averaged surface sensible heat flux, multiplied by -1 because wrf model flux
! has reversed sign convention using gfsio
       if(trim(name)=='shtfl'.and.trim(layer)=='sfc'.and.itr==3) then
!$omp parallel do private(i,j)
         do j=jsta,jend
           do i=1,im
             sfcshx(i,j) = wrt_int_state%WRITE_SUBSET_R(i,j,nfield)
             if (sfcshx(i,j) < spval) then
               sfcshx(i,j) = - sfcshx(i,j)
             else
               sfcshx(i,j) = spval
             endif
           enddo
         enddo
       endif
! time averaged surface latent heat flux, multiplied by -1 because wrf model flux
! has reversed sign vonvention using gfsio
       if(trim(name)=='lhtfl'.and.trim(layer)=='sfc'.and.itr==3) then
!$omp parallel do private(i,j)
         do j=jsta,jend
           do i=1,im
             sfclhx(i,j) = wrt_int_state%WRITE_SUBSET_R(i,j,nfield)
             if (sfclhx(i,j) < spval) then
               sfclhx(i,j) = - sfclhx(i,j)
             else
               sfclhx(i,j) = spval
             endif
           enddo
         enddo
       endif
! time averaged ground heat flux using nemsio
       if(trim(name)=='gflux'.and.trim(layer)=='sfc'.and.itr==3) then
!$omp parallel do private(i,j)
         do j=jsta,jend
           do i=1,im
             subshx(i,j) = wrt_int_state%WRITE_SUBSET_R(i,j,nfield)
             if (subshx(i,j) > spval) subshx(i,j) = spval
           enddo
         enddo
!         print *,'in set_var,gflux=',maxval(subshx(:,jsta:jend)),minval(subshx(:,jsta:jend))
       endif
! time averaged zonal momentum flux using gfsio
       if(trim(name)=='uflx'.and.trim(layer)=='sfc'.and.itr==3) then
!$omp parallel do private(j)
         do j=jsta,jend
           sfcux(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,nfield)
         enddo
       endif
! time averaged meridional momentum flux using nemsio
       if(trim(name)=='vflx'.and.trim(layer)=='sfc'.and.itr==3) then
!$omp parallel do private(j)
         do j=jsta,jend
           sfcvx(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,nfield)
         enddo
       endif
! time averaged zonal gravity wave stress using nemsio
       if(trim(name)=='u-gwd'.and.trim(layer)=='sfc'.and.itr==3) then
!$omp parallel do private(j)
         do j=jsta,jend
           gtaux(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,nfield)
         enddo
       endif
! time averaged meridional gravity wave stress using getgb
       if(trim(name)=='v-gwd'.and.trim(layer)=='sfc'.and.itr==3) then
!$omp parallel do private(j)
         do j=jsta,jend
           gtauy(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,nfield)
         enddo
       endif
! time averaged accumulated potential evaporation
       if(trim(name)=='pevpr'.and.trim(layer)=='sfc'.and.itr==0) then
!$omp parallel do private(i,j)
         do j=jsta,jend
           do i=1,im
             potevp(i,j) = wrt_int_state%WRITE_SUBSET_R(i,j,nfield)
             if (potevp(i,j) > spval) potevp(i,j) = spval
           enddo
         enddo
!         print *,'in set_postvar,pevpr=',maxval(potevp(1:im,jsta:jend)),minval(potevp(1:im,jsta:jend))
       endif
! 10 m u using nemsio
       if(trim(name)=='ugrd'.and.trim(layer)=='10 m above gnd') then
!$omp parallel do private(j)
         do j=jsta,jend
           u10(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,nfield)
         enddo
       endif
! 10 m v using gfsio
       if(trim(name)=='vgrd'.and.trim(layer)=='10 m above gnd') then
!$omp parallel do private(j)
         do j=jsta,jend
           v10(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,nfield)
         enddo
       endif
! vegetation type
       if(trim(name)=='vgtyp'.and.trim(layer)=='sfc') then
!$omp parallel do private(i,j)
         do j=jsta,jend
           do i=1,im
             buf(i,j) = wrt_int_state%WRITE_SUBSET_R(i,j,nfield)
             if (buf(i,j) < spval) then
               ivgtyp(i,j) = nint(buf(i,j))
             else
               ivgtyp(i,j) = 0.0 !need to feed reasonable value to crtm
             endif
           enddo
         enddo
       endif
! soil type 
       if(trim(name)=='sotyp'.and.trim(layer)=='sfc') then
!$omp parallel do private(j)
         do j=jsta,jend
           buf(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,nfield)
         enddo
!$omp parallel do private(i,j)
         do j=jsta,jend
         do i=1,im
           if(buf(i,j)<spval) then
             isltyp(i,j)=nint(buf(i,j))
           else
             isltyp(i,j)=0.  !need to feed reasonable value to crtm
           endif
         enddo
         enddo
!       print *,'in set_postvar,isltyp=',maxval(isltyp(:,jsta:jend)),minval(isltyp(:,jsta:jend))
       endif
! retrieve inst convective cloud top, GFS has cloud top pressure instead of index,
! will need to modify CLDRAD.f to use pressure directly instead of index
       if(trim(name)=='pres'.and.trim(layer)=='convect-cld top') then
!$omp parallel do private(i,j)
         do j=jsta,jend
           do i=1,im
             ptop(i,j) = wrt_int_state%WRITE_SUBSET_R(i,j,nfield)
             if (ptop(i,j) > spval) ptop(i,j) = spval
           enddo
         enddo
       endif
! retrieve inst convective cloud bottom, GFS has cloud top pressure instead of index,
! will need to modify CLDRAD.f to use pressure directly instead of index
       if(trim(name)=='pres'.and.trim(layer)=='convect-cld bot') then
!$omp parallel do private(i,j)
         do j=jsta,jend
           do i=1,im
             pbot(i,j) = wrt_int_state%WRITE_SUBSET_R(i,j,nfield)
             if (pbot(i,j) > spval) pbot(i,j) = spval
           enddo
         enddo
       endif
! retrieve time averaged low cloud top pressure using nemsio
       if(trim(name)=='pres'.and.trim(layer)=='low cld top') then
!$omp parallel do private(i,j)
         do j=jsta,jend
           do i=1,im
             ptopl(i,j) = wrt_int_state%WRITE_SUBSET_R(i,j,nfield)
             if (ptopl(i,j) > spval) ptopl(i,j) = spval
           enddo
         enddo
       endif
! retrieve time averaged low cloud bottom pressure using nemsio
       if(trim(name)=='pres'.and.trim(layer)=='low cld bot') then
!$omp parallel do private(i,j)
         do j=jsta,jend
           do i=1,im
             pbotl(i,j) = wrt_int_state%WRITE_SUBSET_R(i,j,nfield)
             if (pbotl(i,j) > spval) pbotl(i,j) = spval
           enddo
         enddo
       endif
! retrieve time averaged low cloud top temperature using nemsio
       if(trim(name)=='tmp'.and.trim(layer)=='low cld top') then
!$omp parallel do private(i,j)
         do j=jsta,jend
           do i=1,im
             ttopl(i,j) = wrt_int_state%WRITE_SUBSET_R(i,j,nfield)
             if (ttopl(i,j) > spval) ttopl(i,j) = spval
           enddo
         enddo
       endif
! retrieve time averaged middle cloud top pressure using nemsio
       if(trim(name)=='pres'.and.trim(layer)=='mid cld top') then
!$omp parallel do private(i,j)
         do j=jsta,jend
           do i=1,im
             ptopm(i,j) = wrt_int_state%WRITE_SUBSET_R(i,j,nfield)
             if (ptopm(i,j) > spval) ptopm(i,j) = spval
           enddo
         enddo
       endif
! retrieve time averaged middle cloud bottom pressure using  nemsio
       if(trim(name)=='pres'.and.trim(layer)=='mid cld bot') then
!$omp parallel do private(i,j)
         do j=jsta,jend
           do i=1,im
             pbotm(i,j) = wrt_int_state%WRITE_SUBSET_R(i,j,nfield)
             if (pbotm(i,j) > spval) pbotm(i,j) = spval
           enddo
         enddo
       endif
! retrieve time averaged middle cloud top temperature using nemsio
       if(trim(name)=='tmp'.and.trim(layer)=='mid cld top') then
!$omp parallel do private(i,j)
         do j=jsta,jend
           do i=1,im
             ttopm(i,j) = wrt_int_state%WRITE_SUBSET_R(i,j,nfield)
             if (ttopm(i,j) > spval) ttopm(i,j) = spval
           enddo
         enddo
       endif
! retrieve time averaged high cloud top pressure using nemsio *********
       if(trim(name)=='pres'.and.trim(layer)=='high cld top') then
!$omp parallel do private(i,j)
         do j=jsta,jend
           do i=1,im
             ptoph(i,j) = wrt_int_state%WRITE_SUBSET_R(i,j,nfield)
             if (ptoph(i,j) > spval) ptoph(i,j) = spval
           enddo
         enddo
       endif
! retrieve time averaged high cloud bottom pressure using  nemsio
       if(trim(name)=='pres'.and.trim(layer)=='high cld bot') then
!$omp parallel do private(i,j)
         do j=jsta,jend
           do i=1,im
             pboth(i,j) = wrt_int_state%WRITE_SUBSET_R(i,j,nfield)
             if (pboth(i,j) > spval) pboth(i,j) = spval
           enddo
         enddo
       endif
! retrieve time averaged high cloud top temperature using nemsio
       if(trim(name)=='tmp'.and.trim(layer)=='high cld top') then
!$omp parallel do private(i,j)
         do j=jsta,jend
           do i=1,im
             ttoph(i,j) = wrt_int_state%WRITE_SUBSET_R(i,j,nfield)
             if (ttoph(i,j) > spval) ttoph(i,j) = spval
           enddo
         enddo
       endif
! retrieve boundary layer cloud cover using nemsio
       if(trim(name)=='tcdc'.and.trim(layer)=='bndary-layer cld'.and.itr==3) then
!$omp parallel do private(i,j)
         do j=jsta,jend
           do i=1,im
             pblcfr(i,j) = wrt_int_state%WRITE_SUBSET_R(i,j,nfield)
             if (pblcfr(i,j) < spval) then
               pblcfr(i,j) = pblcfr(i,j) * 0.01   ! convert to fraction
             else
               pblcfr(i,j) = spval
             endif
           enddo
         enddo
       endif
! retrieve cloud work function using nemsio
       if(trim(name)=='cwork'.and.trim(layer)=='atmos col'.and.itr==3) then
!$omp parallel do private(i,j)
         do j=jsta,jend
           do i=1,im
             cldwork(i,j) = wrt_int_state%WRITE_SUBSET_R(i,j,nfield)
             if (cldwork(i,j) > spval) cldwork(i,j) = spval
           enddo
         enddo
       endif
! retrieve water runoff using nemsio
       if(trim(name)=='watr'.and.trim(layer)=='sfc'.and.itr==4) then
!$omp parallel do private(i,j)
         do j=jsta,jend
           do i=1,im
             runoff(i,j) = wrt_int_state%WRITE_SUBSET_R(i,j,nfield)
             if (runoff(i,j) > spval) runoff(i,j) = spval
           enddo
         enddo
       endif
! retrieve shelter max temperature using nemsio
       if(trim(name)=='tmax'.and.trim(layer)=='2 m above gnd') then
!$omp parallel do private(i,j)
         do j=jsta,jend
           do i=1,im
             maxtshltr(i,j) = wrt_int_state%WRITE_SUBSET_R(i,j,nfield)
             if (maxtshltr(i,j) > spval) maxtshltr(i,j) = spval
           enddo
         enddo
       endif
! retrieve shelter max temperature using nemsio
       if(trim(name)=='tmin'.and.trim(layer)=='2 m above gnd') then
!$omp parallel do private(i,j)
         do j=jsta,jend
           do i=1,im
             mintshltr(i,j) = wrt_int_state%WRITE_SUBSET_R(i,j,nfield)
             if (mintshltr(i,j) > spval) mintshltr(i,j) = spval
           enddo
         enddo
       endif
! retrieve ice thickness using nemsio
       if(trim(name)=='icetk'.and.trim(layer)=='sfc') then
!$omp parallel do private(i,j)
         do j=jsta,jend
           do i=1,im
             dzice(i,j) = wrt_int_state%WRITE_SUBSET_R(i,j,nfield)
             if (dzice(i,j) > spval) dzice(i,j) = spval
           enddo
         enddo
       endif
! retrieve wilting point using nemsio
       if(trim(name)=='wilt'.and.trim(layer)=='sfc') then
!$omp parallel do private(i,j)
         do j=jsta,jend
           do i=1,im
             smcwlt(i,j) = wrt_int_state%WRITE_SUBSET_R(i,j,nfield)
             if (smcwlt(i,j) > spval) smcwlt(i,j) = spval
           enddo
         enddo
       endif
! retrieve sunshine duration using nemsio
       if(trim(name)=='sunsd'.and.trim(layer)=='sfc') then
!$omp parallel do private(j)
         do j=jsta,jend
           suntime(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,nfield)
         enddo
       endif
! retrieve field capacity using nemsio
       if(trim(name)=='fldcp'.and.trim(layer)=='sfc') then
!$omp parallel do private(i,j)
         do j=jsta,jend
           do i=1,im
             fieldcapa(i,j) = wrt_int_state%WRITE_SUBSET_R(i,j,nfield)
             if (fieldcapa(i,j) > spval) fieldcapa(i,j) = spval
           enddo
         enddo
       endif
!
      endif  ! end of INDX_2D
!
      if(iostatusD3D==0.and. INDX_2D>0)then ! start reading d3d file
!
       INDX_2D2=INDEX(NAMEWL,"_")
       if(INDX_2D2>0) NAME=NAMEWL(1:INDX_2D2-1)
!
! retrieve longwave tendency using getgb
       if(trim(name)==trim(avbl(41))) then
        do l=1,lm
         ll=lm-l+1 !flip 3d fields to count from top down
         do j=jsta,jend
           rlwtt(:,j,ll)=wrt_int_state%WRITE_SUBSET_R(:,j,nfield)
         enddo
        end do
       endif
! retrieve shortwave tendency using getgb
       if(trim(name)==trim(avbl(40))) then
        do l=1,lm
         ll=lm-l+1 !flip 3d fields to count from top down
         do j=jsta,jend
           rswtt(:,j,ll)=wrt_int_state%WRITE_SUBSET_R(:,j,nfield)
         enddo
        end do
       endif
! retrieve vertical diffusion tendency using getgb
       if(trim(name)=='VDIFF TNDY') then
        do l=1,lm
         ll=lm-l+1 !flip 3d fields to count from top down
         do j=jsta,jend
           vdifftt(:,j,ll)=wrt_int_state%WRITE_SUBSET_R(:,j,nfield)
         enddo
        end do
       endif
! retrieve deep convective tendency using getgb
       if(trim(name)==trim(avbl(79))) then
        do l=1,lm
         ll=lm-l+1 !flip 3d fields to count from top down
         do j=jsta,jend
           tcucn(:,j,ll)=wrt_int_state%WRITE_SUBSET_R(:,j,nfield)
         enddo
        end do
       endif
! retrieve shallow convective tendency using getgb
       if(trim(name)=='S VDIFF TNDY') then
        Indexp=358
        do l=1,lm
         ll=lm-l+1 !flip 3d fields to count from top down
         do j=jsta,jend
           tcucns(:,j,ll)=wrt_int_state%WRITE_SUBSET_R(:,j,nfield)
         enddo
        end do
       endif
! retrieve grid scale latent heat tendency using getgb
       if(trim(name)==trim(avbl(78))) then
        do l=1,lm
         ll=lm-l+1 !flip 3d fields to count from top down
         do j=jsta,jend
           train(:,j,ll)=wrt_int_state%WRITE_SUBSET_R(:,j,nfield)
         enddo
        end do
       endif
! retrieve vertical diffusion moistening using getgb
       if(trim(name)=='Vertical diffusion moistening') then
        IndexP=360
        do l=1,lm
         ll=lm-l+1 !flip 3d fields to count from top down
         do j=jsta,jend
           vdiffmois(:,j,ll)=wrt_int_state%WRITE_SUBSET_R(:,j,nfield)
         enddo
        end do
       endif
! retrieve deep convection moistening using getgb
       if(trim(name)=='deep convection moistening') then
        Indexp=361
        do l=1,lm
         ll=lm-l+1 !flip 3d fields to count from top down
         do j=jsta,jend
           dconvmois(:,j,ll)=wrt_int_state%WRITE_SUBSET_R(:,j,nfield)
         enddo
        end do
       endif
! retrieve shallow convection moistening using getgb
       if(trim(name)=='shallow convection moistening') then
        Indexp=362
        do l=1,lm
         ll=lm-l+1 !flip 3d fields to count from top down
         do j=jsta,jend
           sconvmois(:,j,ll)=wrt_int_state%WRITE_SUBSET_R(:,j,nfield)
         enddo
        end do
       endif
! retrieve non-radiation tendency using getgb
       if(trim(name)=='non-radiation tendency') then
        Indexp=363
        do l=1,lm
         ll=lm-l+1 !flip 3d fields to count from top down
         do j=jsta,jend
           nradtt(:,j,ll)=wrt_int_state%WRITE_SUBSET_R(:,j,nfield)    
         end do
        end do
       endif
! retrieve Vertical diffusion of ozone using getgb
       if(trim(name)=='Vertical diffusion of ozone') then
        Indexp=364
        do l=1,lm
         ll=lm-l+1 !flip 3d fields to count from top down
         do j=jsta,jend
           o3vdiff(:,j,ll)=wrt_int_state%WRITE_SUBSET_R(:,j,nfield)
         enddo
        end do
       endif
! retrieve ozone production using getgb
       if(trim(name)=='Ozone production') then
        Indexp=365
        do l=1,lm
         ll=lm-l+1 !flip 3d fields to count from top down
         do j=jsta,jend
           o3prod(:,j,ll)=wrt_int_state%WRITE_SUBSET_R(:,j,nfield)
         enddo
        end do
       endif
! retrieve ozone tendency using getgb
       if(trim(name)=='Ozone tendency') then
        Indexp=366
        do l=1,lm
         ll=lm-l+1 !flip 3d fields to count from top down
         do j=jsta,jend
           o3tndy(:,j,ll)=wrt_int_state%WRITE_SUBSET_R(:,j,nfield)
         enddo
        end do
       endif
! retrieve mass weighted PV using getgb
       if(trim(name)=='Mass weighted PV') then
        Indexp=367
        do l=1,lm
         ll=lm-l+1 !flip 3d fields to count from top down
         do j=jsta,jend
           mwpv(:,j,ll)=wrt_int_state%WRITE_SUBSET_R(:,j,nfield)
         enddo
        end do
       endif
! retrieve OZONE TNDY using getgb
       if(trim(name)=='unknown') then
        Indexp=368
        do l=1,lm
         ll=lm-l+1 !flip 3d fields to count from top down
         do j=jsta,jend
           unknown(:,j,ll)=wrt_int_state%WRITE_SUBSET_R(:,j,nfield)
         enddo
        end do
       endif
! retrieve vertical diffusion zonal acceleration
       if(trim(name)=='VDIFF Z ACCE') then
        Indexp=369
        do l=1,lm
         ll=lm-l+1 !flip 3d fields to count from top down
         do j=jsta,jend
           vdiffzacce(:,j,ll)=wrt_int_state%WRITE_SUBSET_R(:,j,nfield)
         enddo
        end do
       endif
! retrieve gravity drag zonal acceleration
       if(trim(name)=='G DRAG Z ACCE') then
        Indexp=370
        do l=1,lm
         ll=lm-l+1 !flip 3d fields to count from top down
         do j=jsta,jend
           zgdrag(:,j,ll)=wrt_int_state%WRITE_SUBSET_R(:,j,nfield)
         enddo
        end do
       endif
! retrieve convective U momemtum mixing
       if(trim(name)=='CNVCT U M MIX') then
        Indexp=371
        do l=1,lm
         ll=lm-l+1 !flip 3d fields to count from top down
         do j=jsta,jend
           cnvctummixing(:,j,ll)=wrt_int_state%WRITE_SUBSET_R(:,j,nfield)
         enddo
        end do
       endif
! retrieve vertical diffusion meridional acceleration
       if(trim(name)=='VDIFF M ACCE') then
        Indexp=372
        do l=1,lm
         ll=lm-l+1 !flip 3d fields to count from top down
         do j=jsta,jend
           vdiffmacce(:,j,ll)=wrt_int_state%WRITE_SUBSET_R(:,j,nfield)
         enddo
        end do
       endif
! retrieve gravity drag meridional acceleration
       if(trim(name)=='G DRAG M ACCE') then
        Indexp=373
        do l=1,lm
         ll=lm-l+1 !flip 3d fields to count from top down
         do j=jsta,jend
           mgdrag(:,j,ll)=wrt_int_state%WRITE_SUBSET_R(:,j,nfield)
         enddo
        end do
       endif
! retrieve convective V momemtum mixing
       if(trim(name)=='CNVCT V M MIX') then
        Indexp=374
        do l=1,lm
         ll=lm-l+1 !flip 3d fields to count from top down
         do j=jsta,jend
           cnvctvmmixing(:,j,ll)=wrt_int_state%WRITE_SUBSET_R(:,j,nfield)
         enddo
        end do
       endif
! retrieve nonconvective cloud fraction
       if(trim(name)=='N CNVCT CLD FRA') then
        Indexp=375
        do l=1,lm
         ll=lm-l+1 !flip 3d fields to count from top down
         do j=jsta,jend
           ncnvctcfrac(:,j,ll)=wrt_int_state%WRITE_SUBSET_R(:,j,nfield)
         enddo
        end do
       endif
! retrieve convective upward mass flux
       if(trim(name)=='CNVCT U M FLX') then
        Indexp=391
        do l=1,lm
         ll=lm-l+1 !flip 3d fields to count from top down
         do j=jsta,jend
           cnvctumflx(:,j,ll)=wrt_int_state%WRITE_SUBSET_R(:,j,nfield)
         enddo
        end do
       endif
! retrieve convective downward mass flux
       if(trim(name)=='CNVCT D M FLX') then
        Indexp=392
        do l=1,lm
         ll=lm-l+1 !flip 3d fields to count from top down
         do j=jsta,jend
           cnvctdmflx(:,j,ll)=wrt_int_state%WRITE_SUBSET_R(:,j,nfield)
         enddo
        end do
       endif
! retrieve nonconvective detraintment flux
       if(trim(name)=='CNVCT DET M FLX') then
        Indexp=393
        do l=1,lm
         ll=lm-l+1 !flip 3d fields to count from top down
         do j=jsta,jend
           cnvctdetmflx(:,j,ll)=wrt_int_state%WRITE_SUBSET_R(:,j,nfield)
         enddo
        end do
       endif
! retrieve cnvct gravity drag zonal acceleration
       if(trim(name)=='CNVCT G DRAG Z ACCE') then
        Indexp=394
        do l=1,lm
         ll=lm-l+1 !flip 3d fields to count from top down
         do j=jsta,jend
           cnvctzgdrag(:,j,ll)=wrt_int_state%WRITE_SUBSET_R(:,j,nfield)
         enddo
        end do
       endif
! retrieve cnvct gravity drag meridional acceleration
       if(trim(name)=='CNVCT G DRAG M ACCE') then
        Indexp=395
        do l=1,lm
         ll=lm-l+1 !flip 3d fields to count from top down
         do j=jsta,jend
           cnvctmgdrag(:,j,ll)=wrt_int_state%WRITE_SUBSET_R(:,j,nfield)
         enddo
        end do
       endif
!
       end if ! end of d3d file read
!
!for aer file
! Retrieve aer fields if it's listed (GOCART)
      if(trim(wrt_int_state%FILENAME_BASE(NBDL))=='AER.F')then 
!
        INDX_2D2=INDEX(NAMEWL,"_")
        if(INDX_2D2>0) then
          NAME=trim(NAMEWL(1:INDX_2D2-1))
        else
          NAME=TRIM(NAMEWL)
        endif
        INDX_2D3=INDEX(NAMEWL,"_",BACK=.true.)
        if(INDX_2D3>0) then
          layer=NAMEWL(INDX_2D3+1:)
        else
          layer='sfc'
        endiF
        itr=0
        if(INDEX(NAMEWL,"_ave") >0) then
          itr=3
        elseif(INDEX(NAMEWL,"_acc") >0) then
          itr=4
        elseif(INDEX(NAMEWL,"_win") >0) then
          itr=2
        endif

        print *, 'iostatus for aer file, field name=',trim(name),'layer=',trim(layer), &
           'trim(name)==ducmass',trim(name)=='DUSMASS'
!DUEM
        do k=1,nbin_du
          write(ck,'(i3.3)')k
!duem
          VarAerName='DUEM'//ck
          VarAerLevtyp='atmos col'
          if(trim(name)==trim(VarAerName).and.trim(layer)==trim(VarAerLevtyp)) then
!$omp parallel do private(j)
            do j=jsta,jend
              duem(:,j,k)=wrt_int_state%WRITE_SUBSET_R(:,j,nfield)
            enddo
          endif
!dusd
          VarAerName='DUSD'//ck
          VarAerLevtyp='atmos col'
          if(trim(name)==trim(VarAerName).and.trim(layer)==trim(VarAerLevtyp)) then
!$omp parallel do private(j)
            do j=jsta,jend
              dusd(:,j,k)=wrt_int_state%WRITE_SUBSET_R(:,j,nfield)
            enddo
          endif
!dudp
          VarAerName='DUDP'//ck
          VarAerLevtyp='atmos col'
          if(trim(name)==trim(VarAerName).and.trim(layer)==trim(VarAerLevtyp)) then
!$omp parallel do private(j)
            do j=jsta,jend
              dudp(:,j,k)=wrt_int_state%WRITE_SUBSET_R(:,j,nfield)
            enddo
          endif
!duwt
          VarAerName='DUWT'//ck
          VarAerLevtyp='atmos col'
          if(trim(name)==trim(VarAerName).and.trim(layer)==trim(VarAerLevtyp)) then
!$omp parallel do private(j)
            do j=jsta,jend
              duwt(:,j,k)=wrt_int_state%WRITE_SUBSET_R(:,j,nfield)
            enddo
          endif

! enddo aerosol k
        enddo
!DUSMASS
        if(trim(name)=='DUSMASS'.and.trim(layer)=='atmos col') then
!$omp parallel do private(j)
          do j=jsta,jend
            dusmass(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,nfield)
          enddo
          print *,'in aero file,dusmass=',maxval(dusmass(:,jsta:jend)),minval(dusmass(:,jsta:jend))
        endif
!DUCMASS
        if(trim(name)=='DUCMASS'.and.trim(layer)=='atmos col') then
!$omp parallel do private(j)
          do j=jsta,jend
            ducmass(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,nfield)
          enddo
        endif
!DUSMASS25
        if(trim(name)=='DUSMASS25'.and.trim(layer)=='atmos col') then
!$omp parallel do private(j)
          do j=jsta,jend
            dusmass25(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,nfield)
          enddo
        endif
!DUCMASS25
        if(trim(name)=='DUCMASS25'.and.trim(layer)=='atmos col') then
!$omp parallel do private(j)
          do j=jsta,jend
            ducmass25(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,nfield)
          enddo
        endif
!
!end aer file
      endif

      ENDDO
!
! some derived fields, from grid file
!
      extrafield_if: if(trim(wrt_int_state%FILENAME_BASE(NBDL))=='SIG.F') then
! pos east
!
       call collect_loc(gdlat,dummy)
       if(me.eq.0)then
        latstart=nint(dummy(1,1)*gdsdegr)
        latlast=nint(dummy(im,jm)*gdsdegr)
!        print*,'laststart,latlast B bcast= ',latstart,latlast
       endif
       call mpi_bcast(latstart,1,MPI_INTEGER,0,mpi_comm_comp,iret)
       call mpi_bcast(latlast,1,MPI_INTEGER,0,mpi_comm_comp,iret)
       call collect_loc(gdlon,dummy)
       if(me.eq.0)then
        lonstart=nint(dummy(1,1)*gdsdegr)
        lonlast=nint(dummy(im,jm)*gdsdegr)
       endif
       call mpi_bcast(lonstart,1,MPI_INTEGER,0,mpi_comm_comp,iret)
       call mpi_bcast(lonlast,1,MPI_INTEGER,0,mpi_comm_comp,iret)
!      write(6,*)'lonstart,lonlast A calling bcast=',lonstart,lonlast
!
! construct interface pressure from model top (which is zero) and dp from top down
! PDTOP
!        print*,'PT, PDTOP= ',PT,PDTOP
        do l=2,lm
        do j=jsta,jend
          do i=1,im
            pint(i,j,l)=pint(i,j,l-1)+buf3d(i,j,l-1)
!jw real4            ALPINT(I,J,L)=ALOG(PINT(I,J,L))
            ALPINT(I,J,L)=LOG(PINT(I,J,L))
!            if(i==ii .and. j==jj)print*,'sample pint,pmid'                     &
!              ,i,j,l,pint(i,j,l),pmid(i,j,l)
          end do
        end do
!        print *,'in set_postvar,l=',l,'pint=',maxval(pint(1:im,jsta:jend,l)), &
!         minval(pint(1:im,jsta:jend,l))
        end do
!
! GFS probably does not use sigt4, set it to sig*t^4
! GFS probably does not use sigt4, set it to sig*t^4
!$omp parallel do private(i,j,tlmh)
      Do j=jsta,jend
        Do i=1,im
          TLMH       = T(I,J,LM)
          Sigt4(I,j) =  5.67E-8*TLMH*TLMH*TLMH*TLMH
        End do
      End do

! TG is not used, skip it for now
!     allocate(qstl(lm))
      allocate(p2d(im,lm),t2d(im,lm),q2d(im,lm),cw2d(im,lm),          &
               qs2d(im,lm),cfr2d(im,lm))
      do j=jsta,jend
!$omp parallel do private(i,k,es)
        do k=1,lm
          do i=1,im
          p2d(i,k)  = pmid(i,j,k)*0.01
          t2d(i,k)  = t(i,j,k)
          q2d(i,k)  = q(i,j,k)
          cw2d(i,k) = cwm(i,j,k)
          es = fpvsnew(t(i,j,k))
          es = min(es,pmid(i,j,k))
          qs2d(i,k) = con_eps*es/(pmid(i,j,k)+con_epsm1*es)!saturation q for GFS
          enddo
        enddo
        call progcld1                                                 &
!...................................
!  ---  inputs:
             ( p2d,t2d,q2d,qs2d,cw2d,im,lm,0,                         &
!  ---  outputs:
               cfr2d                                                  &
              )
!$omp parallel do private(i,k)
        do k=1,lm
          do i=1,im
            cfr(i,j,k) = cfr2d(i,k)
          enddo
        end do
      end do
      deallocate(p2d,t2d,q2d,qs2d,cw2d,cfr2d)
!
!
!
!htop
      htop=spval
      do j=jsta,jend
        do i=1,im
          if(ptop(i,j) <= 0.0)ptop(i,j)=spval
          if(ptop(i,j) < spval)then
           do l=1,lm
            if(ptop(i,j) <= pmid(i,j,l))then
             htop(i,j)=l
!             if(i==ii .and. j==jj)print*,'sample ptop,pmid pmid-1,pint= ',   &
!                ptop(i,j),pmid(i,j,l),pmid(i,j,l-1),pint(i,j,l),htop(i,j)
             exit
            end if
           end do
          end if
        end do
       end do
!hbot
      hbot=spval
      do j=jsta,jend
        do i=1,im
          if(pbot(i,j) <= 0.0)pbot(i,j)=spval
!         if(.not.lb(i,j))print*,'false bitmask for pbot at '
!     +     ,i,j,pbot(i,j)
          if(pbot(i,j) .lt. spval)then
           do l=lm,1,-1
            if(pbot(i,j) >= pmid(i,j,l))then
             hbot(i,j)=l
!             if(i==ii .and. j==jj)print*,'sample pbot,pmid= ',    &
!                pbot(i,j),pmid(i,j,l),hbot(i,j)
             exit
            end if
           end do
          end if
        end do
       end do

!
!more fields need to be computed
!!!!! COMPUTE Z, GFS integrates Z on mid-layer instead
!!! use GFS contants to see if height becomes more aggreable to GFS pressure grib file
!
       do j = jsta, jend
        do i = 1, im
            ZINT(I,J,LM+1)=FIS(I,J)/con_G
            FI(I,J,1)=FIS(I,J)+T(I,J,LM)                                &
            *(Q(I,J,LM)*con_fvirt+1.0)*con_rd                           &
! using GFS consts
!jw real4            *(ALPINT(I,J,LM+1)-ALOG(PMID(I,J,LM)))
            *(ALPINT(I,J,LM+1)-LOG(PMID(I,J,LM)))
            ZMID(I,J,LM)=FI(I,J,1)/con_G
        end do
       end do
!       print *,'zint=',maxval(zint(1:im,jsta:jend,lp1)),minval(zint(1:im,jsta:jend,lp1)), &
!         'zmid=',maxval(zmid(1:im,jsta:jend,lm)),minval(zmid(1:im,jsta:jend,lm))

! SECOND, INTEGRATE HEIGHT HYDROSTATICLY, GFS integrate height on mid-layer
      DO L=LM-1,1,-1
       do j = jsta, jend
        do i = 1, im
         FI(I,J,2)=0.5*(T(I,J,L)*(Q(I,J,L)*con_fvirt+1.0)               &
                   +T(I,J,L+1)*(Q(I,J,L+1)*con_fvirt+1.0))*con_rd*      &
!         FI(I,J,2)=0.5*(T(I,J,L)*(Q(I,J,L)*D608+1.0)
!     1            +T(I,J,L+1)*(Q(I,J,L+1)*D608+1.0))*RD*
!jw real4                   (ALOG(PMID(I,J,L+1))-ALOG(PMID(I,J,L)))              &
                   (LOG(PMID(I,J,L+1))-LOG(PMID(I,J,L)))              &
                   +FI(I,J,1)
         ZMID(I,J,L)=FI(I,J,2)/con_G

!         if(i.eq.ii.and.j.eq.jj)                                        &
!        print*,'L,sample T,Q,ALPMID(L+1),ALPMID(L),ZMID= '              &
!        ,l,T(I,J,L),Q(I,J,L),LOG(PMID(I,J,L+1)),                       &
!        LOG(PMID(I,J,L)),ZMID(I,J,L)

         FI(I,J,1)=FI(I,J,2)
        ENDDO
       ENDDO
!       print *,'in set_post,l=',l,'zmid=',maxval(zmid(1:im,jsta:jend,l)),minval(zmid(1:im,jsta:jend,l))
      END DO

      DO L=LM,2,-1  ! omit computing model top height because it's infinity
       DO J=JSTA,JEND
        DO I=1,IM
!         ZMID(I,J,L)=(ZINT(I,J,L+1)+ZINT(I,J,L))*0.5  ! ave of z
!jw real4         FACT=(ALPINT(I,J,L)-ALOG(PMID(I,J,L)))/                        &
!jw real4               (ALOG(PMID(I,J,L-1))-ALOG(PMID(I,J,L)))
         FACT=(ALPINT(I,J,L)-LOG(PMID(I,J,L)))/                        &
               (LOG(PMID(I,J,L-1))-LOG(PMID(I,J,L)))
         ZINT(I,J,L)=ZMID(I,J,L)+(ZMID(I,J,L-1)-ZMID(I,J,L))            &
                       *FACT
!         if(i.eq.ii.and.j.eq.jj) print*,'L ZINT= ',l,zint(i,j,l)
        ENDDO
       ENDDO
!       print *,'in set_post,l=',l,'zint=',maxval(zint(1:im,jsta:jend,l)),minval(zint(1:im,jsta:jend,l))
      ENDDO
      deallocate(buf3d,dummy,dummy2)
      deallocate(fi)
!
      endif extrafield_if
!         print *,'in set_postvar,pint=',maxval(pint(1:im,jsta:jend,:)), &
!           minval(pint(1:im,jsta:jend,:)),maxloc(pint(1:im,jsta:jend,:)), &
!           minval(pint(1:im,jsta:jend,:)),'file=',trim(wrt_int_state%FILENAME_BASE(NBDL))

!! -- compute air density RHOMID and remove negative tracer values

      if (gocart_on) then
      do l=1,lm
        do j=jsta,jend
          do i=1,im
           TV=T(I,J,L)*(H1+D608*MAX(Q(I,J,L),QMIN))
           RHOMID(I,J,L)=PMID(I,J,L)/(RD*TV)
           do n = 1,  NBIN_DU
             IF ( dust(i,j,l,n) .LT. SPVAL) THEN
               DUST(i,j,l,n) = MAX(DUST(i,j,l,n), 0.0)
             ENDIF
           enddo
          end do
        end do
      end do
      endif
!
      deallocate(buf)
! generate look up table for lifted parcel calculations

      THL=210.
      PLQ=70000.

      CALL TABLE(PTBL,TTBL,PT,                                     &
                RDQ,RDTH,RDP,RDTHE,PL,THL,QS0,SQS,STHE,THE0)
!
      CALL TABLEQ(TTBLQ,RDPQ,RDTHEQ,PLQ,THL,STHEQ,THE0Q)

!
!
      IF(ME.EQ.0)THEN
        WRITE(6,*)'  SPL (POSTED PRESSURE LEVELS) BELOW: '
        WRITE(6,51) (SPL(L),L=1,LSM)
   50   FORMAT(14(F4.1,1X))
   51   FORMAT(8(F8.1,1X))
      ENDIF
!
!     COMPUTE DERIVED MAP OUTPUT CONSTANTS.
      DO L = 1,LSM
!jw real4         ALSL(L) = ALOG(SPL(L))
         ALSL(L) = LOG(SPL(L))
      END DO
!
!HC WRITE IGDS OUT FOR WEIGHTMAKER TO READ IN AS KGDSIN
        if(me.eq.0)then
!        print*,'writing out igds'
        igdout=110
!        open(igdout,file='griddef.out',form='unformatted'
!     +  ,status='unknown')
        if(maptype .eq. 1)THEN  ! Lambert conformal
          WRITE(igdout)3
          WRITE(6,*)'igd(1)=',3
          WRITE(igdout)im
          WRITE(igdout)jm
          WRITE(igdout)LATSTART
          WRITE(igdout)LONSTART
          WRITE(igdout)8
          WRITE(igdout)CENLON
          WRITE(igdout)DXVAL
          WRITE(igdout)DYVAL
          WRITE(igdout)0
          WRITE(igdout)64
          WRITE(igdout)TRUELAT2
          WRITE(igdout)TRUELAT1
          WRITE(igdout)255
        ELSE IF(MAPTYPE .EQ. 2)THEN  !Polar stereographic
          WRITE(igdout)5
          WRITE(igdout)im
          WRITE(igdout)jm
          WRITE(igdout)LATSTART
          WRITE(igdout)LONSTART
          WRITE(igdout)8
          WRITE(igdout)CENLON
          WRITE(igdout)DXVAL
          WRITE(igdout)DYVAL
          WRITE(igdout)0
          WRITE(igdout)64
          WRITE(igdout)TRUELAT2  !Assume projection at +-90
          WRITE(igdout)TRUELAT1
          WRITE(igdout)255
        ELSE IF(MAPTYPE .EQ. 3)THEN  !Mercator
          WRITE(igdout)1
          WRITE(igdout)im
          WRITE(igdout)jm
          WRITE(igdout)LATSTART
          WRITE(igdout)LONSTART
          WRITE(igdout)8
          WRITE(igdout)latlast
          WRITE(igdout)lonlast
          WRITE(igdout)TRUELAT1
          WRITE(igdout)0
          WRITE(igdout)64
          WRITE(igdout)DXVAL
          WRITE(igdout)DYVAL
          WRITE(igdout)255
        ELSE IF(MAPTYPE.EQ.0 .OR. MAPTYPE.EQ.203)THEN  !A STAGGERED E-GRID
          WRITE(igdout)203
          WRITE(igdout)im
          WRITE(igdout)jm
          WRITE(igdout)LATSTART
          WRITE(igdout)LONSTART
          WRITE(igdout)136
          WRITE(igdout)CENLAT
          WRITE(igdout)CENLON
          WRITE(igdout)DXVAL
          WRITE(igdout)DYVAL
          WRITE(igdout)64
          WRITE(igdout)0
          WRITE(igdout)0
          WRITE(igdout)0
        END IF
        end if
!

      RETURN
!
    end subroutine set_postvars_gfs
!
!     
      SUBROUTINE splat(IDRT,JMAX,ASLAT)
!$$$
      implicit none
      integer(4),intent(in) :: idrt,jmax
      real,intent(out) :: ASLAT(JMAX)
      INTEGER(4),PARAMETER:: KD=SELECTED_REAL_KIND(15,45)
      REAL(8):: PK(JMAX/2),PKM1(JMAX/2),PKM2(JMAX/2)
      REAL(8):: ASLATD(JMAX/2),SP,SPMAX,EPS=10.d0*EPSILON(SP)
      integer,PARAMETER:: JZ=50
      REAL(8) BZ(JZ)
      DATA BZ        / 2.4048255577d0,  5.5200781103d0, &
       8.6537279129d0, 11.7915344391d0, 14.9309177086d0, 18.0710639679d0, &
      21.2116366299d0, 24.3524715308d0, 27.4934791320d0, 30.6346064684d0, &
      33.7758202136d0, 36.9170983537d0, 40.0584257646d0, 43.1997917132d0, &
      46.3411883717d0, 49.4826098974d0, 52.6240518411d0, 55.7655107550d0, &
      58.9069839261d0, 62.0484691902d0, 65.1899648002d0, 68.3314693299d0, &
      71.4729816036d0, 74.6145006437d0, 77.7560256304d0, 80.8975558711d0, &
      84.0390907769d0, 87.1806298436d0, 90.3221726372d0, 93.4637187819d0, &
      96.6052679510d0, 99.7468198587d0, 102.888374254d0, 106.029930916d0, &
      109.171489649d0, 112.313050280d0, 115.454612653d0, 118.596176630d0, &
      121.737742088d0, 124.879308913d0, 128.020877005d0, 131.162446275d0, &
      134.304016638d0, 137.445588020d0, 140.587160352d0, 143.728733573d0, &
      146.870307625d0, 150.011882457d0, 153.153458019d0, 156.295034268d0 /
      REAL(8):: DLT,D1=1.d0
      INTEGER(4):: JHE,JHO,J0=0
      real(8),PARAMETER :: PI=3.14159265358979d0,C=(1.d0-(2.d0/PI)**2)*0.25d0
      real(8) r
      integer jh,js,n,j
!C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!C  GAUSSIAN LATITUDES
      IF(IDRT.EQ.4) THEN
        JH=JMAX/2
        JHE=(JMAX+1)/2
        R=1.d0/SQRT((JMAX+0.5d0)**2+C)
        DO J=1,MIN(JH,JZ)
          ASLATD(J)=COS(BZ(J)*R)
        ENDDO
        DO J=JZ+1,JH
          ASLATD(J)=COS((BZ(JZ)+(J-JZ)*PI)*R)
        ENDDO
        SPMAX=1.d0
        DO WHILE(SPMAX.GT.EPS)
          SPMAX=0.d0
          DO J=1,JH
            PKM1(J)=1.d0
            PK(J)=ASLATD(J)
          ENDDO
          DO N=2,JMAX
            DO J=1,JH
              PKM2(J)=PKM1(J)
              PKM1(J)=PK(J)
              PK(J)=((2*N-1)*ASLATD(J)*PKM1(J)-(N-1)*PKM2(J))/N
            ENDDO
          ENDDO
          DO J=1,JH
            SP=PK(J)*(1.d0-ASLATD(J)**2)/(JMAX*(PKM1(J)-ASLATD(J)*PK(J)))
            ASLATD(J)=ASLATD(J)-SP
            SPMAX=MAX(SPMAX,ABS(SP))
          ENDDO
        ENDDO
!CDIR$ IVDEP
        DO J=1,JH
          ASLAT(J)=ASLATD(J)
          ASLAT(JMAX+1-J)=-ASLAT(J)
        ENDDO
        IF(JHE.GT.JH) THEN
          ASLAT(JHE)=0.d0
        ENDIF
!C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!C  EQUALLY-SPACED LATITUDES INCLUDING POLES
      ELSEIF(IDRT.EQ.0) THEN
        JH=JMAX/2
        JHE=(JMAX+1)/2
        JHO=JHE-1
        DLT=PI/(JMAX-1)
        ASLAT(1)=1.d0
        DO J=2,JH
          ASLAT(J)=COS((J-1)*DLT)
        ENDDO
!CDIR$ IVDEP
        DO J=1,JH
          ASLAT(JMAX+1-J)=-ASLAT(J)
        ENDDO
        IF(JHE.GT.JH) THEN
          ASLAT(JHE)=0.d0
        ENDIF
!C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!C  EQUALLY-SPACED LATITUDES EXCLUDING POLES
      ELSEIF(IDRT.EQ.256) THEN
        JH=JMAX/2
        JHE=(JMAX+1)/2
        JHO=JHE
        DLT=PI/JMAX
        ASLAT(1)=1.d0
        DO J=1,JH
          ASLAT(J)=COS((J-0.5)*DLT)
        ENDDO
!CDIR$ IVDEP
        DO J=1,JH
          ASLAT(JMAX+1-J)=-ASLAT(J)
        ENDDO
        IF(JHE.GT.JH) THEN
          ASLAT(JHE)=0.d0
        ENDIF
      ENDIF
!
      return
!C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     end subroutine splat

!
