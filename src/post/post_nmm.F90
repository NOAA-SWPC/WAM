#include "../../ESMFVersionDefine.h"
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
    subroutine post_run_nmm(wrt_int_state,mype,mpicomp,lead_write,          &
               mygridtype,mymaptype,mynsoil,mynfhr,mynfmin) 

!***  HISTORY
!    28May2013     Lu: Specify iostatusD3D
!    09Oct2015  S Moorthi - adding imp_physics argument to MICROINIT
!-----------------------------------------------------------------------
!*** run post on quilt
!-----------------------------------------------------------------------
!
      use MODULE_WRITE_INTERNAL_STATE
      use CTLBLK_mod, only : komax,ifhr,ifmin,MODELNAME,imp_physics
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      type(WRITE_INTERNAL_STATE),intent(in)    :: wrt_int_state
      integer,intent(in)                       :: mype
      integer,intent(in)                       :: mpicomp
      integer,intent(in)                       :: lead_write
      character(1),intent(in)                  :: mygridtype
      integer,intent(in)                       :: mymaptype
      integer,intent(in)                       :: mynsoil
      integer,intent(in)                       :: mynfhr
      integer,intent(in)                       :: mynfmin
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      integer N,NWTPG,IEOF,LCNTRL
      integer jts,jte
      integer,allocatable  :: jstagrp(:),jendgrp(:)
      integer,save :: kpo,kth,kpv
      real,dimension(komax),save  :: po, th, pv 
      logical,save :: LOG_POSTALCT=.false.

      integer,save :: iostatusD3D=-1
!
      write(0,*)'in post_run start'
!-----------------------------------------------------------------------
!*** set up dimensions
!-----------------------------------------------------------------------
!
      MODELNAME='NMM'
      ifhr=mynfhr
      ifmin=mynfmin
!
      JTS=wrt_int_state%LOCAL_JSTART(wrt_int_state%ID_FTASK_RECV_STA(MYPE))  !<-- Starting J of this write task's subsection
      JTE=wrt_int_state%LOCAL_JEND  (wrt_int_state%ID_FTASK_RECV_END(MYPE))  !<-- Ending J of this write task's subsection
      NWTPG=wrt_int_state%WRITE_TASKS_PER_GROUP
      write(0,*)'in post_run,jts=',jts,'jte=',jte,'nwtpg=',nwtpg, &
        'log_postalct=',log_postalct,'ifhr=',ifhr,'ifmin=',mynfmin
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
          JSTAGRP(N+1)=wrt_int_state%LOCAL_JSTART(wrt_int_state%ID_FTASK_RECV_STA(N+LEAD_WRITE))
          JENDGRP(N+1)=wrt_int_state%LOCAL_JEND  (wrt_int_state%ID_FTASK_RECV_END(N+LEAD_WRITE))
        ENDDO
      write(0,*)'in post_run,jstagrp=',jstagrp,'jendgrp=',jendgrp
!
!-----------------------------------------------------------------------
!*** allocate post variables
!-----------------------------------------------------------------------
!
        call post_alctvars(wrt_int_state%im,wrt_int_state%jm,        &
               wrt_int_state%lm,MYPE,wrt_int_state%WRITE_TASKS_PER_GROUP,   &
               mpicomp,mygridtype,mymaptype,wrt_int_state%post_gribversion, &
               MYNSOIL,         &
               LEAD_WRITE,JTS,JTE,JSTAGRP,JENDGRP)
      write(0,*)'in post_run,aft post_alctvars'
!
!-----------------------------------------------------------------------
!*** read namelist for pv,th,po
!-----------------------------------------------------------------------
!
         call read_postnmlt(kpo,kth,kpv,po,th,pv,wrt_int_state%nlunit, &
                            wrt_int_state%post_namelist)
      write(0,*)'in post_run,aft nmlst po'
!
!-----------------------------------------------------------------------
!*** read namelist for pv,th,po
!-----------------------------------------------------------------------
!
        LOG_POSTALCT=.true.
!
      ENDIF
!       
!-----------------------------------------------------------------------
!*** fill post variables with values from forecast results
!-----------------------------------------------------------------------
!
      call set_postvars_nmm(wrt_int_state,mpicomp,JTS,JTE)
            write(0,*)'af set_postvars'
!
       if(imp_physics==5 .or. imp_physics==85 .or. imp_physics==95) then
         call MICROINIT(imp_physics)
       endif
!
      IEOF=0
      do while( IEOF .eq. 0)
       CALL READCNTRL(kth,IEOF)
          print *,'after readcntrl,IEOF=',IEOF
!      if ( IEOF.eq.0) CALL PROCESS(KTH,KPV,TH(1:KTH),PV(1:KPV))
       if ( IEOF.eq.0) CALL PROCESS(KTH,KPV,TH(1:KTH),PV(1:KPV),iostatusD3D)
          print *,'after readcntrl,IEOF=',IEOF
      enddo
!
!          call de_allocate
      LCNTRL=14
      rewind(LCNTRL)
          print *,'after readcntrl and process'

    end subroutine post_run_nmm
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
    subroutine set_postvars_nmm(wrt_int_state,mpicomp,jts,jte)

!
!***  HISTORY
!     15Jan2013:  Sarah Lu, EL_MYJ changed to EL_PBL to be consistent with
!                 nceppost upgrade

!
!-----------------------------------------------------------------------
!*** set up int_state
!-----------------------------------------------------------------------
!
      use vrbls3d
      use vrbls2d
      use soil
      use masks
      use ctlblk_mod
      use params_mod
      use gridspec_mod
      use lookup_mod
!
      use ESMF
      use MODULE_WRITE_INTERNAL_STATE
!
!-----------------------------------------------------------------------
!
      implicit none
!
      include 'mpif.h'
!
!-----------------------------------------------------------------------
!
      type(WRITE_INTERNAL_STATE),intent(in)    :: wrt_int_state
      integer,intent(in)                       :: mpicomp
      integer,intent(in)                       :: jts,jte
!
!-----------------------------------------------------------------------
!
      integer I,ii,J,jj,L,LL,K,N,N1,N2,NPOSN_1,NPOSN_2, LENGTH
      integer iim1,jm1,im1
      integer NPOS_START,NPOS_END,indx_2d,nfield,ierr,iret
      character(ESMF_MAXSTR) :: NAME
      CHARACTER(3)           :: model_level
      REAL :: FACT
      real degrad
      REAL,dimension(:,:),allocatable :: dummy,vlat,vlon,buf
      REAL,dimension(:,:,:),allocatable :: FI
      REAL,dimension(:),allocatable  :: ETA1, ETA2, DXH
      real, parameter :: G1 = 9.8060226 ! from module_CONSTANTS
!
!-----------------------------------------------------------------------
!***  INTEGER SCALAR/1D HISTORY VARIABLES
!-----------------------------------------------------------------------
!
        N2=0                                                               !<-- Word counter for full string of integer scalar/1D data
!
        DO N=1,wrt_int_state%KOUNT_I1D(1)                                  !<-- Loop through all scalar/1D integer data
!
          NPOSN_1=(N-1)*ESMF_MAXSTR+1
          NPOSN_2=N*ESMF_MAXSTR
          NAME=wrt_int_state%NAMES_I1D_STRING(NPOSN_1:NPOSN_2)             !<-- The variable's na me
          LENGTH=wrt_int_state%LENGTH_DATA_I1D(N)                          !<-- The variable's length in words
!
          DO N1=1,LENGTH
            N2=N2+1
            if(trim(NAME)== 'MP_PHYSICS' .or. trim(NAME)=='MP_PHYSI' )  &
             imp_physics=wrt_int_state%ALL_DATA_I1D(N2)
            if(trim(NAME)== 'SF_SFC_PHYSICS'.or. trim(NAME)== 'SF_SURFA' )  &
             iSF_SURFACE_PHYSICS=wrt_int_state%ALL_DATA_I1D(N2)
            if(trim(NAME)== 'IHRST' )  &
             ihrst=wrt_int_state%ALL_DATA_I1D(N2)
            if(trim(NAME)== 'NPHS' )  &
             NPHS=wrt_int_state%ALL_DATA_I1D(N2)
            if(trim(NAME)== 'NRDSW' )  &
             TRDSW=wrt_int_state%ALL_DATA_I1D(N2)
            if(trim(NAME)== 'NRDLW' )  &
             TRDLW=wrt_int_state%ALL_DATA_I1D(N2)
            if(trim(NAME)== 'NPREC' )  &
             TPREC=wrt_int_state%ALL_DATA_I1D(N2)
            if(trim(NAME)== 'NHEAT' )  &
             THEAT=wrt_int_state%ALL_DATA_I1D(N2)
            if(trim(NAME)== 'NCLOD' )  &
             TCLOD=wrt_int_state%ALL_DATA_I1D(N2)
            if(trim(NAME)== 'NSRFC' )  &
             TSRFC=wrt_int_state%ALL_DATA_I1D(N2)
            if(trim(NAME)== 'IDAT' .and.N1<=size(SDAT)) &
             SDAT(N1)=wrt_int_state%ALL_DATA_I1D(N2)              !<-- Extract the individual data from the data string
          ENDDO
!
       enddo
!SDAT order needs to change
       I=SDAT(1)
       SDAT(1)=SDAT(2)
       SDAT(2)=I
       print *,'imp_physics=',imp_physics,'ihrst=',ihrst,'NPHS=',NPHS, &
         'TRDLW=',TRDLW,'TRDSW=',TRDSW,'TPREC=',TPREC,'THEAT=',THEAT,  &
         'TCLOD=',TCLOD,'TSRFC=',TSRFC,'SDAT=',SDAT,'iSF_SURFACE_PHYSICS=', &
         iSF_SURFACE_PHYSICS
!
!-----------------------------------------------------------------------
!***  REAL SCALAR/1D HISTORY VARIABLES
!-----------------------------------------------------------------------
!
        N2=0                                                               !<-- Word counter for full string of real scalar/1D data
        dt=spval
        dyval=spval
        dxval=spval
        cenlon=spval
        cenlon=spval
        if(.not.allocated(ETA1)) allocate(ETA1(LM+1))
        if(.not.allocated(ETA2)) allocate(ETA2(LM+1))
        if(.not.allocated(DXH)) allocate(DXH(JM))
!
        DO N=1,wrt_int_state%KOUNT_R1D(1)                                  !<-- Loop through all scalar/1D real data
!
          NPOSN_1=(N-1)*ESMF_MAXSTR+1
          NPOSN_2=N*ESMF_MAXSTR
          NAME=wrt_int_state%NAMES_R1D_STRING(NPOSN_1:NPOSN_2)             !<-- The variable's name
          LENGTH=wrt_int_state%LENGTH_DATA_R1D(N)                          !<-- The variable's length
!          print *,'R1D,NAME=',NAME,'LENGTH=',LENGTH,'N=',N,wrt_int_state%KOUNT_R1D(1)
!
          DO N1=1,LENGTH
            N2=N2+1
            if(trim(NAME)== 'DT' )  &
             DT=wrt_int_state%ALL_DATA_R1D(N2)
            if(trim(NAME)== 'PT' )  &
             PT=wrt_int_state%ALL_DATA_R1D(N2)
            if(trim(NAME)== 'PDTOP' )  &
             PDTOP=wrt_int_state%ALL_DATA_R1D(N2)
            if(trim(NAME)== 'DYH' )  &
             DY(:,:)=wrt_int_state%ALL_DATA_R1D(N2)
            if(trim(NAME)== 'DXH' )  &
             DXH(N1)=wrt_int_state%ALL_DATA_R1D(N2)
            if(trim(NAME)== 'DPHD' )  &
             dyval=wrt_int_state%ALL_DATA_R1D(N2)*gdsdegr
            if(trim(NAME)== 'DLMD' )  &
             dxval=wrt_int_state%ALL_DATA_R1D(N2)*gdsdegr
            if(trim(NAME)== 'TPH0D' )  &
             cenlat=nint(wrt_int_state%ALL_DATA_R1D(N2)*gdsdegr)
            if(trim(NAME)== 'TLM0D' )  &
             cenlon=nint(wrt_int_state%ALL_DATA_R1D(N2)*gdsdegr)
            if(trim(NAME)== 'SLDPTH' )  &
             SLDPTH(N1)=wrt_int_state%ALL_DATA_R1D(N2)
            if(trim(NAME)== 'SG1' )  &
             ETA1(N1)=wrt_int_state%ALL_DATA_R1D(N2)
            if(trim(NAME)== 'SG2' )  &
             ETA2(N1)=wrt_int_state%ALL_DATA_R1D(N2)
          ENDDO
!
        ENDDO
!
!flux are averaged, so set:
        DTQ2=dt*NPHS
        TSRFC=TSRFC*dt/3600.
        IF(TSRFC.EQ.0)TSRFC=float(ifhr)  !in case buket does not get emptied
        TRDLW=TRDLW*dt/3600.
        IF(TRDLW.EQ.0)TRDLW=float(ifhr)  !in case buket does not get emptied
        TRDSW=TRDSW*dt/3600.
        IF(TRDSW.EQ.0)TRDSW=float(ifhr)  !in case buket does not get emptied
        THEAT=THEAT*dt/3600.
        IF(THEAT.EQ.0)THEAT=float(ifhr)  !in case buket does not get emptied
        TCLOD=TCLOD*dt/3600.
        IF(TCLOD.EQ.0)TCLOD=float(ifhr)  !in case buket does not get emptied
        TPREC=TPREC*dt/3600.
        IF(TPREC.EQ.0)TPREC=float(ifhr)  !in case buket does not get emptied
!
!for dx
        print *,'in set_postvars,jsta=',jsta,'jend=',jend,'im=',im,'jm=',jm
        do J=jsta,jend
         dx(1:im,j)=dxh(j)
        enddo

!-----------------------------------------------------------------------
!*** set up module variables
!-----------------------------------------------------------------------
! 3-D real var:
!
      tmaxmin=1.
      DEGRAD=90./ASIN(1.)
      write(0,*)'name r2d size=',len(wrt_int_state%NAMES_R2D_STRING)
      allocate(vlat(1:im,jsta_2l:jend_2u),vlon(1:im,jsta_2l:jend_2u))
      allocate(buf(1:im,jsta_2l:jend_2u))
      allocate(dummy(im,jm))
      field_loop_real: DO NFIELD=1,wrt_int_state%KOUNT_R2D(1)
!
        NPOS_START=(NFIELD-1)*ESMF_MAXSTR+1
        NPOS_END=NFIELD*ESMF_MAXSTR
!        print *,'NPOS_START=',NPOS_START,'NPOS_END=',NPOS_END
        NAME=wrt_int_state%NAMES_R2D_STRING(NPOS_START:NPOS_END)
        INDX_2D=index(NAME,"_2D")
!        print *,'in set_postvars,nfield=',nfield,'name=',trim(NAME), &
!          wrt_int_state%WRITE_SUBSET_R(1:2,jsta:jsta+2,NFIELD),  &
!          maxval(wrt_int_state%WRITE_SUBSET_R(1:im,jsta:jend,NFIELD)), &
!          minval(wrt_int_state%WRITE_SUBSET_R(1:im,jsta:jend,NFIELD))

        if (INDX_2D.gt.0) then
           model_level=name(indx_2D-3:indx_2D-1)
           LL=(ichar(model_level(1:1))-48)*100+(ichar(model_level(2:2))-48)*10+ichar(model_level(3:3))-48
           if(name(1:INDX_2d-5).eq.'CLDFRA') then
             do j=jsta,jend
              cfr(:,j,LL)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(name(1:INDX_2d-5).eq.'CW') then
             do j=jsta,jend
              cwm(:,j,LL)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(name(1:INDX_2d-5).eq.'EXCH_H') then
             do j=jsta,jend
              exch_h(:,j,LL)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(name(1:INDX_2d-5).eq.'Q') then
             do j=jsta,jend
              Q(:,j,LL)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
!             print *,'in set_postvars Q=',maxval(Q(1:im,jsta:jend,LL)),   &
!               minval(Q(1:im,jsta:jend,LL)),  &
!               wrt_int_state%WRITE_SUBSET_R(1:5,jsta+3,NFIELD)
           endif
           if(name(1:INDX_2d-5).eq.'Q2') then
             do j=jsta,jend
              Q2(:,j,LL)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
!             print *,'in set_postvars Q2=',maxval(Q2(1:im,jsta:jend,LL)),   &
!               minval(Q2(1:im,jsta:jend,LL)),  &
!               wrt_int_state%WRITE_SUBSET_R(1:5,jsta+3,NFIELD)
           endif
           if(name(1:INDX_2d-5).eq.'PINT') then
             do j=jsta,jend
              PINT(:,j,LL)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
             if(ll/=1) then
               do j=jsta,jend
                 alpint(:,j,ll)=alog(pint(:,j,ll))
               enddo
             elseif(ll==1) then
               do j=jsta,jend
               do i=1,im
                 if(pint(i,j,ll)/=0.) then
                    alpint(:,j,ll)=alog(pint(:,j,ll))
                 else
                    alpint(I,J,ll)=spval
                 endif
               enddo
               enddo
             endif
             if(ll.gt.1) then
               do j = jsta, jend
                 PMID(:,j,ll-1 ) = (PINT(:,J,ll-1)+                              &
                     PINT(:,J,ll))*0.5 ! representative of what model does
               end do
!               print *,'in set_postvars pmid=',maxval(pmid(1:im,jsta:jend,ll-1)), &
!               minval(pmid(1:im,jsta:jend,ll-1))

             endif

!             print *,'in set_postvars,ll=',ll,'pint=',maxval(pint(1:im,jsta:jend,LL)),   &
!               minval(pint(1:im,jsta:jend,LL))
           endif
           if(name(1:INDX_2d-5).eq.'RLWTT') then
             do j=jsta,jend
              RLWTT(:,j,LL)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(name(1:INDX_2d-5).eq.'RSWTT') then
             do j=jsta,jend
              RSWTT(:,j,LL)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(name(1:INDX_2d-5).eq.'T') then
             do j=jsta,jend
              T(:,j,LL)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
!             print *,'in set_postvars T=',maxval(T(1:im,jsta:jend,LL)),   &
!               minval(T(1:im,jsta:jend,LL))
           endif
           if(name(1:INDX_2d-5).eq.'TCUCN') then
             do j=jsta,jend
              TCUCN(:,j,LL)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(name(1:INDX_2d-5).eq.'TRAIN') then
             do j=jsta,jend
              TRAIN(:,j,LL)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(name(1:INDX_2d-5).eq.'U') then
             do j=jsta,jend
              UH(:,j,LL)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
              U(:,j,LL)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
! put u on h point for global nmm
             if(global)then
             buf(:,:)=uh(:,:,ll)
             call exch(buf(1,jsta_2l))
              do j=jsta,jend
               do i=1,im
                im1=i-1
                if(im1<1)im1=im1+im
                jm1=j-1
                if(j==1)then
                 ii=i+im/2
                 iim1=ii-1
                 if(iim1<1)iim1=iim1+im
                 if (ii > im) ii = ii - im
                 uh(i,j,ll)=(buf(i,j)+buf(im1,j)+buf(ii,j)+buf(iim1,j))/4.0
                else
                 uh(i,j,ll)=(buf(i,j)+buf(im1,j)+buf(im1,jm1)+buf(i,jm1))/4.0
                end if
               end do
              end do
             end if ! end of wind interpolation for global NMM
!             print *,'in set_postvars uh=',maxval(uh(1:im,jsta:jend,LL)),   &
!               minval(uh(1:im,jsta:jend,LL)),  &
!               wrt_int_state%WRITE_SUBSET_R(1:5,jsta+3,NFIELD)
           endif
           if(name(1:INDX_2d-5).eq.'V') then
             do j=jsta,jend
              VH(:,j,LL)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
              V(:,j,LL)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
! put u on h point for global nmm
             if(global)then
             buf(:,:)=vh(:,:,ll)
             call exch(buf(1,jsta_2l))
              do j=jsta,jend
               do i=1,im
                im1=i-1
                if(im1<1)im1=im1+im
                jm1=j-1
                if(j==1)then
                 ii=i+im/2
                 iim1=ii-1
                 if(iim1<1)iim1=iim1+im
                 if (ii > im) ii = ii - im
                 vh(i,j,ll)=(buf(i,j)+buf(im1,j)+buf(ii,j)+buf(iim1,j))/4.0
                else
                 vh(i,j,ll)=(buf(i,j)+buf(im1,j)+buf(im1,jm1)+buf(i,jm1))/4.0
                end if
               end do
              end do
             end if ! end of wind interpolation for global NMM
!             print *,'in set_postvars vh=',maxval(vh(1:im,jsta:jend,LL)),   &
!               minval(vh(1:im,jsta:jend,LL)),  &
!               wrt_int_state%WRITE_SUBSET_R(1:5,jsta+3,NFIELD)
           endif
           if(name(1:INDX_2d-5).eq.'W') then
             do j=jsta,jend
              WH(:,j,LL)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
!             print *,'in set_postvars vh=',maxval(vh(1:im,jsta:jend,LL)),   &
!               minval(vh(1:im,jsta:jend,LL)),  &
!               wrt_int_state%WRITE_SUBSET_R(1:5,jsta+3,NFIELD)
           endif

           if(name(1:INDX_2d-5).eq.'XLEN_MIX') then
             do j=jsta,jend
!             EL_MYJ(:,j,LL)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
              EL_PBL(:,j,LL)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(name(1:INDX_2d-5).eq.'F_ICE') then
             do j=jsta,jend
              F_ICE(:,j,LL)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(name(1:INDX_2d-5).eq.'F_RIMEF') then
             do j=jsta,jend
              F_RimeF(:,j,LL)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
!            print *,'ll=',ll,'f_rimef=',maxval(f_rimef(1:im,jsta:jend,ll)), &
!            minval(f_rimef(1:im,jsta:jend,ll))
           endif
           if(name(1:INDX_2d-5).eq.'F_RAIN') then
             do j=jsta,jend
              F_RAIN(:,j,LL)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(name(1:INDX_2d-5).eq.'SH2O') then
             do j=jsta,jend
              sh2o(:,j,LL)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(name(1:INDX_2d-5).eq.'SMC') then
             do j=jsta,jend
              smc(:,j,LL)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(name(1:INDX_2d-5).eq.'STC') then
             do j=jsta,jend
              stc(:,j,LL)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
!
        else

           if(trim(name).eq.'GLAT') then
             do j=jsta,jend
              gdlat(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)*DEGRAD
             enddo
!            print *,'in ste_post,gdlat=',maxval(gdlat(1:im,jstA:jend)),minval(gdlat(1:im,jsta:jend))
           endif
           if(trim(name).eq.'GLON') then
             do j=jsta,jend
              gdlon(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)*DEGRAD
             enddo
!            print *,'in ste_post,gdlon=',maxval(gdlon(1:im,jstA:jend)),minval(gdlon(1:im,jsta:jend))
           endif
           if(trim(name).eq.'PD') then
             do j=jsta,jend
              pd(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
!            print *,'in ste_post,pd=',maxval(pd(1:im,jstA:jend)),minval(pd(1:im,jsta:jend))
           endif
           if(trim(name).eq.'VLAT') then
             do j=jsta,jend
              vlat(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)*DEGRAD
             enddo
           endif
           if(trim(name).eq.'VLON') then
             do j=jsta,jend
              vlon(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)*DEGRAD
             enddo
           endif
           if(trim(name).eq.'ACFRCV') then
             do j=jsta,jend
              ACFRCV(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'ACFRST') then
             do j=jsta,jend
              ACFRST(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'ACPREC') then
             do j=jsta,jend
              ACPREC(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'ACSNOM') then
             do j=jsta,jend
              ACSNOM(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'ACSNOW') then
             do j=jsta,jend
              ACSNOW(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'AKHS_OUT') then
             do j=jsta,jend
              AKHS(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'AKMS_OUT') then
             do j=jsta,jend
              AKMS(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'ALBASE') then
             do j=jsta,jend
              ALBASE(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'ALBEDO') then
             do j=jsta,jend
              ALBEDO(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'ALWIN') then
             do j=jsta,jend
              ALWIN(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'ALWOUT') then
             do j=jsta,jend
              ALWOUT(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'ALWTOA') then
             do j=jsta,jend
              ALWTOA(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'ASWIN') then
             do j=jsta,jend
              ASWIN(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'ASWOUT') then
             do j=jsta,jend
              ASWOUT(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'ASWTOA') then
             do j=jsta,jend
              ASWTOA(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'AVRAIN') then
             AVRAIN=wrt_int_state%WRITE_SUBSET_R(im/2,(jsta+jend)/2,NFIELD)
           endif
           if(trim(name).eq.'AVCNVC') then
             AVCNVC=wrt_int_state%WRITE_SUBSET_R(im/2,(jsta+jend)/2,NFIELD)
           endif
           if(trim(name).eq.'ARDLW') then
             ARDLW=wrt_int_state%WRITE_SUBSET_R(im/2,(jsta+jend)/2,NFIELD)
           endif
           if(trim(name).eq.'ARDSW') then
             ARDSW=wrt_int_state%WRITE_SUBSET_R(im/2,(jsta+jend)/2,NFIELD)
           endif
           if(trim(name).eq.'ASRFC') then
             ASRFC=wrt_int_state%WRITE_SUBSET_R(im/2,(jsta+jend)/2,NFIELD)
           endif
           if(trim(name).eq.'BGROFF') then
             do j=jsta,jend
              BGROFF(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'CFRACH') then
             do j=jsta,jend
              CFRACH(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'CFRACL') then
             do j=jsta,jend
              CFRACL(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'CFRACM') then
             do j=jsta,jend
              CFRACM(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'CLDEFI') then
             do j=jsta,jend
              CLDEFI(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'CMC') then
             do j=jsta,jend
              CMC(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'CNVBOT') then
             do j=jsta,jend
              HBOT(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
             write(0,*) 'HBOT1=',maxval(HBOT(:,jsta:jend)),minval(HBOT(:,jsta:jend))
             where(HBOT<spval)HBOT=float(lm)-hbot+1.0
             write(0,*) 'HBOT2=',maxval(HBOT(:,jsta:jend)),minval(HBOT(:,jsta:jend))
           endif
           if(trim(name).eq.'CNVTOP') then
             do j=jsta,jend
              HTOP(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
             where(htop < spval)htop=float(lm)-htop+1.0
           endif
           if(trim(name).eq.'CPRATE') then
             do j=jsta,jend
              CPRATE(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
!            print *,'in set_post,cpreate=',maxval(cprate(1:im,jsta:jend)), &
!              minval(cprate(1:im,jsta:jend))
           endif
           if(trim(name).eq.'CUPPT') then
             do j=jsta,jend
              CUPPT(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'CUPREC') then
             do j=jsta,jend
              CUPREC(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'LSPA') then
             do j=jsta,jend
              LSPA(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'CZEN') then
             do j=jsta,jend
              CZEN(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'CZMEAN') then
             do j=jsta,jend
              CZMEAN(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'FIS') then
             do j=jsta,jend
              FIS(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
!to be consistent with post
             FIS=FIS/G1
!       print *,'in init_nems,1,fis=',maxval(fis(1:im,jsta:jend)), &
!         minval(fis(1:im,jsta:jend))
             FIS=FIS*G
!       print *,'in init_nems,fis=',maxval(fis(1:im,jsta:jend)), &
!         minval(fis(1:im,jsta:jend))
           endif
           if(trim(name).eq.'GRNFLX') then
             do j=jsta,jend
              GRNFLX(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'PCTSNO') then
             do j=jsta,jend
              PCTSNO(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'HBOTD') then
             do j=jsta,jend
              HBOTD(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
             where(hbotd < spval)hbotd=float(lm)-hbotd+1.0
           endif
           if(trim(name).eq.'HBOTS') then
             do j=jsta,jend
              HBOTS(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
             where(hbots < spval)hbots=float(lm)-hbots+1.0
           endif
           if(trim(name).eq.'HTOPD') then
             do j=jsta,jend
              HTOPD(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
             where(htopd < spval)htopd=float(lm)-htopd+1.0
           endif
           if(trim(name).eq.'HTOPS') then
             do j=jsta,jend
              HTOPS(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
             where(htops < spval)htops=float(lm)-htops+1.0
           endif
           if(trim(name).eq.'MXSNAL') then
             do j=jsta,jend
              MXSNAL(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'PBLH') then
             do j=jsta,jend
              PBLH(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'MIXHT') then
             do j=jsta,jend
              MIXHT(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'PREC') then
             do j=jsta,jend
              PREC(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'PSHLTR') then
             do j=jsta,jend
              PSHLTR(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'POTEVP') then
             do j=jsta,jend
              POTEVP(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'Q10') then
             do j=jsta,jend
              Q10(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'QSH') then
             do j=jsta,jend
              QS(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'QSHLTR') then
             do j=jsta,jend
              QSHLTR(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'T02MAX') then
             do j=jsta,jend
              MAXTSHLTR(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'T02MIN') then
             do j=jsta,jend
              MINTSHLTR(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'RH02MAX') then
             do j=jsta,jend
              MAXRHSHLTR(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'RH02MIN') then
             do j=jsta,jend
              MINRHSHLTR(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'QWBS') then
             do j=jsta,jend
              QWBS(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'QZ0') then
             do j=jsta,jend
              QZ0(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'RADOT') then
             do j=jsta,jend
              RADOT(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'RLWIN') then
             do j=jsta,jend
              RLWIN(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'RLWTOA') then
             do j=jsta,jend
              RLWTOA(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'RSWIN') then
             do j=jsta,jend
              RSWIN(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
!             print *,'in post_nmm, set postvar,rswin=',maxval(rswin(:,jsta:jend)), &
!              minval(rswin(:,jsta:jend))
           endif
           if(trim(name).eq.'RSWINC') then
             do j=jsta,jend
              RSWINC(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'RSWOUT') then
             do j=jsta,jend
              RSWOUT(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'SFCEVP') then
             do j=jsta,jend
              SFCEVP(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'SFCEXC') then
             do j=jsta,jend
              SFCEXC(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'SFCLHX') then
             do j=jsta,jend
              SFCLHX(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'SFCSHX') then
             do j=jsta,jend
              SFCSHX(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'SI') then
             do j=jsta,jend
              SI(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'SICE') then
             do j=jsta,jend
              SICE(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'SIGT4') then
             do j=jsta,jend
              SIGT4(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'SM') then
             do j=jsta,jend
              SM(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'SMSTAV') then
             do j=jsta,jend
              SMSTAV(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'SMSTOT') then
             do j=jsta,jend
              SMSTOT(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'SNO') then
             do j=jsta,jend
              SNO(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'SNOAVG') then
             do j=jsta,jend
              SNOAVG(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'PSFCAVG') then
             do j=jsta,jend
              PSFCAVG(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'T10AVG') then
             do j=jsta,jend
              T10AVG(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'T10') then
             do j=jsta,jend
              T10M(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'AKHSAVG') then
             do j=jsta,jend
              AKHSAVG(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'AKMSAVG') then
             do j=jsta,jend
              AKMSAVG(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'REFDMAX') then
             do j=jsta,jend
              REFD_MAX(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'UPVVELMAX') then
             do j=jsta,jend
              W_UP_MAX(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'DNVVELMAX') then
             do j=jsta,jend
              W_DN_MAX(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'UPHLMAX') then
             do j=jsta,jend
              UP_HELI_MAX(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'SNOPCX') then
             do j=jsta,jend
              SNOPCX(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'SFCUVX') then
             do j=jsta,jend
              SFCUVX(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'SOILTB') then
             do j=jsta,jend
              SOILTB(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'SR') then
             do j=jsta,jend
              SR(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'SSROFF') then
             do j=jsta,jend
              SSROFF(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'SST') then
             do j=jsta,jend
              SST(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
!             print *,'in set_post,sst=',maxval(sst(1:im,jsta:jend)),minval(sst(1:im,jsta:jend))
           endif
           if(trim(name).eq.'SUBSHX') then
             do j=jsta,jend
              SUBSHX(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'TG') then
             do j=jsta,jend
              TG(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'TH10') then
             do j=jsta,jend
              TH10(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'THS') then
             do j=jsta,jend
              THS(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'THZ0') then
             do j=jsta,jend
              THZ0(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'TSHLTR') then
             do j=jsta,jend
              TSHLTR(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'TWBS') then
             do j=jsta,jend
              TWBS(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'U10') then
             do j=jsta,jend
              U10(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'U10MAX') then
             do j=jsta,jend
              U10MAX(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'USTAR') then
             do j=jsta,jend
              USTAR(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'UZ0') then
             do j=jsta,jend
              UZ0(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'V10') then
             do j=jsta,jend
              V10(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'V10MAX') then
             do j=jsta,jend
              V10MAX(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'VEGFRC') then
             do j=jsta,jend
              VEGFRC(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'VZ0') then
             do j=jsta,jend
              VZ0(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif
           if(trim(name).eq.'Z0') then
             do j=jsta,jend
              Z0(:,j)=wrt_int_state%WRITE_SUBSET_R(:,j,NFIELD)
             enddo
           endif

!
        endif


      enddo  field_loop_real
!
      write(0,*)'after field loop real'
!get 2D Integer vars
!
      field_loop_int: DO NFIELD=1,wrt_int_state%KOUNT_I2D(1)
!
        NPOS_START=(NFIELD-1)*ESMF_MAXSTR+1
        NPOS_END=NFIELD*ESMF_MAXSTR
        NAME=wrt_int_state%NAMES_I2D_STRING(NPOS_START:NPOS_END)

           if(trim(name).eq.'ISLTYP') then
             do j=jsta,jend
              ISLTYP(:,j)=wrt_int_state%WRITE_SUBSET_I(:,j,NFIELD)
             enddo
!             print *,'sltyp=',maxval(isltyp(1:im,jsta:jend)),      &
!               minval(isltyp(1:im,jsta:jend))
           endif
           if(trim(name).eq.'IVGTYP') then
             do j=jsta,jend
              IVGTYP(:,j)=wrt_int_state%WRITE_SUBSET_I(:,j,NFIELD)
              buf(:,j)=IVGTYP(:,J)
             enddo
             call collect_loc(buf,dummy)
             if(me==0)novegtype=NINT(maxval(dummy))
             call mpi_bcast(novegtype,1,MPI_INTEGER,0,mpi_comm_comp,iret)
             print*,'novegtype= ',novegtype

!             print *,'vgtyp=',maxval(ivgtyp(1:im,jsta:jend)),      &
!               minval(ivgtyp(1:im,jsta:jend))
           endif
           if(trim(name).eq.'NCFRCV') then
             do j=jsta,jend
              NCFRCV(:,j)=float(wrt_int_state%WRITE_SUBSET_I(:,j,NFIELD))
             enddo
           endif
           if(trim(name).eq.'NCFRST') then
             do j=jsta,jend
              NCFRST(:,j)=float(wrt_int_state%WRITE_SUBSET_I(:,j,NFIELD))
             enddo
           endif
!
!
       enddo field_loop_int
      write(0,*)'after field loop int'


      hbm2=1.0
      write(0,*)'after field loop int'
!
!      CALL EXCH(gdlat(1,JSTA_2L))
!      print *,'after call EXCH,mype=',me
!
!      do j = jsta, jend_m
!        do i = 1, im-1
!          DX ( i, j ) = ERAD*COS(GDLAT(I,J)*DTR)                        &
!            *(GDLON(I+1,J)-GDLON(I,J))*DTR
!          DY ( i, j ) =  ERAD*(GDLAT(I,J)-GDLAT(I,J+1))*DTR  ! like A*DPH
!        end do
!      end do
!
      do j=jsta,jend
      do i=1,im
        F(I,J)=1.454441e-4*sin(gdlat(i,j)*DTR)   ! 2*omeg*sin(phi)
      end do
      end do
!
      call collect_loc(gdlat,dummy)
      if(me.eq.0)then
        ii=(1+im)/2
        jj=(1+jm)/2
        latstart=nint(dummy(1,1)*gdsdegr)
        latlast=nint(dummy(im,jm)*gdsdegr)
        latnw=nint(dummy(1,jm)*gdsdegr)
        latse=nint(dummy(im,1)*gdsdegr)
      end if
      call mpi_bcast(latstart,1,MPI_INTEGER,0,mpi_comm_comp,iret)
      call mpi_bcast(latlast,1,MPI_INTEGER,0,mpi_comm_comp,iret)
!
      if(global)then
       if(gdlon(1,jsta)>0. .and. gdlon(2,jsta)<0.)then
        do j=jsta,jend
         gdlon(1,j)=gdlon(1,j)-360.0
        end do
       end if
       maptype=0 !  for global NMMB on latlon grid
       gridtype='A' ! will put wind on mass point for now to make regular latlon
      end if
      call collect_loc(gdlon,dummy)
      if(me.eq.0)then
        lonstart=nint(dummy(1,1)*1000.)
        lonlast=nint(dummy(im,jm)*1000.)
        lonnw=nint(dummy(1,jm)*1000.)
        lonse=nint(dummy(im,1)*1000.)
      endif
      call mpi_bcast(lonstart,1,MPI_INTEGER,0,mpi_comm_comp,iret)
      call mpi_bcast(lonlast,1,MPI_INTEGER,0,mpi_comm_comp,iret)
!
      call collect_loc(vlat,dummy)
      if(me.eq.0)then
        ii=(1+im)/2
        jj=(1+jm)/2
        latstartv=nint(dummy(1,1)*1000.)
        latlastv=nint(dummy(im,jm)*1000.)
      end if
      call mpi_bcast(latstartv,1,MPI_INTEGER,0,mpi_comm_comp,iret)
      call mpi_bcast(latlastv,1,MPI_INTEGER,0,mpi_comm_comp,iret)
      cenlatv=cenlat
      write(6,*) 'latstartv,cenlatv,latlastv,me A calling bcast=', &
      latstartv,cenlatv,latlastv,me
!
      call collect_loc(vlon,dummy)
      if(me.eq.0)then
        ii=(1+im)/2
        jj=(1+jm)/2
        lonstartv=nint(dummy(1,1)*1000.)
        lonlastv=nint(dummy(im,jm)*1000.)
      end if
      call mpi_bcast(lonstartv,1,MPI_INTEGER,0,mpi_comm_comp,iret)
      call mpi_bcast(lonlastv,1,MPI_INTEGER,0,mpi_comm_comp,iret)
      cenlonv=cenlon
      write(6,*) 'lonstartv,cenlonv,lonlastv,me A calling bcast=', &
      lonstartv,cenlonv,lonlastv,me
      deallocate(vlat,vlon)
!
!get omga
!
      do l = 1, lm
       do j = jsta, jend
        do i = 1, im
            IF(ABS(T(I,J,L)).GT.1.0E-3)                                &
              OMGA(I,J,L) = -WH(I,J,L)*PMID(I,J,L)*G/                   &
                       (RD*T(I,J,L)*(1.+D608*Q(I,J,L)))

        end do
       end do
      end do
      write(0,*)' after OMGA'
!
!get ttnd
!
      ttnd=0.
      where(rlwtt/=spval .and. rswtt/=spval)ttnd=rswtt+rlwtt

!
      if(gridtype=='E')then
       do l = 1, lm
        call exch(PMID(1:IM,JSTA_2L:JEND_2U,L))
        do j = jsta, jend
         do i = 1, im-MOD(J,2)
          IF(J .EQ. 1 .AND. I .LT. IM)THEN   !SOUTHERN BC
           PMIDV(I,J,L)=0.5*(PMID(I,J,L)+PMID(I+1,J,L))
          ELSE IF(J.EQ.JM .AND. I.LT.IM)THEN   !NORTHERN BC
           PMIDV(I,J,L)=0.5*(PMID(I,J,L)+PMID(I+1,J,L))
          ELSE IF(I .EQ. 1 .AND. MOD(J,2) .EQ. 0) THEN   !WESTERN EVEN BC
           PMIDV(I,J,L)=0.5*(PMID(I,J-1,L)+PMID(I,J+1,L))
          ELSE IF(I .EQ. IM .AND. MOD(J,2) .EQ. 0                             &
          .AND. J .LT. JM) THEN   !EASTERN EVEN BC
           PMIDV(I,J,L)=0.5*(PMID(I,J-1,L)+PMID(I,J+1,L))
          ELSE IF (MOD(J,2) .LT. 1) THEN
           PMIDV(I,J,L)=0.25*(PMID(I,J,L)+PMID(I-1,J,L)                       &
             +PMID(I,J+1,L)+PMID(I,J-1,L))
          ELSE
           PMIDV(I,J,L)=0.25*(PMID(I,J,L)+PMID(I+1,J,L)                       &
             +PMID(I,J+1,L)+PMID(I,J-1,L))
          END IF
         end do
        end do
       end do
      else if(gridtype=='B')then
       do l = 1, lm
        call exch(PMID(1:IM,JSTA_2L:JEND_2U,L))
        do j = jsta, jend_m
         do i = 1, im-1
           PMIDV(I,J,L)=0.25*(PMID(I,J,L)+PMID(I+1,J,L)                       &
             +PMID(I,J+1,L)+PMID(I+1,J+1,L))
         end do
        end do
       end do
      end if
      write(0,*)' after PMIDV'

!
!!!!! COMPUTE Z
       allocate(FI(IM,JM,2))
       do j = jsta, jend
        do i = 1, im
            ZINT(I,J,LM+1)=FIS(I,J)/G
        if (I .eq. im/2 .and. J .eq.(jsta+jend)/2 ) then
                   write(6,*) 'G,ZINT: ', G,ZINT(I,J,LM+1)
        endif
            FI(I,J,1)=FIS(I,J)
        end do
       end do
!       print *,'in init_nems,zint=',maxval(zint(1:im,jsta:jend,lm+1)), &
!         minval(zint(1:im,jsta:jend,lm+1))
!
! SECOND, INTEGRATE HEIGHT HYDROSTATICLY
      DO L=LM,1,-1
       do j = jsta, jend
        do i = 1, im
         FI(I,J,2)=HTM(I,J,L)*T(I,J,L)*(Q(I,J,L)*D608+1.0)*RD*                &
                   (ALPINT(I,J,L+1)-ALPINT(I,J,L))+FI(I,J,1)
         ZINT(I,J,L)=FI(I,J,2)/G
!         if(i==im/2.and.j==(jsta+jend)/2)                                              &
!        print*,'L,sample HTM,T,Q,ALPINT(L+1),ALPINT(l),ZINT= '                &
!        ,l,HTM(I,J,L),T(I,J,L),Q(I,J,L),ALPINT(I,J,L+1),                      &
!        ALPINT(I,J,L),ZINT(I,J,L)
         FI(I,J,1)=FI(I,J,2)
        ENDDO
       ENDDO
!       print *,'in set_post,l=',l,'zint=',maxval(zint(1:im,jsta:jend,l)), &
!         minval(zint(1:im,jsta:jend,l))
      END DO
      deallocate(FI)
!
      DO L=1,LM
       do J = JSTA, JEND
          DO I=1,IM
            FACT=(ALOG(PMID(I,J,L))-ALOG(PINT(I,J,L)))/                      &
               (ALOG(PINT(I,J,L+1))-ALOG(PINT(I,J,L)))
            ZMID(I,J,L)=ZINT(I,J,L)+(ZINT(I,J,L+1)-ZINT(I,J,L))*FACT
          ENDDO
        ENDDO
      ENDDO
!
!get tables
       THL=210.
       PLQ=70000.
       CALL TABLE(PTBL,TTBL,PT,                                       &
                RDQ,RDTH,RDP,RDTHE,PL,THL,QS0,SQS,STHE,THE0)

       CALL TABLEQ(TTBLQ,RDPQ,RDTHEQ,PLQ,THL,STHEQ,THE0Q)
      write(0,*)' after TABLEQ, end of set_postvars'

!
    end subroutine set_postvars_nmm
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
