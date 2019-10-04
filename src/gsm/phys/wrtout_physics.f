
      subroutine wrtout_physics(phour,fhour,zhour,idate,
     &                  sl,si,
     &                  sfc_fld, flx_fld, nst_fld, g2d_fld,
     &                  fluxr,
     &                  global_lats_r,lonsperlar,nblck,
     &                  colat1,cfhour1,pl_coeff,sfcf,nstf,flxf,d3df)
!!
!
! May 2009 Jun Wang, modified to use write grid component
! Jan 2010 Sarah Lu, AOD added to flx files
! Feb 2010 Jun Wang, write out restart file
! Jul 2010 S. Moorthi - added nst and other modifications
! Jul 2010 S. Moorthi - added  hchuang  Add flx files output to wrtflx_a
! Jul 2010 Sarah Lu, write out aerosol diag files (for g2d_fld)
! Aug 2010 Sarah Lu, scale the 2d_aer_diag by 1.e6
!                    output time-avg 2d_aer_diag
! Oct 2010 Sarah Lu, add g2d_fld%met
! Oct 2010 Sarah Lu, g2d_fld%met changed from instant to accumulated
! Dec 2010 Sarah Lu, g2d_fld%met contains both instant and time-avg;
!                    wrtaer is called only when gocart is on
! Dec 2010 Jun Wang, change to nemsio library
! Nov 2012 Jun Wang, removing quilting, which is not used
! Nov 2012 Jun Wang, add sfcio opt
! Nov 2012 Jun Wang, add d3d opt
! Nov 2012 Jun Wang, removing quilting, which is not used
! Jan 2013 S. Moorthi, adding sfcf,flxf,d3df to call and related changes
! Feb 2013 Jun Wang, using gribit_gsm from gsm (gridit_gsm)
! May 2013 S. Moorthi, using gribit_gsm from gsm (gridit_gsm)
! Aug 2013 S. Moorthi  Merging with the GFS trunk wrtout(updated with rrtm )
!                      and some cleanup
! Nov 2013 S. Lu       Modifying wrtflx_a and wrtflx_w
! Aug 2015 Xu Li       change nst_fcst to be nstf_name
! Jan 2016 Xu Li       Add nst_only_move, nst_wrt and nst_wrt_nemsio
! Jan 2016 S. Moorthi  bug fix in cloud top temp and pressure diagnostic field
!
! Original wrtout history:
!   program history log:                                               !
!     mmm-yyyy  g. vandengerghe - created program wrtout               !
!      -   --   joe sela     - modified, set lfnhr to false for        !
!                 writing one step output etc.                         !
!      -   --   sarah lu     - added smc(3:4), stc(3:4), slc(1:4),     !
!                 snwdph, canopy, changed 10-200cm to 10-40cm for      !
!                 smc(2), and changed 10-200 to 10-40 for stc(2).      !
!      -   --   jun wang     - modified to add spfhmax/spfhmin.        !
!     nov 2004  xingren wu   - modified to add sea-ice fields.         !
!     oct 2006  helin wei    - modified to add 30 records for land mdl.!
!                 monitoring and analysis purpose.                     !
!     nov 2007  ho-chun huang- code change for gocart, added lggfs3d   !
!                 and wrt3gd_hyb for gocart only g3d files output.     !
!                 also confirmed by helin, correct bug, zorl unit in   !
!                 cm not mm, when callng subuninterprez, array glolal  !
!                 is assign to buff_mult_piecea at the ngrid location, !
!                 then ngrid advanced by 1.  before assign the modified!
!                 value (buffo) to buff_mult_piecea again dial ngrid   !
!                 back by 1 for the correct ngrid index.               !
!      -   --   ho-chun huang- code change add goro lan-sea output     !
!     sep 2008  yu-tai hou   - add sunshine duration time to flux file.!
!         2009  sarah lu     - added 7 clear-sky radiation racords.    !
!      -   --   hc huang     - added control flag wrt_g3d for output   !
!                 g3d calls.                                           !
!      -   --   s. moorthi   - a multi-year continuation to improve    !
!                 the support of requests for adding output fields and !
!                 model upgrades, including changing fields, logical   !
!                 controls, grib data conversions, and incorporating   !
!                 new developments, etc.                               !
!     apr 2012  yu-tai hou   - added 4 sw sfc flux components vis/nir  !
!                 beam/diffused. in subprograms wrtflx_a and wrtflx_w, !
!                 the written out fields have been re-organized that   !
!                 related fields are closer together, and reduced      !
!                 duplications of field labels.                        !
!     jan 2013  s. moorthi   - modified scale factor for spfhmax and   !
!                 spfhmin fields with additional variable ids_iq to    !
!                 avoid conflict with fields dnwd sw/lw (idswf/idlwf). !
!     mar 2013  yu-tai hou   - modified scale factor for sfc-uvb and   !
!                 sfc-csuvb fields with additional variable ids_uvb to !
!                 avoid conflict with field canopy water evap (ievcw). !
!                                                                      !
!!!!!  ==========================================================  !!!!!
!!!!!                       end descriptions                       !!!!!
!!!!!  ==========================================================  !!!!!

      use resol_def,               ONLY: latr, levs, levp1, lonr, nfxr, &
     &                                   ngrids_aer
      use layout1,                 ONLY: me, nodes, lats_node_r,        &
     &                                   nodes_comp
      use namelist_physics_def,    ONLY: gen_coord_hybrid, ldiag3d, 
     &                                   hybrid, fhlwr, fhswr, ens_nam,
     &                                   nstf_name, sfcio_out,iau
      use mpi_def,                 ONLY: liope, info, mpi_comm_all, 
     &                                   mc_comp, mpi_comm_null
      use gfs_physics_sfc_flx_mod, ONLY: Sfc_Var_Data, Flx_Var_Data
      use gfs_physics_nst_var_mod, ONLY: Nst_Var_Data
      use gfs_physics_g2d_mod,     ONLY: G2D_Var_Data
      USE machine,                 ONLY: kind_evod, kind_io8
      implicit none
!!
      TYPE(Sfc_Var_Data)        :: sfc_fld
      TYPE(Flx_Var_Data)        :: flx_fld
      TYPE(Nst_Var_Data)        :: nst_fld
      TYPE(G2D_Var_Data)        :: g2d_fld
      CHARACTER(16)             :: CFHOUR1    ! for ESMF Export State Creation

      character*5          sfcf, nstf, flxf, d3df
      integer              pl_coeff, idate(4), nblck
      real(kind=kind_evod) phour,fhour,zhour
!     real(kind=kind_evod) phour,fhour,zhour, xgf

      real(kind=kind_evod) sl(levs), si(levp1)
      INTEGER              GLOBAL_LATS_R(LATR),   lonsperlar(LATR)
      REAL (KIND=kind_io8) fluxr(nfxr,LONR,LATS_NODE_R)
!!
      logical     lfnhr
      real        colat1
      real(kind=8) t1,t2,t3,t4,t5,ta,tb,tc,td,te,tf,rtc,tx,ty
      real(kind=8) tba,tbb,tbc,tbd
      integer      km,iostat,no3d,ks
      integer      ierr,j,k,l,lenrec,locl,n,node,iret
      integer      nosfc,noflx,nonst,noaer,nfill
      character*16 cosfc,const
      real         timesum
      data         timesum/0./
      integer nfillnum
!!
!!
      character CFHOUR*40,CFORM*40
      integer   jdate(4),ndigyr,ndig,kh,IOPROC
!  ---  locals for iau:
      integer :: jdate_iau(4),idat(8),jdat(8)
      real (kind=kind_evod) :: zhour_iau,fhour_iau
      real :: rinc(5)
!!
      real(kind=kind_io8) slmsk(lonr,lats_node_r)
!!
      real(kind=kind_evod) secphy,secswr,seclwr
!
!-------------------------------------------------------------------------
!     print *,' in wrtout_phyiscs me=',me
      t3  = rtc()
      call mpi_barrier(mpi_comm_all,ierr)
      t4  = rtc()
      tba = t4-t3
      if(nodes_comp < 1 .or. nodes_comp > nodes) then
        print *, '  NODES_COMP UNDEFINED, CANNOT DO I.O '
        call mpi_finalize()
        stop 333
      endif
!
!     ioproc = nodes_comp-1
      ioproc = 0
       
      t1 = rtc()
!!
!     CREATE CFHOUR

      JDATE  = IDATE
      ndigyr = 4
      IF(NDIGYR == 2) THEN
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
!     lfnhr = 3600*abs(fhour-nint(fhour)).le.1.or.phour.eq.0
      lfnhr = 3600*abs(fhour_iau-nint(fhour_iau)) <= 1
      IF (LFNHR) THEN
        KH   = NINT(FHOUR)
        NDIG = MAX(LOG10(KH+0.5)+1.,2.)
        WRITE(CFORM,'("(I",I1,".",I1,")")') NDIG,NDIG
        WRITE(CFHOUR,CFORM) KH
      ELSE
        KS   = NINT(fhour_iau*3600)
        KH   = KS/3600
        KM   = (KS-KH*3600)/60
        KS   = KS-KH*3600-KM*60
        NDIG = MAX(LOG10(KH+0.5)+1.,2.)
        WRITE(CFORM,
     &      '("(I",I1,".",I1,",A1,I2.2,A1,I2.2)")') NDIG,NDIG
        WRITE(CFHOUR,CFORM) KH,':',KM,':',KS
      ENDIF
!
      nfillnum = nfill(ens_nam)
      if (nfillnum > 0) then
        cfhour = cfhour(1:nfill(cfhour)) // ens_nam(1:nfill(ens_nam))
        if (me == 0) print *,' in wrtout_physics cfhour=',cfhour,
     &                                 ' ens_nam=',ens_nam
      else
        cfhour = cfhour(1:nfill(cfhour))
        if (me == 0) print *,' in wrtout_physics cfhour=',cfhour
      endif
!
      nosfc = 62
      noflx = 63
      nonst = 65
      noaer = 66
!!
      t3  = rtc()
      call MPI_BARRIER(mpi_comm_all,ierr)
      t4  = rtc()
      tbd = t4-t3
      t3  = rtc()
      SECPHY = (fhour_iau-zhour_iau)*3600.
      SECSWR = MAX(SECPHY,FHSWR)
      SECLWR = MAX(SECPHY,FHLWR)

!     SECSWR = MAX(SECPHY,FHSWR*3600.)
!     SECLWR = MAX(SECPHY,FHLWR*3600.)
!
!*** BUILD STATE ON EACH NODE ********
! build state on each node.   COMP tasks only
! assemble spectral state first then sfc state,
! then (only if liope)  flux state.
! finally (only if gocart is turned on) aer_diag state
! 
!      print *,'---- start sfc collection section -----'
      t3 = rtc()
      if(mc_comp /= MPI_COMM_NULL) then
        CALL sfc_collect(sfc_fld,global_lats_r,lonsperlar)

        if ( nstf_name(1) > 0 ) then
          call nst_collect(nst_fld,global_lats_r,lonsperlar)
        endif
!
! collect flux grids as was done with sfc grids above.
! but only if liope is true.  If liope is false,
! the fluxes are handled by the original wrtsfc
! predating the I/O task updates.
!
        call   wrtflx_a
     &             (IOPROC,noflx,zhour_iau,fhour_iau,jdate_iau,colat1,
     &              SECSWR,SECLWR,
     &              sfc_fld, flx_fld, fluxr, global_lats_r,lonsperlar,
     &              slmsk)

        if ( ngrids_aer > 0) then
           call   wrtaer
     &             (IOPROC,noaer,ZHOUR_iau,FHOUR_iau,JDATE_iau,
     &              sfc_fld, g2d_fld, global_lats_r, lonsperlar)
        endif

      endif                 ! comp node
      t4 = rtc()
      td = t4-t3
!
!  done with state build
!  NOW STATE IS ASSEMBLED ON EACH NODE.  GET EVERYTHING OFF THE COMPUTE
!  NODES (currently done with a send to the I/O task_
!  send state to I/O task.  All tasks
!
      if(sfcio_out) then
          if (me == 0) print *,'---- start sfc.f section -----'
          call sfc_only_move(ioproc)
          cosfc = sfcf//CFHOUR
          call sfc_wrt(ioproc,nosfc,cosfc,fhour_iau,jdate_iau
     &,                global_lats_r,lonsperlar)
          if ( nstf_name(1) > 0 ) then
            call nst_only_move(ioproc)
            const = nstf//CFHOUR
            call nst_wrt(ioproc,nonst,const,fhour_iau,jdate_iau
     &,                global_lats_r,lonsperlar)
          endif
          call flx_only_move(ioproc)
          cosfc = flxf//CFHOUR

!         print *,'wrtout_physics call wrtsfc to write out ',
!    &    'flx, noflx=',noflx,'cosfc=',trim(cosfc),'ZHOUR=',ZHOUR,
!    &    'FHOUR=',FHOUR,'IDATE=',IDATE,'ioproc=',ioproc

          if(me  == ioproc) then
            call baopenwt(noflx,cosfc,iostat)
!           print *,'after open flx file,',trim(cosfc)
          endif
          call  wrtflx_w(IOPROC,noflx,zhour_iau,fhour_iau,jdate_iau,
     &           colat1,SECSWR,SECLWR,slmsk, global_lats_r,lonsperlar)
      endif          !  sfcio _out
!
      t4 = rtc()
      te = t4-t3
!
!      print *,'---- start diag3d.f section -----'
      if (ldiag3d) then
        if (me == 0) print *,' wrtout_physics ldiag3d on so wrt3d '
        no3d = 64
        if(me == IOPROC)
     &  call BAOPENWT(NO3D,d3df//CFHOUR,iostat)
        if (hybrid .or. gen_coord_hybrid) then
          call WRT3D_hyb(IOPROC,no3d,zhour_iau,fhour_iau,jdate_iau,
     &                   colat1,global_lats_r,lonsperlar,pl_coeff,
     &                   SECSWR,SECLWR,sfc_fld%slmsk,flx_fld%psurf)
        else
          call WRT3D(IOPROC,no3d,zhour_iau,fhour_iau,jdate_iau,colat1,
     &               global_lats_r,lonsperlar,pl_coeff,
     &               SECSWR,SECLWR,sl,si,sfc_fld%slmsk,flx_fld%psurf)
        endif
      endif
!
!      if(me .eq. ioproc)  call wrtlog_physics(phour,fhour,idate)

      tb = rtc()
      tf = tb-t1
!     tf = tb-ta
      t2 = rtc()

       if (me == ioproc) write(0,*)' WRTOUT_PHYSICS TIME=',tf

!     print 1011,tf
!1011 format(' WRTOUT_PHYSICS TIME ',f10.4)
      timesum = timesum+(t2-t1)
!     print 1012,timesum,t2-t1,td,te,tf,t4-t3,tba,tbb,tbc,tbd
 1012 format(
     1 ' WRTOUT_PHYSICS TIME ALL TASKS  ',f10.4,f10.4,
     1 ' state, send, io  iobarr, (beginbarr),
     1 spectbarr,open, openbarr )  ' ,
     1  8f9.4)
!
      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE wrtout_restart_physics(
     &        sfc_fld, nst_fld, fhour,idate,
     &        lats_nodes_r,global_lats_r,lonsperlar,
     &        phy_f3d, phy_f2d, ngptc, nblck, ens_nam)
!!
! Feb 2010 Jun Wang, write out restart file
! Mar 2013 Jun Wang, add idea fields to restart file

      use resol_def,               ONLY: latr, levp1, levs, lonr,
     &                                   ntot2d, ntot3d
      use  namelist_physics_def,   ONLY: nstf_name, lsidea
      use layout1,                 ONLY: me, nodes, lats_node_r
      use gfs_physics_sfc_flx_mod, ONLY: Sfc_Var_Data, Flx_Var_Data
      use gfs_physics_nst_var_mod, ONLY: Nst_Var_Data
      USE machine,                 ONLY: kind_evod, kind_phys
      use idea_composition,        ONLY: pr_idea,gg,prsilvl,amgms
      implicit none
!!
      TYPE(Sfc_Var_Data)        :: sfc_fld
      TYPE(Nst_Var_Data)        :: nst_fld

      real(kind=kind_evod) fhour
      character (len=*)  :: ens_nam
!!
      integer              idate(4), ngptc, nblck
      REAL (KIND=KIND_phys)
     &     phy_f3d(ngptc,levs,ntot3d,nblck,LATS_NODE_R)
     &,    phy_f2d(LONR,LATS_NODE_R,ntot2d)
!!
      real(kind=kind_evod) sl(levs), si(levp1)
!!
      INTEGER lats_nodes_r(nodes), GLOBAL_LATS_R(LATR), lonsperlar(LATR)
      integer IOPROC, IPRINT
!     integer nfill
      character*20 cfile
!!
      IPRINT = 0
!
      cfile = 'SFCR'
!      print *,' cfile=',cfile,'ens_nam=',ens_nam(1:nfill(ens_nam))
!
!      print *,' in rest write fhour=',fhour,
!     &  'idate=',idate,' before para_fixio_w'
!
      IOPROC = nodes-1
      CALL para_fixio_w(ioproc,sfc_fld,cfile,fhour,idate,
     &                  lats_nodes_r,global_lats_r,lonsperlar,
     &                  phy_f3d, phy_f2d, ngptc, nblck, ens_nam,
     &                  lsidea, pr_idea, gg, prsilvl, amgms)
!
      if (nstf_name(1) > 0) then
        cfile = 'NSTR'
        CALL para_nst_w(ioproc,nst_fld,cfile,fhour,idate,
     &                  lats_nodes_r,global_lats_r,lonsperlar,
     &                  ens_nam)
      endif 
!
      return
      end
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE wrtlog_physics(phour,fhour,idate)
      use namelist_physics_def, ONLY: ens_nam
      implicit none

      integer   idate(4),ndigyr,nolog
      integer   ks,kh,km,ndig,nfill
      character CFHOUR*40,CFORM*40
      logical   lfnhr
      real      phour,fhour
!
!     CREATE CFHOUR

!sela set lfnhr to false for writing one step output etc.
      lfnhr = .true.    ! no output
!!mr  lfnhr = .false.   !    output
      lfnhr = 3600*abs(fhour-nint(fhour)).le.1.or.phour.eq.0
      IF(LFNHR) THEN
        KH   = NINT(FHOUR)
        NDIG = MAX(LOG10(KH+0.5)+1.,2.)
        WRITE(CFORM,'("(I",I1,".",I1,")")') NDIG,NDIG
        WRITE(CFHOUR,CFORM) KH
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
      IF(nfill(ens_nam) == 0) THEN
        CFHOUR = CFHOUR(1:nfill(CFHOUR))
      ELSE
        CFHOUR = CFHOUR(1:nfill(CFHOUR)) // ens_nam(1:nfill(ens_nam))
      END IF

      nolog = 99
      OPEN(NOlog,FILE='LOG.F'//CFHOUR,FORM='FORMATTED')
      write(nolog,100)fhour,idate
100   format(' completed mrf fhour=',f10.3,2x,4(i4,2x))
      CLOSE(NOlog)

      RETURN
      END



      SUBROUTINE sfc_collect (sfc_fld,global_lats_r,lonsperlar)
!!
      use resol_def,               ONLY: latr, lonr, ngrids_sfcc, 
     &                                   ngrids_sfcc2d,ngrids_sfcc3d,
     &                                   ngrids_flx, lsoil
      use mod_state,               ONLY:
     &                                   buff_mult_piecea2d,ngrid2d,
     &                                   buff_mult_piecea3d,ngrid3d
      use layout1,                 ONLY: lats_node_r,lats_node_r_max
      use gfs_physics_sfc_flx_mod, ONLY: Sfc_Var_Data
      USE machine,                 ONLY: kind_io8, kind_io4
      implicit none
!!
      TYPE(Sfc_Var_Data)        :: sfc_fld
!
      INTEGER, dimension(latr)  :: GLOBAL_LATS_R, lonsperlar
!!
!!!   real(kind=kind_io4) buff4(lonr,latr,4),bufsave(lonr,lats_node_r)
      real(kind=kind_io8) buffo(lonr,lats_node_r)
      real(kind=kind_io8) buffi(lonr,lats_node_r_max)
      integer             kmsk(lonr,lats_node_r_max)
!     integer kmsk(lonr,lats_node_r_max),kmskcv(lonr,lats_node_r_max)

      integer k

!     integer k, il, ubound, icount, ierr
!!
!     CHARACTER*8 labfix(4)
!     real(kind=kind_io4) yhour
!     integer,save:: version
!     data version/200004/
!     data  icount/0/
!     integer maxlats_comp
!
      ngrid2d = 1
      ngrid3d = 1
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!!
      if(.not. allocated(buff_mult_piecea2d)) then
         allocate
     1 (buff_mult_piecea2d(lonr,lats_node_r_max,1:ngrids_sfcc2d+1),
     1  buff_mult_piecea3d(lonr,lats_node_r_max,1:ngrids_sfcc3d+1))
      endif
!
      kmsk (:,1:lats_node_r) = nint(sfc_fld%slmsk(:,1:lats_node_r))
!
      ngrid2d = 1
      CALL uninterprez(1,kmsk,buffo,sfc_fld%tsea,
     &       global_lats_r,lonsperlar,buff_mult_piecea2d(1,1,ngrid2d))
!
! ngrid=2 here
!
      ngrid3d = 0
      DO k=1,LSOIL
        buffi(:,:) = sfc_fld%SMC(k,:,:)
        ngrid3d = ngrid3d+1
        CALL uninterprez(1,kmsk,buffo,buffi,global_lats_r,lonsperlar,
     &        buff_mult_piecea3d(1,1,ngrid3d))
      ENDDO
!
      ngrid2d = ngrid2d+1
      CALL uninterprez(1,kmsk,buffo,sfc_fld%weasd,
     &   global_lats_r,lonsperlar,buff_mult_piecea2d(1,1,ngrid2d))
!
      DO k=1,LSOIL
        buffi(:,:) = sfc_fld%STC(k,:,:)
!
        ngrid3d = ngrid3d+1
        CALL uninterprez(1,kmsk,buffo,buffi,global_lats_r,lonsperlar,
     &         buff_mult_piecea3d(1,1,ngrid3d))
      ENDDO
!
      ngrid2d = ngrid2d+1
      CALL uninterprez(1,kmsk,buffo,sfc_fld%TG3,
     &       global_lats_r,lonsperlar,buff_mult_piecea2d(1,1,ngrid2d))
!
      ngrid2d=ngrid2d+1
      CALL uninterprez(1,kmsk,buffo,sfc_fld%ZORL,
     &       global_lats_r,lonsperlar,buff_mult_piecea2d(1,1,ngrid2d))
!!
!     where(CV.gt.0.)
!         kmskcv=1
!     elsewhere
!         kmskcv=0
!     endwhere
!
!*********************************************************************
!   Not in version 200501
!     CALL uninterprez(1,kmskcv,buffo,CV,global_lats_r,lonsperlar)
!     CALL uninterprez(1,kmskcv,buffo,CVB,global_lats_r,lonsperlar)
!     CALL uninterprez(1,kmskcv,buffo,CVT,global_lats_r,lonsperlar)
!*********************************************************************
!jws
      ngrid2d = ngrid2d+1
      CALL uninterprez(1,kmsk,buffo,sfc_fld%ALVSF,
     &       global_lats_r,lonsperlar,buff_mult_piecea2d(1,1,ngrid2d))
!
      ngrid2d = ngrid2d+1
      CALL uninterprez(1,kmsk,buffo,sfc_fld%ALVWF,
     &       global_lats_r,lonsperlar,buff_mult_piecea2d(1,1,ngrid2d))
!
      ngrid2d = ngrid2d+1
      CALL uninterprez(1,kmsk,buffo,sfc_fld%ALNSF,
     &       global_lats_r,lonsperlar,buff_mult_piecea2d(1,1,ngrid2d))
!
      ngrid2d = ngrid2d+1
      CALL uninterprez(1,kmsk,buffo,sfc_fld%ALNWF,
     &       global_lats_r,lonsperlar,buff_mult_piecea2d(1,1,ngrid2d))
!
      ngrid2d = ngrid2d+1
      CALL uninterprez(1,kmsk,buffo,sfc_fld%SLMSK,
     &       global_lats_r,lonsperlar,buff_mult_piecea2d(1,1,ngrid2d))
!
      ngrid2d = ngrid2d+1
      CALL uninterprez(1,kmsk,buffo,sfc_fld%VFRAC,
     &       global_lats_r,lonsperlar,buff_mult_piecea2d(1,1,ngrid2d))
!
      ngrid2d = ngrid2d+1
      CALL uninterprez(1,kmsk,buffo,sfc_fld%CANOPY,
     &       global_lats_r,lonsperlar,buff_mult_piecea2d(1,1,ngrid2d))
!
      ngrid2d = ngrid2d+1
      CALL uninterprez(1,kmsk,buffo,sfc_fld%F10M,
     &       global_lats_r,lonsperlar,buff_mult_piecea2d(1,1,ngrid2d))
! T2M
      ngrid2d = ngrid2d+1
      CALL uninterprez(1,kmsk,buffo,sfc_fld%T2M,
     &       global_lats_r,lonsperlar,buff_mult_piecea2d(1,1,ngrid2d))
! Q2M
      ngrid2d = ngrid2d+1
      CALL uninterprez(1,kmsk,buffo,sfc_fld%Q2M,
     &       global_lats_r,lonsperlar,buff_mult_piecea2d(1,1,ngrid2d))
!
      ngrid2d = ngrid2d+1
      CALL uninterprez(1,kmsk,buffo,sfc_fld%VTYPE,
     &       global_lats_r,lonsperlar,buff_mult_piecea2d(1,1,ngrid2d))
!
      ngrid2d = ngrid2d+1
      CALL uninterprez(1,kmsk,buffo,sfc_fld%STYPE,
     &       global_lats_r,lonsperlar,buff_mult_piecea2d(1,1,ngrid2d))

!
      ngrid2d = ngrid2d+1
      CALL uninterprez(1,kmsk,buffo,sfc_fld%FACSF,
     &       global_lats_r,lonsperlar,buff_mult_piecea2d(1,1,ngrid2d))
!
      ngrid2d = ngrid2d+1
      CALL uninterprez(1,kmsk,buffo,sfc_fld%FACWF,
     &       global_lats_r,lonsperlar,buff_mult_piecea2d(1,1,ngrid2d))
!
      ngrid2d = ngrid2d+1
      CALL uninterprez(1,kmsk,buffo,sfc_fld%UUSTAR,
     &       global_lats_r,lonsperlar,buff_mult_piecea2d(1,1,ngrid2d))
!
      ngrid2d = ngrid2d+1
      CALL uninterprez(1,kmsk,buffo,sfc_fld%FFMM,
     &       global_lats_r,lonsperlar,buff_mult_piecea2d(1,1,ngrid2d))
!
      ngrid2d = ngrid2d+1
      CALL uninterprez(1,kmsk,buffo,sfc_fld%FFHH,
     &       global_lats_r,lonsperlar,buff_mult_piecea2d(1,1,ngrid2d))
!
!c-- XW: FOR SEA-ICE Nov04
      ngrid2d = ngrid2d+1
      CALL uninterprez(1,kmsk,buffo,sfc_fld%HICE,
     &       global_lats_r,lonsperlar,buff_mult_piecea2d(1,1,ngrid2d))
      ngrid2d = ngrid2d+1
      CALL uninterprez(1,kmsk,buffo,sfc_fld%FICE,
     &       global_lats_r,lonsperlar,buff_mult_piecea2d(1,1,ngrid2d))
      ngrid2d = ngrid2d+1
      CALL uninterprez(1,kmsk,buffo,sfc_fld%TISFC,
     &       global_lats_r,lonsperlar,buff_mult_piecea2d(1,1,ngrid2d))

!c-- XW: END SEA-ICE Nov04
!
!lu: the addition of 8 Noah-related records starts here ........................
!tprcp
      ngrid2d = ngrid2d+1
      CALL uninterprez(1,kmsk,buffo,sfc_fld%TPRCP,
     &       global_lats_r,lonsperlar,buff_mult_piecea2d(1,1,ngrid2d))
!srflag
      ngrid2d = ngrid2d+1
      CALL uninterprez(1,kmsk,buffo,sfc_fld%SRFLAG,
     &       global_lats_r,lonsperlar,buff_mult_piecea2d(1,1,ngrid2d))
!snwdph
      ngrid2d = ngrid2d+1
      CALL uninterprez(1,kmsk,buffo,sfc_fld%SNWDPH,
     &       global_lats_r,lonsperlar,buff_mult_piecea2d(1,1,ngrid2d))
!slc
!      write(0,*)'in wrt phy, before stc,ngrid2d=',ngrid2d,'ngrid3d=',
!     &   ngrid3d,'lsoil=',lsoil
      DO k=1,LSOIL
        buffi(:,:) = sfc_fld%SLC(k,:,:)
        ngrid3d = ngrid3d+1
        CALL uninterprez(1,kmsk,buffo,buffi,global_lats_r,lonsperlar,
     &         buff_mult_piecea3d(1,1,ngrid3d))
      ENDDO
!shdmin
      ngrid2d = ngrid2d+1
      CALL uninterprez(1,kmsk,buffo,sfc_fld%SHDMIN,
     &       global_lats_r,lonsperlar,buff_mult_piecea2d(1,1,ngrid2d))
!shdmax
      ngrid2d = ngrid2d+1
      CALL uninterprez(1,kmsk,buffo,sfc_fld%SHDMAX,
     &       global_lats_r,lonsperlar,buff_mult_piecea2d(1,1,ngrid2d))
!slope
      ngrid2d = ngrid2d+1
      CALL uninterprez(1,kmsk,buffo,sfc_fld%SLOPE,
     &       global_lats_r,lonsperlar,buff_mult_piecea2d(1,1,ngrid2d))
!snoalb
      ngrid2d = ngrid2d+1
      CALL uninterprez(1,kmsk,buffo,sfc_fld%SNOALB,
     &       global_lats_r,lonsperlar,buff_mult_piecea2d(1,1,ngrid2d))
!lu: the addition of 8 Noah records ends here .........................
!
! Oro
      ngrid2d = ngrid2d+1
      CALL uninterprez(1,kmsk,buffo,sfc_fld%ORO,
     &       global_lats_r,lonsperlar,buff_mult_piecea2d(1,1,ngrid2d))
!
      return
      end

      subroutine sfc_only_move(ioproc)
!
!***********************************************************************
!
      use resol_def, ONLY: ngrids_flx, ngrids_sfcc, lonr,latr
     &                    ,ngrids_sfcc2d,ngrids_sfcc3d,ngrids_sfcc
      use mod_state, ONLY: buff_mult_pieces,buff_mult_piece,
     &                     buff_mult_piecea2d,
     &                     buff_mult_piecea3d, 
     &                     ivar_global_a, ivar_global
      use layout1,   ONLY: nodes, ipt_lats_node_r, lats_node_r, 
     &                     lats_node_r_max, me, nodes_comp
      use mpi_def,   ONLY: mpi_comm_null, mpi_r_io, mc_comp, 
     &                     mpi_integer, mpi_comm_all, liope, 
     &                     info, stat
      implicit none
!
!     integer lats_nodes_r(nodes),ipt,maxfld,ioproc,nproct
      integer ioproc
      integer proc, msgtag, nproc, ierr, illen, nd1
      integer icount
      data icount/0/
      save icount
!     integer maxlats_comp
!  allocate the data structures
!
      if(icount == 0) then
         allocate(ivar_global(10))
         allocate(ivar_global_a(10,nodes))
         ivar_global(1) = ipt_lats_node_r
         ivar_global(2) =  lats_node_r
         ivar_global(3) = lats_node_r_max
         call mpi_gather(ivar_global,10,MPI_INTEGER,
     1       ivar_global_a,10,MPI_INTEGER,ioproc,mc_comp,ierr)
         if(me == ioproc) write(0,*)'in sfc_only_move, ivar_global_a=',
     &                    ivar_global_a(1:3,1:nodes)
         icount = icount+1
      endif
!!
      if(allocated(buff_mult_pieces)) then
          deallocate(buff_mult_pieces)
!     else
!         maxlats_comp=lats_node_r_max
!         if(me .eq. ioproc) then
!           maxlats_comp=ivar_global_a(3,1)
!          endif
      endif
      if(me == ioproc) then
!gwv watch this!!
         allocate (buff_mult_pieces(lonr*latr*ngrids_sfcc))
         buff_mult_pieces = 0.
      endif

      if(.not. allocated(buff_mult_piece)) then
        allocate(buff_mult_piece(lonr*lats_node_r*ngrids_sfcc))
      endif                                                   
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!   SENDLOOP of grids from comp processors to I/O task.  The
!   I/O task may or may not be a comp task also.  The
!   send logic on that task is different for these two cases
!
!  big send
!     if(me .gt. -1) return
!
       buff_mult_piece = 0.
       buff_mult_piece(1:lonr*lats_node_r*ngrids_sfcc2d)=
     1 reshape(buff_mult_piecea2d(1:lonr,1:lats_node_r,1:ngrids_sfcc2d),
     1   (/lonr*lats_node_r*ngrids_sfcc2d/)) 
       buff_mult_piece(lonr*lats_node_r*ngrids_sfcc2d+1:
     1    lonr*lats_node_r*ngrids_sfcc)=
     1 reshape(buff_mult_piecea3d(1:lonr,1:lats_node_r,1:ngrids_sfcc3d),
     1   (/lonr*lats_node_r*ngrids_sfcc3d/) )
!
      IF (ME /= ioproc) THEN      !   Sending the data
         msgtag = me
         illen  = lats_node_r
         CALL mpi_send            !  send the local grid domain
     &(buff_mult_piece,illen*lonr*ngrids_sfcc,MPI_R_IO,ioproc,
     &                  msgtag,MPI_COMM_ALL,info)
      ELSE
        if( MC_COMP /= MPI_COMM_NULL) then
!
! iotask is also a compute task.  send is replaced with direct array copy
!
          if(nodes_comp == 1) then
            buff_mult_pieces(1:lonr*lats_node_r*ngrids_sfcc) =
     1      buff_mult_piece(1:lonr*lats_node_r*ngrids_sfcc)
!                              END COMPUTE TASKS PORTION OF LOGIC
          else
!
!  END COMPUTE TASKS PORTION OF LOGIC
!  receiving part of I/O task, ioproc is the last fcst pe
!
!!      for pes ioproc
!jw        nd1=lonr*lats_node_r*ngrids_sfcc

            nd1 = 0
            DO proc=1,nodes_comp
              illen = ivar_global_a(2,proc)
              if (proc /= ioproc+1) then
                msgtag = proc-1
                CALL mpi_recv(buff_mult_pieces(nd1+1),
     1                        illen*lonr*ngrids_sfcc
     1                       ,MPI_R_IO,proc-1,
     &                        msgtag,MPI_COMM_ALL,stat,info)
              else
                buff_mult_pieces(nd1+1:nd1+lonr*illen*ngrids_sfcc) =
     1           buff_mult_piece(1:lonr*illen*ngrids_sfcc)
              endif
              nd1 = nd1 + illen*lonr*ngrids_sfcc
            enddo
          endif

        endif
      ENDIF                              !end ioproc
!!
      return
      end
!-----------------------------------
      subroutine nst_only_move(ioproc)
!...................................
!
      use resol_def, ONLY: ngrids_nst,lonr,latr
      use mod_state, ONLY: buff_mult_piecenst
      use mod_state, ONLY: buff_mult_pieces,buff_mult_piece,
     &                     ivar_global_a, ivar_global
      use layout1,   ONLY: nodes, ipt_lats_node_r, lats_node_r,
     &                     lats_node_r_max, me, nodes_comp
      use mpi_def,   ONLY: mpi_comm_null, mpi_r_io, mc_comp,
     &                     mpi_integer, mpi_comm_all, liope,
     &                     info, stat
!
      implicit none
!
!  ---  inputs/outputs:
      integer :: ioproc
!  ---  locals:
      integer :: proc, msgtag, ierr, illen, nd1
      integer, save :: icount
      data icount / 0 /
!
!  ---  allocate the data structures

      if (icount == 0) then
        if(.not.allocated(ivar_global)) allocate(ivar_global(10))
        if(.not.allocated(ivar_global_a)) 
     &          allocate(ivar_global_a(10,nodes))
        ivar_global(1) = ipt_lats_node_r
        ivar_global(2) = lats_node_r
        ivar_global(3) = lats_node_r_max
        call mpi_gather(ivar_global,10,MPI_INTEGER,
     &       ivar_global_a,10,MPI_INTEGER,ioproc,mc_comp,ierr)
       if(me == ioproc) write(0,*)'in nst_only_move, ivar_global_a=',
     &                  ivar_global_a(1:3,1:nodes)
        icount = icount + 1
      endif

      if(allocated(buff_mult_pieces)) then
          deallocate(buff_mult_pieces)
      endif
      if(me == ioproc) then
         allocate (buff_mult_pieces(lonr*latr*ngrids_nst))
         buff_mult_pieces = 0.
      endif

      if(.not. allocated(buff_mult_piece)) then
        allocate(buff_mult_piece(lonr*lats_node_r*ngrids_nst))
      endif
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!   SENDLOOP of grids from comp processors to I/O task.  The
!   I/O task may or may not be a comp task also.  The
!   send logic on that task is different for these two cases
!
!  big send
!     if(me .gt. -1) return
!
       buff_mult_piece = 0.
       buff_mult_piece(1:lonr*lats_node_r*ngrids_nst)=
     1 reshape(buff_mult_piecenst(1:lonr,1:lats_node_r,1:ngrids_nst),
     1   (/lonr*lats_node_r*ngrids_nst/))
!
      IF (ME /= ioproc) THEN      !   Sending the data
         msgtag = me
         illen  = lats_node_r
         CALL mpi_send            !  send the local grid domain
     &(buff_mult_piece,illen*lonr*ngrids_nst,MPI_R_IO,ioproc,
     &                  msgtag,MPI_COMM_ALL,info)
      ELSE
        if( MC_COMP /= MPI_COMM_NULL) then
!
! iotask is also a compute task.  send is replaced with direct array copy
!
!
          if(nodes_comp == 1) then
            buff_mult_pieces(1:lonr*lats_node_r*ngrids_nst) =
     1      buff_mult_piece(1:lonr*lats_node_r*ngrids_nst)
!                              END COMPUTE TASKS PORTION OF LOGIC
          else
!
!  END COMPUTE TASKS PORTION OF LOGIC
!  receiving part of I/O task, ioproc is the last fcst pe
!
!!      for pes ioproc
!jw        nd1=lonr*lats_node_r*ngrids_nst

            nd1 = 0
            DO proc=1,nodes_comp
              illen = ivar_global_a(2,proc)
              if (proc /= ioproc+1) then
                msgtag = proc-1
                CALL mpi_recv(buff_mult_pieces(nd1+1),
     &                        illen*lonr*ngrids_nst
     &                       ,MPI_R_IO,proc-1,
     &                        msgtag,MPI_COMM_ALL,stat,info)
              else
                buff_mult_pieces(nd1+1:nd1+lonr*illen*ngrids_nst) =
     &           buff_mult_piece(1:lonr*illen*ngrids_nst)
              endif
              nd1 = nd1 + illen*lonr*ngrids_nst
            enddo
          endif

        endif
      ENDIF                              !end ioproc

      return
!...................................
      end subroutine nst_only_move

!----------------------------------------------------------------------
      SUBROUTINE sfc_wrt(IOPROC,nw,cfile,xhour,idate
     &,                  global_lats_r,lonsperlar)
!!
      use sfcio_module
      use resol_def,    ONLY: lonr, latr, levs,ngrids_sfcc,
     &                        ngrids_sfcc2d,lsoil,ivssfc
      use layout1,      ONLY: me
      USE machine,      ONLY: kind_io8, kind_io4
      implicit none
!!
      integer,intent(in)             :: IOPROC,nw
      character*16,intent(in)        :: cfile
      real(kind=kind_io8),intent(in) :: xhour
      integer,intent(in)             :: idate(4)
      integer,intent(in)             :: GLOBAL_LATS_R(latr),
     &                                  lonsperlar(latr)
!
!-- local vars
      integer ngridss,ngrid,ngrid3

!     CHARACTER*8 labfix(4)
!     real(kind=kind_io4) yhour
!     integer,save:: version
!     data version/200501/
!
      type(sfcio_head) head
      type(sfcio_data) data
      integer iret
      logical first
      save head, first
      data first /.true./
      real(kind=kind_io4), target,allocatable :: buff_mult(:,:,:)
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!    Build surface fields in to buff_mult
!
      if (me == ioproc) then
        if(.not.allocated(buff_mult)) allocate(buff_mult(lonr,latr,
     &                                         ngrids_sfcc))
!
        ngrid = 1
        do ngridss=1,ngrids_sfcc
          call unsplit2z(ngridss,ngrids_sfcc,
     &                   buff_mult(1,1,ngridss),global_lats_r)
        enddo
!    Building surface field is done
!
!
        if (first) then
          head%clabsfc = CHAR(0)//CHAR(0)//CHAR(0)//CHAR(0)//
     &                   CHAR(0)//CHAR(0)//CHAR(0)//CHAR(0)
          head%latb    = latr
          head%lonb    = lonr
          head%ivs     = ivssfc
!         head%irealf  = 1
          head%lsoil   = lsoil
          call sfcio_alhead(head,iret)
          head%lpl     = lonsperlar(1:latr/2)
          if (lsoil == 4) then
            head%zsoil   = (/-0.1,-0.4,-1.0,-2.0/)
          elseif (lsoil == 2) then
            head%zsoil   = (/-0.1,-2.0/)
          endif
          first = .false.
        endif
        head%fhour   = xhour
        head%idate   = idate
!
        PRINT 99,nw,xhour,IDATE
99      FORMAT(1H ,'in fixio nw=',i7,2x,'HOUR=',f8.2,3x,'IDATE=',
     &  4(1X,I4))
!
        ngrid = 1
        ngrid3 = ngrids_sfcc2d+1
!
        data%tsea=>buff_mult(:,:,ngrid)
        ngrid = ngrid+1
        data%smc=>buff_mult(:,:,ngrid3:ngrid3+lsoil-1)
        ngrid3 = ngrid3+lsoil
        data%weasd=>buff_mult(:,:,ngrid)
        ngrid = ngrid+1
        data%stc=>buff_mult(:,:,ngrid3:ngrid3+lsoil-1)
        ngrid3 = ngrid3+lsoil
        data%tg3=>buff_mult(:,:,ngrid)
        ngrid = ngrid+1
        data%zorl=>buff_mult(:,:,ngrid)
        ngrid = ngrid+1
        data%alvsf=>buff_mult(:,:,ngrid)
        ngrid = ngrid+1
        data%alvwf=>buff_mult(:,:,ngrid)
        ngrid=ngrid+1
        data%alnsf=>buff_mult(:,:,ngrid)
        ngrid = ngrid+1
        data%alnwf=>buff_mult(:,:,ngrid)
        ngrid = ngrid+1
        data%slmsk=>buff_mult(:,:,ngrid)
        ngrid = ngrid+1
        data%vfrac=>buff_mult(:,:,ngrid)
        ngrid = ngrid+1
        data%canopy=>buff_mult(:,:,ngrid)
        ngrid = ngrid+1
        data%f10m=>buff_mult(:,:,ngrid)
        ngrid = ngrid+1
        data%t2m=>buff_mult(:,:,ngrid)
        ngrid = ngrid+1
        data%q2m=>buff_mult(:,:,ngrid)
        ngrid = ngrid+1
        data%vtype=>buff_mult(:,:,ngrid)
        ngrid = ngrid+1
        data%stype=>buff_mult(:,:,ngrid)
        ngrid = ngrid+1
        data%facsf=>buff_mult(:,:,ngrid)
        ngrid = ngrid+1
        data%facwf=>buff_mult(:,:,ngrid)
        ngrid = ngrid+1
        data%uustar=>buff_mult(:,:,ngrid)
        ngrid = ngrid+1
        data%ffmm=>buff_mult(:,:,ngrid)
        ngrid = ngrid+1
        data%ffhh=>buff_mult(:,:,ngrid)
        ngrid = ngrid+1
!c-- XW: FOR SEA-ICE Nov04
        data%hice=>buff_mult(:,:,ngrid)
        ngrid = ngrid+1
        data%fice=>buff_mult(:,:,ngrid)
        ngrid = ngrid+1
        data%tisfc=>buff_mult(:,:,ngrid)
        ngrid = ngrid+1
!c-- XW: END SEA-ICE Nov04
!
!lu: the addition of 8 Noah-related records starts here ...............
        data%tprcp=>buff_mult(:,:,ngrid)
        ngrid = ngrid+1
        data%srflag=>buff_mult(:,:,ngrid)
        ngrid = ngrid+1
        data%snwdph=>buff_mult(:,:,ngrid)
        ngrid = ngrid+1
        data%slc=>buff_mult(:,:,ngrid3:ngrid3+lsoil-1)
        ngrid3 = ngrid3+lsoil
        data%shdmin=>buff_mult(:,:,ngrid)
        ngrid = ngrid+1
        data%shdmax=>buff_mult(:,:,ngrid)
        ngrid = ngrid+1
        data%slope=>buff_mult(:,:,ngrid)
        ngrid = ngrid+1
        data%snoalb=>buff_mult(:,:,ngrid)
        ngrid = ngrid+1
!lu: the addition of 8 Noah records ends here .........................
!
        data%orog=>buff_mult(:,:,ngrid)      ! Orography
!
        call sfcio_swohdc(nw,cfile,head,data,iret)
!
!
      endif
      if(allocated(buff_mult)) deallocate(buff_mult)
!

      return
      end subroutine sfc_wrt

!------------------------------------------------------------------------------
      SUBROUTINE sfc_wrt_nemsio(ioproc,cfile,xhour,idate
     &,                  global_lats_r,lonsperlar)
!!
      use nemsio_module,only: nemsio_gfile,nemsio_init,nemsio_open,
     &    nemsio_writerec,nemsio_close
      use resol_def,    ONLY: lonr, latr, levs,ngrids_sfcc,
     &    ncld,ntrac,ntcw,ntoz,lsoil, ivssfc,thermodyn_id,sfcpress_id
      use layout1,      ONLY: me,idrt
      USE machine,      ONLY: kind_io8, kind_io4
!jw
      use gfs_physics_output, only : PHY_INT_STATE_ISCALAR,
     &    PHY_INT_STATE_RSCALAR,
     &    PHY_INT_STATE_1D_I,PHY_INT_STATE_1D_R,
     &    PHY_INT_STATE_2D_R_SFC,PHY_INT_STATE_3D_R
      implicit none
!!
      integer ioproc
      character*16 cfile
      real(kind=kind_io8) xhour
      integer idate(4),k,ngridss
!
      INTEGER GLOBAL_LATS_R(latr), lonsperlar(latr)

      integer i,j,ndim3,N2DR,idate7(7),nrec,kount
      integer nfhour,nfminute,nfsecondd,nfsecondn
      logical  :: outtest
      integer  :: nmetavari,nmetavarr,nmetavarl,nmetaaryi,nmetaaryr

      character(16),allocatable   :: recname(:),reclevtyp(:)
      integer,allocatable         :: reclev(:)
      character(16),allocatable   :: variname(:),varrname(:),
     &                               aryiname(:),aryrname(:)
      integer,allocatable         :: varival(:),aryilen(:),aryival(:,:)
      real(kind_io4),allocatable  :: varrval(:),aryrval(:,:)
      real(kind_io4),allocatable  :: buff_mult(:,:,:),tmp(:)
      type(nemsio_gfile) gfileout
!
!!
!     CHARACTER*8 labfix(4)
!     real(kind=kind_io4) yhour
!     integer,save:: version
!     data version/200501/
!
      integer iret
      logical first
      save first
      save  recname, reclevtyp, reclev
      save nrec,nmetavari,nmetavarr,nmetaaryi,nmetaaryr,
     &     variname,varrname,aryiname,
     &     varival,varrval,aryilen,aryival
      data first /.true./
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!    Build surface fields in to buff_mult
!
      if (me == ioproc) then
!
!
         allocate(buff_mult(lonr,latr,ngrids_sfcc))
         do ngridss=1,ngrids_sfcc
           call unsplit2z(ngridss,ngrids_sfcc,
     &                    buff_mult(1,1,ngridss),global_lats_r)
         enddo
!
!    Building surface field is done
!
         if (first) then
!write out nemsio sfc file:
          nrec=ngrids_sfcc
          kount=size(PHY_INT_STATE_ISCALAR,2)
          do i=1,kount
           if(trim(PHY_INT_STATE_ISCALAR(2,i)).eq.'OGFS_PHY'
     &     .or.trim(PHY_INT_STATE_ISCALAR(2,i)).eq.'OGFS_SFC')
     &        nmetavari=nmetavari+1
          enddo
          allocate(variname(nmetavari),varival(nmetavari))
          j=0
          do i=1,kount
           if(trim(PHY_INT_STATE_ISCALAR(2,i)).eq.'OGFS_PHY' .or.
     &      trim(PHY_INT_STATE_ISCALAR(2,i)).eq.'OGFS_SFC' )then
            j=j+1
            variname(j)=trim(PHY_INT_STATE_ISCALAR(1,i))
            if(trim(variname(j))=='latr') varival(j)=latr
            if(trim(variname(j))=='lonr') varival(j)=lonr
            if(trim(variname(j))=='levs') varival(j)=levs
            if(trim(variname(j))=='ntoz') varival(j)=ntoz
            if(trim(variname(j))=='ntcw') varival(j)=ntcw
            if(trim(variname(j))=='ncld') varival(j)=ncld
            if(trim(variname(j))=='ntrac') varival(j)=ntrac
            if(trim(variname(j))=='thermodyn_id')varival(j)=thermodyn_id
            if(trim(variname(j))=='sfcpress_id') varival(j)=sfcpress_id
            if(trim(variname(j))=='lsoil') varival(j)=lsoil
            if(trim(variname(j))=='idrt') varival(j)=idrt
            if(trim(variname(j))=='ivssfc') varival(j)=ivssfc
           endif
          enddo
!!for real var::
          nmetavarr=0
          do i=1,kount
           if(trim(PHY_INT_STATE_RSCALAR(2,i)).eq.'OGFS_PHY'
     &     .or.trim(PHY_INT_STATE_RSCALAR(2,i)).eq.'OGFS_SFC')
     &     nmetavarr=nmetavarr+1
          enddo
          allocate(varrname(nmetavarr),varrval(nmetavarr))
          j=0
          do i=1,kount
           if(trim(PHY_INT_STATE_RSCALAR(2,i)).eq.'OGFS_PHY'
     &     .or.trim(PHY_INT_STATE_RSCALAR(2,i)).eq.'OGFS_SFC')then
             j=j+1
             varrname(j)=trim(PHY_INT_STATE_RSCALAR(1,i))
             if(trim(varrname(j))=='fhour') varrval(j)=xhour
           endif
          enddo
!!for 1D ary::
          nmetaaryi=0
          do i=1,kount
           if(trim(PHY_INT_STATE_1D_I(2,i)).eq.'OGFS_PHY'
     &     .or.trim(PHY_INT_STATE_1D_I(2,i)).eq.'OGFS_SFC')
     &     nmetaaryi=nmetaaryi+1
          enddo
          allocate(aryiname(nmetaaryi),aryilen(nmetaaryi))
          j=0
          do i=1,kount
           if(trim(PHY_INT_STATE_1D_I(2,i)).eq.'OGFS_PHY'
     &     .or.trim(PHY_INT_STATE_1D_I(2,i)).eq.'OGFS_SFC')then
             j=j+1
             aryiname(j)=trim(PHY_INT_STATE_1D_I(1,i))
             if(aryiname(j)=='IDATE') aryilen(j)=size(idate)
           endif
          enddo
          allocate(aryival(maxval(aryilen),nmetaaryi) )
          aryival(1:aryilen(1),1)=idate(1:aryilen(1))
!
!!for record name, levtyp and lev
          allocate (recname(nrec),reclevtyp(nrec),reclev(nrec))
          N2DR=0
          do i=1,kount
           if(trim(PHY_INT_STATE_2D_R_SFC(2,i)).eq.'OGFS_SFC')then
            N2DR=N2DR+1
            recname(N2DR)=trim(PHY_INT_STATE_2D_R_SFC(1,i))
            reclevtyp(N2DR)=trim(trim(PHY_INT_STATE_2D_R_SFC(3,i)))
            reclev(N2DR)=1
           endif
          enddo
!
          do i=1,kount
           if(trim(PHY_INT_STATE_3D_R(2,i)).eq.'OGFS_SFC')then
            ndim3=0
            if(trim(PHY_INT_STATE_3D_R(4,i)).eq.'lsoil') then
             ndim3=lsoil
            endif
            if(ndim3>0) then
             do j=1,ndim3
              N2DR=N2DR+1
              recname(N2DR)=trim(PHY_INT_STATE_3D_R(1,i))
              reclevtyp(N2DR)=trim(trim(PHY_INT_STATE_3D_R(3,i)) )
              if(trim(PHY_INT_STATE_3D_R(4,i)).eq.'lsoil') then
                reclev(N2DR)=j
              endif
             enddo
            endif
!
           endif
          enddo
!end first
          first=.false.
         endif
     
        idate7=0
        idate7(1)=idate(4)
        idate7(2)=idate(2)
        idate7(3)=idate(3)
        idate7(4)=idate(1)
        idate7(7)=100           !denominator for second
!
        nfhour=int(xhour)
        nfminute=int((xhour-nfhour)*60)
        nfsecondn=int(((xhour-nfhour)*3600-nfminute*60)*100)
        nfsecondd=100
!
        call nemsio_init()
!
        call nemsio_open(gfileout,trim(cfile),'write',
     &    iret = iret,
     &    modelname='GFS',gdatatype='bin4',
     &    idate=idate7,nrec=nrec,
     &    dimx=lonr,dimy=latr,dimz=levs,ncldt=ncld,nmeta=5,
     &    nfhour=nfhour,nfminute=nfminute,nfsecondn=nfsecondn,
     &    nfsecondd=nfsecondd,
     &    extrameta=.true.,nmetavari=nmetavari,
     &    nmetavarr=nmetavarr,
     &    nmetaaryi=nmetaaryi,nmetaaryr=nmetaaryr,
     &    variname=variname,varival=varival,varrname=varrname,
     &    varrval=varrval,
     &    aryiname=aryiname,aryilen=aryilen,aryival=aryival,
     &    ntrac=ntrac,nsoil=lsoil,idrt=idrt,
     &    recname=recname,reclevtyp=reclevtyp,reclev=reclev)
!
        allocate(tmp(lonr*latr))
        do i=1,nrec
         tmp(:)=reshape(buff_mult(:,:,i),(/lonr*latr/) )
         call nemsio_writerec(gfileout,i,tmp,iret=iret)
        enddo
        deallocate(tmp)
        deallocate(buff_mult)
!
        call nemsio_close(gfileout)
!end write pe
      endif
!
      return
      end subroutine sfc_wrt_nemsio
!-----------------------------------
      subroutine nst_wrt(IOPROC,nw,cfile,xhour,idate,
     &                   global_lats_r,lonsperlar)
!...................................
!
      use nstio_module
      use namelist_physics_def, ONLY: lsea
      use resol_def,    ONLY: lonr, latr, ngrids_nst, ivsnst
      use layout1,      ONLY: me
      USE machine,      ONLY: kind_io8, kind_io4
!
      implicit none
!!
      integer,intent(in)             :: IOPROC,nw
      character*16,intent(in)        :: cfile
      real(kind=kind_io8),intent(in) :: xhour
      integer,intent(in)             :: idate(4)
      integer,intent(in)             :: GLOBAL_LATS_R(latr),
     &                                  lonsperlar(latr)
!
!-- local vars
      integer ngridss,ngrid
!
!     integer, save :: version = 200907
      type (nstio_head), save :: head
      type (nstio_data)       :: data
      integer iret
      logical first
      save first
      data first /.true./
      real(kind=kind_io4), target,allocatable :: buff_mult(:,:,:)
!
!  ---  build nsst fields in to buff_mult

      if (me == ioproc) then
        if(.not.allocated(buff_mult)) allocate(buff_mult(lonr,latr,
     &                                         ngrids_nst))
!
        do ngridss=1,ngrids_nst
          call unsplit2z(ngridss,ngrids_nst,
     &                   buff_mult(1,1,ngridss),global_lats_r)
        enddo
!    Building nsst field is done
        if (first) then
          head%clabnst = char(0)//char(0)//char(0)//char(0)//
     &                   char(0)//char(0)//char(0)//char(0)
          head%latb = latr
          head%lonb = lonr
          head%ivn  = ivsnst
          head%lsea = lsea

          call nstio_alhead(head,iret)

          head%lpl  = lonsperlar(1:latr/2)
          first = .false.
        endif
        head%fhour = xhour
        head%idate = idate

        print 99,nw,xhour,idate
 99     format(' in nst_wrt nw=',i7,2x,'hour=',f8.2,3x,'idate=',
     &         4(1x,i4))
        ngrid = 1
        data%slmsk    => buff_mult(:,:,ngrid)
        ngrid = ngrid + 1
        data%xt       => buff_mult(:,:,ngrid)
        ngrid = ngrid + 1
        data%xs       => buff_mult(:,:,ngrid)
        ngrid = ngrid + 1
        data%xu       => buff_mult(:,:,ngrid)
        ngrid = ngrid + 1
        data%xv       => buff_mult(:,:,ngrid)
        ngrid = ngrid + 1
        data%xz       => buff_mult(:,:,ngrid)
        ngrid = ngrid + 1
        data%zm       => buff_mult(:,:,ngrid)
        ngrid = ngrid + 1
        data%xtts     => buff_mult(:,:,ngrid)
        ngrid = ngrid + 1
        data%xzts     => buff_mult(:,:,ngrid)
        ngrid = ngrid + 1
        data%dt_cool  => buff_mult(:,:,ngrid)
        ngrid = ngrid + 1
        data%z_c      => buff_mult(:,:,ngrid)
        ngrid = ngrid + 1
        data%c_0      => buff_mult(:,:,ngrid)
        ngrid = ngrid + 1
        data%c_d      => buff_mult(:,:,ngrid)
        ngrid = ngrid + 1
        data%w_0      => buff_mult(:,:,ngrid)
        ngrid = ngrid + 1
        data%w_d      => buff_mult(:,:,ngrid)
        ngrid = ngrid + 1
        data%d_conv   => buff_mult(:,:,ngrid)
        ngrid = ngrid + 1
        data%ifd      => buff_mult(:,:,ngrid)
        ngrid = ngrid + 1
        data%tref     => buff_mult(:,:,ngrid)
        ngrid = ngrid + 1
        data%qrain    => buff_mult(:,:,ngrid)
        ngrid = ngrid + 1

        call nstio_swohdc(nw,cfile,head,ngrids_nst,data,iret)
      endif
      if(allocated(buff_mult)) deallocate(buff_mult)
!
      return
!...................................
      end subroutine nst_wrt

!------------------------------------------------------------------------------
      SUBROUTINE nst_wrt_nemsio(ioproc,cfile,xhour,idate
     &,                  global_lats_r,lonsperlar)
!!
      use nemsio_module,only: nemsio_gfile,nemsio_init,nemsio_open,
     &    nemsio_writerec,nemsio_close
      use resol_def,    ONLY: lonr, latr, levs, ngrids_nst, ivsnst,
     &    ncld,ntrac,ntcw,ntoz,lsoil,thermodyn_id,sfcpress_id
      use layout1,      ONLY: me,idrt
      USE machine,      ONLY: kind_io8, kind_io4
!jw
      use gfs_physics_output, only : PHY_INT_STATE_ISCALAR,
     &    PHY_INT_STATE_RSCALAR,
     &    PHY_INT_STATE_1D_I,PHY_INT_STATE_1D_R,
     &    PHY_INT_STATE_2D_R_NST
      implicit none
!!
      integer ioproc
      character*16 cfile
      real(kind=kind_io8) xhour
      integer idate(4),k,ngridnst
!
      INTEGER GLOBAL_LATS_R(latr), lonsperlar(latr)

      integer i,j,idate7(7),nrec,kount,N2DR
      integer nfhour,nfminute,nfsecondd,nfsecondn
      integer  :: nmetavari,nmetavarr,nmetavarl,nmetaaryi,nmetaaryr

      character(16),allocatable   :: recname(:),reclevtyp(:)
      integer,allocatable         :: reclev(:)
      character(16),allocatable   :: variname(:),varrname(:),
     &                               aryiname(:),aryrname(:)
      integer,allocatable         :: varival(:),aryilen(:),aryival(:,:)
      real(kind_io4),allocatable  :: varrval(:),aryrval(:,:)
      real(kind_io4),allocatable  :: buff_mult(:,:,:),tmp(:)
      type(nemsio_gfile) gfileout
!
!
      integer iret
      logical first
      save first
      save  recname, reclevtyp, reclev
      save nrec,nmetavari,nmetavarr,nmetaaryi,nmetaaryr,
     &     variname,varrname,aryiname,
     &     varival,varrval,aryilen,aryival
      data first /.true./
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!    Build nsst fields in to buff_mult
!
      if (me == ioproc) then
!
         allocate(buff_mult(lonr,latr,ngrids_nst))
         do ngridnst=1,ngrids_nst
           call unsplit2z(ngridnst,ngrids_nst,
     &                    buff_mult(1,1,ngridnst),global_lats_r)
         enddo
!
!    Building nsst field is done
!
         if (first) then
!write out nemsio nst file:
          nrec=ngrids_nst
          kount=size(PHY_INT_STATE_ISCALAR,2)
          nmetavari=0
          do i=1,kount
           if(trim(PHY_INT_STATE_ISCALAR(2,i)).eq.'OGFS_PHY'
     &     .or.trim(PHY_INT_STATE_ISCALAR(2,i)).eq.'OGFS_NST')
     &        nmetavari=nmetavari+1
          enddo
          allocate(variname(nmetavari),varival(nmetavari))
          j=0
          do i=1,kount
           if(trim(PHY_INT_STATE_ISCALAR(2,i)).eq.'OGFS_PHY' .or.
     &      trim(PHY_INT_STATE_ISCALAR(2,i)).eq.'OGFS_NST' )then
            j=j+1
            variname(j)=trim(PHY_INT_STATE_ISCALAR(1,i))
            if(trim(variname(j))=='latr') varival(j)=latr
            if(trim(variname(j))=='lonr') varival(j)=lonr
            if(trim(variname(j))=='levs') varival(j)=levs
            if(trim(variname(j))=='ntoz') varival(j)=ntoz
            if(trim(variname(j))=='ntcw') varival(j)=ntcw
            if(trim(variname(j))=='ncld') varival(j)=ncld
            if(trim(variname(j))=='ntrac') varival(j)=ntrac
            if(trim(variname(j))=='thermodyn_id')varival(j)=thermodyn_id
            if(trim(variname(j))=='sfcpress_id') varival(j)=sfcpress_id
            if(trim(variname(j))=='lsoil') varival(j)=lsoil
            if(trim(variname(j))=='idrt') varival(j)=idrt
            if(trim(variname(j))=='ivsnst') varival(j)=ivsnst
           endif
          enddo
!!for real var::
          nmetavarr=0
          do i=1,kount
           if(trim(PHY_INT_STATE_RSCALAR(2,i)).eq.'OGFS_PHY'
     &     .or.trim(PHY_INT_STATE_RSCALAR(2,i)).eq.'OGFS_NST')
     &     nmetavarr=nmetavarr+1
          enddo
          allocate(varrname(nmetavarr),varrval(nmetavarr))
          j=0
          do i=1,kount
           if(trim(PHY_INT_STATE_RSCALAR(2,i)).eq.'OGFS_PHY'
     &     .or.trim(PHY_INT_STATE_RSCALAR(2,i)).eq.'OGFS_NST')then
             j=j+1
             varrname(j)=trim(PHY_INT_STATE_RSCALAR(1,i))
             if(trim(varrname(j))=='fhour') varrval(j)=xhour
           endif
          enddo
!!for 1D ary::
          nmetaaryi=0
          do i=1,kount
           if(trim(PHY_INT_STATE_1D_I(2,i)).eq.'OGFS_PHY'
     &     .or.trim(PHY_INT_STATE_1D_I(2,i)).eq.'OGFS_NST')
     &     nmetaaryi=nmetaaryi+1
          enddo
          allocate(aryiname(nmetaaryi),aryilen(nmetaaryi))
          j=0
          do i=1,kount
           if(trim(PHY_INT_STATE_1D_I(2,i)).eq.'OGFS_PHY'
     &     .or.trim(PHY_INT_STATE_1D_I(2,i)).eq.'OGFS_NST')then
             j=j+1
             aryiname(j)=trim(PHY_INT_STATE_1D_I(1,i))
             if(aryiname(j)=='IDAT') aryilen(j)=size(idate)
           endif
          enddo
          allocate(aryival(maxval(aryilen),nmetaaryi) )
          aryival(1:aryilen(1),1)=idate(1:aryilen(1))
!
!!for record name, levtyp and lev
          allocate (recname(nrec),reclevtyp(nrec),reclev(nrec))
          N2DR=0
          do i=1,kount
           if(trim(PHY_INT_STATE_2D_R_NST(2,i)).eq.'OGFS_NST')then
            N2DR=N2DR+1
            recname(N2DR)=trim(PHY_INT_STATE_2D_R_NST(1,i))
           endif
          enddo
!end first
          first=.false.
         endif
     
        idate7=0
        idate7(1)=idate(4)
        idate7(2)=idate(2)
        idate7(3)=idate(3)
        idate7(4)=idate(1)
        idate7(7)=100           !denominator for second
!
        nfhour=int(xhour)
        nfminute=int((xhour-nfhour)*60)
        nfsecondn=int(((xhour-nfhour)*3600-nfminute*60)*100)
        nfsecondd=100
!
        call nemsio_init()
!
        call nemsio_open(gfileout,trim(cfile),'write',
     &    iret = iret,
     &    modelname='GFS',gdatatype='bin4',
     &    idate=idate7,nrec=nrec,
     &    dimx=lonr,dimy=latr,dimz=1,ncldt=ncld,nmeta=5,
     &    nfhour=nfhour,nfminute=nfminute,nfsecondn=nfsecondn,
     &    nfsecondd=nfsecondd,
     &    extrameta=.true.,nmetavari=nmetavari,
     &    nmetavarr=nmetavarr,
     &    nmetaaryi=nmetaaryi,nmetaaryr=nmetaaryr,
     &    variname=variname,varival=varival,varrname=varrname,
     &    varrval=varrval,
     &    aryiname=aryiname,aryilen=aryilen,aryival=aryival,
     &    ntrac=ntrac,nsoil=lsoil,idrt=idrt,
     &    recname=recname,reclevtyp=reclevtyp,reclev=reclev)
!
        allocate(tmp(lonr*latr))
        do i=1,nrec
         tmp(:)=reshape(buff_mult(:,:,i),(/lonr*latr/) )
         call nemsio_writerec(gfileout,i,tmp,iret=iret)
        enddo
        deallocate(tmp)
        deallocate(buff_mult)
!
        call nemsio_close(gfileout)
!end write pe
      endif
!
      return
      end subroutine nst_wrt_nemsio
!-----------------------------------
!------------------------------------------------------------------------------
      SUBROUTINE wrtflx_a(IOPROC,noflx,ZHOUR,FHOUR,IDATE,colat1,
     &                    SECSWR,SECLWR, sfc_fld, flx_fld, fluxr,
     &                    global_lats_r,lonsperlar,slmskful)
!!
!--  revision history
!  May 2013 S. Moorthi: fix weasd, iceth,sndpth,gflxu in flx file
!  Aug 2013 S. Moorthi: merge with gfs trunk version and add sr
!  Nov 2013 Sarah Lu  : remove commented-out lines related to aod
!                       and the lggfs3d option; correct the record
!                       number for 95 to 114

!
!
      use resol_def,               ONLY: lonr, latr, levp1, lsoil, nfxr,
     *                                   ngrids_sfcc, ngrids_flx
      use namelist_physics_def,    ONLY: climate
      use mod_state,               ONLY: buff_mult_piecef
      use layout1,                 ONLY: me, lats_node_r,
     *                                   lats_node_r_max
      use gfs_physics_sfc_flx_mod, ONLY: Sfc_Var_Data, Flx_Var_Data
      USE machine,                 ONLY: kind_io8, kind_io4,grib_undef
      implicit none
!!
!  ---  inputs/outputs:
      integer                   :: ioproc, noflx, idate(4)
      integer                   :: global_lats_r(latr), lonsperlar(latr)

      real (kind=kind_io8)      :: zhour, fhour, colat1, secswr, seclwr
      real (kind=kind_io8)      :: fluxr(nfxr,lonr,lats_node_r)

      TYPE(Sfc_Var_Data)        :: sfc_fld
      TYPE(Flx_Var_Data)        :: flx_fld

      real(kind=kind_io8) slmskful(lonr,lats_node_r)
!!
!*    PARAMETER(NFLD=18)
!     PARAMETER(NFLD=18+6)      ! 550nm AOD added
!     PARAMETER(NFLD=25+6)      ! 550nm AOD added
!  ---  locals:
!     integer, parameter :: nfld = 35
      integer, parameter :: nfld = 29

      integer i,j,k,itop,ibot,k4,l,ngrid2d,len
!     integer i,j,k,itop,ibot,k4,l,nundef,ngrid2d,len
      integer ilpds,iyr,imo,ida,ihr,ifhr,ithr,lg,ierr

      integer, dimension(lonr,lats_node_r_max) :: kmsk, kmsk0, kmskcv,  &
     &                                            kmskgrib

      real (kind=kind_io8) :: rtime, rtimsw, rtimlw, rtimer(nfld)
      real (kind=kind_io8), dimension(lonr,lats_node_r_max) :: slmskloc,&
     &                                                         glolal,  &
     &                                                         buffo

      real (kind=kind_io8) :: rflux(lonr,lats_node_r_max,nfxr)
!
      integer, allocatable :: mskf(:,:)
!
      INTEGER              IDS(255),IENS(5)
      real (kind=kind_io8) SI(LEVP1)
!
!     real (kind=kind_io4)   buff1(lonr,latr)
!     real (kind=kind_io4)   buff1l(lonr*latr)
!jws
!     real(kind=kind_io4) buff_max             ! commented by Moorthi
!jwe
!
!!
      kmsk     = nint(sfc_fld%slmsk)
      kmsk0    = 0
!
      kmskgrib = 0
      ngrid2d  = 1

!     write(0,*)' in wrtflx_a ngrids_flx=',ngrids_flx
!
      CALL uninterprez(1,kmsk,glolal,sfc_fld%slmsk,
     &       global_lats_r,lonsperlar,buff_mult_piecef(1,1,ngrid2d))
      slmskloc = glolal
      slmskful = buff_mult_piecef(1:lonr,1:lats_node_r,ngrid2d)

      if (.not. allocated(mskf)) allocate(mskf(lonr,lats_node_r))
      mskf = nint(slmskful)
!
!     write(0,*)' in wrtflx_a nfxr=',nfxr
      do k=1,nfxr
        do j=1,lats_node_r
          do i=1,lonr
            rflux(i,j,k) = fluxr(k,i,j)
          enddo
        enddo
!       write(0,*)' k=',k,' rflux=',(rflux(90,j,k),j=1,lats_node_r,3)
      enddo
!
      if (fhour > zhour) then
        rtime = 1. / (3600.*(fhour-zhour))
      else
        rtime = 0.
      endif

      if (secswr > 0.) then
        rtimsw = 1. / secswr
      else
        rtimsw = 1.
      endif

      if (seclwr > 0.) then
        rtimlw = 1. / seclwr
      else
        rtimlw = 1.
      endif

!     write(0,*)' secswr=',secswr,' seclwr=',seclwr,' rtimsw=',rtimsw,  &
!    &           ' rtimlw=',rtimlw,' rtime=',rtime

      rtimer     = rtimsw
      rtimer( 1) = rtimlw
      rtimer(20) = rtimlw       ! csulf_toa
      rtimer(22) = rtimlw       ! csdlf_sfc
      rtimer(25) = rtimlw       ! csulf_sfc

!     write(0,*)' rtimer=',rtimer

!..........................................................
!  - zonal component of momentum flux:
!
!$omp parallel do private(j,i)
      do j = 1, lats_node_r
        do i = 1, lonr
          glolal(i,j) = flx_fld%dusfc(i,j)*rtime
        enddo
      enddo
!
      ngrid2d = 1
      CALL uninterprez(2,kmsk0,buffo,glolal,global_lats_r,lonsperlar,
     &                 buff_mult_piecef(1,1,ngrid2d))

!     if (ierr /=0 ) print*,'wrtsfc gribit_gsm ierr=',ierr,'  01) ',
!    & 'Zonal compt of momentum flux (n/m**2) land and sea surface'
!..........................................................
!  - meridional component of momentum flux:
!
!$omp parallel do private(j,i)
      do j = 1, lats_node_r
        do i = 1, lonr
          glolal(i,j) = flx_fld%dvsfc(i,j)*rtime
        enddo
      enddo

      ngrid2d = ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,glolal,global_lats_r,lonsperlar,
     &                 buff_mult_piecef(1,1,ngrid2d))

!     if (ierr /=0 ) print*,'wrtsfc gribit_gsm ierr=',ierr,'  02) ',
!    & 'Merid compt of momentum flux (n/m**2) land and sea surface'
!..........................................................
!  - sensible heat flux:
!
!$omp parallel do private(j,i)
      do j = 1, lats_node_r
        do i = 1, lonr
          glolal(i,j) = flx_fld%dtsfc(i,j)*rtime
        enddo
      enddo

      ngrid2d = ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,glolal,global_lats_r,lonsperlar,
     &                 buff_mult_piecef(1,1,ngrid2d))

!     if (ierr /=0 ) print*,'wrtsfc gribit_gsm ierr=',ierr,'   03) ',
!    & 'Sensible heat flux (w/m**2) land and sea surface'
!..........................................................
!  - latent heat flux:
!
!$omp parallel do private(j,i)
      do j = 1, lats_node_r
        do i = 1, lonr
          glolal(i,j) = flx_fld%dqsfc(i,j)*rtime
        enddo
      enddo

      ngrid2d = ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,glolal,global_lats_r,lonsperlar,
     &                 buff_mult_piecef(1,1,ngrid2d))

!     if (ierr /=0 ) print*,'wrtsfc gribit_gsm ierr=',ierr,'  04) ',
!    & 'Latent heat flux (w/m**2) land and sea surface'
!..........................................................
!  - surface temperature:
!
!     do j=1,lats_node_r
!     do i=1,lonsperlar(j)
!       if (sfc_fld%tsea(i,j) < 10.0) then
!         write(0,*)' IN wrtout tsfc=',sfc_fld%tsea(i,j),' i=',i,' j=',j
!       endif
!     enddo
!     enddo
!     write(0,*)' IN wrtout tseaeq=',(sfc_fld%tsea(90,j),j=44,52)

      ngrid2d = ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,sfc_fld%tsea,
     &           global_lats_r,lonsperlar,buff_mult_piecef(1,1,ngrid2d))

!     if (ierr /=0 ) print*,'wrtsfc gribit_gsm ierr=',ierr,'  05) ',
!    & 'Temperature (k) land and sea surface'
!..........................................................
!  - volumetric soil moist content at layer 10cm and 0cm:
!
!$omp parallel do private(j,i)
      do j = 1, lats_node_r
        do i = 1, lonr
          glolal(i,j) = sfc_fld%smc(1,i,j)
        enddo
      enddo

      ngrid2d = ngrid2d+1
      CALL uninterprez(2,kmsk,buffo,glolal,global_lats_r,lonsperlar,
     &                 buff_mult_piecef(1,1,ngrid2d))

      do j=1,lats_node_r
        do i=1,lonr
          if (mskf(i,j) /= 1) buff_mult_piecef(i,j,ngrid2d) = grib_undef
        enddo
      enddo

!     if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  06) ',
!    & 'Volumetric soil moist content (frac) layer 10cm and 0cm'
!..........................................................
!  - volumetric soil moist content at layer 40cm and 10cm:
!
!$omp parallel do private(j,i)
      do j = 1, lats_node_r
        do i = 1, lonr
          glolal(i,j) = sfc_fld%smc(2,i,j)
        enddo
      enddo

      ngrid2d = ngrid2d+1
      CALL uninterprez(2,kmsk,buffo,glolal,global_lats_r,lonsperlar,
     &                 buff_mult_piecef(1,1,ngrid2d))

      do j=1,lats_node_r
        do i=1,lonr
          if (mskf(i,j) /= 1) buff_mult_piecef(i,j,ngrid2d) = grib_undef
        enddo
	enddo

!     if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  07) ',
!    & 'Volumetric soil moist content (frac) layer 40cm and 10cm'
!..........................................................
!  - temperature at layer 10cm and 0cm:
!
!$omp parallel do private(j,i)
      do j = 1, lats_node_r
        do i = 1, lonr
          glolal(i,j) = sfc_fld%stc(1,i,j)
        enddo
      enddo

      ngrid2d = ngrid2d+1
      CALL uninterprez(2,kmsk,buffo,glolal,global_lats_r,lonsperlar,
     &                 buff_mult_piecef(1,1,ngrid2d))

      do j=1,lats_node_r
        do i=1,lonr
          if (mskf(i,j) /= 1) buff_mult_piecef(i,j,ngrid2d) = grib_undef
        enddo
      enddo

!     nundef   = 0
!     buff_max = 0.
!     do j=1,lats_node_r
!       do i=1,lonr
!         if(buff_mult_piecef(i,j,ngrid2d) /= grib_undef) then
!           if(buff_mult_piecef(i,j,ngrid2d) > buff_max)
!    &                      buff_max = buff_mult_piecef(i,j,ngrid2d)
!           nundef = nundef+1
!         endif
!       enddo
!     enddo
!      write(0,*)'in wrtsfc_a, max stc=',buff_max,' grib_undef=',
!     &   grib_undef,'nundef=',nundef

!     if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  08) ',
!    & 'Temp (k) layer betw two depth below land sfc 10cm and 0cm'

!..........................................................
!  - temperature at layer 40cm and 10cm:
!
!$omp parallel do private(j,i)
      do j = 1, lats_node_r
        do i = 1, lonr
          glolal(i,j) = sfc_fld%stc(2,i,j)
        enddo
      enddo

      ngrid2d = ngrid2d+1
      CALL uninterprez(2,kmsk,buffo,glolal,global_lats_r,lonsperlar,
     &                 buff_mult_piecef(1,1,ngrid2d))
!     where(slmskful /= 1._kind_io8)
!    &     buff_mult_piecef(:,:,ngrid2d) = grib_undef

      do j=1,lats_node_r
        do i=1,lonr
          if (mskf(i,j) /= 1) buff_mult_piecef(i,j,ngrid2d) = grib_undef
        enddo
      enddo

!     if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  09) ',
!    & 'Temp (k) layer betw two depth below land sfc 40cm and 10cm'
!..........................................................
!  - water equivalent of accummulated snow depth:

      ngrid2d = ngrid2d+1
      CALL uninterprez(2,kmsk,buffo,sfc_fld%weasd,
     &           global_lats_r,lonsperlar,buff_mult_piecef(1,1,ngrid2d))

      do j=1,lats_node_r
        do i=1,lonr
          if (mskf(i,j) == 0) buff_mult_piecef(i,j,ngrid2d) = grib_undef
        enddo
      enddo

!     if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  10) ',
!    & 'Water equiv of accum snow depth (kg/m**2) at surface'
!..........................................................
!  - total sky radiation fluxes at toa and surface:
!
      do k = 1, 4
!$omp parallel do private(j,i)
        do j = 1, lats_node_r
          do i = 1, lonr
            glolal(i,j) = rflux(i,j,k)*rtimer(k)
          enddo
        enddo

        ngrid2d = ngrid2d+1
        CALL uninterprez(2,kmsk0,buffo,glolal,global_lats_r,lonsperlar,
     &                   buff_mult_piecef(1,1,ngrid2d))

!       if (ierr/=0 .and. k==1) print*,'wrtsfc gribit_gsm ierr=',ierr,'  ',
!    &   '11) Upward long wave radiation flux (w/m**2) at TOA'
!       if (ierr/=0 .and. k==2) print*,'wrtsfc gribit_gsm ierr=',ierr,'  ',
!    &   '12) Upward solar radiation flux (w/m**2) at TOA'
!       if (ierr/=0 .and. k==3) print*,'wrtsfc gribit_gsm ierr=',ierr,'  ',
!    &   '13) Upward solar radiation flux (w/m**2) at surface'
!       if (ierr/=0 .and. k==4) print*,'wrtsfc gribit_gsm ierr=',ierr,'  ',
!    &   '14) Downward solar radiation flux (w/m**2) at surface'
      enddo
!
!..........................................................
!  - for high, mid, low cloud (cover, pressure, temperature)
!
      lab_do_cloud : do k = 5, 7    ! (high, mid, low clouds)
        k4 = 4 + (k-5)*4

!  - cloud cover (h,m,l):
!
!$omp parallel do private(j,i)
        do j = 1, lats_node_r
          do i = 1, lonr
            glolal(i,j) = rflux(i,j,k)*100.*rtimsw
            if (glolal(i,j) >= 0.5) then
              kmskcv(i,j) = 1
            else
              kmskcv(i,j) = 0
            endif
          enddo
        enddo

        ngrid2d = ngrid2d+1
        CALL uninterprez(2,kmskcv,buffo,glolal,global_lats_r,lonsperlar,
     &                   buff_mult_piecef(1,1,ngrid2d))
!       where(buff_mult_piecef(:,:,ngrid2d) <= 0.5_kind_io4)
!     &       buff_mult_piecef(:,:,ngrid2d) = grib_undef
        kmskgrib = 0
        where(buff_mult_piecef(:,:,ngrid2d)<=0.5_kind_io4) kmskgrib = 1

!       l = k4 + 1
!       if (ierr/=0 .and. k==5) print*,'wrtsfc gribit_gsm ierr=',ierr,'  ',
!    &   '15) Total cloud cover (percent) high cloud layer'
!       if (ierr/=0 .and. k==6) print*,'wrtsfc gribit_gsm ierr=',ierr,'  ',
!    &   '19) Total cloud cover (percent) middle cloud layer'
!       if (ierr/=0 .and. k==7) print*,'wrtsfc gribit_gsm ierr=',ierr,'  ',
!    &   '23) Total cloud cover (percent) low cloud layer'

!  - pressure at cloud top:
!
!$omp parallel do private(j,i)
        do j = 1, lats_node_r
          do i = 1, lonr
            if (rflux(i,j,k) > 0. ) then
              glolal(i,j) = rflux(i,j,k+3)/rflux(i,j,k)
            else
              glolal(i,j) = grib_undef
            endif
          enddo
        enddo

        ngrid2d = ngrid2d+1
        CALL uninterprez(2,kmskcv,buffo,glolal,global_lats_r,lonsperlar,
     &                   buff_mult_piecef(1,1,ngrid2d))
        where(kmskgrib==1) buff_mult_piecef(:,:,ngrid2d) = grib_undef

!       if (ierr/=0 .and. k==5) print*,'wrtsfc gribit_gsm ierr=',ierr,'  ',
!    &   '16) Pressure (pa) high cloud top level'
!       if (ierr/=0 .and. k==6) print*,'wrtsfc gribit_gsm ierr=',ierr,'  ',
!    &   '20) Pressure (pa) middle cloud top level'
!       if (ierr/=0 .and. k==7) print*,'wrtsfc gribit_gsm ierr=',ierr,'  ',
!    &   '24) Pressure (pa) low cloud top level'

!  - pressure at cloud base:
!
!$omp parallel do private(j,i)
        do j = 1, lats_node_r
          do i = 1, lonr
            if (rflux(i,j,k) > 0. ) then
              glolal(i,j) = rflux(i,j,k+6)/rflux(i,j,k)
            else
              glolal(i,j) = grib_undef
            endif
          enddo
        enddo

        ngrid2d = ngrid2d+1
        CALL uninterprez(2,kmskcv,buffo,glolal,global_lats_r,lonsperlar,
     &                   buff_mult_piecef(1,1,ngrid2d))
        where(kmskgrib==1) buff_mult_piecef(:,:,ngrid2d) = grib_undef

!       if (ierr/=0 .and. k==5) print*,'wrtsfc gribit_gsm ierr=',ierr,'  ',
!    &   '17) Pressure (pa) high cloud bottom level'
!       if (ierr/=0 .and. k==6) print*,'wrtsfc gribit_gsm ierr=',ierr,'  ',
!    &   '21) Pressure (pa) middle cloud bottom level'
!       if (ierr/=0 .and. k==7) print*,'wrtsfc gribit_ggsm ierr=',ierr,'  ',
!    &   '25) Pressure (pa) low cloud bottom level'

!  - temperature at cloud top:
!
!$omp parallel do private(j,i)
        do j = 1, lats_node_r
          do i = 1, lonr
            if (rflux(i,j,k) > 0. ) then
              glolal(i,j) = rflux(i,j,k+9)/rflux(i,j,k)
            else
              glolal(i,j) = grib_undef
            endif
          enddo
        enddo

        ngrid2d = ngrid2d+1
        CALL uninterprez(2,kmskcv,buffo,glolal,global_lats_r,lonsperlar,
     &                   buff_mult_piecef(1,1,ngrid2d))
        where(kmskgrib==1) buff_mult_piecef(:,:,ngrid2d) = grib_undef

!       l = k4 + 4
!       if (ierr/=0 .and. k==5) print*,'wrtsfc gribit_gsm ierr=',ierr,'  ',
!    &   '18) Temperature (k) high cloud top level'
!       if (ierr/=0 .and. k==6) print*,'wrtsfc gribiti_gsm ierr=',ierr,'  ',
!    &   '22) Temperature (k) middle cloud top level'
!       if (ierr/=0 .and. k==7) print*,'wrtsfc gribit_gsm ierr=',ierr,'  ',
!    &   '26) Temperature (k) low cloud top level'

      enddo  lab_do_cloud
!
!..........................................................
!  - total cloud amount:
!
!$omp parallel do private(j,i)
      do j = 1, lats_node_r
        do i = 1, lonr
          glolal(i,j) = rflux(i,j,17)*100.*rtimsw
        enddo
      enddo

      ngrid2d = ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,glolal,global_lats_r,lonsperlar,   &
     &                 buff_mult_piecef(1,1,ngrid2d))

!     if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  27) ',
!    & 'Total cloud cover (percent) total atmospheric column'
!..........................................................
!  - boundary layer cloud amount:
!
!$omp parallel do private(j,i)
      do j = 1, lats_node_r
        do i = 1, lonr
          glolal(i,j) = rflux(i,j,18)*100.*rtimsw
        enddo
      enddo

      ngrid2d = ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,glolal,global_lats_r,lonsperlar,
     &                 buff_mult_piecef(1,1,ngrid2d))

!     if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  28) ',
!    & 'Total cloud cover (percent) boundary layer cloud layer'
!..........................................................
!  - surface downward lw fluxes: (use the surface temp adjusted quantities
!    to replace the original on in rec 19 of fluxr)
!
!$omp parallel do private(j,i)
      do j = 1, lats_node_r
        do i = 1, lonr
          glolal(i,j) = flx_fld%dlwsfc(i,j)*rtimlw
        enddo
      enddo

      ngrid2d = ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,glolal,global_lats_r,lonsperlar,
     &                 buff_mult_piecef(1,1,ngrid2d))
!     if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  29) ',
!    & 'Downward long wave radiation flux (w/m**2) land sea surface'
!..........................................................
!  - surface upward lw fluxes: (use the one recalc'ed from surface temp
!    to replace the original one in rec 20 of fluxr)
!
!$omp parallel do private(j,i)
      do j = 1, lats_node_r
        do i = 1, lonr
          glolal(i,j) = flx_fld%ulwsfc(i,j)*rtimlw
        enddo
      enddo

      ngrid2d = ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,glolal,global_lats_r,lonsperlar,
     &                 buff_mult_piecef(1,1,ngrid2d))

!     if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  30) ',
!    & 'Upward long wave radiation flux (w/m**2) land sea surface'
!..........................................................
!
!  - uv-b flux at surface for total sky:
!
!$omp parallel do private(j,i)
      do j=1,LATS_NODE_R
        do i=1,lonr
          glolal(i,j) = rflux(i,j,21)*rtimsw
        enddo
      enddo
      ngrid2d = ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,glolal,global_lats_r,lonsperlar,
     &                 buff_mult_piecef(1,1,ngrid2d))

!     if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  31) ',
!    & 'UV-B downward solar flux (w/m**2) land sea surface'
!..........................................................
!  - uv-b flux at surface for clear sky:

!$omp parallel do private(j,i)
      do j=1,LATS_NODE_R
        do i=1,lonr
          glolal(i,j) = rflux(i,j,22)*rtimsw
        enddo
      enddo
      ngrid2d = ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,glolal,global_lats_r,lonsperlar,
     &                 buff_mult_piecef(1,1,ngrid2d))
!     if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  32) ',
!    & 'Clear sky UV-B downward solar flux (w/m**2) land sea surface'
!     End UV-B fluxes
!..........................................................
!  - incoming solar radiation at toa:

!$omp parallel do private(j,i)
      do j = 1, lats_node_r
        do i = 1, lonr
          glolal(i,j) = rflux(i,j,23)*rtimsw
        enddo
      enddo

      ngrid2d = ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,glolal,global_lats_r,lonsperlar,
     &                 buff_mult_piecef(1,1,ngrid2d))

!     if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  33) ',
!    & 'Downward solar radiation flux (w/m**2) at TOA'
!..........................................................
!  - sw downward surface flux components:
!
      do l = 24, 27
        k = l + 2
!$omp parallel do private(j,i)
        do j = 1, lats_node_r
          do i = 1, lonr
            glolal(i,j) = rflux(i,j,l)*rtimer(k)
          enddo
        enddo

        ngrid2d = ngrid2d+1
        CALL uninterprez(2,kmsk0,buffo,glolal,global_lats_r,lonsperlar,
     &                   buff_mult_piecef(1,1,ngrid2d))

!       if (ierr/=0 .and. l==24) print*,'wrtsfc gribit_gsm ierr=',ierr,
!    &   '   34) Downward sw uv+vis beam radiation flux (w/m**2) sfc '
!       if (ierr/=0 .and. l==25) print*,'wrtsfc gribit_gsm ierr=',ierr,
!    &   '   35) Downward sw uv+vis diffuse radiation flux (w/m**2) sfc'
!       if (ierr/=0 .and. l==26) print*,'wrtsfc gribit_gsm ierr=',ierr,
!    &   '   36) Downward sw nir beam radiation flux (w/m**2) sfc '
!       if (ierr/=0 .and. l==27) print*,'wrtsfc gribit_gsm ierr=',ierr,
!    &   '   37) Downward sw nir diffuse radiation flux (w/m**2) sfc '
      enddo
!...................................................................
!  -  clear sky radiative fluxes at toa and surface:
!
      do l = 28,33
        k = l - 8
!$omp parallel do private(j,i)
        do j = 1, lats_node_r
          do i = 1, lonr
            glolal(i,j) = rflux(i,j,l)*rtimer(k)
          enddo
        enddo

      ngrid2d = ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,glolal,global_lats_r,lonsperlar,
     &                 buff_mult_piecef(1,1,ngrid2d))

!       if (ierr/=0 .and. l==28) print*,'wrtsfc gribit_gsm ierr=',ierr,
!    &   '  38) CS upward long wave radiation flux (w/m**2) at TOA'
!       if (ierr/=0 .and. l==29) print*,'wrtsfc gribit_gsm ierr=',ierr,
!    &   '  39) CS upward solar radiation flux (w/m**2) at TOA'
!       if (ierr/=0 .and. l==30) print*,'wrtsfc gribit_gsm ierr=',ierr,
!    &   '  40) CS downward long radiation flux (w/m**2) at surface'
!       if (ierr/=0 .and. l==31) print*,'wrtsfc gribit_gsm ierr=',ierr,
!    &   '  41) CS upward solar radiation flux (w/m**2) at surface'
!       if (ierr/=0 .and. l==32) print*,'wrtsfc gribit_gsm ierr=',ierr,
!    &   '  42) CS downward solar radiation flux (w/m**2) at surface'
!       if (ierr/=0 .and. l==33) print*,'wrtsfc gribit_gsm ierr=',ierr,
!    &   '  43) CS upward long wave radiation flux (w/m**2) at surface'
      enddo
!...................................................................
!  - surface albedo (derived from radiative fluxes at surface):
!
!$omp parallel do private(j,i)
      do j = 1, lats_node_r
        do i = 1, lonr
          if (rflux(i,j,4) > 0.) then
            glolal(i,j) = max(0.0, 100.0                                &
     &                  * min(1.0,rflux(i,j,3)/rflux(i,j,4)))
          else
            glolal(i,j) = 0.
          endif
        enddo
!       if (j == 46) then
!         write(0,*)' rflux(i,j,3)=',(rflux(i,j,3),i=1,192,10)
!         write(0,*)' rflux(i,j,4)=',(rflux(i,j,4),i=1,192,10)
!         write(0,*)' albedo=',(glolal(i,j),i=1,192,10)
!       endif
      enddo

      ngrid2d = ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,glolal,global_lats_r,lonsperlar,
     &                 buff_mult_piecef(1,1,ngrid2d))

!     write(0,*)' bef albedo in wrtflx_a ngrid2d=',ngrid2d

!     if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  44) ',
!    & 'Albedo (percent) land and sea surface '
!...................................................................
!  - precipitation rate (geshem unit in m, final unit = mm/s = kg/m2/s)

!$omp parallel do private(j,i)
      do j = 1, lats_node_r
        do i = 1, lonr
          glolal(i,j) = flx_fld%geshem(i,j)*1.e3*rtime
        enddo
      enddo

      ngrid2d = ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,glolal,global_lats_r,lonsperlar,
     &                 buff_mult_piecef(1,1,ngrid2d))
!     if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  45) ',
!    & 'Precipitation rate (kg/m**2/s) land and sea surface'
!...................................................................
!  - convective precipitation rate:

!$omp parallel do private(j,i)
      do j = 1, lats_node_r
        do i = 1, lonr
          glolal(i,j) = flx_fld%bengsh(i,j)*1.e3*rtime
        enddo
      enddo

      ngrid2d = ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,glolal,global_lats_r,lonsperlar,
     &                 buff_mult_piecef(1,1,ngrid2d))


!     if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  46) ',
!    & 'Convective precipitation rate (kg/m**2/s) land sea surface'
!...................................................................
!  - ground heat flux:
!
      glolal = flx_fld%gflux*rtime

      ngrid2d = ngrid2d+1
      CALL uninterprez(2,kmsk,buffo,glolal,global_lats_r,lonsperlar,
     &                 buff_mult_piecef(1,1,ngrid2d))

!$omp parallel do private(j,i)
      do j=1,lats_node_r
        do i=1,lonr
          if (mskf(i,j) == 0) buff_mult_piecef(i,j,ngrid2d) = grib_undef
        enddo
      enddo


!     if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  47) ',
!    & 'Ground heat flux (w/m**2) land and sea surface'
!...................................................................
!  - land-sea mask:
!
!     buffo = mod(slmskloc,2._kind_io8)
      ngrid2d = ngrid2d + 1
!$omp parallel do private(j,i)
      do j = 1, lats_node_r
        do i = 1, lonr
          buff_mult_piecef(i,j,ngrid2d) = mod(slmskloc(i,j),2._kind_io8)
        enddo
      enddo

!     if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  48) ',
!    & 'Land-Sea mask (1=land; 0=sea) (integer) land sea surface'
!...................................................................
!  - sea-ice concentration:

      ngrid2d = ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,sfc_fld%fice,
     &       global_lats_r,lonsperlar,buff_mult_piecef(1,1,ngrid2d))

!     if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  49) ',
!    & 'Ice concentration (ice>0; no ice=0) (1/0) land sea surface'
!...................................................................
!  - 10m u wind:

      ngrid2d = ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,flx_fld%u10m,
     &       global_lats_r,lonsperlar,buff_mult_piecef(1,1,ngrid2d))

!     if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  50) ',
!    & 'U wind (m/s) height above ground'
!...................................................................
!  - 10m v wind:

      ngrid2d = ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,flx_fld%v10m,
     &       global_lats_r,lonsperlar,buff_mult_piecef(1,1,ngrid2d))

!     if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  51) ',
!    & 'V wind (m/s) height above ground'
!...................................................................
!  - 2m temperature:

      ngrid2d = ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,sfc_fld%t2m,
     &       global_lats_r,lonsperlar,buff_mult_piecef(1,1,ngrid2d))

!     if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  52) ',
!    & 'Temperature (k) height above ground'

!  - 2m specific humidity:

      ngrid2d = ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,sfc_fld%q2m,
     &       global_lats_r,lonsperlar,buff_mult_piecef(1,1,ngrid2d))

!     if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  53) ',
!    & 'Specific humidity (kg/kg) height above ground'
!...................................................................
!  - surface pressure:

!$omp parallel do private(j,i)
      do j = 1, lats_node_r
        do i = 1, lonr
          glolal(i,j) = flx_fld%psurf(i,j)
        enddo
      enddo
      ngrid2d = ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,glolal,global_lats_r,lonsperlar,
     &     buff_mult_piecef(1,1,ngrid2d))

!     if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  54) ',
!    & 'Pressure (pa) land and sea surface'
!...................................................................
!  - maximum temperature:
!
      ngrid2d = ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,flx_fld%tmpmax,
     &       global_lats_r,lonsperlar,buff_mult_piecef(1,1,ngrid2d))

!     if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  55) ',
!    & 'Maximum temperature (k) height above ground'
!...................................................................
!  - minimum temperature:
!
      ngrid2d = ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,flx_fld%tmpmin,
     &       global_lats_r,lonsperlar,buff_mult_piecef(1,1,ngrid2d))

!     if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  56) ',
!    & 'Minimum temperature (k) height above ground'
!...................................................................
!  - maximum specific humidity

      ngrid2d = ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,flx_fld%spfhmax,
     &       global_lats_r,lonsperlar,buff_mult_piecef(1,1,ngrid2d))

!     if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  57) ',
!    & 'Maximum specific humidity (kg/kg) height above ground'
!...................................................................
!  - minimum specific humidity

      ngrid2d = ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,flx_fld%spfhmin,
     &       global_lats_r,lonsperlar,buff_mult_piecef(1,1,ngrid2d))

!     if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  58) ',
!    & 'Minimum specific humidity (kg/kg) height above ground'
!...................................................................
!  - runoff (accumulative value)

!$omp parallel do private(j,i)
      do j = 1, lats_node_r
        do i = 1, lonr
          glolal(i,j) = flx_fld%runoff(i,j)*1.e3
        enddo
      enddo
      ngrid2d = ngrid2d+1
      CALL uninterprez(2,kmsk,buffo,glolal,global_lats_r,lonsperlar,
     &                 buff_mult_piecef(1,1,ngrid2d))

      do j=1,lats_node_r
        do i=1,lonr
          if (mskf(i,j) == 0) buff_mult_piecef(i,j,ngrid2d) = grib_undef
        enddo
      enddo

!     if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  59) ',
!    & 'Runoff (kg/m**2) land and sea surface'
!...................................................................
!  - potential evaporation rate

!$omp parallel do private(j,i)
      do j = 1, lats_node_r
        do i = 1, lonr
          glolal(i,j) = flx_fld%ep(i,j)*rtime
        enddo
      enddo
      ngrid2d = ngrid2d+1
      CALL uninterprez(2,kmsk,buffo,glolal,global_lats_r,lonsperlar,
     &                 buff_mult_piecef(1,1,ngrid2d))

      do j=1,lats_node_r
        do i=1,lonr
          if (mskf(i,j) == 0) buff_mult_piecef(i,j,ngrid2d) = grib_undef
        enddo
      enddo

!     if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  60) ',
!    & 'Potential evaporation rate (w/m**/) land and sea surface'
!...................................................................
!  - cloud work function
!
!$omp parallel do private(j,i)
      do j = 1, lats_node_r
        do i = 1, lonr
          glolal(i,j) = flx_fld%cldwrk(i,j)*rtime
        enddo
      enddo

      ngrid2d = ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,glolal,global_lats_r,lonsperlar,
     &     buff_mult_piecef(1,1,ngrid2d))

!     if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  61) ',
!    & 'Cloud work function (j/kg) total atmospheric column'
!...................................................................
!  - zonal gravity wave stress
!
!$omp parallel do private(j,i)
      do j = 1, lats_node_r
        do i = 1, lonr
          glolal(i,j) = flx_fld%dugwd(i,j)*rtime
        enddo
      enddo

      ngrid2d = ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,glolal,global_lats_r,lonsperlar,
     &     buff_mult_piecef(1,1,ngrid2d))
!     if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  62) ',
!    & 'Zonal gravity wave stress (n/m**2) land and sea surface'
!...................................................................
!  - meridional gravity wave stress
!
!$omp parallel do private(j,i)
      do j = 1, lats_node_r
        do i = 1, lonr
          glolal(i,j) = flx_fld%dvgwd(i,j)*rtime
        enddo
      enddo

      ngrid2d = ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,glolal,global_lats_r,lonsperlar,
     &     buff_mult_piecef(1,1,ngrid2d))

!     if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  63) ',
!    & 'Meridional gravity wave stress (n/m**2) land sea surface'
!...................................................................
!  - boundary layer height

      ngrid2d = ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,flx_fld%hpbl,
     &       global_lats_r,lonsperlar,buff_mult_piecef(1,1,ngrid2d))

!     if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  64) ',
!    & 'Boundary layer height '
!...................................................................
!  - precipitable water
!
      ngrid2d = ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,flx_fld%pwat,
     &       global_lats_r,lonsperlar,buff_mult_piecef(1,1,ngrid2d))

!     if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  65) ',
!    & 'Precipitable water (kg/m**2) total atmospheric column'
!...................................................................
!  - convective clouds
!    * labeled instantaneous but actually averaged over fhswr seconds
!
!$omp parallel do private(j,i)
      do j = 1, lats_node_r
        do i = 1, lonr
          glolal(i,j) = sfc_fld%cv(i,j)*100.0
          if (glolal(i,j) >= 0.5) then
            kmskcv(i,j) = 1
          else
            kmskcv(i,j) = 0
          endif
        enddo
      enddo

      ngrid2d = ngrid2d+1
      CALL uninterprez(2,kmskcv,buffo,glolal,global_lats_r,lonsperlar,
     &                 buff_mult_piecef(1,1,ngrid2d))
      kmskgrib = 0
      where(buff_mult_piecef(:,:,ngrid2d)<0.5_kind_io8) kmskgrib = 1

!      where(buff_mult_piecef(:,:,ngrid2d)<0.5_kind_io8)
!     &     buff_mult_piecef(:,:,ngrid2d)=grib_undef

!     if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  66) ',
!    & 'Total cloud cover (percent) convective cloud layer'
!.................................................
!  - pressure at convective cloud top
!
!$omp parallel do private(j,i)
      do j = 1, lats_node_r
        do i = 1, lonr
          glolal(i,j) = 0.
          if (sfc_fld%cv(i,j) > 0.) then
            glolal(i,j) = sfc_fld%cvt(i,j)*1.e3
          endif
        enddo
      enddo

      ngrid2d = ngrid2d+1
      CALL uninterprez(2,kmskcv,buffo,glolal,global_lats_r,lonsperlar,
     &                 buff_mult_piecef(1,1,ngrid2d))
      where(kmskgrib == 1) buff_mult_piecef(:,:,ngrid2d) = grib_undef

!     if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  67) ',
!    & 'Pressure (pa) convective cloud top level'
!.................................................
!  - pressure at convective cloud bottom
!
!$omp parallel do private(j,i)
      do j = 1, lats_node_r
        do i = 1, lonr
          glolal(i,j) = 0.
          if (sfc_fld%cv(i,j) > 0.) then
            glolal(i,j) = sfc_fld%cvb(i,j)
          endif
        enddo
      enddo

      ngrid2d = ngrid2d+1
      CALL uninterprez(2,kmskcv,buffo,glolal,global_lats_r,lonsperlar,
     &                 buff_mult_piecef(1,1,ngrid2d))
      where(kmskgrib == 1) buff_mult_piecef(:,:,ngrid2d) = grib_undef

!     if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  68) ',
!    & 'Pressure (pa) convective cloud bottom level'
!.................................................
!  - sea ice thickness

      ngrid2d = ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,sfc_fld%hice,
     &     global_lats_r,lonsperlar,buff_mult_piecef(1,1,ngrid2d))
!     if(ierr.ne.0)print*,'wrtsfc gribit_gsm ierr=',ierr,'  ',

      do j=1,lats_node_r
        do i=1,lonr
          if (mskf(i,j) == 1) buff_mult_piecef(i,j,ngrid2d) = grib_undef
        enddo
      enddo

!     write(0,*)' HICE=', buff_mult_piecef(1:10,2,ngrid2d),
!    &' slmskful=',slmskful(1:10,2),' ngrid2d=',ngrid2d
!     if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  69) ',
!    & 'Sea ice thickness (m) category 1'
!.................................................
!  - volumetric soil moist content (layer 100cm and 40cm)
!

      if (lsoil > 2) then
!$omp parallel do private(j,i)
        do j = 1, lats_node_r
          do i = 1, lonr
            glolal(i,j) = sfc_fld%smc(3,i,j)
          enddo
        enddo

        ngrid2d = ngrid2d+1
        CALL uninterprez(2,kmsk,buffo,glolal,global_lats_r,lonsperlar,
     &                   buff_mult_piecef(1,1,ngrid2d))

      do j=1,lats_node_r
        do i=1,lonr
          if (mskf(i,j) /= 1) buff_mult_piecef(i,j,ngrid2d) = grib_undef
        enddo
      enddo

!       if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  70) ',
!    &   'Volumetric soil moist content (frac) layer 100cm and 40cm'
!..........................................................
!  - volumetric soil moist content (layer 200cm and 100cm)
!
!$omp parallel do private(j,i)
        do j = 1, lats_node_r
          do i = 1, lonr
            glolal(i,j) = sfc_fld%smc(4,i,j)
          enddo
        enddo

        ngrid2d = ngrid2d+1
        CALL uninterprez(2,kmsk,buffo,glolal,global_lats_r,lonsperlar,
     &                   buff_mult_piecef(1,1,ngrid2d))

      do j=1,lats_node_r
        do i=1,lonr
          if (mskf(i,j) /= 1) buff_mult_piecef(i,j,ngrid2d) = grib_undef
        enddo
      enddo

!       if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  71) ',
!    &   'Volumetric soil moist content(frac) layer 200cm and 100cm'
!..........................................................
!  - temperature for layer 100cm and 40cm below sfc   do j=1,LATS_NODE_R

!
!$omp parallel do private(j,i)
        do j = 1, lats_node_r
          do i = 1, lonr
            glolal(i,j) = sfc_fld%stc(3,i,j)
          enddo
        enddo

        ngrid2d = ngrid2d+1
        CALL uninterprez(2,kmsk,buffo,glolal,global_lats_r,lonsperlar,
     &                   buff_mult_piecef(1,1,ngrid2d))

      do j=1,lats_node_r
        do i=1,lonr
          if (mskf(i,j) /= 1) buff_mult_piecef(i,j,ngrid2d) = grib_undef
        enddo
      enddo


!       if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  72) ',
!    &   'Temp (k) betw two depth below land sfc 100cm and 40cm'
!..........................................................
!  - temperature for layer 200cm and 100cm below sfc
!
!$omp parallel do private(j,i)
        do j = 1, lats_node_r
          do i = 1, lonr
            glolal(i,j) = sfc_fld%stc(4,i,j)
          enddo
        enddo

        ngrid2d = ngrid2d+1
        CALL uninterprez(2,kmsk,buffo,glolal,global_lats_r,lonsperlar,
     &                   buff_mult_piecef(1,1,ngrid2d))

      do j=1,lats_node_r
        do i=1,lonr
          if (mskf(i,j) /= 1) buff_mult_piecef(i,j,ngrid2d) = grib_undef
        enddo
      enddo

!       if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  73) ',
!    &   'Temp (k) betw two depth below land sfc 200cm and 100cm'
      endif
!..........................................................
!  - liquid soil moist content layer 10cm and 0cm
!
!$omp parallel do private(j,i)
      do j = 1, lats_node_r
        do i = 1, lonr
          glolal(i,j) = sfc_fld%slc(1,i,j)
        enddo
      enddo

      ngrid2d = ngrid2d+1
      CALL uninterprez(2,kmsk,buffo,glolal,global_lats_r,lonsperlar,
     &                 buff_mult_piecef(1,1,ngrid2d))

      do j=1,lats_node_r
        do i=1,lonr
          if (mskf(i,j) /= 1) buff_mult_piecef(i,j,ngrid2d) = grib_undef
        enddo
      enddo

!     if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  74) ',
!    & 'Liquid soil moist content (frac) layer 10cm and 0cm'
!..........................................................
!  - liquid soil moist content layer 40cm and 10cm
!
!$omp parallel do private(j,i)
      do j = 1, lats_node_r
        do i = 1, lonr
          glolal(i,j) = sfc_fld%slc(2,i,j)
        enddo
      enddo

      ngrid2d = ngrid2d+1
      CALL uninterprez(2,kmsk,buffo,glolal,global_lats_r,lonsperlar,
     &                 buff_mult_piecef(1,1,ngrid2d))

      do j=1,lats_node_r
        do i=1,lonr
          if (mskf(i,j) /= 1) buff_mult_piecef(i,j,ngrid2d) = grib_undef
        enddo
      enddo

!     if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  75) ',
!    & 'Liquid soil moist content (frac) layer 40cm and 10cm'
!..........................................................
!  - liquid soil moist content layer 100cm and 40cm
!
      if (lsoil > 2) then
!$omp parallel do private(j,i)
        do j = 1, lats_node_r
          do i = 1, lonr
            glolal(i,j) = sfc_fld%slc(3,i,j)
          enddo
        enddo

        ngrid2d = ngrid2d+1
        CALL uninterprez(2,kmsk,buffo,glolal,global_lats_r,lonsperlar,
     &                   buff_mult_piecef(1,1,ngrid2d))

        do j=1,lats_node_r
         do i=1,lonr
          if (mskf(i,j) /= 1) buff_mult_piecef(i,j,ngrid2d) = grib_undef
         enddo
        enddo

!       if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  76) ',
!    &   'Liquid soil moist content (frac) layer 100cm and 40cm'
!..........................................................
!  - liquid soil moist content layer 200cm and 100cm
!
!$omp parallel do private(j,i)
        do j = 1, lats_node_r
          do i = 1, lonr
            glolal(i,j) = sfc_fld%slc(4,i,j)
          enddo
        enddo

        ngrid2d = ngrid2d+1
        CALL uninterprez(2,kmsk,buffo,glolal,global_lats_r,lonsperlar,
     &                   buff_mult_piecef(1,1,ngrid2d))

        do j=1,lats_node_r
         do i=1,lonr
          if (mskf(i,j) /= 1) buff_mult_piecef(i,j,ngrid2d) = grib_undef
         enddo
       enddo

!       if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  77) ',
!    &   'Liquid soil moist content (frac) layer 200cm and 100cm'
      endif
!..........................................................
!  - snow depth
!
!$omp parallel do private(j,i)
      do j = 1, lats_node_r
        do i = 1, lonr
          glolal(i,j) = sfc_fld%snwdph(i,j) * 0.001  ! convert from mm to m
        enddo
      enddo

      ngrid2d = ngrid2d+1
      CALL uninterprez(2,kmsk,buffo,glolal,global_lats_r,lonsperlar,
     &                 buff_mult_piecef(1,1,ngrid2d))

      do j=1,lats_node_r
        do i=1,lonr
          if (mskf(i,j) == 0) buff_mult_piecef(i,j,ngrid2d) = grib_undef
        enddo
      enddo

!     if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  78) ',
!    & 'Snow depth (m) land surface'
!..........................................................
!  - canopy water content
!
!     lbm = (slmskful == 1._kind_io8)

      ngrid2d = ngrid2d+1
      CALL uninterprez(2,kmsk,buffo,sfc_fld%canopy,
     &       global_lats_r,lonsperlar,buff_mult_piecef(1,1,ngrid2d))

      do j=1,lats_node_r
        do i=1,lonr
          if (mskf(i,j) /= 1) buff_mult_piecef(i,j,ngrid2d) = grib_undef
        enddo
      enddo

!     if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  79) ',
!    & 'Canopy water content (kg/m^2) land surface'
!..........................................................
!  - the following 30 records are for land mdl use
!  - surface roughness
!
!$omp parallel do private(j,i)
      do j = 1, lats_node_r
        do i = 1, lonr
          glolal(i,j) = sfc_fld%zorl(i,j) * 0.01  ! convert from cm to m
        enddo
      enddo

      ngrid2d = ngrid2d+1
      CALL uninterprez(1,kmsk0,buffo,glolal,global_lats_r,lonsperlar,
     &     buff_mult_piecef(1,1,ngrid2d))

!     if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  80) ',
!    & 'Surface roughness (m)'
!..........................................................
!  - vegetation fraction
!
!$omp parallel do private(j,i)
      do j = 1, lats_node_r
        do i = 1, lonr
          glolal(i,j) = sfc_fld%vfrac(i,j) * 100.0
        enddo
      enddo

      ngrid2d = ngrid2d+1
      CALL uninterprez(1,kmsk,buffo,glolal,global_lats_r,lonsperlar,
     &                 buff_mult_piecef(1,1,ngrid2d))

      do j=1,lats_node_r
        do i=1,lonr
          if (mskf(i,j) /= 1) buff_mult_piecef(i,j,ngrid2d) = grib_undef
        enddo
      enddo

!     if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'   81) ',
!    & 'Vegetation fraction (fractional) land surface'
!..........................................................
!  - vegetation type
!
      ngrid2d = ngrid2d+1
      CALL uninterprez(1,kmsk,glolal,sfc_fld%vtype,
     &       global_lats_r,lonsperlar,buff_mult_piecef(1,1,ngrid2d))
!     buffo = MOD(glolal,2._kind_io8)

      do j=1,lats_node_r
        do i=1,lonr
          if (mskf(i,j) /= 1) buff_mult_piecef(i,j,ngrid2d) = grib_undef
        enddo
      enddo

!     if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  82) ',
!    & 'Vegetation type land surface'
!..........................................................
!  - soil type
!
      ngrid2d = ngrid2d+1
      CALL uninterprez(1,kmsk,glolal,sfc_fld%stype,
     &       global_lats_r,lonsperlar,buff_mult_piecef(1,1,ngrid2d))
!     buffo=MOD(glolal,2._kind_io8)

      do j=1,lats_node_r
        do i=1,lonr
          if (mskf(i,j) /= 1) buff_mult_piecef(i,j,ngrid2d) = grib_undef
        enddo
      enddo

!     if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  83) ',
!    & 'Soil type land surface'
!..........................................................
!  - slope type
!
      ngrid2d = ngrid2d+1
      CALL uninterprez(1,kmsk,glolal,sfc_fld%slope,
     &       global_lats_r,lonsperlar,buff_mult_piecef(1,1,ngrid2d))
!     buffo=MOD(glolal,2._kind_io8)

      do j=1,lats_node_r
        do i=1,lonr
          if (mskf(i,j) /= 1) buff_mult_piecef(i,j,ngrid2d) = grib_undef
        enddo
      enddo

!     if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  84) ',
!    & 'Slope type land surface'
!..........................................................
!  - frictional velocity

      ngrid2d = ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,sfc_fld%uustar,
     &       global_lats_r,lonsperlar,buff_mult_piecef(1,1,ngrid2d))

!     if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  85) ',
!    & 'Frictional velocity (m/s)'
!..........................................................
!  - surface height
!
      ngrid2d = ngrid2d+1
      CALL uninterprez(1,kmsk,buffo,sfc_fld%oro,
     &       global_lats_r,lonsperlar,buff_mult_piecef(1,1,ngrid2d))

!     if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  86) ',
!    & 'Surface height (m)'
!..........................................................
!  - freezing precip flag
!
      ngrid2d = ngrid2d+1
      CALL uninterprez(1,kmsk,buffo,sfc_fld%srflag,
     &       global_lats_r,lonsperlar,buff_mult_piecef(1,1,ngrid2d))

      do j=1,lats_node_r
        do i=1,lonr
          if (mskf(i,j) /= 1) buff_mult_piecef(i,j,ngrid2d) = grib_undef
        enddo
      enddo

!     if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  87) ',
!    & 'Freezing precip flag land surface'
!..........................................................
!  - exchange coefficient ch
!
      ngrid2d = ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,flx_fld%chh,
     &       global_lats_r,lonsperlar,buff_mult_piecef(1,1,ngrid2d))

!     if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  88) ',
!    & 'Exchange coefficient ch(m/s)'
!..........................................................
!  - exchange coefficient cm
!
      ngrid2d = ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,flx_fld%cmm,
     &       global_lats_r,lonsperlar,buff_mult_piecef(1,1,ngrid2d))

!     if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  89) ',
!    & 'Exchange coefficient cm(m/s)'
!..........................................................
!  - potential evaporation rate
!
      ngrid2d = ngrid2d+1
      CALL uninterprez(2,kmsk,buffo,flx_fld%EPI,
     &       global_lats_r,lonsperlar,buff_mult_piecef(1,1,ngrid2d))

      do j=1,lats_node_r
        do i=1,lonr
          if (mskf(i,j) /= 1) buff_mult_piecef(i,j,ngrid2d) = grib_undef
        enddo
      enddo

!     if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  90) ',
!    & 'Potential evaporation rate (w/m**2) land and sea surface'
!..........................................................
!  - downward long wave radiation flux (instantaneous value)
!
      if (.not. climate) then       ! do not output those fields in climate mode

        ngrid2d = ngrid2d+1
        CALL uninterprez(2,kmsk0,buffo,flx_fld%DLWSFCI,
     &         global_lats_r,lonsperlar,buff_mult_piecef(1,1,ngrid2d))

!       if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  91) ',
!    & 'Downward long wave radiation flux (w/m**2) '
!..........................................................
!  - upward long wave radiation flux (instantaneous value)
!
        ngrid2d = ngrid2d+1
        CALL uninterprez(2,kmsk0,buffo,flx_fld%ULWSFCI,
     &         global_lats_r,lonsperlar,buff_mult_piecef(1,1,ngrid2d))

!       if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  92) ',
!    & 'Upward long wave radiation flux (w/m**2)'
!..........................................................
!  - upward short wave radiation flux (instantaneous value)

        ngrid2d = ngrid2d+1
        CALL uninterprez(2,kmsk0,buffo,flx_fld%USWSFCI,
     &         global_lats_r,lonsperlar,buff_mult_piecef(1,1,ngrid2d))
!       if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  93) ',
!    & 'Upward short wave radiation flux (w/m**2)'
!..........................................................
!  - downward short wave radiation flux (instantaneous value)
!
        ngrid2d = ngrid2d+1
        CALL uninterprez(2,kmsk0,buffo,flx_fld%DSWSFCI,
     &         global_lats_r,lonsperlar,buff_mult_piecef(1,1,ngrid2d))

!       if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  94) ',
!    & 'Downward short wave radiation flux (w/m**2)'
!..........................................................
!  - sensible heat flux (instantaneous value)
!
        ngrid2d = ngrid2d+1
        CALL uninterprez(2,kmsk0,buffo,flx_fld%DTSFCI,
     &         global_lats_r,lonsperlar,buff_mult_piecef(1,1,ngrid2d))

!       if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  95) ',
!    & 'Sensible heat flux (w/m**2) land and sea surface'
!..........................................................
!  - latent heat flux (instantaneous value)
!
        ngrid2d = ngrid2d+1
        CALL uninterprez(2,kmsk0,buffo,flx_fld%DQSFCI,
     &         global_lats_r,lonsperlar,buff_mult_piecef(1,1,ngrid2d))

!       if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  96) ',
!    & 'Latent heat flux (w/m**2) land and sea surface'
!..........................................................
!  - ground heat flux (instantaneous value)
!
        ngrid2d = ngrid2d+1
        CALL uninterprez(2,kmsk,buffo,flx_fld%GFLUXI,
     &         global_lats_r,lonsperlar,buff_mult_piecef(1,1,ngrid2d))

        do j=1,lats_node_r
         do i=1,lonr
          if (mskf(i,j) == 0) buff_mult_piecef(i,j,ngrid2d) = grib_undef
         enddo
        enddo

!       if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  97) ',
!    & 'Ground heat flux (w/m**2) land and sea surface'
      endif   ! if_.not._climate
!..........................................................
!  - surface runoff
!
!$omp parallel do private(j,i)
      do j = 1, lats_node_r
        do i = 1, lonr
          glolal(i,j) = flx_fld%runoff(i,j) * 1000.0
        enddo
      enddo

      ngrid2d = ngrid2d+1
      CALL uninterprez(2,kmsk,buffo,glolal,global_lats_r,lonsperlar,
     &       buff_mult_piecef(1,1,ngrid2d))

      do j=1,lats_node_r
        do i=1,lonr
          if (mskf(i,j) /= 1) buff_mult_piecef(i,j,ngrid2d) = grib_undef
        enddo
      enddo

!     if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  98) ',
!    & 'Surface runoff (kg/m^2) land surface'
!..........................................................
!  - lowest model level temp
!
      ngrid2d = ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,flx_fld%t1,
     &       global_lats_r,lonsperlar,buff_mult_piecef(1,1,ngrid2d))

!     if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  99) ',
!    & 'Lowest model level temp (k)'
!..........................................................
!  - lowest model specific humidity
!
      ngrid2d = ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,flx_fld%q1,
     &       global_lats_r,lonsperlar,buff_mult_piecef(1,1,ngrid2d))

!     if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  100) ',
!    & 'Lowest model specific humidity (kg/kg)'
!..........................................................
!  - lowest model u wind
!
      ngrid2d = ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,flx_fld%u1,
     &        global_lats_r,lonsperlar,buff_mult_piecef(1,1,ngrid2d))

!     if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  101) ',
!    & 'Lowest model u wind (m/s)'
!..........................................................
!  - lowest model v wind
!
      ngrid2d = ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,flx_fld%v1,
     &        global_lats_r,lonsperlar,buff_mult_piecef(1,1,ngrid2d))

!     if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  102) ',
!    & 'Lowest model v wind (m/s)'
!..........................................................
!  - lowest model level height
!
      ngrid2d = ngrid2d+1
      CALL uninterprez(2,kmsk,buffo,flx_fld%zlvl,
     &      global_lats_r,lonsperlar,buff_mult_piecef(1,1,ngrid2d))

      do j=1,lats_node_r
        do i=1,lonr
          if (mskf(i,j) /= 1) buff_mult_piecef(i,j,ngrid2d) = grib_undef
        enddo
      enddo

!     if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  103) ',
!    & 'Lowest model level height (m) land surface'
!..........................................................
!  - direct evaporation from bare soil
!
!$omp parallel do private(j,i)
      do j = 1, lats_node_r
        do i = 1, lonr
          glolal(i,j) = flx_fld%evbsa(i,j) * rtime
        enddo
      enddo
      ngrid2d = ngrid2d+1
      CALL uninterprez(2,kmsk,buffo,glolal,global_lats_r,lonsperlar,
     &     buff_mult_piecef(1,1,ngrid2d))

      do j=1,lats_node_r
        do i=1,lonr
          if (mskf(i,j) /= 1) buff_mult_piecef(i,j,ngrid2d) = grib_undef
        enddo
      enddo

!     if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  104) ',
!    & 'Direct evaporation from bare soil(w/m^2) land surface'
!..........................................................
!  - canopy water evaporation
!
!$omp parallel do private(j,i)
      do j = 1, lats_node_r
        do i = 1, lonr
          glolal(i,j) = flx_fld%evcwa(i,j) * rtime
        enddo
      enddo

      ngrid2d = ngrid2d+1
      CALL uninterprez(2,kmsk,buffo,glolal,global_lats_r,lonsperlar,
     &     buff_mult_piecef(1,1,ngrid2d))

      do j=1,lats_node_r
        do i=1,lonr
          if (mskf(i,j) /= 1) buff_mult_piecef(i,j,ngrid2d) = grib_undef
        enddo
      enddo

!     if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  105) ',
!    & 'Canopy water evaporation(w/m^2) land surface'
!..........................................................
!  - transpiration
!
!$omp parallel do private(j,i)
      do j = 1, lats_node_r
        do i = 1, lonr
          glolal(i,j) = flx_fld%transa(i,j) * rtime
        enddo
      enddo

      ngrid2d = ngrid2d+1
      CALL uninterprez(2,kmsk,buffo,glolal,global_lats_r,lonsperlar,
     &     buff_mult_piecef(1,1,ngrid2d))

      do j=1,lats_node_r
        do i=1,lonr
          if (mskf(i,j) /= 1) buff_mult_piecef(i,j,ngrid2d) = grib_undef
        enddo
      enddo

!     if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  106) ',
!    & 'Transpiration (w/m^2) land surface'
!..........................................................
!  - snow sublimation
!
!$omp parallel do private(j,i)
      do j = 1, lats_node_r
        do i = 1, lonr
          glolal(i,j) = flx_fld%sbsnoa(i,j) * rtime
        enddo
      enddo

      ngrid2d = ngrid2d+1
      CALL uninterprez(2,kmsk,buffo,glolal,global_lats_r,lonsperlar,
     &     buff_mult_piecef(1,1,ngrid2d))

      do j=1,lats_node_r
        do i=1,lonr
          if (mskf(i,j) /= 1) buff_mult_piecef(i,j,ngrid2d) = grib_undef
        enddo
      enddo

!     if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  107) ',
!    & 'Snow sublimation (w/m^2) land surface'
!..........................................................
!  - snow cover
!
!$omp parallel do private(j,i)
      do j = 1, lats_node_r
        do i = 1, lonr
          glolal(i,j) = flx_fld%snowca(i,j) * rtime * 100.0
        enddo
      enddo

      ngrid2d = ngrid2d+1
      CALL uninterprez(2,kmsk,buffo,glolal,global_lats_r,lonsperlar,
     &     buff_mult_piecef(1,1,ngrid2d))

      do j=1,lats_node_r
        do i=1,lonr
          if (mskf(i,j) /= 1) buff_mult_piecef(i,j,ngrid2d) = grib_undef
        enddo
      enddo
!     if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  108) ',
!    & 'Snow cover (fraction) land surface'
!..........................................................
!  - total column soil moisture
!
!$omp parallel do private(j,i)
      do j = 1, lats_node_r
        do i = 1, lonr
          glolal(i,j) = flx_fld%soilm(i,j) * 1000.0 ! convert from m to (mm)kg/m^2
!       if (me==42 .and. i == 5 .and. j == 1)
!    & write(0,*)' flx_fld%soilm=',flx_fld%soilm(i,j),' mskf=',
!    & mskf(i,j)

        enddo
      enddo
      ngrid2d = ngrid2d+1
      CALL uninterprez(2,kmsk,buffo,glolal,global_lats_r,lonsperlar,
     &     buff_mult_piecef(1,1,ngrid2d))

      do j=1,lats_node_r
        do i=1,lonr
          if (mskf(i,j) /= 1) buff_mult_piecef(i,j,ngrid2d) = grib_undef
        enddo
      enddo

!     if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  109) ',
!    & 'Total column soil moisture (kg/m^2) land surface'
!..........................................................
!  - snow phase-change heat flux
!
!$omp parallel do private(j,i)
      do j = 1, lats_node_r
        do i = 1, lonr
          glolal(i,j) = flx_fld%snohfa(i,j) * rtime
        enddo
      enddo

      ngrid2d = ngrid2d+1
      CALL uninterprez(2,kmsk,buffo,glolal,global_lats_r,lonsperlar,
     &     buff_mult_piecef(1,1,ngrid2d))

!     if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  110) ',
!    & 'Snow phase-change heat flux [w/m^2] land surface'
!..........................................................
!  - wilting point

!$omp parallel do private(j,i)
      do j = 1, lats_node_r
        do i = 1, lonr
          glolal(i,j) = flx_fld%smcwlt2(i,j)
        enddo
      enddo
      ngrid2d = ngrid2d+1
      CALL uninterprez(2,kmsk,buffo,glolal,global_lats_r,lonsperlar,
     &     buff_mult_piecef(1,1,ngrid2d))

!     if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  111) ',
!    & 'Wilting point [fraction] land surface'
!..........................................................
!  - field capacity
!
!$omp parallel do private(j,i)
      do j = 1, lats_node_r
        do i = 1, lonr
          glolal(i,j) = flx_fld%smcref2(i,j)
        enddo
      enddo
      ngrid2d = ngrid2d+1
      CALL uninterprez(2,kmsk,buffo,glolal,global_lats_r,lonsperlar,
     &     buff_mult_piecef(1,1,ngrid2d))

!     if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  112) ',
!    & 'Field capacity [fraction] land surface'
!..........................................................
!  - accumulated sunshine duration time
!
!$omp parallel do private(j,i)
      do j = 1, lats_node_r
        do i = 1, lonr
          glolal(i,j) = flx_fld%suntim(i,j)
        enddo
      enddo
      ngrid2d = ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,glolal,global_lats_r,lonsperlar,
     &     buff_mult_piecef(1,1,ngrid2d))

!     if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  113) ',
!    & 'Accumulated sunshine duration (sec)'
!..........................................................
!! - frozen precipitation fraction

!$omp parallel do private(j,i)
      do j = 1, lats_node_r
        do i = 1, lonr
          glolal(i,j) = flx_fld%sr(i,j)
        enddo
      enddo

      ngrid2d = ngrid2d+1
      CALL uninterprez(2,kmsk,buffo,glolal,global_lats_r,lonsperlar,
     &                 buff_mult_piecef(1,1,ngrid2d))

!     if (ierr /= 0) print*,'wrtsfc gribit_gsm ierr=',ierr,'  114) ',
!    &  'Frozen precipitation fraction '
!..........................................................
!
!     write(0,*)' at the end in wrtflx_a max_ngrid2d=',ngrid2d

!     if(me == ioproc)
!    &  write(0,*)'(wrtflx_a) GRIB FLUX FILE WRITTEN ',FHOUR,IDATE,noflx
!!
      if (allocated(mskf)) deallocate(mskf)

      RETURN
      END

!!*********************************************************************
!! This routine is added to output 2d aerosol diag fields (Sarah Lu)

      SUBROUTINE wrtaer(IOPROC,noaer,ZHOUR,FHOUR,IDATE,
     &             sfc_fld, g2d_fld,global_lats_r, lonsperlar)
!!
      use resol_def,               ONLY: lonr, latr, ngrids_aer
      use mod_state,               ONLY: buff_mult_pieceg
      use layout1,                 ONLY: me, lats_node_r,lats_node_r_max
      use gfs_physics_sfc_flx_mod, ONLY: Sfc_Var_Data
      use gfs_physics_g2d_mod,     ONLY: G2D_Var_Data
      USE machine,                 ONLY: kind_io8, kind_io4
      implicit none
!!
      TYPE(Sfc_Var_Data)        :: sfc_fld
      TYPE(G2D_Var_Data)        :: g2d_fld
      INTEGER                   GLOBAL_LATS_R(LATR)
      INTEGER                   lonsperlar(LATR)
      integer                   IOPROC
!!
      integer                   i,j,k,l,noaer,ngrid2d,ierr
      real (kind=kind_io8)      rtime
      real (kind=kind_io8)      zhour,fhour

!     real(kind=kind_io8) slmskful(lonr,lats_node_r)
!     real(kind=kind_io8) slmskloc(LONR,LATS_NODE_R)
!
      INTEGER     IDATE(4), IDS(255),IENS(5)
!
      real (kind=kind_io8)   glolal(lonr,LATS_NODE_R_max)
      real (kind=kind_io8)   buffo(lonr,LATS_NODE_R_max)
      integer kmsk  (lonr,lats_node_r_max),kmsk0(lonr,lats_node_r_max)
!
      kmsk=nint(sfc_fld%slmsk)
      kmsk0=0
!
!     ngrid2d=1
!
      IF(FHOUR.GT.ZHOUR) THEN
        RTIME=1./(3600.*(FHOUR-ZHOUR))
      ELSE
        RTIME=0.
      ENDIF
!
!..........................................................
!
      ngrid2d = 0
      if ( g2d_fld%du%nfld > 0 ) then
        do  k = 1, g2d_fld%du%nfld
          glolal = RTIME*1.e6*g2d_fld%du%diag(k)%flds
          ngrid2d=ngrid2d+1
          CALL uninterprez(2,kmsk0,buffo,glolal,global_lats_r,
     &                   lonsperlar,buff_mult_pieceg(1,1,ngrid2d))
        enddo
      endif
!
!..........................................................
!
      if ( g2d_fld%su%nfld > 0 ) then
        do  k = 1, g2d_fld%su%nfld
          glolal = RTIME*1.e6*g2d_fld%su%diag(k)%flds
          ngrid2d=ngrid2d+1
          CALL uninterprez(2,kmsk0,buffo,glolal,global_lats_r,
     &                   lonsperlar,buff_mult_pieceg(1,1,ngrid2d))
        enddo
      endif
!
!..........................................................
!
      if ( g2d_fld%ss%nfld > 0 ) then
        do  k = 1, g2d_fld%ss%nfld
          glolal = RTIME*1.e6*g2d_fld%ss%diag(k)%flds
          ngrid2d=ngrid2d+1
          CALL uninterprez(2,kmsk0,buffo,glolal,global_lats_r,
     &                   lonsperlar,buff_mult_pieceg(1,1,ngrid2d))
        enddo
      endif
!
!..........................................................
!
      if ( g2d_fld%oc%nfld > 0 ) then
        do  k = 1, g2d_fld%oc%nfld
          glolal=RTIME*1.e6*g2d_fld%oc%diag(k)%flds
          ngrid2d=ngrid2d+1
          CALL uninterprez(2,kmsk0,buffo,glolal,global_lats_r,
     &                   lonsperlar,buff_mult_pieceg(1,1,ngrid2d))
        enddo
      endif
!
!..........................................................
!
      if ( g2d_fld%bc%nfld > 0 ) then
        do  k = 1, g2d_fld%bc%nfld
          glolal = RTIME*1.e6*g2d_fld%bc%diag(k)%flds
          ngrid2d=ngrid2d+1
          CALL uninterprez(2,kmsk0,buffo,glolal,global_lats_r,
     &                   lonsperlar,buff_mult_pieceg(1,1,ngrid2d))
        enddo
      endif
!
!..........................................................
! 2d met fields (k=01-10) are time-avg;
! 3d met fields (k=11-24) are instant
! this change makes comparison easier (flx for 2d, sig for 3d)
!
      if ( g2d_fld%met%nfld > 0 ) then
        do  k = 1, g2d_fld%met%nfld
          if (k .le. 10 ) then                      ! time-avg
             glolal=RTIME*g2d_fld%met%diag(k)%flds
          else                                      ! instant
             glolal=g2d_fld%met%diag(k)%flds
          endif
          ngrid2d=ngrid2d+1
          CALL uninterprez(2,kmsk0,buffo,glolal,global_lats_r,
     &                   lonsperlar,buff_mult_pieceg(1,1,ngrid2d))
        enddo
      endif

!!

      if(me.eq.ioproc)
     &   PRINT *,'(wrtaer) GRIB AER FILE WRITTEN ',FHOUR,IDATE,noaer
!!
      RETURN
      END
!!****

       subroutine flx_only_move(ioproc)
!
!***********************************************************************
!
      use resol_def, ONLY: ngrids_flx, ngrids_sfcc, lonr,latr
      use mod_state, ONLY: buff_mult_pieces, buff_mult_piecef,
     &                     ivar_global_a, ivar_global
      use layout1,   ONLY: me, nodes, ipt_lats_node_r, lats_node_r,
     &                     lats_node_r_max, nodes_comp
      use mpi_def,   ONLY: mpi_r_io, stat, mpi_comm_null, info, 
     &                     mc_comp, mpi_integer, mpi_comm_all, liope
      implicit none
!
      integer ipt_lats_node_rl,nodesr
      integer lats_nodes_rl
!      integer lats_nodes_r(nodes),ipt,maxfld,ioproc,nproct
      integer ioproc
      integer proc,j,lat,msgtag,nproc,i,msgtag1,buff,startlat,ierr
      integer illen,nd1
      integer icount
      data icount/0/
      save icount
!     integer maxlats_comp
!     save maxlats_comp,icount
      integer kllen
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      if(icount .eq. 0) then
        if(.not.allocated(ivar_global)) allocate(ivar_global(10))
        if(.not.allocated(ivar_global_a)) 
     &      allocate(ivar_global_a(10,nodes))
        ivar_global(1)=ipt_lats_node_r
        ivar_global(2)= lats_node_r
        ivar_global(3)=lats_node_r_max
        call mpi_gather(ivar_global,10,MPI_INTEGER,
     1    ivar_global_a,10,MPI_INTEGER,ioproc,MPI_COMM_ALL,ierr)
        icount=icount+1
      endif
!!
      if(allocated(buff_mult_pieces)) then
        deallocate(buff_mult_pieces)
!     else
!       maxlats_comp=lats_node_r_max
!       if(me .eq. ioproc) then
!         maxlats_comp=ivar_global_a(3,1)
!       endif
      endif
      if(me .eq. ioproc) then
!gwv watch this!!
!       write(0,*)' In flx_only_move ngrids_flx=',ngrids_flx
          allocate
     1  (buff_mult_pieces(lonr*latr*ngrids_flx))
         buff_mult_pieces=0.
       endif
!
!  big send
       IF (me.ne.ioproc) THEN
!
!         Sending the data
         msgtag=me
         illen=lats_node_r
         kllen=illen*lonr*ngrids_flx
! send the local grid domain
         CALL mpi_send
     &     (buff_mult_piecef(1:lonr,1:lats_node_r,1:ngrids_flx),
     &      kllen,MPI_R_IO,ioproc,msgtag,mc_comp,info)
      ELSE
        if( MC_COMP .ne. MPI_COMM_NULL) then
! iotask is also a compute task.  send is replaced with direct
!  array copy

!     write(0,*)' in flux move buff_mult=',buff_mult_piecef(1:10,2,69)
 
        if(nodes_comp==1) then
          buff_mult_pieces(1:lonr*lats_node_r*ngrids_flx)=
     1   reshape(buff_mult_piecef(1:lonr,1:lats_node_r,1:ngrids_flx),
     1     (/lonr*lats_node_r*ngrids_flx/) )
!     write(0,*)' in flux move buff_mult2=',
!    &          (buff_mult_pieces(68*lonr*latr+lonr+i),i=1,10)
        else

!  END COMPUTE TASKS PORTION OF LOGIC
!  receiving part of I/O task
 
!!
!!     for pes ioproc
        nd1=0
        DO proc=1,nodes_comp
         illen=ivar_global_a(2,proc)
         if (proc.ne.ioproc+1) then
           msgtag=proc-1
           kllen=illen*lonr*ngrids_flx
           CALL mpi_recv
     1       (buff_mult_pieces(nd1+1),kllen,MPI_R_IO,proc-1,
     &                msgtag,mc_comp,stat,info)
!     &                msgtag,MPI_COMM_ALL,stat,info)
         else
           buff_mult_pieces(nd1+1:nd1+lonr*illen*ngrids_flx)=
     1      reshape(buff_mult_piecef(1:lonr,1:illen,1:ngrids_flx),
     1       (/lonr*illen*ngrids_flx/) )
         endif
         nd1=nd1+illen*lonr*ngrids_flx
        enddo
       endif

      endif
!end ioproc
      ENDIF
!
      return
      end
!------------------------------------------------------------------------ 
      SUBROUTINE flx_wrt_nemsio(IOPROC,cfile,ZHOUR,FHOUR,idate
     &,                  global_lats_r,lonsperlar)
!!
      use nemsio_module, only: nemsio_open,nemsio_writerec,nemsio_close
     &  ,nemsio_gfile, nemsio_init,nemsio_finalize
      use resol_def,    ONLY: lonr, latr, levs,ngrids_flx,
     & ncld,ntrac,ntcw,ntoz,lsoil, ivssfc,thermodyn_id,sfcpress_id
      use layout1,      ONLY: me,idrt
      USE machine,      ONLY: kind_io8, kind_io4
!
      use gfs_physics_output, only : PHY_INT_STATE_ISCALAR,
     &    PHY_INT_STATE_RSCALAR,
     &    PHY_INT_STATE_1D_I,PHY_INT_STATE_1D_R,
     &    PHY_INT_STATE_2D_R_FLX
      implicit none
!!
      integer nw,IOPROC
      character*16 cfile,NAME2D
      real(kind=kind_io8) zhour,fhour
      integer idate(4),k,il, ngridss
!
      integer i,j,ndim3,N2DR,INDX,idate7(7),kount,nrec
      integer nfhour,nfminute,nfsecondn,nfsecondd
      logical  :: outtest
      integer ::nmetavari,nmetavarr,nmetavarl,nmetaaryi,nmetaaryr
      character(16),allocatable :: recname(:),reclevtyp(:)
      integer,allocatable :: reclev(:),itr(:)
      character(16),allocatable :: variname(:),varrname(:),
     &    aryiname(:),aryrname(:)
      integer,allocatable :: varival(:),aryilen(:),
     &    aryival(:,:)
      real(kind=kind_io4),allocatable    :: varrval(:)
      real(kind=kind_io4),allocatable    :: buff_mult(:,:,:),tmp(:)
      type(nemsio_gfile) gfileout
!

!!
      CHARACTER*8 labfix(4)
      real(kind=kind_io4) yhour
      integer,save:: version
      data version/200501/
      INTEGER              GLOBAL_LATS_R(latr), lonsperlar(latr)
!
      integer iret
      logical first
      save first
      save  recname, reclevtyp, reclev
      save nrec,nmetavari,nmetavarr,nmetaaryi,nmetaaryr,
     &     variname,varrname,aryiname,
     &     varival,varrval,aryilen,aryival
!jw     &     variname,varrname,aryiname,aryrname,
!jw     &     varival,aryilen,aryrlen,aryival,aryrval,varrval
      data first /.true./
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!    Build surface fields in to buff_mult
!
      if (me == ioproc) then
!
!     write(0,*)' begin of flx_wrt  ngrids_flx=',ngrids_flx

        allocate(buff_mult(lonr,latr,ngrids_flx))
        buff_mult = 0.
        do ngridss=1,ngrids_flx
          print *,' inside flx_wrt calling unsp ngridss=',ngridss
          call unsplit2z(ngridss,ngrids_flx,
     &                   buff_mult(1,1,ngridss),global_lats_r)
        enddo
!    Building surface field is done
!
        if (first) then
!write out nemsio sfc file:
          nrec=ngrids_flx
          kount=size(PHY_INT_STATE_ISCALAR,2)
          do i=1,kount
           if(trim(PHY_INT_STATE_ISCALAR(2,i)).eq.'OGFS_PHY'
     &     .or.trim(PHY_INT_STATE_ISCALAR(2,i)).eq.'OGFS_FLX')
     &        nmetavari=nmetavari+1
          enddo
          allocate(variname(nmetavari),varival(nmetavari))
          j=0
          do i=1,kount
           if(trim(PHY_INT_STATE_ISCALAR(2,i)).eq.'OGFS_PHY' .or.
     &      trim(PHY_INT_STATE_ISCALAR(2,i)).eq.'OGFS_FLX' )then
            j=j+1
            variname(j)=trim(PHY_INT_STATE_ISCALAR(1,i))
            if(trim(variname(j))=='latr') varival(j)=latr
            if(trim(variname(j))=='lonr') varival(j)=lonr
            if(trim(variname(j))=='levs') varival(j)=levs
            if(trim(variname(j))=='ntoz') varival(j)=ntoz
            if(trim(variname(j))=='ntcw') varival(j)=ntcw
            if(trim(variname(j))=='ncld') varival(j)=ncld
            if(trim(variname(j))=='ntrac') varival(j)=ntrac
            if(trim(variname(j))=='thermodyn_id')varival(j)=thermodyn_id
            if(trim(variname(j))=='sfcpress_id') varival(j)=sfcpress_id
            if(trim(variname(j))=='lsoil') varival(j)=lsoil
            if(trim(variname(j))=='idrt') varival(j)=idrt
           endif
          enddo
!!for real var::
          nmetavarr=0
          do i=1,kount
           if(trim(PHY_INT_STATE_RSCALAR(2,i)).eq.'OGFS_PHY'
     &     .or.trim(PHY_INT_STATE_RSCALAR(2,i)).eq.'OGFS_FLX')
     &     nmetavarr=nmetavarr+1
          enddo
          if(nmetavarr>0) then
            allocate(varrname(nmetavarr),varrval(nmetavarr))
            j=0
            do i=1,kount
             if(trim(PHY_INT_STATE_RSCALAR(2,i)).eq.'OGFS_PHY'
     &       .or.trim(PHY_INT_STATE_RSCALAR(2,i)).eq.'OGFS_FLX')then
               j=j+1
               varrname(j)=trim(PHY_INT_STATE_RSCALAR(1,i))
               if(trim(varrname(j))=='fhour') varrval(j)=fhour
               if(trim(varrname(j))=='zhour') varrval(j)=zhour
             endif
            enddo
          endif
!!for 1D ary::
          nmetaaryi=0
          do i=1,kount
           if(trim(PHY_INT_STATE_1D_I(2,i)).eq.'OGFS_PHY'
     &     .or.trim(PHY_INT_STATE_1D_I(2,i)).eq.'OGFS_FLX')
     &     nmetaaryi=nmetaaryi+1
          enddo
          allocate(aryiname(nmetaaryi),aryilen(nmetaaryi))
          j=0
          do i=1,kount
           if(trim(PHY_INT_STATE_1D_I(2,i)).eq.'OGFS_PHY'
     &     .or.trim(PHY_INT_STATE_1D_I(2,i)).eq.'OGFS_FLX')then
             j=j+1
             aryiname(j)=trim(PHY_INT_STATE_1D_I(1,i))
             if(trim(aryiname(j))=='IDATE') aryilen(j)=size(idate)
           endif
          enddo
          allocate(aryival(maxval(aryilen),nmetaaryi) )
          aryival(1:aryilen(1),1)=idate(:)
!!!for 1D real ary::
!          nmetaaryr=0
!          do i=1,kount
!           if(trim(PHY_INT_STATE_1D_R(2,i)).eq.'OGFS_PHY'
!     &     .or.trim(PHY_INT_STATE_1D_R(2,i)).eq.'OGFS_FLX')
!     &     nmetaaryr=nmetaaryr+1
!          enddo
!          allocate(aryrname(nmetaaryr),aryrlen(nmetaaryr))
!          do i=1,kount
!           if(trim(PHY_INT_STATE_1D_R(2,i)).eq.'OGFS_PHY')
!     &     .or.trim(PHY_INT_STATE_1D_R(2,i)).eq.'OGFS_FLX')then
!             aryrname(i)=trim(PHY_INT_STATE_1D_R(1,i))
!             if(i==1) aryrlen(i)=size(ak5)
!             if(i==2) aryrlen(i)=size(bk5)
!             if(i==3) aryrlen(i)=size(ck5)
!           endif
!          enddo
!          allocate(aryrval(maxval(aryrlen),nmetaaryr)
!          aryrval(1:aryrlen(1),1)=ak5(:)
!          aryrval(1:aryrlen(2),2)=bk5(:)
!          aryrval(1:aryrlen(3),2)=ck5(:)
!
!!for record name, levtyp and lev
          allocate (recname(nrec),reclevtyp(nrec),reclev(nrec))
          allocate (itr(nrec))
          N2DR=0
          itr=-99
          do i=1,kount
           if(trim(PHY_INT_STATE_2D_R_FLX(2,i)).eq.'OGFS_FLX')then
            N2DR=N2DR+1
            NAME2D=trim(PHY_INT_STATE_2D_R_FLX(1,i))
            INDX=INDEX(NAME2D,"_")
            if(indx>0) then
              recname(N2DR)=NAME2D(1:INDX-1)
            else
              recname(N2DR)=NAME2D
            endif
!
            reclevtyp(N2DR)=trim(trim(PHY_INT_STATE_2D_R_FLX(3,i)))
            reclev(N2DR)=1
!
!check time average
           if(INDEX(NAME2D,"_ave") >0) then
               itr(N2DR)=3
            elseif(INDEX(NAME2D,"_acc") >0) then
               itr(N2DR)=4
            elseif(INDEX(NAME2D,"_win") >0) then
               itr(N2DR)=2
            endif

           endif
          enddo
!
!end first
          first=.false.
         endif
!
        idate7=0
        idate7(1)=idate(4)
        idate7(2)=idate(2)
        idate7(3)=idate(3)
        idate7(4)=idate(1)
        idate7(7)=100           !denominator for second
!
        nfhour=int(fhour)     
        nfminute=int((fhour-nfhour)*60)
        nfsecondn=int(((fhour-nfhour)*3600-nfminute*60)*100)
        nfsecondd=100
!
        call nemsio_init()
!
        call nemsio_open(gfileout,trim(cfile),'write',
     &    iret = iret,
     &    modelname='GFS',gdatatype='grib',
     &    idate=idate7,nrec=nrec,
     &    dimx=lonr,dimy=latr,dimz=levs,ncldt=ncld,nmeta=5,
     &    nfhour=nfhour,nfminute=nfminute,nfsecondn=nfsecondn,
     &    nfsecondd=nfsecondd,
     &    extrameta=.true.,nmetavari=nmetavari,
     &    nmetavarr=nmetavarr,
     &    nmetaaryi=nmetaaryi,
     &    variname=variname,varival=varival,varrname=varrname,
     &    varrval=varrval,
     &    aryiname=aryiname,aryilen=aryilen,aryival=aryival,
     &    ntrac=ntrac,nsoil=lsoil,idrt=idrt,
     &    recname=recname,reclevtyp=reclevtyp,reclev=reclev)
!
        allocate(tmp(lonr*latr))
        yhour=zhour
        do i=1,nrec
          tmp(:)=reshape(buff_mult(:,:,i),(/lonr*latr/) )
          if(itr(i)==-99) then
            call nemsio_writerec(gfileout,i,tmp,iret=iret)
          else
            call nemsio_writerec(gfileout,i,tmp,iret=iret,itr=itr(i),
     &        zhour=yhour)
          endif
        enddo
        deallocate(tmp)
        deallocate(buff_mult)
!
        call nemsio_close(gfileout)
!end write pe
        call nemsio_finalize()
      endif
!
      print *,' end of flx_wrt '
      return
      end subroutine flx_wrt_nemsio
!------------------------------------------------------------------------
      subroutine wrtflx_w                                               &
!...................................
!  ---  inputs:
     &     ( ioproc,noflx,zhour,fhour,idate,colat1,secswr,seclwr,       &
     &       slmsk,global_lats_r,lonsperlar)
!  ---  outputs: ( none )
! =================   subprogram documentation block   ================
! !
!                                                                       !
!    this program writes out surface flux file in grib form.
!    !
!                                                                       !
!    usage:        call wrtflx_w
!    !
!                                                                       !
!    subprograms called:                                                !
!                  uninterpred, unsplit2d, idsdef, gribit_gsm, wryte    !
!                                                                       !
!    attributes:                                                        !
!      language:   fortran 90                                           !
!                                                                       !
!    external modules referenced:                                       !
!      'module resol_def             in 'resol_def.f'                   !
!      'module mod_state             in 'wrtout.f'                      !
!      'module layout1               in 'layout1.f'                     !
!      'module sig_io                in 'sig_io.f'                      !
!      'module namelist_def          in 'namelist_def.f'                !
!                                                                       !
!    Revision history                                                   !
!    Nov 2013 Sarah Lu  : remove commented-out lines related to aod     !
!                         and the lggfs3d option; correct the record    !
!                         number for 44 to 114                          !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!    input variables:                                            size   !
!      ioproc        : integer, processor id num                   1    !
!      noflx         : integer,                                    1    !
!      zhour         : real, accumulator reset time in hr.         1    !
!      fhour         : real, forecast time in hr.                  1    !
!      idate         : integer, init cond time (h,m,d,y)           4    !
!      colat1        : real, colatitude                            1    !
!      secswr        : real, sw radiation call duration in sec     1    !
!      seclwr        : real, lw radiation call duration in sec     1    !
!      slmsk         : real, sea-lane mask (lonr,lats_node_r)           !
!      global_lats_r : integer, index for global latitudes        latr  !
!      lonsperlar    : integer, num of pts on a given lat circle  latr  !
!                                                                       !
!  ======================  end of definitions  ======================= ! 
!
      use resol_def
      use mod_state
      use layout1
!      use sig_io
      use namelist_physics_def
      USE machine, ONLY: kind_evod, kind_io4,kind_io8
      implicit none
!!
!
!  ---  inputs:
      integer, intent(in) :: ioproc, noflx, idate(4),                   &
     &                       global_lats_r(latr), lonsperlar(latr)

      real (kind=kind_io8), intent(in) :: zhour, fhour, colat1,         &
     &                                    secswr, seclwr
      real (kind=kind_io8), intent(in) :: slmsk(lonr,lats_node_r)

!  ---  outputs: (none)

!  ---  parameters:
      integer, parameter :: nfld = 29
!!
      integer, parameter ::        iprs  =  1, ihgt  =  7, itemp = 11,  &
     &     itmx  = 15, itmn  = 16, iznlw = 33, imerw = 34, isphum= 51,  &
     &     ipwat = 54, ipcpr = 59, isr=   194, isnowd= 65, isnod = 66,  &
     &     icldf = 71, iccldf= 72, islmsk= 81, izorl = 83, ialbdo= 84,  &
     &     istc  = 86, iveg  = 87, irnof = 90, icemsk= 91, isik  = 92,  &
     &     ilhflx=121, ishflx=122, izws  =124, imws  =125, irst  =140,  &
     &     isoilm=144, iep   =145, icldwk=146, izgw  =147, imgw  =148,  &
     &     ighflx=155, icsusw=160, icsdsw=161, icsulw=162, icsdlw=163,  &
     &     iuswfc=160, idswfc=161, iulwfc=162, idlwfc=163, inswfc=164,  &
     &     inlwfc=165, idswvb=166, idswvd=167, idswnb=168, idswnd=169,  &
     &     icmm  =179, isuntm=191, isbs  =198, ievbs =199, ievcw =200,  &
     &     iuvbf =200, iuvbfc=201, idswf =204, idlwf =205, iqmx  =204,  &
     &     iqmn  =205, ichh  =208, itran =210, iuswf =211, iulwf =212,  &
     &     icpcpr=214, ismcwlt=219,ismcref=220,ihpbl =221, islo  =222,  &
     &     icnp  =223, istp  =224, ivtp  =225, isnohf=229, isrf  =235,  &
     &     isnc  =238, iust  =253

      integer, parameter ::        isfc  =  1, itoa  =  8, ielev =105,  &
     &     isglev=109, idbls =111, i2dbls=112, islc  =160, icolmn=200,  &
     &     iblbl =209, ibltl =210, ibllyr=211, ilcbl =212, ilctl =213,  &
     &     ilclyr=214, imcbl =222, imctl =223, imclyr=224, ihcbl =232,  &
     &     ihctl =233, ihclyr=234, icvbl =242, icvtl =243, icvlyr=244 

      integer, parameter ::        ifhour=  1, ifday =  2, inst  = 10,  &
     &                             iwin  =  2, iavg  =  3, iacc  =  4

!  ---  local variables:
      integer :: k, k4, l, il
      integer :: ilpds, iyr, imo, ida, ihr, ifhr, ithr, lg, ierr
      integer :: ids(255), ids_iq, ids_uvb
!     integer, dimension(lonr,lats_node_r) :: kmsk, kmsk0

      real (kind=kind_io8) :: rtimer(nfld), rtime, rtimsw, rtimlw
      real (kind=kind_io8) :: cl1, si(levp1)

      real (kind=kind_io4), dimension(lonr*latr)        :: wrkga

      real (kind=kind_io8), dimension(lonr*latr)        :: slmskful
      real (kind=kind_io8), dimension(lonr,lats_node_r) :: slmskloc

      logical(1) :: lbm(lonr*latr)
      character  :: g(200+lonr*latr*(16+1)/8)

!  ---  label indexes:
      integer, dimension(nfld) :: ipur, itlr
      data  ipur / iulwf , iuswf , iuswf , idswf , icldf , iprs  ,      &
     &             iprs  , itemp , icldf , iprs  , iprs  , itemp ,      &
     &             icldf , iprs  , iprs  , itemp , iuvbf , iuvbfc,      &
     &             idswf , icsulw, icsusw, icsdlw, icsusw, icsdsw,      &
     &             icsulw, idswvb, idswvd, idswnb, idswnd /

      data  itlr / itoa  , itoa  , isfc  , isfc  , ihclyr, ihctl ,      &
     &             ihcbl , ihctl , imclyr, imctl , imcbl , imctl ,      &
     &             ilclyr, ilctl , ilcbl , ilctl , isfc  , isfc  ,      &
     &             itoa  , itoa  , itoa  , isfc  , isfc  , isfc  ,      &
     &             isfc  , isfc  , isfc  , isfc  , isfc   /
!
!
      integer iens(5)
      integer icen,icen2,ienst,iensi
      integer ngridss
      integer i
!!
!===>  begin here
!
!     print *,'in wrtflx_w'

      icen=7 ; icen2=0 ; igen=82 ; ienst=0 ; iensi=0

      ngrid = 0 + ngrids_sfcc + 1 + ngrids_nst
!
      g   = ' '
!
      call unsplit2d_phys_r(ioproc,slmskful,slmsk,global_lats_r)
!
!  ---  set defalt decimal scaling factor array

      ids = 0
      ids_iq      = 5      ! used for iqmx/iqmn due to conflict with idswf/idlwf
      ids_uvb     = 2      ! used for iuvbf/iuvbfc due to conflict with ievcw
      CALL idsdef(1,ids)

!  ---  make adjustment if diff from the defaults
      ids(ihgt)   = 3      ! (007) geopotential height              (def=1)
      ids(itemp)  = 3      ! (011) temperature                      (def=1)
      ids(iznlw)  = 2      ! (033) zonal wind                       (def=1)
      ids(imerw)  = 2      ! (034) meridional wind                  (def=1)
      ids(isphum) = 6      ! (051) specific humidity                (def=5)
      ids(ipcpr)  = 6      ! (059) precipitation rate               (def=6)
      ids(isnowd) = 5      ! (065) water equivalent of snow depth   (def=0)
      ids(isnod)  = 6      ! (066) snow depth/mixed-layer depth     (def=2)
      ids(izorl)  = 4      ! (083) roughness                        (def=5)
      ids(istc)   = 2      ! (086) soil wetness                     (def=0) ! Moorthi
!     ids(istc)   = 4      ! (086) soil wetness                     (def=0)
      ids(iveg)   = 2      ! (087) vegetation/salinity              (def=0)
      ids(irnof)  = 5      ! (090) runoff (def=1)
      ids(icemsk) = 3      ! (091) ice concentration (def=2)
      ids(isik)   = 2      ! (092) ice thickness (def=2)
      ids(isoilm) = 4      ! (144) vol soil moisture content (def=3)
      ids(icmm)   = 4      ! (179) exchange coefficient            (def not set) *table 130
      ids(ievcw)  = 2      ! (200) canopy water evaporation        (def not set)*table 2
!     ids(iuvbf)  = 0      ! (200) sfc dnwd uv-b flux tot sky      (def not set) *table 129
!     ids(iuvbfc) = 0      ! (201) sfc dnwd uv-b flux clr sky      (def not set) *table 129
!     ids(iqmx)   = 5      ! (204) max specific humidity           (def=0) *** table 133
!     ids(iqmn)   = 5      ! (205) min specific humidity           (def=0) *** table 133
      ids(ichh)   = 4      ! (208) exchange coefficient            (def not set)
      ids(icpcpr) = 6      ! (214) convective precipitation rate   (def=6) ***
      ids(ismcwlt)= 4      ! (219) wilting point                   (def=1) *** table 130
      ids(ismcref)= 4      ! (220) field capacity                  (def=4 *** table 130
      ids(icnp)   = 5      ! (223) plant canopy surface water      (def=1)
      ids(isrf)   = 5      ! (235) storm surface runoff            (def=1)
      ids(isnc)   = 3      ! (238) snow cover                      (def=0)
      ids(iust)   = 3      ! (253) friction velocity              (def not set

      ilpds = 28
      if (icen2 == 2) ilpds = 45

      iens(1) = 1
      iens(2) = ienst
      iens(3) = iensi
      iens(4) = 1
      iens(5) = 255

      iyr     = idate(4)
      imo     = idate(2)
      ida     = idate(3)
      ihr     = idate(1)
      ifhr    = nint(zhour)
      ithr    = nint(fhour)

      if (fhour > zhour) then
        rtime = 1./(3600.*(fhour-zhour))
        rtime = 1./(3600.*(fhour-zhour))
      else
        rtime = 0.
      endif

      if (secswr > 0.) then
        rtimsw = 1./secswr
      else
        rtimsw = 1.
      endif

      if (seclwr > 0.) then
        rtimlw = 1./seclwr
      else
        rtimlw = 1.
      endif

!     write(0,*)' secswr=',secswr,' seclwr=',seclwr,' rtimsw=',rtimsw, &
!               ' rtimlw=',rtimlw
      rtimer     = rtimsw
      rtimer( 1) = rtimlw
      rtimer(20) = rtimlw       ! csulf_toa
      rtimer(22) = rtimlw       ! csdlf_sfc
      rtimer(25) = rtimlw       ! csulf_sfc
      cl1        = colat1
!jwang add spfhmax/spfhmin
!     ids(IQMX)   = 5
!     ids(IQMN)   = 5
! UV-B scaling factor, if set up already, comment the next 2 lines out
!     ids(IUVBF)  = 2
!     ids(IUVBFC) = 2
! Ice conentration and thickness scaling factor
!     ids(icemsk) = 3      ! ICE CONCENTRATION ()
!     ids(isik)   = 2      ! ICE THICKNESS (M)
!
!wei added 10/24/2006
!     ids(IZORL)  = 4
!     ids(IHGT)   = 3
!     ids(IVEG)   = 2
!     ids(IUST)   = 3
!     ids(ICHH)   = 4
!     ids(ICMM)   = 4
!     ids(ISRF)   = 5
!     ids(ITEMP)  = 3
!     ids(ISPHUM) = 6
!     ids(IZNLW)  = 2
!     ids(IMERW)  = 2
!     ids(ISNC)   = 3
!     ids(ISTC)   = 4
!     ids(ISOILM) = 4
!     ids(ISNOD)  = 6
!     ids(ISNOWD) = 5
!     ids(ICNP)   = 5
!     ids(IPCPR)  = 6
!     ids(ICPCPR) = 6
!     ids(IRNOF)  = 5
!     ids(ISMCWLT)  = 4
!     ids(ISMCREF)  = 4
!
!     ILPDS = 28
!     IF(ICEN2.EQ.2) ILPDS = 45
!     IENS(1) = 1
!     IENS(2) = IENST
!     IENS(3) = IENSI
!     IENS(4) = 1
!     IENS(5) = 255
!     IYR     = IDATE(4)
!     IMO     = IDATE(2)
!     IDA     = IDATE(3)
!     IHR     = IDATE(1)
!     IFHR    = NINT(ZHOUR)
!     ITHR    = NINT(FHOUR)
!     IF(FHOUR.GT.ZHOUR) THEN
!       RTIME = 1./(3600.*(FHOUR-ZHOUR))
!     ELSE
!       RTIME = 0.
!     ENDIF
!     IF(SECSWR.GT.0.) THEN
!       RTIMSW = 1./SECSWR
!     ELSE
!       RTIMSW=1.
!     ENDIF
!     IF(SECLWR.GT.0.) THEN
!       RTIMLW=1./SECLWR
!     ELSE
!       RTIMLW=1.
!     ENDIF
!     RTIMER=RTIMSW
!     RTIMER(1)=RTIMLW
!*RADFLX*
!     RTIMER(20)=RTIMLW       ! CSULF_TOA
!     RTIMER(22)=RTIMLW       ! CSDLF_SFC
!     RTIMER(25)=RTIMLW       ! CSULF_SFC
!*RADFLX*
!     CL1=colat1
!!
!..........................................................
!  - zonal component of momentum flux:
!
!     glolal = dusfc*rtime

!     print *,' begin of flx_wrt '
!     write(0,*)' begin of flx_wrt  ngrids_flx=',ngrids_flx

      ierr = 0

!
      if (me == ioproc) then
        ngridss = 1
        call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)
!      write(0,*)' wrkga=',maxval(wrkga),minval(wrkga)
        call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,igen,
     &              0,izws,isfc,0,0,iyr,imo,ida,ihr,
     &              ifhour,ifhr,ithr,iavg,0,0,icen2,ids(izws),iens,
     &              0.,0.,0.,0.,0.,0.,g,lg,ierr)

        if (ierr == 0) then
          call wryte(noflx,lg,g)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  01)',
     &     'Zonal compt of momentum flux (n/m**2) land and sea surface'
        endif
!..........................................................
!  - meridional component of momentum flux:
!
!     glolal = dvsfc*rtime

        ngridss = ngridss+1
        call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)
        call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,igen,
     &              0,imws,isfc,0,0,iyr,imo,ida,ihr,
     &              ifhour,ifhr,ithr,iavg,0,0,icen2,ids(imws),iens,
     &              0.,0.,0.,0.,0.,0.,g,lg,ierr)

        if (ierr == 0) then
          call wryte(noflx,lg,g)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  02) ',
     &     'Merid compt of momentum flux (n/m**2) land and sea surface'
        endif
!..........................................................
!  - sensible heat flux:
!
!     glolal = dtsfc*rtime

        ngridss = ngridss+1
        call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

        call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,igen,
     &              0,ishflx,isfc,0,0,iyr,imo,ida,ihr,
     &              ifhour,ifhr,ithr,iavg,0,0,icen2,ids(ishflx),iens,
     &              0.,0.,0.,0.,0.,0.,g,lg,ierr)

        if (ierr == 0) then
          call wryte(noflx,lg,g)
        else
            print*,'wrtsfc gribit ierr=',ierr,'  03) ',
     &       'Sensible heat flux (w/m**2) land and sea surface'
        endif
!..........................................................
!  - latent heat flux:
!
!     glolal = dqsfc*rtime

        ngridss = ngridss+1
        call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

        call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,igen,
     &              0,ilhflx,isfc,0,0,iyr,imo,ida,ihr,
     &              ifhour,ifhr,ithr,iavg,0,0,icen2,ids(ilhflx),iens,
     &              0.,0.,0.,0.,0.,0.,g,lg,ierr)

        if (ierr == 0) then
          call wryte(noflx,lg,g)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  04) ',
     &     'Latent heat flux (w/m**2) land and sea surface'
        endif
!..........................................................
!  - surface temperature:

        ngridss = ngridss+1
        call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

!     do i=1,lonr*latr
!       if (wrkga(i) < 10.0) then
!         write(0,*)' IN wrtout tsfc=',wrkga(i),' i=',i
!       endif
!     enddo
        call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,igen,
     &              0,itemp,isfc,0,0,iyr,imo,ida,ihr,
     &              ifhour,ithr,0,inst,0,0,icen2,ids(itemp),iens,
     &              0.,0.,0.,0.,0.,0.,g,lg,ierr)

        if (ierr == 0) then
          call wryte(noflx,lg,g)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  05) ',
     &     'Temperature (k) land and sea surface'
          stop
        endif
!..........................................................
!  - volumetric soil moist content at layer 10cm and 0cm:
!
!     glolal(:,:) = smc(:,1,:)

        ngridss = ngridss+1
        call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

        lbm = (slmskful == 1._kind_io8)
        call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,igen,
     &              1,isoilm,i2dbls,0,10,iyr,imo,ida,ihr,
     &              ifhour,ithr,0,inst,0,0,icen2,ids(isoilm),iens,
     &              0.,0.,0.,0.,0.,0.,g,lg,ierr)

        if (ierr == 0) then
          call wryte(noflx,lg,g)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  06) ',
     &     'Volumetric soil moist content (frac) layer 10cm and 0cm'
        endif
!..........................................................
!  - volumetric soil moist content at layer 40cm and 10cm:
!
!     glolal(:,:) = smc(:,2,:)

        ngridss = ngridss+1
        call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

        lbm = (slmskful == 1._kind_io8)
        call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,igen,
     &              1,isoilm,i2dbls,10,40,iyr,imo,ida,ihr,
     &              ifhour,ithr,0,inst,0,0,icen2,ids(isoilm),iens,
     &              0.,0.,0.,0.,0.,0.,g,lg,ierr)

        if (ierr == 0) then
          call wryte(noflx,lg,g)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  07) ',
     &     'Volumetric soil moist content (frac) layer 40cm and 10cm'
        endif
!..........................................................
!  - temperature at layer 10cm and 0cm:
!
!     glolal(:,:) = stc(:,1,:)

        ngridss=ngridss+1
        call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

        lbm = (slmskful == 1._kind_io8)
        call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,igen,
     &              1,itemp,i2dbls,0,10,iyr,imo,ida,ihr,
     &              ifhour,ithr,0,inst,0,0,icen2,ids(itemp),iens,
     &              0.,0.,0.,0.,0.,0.,g,lg,ierr)

        if (ierr == 0) then
          call wryte(noflx,lg,g)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  08) ',
     &    'Temp (k) layer betw two depth below land sfc 10cm and 0cm'
        endif
!..........................................................
!  - temperature at layer 40cm and 10cm:
!
!     glolal(:,:) = stc(:,2,:)

        ngridss = ngridss+1
        call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)


        lbm = (slmskful == 1._kind_io8)
        call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,igen,
     &              1,itemp,i2dbls,10,40,iyr,imo,ida,ihr,
     &              ifhour,ithr,0,inst,0,0,icen2,ids(itemp),iens,
     &              0.,0.,0.,0.,0.,0.,g,lg,ierr)

        if (ierr == 0) then
          call wryte(noflx,lg,g)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  09) ',
     &     'Temp (k) layer betw 2 depth below land sfc 40cm and 10cm'
        endif
!..........................................................
!  - water equivalent of accummulated snow depth:

        ngridss = ngridss+1
        call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

        lbm = (slmskful /= 0._kind_io8)
        call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,igen,
     &              1,isnowd,isfc,0,0,iyr,imo,ida,ihr,
     &              ifhour,ithr,0,inst,0,0,icen2,ids(isnowd),iens,
     &              0.,0.,0.,0.,0.,0.,g,lg,ierr)

        if (ierr == 0) then
          call wryte(noflx,lg,g)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  10) ',
     &     'Water equiv of accum snow depth (kg/m**2) land sea surface'
        endif
!..........................................................
!  - total sky radiation fluxes at toa and surface:

        do k = 1, 4

         ngridss = ngridss+1
         call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

         call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,igen,
     &               0,ipur(k),itlr(k),0,0,iyr,imo,ida,ihr,
     &               ifhour,ifhr,ithr,iavg,0,0,icen2,ids(ipur(k)),iens,
     &               0.,0.,0.,0.,0.,0.,g,lg,ierr)

         if (ierr == 0) then
           call wryte(noflx,lg,g)
         else
           if (k == 1) print*,'wrtsfc gribit ierr=',ierr,'  11) ',
     &      'Upward long wave radiation flux (w/m**2) at TOA'
           if (k == 2) print*,'wrtsfc gribit ierr=',ierr,'  12) ',
     &      'Upward solar radiation flux (w/m**2) at TOA'
           if (k == 3) print*,'wrtsfc gribit ierr=',ierr,'  13) ',
     &      'Upward solar radiation flux (w/m**2) at surface'
           if (k == 4) print*,'wrtsfc gribit ierr=',ierr,'  14) ',
     &      'Downward solar radiation flux (w/m**2) at surface'
         endif
        enddo
!..........................................................
!  - for high, mid, low cloud (cover, pressure, temperature)
!
        lab_do_cloud : do k = 5, 7    ! (high, mid, low clouds)
         k4 = 4 + (k-5)*4

!  - cloud cover (h,m,l):

         ngridss = ngridss+1
         call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

         l = k4 + 1
         lbm = (wrkga >= 0.5_kind_io8)
         call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,igen,
     &               0,ipur(l),itlr(l),0,0,iyr,imo,ida,ihr,
     &               ifhour,ifhr,ithr,iavg,0,0,icen2,ids(ipur(l)),iens,
     &               0.,0.,0.,0.,0.,0.,g,lg,ierr)

         if (ierr == 0) then
           call wryte(noflx,lg,g)
         else
           if (k == 5) print*,'wrtsfc gribit ierr=',ierr,'  15) ',
     &      'Total cloud cover (percent) high cloud layer'
           if (k == 6) print*,'wrtsfc gribit ierr=',ierr,'  19) ',
     &      'Total cloud cover (percent) middle cloud layer'
           if (k == 7) print*,'wrtsfc gribit ierr=',ierr,'  23) ',
     &      'Total cloud cover (percent) low cloud layer'
         endif

!  - pressure at cloud top:

         ngridss = ngridss+1
         call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

         l = k4 + 2
         call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,igen,
     &               1,ipur(l),itlr(l),0,0,iyr,imo,ida,ihr,
     &               ifhour,ifhr,ithr,iavg,0,0,icen2,ids(ipur(l)),iens,
     &               0.,0.,0.,0.,0.,0.,g,lg,ierr)

         if (ierr == 0) then
           call wryte(noflx,lg,g)
         else
           if (k == 5) print*,'wrtsfc gribit ierr=',ierr,'  16) ',
     &      'Pressure (pa) high cloud top level'
           if (k == 6) print*,'wrtsfc gribit ierr=',ierr,'  20) ',
     &      'Pressure (pa) middle cloud top level'
           if (k == 7) print*,'wrtsfc gribit ierr=',ierr,'  24) ',
     &      'Pressure (pa) low cloud top level'
         endif

!  - pressure at cloud base:

         ngridss = ngridss+1
         call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)


!..........................................................
!.......  FIX FLUXES FOR APPROX DIURNAL CYCLE
         l = k4 + 3
         call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,igen,
     &               1,ipur(l),itlr(l),0,0,iyr,imo,ida,ihr,
     &               ifhour,ifhr,ithr,iavg,0,0,icen2,ids(ipur(l)),iens,
     &               0.,0.,0.,0.,0.,0.,g,lg,ierr)

         if (ierr == 0) then
           call wryte(noflx,lg,g)
         else
           if (k == 5) print*,'wrtsfc gribit ierr=',ierr,'  17) ',
     &      'Pressure (pa) high cloud bottom level'
           if (k == 6) print*,'wrtsfc gribit ierr=',ierr,'  21) ',
     &      'Pressure (pa) middle cloud bottom level'
           if (k == 7) print*,'wrtsfc gribit ierr=',ierr,'  25) ',
     &      'Pressure (pa) low cloud bottom level'
         endif

!  - temperature at cloud top:

         ngridss = ngridss+1
         call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

         l = k4 + 4
         call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,igen,
     &               1,ipur(l),itlr(l),0,0,iyr,imo,ida,ihr,
     &               ifhour,ifhr,ithr,iavg,0,0,icen2,ids(ipur(l)),iens,
     &               0.,0.,0.,0.,0.,0.,g,lg,ierr)

         if (ierr == 0) then
           call wryte(noflx,lg,g)
         else
           if (k == 5) print*,'wrtsfc gribit ierr=',ierr,'  18) ',
     &      'Temperature (k) high cloud top level'
           if (k == 6) print*,'wrtsfc gribit ierr=',ierr,'  22) ',
     &      'Temperature (k) middle cloud top level'
           if (k == 7) print*,'wrtsfc gribit ierr=',ierr,'  26) ',
     &      'Temperature (k) low cloud top level'
         endif

        enddo  lab_do_cloud

!...................................................................
!  - total cloud amount:

        ngridss = ngridss+1
        call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

        call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,igen,
     &              0,icldf,icolmn,0,0,iyr,imo,ida,ihr,
     &              ifhour,ifhr,ithr,iavg,0,0,icen2,ids(icldf),iens,
     &              0.,0.,0.,0.,0.,0.,g,lg,ierr)

        if (ierr == 0) then
          call wryte(noflx,lg,g)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  27) ',
     &     'Total cloud cover (percent) total atmospheric column'
        endif
!     write(0,*)' gribit 27 ngridss aft tot cloud- ',ngridss
!.................................................
!  - boundary layer cloud amount:

        ngridss=ngridss+1
        call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

        call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,igen,
     &              0,icldf,ibllyr,0,0,iyr,imo,ida,ihr,
     &              ifhour,ifhr,ithr,iavg,0,0,icen2,ids(icldf),iens,
     &              0.,0.,0.,0.,0.,0.,g,lg,ierr)

        if (ierr == 0) then
          call wryte(noflx,lg,g)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  28) ',
     &     'Total cloud cover (percent) boundary layer cloud layer'
        endif
!..........................................................
!  - surface downeard lw fluxes: (use the surface temp adjusted quantities
!    to replace the original on in rec 19 of fluxr)   For UV-B fluxes
!
        ngridss = ngridss+1
        call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

        call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,igen,
     &              0,idlwf,isfc,0,0,iyr,imo,ida,ihr,
     &              ifhour,ifhr,ithr,iavg,0,0,icen2,ids(idlwf),iens,
     &              0.,0.,0.,0.,0.,0.,g,lg,ierr)

        if (ierr == 0) then
          call wryte(noflx,lg,g)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  29) ',
     &     'Downward long wave radiation flux (w/m**2) at surface'
        endif
!..........................................................
!  - surface upward lw fluxes: (use the one recalc'ed from surface temp
!    to replace the original one in rec 20 of fluxr)

        ngridss = ngridss+1
        call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

        call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,igen,
     &              0,iulwf,isfc,0,0,iyr,imo,ida,ihr,
     &              ifhour,ifhr,ithr,iavg,0,0,icen2,ids(iulwf),iens,
     &              0.,0.,0.,0.,0.,0.,g,lg,ierr)

        if (ierr == 0) then
          call wryte(noflx,lg,g)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  30) ',
     &     'Upward long wave radiation flux (w/m**2) at surface'
        endif
!..........................................................
!  - uv-b flux at surface for total sky:

        ngridss = ngridss+1
        call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

        call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,129,icen,
     &             igen,0,iuvbf,isfc,0,0,iyr,imo,ida,ihr,
     &             ifhour,ifhr,ithr,iavg,0,0,icen2,ids_uvb,iens,
     &             0.,0.,0.,0.,0.,0.,g,lg,ierr)

        if (ierr == 0) then
          call wryte(noflx,lg,g)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  31) ',
     &     'UV-B downward solar flux (w/m**2) at surface'
        endif
!..........................................................
!  - uv-b flux at surface for clear sky:

        ngridss = ngridss+1
        call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

        call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,129,icen,
     &             igen,0,iuvbfc,isfc,0,0,iyr,imo,ida,ihr,
     &             ifhour,ifhr,ithr,iavg,0,0,icen2,ids_uvb,iens,
     &             0.,0.,0.,0.,0.,0.,g,lg,ierr)

        if (ierr == 0) then
          call wryte(noflx,lg,g)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  32) ',
     &     'Clear sky UV-B downward solar flux (w/m**2) at surface'
        endif
!..........................................................
!  - incoming solar radiation at toa:

        ngridss = ngridss+1
        call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

        call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,igen,
     &              0,ipur(19),itlr(19),0,0,iyr,imo,ida,ihr,
     &              ifhour,ifhr,ithr,iavg,0,0,icen2,ids(ipur(19)),iens,
     &              0.,0.,0.,0.,0.,0.,g,lg,ierr)

        if (ierr == 0) then
          call wryte(noflx,lg,g)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  33) ',                   &
     &     'Downward solar radiation flux (W/m**2) at TOA'
        endif
!..........................................................
!  - sw downward surface flux components:

        do l = 24, 27
          k = l + 2

          ngridss = ngridss+1
          call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

          call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,
     &               igen,0,ipur(k),itlr(k),0,0,iyr,imo,ida,ihr,
     &               ifhour,ifhr,ithr,iavg,0,0,icen2,ids(ipur(k)),iens,
     &               0.,0.,0.,0.,0.,0.,g,lg,ierr)

          if (ierr == 0) then
            call wryte(noflx,lg,g)
          else
            if (l == 24) print*,'wrtsfc gribit ierr=',ierr,'  34) ',
     &       'Downward sw uv+vis beam radiation flux (w/m**2) sfc '
            if (l == 25) print*,'wrtsfc gribit ierr=',ierr,'  35) ',
     &       'Downward sw uv+vis diffuse radiation flux (w/m**2) sfc'
            if (l == 26) print*,'wrtsfc gribit ierr=',ierr,'  36) ',
     &       'Downward sw nir beam radiation flux (w/m**2) sfc '
            if (l == 27) print*,'wrtsfc gribit ierr=',ierr,'  37) ',
     &       'Downward sw nir diffuse radiation flux (w/m**2) sfc '
          endif
        enddo
!..........................................................
!  -  clear sky radiative fluxes at toa and surface:

        do l = 28, 33
          k = l - 8

          ngridss = ngridss+1
          call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

!         write(0,*)' Clsky=',(wrkga(i),i=1,lonr*latr,lonr),' l=',l
          call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,
     &               igen,0,ipur(k),itlr(k),0,0,iyr,imo,ida,ihr,
     &               ifhour,ifhr,ithr,iavg,0,0,icen2,ids(ipur(k)),iens,
     &               0.,0.,0.,0.,0.,0.,g,lg,ierr)

          if (ierr == 0) then
            call wryte(noflx,lg,g)
          else
            if (l == 28) print*,'wrtsfc gribit ierr=',ierr,'  38) ',
     &       'CS upward long wave radiation flux (w/m**2) at TOA'
            if (l == 29) print*,'wrtsfc gribit ierr=',ierr,'  39) ',
     &       'CS upward solar radiation flux (w/m**2) at TOA'
            if (l == 30) print*,'wrtsfc gribit ierr=',ierr,'  40) ',
     &       'CS downward long radiation flux (w/m**2) at surface'
            if (l == 31) print*,'wrtsfc gribit ierr=',ierr,'  41) ',
     &       'CS upward solar radiation flux (w/m**2)  at surface'
            if (l == 32) print*,'wrtsfc gribit ierr=',ierr,'  42) ',
     &       'CS downward solar radiation flux (w/m**2) at surface'
            if (l == 33) print*,'wrtsfc gribit ierr=',ierr,'  43) ',
     &       'CS upward long wave radiation (w/m**2) at surface'
          endif
        enddo
!...................................................................
!  - surface albedo (derived from radiative fluxes at surface):

        ngridss = ngridss+1
        call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

!      write(0,*)' bef albedo gribit ngridss=',ngridss

!       write(0,*)' albedo_gsm=',(wrkga(i),i=1,lonr*latr,lonr)
        call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,igen,
     &              0,ialbdo,isfc,0,0,iyr,imo,ida,ihr,
     &              ifhour,ifhr,ithr,iavg,0,0,icen2,ids(ialbdo),iens,
     &              0.,0.,0.,0.,0.,0.,g,lg,ierr)

        if (ierr == 0) then
          call wryte(noflx,lg,g)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  44) ',
     &     'Albedo (percent) land and sea surface '
        endif
!...................................................................
!  - precipitation rate (geshem unit in m, final unit = mm/s = kg/m2/s)

        ngridss = ngridss+1
        call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

        call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,igen,
     &              0,ipcpr,isfc,0,0,iyr,imo,ida,ihr,
     &              ifhour,ifhr,ithr,iavg,0,0,icen2,ids(ipcpr),iens,
     &              0.,0.,0.,0.,0.,0.,g,lg,ierr)

        if (ierr == 0) then
          call wryte(noflx,lg,g)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  45) ',
     &     'Precipitation rate (kg/m**2/s) land and sea surface'
        endif
!...................................................................
!  - convective precipitation rate:

        ngridss = ngridss+1
        call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

        call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,igen,
     &              0,icpcpr,isfc,0,0,iyr,imo,ida,ihr,
     &              ifhour,ifhr,ithr,iavg,0,0,icen2,ids(icpcpr),iens,
     &              0.,0.,0.,0.,0.,0.,g,lg,ierr)

        if (ierr == 0) then
          call wryte(noflx,lg,g)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  46) ',
     &     'Convective precipitation rate (kg/m**2/s) at surface'
        endif
!...................................................................
!  - ground heat flux:

        ngridss = ngridss+1
        call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

        lbm = (slmskful /= 0._kind_io8)
        call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,igen,
     &              1,ighflx,isfc,0,0,iyr,imo,ida,ihr,
     &              ifhour,ifhr,ithr,iavg,0,0,icen2,ids(ighflx),iens,
     &              0.,0.,0.,0.,0.,0.,g,lg,ierr)

        if (ierr == 0) then
          call wryte(noflx,lg,g)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  47) ',
     &     'Ground heat flux (w/m**2) land and sea surface'
        endif
!...................................................................
!  - land-sea mask:

        ngridss = ngridss+1
        call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

        call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,igen,
     &              0,islmsk,isfc,0,0,iyr,imo,ida,ihr,
     &              ifhour,ithr,0,inst,0,0,icen2,ids(islmsk),iens,
     &              0.,0.,0.,0.,0.,0.,g,lg,ierr)

        if (ierr == 0) then
          call wryte(noflx,lg,g)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  48) ',
     &     'Land-Sea mask (1=land; 0=sea) (integer) land sea surface'
        endif
!...................................................................
!  - sea-ice concentration:

        ngridss = ngridss+1
        call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

        call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,igen,
     &              0,icemsk,isfc,0,0,iyr,imo,ida,ihr,
     &              ifhour,ithr,0,inst,0,0,icen2,ids(icemsk),iens,
     &              0.,0.,0.,0.,0.,0.,g,lg,ierr)

        if (ierr == 0) then
          call wryte(noflx,lg,g)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  49) ',
     &     'Ice concentration (ice>0; no ice=0) (1/0) land sea surface'
        endif
!...................................................................
!  - 10m u wind:

        ngridss = ngridss+1
        call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

        call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,igen,
     &              0,iznlw,ielev,0,10,iyr,imo,ida,ihr,
     &              ifhour,ithr,0,inst,0,0,icen2,ids(iznlw),iens,
     &              0.,0.,0.,0.,0.,0.,g,lg,ierr)

        if (ierr == 0) then
          call wryte(noflx,lg,g)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  50) ',
     &     'u wind (m/s) height above ground'
        endif
!...................................................................
!  - 10m v wind:

        ngridss = ngridss+1
        call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

        call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,igen,
     &              0,imerw,ielev,0,10,iyr,imo,ida,ihr,
     &              ifhour,ithr,0,inst,0,0,icen2,ids(imerw),iens,
     &              0.,0.,0.,0.,0.,0.,g,lg,ierr)

        if (ierr == 0) then
          call wryte(noflx,lg,g)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  51) ',
     &     'v wind (m/s) height above ground'
        endif
!...................................................................
!  - 2m temperature:

        ngridss = ngridss+1
        call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

        call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,igen,
     &              0,itemp,ielev,0,2,iyr,imo,ida,ihr,
     &              ifhour,ithr,0,inst,0,0,icen2,ids(itemp),iens,
     &              0.,0.,0.,0.,0.,0.,g,lg,ierr)

        if (ierr == 0) then
          call wryte(noflx,lg,g)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  52) ',
     &     'Temperature (k) height above ground'
        endif
!...................................................................
!  - 2m specific humidity:

        ngridss = ngridss+1
        call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

        call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,igen,
     &              0,isphum,ielev,0,2,iyr,imo,ida,ihr,
     &              ifhour,ithr,0,inst,0,0,icen2,ids(isphum),iens,
     &              0.,0.,0.,0.,0.,0.,g,lg,ierr)

        if (ierr == 0) then
          call wryte(noflx,lg,g)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  53) ',
     &     'Specific humidity (kg/kg) height above ground'
        endif
!...................................................................
!  - surface pressure:

        ngridss = ngridss+1
        call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

        call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,igen,
     &              0,iprs,isfc,0,0,iyr,imo,ida,ihr,
     &              ifhour,ithr,0,inst,0,0,icen2,ids(iprs),iens,
     &              0.,0.,0.,0.,0.,0.,g,lg,ierr)

        if (ierr == 0) then
          call wryte(noflx,lg,g)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  54) ',
     &     'Pressure (pa) land and sea surface'
        endif
!...................................................................
!  - maximum temperature:

        ngridss = ngridss+1
        call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

        call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,igen,
     &              0,itmx,ielev,0,2,iyr,imo,ida,ihr,
     &              ifhour,ifhr,ithr,iwin,0,0,icen2,ids(itmx),iens,
     &              0.,0.,0.,0.,0.,0.,g,lg,ierr)

        if (ierr == 0) then
          call wryte(noflx,lg,g)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  55) ',
     &     'Maximum temperature (k) height above ground'
        endif
!...................................................................
!  - minimum temperature:

        ngridss = ngridss+1
        call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

        call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,igen,
     &              0,itmn,ielev,0,2,iyr,imo,ida,ihr,
     &              ifhour,ifhr,ithr,iwin,0,0,icen2,ids(itmn),iens,
     &              0.,0.,0.,0.,0.,0.,g,lg,ierr)

        if (ierr == 0) then
          call wryte(noflx,lg,g)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  56) ',
     &     'Minimum temperature (k) height above ground'
        endif
!...................................................................
!  - maximum specific humidity

        ngridss = ngridss+1
        call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

        call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,133,icen,
     &             igen,0,iqmx,ielev,0,2,iyr,imo,ida,ihr,
     &             ifhour,ifhr,ithr,iwin,0,0,icen2,ids_iq,iens,
     &             0.,0.,0.,0.,0.,0.,g,lg,ierr)

        if (ierr == 0) then
          call wryte(noflx,lg,g)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  57) ',
     &     'Maximum specific humidity (kg/kg) height above ground'
        endif
!...................................................................
!  - minimum specific humidity

        ngridss = ngridss+1
        call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

       call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,133,icen,igen,
     &             0,iqmn,ielev,0,2,iyr,imo,ida,ihr,
     &             ifhour,ifhr,ithr,iwin,0,0,icen2,ids_iq,iens,
     &             0.,0.,0.,0.,0.,0.,g,lg,ierr)

       if (ierr == 0) then
         call wryte(noflx,lg,g)
       else
         print*,'wrtsfc gribit ierr=',ierr,'  58) ',
     &    'Minimum specific humidity (kg/kg) height above ground'
       endif
!...................................................................
!  - runoff, the output unit of runoff is kg/m2 (accumulative value)

        ngridss = ngridss+1
        call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

        lbm = (slmskful /= 0._kind_io8)
        call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,igen,
     &              1,irnof,isfc,0,0,iyr,imo,ida,ihr,
     &              ifhour,ifhr,ithr,iacc,0,0,icen2,ids(irnof),iens,
     &              0.,0.,0.,0.,0.,0.,g,lg,ierr)

        if (ierr == 0) then
          call wryte(noflx,lg,g)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  59) ',
     &     'Runoff (kg/m**2) land and sea surface'
        endif
!...................................................................
!  - potential evaporation rate

        ngridss = ngridss+1
        call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

        lbm = (slmskful /= 0._kind_io8)
        call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,igen,
     &              1,iep,isfc,0,0,iyr,imo,ida,ihr,
     &              ifhour,ifhr,ithr,iavg,0,0,icen2,ids(iep),iens,
     &              0.,0.,0.,0.,0.,0.,g,lg,ierr)

        if (ierr == 0) then
          call wryte(noflx,lg,g)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  60) ',
     &     'Potential evaporation rate (w/m**/) land and sea surface'
        endif
!...................................................................
!  - cloud work function

        ngridss = ngridss+1
        call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

        call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,igen,
     &              0,icldwk,icolmn,0,0,iyr,imo,ida,ihr,
     &              ifhour,ifhr,ithr,iavg,0,0,icen2,ids(icldwk),iens,
     &              0.,0.,0.,0.,0.,0.,g,lg,ierr)

        if (ierr == 0) then
          call wryte(noflx,lg,g)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  61) ',
     &     'Cloud work function (j/kg) total atmospheric column'
        endif
!...................................................................
!  - zonal gravity wave stress

        ngridss = ngridss+1
        call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

        call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,igen,
     &              0,izgw,isfc,0,0,iyr,imo,ida,ihr,
     &              ifhour,ifhr,ithr,iavg,0,0,icen2,ids(izgw),iens,
     &              0.,0.,0.,0.,0.,0.,g,lg,ierr)

        if (ierr == 0) then
          call wryte(noflx,lg,g)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  62) ',
     &     'Zonal gravity wave stress (n/m**2) land and sea surface'
        endif
!...................................................................
!  - meridional gravity wave stress

        ngridss = ngridss+1
        call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

        call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,igen,
     &              0,imgw,isfc,0,0,iyr,imo,ida,ihr,
     &              ifhour,ifhr,ithr,iavg,0,0,icen2,ids(imgw),iens,
     &              0.,0.,0.,0.,0.,0.,g,lg,ierr)

        if (ierr == 0) then
          call wryte(noflx,lg,g)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  63) ',
     &     'Meridional gravity wave stress (n/m**2) land sea surface'
        endif

!...................................................................
!  - boundary layer height

        ngridss = ngridss+1
        call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

        call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,igen,
     &              0,ihpbl,isfc,0,0,iyr,imo,ida,ihr,
     &              ifhour,ithr,0,inst,0,0,icen2,ids(ihpbl),iens,
     &              0.,0.,0.,0.,0.,0.,g,lg,ierr)

        if (ierr == 0) then
          call wryte(noflx,lg,g)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  64) ',
     &     'Boundary layer height'
        endif
!...................................................................
!  - precipitable water

        ngridss = ngridss+1
        call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

        call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,igen,
     &              0,ipwat,icolmn,0,0,iyr,imo,ida,ihr,
     &              ifhour,ithr,0,inst,0,0,icen2,ids(ipwat),iens,
     &              0.,0.,0.,0.,0.,0.,g,lg,ierr)

        if (ierr == 0) then
          call wryte(noflx,lg,g)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  65) ',
     &     'Precipitable water (kg/m**2) total atmospheric column'
        endif
!...................................................................
!  - convective clouds
!    * labeled instantaneous but actually averaged over fhswr seconds

        ngridss = ngridss+1
        call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

        lbm = (wrkga >= 0.5_kind_io8)
        call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,igen,
     &              0,icldf,icvlyr,0,0,iyr,imo,ida,ihr,
     &              ifhour,ithr,0,inst,0,0,icen2,ids(icldf),iens,
     &              0.,0.,0.,0.,0.,0.,g,lg,ierr)

        if (ierr == 0) then
          call wryte(noflx,lg,g)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  66) ',
     &     'Total cloud cover (percent) convective cloud layer'
        endif
!.................................................
!  - pressure at convective cloud top

        ngridss = ngridss+1
        call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

        call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,igen,
     &              1,iprs,icvtl,0,0,iyr,imo,ida,ihr,
     &              ifhour,ithr,0,inst,0,0,icen2,ids(iprs),iens,
     &              0.,0.,0.,0.,0.,0.,g,lg,ierr)

        if (ierr == 0) then
          call wryte(noflx,lg,g)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  67) ',
     &     'Pressure (pa) convective cloud top level'
        endif
!.................................................
!  - pressure at convective cloud bottom

        ngridss = ngridss+1
        call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

        call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,igen,
     &              1,iprs,icvbl,0,0,iyr,imo,ida,ihr,
     &              ifhour,ithr,0,inst,0,0,icen2,ids(iprs),iens,
     &              0.,0.,0.,0.,0.,0.,g,lg,ierr)

        if (ierr == 0) then
          call wryte(noflx,lg,g)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  68) ',
     &     'Pressure (pa) convective cloud bottom level'
        endif
!.................................................
!  - sea ice thickness

        ngridss = ngridss+1
        call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

        lbm = (slmskful /= 1._kind_io8)
!       lbm = (slmskful == 2._kind_io8)
!     write(0,*)' calling gribit_gsm for sea ice thickness'
!     write(0,*)' sea ice=',wrkga(lonr+1:lonr+10),' ngridss=',ngridss
        call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,igen,
     &              1,isik,isfc,0,0,iyr,imo,ida,ihr,
     &              ifhour,ithr,0,inst,0,0,icen2,ids(isik),iens,
     &              0.,0.,0.,0.,0.,0.,g,lg,ierr)

        if (ierr == 0) then
          call wryte(noflx,lg,g)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  69) ',
     &     'Sea ice thickness (m) category 1'
        endif
!.................................................
!  - volumetric soil moist content (layer 100cm and 40cm)

        if (lsoil > 2) then

          ngridss = ngridss+1
          call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

!     write(0,*)' soil moist=',wrkga(lonr+1:lonr+10),' ngridss=',ngridss

          lbm = (slmskful == 1._kind_io8)
          call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,
     &               igen,1,isoilm,i2dbls,40,100,iyr,imo,ida,ihr,
     &               ifhour,ithr,0,inst,0,0,icen2,ids(isoilm),iens,
     &               0.,0.,0.,0.,0.,0.,g,lg,ierr)

          if (ierr == 0) then
            call wryte(noflx,lg,g)
          else
            print*,'wrtsfc gribit ierr=',ierr,'  70) ',
     &       '71) Volumetric soil moist content (frac) 100cm and 40cm'
          endif
!..........................................................
!  - volumetric soil moist content (layer 200cm and 100cm)

          ngridss = ngridss+1
          call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

          lbm = (slmskful == 1._kind_io8)
          call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,
     &               igen,1,isoilm,i2dbls,100,200,iyr,imo,ida,ihr,
     &               ifhour,ithr,0,inst,0,0,icen2,ids(isoilm),iens,
     &               0.,0.,0.,0.,0.,0.,g,lg,ierr)

          if (ierr == 0) then
            call wryte(noflx,lg,g)
          else
            print*,'wrtsfc gribit ierr=',ierr,'  71) ',
     &       'Volumetric soil moist content (frac) 200cm and 100cm'
          endif
!..........................................................
!  - temperature for layer 100cm and 40cm below sfc

          ngridss = ngridss+1
          call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

          lbm = (slmskful == 1._kind_io8)
          call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,
     &                igen,1,itemp,i2dbls,40,100,iyr,imo,ida,ihr,
     &                ifhour,ithr,0,inst,0,0,icen2,ids(itemp),iens,
     &                0.,0.,0.,0.,0.,0.,g,lg,ierr)

          if (ierr == 0) then
            call wryte(noflx,lg,g)
          else
            print*,'wrtsfc gribit ierr=',ierr,'  72) ',
     &       'Temp (k) layer betw two depth below sfc 100cm and 40cm'
          endif
!..........................................................
!  - temperature for layer 200cm and 100cm below sfc
!
          ngridss = ngridss+1
          call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

          lbm = (slmskful == 1._kind_io8)
          call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,
     &                igen,1,itemp,i2dbls,100,200,iyr,imo,ida,ihr,
     &                ifhour,ithr,0,inst,0,0,icen2,ids(itemp),iens,
     &                0.,0.,0.,0.,0.,0.,g,lg,ierr)

          if (ierr == 0) then
            call wryte(noflx,lg,g)
          else
            print*,'wrtsfc gribit ierr=',ierr,'  73) ',
     &       'Temp (k) layer betw two depth below sfc 200cm and 100cm'
          endif
        endif   ! end_if_lsoil
!..........................................................
!  - liquid soil moist content layer 40cm and 10cm
!
        ngridss = ngridss+1
        call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

        lbm = (slmskful == 1._kind_io8)
        call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,130,icen,
     &             igen,1,islc,i2dbls,0,10,iyr,imo,ida,ihr,
     &             ifhour,ithr,0,inst,0,0,icen2,ids(isoilm),iens,
     &             0.,0.,0.,0.,0.,0.,g,lg,ierr)

        if (ierr == 0) then
          call wryte(noflx,lg,g)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  74) ',
     &     'Liquid soil moist content (frac) layer 10cm and 0cm'
        endif
!..........................................................
!  - liquid soil moist content layer 40cm and 10cm

        ngridss = ngridss+1
        call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

        lbm = (slmskful == 1._kind_io8)
        call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,130,icen,
     &             igen,1,islc,i2dbls,10,40,iyr,imo,ida,ihr,
     &             ifhour,ithr,0,inst,0,0,icen2,ids(isoilm),iens,
     &             0.,0.,0.,0.,0.,0.,g,lg,ierr)

        if (ierr == 0) then
          call wryte(noflx,lg,g)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  75) ',
     &     'Liquid soil moist content (frac) layer 40cm and 10cm'
        endif
!..........................................................
!  - liquid soil moist content layer 100cm and 40cm
!
        if (lsoil > 2) then

          ngridss = ngridss+1
          call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

          lbm = (slmskful == 1._kind_io8)
          call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,130,icen,
     &               igen,1,islc,i2dbls,40,100,iyr,imo,ida,ihr,
     &               ifhour,ithr,0,inst,0,0,icen2,ids(isoilm),iens,
     &               0.,0.,0.,0.,0.,0.,g,lg,ierr)

          if (ierr == 0) then
            call wryte(noflx,lg,g)
          else
            print*,'wrtsfc gribit ierr=',ierr,'  76) ',
     &       'Liquid soil moist content (frac) layer 100cm and 40cm'
          endif
!..........................................................
!  - liquid soil moist content layer 200cm and 100cm

          ngridss = ngridss+1
          call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

          lbm = (slmskful == 1._kind_io8)
          call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,130,icen,
     &               igen,1,islc,i2dbls,100,200,iyr,imo,ida,ihr,
     &               ifhour,ithr,0,inst,0,0,icen2,ids(isoilm),iens,
     &               0.,0.,0.,0.,0.,0.,g,lg,ierr)

          if (ierr == 0) then
            call wryte(noflx,lg,g)
          else
            print*,'wrtsfc gribit ierr=',ierr,'  77) ',
     &       'Liquid soil moist content (frac) layer 200cm and 100cm'
          endif
        endif
!..........................................................
!  - snow depth

        ngridss = ngridss+1
        call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

!       lbm = (slmskful == 1._kind_io8)
        lbm = (slmskful == 1._kind_io8 .or. slmskful == 2._kind_io8)
        call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,igen,
     &              1,isnod,isfc,0,0,iyr,imo,ida,ihr,
     &              ifhour,ithr,0,inst,0,0,icen2,ids(isnod),iens,
     &              0.,0.,0.,0.,0.,0.,g,lg,ierr)

        if (ierr == 0) then
          call wryte(noflx,lg,g)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  78) ',
     &     'Snow depth (m) land surface'
        endif
!..........................................................
!  - canopy water content
!
        ngridss = ngridss+1
        call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

        lbm = (slmskful == 1._kind_io8)
        call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,igen,
     &              1,icnp,isfc,0,0,iyr,imo,ida,ihr,
     &              ifhour,ithr,0,inst,0,0,icen2,ids(icnp),iens,
     &              0.,0.,0.,0.,0.,0.,g,lg,ierr)

        if (ierr == 0) then
          call wryte(noflx,lg,g)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  79) ',
     &     'Canopy water content (kg/m^2) land surface'
        endif
!..........................................................
!  - the following 30 records are for land mdl use
!  - surface roughness
!
        ngridss = ngridss+1
        call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

        call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,igen,
     &              0,izorl,isfc,0,0,iyr,imo,ida,ihr,
     &              ifhour,ithr,0,inst,0,0,icen2,ids(izorl),iens,
     &              0.,0.,0.,0.,0.,0.,g,lg,ierr)

        if (ierr == 0) then
          call wryte(noflx,lg,g)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  80) ',
     &     'Surface roughness (m)'
        endif
!..........................................................
!  - vegetation fraction
!
        ngridss = ngridss+1
        call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

        lbm = (slmskful == 1._kind_io8)
        call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,igen,
     &              1,iveg,isfc,0,0,iyr,imo,ida,ihr,
     &              ifhour,ithr,0,inst,0,0,icen2,ids(iveg),iens,
     &              0.,0.,0.,0.,0.,0.,g,lg,ierr)

        if (ierr == 0) then
          call wryte(noflx,lg,g)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  81) ',
     &     'Vegetation fraction (fractional) land surface'
        endif
!..........................................................
!  - vegetation type
!
        ngridss = ngridss+1
        call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

        lbm = (slmskful == 1._kind_io8)
        call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,igen,
     &              1,ivtp,isfc,0,0,iyr,imo,ida,ihr,
     &              ifhour,ithr,0,inst,0,0,icen2,ids(ivtp),iens,
     &              0.,0.,0.,0.,0.,0.,g,lg,ierr)

        if (ierr == 0) then
          call wryte(noflx,lg,g)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  82) ',
     &     'Vegetation type land surface'
        endif
!..........................................................
!  - soil type
!
        ngridss = ngridss+1
        call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

        lbm = (slmskful == 1._kind_io8)
        call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,igen,
     &              1,istp,isfc,0,0,iyr,imo,ida,ihr,
     &              ifhour,ithr,0,inst,0,0,icen2,ids(istp),iens,
     &              0.,0.,0.,0.,0.,0.,g,lg,ierr)

        if (ierr == 0) then
          call wryte(noflx,lg,g)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  83) ',
     &     'Soil type land surface'
        endif
!..........................................................
!  - slope type
!
        ngridss = ngridss+1
        call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

        lbm = (slmskful == 1._kind_io8)
        call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,130,icen,
     &              igen,1,islo,isfc,0,0,iyr,imo,ida,ihr,
     &              ifhour,ithr,0,inst,0,0,icen2,ids(islo),iens,
     &              0.,0.,0.,0.,0.,0.,g,lg,ierr)

        if (ierr == 0) then
          call wryte(noflx,lg,g)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  84) ',
     &     'Slope type land surface'
        endif
!..........................................................
!  - frictional velocity
!
        ngridss = ngridss+1
        call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

        call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,igen,
     &              0,iust,isfc,0,0,iyr,imo,ida,ihr,
     &              ifhour,ithr,0,inst,0,0,icen2,ids(iust),iens,
     &              0.,0.,0.,0.,0.,0.,g,lg,ierr)

        if (ierr == 0) then
          call wryte(noflx,lg,g)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  85) ',
     &     'Frictional velocity (m/s)'
        endif
!..........................................................
!  - surface height
!
        ngridss = ngridss+1
        call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

        call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,igen,
     &              0,ihgt,isfc,0,0,iyr,imo,ida,ihr,
     &              ifhour,ithr,0,inst,0,0,icen2,ids(ihgt),iens,
     &              0.,0.,0.,0.,0.,0.,g,lg,ierr)

        if (ierr == 0) then
          call wryte(noflx,lg,g)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  86) ',
     &     'Surface height (m)'
        endif
!..........................................................
!  - freezing precip flag
!
        ngridss = ngridss+1
        call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

        lbm = (slmskful == 1._kind_io8)
        call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,igen,
     &              1,irst,isfc,0,0,iyr,imo,ida,ihr,
     &              ifhour,ithr,0,inst,0,0,icen2,ids(irst),iens,
     &              0.,0.,0.,0.,0.,0.,g,lg,ierr)

        if (ierr == 0) then
          call wryte(noflx,lg,g)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  87) ',
     &     'Freezing precip flag land surface'
        endif
!..........................................................
!  - exchange coefficient ch
!
        ngridss = ngridss+1
        call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

        call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,igen,
     &              0,ichh,isfc,0,0,iyr,imo,ida,ihr,
     &              ifhour,ithr,0,inst,0,0,icen2,ids(ichh),iens,
     &              0.,0.,0.,0.,0.,0.,g,lg,ierr)

        if (ierr == 0) then
          call wryte(noflx,lg,g)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  88) ',
     &     'Exchange coefficient CH(m/s)'
        endif
!..........................................................
!  - exchange coefficient cm
!
        ngridss = ngridss+1
        call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

        call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,130,icen,
     &              igen,0,icmm,isfc,0,0,iyr,imo,ida,ihr,
     &              ifhour,ithr,0,inst,0,0,icen2,ids(icmm),iens,
     &              0.,0.,0.,0.,0.,0.,g,lg,ierr)

        if (ierr == 0) then
          call wryte(noflx,lg,g)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  89) ',
     &     'Exchange coefficient CM(m/s)'
        endif
!..........................................................
!  - potential evaporation rate
!
        ngridss = ngridss+1
        call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

        lbm = (slmskful == 1._kind_io8)
        call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,igen,
     &              1,iep,isfc,0,0,iyr,imo,ida,ihr,
     &              ifhour,ithr,0,inst,0,0,icen2,ids(iep),iens,
     &              0.,0.,0.,0.,0.,0.,g,lg,ierr)

        if (ierr == 0) then
          call wryte(noflx,lg,g)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  90) ',
     &     'Potential evaporation rate (w/m**2) land and sea surface'
        endif

        if (.not. climate) then     ! do not output those fields in climate mode
!..........................................................
!  - downward long wave radiation flux (instantaneous value)
!
          ngridss = ngridss+1
          call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

          call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,
     &                igen,0,idlwf,isfc,0,0,iyr,imo,ida,ihr,
     &                ifhour,ithr,0,inst,0,0,icen2,ids(idlwf),iens,
     &                0.,0.,0.,0.,0.,0.,g,lg,ierr)

          if (ierr == 0) then
            call wryte(noflx,lg,g)
          else
            print*,'wrtsfc gribit ierr=',ierr,'  91) ',
     &       'Downward long wave radiation flux (w/m**2)'
          endif
!..........................................................
!  - upward long wave radiation flux (instantaneous value)
!
          ngridss = ngridss+1
          call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

          call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,
     &                igen,0,iulwf,isfc,0,0,iyr,imo,ida,ihr,
     &                ifhour,ithr,0,inst,0,0,icen2,ids(iulwf),iens,
     &                0.,0.,0.,0.,0.,0.,g,lg,ierr)

          if (ierr == 0) then
            call wryte(noflx,lg,g)
          else
            print*,'wrtsfc gribit ierr=',ierr,'  92) ',
     &       'Upward long wave radiation flux (w/m**2)'
          endif
!..........................................................
!  - upward short wave radiation flux (instantaneous value)
!
          ngridss = ngridss+1
          call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

          call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,
     &                igen,0,iuswf,isfc,0,0,iyr,imo,ida,ihr,
     &                ifhour,ithr,0,inst,0,0,icen2,ids(iuswf),iens,
     &                0.,0.,0.,0.,0.,0.,g,lg,ierr)

          if (ierr == 0) then
            call wryte(noflx,lg,g)
          else
            print*,'wrtsfc gribit ierr=',ierr,'  93) ',
     &       'Upward short wave radiation flux (w/m**2)'
          endif
!..........................................................
!  - downward short wave radiation flux (instantaneous value)
!
          ngridss = ngridss+1
          call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

          call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,
     &                igen,0,idswf,isfc,0,0,iyr,imo,ida,ihr,
     &                ifhour,ithr,0,inst,0,0,icen2,ids(idswf),iens,
     &                0.,0.,0.,0.,0.,0.,g,lg,ierr)

          if (ierr == 0) then
            call wryte(noflx,lg,g)
          else
            print*,'wrtsfc gribit ierr=',ierr,'  94) ',
     &       'Downward short wave radiation flux (w/m**2)'
          endif
!..........................................................
!  - sensible heat flux (instantaneous value)
!
          ngridss = ngridss+1
          call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

          call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,
     &                igen,0,ishflx,isfc,0,0,iyr,imo,ida,ihr,
     &                ifhour,ithr,0,inst,0,0,icen2,ids(ishflx),iens,
     &                0.,0.,0.,0.,0.,0.,g,lg,ierr)

          if (ierr == 0) then
            call wryte(noflx,lg,g)
          else
            print*,'wrtsfc gribit ierr=',ierr,'  95) ',
     &       'Sensible heat flux (w/m**2) land and sea surface'
          endif
!..........................................................
!  - latent heat flux (instantaneous value)
!
          ngridss = ngridss+1
          call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

          call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,
     &                igen,0,ilhflx,isfc,0,0,iyr,imo,ida,ihr,
     &                ifhour,ithr,0,inst,0,0,icen2,ids(ilhflx),iens,
     &                0.,0.,0.,0.,0.,0.,g,lg,ierr)

          if (ierr == 0) then
            call wryte(noflx,lg,g)
          else
            print*,'wrtsfc gribit ierr=',ierr,'  96) ',
     &       'Latent heat flux (w/m**2) land and sea surface'
          endif
!..........................................................
!  - ground heat flux (instantaneous value)
!
          ngridss = ngridss+1
          call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

          lbm = (slmskful /= 0._kind_io8)
          call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,
     &                igen,1,ighflx,isfc,0,0,iyr,imo,ida,ihr,
     &                ifhour,ithr,0,inst,0,0,icen2,ids(ighflx),iens,
     &                0.,0.,0.,0.,0.,0.,g,lg,ierr)

          if (ierr == 0) then
            call wryte(noflx,lg,g)
          else
            print*,'wrtsfc gribit ierr=',ierr,'  97) ',
     &       'Ground heat flux (w/m**2) land and sea surface'
          endif
        endif             ! end of if(.not climate)
!..........................................................
!  - surface runoff
!
        ngridss = ngridss+1
        call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

        lbm = (slmskful == 1._kind_io8)
        call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,igen,
     &              1,isrf,isfc,0,0,iyr,imo,ida,ihr,
     &              ifhour,ifhr,ithr,iacc,0,0,icen2,ids(isrf),iens,
     &              0.,0.,0.,0.,0.,0.,g,lg,ierr)

        if (ierr == 0) then
          call wryte(noflx,lg,g)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  98) ',
     &     'Surface runoff (kg/m^2) land surface'
        endif
!..........................................................
!  - lowest model level temp
!
        ngridss = ngridss+1
        call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

        call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,igen,
     &              0,itemp,isglev,1,1,iyr,imo,ida,ihr,
     &              ifhour,ithr,0,inst,0,0,icen2,ids(itemp),iens,
     &              0.,0.,0.,0.,0.,0.,g,lg,ierr)

        if (ierr == 0) then
          call wryte(noflx,lg,g)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  99) ',
     &     ' Lowest model level temp (k)'
        endif
!..........................................................
!  - lowest model specific humidity
!
        ngridss = ngridss+1
        call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

        call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,igen,
     &              0,isphum,isglev,1,1,iyr,imo,ida,ihr,
     &              ifhour,ithr,0,inst,0,0,icen2,ids(isphum),iens,
     &              0.,0.,0.,0.,0.,0.,g,lg,ierr)

        if (ierr == 0) then
          call wryte(noflx,lg,g)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  100) ',
     &     'Lowest model specific humidity (kg/kg)'
        endif
!..........................................................
!  - lowest model u wind
!
        ngridss = ngridss+1
        call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

        call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,igen,
     &              0,iznlw,isglev,1,1,iyr,imo,ida,ihr,
     &              ifhour,ithr,0,inst,0,0,icen2,ids(iznlw),iens,
     &              0.,0.,0.,0.,0.,0.,g,lg,ierr)

        if (ierr == 0) then
          call wryte(noflx,lg,g)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  101) ',
     &     'Lowest model u wind (m/s)'
        endif
!..........................................................
!  - lowest model v wind
!
        ngridss = ngridss+1
        call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

        call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,igen,
     &              0,imerw,isglev,1,1,iyr,imo,ida,ihr,
     &              ifhour,ithr,0,inst,0,0,icen2,ids(imerw),iens,
     &              0.,0.,0.,0.,0.,0.,g,lg,ierr)

        if (ierr == 0) then
          call wryte(noflx,lg,g)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  102) ',
     &     'Lowest model v wind (m/s)'
        endif
!..........................................................
!  - lowest model level height
!
        ngridss = ngridss+1
        call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

        lbm = (slmskful == 1._kind_io8)
        call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,igen,
     &              1,ihgt,isglev,1,1,iyr,imo,ida,ihr,
     &              ifhour,ithr,0,inst,0,0,icen2,ids(ihgt),iens,
     &              0.,0.,0.,0.,0.,0.,g,lg,ierr)

        if (ierr == 0) then
          call wryte(noflx,lg,g)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  103) ',
     &     'Lowest model level height (m) land surface'
        endif
!..........................................................
!  - direct evaporation from bare soil
!
        ngridss = ngridss+1
        call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

        lbm = (slmskful == 1._kind_io8)
        call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,igen,
     &              1,ievbs,isfc,0,0,iyr,imo,ida,ihr,
     &              ifhour,ifhr,ithr,iavg,0,0,icen2,ids(ievbs),iens,
     &              0.,0.,0.,0.,0.,0.,g,lg,ierr)

        if (ierr == 0) then
          call wryte(noflx,lg,g)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  104) ',
     &     'Direct evaporation from bare soil(w/m^2) land surface'
        endif
!..........................................................
!  - canopy water evaporation
!
        ngridss = ngridss+1
        call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

        lbm = (slmskful == 1._kind_io8)
        call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,igen,
     &              1,ievcw,isfc,0,0,iyr,imo,ida,ihr,
     &              ifhour,ifhr,ithr,iavg,0,0,icen2,ids(ievbs),iens,
     &              0.,0.,0.,0.,0.,0.,g,lg,ierr)

        if (ierr == 0) then
          call wryte(noflx,lg,g)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  105) ',
     &     'Canopy water evaporation(w/m^2) land surface'
        endif
!..........................................................
!  - transpiration
!
        ngridss = ngridss+1
        call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

        lbm = (slmskful == 1._kind_io8)
        call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,igen,
     &              1,itran,isfc,0,0,iyr,imo,ida,ihr,
     &              ifhour,ifhr,ithr,iavg,0,0,icen2,ids(itran),iens,
     &              0.,0.,0.,0.,0.,0.,g,lg,ierr)

        if (ierr == 0) then
          call wryte(noflx,lg,g)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  106) ',
     &     'Transpiration (w/m^2) land surface'
        endif
!..........................................................
!  - snow sublimation
!
        ngridss = ngridss+1
        call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

        lbm = (slmskful == 1._kind_io8)
        call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,130,icen,
     &              igen,1,isbs,isfc,0,0,iyr,imo,ida,ihr,
     &              ifhour,ifhr,ithr,iavg,0,0,icen2,ids(isbs),iens,
     &              0.,0.,0.,0.,0.,0.,g,lg,ierr)

        if (ierr == 0) then
          call wryte(noflx,lg,g)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  107) ',
     &     'Snow sublimation (w/m^2) land surface'
        endif
!..........................................................
!  - snow cover
!
        ngridss = ngridss+1
        call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

        lbm = (slmskful == 1._kind_io8)
        call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,igen,
     &              1,isnc,isfc,0,0,iyr,imo,ida,ihr,
     &              ifhour,ifhr,ithr,iavg,0,0,icen2,ids(isnc),iens,
     &              0.,0.,0.,0.,0.,0.,g,lg,ierr)

        if (ierr == 0) then
          call wryte(noflx,lg,g)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  108) ',
     &     'Snow cover (fraction) land surface'
        endif
!..........................................................
!  - total column soil moisture
!
        ngridss = ngridss+1
        call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

        lbm = (slmskful == 1._kind_io8)
        call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,igen,
     &              1,istc,i2dbls,0,200,iyr,imo,ida,ihr,
     &              ifhour,ithr,0,inst,0,0,icen2,ids(istc),iens,
     &              0.,0.,0.,0.,0.,0.,g,lg,ierr)

        if (ierr == 0) then
          call wryte(noflx,lg,g)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  109) ',
     &     'Total column soil moisture (kg/m^2) land surface'
        endif
!..........................................................
!  - snow phase-change heat flux
!
        ngridss = ngridss+1
        call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

        lbm = (slmskful == 1._kind_io8)
        call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,130,icen,
     &              igen,1,isnohf,isfc,0,0,iyr,imo,ida,ihr,
     &              ifhour,ifhr,ithr,iavg,0,0,icen2,ids(isnohf),iens,
     &              0.,0.,0.,0.,0.,0.,g,lg,ierr)

        if (ierr == 0) then
          call wryte(noflx,lg,g)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  110) ',
     &     'Snow phase-change heat flux [w/m^2] land surface'
        endif
!..........................................................
!  - wilting point
!
        ngridss = ngridss+1
        call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

        lbm = (slmskful == 1._kind_io8)
        call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,130,icen,
     &              igen,1,ismcwlt,isfc,0,0,iyr,imo,ida,ihr,
     &              ifhour,ithr,0,inst,0,0,icen2,ids(ismcwlt),iens,
     &              0.,0.,0.,0.,0.,0.,g,lg,ierr)

        if (ierr == 0) then
          call wryte(noflx,lg,g)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  111) ',
     &     'Wilting point [fraction] land surface'
        endif
!..........................................................
!  - field capacity
!
        ngridss = ngridss+1
        call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

        lbm = (slmskful == 1._kind_io8)
        call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,130,icen,
     &              igen,1,ismcref,isfc,0,0,iyr,imo,ida,ihr,
     &              ifhour,ithr,0,inst,0,0,icen2,ids(ismcref),iens,
     &              0.,0.,0.,0.,0.,0.,g,lg,ierr)

        if (ierr == 0) then
          call wryte(noflx,lg,g)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  112) ',
     &     'Field capacity [fraction] land surface'
        endif
!..........................................................
!  - accumulated sunshine duration time
!
        ngridss = ngridss+1
        call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

        call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,133,icen,
     &              igen,0,isuntm,isfc,0,0,iyr,imo,ida,ihr,
     &              ifhour,ifhr,ithr,iacc,0,0,icen2,ids(isuntm),iens,
     &              0.,0.,0.,0.,0.,0.,g,lg,ierr)

        if (ierr == 0) then
          call wryte(noflx,lg,g)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  113) ',
     &     'Accumulated sunshine duration time (sec)'
        endif
!
!..........................................................
!  - frozen precipitation fraction
!
        ngridss = ngridss+1
        call unsplit2z(ngridss,ngrids_flx,wrkga,global_lats_r)

        call gribit_gsm(wrkga,lbm,4,lonr,latr,16,cl1,ilpds,2,icen,igen,
     &              0,isr,isfc,0,0,iyr,imo,ida,ihr,
     &              ifhour,ithr,0,inst,0,0,icen2,1,iens,
     &              0.,0.,0.,0.,0.,0.,g,lg,ierr)

        if (ierr == 0) then
          call wryte(noflx,lg,g)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  114) ',
     &     'Frozen precipitation (fraction)'
        endif

!       write(0,*)' at the end in wrtflx_w max_ngridss=',ngridss
!
        PRINT *,'GRIB FLUX FILE WRITTEN ',FHOUR,IDATE,noflx

      endif           ! if (me == ioproc) then
!
      RETURN

      end subroutine wrtflx_w

!
!-------------------------------------------------------------------------
!
!      INTEGER FUNCTION nfill(C)
!      implicit none
!      integer j
!      CHARACTER*(*) C
!      NFILL=LEN(C)
!      DO J=1,NFILL
!        IF(C(J:J).EQ.' ') THEN
!          NFILL=J-1
!          RETURN
!        ENDIF
!      ENDDO
!      RETURN
!      END
 
 
      SUBROUTINE nst_collect (nst_fld,global_lats_r,lonsperlar)
!!
      use resol_def,               ONLY: latr, lonr,ngrids_nst
      use mod_state,               ONLY: buff_mult_piecenst,ngridnst
      use layout1,                 ONLY: lats_node_r,lats_node_r_max
      use gfs_physics_nst_var_mod, ONLY: Nst_Var_Data
      USE machine,                 ONLY: kind_io8, kind_io4
      implicit none
!!
      TYPE(Nst_Var_Data)        :: nst_fld
!
      INTEGER             GLOBAL_LATS_R(latr)
      INTEGER             lonsperlar(latr)
!!
      real(kind=kind_io8) buffo(lonr,lats_node_r)
      integer             kmsk(lonr,lats_node_r_max)
!!
!
      if(.not. allocated(buff_mult_piecenst)) then
       allocate(buff_mult_piecenst(lonr,lats_node_r_max,1:ngrids_nst+1))
      endif
!
      kmsk = nint(nst_fld%slmsk)
!
!-- slmsk
      ngridnst = 1
      CALL uninterprez(1,kmsk,buffo,nst_fld%slmsk,
     &       global_lats_r,lonsperlar,buff_mult_piecenst(1,1,ngridnst))
!
!-- xt
      ngridnst = ngridnst+1
      CALL uninterprez(1,kmsk,buffo,nst_fld%xt,
     &       global_lats_r,lonsperlar,buff_mult_piecenst(1,1,ngridnst))
!
!-- xs
      ngridnst = ngridnst+1
      CALL uninterprez(1,kmsk,buffo,nst_fld%xs,
     &       global_lats_r,lonsperlar,buff_mult_piecenst(1,1,ngridnst))
!
!-- xu
      ngridnst = ngridnst+1
      CALL uninterprez(1,kmsk,buffo,nst_fld%xu,
     &       global_lats_r,lonsperlar,buff_mult_piecenst(1,1,ngridnst))
!
!-- xv
      ngridnst = ngridnst+1
      CALL uninterprez(1,kmsk,buffo,nst_fld%xv,
     &       global_lats_r,lonsperlar,buff_mult_piecenst(1,1,ngridnst))
!
!-- 6 xz
      ngridnst = ngridnst+1
      CALL uninterprez(1,kmsk,buffo,nst_fld%xz,
     &       global_lats_r,lonsperlar,buff_mult_piecenst(1,1,ngridnst))
!
!-- zm
      ngridnst = ngridnst+1
      CALL uninterprez(1,kmsk,buffo,nst_fld%zm,
     &       global_lats_r,lonsperlar,buff_mult_piecenst(1,1,ngridnst))
!
!-- xtts
      ngridnst = ngridnst+1
      CALL uninterprez(1,kmsk,buffo,nst_fld%xtts,
     &       global_lats_r,lonsperlar,buff_mult_piecenst(1,1,ngridnst))
!
!-- xzts
      ngridnst = ngridnst+1
      CALL uninterprez(1,kmsk,buffo,nst_fld%xzts,
     &       global_lats_r,lonsperlar,buff_mult_piecenst(1,1,ngridnst))
!
!-- 10 dt_cool
      ngridnst = ngridnst+1
      CALL uninterprez(1,kmsk,buffo,nst_fld%dt_cool,
     &       global_lats_r,lonsperlar,buff_mult_piecenst(1,1,ngridnst))
!
!-- z_c
      ngridnst = ngridnst+1
      CALL uninterprez(1,kmsk,buffo,nst_fld%z_c,
     &       global_lats_r,lonsperlar,buff_mult_piecenst(1,1,ngridnst))
!
!-- c_0
      ngridnst = ngridnst+1
      CALL uninterprez(1,kmsk,buffo,nst_fld%c_0,
     &       global_lats_r,lonsperlar,buff_mult_piecenst(1,1,ngridnst))
!
!-- c_d
      ngridnst = ngridnst+1
      CALL uninterprez(1,kmsk,buffo,nst_fld%c_d,
     &       global_lats_r,lonsperlar,buff_mult_piecenst(1,1,ngridnst))
!
!-- w_0
      ngridnst = ngridnst+1
      CALL uninterprez(1,kmsk,buffo,nst_fld%w_0,
     &       global_lats_r,lonsperlar,buff_mult_piecenst(1,1,ngridnst))
!
!-- w_d
      ngridnst = ngridnst+1
      CALL uninterprez(1,kmsk,buffo,nst_fld%w_d,
     &       global_lats_r,lonsperlar,buff_mult_piecenst(1,1,ngridnst))
!
!-- d_conv
      ngridnst = ngridnst+1
      CALL uninterprez(1,kmsk,buffo,nst_fld%d_conv,
     &       global_lats_r,lonsperlar,buff_mult_piecenst(1,1,ngridnst))
!
!-- ifd
      ngridnst = ngridnst+1
      CALL uninterprez(1,kmsk,buffo,nst_fld%ifd,
     &       global_lats_r,lonsperlar,buff_mult_piecenst(1,1,ngridnst))
!
!-- tref
      ngridnst = ngridnst+1
      CALL uninterprez(1,kmsk,buffo,nst_fld%tref,
     &       global_lats_r,lonsperlar,buff_mult_piecenst(1,1,ngridnst))
!
!-- qrain
      ngridnst = ngridnst+1
      CALL uninterprez(1,kmsk,buffo,nst_fld%qrain,
     &       global_lats_r,lonsperlar,buff_mult_piecenst(1,1,ngridnst))
!

      return
      end subroutine nst_collect
