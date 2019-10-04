!-----------------------------------------------------------------------
                        module module_INIT_READ_NEMSIO
!-----------------------------------------------------------------------
use mpi
use esmf
use module_kinds
use module_dm_parallel,only : ids,ide,jds,jde &
                             ,ims,ime,jms,jme &
                             ,its,ite,jts,jte &
                             ,its_h2,ite_h2,jts_h2,jte_h2 &
                             ,lm &
                             ,mype_share,npes,num_pts_max &
                             ,mpi_comm_comp &
                             ,dstrb
use module_exchange
use module_constants
use module_solver_internal_state,only: solver_internal_state
use nemsio_module_mpi
!
!-----------------------------------------------------------------------
!
      implicit none
!
      private
!
      public read_nemsio
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      contains
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
                        subroutine read_nemsio &
      (INT_STATE,my_domain_id,rc)
!
!-----------------------------------------------------------------------
!
implicit none
!
!------------------------
!***  Argument variables
!------------------------
!
type(solver_internal_state),pointer :: int_state

integer(kind=kint),intent(in) :: &
 my_domain_id

integer(kind=kint),intent(out) :: &
 rc
!
!---------------------
!***  Local variables
!---------------------
!
integer(kind=kint) :: &
 i &                         ! index in x direction
,iend &
,ierr &
,irtn &
,j &                         ! index in y direction
,jend &
,k &                         ! index
,kount &
,ks &                        ! tracer index
,l &                         ! index in p direction
,length &
,n &
,nl &
,nv

integer(kind=kint) :: &
 iyear_fcst           &
,imonth_fcst          &
,iday_fcst            &
,ihour_fcst
 
integer(kind=kint) :: &
 lb &
,ub &
,nvars_bc_2d_h &
,nvars_bc_3d_h &
,nvars_bc_4d_h &
,nvars_bc_2d_v &
,nvars_bc_3d_v 

integer(kind=kint) :: idate(7)
integer(kind=kint) :: fcstdate(7)
character(3)       ::tn
 
real(kind=kfpt):: &
 tend,tend_max,arg

real(kind=kfpt),dimension(:),allocatable :: &
 all_bc_data &
,tmp, psint

logical(kind=klog) :: opened

type(nemsio_gfile) :: gfile

integer(kind=kint) :: &
 mype

character(64):: &
 infile

integer(kind=kint):: &
 ihrend &                    ! maximum forecast length, hours
,ntsd &
,ntstm_max &
,nfcst

integer nrec,mysize,myrank
integer fldsize,fldst,js,recn
character(16),allocatable :: recname(:), reclevtyp(:)
integer,allocatable       :: reclev(:)
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      mype=mype_share
!
!-----------------------------------------------------------------------
!***  Initialize nemsio_module
!-----------------------------------------------------------------------
!
      call nemsio_init()
!
!-----------------------------------------------------------------------
!
      read_blocks: if(.not.int_state%RESTART) then                      ! cold start
!
!-----------------------------------------------------------------------
!
        write(infile,'(a,i2.2,a)')'input_domain_',my_domain_id,'_nemsio'
        call nemsio_open(gfile,infile,'read',mpi_comm_comp,iret=ierr)
        if(ierr/=0)then
          write(0,*)' Unable to open file ',trim(infile),' in READ_NEMSIO'
          rc = ierr
          return
        endif
!
        call nemsio_getfilehead(gfile,nrec=nrec,iret=ierr)
!
        allocate(recname(nrec),reclevtyp(nrec),reclev(nrec))
        call nemsio_getfilehead(gfile,recname=recname,reclevtyp=reclevtyp,   &
                                reclev=reclev,iret=ierr)
!
!-----------------------------------------------------------------------
!***  Get run,idat,ihrst,ihrend,ntsd
!-----------------------------------------------------------------------
!
        call nemsio_getheadvar(gfile,'RUN',int_state%RUN,ierr)
        call nemsio_getheadvar(gfile,'IDAT',int_state%IDAT,ierr)
        call nemsio_getheadvar(gfile,'IHRST',int_state%IHRST,ierr)
        call nemsio_getheadvar(gfile,'IHREND',int_state%IHREND,ierr)
        call nemsio_getheadvar(gfile,'NTSD',ntsd,ierr)
!
!-----------------------------------------------------------------------
!***  Print the time information.
!-----------------------------------------------------------------------
!
        if(mype==0)then
          write(0,*) 'run, idat,ntsd: ', int_state%RUN, int_state%IDAT, NTSD
          write(0,*)' Start year =',int_state%IDAT(3)
          write(0,*)' Start month=',int_state%IDAT(2)
          write(0,*)' Start day  =',int_state%IDAT(1)
          write(0,*)' Start hour =',int_state%IHRST
          write(0,*)' Timestep   =',int_state%DT
          write(0,*)' Steps/hour =',3600./int_state%DT
          if(.not.int_state%GLOBAL)write(0,*)' Max fcst hours=',int_state%IHREND
        endif
!
!-----------------------------------------------------------------------
!***  Get SW corner of nest domains on their parent grid.
!-----------------------------------------------------------------------
!
        call nemsio_getheadvar(gfile,'IPARSTRT',int_state%I_PAR_STA,ierr)
        call nemsio_getheadvar(gfile,'JPARSTRT',int_state%J_PAR_STA,ierr)
!
!-----------------------------------------------------------------------
!
        call nemsio_getheadvar(gfile,'PT',int_state%PT,ierr)
        call nemsio_getheadvar(gfile,'PDTOP',int_state%PDTOP,ierr)
        call nemsio_getheadvar(gfile,'LPT2',int_state%LPT2,ierr)
        call nemsio_getheadvar(gfile,'SGM',int_state%SGM,ierr)
        call nemsio_getheadvar(gfile,'SG1',int_state%SG1,ierr)
        call nemsio_getheadvar(gfile,'DSG1',int_state%DSG1,ierr)
        call nemsio_getheadvar(gfile,'SGML1',int_state%SGML1,ierr)
        call nemsio_getheadvar(gfile,'SG2',int_state%SG2,ierr)
        call nemsio_getheadvar(gfile,'DSG2',int_state%DSG2,ierr)
        call nemsio_getheadvar(gfile,'SGML2',int_state%SGML2,ierr)
!
        fldsize=(jte-jts+1)*(ite-its+1)
        allocate(tmp((ite-its+1)*(jte-jts+1)*nrec),stat=i)
!
        call nemsio_denseread(gfile,its,ite,jts,jte,tmp,iret=ierr)
        if(ierr/=0) then
          write(0,*)'WRONG: Could not read all the fields in the file!'
          rc = ierr
          return
        endif
!-----------------------------------------------------------------------
!
!-- fis
        int_state%fis=0.
        call getrecn(recname,reclevtyp,reclev,nrec,'fis','sfc',1,recn)
        if(recn/=0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              int_state%fis(i,j)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
        call halo_exch(int_state%fis,1,3,3)
!-----------------------------------------------------------------------
!
!-- stdh
        int_state%stdh=0.
        call getrecn(recname,reclevtyp,reclev,nrec,'stdh','sfc',1,recn)
        if(recn>0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              int_state%stdh(i,j)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
        call halo_exch(int_state%stdh,1,3,3)
!-----------------------------------------------------------------------
!
!-- sm
        int_state%sm=0.
        call getrecn(recname,reclevtyp,reclev,nrec,'sm','sfc',1,recn)
        if(recn/=0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              int_state%sm(i,j)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
        call halo_exch(int_state%sm,1,2,2)
!-----------------------------------------------------------------------
!
!-- dpres
        int_state%pd=0.
        call getrecn(recname,reclevtyp,reclev,nrec,'dpres','hybrid sig lev',1,recn)
        if(recn/=0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              int_state%pd(i,j)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
        call halo_exch(int_state%pd,1,2,2)
!-----------------------------------------------------------------------
!
!-- ugrd
      int_state%u=0.
      do l=1,lm
        call getrecn(recname,reclevtyp,reclev,nrec,'ugrd','mid layer',l,recn)
        if(recn/=0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              int_state%u(i,j,l)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
      enddo
      call halo_exch(int_state%u,lm,2,2)
!-----------------------------------------------------------------------
!
!--vgrd
        int_state%v=0.
        do l=1,lm
         call getrecn(recname,reclevtyp,reclev,nrec,'vgrd','mid layer',l,recn)
         if(recn/=0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              int_state%v(i,j,l)=tmp(i-its+1+js+fldst)
            enddo
          enddo
         endif
        enddo
        call halo_exch(int_state%v,lm,2,2)
!-----------------------------------------------------------------------
!--tmp
        int_state%t=0.
        do l=1,lm
         call getrecn(recname,reclevtyp,reclev,nrec,'tmp','mid layer',l,recn)
         if(recn/=0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              int_state%t(i,j,l)=tmp(i-its+1+js+fldst)
            enddo
          enddo
         endif
        enddo
        call halo_exch(int_state%t,lm,2,2)
!-----------------------------------------------------------------------
!
!--spfh
        int_state%q=0.
        do l=1,lm
         call getrecn(recname,reclevtyp,reclev,nrec,'spfh','mid layer',l,recn)
         if(recn/=0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              int_state%q(i,j,l)=tmp(i-its+1+js+fldst)
            enddo
          enddo
         endif
        enddo
        call halo_exch(int_state%q,lm,2,2)
!-----------------------------------------------------------------------
!
!-- clwmr
        int_state%cw=0.
        do l=1,lm
         call getrecn(recname,reclevtyp,reclev,nrec,'clwmr','mid layer',l,recn)
         if(recn/=0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              int_state%cw(i,j,l)=tmp(i-its+1+js+fldst)
            enddo
          enddo
         endif
        enddo
        call halo_exch(int_state%cw,lm,2,2)
!
!-----------------------------------------------------------------------
!
        tend_max=real(int_state%IHREND)
        ntstm_max=nint(tend_max*3600./int_state%DT)+1
        tend=real(int_state%nhours_fcst)
        int_state%NTSTM=nint(tend*3600./int_state%DT)+1
        if(.not.int_state%global)then
          if(mype==0)then
            write(0,*)' Max runtime is ',tend_max,' hours'
          endif
        endif
        if(mype==0)then
          write(0,*)' Requested runtime is ',tend,' hours'
          write(0,*)' NTSTM=',int_state%NTSTM
        endif
        if(int_state%NTSTM>ntstm_max.and..not.int_state%global)then
          if(mype==0)then
            write(0,*)' Requested fcst length exceeds maximum'
            write(0,*)' Resetting to maximum'
          endif
          int_state%NTSTM=min(int_state%NTSTM,ntstm_max)
        endif
!
        int_state%ihr=nint(ntsd*int_state%DT/3600.)
!
!-----------------------------------------------------------------------
        do l=1,lm
          int_state%pdsg1(l)=int_state%dsg1(l)*int_state%pdtop
          int_state%psgml1(l)=int_state%sgml1(l)*int_state%pdtop+int_state%pt
        enddo
!
        do l=1,lm+1
          int_state%psg1(l)=int_state%sg1(l)*int_state%pdtop+int_state%pt
        enddo
!
!-- Compute pressure dependent floor values of EPSL & EPSQ2
!
        ALLOCATE(PSINT(1:LM+1))
        DO K=1,LM+1
          PSINT(K)=int_state%SG2(k)*101325.+int_state%PSG1(K)
        ENDDO
  !
        DO K=2,LM
          ARG=(PSINT(K)-PSINT(2))/(PSINT(LM)-PSINT(2))*PI
!         int_state%EPSQ2(K-1)=(1.+COS(ARG))*0.09+0.02
          int_state%EPSQ2(K-1)=0.02
          int_state%EPSL(K-1)=SQRT(int_state%EPSQ2(K-1)*0.5)
        ENDDO   
        int_state%EPSQ2(LM)=int_state%EPSQ2(LM-1)
        DEALLOCATE(PSINT)

!write(0,*) 'epsl=',epsl
!write(0,*) 'epsq2=',epsq2
!-----------------------------------------------------------------------
        do j=jts,jte
          do i=its,ite
            int_state%pdo(i,j)=int_state%pd(i,j)
          enddo
        enddo
        call halo_exch(int_state%pdo,1,2,2)
!
        do l=1,lm
          do j=jts,jte
            do i=its,ite
              int_state%up(i,j,l)=int_state%u(i,j,l)
              int_state%vp(i,j,l)=int_state%v(i,j,l)
              int_state%tp(i,j,l)=int_state%t(i,j,l)
            enddo
          enddo
        enddo
        call halo_exch(int_state%tp,lm,int_state%up,lm,int_state%vp,lm,2,2)
!-----------------------------------------------------------------------
        do l=1,lm
          do j=jms,jme
            do i=ims,ime
              int_state%q2(i,j,l)=0.02
              int_state%o3(i,j,l)=0.
              if(i.ge.ide  /2+1- 6.and.i.le.ide  /2+1+ 6.and. &
                 j.ge.jde*3/4+1- 6.and.j.le.jde*3/4+1+ 6.) then !global
!                 j.ge.jde  /2+1- 6.and.j.le.jde  /2+1+ 6.) then !regional
                int_state%o3(i,j,l)=10.
              endif
              int_state%dwdt(i,j,l)=1.
              int_state%w(i,j,l)=0.
            enddo
          enddo
        enddo
        call halo_exch(int_state%dwdt,lm,2,2)
!
        do j=jts,jte
          do i=its,ite
            int_state%pint(i,j,1)=int_state%pt
          enddo
        enddo
!
        do l=1,lm
          do j=jts,jte
            do i=its,ite
              int_state%pint(i,j,l+1)=int_state%pint(i,j,l)+int_state%dsg2(l)*int_state%pd(i,j)+int_state%pdsg1(l)
            enddo
          enddo
        enddo
        call halo_exch(int_state%pint,lm+1,2,2)
!
        call halo_exch(int_state%q2,lm,int_state%o3,lm,2,2)
        do l=1,lm
          do j=jms,jme
            do i=ims,ime
              int_state%tracers_prev(i,j,l,int_state%indx_q )=sqrt(max(int_state%q (i,j,l),0.))
              int_state%tracers_prev(i,j,l,int_state%indx_cw)=sqrt(max(int_state%cw(i,j,l),0.))
              int_state%tracers_prev(i,j,l,int_state%indx_q2)=sqrt(max(int_state%q2(i,j,l),0.))
            enddo
          enddo
        enddo
!
!-----------------------------------------------------------------------
!---reading surface data------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  ALBEDO
!-----------------------------------------------------------------------
!
      call getrecn(recname,reclevtyp,reclev,nrec,'albedo','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%ALBEDO(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!
!-----------------------------------------------------------------------
!***  ALBASE
!-----------------------------------------------------------------------
!
      call getrecn(recname,reclevtyp,reclev,nrec,'albase','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%ALBASE(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!
!-----------------------------------------------------------------------
!*** EPSR
!-----------------------------------------------------------------------
!
      call getrecn(recname,reclevtyp,reclev,nrec,'epsr','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%EPSR(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!
!-----------------------------------------------------------------------
!*** MXSNAL
!-----------------------------------------------------------------------
!
      call getrecn(recname,reclevtyp,reclev,nrec,'mxsnal','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%MXSNAL(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!
!-----------------------------------------------------------------------
!*** TSKIN
!-----------------------------------------------------------------------
!
      call getrecn(recname,reclevtyp,reclev,nrec,'tskin','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%TSKIN(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!
!-----------------------------------------------------------------------
!*** SST
!-----------------------------------------------------------------------
!
      call getrecn(recname,reclevtyp,reclev,nrec,'tsea','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%SST(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!
!-----------------------------------------------------------------------
!***  SNO and SNOWC
!-----------------------------------------------------------------------
!
      call getrecn(recname,reclevtyp,reclev,nrec,'sno','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%SNO(i,j)=tmp(i-its+1+js+fldst)
            if(int_state%SNO(i,j).gt.0.) then          !2013
               int_state%SNOWC(i,j) = 0.98
            else
               int_state%SNOWC(i,j) = 0.0
            endif
          enddo
        enddo
      endif
!
!-----------------------------------------------------------------------
!***  SI
!-----------------------------------------------------------------------
!
      call getrecn(recname,reclevtyp,reclev,nrec,'si','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%SI(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!
!-----------------------------------------------------------------------
!***  SICE
!-----------------------------------------------------------------------
!
      call getrecn(recname,reclevtyp,reclev,nrec,'sice','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%SICE(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!      write(0,*)'aft sice, sice=',maxval(int_state%SICE),minval(int_state%SICE)
!
!-----------------------------------------------------------------------
!***  TG
!-----------------------------------------------------------------------
!
      call getrecn(recname,reclevtyp,reclev,nrec,'tg','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%TG(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!
!-----------------------------------------------------------------------
!***  CMC
!-----------------------------------------------------------------------
!
      call getrecn(recname,reclevtyp,reclev,nrec,'cmc','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%CMC(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!
!-----------------------------------------------------------------------
!***  SR
!-----------------------------------------------------------------------
!
      call getrecn(recname,reclevtyp,reclev,nrec,'sr','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%SR(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!
!-----------------------------------------------------------------------
!***  USTAR
!-----------------------------------------------------------------------
!
      call getrecn(recname,reclevtyp,reclev,nrec,'ustar','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%ustar(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!
!-----------------------------------------------------------------------
!***  Z0
!-----------------------------------------------------------------------
!
      call getrecn(recname,reclevtyp,reclev,nrec,'zorl','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%Z0(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
      CALL HALO_EXCH(int_state%Z0,1,3,3)
!
!-----------------------------------------------------------------------
!***  Z0BASE
!-----------------------------------------------------------------------
!
      call getrecn(recname,reclevtyp,reclev,nrec,'z0base','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%Z0BASE(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
      CALL HALO_EXCH(int_state%Z0BASE,1,3,3)
!
!-----------------------------------------------------------------------
!***  STC, SMC, SH2O
!-----------------------------------------------------------------------
!
      DO L=1,int_state%NSOIL
!
        call getrecn(recname,reclevtyp,reclev,nrec,'stc','soil layer',l,recn)
        if(recn>0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              int_state%STC(i,j,l)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
!
        call getrecn(recname,reclevtyp,reclev,nrec,'smc','soil layer',l,recn)
        if(recn>0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              int_state%SMC(i,j,l)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
!
        call getrecn(recname,reclevtyp,reclev,nrec,'sh2o','soil layer',l,recn)
        if(recn>0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              int_state%SH2O(i,j,l)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
!
      ENDDO
!
!-----------------------------------------------------------------------
!***  ISLTYP
!-----------------------------------------------------------------------
!
      call getrecn(recname,reclevtyp,reclev,nrec,'sltyp','sfc',1,recn)
      int_state%ISLTYP=0
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%ISLTYP(i,j)=INT(tmp(i-its+1+js+fldst))
          enddo
        enddo
      endif
!
!-----------------------------------------------------------------------
!***  IVGTYP
!-----------------------------------------------------------------------
!
      call getrecn(recname,reclevtyp,reclev,nrec,'vgtyp','sfc',1,recn)
      int_state%IVGTYP=0
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%IVGTYP(i,j)=INT(tmp(i-its+1+js+fldst))
          enddo
        enddo
      endif
!
!-----------------------------------------------------------------------
!***  VEGFRC
!-----------------------------------------------------------------------
!
      call getrecn(recname,reclevtyp,reclev,nrec,'vegfrc','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%VEGFRC(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!
!! FOR Hurricane application, U10, V10 at initial time are needed for running tracker
! Weiguo Wang 2014-06-22
 !        write(0,*)'read v10'
 !        write(0,*)'int_state%RUN_TC=',int_state%RUN_TC 
     if (int_state%RUN_TC) then
!-----------------------------------------------------------------------
!***  U10
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'u10','10 m above gnd',1,recn)
       int_state%u10=0.0
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%U10(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
      call halo_exch(int_state%U10,1,2,2)
     !    write(0,*)'read u10'
     !    write(0,*)int_state%U10(1:10,5)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  V10
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'v10','10 m above gnd',1,recn)
       int_state%v10=0.0
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%V10(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
      call halo_exch(int_state%v10,1,2,2)
      !   write(0,*)'read v10'
      !   write(0,*)int_state%v10(1:10,5)
     endif   !  if hurricane
!-----------------------------------------------------------------------
!
        call nemsio_getheadvar(gfile,'dlmd',int_state%dlmd,ierr)
        call nemsio_getheadvar(gfile,'dphd',int_state%dphd,ierr)
        call nemsio_getheadvar(gfile,'wbd',int_state%wbd,ierr)
        call nemsio_getheadvar(gfile,'sbd',int_state%sbd,ierr)
        call nemsio_getheadvar(gfile,'tlm0d',int_state%tlm0d,ierr)
        call nemsio_getheadvar(gfile,'tph0d',int_state%tph0d,ierr)
!
        call mpi_bcast(int_state%pt,1,mpi_real,0,mpi_comm_comp,irtn)
        call mpi_bcast(int_state%dlmd,1,mpi_real,0,mpi_comm_comp,irtn)
        call mpi_bcast(int_state%dphd,1,mpi_real,0,mpi_comm_comp,irtn)
        call mpi_bcast(int_state%wbd,1,mpi_real,0,mpi_comm_comp,irtn)
        call mpi_bcast(int_state%sbd,1,mpi_real,0,mpi_comm_comp,irtn)
        call mpi_bcast(int_state%tlm0d,1,mpi_real,0,mpi_comm_comp,irtn)
        call mpi_bcast(int_state%tph0d,1,mpi_real,0,mpi_comm_comp,irtn)
!
!-----------------------------------------------------------------------
!
        call nemsio_close(gfile,iret=ierr)
!
!-----------------------------------------------------------------------
      else read_blocks                          ! restart
!-----------------------------------------------------------------------
!
        write(infile,'(a,i2.2,a)')'restart_file_',my_domain_id,'_nemsio'
        call nemsio_open(gfile,infile,'read',mpi_comm_comp,iret=ierr)
        if(ierr/=0)then
          write(0,*)' Unable to open ',trim(infile),' in READ_NEMSIO'
          rc = ierr
          return
        endif
!
!-----------------------------------------------------------------------
!***  Read from restart file: Integer scalars
!-----------------------------------------------------------------------
        call nemsio_getheadvar(gfile,'FCSTDATE',FCSTDATE,ierr)
        iyear_fcst=FCSTDATE(1)
        imonth_fcst=FCSTDATE(2)
        iday_fcst=FCSTDATE(3)
        ihour_fcst=FCSTDATE(4)
        call nemsio_getheadvar(gfile,'IHRST',int_state%IHRST,ierr)
        call nemsio_getheadvar(gfile,'LPT2',int_state%LPT2,ierr)
        call nemsio_getheadvar(gfile,'I_PAR_STA',int_state%I_PAR_STA,ierr)
        call nemsio_getheadvar(gfile,'J_PAR_STA',int_state%J_PAR_STA,ierr)
        call nemsio_getheadvar(gfile,'LAST_STEP_MOVED',int_state%LAST_STEP_MOVED,ierr)
        call nemsio_getheadvar(gfile,'NMTS',int_state%NMTS,ierr)
!
        call nemsio_getheadvar(gfile,'TRACK_N_OLD',int_state%TRACK_N_OLD,ierr)
        call nemsio_getheadvar(gfile,'TRACKER_IFIX',int_state%TRACKER_IFIX,ierr)
        call nemsio_getheadvar(gfile,'TRACKER_JFIX',int_state%TRACKER_JFIX,ierr)
        call nemsio_getheadvar(gfile,'TRACKER_GAVE_UP',int_state%TRACKER_GAVE_UP,ierr)
        call nemsio_getheadvar(gfile,'TRACK_HAVE_GUESS',int_state%TRACK_HAVE_GUESS,ierr)
        call nemsio_getheadvar(gfile,'TRACKER_HAVEFIX',int_state%TRACKER_HAVEFIX,ierr)
!
!-----------------------------------------------------------------------
!***  Read from restart file: Integer 1D arrays
!-----------------------------------------------------------------------
        call nemsio_getheadvar(gfile,'IDAT',int_state%idat,ierr)
        call nemsio_getheadvar(gfile,'NTSCM',int_state%ntscm,ierr)
!
        if(mype==0)then
          write(0,*)'**** in read_nemsio *****************'
          write(0,*)' Restart year =',iyear_fcst
          write(0,*)' Restart month=',imonth_fcst
          write(0,*)' Restart day  =',iday_fcst
          write(0,*)' Restart hour =',ihour_fcst
          write(0,*)' Original start year =',int_state%idat(3)
          write(0,*)' Original start month=',int_state%idat(2)
          write(0,*)' Original start day  =',int_state%idat(1)
          write(0,*)' Original start hour =',int_state%ihrst
          write(0,*)' Timestep   =',int_state%dt
          write(0,*)' Steps/hour =',3600./int_state%dt
          write(0,*)'*************************************'
        endif
!
!-----------------------------------------------------------------------
!***  Read from restart file: Real scalars
!-----------------------------------------------------------------------
        call nemsio_getheadvar(gfile,'PDTOP',int_state%pdtop,ierr)
        call nemsio_getheadvar(gfile,'PT',int_state%pt,ierr)
        call nemsio_getheadvar(gfile,'TLM0D',int_state%tlm0d,ierr)
        call nemsio_getheadvar(gfile,'TPH0D',int_state%tph0d,ierr)
        call nemsio_getheadvar(gfile,'DPHD',int_state%dphd,ierr)
        call nemsio_getheadvar(gfile,'DLMD',int_state%dlmd,ierr)
        call nemsio_getheadvar(gfile,'SBD',int_state%sbd,ierr)
        call nemsio_getheadvar(gfile,'WBD',int_state%wbd,ierr)
!
        call nemsio_getheadvar(gfile,'TRACK_LAST_HOUR',int_state%TRACK_LAST_HOUR,ierr)
        call nemsio_getheadvar(gfile,'TRACK_GUESS_LAT',int_state%TRACK_GUESS_LAT,ierr)
        call nemsio_getheadvar(gfile,'TRACK_GUESS_LON',int_state%TRACK_GUESS_LON,ierr)
        call nemsio_getheadvar(gfile,'TRACK_EDGE_DIST',int_state%TRACK_EDGE_DIST,ierr)
        call nemsio_getheadvar(gfile,'TRACK_STDERR_M1',int_state%TRACK_STDERR_M1,ierr)
        call nemsio_getheadvar(gfile,'TRACK_STDERR_M2',int_state%TRACK_STDERR_M2,ierr)
        call nemsio_getheadvar(gfile,'TRACK_STDERR_M3',int_state%TRACK_STDERR_M3,ierr)
        call nemsio_getheadvar(gfile,'TRACKER_FIXLAT',int_state%TRACKER_FIXLAT,ierr)
        call nemsio_getheadvar(gfile,'TRACKER_FIXLON',int_state%TRACKER_FIXLON,ierr)
        call nemsio_getheadvar(gfile,'TRACKER_RMW',int_state%TRACKER_RMW,ierr)
        call nemsio_getheadvar(gfile,'TRACKER_PMIN',int_state%TRACKER_PMIN,ierr)
        call nemsio_getheadvar(gfile,'TRACKER_VMAX',int_state%TRACKER_VMAX,ierr)
!
!-----------------------------------------------------------------------
!***  Read from restart file: Real 1D arrays
!-----------------------------------------------------------------------
        call nemsio_getheadvar(gfile,'SG1',int_state%sg1,ierr)
        call nemsio_getheadvar(gfile,'SG2',int_state%sg2,ierr)
        call nemsio_getheadvar(gfile,'DSG1',int_state%dsg1,ierr)
        call nemsio_getheadvar(gfile,'DSG2',int_state%dsg2,ierr)
        call nemsio_getheadvar(gfile,'SGML1',int_state%sgml1,ierr)
        call nemsio_getheadvar(gfile,'SGML2',int_state%sgml2,ierr)
        call nemsio_getheadvar(gfile,'SGM',int_state%sgm,ierr)
        call nemsio_getheadvar(gfile,'EPSL',int_state%EPSL,ierr)
        call nemsio_getheadvar(gfile,'EPSQ2',int_state%EPSQ2,ierr)
        CALL NEMSIO_GETHEADVAR(gfile,'SLDPTH',int_state%SLDPTH,iret=irtn)
        CALL NEMSIO_GETHEADVAR(gfile,'MP_RESTART',int_state%MP_RESTART_STATE,iret=irtn)
        CALL NEMSIO_GETHEADVAR(gfile,'TBPVS_STAT',int_state%TBPVS_STATE,iret=irtn)
        CALL NEMSIO_GETHEADVAR(gfile,'TBPVS0_STA',int_state%TBPVS0_STATE,iret=irtn)
!
        DO L=1,LM
          int_state%PDSG1(L)=int_state%DSG1(L)*int_state%PDTOP
          int_state%PSGML1(L)=int_state%SGML1(L)*int_state%PDTOP+int_state%PT
        ENDDO
!
        DO L=1,LM+1
          int_state%PSG1(L)=int_state%SG1(L)*int_state%PDTOP+int_state%PT
        ENDDO
!
        CALL NEMSIO_GETHEADVAR(gfile,'TRACK_OLD_NTSD',int_state%TRACK_OLD_NTSD,iret=irtn)
        CALL NEMSIO_GETHEADVAR(gfile,'TRACK_OLD_LAT',int_state%TRACK_OLD_LAT,iret=irtn)
        CALL NEMSIO_GETHEADVAR(gfile,'TRACK_OLD_LON',int_state%TRACK_OLD_LON,iret=irtn)
!
!-----------------------------------------------------------------------
!***  Read in the full-domain 1-D string of boundary data.
!***  Each task isolates its own piece of that data.
!-----------------------------------------------------------------------
!
        length=(int_state%nlev_h*int_state%lnsh                         &  !<-- Total # of words
               +int_state%nlev_v*int_state%lnsv)                        &  !    in full-domain
               *2*2*(ide-ids+jde-jds+2)                                    !    boundary arrays.
        allocate(all_bc_data(1:length))
!
        call nemsio_getheadvar(gfile,'ALL_BC_DATA',all_bc_data,ierr)
!
!-----------------------------------------------------------------------
!
        nvars_bc_2d_h=int_state%nvars_bc_2d_h
        nvars_bc_3d_h=int_state%nvars_bc_3d_h
        nvars_bc_4d_h=int_state%nvars_bc_4d_h
        nvars_bc_2d_v=int_state%nvars_bc_2d_v
        nvars_bc_3d_v=int_state%nvars_bc_3d_v
!
        kount=0
!
!-----------------------------------------------------------------------
!
        iend=min(ite_h2,ide)
!
        if(nvars_bc_2d_h>0)then
          do nv=1,nvars_bc_2d_h
            do n=1,2
            do j=1,int_state%lnsh
            do i=ids,ide
              kount=kount+1
              if(jts==jds.and.i>=its_h2.and.i<=iend)then                   !<-- South boundary tasks extract 2-D BC H-pt data
                int_state%bnd_vars_h%var_2d(nv)%south(i,j,n)=           &
                  all_bc_data(kount)
              endif
            enddo
            enddo
            enddo
          enddo
        endif
!
        if(nvars_bc_3d_h>0)then
          do nv=1,nvars_bc_3d_h
            do n=1,2
            do l=1,lm
            do j=1,int_state%lnsh
            do i=ids,ide
              kount=kount+1
              if(jts==jds.and.i>=its_h2.and.i<=iend)then                   !<-- South boundary tasks extract 3-D BC H-pt data
                int_state%bnd_vars_h%var_3d(nv)%south(i,j,l,n)=         &
                  all_bc_data(kount)
              endif
            enddo
            enddo
            enddo
            enddo
          enddo
        endif
!
        if(nvars_bc_4d_h>0)then
          do nv=1,nvars_bc_4d_h
            lb=int_state%lbnd_4d(nv)
            ub=int_state%ubnd_4d(nv)
            do nl=lb,ub
            do n=1,2
            do l=1,lm
            do j=1,int_state%lnsh
            do i=ids,ide
              kount=kount+1
              if(jts==jds.and.i>=its_h2.and.i<=iend)then                   !<-- South boundary tasks extract 4-D BC H-pt data
                int_state%bnd_vars_h%var_4d(nv)%south(i,j,l,n,nl)=      &
                  all_bc_data(kount)
              endif
            enddo
            enddo
            enddo
            enddo
            enddo
          enddo
        endif
!
        if(nvars_bc_2d_v>0)then
          do nv=1,nvars_bc_2d_v
            do n=1,2
            do j=1,int_state%lnsv
            do i=ids,ide
              kount=kount+1
              if(jts==jds.and.i>=its_h2.and.i<=iend)then                   !<-- South boundary tasks extract 2-D BC V-pt data
                int_state%bnd_vars_v%var_2d(nv)%south(i,j,n)=           &
                  all_bc_data(kount)
              endif
            enddo
            enddo
            enddo
          enddo
        endif
!
        if(nvars_bc_3d_v>0)then
          do nv=1,nvars_bc_3d_v
            do n=1,2
            do l=1,lm
            do j=1,int_state%lnsv
            do i=ids,ide
              kount=kount+1
              if(jts==jds.and.i>=its_h2.and.i<=iend)then                   !<-- South boundary tasks extract 3-D BC V-pt data
                int_state%bnd_vars_v%var_3d(nv)%south(i,j,l,n)=         &
                  all_bc_data(kount)
              endif
            enddo
            enddo
            enddo
            enddo
          enddo
        endif
!
        if(nvars_bc_2d_h>0)then
          do nv=1,nvars_bc_2d_h
            do n=1,2
            do j=1,int_state%lnsh
            do i=ids,ide
              kount=kount+1
              if(jte==jde.and.i>=its_h2.and.i<=iend)then                   !<-- North boundary tasks extract 2-D BC H-pt data
                int_state%bnd_vars_h%var_2d(nv)%north(i,j,n)=           &
                  all_bc_data(kount)
              endif
            enddo
            enddo
            enddo
          enddo
        endif
!
        if(nvars_bc_3d_h>0)then
          do nv=1,nvars_bc_3d_h
            do n=1,2
            do l=1,lm
            do j=1,int_state%lnsh
            do i=ids,ide
              kount=kount+1
              if(jte==jde.and.i>=its_h2.and.i<=iend)then                   !<-- North boundary tasks extract 3-D BC H-pt data
                int_state%bnd_vars_h%var_3d(nv)%north(i,j,l,n)=         &
                  all_bc_data(kount)
              endif
            enddo
            enddo
            enddo
            enddo
          enddo
        endif
!
        if(nvars_bc_4d_h>0)then
          do nv=1,nvars_bc_4d_h
            lb=int_state%lbnd_4d(nv)
            ub=int_state%ubnd_4d(nv)
            do nl=lb,ub
            do n=1,2
            do l=1,lm
            do j=1,int_state%lnsh
            do i=ids,ide
              kount=kount+1
              if(jte==jde.and.i>=its_h2.and.i<=iend)then                   !<-- North boundary tasks extract 4-D BC H-pt data
                int_state%bnd_vars_h%var_4d(nv)%north(i,j,l,n,nl)=      &
                  all_bc_data(kount)
              endif
            enddo
            enddo
            enddo
            enddo
            enddo
          enddo
        endif
!
        if(nvars_bc_2d_v>0)then
          do nv=1,nvars_bc_2d_v
            do n=1,2
            do j=1,int_state%lnsv
            do i=ids,ide
              kount=kount+1
              if(jte==jde.and.i>=its_h2.and.i<=iend)then                   !<-- North boundary tasks extract 2-D BC V-pt data
                int_state%bnd_vars_v%var_2d(nv)%north(i,j,n)=           &
                  all_bc_data(kount)
              endif
            enddo
            enddo
            enddo
          enddo
        endif
!
        if(nvars_bc_3d_v>0)then
          do nv=1,nvars_bc_3d_v
            do n=1,2
            do l=1,lm
            do j=1,int_state%lnsv
            do i=ids,ide
              kount=kount+1
              if(jte==jde.and.i>=its_h2.and.i<=iend)then                   !<-- North boundary tasks extract 3-D BC V-pt data
                int_state%bnd_vars_v%var_3d(nv)%north(i,j,l,n)=         &
                  all_bc_data(kount)
              endif
            enddo
            enddo
            enddo
            enddo
          enddo
        endif
!
        jend=min(jte_h2,jde)
!
        if(nvars_bc_2d_h>0)then
          do nv=1,nvars_bc_2d_h
            do n=1,2
            do j=jds,jde
            do i=1,int_state%lnsh
              kount=kount+1
              if(its==ids.and.j>=jts_h2.and.j<=jend)then                   !<-- West boundary tasks extract 2-D BC H-pt data
                int_state%bnd_vars_h%var_2d(nv)%west(i,j,n)=            &
                  all_bc_data(kount)
              endif
            enddo
            enddo
            enddo
          enddo
        endif
!
        if(nvars_bc_3d_h>0)then
          do nv=1,nvars_bc_3d_h
            do n=1,2
            do l=1,lm
            do j=jds,jde
            do i=1,int_state%lnsh
              kount=kount+1
              if(its==ids.and.j>=jts_h2.and.j<=jend)then                   !<-- West boundary tasks extract 3-D BC H-pt data
                int_state%bnd_vars_h%var_3d(nv)%west(i,j,l,n)=          &
                  all_bc_data(kount)
              endif
            enddo
            enddo
            enddo
            enddo
          enddo
        endif
!
        if(nvars_bc_4d_h>0)then
          do nv=1,nvars_bc_4d_h
            lb=int_state%lbnd_4d(nv)
            ub=int_state%ubnd_4d(nv)
            do nl=lb,ub
            do n=1,2
            do l=1,lm
            do j=jds,jde
            do i=1,int_state%lnsh
              kount=kount+1
              if(its==ids.and.j>=jts_h2.and.j<=jend)then                   !<-- West boundary tasks extract 4-D BC H-pt data
                int_state%bnd_vars_h%var_4d(nv)%west(i,j,l,n,nl)=       &
                  all_bc_data(kount)
              endif
            enddo
            enddo
            enddo
            enddo
            enddo
          enddo
        endif
!
        if(nvars_bc_2d_v>0)then
          do nv=1,nvars_bc_2d_v
            do n=1,2
            do j=jds,jde
            do i=1,int_state%lnsv
              kount=kount+1
              if(its==ids.and.j>=jts_h2.and.j<=jend)then                   !<-- West boundary tasks extract 2-D BC V-pt data
                int_state%bnd_vars_v%var_2d(nv)%west(i,j,n)=            &
                  all_bc_data(kount)
              endif
            enddo
            enddo
            enddo
          enddo
        endif
!
        if(nvars_bc_3d_v>0)then
          do nv=1,nvars_bc_3d_v
            do n=1,2
            do l=1,lm
            do j=jds,jde
            do i=1,int_state%lnsv
              kount=kount+1
              if(its==ids.and.j>=jts_h2.and.j<=jend)then                   !<-- West boundary tasks extract 3-D BC V-pt data
                int_state%bnd_vars_v%var_3d(nv)%west(i,j,l,n)=          &
                  all_bc_data(kount)
              endif
            enddo
            enddo
            enddo
            enddo
          enddo
        endif
!
        if(nvars_bc_2d_h>0)then
          do nv=1,nvars_bc_2d_h
            do n=1,2
            do j=jds,jde
            do i=1,int_state%lnsh
              kount=kount+1
              if(ite==ide.and.j>=jts_h2.and.j<=jend)then                   !<-- East boundary tasks extract 2-D BC H-pt data
                int_state%bnd_vars_h%var_2d(nv)%east(i,j,n)=            &
                  all_bc_data(kount)
              endif
            enddo
            enddo
            enddo
          enddo
        endif
!
        if(nvars_bc_3d_h>0)then
          do nv=1,nvars_bc_3d_h
            do n=1,2
            do l=1,lm
            do j=jds,jde
            do i=1,int_state%lnsh
              kount=kount+1
              if(ite==ide.and.j>=jts_h2.and.j<=jend)then                   !<-- East boundary tasks extract 3-D BC H-pt data
                int_state%bnd_vars_h%var_3d(nv)%east(i,j,l,n)=          &
                  all_bc_data(kount)
              endif
            enddo
            enddo
            enddo
            enddo
          enddo
        endif
!
        if(nvars_bc_4d_h>0)then
          do nv=1,nvars_bc_4d_h
            lb=int_state%lbnd_4d(nv)
            ub=int_state%ubnd_4d(nv)
            do nl=lb,ub
            do n=1,2
            do l=1,lm
            do j=jds,jde
            do i=1,int_state%lnsh
              kount=kount+1
              if(ite==ide.and.j>=jts_h2.and.j<=jend)then                   !<-- East boundary tasks extract 4-D BC H-pt data
                int_state%bnd_vars_h%var_4d(nv)%east(i,j,l,n,nl)=       &
                  all_bc_data(kount)
              endif
            enddo
            enddo
            enddo
            enddo
            enddo
          enddo
        endif
!
        if(nvars_bc_2d_v>0)then
          do nv=1,nvars_bc_2d_v
            do n=1,2
            do j=jds,jde
            do i=1,int_state%lnsv
              kount=kount+1
              if(ite==ide.and.j>=jts_h2.and.j<=jend)then                   !<-- East boundary tasks extract 2-D BC V-pt data
                int_state%bnd_vars_v%var_2d(nv)%east(i,j,n)=            &
                  all_bc_data(kount)
              endif
            enddo
            enddo
            enddo
          enddo
        endif
!
        if(nvars_bc_3d_v>0)then
          do nv=1,nvars_bc_3d_v
            do n=1,2
            do l=1,lm
            do j=jds,jde
            do i=1,int_state%lnsv
              kount=kount+1
              if(ite==ide.and.j>=jts_h2.and.j<=jend)then                   !<-- East boundary tasks extract 3-D BC V-pt data
                int_state%bnd_vars_v%var_3d(nv)%east(i,j,l,n)=          &
                  all_bc_data(kount)
              endif
            enddo
            enddo
            enddo
            enddo
          enddo
        endif
!
        deallocate(all_bc_data)
!
!-----------------------------------------------------------------------
!***  Read from restart file: Logical
!-----------------------------------------------------------------------
        call nemsio_getheadvar(gfile,'RUN',int_state%run,ierr)
!-----------------------------------------------------------------------
!***  Read from restart file: Real 2D arrays
!-----------------------------------------------------------------------
!
      fldsize=(jte-jts+1)*(ite-its+1)
!
      call nemsio_getfilehead(gfile,nrec=nrec,iret=ierr)
      allocate(recname(nrec),reclevtyp(nrec),reclev(nrec))
      call nemsio_getfilehead(gfile,recname=recname,reclevtyp=reclevtyp,   &
                              reclev=reclev,iret=ierr)
!
      allocate(tmp(fldsize*nrec))
      tmp=0.
      call nemsio_denseread(gfile,its,ite,jts,jte,tmp,iret=ierr)
!
!-----------------------------------------------------------------------
!***  close nemsio file
!-----------------------------------------------------------------------
!
      CALL NEMSIO_CLOSE(GFILE)
!-----------------------------------------------------------------------
!
      CALL NEMSIO_FINALIZE()
!
!-----------------------------------------------------------------------
!***  assign data: Integer 2D arrays
!-----------------------------------------------------------------------
!
      call getrecn(recname,reclevtyp,reclev,nrec,'sltyp','sfc',1,recn)
      int_state%ISLTYP=0
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%ISLTYP(i,j)=NINT(tmp(i-its+1+js+fldst))
          enddo
        enddo
      endif
!
      call getrecn(recname,reclevtyp,reclev,nrec,'vgtyp','sfc',1,recn)
      int_state%IVGTYP=0
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%IVGTYP(i,j)=NINT(tmp(i-its+1+js+fldst))
          enddo
        enddo
      endif
!
      call getrecn(recname,reclevtyp,reclev,nrec,'cfrcv','sfc',1,recn)
      int_state%NCFRCV=0
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%NCFRCV(i,j)=NINT(tmp(i-its+1+js+fldst))
          enddo
        enddo
      endif
!
      call getrecn(recname,reclevtyp,reclev,nrec,'cfrst','sfc',1,recn)
      int_state%NCFRST=0
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%NCFRST(i,j)=NINT(tmp(i-its+1+js+fldst))
          enddo
        enddo
      endif
!
      call getrecn(recname,reclevtyp,reclev,nrec,'tracker_fixes','sfc',1,recn)
      int_state%TRACKER_FIXES=0
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%TRACKER_FIXES(i,j)=NINT(tmp(i-its+1+js+fldst))
          enddo
        enddo
      endif
!
!-----------------------------------------------------------------------
!***  Assign data: Real 2D arrays
!-----------------------------------------------------------------------
!-- fis
      int_state%fis=0.
      call getrecn(recname,reclevtyp,reclev,nrec,'fis','sfc',1,recn)
      if(recn/=0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%fis(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
      call halo_exch(int_state%fis,1,2,2)
!-----------------------------------------------------------------------
!-- glat
      int_state%glat=0.
      call getrecn(recname,reclevtyp,reclev,nrec,'glat','sfc',1,recn)
      if(recn/=0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%glat(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
      call halo_exch(int_state%glat,1,2,2)
!-----------------------------------------------------------------------
!-- glon
      int_state%glon=0.
      call getrecn(recname,reclevtyp,reclev,nrec,'glon','sfc',1,recn)
      if(recn/=0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%glon(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
      call halo_exch(int_state%glon,1,2,2)
!-----------------------------------------------------------------------
!-- vlat
      int_state%vlat=0.
      call getrecn(recname,reclevtyp,reclev,nrec,'vlat','sfc',1,recn)
      if(recn/=0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%vlat(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
      call halo_exch(int_state%vlat,1,2,2)
!-----------------------------------------------------------------------
!-- vlon
      int_state%vlon=0.
      call getrecn(recname,reclevtyp,reclev,nrec,'vlon','sfc',1,recn)
      if(recn/=0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%vlon(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
      call halo_exch(int_state%vlon,1,2,2)
!-----------------------------------------------------------------------
!-- hdacx
      int_state%hdacx=0.
      call getrecn(recname,reclevtyp,reclev,nrec,'hdacx','sfc',1,recn)
      if(recn/=0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%hdacx(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
      call halo_exch(int_state%hdacx,1,2,2)
!-----------------------------------------------------------------------
!-- hdacy
      int_state%hdacy=0.
      call getrecn(recname,reclevtyp,reclev,nrec,'hdacy','sfc',1,recn)
      if(recn/=0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%hdacy(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
      call halo_exch(int_state%hdacy,1,2,2)
!-----------------------------------------------------------------------
!-- f
      int_state%f=0.
      call getrecn(recname,reclevtyp,reclev,nrec,'f','sfc',1,recn)
      if(recn/=0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%f(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
      call halo_exch(int_state%f,1,2,2)
!-----------------------------------------------------------------------
!-- hdacvx
      int_state%hdacvx=0.
      call getrecn(recname,reclevtyp,reclev,nrec,'hdacvx','sfc',1,recn)
      if(recn/=0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%hdacvx(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
      call halo_exch(int_state%hdacvx,1,2,2)
!-----------------------------------------------------------------------
!-- hdacvy
      int_state%hdacvy=0.
      call getrecn(recname,reclevtyp,reclev,nrec,'hdacvy','sfc',1,recn)
      if(recn/=0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%hdacvy(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
      call halo_exch(int_state%hdacvy,1,2,2)
!-----------------------------------------------------------------------
!--pd
      int_state%pd=0.
      call getrecn(recname,reclevtyp,reclev,nrec,'dpres','hybrid sig lev',1,recn)
      if(recn/=0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%pd(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
      call halo_exch(int_state%pd,1,2,2)
!-----------------------------------------------------------------------
!-- pdo
      int_state%pdo=0.
      call getrecn(recname,reclevtyp,reclev,nrec,'pdo','sfc',1,recn)
      if(recn/=0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%pdo(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
      call halo_exch(int_state%pdo,1,2,2)
!
!-----------------------------------------------------------------------
!***  Read from restart file: Real 2D arrays (contd.)
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  ACFRCV
!-----------------------------------------------------------------------
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%ACFRCV(I,J)=0.
      ENDDO
      ENDDO
      call getrecn(recname,reclevtyp,reclev,nrec,'acfrcv','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%ACFRCV(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  ACFRST
!-----------------------------------------------------------------------
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%ACFRST(I,J)=0.
      ENDDO
      ENDDO
      call getrecn(recname,reclevtyp,reclev,nrec,'acfrst','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%ACFRST(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  ACPREC
!-----------------------------------------------------------------------
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%ACPREC(I,J)=0.
      ENDDO
      ENDDO
      call getrecn(recname,reclevtyp,reclev,nrec,'acprec','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%ACPREC(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  ACPREC_TOT
!-----------------------------------------------------------------------
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%ACPREC_TOT(I,J)=0.
      ENDDO
      ENDDO
      call getrecn(recname,reclevtyp,reclev,nrec,'acprec_tot','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%ACPREC_TOT(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  ACSNOM
!-----------------------------------------------------------------------
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%ACSNOM(I,J)=0.
      ENDDO
      ENDDO
      call getrecn(recname,reclevtyp,reclev,nrec,'acsnom','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%ACSNOM(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  ACSNOW
!-----------------------------------------------------------------------
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%ACSNOW(I,J)=0.
      ENDDO
      ENDDO
      call getrecn(recname,reclevtyp,reclev,nrec,'acsnow','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%ACSNOW(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  AKHS_OUT
!-----------------------------------------------------------------------
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%AKHS_OUT(I,J)=0.
      ENDDO
      ENDDO
      call getrecn(recname,reclevtyp,reclev,nrec,'akhs_out','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%AKHS_OUT(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  AKMS_OUT
!-----------------------------------------------------------------------
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%AKMS_OUT(I,J)=0.
      ENDDO
      ENDDO
      call getrecn(recname,reclevtyp,reclev,nrec,'akms_out','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%AKMS_OUT(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  ALBASE
!-----------------------------------------------------------------------
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%ALBASE(I,J)=0.
      ENDDO
      ENDDO
      call getrecn(recname,reclevtyp,reclev,nrec,'albase','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%ALBASE(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  ALBEDO
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'albedo','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%ALBEDO(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  ALWIN
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'alwin','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%ALWIN(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  ALWOUT
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'alwout','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%ALWOUT(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  ALWTOA
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'alwtoa','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%ALWTOA(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  ASWIN
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'aswin','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%ASWIN(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  ASWOUT
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'aswout','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%ASWOUT(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  ASWTOA
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'aswtoa','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%ASWTOA(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  BGROFF
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'bgroff','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%BGROFF(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  CFRACH
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'cfrach','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%CFRACH(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  CFRACL
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'cfracl','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%CFRACL(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  CFRACM
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'cfracm','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%CFRACM(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  CLDEFI
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'cldefi','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%CLDEFI(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  CMC
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'cmc','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%CMC(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  CNVBOT
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'cnvbot','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%CNVBOT(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  CNVTOP
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'cnvtop','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%CNVTOP(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  CPRATE
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'cprate','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%CPRATE(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  CUPPT
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'cuppt','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%CUPPT(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  CUPREC
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'cuprec','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%CUPREC(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  CZEN
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'czen','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%CZEN(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  CZMEAN
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'czmean','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%CZMEAN(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  EPSR
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'epsr','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%EPSR(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  GRNFLX
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'grnflx','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%GRNFLX(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  HBOTD
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'hbotd','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%HBOTD(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  HBOTS
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'hbots','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%HBOTS(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  HTOPD
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'htopd','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%HTOPD(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  HTOPS
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'htops','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%HTOPS(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  SNOW ALBEDO
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'mxsnal','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%MXSNAL(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  PBLH
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'pblh','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%PBLH(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  POTEVP
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'potevp','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%POTEVP(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  PREC
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'prec','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%PREC(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  PSHLTR
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'pshltr','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%PSHLTR(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  Q10
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'q10','10 m above gnd',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%Q10(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  QSH
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'qsh','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%QSH(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  QSHLTR
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'qshltr','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%QSHLTR(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  QWBS
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'qwbs','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%QWBS(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  QZ0
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'qz0','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%QZ0(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  RADOT
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'radot','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%RADOT(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  RLWIN
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'rlwin','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%RLWIN(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  RLWTOA
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'rlwtoa','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%RLWTOA(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  RSWIN
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'rswin','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%RSWIN(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  RSWINC
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'rswinc','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%rswinc(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  RSWOUT
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'rswout','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%RSWOUT(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  SFCEVP
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'sfcevp','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%SFCEVP(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  SFCEXC
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'sfcexc','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%SFCEXC(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  SFCLHX
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'sfclhx','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%SFCLHX(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  SFCSHX
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'sfcshx','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%SFCSHX(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  SI
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'si','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%SI(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  SICE
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'sice','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%SICE(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  SIGT4
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'sigt4','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%SIGT4(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  SM (Seamask)
!-----------------------------------------------------------------------
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%SM(I,J)=0.
      ENDDO
      ENDDO
      call getrecn(recname,reclevtyp,reclev,nrec,'sm','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%SM(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  SMSTAV
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'smstav','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%SMSTAV(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  SMSTOT
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'smstot','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%SMSTOT(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  SNO
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'sno','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%SNO(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!
!-----------------------------------------------------------------------
!***  SNOWC
!-----------------------------------------------------------------------
!
      call getrecn(recname,reclevtyp,reclev,nrec,'snowc','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%SNOWC(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!
!-----------------------------------------------------------------------
!***  SNOPCX
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'snopcx','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%SNOPCX(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  SOILTB
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'soiltb','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%SOILTB(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  SR
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'sr','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%SR(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  SSROFF
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'ssroff','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%SSROFF(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  SST
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'tsea','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%SST(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  SUBSHX
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'subshx','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%SUBSHX(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  TG
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'tg','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%TG(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  TH10
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'th10','10 m above gnd',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%TH10(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  THS
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'ths','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%THS(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  THZ0
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'thz0','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%THZ0(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  TSHLTR
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'tshltr','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%TSHLTR(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  TWBS
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'twbs','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%TWBS(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  U10
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'u10','10 m above gnd',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%U10(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  USTAR
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'uustar','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%USTAR(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  UZ0
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'uz0','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%UZ0(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
      CALL HALO_EXCH(int_state%UZ0,1,3,3)
!-----------------------------------------------------------------------
!***  V10
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'v10','10 m above gnd',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%V10(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  VEGFRC
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'vegfrc','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%VEGFRC(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  VZ0
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'vz0','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%VZ0(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
      CALL HALO_EXCH(int_state%VZ0,1,3,3)
!-----------------------------------------------------------------------
!***  Z0
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'zorl','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%Z0(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
      CALL HALO_EXCH(int_state%Z0,1,3,3)
!-----------------------------------------------------------------------
!***  TSKIN
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'tskin','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%TSKIN(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  AKHS
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'akhs','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%AKHS(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  AKMS
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'akms','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%AKMS(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  HBOT
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'hbot','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%HBOT(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  HTOP
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'htop','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%HTOP(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  RSWTOA
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'rswtoa','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%RSWTOA(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  POTFLX
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'potflx','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%POTFLX(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  RMOL
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'rmol','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%RMOL(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  T2
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'t2','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%T2(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  Z0BASE
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'z0base','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%z0base(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  TLMIN
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'tlmin','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%TLMIN(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  TLMAX
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'tlmax','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%TLMAX(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  ACUTIM
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'acutim','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%ACUTIM(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  APHTIM
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'aphtim','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%APHTIM(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  ARDLW
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'ardlw','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%ARDLW(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  ARDSW
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'ardsw','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%ARDSW(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  ASRFC
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'asrfc','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%ASRFC(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  AVRAIN
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'avrain','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%AVRAIN(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  AVCNVC
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'avcnvc','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%AVCNVC(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  M10RV
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'m10rv','10 m above gnd',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%M10RV(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  M10WIND
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'m10wind','10 m above gnd',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%M10WIND(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  MEMBRANE_MSLP
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'membrane_mslp','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%MEMBRANE_MSLP(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  P500U
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'p500u','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%P500U(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  P500V
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'p500v','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%P500V(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  P700RV
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'p700rv','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%P700RV(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  P700U
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'p700u','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%P700U(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  P700V
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'p700v','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%P700V(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  P700WIND
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'p700wind','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%P700WIND(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  P700Z
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'p700z','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%P700Z(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  P850RV
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'p850rv','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%P850RV(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  P850U
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'p850u','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%P850U(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  P850V
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'p850v','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%P850V(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  P850WIND
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'p850wind','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%P850WIND(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  P850Z
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'p850z','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%P850Z(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  SM10RV
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'sm10rv','10 m above gnd',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%SM10RV(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  SM10WIND
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'sm10wind','10 m above gnd',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%SM10WIND(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  SMSLP
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'smslp','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%SMSLP(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  SP700RV
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'sp700rv','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%SP700RV(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  SP700WIND
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'sp700wind','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%SP700WIND(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  SP700Z
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'sp700z','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%SP700Z(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  SP850RV
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'sp850rv','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%SP850RV(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  SP850WIND
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'sp850wind','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%SP850WIND(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  SP850Z
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'sp850z','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%SP850Z(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  TRACKER_ANGLE
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'tracker_angle','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%tracker_angle(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  TRACKER_DISTSQ
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'tracker_distsq','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%tracker_distsq(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  Read from restart file: Real 3D arrays (only DYN)
!-----------------------------------------------------------------------
!w
      int_state%w=0.
      do l=1,lm
        call getrecn(recname,reclevtyp,reclev,nrec,'vvel','mid layer',l,recn)
        if(recn/=0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              int_state%w(i,j,l)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
      enddo
      call halo_exch(int_state%w,lm,2,2)
!
!-----------------------------------------------------------------------
!-- dwdt
      int_state%dwdt=0.
      do l=1,lm
        call getrecn(recname,reclevtyp,reclev,nrec,'dwdt','mid layer',l,recn)
        if(recn/=0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              int_state%dwdt(i,j,l)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
      enddo
      call halo_exch(int_state%dwdt,lm,2,2)
!-----------------------------------------------------------------------
!-- pres
      int_state%pint=0.
      do l=1,lm+1
        call getrecn(recname,reclevtyp,reclev,nrec,'pres','layer',l,recn)
        if(recn/=0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              int_state%pint(i,j,l)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
      enddo
      call halo_exch(int_state%pint,lm+1,2,2)
!-----------------------------------------------------------------------
!-- omgalf
      int_state%omgalf=0.
      do l=1,lm
        call getrecn(recname,reclevtyp,reclev,nrec,'omgalf','mid layer',l,recn)
        if(recn/=0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              int_state%omgalf(i,j,l)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
      enddo
      call halo_exch(int_state%omgalf,lm,2,2)
!-----------------------------------------------------------------------
!-- o3mr
      int_state%o3=0.
      do l=1,lm
        call getrecn(recname,reclevtyp,reclev,nrec,'o3mr','mid layer',l,recn)
        if(recn/=0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              int_state%o3(i,j,l)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
      enddo
      call halo_exch(int_state%o3,lm,2,2)
!-----------------------------------------------------------------------
!-- div
      int_state%div=0.
      do l=1,lm
        call getrecn(recname,reclevtyp,reclev,nrec,'div','mid layer',l,recn)
        if(recn/=0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              int_state%div(i,j,l)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
      enddo
      call halo_exch(int_state%div,lm,2,2)
!-----------------------------------------------------------------------
!-- tcu
      int_state%tcu=0.
      do l=1,lm
        call getrecn(recname,reclevtyp,reclev,nrec,'tcu','mid layer',l,recn)
        if(recn/=0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              int_state%tcu(i,j,l)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
      enddo
      call halo_exch(int_state%tcu,lm,2,2)
!-----------------------------------------------------------------------
!-- tcv
      int_state%tcv=0.
      do l=1,lm
        call getrecn(recname,reclevtyp,reclev,nrec,'tcv','mid layer',l,recn)
        if(recn/=0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              int_state%tcv(i,j,l)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
      enddo
      call halo_exch(int_state%tcv,lm,2,2)
!-----------------------------------------------------------------------
!-- tct
      int_state%tct=0.
      do l=1,lm
        call getrecn(recname,reclevtyp,reclev,nrec,'tct','mid layer',l,recn)
        if(recn/=0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              int_state%tct(i,j,l)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
      enddo
      call halo_exch(int_state%tct,lm,2,2)
!-----------------------------------------------------------------------
!--tp
      int_state%tp=0.
      do l=1,lm
        call getrecn(recname,reclevtyp,reclev,nrec,'tp','mid layer',l,recn)
        if(recn/=0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              int_state%tp(i,j,l)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
      enddo
      call halo_exch(int_state%tp,lm,2,2)
!-----------------------------------------------------------------------
!-- up
      int_state%up=0.
      do l=1,lm
        call getrecn(recname,reclevtyp,reclev,nrec,'up','mid layer',l,recn)
        if(recn/=0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              int_state%up(i,j,l)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
      enddo
      call halo_exch(int_state%up,lm,2,2)
!-----------------------------------------------------------------------
!-- vp
      int_state%vp=0.
      do l=1,lm
        call getrecn(recname,reclevtyp,reclev,nrec,'vp','mid layer',l,recn)
        if(recn/=0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              int_state%vp(i,j,l)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
      enddo
      call halo_exch(int_state%vp,lm,2,2)
!-----------------------------------------------------------------------
!-- z
      int_state%z=0.
      do l=1,lm
        call getrecn(recname,reclevtyp,reclev,nrec,'z','mid layer',l,recn)
        if(recn/=0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              int_state%z(i,j,l)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
      enddo
      call halo_exch(int_state%z,lm,2,2)
!-----------------------------------------------------------------------
!
      DO K=1,LM
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%Told(I,J,K)=0.
        ENDDO
        ENDDO
        call getrecn(recname,reclevtyp,reclev,nrec,'told','mid layer',k,recn)
        if(recn>0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              int_state%Told(i,j,k)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
!
      ENDDO
!-----------------------------------------------------------------------
!
      DO K=1,LM
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%Tadj(I,J,K)=0.
        ENDDO
        ENDDO
        call getrecn(recname,reclevtyp,reclev,nrec,'tadj','mid layer',k,recn)
        if(recn>0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              int_state%Tadj(i,j,k)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
!
      ENDDO
!-----------------------------------------------------------------------
!
      DO K=1,LM
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%CLDFRA(I,J,K)=0.
        ENDDO
        ENDDO
        call getrecn(recname,reclevtyp,reclev,nrec,'cldfra','mid layer',k,recn)
        if(recn>0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              int_state%CLDFRA(i,j,k)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
!
      ENDDO
!-----------------------------------------------------------------------
!
      DO K=1,LM
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%Q2(I,J,K)=0.
        ENDDO
        ENDDO
        call getrecn(recname,reclevtyp,reclev,nrec,'q2','mid layer',k,recn)
        if(recn>0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              int_state%Q2(i,j,k)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif

      ENDDO
      CALL HALO_EXCH(int_state%Q2,LM,2,2)
!-----------------------------------------------------------------------
!
      DO K=1,LM
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%RLWTT(I,J,K)=0.
        ENDDO
        ENDDO
        call getrecn(recname,reclevtyp,reclev,nrec,'rlwtt','mid layer',k,recn)
        if(recn>0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              int_state%RLWTT(i,j,k)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
!
      ENDDO
!-----------------------------------------------------------------------
!
      DO K=1,LM
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%RSWTT(I,J,K)=0.
        ENDDO
        ENDDO
        call getrecn(recname,reclevtyp,reclev,nrec,'rswtt','mid layer',k,recn)
        if(recn>0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              int_state%RSWTT(i,j,k)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
!
      ENDDO
!-----------------------------------------------------------------------
!
      DO K=1,LM
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%T(I,J,K)=0.
        ENDDO
        ENDDO
        call getrecn(recname,reclevtyp,reclev,nrec,'tmp','mid layer',k,recn)
        if(recn>0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              int_state%T(i,j,k)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
!
      ENDDO
      CALL HALO_EXCH(int_state%T,LM,2,2)
!-----------------------------------------------------------------------
!
      DO K=1,LM
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%TCUCN(I,J,K)=0.
        ENDDO
        ENDDO
        call getrecn(recname,reclevtyp,reclev,nrec,'tcucn','mid layer',k,recn)
        if(recn>0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              int_state%TCUCN(i,j,k)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
!
      ENDDO
!-----------------------------------------------------------------------
!
      DO K=1,LM
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%TRAIN(I,J,K)=0.
        ENDDO
        ENDDO
        call getrecn(recname,reclevtyp,reclev,nrec,'train','mid layer',k,recn)
        if(recn>0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              int_state%TRAIN(i,j,k)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
!
      ENDDO
!-----------------------------------------------------------------------
!
      DO K=1,LM
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%U(I,J,K)=0.
        ENDDO
        ENDDO
        call getrecn(recname,reclevtyp,reclev,nrec,'ugrd','mid layer',k,recn)
        if(recn>0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              int_state%U(i,j,k)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
!
      ENDDO
      CALL HALO_EXCH(int_state%U,LM,2,2)
!-----------------------------------------------------------------------
!
      DO K=1,LM
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%V(I,J,K)=0.
        ENDDO
        ENDDO
        call getrecn(recname,reclevtyp,reclev,nrec,'vgrd','mid layer',k,recn)
        if(recn>0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              int_state%V(i,j,k)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
!
      ENDDO
      CALL HALO_EXCH(int_state%V,LM,2,2)
!-----------------------------------------------------------------------
!
      DO K=1,LM
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%XLEN_MIX(I,J,K)=0.
        ENDDO
        ENDDO
        call getrecn(recname,reclevtyp,reclev,nrec,'xlen_mix','mid layer',k,recn)
        if(recn>0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              int_state%XLEN_MIX(i,j,k)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
!
      ENDDO
!-----------------------------------------------------------------------
!
      DO K=1,LM
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%F_ICE(I,J,K)=0.
        ENDDO
        ENDDO
        call getrecn(recname,reclevtyp,reclev,nrec,'f_ice','mid layer',k,recn)
        if(recn>0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              int_state%F_ICE(i,j,k)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
!
      ENDDO
!-----------------------------------------------------------------------
!
      DO K=1,LM
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%F_RIMEF(I,J,K)=0.
        ENDDO
        ENDDO
        call getrecn(recname,reclevtyp,reclev,nrec,'f_rimef','mid layer',k,recn)
        if(recn>0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              int_state%F_RIMEF(i,j,k)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
!
      ENDDO
!-----------------------------------------------------------------------
!
      DO K=1,LM
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%F_RAIN(I,J,K)=0.
        ENDDO
        ENDDO
        call getrecn(recname,reclevtyp,reclev,nrec,'f_rain','mid layer',k,recn)
        if(recn>0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              int_state%F_RAIN(i,j,k)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
!
      ENDDO
!-----------------------------------------------------------------------
!
      DO K=1,LM
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%REFL_10CM(I,J,K)=DBZmin
        ENDDO
        ENDDO
        call getrecn(recname,reclevtyp,reclev,nrec,'refl_10cm','mid layer',k,recn)
        if(recn>0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              int_state%REFL_10CM(i,j,k)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
!
      ENDDO
!-----------------------------------------------------------------------
!***  Radar-derived T tendencies from GSI analysis
!-----------------------------------------------------------------------
!
      if(int_state%RADAR_INIT==1) then
      DO K=1,LM
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%DFI_TTEN(I,J,K)=0.
        ENDDO
        ENDDO
        call getrecn(recname,reclevtyp,reclev,nrec,'dfi_tten','mid layer',k,recn)
        if(recn>0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              if(tmp(i-its+1+js+fldst)>0.0 .and. tmp(i-its+1+js+fldst)<0.1)then
                int_state%DFI_TTEN(i,j,k)=tmp(i-its+1+js+fldst)*1.00
              end if
            enddo
          enddo
        endif
!
      ENDDO
      end if
!
!-----------------------------------------------------------------------
!***  SH2O, SMC, STC
!-----------------------------------------------------------------------
!
      DO K=1,int_state%NSOIL
!
        call getrecn(recname,reclevtyp,reclev,nrec,'sh2o','soil layer',k,recn)
        if(recn>0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              int_state%SH2O(i,j,k)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif

        call getrecn(recname,reclevtyp,reclev,nrec,'smc','soil layer',k,recn)
        if(recn>0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              int_state%SMC(i,j,k)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
!
        call getrecn(recname,reclevtyp,reclev,nrec,'stc','soil layer',k,recn)
        if(recn>0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              int_state%STC(i,j,k)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif

      ENDDO


!-- tracers_prev
      do n=1,int_state%num_tracers_total
        int_state%tracers_prev(:,:,:,n)=0.
        write(tn,'(I3.3)')n
        do l=1,lm
          call getrecn(recname,reclevtyp,reclev,nrec,'tracers_prev_'//tn,  &
               'mid layer',l,recn)
          if(recn/=0) then
            fldst=(recn-1)*fldsize
            do j=jts,jte
              js=(j-jts)*(ite-its+1)
              do i=its,ite
                int_state%tracers_prev(i,j,l,n)=tmp(i-its+1+js+fldst)
              enddo
            enddo
          endif
        enddo
      enddo
      call halo_exch(int_state%tracers_prev,lm,int_state%num_tracers_total,1,2,2)
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Read from restart file: Real 2D arrays (contd.)
!-----------------------------------------------------------------------
!
!-- sice
      int_state%sice=0.
      call getrecn(recname,reclevtyp,reclev,nrec,'sice','sfc',1,recn)
      if(recn/=0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%sice(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
      call halo_exch(int_state%sice,1,2,2)
!-----------------------------------------------------------------------
!-- sm
      int_state%sm=0.
      call getrecn(recname,reclevtyp,reclev,nrec,'sm','sfc',1,recn)
      if(recn/=0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%sm(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
      call halo_exch(int_state%sm,1,2,2)
!
!-----------------------------------------------------------------------
!***  Read from restart file: Real 3D arrays
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!-- cw
      int_state%cw=0.
      do l=1,lm
        call getrecn(recname,reclevtyp,reclev,nrec,'clwmr','mid layer',l,recn)
        if(recn/=0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              int_state%cw(i,j,l)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
      enddo
      call halo_exch(int_state%cw,lm,2,2)
!      write(0,*)'in init restart after2,clwmr =',maxval(cw),minval(cw)
!-----------------------------------------------------------------------
!-- spfh
      int_state%q=0.
      do l=1,lm
        call getrecn(recname,reclevtyp,reclev,nrec,'spfh','mid layer',l,recn)
        if(recn/=0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              int_state%q(i,j,l)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
      enddo
      call halo_exch(int_state%q,lm,2,2)
!-----------------------------------------------------------------------
!-- q2
      int_state%q2=0.
      do l=1,lm
        call getrecn(recname,reclevtyp,reclev,nrec,'q2','mid layer',l,recn)
        if(recn/=0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              int_state%q2(i,j,l)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
      enddo
      call halo_exch(int_state%q2,lm,2,2)
!-----------------------------------------------------------------------
!-- t
      int_state%t=0.
      do l=1,lm
        call getrecn(recname,reclevtyp,reclev,nrec,'tmp','mid layer',l,recn)
        if(recn/=0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              int_state%t(i,j,l)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
      enddo
      call halo_exch(int_state%t,lm,2,2)
!-----------------------------------------------------------------------
!-- u
      int_state%u=0.
      do l=1,lm
        call getrecn(recname,reclevtyp,reclev,nrec,'ugrd','mid layer',l,recn)
        if(recn/=0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              int_state%u(i,j,l)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
      enddo
      call halo_exch(int_state%u,lm,2,2)
!-----------------------------------------------------------------------
!-- v
      int_state%v=0.
      do l=1,lm
        call getrecn(recname,reclevtyp,reclev,nrec,'vgrd','mid layer',l,recn)
        if(recn/=0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              int_state%v(i,j,l)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
      enddo
      call halo_exch(int_state%v,lm,2,2)
!-----------------------------------------------------------------------
!-- rest of tracers
      int_state%tracers(:,:,:,int_state%indx_q2+1:int_state%num_tracers_total)=0.
      do n=int_state%indx_q2+1,int_state%num_tracers_total                     !<-- The first 'indx_q2' arrays are unallocated pointers
        write(tn,'(I3.3)')n
        do l=1,lm
          call getrecn(recname,reclevtyp,reclev,nrec,'tracers_'//tn,     &
                       'mid layer',l,recn)
          if(recn/=0) then
            fldst=(recn-1)*fldsize
            do j=jts,jte
              js=(j-jts)*(ite-its+1)
              do i=its,ite
                int_state%tracers(i,j,l,n)=tmp(i-its+1+js+fldst)
              enddo
            enddo
          endif
        enddo
      enddo
      call halo_exch(int_state%tracers,lm,int_state%num_tracers_total,1,2,2)
!-----------------------------------------------------------------------
!
        tend_max=real(int_state%ihrend)
        ntstm_max=nint(tend_max*3600./int_state%dt)+1
        tend=real(int_state%nhours_fcst)
        int_state%ntstm=nint(tend*3600./int_state%dt)+1
        if(.not.int_state%global)then
          if(mype==0)then
            write(0,*)' Max runtime is ',tend_max,' hours'
          endif
        endif
        if(mype==0)then
          write(0,*)' Requested runtime is ',tend,' hours'
          write(0,*)' NTSTM=',int_state%ntstm
        endif
        if(int_state%ntstm>ntstm_max.and..not.int_state%global)then
          if(mype==0)then
            write(0,*)' Requested fcst length exceeds maximum'
            write(0,*)' Resetting to maximum'
          endif
          int_state%ntstm=min(int_state%ntstm,ntstm_max)
        endif
!
        ntsd=0
        int_state%ihr=nint(ntsd*int_state%dt/3600.)
!-----------------------------------------------------------------------
        do l=1,lm
          int_state%pdsg1(l)=int_state%dsg1(l)*int_state%pdtop
          int_state%psgml1(l)=int_state%sgml1(l)*int_state%pdtop+int_state%pt
        enddo
!
        do l=1,lm+1
          int_state%psg1(l)=int_state%sg1(l)*int_state%pdtop+int_state%pt
        enddo
!-----------------------------------------------------------------------
      endif  read_blocks                        ! cold start /restart
!-----------------------------------------------------------------------
!
      deallocate(tmp)
!
!-----------------------------------------------------------------------
!
      end subroutine read_nemsio
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE getrecn(recname,reclevtyp,reclev,nrec,fldname,          &
                         fldlevtyp,fldlev,recn)
!-----------------------------------------------------------------------
!-- this subroutine searches the field list to find out a specific field,
!-- and return the field number for that field
!-----------------------------------------------------------------------
!
        implicit none
!
        integer,intent(in)      :: nrec
        character(*),intent(in) :: recname(nrec)
        character(*),intent(in) :: reclevtyp(nrec)
        integer,intent(in)      :: reclev(nrec)
        character(*),intent(in) :: fldname
        character(*),intent(in) :: fldlevtyp
        integer,intent(in)      :: fldlev
        integer,intent(out)     :: recn
!
        integer i
!
        recn=0
        do i=1,nrec
          if(trim(recname(i))==trim(fldname).and.                        &
            trim(reclevtyp(i))==trim(fldlevtyp) .and.                    &
            reclev(i)==fldlev) then
            recn=i
            return
          endif
        enddo
!
        if(recn==0) print *,'WARNING: field ',trim(fldname),' ',         &
          trim(fldlevtyp),' ',fldlev,' is not in the nemsio file!'
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE getrecn
!
!-----------------------------------------------------------------------
!
      end module module_INIT_READ_NEMSIO
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
