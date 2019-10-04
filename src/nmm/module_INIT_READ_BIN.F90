!-----------------------------------------------------------------------
                        module module_INIT_READ_BIN
!-----------------------------------------------------------------------
use module_kinds
use module_dm_parallel,only : dstrb,idstrb
use module_exchange
use module_constants
use module_solver_internal_state,only: solver_internal_state
use module_microphysics_nmm
!
!-----------------------------------------------------------------------
!
      implicit none
!
      private
!
      public :: read_binary,physics_read_gwd
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
                        subroutine read_binary &
      (INT_STATE &
      ,my_domain_id &
      ,mpi_comm_comp &
      ,mype &
      ,its,ite,jts,jte &
      ,ims,ime,jms,jme &
      ,ids,ide,jds,jde &
      ,its_h2,ite_h2,jts_h2,jte_h2 &
      ,lm &
      ,rc)
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
 its,ite &
,ims,ime &
,ids,ide &
,its_h2,ite_h2 &
,jts,jte &
,jms,jme &
,jds,jde &
,jts_h2,jte_h2 &
,lm &
,mpi_comm_comp &
,my_domain_id &
,mype 

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
,irtn &
,j &                         ! index in y direction
,jend &
,k &                         ! index
,kount &
,ks &                        ! tracer index
,l &                         ! index in p direction
,lb &
,length &
,n &
,nl &
,nv &
,nvars_bc_2d_h &
,nvars_bc_3d_h &
,nvars_bc_4d_h &
,nvars_bc_2d_v &
,nvars_bc_3d_v &
,ub

integer(kind=kint) :: &      ! dimensions from input file
 im &
,jm &
,lmm &
,lnsh

integer(kind=kint) :: &
 iyear_fcst           &
,imonth_fcst          &
,iday_fcst            &
,ihour_fcst

real(kind=kfpt):: &
 tend,tend_max,arg

real(kind=kfpt),dimension(int_state%NSOIL) :: &
 soil1din

real(kind=kfpt),dimension(:),allocatable :: &
 all_bc_data,psint

real(kind=kfpt),allocatable,dimension(:,:) :: &
 temp1

integer(kind=kfpt),allocatable,dimension(:,:) :: &
 itemp

real(kind=kfpt),allocatable,dimension(:,:,:) :: &
 tempsoil

logical(kind=klog) :: opened

integer(kind=kint) :: &
 ierr 

character(64):: &
 infile

integer(kind=kint):: &
 ntsd &
,ntstm_max &
,nfcst

!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      allocate(temp1(ids:ide,jds:jde),stat=i)
!
!-----------------------------------------------------------------------
!
      select_unit: do n=51,59
        inquire(n,opened=opened)
        if(.not.opened)then
          nfcst=n
          exit select_unit
        endif
      enddo select_unit
!
!-----------------------------------------------------------------------
!
      read_blocks: if(.not.int_state%RESTART) then                     ! cold start
!
!-----------------------------------------------------------------------
!
        write(infile,'(a,i2.2)')'input_domain_',my_domain_id
        open(unit=nfcst,file=infile,status='old',form='unformatted'     &
            ,iostat=ierr)
        if(ierr/=0)then
          write(0,*)' Unable to open ',trim(infile),' in READ_BINARY'
          rc = ierr
          return
        endif
!
!-----------------------------------------------------------------------
!
        read (nfcst) int_state%RUN, &
                     int_state%IDAT, &
                     int_state%IHRST, &
                     int_state%IHREND, &
                     NTSD
        read (nfcst) int_state%PT, &
                     int_state%PDTOP, &
                     int_state%LPT2, &
                     int_state%SGM, &
                     int_state%SG1, &
                     int_state%DSG1, &
                     int_state%SGML1, &
                     int_state%SG2, &
                     int_state%DSG2, &
                     int_state%SGML2
        read (nfcst) int_state%I_PAR_STA, &
                     int_state%J_PAR_STA
        read (nfcst) int_state%DLMD, &
                     int_state%DPHD, &
                     int_state%WBD, &
                     int_state%SBD, &
                     int_state%TLM0D, &
                     int_state%TPH0D
        read (nfcst) im,jm,lmm,lnsh
!
!-----------------------------------------------------------------------
!***  Print the time & domain information.
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
          write(0,*) 'nmm_dyn reads of PT, PDTOP: ',int_state%PT,int_state%PDTOP
          write(0,*) 'nmm_dyn reads of I_PAR_STA, J_PAR_STA: ',int_state%I_PAR_STA,int_state%J_PAR_STA
          write(0,*) 'nmm_dyn reads of TLM0D, TPH0D: ',int_state%TLM0D,int_state%TPH0D
          write(0,*) 'nmm_dyn reads of DLMD, DPHD: ',int_state%DLMD,int_state%DPHD
          write(0,*) 'nmm_dyn reads of WBD, SBD: ',int_state%WBD,int_state%SBD
          write(0,*) 'nmm_dyn reads of IM, JM, LM, LNSH: ',im,jm,lmm,lnsh
        endif
!
!-----------------------------------------------------------------------
!
      DO L=1,LM+1
        int_state%PSG1(L)=int_state%SG1(L)*int_state%PDTOP+int_state%PT
      ENDDO
      DO L=1,LM
        int_state%PDSG1(L)=int_state%DSG1(L)*int_state%PDTOP
        int_state%PSGML1(L)=int_state%SGML1(L)*int_state%PDTOP+int_state%PT
      ENDDO
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
!       int_state%EPSQ2(K-1)=(1.+COS(ARG))*0.09+0.02
        int_state%EPSQ2(K-1)=0.02
        int_state%EPSL(K-1)=SQRT(int_state%EPSQ2(K-1)*0.5)
      ENDDO   
      int_state%EPSQ2(LM)=int_state%EPSQ2(LM-1)
      DEALLOCATE(PSINT)

!write(0,*) 'epsl=',epsl
!write(0,*) 'epsq2=',epsq2

!
!-----------------------------------------------------------------------
!***  Proceed with getting fields from input file.
!***  NOTE: Five records were already read at the top of this routine.
!-----------------------------------------------------------------------
!
        if(mype==0)then
          read(nfcst)temp1
        endif
        do j=jms,jme
        do i=ims,ime
          int_state%fis(i,j)=0.
        enddo
        enddo
        call dstrb(temp1,int_state%fis,1,1,1,1,1,mype,mpi_comm_comp)
        call halo_exch(int_state%fis,1,2,2)
!-----------------------------------------------------------------------
!
        if(mype==0)then
          read(nfcst)temp1
        endif
        do j=jms,jme
        do i=ims,ime
          int_state%stdh(i,j)=0.
        enddo
        enddo
        call dstrb(temp1,int_state%stdh,1,1,1,1,1,mype,mpi_comm_comp)
        call halo_exch(int_state%stdh,1,2,2)
!-----------------------------------------------------------------------
!
        if(mype==0)then
          read(nfcst)temp1
        endif
        do j=jms,jme
        do i=ims,ime
          int_state%sm(i,j)=0.
        enddo
        enddo
        call dstrb(temp1,int_state%sm,1,1,1,1,1,mype,mpi_comm_comp)
        call halo_exch(int_state%sm,1,2,2)
!-----------------------------------------------------------------------
!
        if(mype==0)then
          read(nfcst)temp1
        endif
        do j=jms,jme
        do i=ims,ime
          int_state%pd(i,j)=0.
        enddo
        enddo
        call dstrb(temp1,int_state%pd,1,1,1,1,1,mype,mpi_comm_comp)
        call halo_exch(int_state%pd,1,2,2)
!-----------------------------------------------------------------------
!
        call mpi_barrier(mpi_comm_comp,irtn)
        do l=1,lm
          if(mype==0)then
            read(nfcst)temp1
          endif
          do j=jms,jme
          do i=ims,ime
            int_state%u(i,j,l)=0.
          enddo
          enddo
          call dstrb(temp1,int_state%u,1,1,1,lm,l,mype,mpi_comm_comp)
        enddo
        call halo_exch(int_state%u,lm,2,2)
!-----------------------------------------------------------------------
!
        do l=1,lm
          if(mype==0)then
            read(nfcst)temp1
          endif
          do j=jms,jme
          do i=ims,ime
            int_state%v(i,j,l)=0.
          enddo
          enddo
          call dstrb(temp1,int_state%v,1,1,1,lm,l,mype,mpi_comm_comp)
        enddo
        call halo_exch(int_state%v,lm,2,2)
!-----------------------------------------------------------------------
!
        do l=1,lm
          if(mype==0)then
            read(nfcst)temp1
            write(0,*) 'L, T extremes: ', L, minval(temp1),maxval(temp1)
          endif
          do j=jms,jme
          do i=ims,ime
            int_state%t(i,j,l)=0.
          enddo
          enddo
          call dstrb(temp1,int_state%t,1,1,1,lm,l,mype,mpi_comm_comp)
        enddo
        call halo_exch(int_state%t,lm,2,2)
!-----------------------------------------------------------------------
!
        do l=1,lm
          if(mype==0)then
            read(nfcst)temp1
          endif
          do j=jms,jme
          do i=ims,ime
            int_state%q(i,j,l)=0.
          enddo
          enddo
          call dstrb(temp1,int_state%q,1,1,1,lm,l,mype,mpi_comm_comp)
        enddo
        call halo_exch(int_state%q,lm,2,2)
!
!-----------------------------------------------------------------------
!
        do l=1,lm
          if(mype==0)then
            read(nfcst)temp1
          endif
          do j=jms,jme
          do i=ims,ime
            int_state%cw(i,j,l)=0.
          enddo
          enddo
          call dstrb(temp1,int_state%cw,1,1,1,lm,l,mype,mpi_comm_comp)
        enddo
        call halo_exch(int_state%cw,lm,2,2)
!
!-----------------------------------------------------------------------
!
        do l=1,lm
          if(mype==0)then
            read(nfcst)temp1   ! O3
          endif
          do j=jms,jme
          do i=ims,ime
            int_state%o3(i,j,l)=0.
          enddo
          enddo
          call dstrb(temp1,int_state%o3,1,1,1,lm,l,mype,mpi_comm_comp)
        enddo
        call halo_exch(int_state%o3,lm,2,2)
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
              int_state%pint(i,j,l+1)=int_state%PINT(i,j,l)+int_state%DSG2(l)*int_state%PD(i,j)+int_state%PDSG1(l)
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
        if(mype==0)then
          read(nfcst)temp1  ! ALBEDO
        endif
      CALL DSTRB(TEMP1,int_state%ALBEDO,1,1,1,1,1,mype,mpi_comm_comp)
!
        if(mype==0)then
          read(nfcst)temp1  ! ALBASE
        endif
      CALL DSTRB(TEMP1,int_state%ALBASE,1,1,1,1,1,mype,mpi_comm_comp)
!
        if(mype==0)then
          read(nfcst)temp1  ! EPSR
        endif
      CALL DSTRB(TEMP1,int_state%EPSR,1,1,1,1,1,mype,mpi_comm_comp)
!
        if(mype==0)then
          read(nfcst)temp1  ! MXSNAL
        endif
      CALL DSTRB(TEMP1,int_state%MXSNAL,1,1,1,1,1,mype,mpi_comm_comp)
!
        if(mype==0)then
          read(nfcst)temp1  ! TSKIN
        endif
      CALL DSTRB(TEMP1,int_state%TSKIN,1,1,1,1,1,mype,mpi_comm_comp)
!
        if(mype==0)then
          read(nfcst)temp1  ! SST
        endif
      CALL DSTRB(TEMP1,int_state%SST,1,1,1,1,1,mype,mpi_comm_comp)
!
        if(mype==0)then
          read(nfcst)temp1  ! SNO and SNOWC
        endif
      CALL DSTRB(TEMP1,int_state%SNO,1,1,1,1,1,mype,mpi_comm_comp)

      DO J=JMS,JME
      DO I=IMS,IME
        if(int_state%SNO(I,J).gt.0.) then
           int_state%SNOWC(I,J)=0.98
        else
           int_state%SNOWC(I,J)=0.0
        endif
      ENDDO
      ENDDO
!
        if(mype==0)then
          read(nfcst)temp1  ! SI
        endif
      CALL DSTRB(TEMP1,int_state%SI,1,1,1,1,1,mype,mpi_comm_comp)
!
        if(mype==0)then
          read(nfcst)temp1  ! SICE
        endif
      CALL DSTRB(TEMP1,int_state%SICE,1,1,1,1,1,mype,mpi_comm_comp)
!
        if(mype==0)then
          read(nfcst) temp1  ! TG
        endif
      CALL DSTRB(TEMP1,int_state%TG,1,1,1,1,1,mype,mpi_comm_comp)
!
        if(mype==0)then
          read(nfcst) temp1  ! CMC
        endif
      CALL DSTRB(TEMP1,int_state%CMC,1,1,1,1,1,mype,mpi_comm_comp)
!
        if(mype==0)then
          read(nfcst) temp1  ! SR
        endif
      CALL DSTRB(TEMP1,int_state%SR,1,1,1,1,1,mype,mpi_comm_comp)
!
        if(mype==0)then
          read(nfcst) temp1  ! USTAR
        endif
      CALL DSTRB(TEMP1,int_state%USTAR,1,1,1,1,1,mype,mpi_comm_comp)
!
        if(mype==0)then
          read(nfcst) temp1  ! Z0
        endif
      CALL DSTRB(TEMP1,int_state%Z0,1,1,1,1,1,mype,mpi_comm_comp)
      CALL HALO_EXCH(int_state%Z0,1,3,3)
!
        if(mype==0)then
          read(nfcst) temp1  ! Z0BASE
        endif
      CALL DSTRB(TEMP1,int_state%Z0BASE,1,1,1,1,1,mype,mpi_comm_comp)
      CALL HALO_EXCH(int_state%Z0BASE,1,3,3)
!
      ALLOCATE(TEMPSOIL(1:int_state%NSOIL,IDS:IDE,JDS:JDE),STAT=I)
!
        if(mype==0)then
          read(nfcst) TEMPSOIL  ! STC
        endif
      CALL DSTRB(TEMPSOIL(1,IDS:IDE,JDS:JDE),int_state%STC(IMS:IME,JMS:JME,1),1,1,1,1,1,mype,mpi_comm_comp)
      CALL DSTRB(TEMPSOIL(2,IDS:IDE,JDS:JDE),int_state%STC(IMS:IME,JMS:JME,2),1,1,1,1,1,mype,mpi_comm_comp)
      CALL DSTRB(TEMPSOIL(3,IDS:IDE,JDS:JDE),int_state%STC(IMS:IME,JMS:JME,3),1,1,1,1,1,mype,mpi_comm_comp)
      CALL DSTRB(TEMPSOIL(4,IDS:IDE,JDS:JDE),int_state%STC(IMS:IME,JMS:JME,4),1,1,1,1,1,mype,mpi_comm_comp)
!
        if(mype==0)then
          read(nfcst) TEMPSOIL  ! SMC
        endif
      CALL DSTRB(TEMPSOIL(1,:,:),int_state%SMC(:,:,1),1,1,1,1,1,mype,mpi_comm_comp)
      CALL DSTRB(TEMPSOIL(2,:,:),int_state%SMC(:,:,2),1,1,1,1,1,mype,mpi_comm_comp)
      CALL DSTRB(TEMPSOIL(3,:,:),int_state%SMC(:,:,3),1,1,1,1,1,mype,mpi_comm_comp)
      CALL DSTRB(TEMPSOIL(4,:,:),int_state%SMC(:,:,4),1,1,1,1,1,mype,mpi_comm_comp)
!
        if(mype==0)then
          read(nfcst) TEMPSOIL  ! SH2O
        endif
      CALL DSTRB(TEMPSOIL(1,:,:),int_state%SH2O(:,:,1),1,1,1,1,1,mype,mpi_comm_comp)
      CALL DSTRB(TEMPSOIL(2,:,:),int_state%SH2O(:,:,2),1,1,1,1,1,mype,mpi_comm_comp)
      CALL DSTRB(TEMPSOIL(3,:,:),int_state%SH2O(:,:,3),1,1,1,1,1,mype,mpi_comm_comp)
      CALL DSTRB(TEMPSOIL(4,:,:),int_state%SH2O(:,:,4),1,1,1,1,1,mype,mpi_comm_comp)
!
      DEALLOCATE(TEMPSOIL)
      ALLOCATE(ITEMP(IDS:IDE,JDS:JDE),STAT=I)
!
        if(mype==0)then
          read(nfcst) ITEMP  ! ISLTYP
        endif
      CALL IDSTRB(ITEMP,int_state%ISLTYP,mype,mpi_comm_comp)
!
        if(mype==0)then
          read(nfcst) ITEMP  ! IVGTYP
        endif
      CALL IDSTRB(ITEMP,int_state%IVGTYP,mype,mpi_comm_comp)
!
      DEALLOCATE(ITEMP)
!
        if(mype==0)then
          read(nfcst) temp1  ! VEGFRC
        endif
      CALL DSTRB(TEMP1,int_state%VEGFRC,1,1,1,1,1,mype,mpi_comm_comp)
!
        if(mype==0)then
          read(nfcst) SOIL1DIN  ! DZSOIL
        endif
!
        if(mype==0)then
          read(nfcst) SOIL1DIN  ! SLDPTH
        endif
!
!       if(mype==0)then               ! here will be 14 orography fields for GWD
!         do n=1,14
!           read(nfcst) temp1
!         enddo
!       endif
!
!-----------------------------------------------------------------------
!
        close(nfcst)
!
!-----------------------------------------------------------------------
      else  read_blocks                         ! restart
!-----------------------------------------------------------------------
!
        write(infile,'(a,i2.2)')'restart_file_',my_domain_id
        open(unit=nfcst,file=infile,status='old',form='unformatted'     &
            ,iostat=ierr)
        if(ierr/=0)then
          write(0,*)' Unable to open ',trim(infile),' in READ_BINARY'
          rc = ierr
          return
        endif
!
!-----------------------------------------------------------------------
!***  Read from restart file: Integer scalars
!-----------------------------------------------------------------------
        read(nfcst) iyear_fcst
        read(nfcst) imonth_fcst
        read(nfcst) iday_fcst
        read(nfcst) ihour_fcst
        read(nfcst) !iminute_fcst
        read(nfcst) ! second_fcst
        read(nfcst) ! ntsd
        read(nfcst) ! im
        read(nfcst) ! jm
        read(nfcst) ! lm
        read(nfcst) int_state%IHRST
        read(nfcst) int_state%I_PAR_STA
        read(nfcst) int_state%J_PAR_STA
        read(nfcst) int_state%LAST_STEP_MOVED
        read(nfcst) int_state%LPT2
        read(nfcst) ! nsoil
        read(nfcst) ! nphs
        read(nfcst) ! nclod
        read(nfcst) ! nheat
        read(nfcst) ! nmts
        read(nfcst) ! nprec
        read(nfcst) ! nrdlw
        read(nfcst) ! nrdsw
        read(nfcst) ! nsrfc
!
        read(nfcst) ! cu_physics
        read(nfcst) ! mp_physics
        read(nfcst) ! lsm_physics
!
        read(nfcst) int_state%NTRACK
        read(nfcst) int_state%TRACK_HAVE_GUESS
        read(nfcst) int_state%TRACK_N_OLD
        read(nfcst) int_state%TRACKER_HAVEFIX
        read(nfcst) int_state%TRACKER_GAVE_UP
        read(nfcst) int_state%TRACKER_IFIX
        read(nfcst) int_state%TRACKER_JFIX
!
!-----------------------------------------------------------------------
!***  Read from restart file: Integer 1D arrays
!-----------------------------------------------------------------------
        read(nfcst) ! int_state%NTSCM
        read(nfcst) int_state%IDAT
!
        if(mype==0)then
          write(0,*)'**** read in core *******************'
          write(0,*)' Restart year =',iyear_fcst
          write(0,*)' Restart month=',imonth_fcst
          write(0,*)' Restart day  =',iday_fcst
          write(0,*)' Restart hour =',ihour_fcst
          write(0,*)' Original start year =',int_state%IDAT(3)
          write(0,*)' Original start month=',int_state%IDAT(2)
          write(0,*)' Original start day  =',int_state%IDAT(1)
          write(0,*)' Original start hour =',int_state%IHRST
          write(0,*)' Timestep   =',int_state%DT
          write(0,*)' Steps/hour =',3600./int_state%DT
          write(0,*)'*************************************'
        endif
!
        read(nfcst) int_state%TRACK_OLD_NTSD
!
!-----------------------------------------------------------------------
!***  Read from restart file: Real scalars
!-----------------------------------------------------------------------
        read(nfcst) ! dt
        read(nfcst) ! dyh
        read(nfcst) int_state%PDTOP
        read(nfcst) int_state%PT
        read(nfcst) int_state%TLM0D
        read(nfcst) int_state%TPH0D
        read(nfcst) ! tstart
        read(nfcst) int_state%DPHD
        read(nfcst) int_state%DLMD
        read(nfcst) int_state%SBD
        read(nfcst) int_state%WBD
!
        read(nfcst) int_state%TRACK_LAST_HOUR
        read(nfcst) int_state%TRACK_GUESS_LAT
        read(nfcst) int_state%TRACK_GUESS_LON
        read(nfcst) int_state%TRACK_EDGE_DIST
        read(nfcst) int_state%TRACK_STDERR_M1
        read(nfcst) int_state%TRACK_STDERR_M2
        read(nfcst) int_state%TRACK_STDERR_M3
        read(nfcst) int_state%TRACKER_FIXLAT
        read(nfcst) int_state%TRACKER_FIXLON
        read(nfcst) int_state%TRACKER_RMW
        read(nfcst) int_state%TRACKER_PMIN
        read(nfcst) int_state%TRACKER_VMAX
!-----------------------------------------------------------------------
!***  Read from restart file: Real 1D arrays
!-----------------------------------------------------------------------
        read(nfcst) ! dxh
        read(nfcst) int_state%SG1
        read(nfcst) int_state%SG2
        read(nfcst) int_state%DSG1
        read(nfcst) int_state%DSG2
        read(nfcst) int_state%SGML1
        read(nfcst) int_state%SGML2
        read(nfcst) int_state%SGM
        read(nfcst) int_state%EPSL
        read(nfcst) int_state%EPSQ2
        read(nfcst) int_state%SLDPTH
        read(nfcst) int_state%MP_RESTART_STATE
        read(nfcst) int_state%TBPVS_STATE
        read(nfcst) int_state%TBPVS0_STATE
!
        read(nfcst) int_state%TRACK_OLD_LAT
        read(nfcst) int_state%TRACK_OLD_LON
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
!-----------------------------------------------------------------------
!***  Read in the full-domain 1-D datastring of boundary winds.
!***  Each task isolates its own piece of that data.
!-----------------------------------------------------------------------
!
        length=(int_state%nlev_h*int_state%lnsh                         &  !<-- Total # of words
               +int_state%nlev_v*int_state%lnsv)                        &  !    in full-domain
               *2*2*(ide-ids+jde-jds+2)                                    !    boundary arrays.
        allocate(all_bc_data(1:length))
!
        read(nfcst) all_bc_data
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
        read(nfcst) ! global
        read(nfcst) int_state%RUN
        read(nfcst) ! adiabatic
!
!-----------------------------------------------------------------------
!***  Read from restart file: Integer 2D arrays
!-----------------------------------------------------------------------
!
      ALLOCATE(ITEMP(IDS:IDE,JDS:JDE),STAT=I)
!
      IF(MYPE==0)THEN
        READ(NFCST) ITEMP
      ENDIF
      CALL IDSTRB(ITEMP,int_state%ISLTYP,MYPE,MPI_COMM_COMP)
!
      IF(MYPE==0)THEN
        READ(NFCST) ITEMP
      ENDIF
      CALL IDSTRB(ITEMP,int_state%IVGTYP,MYPE,MPI_COMM_COMP)
!
      IF(MYPE==0)THEN
        READ(NFCST) ITEMP
      ENDIF
      CALL IDSTRB(ITEMP,int_state%NCFRCV,MYPE,MPI_COMM_COMP)
!
      IF(MYPE==0)THEN
        READ(NFCST) ITEMP
      ENDIF
      CALL IDSTRB(ITEMP,int_state%NCFRST,MYPE,MPI_COMM_COMP)
!
      IF(MYPE==0)THEN
        READ(NFCST) ITEMP
      ENDIF
      CALL IDSTRB(ITEMP,int_state%TRACKER_FIXES,MYPE,MPI_COMM_COMP)
!
      DEALLOCATE(ITEMP)
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  Read from restart file by alphabetical order (new in ESMF6)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  ACFRCV
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%ACFRCV,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  ACFRST
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%ACFRST,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  ACPREC
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%ACPREC,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  ACPREC_TOT
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%ACPREC_TOT,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  ACSNOM
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%ACSNOM,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  ACSNOW
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%ACSNOW,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  ACUTIM
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%ACUTIM,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  AKHS
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%AKHS,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  AKHS_OUT
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%AKHS_OUT,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  AKMS
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%AKMS,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  AKMS_OUT
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%AKMS_OUT,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  ALBASE
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%ALBASE,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  ALBEDO
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%ALBEDO,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  ALWIN
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%ALWIN,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  ALWOUT
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%ALWOUT,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  ALWTOA
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%ALWTOA,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  APHTIM
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%APHTIM,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  ARDLW
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%ARDLW,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  ARDSW
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%ARDSW,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  ASRFC
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%ASRFC,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  ASWIN
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%ASWIN,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  ASWOUT
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%ASWOUT,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  ASWTOA
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%ASWTOA,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  AVCNVC
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%AVCNVC,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  AVRAIN
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%AVRAIN,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  BGROFF
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%BGROFF,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  CFRACH
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%CFRACH,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  CFRACL
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%CFRACL,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  CFRACM
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%CFRACM,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  CLDEFI
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%CLDEFI,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  CLDFRA
!-----------------------------------------------------------------------
      DO L=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1
        ENDIF
        CALL DSTRB(TEMP1,INT_STATE%CLDFRA,1,1,1,LM,L,MYPE,MPI_COMM_COMP)
      ENDDO
!-----------------------------------------------------------------------
!***  CMC
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%CMC,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  CNVBOT
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%CNVBOT,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  CNVTOP
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%CNVTOP,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  CPRATE
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%CPRATE,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  CUPPT
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%CUPPT,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  CUPREC
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%CUPREC,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  CW
!-----------------------------------------------------------------------
      DO L=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1
        ENDIF
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%CW(I,J,L)=0.
        ENDDO
        ENDDO
        CALL DSTRB(TEMP1,int_state%CW,1,1,1,LM,L,MYPE,MPI_COMM_COMP)
      ENDDO
      CALL HALO_EXCH(int_state%CW,LM,2,2)
!-----------------------------------------------------------------------
!***  CZEN
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%CZEN,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  CZMEAN
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%CZMEAN,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  DFI_TTEN
!-----------------------------------------------------------------------
      DO L=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1
        ENDIF
        CALL DSTRB(TEMP1,int_state%DFI_TTEN,1,1,1,LM,L,MYPE,MPI_COMM_COMP)
      ENDDO
      CALL HALO_EXCH(int_state%DFI_TTEN,LM,2,2)
!-----------------------------------------------------------------------
!***  DIV
!-----------------------------------------------------------------------
      DO L=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1
        ENDIF
        CALL DSTRB(TEMP1,int_state%DIV,1,1,1,LM,L,MYPE,MPI_COMM_COMP)
      ENDDO
      CALL HALO_EXCH(int_state%DIV,LM,2,2)
!-----------------------------------------------------------------------
!***  DWDT
!-----------------------------------------------------------------------
      DO L=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1
        ENDIF
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%DWDT(I,J,L)=0.
        ENDDO
        ENDDO
        CALL DSTRB(TEMP1,int_state%DWDT,1,1,1,LM,L,MYPE,MPI_COMM_COMP)
      ENDDO
      CALL HALO_EXCH(int_state%DWDT,LM,2,2)
!-----------------------------------------------------------------------
!***  EPSR
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%EPSR,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  F
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%F(I,J)=0.
      ENDDO
      ENDDO
      CALL DSTRB(TEMP1,int_state%F,1,1,1,1,1,MYPE,MPI_COMM_COMP)
      CALL HALO_EXCH(int_state%F,1,2,2)
!-----------------------------------------------------------------------
!***  FIS
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%FIS(I,J)=0.
      ENDDO
      ENDDO
      CALL DSTRB(TEMP1,int_state%FIS,1,1,1,1,1,MYPE,MPI_COMM_COMP)
      CALL HALO_EXCH(int_state%FIS,1,2,2)
!-----------------------------------------------------------------------
!***  F_ICE
!-----------------------------------------------------------------------
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1
        ENDIF
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%F_ICE(I,J,K)=0.
        ENDDO
        ENDDO
        CALL DSTRB(TEMP1,int_state%F_ICE,1,1,1,LM,K,MYPE,MPI_COMM_COMP)
      ENDDO
!-----------------------------------------------------------------------
!***  F_RAIN
!-----------------------------------------------------------------------
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1
        ENDIF
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%F_RAIN(I,J,K)=0.
        ENDDO
        ENDDO
        CALL DSTRB(TEMP1,int_state%F_RAIN,1,1,1,LM,K,MYPE,MPI_COMM_COMP)
      ENDDO
!-----------------------------------------------------------------------
!***  F_RIMEF
!-----------------------------------------------------------------------
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1
        ENDIF
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%F_RIMEF(I,J,K)=0.
        ENDDO
        ENDDO
        CALL DSTRB(TEMP1,int_state%F_RIMEF,1,1,1,LM,K,MYPE,MPI_COMM_COMP)
      ENDDO
!-----------------------------------------------------------------------
!***  GLAT
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%GLAT(I,J)=0.
      ENDDO
      ENDDO
      CALL DSTRB(TEMP1,int_state%GLAT,1,1,1,1,1,MYPE,MPI_COMM_COMP)
      CALL HALO_EXCH(int_state%GLAT,1,2,2)
!-----------------------------------------------------------------------
!***  GLON
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%GLON(I,J)=0.
      ENDDO
      ENDDO
      CALL DSTRB(TEMP1,int_state%GLON,1,1,1,1,1,MYPE,MPI_COMM_COMP)
      CALL HALO_EXCH(int_state%GLON,1,2,2)
!-----------------------------------------------------------------------
!***  GRNFLX
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%GRNFLX,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  HBOT
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%HBOT,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  HBOTD
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%HBOTD,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  HBOTS
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%HBOTS,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  HDACVX
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%HDACVX(I,J)=0.
      ENDDO
      ENDDO
      CALL DSTRB(TEMP1,int_state%HDACVX,1,1,1,1,1,MYPE,MPI_COMM_COMP)
      CALL HALO_EXCH(int_state%HDACVX,1,2,2)
!-----------------------------------------------------------------------
!***  HDACVY
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%HDACVY(I,J)=0.
      ENDDO
      ENDDO
      CALL DSTRB(TEMP1,int_state%HDACVY,1,1,1,1,1,MYPE,MPI_COMM_COMP)
      CALL HALO_EXCH(int_state%HDACVY,1,2,2)
!-----------------------------------------------------------------------
!***  HDACX
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%HDACX(I,J)=0.
      ENDDO
      ENDDO
      CALL DSTRB(TEMP1,int_state%HDACX,1,1,1,1,1,MYPE,MPI_COMM_COMP)
      CALL HALO_EXCH(int_state%HDACX,1,2,2)
!-----------------------------------------------------------------------
!***  HDACY
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%HDACY(I,J)=0.
      ENDDO
      ENDDO
      CALL DSTRB(TEMP1,int_state%HDACY,1,1,1,1,1,MYPE,MPI_COMM_COMP)
      CALL HALO_EXCH(int_state%HDACY,1,2,2)
!-----------------------------------------------------------------------
!***  HTOP
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%HTOP,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  HTOPD
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%HTOPD,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  HTOPS
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%HTOPS,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  M10RV
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%M10RV,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  M10WIND
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%M10WIND,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  MEMBRANE_MSLP
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%MEMBRANE_MSLP,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  MXSNAL (SNOW ALBEDO)
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%MXSNAL,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  O3
!-----------------------------------------------------------------------
      DO L=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1
        ENDIF
        CALL DSTRB(TEMP1,int_state%O3,1,1,1,LM,L,MYPE,MPI_COMM_COMP)
      ENDDO
      CALL HALO_EXCH(int_state%O3,LM,2,2)
!-----------------------------------------------------------------------
!***  OMGALF
!-----------------------------------------------------------------------
      DO L=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1
        ENDIF
        CALL DSTRB(TEMP1,int_state%OMGALF,1,1,1,LM,L,MYPE,MPI_COMM_COMP)
      ENDDO
      CALL HALO_EXCH(int_state%OMGALF,LM,2,2)
!-----------------------------------------------------------------------
!***  P500U
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%P500U,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  P500V
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%P500V,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  P700RV
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%P700RV,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  P700U
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%P700U,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  P700V
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%P700V,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  P700WIND
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%P700WIND,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  P700Z
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%P700Z,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  P850RV
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%P850RV,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  P850U
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%P850U,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  P850V
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%P850V,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  P850WIND
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%P850WIND,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  P850Z
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%P850Z,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  PBLH
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%PBLH,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  PD
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%PD(I,J)=0.
      ENDDO
      ENDDO
      CALL DSTRB(TEMP1,int_state%PD,1,1,1,1,1,MYPE,MPI_COMM_COMP)
      CALL HALO_EXCH(int_state%PD,1,2,2)
!-----------------------------------------------------------------------
!***  PDO
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%PDO(I,J)=0.
      ENDDO
      ENDDO
      CALL DSTRB(TEMP1,int_state%PDO,1,1,1,1,1,MYPE,MPI_COMM_COMP)
      CALL HALO_EXCH(int_state%PDO,1,2,2)
!-----------------------------------------------------------------------
!***  PINT
!-----------------------------------------------------------------------
      DO L=1,LM+1
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1
        ENDIF
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%PINT(I,J,L)=0.
        ENDDO
        ENDDO
        CALL DSTRB(TEMP1,int_state%PINT,1,1,1,LM+1,L,MYPE,MPI_COMM_COMP)
      ENDDO
      CALL HALO_EXCH(int_state%PINT,LM+1,2,2)
!-----------------------------------------------------------------------
!***  POTEVP
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%POTEVP,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  POTFLX
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%POTFLX,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  PREC
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%PREC,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  PSGDT
!-----------------------------------------------------------------------
      DO L=1,LM-1
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1
        ENDIF
      ENDDO
!-----------------------------------------------------------------------
!***  PSHLTR
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%PSHLTR,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  Q10
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%Q10,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  Q2
!-----------------------------------------------------------------------
      DO L=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1
        ENDIF
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%Q2(I,J,L)=0.
        ENDDO
        ENDDO
        CALL DSTRB(TEMP1,int_state%Q2,1,1,1,LM,L,MYPE,MPI_COMM_COMP)
      ENDDO
      CALL HALO_EXCH(int_state%Q2,LM,2,2)
!-----------------------------------------------------------------------
!***  QSH
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%QSH,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  QSHLTR
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%QSHLTR,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  QWBS
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%QWBS,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  QZ0
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%QZ0,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  Q
!-----------------------------------------------------------------------
      DO L=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1
        ENDIF
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%Q(I,J,L)=0.
        ENDDO
        ENDDO
        CALL DSTRB(TEMP1,int_state%Q,1,1,1,LM,L,MYPE,MPI_COMM_COMP)
      ENDDO
      CALL HALO_EXCH(int_state%Q,LM,2,2)
!-----------------------------------------------------------------------
!***  RADOT
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%RADOT,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  RLWIN
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%RLWIN,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  RLWTOA
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%RLWTOA,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  RLWTT
!-----------------------------------------------------------------------
      DO L=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1
        ENDIF
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%RLWTT(I,J,L)=0.
        ENDDO
        ENDDO
        CALL DSTRB(TEMP1,int_state%RLWTT,1,1,1,LM,L,MYPE,MPI_COMM_COMP)
      ENDDO
!-----------------------------------------------------------------------
!***  RMOL
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%RMOL,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  RSWIN
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%RSWIN,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  RSWINC
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%RSWINC,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  RSWOUT
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%RSWOUT,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  RSWTOA
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%RSWTOA,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  RSWTT
!-----------------------------------------------------------------------
      DO L=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1
        ENDIF
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%RSWTT(I,J,L)=0.
        ENDDO
        ENDDO
        CALL DSTRB(TEMP1,int_state%RSWTT,1,1,1,LM,L,MYPE,MPI_COMM_COMP)
      ENDDO
!-----------------------------------------------------------------------
!***  SFCEVP
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SFCEVP,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  SFCEXC
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SFCEXC,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  SFCLHX
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SFCLHX,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  SFCSHX
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SFCSHX,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  SH2O
!-----------------------------------------------------------------------
      DO K=1,int_state%NSOIL
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1
        ENDIF
        CALL DSTRB(TEMP1,int_state%SH2O,1,1,1,int_state%NSOIL,K        &
                  ,MYPE,MPI_COMM_COMP)
      ENDDO
!-----------------------------------------------------------------------
!***  SI
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SI,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  SICE
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SICE,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  SIGT4
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SIGT4,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  SM (Seamask)
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%SM(I,J)=0.
      ENDDO
      ENDDO
      CALL DSTRB(TEMP1,int_state%SM,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  SM10RV
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SM10RV,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  SM10WIND
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SM10WIND,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  SMC
!-----------------------------------------------------------------------
      DO K=1,int_state%NSOIL
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1
        ENDIF
        CALL DSTRB(TEMP1,int_state%SMC,1,1,1,int_state%NSOIL,K        &
                  ,MYPE,MPI_COMM_COMP)
      ENDDO
!-----------------------------------------------------------------------
!***  SMSLP
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SMSLP,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  SMSTAV
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SMSTAV,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  SMSTOT
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SMSTOT,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  SNO
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SNO,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  SNOPCX
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SNOPCX,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  SNOWC
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SNOWC,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  SOILTB
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SOILTB,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  SP700RV
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SP700RV,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  SP700WIND
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SP700WIND,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  SP700Z
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SP700Z,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  SP850RV
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SP850RV,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  SP850WIND
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SP850WIND,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  SP850Z
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SP850Z,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  SR
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SR,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  SSROFF
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SSROFF,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  SST
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SST,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  STC
!-----------------------------------------------------------------------
      DO K=1,int_state%NSOIL
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1
        ENDIF
        CALL DSTRB(TEMP1,int_state%STC,1,1,1,int_state%NSOIL,K          &
                  ,MYPE,MPI_COMM_COMP)
      ENDDO
!-----------------------------------------------------------------------
!***  SUBSHX
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SUBSHX,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  T2
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%T2,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  TCT
!-----------------------------------------------------------------------
      DO L=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1
        ENDIF
        CALL DSTRB(TEMP1,int_state%TCT,1,1,1,LM,L,MYPE,MPI_COMM_COMP)
      ENDDO
      CALL HALO_EXCH(int_state%TCT,LM,2,2)
!-----------------------------------------------------------------------
!***  TCUCN
!-----------------------------------------------------------------------
      DO L=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1
        ENDIF
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%TCUCN(I,J,L)=0.
        ENDDO
        ENDDO
        CALL DSTRB(TEMP1,int_state%TCUCN,1,1,1,LM,L,MYPE,MPI_COMM_COMP)
      ENDDO
!-----------------------------------------------------------------------
!***  TCU
!-----------------------------------------------------------------------
      DO L=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1
        ENDIF
        CALL DSTRB(TEMP1,int_state%TCU,1,1,1,LM,L,MYPE,MPI_COMM_COMP)
      ENDDO
      CALL HALO_EXCH(int_state%TCU,LM,2,2)
!-----------------------------------------------------------------------
!***  TCV
!-----------------------------------------------------------------------
      DO L=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1
        ENDIF
        CALL DSTRB(TEMP1,int_state%tcv,1,1,1,LM,L,MYPE,MPI_COMM_COMP)
      ENDDO
      CALL HALO_EXCH(int_state%TCV,LM,2,2)
!-----------------------------------------------------------------------
!***  TG
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%TG,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  TH10
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%TH10,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  THS
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%THS,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  THZ0
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%THZ0,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  TLMAX
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%TLMAX,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  TLMIN
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%TLMIN,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  TP
!-----------------------------------------------------------------------
      DO L=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1
        ENDIF
        CALL DSTRB(TEMP1,int_state%TP,1,1,1,LM,L,MYPE,MPI_COMM_COMP)
      ENDDO
      CALL HALO_EXCH(int_state%TP,LM,2,2)
!-----------------------------------------------------------------------
!***  TRACERS
!-----------------------------------------------------------------------
      DO N=int_state%INDX_Q2+1,int_state%NUM_TRACERS_TOTAL  !<-- The first 'indx_q2' arrays are unallocated pointers
        DO L=1,LM
          IF(MYPE==0)THEN
            READ(NFCST)TEMP1
          ENDIF
          DO J=JMS,JME
          DO I=IMS,IME
            int_state%TRACERS(I,J,L,N)=0.
          ENDDO
          ENDDO
          CALL DSTRB(TEMP1,int_state%TRACERS(:,:,:,N),1,1,1,LM,L      &
                    ,MYPE,MPI_COMM_COMP)
        ENDDO
      ENDDO
      CALL HALO_EXCH(int_state%TRACERS,LM,int_state%NUM_TRACERS_TOTAL,1,2,2)
!-----------------------------------------------------------------------
!***  TRACERS_PREV
!-----------------------------------------------------------------------
      DO N=1,int_state%NUM_TRACERS_TOTAL
        DO L=1,LM
          IF(MYPE==0)THEN
            READ(NFCST)TEMP1
          ENDIF
          CALL DSTRB(TEMP1,int_state%TRACERS_PREV(:,:,:,N),1,1,1,LM,L &
                    ,MYPE,MPI_COMM_COMP)
        ENDDO
      ENDDO
      CALL HALO_EXCH(int_state%TRACERS_PREV,LM,int_state%NUM_TRACERS_TOTAL,1,2,2)
!-----------------------------------------------------------------------
!***  TRACKER_ANGLE
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%TRACKER_ANGLE,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  TRACKER_DISTSQ
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%TRACKER_DISTSQ,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  TRAIN
!-----------------------------------------------------------------------
      DO L=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1
        ENDIF
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%TRAIN(I,J,L)=0.
        ENDDO
        ENDDO
        CALL DSTRB(TEMP1,int_state%TRAIN,1,1,1,LM,L,MYPE,MPI_COMM_COMP)
      ENDDO
!-----------------------------------------------------------------------
!***  TSHLTR
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%TSHLTR,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  TSKIN
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%TSKIN,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  TWBS
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%TWBS,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  T
!-----------------------------------------------------------------------
      DO L=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1
          WRITE(0,*) 'L, T extremes: ', L, MINVAL(TEMP1),MAXVAL(TEMP1)
        ENDIF
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%T(I,J,L)=0.
        ENDDO
        ENDDO
        CALL DSTRB(TEMP1,int_state%T,1,1,1,LM,L,MYPE,MPI_COMM_COMP)
      ENDDO
      CALL HALO_EXCH(int_state%T,LM,2,2)
!-----------------------------------------------------------------------
!***  TADJ
!-----------------------------------------------------------------------
      DO L=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1
        ENDIF
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%TADJ(I,J,L)=0.
        ENDDO
        ENDDO
        CALL DSTRB(TEMP1,int_state%TADJ,1,1,1,LM,L,MYPE,MPI_COMM_COMP)
      ENDDO
      CALL HALO_EXCH(int_state%TADJ,LM,2,2)
!-----------------------------------------------------------------------
!***  TOLD
!-----------------------------------------------------------------------
      DO L=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1
        ENDIF
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%TOLD(I,J,L)=0.
        ENDDO
        ENDDO
        CALL DSTRB(TEMP1,int_state%TOLD,1,1,1,LM,L,MYPE,MPI_COMM_COMP)
      ENDDO
      CALL HALO_EXCH(int_state%TOLD,LM,2,2)
!-----------------------------------------------------------------------
!***  U10
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%U10,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  UP
!-----------------------------------------------------------------------
      DO L=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1
        ENDIF
        CALL DSTRB(TEMP1,int_state%UP,1,1,1,LM,L,MYPE,MPI_COMM_COMP)
      ENDDO
      CALL HALO_EXCH(int_state%UP,LM,2,2)
!-----------------------------------------------------------------------
!***  USTAR
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%USTAR,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  UZ0
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%UZ0,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  U
!-----------------------------------------------------------------------
      DO L=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1
        ENDIF
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%U(I,J,L)=0.
        ENDDO
        ENDDO
        CALL DSTRB(TEMP1,int_state%U,1,1,1,LM,L,MYPE,MPI_COMM_COMP)
      ENDDO
      CALL HALO_EXCH(int_state%U,LM,2,2)
!-----------------------------------------------------------------------
!***  V10
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%V10,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  VEGFRC
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST) TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%VEGFRC,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  VLAT
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%VLAT(I,J)=0.
      ENDDO
      ENDDO
      CALL DSTRB(TEMP1,int_state%VLAT,1,1,1,1,1,MYPE,MPI_COMM_COMP)
      CALL HALO_EXCH(int_state%VLAT,1,2,2)
!-----------------------------------------------------------------------
!***  VLON
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%VLON(I,J)=0.
      ENDDO
      ENDDO
      CALL DSTRB(TEMP1,int_state%VLON,1,1,1,1,1,MYPE,MPI_COMM_COMP)
      CALL HALO_EXCH(int_state%VLON,1,2,2)
!-----------------------------------------------------------------------
!***  VP
!-----------------------------------------------------------------------
      DO L=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1
        ENDIF
        CALL DSTRB(TEMP1,int_state%VP,1,1,1,LM,L,MYPE,MPI_COMM_COMP)
      ENDDO
      CALL HALO_EXCH(int_state%VP,LM,2,2)
!-----------------------------------------------------------------------
!***  VZ0
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST) TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%VZ0,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  V
!-----------------------------------------------------------------------
      DO L=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1
        ENDIF
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%V(I,J,L)=0.
        ENDDO
        ENDDO
        CALL DSTRB(TEMP1,int_state%V,1,1,1,LM,L,MYPE,MPI_COMM_COMP)
      ENDDO
      CALL HALO_EXCH(int_state%V,LM,2,2)
!-----------------------------------------------------------------------
!***  W
!-----------------------------------------------------------------------
      int_state%W=0.0
      DO L=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1
        ENDIF
        CALL DSTRB(TEMP1,int_state%W,1,1,1,LM,L,MYPE,MPI_COMM_COMP)
      ENDDO
      CALL HALO_EXCH(int_state%W,LM,2,2)
!-----------------------------------------------------------------------
!***  XLEN_MIX
!-----------------------------------------------------------------------
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1
        ENDIF
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%XLEN_MIX(I,J,K)=0.
        ENDDO
        ENDDO
        CALL DSTRB(TEMP1,int_state%XLEN_MIX,1,1,1,LM,K,MYPE,MPI_COMM_COMP)
      ENDDO
!-----------------------------------------------------------------------
!***  Z0
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%Z0,1,1,1,1,1,MYPE,MPI_COMM_COMP)
      CALL HALO_EXCH(int_state%Z0,1,3,3)
!-----------------------------------------------------------------------
!***  Z0BASE
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%Z0BASE,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!***  Z
!-----------------------------------------------------------------------
      DO L=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1
        ENDIF
        CALL DSTRB(TEMP1,int_state%Z,1,1,1,LM,L,MYPE,MPI_COMM_COMP)
      ENDDO
      CALL HALO_EXCH(int_state%Z,LM,2,2)
!-----------------------------------------------------------------------
!
        close(nfcst)
!
!-----------------------------------------------------------------------
!
        tend_max=real(int_state%IHREND)
        ntstm_max=nint(tend_max*3600./int_state%DT)+1
        tend=real(int_state%nhours_fcst)
        int_state%NTSTM=nint(tend*3600./int_state%DT)+1
        if(.not.int_state%GLOBAL)then
          if(mype==0)then
            write(0,*)' Max runtime is ',tend_max,' hours'
          endif
        endif
        if(mype==0)then
          write(0,*)' Requested runtime is ',tend,' hours'
          write(0,*)' NTSTM=',int_state%NTSTM
        endif
        if(int_state%NTSTM>ntstm_max.and..not.int_state%GLOBAL)then
          if(mype==0)then
            write(0,*)' Requested fcst length exceeds maximum'
            write(0,*)' Resetting to maximum'
          endif
          int_state%NTSTM=min(int_state%NTSTM,ntstm_max)
        endif
!
        ntsd=0
        int_state%IHR=nint(ntsd*int_state%DT/3600.)
!-----------------------------------------------------------------------
        do l=1,lm
          int_state%PDSG1(l)=int_state%DSG1(l)*int_state%PDTOP
          int_state%PSGML1(l)=int_state%SGML1(l)*int_state%PDTOP+int_state%PT
        enddo
!
        do l=1,lm+1
          int_state%PSG1(l)=int_state%SG1(l)*int_state%PDTOP+int_state%PT
        enddo
!-----------------------------------------------------------------------
      endif  read_blocks                        ! cold start /restart
!-----------------------------------------------------------------------
!
      deallocate(temp1)
!
!-----------------------------------------------------------------------
!
      end subroutine read_binary
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
      SUBROUTINE PHYSICS_READ_GWD(INFILE,NGWD,INT_STATE                &
                                 ,MYPE,MPI_COMM_COMP                   &
                                 ,IDS,IDE,JDS,JDE,RC)
!----------------------------------------------------------------------
!
!------------------------
!***  Argument variables
!------------------------
!
      INTEGER,INTENT(IN) :: NGWD,MYPE,MPI_COMM_COMP
      INTEGER,INTENT(IN) :: IDS,IDE,JDS,JDE
!
      CHARACTER(LEN=*),INTENT(IN) :: INFILE
!
      TYPE(SOLVER_INTERNAL_STATE),POINTER,INTENT(INOUT) :: INT_STATE     !<-- The physics internal state
!
      INTEGER,INTENT(OUT) :: RC
!
!---------------------
!***  Local variables
!---------------------
!
      INTEGER :: IERR
!
      REAL,DIMENSION(:,:),ALLOCATABLE :: TEMP_GWD
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC = 0
!
      ALLOCATE(TEMP_GWD(IDS:IDE,JDS:JDE))      
!
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        OPEN(unit=NGWD,file=INFILE,status='old',form='unformatted'      &
            ,iostat=IERR)
        IF(IERR/=0)THEN
          WRITE(0,*)' Unable to open file ',TRIM(INFILE)                &
                   ,' in PHYSICS_READ_GWD'
          RC = IERR
          RETURN
        ENDIF
      ENDIF
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NGWD)TEMP_GWD
      ENDIF
!
      CALL DSTRB(TEMP_GWD,int_state%HSTDV,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NGWD)TEMP_GWD
      ENDIF
!
      CALL DSTRB(TEMP_GWD,int_state%HCNVX,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NGWD)TEMP_GWD
      ENDIF
!
      CALL DSTRB(TEMP_GWD,int_state%HASYW,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NGWD)TEMP_GWD
      ENDIF
!
      CALL DSTRB(TEMP_GWD,int_state%HASYS,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NGWD)TEMP_GWD
      ENDIF
!
      CALL DSTRB(TEMP_GWD,int_state%HASYSW,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NGWD)TEMP_GWD
      ENDIF
!
      CALL DSTRB(TEMP_GWD,int_state%HASYNW,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NGWD)TEMP_GWD
      ENDIF
!
      CALL DSTRB(TEMP_GWD,int_state%HLENW,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NGWD)TEMP_GWD
      ENDIF
!
      CALL DSTRB(TEMP_GWD,int_state%HLENS,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NGWD)TEMP_GWD
      ENDIF
!
      CALL DSTRB(TEMP_GWD,int_state%HLENSW,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NGWD)TEMP_GWD
      ENDIF
!
      CALL DSTRB(TEMP_GWD,int_state%HLENNW,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NGWD)TEMP_GWD
      ENDIF
!
      CALL DSTRB(TEMP_GWD,int_state%HANGL,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NGWD)TEMP_GWD
      ENDIF
!
      CALL DSTRB(TEMP_GWD,int_state%HANIS,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NGWD)TEMP_GWD
      ENDIF
!
      CALL DSTRB(TEMP_GWD,int_state%HSLOP,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NGWD)TEMP_GWD
      ENDIF
!
      CALL DSTRB(TEMP_GWD,int_state%HZMAX,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!-----------------------------------------------------------------------
!
      IF(MYPE==0)THEN
        CLOSE(NGWD)
      ENDIF
!
      DEALLOCATE(TEMP_GWD)
!-----------------------------------------------------------------------
!
      END SUBROUTINE PHYSICS_READ_GWD 
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      end module module_INIT_READ_BIN
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------

