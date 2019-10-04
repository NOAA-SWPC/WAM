!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
    subroutine post_alctvars(imi,jmi,lmi,mypei,nwtlpes,mpicomp,mygridtype,  &
               mymaptype,post_gribversion,mynsoil,lead_write,jts,jte,jtsgrp,jtegrp)
!
!
!   revision history:
!    May 2011 Jun Wang: generate code from MPI_FIRST.f in post trunk
!    Jan 2012 Lu/Wang:  add gocart variables
!    Feb 2012 Jun Wang: add post_finalize
!    May 2013 Lu: add allocate_all
!
!
!-----------------------------------------------------------------------
!*** allocate post variables
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
!
!-----------------------------------------------------------------------
!
      implicit none
!
      include 'mpif.h'
!
!-----------------------------------------------------------------------
!
      integer,intent(in)            :: imi,jmi,lmi,mypei,nwtlpes,mpicomp
      character(1),intent(in)       :: mygridtype
      character(5),intent(in)       :: post_gribversion
      integer,intent(in)            :: mymaptype,mynsoil
      integer,intent(in)            :: lead_write
      integer,intent(in)            :: jts,jte
      integer,intent(in)            :: jtsgrp(nwtlpes),jtegrp(nwtlpes)
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      integer ii,jj,i,j,l,NPOS_END,NPOS_START,indx_2d,indx_num,nfield
      integer ierr,last_write_task
      integer indx,ll,mype
!
      REAL,allocatable :: SLDPTH2(:),dxh(:),dummy(:,:)
      REAL FACT,tsph,tstart
      REAL RINC(5)
!
!-----------------------------------------------------------------------
!*** get dims from int_state
!-----------------------------------------------------------------------
!
      print *,'in post_alctvars,im=',imi,'jm=',jmi,'lm=',lmi,'grib=',post_gribversion
      im=imi
      jm=jmi
      lm=lmi
      im_jm=im*jm
      lp1=lm + 1
      grib=trim(post_gribversion)
! set ndegr
      if(grib=='grib1') then
        gdsdegr=1000.
      else if (grib=='grib2') then
        gdsdegr=1000000.
      endif
      IOFORM='grib'
      mype=mypei
      me=mype-lead_write
      last_write_task=lead_write+nwtlpes-1
      MPI_COMM_COMP=mpicomp
      num_servers=nwtlpes
      NUM_PROCS=nwtlpes
      NUM_SERVERS=0
      GRIDTYPE=mygridtype
      MAPTYPE=mymaptype
      NSOIL=mynsoil
      print *,'grib=',grib,'ioform=',ioform,'mype=',mype,'me=',me, &
         'lead_write=',lead_write,'last_write_task=',last_write_task, &
         'num_servers=',num_servers,'NUM_PROCS=',NUM_PROCS,'GRIDTYPE=', &
          GRIDTYPE,'maptype=',maptype,'nsoil=',nsoil,'gdsdegr=',gdsdegr
!
      allocate(dxh(jm))
      allocate(SLDPTH2(nsoil))
!
!-----------------------------------------------------------------------
!***  ALLOCATE THE ARRAYS OF THE POST.
!-----------------------------------------------------------------------
!
      jsta=jts
      jend=jte
      jsta_m  = jsta
      jsta_m2 = jsta
      jend_m  = jend
      jend_m2 = jend
      if ( mype .eq. lead_write ) then
         jsta_m  = 2
         jsta_m2 = 3
      end if
      if ( mype .eq. last_write_task ) then
         jend_m  = jm - 1
         jend_m2 = jm - 2
      end if
!** neighbors
!jw      iup = me + 1
!jw      idn = me - 1
      iup = mype + 1 - lead_write
      idn = mype - 1 - lead_write
      if ( mype .eq. lead_write ) then
         idn = MPI_PROC_NULL
      end if
      if ( mype .eq. last_write_task ) then
         iup = MPI_PROC_NULL
      end if
      print *,'lead_write_task=',lead_write,'last taks=',last_write_task, &
        'idn=',idn,'iup=',iup,'MPI_PROC_NULL=',MPI_PROC_NULL,'jsta=',jsta,'jend=',jend
!
!     counts, disps for gatherv and scatterv
!
      do i = 1, NUM_PROCS
       icnt(i-1) = (jtegrp(i)-jtsgrp(i)+1)*im
       idsp(i-1) = (jtsgrp(i)-1)*im
       if ( mype .eq. lead_write ) then
           print *, ' i, icnt(i),idsp(i) = ',i-1,icnt(i-1),idsp(i-1)
       end if
      enddo
!
!     extraction limits -- set to two rows
!
      jsta_2l = max(jsta - 2,  1 )
      jend_2u = min(jend + 2, jm )
! special for c-grid v
      jvend_2u = min(jend + 2, jm+1 )
      print *,'im=',im,'jsta_2l=',jsta_2l,'jend_2u=',jend_2u,'lm=',lm
!
!
! SETS UP MESSAGE PASSING INFO

      call ALLOCATE_ALL()

!***
! LMH always = LM for sigma-type vert coord
! LMV always = LM for sigma-type vert coord

       do j = jsta_2l, jend_2u
        do i = 1, im
            LMV ( i, j ) = lm
            LMH ( i, j ) = lm
        end do
       end do
!
! HTM VTM all 1 for sigma-type vert coord

      do l = 1, lm
       do j = jsta_2l, jend_2u
        do i = 1, im
            HTM ( i, j, l ) = 1.0
            VTM ( i, j, l ) = 1.0
        end do
       end do
      end do
    end subroutine post_alctvars
!
!---------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!---------------------------------------------------------------------
!
  subroutine read_postnmlt(kpo,kth,kpv,po,th,pv,nlunit,post_namelist)
!
      use ctlblk_mod, only : komax,fileNameD3D,lsm,lsmp1,SPL,SPLDEF,  &
                              lsmdef,ALSL,me,d3d_on,gocart_on
!
      implicit none
!---
      character (len=*), intent(in) :: post_namelist
      integer :: kpo,kth,kpv,nlunit
      real,dimension(komax) :: po,th,pv
      namelist/nampgb/kpo,po,kth,th,kpv,pv,d3d_on,gocart_on
      integer l,k,iret
!---------------------------------------------------------------------
!
      print *,'in read_postnmlt'
!
! set default for kpo, kth, th, kpv, pv
      kpo = 0
      po  = 0
      kth = 1
      th  = (/320.,(0.,k=kth+1,komax)/) ! isentropic level to output
      kpv = 8
      pv  = (/0.5,-0.5,1.0,-1.0,1.5,-1.5,2.0,-2.0,(0.,k=kpv+1,komax)/)
      d3d_on    = .false.
      gocart_on = .false.
!
      if (me == 0) print *,' nlunit=',nlunit,' post_namelist=', &
     &                      post_namelist
!     read(5,nampgb,iostat=iret,end=118)
      if (nlunit > 0) then
        open (unit=nlunit,file=post_namelist)
        rewind(nlunit)
        read(nlunit,nampgb,iostat=iret,end=118)
      endif
 118  continue
      if (me == 0) then
        print*,'komax,iret for nampgb= ',komax,iret
        print*,'komax,kpo,kth,th,kpv,pv= ',komax,kpo            &
     &  ,kth,th(1:kth),kpv,pv(1:kpv),' gocart_on=',gocart_on
       endif
       fileNameD3D = '/dev/null'
!
!119  continue
! set up pressure level from POSTGPVARS or DEFAULT
      if(kpo == 0)then
! use default pressure levels
        print*,'using default pressure levels,spldef=',(spldef(l),l=1,lsmdef)
        lsm = lsmdef
        do l=1,lsm
         spl(l) = spldef(l)
        end do
      else
! use POSTGPVARS
        print*,'using pressure levels from POSTGPVARS'
        lsm = kpo
        if(po(lsm)<po(1))then ! post logic assumes asscending
         do l=1,lsm
          spl(l) = po(lsm-l+1)*100.
         end do
        else
         do l=1,lsm
          spl(l) = po(l)*100.
         end do
        end if
      end if
      print*,'LSM, SPL = ',lsm,spl(1:lsm)
      lsmp1 = lsm + 1
!
!     COMPUTE DERIVED MAP OUTPUT CONSTANTS.
      DO L = 1,LSM
!jw real4         ALSL(L) = ALOG(SPL(L))
         ALSL(L) = LOG(SPL(L))
      END DO
      write(0,*)' after ALSL'
!
1000  continue

      end subroutine read_postnmlt
!
!---------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!---------------------------------------------------------------------
!
    subroutine post_finalize(post_gribversion)
!
      use grib2_module, only : grib_info_finalize
!
      character(*),intent(in) :: post_gribversion
!
      IF(trim(POST_GRIBVERSION)=='grib2') then
         call  grib_info_finalize()
      ENDIF
!
      call de_allocate
!
    end subroutine post_finalize

