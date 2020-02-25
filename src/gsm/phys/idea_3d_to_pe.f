;==============================================
;
; spread/gain subroutines for WAM reduced gaus-grids
;
; copies of   split2d_phys
;           unsplit2d_phys
;           interpred_phys
;         uninterpred_phys
;              intlon_phys
;
;        spread3d_wam (split2d_wam  +  interpred_wam)
;        gain3d_wam   (unsplit2d_wam+uninterpred_wam)
;
;==============================================
	  idea_lat62_gaus.o              \
	  idea_iau_split.o              \
	  idea_ncout_gain3d.o           \
          idea_write_wamday.o           \
	  idea_ncout_phys.o              \
          idea_das_smlt.o               \
          idea_saber.o                  \







!  Read T   -- this is real temperature
!  All PEs read  records ???
!        etime = timef()
!        print *,'after nemsioheader,time=',etime-stime
      allocate (nemsio_data(lonf*latg))
      do k=1,levs
       if(w3rlkind==8) then
        call nemsio_readrecv(gfile_in,'tmp','mid layer',k,nemsio_data,
     &     iret=iret)
       elseif(w3rlkind==4) then
        call nemsio_readrecvw34(gfile_in,'tmp','mid layer',k,
     &     nemsio_data,iret=iret)
       endif
        call split2d(nemsio_data,buffo,global_lats_a)
        CALL interpred(1,kmsk,buffo,ttg(1,1,k),global_lats_a,lonsperlat)
      enddo
c**********************************************************************
c
      subroutine split2d(x,xl,global_lats_a)
c
c***********************************************************************
c
      use gfs_dyn_resol_def, ONLY: lonf, latg
      use gfs_dyn_layout1, ONLY: lats_node_a, ipt_lats_node_a, me
      use gfs_dyn_mpi_def, ONLY: MPI_R_IO, MPI_COMM_ALL, info
      USE gfs_dyn_MACHINE, ONLY: kind_io4, kind_io8
      implicit none
!!
      real(kind=kind_io4) x(lonf,latg)
      real (kind=kind_io8) xl(lonf,lats_node_a)
      real(kind=kind_io4) tmp(lonf,latg)
      integer global_lats_a(latg)
      integer nprocf,nodesr
!     integer maxfld,nprocf,nodesr
!     integer proc,j,lat,msgtag,nproc,i,msgtag1,buff,startlat,ierr
      integer proc,j,lat,nproc,i,buff,startlat,ierr
      integer ifld/0/
      save ifld
      real t1,t2,t3,t4,timef,ta,tb

!!
      XL=0.
    
      IF (me==0) THEN
!
!         Sending the data
!         ----------------
!-- do not need to send data, all processores read the data
         tmp=0.
         do j=1,latg
            do i=1,lonf
              tmp(i,j)=X(i,j)
            enddo
         enddo
      ENDIF
	
      call mpi_bcast
     1 (tmp,lonf*latg,MPI_R_IO,0,MPI_COMM_ALL,info)

       call mpi_barrier(mpi_comm_all,info)
!
!-- get subdomain of data from 2D-global
!
        do j=1,lats_node_a
           lat=global_lats_a(ipt_lats_node_a-1+j)
           do i=1,lonf
              xl(i,j)=tmp(i,lat)
           enddo
        enddo
      return
      end
!
!**********************************************************************

      subroutine interpred(iord,kmsk,f,fi,global_lats_a,lonsperlat)
!!
      use gfs_dyn_resol_def
      use gfs_dyn_layout1
      implicit none
!!
      integer              global_lats_a(latg)
      integer,intent(in):: iord
      integer,intent(in):: kmsk(lonf,lats_node_a)
      integer,intent(in):: lonsperlat(latg)
      real(kind=kind_io8),intent(in):: f(lonf,lats_node_a)
      real(kind=kind_io8),intent(out):: fi(lonf,lats_node_a)
      integer j,lons,lat
!!
      do j=1,lats_node_a
          lat=global_lats_a(ipt_lats_node_a-1+j)
          lons=lonsperlat(lat)
          if(lons.ne.lonf) then
            call intlon(iord,1,1,lonf,lons,
     &                  kmsk(1,j),f(1,j),fi(1,j))
cjfe        fi(lons+1:lonf,j)=-9999.e9
            fi(lons+1:lonf,j)=0.
          else
            fi(:,j)=f(:,j)
          endif

        enddo
      end subroutine
c
c***********************************************************************
cdynamics
      subroutine intlon(iord,imon,imsk,m1,m2,k1,f1,f2)
      use gfs_dyn_machine
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
!======================================================================
! from read-fix

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
      real(kind=kind_io8),intent(out):: fi(lonr,lats_node_r)   ! first elements 1:lons are taken
      integer j,lons,lat
!!
      do j=1,lats_node_r
          lat=global_lats_r(ipt_lats_node_r-1+j)
          lons=lonsperlar(lat)
          if(lons.ne.lonr) then
            call intlon_phys(iord,1,1,lonr,lons,
     &                  kmsk(1,j),f(1,j),fi(1,j))
cjfe        fi(lons+1:lonr,j)=-9999.e9
            fi(lons+1:lonr,j)=0.
          else
            fi(:,j)=f(:,j)
          endif
        enddo
      end subroutine interpred_phys     !(iord,kmsk,f,fi,global_lats_r,lonsperlar)
c=========================================
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
      real(kind=kind_io8),intent(out):: f(lonr,lats_node_r)     !filled array lonr
      real(kind=kind_io8),intent(in) :: fi(lonr,lats_node_r)    ! input (1:lons) =/=0
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
!***********************************************************************
!  for spread (fi, f) => (f1, f2) and f2-out everytime
!
      subroutine intlon_phys(iord,imon,imsk,m1, m2,  k1, f1, f2)
!           call intlon_phys(iord, 1,  1, lonr,lons,kmsk f   fi=f
!
! m1   >  m2   for interpred (lonr, lons)
! lonr  lons
! m1   <  m2   for uninterpred (lons, lonr)
! lons   lonr
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
