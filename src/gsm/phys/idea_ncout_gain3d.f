      SUBROUTINE IDEA_NCOUT_GAIN3D
     & (ioproc, t3d, levs, t3dpe, GLOBAL_LATS_R,LONSPERLAR)
!--------------------------------------------------------
!
      use resol_def,   ONLY: latr, lonr
      use layout1,     ONLY: me, lats_node_r, lats_node_r_max,
     &                       ipt_lats_node_r, nodes
      use mpi_def,     ONLY: info, mpi_comm_all, liope, mpi_r_io,
     &                       stat,mpi_r_io_r
      USE machine,     ONLY: kind_io4, kind_io8
     
!--------------------------------------------------------
      IMPLICIT NONE
!     
      integer :: ioproc
      integer :: levs
      integer :: global_lats_r(latr), lonsperlar(latr)
!!
      real (kind=kind_io8) t3d(lonr, latr,         levs)     !out
      real (kind=kind_io8) t3dpe(lonr,lats_node_r, levs)
      real (kind=kind_io8) xl(lonr,lats_node_r)  ! for in-copy
      real (kind=kind_io8) xf(lonr,lats_node_r)  ! for in-copy
      real (kind=kind_io8)  x(lonr,latr)         !out
      real (kind=kind_io8) tmp(lonr,latr+2)

      
      integer :: ipt_lats_node_rl,nodesr
      integer :: lats_nodes_rl
      integer :: maxfld, nproct
      integer ::  proc,j,lat,msgtag,nproc,i,buff,startlat,ierr
!    
      integer ifldu/0/
!      save ifldu
      integer illen,ncc
      integer :: K
      data ncc/0/
!
! t3dpe => xl => xf => X => t3d
!  
      DO k=1, levs
      xl(1:lonr,1:lats_node_r) = t3dpe(1:lonr,1:lats_node_r,K)
!
       CALL wam_fill(xl,xf,global_lats_r,lonsperlar)
       CALL wam_unsplit2d(ioproc,x,xf,global_lats_r)
!$omp parallel do private(i,j)
        do j=1,latr
           do i=1,lonr
              t3d(i, j,k) = X(i, j)     !fill-in is needed only for IOPROC
          enddo
        enddo
       call mpi_barrier(mpi_comm_all,ierr)
       enddo

 
       RETURN
!
      END SUBROUTINE IDEA_NCOUT_GAIN3D
!
      subroutine wam_unsplit2d(ioproc,x,xl,global_lats_r)
!
! I cannot get the ESSENCE: mpi_r_io_r => r8
!
      use resol_def,   ONLY: latr, lonr
      use layout1,     ONLY: me, lats_node_r, lats_node_r_max,
     &                       ipt_lats_node_r, nodes
      use mpi_def,     ONLY: info, mpi_comm_all, liope, mpi_r_io,
     &                       stat, mpi_r_io_r
      USE machine,     ONLY: kind_io4, kind_io8
      implicit none
!!
      real(kind=kind_io8)  x(lonr,latr)         !out
      real (kind=kind_io8) xl(lonr,lats_node_r) ! in

      real(kind=kind_io8) tmp(lonr,latr+2)


      integer global_lats_r(latr),ipt_lats_node_rl,nodesr
      integer lats_nodes_rl
      integer maxfld,ioproc,nproct
      integer proc,j,lat,msgtag,nproc,i,buff,startlat,ierr
!    
      integer ifldu/0/
      save ifldu
      integer illen,ncc
      data ncc/0/
!
!
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
!         msgtag = 1000 + (me+1)*nodesr        + ifldu
         msgtag = me + 1
!       write(6,*)' calling mpi_send for'
!     &,' pe=',me,' msgtag=',msgtag

         call MPI_SEND(tmp(lonr,latr+1),1,MPI_R_IO_R,ioproc,
     &                 msgtag,MPI_COMM_ALL,info)

          call MPI_SEND(tmp(lonr,latr+2),1,MPI_R_IO_R,ioproc,
     &                  msgtag,MPI_COMM_ALL,info)

         illen = tmp(lonr,latr+2)     ! lats_node_r
! send the local grid domain
         CALL mpi_send(tmp(1,1),illen*lonr,MPI_R_IO_R,ioproc,
     &                  msgtag,MPI_COMM_ALL,info)
!       write(6,*)' after mpi_send for'
!     &,' pe=',me,' msgtag=',msgtag
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
            CALL mpi_recv(tmp(lonr,latr+1),1,MPI_R_IO_R,proc-1,
     &                    msgtag,MPI_COMM_ALL,stat,info) 
!          write(6,*)' after mpi_recv1'
!     &,' pe=',me,' proc=',proc, tmp(lonr,latr+1)
            CALL mpi_recv(tmp(lonr,latr+2),1,MPI_R_IO_R,proc-1,
     &                    msgtag,MPI_COMM_ALL,stat,info)
!         write(6,*)' after mpi_recv2'
!     &,' pe=',me,' proc=', proc, tmp(lonr,latr+2)
            illen = int(tmp(lonr,latr+2))
            CALL mpi_recv(tmp(1,1),illen*lonr ,MPI_R_IO_R,proc-1,
     &                     msgtag,MPI_COMM_ALL,stat,info)
!         write(6,*)' after mpi_recv3'
!     &,' pe=',me,' proc=', proc, maxval(tmp(1:lonr,1:illen))
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
!
! Got it:  precise values for ipt_lats_node_rl  & lats_nodes_rl  
!                  tmp(i,j) = XL(i,j)
!
      return
      end subroutine wam_unsplit2d
!
      subroutine wam_fill(fin,fout,global_lats_r,lonsperlar)
!
      use resol_def,   ONLY: latr, lonr
      use layout1,     ONLY: lats_node_r, ipt_lats_node_r
      USE machine,     ONLY: kind_io8
      implicit none
!
      integer              global_lats_r(latr)

      integer,intent(in):: lonsperlar(latr)
      real(kind=kind_io8),intent(in) :: fin(lonr,lats_node_r)
      real(kind=kind_io8),intent(out):: fout(lonr,lats_node_r)

      integer j,lons,lat
!
      do j=1,lats_node_r

        lat  = global_lats_r(ipt_lats_node_r-1+j)
        lons = lonsperlar(lat)
        if(lons .lt. lonr) then
         call wam_intlon_phys(lons,lonr,fin(1,j),fout(1,j))
!
!          fin(1:lons,j) = fi(1:lons,j)
!          fout(1:lonr,j)=/=0.
        else
          fout(:,j) = fin(:,j)
        endif
      enddo
      end subroutine wam_fill
!
      subroutine wam_spread(fin,fout,global_lats_r,lonsperlar)
!
      use resol_def,   ONLY: latr, lonr
      use layout1,     ONLY: lats_node_r, ipt_lats_node_r
      USE machine,     ONLY: kind_io8
      implicit none
!
      integer ,intent(in)::global_lats_r(latr)

      integer,intent(in):: lonsperlar(latr)

      real(kind=kind_io8),intent(out):: fout(lonr,lats_node_r)

      real(kind=kind_io8),intent(in) :: fin(lonr,lats_node_r)
      integer j,lons,lat
!!
      do j=1,lats_node_r

        lat  = global_lats_r(ipt_lats_node_r-1+j)
        lons = lonsperlar(lat)
        if(lons .ne. lonr) then
          call wam_intlon_phys(lonr,lons,fin(1,j),fout(1,j))
!
!          f(1:lons,j) = fi(1:lons,j)
          fout(lons+1:lonr,j)=0.
        else
          fout(:,j) = fin(:,j)
        endif
      enddo
      end subroutine wam_spread
!***********************************************************************
      subroutine wam_intlon_phys(     m1,  m2,   f1,f2)
!                intlon_phys(  gain  lons<lonr & lonr > lons -spread)
!
! f2(m2) -output   spread( m2 < m1) gain (m2 > m1  and repeats values)
!
      use machine, ONLY: kind_io8
      implicit none
      integer,intent(in)              ::   m1,m2
   
      real (kind=kind_io8),intent(in) ::  f1(m1)
      real (kind=kind_io8),intent(out)::  f2(m2)
      integer i2,in,il,ir

      real (kind=kind_io8) r,x1

      r=real(m1)/real(m2)       !  lons/lonr 50/100

      do i2=1,m2                !       lonr
         x1=(i2-1)*r            !  shrink index i2=m1+1
         in=mod(nint(x1),m1)+1  !  copy i2 <=> "in = 
         f2(i2)=f1(in)         
      enddo

      end subroutine wam_intlon_phys
!***********************************************************************
