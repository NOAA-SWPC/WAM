      subroutine twrites_rst(fname,IOPROC,fhour,idate,
     &           si,ls_nodes,max_ls_nodes,step,
     &           trie_ls,trio_ls,trie_ls_rqt,trio_ls_rqt)
!
!-------------------------------------------------------------------
!*** program log
!*** Dec, 2009 Jun Wang:  write spectral variables for restart
!*** Dec, 2010 Jun Wang:  change to nemsio library
!*** Feb, 2011 Henry Juang: use generic argument for spectral fit to mass_dp and ndsl
!    Jun 26 2014 S. Moorthi - added mpigathe8
!-------------------------------------------------------------------
!
      use gfs_dyn_resol_def
      use gfs_dyn_layout1
      use gfs_dyn_coordinate_def	
      use namelist_dynamics_def
      use gfs_dyn_mpi_def
      use nemsio_module
!
      implicit none
!
      character(*),intent(in) :: fname
      integer,intent(in) :: ioproc
      real(kind=kind_evod),intent(in) :: fhour
      integer,intent(in) :: idate(4),step
!
      real(kind=kind_evod),intent(in):: si(levp1)
      integer,intent(in)             :: ls_nodes(ls_dim,nodes)
      integer,intent(in)             :: max_ls_nodes(nodes)
!
      real(kind=kind_evod),intent(in) :: trie_ls(len_trie_ls,2,lotls)
      real(kind=kind_evod),intent(in) :: trio_ls(len_trio_ls,2,lotls)
      real(kind=kind_evod),intent(in) :: trie_ls_rqt(len_trie_ls,2,levh)
      real(kind=kind_evod),intent(in) :: trio_ls_rqt(len_trio_ls,2,levh)
!
!local variables:
      REAL(kind=8) t1,t2,t3,t4,t5,t6,ta,tb,rtc
!
      integer              ierr,j,k,l,lenrec,locl,n,node,idate7(7)
!
      integer              indjoff
      integer              indev
      integer              indod
      integer              indev1,indev2
      integer              indod1,indod2
      integer              lots_rst
!
      real(kind=kind_ior), target ::   buf(lnt2)
      real(kind=kind_ior),allocatable ::   Z_R(:)
      real(kind=kind_ior)   tmps(4+nodes+jcap1*nodes)
      real(kind=kind_ior)   tmpr(3+nodes+jcap1*(nodes-1))
!
      type(nemsio_gfile) gfile
      integer nmetavarr8
      character(16),allocatable :: recname(:),reclevtyp(:)
      integer,allocatable :: reclev(:)
      character(16),allocatable :: varr8name(:)
      real(kind_ior),allocatable :: varr8val(:)
      integer iret, idvt
      integer il,ilen,i,msgtag,ls_diml,nrec
      integer nfhour,nfminute,nfsecondn,nfsecondd,nframe,jrec,nmeta
      logical first
      save first,nmetavarr8,recname,reclevtyp, 
     &     reclev,varr8name,varr8val,nmeta,nrec
      save Z_R
      data first /.true./
!
      integer              indlsev,jbasev
      integer              indlsod,jbasod
      integer              joff, jj
!
      include 'function2'
!
      joff(n,l)=(jcap1)*(jcap2)-(jcap1-l)*(jcap2-l)+2*(n-l)
!
      real(kind=kind_mpi_r),allocatable :: trieo_ls_node (:,:,:)
      real(kind=kind_mpi_r),allocatable :: trieo_ls_nodes(:,:,:,:)
!
      real(kind=kind_mpi_r),allocatable :: trieo_gz_node (:,:)
      real(kind=kind_mpi_r),allocatable :: trieo_gz_nodes(:,:,:)
!
      integer      lan,lat,iblk,lons_lat,lon,NJEFF,nn,lv
      integer      kwq,kwdp,kwte,kwdi,kwze,kwrq
      integer      kkq,kkdp,kkte,kkdi,kkze,kkrq
!
!---------------------------------------------------------------------
!
!      print *,' enter twrites_rst ' 

      kwq = 1
      kwdp = kwq + 1
      kwte = kwdp + levs
      kwdi = kwte + levs
      kwze = kwdi + levs
      kwrq = kwze + levs
      lots_rst = 4*levs + levh + 1

      call mpi_comm_size(MPI_COMM_ALL,i,ierr)
!
!-- allocate
      allocate ( trieo_ls_node  ( len_trie_ls_max+len_trio_ls_max,
     x                            2, lots_rst ) )
!
!-- compute gze_ls only once
      if(first ) then
        allocate(trieo_gz_node(len_trie_ls_max+len_trio_ls_max,2))
        do j=1,len_trie_ls
          trieo_gz_node(j,1) = trie_ls(j,1,p_gz)
          trieo_gz_node(j,2) = trie_ls(j,2,p_gz)
        enddo
        do j=1,len_trio_ls
          trieo_gz_node(j+len_trie_ls_max,1) = trio_ls(j,1,p_gz)
          trieo_gz_node(j+len_trie_ls_max,2) = trio_ls(j,2,p_gz)
        enddo
      endif
!
!-- collect data
      if( step==-1 ) then	! time step n-1
        kkq  = p_qm
        kkdp = p_dpm
        kkte = p_tem
        kkdi = p_dim
        kkze = p_zem
        kkrq = p_rm
      else if( step==0 ) then	! time step n
        kkq  = p_q
        kkdp = p_dp
        kkte = p_te
        kkdi = p_di
        kkze = p_ze
        kkrq = p_rq
      else
        print *,' **** error in twrites_rst: unknown step=',step
      endif
        
      do j=1,len_trie_ls
        trieo_ls_node(j,1,kwq) = trie_ls(j,1,kkq)
        trieo_ls_node(j,2,kwq) = trie_ls(j,2,kkq)
      enddo
!
      do j=1,len_trio_ls
        trieo_ls_node(j+len_trie_ls_max,1,kwq) = trio_ls(j,1,kkq)
        trieo_ls_node(j+len_trie_ls_max,2,kwq) = trio_ls(j,2,kkq)
      enddo
!
      do k=1,levs
        do j=1,len_trie_ls
          trieo_ls_node(j,1,kwdp+k-1) = trie_ls(j,1,kkdp+k-1)
          trieo_ls_node(j,2,kwdp+k-1) = trie_ls(j,2,kkdp+k-1)
          trieo_ls_node(j,1,kwte+k-1) = trie_ls(j,1,kkte+k-1)
          trieo_ls_node(j,2,kwte+k-1) = trie_ls(j,2,kkte+k-1)
          trieo_ls_node(j,1,kwdi+k-1) = trie_ls(j,1,kkdi+k-1)
          trieo_ls_node(j,2,kwdi+k-1) = trie_ls(j,2,kkdi+k-1)
          trieo_ls_node(j,1,kwze+k-1) = trie_ls(j,1,kkze+k-1)
          trieo_ls_node(j,2,kwze+k-1) = trie_ls(j,2,kkze+k-1)
        enddo
        do j=1,len_trio_ls
          jj = j+len_trie_ls_max
          trieo_ls_node(jj,1,kwdp+k-1) = trio_ls(j,1,kkdp+k-1)
          trieo_ls_node(jj,2,kwdp+k-1) = trio_ls(j,2,kkdp+k-1)
          trieo_ls_node(jj,1,kwte+k-1) = trio_ls(j,1,kkte+k-1)
          trieo_ls_node(jj,2,kwte+k-1) = trio_ls(j,2,kkte+k-1)
          trieo_ls_node(jj,1,kwdi+k-1) = trio_ls(j,1,kkdi+k-1)
          trieo_ls_node(jj,2,kwdi+k-1) = trio_ls(j,2,kkdi+k-1)
          trieo_ls_node(jj,1,kwze+k-1) = trio_ls(j,1,kkze+k-1)
          trieo_ls_node(jj,2,kwze+k-1) = trio_ls(j,2,kkze+k-1)
        enddo
      enddo
!
      if( .not. ndslfv ) then
      do k=1,levh
        do j=1,len_trie_ls
          trieo_ls_node(j,1,kwrq+k-1) = trie_ls(j,1,kkrq+k-1)
          trieo_ls_node(j,2,kwrq+k-1) = trie_ls(j,2,kkrq+k-1)
        enddo
        do j=1,len_trio_ls
          jj = j+len_trie_ls_max
          trieo_ls_node(jj,1,kwrq+k-1) = trio_ls(j,1,kkrq+k-1)
          trieo_ls_node(jj,2,kwrq+k-1) = trio_ls(j,2,kkrq+k-1)
        enddo
      enddo
      else
      do k=1,levh
        do j=1,len_trie_ls
          trieo_ls_node(j,1,kwrq+k-1) = trie_ls_rqt(j,1,k)
          trieo_ls_node(j,2,kwrq+k-1) = trie_ls_rqt(j,2,k)
        enddo
        do j=1,len_trio_ls
          jj = j+len_trie_ls_max
          trieo_ls_node(jj,1,kwrq+k-1) = trio_ls_rqt(j,1,k)
          trieo_ls_node(jj,2,kwrq+k-1) = trio_ls_rqt(j,2,k)
        enddo
      enddo
      endif
!
!-- collect data to ioproc
!WY bug fix.
!-----------
      if ( me .eq. ioproc ) then
        write(0,*)'ALLOC PARMS TWRITE ',len_trie_ls_max+len_trio_ls_max,
     &      2, lots_rst, nodes,1
 
         allocate ( trieo_ls_nodes ( len_trie_ls_max+len_trio_ls_max,
     &                               2, lots_rst, nodes ),
     &              stat=ierr )
         if(first) then
           allocate ( trieo_gz_nodes ( len_trie_ls_max+len_trio_ls_max,
     &                               2, nodes ),stat=ierr )
         endif
       else
         allocate (trieo_ls_nodes(1, 1, 1, 1), stat = ierr)
         if(first) then
           allocate (trieo_gz_nodes(1, 1, 1), stat = ierr )
         endif
      endif
      if (ierr .ne. 0) then
        write (0,*) ' GWX trieo_ls_nodes allocate failed'
        call mpi_abort(mpi_comm_all,ierr,i)
       endif
!
      call mpi_barrier(MPI_COMM_ALL,ierr)
      if(nodes >1 )then

        lenrec = (len_trie_ls_max+len_trio_ls_max) * 2 * lots_rst
!
        t1=rtc()
!       if (jcap < 1200 .and. ixgr <= 0) then
        if (jcap < 1200) then
          call mpi_gather(trieo_ls_node , lenrec, MPI_R_MPI_R,
     &                    trieo_ls_nodes, lenrec, MPI_R_MPI_R,
     &                    ioproc, MPI_COMM_ALL, ierr)
        else
          call mpi_gathe8(trieo_ls_node , lenrec, MPI_R_MPI_R,
     &                    trieo_ls_nodes, lenrec, MPI_R_MPI_R,
     &                    ioproc, MPI_COMM_ALL, ierr)
        endif
!       call mpi_gather( trieo_ls_node , lenrec, MPI_R_MPI_R,
!    x                  trieo_ls_nodes, lenrec, MPI_R_MPI_R,
!    x                  ioproc, MPI_COMM_ALL, ierr)
        call mpi_barrier(MPI_COMM_ALL,ierr)
        if(first) then
          lenrec=(len_trie_ls_max+len_trio_ls_max) * 2 
          call mpi_gather( trieo_gz_node , lenrec, MPI_R_MPI_R,
     &                     trieo_gz_nodes, lenrec, MPI_R_MPI_R,
     &                     ioproc, MPI_COMM_ALL, ierr)
        endif
      else
       trieo_ls_nodes(:,:,:,1)=trieo_ls_node(:,:,:)
       if(first) then
        trieo_gz_nodes(:,:,1)=trieo_gz_node(:,:)
       endif
      endif
      deallocate ( trieo_ls_node  )
      if(first) then
        deallocate ( trieo_gz_node  )
        if ( me .eq. ioproc ) then
!
         allocate(Z_R(lnt2) )
!
         do node=1,nodes
          jbasev=0
          do locl=1,max_ls_nodes(node)
            l=ls_nodes(locl,node)
            indev1 = indlsev(L,L)
            if (mod(L,2).eq.mod(jcap+1,2)) then
              indev2 = indlsev(jcap-1,L)
            else
              indev2 = indlsev(jcap  ,L)
            endif
            indjoff=joff(l,l)
            do indev = indev1 , indev2
              Z_R(indjoff+1) = trieo_gz_nodes(indev,1,node)
              Z_R(indjoff+2) = trieo_gz_nodes(indev,2,node)
              indjoff=indjoff+4
            end do
            jbasev=jbasev+(jcap+3-l)/2
           end do
!
           jbasod=len_trie_ls_max
           do locl=1,max_ls_nodes(node)
              l=ls_nodes(locl,node)
!
              if ( l .ne. jcap ) then  ! fix for out of bound error
                indod1 = indlsod(L+1,L)
                if (mod(L,2).eq.mod(jcap+1,2)) then
                  indod2 = indlsod(jcap  ,L)
                else
                  indod2 = indlsod(jcap-1,L)
                endif
                indjoff=joff(l+1,l)
                do indod = indod1 , indod2
                  Z_R(indjoff+1) = trieo_gz_nodes(indod,1,node)
                  Z_R(indjoff+2) = trieo_gz_nodes(indod,2,node)
                  indjoff=indjoff+4
                end do
                jbasod=jbasod+(jcap+2-l)/2
              endif  ! fix for out of bound error
            end do
         end do
        endif    !end ioproc
        deallocate(trieo_gz_nodes)
      endif
!
      t2=rtc()
      call mpi_barrier(MPI_COMM_ALL,ierr)
!
!-- write out data
      IF (me.eq.ioproc) THEN
 
        print *,' in TWRITES fhour=',fhour
!
        if (first) then
!
          nmeta=5
          nrec=1+lots_rst
          ntrac=levh/levs
!          print *,'write spec rst, nrec=',nrec,'ntrac=',ntrac
!
!-- field infomation
          allocate(recname(nrec),reclevtyp(nrec),reclev(nrec))
          recname(1)='gz'
          recname(2)='pres'
          recname(       3:  levs+2)='dpres'
          recname(  levs+3:2*levs+2)='tmp'
          recname(2*levs+3:3*levs+2)='di'
          recname(3*levs+3:4*levs+2)='ze'
          recname(4*levs+3:4*levs+levh+2)='rq'
          reclevtyp(1)='sfc'
          reclevtyp(2)='sfc'
          reclevtyp(3:4*levs+2)='mid layer'
          reclevtyp(4*levs+3:4*levs+levh+2)='tracer layer'
          reclev(1)=1
          reclev(2)=1
          do i=1,levs
            reclev(i+2)=i
            reclev(i+2+  levs)=i
            reclev(i+2+2*levs)=i
            reclev(i+2+3*levs)=i
          enddo
          do i=1,levh
            reclev(i+2+4*levs)=i
          enddo
!
          nmetavarr8=1
          allocate(varr8name(nmetavarr8),varr8val(nmetavarr8))
          varr8name(1:nmetavarr8)=(/'fhour  '/)
!
!endof first
        endif
!
        idate7=0.;idate7(7)=1
        idate7(1)=idate(4)
        idate7(2:3)=idate(2:3)
        idate7(4)=idate(1)
        nfhour=int(fhour)
        nfminute=int((fhour-nfhour)*60)
        nfsecondn=int(((fhour-nfhour)*60-nfminute)*60)
        nfsecondd=1
!
        varr8val(1)=fhour
!
        call nemsio_init()
!
        call nemsio_open(gfile,fname,'write',iret,modelname='GFS',   
     &    gdatatype='bin8',idate=idate7,nfhour=nfhour,nfminute=nfminute,
     &    nfsecondn=nfsecondn,nfsecondd=nfsecondd,dimx=lnt2,dimy=1,
     &    dimz=levs,nmeta=nmeta,nrec=nrec,jcap=jcap,ntrac=ntrac,
     &    recname=recname,reclevtyp=reclevtyp,reclev=reclev,
     &    extrameta=.true.,
     &    nmetavarr8=nmetavarr8,varr8name=varr8name,varr8val=varr8val)
!        print *,'after nemsio_open,',trim(fname),'iret=',iret
!
!--- write out data
!
       jrec=1
       call nemsio_writerec(gfile,1,Z_R,iret=iret)
!       print *,'after snemsio_writerec gz,',maxval(Z_r),minval(Z_R),
!     &   'iret=',iret
!
       do k=1,lots_rst
         jrec=k+1
         do node=1,nodes
!
           jbasev=0
           do locl=1,max_ls_nodes(node)
              l=ls_nodes(locl,node)
              indev1 = indlsev(L,L)
              if (mod(L,2).eq.mod(jcap+1,2)) then
                indev2 = indlsev(jcap-1,L)
              else
                indev2 = indlsev(jcap  ,L)
              endif
              indjoff=joff(l,l)
              do indev = indev1 , indev2
                buf(indjoff+1) = trieo_ls_nodes(indev,1,k,node)
                buf(indjoff+2) = trieo_ls_nodes(indev,2,k,node)
                indjoff=indjoff+4
              end do
              jbasev=jbasev+(jcap+3-l)/2
            end do
!
            jbasod=len_trie_ls_max
            do locl=1,max_ls_nodes(node)
              l=ls_nodes(locl,node)
!
              if ( l .ne. jcap ) then  ! fix for out of bound error
                indod1 = indlsod(L+1,L)
                if (mod(L,2).eq.mod(jcap+1,2)) then
                  indod2 = indlsod(jcap  ,L)
                else
                  indod2 = indlsod(jcap-1,L)
                endif
                indjoff=joff(l+1,l)
                do indod = indod1 , indod2
                  buf(indjoff+1) = trieo_ls_nodes(indod,1,k,node)
                  buf(indjoff+2) = trieo_ls_nodes(indod,2,k,node)
                  indjoff=indjoff+4
                end do

                jbasod=jbasod+(jcap+2-l)/2
              endif  ! fix for out of bound error
            end do
          end do
!
          call nemsio_writerec(gfile,jrec,buf,iret=iret)
!        print *,'after snemsio_writerec,jrec=',jrec,'buf=',maxval(buf),
!     &     minval(buf),'iret=',iret
        end do
!
        t4=rtc  ()
!sela print *, ' DISK TIME FOR SIG TWRITEO WRT ',t4-t3
!
!WY bug fix.
!-----------
!        deallocate ( trieo_ls_nodes )
!
        call nemsio_close(gfile)
        call nemsio_finalize()
!
!
      endif   !me.eq.ioproc
      deallocate ( trieo_ls_nodes )
!!
      if(first) then
          first = .false.
      endif
      call mpi_barrier(MPI_COMM_ALL,ierr)

      return
      end
