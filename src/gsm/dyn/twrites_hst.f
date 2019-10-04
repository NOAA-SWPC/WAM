      subroutine twrites_hst(fname,IOPROC,fhour,idate,
     &           ls_nodes,max_ls_nodes,trie_ls,trio_ls,
     &           trie_ls_rqt,trio_ls_rqt,
     &           thermodyn_id_out,sfcpress_id_out,pdryini)
!
!-------------------------------------------------------------------
!*** program log
!*** Sep, 2012 Jun Wang:  modified from twrites_rst
!    Jun 26 2014 S. Moorthi - added mpigathe8
!-------------------------------------------------------------------
!
      use gfs_dyn_resol_def
      use gfs_dyn_layout1
      use gfs_dyn_coordinate_def
      use namelist_dynamics_def
      use gfs_dyn_mpi_def
      use gfs_dyn_io_header, only : z,idvt,icen2,ienst,iensi,lonr,latr
!     use gfs_dyn_coordinate_def, only: vertcoord_id,ak5,bk5,ck5
      use gfs_dyn_tracer_const, only : cpi,ri
      use sigio_module
      use sigio_r_module
!
      implicit none
!
      character(*),intent(in)         :: fname
      integer,intent(in)              :: ioproc
      real(kind=kind_evod),intent(in) :: fhour
      integer,intent(in)              :: idate(4)
!
      integer,intent(in)              :: ls_nodes(ls_dim,nodes)
      integer,intent(in)              :: max_ls_nodes(nodes)
      REAL (KIND=KIND_grid),intent(in):: pdryini
      integer,intent(in)              :: thermodyn_id_out
      integer,intent(in)              :: sfcpress_id_out
!
      real(kind=kind_evod),intent(in) :: trie_ls(len_trie_ls,2,lotls)
      real(kind=kind_evod),intent(in) :: trio_ls(len_trio_ls,2,lotls)
      real(kind=kind_evod),intent(in) :: trie_ls_rqt(len_trie_ls,2,levh)
      real(kind=kind_evod),intent(in) :: trio_ls_rqt(len_trio_ls,2,levh)
!
!local variables:
      REAL(kind=8) t1,t2,t3,t4,t5,t6,ta,tb,rtc
!
      integer              i,ierr,j,k,l,lenrec,locl,n,node,idate7(7)
      integer              nsig,   indjoff, indev, indod, indev1, indev2
     &,                    indod1, indod2, lots_hst
!
      real(kind=kind_io4), target ::   buf(lnt2)
!
!sigio
      type(sigio_head) head
      type(sigio_dati) dati
      integer iret, num_dta
      logical first
      save head, first
      data first /.true./
!
      integer              indlsev, jbasev, indlsod, jbasod, joff, jj
!
      include 'function2'
!
      joff(n,l) = (jcap1)*(jcap2) - (jcap1-l)*(jcap2-l) + 2*(n-l)
!
      real(kind=kind_mpi_r),allocatable :: trieo_ls_node (:,:,:)
      real(kind=kind_mpi_r),allocatable :: trieo_ls_nodes(:,:,:,:)
!
      real(kind=kind_mpi_r),allocatable :: trieo_gz_node (:,:)
      real(kind=kind_mpi_r),allocatable :: trieo_gz_nodes(:,:,:)
!
      integer      lan,lat,iblk,lons_lat,lon,NJEFF,nn,lv
      integer      kwq,kwdp,kwte,kwdi,kwrq
      integer      kkq,kkdp,kkte,kkdi,kkze,kkrq
!
!---------------------------------------------------------------------
!
!     if (me == 0)  print *,' enter twrites_hst ' 

      kwq = 1
      kwte = kwq + 1
      kwdi = kwte + levs
      kwrq = kwdi + 2*levs
      lots_hst = 3*levs + levh + 1

      call mpi_comm_size(MPI_COMM_ALL,i,ierr)
!
!-- allocate
      allocate ( trieo_ls_node  ( len_trie_ls_max+len_trio_ls_max,
     &                            2, lots_hst ) )
!
!-- get p_ze only once from input
!      if(first ) then
!        allocate(trieo_gz_node(len_trie_ls_max+len_trio_ls_max,2))
!        do j=1,len_trie_ls
!          trieo_gz_node(j,1) = trie_ls(j,1,p_gz)
!          trieo_gz_node(j,2) = trie_ls(j,2,p_gz)
!        enddo
!        do j=1,len_trio_ls
!          trieo_gz_node(j+len_trie_ls_max,1) = trio_ls(j,1,p_gz)
!          trieo_gz_node(j+len_trie_ls_max,2) = trio_ls(j,2,p_gz)
!        enddo
!      endif
!
!-- collect data
      kkq  = p_q
      kkte = p_te
      kkdi = p_di
      kkze = p_ze
        
!     write(0,*)' in twrite_hst p_q=',p_q,' p_te=',p_te,
!    &' p_di=',p_di,' p_ze=',p_ze,' p_rq=',p_rq
!    &,'trie=',trie_ls(1,1,kkq)

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
!     write(0,*)' in hist kkte=',kkte,'temp=',trie_ls(1,1,kkte+k-1)
!    &,trie_ls(1,1,kkdi+k-1),trie_ls(1,1,kkze+k-1)
        do j=1,len_trie_ls
          trieo_ls_node(j,1,kwte+k-1)   = trie_ls(j,1,kkte+k-1)
          trieo_ls_node(j,2,kwte+k-1)   = trie_ls(j,2,kkte+k-1)
          trieo_ls_node(j,1,kwdi+2*k-2) = trie_ls(j,1,kkdi+k-1)
          trieo_ls_node(j,2,kwdi+2*k-2) = trie_ls(j,2,kkdi+k-1)
          trieo_ls_node(j,1,kwdi+2*k-1) = trie_ls(j,1,kkze+k-1)
          trieo_ls_node(j,2,kwdi+2*k-1) = trie_ls(j,2,kkze+k-1)
        enddo
        do j=1,len_trio_ls
          jj = j+len_trie_ls_max
          trieo_ls_node(jj,1,kwte+k-1)   = trio_ls(j,1,kkte+k-1)
          trieo_ls_node(jj,2,kwte+k-1)   = trio_ls(j,2,kkte+k-1)
          trieo_ls_node(jj,1,kwdi+2*k-2) = trio_ls(j,1,kkdi+k-1)
          trieo_ls_node(jj,2,kwdi+2*k-2) = trio_ls(j,2,kkdi+k-1)
          trieo_ls_node(jj,1,kwdi+2*k-1) = trio_ls(j,1,kkze+k-1)
          trieo_ls_node(jj,2,kwdi+2*k-1) = trio_ls(j,2,kkze+k-1)
        enddo
      enddo
!
      if( .not. ndslfv ) then
        kkrq = p_rq
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
      if ( me  == ioproc ) then

!       write(0,*)'ALLOC PARMS TWRITE ',len_trie_ls_max+len_trio_ls_max,
!     &      2, lots_hst, nodes,1
 
         allocate ( trieo_ls_nodes ( len_trie_ls_max+len_trio_ls_max,
     &                               2, lots_hst, nodes ), stat=ierr )
      else
         allocate ( trieo_ls_nodes (1,1,1,1), stat=ierr )
      endif
      if (ierr /= 0) then
        write (0,*) ' GWX trieo_ls_nodes allocate failed'
        call mpi_abort(mpi_comm_all,ierr,i)
      endif
!
      call mpi_barrier(MPI_COMM_ALL,ierr)
!
!     if (me == 0) print *,'in twrites,collect data'   
      if(nodes >1 )then

        lenrec = (len_trie_ls_max+len_trio_ls_max) * 2 * lots_hst
!
        t1=rtc()
!       call mpi_gather( trieo_ls_node , lenrec, MPI_R_MPI_R,
!    &                   trieo_ls_nodes, lenrec, MPI_R_MPI_R,
!    &                   ioproc, MPI_COMM_ALL, ierr)
!     if (jcap < 1200 .and. ixgr <= 0) then
      if (jcap < 1200) then
        call mpi_gather(trieo_ls_node , lenrec, MPI_R_MPI_R,
     &                  trieo_ls_nodes, lenrec, MPI_R_MPI_R,
     &                  ioproc, MPI_COMM_ALL, ierr)
      else
        call mpi_gathe8(trieo_ls_node , lenrec, MPI_R_MPI_R,
     &                  trieo_ls_nodes, lenrec, MPI_R_MPI_R,
     &                  ioproc, MPI_COMM_ALL, ierr)
      endif
        call mpi_barrier(MPI_COMM_ALL,ierr)
      else
        trieo_ls_nodes(:,:,:,1) = trieo_ls_node(:,:,:)
      endif
      deallocate ( trieo_ls_node  )
!
!     t2=rtc()

      call mpi_barrier(MPI_COMM_ALL,ierr)
!
!-- write out data
!     if (me == 0) print *,'in twrites,write out data,ioprec=',ioproc   

      if (me == ioproc) then
 
        nsig=61
!       print *,' in TWRITES bf sigio_rwopen=',trim(fname)
        call sigio_rwopen(nsig,trim(fname),iret)
!       print *,' in TWRITES fhour=',fhour,'iret=',iret
!
        if (first) then
!
          first = .false.
!
!change to sigio output
!
          head%ivs = 198410
          if (levs > 99) head%ivs = 200509
          head%levs    = levs
          head%ntrac   = ntrac
          if (gen_coord_hybrid) then                            ! hmhj
            head%ivs = 200509
            head%nvcoord = 3                                    ! hmhj
          else if (hybrid) then                                 ! hmhj
            head%nvcoord = 2
          else
            head%nvcoord = 1
          endif
          idvm = thermodyn_id_out*10 + sfcpress_id_out
!        print *,' before alhead idvm=',idvm,thermodyn_id_out,
!    &      sfcpress_id_out
          call sigio_alhead(head,iret,levs,head%nvcoord,ntrac,idvm)
          head%nxgr = 0
          head%nxss = 0
        endif
        rewind(nsig)
        head%clabsig=char(0)//char(0)//char(0)//char(0)//
     &               char(0)//char(0)//char(0)//char(0)
!
        head%fhour   = fhour
        head%idate   = idate
        head%jcap    = jcap
        head%latb    = latr
        head%lonb    = lonr
        head%itrun   = 1
        head%iorder  = 2
        head%irealf  = 1      ! for real * 4 data
        head%igen    = igen
        head%latf    = latg
        head%latr    = latr
        head%lonf    = lonf
        head%lonr    = lonr
        head%icen2   = icen2
        head%iens(1) = ienst
        head%iens(2) = iensi
        head%idpp    = 0
        head%idsl    = idsl ! idsl=2 for middle of layer
        head%idvm    = 0
        head%idvt    = idvt
        head%idrun   = 0
        head%idusr   = 0
        head%pdryini = pdryini
        head%ncldt   = ncld
        ixgr         = 0              !!! Set by Moorthi
        head%ixgr    = ixgr
!!!     head%ldata(:) = lnt2          !!! Commented by Moorthi
!!
        if (gen_coord_hybrid) then                              ! hmhj
          head%idvc    = 3                                      ! hmhj
          head%idvm    = thermodyn_id*10+sfcpress_id            ! hmhj
          head%idsl    = 2                                      ! hmhj
          do k=1,levp1                                          ! hmhj
            head%vcoord(k,1) = ak5(k)*1000.                     ! hmhj
            head%vcoord(k,2) = bk5(k)                           ! hmhj
            head%vcoord(k,3) = ck5(k)*1000.                     ! hmhj
          enddo                                                 ! hmhj
          if (thermodyn_id == 3) then
            head%cpi(1:ntrac+1) = cpi(0:ntrac)
            head%ri(1:ntrac+1)  = ri(0:ntrac)
          endif
!      print *,' thermodyn_id=',thermodyn_id_out,
!     & ' cpi=',head%cpi(0:ntrac+1)
!     &,' ri=',head%ri(1:ntrac+1)
        else if (hybrid) then                                   ! hmhj
          head%idvc    = 2    ! for hybrid vertical coord.
          head%idsl    = 1                                      ! hmhj
          do k=1,levp1
            head%vcoord(k,1) = ak5(levp1+1-k)*1000.
            head%vcoord(k,2) = bk5(levp1+1-k)

!        print 190,k,head%vcoord(k,1),head%vcoord(k,2)
!190   format('in twrite k=',i2,'  ak5r4=',f13.6,'  bk5r4=',e13.5)
          enddo
        else
          head%idsl    = 1                                      ! hmhj
          head%idvc    = 1    ! for sigma vertical coord. (default)
          head%vcoord(:,1) = ak5(:)
        endif
!
        call sigio_rwhead(nsig,head,iret)
!write out z

        buf    = Z
        dati%i = 1
        dati%f => buf
        call sigio_rwdati(nsig,head,dati,iret)
!
        do k=1,lots_hst
          do node=1,nodes
!
            jbasev=0
            do locl=1,max_ls_nodes(node)
              l = ls_nodes(locl,node)
              indjoff = joff(l,l)
              do indev = indlsev(l,l) , indlsev(jcap-mod(l,2),l)
                buf(indjoff+1) = trieo_ls_nodes(indev,1,k,node)
                buf(indjoff+2) = trieo_ls_nodes(indev,2,k,node)
                indjoff = indjoff + 4
              end do
              jbasev = jbasev + (jcap+3-l)/2
            end do
!
            jbasod = len_trie_ls_max
            do locl=1,max_ls_nodes(node)
              l = ls_nodes(locl,node)
              if ( l .ne. jcap ) then  ! fix for out of bound error
                indjoff = joff(l+1,l)
                do indod = indlsod(l+1,l) , indlsod(jcap-mod(l+1,2),l)
                  buf(indjoff+1) = trieo_ls_nodes(indod,1,k,node)
                  buf(indjoff+2) = trieo_ls_nodes(indod,2,k,node)
                  indjoff = indjoff + 4
                end do
                jbasod = jbasod + (jcap+2-l)/2
              endif  ! fix for out of bound error
            end do
          end do
!
          dati%i = k+1
          dati%f => buf
          call sigio_rwdati(nsig,head,dati,iret)
        end do
!
        print *,' twriteeo fhour=',fhour,' idate=',idate,' nsig=',nsig
!
      endif   !me.eq.ioproc
      deallocate ( trieo_ls_nodes )
!!
      call mpi_barrier(MPI_COMM_ALL,ierr)

      return
      end
