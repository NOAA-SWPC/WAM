      SUBROUTINE grid_collect(zsg,psg,uug,vvg,teg,rqg,dpg,
     &                        global_lats_a,lonsperlat)
!!
!! Revision history:
!  2007           Henry Juang, original code
!  2008           Jun Wang  modified buff for write grid component
!  Nov 23 2009    Sarah Lu, comment out 4D tracer
!  Sep 08 2010    Jun Wang  change gfsio to nemsio
!  Dec 16 2010    Jun Wang  change to nemsio library
!  Feb 20 2011    Hann-Ming Henry Juang add code for NDSL
!  Sep 24 2014    S Moorthi - some cleanup and optimization
!  Feb 04 2015    S. Moorthi - threading and optimization
!

      use gfs_dyn_resol_def
      use gfs_dyn_write_state,    only: buff_mult_pieceg,ngrid,ngridg
      use gfs_dyn_layout1
      use gfs_dyn_mpi_def
      use namelist_dynamics_def,  only: gen_coord_hybrid,iau
      use gfs_dyn_coordinate_def, only: ak5, bk5
      use gfs_dyn_physcons, rk => con_rocp

      implicit none

      real (kind=kind_grid), parameter :: rk1 = rk + 1.0, rkr = 1.0/rk
     &,                                   p0=100000.0, p0i=1.0/p0
      real (kind=kind_io4), parameter  :: zero4=0.0
!
      integer, dimension(latg) :: global_lats_a, lonsperlat
!
      real(kind=kind_grid), dimension(lonf,lats_node_a)      ::  zsg,psg
      real(kind=kind_grid), dimension(lonf,lats_node_a,levs) ::  uug,vvg
     &,                                                          teg,dpg
      real(kind=kind_grid), dimension(lonf,lats_node_a,levh) ::  rqg
!
!     real(kind=kind_io4), dimension(lonf,lats_node_a) :: pup,  pdn,
!    &                                                    pupk, pdnk
      real(kind=kind_io4), dimension(lonf,lats_node_a) :: pdn, plyr,pdnk
      real(kind=kind_io4)                              :: pup, pupk
!
      real(kind=kind_io8), dimension(lonf,lats_node_a) :: buffo, buffi
      integer, dimension(lonf,lats_node_a)             :: kmsk
      integer i, j, k, kk, il, icount, ierr, ll
!     integer i, j, k, kk, il, icount, ierr, maxlats_comp
!
      data  icount/0/
!
      ngridg = 1
      if(.not. allocated(buff_mult_pieceg)) then
         allocate(buff_mult_pieceg(lonf,lats_node_a_max,ngrids_gg))
      endif
!
!$omp parallel do private(i,j)
      do j=1,lats_node_a
        do i=1,lonf
          kmsk(i,j)  = 0
          buffi(i,j) = zsg(i,j)
        enddo
      enddo
!
      CALL uninterpreg(1,kmsk,buffo,buffi,global_lats_a,lonsperlat,
     &                 buff_mult_pieceg(1,1,1) )

!      write(1000+me,*)'in grid collect, buff_zsg=',
!      write(0,*)'in grid collect, buff_zsg=',' me=',me,
!    &  maxval(buff_mult_pieceg(1:lonf,1:lats_node_a,1)),
!    &  minval(buff_mult_pieceg(1:lonf,1:lats_node_a,1))
!
!$omp parallel do private(i,j)
      do j=1,lats_node_a
        do i=1,lonf
          buffi(i,j) = psg(i,j)
        enddo
      enddo
      CALL uninterpreg(1,kmsk,buffo,buffi,global_lats_a,lonsperlat,
     &                 buff_mult_pieceg(1,1,2) )
!      write(0,*)'in grid collect, buff_psg=',' me=',me,
!    & maxval(buff_mult_pieceg(1:lonf,1:lats_node_a,2)),
!    & minval(buff_mult_pieceg(1:lonf,1:lats_node_a,2))
!      print *,'in grid collect, buff_psg(135,23)=',
!     &  buff_mult_pieceg(135,23,2)
!
      do k=1,levs
!$omp parallel do private(i,j)
        do j=1,lats_node_a
          do i=1,lonf
            buffi(i,j) = dpg(i,j,k)
          enddo
        enddo
        CALL uninterpreg(1,kmsk,buffo,buffi,global_lats_a,lonsperlat,
     &                   buff_mult_pieceg(1,1,2+k) )
!      write(0,*)'in grid collect, buff_dpg=',' me=',me,
!    & maxval(buff_mult_pieceg(1:lonf,1:lats_node_a,2+k)),
!    & minval(buff_mult_pieceg(1:lonf,1:lats_node_a,2+k))
      enddo
!!
!***  write out the layer mean pressure
!$omp parallel do private(i,j)
      do j=1,lats_node_a
        do i=1,lonf
          pdn(i,j)  = buff_mult_pieceg(i,j,2)
          pdnk(i,j) = (pdn(i,j)*p0i) ** rk
        enddo
      enddo
      if (gen_coord_hybrid) then
        do k=1,levs
          kk = 2+levs+k
!$omp parallel do private(i,j,pup)
          do j=1,lats_node_a
            do i=1,lonf
              pup  = max(pdn(i,j)-buff_mult_pieceg(i,j,2+k),zero4)
              buff_mult_pieceg(i,j,2+levs+k) = 0.5*(pup + pdn(i,j))
              pdn(i,j)  = pup
            enddo
          enddo
        enddo
      else
        do k=1,levs
          kk = 2+levs+k
!!$omp parallel do private(i,j)
!         do j=1,lats_node_a
!           do i=1,lonf
!             pup(i,j)  = max(pdn(i,j)-buff_mult_pieceg(i,j,2+k),zero4)
!             if (pupk(i,j) > 1.0e-6) then
!               pupk(i,j) = (pup(i,j)*p0i) ** rk
!             else
!               pupk(i,j) = zero4
!             endif
!             buff_mult_pieceg(i,j,kk) = p0*((pdnk(i,j)*pdn(i,j) -
!     &          pupk(i,j)*pup(i,j)) /(rk1*(pdn(i,j)-pup(i,j)))) ** rkr
!             pdn(i,j)  = pup(i,j)
!             pdnk(i,j) = pupk(i,j)
!           enddo
!         enddo
!       enddo
          ll = levs + 1 - k
!$omp parallel do private(i,j,pup,pupk)
          do j=1,lats_node_a
            do i=1,lonf
              pup  = ak5(ll) + bk5(ll)*buff_mult_pieceg(i,j,2)
              if (pup > 1.0e-6) then
                pupk = (pup*p0i) ** rk
              else
                pupk = zero4
              endif
              buff_mult_pieceg(i,j,kk) = p0*((pdnk(i,j)*pdn(i,j) -
     &           pupk*pup) /(rk1*(pdn(i,j)-pup))) ** rkr
              pdn(i,j)  = pup
              pdnk(i,j) = pupk
            enddo
          enddo
!      write(0,*)'in grid collect, buff_pg=',' me=',me,
!    & maxval(buff_mult_pieceg(1:lonf,1:lats_node_a,kk)),
!    & minval(buff_mult_pieceg(1:lonf,1:lats_node_a,kk))
!    &,' k=',k,'ll=',ll,' akbk=',ak5(ll),bk5(ll)
!
        enddo
      endif
!
      do k=1,levs
!$omp parallel do private(i,j)
        do j=1,lats_node_a
          do i=1,lonf
            buffi(i,j) = uug(i,j,k)
          enddo
        enddo
        CALL uninterpreg(1,kmsk,buffo,buffi,global_lats_a,lonsperlat,
     &                   buff_mult_pieceg(1,1,2+2*levs+k) )
!      write(0,*)'in grid collect, buff_uug=',' me=',me,
!    & maxval(buff_mult_pieceg(1:lonf,1:lats_node_a,2+2*levs+k)),
!    & minval(buff_mult_pieceg(1:lonf,1:lats_node_a,2+2*levs+k))
      enddo
!
      do k=1,levs
!$omp parallel do private(i,j)
        do j=1,lats_node_a
          do i=1,lonf
            buffi(i,j) = vvg(i,j,k)
          enddo
        enddo
        CALL uninterpreg(1,kmsk,buffo,buffi,global_lats_a,lonsperlat,
     &                   buff_mult_pieceg(1,1,2+3*levs+k) )
!      write(0,*)'in grid collect, buff_vvg=',' me=',me,
!    & maxval(buff_mult_pieceg(1:lonf,1:lats_node_a,2+3*levs+k)),
!    & minval(buff_mult_pieceg(1:lonf,1:lats_node_a,2+3*levs+k))
      enddo
!
      do k=1,levs
!$omp parallel do private(i,j)
        do j=1,lats_node_a
          do i=1,lonf
            buffi(i,j) = teg(i,j,k)
          enddo
        enddo
        CALL uninterpreg(1,kmsk,buffo,buffi,global_lats_a,lonsperlat,
     &                   buff_mult_pieceg(1,1,2+4*levs+k) )
!      write(0,*)'in grid collect, buff_teg=',' me=',me,
!    & maxval(buff_mult_pieceg(1:lonf,1:lats_node_a,2+4*levs+k)),
!    & minval(buff_mult_pieceg(1:lonf,1:lats_node_a,2+4*levs+k))
      enddo
!
      if (levh > 0) then
        do k=1,levh
!$omp parallel do private(i,j)
          do j=1,lats_node_a
            do i=1,lonf
              buffi(i,j) = rqg(i,j,k)
!             if (abs(buffi(i,j)) < 1.0e-15) buffi(i,j) = 0.0
            enddo
          enddo
          CALL uninterpreg(1,kmsk,buffo,buffi,global_lats_a,lonsperlat,
     &                     buff_mult_pieceg(1,1,2+5*levs+k) )

!      write(0,*)'in grid collect, buff_rqg=',' me=',me,
!    & maxval(buff_mult_pieceg(1:lonf,1:lats_node_a,2+5*levs+k)),
!    & minval(buff_mult_pieceg(1:lonf,1:lats_node_a,2+5*levs+k))
        enddo
      endif
!
      return
      end

       subroutine atmgg_move(ioproc)
c
c***********************************************************************
c
      use gfs_dyn_resol_def
      use gfs_dyn_write_state
      use gfs_dyn_layout1
      use gfs_dyn_mpi_def
      implicit none
!
      integer ipt_lats_node_al,nodesr,lats_nodes_al,ioproc
!     integer lats_nodes_a(nodes),ipt,maxfld,ioproc,nproct
      integer proc,j,lat,msgtag,nproc,i,msgtag1,buff,startlat,ierr
      integer illen,nd1,icount
!     integer illen,nd1,icount,maxlats_comp
      data icount/0/
!  allocate the data structures
!
      if(icount == 0) then
         allocate(ivarg_global(10))
         allocate(ivarg_global_a(10,nodes))
         ivarg_global(1) = ipt_lats_node_a
         ivarg_global(2) = lats_node_a
         ivarg_global(3) = lats_node_a_max
         call mpi_gather(ivarg_global,10,MPI_INTEGER,
     &       ivarg_global_a,10,MPI_INTEGER,ioproc,MPI_COMM_ALL,ierr)
         icount = icount + 1
      endif
!     print *,' icount=',icount
!!
      if(.not. allocated(buff_mult_piecesg)) then
!        maxlats_comp = lats_node_a_max
         if(me == ioproc) then
!          maxlats_comp = ivarg_global_a(3,1)
           allocate(buff_mult_piecesg(lonf*latg*ngrids_gg))
        endif
      endif


!
!   SENDLOOP of grids from comp processors to I/O task.  The
!   I/O task may or may not be a comp task also.  The
!   send logic on that task is different for these two cases
!
!  big send
!     if(me .gt. -1) return
!
!
      IF (ME /= ioproc) THEN      !   Sending the data
         msgtag = me
         illen  = lats_node_a
         CALL mpi_send            !  send the local grid domain
     &(buff_mult_pieceg,illen*lonf*ngrids_gg,MPI_R_IO,ioproc,
     &                  msgtag,MPI_COMM_ALL,info)
      ELSE
        if( MC_COMP /= MPI_COMM_NULL) then
!
c iotask is also a compute task.  send is replaced with direct
c  array copy
          if(nodes_comp == 1) then
            buff_mult_piecesg(1:lonf*lats_node_a*ngrids_gg)=
     &     reshape(buff_mult_pieceg(1:lonf,1:lats_node_a,1:ngrids_gg),
     &       (/lonf*lats_node_a*ngrids_gg/) )
!                              END COMPUTE TASKS PORTION OF LOGIC
          else
!
c  END COMPUTE TASKS PORTION OF LOGIC
c  receiving part of I/O task
!
!!      for pes ioproc
            nd1 = 0
            DO proc=1,nodes_comp
              illen = ivarg_global_a(2,proc)
              if (proc /= ioproc+1) then
                msgtag = proc-1
!               print *,' pux target ',ubound(buff_mult_piecesg)
                CALL mpi_recv(buff_mult_piecesg(nd1+1),
     &                        illen*lonf*ngrids_gg
     &                       ,MPI_R_IO,proc-1,
     &                        msgtag,MPI_COMM_ALL,stat,info)
              else
                buff_mult_piecesg(nd1+1:nd1+illen*lonf*ngrids_gg)=
     &          reshape(buff_mult_pieceg(1:lonf,1:illen,1:ngrids_gg),
     &                  (/lonf*illen*ngrids_gg/) )
              endif
              nd1 = nd1 + illen*lonf*ngrids_gg
            enddo
          endif 
        endif
      ENDIF
      call mpi_barrier(mpi_comm_all,ierr)
!!
      return
      end
      SUBROUTINE atmgg_wrt(IOPROC,cfile,xhour,idate
     &,                  global_lats_a,lonsperlat,pdryini)
!!
      use gfs_dyn_resol_def
      use namelist_dynamics_def
      use gfs_dyn_vert_def
      use gfs_dyn_coordinate_def
      use gfs_dyn_io_header
      use gfs_dyn_write_state
      use gfs_dyn_layout1
      use gfs_dyn_mpi_def
      use gfs_dyn_physcons, rk => con_rocp
      use gfs_dyn_tracer_const, only : cpi,ri
      use gfs_dynamics_output, only : DYN_INT_STATE_ISCALAR,
     & DYN_INT_STATE_RSCALAR,DYN_INT_STATE_LSCALAR, DYN_INT_STATE_1D_I,
     & DYN_INT_STATE_2D_I,DYN_INT_STATE_1D_R,DYN_INT_STATE_2D_R,
     & DYN_INT_STATE_3D_R_DIAB,DYN_INT_STATE_3D_R_ADIAB,
     & DYN_INT_STATE_4D_R
      use nemsio_module
!
      implicit none
!!
      real (kind=kind_grid), parameter :: rk1 = rk + 1.0, rkr = 1.0/rk
     &,                                   p0=100000.0, p0i=1.0/p0
      real (kind=kind_io4), parameter  :: zero4=0.0
      integer IOPROC
      character*40 cfile, tracer
      real(kind=kind_grid) xhour, pdryini
      integer idate(4),k,il, ngridgg, nt,idate7(7)
!
      INTEGER, dimension(latg) :: global_lats_a, lonsperlat
!!
!     real(kind=kind_evod) gencode,ppid,realform
      real(kind=kind_io4) yhour, pdryini4
      real(kind=kind_io4), allocatable          :: vcoord4(:,:,:)
      real(kind=kind_io4), dimension(lonf*latg) :: pup,  pdn, plyr
     &,                                            pupk, pdnk
!     real(kind=kind_io8), allocatable :: buff(:)
!     real(kind=kind_io4) yhour, pdryini4, vcoord4(levp1,3)
      integer iret, ks,iorder_l,i
!     integer iret, ks,irealf,iorder, idusr
      character * 16, allocatable :: recname(:), reclevtyp(:)
      integer,       allocatable :: reclev(:)
!
      integer j,ndim3,N2DR,kount,nrec,nfhour,nfminute,nfsecondn
     &,      nfsecondd,nmetavari,nmetavarr,nmetavarl,nmetaaryi,nmetaaryr
      character(16),allocatable :: variname(:),varrname(:),
     &                             varlname(:),aryiname(:),aryrname(:)
      integer,allocatable       :: varival(:),aryilen(:),aryrlen(:),
     &                             aryival(:,:)
      real(kind=kind_io4),allocatable  :: varrval(:),aryrval(:,:)
      logical,allocatable              :: varlval(:)
      REAL(KIND=KIND_io4) ,allocatable :: buff_multg(:,:)
      character * 16                   :: recname1, reclevtyp1
      integer reclev1
!  ---  locals for iau:
      integer :: jdate_iau(4),idat(8),jdat(8)
      real (kind=kind_evod) :: fhour_iau
      real :: rinc(5)

      real(kind=kind_io4),allocatable  :: tmp(:)
      type(nemsio_gfile)               :: gfileout
!
      logical first
      save first, recname, reclevtyp, reclev, vcoord4
      save nrec,nmetavari,nmetavarr,nmetavarl,nmetaaryi,nmetaaryr,
     &     variname,varrname,varlname,aryiname,aryrname,
     &     varival,aryilen,aryrlen,aryival,aryrval,varrval,varlval
      data first /.true./
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!    Build upper air fields in to buff_mult
!
      if (me == ioproc) then
!
        allocate(buff_multg(lonf*latg,ngrids_gg))
        do ngridgg=1,ngrids_gg
          call unsplit2g(ioproc,ngridgg,ngrids_gg,buff_multg(1,ngridgg),
     &        global_lats_a)
        enddo
!     print *,' finished ngrid loop ngrids_gg=',ngrids_gg
!    Building upper air  field is done
!
        if (first) then
          first = .false.
!***
!--- get names from module output 
          nrec=ngrids_gg
          kount=size(DYN_INT_STATE_ISCALAR,2)
!for integer var::
          nmetavari=0
          do i=1,kount
           if(trim(DYN_INT_STATE_ISCALAR(2,i)).eq.'OGFS_SIG') 
     &     nmetavari=nmetavari+1
          enddo
          if(nmetavari>0) then
            allocate(variname(nmetavari),varival(nmetavari))
            j=0
            do i=1,kount
             if(trim(DYN_INT_STATE_ISCALAR(2,i)).eq.'OGFS_SIG')then
               j=j+1
               variname(j)=trim(DYN_INT_STATE_ISCALAR(1,i))
               if(variname(j)=='latg')         varival(j)=latg
               if(variname(j)=='lonf')         varival(j)=lonf
               if(variname(j)=='levs')         varival(j)=levs
               if(variname(j)=='jcap')         varival(j)=jcap
               if(variname(j)=='ntoz')         varival(j)=ntoz
               if(variname(j)=='ntcw')         varival(j)=ntcw
               if(variname(j)=='ncld')         varival(j)=ncld
               if(variname(j)=='ntrac')        varival(j)=ntrac
               if(variname(j)=='vertcoord_id') varival(j)=vertcoord_id
               if(variname(j)=='thermodyn_id') varival(j)=thermodyn_id
               if(variname(j)=='sfcpress_id')  varival(j)=sfcpress_id
               if(variname(j)=='ienst')        varival(j)=ienst
               if(variname(j)=='iensi')        varival(j)=iensi
               if(variname(j)=='itrun')        varival(j)=itrun
               if(variname(j)=='icen2')        varival(j)=icen2
             endif
            enddo
          endif
!for real var::
          nmetavarr=0
          do i=1,kount
           if(trim(DYN_INT_STATE_RSCALAR(2,i)).eq.'OGFS_SIG')
     &     nmetavarr=nmetavarr+1
          enddo
          if(nmetavarr>0) then
            allocate(varrname(nmetavarr),varrval(nmetavarr))
            j=0
            do i=1,kount
             if(trim(DYN_INT_STATE_RSCALAR(2,i)).eq.'OGFS_SIG')then
              j=j+1
              varrname(j)=trim(DYN_INT_STATE_RSCALAR(1,i))
              if(varrname(j)=='pdryini') varrval(j)=pdryini
              if(varrname(j)=='fhour')   varrval(j)=xhour
             endif
            enddo
         endif

!for logical var::
          nmetavarl=0
          do i=1,kount
           if(trim(DYN_INT_STATE_LSCALAR(2,i)).eq.'OGFS_SIG')
     &     nmetavarl=nmetavarl+1
          enddo
          if(nmetavarl>0) then
            allocate(varlname(nmetavarl),varlval(nmetavarl))
            j=0
            do i=1,kount
            if(trim(DYN_INT_STATE_LSCALAR(2,i)).eq.'OGFS_SIG')then
              j=j+1
              varlname(j)=trim(DYN_INT_STATE_LSCALAR(1,i))
              if(varlname(j)=='HYBRID') varlval(j)=hybrid
              if(varlname(j)=='GEN_COORD_HYBRID') 
     &         varlval(j)=gen_coord_hybrid
              if(varlname(j)=='ADIABATIC') varlval(j)=adiabatic
             endif
            enddo
         endif
!for 1D integer array
          nmetaaryi=0
          do i=1,kount
           if(trim(DYN_INT_STATE_1D_I(2,i)).eq.'OGFS_SIG')
     &     nmetaaryi=nmetaaryi+1
          enddo
!  Add fcst_date into ARYI
!jw          nmetaaryi=nmetaaryi+1
          allocate(aryiname(nmetaaryi),aryilen(nmetaaryi))
          j=0
          do i=1,kount
           if(trim(DYN_INT_STATE_1D_I(2,i)).eq.'OGFS_SIG') then
            j=j+1
            aryiname(j)=trim(DYN_INT_STATE_1D_I(1,i))
            if(aryiname(j)=='IDATE') aryilen(j)=size(idate)
           endif
          enddo
!jw          aryiname(nmetaaryi)='FCSTDATE'
!jw          aryilen(nmetaaryi)=7
          allocate(aryival(maxval(aryilen),nmetaaryi))
          aryival(1:aryilen(1),1)=idate(1:aryilen(1))
!jw          aryival(1:7,nmetaaryi)=FCSTDATE(1:7)
!for 1D real array
          nmetaaryr=0
          do i=1,kount
           if(trim(DYN_INT_STATE_1D_R(2,i)).eq.'OGFS_SIG')
     &     nmetaaryr=nmetaaryr+1
          enddo
          if(nmetaaryr>0) then
            allocate(aryrname(nmetaaryr),aryrlen(nmetaaryr))
            j=0
            do i=1,kount
             if(trim(DYN_INT_STATE_1D_R(2,i)).eq.'OGFS_SIG') then
              j=j+1
              aryrname(j)=trim(DYN_INT_STATE_1D_R(1,i))
              if(aryrname(j)=='AK5') aryrlen(j)=size(ak5)
              if(aryrname(j)=='BK5') aryrlen(j)=size(bk5)
              if(aryrname(j)=='CK5') aryrlen(j)=size(ck5)
              if(aryrname(j)=='SI') aryrlen(j)=size(si)
              if(aryrname(j)=='CPI') aryrlen(j)=ntrac+1
              if(aryrname(j)=='RI') aryrlen(j)=ntrac+1
             endif
            enddo
            allocate(aryrval(maxval(aryrlen),nmetaaryr))
            aryrval(1:aryrlen(1),1)=ak5(1:aryrlen(1))
            aryrval(1:aryrlen(2),2)=bk5(1:aryrlen(2))
            aryrval(1:aryrlen(3),3)=ck5(1:aryrlen(3))
            aryrval(1:aryrlen(4),4)=si(1:aryrlen(4))
            aryrval(1:aryrlen(5),5)=cpi(0:aryrlen(5)-1)
            aryrval(1:aryrlen(6),6)=ri(0:aryrlen(6)-1)
          endif
!
!for record name, levtyp and lev          
          allocate (recname(nrec),reclevtyp(nrec),reclev(nrec))
          N2DR=0
          do i=1,kount
           if(trim(DYN_INT_STATE_2D_R(2,i)).eq.'OGFS_SIG') then
            recname(i)=trim(DYN_INT_STATE_2D_R(1,i))
            reclevtyp(i)='sfc'
            reclev(i)=1
            N2DR=N2DR+1
           endif
          enddo
          do j=1,Kount
           if(trim(DYN_INT_STATE_3D_R_DIAB(2,j)).eq.'OGFS_SIG') then
            if(trim(DYN_INT_STATE_3D_R_DIAB(3,j)).eq.'levs') then
             do i=1,levs
               recname(N2DR+1)=trim(DYN_INT_STATE_3D_R_DIAB(1,j))
               reclevtyp(N2DR+1)='mid layer'
               reclev(N2DR+1)=i
               N2DR=N2DR+1
             enddo
            endif
           endif
          enddo

!! generalized tracers:
!  comment out the DYN_INT_STATE_4D_R, the met+chem tracers are included
!  in DYN_INT_STATE_3D_R_DIAB

!          if(ntrac>ntcw) then
!           do j=1,Kount
!            if(trim(DYN_INT_STATE_4D_R(2,j)).eq.'OGFS_SIG') then
!             if(trim(DYN_INT_STATE_4D_R(3,j))=='levs')then
!              NDIM3=levs
!             elseif(trim(DYN_INT_STATE_4D_R(3,j))=='levsp1')then
!              NDIM3=levs+1
!             endif
!             write(tracer,'("_",i2)') j

!             do i=1,NDIM3
!                recname(N2DR+1)=trim(DYN_INT_STATE_4D_R(1,i))//TRACER
!                reclevtyp(N2DR+1)='mid layer'
!                reclev(N2DR+1)=i
!                N2DR=N2DR+1
!             enddo
!            Endif
!           enddo
!          endif
!
          idpp  = 0
          idusr = 0
          idrun = 0
          ALLOCATE(VCOORD4(levs+1,3,2))
          vcoord4=0.
!for output:
!         idvm    = thermodyn_id*10 + sfcpress_id  
          idvm    = 22                                 ! 1:  ln(ps) 2:ps   ! hmhj
                                                       ! 1: Tv, 2: T, 3:Th
          if (gen_coord_hybrid) then                                      ! hmhj
            idvc    = 3
            idsl    = 2    ! idsl=2 for middle of layer                   ! hmhj
            nvcoord = 3
!            allocate (vcoord4(levp1,nvcoord))
            do k=1,levp1                                                  ! hmhj
              vcoord4(k,1,1) = real(ak5(k),4)*1000.                       ! hmhj
              vcoord4(k,2,1) = bk5(k)                                     ! hmhj
              vcoord4(k,3,1) = ck5(k)*1000.                               ! hmhj
            enddo                                                         ! hmhj
            write(0,*)'in atmgg_wrt,vcoord4(1,1)=',vcoord4(1:levp1,1,1)
          else if (hybrid) then
            idvc    = 2                        ! for hybrid vertical coord.
            idsl    = 1    ! idsl=1 for Phillips 
            nvcoord = 2
!            allocate (vcoord4(levp1,nvcoord))
            do k=1,levp1
              vcoord4(k,1,1) = real(ak5(levp1+1-k),4)*1000.
              vcoord4(k,2,1) = bk5(levp1+1-k)
!              print 190,k,vcoord4(k,1),vcoord4(k,2)
190           format('in gfsio k=',i2,'  ak5r4=',f13.6,'  bk5r4=',e13.5)
            enddo
          else
            idvc    = 1    ! for sigma vertical coord. (default)
            idsl    = 1    ! idsl=1 for Phillips 
            nvcoord = 1
            vcoord4(:,1,1) = si (:)
          endif
!end first
        endif
!
      print*,'in grid collect',iau,xhour
      if ( iau .and. xhour >= 6.0 ) then
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
        fhour_iau = xhour - 6.0
      else
        jdate_iau = idate
        fhour_iau = xhour
      endif
 
        pdryini4 = pdryini
        iorder_l = 2
        irealf   = 2
        yhour    = fhour_iau
        idvt    = (ntoz-1) + 10 * (ntcw-1)
        idate7=0
        idate7(1)=jdate_iau(4)
        idate7(2)=jdate_iau(2)
        idate7(3)=jdate_iau(3)
        idate7(4)=jdate_iau(1)
        idate7(7)=100
!
        nfhour=int(yhour)
        nfminute=int((yhour-nfhour)*60.)
        nfsecondn=int(((yhour-nfhour)*3600.-nfminute*60)*100)
        nfsecondd=100
!
!      write(0,*)' calling nemsio_open lonf=',lonf,' latg=',latg
!     &,' idate=',idate,' yhour=',yhour
!
        call nemsio_open(gfileout,trim(cfile),'write',iret=iret,
     &    modelname='GFS',gdatatype='grib',version=ivsupa,
     &    idate=idate7,nrec=nrec,
     &    nfhour=nfhour,nfminute=nfminute,nfsecondn=nfsecondn,
     &    nfsecondd=nfsecondd,
     &    dimx=lonf,dimy=latg,dimz=levs,jcap=jcap,ncldt=ncld,
     &    idsl=idsl,idvc=idvc,idvm=idvm,
     &    ntrac=ntrac,vcoord=vcoord4,cpi=real(cpi,4),ri=real(ri,4),
     &    extrameta=.true.,nmetavari=nmetavari,
     &    nmetavarr=nmetavarr,nmetavarl=nmetavarl,
     &    nmetaaryi=nmetaaryi,nmetaaryr=nmetaaryr,
     &    variname=variname,varival=varival,varrname=varrname,
     &    varrval=varrval,varlname=varlname,varlval=varlval,
     &    aryiname=aryiname,aryilen=aryilen,aryival=aryival,
     &    aryrname=aryrname,aryrlen=aryrlen,aryrval=aryrval,
     &    recname=recname,reclevtyp=reclevtyp,
     &    reclev=reclev)
!
      print *,' after calling nemfio_open iret=',iret
!     if (yhour .gt. 5.99) call mpi_quit(3333)
!
!     print *,' buff_multg=',buff_multg(lonb*latb/2,:)
      allocate(tmp(lonf*latg) )
       do k=1,nrec
         call nemsio_getrechead(gfileout,k,recname1,reclevtyp1,
     &     reclev1,iret)
         tmp(1:lonf*latg)=buff_multg(1:lonf*latg,k)
         call nemsio_writerec(gfileout,k,tmp,iret=iret)
       enddo
       deallocate(tmp)
       deallocate(buff_multg)
!     print *,' return code before closing iret=',iret

      call  nemsio_close(gfileout,iret)
!
!endif ioproc
      endif
!     print *,' return code after closing iret=',iret
!     if (allocated(vcoord4)) deallocate(vcoord4)
!     print *,' after all atmgg writes iret=',iret
      return
      end
!
!

      subroutine uninterpreg(iord,kmsk,f,fi,global_lats_a,lonsperlat, 
     &                       buff_mult)
!!
      use gfs_dyn_resol_def
      use gfs_dyn_layout1
      implicit none
!!
      integer,intent(in)                  :: iord
      integer,intent(in)                  :: kmsk(lonf,lats_node_a)
      integer,intent(in), dimension(latg) :: global_lats_a, lonsperlat
      real(kind=kind_io8),intent(out)     :: f(lonf,lats_node_a)
      real(kind=kind_io8),intent(in)      :: fi(lonf,lats_node_a)

      real(kind=kind_io4),intent(inout)::buff_mult(lonf,lats_node_a_max)
!      real(kind=4) f4(lonf,lats_node_a)
      integer i,j,lons,lat
!!
      do j=1,lats_node_a
         lat  = global_lats_a(ipt_lats_node_a-1+j)
         lons = lonsperlat(lat)
         if(lons /= lonf) then
           call intlon(iord,1,1,lons,lonf,
     &                 kmsk(1,j),fi(1,j),f(1,j))
!          f4(:,j)=fi(:,j)
         else
            f(:,j)=fi(:,j)
!           f4(:,j)=fi(:,j)
         endif
      enddo
!     print *,' ngridg=',ngridg
!$omp parallel do private(i,j)
      do j=1,lats_node_a
        do i=1,lonf
!jw          buff_mult_pieceg(i,ngridg,j) = f (i,j)
          buff_mult(i,j) = f (i,j)
        end do
      end do
!jw      ngridg=ngridg+1
      end subroutine
       subroutine unsplit2g(ioproc,ngridx,ngridt,x,global_lats_a)
c
c***********************************************************************
c
      use gfs_dyn_resol_def
      use gfs_dyn_write_state
      use gfs_dyn_layout1
      use gfs_dyn_mpi_def
      implicit none
!!
      real(kind=kind_io4) x(lonf*latg)
      integer global_lats_a(latg),ipt_lats_node_al,nodesr
      integer lats_nodes_al
      integer maxfld,ioproc,nproct,ngridx,ngridt
      integer proc,j,lat,msgtag,nproc,i,msgtag1,buff,startlat,ierr
      integer ifldu/0/
      save ifldu
      integer illen,nd1,nd2
       character*8 cna
c
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c
!!
!     write(cna,985)600+ngridg
!985   format('fort.',i3)
      X = 0.
!     maxfld=50
      ifldu=ifldu+1
!!
      IF (me == ioproc) THEN

!!     for pes ioproc
        nproct = nodes_comp
!       print *,' NGRIDG=',ngridg,' ifldu=',ifldu,' nproct=',nproct
!for all the other pes
        nd1 = 0
        DO proc=1,nproct
          ipt_lats_node_al = ivarg_global_a(1,proc)
          lats_nodes_al    = ivarg_global_a(2,proc)
          nd2=lats_nodes_al*lonf*(ngridx-1)
          do j=1,lats_nodes_al
            lat = global_lats_a(ipt_lats_node_al-1+j)
            do i=1,lonf
             x(i+(lat-1)*lonf) = buff_mult_piecesg(nd1+nd2+i+(j-1)*lonf)
            enddo
          enddo
          nd1=nd1+lats_nodes_al*lonf*ngridt
        enddo
!
      ENDIF
!!
      return
      end
