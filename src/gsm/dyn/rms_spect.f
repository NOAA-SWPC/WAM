      subroutine rms_spect(qe_ls,xe_ls,ye_ls,we_ls,re_ls,
     &                     qo_ls,xo_ls,yo_ls,wo_ls,ro_ls,
     &                     ls_nodes,max_ls_nodes)
!
!    September 24, 2004
!    version of rms_spect with ntrac generalized from 3.
!
!    May 2, 2003
!    version of rms_spect with ntrac re_ls ro_ls.
!
!    April 10, 2003
!    version of rms_spect with gather for each level.
!
!    February 24, 2015 - S Moorthi added gg_tracer option and cleaned up
!

!
      use gfs_dyn_resol_def
      use gfs_dyn_layout1
      use namelist_dynamics_def                                         ! hmhj
      use gfs_dyn_mpi_def
      implicit none
!
      real(kind=kind_evod), dimension(len_trie_ls,2) :: qe_ls 
      real(kind=kind_evod), dimension(len_trie_ls,2,levs) :: xe_ls,ye_ls
     &,                                                      we_ls
      real(kind=kind_evod), dimension(len_trie_ls,2,levs,ntrac) :: re_ls
!
      real(kind=kind_evod), dimension(len_trio_ls,2) :: qo_ls
      real(kind=kind_evod), dimension(len_trio_ls,2,levs) :: xo_ls,yo_ls
     &,                                                      wo_ls
      real(kind=kind_evod), dimension(len_trio_ls,2,levs,ntrac) :: ro_ls
!
!mjr  real(kind=kind_evod) del(levs)
!
      integer              ls_nodes(ls_dim,nodes),  max_ls_nodes(nodes)
!
      integer              ierr,j,k,l,lenrec,lev,locl,n,node
     &,                    indev, indod,indlsev,jbasev,indlsod,jbasod
     &,                    kmq, kmx, kmy, kmw, kmr, ktot, len_tot, kk
!
      real(kind=kind_evod) fgbar(3+ntrac+1), fgbar_sav(3,levs)
     &,                    fgbar_ntrac(ntrac,levs)
!
!mjr  real(kind=kind_evod) vx
!mjr  real(kind=kind_evod) vy
!mjr  real(kind=kind_evod) vw
!mjr  real(kind=kind_evod) vr
!
      include 'function2'
!
      real(kind=kind_mpi),allocatable :: trieo_ls_node (:,:,:)
      real(kind=kind_mpi),allocatable :: trieo_ls_nodes(:,:,:,:)
!
      integer       node_l_all(0:jcap)
      integer     jbasev_l_all(0:jcap)
      integer     jbasod_l_all(0:jcap)
!
      real(kind=kind_evod), parameter ::  cons0=0.d0, cons0p5=0.5d0
!
      kmq = 1      !  qe/o_ls
!
      kmx = 1+kmq  !  xe/o_ls
      kmy = 2+kmq  !  ye/o_ls
      kmw = 3+kmq  !  we/o_ls
      kmr = 4+kmq  !  re/o_ls
      if (gg_tracers) then
        ktot = 3 + kmq
      else
        ktot = 3 + kmq + ntrac
      endif
      len_tot = len_trie_ls_max+len_trio_ls_max
!
      allocate (trieo_ls_node(len_tot,2,ktot))
      trieo_ls_node = 0.0
!
      if ( me == 0 ) then
        allocate (trieo_ls_nodes(len_tot,2,ktot,nodes))
      else
         allocate (trieo_ls_nodes( 2, 2, 2, 2 ))
      endif
!
      do j=1,len_trie_ls
         trieo_ls_node(j,1,kmq) = qe_ls(j,1)
         trieo_ls_node(j,2,kmq) = qe_ls(j,2)
      enddo
!
      do j=1,len_trio_ls
         trieo_ls_node(j+len_trie_ls_max,1,kmq) = qo_ls(j,1)
         trieo_ls_node(j+len_trie_ls_max,2,kmq) = qo_ls(j,2)
      enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do lev=1,levs  !                              begin big lev loop
!
        if ( lev == 2 ) then
!
          kmq = 0      !  qe/o_ls
          kmx = 1+kmq  !  xe/o_ls
          kmy = 2+kmq  !  ye/o_ls
          kmw = 3+kmq  !  we/o_ls
          kmr = 4+kmq  !  re/o_ls
          if (gg_tracers) then
            ktot = 3 + kmq
          else
            ktot = 3 + kmq + ntrac
          endif
!
          if ( me == 0 ) then
            deallocate (trieo_ls_nodes)
              allocate (trieo_ls_nodes (len_tot,2,ktot,nodes))
          endif
!
        endif
!
        do j=1,len_trie_ls
!
          trieo_ls_node(j,1,kmx) = xe_ls(j,1,lev)
          trieo_ls_node(j,2,kmx) = xe_ls(j,2,lev)
!
          trieo_ls_node(j,1,kmy) = ye_ls(j,1,lev)
          trieo_ls_node(j,2,kmy) = ye_ls(j,2,lev)
!
          trieo_ls_node(j,1,kmw) = we_ls(j,1,lev)
          trieo_ls_node(j,2,kmw) = we_ls(j,2,lev)
!
          if (.not. gg_tracers) then
            do k=1,ntrac
              trieo_ls_node(j,1,kmr+k-1) = re_ls(j,1,lev,k)
              trieo_ls_node(j,2,kmr+k-1) = re_ls(j,2,lev,k)
            enddo
          endif
!
        enddo
!
        do j=1,len_trio_ls
!
          trieo_ls_node(j+len_trie_ls_max,1,kmx) = xo_ls(j,1,lev)
          trieo_ls_node(j+len_trie_ls_max,2,kmx) = xo_ls(j,2,lev)
!
          trieo_ls_node(j+len_trie_ls_max,1,kmy) = yo_ls(j,1,lev)
          trieo_ls_node(j+len_trie_ls_max,2,kmy) = yo_ls(j,2,lev)
!
          trieo_ls_node(j+len_trie_ls_max,1,kmw) = wo_ls(j,1,lev)
          trieo_ls_node(j+len_trie_ls_max,2,kmw) = wo_ls(j,2,lev)
!
          if (.not. gg_tracers) then
            do k=1,ntrac
              kk = kmr+k-1
              trieo_ls_node(j+len_trie_ls_max,1,kk) = ro_ls(j,1,lev,k)
              trieo_ls_node(j+len_trie_ls_max,2,kk) = ro_ls(j,2,lev,k)
            enddo
          endif
!
        enddo
!
        lenrec = (len_trie_ls_max+len_trio_ls_max) * 2 * ktot
!
!
        call mpi_gather( trieo_ls_node , lenrec, mpi_r_mpi,
     &                   trieo_ls_nodes, lenrec, mpi_r_mpi,
     &                   0, mc_comp, ierr)
!
!
        if ( me /= 0 ) cycle
!
        do node=1,nodes
!
          jbasev = 0
          jbasod = len_trie_ls_max
!cmr      do l = 0, jcap
          do locl=1,max_ls_nodes(node)
            l = ls_nodes(locl,node)
!
              node_l_all(l) = node
            jbasev_l_all(l) = jbasev
            jbasod_l_all(l) = jbasod
!
            jbasev = jbasev + (jcap+3-l)/2
            jbasod = jbasod + (jcap+2-l)/2
          end do
!
        end do
!
        do k=1,ktot
!
!cmr      fgbar(k) = 0.
          fgbar(k) = cons0
!
          l=0
          node   =   node_l_all(l)
          jbasev = jbasev_l_all(l)
          jbasod = jbasod_l_all(l)
!
          indev  = indlsev(0,l)
          indod  = indlsod(1,l)
          do n=0, jcap
            if(mod(n+l,2) == 0) then
               fgbar(k) = fgbar(k) + trieo_ls_nodes(indev,1,k,node)
     &                             * trieo_ls_nodes(indev,1,k,node)
               indev = indev + 1
            else
               fgbar(k) = fgbar(k) + trieo_ls_nodes(indod,1,k,node)
     &                             * trieo_ls_nodes(indod,1,k,node)
               indod = indod + 1
            endif
          end do
!
          indev = indlsev(0,l)
          indod = indlsod(1,l)
          do n=0, jcap
            if(mod(n+l,2) == 0) then
               fgbar(k) = fgbar(k) + trieo_ls_nodes(indev,2,k,node)
     &                             * trieo_ls_nodes(indev,2,k,node)
               indev = indev + 1
            else
               fgbar(k) = fgbar(k) + trieo_ls_nodes(indod,2,k,node)
     &                             * trieo_ls_nodes(indod,2,k,node)
               indod = indod + 1
            endif
          end do
!
!cmr      fgbar(k)=fgbar(k)*0.5
          fgbar(k)=fgbar(k)*cons0p5
!
          do l=1, jcap
            node   =   node_l_all(l)
            jbasev = jbasev_l_all(l)
            jbasod = jbasod_l_all(l)
!
            indev  = indlsev(l  ,l)
            indod  = indlsod(l+1,l)
            do n=l, jcap
               if(mod(n+l,2) == 0) then
                  fgbar(k) = fgbar(k) + trieo_ls_nodes(indev,1,k,node)
     &                                * trieo_ls_nodes(indev,1,k,node)
                  indev = indev + 1
               else
                  fgbar(k) = fgbar(k) + trieo_ls_nodes(indod,1,k,node)
     &                                * trieo_ls_nodes(indod,1,k,node)
                  indod = indod + 1
               endif
            end do
!
            indev = indlsev(l  ,l)
            indod = indlsod(l+1,l)
            do n=l, jcap
               if(mod(n+l,2) == 0) then
                  fgbar(k) = fgbar(k) + trieo_ls_nodes(indev,2,k,node)
     &                                * trieo_ls_nodes(indev,2,k,node)
                  indev = indev + 1
               else
                  fgbar(k) = fgbar(k) + trieo_ls_nodes(indod,2,k,node)
     &                                * trieo_ls_nodes(indod,2,k,node)
                  indod = indod + 1
               endif
            end do
!
          end do
!
          fgbar(k) = sqrt(fgbar(k))
!
        end do
!
!cmr  vx=0.e0      !constant
!mjr  vx=cons0     !constant
!cmr  vy=0.e0      !constant
!mjr  vy=cons0     !constant
!cmr  vw=0.e0      !constant
!mjr  vw=cons0     !constant
!cmr  vr=0.e0      !constant
!mjr  vr=cons0     !constant
!
!mjr  do k=1,levs
!mjr     vx=vx+fgbar(kmx+k-1)*del(k)
!mjr     vy=vy+fgbar(kmy+k-1)*del(k)
!mjr     vw=vw+fgbar(kmw+k-1)*del(k)
!mjr     vr=vr+fgbar(kmr+k-1)*del(k)
!mjr  end do
!
        if ( lev == 1 ) then
!
          if ( gen_coord_hybrid ) then					! hmhj
            print 101,fgbar(kmq)						! hmhj
  101       format(/ 1x,3x, ' rms_spect   ps= ',1(es17.10,1x) /)		! hmhj
          else								! hmhj
            print 100,fgbar(kmq)
  100       format(/ 1x,3x, ' rms_spect   rms_ln(ps)= ',1(es17.10,1x) /)
          endif								! hmhj
!
          print 50
   50     format(1x,3x,
     &           ' rms_spect      div',
     &           '               vort',
     &           '               temp')
!
        endif
!
        fgbar_sav(1,lev) = fgbar(kmx)
        fgbar_sav(2,lev) = fgbar(kmw)
        fgbar_sav(3,lev) = fgbar(kmy)
!mjr  fgbar_sav(4,lev) = fgbar(kmr)
!mjr  fgbar_sav(5,lev) = fgbar(kmr+1)
!mjr  fgbar_sav(6,lev) = fgbar(kmr+2)
!
        if (.not. gg_tracers) then
          do k=1,ntrac
            fgbar_ntrac(k,lev) = fgbar(kmr+k-1)
          enddo
        endif
!
        print 200,lev, fgbar(kmx), fgbar(kmw), fgbar(kmy)
!cmr &          , ( fgbar(kmr+k-1), k=1,ntrac )
  200   format(1x,i3,6(2x,es17.10))
!
      enddo                               ! end big lev loop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (.not. gg_tracers) then
        if ( me == 0 ) then               ! begin ntrac print
!
          print 52
   52     format(/ 1x,3x,
     &           ' rms_spect mixratio',
     &           '              ozone',
     &           '           cld_liqo')
! 
          do lev=1,levs
            print 200, lev, (fgbar_ntrac(k,lev),
     &                         k=1,min(3,ntrac))
          enddo
!
!----------------------------------------------------------------------
!
!         do j=4,ntrac,4
!
!           print 53, (k, ntrac, k=j,min(j+3,ntrac))
!  53       format(/ 1x,3x, 6(2x,i3,' of ntrac =',i3:))
!
!           do lev=1,levs
!              print 200, lev, (fgbar_ntrac(k,lev),
!    &                            k=j,min(j+3,ntrac))
!           enddo
!         enddo
!
!----------------------------------------------------------------------
!
          print 54
   54     format(/ ' fin rms_spect ' /)
!
        endif  !  end ntrac print
      endif
!
      deallocate ( trieo_ls_node  )
      deallocate ( trieo_ls_nodes )
!
      return
      end
