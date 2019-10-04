      subroutine triseori(trisca,triev,triod,levels,ls_node)
!
!02-11-2015  S. Moorthi - some cleanup
! 
      use gfs_dyn_resol_def
      use gfs_dyn_layout1
      implicit none
      integer              levels

!     real(kind=kind_evod) trisca(lnt22,levels)
      real(kind=kind_io8)   trisca(lnt2,levels)
      real(kind=kind_evod)  triev(len_trie_ls,2,levels)
      real(kind=kind_evod)  triod(len_trio_ls,2,levels)
!
      integer              ls_node(ls_dim,3)
!
!cmr  ls_node(1,1) ... ls_node(ls_max_node,1) : values of L
!cmr  ls_node(1,2) ... ls_node(ls_max_node,2) : values of jbasev
!cmr  ls_node(1,3) ... ls_node(ls_max_node,3) : values of jbasod
cc
!cmr    copy real even elements of scalar complex even-odd triangles
!cmr    to   real even triangles.  set top rows to zeros.
!cmr    copy imaginary even elements of scalar complex even-odd triangles
!cmr    to   imaginary even triangles.  set top rows to zeros.
!cmr    copy real odd elements of scalar complex even-odd triangles
!cmr    to   real odd triangles.  set top rows to zeros.
!cmr    copy imaginary odd elements of scalar complex even-odd triangles
!cmr    to   imaginary odd triangles.  set top rows to zeros.
!
!    local scalars
!    -------------
!
      integer   indsca, indev, indod, indev1,indev2, indod1, indod2
     &,         n, l, locl, k
!
!    statement functions
!    -------------------
!
!    offsca is scalar complex even-odd triangle offset in words.
!
      integer   offsca
!
      offsca(n,l) = (jcap+1)*(jcap+2) - (jcap-l+1)*(jcap-l+2) + 2*(n-l)
!
      integer   indlsev, jbasev, indlsod, jbasod
!
      include 'function2'
!
!
      real(kind=kind_evod), parameter :: cons0=0.d0
!
!......................................................................
!
      do k = 1, levels
!
         do locl=1,ls_max_node
                 l=ls_node(locl,1)
            jbasev=ls_node(locl,2)
            indev1 = indlsev(L,L)
            if (mod(L,2) == mod(jcap+1,2)) then
               indev2 = indlsev(jcap-1,L)
            else
               indev2 = indlsev(jcap  ,L)
            endif
!          copy the even (n-l) terms for each level
            indsca = offsca(l,l)
            do indev = indev1 , indev2
               triev(indev,1,k) = trisca(indsca+1,k) ! real part
               triev(indev,2,k) = trisca(indsca+2,k) ! imaginary part
               indsca = indsca+4
            end do
         end do
!
!......................................................................
!
         do locl=1,ls_max_node
                 l=ls_node(locl,1)
           if (l /= jcap) then
              jbasod = ls_node(locl,3)
              indod1 = indlsod(L+1,L)
              if (mod(L,2) == mod(jcap+1,2)) then
                 indod2 = indlsod(jcap  ,L)
              else
                 indod2 = indlsod(jcap-1,L)
              endif
!          copy the odd (n-l) terms for each level
              indsca   = offsca(l+1,l)
              do indod = indod1 , indod2
                 triod(indod,1,k) = trisca(indsca+1,k) ! real part
                 triod(indod,2,k) = trisca(indsca+2,k) ! imaginary part
                 indsca = indsca+4
              end do
 
           endif
         end do
!
!......................................................................
!
         do locl=1,ls_max_node
                    l = ls_node(locl,1)
               jbasev = ls_node(locl,2)
               jbasod = ls_node(locl,3)
            if (mod(L,2) == mod(jcap+1,2)) then
!             set the even (n-l) terms of the top row to zero
               triev(indlsev(jcap+1,l),1,k) = cons0 ! real part
               triev(indlsev(jcap+1,l),2,k) = cons0 ! imaginary part
            else
!             set the  odd (n-l) terms of the top row to zero
               triod(indlsod(jcap+1,l),1,k) = cons0 ! real part
               triod(indlsod(jcap+1,l),2,k) = cons0 ! imaginary part
            endif
         end do
!
      end do
!
      return
      end
