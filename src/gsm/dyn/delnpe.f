      subroutine delnpe(qe,odphi,edlam,epse,epso,ls_node)
!
! 02-11-2015 s moorthi - some cleanup
!
      use gfs_dyn_resol_def
      use gfs_dyn_layout1
      use gfs_dyn_physcons, rerth => con_rerth
      implicit none
!
!    input q is in ibm triang. order
!    output  is in ibm triang. order
!
      real(kind=kind_evod), dimension(len_trie_ls,2) :: qe, edlam
      real(kind=kind_evod), dimension(len_trio_ls,2) :: odphi
!
      real(kind=kind_evod)  epse(len_trie_ls), epso(len_trio_ls)
!
      integer              ls_node(ls_dim,3)
!
!cmr  ls_node(1,1) ... ls_node(ls_max_node,1) : values of L
!cmr  ls_node(1,2) ... ls_node(ls_max_node,2) : values of jbasev
!cmr  ls_node(1,3) ... ls_node(ls_max_node,3) : values of jbasod
!
      real(kind=kind_evod) aa,r1mn,rl,rnp2
      real(kind=kind_evod), parameter :: cons1=1.0d0, cons2=2.0d0
!
      integer              indlsev, jbasev, indlsod, jbasod, l, locl, n
     &,                    indev,indev1,indev2, indod,indod1,indod2
     &,                    inddif, ii
!
      include 'function2'
!
!......................................................................
!
!     CALL countperf(0,13,0.)
!!
      do locl=1,ls_max_node
              l=ls_node(locl,1)
         jbasev=ls_node(locl,2)
         indev1 = indlsev(L,L)
         if (mod(L,2) == mod(jcap+1,2)) then
            indev2 = indlsev(jcap+1,L)
         else
            indev2 = indlsev(jcap  ,L)
         endif
         rl=l
         do indev = indev1 , indev2
!          dlam(l,n)= i*l*q(l,n)
!
            edlam(indev,1) = -rl * qe(indev,2)
            edlam(indev,2) =  rl * qe(indev,1)
!
         enddo
!
      enddo
!
!......................................................................
!
      do locl=1,ls_max_node
             l=ls_node(locl,1)
        jbasev=ls_node(locl,2)
        jbasod=ls_node(locl,3)
        indev1 = indlsev(L,L)
        if (mod(L,2) == mod(jcap+1,2)) then
          indev2 = indlsev(jcap-1,L)
        else
          indev2 = indlsev(jcap  ,L)
        endif
        indod1 = indlsod(l+1,l)
        inddif = indev1 - indod1
!
        r1mn = -l
        do indev = indev1 , indev2
           ii = indev-inddif
           odphi(ii,1) = r1mn * epso(ii) * qe(indev,1)
           odphi(ii,2) = r1mn * epso(ii) * qe(indev,2)
           r1mn = r1mn - cons2
        enddo
!
      enddo
!
!......................................................................
!
      do locl=1,ls_max_node
             l=ls_node(locl,1)
        jbasev=ls_node(locl,2)
        jbasod=ls_node(locl,3)
        indev1 = indlsev(L,L) + 1
        if (mod(L,2) == mod(jcap+1,2)) then
          indev2 = indlsev(jcap-1,L)
        else
          indev2 = indlsev(jcap  ,L)
        endif
        indod1 = indlsod(l+1,l)
        inddif = indev1 - indod1
!
        rnp2 = l+3
        do indev = indev1 , indev2
          ii = indev-inddif
          odphi(ii,1) = odphi(ii,1) + rnp2 * epse(indev) * qe(indev,1)
          odphi(ii,2) = odphi(ii,2) + rnp2 * epse(indev) * qe(indev,2)
          rnp2 = rnp2 + cons2
        enddo
!
      enddo
!
!......................................................................
!
      aa = cons1 / rerth
!
      do locl=1,ls_max_node
             l=ls_node(locl,1)
        jbasev=ls_node(locl,2)
        jbasod=ls_node(locl,3)
        indev1 = indlsev(L,L)
        indod1 = indlsod(L+1,L)
        if (mod(L,2) == mod(jcap+1,2)) then
          indev2 = indlsev(jcap+1,L)
          indod2 = indlsod(jcap  ,L)
        else
          indev2 = indlsev(jcap  ,L)
          indod2 = indlsod(jcap+1,L)
        endif

        DO indev = indev1 , indev2
          edlam(indev,1) = edlam(indev,1) * aa
          edlam(indev,2) = edlam(indev,2) * aa
        enddo
!
        do indod = indod1 , indod2
          odphi(indod,1) = odphi(indod,1) * aa
          odphi(indod,2) = odphi(indod,2) * aa
        enddo
!
      enddo
!
!     CALL countperf(1,13,0.)
!!
      return
      end
