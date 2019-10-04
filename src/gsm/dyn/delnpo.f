      subroutine delnpo(qo,edphi,odlam,epse,epso,ls_node)
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
      real(kind=kind_evod), dimension(len_trio_ls,2) :: qo, odlam
      real(kind=kind_evod), dimension(len_trie_ls,2) :: edphi
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
!
      real(kind=kind_evod), parameter :: cons0=0.0d0, cons1=1.0d0,
     &                                   cons2=2.0d0
!
      integer              indlsev,jbasev,indlsod,jbasod,l,locl,n
     &,                    indev,indev1,indev2,indod,indod1,indod2
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
        jbasod=ls_node(locl,3)
        indod1 = indlsod(L+1,L)
        if (mod(L,2) == mod(jcap+1,2)) then
          indod2 = indlsod(jcap  ,L)
        else
          indod2 = indlsod(jcap+1,L)
        endif
!
        rl = l
        do indod = indod1 , indod2
!         dlam(l,n)= i*l*q(l,n)
          odlam(indod,1) = -rl * qo(indod,2)
          odlam(indod,2) =  rl * qo(indod,1)
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
          indev2 = indlsev(jcap+1,L)
        else
          indev2 = indlsev(jcap  ,L)
        endif
        indod1 = indlsod(l+1,l)
        inddif = indev1 - indod1
!
        r1mn = -l - 1
        do indev = indev1 , indev2
          ii = indev - inddif
          edphi(indev,1) = r1mn * epse(indev) * qo(ii,1)
          edphi(indev,2) = r1mn * epse(indev) * qo(ii,2)
          r1mn           = r1mn - cons2
        enddo
!
      enddo
!
!......................................................................
!
      do locl=1,ls_max_node
             l = ls_node(locl,1)
!       jbasev = ls_node(locl,2)
        edphi(indlsev(l,l),1) = cons0
        edphi(indlsev(l,l),2) = cons0
      enddo
!
!......................................................................
!
      do locl=1,ls_max_node
             l = ls_node(locl,1)
        jbasev = ls_node(locl,2)
        jbasod = ls_node(locl,3)
        indev1 = indlsev(L,L)
        if (mod(L,2) == mod(jcap+1,2)) then
          indev2 = indlsev(jcap+1,L) - 1
        else
          indev2 = indlsev(jcap  ,L) - 1
        endif
        indod1 = indlsod(l+1,l)
        inddif = indev1 - indod1
!
        rnp2 = l+2
        do indev = indev1 , indev2
          ii = indev - inddif
          edphi(indev,1) = edphi(indev,1) + rnp2 * epso(ii) * qo(ii,1)
          edphi(indev,2) = edphi(indev,2) + rnp2 * epso(ii) * qo(ii,2)
          rnp2           = rnp2 + cons2
        enddo
!
      enddo
!
!......................................................................
!
      aa = cons1 / rerth
!
!
      do locl=1,ls_max_node
             l = ls_node(locl,1)
        jbasev = ls_node(locl,2)
        jbasod = ls_node(locl,3)
        indev1 = indlsev(L,L)
        indod1 = indlsod(L+1,L)
        if (mod(L,2) == mod(jcap+1,2)) then
          indev2 = indlsev(jcap+1,L)
          indod2 = indlsod(jcap  ,L)
        else
          indev2 = indlsev(jcap  ,L)
          indod2 = indlsod(jcap+1,L)
        endif

        do indod = indod1 , indod2
          odlam(indod,1) = odlam(indod,1) * aa
          odlam(indod,2) = odlam(indod,2) * aa
        enddo
!
        do indev = indev1 , indev2
          edphi(indev,1) = edphi(indev,1) * aa
          edphi(indev,2) = edphi(indev,2) * aa
        enddo
!
      enddo
!
!     CALL countperf(1,13,0.)
!!
      return
      end
