      subroutine filtereo
     x     (teme,tee,ye,dime,die,xe,zeme,zee,we,rme,rqe,rte,
     x      dpme,dpe,dpne,
     x      temo,teo,yo,dimo,dio,xo,zemo,zeo,wo,rmo,rqo,rto,
     x      dpmo,dpo,dpno,
     x      filta,ls_node)
!
! program log
! 2011 02 20 : henry juang, add dp for mass_dp and ndsl advection
cc
      use gfs_dyn_resol_def
      use gfs_dyn_layout1
      implicit none
cc
      real(kind=kind_evod)
     1  dpe(len_trie_ls,2,levs),
     1 dpme(len_trie_ls,2,levs),
     1 dpne(len_trie_ls,2,levs)
c
      real(kind=kind_evod)
     1  tee(len_trie_ls,2,levs),  die(len_trie_ls,2,levs),
     1 teme(len_trie_ls,2,levs), dime(len_trie_ls,2,levs),
     1   ye(len_trie_ls,2,levs),   xe(len_trie_ls,2,levs),
     1  rqe(len_trie_ls,2,levh),  rte(len_trie_ls,2,levh)
cc
      real(kind=kind_evod)
     1  zee(len_trie_ls,2,levs),
     1 zeme(len_trie_ls,2,levs),
     1   we(len_trie_ls,2,levs),
     1  rme(len_trie_ls,2,levh)
cc
      real(kind=kind_evod)
     1  dpo(len_trio_ls,2,levs),
     1 dpmo(len_trio_ls,2,levs),
     1 dpno(len_trio_ls,2,levs)
c
      real(kind=kind_evod)
     1  teo(len_trio_ls,2,levs),  dio(len_trio_ls,2,levs),
     1 temo(len_trio_ls,2,levs), dimo(len_trio_ls,2,levs),
     1   yo(len_trio_ls,2,levs),   xo(len_trio_ls,2,levs),
     1  rqo(len_trio_ls,2,levh),  rto(len_trio_ls,2,levh)
cc
      real(kind=kind_evod)
     1  zeo(len_trio_ls,2,levs),
     1 zemo(len_trio_ls,2,levs),
     1   wo(len_trio_ls,2,levs),
     1  rmo(len_trio_ls,2,levh)
cc
      real(kind=kind_evod) filta
cc
      integer              ls_node(ls_dim,3)
cc
cc
      integer              k,l,locl,n
cc
      integer              indev
      integer              indod
      integer              indev1,indev2
      integer              indod1,indod2
      real(kind=kind_evod) filtb
cc
      real(kind=kind_evod) cons0p5,cons1     !constant
cc
      integer              indlsev,jbasev
      integer              indlsod,jbasod
cc
      include 'function2'
cc
cc
      CALL countperf(0,13,0.)
!!
      cons0p5 = 0.5d0                        !constant
      cons1   = 1.d0                         !constant
cc
cc
      filtb = (cons1-filta)*cons0p5          !constant
cc
cc
!$omp parallel do private(k,jbasev,jbasod,locl,l,indev,indod)
!$omp+private(indev1,indev2,indod1,indod2)
      do k=1,levs
cc
         do locl=1,ls_max_node
                 l=ls_node(locl,1)
            jbasev=ls_node(locl,2)
            indev1 = jbasev+(L-L)/2+1
            if (mod(L,2).eq.mod(jcap+1,2)) then
               indev2 = jbasev+(jcap+1-L)/2+1
            else
               indev2 = jbasev+(jcap  -L)/2+1
            endif
            do indev = indev1 , indev2
!hmhj no need to have filter for dp
               dpme(indev,1,k) =         dpe (indev,1,k)
               dpme(indev,2,k) =         dpe (indev,2,k)
!hmhj         dpme(indev,1,k) = filta * dpe (indev,1,k)
!hmhjx                        + filtb *(dpme(indev,1,k)+dpne(indev,1,k))
!hmhj         dpme(indev,2,k) = filta * dpe (indev,2,k)
!hmhjx                        + filtb *(dpme(indev,2,k)+dpne(indev,2,k))
               teme(indev,1,k) = filta * tee (indev,1,k)
     x                         + filtb *(teme(indev,1,k)+ye(indev,1,k))
               teme(indev,2,k) = filta * tee (indev,2,k)
     x                         + filtb *(teme(indev,2,k)+ye(indev,2,k))
               dime(indev,1,k) = filta * die (indev,1,k)
     x                         + filtb *(dime(indev,1,k)+xe(indev,1,k))
               dime(indev,2,k) = filta * die (indev,2,k)
     x                         + filtb *(dime(indev,2,k)+xe(indev,2,k))
               zeme(indev,1,k) = filta * zee (indev,1,k)
     x                         + filtb *(zeme(indev,1,k)+we(indev,1,k))
               zeme(indev,2,k) = filta * zee (indev,2,k)
     x                         + filtb *(zeme(indev,2,k)+we(indev,2,k))
                dpe(indev,1,k) =         dpne(indev,1,k)
                dpe(indev,2,k) =         dpne(indev,2,k)
                tee(indev,1,k) =           ye(indev,1,k)
                tee(indev,2,k) =           ye(indev,2,k)
                die(indev,1,k) =           xe(indev,1,k)
                die(indev,2,k) =           xe(indev,2,k)
                zee(indev,1,k) =           we(indev,1,k)
                zee(indev,2,k) =           we(indev,2,k)
cc
            enddo
         enddo
cc
cc......................................................................
cc
         do locl=1,ls_max_node
                 l=ls_node(locl,1)
            jbasod=ls_node(locl,3)
            indod1 = jbasod+(L+1-L)/2+1
            if (mod(L,2).eq.mod(jcap+1,2)) then
               indod2 = jbasod+(jcap  -L)/2+1
            else
               indod2 = jbasod+(jcap+1-L)/2+1
            endif
            do indod = indod1 , indod2
!hmhj no need to have filter for dp
               dpmo(indod,1,k) =         dpo (indod,1,k)
               dpmo(indod,2,k) =         dpo (indod,2,k)
!hmhj         dpmo(indod,1,k) = filta * dpo (indod,1,k)
!hmhjx                        + filtb *(dpmo(indod,1,k)+dpno(indod,1,k))
!hmhj         dpmo(indod,2,k) = filta * dpo (indod,2,k)
!hmhjx                        + filtb *(dpmo(indod,2,k)+dpno(indod,2,k))
               temo(indod,1,k) = filta * teo (indod,1,k)
     x                         + filtb *(temo(indod,1,k)+yo(indod,1,k))
               temo(indod,2,k) = filta * teo (indod,2,k)
     x                         + filtb *(temo(indod,2,k)+yo(indod,2,k))
               dimo(indod,1,k) = filta * dio (indod,1,k)
     x                         + filtb *(dimo(indod,1,k)+xo(indod,1,k))
               dimo(indod,2,k) = filta * dio (indod,2,k)
     x                         + filtb *(dimo(indod,2,k)+xo(indod,2,k))
               zemo(indod,1,k) = filta * zeo (indod,1,k)
     x                         + filtb *(zemo(indod,1,k)+wo(indod,1,k))
               zemo(indod,2,k) = filta * zeo (indod,2,k)
     x                         + filtb *(zemo(indod,2,k)+wo(indod,2,k))
                dpo(indod,1,k) =         dpno(indod,1,k)
                dpo(indod,2,k) =         dpno(indod,2,k)
                teo(indod,1,k) =           yo(indod,1,k)
                teo(indod,2,k) =           yo(indod,2,k)
                dio(indod,1,k) =           xo(indod,1,k)
                dio(indod,2,k) =           xo(indod,2,k)
                zeo(indod,1,k) =           wo(indod,1,k)
                zeo(indod,2,k) =           wo(indod,2,k)
cc
            enddo
         enddo
cc
      enddo
cc
cc......................................................................
cc
!$omp parallel do private(k,jbasev,jbasod,locl,l,indev,indod)
!$omp+private(indev1,indev2,indod1,indod2)
      do k=1,levh
cc
         do locl=1,ls_max_node
                 l=ls_node(locl,1)
            jbasev=ls_node(locl,2)
            indev1 = jbasev+(L-L)/2+1
            if (mod(L,2).eq.mod(jcap+1,2)) then
               indev2 = jbasev+(jcap+1-L)/2+1
            else
               indev2 = jbasev+(jcap  -L)/2+1
            endif
            do indev = indev1 , indev2
cc
               rme(indev,1,k) = filta * rqe(indev,1,k)
     x                        + filtb *(rme(indev,1,k)+rte(indev,1,k))
               rme(indev,2,k) = filta * rqe(indev,2,k)
     x                        + filtb *(rme(indev,2,k)+rte(indev,2,k))
               rqe(indev,1,k) =         rte(indev,1,k)
               rqe(indev,2,k) =         rte(indev,2,k)
cc
            enddo
         enddo
cc
cc......................................................................
cc
         do locl=1,ls_max_node
                 l=ls_node(locl,1)
            jbasod=ls_node(locl,3)
            indod1 = jbasod+(L+1-L)/2+1
            if (mod(L,2).eq.mod(jcap+1,2)) then
               indod2 = jbasod+(jcap  -L)/2+1
            else
               indod2 = jbasod+(jcap+1-L)/2+1
            endif
            do indod = indod1 , indod2
cc
               rmo(indod,1,k) = filta * rqo(indod,1,k)
     x                        + filtb *(rmo(indod,1,k)+rto(indod,1,k))
               rmo(indod,2,k) = filta * rqo(indod,2,k)
     x                        + filtb *(rmo(indod,2,k)+rto(indod,2,k))
               rqo(indod,1,k) =         rto(indod,1,k)
               rqo(indod,2,k) =         rto(indod,2,k)
cc
            enddo
         enddo
cc
      enddo
cc
      CALL countperf(1,13,0.)
!!
      return
      end
