      subroutine matinv(a,m,n,d,p,r)
!
! program log
! 2011 02 20 : henry juang, need for mass_dp and ndsl advection
!
      use gfs_dyn_machine
      implicit none
      integer              i,j,k,l,m,n
      real(kind=kind_evod) a(m,n,n),d(m),p(m),r(m)
      real(kind=kind_evod) cons0,cons1     !constant
      cons0 = 0.d0     !constant
      cons1 = 1.d0     !constant
      do 200 l=1,m
      d(l)=cons1     !constant
  200 continue
      do 100 k=1,n
      do 250 l=1,m
      p(l)=a(l,k,k)
  250 continue
      do 300 l=1,m
      r(l)=-cons1/p(l)     !constant
  300 continue
      do 350 l=1,m
      a(l,k,k)=cons0       !constant
  350 continue
      do  20 i=1,n
      do 400 l=1,m
      a(l,i,k)=a(l,i,k)*r(l)
  400 continue
   20 continue
      do 60 i=1,n
      if(i.eq.k) go to 60
      do  40 j=1,n
      do 450 l=1,m
      a(l,i,j)=a(l,i,k)*a(l,k,j)+a(l,i,j)
  450 continue
   40 continue
   60 continue
      do 600 l=1,m
      r(l)=-r(l)
  600 continue
      do  80 j=1,n
      do 650 l=1,m
      a(l,k,j)=a(l,k,j)*r(l)
  650 continue
   80 continue
      do 700 l=1,m
      d(l)=d(l)*p(l)
  700 continue
      do 750 l=1,m
      a(l,k,k)=r(l)
  750 continue
  100 continue
      return
      end
