      module reduce_lons_grid_module
      use machine , only : kind_dbl_prec
      implicit none

      real(kind=kind_dbl_prec) , allocatable :: colrad(:)
      real(kind=kind_dbl_prec) , allocatable :: dx(:),yy(:)
      real(kind=kind_dbl_prec) , allocatable :: deps(:),rdeps(:)
      integer , allocatable :: indxmv(:)

      integer jcap,jcap1,jcap2,lnt2,lnut2,latrh

      contains
! -------------------------------------------------------------------   
      subroutine reduce_grid (numreduce,jcapi,latr,lons_lat)
      use machine , only : kind_dbl_prec
      implicit none
      integer numreduce,jcapi,latr
      integer lons_lat(latr)

      real(kind=kind_dbl_prec) , allocatable :: qtt(:)
      real(kind=kind_dbl_prec) qmaxall,qttcut
      integer lat,n

      jcap=jcapi
      latrh=latr/2
      jcap1=jcap+1
      jcap2=jcap+2
      lnt2=jcap1*jcap2
      lnut2=lnt2+jcap1*2
c
      allocate(   qtt(   lnt2) )
      allocate(colrad(  latrh) )
      allocate(  deps(  lnut2) )
      allocate( rdeps(  lnut2) )
      allocate(    dx(jcap1*2) )
      allocate(    yy(  jcap1) )
      allocate(indxmv(  lnut2) )

      if( numreduce.gt.0 ) then
        call reduce_glats()
        call reduce_gpln2i()
        qmaxall=0.0
        do lat=1,latrh
          call reduce_pln2i (lat,qtt)
          do n=1,lnt2
            qmaxall=max(qmaxall,abs(qtt(n)))
          enddo
        enddo
! turn on the safe limit
!       numreduce=max(numreduce,2)
        qttcut=qmaxall/(10.**numreduce)
        do lat=1,latrh
          call reduce_pln2i (lat,qtt)
          call reduce_per_lat (qtt,qttcut,lons_lat(lat))
        enddo
        do lat=1,latrh
          lons_lat(latr-lat+1)=lons_lat(lat)
        enddo
      endif

      deallocate(   qtt)
      deallocate(colrad)
      deallocate(  deps)
      deallocate( rdeps)
      deallocate(    dx)
      deallocate(    yy)
      deallocate(indxmv)
      return
      end subroutine reduce_grid

! -------------------------------------------------------------
      subroutine reduce_per_lat(qtt,qttcut,lonrd)
      use machine , only : kind_dbl_prec
      implicit none
      real(kind=kind_dbl_prec) qtt(lnt2),qttcut
      integer lonrd
c
      integer ind,mwave,need,m,n,ireal,imagi
      integer lcapd,lonri,lonro,lonr,lonrf
      integer ii,j,k,l,jtime,ktime,ltime

      ind=0
      mwave=0
      do m=0,jcap
        need=0
        do n=m,jcap
          ind=ind+1
          ireal=2*ind-1
          imagi=2*ind
          if(abs(qtt(ireal)).ge.qttcut.and.
     &       abs(qtt(imagi)).ge.qttcut) then
             need=1
          endif
        enddo
        mwave=mwave+need
        if( need.eq.0 ) go to 123
      enddo
 123  lcapd=mwave
c
      lonri=3*(lcapd-1)+1
      lonrf=lonri+mod(lonri,2)
 100  continue
        lonr=lonrf
        lonro=lonr
        jtime=nint(log(real(lonr))/log(2.))
        jtime=min(jtime,25)			! consider DCRFT
        do j=1,jtime
          if( mod(lonr,2).eq.0 ) then
            lonr=lonr/2
            if( lonr.eq.1 ) go to 200
          endif
        enddo
        ktime=nint(log(real(lonr))/log(3.))
        ktime=min(ktime,2)			! consider DCRFT
        do k=1,ktime
          if( mod(lonr,3).eq.0 ) then
            lonr=lonr/3
            if( lonr.eq.1 ) go to 200
          endif
        enddo
        ltime=nint(log(real(lonr))/log(5.))
        ltime=min(ltime,1)			! consider DCRFT
        do l=1,ltime
          if( mod(lonr,5).eq.0 ) then
            lonr=lonr/5
            if( lonr.eq.1 ) go to 200
          endif
        enddo
        if( mod(lonr,7).eq.0 ) then		! consider DCRFT
          lonr=lonr/7				! consider DCRFT
          if( lonr.eq.1 ) go to 200		! consider DCRFT
        endif					! consider DCRFT
        if( mod(lonr,11).eq.0 ) then		! consider DCRFT
          lonr=lonr/11				! consider DCRFT
          if( lonr.eq.1 ) go to 200		! consider DCRFT
        endif					! consider DCRFT
      lonrf=lonrf+2
      go to 100
 200  lonrd=lonro
! give more points at polar area
      if( jcap.ge.190 ) then
        lonrd=max(lonrd,64)
      else
        lonrd=max(lonrd,30)
      endif
c
      return
      end subroutine reduce_per_lat
! -------------------------------------------------------------
      subroutine reduce_pln2i(lat,qtt)
      use machine , only : kind_dbl_prec
      implicit none
      real(kind=kind_dbl_prec)  qtt(jcap1*jcap2)
      integer   lat
cc
      real(kind=kind_dbl_prec) x(jcap1),dpln(lnut2),qlnv(lnut2)
      real(kind=kind_dbl_prec) sinlat,cos2,prod
      integer twoj1,ll,lplus,lpt,lpv,lp2,lp1,lp0,len,n,i
cc
      twoj1=2*jcap1
c
      sinlat = cos(colrad(lat))
      cos2   = 1.0 - sinlat * sinlat
      prod   = 1.0
      do ll=1,jcap1
        x(ll) = 0.5*prod
        prod = prod*cos2*yy(ll)
      enddo
      do ll=1,jcap1
        x(ll) = sqrt(x(ll))
      enddo
      do ll=1,jcap1
        dpln(2*ll-1) = x(ll)
        dpln(2*ll  ) = x(ll)
      enddo
c
      lplus = twoj1
      do ll=1,twoj1
        dpln(ll+lplus) = dx(ll) * sinlat * dpln(ll)
      enddo
      lp2 = 0
      lp1 =     twoj1
      lp0 = 2 * twoj1
      len =     twoj1 - 2
      do n=3,jcap2
        do ll=1,len
             dpln(ll+lp0) = (sinlat * dpln(ll+lp1)
     1               - deps(ll+lp1) * dpln(ll+lp2)) * rdeps(ll+lp0)
        enddo
        lp2 = lp1
        lp1 = lp0
        lp0 = lp0 + len
        len = len - 2
      enddo
cc
cc    transpose vector dpln array from cray order to ibm order.
      do i=1,lnut2
        qlnv(indxmv(i)) = dpln(i)
      enddo
cc
      lpv = 0
      lpt = 0
      len = twoj1
      do n=1,jcap1
        do ll=1,len
          qtt(ll+lpt) = qlnv(ll+lpv)
        enddo
        lpv = lpv + len + 2
        lpt = lpt + len
        len = len - 2
      enddo
c
      return
      end subroutine reduce_pln2i
! -------------------------------------------------------------
      subroutine reduce_gpln2i()
      use machine , only : kind_dbl_prec
      implicit none
c
      real(kind=kind_dbl_prec) , allocatable ::  x(:)
      integer ll,lplus,len,inde,l,n,nn
      integer lnut,i,twoj1,ibegin
      integer lln,indx
cc
      allocate( x(jcap1) )
c
      do ll=1,jcap1
        rdeps(ll) = 0.0
      enddo
      lplus = jcap1
      len   = jcap1
      do inde=2,jcap2
        do ll=1,len
          l = ll - 1
          n = l + inde - 1
          rdeps(ll+lplus) = (n*n - l*l) / (4.0 * n*n - 1.0)
        enddo
        lplus = lplus + len
        len = len - 1
      enddo
c
      lnut=jcap1*jcap2/2+jcap1
      do i=jcap2,lnut
        rdeps(i) = sqrt(rdeps(i))
      enddo
      do i=1,lnut
        deps(2*i-1) = rdeps(i)
        deps(2*i  ) = rdeps(i)
      enddo
c
      twoj1=jcap1*2
      ibegin = twoj1 + 1
      do i=ibegin,lnut2
        rdeps(i) = 1.0/deps(i)
      enddo
      do i=1,jcap1
        x(i) = i*2+1
      enddo
      do i=1,jcap1
        yy(i) = x(i)/(x(i)-1.)
      enddo
      do i=1,jcap1
        x(i) = sqrt(x(i))
      enddo
      do i=1,jcap1
        dx(2*i-1) = x(i)
        dx(2*i  ) = x(i)
      enddo
cc
cc    set index array for transposing vector array
cc    from cray order to ibm order.
      l=0
      do nn=1,jcap2
        lln=min0(jcap2-nn+1,jcap1)
        do ll=1,lln
          indx=((jcap1+2)*(ll-1)-(ll-1)*ll/2+nn)*2
          l=l+2
          indxmv(l-1)=indx-1
          indxmv(l  )=indx
        enddo
      enddo
c
      deallocate( x )
cc
      return
      end subroutine reduce_gpln2i
! -------------------------------------------------------------
      subroutine reduce_glats()
      use machine , only : kind_dbl_prec
      implicit none
c
      real(kind=kind_dbl_prec) eps,si,rl2,scale,dradz,rad,drad,p2,p1
      real(kind=kind_dbl_prec) x,w,sn,rc,pi
      integer l2,k1,k,iter
c
      eps=1.e-12
      si = 1.0
      l2=2*latrh
      rl2=l2
      scale = 2.0/(rl2*rl2)
      k1=l2-1
      pi = atan(si)*4.0 
      dradz = pi / 360. / 10.
      rad = 0.0
      do 1000 k=1,latrh
        iter=0
        drad=dradz
1       call reduce_poly(l2,rad,p2)
2       p1 =p2
        iter=iter+1
        rad=rad+drad
        call reduce_poly(l2,rad,p2)
        if(sign(si,p1).eq.sign(si,p2)) go to 2
        if(drad.lt.eps)go to 3
        rad=rad-drad
        drad = drad * 0.25
        go to 1
3       continue
        colrad(k)=rad
!       print *,' k colat ',k,colrad(k)
1000  continue
      return
      end subroutine reduce_glats
! ------------------------------------------------------------------
      subroutine reduce_poly(n,rad,p)
      use machine , only : kind_dbl_prec
      implicit none
      integer n
      real(kind=kind_dbl_prec) rad,p
      real(kind=kind_dbl_prec) x,y1,y2,g,y3
      integer i
      x = cos(rad)
      y1 = 1.0
      y2=x
      do i=2,n
      g=x*y2
      y3=g-y1+g-(g-y1)/float(i)
      y1=y2
      y2=y3
      enddo
      p=y3
      return
      end subroutine reduce_poly

      end module reduce_lons_grid_module
