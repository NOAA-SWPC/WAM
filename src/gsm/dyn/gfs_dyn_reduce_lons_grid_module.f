      module gfs_dyn_reduce_lons_grid_module
!
! This moudule is based on Juang 2004 paper about reduced spherical
! transform which has both reduced Gaussian grid and reduced Legendre
! transform
!
! program lot:
! 20110211:	Henry Juang upgraded to have linear grid option for 
!               non-iterating dimensional-split semi-Lagrangian NDSL
! 20130202:	Henry Juang add x number and increse precision for glats
!

      use gfs_dyn_machine , only : kind_dbl_prec, kind_qdt_prec
      implicit none

      real(kind=kind_dbl_prec) , allocatable :: colrad(:)
      real(kind=kind_dbl_prec) , allocatable :: dx(:),yy(:)
      real(kind=kind_dbl_prec) , allocatable :: deps(:),rdeps(:)
      integer , allocatable :: indxmv(:)

      integer jcap1,jcap2,lnt2,lnut2,latgh

      contains
! -------------------------------------------------------------------   
      subroutine gfs_dyn_reduce_grid (numreduce,jcap,lonf,latg,lons_lat)
      implicit none
      integer numreduce,jcap,lonf,latg
      integer lons_lat(latg)

      real(kind=kind_dbl_prec) , allocatable :: qtt(:)
      real(kind=kind_dbl_prec) qmaxall,qttcut
      integer lat,n,is

      latgh=latg/2
      jcap1=jcap+1
      jcap2=jcap+2
      lnt2=jcap1*jcap2
      lnut2=lnt2+jcap1*2
c
      allocate(   qtt(   lnt2) )
      allocate(colrad(  latgh) )
      allocate(  deps(  lnut2) )
      allocate( rdeps(  lnut2) )
      allocate(    dx(jcap1*2) )
      allocate(    yy(  jcap1) )
      allocate(indxmv(  lnut2) )

      if( numreduce.gt.0 ) then
        call gfs_dyn_reduce_glats()
        call gfs_dyn_reduce_gpln2i()
        qmaxall=0.0
        do lat=1,latgh
          call gfs_dyn_reduce_pln2i (lat,qtt)
          do n=1,lnt2
            qmaxall=max(qmaxall,abs(qtt(n)))
          enddo
        enddo
! turn on the safe limit
        numreduce=max(numreduce,4)
        qttcut=qmaxall/(10.**numreduce)
        is = nint( (lonf-1.)/jcap )
        do lat=1,latgh
          call gfs_dyn_reduce_pln2i (lat,qtt)
          call gfs_dyn_reduce_per_lat (qtt,qttcut,is,lons_lat(lat))
        enddo
        do lat=1,latgh
          lons_lat(latg-lat+1)=lons_lat(lat)
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
      end subroutine gfs_dyn_reduce_grid

! -------------------------------------------------------------
      subroutine gfs_dyn_reduce_per_lat(qtt,qttcut,is,lonfd)
      implicit none
      real(kind=kind_dbl_prec) qtt(lnt2),qttcut
      integer lonfd,is
c
      integer ind,mwave,need,m,n,ireal,imagi
      integer lcapd,lonfi,lonfo,lonf,lonff
      integer ii,j,k,l,jtime,ktime,ltime

      ind=0
      mwave=0
      do m=1,jcap1
        need=0
        do n=m,jcap1
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
      lonfi=is*(lcapd-1)+1
      lonff=lonfi+mod(lonfi,2)
 100  continue
        lonf=lonff
        lonfo=lonf
        jtime=nint(log(real(lonf))/log(2.))
        jtime=min(jtime,25)			! consider DCRFT
        do j=1,jtime
          if( mod(lonf,2).eq.0 ) then
            lonf=lonf/2
            if( lonf.eq.1 ) go to 200
          endif
        enddo
        ktime=nint(log(real(lonf))/log(3.))
        ktime=min(ktime,2)			! consider DCRFT
        do k=1,ktime
          if( mod(lonf,3).eq.0 ) then
            lonf=lonf/3
            if( lonf.eq.1 ) go to 200
          endif
        enddo
        ltime=nint(log(real(lonf))/log(5.))
        ltime=min(ltime,1)			! consider DCRFT
        do l=1,ltime
          if( mod(lonf,5).eq.0 ) then
            lonf=lonf/5
            if( lonf.eq.1 ) go to 200
          endif
        enddo
        if( mod(lonf,7).eq.0 ) then		! consider DCRFT
          lonf=lonf/7				! consider DCRFT
          if( lonf.eq.1 ) go to 200		! consider DCRFT
        endif					! consider DCRFT
        if( mod(lonf,11).eq.0 ) then		! consider DCRFT
          lonf=lonf/11				! consider DCRFT
          if( lonf.eq.1 ) go to 200		! consider DCRFT
        endif					! consider DCRFT
      lonff=lonff+2
      go to 100
 200  lonfd=lonfo
!hmhj give more points at polar area
!hmhj if( jcap.ge.190 ) then
!hmhj   lonfd=max(lonfd,64)
!hmhj else
!hmhj   lonfd=max(lonfd,30)
!hmhj endif
c
      return
      end subroutine gfs_dyn_reduce_per_lat
! -------------------------------------------------------------
      subroutine gfs_dyn_reduce_pln2i(lat,qtt)
      implicit none
!
! define x number constant start
      integer,  parameter :: in_f = 960 , in_h = in_f/2
      real(kind=kind_dbl_prec), parameter :: bb_f = 2.d0 ** ( in_f )
      real(kind=kind_dbl_prec), parameter :: bs_f = 2.d0 ** (-in_f )
      real(kind=kind_dbl_prec), parameter :: bb_h = 2.d0 ** ( in_h )
      real(kind=kind_dbl_prec), parameter :: bs_h = 2.d0 ** (-in_h )
! define x number constant end

      real(kind=kind_dbl_prec)  qtt(jcap1*jcap2)
      integer   lat
cc
      real(kind=kind_dbl_prec) x(jcap1),dpln(lnut2),qlnv(lnut2)
      real(kind=kind_dbl_prec) sinlat,cos2,prod,aa,bb,w
      integer twoj1,ll,lplus,lpt,lpv,lp2,lp1,lp0,len,n,i
      integer iprod,id,ix(jcap1),idpln(lnut2)
cc
      twoj1=2*jcap1
c
      sinlat = cos(colrad(lat))
      cos2   = 1.0 - sinlat * sinlat

      prod   = 1.0
      iprod  = 0
      do ll=1,jcap1
        x(ll) = sqrt(0.5*prod)
        ix(ll) = iprod/2
        if( mod(iprod,2).ne.0 ) then
          if( iprod.gt.0 ) then
            x(ll) = x(ll) * bb_h
          else
            x(ll) = x(ll) * bs_h
          endif
        endif
        prod = prod*cos2*yy(ll)
        w = abs(prod)
        if( w.ge.bb_h ) then
          prod = prod * bs_f
          iprod = iprod + 1
        elseif( w.lt.bs_h ) then
          prod = prod * bb_f
          iprod = iprod - 1
        endif
      enddo

      do ll=1,jcap1
        w = abs(x(ll))
        if( w.ge.bb_h ) then
          x(ll) = x(ll) * bs_f
          ix(ll) = ix(ll) + 1
        elseif( w.lt.bs_h ) then
          x(ll) = x(ll) * bb_f
          ix(ll) = ix(ll) - 1
        endif
      enddo

      do ll=1,jcap1
        dpln(2*ll-1) = x(ll)
        dpln(2*ll  ) = x(ll)
        idpln(2*ll-1) = ix(ll)
        idpln(2*ll  ) = ix(ll)
      enddo
c
      lplus = twoj1
      do ll=1,twoj1
        dpln(ll+lplus) = dx(ll) * sinlat * dpln(ll)
        idpln(ll+lplus) = idpln(ll)
        w = abs(dpln(ll+lplus))
        if( w.ge.bb_h ) then
          dpln(ll+lplus) = dpln(ll+lplus) * bs_f
          idpln(ll+lplus) = idpln(ll+lplus) + 1
        elseif( w.lt.bs_h ) then
          dpln(ll+lplus) = dpln(ll+lplus) * bb_f
          idpln(ll+lplus) = idpln(ll+lplus) - 1
        endif
      enddo
      lp2 = 0
      lp1 =     twoj1
      lp0 = 2 * twoj1
      len =     twoj1 - 2
      do n=3,jcap2
        do ll=1,len
!x           dpln(ll+lp0) = (sinlat * dpln(ll+lp1)
!x   1               - deps(ll+lp1) * dpln(ll+lp2)) * rdeps(ll+lp0)
             aa = sinlat * rdeps(ll+lp0)
             bb = deps(ll+lp1) * rdeps(ll+lp0)
             id = idpln(ll+lp1)-idpln(ll+lp2)
             if( id.eq.0 ) then
               dpln(ll+lp0) = aa*dpln(ll+lp1) - bb*dpln(ll+lp2)
               idpln(ll+lp0) = idpln(ll+lp2)
             elseif( id.eq.1 ) then
               dpln(ll+lp0) = aa*dpln(ll+lp1) - bb*dpln(ll+lp2)*bs_f
               idpln(ll+lp0) = idpln(ll+lp1)
             elseif( id.eq.-1 ) then
               dpln(ll+lp0) = aa*dpln(ll+lp1)*bs_f - bb*dpln(ll+lp2)
               idpln(ll+lp0) = idpln(ll+lp2)
             elseif( id.gt.1 ) then
               dpln(ll+lp0) = aa*dpln(ll+lp1)
               idpln(ll+lp0) = idpln(ll+lp1)
             else
               dpln(ll+lp0) = - bb*dpln(ll+lp2)
               idpln(ll+lp0) = idpln(ll+lp2)
             endif
             w = abs(dpln(ll+lp0))
             if( w.ge.bb_h ) then
               dpln(ll+lp0) = dpln(ll+lp0) * bs_f
               idpln(ll+lp0) = idpln(ll+lp0) + 1
             elseif( w.lt.bs_h ) then
               dpln(ll+lp0) = dpln(ll+lp0) * bb_f
               idpln(ll+lp0) = idpln(ll+lp0) - 1
             endif
        enddo
        lp2 = lp1
        lp1 = lp0
        lp0 = lp0 + len
        len = len - 2
      enddo
cc
cc    transpose vector dpln array from cray order to ibm order.
      do i=1,lnut2
!x      qlnv(indxmv(i)) = dpln(i)
        if( idpln(i).eq.0 ) then
          qlnv(indxmv(i)) = dpln(i)
        elseif( idpln(i).eq.-1 ) then
          qlnv(indxmv(i)) = dpln(i) * bs_f
        elseif( idpln(i).lt.-1 ) then
          qlnv(indxmv(i)) = 0.0
        else
          qlnv(indxmv(i)) = dpln(i) * bb_f
        endif
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
      end subroutine gfs_dyn_reduce_pln2i
! -------------------------------------------------------------
      subroutine gfs_dyn_reduce_gpln2i()
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
      end subroutine gfs_dyn_reduce_gpln2i
! -------------------------------------------------------------
      subroutine gfs_dyn_reduce_glats()
      implicit none
c
      real(kind=kind_qdt_prec) eps,si,rl2,scale,dradz,rad,drad,p2,p1
      real(kind=kind_qdt_prec) x,w,sn,rc,pi
      integer l2,k1,k,iter
c
!x    eps=1.e-12
      eps=1.e-20
      si = 1.0
      l2=2*latgh
      rl2=l2
      scale = 2.0/(rl2*rl2)
      k1=l2-1
      pi = atan(si)*4.0 
      dradz = pi / float(latgh) / 200.
      rad = 0.0
      do 1000 k=1,latgh
        iter=0
        drad=dradz
1       call gfs_dyn_reduce_poly(l2,rad,p2)
2       p1 =p2
        iter=iter+1
        rad=rad+drad
        call gfs_dyn_reduce_poly(l2,rad,p2)
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
      end subroutine gfs_dyn_reduce_glats
! ------------------------------------------------------------------
      subroutine gfs_dyn_reduce_poly(n,rad,p)
      implicit none
      integer n
      real(kind=kind_qdt_prec) rad,p
      real(kind=kind_qdt_prec) x,y1,y2,g,y3
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
      end subroutine gfs_dyn_reduce_poly

      end module gfs_dyn_reduce_lons_grid_module
