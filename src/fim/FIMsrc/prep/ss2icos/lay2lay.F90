      subroutine lay2lay(kkold,pold,thold,				&
                         varin1,varin2,varin3,varin4,varin5,		&
                         kknew,pnew,thnew,				&
                         varou1,varou2,varou3,varou4,varou5,		&
                         targt,nip,PrintDiags)
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!   Layer-to-layer conversion routine, obtained by combining routines
!   lay2lev4 (layer to level) and lin2stp (linear to step fct)
!
!   (Lay2lev based on Bleck, R., 1984: Vertical Coordinate Transformation
!   of Vertically Discretized Atmospheric Fields. Mon. Wea. Rev., 112,
!   2537-2541)
!
!   Specifically, lay2lay ....

!     (1) fits continuous, piecewise linear curves to input stairstep
!         profiles thold(pold), varin1(pold), varin2(pold), ...;
!     (2) creates a new stairstep profile of theta by piecewise
!         integrating over the continous theta curve such that the new
!         steps theta(pnew) match prescribed 'target' values;
!     (3) modifies the new pressure -pnew- where necessary to satisfy
!         minimum layer thickness constraints;
!     (4) integrates the continuous curves for thold, varin1, varin2,...
!         over pressure intervals obtained in (3) to produce new
!         stairstep profiles varout1(pnew), varout2(pnew), ...;
!
!   The routine is presently configured to conserve the height of the
!   air column. It can easily be modified to conserve column thermal
!   energy instead.
!
!   Rainer Bleck	2008
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
! --- More details:
!
! --- 'Old' data consisting of 3-dim. arrays of
!
! ---     (a) 'kkold+1' interface pressure (or Exner fct) values (pold)
! ---     (b) 'kkold' layer averages of pot.temperature (thold)
! ---     (b) 'kkold' layer averages of dep.variables (varin1,varin2,..)
!
! --- ...are transformed from piecewise constant (stairstep-type)
! --- functions of -p- to continuous, piecewise linear functions of -p-
! --- preserving the column mean of the input data. This intrinsically
! --- nonunique curve fitting problem is rendered unique by minimizing
! --- zigzags in the curves, i.e., imposing penalties for large 2nd
! --- derivatives. Additional penalties are imposed for large overshoots
! --- in the center portion of each layer. (Sharply raising the latter
! --- penalties will force the algorithm to generate curves resembling
! --- the input stairstep curves.)
!
! --- The resulting line segments are integrated to form layers
! --- constrained to yield prescribed values of pot.temperature (targt).
! --- New values of intfc pressure, theta, and dep.variables are returned
! --- in pnew(:,:,kknew+1), thnew(:,:,kknew), varou1(:,:,kknew), ....
!
      use module_hybgen
      use module_constants,only: p1000,rd,cp
      use module_control, only: yyyymmddhhmm,PrintIpnDiag,ptop
!
      implicit none
      integer,intent(IN) :: nip,kkold,kknew
      logical,intent(IN) ::  PrintDiags
!SMS$DISTRIBUTE(dh,NIP) BEGIN
      real,intent(IN)  :: pold(kkold+1,nip),thold (kkold,nip),		&
                            varin1(kkold,nip),varin2(kkold,nip),	&
                            varin3(kkold,nip),varin4(kkold,nip),	&
                            varin5(kkold,nip),				&
                            targt (kknew,nip)
      real,intent(OUT) :: pnew(kknew+1,nip),thnew (kknew,nip),		&
                          varou1(kknew,nip),varou2(kknew,nip),		&
                          varou3(kknew,nip),varou4(kknew,nip),		&
                          varou5(kknew,nip)
!SMS$DISTRIBUTE END
!
      integer,parameter :: nsize = 321		!  must be 5*kkold+1 or larger
      integer ipn,k,l,n,indx(nsize),last
      real matrix(nsize,nsize),pint(nsize),solth(nsize),		&
           solu1(nsize),solu2(nsize),solu3(nsize),solu4(nsize),		&
           solu5(nsize),excol(kknew+1),thcol(kknew),			&
           v1col(kknew),v2col(kknew),v3col(kknew),v4col(kknew),		&
           v5col(kknew),tarcol(kknew),prcol(kknew+1),			&
           exold(kknew+1),exwrk(kknew+1),unusd(nsize),			&
           oddev,p2ex,ex2p,arg
      real pr_extd(nsize+3),th_extd(nsize+3),tg_extd(kknew+1)
!
      real,parameter :: penlty = 0.		! penalty for midlyr overshoots
      real,parameter :: flag = -.03125		! missing data
      logical realyr,vrbos
!
      ex2p(arg)=p1000*(arg/cp)**(cp/rd)		!  convert Pi => p
!     p2ex(arg)=cp*(arg/p1000)**(rd/cp)		!  convert p  => Pi

      print *,'entering lay2lay...'
      if (nsize.lt.5*kkold+1) stop '(nsize too small in subr.lay2lev)'
!
      n=kkold				!  number of input layers
!
      do 30 k=1,nsize
      solth(k)=0.
      solu1(k)=0.
      solu2(k)=0.
      solu3(k)=0.
      solu4(k)=0.
      solu5(k)=0.
      do 30 l=1,nsize
  30  matrix(k,l)=0.
!
      do 18 k=1,n
! --- upper left quadrant:
      matrix(k,4*k-3)=1.
      matrix(k,4*k-2)=2.
      matrix(k,4*k-1)=2.
      matrix(k,4*k  )=2.
      matrix(k,4*k+1)=1.
! --- lower right quadrant:
      matrix(n+4*k-3,4*n+1+k)=1.
      matrix(n+4*k-2,4*n+1+k)=2.
      matrix(n+4*k-1,4*n+1+k)=2.
      matrix(n+4*k  ,4*n+1+k)=2.
  18  matrix(n+4*k+1,4*n+1+k)=1.
!
! --- lower left quadrant:
      do 19 k=2,4*n
      matrix(n+k+1,k-1)=1.
      matrix(n+k-1,k+1)=1.
  19  matrix(n+k  ,k  )=6.
!
      do 20 k=2,4*n-1
      matrix(n+k+1,k  )=-4.
  20  matrix(n+k  ,k+1)=-4.
!
      matrix(  n+1,    1)=1.
      matrix(5*n+1,4*n+1)=1.
      matrix(5*n  ,4*n+1)=-2.
      matrix(5*n+1,4*n  )=-2.
      matrix(  n+1,    2)=-2.
      matrix(  n+2,    1)=-2.
      matrix(5*n  ,4*n  )=5.
      matrix(  n+2,    2)=5.
!
! --- penalize overshoots at layer midpoints
!
      do 15 k=1,n
      matrix(n+4*k-2,4*k-2)=matrix(n+4*k-2,4*k-2)+penlty*.25
      matrix(n+4*k-1,4*k-1)=matrix(n+4*k-1,4*k-1)+penlty
      matrix(n+4*k  ,4*k  )=matrix(n+4*k  ,4*k  )+penlty*.25
  15  continue
!
! --- decompose matrix
      call ludcmp(matrix,5*n+1,nsize,indx,oddev)
!
!! Cray modification (thanks to Pete Johnsen):
!!      Make multi-line OMP directive conform to the
!!      OpenMP Application Program Interface Version 2.5 May 2005,
!!      Sect. 2.1.2 Free Source Form Directives, according to which:
!!         Continued directive lines must have an ampersand as the
!!         last nonblank character on the line, prior to any comment
!!         placed inside the directive. Continuation directive lines
!!         can have an ampersand after the directive sentinel with
!!         optional white space before and after the ampersand.
!TBH:  NOTE that this file has not been tested with OpenMP in quite a while 
!TBH:  so the directive below mauy require adjustment.  
!!$OMP PARALLEL DO PRIVATE(solth,solu1,solu2,solu3,solu4,solu5,pint,     &
!!$OMP&            kold,knew,dp1,dp2,realyr,vrbos , ipnGlob,mype,last,   &
!!$OMP&            tarcol,pr_extd,th_extd,tg_extd,excol,exold,           &
!!$OMP&            thcol,prcol,dpcol,v1col,v2col,v3col,v4col,v5col)      &
!!$OMP&            SHARED(matrix,indx)
!SMS$PARALLEL (dh,ipn) BEGIN
      do 1 ipn=1,nip
      vrbos=ipn.eq.PrintIpnDiag
!
! --- open file for saving step-by-step details of the interpolation procedure
      if (vrbos) then
!SMS$ignore begin
       print  '(3a,i8,a)','store data for date = ',yyyymmddhhmm,	&
        ', ipn =',ipn,' in file "stairstep_details" for offln plotting'
       open (31,file='stairstep_details',form='formatted')
       write (31,'(a10,i8,4i5)') yyyymmddhhmm,ipn,kkold,4*n+1,kknew,kknew
!SMS$ignore end
      end if
!
      do 10 k=1,kknew
      varou1(k,ipn)=flag
      varou2(k,ipn)=flag
      varou3(k,ipn)=flag
      varou4(k,ipn)=flag
      varou5(k,ipn)=flag
  10  continue
!
      do 2 k=1,n
      solth(k)=8.*thold (k,ipn)
      solu1(k)=8.*varin1(k,ipn)
      solu2(k)=8.*varin2(k,ipn)
      solu3(k)=8.*varin3(k,ipn)
      solu4(k)=8.*varin4(k,ipn)
      solu5(k)=8.*varin5(k,ipn)
   2  continue
!
      do 31 k=n+1,nsize
      solth(k)=0.
      solu1(k)=0.
      solu2(k)=0.
      solu3(k)=0.
      solu4(k)=0.
      solu5(k)=0.
  31  continue
!
! --- penalize overshoots at layer midpoints
!
      do 16 k=1,n
      solth(n+4*k-2)=thold (k,ipn)*penlty*.25
      solth(n+4*k-1)=thold (k,ipn)*penlty
      solth(n+4*k  )=thold (k,ipn)*penlty*.25
!
      solu1(n+4*k-2)=varin1(k,ipn)*penlty*.25
      solu1(n+4*k-1)=varin1(k,ipn)*penlty
      solu1(n+4*k  )=varin1(k,ipn)*penlty*.25
!
      solu2(n+4*k-2)=varin2(k,ipn)*penlty*.25
      solu2(n+4*k-1)=varin2(k,ipn)*penlty
      solu2(n+4*k  )=varin2(k,ipn)*penlty*.25
!
      solu3(n+4*k-2)=varin3(k,ipn)*penlty*.25
      solu3(n+4*k-1)=varin3(k,ipn)*penlty
      solu3(n+4*k  )=varin3(k,ipn)*penlty*.25
!
      solu4(n+4*k-2)=varin4(k,ipn)*penlty*.25
      solu4(n+4*k-1)=varin4(k,ipn)*penlty
      solu4(n+4*k  )=varin4(k,ipn)*penlty*.25
!
      solu5(n+4*k-2)=varin5(k,ipn)*penlty*.25
      solu5(n+4*k-1)=varin5(k,ipn)*penlty
      solu5(n+4*k  )=varin5(k,ipn)*penlty*.25
!
  16  continue
!
! --- replace values in massless layers by values from nearest 'real' layer
!
!cc      realyr=.false.
!cc      do 8 k=1,n
!cc      if (pold(ipn,k+1).gt.pold(ipn,k)+.01) realyr=.true.
!cc      if (k.eq.1) go to 8
!cc      if (realyr .and. pold(ipn,k).ge.pold(ipn,k+1)-.01) then
!cc        solth(k)=solth(k-1)
!cc        solu1(k)=solu1(k-1)
!cc        solu2(k)=solu2(k-1)
!cc        solu3(k)=solu3(k-1)
!cc      end if
!cc   8  continue
!
!cc      realyr=.false.
!cc      do 9 k=n,1,-1
!cc      if (pold(ipn,k+1).gt.pold(ipn,k)+.01) realyr=.true.
!cc      if (k.eq.n) go to 9
!cc      if (realyr .and. pold(ipn,k).ge.pold(ipn,k+1)-.01) then
!cc       solth(k)=solth(k+1)
!cc       solu1(k)=solu1(k+1)
!cc       solu2(k)=solu2(k+1)
!cc       solu3(k)=solu3(k+1)
!cc      end if
!cc   9  continue
!
      if (vrbos) then
!SMS$ignore begin
       write (*,101) ipn,'thold  inpt',(pold(k,ipn),.125*solth(k),k=1,n),&
         pold(n+1,ipn)
!!     write (*,101) ipn,'vrbl 1 inpt',(pold(k,ipn),.125*solu1(k),k=1,n)
!!     write (*,101) ipn,'vrbl 2 inpt',(pold(k,ipn),.125*solu2(k),k=1,n)
!!     write (*,101) ipn,'vrbl 3 inpt',(pold(k,ipn),.125*solu3(k),k=1,n)
!!     write (*,101) ipn,'vrbl 4 inpt',(pold(k,ipn),.125*solu4(k),k=1,n)
!!     write (*,101) ipn,'vrbl 5 inpt',(pold(k,ipn),.125*solu5(k),k=1,n)
       write (31,*) 'input profile'
       write (31,100) (pold(k,ipn),.125*solth(k),k=1,n)
 100   format (2es15.7)
!SMS$ignore end
      end if
 101  format (i10,3x,a,' profile:'/(4(f9.1,es10.3)))
!
      call lubksb(matrix,5*n+1,nsize,indx,solth)
      call lubksb(matrix,5*n+1,nsize,indx,solu1)
      call lubksb(matrix,5*n+1,nsize,indx,solu2)
      call lubksb(matrix,5*n+1,nsize,indx,solu3)
      call lubksb(matrix,5*n+1,nsize,indx,solu4)
      call lubksb(matrix,5*n+1,nsize,indx,solu5)
!
! --- assign pressure values to the end points of the 4*n line segments
!
      do 4 k=1,n
   4  pint(4*k-3)=pold(k,ipn      )
      pint(4*n+1)=pold(kkold+1,ipn)
      do 5 k=1,n
      pint(4*k-2)=.75*pint(4*k-3)+.25*pint(4*k+1)
      pint(4*k-1)=.50*pint(4*k-3)+.50*pint(4*k+1)
   5  pint(4*k  )=.25*pint(4*k-3)+.75*pint(4*k+1)
!
      if (vrbos) then
!SMS$ignore begin
       write (*,101) ipn,'theta  pcwise lin.',(pint(k),solth(k),k=1,4*n+1)
!!     write (*,101) ipn,'vrbl 1 pcwise lin.',(pint(k),solu1(k),k=1,4*n+1)
!!     write (*,101) ipn,'vrbl 2 pcwise lin.',(pint(k),solu2(k),k=1,4*n+1)
!!     write (*,101) ipn,'vrbl 3 pcwise lin.',(pint(k),solu3(k),k=1,4*n+1)
!!     write (*,101) ipn,'vrbl 4 pcwise lin.',(pint(k),solu4(k),k=1,4*n+1)
!!     write (*,101) ipn,'vrbl 5 pcwise lin.',(pint(k),solu5(k),k=1,5*n+1)
       write (31,*) 'piecewise linear fit'
       write (31,100) (pint(k),solth(k),k=1,4*n+1)
!SMS$ignore end
      end if
!
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
! --- Step 1 completed. Now break -solth- into stairsteps whose 'risers'
! --- are at prescribed -targt- values
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!
      tarcol(:)=targt(:,ipn)
!
! --- extend fitted curve at bottom and top to match range of target values
! --- extend list of target values to match range of fitted curve
!
      pr_extd(2:4*n+2)=pint(1:4*n+1)
      pr_extd(    1)=pint(    1)
      pr_extd(4*n+3)=pint(4*n+1)
!
      th_extd(2:4*n+2)=solth(1:4*n+1)
      th_extd(    1)=min(tarcol(    1),solth(    1))
      th_extd(4*n+3)=max(tarcol(kknew),solth(4*n+1))
!
      tg_extd(1:kknew)=tarcol
      tg_extd(      1)=min(tarcol(    1),solth(    1))
      tg_extd(kknew+1)=max(tarcol(kknew),solth(4*n+1))
!
      do k=2,4*n+3
        th_extd(k)=max(th_extd(k),th_extd(k-1))		! remove superadiabats
      end do
!
      call lin2stp(th_extd,pr_extd,4*n+3,tg_extd,excol,kknew,		&
                   vrbos,ipn)
      do k=kknew,2,-1
       excol(k+1)=excol(k)
      end do
      excol(1)=pint(1)
      thcol(:)=tarcol(:)
!
      if (vrbos) then
!SMS$ignore begin
       write (31,*) 'output profile before hybridization'
       write (31,100) (excol(k),thcol(k),k=1,kknew)
!SMS$ignore end
      end if
!
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
! --- Step 2 completed. Now hybridize the grid
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!
! --- FIM atmosphere ends at pressure -ptop-
      excol(kknew+1)=cp*(ptop/p1000)**(rd/cp)
      do 14 k=kknew,kknew/2,-1
   14 excol(k)=max(excol(k),excol(k+1))
!
! --- inflate massless layers to create hybrid-isentropic grid.
! --- it is assumed here that -pint- represents Exner fctn, *not* pressure
!
      prcol(  1)=ex2p(excol(  1))
      do 13 k=1,kknew
      excol(k+1)=min(excol(k),excol(k+1))		! remove superadiabats
      prcol(k+1)=ex2p(excol(k+1))
  13  continue
      exold(:)=excol(:)
!
! --- for consistency between initialization and time integration,
! --- regrid_1d (part of the FIM grid generator) is used to here
! --- to inflate massless layers. however, vertical regridding of
! --- dependent variables will be done separately by lin2stp.
!
!SMS$ignore begin
      if (vrbos) print *,'lay2lay calling regrid_1d/remap_1d ...'
!SMS$ignore end
!
      call regrid_1d(0,tarcol,thcol,excol,prcol,vrbos,ipn,PrintDiags)
      exwrk(:)=exold(:)
      call remap_1d(0,tarcol,1,thcol,unusd,unusd,exwrk,excol,		&
                unusd,prcol,vrbos,ipn,PrintDiags)
      excol(:)=exwrk(:)
!
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
! --- Step 3 completed. Now interpolate input variables to the new grid.
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!
      do k=1,4*n+1
       pint(k)=max(pint(k),excol(kknew+1))
      end do
! 
! --- do full-column regridding for all variables except theta
!
      call lin2stp(pint,solu1,4*n+1,excol,v1col,kknew,vrbos,ipn)
      call lin2stp(pint,solu2,4*n+1,excol,v2col,kknew,vrbos,ipn)
      call lin2stp(pint,solu3,4*n+1,excol,v3col,kknew,vrbos,ipn)
      call lin2stp(pint,solu4,4*n+1,excol,v4col,kknew,vrbos,ipn)
      call lin2stp(pint,solu5,4*n+1,excol,v5col,kknew,vrbos,ipn)
!
! --- regrid theta in inflated layers only. this keeps theta values in
! --- isentropic region on target but may violate the integral constraint
!
      do k=2,kkold+1
       if (excol(k).eq.exold(k)) go to 21
      end do
      stop '(lay2lay error)'
  21  last=k-1
      do k=1,4*n+1
       pint(k)=max(pint(k),excol(last+1))
      end do
      call lin2stp(pint,solth,4*n+1,excol,thcol,last ,vrbos,ipn)
!
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
! --- Steps 4 completed. Finish up.
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!
      do 6 k=1,kknew
      pnew  (k,ipn)=excol(k)
      thnew (k,ipn)=thcol(k)
      varou1(k,ipn)=v1col(k)
      varou2(k,ipn)=v2col(k)
      varou3(k,ipn)=v3col(k)
      varou4(k,ipn)=v4col(k)
      varou5(k,ipn)=v5col(k)
   6  continue
      pnew(kknew+1,ipn)=excol(kknew+1)
!
      if (vrbos) then
!SMS$ignore begin
       write (*,101) ipn,'theta  outp',(excol(k),thcol(k),k=1,kknew) &
        ,excol(kknew+1)
!!     write (*,101) ipn,'vrbl 1 outp',(excol(k),v1col(k),k=1,kknew) &
!!      ,excol(kknew+1)
!!     write (*,101) ipn,'vrbl 2 outp',(excol(k),v2col(k),k=1,kknew) &
!!      ,excol(kknew+1)
!!     write (*,101) ipn,'vrbl 3 outp',(excol(k),v3col(k),k=1,kknew) &
!!      ,excol(kknew+1)
!!     write (*,101) ipn,'vrbl 4 outp',(excol(k),v4col(k),k=1,kknew) &
!!      ,excol(kknew+1)
!!     write (*,101) ipn,'vrbl 5 outp',(excol(k),v5col(k),k=1,kknew) &
!!      ,excol(kknew+1)
       write (31,*) 'output profile after hybridization'
       write (31,100) (excol(k),thcol(k),k=1,kknew)
       close (31)
!SMS$ignore end
      end if
!
   1  continue
!SMS$PARALLEL END
!
      print *,'... exiting lay2lay'
      return
      end subroutine lay2lay
!
!
      SUBROUTINE lubksb(a,n,np,indx,b)
      INTEGER n,np,indx(n)
      REAL a(np,np),b(n)
      INTEGER i,ii,j,ll
      REAL sum
      ii=0
!
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        do 13 j=i+1,n
          sum=sum-a(i,j)*b(j)
13      continue
        b(i)=sum/a(i,i)
14    continue
      return
      END
!  (C) Copr. 1986-92 Numerical Recipes Software 'W3.
!
!
      SUBROUTINE ludcmp(a,n,np,indx,d)
      INTEGER n,np,indx(n),NMAX
      REAL d,a(np,np),TINY
      PARAMETER (NMAX=500,TINY=1.0e-20)
      INTEGER i,imax,j,k
      REAL aamax,dum,sum,vv(NMAX)
      d=1.
      do 12 i=1,n
        aamax=0.
        do 11 j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11      continue
        if (aamax.eq.0.) pause 'singular matrix in ludcmp'
        vv(i)=1./aamax
12    continue
      do 19 j=1,n
        do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
14      continue
        aamax=0.
        do 16 i=j,n
          sum=a(i,j)
          do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
15        continue
          a(i,j)=sum
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0.)a(j,j)=TINY
        if(j.ne.n)then
          dum=1./a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
      return
      END
!  (C) Copr. 1986-92 Numerical Recipes Software 'W3.
