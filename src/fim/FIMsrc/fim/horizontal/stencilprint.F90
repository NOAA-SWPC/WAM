module stencilprint
contains
! -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
! This routine displays a cluster of icos values centered on grid cell
! 'PrintIpnDiag' (specified in FIMnamelist). The user-selectable
! parameter 'num_rings' specifies the number of concentric rings of
! data points displayed around the center point. Values are printed in
! their approximate geographic location relative to the center point,
! scaled to an 80-col screen with south pointing down.
!
! If 'field' has a vertical dimension (i.e. kdm > 1), data are printed
! for levels kfrst ... klast in steps of kstep.
!
! If 'field' has no vertical dimension, call stencl with kdm=1. 
!
!     R.Bleck              October 2009
! -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

   subroutine stencl(field,kdm,sclfac,what)

   use module_control,  only: nip, PrintIpnDiag,			&
                              stencl_frst,stencl_last,stencl_step	
   use module_constants,only: deg_lat, deg_lon, nprox, prox, pi

   implicit none
   real     ,intent(IN) :: sclfac	! scale factor for F format printing
   integer  ,intent(IN) :: kdm
   character,intent(IN) :: what*(*)
!SMS$DISTRIBUTE (dh,nip) BEGIN
   real,     intent(IN) :: field(kdm,nip)
!SMS$DISTRIBUTE END

   integer,parameter :: num_rings=2	! number of concentric rings
!  integer,parameter :: num_rings=3	! number of concentric rings

   integer,parameter :: maxpts=1+6*(2**num_rings-1)
   integer,parameter :: idm=75*num_rings/max(num_rings,3),jdm=idm/3
   character map(idm,jdm)*1,string*7
   integer edg,i,j,k,m,n,ntot,ipn,ipx,nxtring,nr,nold,nnew,cumutot
   logical newcell
   real*8 xmin,xmax,ymin,ymax,diflon,arg
   real*8 x(maxpts),y(maxpts),dist(maxpts),angl(maxpts),		&
          valu(kdm,maxpts),alat(maxpts),alon(maxpts),			&
          cumulat(maxpts),cumulon(maxpts),cumuval(kdm,maxpts)
   integer cell(maxpts),iloc(maxpts),jloc(maxpts),cumucell(maxpts),	&
           cell1(maxpts)
   integer :: kfrst,kstep,klast

#if ( defined NEED_SINDCOSD )
!JR Define required statement functions for situations where
!JR compiler doesn't support them (e.g. gfortran)
   real*8 :: val, sind, cosd, asind, acosd
   sind(val) = sin(val*pi/180.)
   cosd(val) = cos(val*pi/180.)

   asind(val) = (180./pi)*asin(val)
   acosd(val) = (180./pi)*acos(val)
#endif

   if (PrintIpnDiag.le.0) return
!  print '(3a,i8)','entering stencl to print ',what,' at ipn =',	&
!    PrintIpnDiag
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
! --- layers printed are: kfrst, kfrst+kstep, kfrst+2*kstep, ..., klast
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
   if (kdm.eq.1) then
     kfrst=1
     klast=1
     kstep=1
   else
     if (stencl_step.eq.0) return
     kfrst=stencl_frst
     klast=stencl_last
     kstep=stencl_step
   end if

! --- build up inventory of cell indices, starting with center cell

!SMS$PARALLEL (dh,ipn) BEGIN
   alat(:)=-999.
   alon(:)=-999.
   valu(:,:)=-1.e33
   do ipn=1,nip
!SMS$IGNORE BEGIN
    if (ipn.eq.PrintIpnDiag) then
     alat(1)=deg_lat(ipn)
     alon(1)=deg_lon(ipn)
     valu(kfrst:klast,1)=field(kfrst:klast,ipn)
!    print '(a,i3,2f9.2)','lat/lon/ipn of cell',0,alat(1),alon(1)
    end if
!SMS$IGNORE END
   end do
!SMS$REDUCE (ntot,alat,alon,valu,MAX)
!SMS$PARALLEL END

   cumucell(1)=PrintIpnDiag
   cumulat(1)=alat(1)
   cumulon(1)=alon(1)
   cumuval(kfrst:klast,1)=valu(kfrst:klast,1)
   cumutot=1
   nold=0

! --- expand cluster. each new ring is union of neighbors of previous ring

   do nxtring=1,num_rings

! --- neighbors of cells 1...nold have already been located.
! --- new ring consists of neighbors of cells (nold+1)...cumutot

    nnew=cumutot
    do nr=nold+1,nnew
!   print '(3(a,i3))','starting ring',nxtring,				&
!    '. looking for neighbors of cells',nold+1,' -',nnew

! --- find lat/lon of cells in ring 'nxtring'

!SMS$PARALLEL (dh,ipn) BEGIN
!SMS$EXCHANGE(field)
     ntot=-99
     alat(:)=-999.
     alon(:)=-999.
     valu(:,:)=-1.e33
     do ipn=1,nip
!SMS$IGNORE BEGIN
      if (ipn.eq.cumucell(nr)) then
       ntot=0
       do edg=1,nprox(ipn)
        ipx=prox(edg,ipn)
        ntot=ntot+1
        alat(ntot)=deg_lat(ipx)
        alon(ntot)=deg_lon(ipx)
        valu(kfrst:klast,ntot)=field(kfrst:klast,ipx)
!       print '(a,i9,2f9.2)','lat/lon/ipn of neighbor',ipx,	&
!         alat(ntot),alon(ntot)
       end do
      end if
!SMS$IGNORE END
     end do
!SMS$REDUCE (ntot,alat,alon,valu,MAX)
!SMS$PARALLEL END

! --- knowing lat/lon of cells in ring 'nxtring', find their indices

     do n=1,ntot
!     print '(a,i3,a,2f9.2)','n =',n,'   looking for cell at lat/lon',	&
!       alat(n),alon(n)

!SMS$PARALLEL (dh,ipn) BEGIN
      cell1(:)=-99
      do ipn=1,nip
!SMS$IGNORE BEGIN
       if (deg_lat(ipn).eq.alat(n) .and. deg_lon(ipn).eq.alon(n)) then
!       print '(a,i9,a,2f9.2)','cell',ipn,' has desired lat/lon',	&
!         alat(n),alon(n)
        cell1(n)=ipn
       end if
!SMS$IGNORE END
      end do
!SMS$REDUCE (cell1,MAX)
!SMS$PARALLEL END
      cell(n)=cell1(n)
     end do

!    print '(a,i4,a/(i9,2f9.2,f9.1))','considering',ntot,' cells', 	&
!      (cell(n),alat(n),alon(n),valu(kfrst,n),n=1,ntot)

! --- eliminate duplicates

     do n=1,ntot
      newcell=.true.
      do m=1,cumutot
!      print *,'comparing',n,cell(n),m,cumucell(m)
       if (cell(n).eq.cumucell(m)) then
        newcell=.false.
!       print *,'rejecting',cell(n)
        exit
       end if
      end do
      if (newcell) then
!      print *,'accepting',cell(n)
       cumutot=cumutot+1
       if (cumutot.gt.maxpts) stop '(number of cells > maxpts)'
       cumucell(cumutot)=cell(n)
       cumulat(cumutot)=alat(n)
       cumulon(cumutot)=alon(n)
       cumuval(kfrst:klast,cumutot)=valu(kfrst:klast,n)
      end if
     end do
!    print '(3(a,i4)/(i9,2f9.2,f9.1))','total cell count is',cumutot,	&
!     ' after locating neighbors of ',nr,' cells in ring',nxtring-1,	&
!      (cumucell(n),cumulat(n),cumulon(n),cumuval(kfrst,n),n=1,cumutot)
    end do			! nr loop
!   print *,'ring',nxtring,'  completed'
    nold=nnew
   end do			! nxtring

   ntot=cumutot
   do n=1,ntot
    cell(n)=cumucell(n)
    alat(n)=cumulat(n)
    alon(n)=cumulon(n)
    valu(kfrst:klast,n)=cumuval(kfrst:klast,n)
   end do

!  print '(a/(i4,i9,2f9.2,f9.1))',					&
!   'all cluster points have now been identified:',			&
!    (n,cell(n),alat(n),alon(n),valu(kfrst,n),n=1,ntot)

   xmin= 1.e33
   xmax=-1.e33
   ymin= 1.e33
   ymax=-1.e33

   map(:,:)=' '

! --- dist = distance (deg) between cell 'PrintIpnDiag' and cell 'n'.
! --- angl = compass heading from cell 'PrintIpnDiag' to cell 'n'.

   x(1)=0.
   y(1)=0.
   alat(1) = max(-89.999_8, min(alat(1),89.999_8))  ! avoid centering on poles
   alat(:) = 90. - alat(:)                          ! convert to colatitude
!  print 100
100 format ('    cell    lat     lon    dist    angl      x       y')
!  print 101,cell(1),alat(1),alon(1),0.,0.,x(1),y(1)
101 format (i8,2f8.2,f8.3,f8.2,2f8.3)
   do n=2,ntot
    diflon=alon(n)-alon(1)
    dist(n)=acosd(cosd(alat(1))*cosd(alat(n)) &
                 +sind(alat(1))*sind(alat(n))*cosd(diflon))
    arg = max(-1._8, min(1._8, sind(alat(n))*sind(diflon)/sind(dist(n))))
    angl(n)=asind(arg)
! --- determine quadrant by plugging angl=90 into law of cosines
    if (cosd(alat(1))*cosd(dist(n)).gt.cosd(alat(n))) then
      angl(n)= 180.-angl(n)
    end if
    x(n)=dist(n)*cosd(angl(n))
    y(n)=dist(n)*sind(angl(n))
!   print 101,cell(n),alat(n),alon(n),dist(n),angl(n),x(n),y(n)
    xmin=min(xmin,x(n))
    ymin=min(ymin,y(n))
    xmax=max(xmax,x(n))
    ymax=max(ymax,y(n))
   end do

!  print *,'lat/lon limits:',xmin,xmax,ymin,ymax

! --- print cluster of icos cell identifiers

   do n=1,ntot
    i=1.5+(idm-8.)*(y(n)-ymin)/(ymax-ymin)
    j=1.5+(jdm-1.)*(x(n)-xmin)/(xmax-xmin)
    iloc(n)=i
    jloc(n)=j
    write (string,'(i7)') cell(n)
    if (cell(n).lt.100000) write (string,'(i6,1x)') cell(n)
    if (cell(n).lt.1000  ) write (string,'(i5,2x)') cell(n)
    if (cell(n).lt.10    ) write (string,'(i4,3x)') cell(n)
    do m=1,7
     map(i+m,j)=string(m:m)
    end do
   end do
   print *,             &
   '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
   print *,'shown below: grid indices'
   do j=jdm,1,-1
    print '(2x,75a1)',(map(i,j),i=1,idm)
   end do
   print *,             &
   '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'

! --- print cluster of icos cell values

   do k=kfrst,klast,kstep
    do n=1,ntot
     i=iloc(n)
     j=jloc(n)
     write (string,'(f7.1)') valu(k,n)*sclfac
     do m=1,7
      map(i+m,j)=string(m:m)
     end do
    end do
    if (kdm.eq.1) then
     print '(2a)','shown below: ',trim(what)
    else
     print '(3a,i4)','shown below: ',trim(what),', k =',k
    end if
    do j=jdm,1,-1
     print '(2x,75a1)',(map(i,j),i=1,idm)
    end do
   print *,             &
   '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
   end do			! vertical loop

!  print *,'exiting stencl, PrintIpnDiag =',PrintIpnDiag
!!SMS$BARRIER
   return
   end subroutine stencl
end module stencilprint
