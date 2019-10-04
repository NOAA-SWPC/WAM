module stenedgprint
contains
! -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
! stenedg is an offshoot of subroutine stencl. besides displaying grid
! point values associated with cell centers, it displays a second field
! of variables defined on cell edges. Specifically, stenedg displays
!
! (a) 'fld' at the center point (ipn=PrintIpnDiag),
! (b) 'fld' at the 5 or 6 icos points surrounding PrintIpnDiag,
! (c) 'fld_edg' on the 5 or 6 edges of cell PrintIpnDiag.
!
! data are printed at their approximate geographic location relative
! to the center point, with south pointing down.
!
!     R.Bleck              October 2009
! -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

   subroutine stenedg(fld,fld_edg,kdm,what)

   use module_control,  only: nip,npp,PrintIpnDiag,			&
                              stencl_frst,stencl_last,stencl_step	
   use module_constants,only: deg_lat, deg_lon, nprox, prox, pi
   implicit none
!  real     ,intent(IN) :: sclfac,scl_edg	! scale factors for F fmt printg
   integer  ,intent(IN) :: kdm
   character,intent(IN) :: what*(*)
!SMS$DISTRIBUTE (dh,nip) BEGIN
   real,     intent(IN) :: fld    (kdm    ,nip)	! values in cell centers
   real,     intent(IN) :: fld_edg(kdm,npp,nip)	! values along cell edges
!SMS$DISTRIBUTE END

   integer,parameter :: maxpts=6
!  integer,parameter :: idm=50,jdm=idm/3		! F format option
   integer,parameter :: idm=60,jdm=idm/3		! E format option
   character map(idm,jdm)*1,string*9
   integer edg,i,j,k,m,n,ntot,ipn,ipx,jmin,jmax
   real*8 xmin,xmax,ymin,ymax,diflon,arg
   real*8 alat(0:maxpts),alon(0:maxpts),				&
          dist(maxpts),angl(maxpts),x(0:maxpts),y(0:maxpts),		&
          valu(kdm,0:maxpts),valu_edg(kdm,maxpts)
   integer cell(0:maxpts),iloc(0:maxpts),jloc(0:maxpts),cell1(0:maxpts)
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
!  print '(3a,i8)','entering stenedg to print ',what,' at ipn =',	&
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
     alat(0)=deg_lat(ipn)
     alon(0)=deg_lon(ipn)
     valu(kfrst:klast,0)=fld(kfrst:klast,ipn)
!    print '(a,i3,2f9.2,i9)','lat/lon/ipn of cell',0,alat(0),alon(0),ipn
    end if
!SMS$IGNORE END
   end do
!SMS$REDUCE (alat,alon,valu,MAX)
!SMS$PARALLEL END

! --- find lat/lon of first (and only) ring of cells

!SMS$PARALLEL (dh,ipn) BEGIN
!SMS$EXCHANGE(fld)
   ntot=-99
   alat    (1:maxpts)=-999.
   alon    (1:maxpts)=-999.
   valu    (:,1:maxpts)=-1.e33
   valu_edg(:,1:maxpts)=-1.e33
   do ipn=1,nip
!SMS$IGNORE BEGIN
    if (ipn.eq.PrintIpnDiag) then
     ntot=0
     do edg=1,nprox(ipn)
      ipx=prox(edg,ipn)
      ntot=ntot+1
      alat(ntot)=deg_lat(ipx)
      alon(ntot)=deg_lon(ipx)
      valu(kfrst:klast,ntot)=fld(kfrst:klast,ipx)
!     print '(a,i3,2f9.2)','lat/lon of cell',ntot,alat(ntot),alon(ntot)
      valu_edg(kfrst:klast,edg)=fld_edg(kfrst:klast,edg,ipn)
     end do
!    print '(a,i5,es11.2,2f8.2/(i25,es11.2,2f8.2))',			&
!     'edge values,lat/lon:',(edg,valu_edg(kfrst,edg),alat(edg),	&
!      alon(edg),edg=1,nprox(ipn))
    end if
!SMS$IGNORE END
   end do
!SMS$REDUCE (ntot,alat,alon,valu,valu_edg,MAX)
!SMS$PARALLEL END

! --- knowing lat/lon of cells, find their indices

   do n=0,ntot
!   print '(a,i3,a,2f9.2)','n =',n,'   looking for cell at lat/lon',	&
!     alat(n),alon(n)

!SMS$PARALLEL (dh,ipn) BEGIN
    cell1(:)=-99
    do ipn=1,nip
!SMS$IGNORE BEGIN
     if (deg_lat(ipn).eq.alat(n) .and. deg_lon(ipn).eq.alon(n)) then
!     print *,'cell',ipn,' has desired lat/lon'
      cell1(n)=ipn
     end if
!SMS$IGNORE END
    end do
!SMS$REDUCE (cell1,MAX)
!SMS$PARALLEL END
   cell(n)=cell1(n)
   end do

!  print '(a,i9,a/(i9,2f9.2,f9.1))','adding',ntot,' cells', 		&
!    (cell(n),alat(n),alon(n),valu(kfrst,n),n=1,ntot)

! --- all cluster points have now been identified

   xmin= 1.e33
   xmax=-1.e33
   ymin= 1.e33
   ymax=-1.e33

   map(:,:)=' '
! --- dist = distance (deg) between cell 'PrintIpnDiag' and cell 'n'.
! --- angl = compass heading from cell 'PrintIpnDiag' to cell 'n'.

   x(0)=0.
   y(0)=0.
   alat(0) = max(-89.999_8, min(alat(0), 89.999_8)) ! avoid centering on poles
   alat(:) = 90. - alat(:)                          ! convert to colatitude
!  print 100
100 format ('    cell    lat     lon    dist    angl      x       y')
!  print 101,cell(0),alat(0),alon(0),0.,0.,x(0),y(0)
101 format (i8,2f8.2,f8.3,f8.2,2f8.3)
   do n=1,ntot
    diflon=alon(n)-alon(0)
    dist(n)=acosd(cosd(alat(0))*cosd(alat(n)) &
                 +sind(alat(0))*sind(alat(n))*cosd(diflon))
    arg = max(-1._8, min(1._8, sind(alat(n))*sind(diflon)/sind(dist(n))))
    angl(n)=asind(arg)
! --- determine quadrant by plugging angl=90 into law of cosines
    if (cosd(alat(0))*cosd(dist(n)).gt.cosd(alat(n))) then
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

   jmin= 999
   jmax=-999
   do n=0,ntot
!   i=1.5+(idm-8.)*(y(n)-ymin)/(ymax-ymin)		! F format option
    i=1.5+(idm-10.)*(y(n)-ymin)/(ymax-ymin)		! E format option
    j=1.5+(jdm-1.)*(x(n)-xmin)/(xmax-xmin)
    iloc(n)=i
    jloc(n)=j

! --- to save space, crowd entries around mid point
!   i=1.5+(idm-8.)*(.54*y(n)-ymin)/(ymax-ymin)		! F format option
    i=1.5+(idm-10.)*(.54*y(n)-ymin)/(ymax-ymin)		! E format option
    j=1.5+(jdm-1.)*(.54*x(n)-xmin)/(xmax-xmin)
    jmin=min(jmin,j)
    jmax=max(jmax,j)

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
   print *,'shown below: grid indices (compressed stencil)'
   do j=jmax,jmin,-1
    print '(2x,75a1)',(map(i,j),i=1,idm)
   end do
   print *,             &
   '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'

! --- print cluster of icos cell values

   do k=kfrst,klast,kstep
    do n=0,ntot
     i=iloc(n)
     j=jloc(n)
!    write (string,'(f7.1)') valu(k,n)*sclfac		! F format option
     write (string,'(es9.2)') valu(k,n)			! E format option
!    do m=1,7						! F format option
     do m=1,9						! E format option
      map(i+m,j)=string(m:m)
     end do
    end do

! --- print values on icos cell edges

    do n=1,ntot
!    i=1.5+(idm-8.)*(.54*y(n)-ymin)/(ymax-ymin)		! F format option
     i=1.5+(idm-10.)*(.54*y(n)-ymin)/(ymax-ymin)	! E format option
     j=1.5+(jdm-1.)*(.54*x(n)-xmin)/(xmax-xmin)
!    write (string,'(f7.1)') valu_edg(k,n)*scl_edg	! F format option
     write (string,'(es9.2)') valu_edg(k,n)		! E format option
!    do m=1,7						! F format option
     do m=1,9						! E format option
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

!  print *,'exiting stenedg, PrintIpnDiag =',PrintIpnDiag
!!SMS$BARRIER
   return
   end subroutine stenedg
end module stenedgprint
