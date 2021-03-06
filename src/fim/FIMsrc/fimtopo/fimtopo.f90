
module const

real, parameter :: pi=3.1415926535897931
real, parameter :: pi4=pi/4.0
real, parameter :: pi2=pi/2.0
real, parameter :: deg2rad=pi/180.0
real, parameter :: rad2deg=1.0/deg2rad

real, parameter :: rearth=6371.0
real, parameter :: earthcircum=2.0*pi*rearth
real, parameter :: earthomega=7.292e-5

real, parameter :: km2nm=60.0/(2*pi*rearth/360.0)
real, parameter :: nm2km=1.0/km2nm
real, parameter :: deglat2km=((2.0*pi*rearth)/360.0)
real, parameter :: deglat2nm=60.0
real, parameter :: km2deglat=1.0/deglat2km
real, parameter :: nm2deglat=1.0/deglat2nm
real, parameter :: knots2ms=1000.0/(km2nm*3600.0)
real, parameter :: yms2knots=1.0/knots2ms
real, parameter :: epsilonm5=1.0e-5
real, parameter :: gravity=9.80665

end module const

module libmf

contains

subroutine stat2(a,m,n,amin,amax,amean,avar,asigma)

  real(kind=4) :: a(m,n)

  amean = 0.0
  amin = 9.9e25
  amax = -9.9e25
  avar = 0.0
  asigma = 0.0
  do i=1,m
     do j=1,n
        if(a(i,j).lt.amin) amin=a(i,j)
        if(a(i,j).gt.amax) amax=a(i,j)
        amean=amean+a(i,j)
     end do
  end do
  amean = amean/m*n
  do i=1,m
     do j=1,n
        avar = avar + (a(i,j)-amean)**2
     end do
  end do
  avar = avar/(m*n-1)
  asigma=sqrt(avar)
  return
  
end subroutine stat2

subroutine qprntn(a,qtitle,ibeg,jbeg,m,n,iskip,iunit)
      
!
!**********	12 APR 91 this version outputs to iunit 
!**********	using write on the Cray Y/MP 
!
!***************************************************************
!***************************************************************
!*****                                                     *****
!*****       qprint output routine (!orrected 4/26/86)     *****
!*****                                                     *****
!***************************************************************
!***************************************************************
!
! a= fwa of m x n array
! qtitle - title
! ibeg,jbeg=lower left corner coords to be printed
! up to 43 x 83 points printed
!
      real(kind=4) a(m,n),ix(81)
      real(kind=4) xm
      character qtitle*24
!
!  determine grid limits
!
      if(iskip.eq.0) iskip=1
      iend=min0(ibeg+79*iskip,m)
      jend=min0(jbeg+79*iskip,n)

      half=0.5
!
   24 continue
!
!  index backwards checking for max
!
   11 xm=0.
      jendsc=min0(jend,n)
      do j=jbeg,jendsc,iskip
      jend_qp = j
      do i=ibeg,iend,iskip
        xm=max(xm*1.d0,abs(a(i,j)))
      end do
      end do
!
!  determine scaling factor limits
!
      if(xm.lt.1.0e-32.or.xm.eq.0.0) xm=99.0
      xm=alog10(99.0/xm)
      kp=xm
      if(xm.lt.0.0)kp=kp-1
!
!  print scaling constants
!
   12 write(iunit,1) qtitle,kp,iskip,(i,i=ibeg,iend,2*iskip)

    1 format('0',a,'   k=',i3,' iskip=',i2,/,' ',41i6) 
      fk=10.0**kp
!
!  quickprint field
!
      do  jli=jend_qp,jbeg,-iskip
         ii= 0
         if(kp.eq.0) then 
          do i=ibeg,iend,iskip
            ii=ii+1
            ix(ii)=a(i,jli)+sign(half,a(i,jli))
          end do
        else
          do i=ibeg,iend,iskip
            ii=ii+1
            ix(ii)=a(i,jli)*fk+sign(half,a(i,jli))
          end do
        end if
        write(iunit,'(i4,81i3)') jli,(ix(i),i=1,ii),jli
     enddo
     return
     
end subroutine qprntn

subroutine smth2d(a,ni,nj,ib,ie,jb,je, &
     anu,npass,nnu,ioresp,io,iskip,dx,b,undef)

!...      routine to smooth a 2-d field at subsection of interior points
!...      using a noncomplex shuman (1957) smoother-desmoother 

      real (kind=4) a(ni,nj),b(ni,nj),anu(nnu)
      real (kind=4) pi,rlambda,dx,undef

      logical ioresp,io
      character qtitle*24
      
!...     output unsmoothed field if io.ne.0

      if(io) then
        call stat2(a,ni,nj,amin,amax,amean,avar,asd)
        write(6,12) amean,amin,amax,avar,asd
12      format(' ',/,' ',' input field mean = ',1pe13.4,/ &
            ' ','             amin = ',1pe13.4,/ &
            ' ','             amax = ',1pe13.4,/ &
            ' ','         variance = ',1pe13.4,/ &
            ' ','         stnd dev = ',1pe13.4,//) 
        qtitle='raw field                '
        call qprntn(a,qtitle,1,1,ni,nj,iskip,6)

      end if

!         mmmmmmmmmmmmmmmmm main loops, npass, the nus 

      do nn=1,npass

        do l=1,nnu

          do i=ib,ie
            do j=jb,je

              if( &
                  a(i,j).eq.undef.or. &
                  a(i+1,j).eq.undef.or. &
                  a(i-1,j).eq.undef.or. &
                  a(i,j-1).eq.undef.or. &
                  a(i,j+1).eq.undef.or. &
                  a(i+1,j+1).eq.undef.or. &
                  a(i+1,j-1).eq.undef.or. &
                  a(i-1,j-1).eq.undef.or. &
                  a(i-1,j+1).eq.undef &
                  ) then
                b(i,j)=a(i,j)

              else

                b(i,j)=a(i,j)*(1.0-anu(l))**2 &
                    + 0.5*anu(l)*(1.0-anu(l))* &
                    (a(i+1,j)+a(i-1,j)+a(i,j+1)+a(i,j-1)) &
                    + 0.25*(anu(l)**2)* &
                    (a(i-1,j-1)+a(i-1,j+1)+a(i+1,j-1)+a(i+1,j+1))
             endif
            end do
          end do

          do i=ib,ie
            do j=jb,je
              a(i,j)=b(i,j)
            end do
          end do

        end do 

      end do

      if(ioresp) then

        write(6,200) npass,nnu 
200     format(' ',//,' ','smoothing function analysis'/ &
         ' ',5x,'number of passes = ',i2/ &
         ' ',5x,'number of elements per pass = ',i2)
        do k=1,nnu
          write(6,201) k,anu(k)
201       format(' ',7x,'k = ',i2, &
       '  smoothing coefficient nu = ',f6.3)
        end do

        pi=4.0*atan(1.0)

        do i=2,ni
          b(i,1)=float(i)
          b(i,2)=1.0
          do mm=1,nnu
            b(i,2)=b(i,2)*(1.0-anu(mm)*(1.0-cos(2.0*pi/float(i))))
          end do

          b(i,2)=b(i,2)**npass
        end do

!
        write(6,222)
222     format(' ','response function as a function of wavelength ', &
         'in grid units*dx',//, &
         ' ','  lambda  response  ',//)
!
        do i=2,ni
          rlambda=dx*i
          write(6,225) rlambda,b(i,2)
225       format(' ',f7.1,3x,f6.3)
        end do

      end if 

      return

end subroutine smth2d

end module libmf




program fimtopo

  use const
  use libmf

  character topodatfile*120,topoglvldir*120

  real(kind=4), allocatable :: z0out(:)

  namelist /TOPOnamelist/ toponpass,toposmoothfact,&
       topodatfile,topoglvldir

  integer        :: glvl                                   ! The grid level defined in the Makefile
  integer        :: SubdivNum(20)                          ! Subdivision numbers for each recursive refinement(2: bisection, 3: trisection, etc.)
  integer        :: nvl                                    ! Number of vertical native levels

  namelist /CNTLnamelist/ glvl, SubdivNum, nvl

  OPEN (10,file="FIMnamelist")
  READ (10,NML=CNTLnamelist)
  close(10)


  OPEN(12,file="fimtopo.nl")
  read(12, nml=TOPOnamelist)
  write(*,nml=TOPOnamelist)
  close(12)


  !
  ! fim properties
  !

  nip=10*(2**glvl)**2+2     ! # of icosahedral point
  if(glvl <= 9) then
     gridscale=15.0*(10-glvl) 
  else
     print*,'EEE invalid glvl must be <= 9'
     stop 'bad glvl'
  endif

  allocate (z0out(nip),stat=irc)

  call mktopo(topodatfile,topoglvldir,glvl,nip,gridscale,toposmoothfact,toponpass, &
       z0out)

  deallocate (z0out)

  stop

!!$  !dddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd
!!$  !
!!$  ! diagnostic code
!!$  !
!!$
!!$  allocate (z0(nip),stat=irc)
!!$
!!$  compdata='wrf'
!!$  compdata='gfs'
!!$
!!$
!!$  open(14,file='g'//cglvl//'z0_'//compdata//'.dat',form='unformatted',status='old',err=801)
!!$  read(14) z0
!!$  close(14) 
!!$  
!!$  open(iunitobs,file='g'//cglvl//'z0'//version//'.obs',form='unformatted',status='unknown',err=812)
!!$
!!$  z0mnc=0.0
!!$  z0maxc=-1e20
!!$  z0minc=1e20
!!$ 
!!$  nmc=0
!!$  do i=1,nip
!!$     if(z0(i) > -1e10 .and. z0(i) < 1e20) then
!!$        nmc=nmc+1
!!$        z0mnc=z0mnc+z0(i)
!!$     endif
!!$
!!$     if(z0(i) < z0minc) z0minc=z0(i)
!!$     if(z0(i) > z0maxc) z0maxc=z0(i)
!!$  enddo
!!$  
!!$  z0mnc=z0mnc/nmc
!!$  
!!$
!!$  dlon=360.0/(ni)
!!$  dlat=180.0/(nj)
!!$  
!!$  print*,'dlon: ',dlon,' dlat: ',dlat,' elon: ',blon+(ni-1)*dlon,' elat: ',blat+(nj-1)*dlat
!!$
!!$  radinfj=(gridscale*smoothfact*km2deglat)/dlat
!!$  radinfj=radinfj*0.5
!!$
!!$  do i=1,nip
!!$     rlat=lat(i)*rad2deg
!!$     rlon=lon(i)*rad2deg
!!$     if(rlon > 180.0) rlon=rlon-360.0
!!$     z0in=z0(i)
!!$     call anltopo(z0(i),topo,rlat,rlon,ni,nj,z0out(i))
!!$     !write(*,'("i:",i6,2x,4f12.3,2x,"ib,ie: ",2(i5,1x),"jb,je: ",2(i5,1x))') i,rlat,rlon,ri,rj,ib,ie,jb,je
!!$  enddo
!!$
!!$
!!$  z0mn=0.0
!!$  z0max=-1e20
!!$  z0min=1e20
!!$
!!$  z0dmn=0.0
!!$  z0dmna=0.0
!!$  z0dmax=-1e20
!!$  z0dmin=1e20
!!$ 
!!$  nm=0
!!$  do i=1,nip
!!$
!!$     if(z0(i) > -1e10 .and. z0(i) < 1e20) then
!!$        dz0=z0(i)-z0out(i)
!!$        nm=nm+1
!!$        z0mn=z0mn+z0out(i)
!!$        z0dmn=z0dmn+dz0
!!$        z0dmna=z0dmna+abs(dz0)
!!$
!!$        if(abs(dz0) > 500.0) then
!!$           write(*,'(a,2x,f9.2,2x,f7.2,1x,f7.2)') 'BBBBBBBBBBBBBBB dz0: ',dz0,lat(i)*rad2deg,lon(i)*rad2deg
!!$        endif
!!$           
!!$        if(dz0 < z0dmin) z0dmin=dz0
!!$        if(dz0 > z0dmax) z0dmax=dz0
!!$        if(z0out(i) < z0min) z0min=z0out(i)
!!$        if(z0out(i) > z0max) z0max=z0out(i)
!!$     endif
!!$
!!$     write(stid,'(a1,i7.7)') 's',i
!!$     rlat=lat(i)*rad2deg
!!$     rlon=lon(i)*rad2deg
!!$     tim = 0.0 
!!$     nlev = 1 
!!$     nflag = 1
!!$     write (iunitobs) stid,rlat,rlon,tim,nlev,nflag 
!!$     write (iunitobs) z0out(i)
!!$
!!$  enddo
!!$
!!$  z0mn=z0mn/nm
!!$  z0dmn=z0dmn/nm
!!$  z0dmna=z0dmna/nm
!!$  
!!$  print*,'               cglvl: ',cglvl
!!$  print*,'           compdata: ',compdata
!!$  print*,'          gridscale: ',gridscale
!!$  print*,'         smoothfact: ',smoothfact
!!$  print*,' radinf in dj units: ',radinfj
!!$  print*,'                nip: ',nip
!!$  print*,'             nmc,nm: ',nmc,nm
!!$  print*
!!$  print*,'            version: ',version
!!$  print*
!!$  print*,'         z0mnc,z0mn: ',z0mnc,z0mn
!!$  print*,'       z0maxc,z0max: ',z0maxc,z0max
!!$  print*,'       z0minc,z0max: ',z0minc,z0min
!!$  print* 
!!$  print*,'              z0dmn: ',z0dmn
!!$  print*,'             z0dmna: ',z0dmna
!!$  print*,'             z0dmax: ',z0dmax
!!$  print*,'             z0dmin: ',z0dmin
!!$
!!$
!!$
!!$  !
!!$  ! finalize .obs
!!$  !
!!$
!!$  nlev = 0 
!!$  write (iunitobs) stid,rlat,rlon,tim,nlev,nflag 
!!$  close(iunitobs)
!!$
!!$

  goto 900
800 continue
  print*,'error in open of g?glvl.dat'
  
801 continue
  print*,'error in open of ',topodatfile
  
  
802 continue
  print*,'read end'
  
803 continue
  print*,'read err'
  
900 continue
  stop

end program fimtopo




subroutine mktopo(topodatfile,glvldir,glvl,nip,gridscale,smoothfact,npass, &
     z0out)

  use const
  use libmf

  parameter(ni=4001,nj=2000,nnu=2)

  character(120)           :: topodatfile,glvldir,tpath,gpath
  character(16)            :: header

  real(kind=4), allocatable :: lat(:),lon(:)

  real(4) topo(ni,nj),dum(ni,nj),anu(nnu),z0out(nip)

  integer        :: glvl

  logical iosmth2dresp,iosmth2d,diag

  blon=-179.9550
  blat=-89.9550

  dlon=360.0/(ni)
  dlat=180.0/(nj)
  
  print*,'dlon: ',dlon,' dlat: ',dlat,' elon: ',blon+(ni-1)*dlon,' elat: ',blat+(nj-1)*dlat

  !
  ! allocate arrays
  !

  allocate (lat(nip),lon(nip),stat=irc)

  print*, 'allocate irc = ',irc,' nip: ',nip,gridscale,smoothfact

  iunittopo=10
  iunitglvl=12

  
  gpath=trim(glvldir)//'glvl.dat'
  tpath=trim(topodatfile)

  print*,'iiiiiiiiii ',gpath
  print*,'tttttttttt ',tpath

  !
  !  read in the glvl.dat for lat/lons
  !
  open(iunitglvl,file=gpath,form='unformatted',status='old',err=800)
  
  read(iunitglvl) header
  print*,'h1 ',header
  read(iunitglvl) header

  print*,'h2 ',header
  read(iunitglvl) lat
  read(iunitglvl) lon
  close(iunitglvl)

  print*,'HHHHHHHHHHHHHHHHH ',header

  !
  !  read in 5' wrf topo file
  ! 

  open(iunittopo,file=tpath,form='unformatted',status='old',err=801)
  read(iunittopo) topo

  !
  ! smooth topo
  !

  if(npass > 0) then

     ib=1
     ie=ni
     jb=1
     je=nj
     anu(1)=0.5
     anu(2)=0.5
     undef=1e20
     iskip=1
     dx=gridscale
     
     call smth2d(topo,ni,nj,ib,ie,jb,je, &
          anu,npass,nnu,iosmth2dresp,iosmth2d,iskip,dx,dum,undef)
  endif

  !
  ! analyze the topo to the icos grid
  !

  radinf=gridscale*smoothfact*0.5
  radinfj=(radinf*km2deglat)/dlat

  print*,'nnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnn nip:',nip,rad2deg
  do i=1,nip
     rlat=lat(i)*rad2deg
     rlon=lon(i)*rad2deg
     if(rlon > 180.0) rlon=rlon-360.0
     z0comp=0.0
     call anltopo(z0comp,topo,rlat,rlon,ni,nj,z0out(i), &
     blat,dlat,blon,dlon,radinf,radinfj)
  enddo


  z0outmean=0.0
  z0outmax=-1e20
  z0outmin=1e20
 
  nm=0
  do i=1,nip
     if(z0out(i) > -1e10 .and. z0out(i) < 1e20) then
        nm=nm+1
        z0outmean=z0outmean+z0out(i)
     endif

     if(z0out(i) < z0outmin) z0outmin=z0out(i)
     if(z0out(i) > z0outmax) z0outmax=z0out(i)
  enddo
  
  z0outmean=z0outmean/nm
  
  print*,'               glvl: ',glvl
  print*,'          gridscale: ',gridscale
  print*,'         smoothfact: ',smoothfact
  print*,' radinf in dj units: ',radinfj
  print*,'                nip: ',nip
  print*,'                 nm: ',nm
  print*
  print*,'            version: ',version
  print*
  print*,'          z0outmean: ',z0outmean
  print*,'           z0outmax: ',z0outmax
  print*,'           z0outmin: ',z0outmin
  print*

  deallocate(lat,lon)

  goto 900
800 continue
  print*,'error in open of g?glvl.dat'
  
801 continue
  print*,'error in open of ',topodatfile
  
  
900 continue
  stop
  return

  
end subroutine mktopo



subroutine anltopo(z0test,topo,rlat,rlon,ni,nj,z0out, &
     blat,dlat,blon,dlon,radinf,radinfj)

  use const

  real(4) topo(ni,nj),z0s(ni*20)

  integer verb

  verb=0

  ri=(rlon-blon)/dlon+1.0
  rj=(rlat-blat)/dlat+1.0
  
  rlatfact=cos(rlat*deg2rad)

  radinfi=0.0
  if(rlatfact > 0.0) radinfi=radinfj/rlatfact
  
  rjb=rj-radinfj-1.0
  rje=rj+radinfj+1.0

  if(rjb < 1) rjb=1.0
  if(rje > nj) rje=nj
  
  if(radinfi == 0) then
     rib=1
     rie=ni
  else
     rib=ri-radinfi-1.0
     rie=ri+radinfi+1.0
  endif

  if(rlatfact == 0.0) then
     rib=1.0
     rie=ni/2
  else
     if(rib < 1) rib=1.0
     if(rie > ni) rie=ni
  endif
  
  ib=nint(rib)
  ie=nint(rie)

  jb=nint(rjb)
  je=nint(rje)

  if(ib < 1) ib=1
  if(ie > ni) ie=ni

  if(jb < 1) jb=1
  if(je > nj) je=nj


  z0bar=0.0
  z0rms=0.0
  nz0=0


  do ii=ib,ie
     do jj=jb,je
        tlat=blat+(jj-1)*dlat
        tlon=blon+(ii-1)*dlon
        dy=(tlat-rlat)
        dx=(tlon-rlon)*rlatfact
        dist=sqrt(dx*dx+dy*dy)*deglat2km
        if(dist < distmin) then
           distmin=dist
           z0min=topo(ii,jj)
        endif

        if(dist <= radinf) then
           nz0=nz0+1
           z0s(nz0)=topo(ii,jj)
        endif

     end do
  end do

  do n=1,nz0
     z0bar=z0bar+z0s(n)
  enddo
  
  if(nz0 > 0) then
     z0bar=z0bar/nz0
  else
     print*,'in analtopo no points! nz0 = 0',rlat,rlon
     stop 'no points'
  endif

  if(verb == 1) then
     dz0=z0test-z0bar
     write(*,'("final",2x,2(f7.2,1x),2x,i6,1x,2(f7.2,1x),2x,(2(f7.2,1x)),2x,4(i5,1x))') &
          rlat,rlon,nz0,z0bar,dz0,distmin,z0min,ib,(ie-ib),jb,(je-jb)
  endif

  z0out=z0bar

  return
  
end subroutine anltopo
