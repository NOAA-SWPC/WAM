
! Thanks to Pete Johnsen of Cray for the following optimization which 
! allows 10km FIM to run on as many as 33,000 cores.  This optimization has 
! general utility.  See "pjj" for details.  
!
!pjj
! Cray XT code to reduce amount of memory used by icos grid
! variables.  This splits the 10 variables into 5 sets of 2
! which reuses available memory.

subroutine ss2icos(nvp,sanlFile,us3d,vs3d,dp3d,mp3d,pr3d,ex3d,ph3d,tr3d,gfsltln_file)

! read spherical data (GFS) and perform 2-step transform:
! --- (1) horizontal transform from spherical to icos grid
! --- (2) vertical transform from sigma to hybrid-isentropic coord.

!SMS$ignore begin
  use sigio_module
!SMS$ignore end
  use module_control,only: glvl,nvl,nvlp1,nip,ntra,ntrb,curve,		&
                           NumCacheBLocksPerPE,PrintIpnDiag,		&
                           PrintDiagProgVars,PrintDiagNoise,PrintDiags, &
                           alt_topo,pure_sig
  use module_constants,only: grvity,sigak,sigbk
!SMS$ignore begin
  USE slint, ONLY: bilinear_init
!SMS$ignore end
  use stencilprint
  implicit none

  integer          ,intent(IN)     :: nvp
  CHARACTER(len=80),intent(IN)     :: sanlFile
  CHARACTER(len=80),intent(IN)     :: gfsltln_file

!SMS$DISTRIBUTE(dh,NIP) BEGIN
  real,intent(out) :: us3d(nvl  ,nip)	 ! zonal wind (m/s), layer
  real,intent(out) :: vs3d(nvl  ,nip)	 ! meridional wind (m/s), layer
  real,intent(out) :: dp3d(nvl  ,nip)	 ! del p between coord levels (pascals)
  real,intent(out) :: mp3d(nvl  ,nip)	 ! Montgomery Potential (m^2/s^2)
  real,intent(out) :: pr3d(nvlp1,nip)	 ! pressure (pascal)
  real,intent(out) :: ex3d(nvlp1,nip)	 ! exner function
  real,intent(out) :: ph3d(nvlp1,nip)	 ! geopotential (=gz), m^2/s^2
  real,intent(out) :: tr3d(nvl  ,nip,ntra+ntrb)! 1=pot.temp, 2=water vapor, 3=cloud condensate, 4=ozone

  real(4)          :: hs_lev(      nip)    ! surface height (m)
  real(4)          :: ps_lev(      nip)    ! surface pressure in pascals
  real(4)          :: t_lyr (nvp  ,nip)    ! temperature in Kelvins
  real(4)          :: qv_lyr(nvp  ,nip)    ! specific humidity
  real(4)          :: qc_lyr(nvp  ,nip)    ! cloud condensate
  real(4)          :: u_lyr (nvp  ,nip)    ! zonal velocity
  real(4)          :: v_lyr (nvp  ,nip)    ! meridional velocity
  real(4)          :: o3_lyr(nvp  ,nip)    ! ozone mixing ratio
  real(4)          :: p_lev (nvp+1,nip)    ! interface pressure
  real(4)          :: z_lev (nvp+1,nip)    ! geopotential height 
!SMS$DISTRIBUTE END

  integer :: ipn
  integer(sigio_intkind),parameter :: lusig=82
  type(sigio_head)                 :: head
  type(sigio_data)                 :: data

  real(4) :: sig_lyr(nvp  ) ! sig at layer midpoints and interfaces.
  real(4) :: sig_lev(nvp+1) ! sig at layer midpoints and interfaces.

  integer(sigio_intkind)           :: iret

!  integer                          :: n
! pjj/cray - from ss2g1
   integer :: imax, jmax
! Storage routines for GFS input -- may be used multiply as indicated by '/'
  real(4), allocatable :: &
    f1(:),                &       ! surface pressure
    f2(:),                &       ! surface height (m)
    g1(:,:),              &       ! layer virt.temp/u_wind/ozone/condensate
    g2(:,:),              &       ! layer specif.hum/v_wind
    g3(:,:),              &       ! interface geopotential
    pi(:,:),              &       ! interface pressure
    pl(:,:),              &       ! layer midpoint pressure
    dl(:,:)                       ! sig layer thickness

!SMS$SERIAL (<hs_lev,ps_lev,OUT> : default=ignore)  BEGIN
  call sigio_srohdc(lusig,sanlFile,head,data,iret)
  if (iret .ne. 0) then
     call errmsg('ss2icos: error reading '//sanlFile)
!    call errexit(2)		!  errexit doesn't call MPI_ABORT -> use STOP
     STOP 
  endif
  if(nvp /= head%levs) then
     call errmsg('ss2icos: nvp differs from head%levs')
     print '(a,2i5)','nvp,head%levs =',nvp,head%levs
!    call errexit(2)		!  errexit doesn't call MPI_ABORT -> use STOP
     STOP 
  endif    
  if (pure_sig .and. nvl /= head%levs) then
     call errmsg('ss2icos: in "pure_sig" mode, nvl must match head%levs')
     print '(a,2i5)','nvl,head%levs =',nvl,head%levs
!    call errexit(2)		!  errexit doesn't call MPI_ABORT -> use STOP
     STOP 
  end if

! pjj/cray - from ss2gg1
! grid_file2 = "glvl.dat"
!  imax = 1152
!  jmax = 576

  imax=head%lonb
  jmax=head%latb


!SMS$ignore begin
  OPEN                (66,file="glvl.dat",status='old',form='unformatted')
  call TestGlvlHeader (66,     "glvl.dat",'ss2icos',glvl)
  call TestCurveHeader(66,     "glvl.dat",'ss2icos',curve)
  CALL bilinear_init(gfsltln_file, imax*jmax, 66, nip)
  close(66)
!SMS$ignore end

  ! get surface height & pres (hs_lev,ps_lev) 
  ! interpolated to icos grid.

  allocate(f1(imax*jmax), f2(imax*jmax))

  call ss2gg_xt1(4,imax, jmax, head,data,nip,glvl,curve,nvp,	&
                 hs_lev,ps_lev, f1, f2)

  ! write out GrADS control file.
  call ss2gg2(4,imax, jmax, head,'siganl.ieee',sig_lev,sig_lyr)
!SMS$SERIAL end

  allocate (sigak(nvp+1),sigbk(nvp+1))

!SMS$SERIAL (<z_lev,p_lev,sigak,sigbk,OUT> : default=ignore)  BEGIN

  ! get temp (t_lyr), pres (p_lev), 
  ! spec hum (qv_lyr), cloud condensate (qc_lyr),
  ! geopot height (z_lev) winds (u_lyr,v_lyr) and
  ! ozone (o3_lyr) interpolated to icos grid.

  allocate(                      &  
    g1(imax*jmax,head%levs),     &    ! layer virt.temp/u_wind/ozone/condensate
    g2(imax*jmax,head%levs),     &    ! layer specif.hum/v_wind
    g3(imax*jmax,head%levs+1),   &    ! interface geopotential
    pi(imax*jmax,head%levs+1),   &    ! interface pressure
    pl(imax*jmax,head%levs),     &    ! layer midpoint pressure
    dl(imax*jmax,head%levs) )         ! sig layer thickness

  call ss2gg_xt2(4,imax, jmax, head,data,nip,glvl,curve,nvp,    &
                 sigak,sigbk,z_lev,p_lev, f1,f2,g1,g2,g3,pi,pl,dl )

!SMS$SERIAL end

!SMS$SERIAL (<t_lyr,qv_lyr,OUT> : default=ignore)  BEGIN
  call ss2gg_xt3(4,imax, jmax, head,data,nip,glvl,curve,nvp,    &
                 t_lyr,qv_lyr, g1,g2 )
!SMS$SERIAL end

!SMS$SERIAL (<u_lyr,v_lyr,OUT> : default=ignore)  BEGIN
  call ss2gg_xt4(4,imax, jmax, head,data,nip,glvl,curve,nvp,    &
                 u_lyr,v_lyr, g1,g2 )

! --- reverse velocity vectors at poles
       u_lyr(:,1  )=-u_lyr(:,1  )
       v_lyr(:,1  )=-v_lyr(:,1  )
       u_lyr(:,nip)=-u_lyr(:,nip)
       v_lyr(:,nip)=-v_lyr(:,nip)
!SMS$SERIAL end

!SMS$SERIAL (<o3_lyr,qc_lyr,OUT> : default=ignore)  BEGIN
  call ss2gg_xt5(4,imax, jmax, head,data,nip,glvl,curve,nvp,    &
                 o3_lyr,qc_lyr, g1)
  call sigio_axdata(data,iret)     ! deallocate array
! also done with these
  deallocate (f1)
  deallocate (f2)
  deallocate (g1)
  deallocate (g2)
  deallocate (g3)
  deallocate (pi)
  deallocate (pl)
  deallocate (dl)
!SMS$SERIAL end

  ! horizontal interpolation done.

  if (alt_topo) then
     
!SMS$SERIAL (<hs_lev,INOUT> : default=ignore)  BEGIN
  ! get topo on icos grid
     !call rdtopo(hs_lev,nip)
     call mktopo(hs_lev,nip)
!SMS$SERIAL END

  ! correct zg for new topo
     call ss2ggtopo(nip,nvp,	&
          hs_lev,ps_lev,z_lev,p_lev,t_lyr,qv_lyr,	&
          u_lyr,v_lyr,o3_lyr,qc_lyr)

  end if

!SMS$PARALLEL (dh,ipn) BEGIN
  do ipn=1,nip
    hs_lev(ipn)=hs_lev(ipn)*grvity	! surface height => surface geopot
  enddo
!SMS$PARALLEL END

  call stencl(hs_lev,1,1.,'surface height (m)')

  ! Now do vertical interpolation.

  call fimini(nvp,hs_lev,ps_lev,z_lev,p_lev,t_lyr,qv_lyr,	&
              u_lyr,v_lyr,o3_lyr,qc_lyr,us3d,vs3d,		&
              dp3d,mp3d,pr3d,ex3d,ph3d,tr3d)
  return 
end subroutine ss2icos




! pjj/cray - first set 
subroutine ss2gg_xt1(idrt,imax,jmax,head,data,nip,glvl,curve,nvp,	&
                     hs_lev,ps_lev, f1, f2)

  use module_constants,only: rd,cp,qvmin,grvity,p1000
!SMS$ignore begin
  use sigio_module
  use physcons
  USE slint, ONLY:bl_int
!SMS$ignore end
  implicit none

  integer,intent(in):: idrt,imax,jmax,nip,glvl,curve
  type(sigio_head),intent(in):: head
  type(sigio_data),intent(in):: data
! Storage routines for GFS input -- may be used multiply as indicated by '/'
  real(4) f1(imax*jmax),		&	! surface pressure
          f2(imax*jmax)				! surface height (m)

  integer j,k,k1,ipn,nvp,kgrnd(nip)
  real(4) icos2d(nip),icos2d1(nip)
  real(4),intent(OUT) :: hs_lev(nip),ps_lev(nip)
  real exlo,exup,th_lyr,pkap,zold(nip),znew

! perform spherical transform on surface height (f2) field
  call sptez(0,head%jcap,idrt,imax,jmax,data%hs,f2,1)

! interpolate surface height to icos grid
  CALL bl_int (f2, hs_lev)

! perform spherical transform on surface pressure (f1) field
  call sptez(0,head%jcap,idrt,imax,jmax,data%ps,f1,1)
  f1=exp(f1)*1.e3		! convert ln(ps) in centibars to ps in Pa.

! interpolate surface pressure to icos grid
  CALL bl_int (f1, ps_lev)	! unit in pascal

  print 100,'min,max of srf.height on spherical grid:',			&
      minval(f2),maxval(f2)
  print 100,'min,max of srf.press  on spherical grid:',			&
      minval(f1),maxval(f1)
  print 100,'min,max of srf.height on icos grid:',			&
      minval(hs_lev),maxval(hs_lev)
  print 100,'min,max of srf.press  on icos grid:',			&
      minval(ps_lev),maxval(ps_lev)
 100 format (a,2f13.2)

  return
end subroutine ss2gg_xt1



subroutine ss2gg_xt2(idrt,imax,jmax,head,data,nip,glvl,curve,nvp,	&
                     sigak,sigbk,z_lev,p_lev,f1,f2,g1,g2,g3,pi,pl,dl)
  use module_constants,only: rd,cp,qvmin,grvity,p1000
!SMS$ignore begin
  use sigio_module
  use physcons
  USE slint, ONLY:bl_int
!SMS$ignore end
  implicit none

  integer,intent(in):: idrt,imax,jmax,nip,glvl,curve
  type(sigio_head),intent(inout):: head
  type(sigio_data),intent(in):: data
! Storage routines for GFS input -- may be used multiply as indicated by '/'
  real(4) f1(imax*jmax),		&	! surface pressure
          f2(imax*jmax),		&	! surface height (m)
          g1(imax*jmax,head%levs),	&	! layer virt.temp/u_wind/ozone/condensate
          g2(imax*jmax,head%levs),	&	! layer specif.hum/v_wind
          g3(imax*jmax,head%levs+1),	&	! interface geopotential
          pi(imax*jmax,head%levs+1),	&	! interface pressure
          pl(imax*jmax,head%levs),	&	! layer midpoint pressure
          dl(imax*jmax,head%levs)		! sig layer thickness
  real,intent(OUT) :: sigak(head%levs+1),sigbk(head%levs+1)

  integer j,k,ipn,nvp
  real(4) icos2d(nip),icos2d1(nip), exlo, exup
  real(4),intent(OUT) :: z_lev(nvp+1,nip)
  real(4),intent(OUT) :: p_lev(nvp+1,nip)

! perform spherical transform on virt.temp. (g1) and specif.hum. (g2) field
  call sptezm(0,head%jcap,idrt,imax,jmax,head%levs,data%t,g1,1)
  call sptezm(0,head%jcap,idrt,imax,jmax,head%levs,data%q,g2,1)

! calculate pressure at sigma layer midlevels (pl) and
! interfaces (pi).  dl is pressure drop across layer.

  call modpr(imax*jmax,imax*jmax,head%levs,head%idvc,head%idsl,&
             head%si,head%ak,head%bk,f1,pl,dl)
  if(head%idvc.eq.2) then
    write (*,'(/a)')					&
    'GFS intfc.prs is defined as p(k) = ak(k) + bk(k) * surf.prs'

! avoid zero pressure at top
    head%ak(head%levs+1)=max(head%ak(head%levs+1),.2*head%ak(head%levs))
    head%bk(head%levs+1)=0.

     do k=1,head%levs+1
       pi(:,k)=head%ak(k)+head%bk(k)*f1
       sigak(k)=head%ak(k)
       sigbk(k)=head%bk(k)
     enddo
  else
    write (*,'(/a)')					&
    'GFS intfc.prs is defined as p(k) = bk(k) * surf.prs'
     do k=1,head%levs+1
       pi(:,k)=head%si(k)*f1
       sigak(k)=0.
       sigbk(k)=head%si(k)
     enddo
  endif
  write (*,'(a/(5f14.6))') 'ak array:',(head%ak(k),k=1,head%levs+1)
  write (*,'(a/(5f14.6))') 'bk array:',(head%bk(k),k=1,head%levs+1)

! Compute geopotential (g3) on interfaces (still on spherical grid!)

  do j=1,imax*jmax
   exup=cp*(pi(j,1)/p1000)**(rd/cp)
   g3(j,1)=f2(j)*grvity		! f2 = surface height
   do k=1,head%levs
    exlo=exup
    exup=cp*(pi(j,k+1)/p1000)**(rd/cp)
    g3(j,k+1)=g3(j,k)+(exlo-exup)*g1(j,k)*(p1000/pl(j,k))**(rd/cp)
   end do
  end do

! interpolate geopotential and intfc.pressure to icos grid

  do k=1,head%levs+1
     CALL bl_int (g3(1,k), icos2d)
     do ipn=1,nip
      z_lev(k,ipn)=icos2d(ipn)
     end do
     CALL bl_int (pi(1,k), icos2d)
     do ipn=1,nip
      p_lev(k,ipn)=icos2d(ipn)
     end do
  enddo

  return
end subroutine ss2gg_xt2



subroutine ss2gg_xt3(idrt,imax,jmax,head,data,nip,glvl,curve,nvp,	&
                     t_lyr,qv_lyr, g1,g2)
  use module_constants,only: rd,cp,qvmin,grvity,p1000
!SMS$ignore begin
  use sigio_module
  use physcons
  USE slint, ONLY:bl_int
!SMS$ignore end
  implicit none

  integer,intent(in):: idrt,imax,jmax,nip,glvl,curve
  type(sigio_head),intent(in):: head
  type(sigio_data),intent(in):: data
! Storage routines for GFS input -- may be used multiply as indicated by '/'
  real(4) g1(imax*jmax,head%levs),	&	! layer virt.temp/u_wind/ozone/condensate
          g2(imax*jmax,head%levs)	! layer specif.hum/v_wind

  integer j,k,ipn,nvp
  real(4) icos2d(nip),icos2d1(nip), exlo, exup
  real(4),intent(OUT) :: t_lyr(nvp,nip),qv_lyr(nvp,nip)

! interpolate virt.temp and specif.hum to icos grid

  do k=1,head%levs
     CALL bl_int (g1(1,k), icos2d)
     do ipn=1,nip
      t_lyr(k,ipn)=icos2d(ipn)
     end do
     CALL bl_int (g2(1,k), icos2d)
     do ipn=1,nip
      qv_lyr(k,ipn)=max( icos2d(ipn), qvmin )
     end do
  enddo

  return
end subroutine ss2gg_xt3



subroutine ss2gg_xt4(idrt,imax,jmax,head,data,nip,glvl,curve,nvp,	&
                     u_lyr,v_lyr, g1,g2)
  use module_constants,only: rd,cp,qvmin,grvity,p1000
!SMS$ignore begin
  use sigio_module
  use physcons
  USE slint, ONLY:bl_int
!SMS$ignore end
  implicit none

  integer,intent(in):: idrt,imax,jmax,nip,glvl,curve
  type(sigio_head),intent(in):: head
  type(sigio_data),intent(in):: data
! Storage routines for GFS input -- may be used multiply as indicated by '/'
  real(4) g1(imax*jmax,head%levs),	&	! layer virt.temp/u_wind/ozone/condensate
          g2(imax*jmax,head%levs)	! layer specif.hum/v_wind

  integer j,k,ipn,nvp
  real(4) icos2d(nip),icos2d1(nip), exlo, exup
  real(4),intent(OUT) :: u_lyr(nvp,nip),v_lyr(nvp,nip)

! perform spherical transform on u wind (g1) and v wind (g2) field
  call sptezmv(0,head%jcap,idrt,imax,jmax,head%levs,data%d,data%z,g1,g2,1)

! interpolate u,v to icos grid
  do k=1,head%levs
       CALL bl_int (g1(1,k), icos2d)
       CALL bl_int (g2(1,k), icos2d1)
       do ipn=1,nip
        u_lyr(k,ipn)=icos2d(ipn)
        v_lyr(k,ipn)=icos2d1(ipn)
       end do
  enddo
 
  return
end subroutine ss2gg_xt4



subroutine ss2gg_xt5(idrt,imax,jmax,head,data,nip,glvl,curve,nvp,	&
                     o3_lyr,qc_lyr, g1 )
  use module_constants,only: rd,cp,qvmin,grvity,p1000
!SMS$ignore begin
  use sigio_module
  use physcons
  USE slint, ONLY:bl_int
!SMS$ignore end
  implicit none

  integer,intent(in):: idrt,imax,jmax,nip,glvl,curve
  type(sigio_head),intent(in):: head
  type(sigio_data),intent(in):: data
! Storage routines for GFS input -- may be used multiply as indicated by '/'
  real(4) :: g1(imax*jmax,head%levs)	! layer virt.temp/u_wind/ozone/condensate

  integer j,k,ipn,nvp
! CHARACTER(len=80)   :: grid_file1,grid_file2
  real(4) icos2d(nip),icos2d1(nip), exlo, exup
  real(4),intent(OUT) :: qc_lyr(nvp,nip)
  real(4),intent(OUT) :: o3_lyr(nvp,nip)

! perform spherical transform on ozone (g1) field
  call sptezm(0,head%jcap,idrt,imax,jmax,head%levs,data%q(1,1,2),g1,1)

! interpolate ozone to icos grid
  do k=1,head%levs
       CALL bl_int (g1(1,k), icos2d)
       do ipn=1,nip
        o3_lyr(k,ipn)=max( icos2d(ipn), 0. )
       end do
  enddo

! perform spherical transform on cloud condensate (g1) field
  call sptezm(0,head%jcap,idrt,imax,jmax,head%levs,data%q(1,1,3),g1,1)

! interpolate cloud condensate to icos grid
  do k=1,head%levs
       CALL bl_int (g1(1,k), icos2d)
       do ipn=1,nip
        qc_lyr(k,ipn)=max( icos2d(ipn), 0. )
       end do
  enddo

  return
end subroutine ss2gg_xt5



subroutine ss2gg2(idrt,imax,jmax,head,cfggg,si,sl)
!SMS$ignore begin
  use sigio_module
!SMS$ignore end
  implicit none
  integer,         intent(in) :: idrt,imax,jmax
  type(sigio_head),intent(in) :: head
  character*(*) cfggg
  real(4) slat(jmax),wlat(jmax)
  integer idat(8),jdat(8),jhr
  real(8) rincin(5) !r8 because w3movdat is in libcol where it is r8
  integer idatin(8)
  character*10 cdat(8)
  integer n,luctl,iret
  real  ::  ps(1) = 1.e5
  real(4) sl(head%levs),dl(head%levs),si(head%levs+1)

  luctl=12
  open(luctl,file='siganl.ctl',status='replace',iostat=iret)

  rincin    = 0.
  rincin(2) = head%fhour
  idatin    = 0
  idatin(1) = head%idate(4)
  idatin(2) = head%idate(2)
  idatin(3) = head%idate(3)
  idatin(5) = head%idate(1)
  call w3movdat(rincin,idatin,idat)
  call w3pradat(idat,cdat)
  jhr=12 ! ?? jsw - what is this for ??

  call modpr(1,1,head%levs,head%idvc,head%idsl,&
             head%si,head%ak,head%bk,ps,sl,dl)
  sl=sl/1.e5
  dl=dl/1.e5
  do n=1,head%levs+1
     si(n)=(head%ak(n) + head%bk(n)*ps(1))/1.e5 
  enddo
  if(cfggg(1:1).eq.'/') then
    write(luctl,'("dset ",a)') cfggg
  else
    write(luctl,'("dset ^",a)') cfggg
  endif
  write(luctl,'("options yrev")')
  write(luctl,'("undef -9.99E+33")')
  write(luctl,'("title ss2icos")')
  write(luctl,'("xdef",i6," linear",2f12.6)') imax,0.d0,360.d0/imax
  if(idrt.eq.0) then
    write(luctl,'("ydef",i6," linear",2f12.6)')&
     jmax,-90.d0,180.d0/(jmax-1)
  elseif(idrt.eq.256) then
    write(luctl,'("ydef",i6," linear",2f12.6)')&
     jmax,-90.d0*(jmax-1)/jmax,180.d0/jmax
  elseif(idrt.eq.4) then
    call splat(idrt,jmax,slat,wlat)
    write(luctl,'("ydef",i6," levels")') jmax
    write(luctl,'(5f12.6)') 180.d0/acos(-1.d0)*asin(dble(slat(jmax:1:-1)))
  endif
  write(luctl,'("zdef",i6," levels")') head%levs
  write(luctl,'(5f12.6)') sl
  write(luctl,'("tdef",i6," linear ",i2.2,"Z",i2.2,a3,i4.4,1x,i6,"hr")')&
   1,idat(5),idat(3),cdat(2)(1:3),idat(1),jhr
  write(luctl,'("vars",i6)') 10+head%ntrac
  write(luctl,'("HS  ",i3," 99 surface orography (m)")') 1
  write(luctl,'("PS  ",i3," 99 surface pressure (Pa)")') 1
  write(luctl,'("P   ",i3," 99 pressure (Pa)")') head%levs
  write(luctl,'("DP  ",i3," 99 delta pressure (Pa)")') head%levs
  write(luctl,'("T   ",i3," 99 temperature (K)")') head%levs
  write(luctl,'("Q   ",i3," 99 specific humidity (kg/kg)")') head%levs
  write(luctl,'("RH  ",i3," 99 relative humidity (%)")') head%levs
  write(luctl,'("U   ",i3," 99 zonal wind (m/s)")') head%levs
  write(luctl,'("V   ",i3," 99 meridional wind (m/s)")') head%levs
  write(luctl,'("DIV ",i3," 99 divergence (m/s**2)")') head%levs
  write(luctl,'("VOR ",i3," 99 vorticity (m/s**2)")') head%levs
  do n=2,min(head%ntrac,9)
    write(luctl,'("Q",i1,2x,i3," 99 tracer ",i1," (kg/kg)")') n,head%levs,n
  enddo
  do n=10,head%ntrac
    write(luctl,'("Q",i2,1x,i3," 99 tracer ",i2," (kg/kg)")') n,head%levs,n
  enddo
  write(luctl,'("endvars")')
  close (luctl)
end subroutine ss2gg2


subroutine modpr(im,ix,km,idvc,idsl,si,ak,bk,ps,pm,pd)
!$$$  subprogram documentation block
!
! subprogram:    modpr       compute model pressures
!   prgmmr: iredell          org: w/nmc23     date: 92-10-31
!
! abstract: compute model pressures.
!
! program history log:
! 2001-07-25  mark iredell
!
! usage:    call modpr(im,ix,km,idvc,idsl,si,ak,bk,ps,pm,pd)
!   input argument list:
!     im           integer number of points to compute
!     ix           integer first dimension
!     km           integer number of levels
!     idvc         integer vertical coordinate id
!                  (1 for sigma and 2 for hybrid)
!     idsl         integer type of sigma structure
!                  (1 for phillips or 2 for mean)
!     si           real (km+1) sigma interface values (idvc=1)
!     ak           real (km+1) hybrid interface a (idvc=2)
!     bk           real (km+1) hybrid interface b (idvc=2)
!     ps           real (ix) surface pressure (pa)
!   output argument list:
!     pm           real (ix,km) mid-layer pressure (pa)
!     pd           real (ix,km) delta pressure (pa)
!
! attributes:
!   language: fortran
!
!$$$
  implicit none
  integer,intent(in):: im,ix,km,idvc,idsl
  real,intent(in):: si(km+1),ak(km+1),bk(km+1),ps(im)
  real,intent(out):: pm(ix,km),pd(ix,km)
  real,parameter:: rocp=287.05/1004.6,rocp1=rocp+1,rocpr=1/rocp
  real pid,piu
  integer i,k
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  do k=1,km
    do i=1,im
      if(idvc.eq.2) then
        pid=ak(k)+bk(k)*ps(i)
        piu=ak(k+1)+bk(k+1)*ps(i)
      else
        pid=si(k)*ps(i)
        piu=si(k+1)*ps(i)
      endif
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if(idsl.eq.2) then
        pm(i,k)=(pid+piu)/2
      else
        pm(i,k)=((pid**rocp1-piu**rocp1)/(rocp1*(pid-piu)))**rocpr
      endif
      pd(i,k)=pid-piu
    enddo
  enddo
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
end subroutine modpr



subroutine ss2ggtopo(nip,nvp,	&
     hs_lev,ps_lev,z_lev,p_lev,t_lyr,qv_lyr,	&
     u_lyr,v_lyr,o3_lyr,qc_lyr)

  use module_constants,only: rd,cp,qvmin,grvity,p1000

  implicit none

!SMS$DISTRIBUTE(dh,nip) BEGIN
  integer j,k,k1,ipn,nip,nvp,kgrnd(nip)

  real(4) :: hs_lev(nip),ps_lev(nip),z_lev(nvp+1,nip),	&
       p_lev(nvp+1,nip),t_lyr(nvp,nip),		&
       qv_lyr(nvp,nip),qc_lyr(nvp,nip),		&
       u_lyr(nvp,nip),v_lyr(nvp,nip),o3_lyr(nvp,nip)

  real exlo,exup,th_lyr,pkap,zold(nip),znew
  real hmax,hmin,pmax,pmin
!SMS$DISTRIBUTE END

!SMS$PARALLEL (dh,ipn) BEGIN
   print *,'switch to non-GFS surface height (topo dat file) ....'
   do ipn=1,nip			! horiz. loop
    do k=1,nvp			! vert. loop
     if (z_lev(k+1,ipn).gt.hs_lev(ipn)*grvity) then

! --- level k+1 is above ground. integrate hydrostat.eqn down from there.
! --- sequence of operations (layer k is sandwiched between interfaces k,k+1):
! --- (a) get old midlayer p^kappa from (partial p^(1+kappa))/(partial p)
! --- (b) get old theta from
! ---     partial phi_old / partial pi_old = -theta_old
! --- (c) set theta_new = theta_old (not optimal, but tolerable for now)
! --- (d) get new bottom pressure from
! ---     partial phi_new / partial pi_new = -theta_new
! --- (e) get new surf.temp. from theta_new and new bottom pressure

      exup=cp*(p_lev(k+1,ipn)/p1000)**(rd/cp)
      exlo=cp*(p_lev(k  ,ipn)/p1000)**(rd/cp)
!      pkap=(exlo*p_lev(k,ipn)-exup*p_lev(k+1,ipn))/		&
!       ((rd+cp)*(p_lev(k,ipn)-     p_lev(k+1,ipn)))
      th_lyr=(z_lev(k+1,ipn)-z_lev(k,ipn))/(exlo-exup)
      exlo=exup+(z_lev(k+1,ipn)-hs_lev(ipn)*grvity)/th_lyr
      p_lev(1,ipn)=p1000*(exlo/cp)**(cp/rd)		! new srf.pres.
      ps_lev(ipn)=p_lev(1,ipn)
      if (exlo.gt.exup+.01) then
       pkap=(exlo*p_lev(1,ipn)-exup*p_lev(k+1,ipn))/		&
        ((rd+cp)*(p_lev(1,ipn)-     p_lev(k+1,ipn)))
      else
       pkap=.5*(exlo+exup)/cp
      end if
!      pkap=(exlo*p_lev(1,ipn)-exup*p_lev(k+1,ipn))/		&
!       ((rd+cp)*(p_lev(1,ipn)-     p_lev(k+1,ipn)))
      t_lyr(1,ipn)=th_lyr*pkap				! new srf.temp.
      z_lev(1,ipn)=hs_lev(ipn)*grvity			! new srf.geopot.
      do k1=2,k
       p_lev(k1,ipn)=p_lev(1,ipn)
       t_lyr(k1,ipn)=t_lyr(1,ipn)
       z_lev(k1,ipn)=z_lev(1,ipn)
      end do
      kgrnd(ipn)=k
      zold(ipn)=z_lev(k+1,ipn)
      exit
     end if
    end do                      ! vert. loop
   end do                       ! horiz. loop

   hmin = minval(hs_lev(1:nip))
   hmax = maxval(hs_lev(1:nip))
   pmin = minval(ps_lev(1:nip))
   pmax = maxval(ps_lev(1:nip))
!SMS$REDUCE(hmax,pmax,max)
!SMS$REDUCE(hmin,pmin,min)

   print 100,'min,max of new srf.height on icos grid:',hmin,hmax
   print 100,'min,max of new srf.press  on icos grid:',pmin,pmax
 100 format (a,2f13.2)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! optional: re-compute geopotential on interfaces to check for errors
   do ipn=1,nip
    znew=z_lev(1,ipn)
    exup=cp*(p_lev(1,ipn)/p1000)**(rd/cp)
    do k=1,kgrnd(ipn)
     exlo=exup
     exup=cp*(p_lev(k+1,ipn)/p1000)**(rd/cp)
      if (exlo.gt.exup+.01) then
       pkap=(exlo*p_lev(1,ipn)-exup*p_lev(k+1,ipn))/		&
        ((rd+cp)*(p_lev(1,ipn)-     p_lev(k+1,ipn)))
      else
       pkap=.5*(exlo+exup)/cp
      end if
     !pkap=(exlo*p_lev(k,ipn)-exup*p_lev(k+1,ipn))/			&
     ! ((rd+cp)*(p_lev(k,ipn)-     p_lev(k+1,ipn)))
     znew=znew+(exlo-exup)*t_lyr(k,ipn)/pkap
    end do
    k=kgrnd(ipn)
    if (abs(znew-zold(ipn)).gt.1.)					&
       print '(a,2i7,f11.1,f9.1)',					&
       'height discrepancy at ipn,kgrnd =',ipn,k,zold(ipn),znew
   end do
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!SMS$PARALLEL END

  return
end subroutine ss2ggtopo





