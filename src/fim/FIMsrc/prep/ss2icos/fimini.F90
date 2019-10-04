
!*********************************************************************
!     fimini    
!       Initialization subroutine to restep/remap GFS data for FIM  
!       initial condition. 
!     
!       R. Bleck                April,2007       
!       J. Lee                  July,2007
!       R. Bleck                Aug. 2008
!*********************************************************************

  subroutine fimini(nlyr,ter_in,psf_in,gz_in,p_in,t_in,q_in,		&
                    u_in,v_in,o3_in,qc_in,				&
                    us3d,vs3d,dp3d,mp3d,pr3d,ex3d,ph3d,tr3d)

  use module_control  ,only: nip,nvl,nvlp1,glvl,kbl,ntra,ntrb,npp,	&
                             PrintDiags,PrintIpnDiag,pure_sig,EnKFAnl
  use module_constants,only: p1000,rd,cp,qvmin,deg_lat
  use findmaxmin2
  use findmaxmin3

  implicit none
  integer,intent(IN) :: nlyr		   ! number of input layers
!SMS$DISTRIBUTE(dh,NIP) BEGIN
  real,intent(IN ) :: ter_in(      nip)	   ! surface geopotential
  real,intent(IN ) :: psf_in(      nip)	   ! surface pressure
  real,intent(IN ) :: gz_in(nlyr+1,nip)	   ! geopotential
  real,intent(IN ) :: p_in (nlyr+1,nip)	   ! interface pressure
  real,intent(IN ) :: t_in (nlyr  ,nip)	   ! layer temperature
  real,intent(IN ) :: q_in (nlyr  ,nip)	   ! layer rel humid
  real,intent(IN ) :: qc_in (nlyr  ,nip)   ! layer cloud condensate
  real,intent(IN ) :: u_in (nlyr  ,nip)	   ! layer u wind
  real,intent(IN ) :: v_in (nlyr  ,nip)	   ! layer v wind
  real,intent(IN ) :: o3_in(nlyr  ,nip)	   ! layer ozone
  real,intent(OUT) :: us3d (nvl   ,nip)	   ! zonal wind (m/s), layer
  real,intent(OUT) :: vs3d (nvl   ,nip)	   ! meridional wind (m/s), layer
  real,intent(OUT) :: dp3d (nvl   ,nip)	   ! layer thickness (pascal)
  real,intent(OUT) :: mp3d (nvl   ,nip)	   ! Montgomery Potential (m^2/s^2)
  real,intent(OUT) :: pr3d (nvlp1 ,nip)	   ! pressure (pascal)
  real,intent(OUT) :: ex3d (nvlp1 ,nip)	   ! exner function
  real,intent(OUT) :: ph3d (nvlp1 ,nip)    ! geopotential (=gz), m^2/s^2
  real,intent(OUT) :: tr3d (nvl   ,nip,ntra+ntrb)! 1=pot.temp
                                           ! 2=water vapor
                                           ! 3=cloud water/condensate
                                           ! 4=ozone

  real, dimension(nvl  ,nip)     :: uswrk,vswrk,thwrk,qvwrk,o3wrk,qcwrk
  real, dimension(nvlp1,nip)     :: exwrk

! --- variables on spherical grid in sigma layers:
  real :: th_in (nlyr,nip)
  real :: targ_in(nvl,nip)

! --- variables on spherical grid at sigma levels:
  real :: exlev_in(nlyr+1,nip)

!SMS$DISTRIBUTE END

  integer :: ipn
  integer :: k,k1,k2

  real :: theta_lyrs(nvl)
  real :: ex_top

  character(len=16)  :: string
  logical :: vrbos                ! switch for 'verbose' mode

!<><><><><><><><><><><><><><>><><><><><><><><><><><><><><><><><><><><><><>
!     e x t e r n a l   g r i d d e d   d a t a    i n g e s t
!        f o r   F I M   i n i t i a l   c o n d i t i o n s
!<><><><><><><><><><><><><><>><><><><><><><><><><><><><><><><><><><><><><>

!SMS$PARALLEL (dh,ipn) BEGIN

  open(10, file="theta_coor.txt",form="formatted",action="read")
  read(10,*) theta_lyrs(1:nvl)
  close(10)

! --- diagnostic output
  if (PrintDiags) then
!SMS$SERIAL  BEGIN
   print 102,'terrain hgt min,max:',minval(ter_in),maxval(ter_in)
   print 102,'surface prs min,max:',minval(psf_in),maxval(psf_in)
   print 102,'geopot      min,max:',minval(gz_in),maxval(gz_in)
   print 102,'pressure    min,max:',minval(p_in),maxval(p_in)
   print 102,'temperature min,max:',minval(t_in),maxval(t_in)
   print 102,'moisture    min,max:',minval(q_in),maxval(q_in)
   print 102,'condensate  min,max:',minval(qc_in),maxval(qc_in)
   print 102,'u velocity  min,max:',minval(u_in),maxval(u_in)
   print 102,'v velocity  min,max:',minval(v_in),maxval(v_in)
   print 102,'   O3       min,max:',minval(o3_in),maxval(o3_in)
!SMS$SERIAL  END
  end if
102  format (a,2es15.5)

! --- compute Exner function and potential temperature in input layers

  do ipn = 1, nip		! horizontal loop
   vrbos=ipn.eq.PrintIpnDiag

   do k = 1, nlyr+1
    exlev_in(k,ipn)=cp*(p_in(k,ipn)/p1000)**(rd/cp)
   end do

! --- infer virt.pot.temperature from geopotential and Exner fct
   do k=1,nlyr
    if (exlev_in(k,ipn).gt.exlev_in(k+1,ipn)+.01) then
     th_in(k,ipn)=(gz_in(k+1,ipn)-gz_in(k,ipn))			&
      /(exlev_in(k,ipn)-exlev_in(k+1,ipn))
    else
     th_in(k,ipn)=t_in(k,ipn)*(p1000/p_in(k,ipn))**(rd/cp)
    end if
   end do

   if (vrbos) then 
!SMS$ignore begin
    print '(a,i8/7x,a)','fimini  input at ipn =',ipn,			&
     'geopot      pres    exner     temp    theta     moist      u      v'
    print 103,0,ter_in(ipn)
    print 103,(k,gz_in(k,ipn),p_in(k,ipn),exlev_in(k,ipn),		&
     t_in(k,ipn),th_in(k,ipn),q_in(k,ipn),u_in(k,ipn), v_in(k,ipn),	&
      k=1,nlyr),nlyr+1,gz_in(nlyr+1,ipn),p_in(nlyr+1,ipn),		&
       exlev_in(nlyr+1,ipn)
103 format (i4,2f10.1,3f9.2,es10.2,2f7.2)
!SMS$ignore end
   end if

  end do			! horizontal loop

  if (PrintDiags) then
   do k=1,nlyr
    write (string,'(a,i3)') 'theta lyr',k
    call findmxmn2(th_in,nlyr,nip,k,string)
    write (string,'(a,i3)') 'u-vel lyr',k
    call findmxmn2(u_in,nlyr,nip,k,string)
    write (string,'(a,i3)') 'v-vel lyr',k
    call findmxmn2(v_in,nlyr,nip,k,string)
    write (string,'(a,i3)') 'moist lyr',k
    call findmxmn2(q_in,nlyr,nip,k,string)
   end do
  end if

!<><><><><><><><><><><><><><>><><><><><><><><><><><><><><><><><><><><><><><>
!     c o n v e r t   t o   h y b r i d - i s e n t r o p i c   l a y e r s
!<><><><><><><><><><><><><><>><><><><><><><><><><><><><><><><><><><><><><><>

  if (pure_sig .or. EnKFAnl) then
   us3d(:,:)   = u_in     (:,:)
   vs3d(:,:)   = v_in     (:,:)
   ex3d(:,:)   = exlev_in (:,:)
   tr3d(:,:,1) = th_in    (:,:)			! virt.pot.temperature
   tr3d(:,:,2) = max(q_in (:,:),qvmin)		! water vapor
   tr3d(:,:,3) = max(qc_in(:,:),0.   )		! liquid water/condensate
   tr3d(:,:,4) = o3_in    (:,:)			! ozone

  else				! grid is hybrid-isentropic
   do ipn=1,nip
    targ_in(:,ipn)=theta_lyrs
   end do
 
! --- convert input stairstep profiles into continuous line segments.
! --- integrate those segments over new layers consistent with chosen
! --- coordinate ('target') values, subject to min thknss constraints.

   call lay2lay(nlyr,exlev_in,th_in,q_in, u_in, v_in, o3_in,qc_in,	&
                nvl ,exwrk   ,thwrk,qvwrk,uswrk,vswrk,o3wrk,qcwrk,	&
                targ_in,nip, PrintDiags)

   us3d(:,:)   = uswrk(:,:)
   vs3d(:,:)   = vswrk(:,:)
   ex3d(:,:)   = exwrk(:,:)
   tr3d(:,:,1) = thwrk(:,:)			! virt.pot.temperature
   tr3d(:,:,2) = max(qvwrk(:,:),qvmin)		! water vapor
   tr3d(:,:,3) = max(qcwrk(:,:),0.   )		! liquid water/condensate
   tr3d(:,:,4) = o3wrk(:,:)			! ozone

! --- add nontrivial content to class B and remaining class A tracer arrays.
! --- for testing purposes, use initl pressure and initial latitude as tracers.
! --- copy class B tracers into class A to provide reference tracer fields.

   do ipn = 1,nip				! horizontal loop
    do k=1,nvl
     if (ntrb.ge.1) then
! --- choose initial pressure as tracer 1
      tr3d(k,ipn,ntra+1)=p1000*(exwrk(k,ipn)/cp)**3.5	! => class B
      if (ntra.ge.5) tr3d(k,ipn,5)=tr3d(k,ipn,ntra+1)	! duplicate in class A
     end if
   
     if (ntrb.ge.2) then
! --- choose initial latitude (+90 for pos.def.) as tracer 2
      tr3d(k,ipn,ntra+2)=deg_lat(ipn) + 90.		! => class B
      if (ntra.ge.6) tr3d(k,ipn,6)=tr3d(k,ipn,ntra+2)	! duplicate in class A
     end if
    end do
   end do
  end if			! pure_sig:  true or false

  do ipn = 1,nip		! horizontal loop
   vrbos=ipn.eq.PrintIpnDiag

   pr3d(nvlp1,ipn)=p1000*(ex3d(nvlp1,ipn)/cp)**(cp/rd)
   do k = nvl,1,-1
    pr3d(k,ipn)=p1000*(ex3d(k,ipn)/cp)**(cp/rd)
    dp3d(k,ipn)=pr3d(k,ipn)-pr3d(k+1,ipn)
   end do

! --- solve hydrostatic eqn for geo- and montgomery potential
   mp3d(1,ipn)=ter_in(ipn)+ex3d(1,ipn)*tr3d(1,ipn,1)	! first layer
   ph3d(1,ipn)=ter_in(ipn)
   do k = 2,nvl
    mp3d(k,ipn)=mp3d(k-1,ipn)+ex3d(k,ipn)*(tr3d(k,ipn,1)-tr3d(k-1,ipn,1))
    ph3d(k,ipn)=mp3d(k,ipn)-ex3d(k,ipn)*tr3d(k,ipn,1)
   end do
   ph3d(nvlp1,ipn)=mp3d(nvl,ipn)-ex3d(nvlp1,ipn)*tr3d(nvl,ipn,1)

   if (vrbos) then
!SMS$ignore begin
    write (6,99) ipn,'o u t p u t     p r o f i l e :'
99  format ('ipn =',i8,'   fimini    ',a/				&
     '(5-line groups: pres, exn.fct, geopot(km), theta, montg.pot/1000)')
    do k2=1,nvl,10
      print '(   -2p,11f7.1)',(pr3d(k1,ipn)  ,k1=k2,min(nvlp1,k2+10) )
      print '(       11f7.1)',(ex3d(k1,ipn)  ,k1=k2,min(nvlp1,k2+10) )
      print '(   -4p,11f7.3)',(ph3d(k1,ipn)  ,k1=k2,min(nvlp1,k2+10) )
      print '(4x,    10f7.2)',(tr3d(k1,ipn,1),k1=k2,min(nvl  ,k2+9) )
      print '(4x,-3p,10f7.2)',(mp3d(k1,ipn)  ,k1=k2,min(nvl  ,k2+9) )
      print *
    end do
!SMS$ignore end
   end if

  end do			! horizontal loop
!SMS$PARALLEL END
 
! --- show max/min of all tracer fields in a single mid-range layer
  k1=nvl/2
  do k2=1,ntra
   write (string,'(a,i2,a,i3)') 'trc A',k2,' lyr',k1
   call findmxmn3(tr3d,nvl,nip,ntra+ntrb,k1,k2,string)
  end do
  do k2=1,ntrb
   write (string,'(a,i2,a,i3)') 'trc B',k2,' lyr',k1
   call findmxmn3(tr3d,nvl,nip,ntra+ntrb,k1,ntra+k2,string)
  end do

  if(PrintDiags) then
!SMS$SERIAL  BEGIN
   print 102,'minmax u      ',minval(us3d(1,:)),maxval(us3d(nvl,:))
   print 102,'minmax v      ',minval(vs3d(1,:)),maxval(vs3d(nvl,:))
   print 102,'minmax dp     ',minval(dp3d(1,:)),maxval(dp3d(nvl,:))
   print 102,'minmax ph     ',minval(ph3d(1,:)),maxval(ph3d(nvl+1,:))
   print 102,'minmax pres   ',minval(pr3d(nvl+1,:)),maxval(pr3d(1,:))
   print 102,'minmax pi     ',minval(ex3d(nvl+1,:)),maxval(ex3d(1,:))
   print 102,'minmax mp     ',minval(mp3d(1,:)),maxval(mp3d(nvl,:))
   print 102,'minmax tr(1)  ',minval(tr3d(:,:,1)),maxval(tr3d(:,:,1))
   print 102,'minmax tr(2)  ',minval(tr3d(:,:,2)),maxval(tr3d(:,:,2))
   print 102,'minmax tr(3)  ',minval(tr3d(:,:,3)),maxval(tr3d(:,:,3))
   print 102,'minmax tr(4)  ',minval(tr3d(:,:,4)),maxval(tr3d(:,:,4))
   print 102,'max ter&ph1 ',maxval(ter_in(:)),maxval(ph3d(1,:))
!SMS$SERIAL  END
  end if

  print *,'... exiting fimini'
  return 
  end subroutine fimini
