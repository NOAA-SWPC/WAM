module module_hybgen
use findmaxmin2
contains
!***********************************************************************
!     hybgen
!       Adapted from 2nd generation restep-based grid generator in HYCOM
!       R. Bleck     Mar 2006
!       R. Bleck     Mar 2008       arbitrary number of tracers
!       R. Bleck     Jul 2009       inflation of top layers
!       R. Bleck     Nov 2009       merged with hybgen_sig
!       S. Sun       Nov 2009       added option to smooth th-dot
!       S. Sun       Nov 2009       added lateral intfc smoothing option
!***********************************************************************

  subroutine hybgen (its,   &
            targt,          &	! target potential temperature (1-D)
            us3d,vs3d,      &	! horiz.velocity (3-D)
            tr3d,           &	! mass field tracers (ntra 3-D fields)
            sdot,           &	! intfc displacement sdot*(dp/ds)*dt (3-D)
            ex3d,dp3d,pr3d  )	! Exner function, layer thknss, pressure (3-D)
 
  use module_control  ,only: nvl,nvlp1,ntra,kbl,nip,dt,			&
                             PrintIpnDiag,PrintDiags,			& 
                             intfc_smooth,pure_sig
  use module_constants,only: p1000,rd,cp,sigak,sigbk
  use module_variables,only: worka,workb
  use module_dffusn_lev
  use module_dffusn_lyr
  implicit none

! Type and dimension external variables:

  integer,intent(IN)     :: its			! model time step
  real   ,intent(IN)     :: targt(nvl)	        ! target pot.temp.
!SMS$DISTRIBUTE (dh,nip) BEGIN
  real   ,intent(INOUT)  :: tr3d(nvl,nip,ntra)	! mass field tracers
  real   ,intent(INOUT)  :: ex3d(nvlp1,nip)	! Exner function
  real   ,intent(INOUT)  :: dp3d(nvl  ,nip)	! layer thickness
  real   ,intent(INOUT)  :: pr3d(nvlp1,nip)	! pressure
  real   ,intent(OUT)    :: sdot(nvlp1,nip)	! sdot*(dp/ds)
  real   ,intent(INOUT)  :: us3d(nvl  ,nip)	! u velocity
  real   ,intent(INOUT)  :: vs3d(nvl  ,nip)	! v velocity

! Type and dimension of local variables:

  real    :: exsmo3d(nvlp1,nip)			! 3d array for smoothing
  real    :: exdif3d(nvl,nip)
!SMS$DISTRIBUTE END
  integer :: ipn,k,k1,k2
  logical :: vrbos				! switch for 'verbose' mode
  real    :: ucol(nvl),vcol(nvl)		! velocity column vectors
  real    :: dpcol(nvl),prcol(nvlp1)		! thknss, pres column vectors
  real    :: excol(nvlp1)			! Exner fcn column vector
  real    :: exsmo(nvlp1)			! smoothed Exner fcn col.vector
  real    :: trcol(nvl,ntra)			! tracer column vectors
  real    :: thcol(nvl)				! theta col.vector for remap
  real    :: valmin
  logical :: thsmoo = .false.			! use smoothed thdot in thcol
  character :: string*20
  real    :: smoo_coeff(nvlp1),taper

! taper(k1,k2)=float(k2-k1)**2/(19.+float(k2-k1)**2)	! range: 0.05...1
  taper(k1,k2)=float(k2-k1)**2/( 9.+float(k2-k1)**2)	! range: 0.1....1
! taper(k1,k2)=float(k2-k1)**2/( 4.+float(k2-k1)**2)	! range: 0.2....1
! taper(k1,k2)=float(k2-k1)**4/(19.+float(k2-k1)**4)	! range: 0.05...1
! taper(k1,k2)=float(k2-k1)**4/( 9.+float(k2-k1)**4)	! range: 0.1....1
! taper(k1,k2)=float(k2-k1)**4/( 4.+float(k2-k1)**4)	! range: 0.2....1

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! --- if thsmoo=true, smooth diabatic theta increment (thdot) laterally.
! --- the smoothed field is used   o n l y   to determine coordinate movement
! --- (regridding). the regular theta field (tracer 1) is not smoothed.

  if (pure_sig) then
   thsmoo      =.false.
   intfc_smooth=0.
  end if

  if (thsmoo .and. its.gt.0) then
!SMS$PARALLEL(dh, ipn) BEGIN
   workb(:,:)=tr3d(:,:,1)-worka(:,:)	! worka = theta before physics call
!SMS$PARALLEL END

   call dffusn_lyr(workb,dp3d,dt*20.)		! lateral thdot smoothing

!SMS$PARALLEL(dh, ipn) BEGIN
   worka(:,:)=worka(:,:)+workb(:,:)
!SMS$PARALLEL END
  else						! thsmoo = false
!SMS$PARALLEL(dh, ipn) BEGIN
   worka(:,:)=tr3d(:,:,1)
!SMS$PARALLEL END
  end if					! thsmoo
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
!SMS$PARALLEL(dh, ipn) BEGIN
!sms$compare_var(ex3d, "hybgen.F90 - ex3d7 ")

  do ipn=1,nip					! loop over icos grid
   vrbos=ipn.eq.PrintIpnDiag

!   if (vrbos .and. pure_sig) then
!!SMS$IGNORE BEGIN
!    write (6,'(a/(5f14.6))') '(hybgen) sigak array:',sigak
!    write (6,'(a/(5f14.6))') '(hybgen) sigbk array:',sigbk
!!SMS$IGNORE END
!   end if

! --- call single-column version of hybgen
  
! --- interface smoothing requires separation of regridding and remappping.
! --- step 1: subr. regrid_1d does the regridding 
! --- step 2: subr. dffusn_lev does the (optional) interface smoothing
! --- step 3: subr. remap_1d does the remapping, i.e., it vertically advects
! ---         all prognostic variables (after re-inflating layers that may
! ---         have become too thin during smoothing)

   thcol(:) = worka(:,ipn)
   excol(:) = ex3d(:,ipn)
   prcol(:) = pr3d(:,ipn)

   call regrid_1d(its,targt,thcol,excol,prcol,vrbos,ipn,PrintDiags)

   exsmo3d(:,ipn) = excol(:)
  enddo						! loop over icos grid
!SMS$PARALLEL END

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (intfc_smooth.gt.0.) then
! --- make smoothing coefficient a function of layer index
   smoo_coeff(1:kbl+1)=0.		! don't modify lowest -kbl- layers
   do k=kbl+2,nvl
    smoo_coeff(k)=dt*intfc_smooth*taper(k,kbl+1)*taper(k,nvlp1)
   end do

   if (intfc_smooth.gt.0.) call dffusn_lev(exsmo3d,smoo_coeff,nvlp1,kbl+2,nvl)

! --- check occasionally whether smoothing is generating neg.lyr.thknss
   if (mod(its,100).eq.0) then
    do k=nvl/3,2*nvl/3
!SMS$PARALLEL(dh, ipn) BEGIN
     do ipn=1,nip
      exdif3d(k,ipn)=exsmo3d(k,ipn)-exsmo3d(k+1,ipn)
     end do
     valmin=minval(exdif3d(k,:))
!SMS$REDUCE(valmin,MIN)
!SMS$PARALLEL END
     if (valmin.lt.0.) then
      write (string,'(a,i2)') 'exdif aftr smoo k=',k
      call findmxmn2(exdif3d,nvl,nip,k,string)
     end if
    end do					! k loop
   end if
  end if					! intfc_smooth
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!SMS$PARALLEL(dh, ipn) BEGIN
  do ipn=1,nip					! loop over icos grid
   vrbos=ipn.eq.PrintIpnDiag

   trcol(:,:) = tr3d(:,ipn,:)
   excol(:)   = ex3d(:,ipn)
   dpcol(:)   = dp3d(:,ipn)
   prcol(:)   = pr3d(:,ipn)
   ucol (:)   = us3d(:,ipn)
   vcol (:)   = vs3d(:,ipn)
   exsmo(:)   = exsmo3d(:,ipn)

   call remap_1d(its,targt,ntra,trcol,ucol,vcol,excol,exsmo,		&
                dpcol,prcol,vrbos,ipn,PrintDiags)

   tr3d(:,ipn,:) = trcol(:,:)
   sdot(:,ipn)   = excol(:)-ex3d(:,ipn)
   ex3d(:,ipn)   = excol(:)
   dp3d(:,ipn)   = dpcol(:)
   pr3d(:,ipn)   = prcol(:)
   us3d(:,ipn)   = ucol (:)
   vs3d(:,ipn)   = vcol (:)

  end do					! loop over icos grid

!!sms$compare_var(ex3d, "hybgen.F90 - ex3d8 ")
!SMS$PARALLEL END

  return
  end subroutine hybgen


  subroutine regrid_1d(its,targt,thcol,excol,prcol,vrbos,ipn,PrintDiags)

  use module_control  ,only: nvl,nvlp1,kbl,nip,dt,thktop,pure_sig,ptop
  use module_constants,only: p1000,rd,cp,dpsig,grvity,deg_lat,deg_lon,	&
                             sigak,sigbk
  implicit none

  real,intent(IN)    :: prcol(nvlp1)		! must be consistent with excol
  real,intent(INOUT) :: excol(nvlp1)		! must be consistent with prcol
  real,intent(IN)    :: targt(nvl)		! target pot.temp.
  real,intent(INOUT) :: thcol(nvl)		! actual pot.temp.
  integer,intent(IN) :: its,ipn
  logical,intent(IN) :: vrbos,PrintDiags

! Type and dimension of local variables:

  integer :: k,k1,k2,iter		! layer/level indices
  real    :: thnew(nvl)			! new pot.temp.
  real    :: exnew(nvlp1)		! new Exner function
  real    :: exwrk(nvlp1)		! intermediate column values
  real    :: heatfx(nvl)		! turbulent vertical heat flux
  real    :: arg,wgt1,wgt2,try,coeff,cnv
  real    :: ex2p,dex2dp,p2ex,dp2dex
  real    :: eqlb_slak
  logical :: event
  real,parameter :: dffudt=.1		! therm.diffu.coeff x time step [m^2]
  real,parameter :: thin=rd/p1000	! approx. 1 Pa in Exner fcn units
  real,parameter :: tolrnce=0.001	! in degrees

  ex2p(arg)=p1000*(arg/cp)**(cp/rd)	!  convert Pi => p
  p2ex(arg)=cp*(arg/p1000)**(rd/cp)	!  convert p  => Pi

  eqlb_slak=.1*float(its)/(10.+float(its))

!SMS$IGNORE BEGIN
  if (vrbos) then
    write (6,99) its,ipn,deg_lat(ipn),deg_lon(ipn),'i n p u t  profile:'
    do k2=1,nvl,10
     write (6,100) (prcol(k1),k1=k2,min(nvlp1,k2+10) )
     write (6,102) (excol(k1),k1=k2,min(nvlp1,k2+10) )
     write (6,101) (thcol(k1),k1=k2,min(nvl  ,k2+9 ) )
     write (6,101) (targt(k1),k1=k2,min(nvl  ,k2+9 ) )
     write (6,101)
    end do
  end if
 99  format ('its,ipn=',i6,i8,'  lat/lon=',2f7.1,'  hybgen  ',a/	&
       '(4-line groups: pressure, exn.fcn, theta, target)')
 100 format (-2p,11f7.1)
 102 format (   11f7.1)
 101 format (5x,10f7.2)
!SMS$IGNORE END

  do k=1,nvl
   if (excol(k).lt.excol(k+1)-thin) then
!SMS$IGNORE BEGIN
    write (6,99) its,ipn,deg_lat(ipn),deg_lon(ipn),'i n p u t  profile:'
    do k2=1,nvl,10
     write (6,100) (prcol(k1),k1=k2,min(nvlp1,k2+10) )
     write (6,102) (excol(k1),k1=k2,min(nvlp1,k2+10) )
     write (6,101) (thcol(k1),k1=k2,min(nvl  ,k2+9 ) )
     write (6,101) (targt(k1),k1=k2,min(nvl  ,k2+9 ) )
     write (6,101)
    end do
    write (6,103) 'nonmonotonic Exner function on input at ipn,k =',	&
      ipn,k,excol(k),excol(k+1)
 103 format (a,i8,i4,2f9.2)
!SMS$IGNORE END
!   stop '(error: non-monotonic Exner fcn)'
    excol(k+1)=excol(k)
   end if
  end do

  if (pure_sig) then

!  if (vrbos) then
!   write (6,'(a/(5f14.6))') '(regrid_1d) sigak array:',sigak
!   write (6,'(a/(5f14.6))') '(regrid_1d) sigbk array:',sigbk
!  end if

! --- ><><><><><><><><><><><><><><><><><><>
! ---   restore to sigma coordinate grid
! --- ><><><><><><><><><><><><><><><><><><>

   exnew(    1)=excol(    1)
   exnew(nvlp1)=excol(nvlp1)
   do k=2,nvl
! --- sigak,sigbk define the sigma-p levels used in the GFS model.
    exnew(k)=p2ex(sigak(k)+sigbk(k)*prcol(1))
   end do

   if (vrbos) then
!SMS$IGNORE BEGIN
    write (6,98) its,ipn,'o u t p u t     p r e s s u r e'
 98 format ('time step',i6,'  ipn=',i8,'   hybgen    ',a)
    do k2=1,nvl,10
     write (6,100) (ex2p(exnew(k1)),k1=k2,min(nvlp1,k2+10) )
    end do
!SMS$IGNORE END
   end if
  else

! --- ><><><><><><><><><><><><><><><><><><>
! ---    restore to hybrid-isentropic grid
! --- ><><><><><><><><><><><><><><><><><><>

   ! --- eliminate static instabilities

   do k=nvl-1,1,-1
    thcol(k)=min(thcol(k),thcol(k+1))
   end do
   exwrk(:)=excol(:)

! --- 'thcol' profile is now stably (or neutrally) stratified.

! --- check whether theta values exceed 'targt' range at bottom of column.
! --- if so, homogenize layers to eliminate this condition

!?     do  k=1,nvl-1
!?      if (thcol(k).ge.targt(1)-tolrnce) then
!?       exit				!  no action required
!?      else
!? ! --- average of layers 1...k colder than lowest target
!? 
!? !SMS$IGNORE BEGIN
!?       if (vrbos) then
!?        write (*,107) ipn,'lyrs',1,k+1,' bfore range limiting:',	&
!?         (exwrk(k2),thcol(k2),k2=1,k+1),exwrk(k+2)
!?       end if
!? !SMS$IGNORE END
!? 
!? ! --- compute pot.temp. obtained by homogenizing layers 1,...,k+1
!?       wgt1=max(thin,exwrk(k  )-exwrk(k+1))
!?       wgt2=max(  0.,exwrk(k+1)-exwrk(k+2))
!?       try=(thcol(k)*wgt1+thcol(k+1)*wgt2)/(wgt1+wgt2)
!?       if (try.lt.targt(1)+tolrnce) then
!? ! --- average of layers 1,...,k+1 still too cold. continue adding layers
!?        exwrk(k+1)=exwrk(1)
!?        do k1=1,k+1
!?         thcol(k1)=try
!?        end do
!? 
!? !SMS$IGNORE BEGIN
!?        if (vrbos) then
!?         write (*,107) ipn,'lyrs',1,k+1,' after range limiting:',	&
!?          (exwrk(k2),thcol(k2),k2=1,k+1),exwrk(k+2)
!?        end if
!? !SMS$IGNORE END
!? 
!?       else
!? ! --- adding all of layer k+1 is overkill; entrain only part of lyr k+1
!?        exwrk(k+1)=min(exwrk(k  ),max(exwrk(k+2),			&
!?                      (exwrk(k  )*(thcol(k  )-targt(1))		&
!?                  +    exwrk(k+1)*(thcol(k+1)-thcol(k)))		&
!?                  /               (thcol(k+1)-targt(1))))
!?        do k1=1,k
!?         thcol(k1)=targt(1)
!?        end do
!? 
!? !SMS$IGNORE BEGIN
!?        if (vrbos) then
!?         write (*,107) ipn,'lyrs',1,k+1,' after range limiting:',	&
!?          (exwrk(k2),thcol(k2),k2=1,k+1),exwrk(k+2)
!?        end if
!? !SMS$IGNORE END
!? 
!?        exit				!  range limiting completed
!?       end if
!?      end if
!?     end do

! --- check whether theta values exceed 'targt' range at top of column.
! --- if so, homogenize layers to eliminate this condition

    do  k=nvl,2,-1
     if (thcol(k).le.targt(nvl)+tolrnce) then
      exit				!  no action required
     else
! --- average of layers k...nvl warmer than highest target

      if (vrbos) then
!SMS$IGNORE BEGIN
       write (*,107) ipn,'lyrs',k-1,nvl,' bfore range limiting:',	&
        (exwrk(k2),thcol(k2),k2=k-1,nvl),exwrk(nvlp1)
!SMS$IGNORE END
      end if

! --- compute pot.temp. obtained by homogenizing layers k-1,...,nvl
      wgt1=max(  0.,exwrk(k-1)-exwrk(k  ))
      wgt2=max(thin,exwrk(k  )-exwrk(k+1))
      try=(thcol(k-1)*wgt1+thcol(k)*wgt2)/(wgt1+wgt2)
      if (try.gt.targt(nvl)-tolrnce) then
!---average of layers k-1,...,nvl still too warm. continue adding layers
       exwrk(k)=exwrk(nvlp1)
       do k1=k-1,nvl
        thcol(k1)=try
       end do

       if (vrbos) then
!SMS$IGNORE BEGIN
        write (*,107) ipn,'lyrs',k-1,nvl,' after range limiting:',&
         (exwrk(k2),thcol(k2),k2=k-1,nvl),exwrk(nvlp1)
!SMS$IGNORE END
       end if

      else
! --- adding all of layer k-1 is overkill; entrain only part of lyr k-1
       exwrk(k)=min(exwrk(k-1),max(exwrk(k+1),				&
                   (exwrk(k  )*(thcol(k  )-thcol(k-1))			&
               +    exwrk(k+1)*(targt(nvl)-thcol(k  )))			&
               /               (targt(nvl)-thcol(k-1))))
       do k1=k,nvl
        thcol(k1)=targt(nvl)
       end do

       if (vrbos) then
!SMS$IGNORE BEGIN
        write (*,107) ipn,'lyrs',k-1,nvl,' after range limiting:',	&
         (exwrk(k2),thcol(k2),k2=k-1,nvl),exwrk(nvlp1)
!SMS$IGNORE END
       end if

       exit				!  range limiting completed
      end if
     end if
    end do
 107 format (i8,2x,a,i3,'-',i2,a,20(f7.1,f8.3))

! --- now convert column to purely isentropic coordinates, i.e.,
! --- find pressure levels where all pot.temps are on target

   call restp_1d(thcol,exwrk,nvl,thnew,exnew,targt,nvl,vrbos,ipn)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! --- invoke heat diffusion (McDougall-Dewar) to inflate thin layers

  heatfx(nvl)=0.
  heatfx(  1)=0.
  do iter=1,3				! apply the scheme repeatedly
   exwrk(:)=exnew(:)
   do k=2,nvl-1
    heatfx(k)=0.    
    if (exwrk(k  ).lt.exwrk(    1)-.01 .and.				&
        exwrk(k+1).gt.exwrk(nvlp1)+.01) then
     cnv=grvity/targt(k)				! Exner fcn units/meter
     heatfx(k)=dffudt*cnv*cnv*.5*(targt(k+1)-targt(k-1))		&
      /max(.03*cnv,exwrk(k)-exwrk(k+1))
    end if
   end do
   do k=1,nvl-1
    if (exwrk(k+1).lt.exwrk(    1)-.01 .and.				&
        exwrk(k+1).gt.exwrk(nvlp1)+.01)					&
     exnew(k+1)=max(exwrk(k+2),min(exwrk(k),				&
       exwrk(k+1)+(heatfx(k+1)-heatfx(k))/(targt(k+1)-targt(k))))
   end do
   event=.false.
   do k=2,nvlp1
    if (exnew(k).gt.exnew(k-1)+thin) then
     event=.true.
!SMS$IGNORE BEGIN
     print '(a,i8,i4,a,2F8.1)','dp<0 due to heat diffusion at ipn,k=',	&
      ipn,k,'  lat/lon =',deg_lat(ipn),deg_lon(ipn)
!SMS$IGNORE END
    end if
   end do				! iter

   if (vrbos .or. event) then
!SMS$IGNORE BEGIN
    print '(i8,a,i2,a)',ipn,' heat diffusion, iter',iter,		&
     '  Ex.fcn (3-line groups: old,new,dif x 10^4)'
    do k2=1,nvl,10
     write (6,108) (exwrk(k1),k1=k2,min(nvlp1,k2+9) )
     write (6,108) (exnew(k1),k1=k2,min(nvlp1,k2+9) )
     write (6,109) (int(1.e4*(exnew(k1)-exwrk(k1))),k1=k2,min(nvlp1,k2+9) )
     write (6,108)
    end do
 108 format (10f8.2)
 109 format (10i8)
!SMS$IGNORE END
   end if

   if (.not.event) exit
  end do			! iter
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

! --- suppress comput.mode causing gradual depletion of alternate layers

   call equilb(thnew,exnew,nvl,eqlb_slak,vrbos,ipn,PrintDiags)

! --- inflate massless layers

   call inflate(ipn,prcol(1),exnew)

  end if			!  pure sigma or hybrid-isentropic option

  if (vrbos) then
!SMS$IGNORE BEGIN
   write (6,104) ipn,'  hybgen: old Exner fcn (excol)',excol
   write (6,104) ipn,'  hybgen: new Exner fcn (exnew)',exnew
   write (6,104) ipn,'  hybgen: Exner fcn tndcy (displ)',exnew-excol
 104 format (i8,a/(8f9.2))
!SMS$IGNORE END
  end if

  excol(:)=exnew(:)
  return
  end subroutine regrid_1d
end module module_hybgen	! SMS doesn't like multiple routines in module


  subroutine remap_1d(its,targt,ntra,trcol,ucol,vcol,excol,exsmo,	&
                      dpcol,prcol,vrbos,ipn,PrintDiags)

  use module_control  ,only: nvl,nvlp1,kbl,nip,dt,thktop,pure_sig,	&
                             intfc_smooth,slak
  use module_constants,only: p1000,rd,cp,dpsig,deg_lat,deg_lon
  implicit none

  integer,intent(IN) :: ntra		! number of tracer fields
  real,intent(INOUT) :: trcol(nvl,ntra),ucol(nvl),vcol(nvl),		&
                        excol(nvlp1),dpcol(nvl),prcol(nvlp1)
  real,intent(IN)    :: targt(nvl)	! target pot.temp.
  real,intent(IN)    :: exsmo(nvlp1)	! Exner fcn after regridding & smoothing
  integer,intent(IN) :: its,ipn
  logical,intent(IN) :: vrbos,PrintDiags

! Type and dimension of local variables:

  integer :: k,k1,k2,n			! layer/level indices
  real    :: trnew(nvl,ntra)		! new tracers (1=pot.temp.)
  real    :: prnew(nvlp1)		! new pres & lyr thickness
  real    :: exnew(nvlp1)		! new Exner fcn
  real    :: unew(nvl),vnew(nvl)	! new velocities
  real    :: exwrk(nvlp1)		! intermediate column values
  real    :: thwrk(nvlp1)		! intermediate column values
  real    :: pk1col(nvlp1)		! vert.coord. used for theta remap
  real    :: pk1new(nvlp1)		! vert.coord. used for theta remap

  real    :: arg,colin,clout,qq
  real    :: ex2p
  real    :: dplo,dpup,devilo,deviup,tha,thb
  real    :: kappa(nvl),avgkap(nvl)     ! variables used in kappa diagno
  real    :: exdif(nvlp1),rmsdsp        ! variables used in kappa diagno
  real    :: exsav(nvlp1)               ! variables used in kappa diagno

  logical,parameter :: kappa_diag=.FALSE.	! vertical diffusivity diagno
  real   ,parameter :: small=1.e-6
  real   ,parameter :: acurcy=1.3e-6	! for 32-bit word length

  integer,parameter :: cnsv=1		! cnsv=1: conserve pot+intern.energy
! integer,parameter :: cnsv=2		! cnsv=2: conserve column height

  ex2p(arg)=p1000*(arg/cp)**(cp/rd)	!  convert Pi => p

  if (cnsv.eq.1) pk1col=excol*prcol	!  pk1 = p^{kappa+1)
  if (cnsv.eq.2) pk1col=excol		!  pk1 = p^k (Exner fcn)

  exnew(:)=exsmo(:)

  if (.not.pure_sig) then
   if (intfc_smooth.gt.0.) call inflate(ipn,prcol(1),exnew)
   do k=1,nvlp1
! --- (optional:) retard restoration to reduce overshooting
    if (its.gt.0) exnew(k)=slak*exnew(k)+(1.-slak)*excol(k)
   end do
  end if

! --- update pressure, layer thickness, pressure^(kappa+1)

  do k=nvlp1,1,-1
   prnew(k)=ex2p(exnew(k))
   if (k.le.nvl) dpcol(k)=prnew(k)-prnew(k+1)
  end do

  if (cnsv.eq.1) pk1new=exnew*prnew	!  pk1 = p^(kappa+1)
  if (cnsv.eq.2) pk1new=exnew		!  pk1 = p^k (Exner fcn)

  prnew (nvlp1)=prcol (nvlp1)		! safeguard against roundoff error
  pk1new(nvlp1)=pk1col(nvlp1)		! safeguard against roundoff error
  prnew (1)=prcol (1)			! safeguard against roundoff error
  pk1new(1)=pk1col(1)			! safeguard against roundoff error

! --- interface movement spawns vertical advection of dependent variables.
! --- we have 3 advection choices: PCM,PLM,PPM. use in any combination by
! --- selectively activating lines below. (one subr. call per variable!)

!!call pcmadv(prcol, ucol,nvl,prnew, unew,nvl,.false.,ipn)	!  u vel.
  call plmadv(prcol, ucol,nvl,prnew, unew,nvl,.false.,ipn)	!  u vel.
!!call ppmadv(prcol, ucol,    prnew, unew,nvl,.false.,ipn)	!  u vel.

!!call pcmadv(prcol, vcol,nvl,prnew, vnew,nvl,.false.,ipn)	!  v vel.
  call plmadv(prcol, vcol,nvl,prnew, vnew,nvl,.false.,ipn)	!  v vel.
!!call ppmadv(prcol, vcol,    prnew, vnew,nvl,.false.,ipn)	!  v vel.

  do n=2,ntra			! all tracers except pot.temp
!! call pcmadv(prcol,trcol(1,n),nvl,prnew,trnew(1,n),nvl,.false.,ipn)
   call plmadv(prcol,trcol(1,n),nvl,prnew,trnew(1,n),nvl,.false.,ipn)
!! call ppmadv(prcol,qvcol(1,n),    prnew,trnew(1,n),nvl,.false.,ipn)
  end do

! --- now advect pot.temp (k=1)
!!call pcmadv(pk1col,trcol,nvl,pk1new,trnew,nvl,vrbos,ipn)	!  theta
  call plmadv(pk1col,trcol,nvl,pk1new,trnew,nvl,vrbos,ipn)	!  theta
!!call ppmadv(pk1col,trcol,    pk1new,trnew,nvl,vrbos,ipn)	!  theta

  if (.not.pure_sig) then
! --- redistribute theta among neighboring layers to help them stay on target.
! --- this is to counteract a comput.mode associated with vertical advection

  do k=nvl,3,-1
   dplo=max(pk1new(k-1)-pk1new(k  ),pk1new(1)*small)
   dpup=max(pk1new(k  )-pk1new(k+1),pk1new(1)*small)
   tha=trnew(k-1,1)
   thb=trnew(k  ,1)
   devilo=(tha-targt(k-1))*dplo
   deviup=(thb-targt(k  ))*dpup
   if (deviup.gt.0. .and. devilo.lt.0.) then
    trnew(k-1,1)=tha+min(deviup,-devilo)/dplo
    trnew(k  ,1)=thb-min(deviup,-devilo)/dpup
   else									&
   if (deviup.lt.0. .and. devilo.gt.0.) then
    trnew(k-1,1)=tha-min(devilo,-deviup)/dplo
    trnew(k  ,1)=thb+min(devilo,-deviup)/dpup
   end if

   if (vrbos .and. deviup*devilo.lt.0.) then   
!SMS$IGNORE BEGIN
    write (6,'(a,i8,i4,2(3x,a,2f9.4))') 'ipn,k =',ipn,k,		&
     'targ dev''n',tha-targt(k-1),thb-targt(k),'cut to',		&
      trnew(k-1,1)-targt(k-1),trnew(k,1)-targt(k)
!SMS$IGNORE END
   end if
  end do

! --- fill massless cells with data from mass-containing layer below

  if (thktop.eq.0.) then
   do k=3,nvl
    qq=1./(dpcol(k)+small)
    ucol(k)=(ucol(k)*dpcol(k)+ucol(k-1)*small)*qq
    vcol(k)=(vcol(k)*dpcol(k)+vcol(k-1)*small)*qq
    trnew(k,2:ntra)=(trnew(k,2:ntra)*dpcol(k)+trnew(k-1,2:ntra)*small)*qq
   end do
  end if

  end if			! hybrid-isentropic option

! --- column integrals (colin/clout) are for diagnostic purposes only
  colin=0.
  clout=0.
  do k=1,nvl
   colin=colin+trcol(k,1)*(pk1col(k)-pk1col(k+1))
   clout=clout+trnew(k,1)*(pk1new(k)-pk1new(k+1))
  end do

  if (abs(clout-colin).gt.acurcy*abs(colin)) then
!SMS$IGNORE BEGIN
    write (6,106) ipn,'hybgen - column intgl.error',		&
     colin,clout,(clout-colin)/colin
 106 format (i8,3x,a,2es14.6,es9.1)
  end if
!SMS$IGNORE END

  !--------------------------------------------------------------------
  ! --- vertical diffusivity diagnostics (optional)

  rmsdsp=0.
  do k=1,nvl
    rmsdsp=rmsdsp+(exnew(k+1)-excol(k+1))**2
  end do
  rmsdsp=sqrt(rmsdsp/float(nvl))

  if (kappa_diag .and. rmsdsp.gt..01) then
    thwrk(:)=trnew(:,1)
    exwrk(:)=exnew(:)
    call restp_1d(thwrk,exwrk,nvl,trnew,exnew,targt,nvl,vrbos,ipn)
    exdif(    1)=0.
    exdif(nvlp1)=0.
    do k=2,nvl
      exdif(k)=(exnew(k)-exsav(k))/dt
    end do

    call diagkp(nvl,exsav,exdif,targt,kappa,.false.,ipn)
!!!    call diagkp(nvl,exsav,exdif,targt,kappa,.false.,ipn)
    if (minval(kappa).lt.-1.e-2) then
!!!    if (maxval(kappa).gt. 1.e-2) then
     call diagkp(nvl,exsav,exdif,targt,kappa,.true.,ipn)

!SMS$IGNORE BEGIN
     write (6,'(a,2i5,5x,a/12(f6.1))') 'its,ipn =',its,ipn,		&
      'vert.diffusivity (cm^2/s):',(1.e4*kappa(k),k=2,nvl-1)
!!!     print *,'maxval:',maxval(kappa)
     print *,'minval:',minval(kappa)
!SMS$IGNORE END

    end if
    avgkap(:)=avgkap(:)+kappa(:)
  end if			!  diffusivity diagnostics

!SMS$IGNORE BEGIN
  if (kappa_diag .and. mod(its,15).eq.0)				&
   write (6,'(a,i5,5x,a/12(f6.1))') 'its =',its,			&
    'avg.vert.diffusivity (cm^2/s):',(1.e4*avgkap(k)/nip,k=2,nvl-1)
!SMS$IGNORE END
  !--------------------------------------------------------------------

  ucol(:)=unew(:)
  vcol(:)=vnew(:)
  trcol(:,:)=trnew(:,:)
  excol(:)=exnew(:)
  prcol(:)=prnew(:)

  do k=1,nvl
   if (excol(k).lt.excol(k+1)) then
!SMS$IGNORE BEGIN
    write (6,99) its,ipn,deg_lat(ipn),deg_lon(ipn),'o u t p u t  profile:'
    do k2=1,nvl,10
     write (6,100) (prcol(k1)  ,k1=k2,min(nvlp1,k2+10) )
     write (6,102) (excol(k1)  ,k1=k2,min(nvlp1,k2+10) )
     write (6,101) (trcol(k1,1),k1=k2,min(nvl  ,k2+9 ) )
     write (6,101) (targt(k1)  ,k1=k2,min(nvl  ,k2+9 ) )
     write (6,101)
    end do
    write (6,103) 'nonmonotonic Exner function on return at ipn,k =',	&
     ipn,k,excol(k),excol(k+1)
 103 format (a,i8,i4,2f9.2)
!SMS$IGNORE END
!   stop '(error: non-monotonic Exner fcn)'
    excol(k+1)=excol(k)
   end if
  end do

  if (vrbos) then
!SMS$IGNORE BEGIN
    write (6,99) its,ipn,deg_lat(ipn),deg_lon(ipn),'o u t p u t  profile:'
   do k2=1,nvl,10
    write (6,100) (prcol(k1)  ,k1=k2,min(nvlp1,k2+10) )
    write (6,102) (excol(k1)  ,k1=k2,min(nvlp1,k2+10) )
    write (6,101) (trcol(k1,1),k1=k2,min(nvl  ,k2+9 ) )
    write (6,101) (targt(k1)  ,k1=k2,min(nvl  ,k2+9 ) )
    write (6,101)
   end do
 99  format ('its,ipn=',i6,i8,'  lat/lon=',2f7.1,'  hybgen  ',a/	&
       '(4-line groups: pressure, exn.fcn, theta, target)')
 100 format (-2p,11f7.1)
 102 format (   11f7.1)
 101 format (5x,10f7.2)
!SMS$IGNORE END
  end if

  return
  end subroutine remap_1d


  subroutine restp_1d(thold,pkold,kold,thnew,pknew,targt,knew,vrbos,ipn)
  
! --- convert a stairstep (i.e., piecewise constant) theta profile into a
! --- stairstep profile constrained to have prescribed theta ('targt') steps.
  
! --- input  variables: thold,pkold,targt,kold,knew,vrbos
! --- output variables: thnew,pknew
  
  use module_constants, only: cp,rd,p1000
  
  implicit none
  integer,intent(IN) :: kold,knew,ipn
  real,intent(IN)    :: thold(kold),pkold(kold+1),targt(knew)
  real,intent(OUT)   :: thnew(knew),pknew(knew+1)
  logical,intent(IN) :: vrbos

  integer k,ko
  real oldth(kold)
  real cloutt,colint,pinteg,tha,thb,ex2p,arg
  real, parameter :: acurcy=1.e-6
  
  ex2p(arg)=p1000*(arg/cp)**(cp/rd)      !  convert Pi => p (mb)
  
  if (vrbos) then
!SMS$IGNORE BEGIN
    write (6,101) ipn,							&
    'restp1 -- input profile:   theta    thknss    press     p^kap',	&
    ex2p(pkold(1)),pkold(1),						&
    (k,thold(k),ex2p(pkold(k))-ex2p(pkold(k+1)),ex2p(pkold(k+1)),	&
    pkold(k+1),k=1,kold)
 101 format (i8,4x,a/54x,f9.1,f11.3/(33x,i3,f9.3,2f9.1,f11.3))
!SMS$IGNORE END
  end if
  
! --- remove theta inversions from input profile
  oldth(kold)=thold(kold)
  do k=kold,2,-1
    oldth(k-1)=min(oldth(k),thold(k-1))
  end do
  
  thnew(:)=targt(:)
  thnew(   1)=min(oldth(1),oldth(kold),targt(   1))
  thnew(knew)=max(oldth(1),oldth(kold),targt(knew))
  pknew(     1)=pkold(     1)
  pknew(knew+1)=pkold(kold+1)
  
! --- column integrals (colin/clout) are computed for diagnostic purposes only
  
  cloutt=0.
  colint=0.
  do k=1,kold
    colint=colint+oldth(k)*(pkold(k+1)-pkold(k))
  end do
  
! --- find interface pknew(k+1) separating layers k and k+1 by requiring
! --- that integral over pk*d(theta) from thnew(k) to thnew(k+1) be preserved.
  
  ko=1
  do k=1,knew-1
    pinteg=0.
    thb=thnew(k)
  5 tha=thb
    thb=min(thnew(k+1),max(thnew(k),oldth(ko)))
    pinteg=pinteg+pkold(ko)*(thb-tha)
    if (oldth(ko) < thnew(k+1)) then
      if (ko.lt.kold) then
        ko=ko+1
        go to 5
      end if
      tha=thb
      thb=thnew(k+1)
      pinteg=pinteg+pkold(kold+1)*(thb-tha)
    end if
    pknew(k+1)=pknew(k)
    if (thnew(k+1) > thnew(k)) pknew(k+1)=pinteg/(thnew(k+1)-thnew(k))
    cloutt=cloutt+thnew(k)*(pknew(k+1)-pknew(k))
  enddo
  
  cloutt=cloutt+thnew(knew)*(pknew(knew+1)-pknew(knew))
  if (abs(cloutt-colint).gt.acurcy*abs(colint)) then
!SMS$IGNORE BEGIN
    write (6,100) ipn,'restp1 - column intgl.error',			&
    colint,cloutt,(cloutt-colint)/colint
 100 format (i8,3x,a,2es14.6,es9.1)
!SMS$IGNORE END
  end if
  
  if (vrbos) then
!SMS$IGNORE BEGIN
    write (6,101) ipn,					     		&
    'restp1 -- outpt profile:   theta    thknss    press     p^kap',	&
    ex2p(pknew(1)),pknew(1),						&
    (k,thnew(k),ex2p(pknew(k))-ex2p(pknew(k+1)),ex2p(pknew(k+1)),	&
    pknew(k+1),k=1,knew)
!SMS$IGNORE END
  end if
  
  return
  end subroutine restp_1d


  subroutine inflate(ipn,psurf,exner)

! --- <><><><><><<><><><><><><><><><><><><><><><><>
! --- hybridization (inflation of massless layers)
! --- <><><><><><<><><><><><><><><><><><><><><><><>

  use module_constants, only: cp,rd,p1000,dpsig
  use module_control  ,only: nvl,nvlp1,kbl,PrintIpnDiag,thktop,ptop
 
  implicit none
  integer, intent(IN)       :: ipn		! icos index
  real   , intent(IN)       :: psurf		! surface pressure
  real   , intent(INOUT)    :: exner(nvlp1)	! Exner function

  real    :: dp0,dp1,dpsum,arg,arg1,exold,q,extop
  real    :: dex2dp,dp2dex,shrink,top,cushn,taper
  integer :: last,k,k1,k2
  logical :: vrbos			!  switch for 'verbose' mode
  real,parameter :: thin=rd/p1000	! approx. 1 Pa in Exner fcn units

! taper(k1,k2)=1./float(k2-k1)
  taper(k1,k2)=1./(1.+(.25*(k2-k1))**2)
! taper(k1,k2)=1./(1.+(.5*(k2-k1))**2)
! dex2dp(arg)=p1000*arg/rd		!  convert d Pi => dp (for p ~ p1000)
  dp2dex(arg)=rd*arg/p1000		!  convert dp => d Pi (for p ~ p1000)

! --- simplified cushion function suitable for nonnegative first argument
  cushn(arg,arg1)=.25*min(arg,2.*arg1)**2/arg1+max(arg,2.*arg1)-arg1
! cushn(arg,arg1)=arg1			! no cushn

  vrbos=ipn.eq.PrintIpnDiag

  dpsum=0.
  extop=cp*(ptop/p1000)**(rd/cp)
  top=.2		! sigma layers are shrunk in proportion to (psurf-top) 
  shrink=min(1.,(psurf/p1000-top)/(1.-top))
  last=1

! --- step 1: inflate layers from bottom up

  do k=1,nvl-1

! --- set lower limits for layer thknss (dp0)
! --- and upper limits for upper intfc pressure (dpsum)

   dp0=dp2dex(dpsig(k))*shrink	!  convert to Exner fcn, shrink abv mountains
   dpsum=dpsum+dp0		!  minimum cumulative distance to ground
   exold=exner(k+1)
!  if (exner(k)-dp0.lt.extop) then
!!SMS$IGNORE BEGIN
!   print '(2a,i8,i4,a,f9.2)','warning: stack of min.thknss layers',	&
!    ' extends past top of atmosphere at ipn,k=',ipn,k,'  extop=',extop
!!SMS$IGNORE END
!  end if

   if (k.le.kbl) then
! --- maintain -kbl- fixed-depth layers near surface
    dp1=dp0
    last=k+1
    exner(k+1)=exner(k)-dp1

   else				!  k > kbl
    if (k.gt.last) then		!  out of sigma domain -> reduce dp0
     dp0=dp0/min(20.,1.25*2.**(k-last))
!!   dp0=dp0/min(100.,12.5*2.**(k-last))
!!   dp0=dp0/min(100.,12.5*float(k-last)**2)
!!   dp0=dp0/min(10. ,1.25*float(k-last)**2)
! next line sets 1 hPa as minimum thickness in isentropically resolved layers
!    26 Oct 2010 - Rainer, Stan
     dp0 = max(dp0,dp2dex(100.))
    end if
    dp1=cushn(max(0.,exner(k)-min(exold,exner(1)-dpsum)),dp0)
    if (exner(k)-dp1.lt.exold-thin) then	!  inflation required
     if (k.eq.last) last=k+1			!  still in sigma domain
     exner(k+1)=max(extop,exner(k)-dp1)
    end if
   end if			!  k < or > kbl

   if (vrbos .and. exold.ne.exner(k+1)) then
!SMS$IGNORE BEGIN
    write (6,105) ipn,k+1,exold,exner(k+1),			&
     max(0.,exner(k)-min(exold,exner(1)-dpsum)),dp0,dp1
 105 format (i8,'  k=',i3,'  exn.',f8.2,'  =>',f8.2,            &
      '  arg1/2,cush =',3f8.2)
!SMS$IGNORE END
   end if

  end do			!  k loop

! --- step 2: inflate layers from top down

  if (thktop.gt.0.) then
   do k=nvl,2,-1
    q=exner(k)/cp
    q=q*q*sqrt(q)		!  (p/p0)**(5/7) = (Exn/cp)**(5/2)
    dp0=dp2dex(thktop)*taper(k,nvl)/q
    exold=exner(k)
    dp1=cushn(max(0.,exold-exner(k+1)),dp0)
    if (exner(k+1)+dp1.gt.exold+thin) then
     exner(k)=exner(k+1)+dp1

     if (vrbos) then
!SMS$IGNORE BEGIN
      write (6,105) ipn,k,exold,exner(k),max(0.,exold-exner(k+1)),dp0,dp1
!SMS$IGNORE END
     end if
    
    end if
    if (exner(k).lt.exner(k-1)-dp0) exit
   end do			!  k loop
  end if			!  thktop > 0

  return
  end subroutine inflate


    subroutine ppmadv(xold,fldold,xnew,fldnew,kk,vrbos,ipn)

!-- PPM-based 1-dim transport routine, extracted from HYCOM's fct3d.f.
!-- Note: |xold-xnew| can exceed cell width (i.e., no CFL constraints)

!-- xold/new	- old/new cell boundaries
!-- fldold/new	- mixing ratio of dep.variable before and after transport
!-- kk		- number of layers
!-- vrbos	- if .true., print diagnostic messages for grid point -ipn-

    implicit none
    integer, intent(IN) :: kk,ipn
    real,    intent(IN) :: xold(kk+1),xnew(kk+1),fldold(kk)
    real,   intent(OUT) :: fldnew(kk)
    logical, intent(IN) :: vrbos	!  switch for 'verbose' mode

    real zold(kk+1),znew(kk+1),delx(kk+1),delz(kk+1),fco(kk),fcn(kk),	&
         vertfx(kk+1),vertdv(kk)
    real a(kk),b(kk),c(kk),dx,fcdx,yl,yr
    real amount,bfore,after,dpth,scale,slab,dslab
    integer k,lyr
    real, parameter :: athird=1./3.
    real, parameter :: small=1.e-9
    real, parameter :: acurcy=1.e-5

    delx=xnew-xold
!-- make sure -x- increases with -k-. change sign of -xold/new- if necessary
    if (xold(1).lt.xold(kk+1)) then
      zold=xold
      znew=xnew
    else
      zold=-xold
      znew=-xnew
    end if
    delz=znew-zold

!SMS$IGNORE BEGIN
    if (vrbos)								&
     write (*,100) ipn,'entering ppmadv: old_p^kap d(p^kap)  variable',	&
      (k,xold(k),delx(k),fldold(k),k=1,kk),kk+1,xold(kk+1),delx(kk+1)
 100 format (i8,3x,a/(23x,i3,2f10.3,es11.3))
!SMS$IGNORE END

!-- deduce old and new cell width from -zold,znew-
    do 15 k=1,kk
    fco(k)=max(0.,zold(k+1)-zold(k))
15  fcn(k)=max(0.,znew(k+1)-znew(k))

    bfore=0.
    scale=0.
    dpth=0.
    do k=1,kk
      bfore=bfore+fldold(k)*fco(k)
      dpth=dpth+fco(k)
      scale=scale+abs(fldold(k))
    end do
    fldnew=fldold

!-- start by filling zero-width cells with data from neighboring cells

    do 17 k=kk-1,1,-1
17  fldnew(k)=(fldnew(k)*fco(k)+fldnew(k+1)*small)			&
             /(          fco(k)+            small)
    do 18 k=2,kk
18  fldnew(k)=(fldnew(k)*fco(k)+fldnew(k-1)*small)			&
             /(          fco(k)+            small)

!-- fit 0th, 1st, or 2nd deg. polynomial to -fldnew- in each cell
    a(1 )=fldnew(1 )
    b(1 )=0.
    c(1 )=0.
    a(kk)=fldnew(kk)
    b(kk)=0.
    c(kk)=0.

    do 16 k=2,kk-1
!-- uncomment one of the following 3 options to activate pcm,plm,ppm resp.
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!-- piecewise constant method:
!cc    a(k)=fldnew(k)
!cc    b(k)=0.
!cc    c(k)=0.
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!-- piecewise linear method:
!-- fit linear function a+bx to tracer in each cell (-.5 < x < +.5)
!cc    a(k)=fldnew(k)
!cc    b(k)=0.
!cc    if (fldnew(k).le.min(fldnew(k-1),fldnew(k+1)) .or.		&
!cc        fldnew(k).ge.max(fldnew(k-1),fldnew(k+1))) then
!cc      b(k)=0.
!cc    else if ((fldnew(k+1)-fldnew(k-1))*(fldnew(k-1)+fldnew(k+1)	&
!cc      -2.*fldnew(k)).gt.0.) then
!cc      b(k)=fldnew(k)-fldnew(k-1)
!cc    else
!cc      b(k)=fldnew(k+1)-fldnew(k)
!cc    end if
!cc    c(k)=0.
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!-- piecewise parabolic method:
!-- construct parabola  a+bx+cx^2  whose integral over [-.5,+.5] equals
!-- fldnew(k) and which passes though points yl,yr at [-.5,+.5] resp.
!!    yl=.5*(fldnew(k-1)+fldnew(k))
!!    yr=.5*(fldnew(k+1)+fldnew(k))
    yl=(max(small,fco(k-1))*fldnew(k)+max(small,fco(k))*fldnew(k-1))/	&
       (max(small,fco(k-1))          +max(small,fco(k)))
    yr=(max(small,fco(k+1))*fldnew(k)+max(small,fco(k))*fldnew(k+1))/	&
       (max(small,fco(k+1))          +max(small,fco(k)))
    a(k)=1.5*fldnew(k)-.25*(yl+yr)
    b(k)=yr-yl
    c(k)=6.*(.5*(yl+yr)-fldnew(k))
    if (abs(yr-yl) .lt. 6.*abs(.5*(yl+yr)-fldnew(k))) then
!-- apex of parabola lies inside interval [-.5,+.5], implying an over-
!-- or undershoot situation. change curve to prevent over/undershoots.
      if (abs(yr-yl) .gt. 2.*abs(.5*(yl+yr)-fldnew(k))) then
!-- put apex of parabola on edge of interval [-.5,+.5]
        if ((yr-yl)*(.5*(yl+yr)-fldnew(k)) .gt. 0.) then
!-- apex at x=-.5
          a(k)=.25*(3.*fldnew(k)+yl)
          c(k)=3.*(fldnew(k)-yl)
          b(k)=c(k)
        else
!-- apex at x=+.5
          a(k)=.25*(3.*fldnew(k)+yr)
          c(k)=3.*(fldnew(k)-yr)
          b(k)=-c(k)
        end if
      else			!  -1/6 < x < +1/6
!-- moving apex won't help. replace parabola by constant.
        a(k)=fldnew(k)
        b(k)=0.
        c(k)=0.
      end if
    end if
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
16  continue

!-- get flux by summing -fldnew- over upstream slab of thickness -delz-

    do 22 k=2,kk
    slab=0.
    amount=0.
    vertfx(k)=0.
    if (delz(k).gt.0.) then			! interface moves in +k dir.
      lyr=k-1
24    lyr=lyr+1
      if (slab.ge.delz(k)) goto 23
      if (fco(lyr).gt.0.) then
        dslab=min(slab+fco(lyr), delz(k))	&
             -min(slab         , delz(k))
        dx=dslab/fco(lyr)
        fcdx=a(lyr)				&
            +b(lyr)*.5*(dx-1.)			& !  not needed in pcm
            +c(lyr)*(.25-dx*(.5-dx*athird))	  !  not needed in pcm,plm
        amount=amount+fcdx*dslab
        slab=slab+dslab
      end if
      if (lyr.lt.kk) go to 24
    else if (delz(k).lt.0.) then		! interface moves in -k dir.
      lyr=k
25    lyr=lyr-1
      if (slab.ge.-delz(k)) goto 23
      if (fco(lyr).gt.0.) then
        dslab=min(slab+fco(lyr),-delz(k))	&
             -min(slab         ,-delz(k))
        dx=dslab/fco(lyr)
        fcdx=a(lyr)				&
            +b(lyr)*.5*(1.-dx)			& !  not needed in pcm
            +c(lyr)*(.25-dx*(.5-dx*athird))	  !  not needed in pcm,plm
        amount=amount+fcdx*dslab
        slab=slab+dslab
      end if
      if (lyr.gt.2) go to 25
    end if
23  if (slab.ne.0.) vertfx(k)=-delz(k)*amount/slab
22  continue

    vertfx(   1)=0.			!  don't allow flux through lower bdry
    vertfx(kk+1)=0.			!  don't allow flux through upper bdry
    do 26 k=1,kk
26  vertdv(k)=vertfx(k+1)-vertfx(k)

!SMS$IGNORE BEGIN
    if (vrbos) write (*,'(a/(i3,4es12.3))')				&
     'ppmadv:   flux  flx.div/thk    old_thk     new_thk',		&
      (k,vertfx(k),vertdv(k)/max(small,fcn(k)),fco(k),fcn(k),k=1,kk),	&
       kk+1,vertfx(kk+1)
!SMS$IGNORE END

    do 4 k=1,kk
    amount=fldnew(k)*fco(k)-vertdv(k)
4   fldnew(k)=(fldnew(k)*small+amount)/(small+fcn(k))

    after=0.
    do k=1,kk
      after=after+fldnew(k)*fcn(k)
    end do

!SMS$IGNORE BEGIN
    if (abs(bfore-after)*kk.gt.acurcy*scale*dpth) then
      write (*,104) ipn,'ppmadv - bad column intgl.:',bfore,after
    end if
 104 format (i8,3x,a,2es15.7)

    if (vrbos)								&
     write (*,100) ipn,'exiting ppmadv:  d(p^kap) new_p^kap  variable',	&
      (k,delx(k),xnew(k),fldnew(k),k=1,kk),kk+1,delx(kk+1),xnew(kk+1)
!SMS$IGNORE END

    return
    end subroutine ppmadv


    subroutine plmadv(xold_r4,yold_r4,kold,xnew_r4,ynew_r4,knew,vrbos,ipn)
!
! --- this version performs all calculations in double precision (real*8).
! --- data are passed in and out as 'real'.
!
! --- consider two stepwise constant functions -yold,ynew- whose
! --- discontinuities are at abscissa values -xold,xnew- respectively.
! --- treat -ynew- as unknown. solve for -ynew- under the condition that
! --- the integral over y*dx is preserved (integration based on PLM).
!
    implicit none
    integer,intent(IN) :: kold,knew
    integer,intent(IN) :: ipn		!  current location in horiz.grid
    real,intent(IN)    :: yold_r4(kold),xold_r4(kold+1),xnew_r4(knew+1)
    real,intent(OUT)   :: ynew_r4(knew)
    logical,intent(IN) :: vrbos		!  if true, print diagnostics
    integer k,ko,n,kstart
    real*8 yold(kold),xold(kold+1),xnew(knew+1),ynew(knew),		&
         zold(kold+1),znew(knew+1),colin,clout,slope,wgta,wgtb,wgtc,	&
         yinteg,ylft(kold),yrgt(kold),xlo,xhi,ra,rb,ya,yb,q,scale,	&
         yrka,ylk,yrk,ylkb
    real*8 plmslp
    external plmslp
    logical at_top
    real,parameter    :: onemu=1.e-6, acurcy=1.e-11, flag=-999.
    integer,parameter :: iter=1
!
    xold=xold_r4
    yold=yold_r4
    xnew=xnew_r4

!SMS$IGNORE BEGIN
    if (vrbos)								&
     write (*,101) ipn,'old',(k,xold(k),yold(k),k=1,kold),		&
      kold+1,xold(kold+1)
 101 format (i8,'  plmadv -- ',a,' profile:',5x,'x',11x,'y'/		&
        (i27,f14.1,es12.4))
!SMS$IGNORE END
 
    if (xold_r4(     1).ne.xnew_r4(     1)) write (*,102) ipn,		&
     ' plmadv warning: bottom xold,xnew differ:',			&
      xold_r4(     1),xnew_r4(     1)
    if (xold_r4(kold+1).ne.xnew_r4(knew+1)) write (*,102) ipn,		&
     ' plmadv warning: top xold,xnew differ:',				&
      xold_r4(kold+1),xnew_r4(knew+1)
 102 format (i8,a,2f10.2)
    xold(1)=max(xold(1),xnew(1))
    xnew(1)=max(xold(1),xnew(1))
 
!-- make sure -x- increases with -k-. change sign of -xold/new- if necessary
    if (xold(1).lt.xold(kold+1)) then
      zold=xold
      znew=xnew
    else
      zold=-xold
      znew=-xnew
    end if
!
! --- column integrals (colin/clout) are computed for diagnostic purposes only
    scale=0.
    colin=0.
    clout=0.
    do 3 k=1,kold
    scale=scale+abs(yold(k))
 3  colin=colin+yold(k)*(zold(k+1)-zold(k))
!
! --- replace each flat segment of stairstep curve by
! --- a slanting segment, using PLM-type limiters.
!
    ylft(   1)=yold(   1)
    yrgt(   1)=yold(   1)
    ylft(kold)=yold(kold)
    yrgt(kold)=yold(kold)
    do 6 n=1,iter			!  iterate to optimize limiters
    do 2 k=2,kold-1
    if (n.eq.1) then
      yrka=yold(k-1)
      if (zold(k  ).eq.zold(     1)) yrka=yold(k)
      ylk=yold(k)
      yrk=yold(k)
      ylkb=yold(k+1)
      if (zold(k+1).eq.zold(kold+1)) ylkb=yold(k)
    else
      yrka=yrgt(k-1)
      if (zold(k  ).eq.zold(     1)) yrka=yold(k)
      ylk=ylft(k)
      yrk=yrgt(k)
      ylkb=ylft(k+1)
      if (zold(k+1).eq.zold(kold+1)) ylkb=yold(k)
    end if
    wgta=max(DBLE(onemu),zold(k  )-zold(k-1))
    wgtb=max(DBLE(onemu),zold(k+1)-zold(k  ))
    wgtc=max(DBLE(onemu),zold(k+2)-zold(k+1))
    if (k.eq.     1) wgta=onemu
    if (k.eq.kold-1) wgtc=onemu
    slope=plmslp((wgtb*yrka+wgta*ylk)/(wgtb+wgta),			&
         yold(k),(wgtb*ylkb+wgtc*yrk)/(wgtb+wgtc))
    ylft(k)=yold(k)-slope
 2  yrgt(k)=yold(k)+slope

!!SMS$IGNORE BEGIN
!   if (vrbos) print '(a,i2,5x,a,14x,a,8x,a/(i3,es12.4,5x,2es12.4))',	&
!    'iter',n,'y','ylft','yrgt',(ko,yold(ko),ylft(ko),yrgt(ko),		&
!     ko=1,kold)
!!SMS$IGNORE END
 6  continue
!
! --- y in k-th interval now varies from ylft at zold(k) to yrgt at zold(k+1).
! --- find ynew(k) by requiring 
! --- that the integral over y*dx from znew(k) to znew(k+1) be preserved.
!
    at_top=.true.
    kstart=1
    do 4 k=1,knew
    yinteg=0.
    xlo=znew(k  )
    xhi=znew(k+1)
!cc    if (vrbos) print '(a,2f9.3)','xlo,xhi =',xlo,xhi
    if (xhi.gt.xlo) then
      at_top=.false.
      do 5 ko=kstart,kold
      if (zold(ko+1).le.xlo) then
        kstart=ko+1
        go to 5
      end if
      if (zold(ko  ).ge.xhi) go to 1
! --- integrate over sloping portions of y(x) curve:
      ra=max(xlo,min(xhi,zold(ko  )))
      rb=max(xlo,min(xhi,zold(ko+1)))
      ya=ylft(k)
      yb=yrgt(k)
      wgta=flag
      wgtb=flag
      if (zold(ko+1).ne.zold(ko)) then
        if (ra.ge.zold(ko).and.ra.le.zold(ko+1)) then
          wgta=(zold(ko+1)-ra)/(zold(ko+1)-zold(ko))
          ya=ylft(ko)*wgta+yrgt(ko)*(1.-wgta)
        end if
        if (rb.ge.zold(ko).and.rb.le.zold(ko+1)) then
          wgtb=(zold(ko+1)-rb)/(zold(ko+1)-zold(ko))
          yb=ylft(ko)*wgtb+yrgt(ko)*(1.-wgtb)
        end if
      end if
      yinteg=yinteg+.5*(ya+yb)*(rb-ra)
!cc      if (vrbos) print '(2i4,4f9.3,3f11.1)',				&
!cc        k,ko,ra,rb,wgta,wgtb,ya,yb,yinteg
 5    continue
      yinteg=yinteg+yb*(xhi-rb)
!cc      if (vrbos) print '(2i4,4f9.3,3f11.1)',				&
!cc        k,0,rb,xhi,wgta,wgtb,yb,yb,yinteg
 1    ynew(k)=yinteg/(xhi-xlo)
    else if (at_top) then
      ynew(k)=yold(   1)
    else                              !  at end
      ynew(k)=yold(kold)
    end if
!cc    if (vrbos) print '(a,f11.1)','ynew =',ynew(k)
    clout=clout+ynew(k)*(znew(k+1)-znew(k))
 4  continue
!
!SMS$IGNORE BEGIN
    if (abs(clout-colin).gt.acurcy*scale*(zold(kold+1)-zold(1)))	&
      write (*,100) ipn,' plmadv - column intgl.error',			&
       colin,clout,(clout-colin)/colin
 100 format (i8,a,2es14.6,es9.1)
    if (vrbos)								&
     write (*,101) ipn,'new',(k,xnew(k),ynew(k),k=1,knew),		&
      knew+1,xnew(knew+1)
!SMS$IGNORE END
    ynew_r4=ynew
    return
    end subroutine plmadv


    real*8 function plmslp(ylft,ymid,yrgt)
!
! --- get slope at point 'ymid' for piecewise linear interpolation
    implicit none
    real*8,intent(IN) :: ylft,ymid,yrgt
!
    if (ymid.le.min(ylft,yrgt) .or.					&
        ymid.ge.max(ylft,yrgt)) then
      plmslp=0.
    else if ((yrgt-ylft)*(ylft+yrgt-2.*ymid).gt.0.) then
      plmslp=ymid-ylft
    else
      plmslp=yrgt-ymid
    end if
    return
    end function plmslp


    subroutine pcmadv(xold_r4,yold_r4,kold,xnew_r4,ynew_r4,knew,vrbos,ipn)
!
! --- this version performs all calculations in double precision (real*8).
! --- data are passed in and out as 'real'.
!
! --- consider two stepwise constant functions -yold,ynew- whose
! --- discontinuities are at abscissa values -xold,xnew- respectively.
! --- treat -ynew- as unknown. solve for -ynew- under the condition that
! --- the integral over y*dx is preserved (integration based on PCM).
!
    implicit none
    integer,intent(IN) :: kold,knew
    integer,intent(IN) :: ipn		!  current location in horiz.grid
    real,intent(IN)    :: yold_r4(kold),xold_r4(kold+1),xnew_r4(knew+1)
    real,intent(OUT)   :: ynew_r4(knew)
    logical,intent(IN) :: vrbos		!  if true, print diagnostics
    integer k,ko,n,kstart
    real*8 yold(kold),xold(kold+1),xnew(knew+1),ynew(knew),		&
         zold(kold+1),znew(knew+1),colin,clout,yinteg,xlo,xhi,ra,rb,	&
         scale
    logical at_top
    real,parameter     :: acurcy=1.e-6
!
    xold=xold_r4
    yold=yold_r4
    xnew=xnew_r4

!SMS$IGNORE BEGIN
    if (vrbos)								&
     write (*,101) ipn,'old',(k,xold(k),yold(k),k=1,kold),		&
      kold+1,xold(kold+1)
 101 format (i8,'  pcmadv -- ',a,' profile:',5x,'x',11x,'y'/		&
        (i27,f14.1,es12.4))
!SMS$IGNORE END
!
!   if (xold_r4(1).ne.xnew_r4(1)) write (*,*) ipn,			&
!    ' plmadv warning: bottom xold,xnew differ:',			&
!     xold_r4(1),xnew_r4(1)
!   if (xold_r4(kold+1).ne.xnew_r4(knew+1)) write (*,*) ipn,		&
!    ' plmadv warning: top xold,xnew differ:',				&
!     xold_r4(kold+1),xnew_r4(knew+1)
!
!-- make sure -x- increases with -k-. change sign of -xold/new- if necessary
    if (xold(1).lt.xold(kold+1)) then
      zold=xold
      znew=xnew
    else
      zold=-xold
      znew=-xnew
    end if
!
! --- column integrals (colin/clout) are computed for diagnostic purposes only
    scale=0.
    colin=0.
    clout=0.
    do 3 k=1,kold
    scale=scale+abs(yold(k))
 3  colin=colin+yold(k)*(zold(k+1)-zold(k))
!
! --- find ynew(k) by requiring 
! --- that the integral over y*dx from znew(k) to znew(k+1) be preserved.
!
    at_top=.true.
    kstart=1
    do 4 k=1,knew
    yinteg=0.
    xlo=znew(k  )
    xhi=znew(k+1)
!cc    if (vrbos) print '(a,2f9.3)','xlo,xhi =',xlo,xhi
    if (xhi.gt.xlo) then
      at_top=.false.
      do 5 ko=kstart,kold
      if (zold(ko+1).le.xlo) then
        kstart=ko+1
        go to 5
      end if
      if (zold(ko  ).ge.xhi) go to 1
      ra=max(xlo,min(xhi,zold(ko  )))
      rb=max(xlo,min(xhi,zold(ko+1)))
      yinteg=yinteg+yold(ko)*(rb-ra)
 5    continue
 1    ynew(k)=yinteg/(xhi-xlo)
    else if (at_top) then
      ynew(k)=yold(   1)
    else                              !  at end
      ynew(k)=yold(kold)
    end if
!cc    if (vrbos) print '(a,f11.1)','ynew =',ynew(k)
    clout=clout+ynew(k)*(znew(k+1)-znew(k))
 4  continue
!
!SMS$IGNORE BEGIN
    if (abs(clout-colin).gt.acurcy*scale*(zold(kold+1)-zold(1)))	&
      write (*,100) ipn,' pcmadv - column intgl.error',			&
       colin,clout,(clout-colin)/colin
 100 format (i8,a,2es14.6,es9.1)
    if (vrbos)								&
     write (*,101) ipn,'new',(k,xnew(k),ynew(k),k=1,knew),		&
      knew+1,xnew(knew+1)
!SMS$IGNORE END
    ynew_r4=ynew
    return
    end subroutine pcmadv


   subroutine equilb(thet,pres,kk,slak,vrbos,ipn,PrintDiags)

! --- expand thin layers at the expense of thick layers above and below.
! --- do this without changing theta

   use module_constants,only: rd,p1000

   implicit none
   integer,intent(IN) :: kk,ipn		! no. of layers, test point location
   logical,intent(IN) :: vrbos		! if true, print results at test point
   real,intent(IN)    :: thet(kk)	! pot.temp. in layers
   real,intent(INOUT) :: pres(kk+1)	! Exner fcn on interfaces
   real,intent(IN)    :: slak		! retardation coefficient
   logical,intent(IN) :: PrintDiags
   integer k,k1,k2,ncount
   real dp1,dp2,dp3,dp4,dp5,th1,th2,th3,th4,th5,dis3,dis4,	&
        ratio,goal,pnew(kk+1)
   logical event
   real,parameter :: thin=rd/p1000	! approx. 1 Pa in Exner fcn units

   if (vrbos) then
!SMS$IGNORE BEGIN
    write (6,99) ipn,'  equilb input profile:'
    do k2=1,kk,10
     write (6,100) (pres(k1),k1=k2,min(kk+1,k2+10) )
     write (6,101) (thet(k1),k1=k2,min(kk  ,k2+9 ) )
     write (6,100)
    end do
 99  format ('ipn=',i8,a,i9,a)
 100 format (11f7.1)
 102 format (11i7)
 101 format (5x,10f7.2)
!SMS$IGNORE END
   end if

! --- scenario 1: sequence of 5 thin-thick-thin-thick-thin layers

   pnew(:)=pres(:)
   ncount=0
   do 1 k=3,kk-2
    if (pnew(k-2).gt.pnew(1)-thin) go to 1
    dp1=pnew(k-2)-pnew(k-1)
    dp2=pnew(k-1)-pnew(k  )
    dp3=pnew(k  )-pnew(k+1)
    dp4=pnew(k+1)-pnew(k+2)
    dp5=pnew(k+2)-pnew(k+3)
    th1=thet(k-2)
    th2=thet(k-1)
    th3=thet(k  )
    th4=thet(k+1)
    th5=thet(k+2)
! --- look for small dp1,dp3,dp5 in combination with large dp2,dp4
    if (dp2.gt.dp1 .and. dp4.gt.dp5) then
     goal=.5*(dp3+min(dp2,dp4))		!  desired thknss of lyr 3
     if (dp2.gt.goal .and. dp4.gt.goal) then
! --- thin-thick-thin-thick-thin combination found -> inflate lyr 3
      dis3=min(dp2-goal,goal-dp3) * slak
      dis4=min(dp4-goal,goal-dp3) * slak
      if (th3.gt.th2 .and. th4.gt.th3) then
! --- theta conservation requires  dis3*(th3-th2)+dis4*(th4-th3)=0
       ratio=(th4-th3)/(th3-th2)
       if (ratio.gt.1.) then
        dis4=dis3/ratio
       else
        dis3=dis4*ratio
       end if
      end if
! --- ready to expand middle layer
      pnew(k  )=pnew(k  )+dis3
      pnew(k+1)=pnew(k+1)-dis4

      if (vrbos) then
!SMS$IGNORE BEGIN
       write (6,'(a,5f8.2)') 'thknss quintuplet',dp1,dp2,		&
        dp3,dp4,dp5, '          becomes',dp1,pnew(k-1)-pnew(k),		&
         pnew(k)-pnew(k+1),pnew(k+1)-pnew(k+2),dp5,			&
          'orig dis3,dis4, mod dis3,dis4,ratio =',			&
           min(dp2-goal,goal-dp3),min(dp4-goal,goal-dp3),		&
            dis3,dis4,ratio
!SMS$IGNORE END
      end if
      ncount=ncount+1
     end if
    end if
 1 continue

! --- scenario 2: sequence of 3 thin-thick-thin layers

!  do 2 k=2,kk-1
!   if (pnew(k).gt.pnew(1)-thin) go to 2
!   dp2=pnew(k-1)-pnew(k  )
!   dp3=pnew(k  )-pnew(k+1)
!   dp4=pnew(k+1)-pnew(k+2)
!   th2=thet(k-1)
!   th3=thet(k  )
!   th4=thet(k+1)
! --- look for small dp2,dp4 in combination with large dp3
!   if (dp3.gt.dp2 .and. dp3.gt.dp4) then
!    goal=.5*(dp3+max(dp2,dp4))		!  desired thknss of lyr 3
!    if (dp2.lt.goal .and. dp4.lt.goal) then
! --- thin-thick-thin combination found -> deflate lyr 3
!     dis3=min(goal-dp2,dp3-goal) * slak
!     dis4=min(goal-dp4,dp3-goal) * slak
!     if (th3.gt.th2 .and. th4.gt.th3) then
! --- theta conservation requires  dis3*(th3-th2)+dis4*(th4-th3)=0
!      ratio=(th4-th3)/(th3-th2)
!      if (ratio.gt.1.) then
!       dis4=dis3/ratio
!      else
!       dis3=dis4*ratio
!      end if
!     end if
! --- ready to shrink middle layer
!     pnew(k  )=pnew(k  )-dis3
!     pnew(k+1)=pnew(k+1)+dis4

!!SMS$IGNORE BEGIN
!     if (vrbos) then
!      write (6,'(a,i3,a,3f9.3/a,3f9.3)') 'k=',k,' thknss triple',	&
!       dp2,dp3,dp4,'            becomes',pnew(k-1)-pnew(k),		&
!        pnew(k)-pnew(k+1),pnew(k+1)-pnew(k+2)
!      write (6,'(a,3es11.3)') 'pres.chgs.,ratio=',-dis3,dis4,ratio
!     end if
!!SMS$IGNORE END

!     ncount=ncount+1
!    end if
!   end if
!2 continue

   event = ncount>15			!  find interesting cases

   do k=1,kk
    if (pnew(k+1).gt.pnew(k)+thin) then
     event=.true.
!SMS$IGNORE BEGIN
     write (*,'(a,i3)')						&
     'error: nonmonotonic pressure on return from equilb, k=',k
!SMS$IGNORE END
    end if
   end do

   if (event .and. PrintDiags) then
!SMS$IGNORE BEGIN
    write (6,99) ipn,'  equilb input profile:'
    do k2=1,kk,10
     write (6,100) (pres(k1),k1=k2,min(kk+1,k2+10) )
     write (6,101) (thet(k1),k1=k2,min(kk  ,k2+9 ) )
     write (6,100)
    end do
!SMS$IGNORE END
   end if
   if ((event .and. PrintDiags) .or. vrbos) then
!SMS$IGNORE BEGIN
    write (6,99) ipn,'  equilb output profile:',ncount,' inflations'
    do k2=1,kk,10
     write (6,100) (pnew(k1),k1=k2,min(kk+1,k2+10) )
     write (6,102) (int(1000.*(pnew(k1)-pres(k1))),k1=k2,min(kk+1,k2+10) )
     write (6,101) (thet(k1),k1=k2,min(kk  ,k2+9 ) )
     write (6,100)
    end do
!SMS$IGNORE END
   end if

   pres(:)=pnew(:)
   return
end subroutine equilb


subroutine diagkp(kk,z,zdot,theta,kappa,vrbos,ipn)

! --- diagnose vertical mixing coefficient resulting from manipulation
! --- of the air column (such as vertical advection, regridding,...)

! --- approach: transform diffusion eqn for theta
! ---  (d theta/dt)_z = (d/dz) [kappa d theta/dz]
! --- into
! ---  (dz/dt)_theta = -(d/d theta) [kappa d theta/dz]
! --- ('d' = 'partial'). represent kappa in layers as averages of kappa
! --- on interfaces. solve tridiagonal system for interface kappas.
! --- kappa values returned are   l a y e r   averages.
! --- kappa is assumed zero in first and last layer (ref. McDougall & Dewar).

! --- input variables:
! --- z(kk+1)    -  interface depths (meters, massless layers allowed)
! --- zdot(kk+1) -  rate of interface movement (m/s)
! --- theta(kk)  -  buoyancy (independent variable!)

! --- output: kappa (m^2/s) in layers whose thickness exceeds 'thresh'

   implicit none
   integer i,j,k,kk,kp
   real z(kk+1),theta(kk),zdot(kk+1),kappa(kk),zmid(kk),	&
        thmid(kk),zdotm(kk),a(kk+1,kk+1),a1(kk+1,kk+1),		&
        b(kk+1),b1(kk+1),d,dpidz
   integer indx(kk+1),klist(kk),ipn
   logical vrbos
   real,parameter :: thresh=1.e-3

! --- weed out massless layers
   kp=0
   do 5 k=1,kk
   klist(k)=0
   if (z(k+1).lt.z(k)-thresh) then	!  atmospheric case: z = Exn.fcn, decreasing with k
     kp=kp+1
     klist(k)=kp
     zmid(kp)=.5*(z(k)+z(k+1))		! mid layer depth
     thmid(kp)=theta(k)
     zdotm(kp)=zdot(k)
   end if
5  continue

   if (vrbos) then
!SMS$IGNORE BEGIN
     write (*,'(a,i8,a/(5(f9.1,f6.1)))') 'ipn =',ipn,			&
      '  input profile:',(z(k),theta(k),k=1,kk),z(kk+1)
     write (*,'(2(i5,a))') kp,' non-massless layers  =>',kp-3,' unknowns'
!SMS$IGNORE END
   end if

   a=0.
   do 10 k=3,kp-1
   a(k,k+1)= (thmid(k+1)-thmid(k-1))/(zmid(k+1)-zmid(k-1))
   a(k,k-1)=-(thmid(k  )-thmid(k-2))/(zmid(k  )-zmid(k-2))
   a(k,k  )=a(k,k+1)+a(k,k-1)
10 continue

   do 11 k=3,kp-1
11 b(k)=2.*zdotm(k)*(thmid(k-1)-thmid(k))

   if (vrbos) then
!SMS$IGNORE BEGIN
     write(*,*) 'ipn =',ipn,'  input matrix * 10^3'
     do k=3,kp-1
       write(*,'(3p,20f8.2)') (a(k,j),j=3,kp-1)
     end do
!SMS$IGNORE END
   end if

   a1=a
   b1=b

!SMS$IGNORE BEGIN
   if (vrbos) write(*,*) 'ipn =',ipn,'  input rhs * 10^3'
   if (vrbos) write(*,'(3p,20f8.2)')  (b(j),j=3,kp-1)
!SMS$IGNORE END

   call ludcmp(a(3,3),kp-3,kk+1,indx,d)
   call lubksb(a(3,3),kp-3,kk+1,indx,b(3))

   if (vrbos) then
!SMS$IGNORE BEGIN
! --- did the equation solver return a credible solution?
     do 22 k=3,kp-1
     b1(k)=0.
     do 22 i=3,kp-1
22   b1(k)=b1(k)+a1(k,i)*b(i)
     write(*,*) 'ipn =',ipn,'  matrix * kappa * 10^3'
     write(*,'(3p,20f8.2)')  (b1(j),j=3,kp-1)

     do k=2,kk-1
       j=klist(k)
       if (j.gt.1 .and. j.lt.kp) then
         dpidz=9.806/thmid(j)			!  exn.fcn => z conversion
         write(*,'(a,i4,6(a,f8.2))')	 				&
!!       write(*,'(a,i4,2(a,f8.1),a,es11.2,3(a,f8.1))')	 		&
         'k=',k,' z_uppr=',z(k),' zmid=',zmid(j),			&
         ' zdot(x10^3)=',1.e3*zdot(k),					&
         ' k_uppr(cm^2/s)=',1.e4*b(j)*dpidz**2,				&
         ' k_mid=',.5e4*(b(j)+b(j+1))*dpidz**2
       end if
     end do
!SMS$IGNORE END
   end if

! --- set diffusivity=0 in 1st and last layer
   b(   1)=0.
   b(   2)=0.
   b(kp  )=0.
   b(kp+1)=0.
   kappa=0.
   do k=2,kk-1
     dpidz=9.806/theta(k)			!  atmosphere
     j=klist(k)
     if (j.gt.1 .and. j.lt.kp) kappa(k)=.5*(b(j)+b(j+1))*dpidz**2
   end do
   return
end subroutine diagkp


!! SUBROUTINE ludcmp(a,n,np,indx,d)
!!       INTEGER n,np,indx(n),NMAX
!!       REAL d,a(np,np),TINY
!!       PARAMETER (NMAX=500,TINY=1.0e-20)
!!       INTEGER i,imax,j,k
!!       REAL aamax,dum,sum,vv(NMAX)
!!       d=1.
!! 
!!       do 12 i=1,n
!!         aamax=0.
!!         do 11 j=1,n
!!           if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
!! 11      continue
!!         if (aamax.eq.0.) pause 'singular matrix in ludcmp'
!!         vv(i)=1./aamax
!! 12    continue
!!       do 19 j=1,n
!!         do 14 i=1,j-1
!!           sum=a(i,j)
!!           do 13 k=1,i-1
!!             sum=sum-a(i,k)*a(k,j)
!! 13        continue
!!           a(i,j)=sum
!! 14      continue
!!         aamax=0.
!!         do 16 i=j,n
!!           sum=a(i,j)
!!           do 15 k=1,j-1
!!             sum=sum-a(i,k)*a(k,j)
!! 15        continue
!!           a(i,j)=sum
!!           dum=vv(i)*abs(sum)
!!           if (dum.ge.aamax) then
!!             imax=i
!!             aamax=dum
!!           endif
!! 16      continue
!!         if (j.ne.imax)then
!!           do 17 k=1,n
!!             dum=a(imax,k)
!!             a(imax,k)=a(j,k)
!!             a(j,k)=dum
!! 17        continue
!!           d=-d
!!           vv(imax)=vv(j)
!!         endif
!!         indx(j)=imax
!!         if(a(j,j).eq.0.)a(j,j)=TINY
!!         if(j.ne.n)then
!!           dum=1./a(j,j)
!!           do 18 i=j+1,n
!!             a(i,j)=a(i,j)*dum
!! 18        continue
!!         endif
!! 19    continue
!!       return
!! end subroutine ludcmp
!! !  (C) Copr. 1986-92 Numerical Recipes Software 'W3.
!! 
!! 
!! SUBROUTINE lubksb(a,n,np,indx,b)
!!       INTEGER n,np,indx(n)
!!       REAL a(np,np),b(n)
!!       INTEGER i,ii,j,ll
!!       REAL sum
!!       ii=0
!!       do 12 i=1,n
!!         ll=indx(i)
!!         sum=b(ll)
!!         b(ll)=b(i)
!!         if (ii.ne.0)then
!!           do 11 j=ii,i-1
!!             sum=sum-a(i,j)*b(j)
!! 11        continue
!!         else if (sum.ne.0.) then
!!           ii=i
!!         endif
!!         b(i)=sum
!! 12    continue
!!       do 14 i=n,1,-1
!!         sum=b(i)
!!         do 13 j=i+1,n
!!           sum=sum-a(i,j)*b(j)
!! 13      continue
!!         b(i)=sum/a(i,i)
!! 14    continue
!!       return
!! end subroutine lubksb
!  (C) Copr. 1986-92 Numerical Recipes Software 'W3.
!end module module_hybgen	! SMS doesn't like multiple routines in module
