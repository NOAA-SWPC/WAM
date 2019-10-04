module module_op_diag
contains
!*********************************************************************
!	op_diag
!	Calculate derived variables for output from FIM global model
!	S. Benjamin    Feb 2008
!       S. Benjamin  - Apr 2008
!                      -  mods for use of theta-v in th3d prog variable
!                         instead of previous non-virt theta in th3d
!       S. Benjamin  - May 2008
!                      - calculation of isobaric 25-mb fields 
!                         on native icosahedral horizontal grid
!                         added in array g3p
!                        Will allow improved isobaric fields to be
!                         output from FIMpost, with improved reduction
!                         from full multivariate calculations using
!                         all fields (which are available in FIM but not FIMpost)
!*********************************************************************

subroutine op_diag(         &
  its,nts,                 &	! index time step, final timestep
  us3d,vs3d,dp3d,          &	! west wind, south wind, delta pres 
  pr3d,ex3d,mp3d,          &	! pressure, Exner, mont pot, 
  tr,      vor,ws3d,       &    ! tracers
  ph3d,rn2d,rc2d,          &    ! geopotential,accumulated precipitation/rainfall
  ts2d,us2d,hf2d,qf2d,sw2d,lw2d,st3d,sm3d  &
!   Below are output variables from op_diag &
  ,tk3d, rh3d, td3d, pw2d, mslp, g3p, g2d,time &
    )


use module_constants
use module_control,only: nvl,nvlp1,nip,ntra,ntrb,nvlp,nvarp,nvar2d,  &
                         dt,glvl,curve,yyyymmddhhmm,ArchvTimeUnit,pres_pa
use restart, only: write_restart

USE MACHINE  , ONLY : kind_rad
USE FUNCPHYS , ONLY : fpvs, fpvsl 
use module_outqv           ,only: outqv
use module_outqv_mn        ,only: outqv_mn
use module_outqv_mn_lat    ,only: outqv_mn_lat
use module_outqv_mn_lat_abs,only: outqv_mn_lat_abs
use module_out4d_mn        ,only: out4d_mn

!USE FUNCPHYS , ONLY : fpvs, ftdp
implicit none

! External variable declarations:
integer,intent(IN)    :: its,nts
!SMS$DISTRIBUTE (dh,nip) BEGIN
real   ,intent(IN)    :: us3d(nvl  ,nip),vs3d(nvl,nip),dp3d(nvl,nip)
real   ,intent(IN)    :: pr3d(nvlp1,nip)
real   ,intent(IN)    :: mp3d(nvl  ,nip)
real   ,intent(IN)    :: vor (nvl  ,nip) !not used
real   ,intent(IN)    :: ws3d(nvl  ,nip)
real   ,intent(IN)    :: ph3d(nvlp1,nip)
real   ,intent(IN)    :: ex3d(nvlp1,nip)
real   ,intent(IN)    :: tr(nvl    ,nip,ntra+ntrb)
real   ,intent(OUT) :: rh3d(nvl  ,nip)
real   ,intent(OUT) :: tk3d(nvl  ,nip)
real   ,intent(OUT) :: g3p (nvlp ,nip, nvarp)
real   ,intent(OUT) :: g2d (nip  ,nvar2d)

!       nvar2d indices
!       1  -  height cloud base (in meters above sea level (ASL))
!       2  -  pressure cloud base (Pa)
!       3  -  height cloud top (m ASL)
!       4  -  pressure cloud top (Pa)
!       5  -  column relative humidity with respect to saturated column precipitable water (0-1)

real                  :: cloudbase(nip)
real   ,intent(INOUT) :: td3d(nvl,nip)
real   ,intent(INOUT) :: pw2d(nip)
real   ,intent(INOUT) :: mslp(nip)

real   ,intent(IN)    :: rn2d(nip),rc2d(nip)
real   ,intent(IN)    :: ts2d(nip),us2d(nip),hf2d(nip),qf2d(nip),sw2d(nip),lw2d(nip),st3d(4,nip),sm3d(4,nip)
real                  :: th3d(nvl,nip),qv3d(nvl,nip),qw3d(nvl,nip)
real                  :: dummy(nvl,nip)
!  pkap = (p/p0)**(R/Cp)  i.e., p(in bars) to the kappa (Rd/Cp) power
!  pkap = Exner / Cp
real                  :: pkap (nvl), hgt(nvl), pres(nvl), tvirt(nvl), thv(nvl), thnv(nvl)
real                  :: pklv (nvlp1)
real                  :: pkapp (nvlp), prsp(nvlp),tsfc
integer               :: nsmooth, isn, ism
!real                 :: oz3d(nvl,nip)
real                  :: rhpw(nip), satpw(nip)
!SMS$DISTRIBUTE END

integer,intent(in)   :: time
integer, parameter   :: lunout=40,lun3d=50,lun2d=60  ! Logical units for output
character(80)        :: filename
character(80)        :: header(10)
character(4)         :: var(5)
integer              :: ivl,ipn,nop=0,nrs=0
integer              :: nv,ivlp, ivar
real                 :: drh, du, dv, rh, fact
real                 :: maxqv3d,minqv3d,aveqv3d,maxdp3d,mindp3d,avedp3d
real                 :: esw,qsw,es,esln
real                 :: thet,dthet,dpkap,dhgt
real                 :: gam, gamd, gams, ex1, ex1s, ex1sinv,tt1,t6
real                 :: tsfcnew, thbar
real                 :: cloud_def_p, zcldbase, pcldbase, zcldtop, pcldtop, watericemax

integer              :: LB

real (kind=kind_rad) :: tk

! 6.5 K/km - Standard lapse rate
        DATA     GAMs    /0.0065/
! 10  K/km - Dry adiabatic lapse rate
        DATA     GAMD    /0.0100/
        DATA     cloud_def_p /0.00001/
!       DATA     cloud_def_p /0.000001/

!  Variable names for outqv* printout
   DATA var /'Hgt ','Temp','RH  ','Uwnd','Vwnd' /

nsmooth = 1 

!there is an SMS problem writing tr(:,:,1)
  th3d(:,:)=tr(:,:,1)
  qv3d(:,:)=tr(:,:,2)
  qw3d(:,:)=tr(:,:,3)
! oz3d(:,:)=tr(:,:,4)
  header(1:10) = '                                                                                '
! Calculate isobaric levels every 25 hPa
  do ivlp = 1,nvlp
    ! prsp(ivlp) = p1000 -(float(ivlp-1) * 2500.)     ! in Pa
    prsp(ivlp) = float(pres_pa(ivlp))     ! in Pa
!   pkapp(ivlp) = (1.- float(ivlp-1)* 0.025)**(rd/cp)
    pkapp(ivlp) = (prsp(ivlp)/100000.)**(rd/cp)
    print *, "in op_diag calc levels: prsp(",ivlp,"): ",prsp(ivlp)," pkapp(",ivlp,"): ", pkapp(ivlp)
  end do

  ex1s   = rd*gams/grvity
  ex1sinv= 1./ex1s

!SMS$PARALLEL (dh,ipn) BEGIN

!............................................................
!  Beginning of horizontal loop
!............................................................
  do ipn=1,nip
    pw2d(ipn) = 0.
    satpw (ipn)     = 0.

!  Initialize isobaric fields as -999999. = -spval_p
      g3p(:,ipn,:) = -spval_p
      g2d(ipn,:)   = -spval_p

    do ivl=1,nvl
!   PW calculation here uses water vapor and condensate.
      pw2d(ipn) = pw2d(ipn)+dp3d(ivl,ipn)*(tr(ivl,ipn,2)+tr(ivl,ipn,3))/grvity
!     pw2d(ipn) = pw2d(ipn)+dp3d(ivl,ipn)*(tr(ivl,ipn,2))/grvity

!rb   pres(ivl) = 0.5*(pr3d(ivl,ipn)+pr3d(ivl+1,ipn))
!rb   hgt1(ivl)  = 0.5*(ph3d(ivl,ipn)+ph3d(ivl+1,ipn)) / 9.8
!rb   pkap(ivl) = (pres(ivl)/p1000)**(rd/cp)

      if (pr3d(ivl,ipn).gt.pr3d(ivl+1,ipn)+0.1) then
        pkap(ivl) = (ex3d(ivl  ,ipn)*pr3d(ivl  ,ipn)                     &
                    -ex3d(ivl+1,ipn)*pr3d(ivl+1,ipn))/                   &
                   ((cp+rd)*(pr3d(ivl,ipn)-pr3d(ivl+1,ipn)))
      else
        pkap(ivl)=0.5*(ex3d(ivl,ipn)+ex3d(ivl+1,ipn))/cp
      end if
      pres(ivl)=p1000*pkap(ivl)**(cp/rd)
      hgt(ivl)=(ph3d(ivl,ipn)                                           &
       +(ex3d(ivl,ipn)-cp*pkap(ivl))*th3d(ivl,ipn)) / 9.8

      tvirt(ivl) = th3d(ivl,ipn) * pkap(ivl)                ! virtual temp
!   temperature
      tk3d(ivl,ipn)=tvirt(ivl) / (1.+0.6078*qv3d(ivl,ipn)) ! temperature
      thv  (ivl) = th3d(ivl,ipn)                           ! virtual pot temp
      thnv (ivl) = th3d(ivl,ipn)/(1.+0.6078*qv3d(ivl,ipn)) ! non-virt pot temp 
      tk  = tk3d(ivl,ipn)
!  sat vapor pressure w.r.t. water/ice combination
!     esw=fpvs(tk)
!  sat vapor pressure w.r.t. water (liquid) only
      esw=fpvsl(tk)
      esw=min (real(fpvsl(tk)) , pres(ivl)) ! qsw <= 1, as it has to be

      qsw=0.62197*esw / (pres(ivl)-0.37803*esw)
      qsw = max(1.e-8,min(qsw,0.1))
      dummy(ivl,ipn) = qsw
      rh3d(ivl,ipn)=100.*max(0., min(1., tr(ivl,ipn,2)/qsw))
      ! rh for 0-100, as done in GRIB for other models
!   PW in saturated column
      satpw(ipn) = satpw(ipn)+(dp3d(ivl,ipn)*qsw/grvity)

      es  = pres(ivl)*(tr(ivl,ipn,2)+1.e-8)/(0.62197+(tr(ivl,ipn,2)+1.e-8))
      esln  = log(es)
      td3d(ivl,ipn)= (35.86*esln-4947.2325)/(esln-23.6837)
    end do

    rhpw(ipn) =  pw2d(ipn) / satpw(ipn)

!   Exner-like variables
!   ===================
!  pkap = (p/p0)**(R/Cp)  i.e., p(in bars) to the kappa (Rd/Cp) power
!  pkap = Exner / Cp
!   ===================
!   pklv  - Exner fn / Cp (0-1) on nvl+1 LEVELS (not layer midpoints)
!   pkap  -  "      "        "     nvl   LAYER midpoints
!   pkapp -  "      "        "  on nvlp  ISOBARIC levels

!   temperature variables
!   ===================
!   th3d, thv - virtual  potential temperature on nvl LAYER midpoints
!   thnv      - NON-virt     "         "        " "   "      "
!   tvirt     - virtual           temperature on nvl LAYER midpoint
!   tk3d      - nonvirtual           "        "   "   "      "
!

    do ivl=1,nvlp1
      pklv(ivl)= ex3d(ivl,ipn)/cp
    end do

!   Calculate isobaric values

!    -- Step 1
!   Calculate lapse rate (in virtual temp!) near surface:
!      Level 6 should be about 60 hPa above surface.
!      To avoid excessive effect (warm or cold) from lowest level,
!        use theta(ivl=3) with pres(ivl=1) to obtain
!        estimated temp at level 1 (also used in RUC post - hb2p.f)

      tt1=thnv (3 )*pkap(1)
      t6 =thnv (6 )*pkap(6)
      gam = (tt1-t6)/(hgt(6)-hgt(1))
      gam = min (gamd, max(gam,gams))
      ex1 = rd*gam/grvity
      tsfc=thnv(3 )*(pr3d(1,ipn)/p1000)**(rd/cp)  ! tsfc only for reduction here

!   tsfc      - virtual temp at lowest LEVEL
!   tt1       - virtual temp at lowest LAYER midpoint

! Step 2 - obtain isobaric variables (except height) IF terrain elevation allows
    do ivlp = 1,nvlp
      do ivl=2,nvl
        dthet = thnv (ivl) - thnv (ivl-1)
        dpkap  = pkap(ivl)  - pkap(ivl-1)
        drh   = rh3d(ivl,ipn) - rh3d(ivl-1,ipn)
        du    = us3d(ivl,ipn) - us3d(ivl-1,ipn)
        dv    = vs3d(ivl,ipn) - vs3d(ivl-1,ipn)
       if (pkap(ivl) .lt. pkapp(ivlp) .and.pkap(ivl-1).ge.pkapp(ivlp) ) then
        thet  = thnv (ivl-1) &
                 + (pkapp(ivlp)-pkap(ivl-1)) * dthet/dpkap
        rh    = rh3d(ivl-1,ipn) &
                 + (pkapp(ivlp)-pkap(ivl-1)) * drh  /dpkap
        g3p(ivlp,ipn,4)    = us3d(ivl-1,ipn) &
                 + (pkapp(ivlp)-pkap(ivl-1)) * du   /dpkap
        g3p(ivlp,ipn,5)    = vs3d(ivl-1,ipn) &
                 + (pkapp(ivlp)-pkap(ivl-1)) * dv   /dpkap
        g3p(ivlp,ipn,2) = thet*pkapp(ivlp)
        g3p(ivlp,ipn,3) = rh             
        if(ntrb.gt.0)then
          do nv=ntra+1,ntra+ntrb
             dv    = tr(ivl,ipn,nv) - tr(ivl-1,ipn,nv)
             g3p(ivlp,ipn,5+nv-ntra-1)    = tr(ivl-1,ipn,nv) &
                   + (pkapp(ivlp)-pkap(ivl-1)) * dv   /dpkap
          enddo
        endif
       end if
      end do  ! ivl

    end do    ! ivlp

! Step 3 - Now, obtain isobaric heights IF terrain elevation allows
    do 10 ivlp = 1,nvlp
      do ivl=2,nvlp1
        if (pklv(ivl) .le. pkapp(ivlp) .and.pklv(ivl-1).gt.pkapp(ivlp) ) then
        dpkap  = pklv(ivl)  - pklv(ivl-1)
        dhgt  = (ph3d(ivl,ipn) - ph3d(ivl-1,ipn))/9.8
        g3p(ivlp,ipn,1) = ph3d(ivl-1,ipn)/9.8 + (pkapp(ivlp)-pklv(ivl-1)) * dhgt/dpkap 
        go to 10
        end if
      end do
10  continue


! Step 4 - REDUCE from atmosphere above to obtain heights/temps below terrain surface
!   Set other isobaric variables also below terrain surface
    do ivlp = nvlp-1,1,-1
      if (pkap(1) < pkapp(ivlp)) then
         g3p(ivlp,ipn,2) = tt1*(prsp(ivlp)/pres(1))**ex1s
         g3p(ivlp,ipn,1) = ph3d(1,ipn)/9.8 - (g3p(ivlp,ipn,2)-tsfc)/gams
         g3p(ivlp,ipn,3) = g3p(ivlp+1,ipn,3)
         g3p(ivlp,ipn,4) = g3p(ivlp+1,ipn,4)
         g3p(ivlp,ipn,5) = g3p(ivlp+1,ipn,5)
        if(ntrb.gt.0)then
         do nv=ntra+1,ntra+ntrb
           g3p(ivlp,ipn,nv)    = g3p(ivlp+1,ipn,nv)
         enddo
        endif
      end if
    end do

!   do ivlp = nvlp-3, nvlp
! Step 5 - EXTRAPOLATE to obtain height/temp above top native level
!   Set other isobaric variables also above top native level
!     if (g3p(ivlp,ipn,2) < 10.) then
!        g3p(ivlp,ipn,2) = tk3d(nvl,ipn)     !  isothermal lapse rate
!        g3p(ivlp,ipn,1) = hgt(nvl) + rd*tk3d(nvl,ipn)/grvity &
!                * alog(pres(nvl)/prsp(ivlp))
!        g3p(ivlp,ipn,3) = g3p(ivlp-1,ipn,3)
!        g3p(ivlp,ipn,4) = g3p(ivlp-1,ipn,4)
!        g3p(ivlp,ipn,5) = g3p(ivlp-1,ipn,5)
!       if(ntrb.gt.0)then
!        do nv=ntra+2,ntra+ntrb
!          g3p(ivlp,ipn,nv)    = g3p(ivlp-1,ipn,nv)
!        enddo
!       endif
!     end if
!   end do

!   calculate mean pot temp from surface (actually, theta at level 2) down to sea level
!     using 0.8e-3 lapse rate
    thbar = th3d(2,ipn) - ph3d(1,ipn)* 0.5 * 0.8e-3 / grvity
! Calculate sea-level pressure reduction
    mslp(ipn) = p1000 * ((pr3d(1,ipn)/p1000)**(rd/cp) + ph3d(1,ipn)/(cp*thbar))**(cp/rd) 

!===============================================================
!  Cloud base/top diagnosis
!===============================================================

! g2d(ipn,1) = -spval_p
 watericemax = 0.
 zcldbase    = -spval_p
 zcldtop     = -spval_p
 pcldbase    = 0.
 pcldtop     = 0.

 do ivl = 1,nvl
   watericemax = max(watericemax, tr(ivl,ipn,3) )
 end do

 if (watericemax < cloud_def_p) go to 500

!     At surface?
 if (tr(1,ipn,3) > cloud_def_p) then
   zcldbase = hgt (1)
   pcldbase = pres(1)
   go to 400
 end if

!    Cloud aloft?
   do ivl = 3,nvl
     if (tr(ivl,ipn,3) > cloud_def_p) then
       zcldbase = hgt(ivl-1) + (tr(ivl,ipn,3)-cloud_def_p)   &
                *(hgt(ivl  ) - hgt(ivl-1))                    &
           / max(cloud_def_p,(tr (ivl, ipn,3)-tr(ivl-1,ipn,3)))
!      zcldbase = hgt(ivl)

       pcldbase = pres(ivl-1) + (tr(ivl,ipn,3)-cloud_def_p)  &
                *(pres(ivl  ) - pres(ivl-1))             &
          / max(cloud_def_p,(tr (ivl, ipn,3)-tr(ivl-1,ipn,3)))
!      pcldbase = pres(ivl)
       go to 400
     end if
   end do

400  continue

  do ivl = nvl-2, 2, -1
     if (tr(ivl,ipn,3) .gt. cloud_def_p) then
       zcldtop = hgt(ivl) + (tr(ivl,ipn,3)-cloud_def_p)  &
                *(hgt(ivl+1) - hgt(ivl))             &
             / max(cloud_def_p,(tr (ivl, ipn,3)-tr(ivl+1,ipn,3)))
!      zcldtop = hgt(ivl)
       pcldtop = pres(ivl) + (tr(ivl,ipn,3)-cloud_def_p)  &
                *(pres(ivl+1) - pres(ivl  ))             &
          / max(cloud_def_p,(tr (ivl, ipn,3)-tr(ivl+1,ipn,3)))
!      pcldtop = pres(ivl)
       go to 500
     end if
  end do

500  continue

  cloudbase(ipn) = zcldbase
  g2d(ipn,1) = zcldbase
  g2d(ipn,2) = pcldbase
  g2d(ipn,3) = zcldtop
  g2d(ipn,4) = pcldtop
  g2d(ipn,5) = rhpw(ipn)
  g2d(ipn,6) = satpw(ipn)


!............................................................
  end do     ! ipn loop - primary horizontal loop
!............................................................

!SMS$PARALLEL END


    LB = LBOUND(g2d,1)
    write (6,'(a,i9,1x,2a,i8)') ' MAPS SLP at ',time,ArchvTimeUnit,', time step=',its
    call outqv_mn    (  mslp ,1,nip,its ,0.01)
    write (6,'(a,i9,1x,2a,i8)') ' MAPS SLP  - max/min  at ',time,ArchvTimeUnit,', time step=',its
    call outqv       (  mslp ,deg_lat,deg_lon,nip,1,0.01)
    write (6,'(a,i9,1x,2a,i8)') ' Precipitable water at ',time,ArchvTimeUnit,', time step=',its
    call outqv_mn    (  pw2d ,1,nip,its ,1.)
    write (6,'(a,i9,1x,2a,i8)') ' Precipitable water  - max/min  at ',time,ArchvTimeUnit,', time step=',its
    call outqv       (  pw2d ,deg_lat,deg_lon,nip,1 ,1.)

    write (6,'(a,i9,1x,2a,i8)') ' SatPW  at ',time,ArchvTimeUnit,', time step=',its
    call outqv_mn    (  g2d(LB,6) ,1,nip,its ,1.)
    write (6,'(a,i9,1x,2a,i8)') ' SatPW  - max/min  at ',time,ArchvTimeUnit,', time step=',its
    call outqv       (  g2d(LB,6) ,deg_lat,deg_lon,nip,1,1.)

    write (6,'(a,i9,1x,2a,i8)') ' RH w.r.t. PW  at ',time,ArchvTimeUnit,', time step=',its
    call outqv_mn    (  g2d(LB,5) ,1,nip,its ,1.)
    write (6,'(a,i9,1x,2a,i8)') ' RH w.r.t. PW  - max/min  at ',time,ArchvTimeUnit,', time step=',its
    call outqv       (  g2d(LB,5) ,deg_lat,deg_lon,nip,1,1.)

    write (6,'(a,i9,1x,2a,i8)') ' Temperature at ',time,ArchvTimeUnit,', time step=',its
    call outqv_mn    (  tk3d ,nvl,nip,its ,1.)
    write (6,'(a,i9,1x,2a,i8)') ' Temperature  -max/min  at ',time,ArchvTimeUnit,', time step=',its
    call outqv       (  tk3d ,deg_lat,deg_lon,nip,nvl,1.)
    write (6,'(a,i9,1x,2a,i8)') ' RH at ',time,ArchvTimeUnit,', time step=',its
    call outqv_mn    (  rh3d ,nvl,nip,its ,1.)
    write (6,'(a,i9,1x,2a,i8)') ' RH  - max/min  at ',time,ArchvTimeUnit,', time step=',its
    call outqv       (  rh3d ,deg_lat,deg_lon,nip,nvl,1.)

    write (6,'(a,i9,1x,2a,i8)') ' qsw at ',time,ArchvTimeUnit,', time step=',its
    call outqv_mn    (  dummy ,nvl,nip,its ,1.)
    write (6,'(a,i9,1x,2a,i8)') ' qsw  - max/min  at ',time,ArchvTimeUnit,', time step=',its
    call outqv       (  dummy ,deg_lat,deg_lon,nip,nvl,1.)

    write (6,'(a,i9,1x,2a,i8)') ' Cloudbase  - max/min  at ',time,ArchvTimeUnit,', time step=',its
    call outqv       (  cloudbase ,deg_lat,deg_lon,nip,1,1.)
    write (6,'(a,i9,1x,2a,i8)') ' Cloud base height  - max/min  at ',time,ArchvTimeUnit,', time step=',its
    call outqv       (  g2d(LB,1) ,deg_lat,deg_lon,nip,1,1.)
    write (6,'(a,i9,1x,2a,i8)') ' Cloud base pressure  - max/min  at ',time,ArchvTimeUnit,', time step=',its
    call outqv       (  g2d(LB,2) ,deg_lat,deg_lon,nip,1,1.)
    write (6,'(a,i9,1x,2a,i8)') ' Cloud top height  - max/min  at ',time,ArchvTimeUnit,', time step=',its
    call outqv       (  g2d(LB,3) ,deg_lat,deg_lon,nip,1,1.)
    write (6,'(a,i9,1x,2a,i8)') ' Cloud top pressure  - max/min  at ',time,ArchvTimeUnit,', time step=',its
    call outqv       (  g2d(LB,4) ,deg_lat,deg_lon,nip,1,1.)

    write (6,'(a,i9,1x,2a,i8)') ' Soil moisture at ',time,ArchvTimeUnit,', time step=',its
    call outqv_mn    (  sm3d ,4,nip,its ,1.)
    write (6,'(a,i9,1x,2a,i8)') ' Soil moisture  - max/min  at ',time,ArchvTimeUnit,', time step=',its
    call outqv       (  sm3d ,deg_lat,deg_lon,nip,4,1.)
    
    fact = 1000.
    write (6,'(a,i9,1x,2a,i8)') ' Sensible heat flux at ',time,ArchvTimeUnit,', time step=',its
    call outqv_mn (hf2d,1,nip,its,1.)
    write (6,'(a,i9,1x,2a,i8)') ' Sensible heat flux at ',time,ArchvTimeUnit,', time step=',its
    call outqv_mn_lat (hf2d,deg_lat,deg_lon,30.,1,nip,its,fact)
    write (6,'(a,i9,1x,2a,i8)') ' Latent heat flux at ',time,ArchvTimeUnit,', time step=',its
    call outqv_mn (qf2d,1,nip,its,1.)
    write (6,'(a,i9,1x,2a,i8)') ' Latent heat flux at ',time,ArchvTimeUnit,', time step=',its
    call outqv_mn_lat (qf2d,deg_lat,deg_lon,30.,1,nip,its,fact)

  LB = LBOUND(g3p,2)
  do ivar=1,5
 
    write (6,'(2a,i9,1x,2a,i8)') var(ivar),'-isobaric at ',time,ArchvTimeUnit,', time step=',its
    call outqv_mn    (  g3p(1,LB,ivar),nvlp,nip,its,1.)
    write (6,'(2a,i9,1x,2a,i8)') var(ivar),'-isobaric  - max/min  at ',time,ArchvTimeUnit,', time step=',its
    call outqv       (  g3p(1,LB,ivar) ,deg_lat,deg_lon,nip,nvlp,1.)
  end do

  
!SMS$PARALLEL (dh,ipn) BEGIN
    do ipn = 1,nip
      do ivl=1,nvl
        dummy (ivl,ipn) = tk3d(ivl,ipn) - td3d(ivl,ipn)
      end do
    end do
!SMS$PARALLEL END
    write (6,'(a,i9,1x,2a,i8)') ' T-Td at ',time,ArchvTimeUnit,', time step=',its
    call outqv_mn    (  dummy ,nvl,nip,its ,1.)
    write (6,'(a,i9,1x,2a,i8)') ' T-Td  - max/min  at ',time,ArchvTimeUnit,', time step=',its
    call outqv       (  dummy ,deg_lat,deg_lon,nip,nvl,1.)

 

return
end subroutine op_diag
end module module_op_diag
