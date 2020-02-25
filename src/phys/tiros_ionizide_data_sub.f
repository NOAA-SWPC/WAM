!
      subroutine tiros_ionize_data
     & (pres, lev1,levs,z,emaps,cmaps,djspectra,
     & grav,on,o2n, n2n,tn,gm_lat,essa1,tiros_activity_level,GW,
     & eden_aurora1D,me)
!    &   ,eflux,ch)     
! 
!vay-2015: pass den = rho from the TOP-level program take-out comput/arrays ntot & meanmass
!          version with data/xxx/-statements.......eden_aurora1D   3-density due to aurora
!
!         output: of tiros_ionize_data
!
!      use IDEA_ION_INPUT, only : EMAPS1 => EMAPS, CMAPS1 => CMAPS 
!      use IDEA_ION_INPUT, only : DJSPECTRA1 => DJSPECTRA
      use idea_mpi_def,      only : mpi_id
      use idea_composition,  only :  DTR, ELCH, R_2_d, PI
!     use tirosdata
      implicit none
!
      INTEGER :: j, i, m, l, tiros_activity_level, iband
      INTEGER :: levs, lev1
      INTEGER, parameter :: jmaxwell = 6

      real, intent(in) :: pres(levs)
      real :: pres1(levs)
      integer, intent(in) :: me
      real, intent(in) :: Z(levs),GRAV(levs),ON(levs),
     &        o2n(levs),n2n(levs)
     
      real :: gl, mlt, GW, gm_lat
      real :: bz, gscon, amu, e0
!
      real  :: emaps(21,20,7),cmaps(21,20,7),djspectra(15,21)
      real ::  NTOT(levs),meanmass(levs)
      REAL :: eden_aurora1D(levs)
!      goto 777
      
      real ::  QIONT(levs),RATIO(21),RLAM(21)
     &,den(levs),dl_lower,dl_upper,qiont_lower,qiont_upper
     &,tn(levs),mo,mo2,mn2,alpha
     &,rno,RANGE_en,pr,ratioz,rlamz,mh,mhe,q
     &,qiont_O(levs),qiont_O2(levs),qiont_N2(levs)
      real :: width(15),en(15),TE11(21),TE15(21),width_maxwell
      real :: ionchr(21),ratio_ch,en_maxwell(jmaxwell),dl(jmaxwell)
     &,qion_maxwell(levs),lognpres(8),
     &ion_recomb(8),logpres,rr
c
      
      
      real :: ch , chi , dfac , diff , dprof , ed ,
     &     eflux , essa1 , qdmsp , QT(levs) ,
     &     ri , rj , th , swbz , offset, THMagd
      INTEGER i1 , i2 , j1 , j2 , k , kk , ld , n , nn , jj , jjj
      data en/.37,.6,.92,1.37,2.01,2.91,4.19,6.,8.56,12.18,
     &17.3,24.49,36.66,54.77,81.82/
      data width/.158,.315,.315,.63,.631,1.261,1.26,2.522,
     &2.522,5.043,5.043,10.,14.81,22.13,33.06/
      data RLAM/1.49,1.52,1.51,1.48,1.43,1.37,1.30,1.22,
     &1.12,1.01,0.895,0.785,0.650,0.540,0.415,0.320,0.225,
     &0.14,0.08,0.04,0.0/
c
      DATA ionchr/.378 , .458 , .616 , .773 , .913 , 1.088 , 1.403 ,
     &     1.718 , 2.033 , 2.349 , 2.979 , 3.610 , 4.250 , 4.780 ,
     &     6.130 , 7.392 , 8.653 , 9.914 , 12.436 , 14.957 , 17.479/
      data lognpres/-3.425,-4.835,-5.918,-7.066,-7.784,-8.366,-9.314,
     &-10.507/
      data ion_recomb/3.20e-13,3.20e-13,2.75e-13,1.45e-13,1.13e-13,
     &8.30e-14,3.70e-14,2.00e-14/
!      DATA dprof/4.23E19 , 5.62E19 , 5.77E19 , 5.70E19 , 1.04E19 ,
!     &     1.03E20 , 1.22E20 , 1.23E20 , 0.00E19 , 8.36E19 , 2.37E20 ,
!     &     2.61E20 , 0.00E19 , 0.00E18 , 3.07E20 , 5.26E20 , 0.00E19 ,
!     &     0.00E18 , 0.00E18 , 8.57E20/
!      DATA qdmsp/15*0.0/
      if (me .eq. 0) then
!!      print *,'eden_aur=', maxval(eden_aurora1D), minval(eden_aurora1D) 
!      print *,'start_tn-adt=',maxval(tn),minval(tn)   
       print *,'levslev1=', levs,lev1
      endif
      do i=1,levs
      qiont(i) = 0.0
      qion_maxwell(i) = 0.0
      qiont_O(i) = 0.0
      qiont_O2(i) = 0.0
      qiont_N2(i) = 0.0
      eden_aurora1D(i) = 0.0
      enddo
!      goto 777
      do 16 m=1,21
   16 ratio(m) = (m-1)*0.05
      do iband=1,21
      te15(iband)=0.0
      te11(iband)=0.0
! the ionization rates will need to be normalized to TE11, which is
! the energy flux between 300eV and 20keV, which is provided by the
! TIROS energy influx maps emaps, rather than the energy from
! 300eV to 100keV, which is what the spectra were normalized to
!
! check the energy influx is normalized to 1 erg/cm2/s
      do 17 m=1,15
   17 TE15(iband)=TE15(iband)+djspectra(m,iband)*en(m)*width(m)*1.6E-06
! normalize with the energy influx 300eV to 20keV
      do 18 m=1,11
   18 TE11(iband)=TE11(iband)+djspectra(m,iband)*en(m)*width(m)*1.6E-06
      if(me .eq. 0) then
!      print *, iband, TE11(iband), TE15(iband)
      endif
      enddo
      bz = 1.38e-23
      gscon = 8.314e3
      mo = 16.
      mo2 = 32.
      mn2 = 28.
      mh = 1.
      mhe = 4.
      amu = 1.661e-27
      E0=0.035
      WIDTH_maxwell=0.050
      do j = 1,jmaxwell
      en_maxwell(j) = j*0.05 - 0.025
      enddo
! initialize qiont, etc
!      do i=1,levs
!      qiont(i) = 0.0
!      qion_maxwell(i) = 0.0
!      qiont_O(i) = 0.0
!      qiont_O2(i) = 0.0
!      qiont_N2(i) = 0.0
!      eden_aurora1D(i) = 0.0
!      enddo
! convert magnetic latitude from radians to degrees
      thmagd = gm_lat * R_2_D
!     print *, 'gm_lat   thmagd  essa1', gm_lat, thmagd, essa1
      th = abs(thmagd) - 50.
      if(abs(thmagd).le.50.) goto 200
! calculate magnetic hour angle from noon in gregrees
!      essa1 = (mlt + 12.)*15.
! now passed essa1 directly
      IF ( essa1.GE.360.0 ) THEN
          essa1 = essa1 - 360.0
      ELSEIF ( essa1.LT.0.0 ) THEN
          essa1 = essa1 + 360.
      ENDIF
cc  **
      l = tiros_activity_level - 2
      IF ( l.LT.1 ) l = 1
      IF ( l.GT.7 ) l = 7
! define dfac to scale qiont later with large GW and tiros_activity_level
! Added by Zhuxiao.Li
      dfac = 1.0
      IF (tiros_activity_level.gt.9 .and. GW.gt.96.0)
     &   dfac = GW/96.0
 
cc  **
      ri = essa1/18.0 + 11.
      i1 = ri                   ! i1 =int(ri) ?
      ri = ri - i1
      IF ( i1.GT.20 ) i1 = i1 - 20
      i2 = i1 + 1
      IF ( i2.GT.20 ) i2 = i2 - 20
      rj = th/2. + 1.
      j1 = rj                   ! j1 =int(rj) ?
      rj = rj - j1
      j2 = j1 + 1
!
      eflux = rj*ri*EMAps(j2,i2,l) + (1.-rj)*ri*EMAps(j1,i2,l)
     &        + rj*(1.-ri)*EMAps(j2,i1,l) + (1.-rj)*(1.-ri)
     &        *EMAps(j1,i1,l)
      eflux = 10.**(eflux)/1000.
!             print *, 'eflux   ', eflux
!
      ch = rj*ri*CMAps(j2,i2,l) + (1.-rj)*ri*CMAps(j1,i2,l) + rj*(1.-ri)
     &     *CMAps(j2,i1,l) + (1.-rj)*(1.-ri)*CMAps(j1,i1,l)
!      print *, 'ch   ', ch
! validation tests:
! to compare with figure 4 or 5 F-R and Evans 1987
! set ch to 5 different mean energies and
! set eflux to 1.0 mW/m2
!       ch = 2.98
!       eflux=1.0
! a useful thing to compare is the ionization rate profile for ch=2.98
! with the equivalent output assuming a Maxwellian spectrum using the
! other code ionize_ipe_3 with the same mean energy of 2.98
! in this case the profiles are similar, at other other values of ch
! they can be quite different.
!
!      print *, 'set for test ch eflux', ch, eflux
!
      IF ( ch.LT.0.378 ) ch = 0.379
!      IF ( ch.GT.17.479 ) WRITE (6,99001) ch
!
      DO 300 kk = 2 , 21
         IF ( ch.LE.ionchr(kk) ) THEN
            k = kk - 1
            GOTO 400
         ENDIF
 300  CONTINUE
!      if (me .eq. 0) then
!      print *,'k,kk,ch=',k,kk,ch
!      endif
 400  continue 
      if (me .eq. 0) then
!      print *,'400k,kk,ch=',k,kk,ch
      endif
      chi = ch - ionchr(k)
      diff = ionchr(kk) - ionchr(k)
      ratio_ch = chi/diff
!      if(ratio_ch.gt.1.) print *, 'ratio_ch out of bounds', ratio_ch
c
c
99001 FORMAT ('  ch value outof bound in tiros',f10.6)
!!! loop through ipe height levels

!      goto 777

      do 2000 i=lev1,levs
! stop at 1000km altitude
      if(z(i)*1.e-3 .gt. 1000.) goto 200                ! zwam < 1000, vay       
! set up neutral parameters
! ntot  total number density m-3
! pres  pressure Pa
! meanmass amu
! den neutral mass density kg/m3
!      grav(i)=-gr(i)/100.
!
! rewrite
!............ pres =constant......den = rho !!! all computed before
!             should be passed to "tiros", neutrals are fixed here

      ntot(i) = on(i)+o2n(i)+n2n(i)
!      pres1(i) = ntot(i)*bz*tn(i)
      meanmass(i) = (on(i)*mo+o2n(i)*mo2+n2n(i)*mn2)/ntot(i)
      den(i) = pres(i)*meanmass(i)/(gscon*tn(i))
!
! calculate ion recombination rate
!      data lognpres/-3.425,-4.835,-5.918,-7.066,-7.784,-8.366,-9.314,-10.507/
      logpres = log(pres(i))

      if (logpres.ge.-4.835)then
      rr=3.20e-13
      goto 450
      endif
      if (logpres.le.-10.507)then
      rr=2.00e-14
      goto 450
      endif
      do jjj=3,8
      if(lognpres(jjj).le.logpres)then
      jj=jjj-1
      goto 451
      endif
      enddo
  451 continue

!err      rr = ion_recomb(jj)-logpres*(ion_recomb(jjj)-ion_recomb(jj))/
!err     &(lognpres(jjj)-lognpres(jj))
!update Vay-2016/10
!NEMS.x             0000000000ADF0A4  tiros_ionize_data         449/450  idea_ion_empirmodels.f
!ion_recomb(jjj =>  ion_recomb(jj)
!
      rr = ion_recomb(jjj)-(lognpres(jjj)-logpres)*
     & (ion_recomb(jjj)-ion_recomb(jj))/(lognpres(jjj)-lognpres(jj))
  450 continue
!
!

 3000 format(1x,i5,2f10.2,1p3e9.2)
! loop through all energy bands
      DO 10 L=1,15
!      DL(L)=RNO*en(l)*EXP(-EN(L)/ALPHA)
      DL_lower = djspectra(L,K)
      DL_upper = djspectra(L,KK)
      RANGE_en=4.57E-05*EN(L)**1.75
      PR=RANGE_en*grav(i)
      RATIOZ=PRES(i)/PR
      IF(RATIOZ.GT.1.0) GOTO 20
!      DO 12 M=1,21
      DO 12 M=2,21
      IF(RATIOZ.GT.RATIO(M)) GOTO 12
      RLAMZ=RLAM(M-1)+(RATIOZ-RATIO(M-1))*(RLAM(M)-RLAM(M-1))/
     &(RATIO(M)-RATIO(M-1))
      GOTO 13
   12 CONTINUE
   13 CONTINUE
      GOTO 21
   20 RLAMZ=0.0
   21 CONTINUE
      QIONt_lower=den(i)*EN(L)*RLAMZ*DL_lower*WIDTH(l)*1.E7/RANGE_en/E0
      QIONt_upper=den(i)*EN(L)*RLAMZ*DL_upper*WIDTH(l)*1.E7/RANGE_en/E0
      QIONt(i)=qiont(i)+(ratio_ch*qiont_upper+(1.-ratio_ch)*qiont_lower)
     &*eflux/(ratio_ch*te11(kk)+(1.-ratio_ch)*te11(k))
   11 CONTINUE
   10 continue
! add 0 - 300eV as Maxwellian with ch, and eflux
      alpha = ch/2.
      RNO=eflux*6.24E12/2./ALPHA**3
      DO 110 L=1,jmaxwell
      DL(L)=RNO*en_maxwell(l)*EXP(-EN_maxwell(L)/ALPHA)
      RANGE_en=4.57E-05*EN_maxwell(L)**1.75
      PR=RANGE_en*grav(i)
      RATIOZ=PRES(i)/PR
      IF(RATIOZ.GT.1.0) GOTO 120
!      DO 112 M=1,21
      DO 112 M=2,21
      IF(RATIOZ.GT.RATIO(M)) GOTO 112
      RLAMZ=RLAM(M-1)+(RATIOZ-RATIO(M-1))*(RLAM(M)-RLAM(M-1))/
     &(RATIO(M)-RATIO(M-1))
      GOTO 113
  112 CONTINUE
  113 CONTINUE
      GOTO 121
  120 RLAMZ=0.0
  121 CONTINUE
      qion_maxwell(i)=qion_maxwell(i)+den(i)*en_maxwell(L)*RLAMZ*DL(L)
     &*WIDTH_maxwell/RANGE_en/E0
  110 CONTINUE
      qiont(i) = qiont(i) + qion_maxwell(i)

!  the qiont scaled by dfac with large GW and tiros_activity level
!  added by Zhuxiao
       qiont(i) = qiont(i)*dfac
!
!vay-2016, extra security >0 and non-zero
!
      if (rr.gt.0.) then 
         eden_aurora1D(i)=sqrt(qiont(i)/rr)
      else
         eden_aurora1D(i)=0.
      endif
!
      q=qiont(i)/(0.92*n2n(i)+1.5*o2n(i)+0.56*on(i))
      qiont_O(i)=(0.5*o2n(i)+0.56*on(i))*q
      qiont_O2(i)=o2n(i)*q
      qiont_N2(i)=0.92*n2n(i)*q
 2000 continue
  777 continue   
!      do i=1,levs

!      eden_aurora1D(i) = 0.0
!      enddo


  200 continue

      if (me .eq. 0) then
      print *,'eden_aur=', maxval(eden_aurora1D), minval(eden_aurora1D) 
      print *,'tn-adt=',maxval(tn),minval(tn)   
      endif


      RETURN
      END subroutine tiros_ionize_data
