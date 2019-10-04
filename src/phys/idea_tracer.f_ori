! Apr 06 2012    Henry Juang, initial implement for nems
! Oct 20 2015    Weiyu Yang,  add f10.7 inputted data.
!
! Oct 15 2016    VAY, take out Weiyu f10.7 and put them as a parameter
! Nov    2016    Correct implementation of JO2 and simplified O2-O chemistry
! Dec-Jan 16/17  Addition of ozone "diurnal" chemistry and extra reactions
!-----------------------------------------------------------------------

      module idea_tracer_mod
!-----------------------------------------------------------------------
! hold jprofile-old; oh-ho2 global profiles
!--------------------------------------------

      implicit none
!
      real, allocatable::  jj(:)
      real, allocatable::  oh(:), ho2(:)
      real, allocatable:: vmr_glob(:,:)     !(levs, nvmr)
      end module idea_tracer_mod
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine idea_tracer_init(levs)
!
!
!
      use idea_tracer_mod,    only :   jj, oh, ho2, vmr_glob
      use idea_tracers_input, only :
     &    jprofile, hprofile, init_tracer_constants,  WAM_GLOBAL_TRACERS
      use IDEA_MPI_def,       only : mpi_WAM_quit

      implicit none
      real, parameter :: f107=100.  ! just any values of f107 to init old 1D-JJs
                                    ! now updated with time...JO2-3D
      integer, parameter :: nvmr=15 !  15-global vertical arrays returned by WAM_GLOBAL_TRACERS
      integer, intent(in):: levs    !number of pres levels
!
      allocate (jj(levs))
      allocate (oh(levs), ho2(levs))
      allocate (vmr_glob(levs, nvmr))

      call jprofile(levs,f107,jj)    ! old 1D-Jo2 profile for day and night
      call hprofile(levs,oh,ho2)     ! 1D [oh-ho2]  profiles can be replaced on monthly var. zonal fields
      call init_tracer_constants     ! comments on the chemical mechanisms with reaction constants
!
      call WAM_GLOBAL_TRACERS(levs, nvmr, vmr_glob)
      print *, ' VAY WAM_GLOBAL_TRACERS INIT'
!
!     call  mpi_WAM_quit(iret, message)
!
      return
      end subroutine idea_tracer_init
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine idea_tracer(im,ix,levs,ntrac,ntrac_i,grav,prsi,prsl,   
     &  adt,q,dtp,n1,n2,ozn, n3,n,rho,am, am29, 
     &  cospass,dayno,zg,f107,f107d,me)
!
      use physcons, only          : avgd => con_avgd             
!     &                     amo3 => con_amo3,amh2o => con_amw
      use idea_composition, only  : bz,amo,amn2, amo2, amo3, amh2o
      use idea_composition, only  : rbz, rmo, rmo2, rmn2, rmh2o,rmo3
!
!      use idea_tracers_input,only : mo, mo2, mn2, mh2o,mo3
!
      use idea_tracer_mod,   only : jj
      implicit none
! Argument
      integer, intent(in) :: me              ! current PE
      integer, intent(in) :: im              ! number of data points in up,dudt(first dim)
      integer, intent(in) :: ix              ! max data points in fields
      integer, intent(in) :: levs            ! number of pressure levels
      integer, intent(in) :: ntrac           ! number of tracer (total)
      integer, intent(in) :: ntrac_i         ! number of tracer add by IDEA

      integer, intent(in) :: dayno           ! calendar day
      real, intent(in)    :: cospass(im)     ! cos zenith angle
      real, intent(in)    :: zg(ix,levs)     !layer height (m)
      real, intent(in)    :: f107, f107d     ! variable F107 to recomput JJ
!
      real, intent(in)    :: prsi(ix,levs+1) ! interface pressure in KPa
      real, intent(in)    :: prsl(ix,levs)   ! layer pressure in KPa
      real, intent(in)    :: grav(ix,levs)   ! (m/s2)
      real, intent(in)    :: adt(ix,levs)    ! input  temp at dt=0
      real, intent(in)    :: dtp             ! time step in second
      real, intent(inout) :: q(ix,levs,ntrac)   ! input output tracer
      real, intent(out)   :: n1(ix,levs)     ! number density of o (/cm3)
      real, intent(out)   :: n2(ix,levs)     ! number density of o2 (/cm3)
      real, intent(out)   :: ozn(ix,levs)    ! number density of o2 (/cm3)
      real, intent(out)   :: n3(ix,levs)     ! number density of n2 (/cm3)
      real, intent(out)   :: n(ix,levs)      ! total number density (/cm3)
      real, intent(out)   :: rho(ix,levs)    ! density of  (kg/m3)
      real, intent(out)   :: am(ix,levs)     ! avg mass of mix  (kg)
      real, intent(out)   :: am29(ix,levs)   ! avg mass of air from 28.84 => 16.
! local argument
      real ::  dq1(ix,levs,ntrac_i),dq2(ix,levs,ntrac_i)
      real ::  mh2o,mo3, mo,mo2,mn2        
      real ::  qin(ix,levs,ntrac_i)   !, qsumo(ix,levs)
      real ::  q3(ix,levs), qsumo3
      integer i,k,in
! 
       real :: Jrates_O2(ix, levs)            !  replacement of  jj-2014 
       real :: Jrates_O3(ix, levs)            !  O3-photorates
!   
      do in=1,ntrac_i
        do i=1,im
          do k=1,levs
            qin(i,k,in)=max(q(i,k,ntrac-ntrac_i+in), 1.e-36)

          enddo
        enddo
      enddo

! change unit from g/mol to kg
!      mo=amo*1.e-3/avgd
!      mo2=amo2*1.e-3/avgd
!      mn2=amn2*1.e-3/avgd
!      mh2o=amh2o*1.e-3/avgd
!      mo3=amo3*1.e-3/avgd
! at layer , here n,n1,n2 unit is /m3 , rho is in kg/m3
      do i=1,im
        do k=1,levs
!
!          am(i,k)=1./(qin(i,k,1)/mo+qin(i,k,2)/mo2+q(i,k,1)/mh2o+       
!     & q(i,k,2)/mo3+(1.-qin(i,k,1)-qin(i,k,2)-qsumo(i,k))/mn2)
!
          q3(i,k)=    max(q(i,k, 2), 1.e-36)                        !       ozone
          qsumo3= max(q(i,k, 1), 1.e-36)+q3(i,k)                    !       O3+H2O
        am(i,k)=1./(qin(i,k,1)*rmo+qin(i,k,2)*rmo2+q(i,k,1)*rmh2o+       
     &    q(i,k,2)*rmo3+(1.-qin(i,k,1)-qin(i,k,2)-qsumo3)*rmn2)
!
! am = mu29*1.e-3/avgd
!
          n(i,k)=rbz*prsl(i,k)/adt(i,k)
          rho(i,k)=am(i,k)*n(i,k)
          n1(i,k)=qin(i,k,1)*rho(i,k)*rmo
          n2(i,k)=qin(i,k,2)*rho(i,k)*rmo2
          ozn(i,k)= q(i,k,2)*rho(i,k)*rmo3 
      
          n3(i,k) = n(i,k) -n1(i,k) -n2(i,k) -ozn(i,k)
        enddo
      enddo
       am29 = am * 1.e3 * avgd
!       if ( me == 0) print *, maxval(am29), minval(am29), 'VAY29-am'
!======================================================================
!     Molecular "only" diffusion of major species
!======================================================================   
! 
      call idea_tracer_m(im,ix,levs,ntrac_i,grav,prsi,prsl,adt,dtp,     
     &     qin,am,dq1)
!
! do we need sequential update for n1 & n2 & n3 after Molecular diffusion ?
! order Chemistry=> Difusion or Diffusion => Chemistry
!
!
! chemistry O-O2 with updated  Jrates_O2 that replaces JJ-1D
!
      call idea_dissociation_jo2(im,ix,levs,adt,cospass,n1,n2, ozn, n3,
     &          dayno,zg,grav, f107, f107d, Jrates_O2)

!======================================================================
!     Oxygen chemistry O-O2
!======================================================================   
      call idea_tracer_c(im,ix,levs,ntrac_i,adt,dtp,Jrates_O2,
     &     n1,n2,n,rho, qin,dq2)

!=====================================================================
! Research chemical mechanisms & eddy diffusion
!=====================================================================
! Oxygen chemistry O-O2-O3 with prescribed H-OH-HO2
!=====================================================================
!
!      call idea_dissociation_jo3(im,ix,levs,adt,cospass,n1,n2, ozn, n3,
!     &          dayno,zg,grav, f107, f107d, Jrates_O2, Jrates_O3)
!
!      call idea_tracer_cozn(im,ix,levs,ntrac_i,adt,dtp, Jrates_O2,
!     &     Jrates_O3, n1,n2, ozn, n,rho, qin,dq2, q3)
!
!=====================================================================
! Oxygen chemistry O-O2-O3 with prescribed H-OH-HO2 and turbulent MLT
!        diffusion
!=====================================================================
!      call idea_tracer_cdif(im,ix,levs,ntrac_i,adt,dtp,Jrates_O2,
!     &     Jrates_O3, n1,n2,ozn, n,rho, qin,dq2, q3,
!     &     grav,prsi,prsl, am)
!=====================================================================
!  upgrades only for mmr of [O & O2] not for n1 & n2
      do in=1,ntrac_i
        do i=1,im
          do k=1,levs
            q(i,k,in+ntrac-ntrac_i)=q(i,k,in+ntrac-ntrac_i)+            
     &        dq1(i,k,in)  +dq2(i,k,in)
            q(i,k,in+ntrac-ntrac_i)=max(q(i,k,in+ntrac-ntrac_i),1.e-36)            
!
! here upgrades for n1 & n2 needs to be done
!
          enddo
        enddo
      enddo
!======================================================================
! 
!     Here we need to update "n1, n2, ozn, amu, rho, n3"
!
!     nair = P/kT =  constant
!======================================================================
      do i=1,im
        do k=1,levs
!
! update am
!
        q(i,k,2) = max(q3(i,k), 1.e-36)  ! ozone
        qsumo3 = q(i,k,2) + q(i,k,1) 
       am(i,k)=1./(q(i,k,4)*rmo+q(i,k,5)*rmo2+q(i,k,1)*rmh2o+       
     & q(i,k,2)*rmo3+(1.-q(i,k,4)-q(i,k,5)-qsumo3)*rmn2)
!
! new..... n1, n2, n3, ozn, rho for Radiation
!
        rho(i,k)=am(i,k)*n(i,k)
        n1(i,k) =q(i,k,4)*rho(i,k)*rmo
        n2(i,k) =q(i,k,5)*rho(i,k)*rmo2
        ozn(i,k) =q(i,k,2)*rho(i,k)*rmo3
        n3(i,k) =n(i,k)-n1(i,k)-n2(i,k)
        
        enddo
      enddo
         am29 = am * 1.e3 * avgd
!         if ( me == 0) print *, maxval(am29), minval(am29), 'VAY-am29C'
      return
      end
!====================================================================
! BK, Diffusion of Major Species (O-O2-N2), RAA-2006/14 version 
!====================================================================
      subroutine idea_tracer_m(im,ix,levs,ntrac_i,grav,prsi,prsl,adt,   
     &dtp,qin,am,dq)
!-----------------------------------------------------------------------
!
! calaulate tracer changes caused by molecular diffusion
!
!-----------------------------------------------------------------------
      use physcons,  only :rgas=>con_rgas,            
     &               avgd => con_avgd
      use machine, only : kind_phys
      use idea_composition, only:  amo, amo2, amn2, bz, rbz
      implicit none
! Argument
      integer, intent(in) :: im    ! number of data points in up,dudt(first dim)
      integer, intent(in) :: ix    ! max data points in fields
      integer, intent(in) :: levs  ! number of pressure levels
      integer, intent(in) :: ntrac_i ! number of tracer add by IDEA
      real,    intent(in) :: dtp   ! time step in second
      real, intent(in)    :: prsi(ix,levs+1) ! interface pressure in KPa
      real, intent(in)    :: prsl(ix,levs)   ! layer pressure in KPa
      real, intent(in)    :: grav(ix,levs)   ! (m/s2)
      real, intent(in) :: adt(ix,levs)   ! input  temp at dt=0
      real, intent(in) :: qin(ix,levs,ntrac_i)   ! input tracer
      real, intent(in)   :: am(ix,levs)   ! avg mass of mix  (kg)
      real, intent(out):: dq(ix,levs,ntrac_i) ! output tracer changes
!local  variables
      real n1_i(levs+1),n2_i(levs+1),n3_i(levs+1),n_i(levs+1)
      real t_i(levs+1),am_i(levs+1),qout(ix,levs,ntrac_i)
      real beta(2,2,levs+1),a(2,2,levs),b(2,2,levs),c(2,2,levs)
      real ggg(2,2),ee(2,2,levs+1),f(2,levs+1),                         
     &     d12,d13,d23,a12,a13,a23,s12,s13,s23,mo,mo2,mn2,              
     &     dp1(levs),dp1_i(levs+1)
      real partb_i(levs+1),parta(levs),hold1,dtp1,hold2
      integer k,i,kk,kk1,in
! change unit from g/mol to kg
      mo=amo*1.e-3/avgd
      mo2=amo2*1.e-3/avgd
      mn2=amn2*1.e-3/avgd
! some constants
      a12=9.69e18
      a13=9.69e18
      a23=8.3e18
c
      s12=0.774
      s13=0.774
      s23=0.724
! set boundary
      beta(1:2,1:2,1)=0.
      beta(1:2,1:2,levs+1)=0.
      a(1:2,1:2,1)=0.
      c(1:2,1:2,levs)=0.
      ee(1:2,1:2,levs+1)=0.
      f(1:2,levs+1)=0.
!
      dtp1=1./dtp
      t_i=0.
      am_i=0.
      n_i=0.
      n1_i=0.
      n2_i=0.
      n3_i=0.
!
! for each longitude
!
      do i=1,im
! calculate temp in interface pressure levels
! get compositions at interface pressure levels
        do k=2,levs
          t_i(k)=(adt(i,k-1)+adt(i,k))*.5
          am_i(k)=.5*(am(i,k-1)+am(i,k))
          n_i(k)=prsi(i,k)/bz/t_i(k)             ! get out from 1000 .....of kPa
          n1_i(k)=.5*(qin(i,k,1)+qin(i,k-1,1))*am_i(k)*n_i(k)/mo
          n2_i(k)=.5*(qin(i,k,2)+qin(i,k-1,2))*am_i(k)*n_i(k)/mo2
          n3_i(k)=n_i(k)-n1_i(k)-n2_i(k)
        enddo
       if(i.eq.6) then
! 
       endif
!
! calculate beta at interface pressure
        do k=2,levs
          d12=a12*t_i(k)**(s12)
          d13=a13*t_i(k)**(s13)
          d23=a23*t_i(k)**(s23)
          hold1=1./(n1_i(k)*d23+n2_i(k)*d13+n3_i(k)*d12)
          beta(1,1,k)=hold1*d13*mo*(n1_i(k)*mn2*d23+                    
     &            (n2_i(k)*mo2+n3_i(k)*mn2)*d12)
          beta(2,2,k)=hold1*d23*mo2*(n2_i(k)*mn2*d13+                   
     &            (n1_i(k)*mo+n3_i(k)*mn2)*d12)
          beta(1,2,k)=hold1*d23*mo*n1_i(k)*(mn2*d13-mo2*d12)
          beta(2,1,k)=hold1*d13*mo2*n2_i(k)*(mn2*d23-mo*d12)
!      if(i.eq.6) print*,'www6-beta',i,k,beta(1,1,k),hold1,n1_i(k),     
!    & n2_i(k),n3_i(k),d12,d13,d23,t_i(k),mo2,mn2
        enddo
!      if(i.eq.6) print*,'www6-beta',i,beta(1,1,2:levs)
! solve tridiagonal problem
        do k=1,levs
          dp1(k)=1./(prsi(i,k)-prsi(i,k+1))
          parta(k)=dtp*grav(i,k)*dp1(k)/bz
        enddo
        do k=2,levs
          dp1_i(k)=1./(prsl(i,k-1)-prsl(i,k))
          partb_i(k)=.5*(grav(i,k)+grav(i,k-1))/t_i(k)
        enddo
        do k=2,levs
          hold1=parta(k)*partb_i(k)
          hold2=am(i,k-1)*prsl(i,k-1)*dp1_i(k)
          a(1,1,k)=hold1*beta(1,1,k)*(hold2/mo-.5)
          a(1,2,k)=hold1*beta(1,2,k)*(hold2/mo2-.5)
          a(2,1,k)=hold1*beta(2,1,k)*(hold2/mo-.5)
          a(2,2,k)=hold1*beta(2,2,k)*(hold2/mo2-.5)
         enddo
!      print*,'www6-a',i,a(1:2,1:2,levs-3:levs)
        do k=1,levs-1
          hold1=parta(k)*partb_i(k+1)
          hold2=am(i,k+1)*prsl(i,k+1)*dp1_i(k+1)
          c(1,1,k)=hold1*beta(1,1,k+1)*(hold2/mo+.5)
          c(1,2,k)=hold1*beta(1,2,k+1)*(hold2/mo2+.5)
          c(2,1,k)=hold1*beta(2,1,k+1)*(hold2/mo+.5)
          c(2,2,k)=hold1*beta(2,2,k+1)*(hold2/mo2+.5)
         enddo
        do k=2,levs-1
          hold1=am(i,k)*prsl(i,k)*dp1_i(k+1)
          hold2=am(i,k)*prsl(i,k)*dp1_i(k)
      b(1,1,k)=1.+parta(k)*(partb_i(k+1)*beta(1,1,k+1)*(hold1/mo-.5)    
     &                    +partb_i(k)*beta(1,1,k)*(hold2/mo+.5))
      b(2,2,k)=1.+parta(k)*(partb_i(k+1)*beta(2,2,k+1)*(hold1/mo2-.5)   
     &                    +partb_i(k)*beta(2,2,k)*(hold2/mo2+.5))
      b(1,2,k)=parta(k)*(partb_i(k+1)*beta(1,2,k+1)*(hold1/mo2-.5)      
     &                    +partb_i(k)*beta(1,2,k)*(hold2/mo2+.5))
      b(2,1,k)=parta(k)*(partb_i(k+1)*beta(2,1,k+1)*(hold1/mo-.5)       
     &                    +partb_i(k)*beta(2,1,k)*(hold2/mo+.5))
        enddo
          hold1=am(i,1)*prsl(i,1)*dp1_i(2)
      b(1,1,1)=1.+parta(1)*partb_i(2)*beta(1,1,2)*(hold1/mo-.5)
      b(2,2,1)=1.+parta(1)*partb_i(2)*beta(2,2,2)*(hold1/mo2-.5)
      b(1,2,1)=parta(1)*partb_i(2)*beta(1,2,2)*(hold1/mo2-.5)
      b(2,1,1)=parta(1)*partb_i(2)*beta(2,1,2)*(hold1/mo-.5)
          hold2=am(i,levs)*prsl(i,levs)*dp1_i(levs)
      b(1,1,levs)=1.+parta(levs)*partb_i(levs)*beta(1,1,levs)*          
     &(hold2/mo+.5)
      b(2,2,levs)=1.+parta(levs)*partb_i(levs)*beta(2,2,levs)*          
     &(hold2/mo2+.5)
      b(1,2,levs)=parta(levs)*partb_i(levs)*beta(1,2,levs)*             
     &(hold2/mo2+.5)
      b(2,1,levs)=parta(levs)*partb_i(levs)*beta(2,1,levs)*             
     &(hold2/mo+.5)
       do k=levs,1,-1
         ggg(1,1)=b(2,2,k)-c(2,1,k)*ee(1,2,k+1)-c(2,2,k)*ee(2,2,k+1)
         ggg(2,2)=b(1,1,k)-c(1,1,k)*ee(1,1,k+1)-c(1,2,k)*ee(2,1,k+1)
         ggg(1,2)=-1.*b(1,2,k)+c(1,1,k)*ee(1,2,k+1)+c(1,2,k)*ee(2,2,k+1)
         ggg(2,1)=-1.*b(2,1,k)+c(2,1,k)*ee(1,1,k+1)+c(2,2,k)*ee(2,1,k+1)
         hold1=1./(ggg(1,1)*ggg(2,2)-ggg(1,2)*ggg(2,1))
         ggg=ggg*hold1
         ee(1,1,k)=ggg(1,1)*a(1,1,k)+ggg(1,2)*a(2,1,k)       
         ee(1,2,k)=ggg(1,1)*a(1,2,k)+ggg(1,2)*a(2,2,k)       
         ee(2,1,k)=ggg(2,1)*a(1,1,k)+ggg(2,2)*a(2,1,k)       
         ee(2,2,k)=ggg(2,1)*a(1,2,k)+ggg(2,2)*a(2,2,k)       
      f(1,k)=ggg(1,1)*(qin(i,k,1)+c(1,1,k)*f(1,k+1)                      
     &+c(1,2,k)*f(2,k+1))+ggg(1,2)*(qin(i,k,2)+c(2,1,k)*f(1,k+1)         
     &+c(2,2,k)*f(2,k+1))
      f(2,k)=ggg(2,1)*(qin(i,k,1)+c(1,1,k)*f(1,k+1)                      
     &+c(1,2,k)*f(2,k+1))+ggg(2,2)*(qin(i,k,2)+c(2,1,k)*f(1,k+1)         
     &+c(2,2,k)*f(2,k+1))
        enddo
        do in=1,ntrac_i
          qout(i,1,in)=f(in,1)
          dq(i,1,in)=qout(i,1,in)-qin(i,1,in)
        enddo
        do k=2,levs
          qout(i,k,1)=ee(1,1,k)*qout(i,k-1,1)+ee(1,2,k)*qout(i,k-1,2)+  
     &              f(1,k)
          qout(i,k,2)=ee(2,1,k)*qout(i,k-1,1)+ee(2,2,k)*qout(i,k-1,2)+  
     &              f(2,k)
          do in=1,ntrac_i
            dq(i,k,in)=qout(i,k,in)-qin(i,k,in)
          enddo
        enddo
      enddo !i
      return
      end subroutine idea_tracer_m


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!   
      subroutine idea_tracer_c(im,ix,levs,ntrac_i,adt,dtp,jo2,n1,n2,     
     &        n,rho,qin,dq)
!-----------------------------------------------------------------------
!
! calculate tracer changes caused by chemistry reaction
!    [h2o o3 n2 o o2 h oh ho2 no n co2 ndens tg mu zkm]  
!-----------------------------------------------------------------------
      use physcons, only : rgas=>con_rgas
      use physcons, only : avgd => con_avgd
      use machine,  only : kind_phys
      use idea_composition, only : amo, amn2, amo2
!      use idea_tracers_input,only : mo, mo2, mn2
      use idea_tracer_mod,  only : oh, ho2, vmr_glob
!
!==========================================================   
!         1          5    6  7  8
!"names: h2o o3 n2 o o2  h oh ho2  no n co2 ndens tg mu zkm"
!         1/m3
!==========================================================
!
      implicit none
!
! Argument
      integer, intent(in) :: im          ! number of data points in up,dudt(first dim)
      integer, intent(in) :: ix          ! max data points in fields
      integer, intent(in) :: levs        ! number of pressure levels
      integer, intent(in) :: ntrac_i     ! number of tracer add by IDEA
      real,    intent(in) :: dtp         ! time step in second
      real, intent(in) :: adt(ix,levs)   ! input  temp at dt=0

      real, intent(in) :: qin(ix,levs,ntrac_i)   ! input tracer

!!2013      real, intent(in) :: jj(levs)  ! input photo diss rate
      real, intent(in) :: jo2(ix, levs)   ! input photo diss rate
      real, intent(in) :: n1(ix,levs)     ! number density of o
      real, intent(in) :: n2(ix,levs)     ! number density of o2
      real, intent(in) :: n(ix,levs)      ! number density of mixture
      real, intent(in) :: rho(ix,levs)    ! density of mixture
      real, intent(out):: dq(ix,levs,ntrac_i) ! output
! Local variables
      real, dimension(levs) :: noh, nh, nho2, natom, ndens
!
      real k1,k2,p1,p2,L1,L2
      real :: k3, k4, k5
      real :: mo,mo2,mn2
      real :: qout(ix,levs,ntrac_i)
      integer k,i
      real :: q1new, q2new
      real :: koom, koomt
!
! [h2o o3 n2 o o2 h oh ho2 no n co2 ndens tg mu zkm]  
!   1           5             10      12        15
!
!      nh(1:levs)   = vmr_glob(1:levs,6)
!      noh(1:levs)  = vmr_glob(1:levs,7)
!      nho2(1:levs)  = vmr_glob(1:levs,8)
!
!
! fill in noh, nh, nho2, natom in 1/m3 or vmr
!
      mo=amo*1.e-3/avgd
      mo2=amo2*1.e-3/avgd
      mn2=amn2*1.e-3/avgd
!
        k3=4.2e-17
        k4=1.0e-26
        k5=3.5e-17
        k2=0.
        koom =6.e-34*1.e-12  ! 1/m3/m3 n*n1
!-------------------------      should be in tracer_init
      do k=1,levs
      do i=1,im

!========================================================================
! O+O+M  =k1 => O2 +M
! O2 +hv =JJ => O + O
! O + O+ M  =k1 => O2+M 
! O + O2    =k2=>  O3
! O + O     =k4=>  O2
!
! O +OH  =k3=> O2 +H          rate = 1.80E-11*exp(    180./t)   (165      
! O +HO2 =grk5=> O2 +OH       rate = 3.00E-11*exp(    200./t)
!
!  n1 = sqrt(J2*[n2/N]/k1) .......with J3
!  n1 = sqrt(J2/k1*J3/Ko)/sqrt(N)
!===================================================
! get coefficent array o o2 n2
!----------------------------------------------------
!        k1=2.76e-46*exp(710./adt(i,k))             ! JPL-16
        k1=4.7e-45*(300./adt(i,k))**2               ! OLD-WAM & BK-73
        koomt = koom*(300./adt(i,k))**2.4
!
!
        p1=2.*Jo2(i,k)*n2(i,k) * mo/rho(i,k)        !O-production mmr/density*
!
       L1=2.*k1*n1(i,k)*n(i,k)+k2*n2(i,k)*n(i,k)    ! O-loss
     &   +k3*oh(k)+k5*ho2(k)                        ! VAY-2017 Noh & Nho2
!
       L1 = L1+ koomt*n(i,k)*n2(i,k)

        p2= (k1*n1(i,k)*n1(i,k)*n(i,k))
     &       * mo2/rho(i,k)  
!
        L2 = koomt*n(i,k)*n1(i,k) + Jo2(i,k)
!
        q1new=(qin(i,k,1)+p1*dtp)/(1.+L1*dtp)
        q2new=(qin(i,k,2)+p2*dtp)/(1.+L2*dtp)

        dq(i,k,1)=q1new-qin(i,k,1)
        dq(i,k,2)= -dq(i,k,1)     ! TFR-version loss O  =  gain O2 & vise versa, not important
!        dq(i,k,2)=q2new-qin(i,k,2)        
! 
        enddo
      enddo

      RETURN
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! wam-2014, below
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      do k=1,levs
      do i=1,im
! get coefficent array o o2 n2
        k1=4.7e-45*(300./adt(i,k))**2
        k2=6.e-46*(300./adt(i,k))**(2.4)
        p1=2.*jo2(i,k)*n2(i,k)*mo/rho(i,k)
        p2=k1*n1(i,k)**2*n(i,k)*mo2/rho(i,k)
        L1=2.*k1*n1(i,k)*n(i,k)+k2*n2(i,k)*n(i,k)
        L2=k2*n1(i,k)*n(i,k)+jo2(i,k)
        qout(i,k,1)=(qin(i,k,1)+p1*dtp)/(1.+L1*dtp)
        qout(i,k,2)=(qin(i,k,2)+p2*dtp)/(1.+L2*dtp)
        dq(i,k,1)=qout(i,k,1)-qin(i,k,1)
        dq(i,k,2)=qout(i,k,2)-qin(i,k,2)
      enddo
      enddo
      return
      end subroutine idea_tracer_c


!-----------------------------------------------------------------------
!
! Update Ox-chemistry only above the certain level where JO2 & JO3 
!      can be estimated + below make sure that 
!             "no"-updates for O-O2-O3
!
!------------------------------------------------------------------------
      subroutine idea_tracer_cozn(im,ix,levs,ntrac_i,adt,dtp,
     & jo2,jo3, n1,n2,n3, n,rho,qin,dq, q3)
!-----------------------------------------------------------------------
!
! calaulate tracer changes caused by chemistry reaction
!
!-----------------------------------------------------------------------
      use physcons, only : rgas=>con_rgas
      use physcons, only : avgd => con_avgd
      use machine,  only : kind_phys
      use idea_solar      , only : nps
      use idea_composition, only : amo, amn2, amo2, amo3
      use idea_tracer_mod,  only : oh, ho2, vmr_glob     ! vmr_glob GLOBAL profiles of species in 1/m3
!      use idea_tracers_input,only : mo, mo2, mn2
      use idea_tracer_mod,  only : oh, ho2
      implicit none
!
! Argument
      integer, intent(in) :: im          ! number of data points in up,dudt(first dim)
      integer, intent(in) :: ix          ! max data points in fields
      integer, intent(in) :: levs        ! number of pressure levels
      integer, intent(in) :: ntrac_i     ! number of tracer add by IDEA
      real,    intent(in) :: dtp         ! time step in second
      real, intent(in) :: adt(ix,levs)   ! input  temp at dt=0

      real, intent(in) :: qin(ix,levs,ntrac_i)   ! input tracer

!!2014      real, intent(in) :: jj(levs)         ! input photo diss rate
      real, intent(in) :: jo2(ix, levs)          ! input photo diss rate
      real, intent(in) :: jo3(ix, levs)          ! input photo diss rate
      real, intent(in) :: n1(ix,levs)            ! number density of o
      real, intent(in) :: n2(ix,levs)            ! number density of o2
      real, intent(in) :: n3(ix,levs)            ! number density of o2
      real, intent(in) :: n(ix,levs)             ! number density of mixture
      real, intent(in) :: rho(ix,levs)           ! density of mixture
!
!
      real, intent(out):: dq(ix,levs,ntrac_i) ! output
      real, intent(inout):: q3(ix,levs)
!
! Local variables
      integer            ::    k,i
      real :: dq3
      real :: p1,p2,p3
      real :: L1,L2, L3
      real :: k1,k2,k3, k4, k5, ko3o, ko3h
      real :: mo,mo2,mn2, mo3
      real, dimension(levs) :: noh, nh, nho2, natom, ndens
!      real :: qout(ix,levs,ntrac_i)     
      real :: q1new, q2new
      real :: koom, koomt
      real :: rt1            ! 1/T
!
      mo=amo*1.e-3/avgd
      mo2=2.*mo              !amo2*1.e-3/avgd
      mo3=3.*mo 
      mn2=amn2*1.e-3/avgd

!
        k3=4.2e-17           ! O+OH  2.2e-17*exp(120/T)
        k4=1.0e-26           ! O+O = O2 artificila LOSS of "O"     ????
        k5=3.5e-17           ! O+HO2 3.5e-17 + 7.3E-17 +3.3e-17*exp(200/T)
        k2=0.
        koom =6.e-34*1.e-12  ! O+O+M   1/m3/m3 n*n1
!-------------------------      should be in tracer_init or "idea_chem_init"
!
! [h2o o3 n2 o o2 h oh ho2 no n co2 ndens tg mu zkm]  
!   1           5             10      12        15
!
      nh(1:levs)   = vmr_glob(1:levs,6)
      noh(1:levs)  = vmr_glob(1:levs,7)
      nho2(1:levs)  = vmr_glob(1:levs,8)
!
!fake don't use for extra-checks
!
      natom(1:levs)  = vmr_glob(1:levs,4)
      ndens(1:levs)  = vmr_glob(1:levs,12)
!
! fill in noh, nh, nho2, natom in 1/m3 or vmr                            
      do i=1,im 
         do k=nps,levs
          rt1 =1./adt(i,k)
!-----------------------------------------------------------------------
! O3 + hV =J3=> O2 +O
! O3  + O =Ko=>  2*O2 
! O3  + H =Kh=> OH +O2 
! O + O2 + M ->koom-> O3 + M
!        rxt(:,k,usr_O_O2_ndx) = 6.e-34_r8 * (300./tp(:))**2.4_r8
![O_O3,cph=392.19]     O + O3 -> 2*O2    ; 8.00e-12, -2060.
![H_O3,cph=194.71]   H + O3 -> OH + O2   ; 1.40e-10,   -470.
! O3_aprox_day =   n1*n2*M*koom/J3
! O3_aprox_night = n1*n2*M*koom/[Kh*H+Ko*O]
!   Can be used as correction of O3 at Z > 50 km, needs J3 & [H]-empirical
!
!  See /scratch3/NCEPDEV/swpc/save/Valery.Yudin/NCAR/CESM3_b13/chemistry/pp_waccm_mozart
!
!========================================================================
! O+O+M  =k1 => O2 +M
! O2 +hv =JJ => O + O
! O + O+ M  =k1  => O2+M 
! O + O2+M  =k2  =>  O3
! O + O3    =ko3o=> 2O2
!
! O +OH  =k3=> O2 +H          rate = 1.80E-11*exp(    180./t)   (165      
! O +HO2 =grk5=> O2 +OH       rate = 3.00E-11*exp(    200./t)
!
!  n1 = sqrt(J2*[n2/N]/k1) .......with J3
!  n1 = sqrt(J2/k1*J3/Ko)/sqrt(N)
!
!----------------------------------------------------    
!JPL        k1=2.76e-46*exp(710.*rt1)
        k1    = 4.70e-45*(300.*rt1)**2.                         ! BS & BK
        koomt = koom*(300.*rt1)**2.4                            ! O+O+M    AKS-2005
        ko3o  = 8.0e-18*exp(-2060.*rt1)                         ! x 1e-6  1/m3 vs 1/cm3 AKS-2005
        ko3h  = 1.4e-16*exp(- 270.*rt1)                         ! x 1e-6    AKS-2005
!
!
        p1=(2.*Jo2(i,k)*n2(i,k)+Jo3(i,k)*n3(i,k))
     &       * mo/rho(i,k)                                      !O-production mmr = rhoi/density*
        p2= (k1*n1(i,k)*n1(i,k)*n(i,k)+Jo3(i,k)*n3(i,k))
     &       * mo2/rho(i,k)
  
!
       L2 = koomt*n(i,k)*n1(i,k) + Jo2(i,k)
       L1=  2.*k1*n1(i,k)*n(i,k)                                ! +k2*n2(i,k)*n(i,k)    ! O-loss
     &      +k3*Noh(k)+k5*Nho2(k)                               ! + O-loss by HOx
!
       L1 = L1+ koomt*n(i,k)*n2(i,k) +ko3o*n3(i,k)*0.           ! + O3-related losses O+O2+M=> O3 O3+O=>
!
!        q1new=(qin(i,k,1)+p1*dtp)/(1.+L1*dtp)
!        dq(i,k,1)=  q1new-qin(i,k,1)
!        dq(i,k,2)= -dq(i,k,1)

        q2new=(qin(i,k,2)+p2*dtp)/(1.+L2*dtp)
        dq(i,k,2)=q2new-qin(i,k,2) 
        dq(i,k,1)= -dq(i,k,2) 
                                                                 ! "loss O"  =  "gain O2" & vise versa       
!
!
! update ozone with new n2 & n1
!
        L3 = Jo3(i,k) + ko3o*n1(i,k)+ko3h*Nh(k)                 !  +ko3h*NH 
        P3 = koomt*n1(i,k)*n2(i,k)*n(i,k)*mo3/rho(i,k)          !  production only "oxygen" cycle
        dq3 =dtp*p3/(1.+L3*dtp) 
        q3(i,k)= q3(i,k)+ dq3                                   !  bkg from NRL "+" P3/L3
!
!        dq(i,k,2)= -dq(i,k,1)- 2./3.*dq3                       !  O3+O = 2O2 =>  dq1 +2./3.*dq3 = -dq2
        enddo
!=================================================================================
! Add eddy + Molecular Diffusion diffusion for "2"-tracers Ox=O3+O & O2
!
! above ~80 km let's make NRL-scheme => zero Production and Loss of O3
!
!     and elaborate O3 from the Zonal Mean NRL Chem-ry  to "Diurnal Variations"
!
! Add for each tracer O3-O2-O Molec+Eddy Diffusion !
!     Check Jo3/Jo2 in standalone.... Extend J's => surface
!=================================================================================
      enddo

      RETURN

      end subroutine idea_tracer_cozn


