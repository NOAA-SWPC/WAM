
                        module module_control
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
use mpi
!
use module_kinds
use module_exchange
use module_constants
use module_derived_types,only: bc_h_all,bc_v_all
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer(kind=kint),parameter :: num_domains_max=99
!
!-----------------------------------------------------------------------
!
      public NMMB_Finalize
      public grid_consts
      public num_domains_max
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!---look-up tables------------------------------------------------------
!-----------------------------------------------------------------------
integer(kind=kint),parameter :: &
 itb=201 &                   ! convection tables, dimension 1
,jtb=601 &                   ! convection tables, dimension 2
,kexm=10001 &                ! size of exponentiation table
,kztm=10001                  ! size of surface layer stability function table

real(kind=kfpt),parameter :: &
 ph=105000. &                ! upper bound of pressure range
,thh=350. &                  ! upper bound of potential temperature range
,thl=200.                    ! upper bound of potential temperature range

integer(kind=kint):: &
 kexm2 &                     ! internal points of exponentiation table
,kztm2                       ! internal pts. of the stability function table

real(kind=kfpt) :: &
 dex &                       ! exponentiation table step
,dzeta1 &                    ! sea table z/L step
,dzeta2 &                    ! land table z/L step 
,fh01 &                      ! prandtl number for sea stability functions
,fh02 &                      ! prandtl number for land stability functions
,pl &                        ! lower bound of pressure range
,rdp &                       ! scaling factor for pressure
,rdq &                       ! scaling factor for humidity
,rdth &                      ! scaling factor for potential temperature
,rdthe &                     ! scaling factor for equivalent pot. temperature
,rdex &                      ! exponentiation table scaling factor
,xmax &                      ! upper bound for exponent in the table
,xmin &                      ! lower bound for exponent in the table
,ztmax1 &                    ! upper bound for z/L for sea stab. functions
,ztmin1 &                    ! lower bound for z/L for sea stab. functions
,ztmax2 &                    ! upper bound for z/L for land stab. functions
,ztmin2                      ! lower bound for z/L for land stab. functions
 
real(kind=kfpt),dimension(1:itb):: &
 sthe &                      ! range for equivalent potential temperature
,the0                        ! base for equivalent potential temperature           

real(kind=kfpt),dimension(1:jtb):: &
 qs0 &                       ! base for saturation specific humidity
,sqs                         ! range for saturation specific humidity

real(kind=kfpt),dimension(1:kexm):: &
 expf                        ! exponentiation table

real(kind=kfpt),dimension(1:kztm):: &
 psih1 &                     ! sea heat stability function
,psim1 &                     ! sea momentum stability function
,psih2 &                     ! land heat stability function
,psim2                       ! land momentum stability function

real(kind=kfpt),dimension(1:itb,1:jtb):: &
 ptbl                        ! saturation pressure table

real(kind=kfpt),dimension(1:jtb,1:itb):: &
 ttbl                        ! temperature table
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!---miscellaneous control parameters------------------------------------
!-----------------------------------------------------------------------
character(64):: &
 infile
 
logical(kind=klog):: &
 readbc &                    ! read regional boundary conditions
,runbc                       ! boundary data ready, start run

integer(kind=kint):: &
 ihr &                       ! current forecast hour
,ihrbc &                     ! boundary condition hour
,ihrstbc &                   ! boundary conditions starting time
,nbc                         ! boundary data logical unit

integer(kind=kint),dimension(1:3):: &
 idatbc(3)                   ! date of boundary data, day, month, year

real(kind=kfpt):: &
 bofac                       ! amplification of diffusion along bndrs.

integer(kind=kint):: &
 lnsbc                       ! # of boundary lines with enhanced diffusion
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
       contains
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
      subroutine consts &
      (global &
      ,dt &
      ,smag2,codamp,wcor &
      ,pt &
      ,tph0d,tlm0d &
      ,sbd,wbd &
      ,dphd,dlmd &
      ,dxh,rdxh &
      ,dxv,rdxv &
      ,dyh,rdyh &
      ,dyv,rdyv &
      ,ddv,rddv &
      ,ddmpu,ddmpv &
      ,ef4t,wpdar &
      ,fcp,fdiv &
      ,curv,f &
      ,fad,fah &
      ,dare,rare &
      ,glat,glon &
      ,glat_sw,glon_sw &
      ,vlat,vlon &
      ,hdacx,hdacy &
      ,hdacvx,hdacvy &
      ,lnsh,lnsad &
      ,adv_standard,adv_upstream &
      ,e_bdy,n_bdy,s_bdy,w_bdy &
      ,nboco,tboco &
      ,my_domain_id,mype &
      ,its,ite,jts,jte &
      ,ims,ime,jms,jme &
      ,ids,ide,jds,jde )
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
integer(kind=kint),intent(in) :: &
 ids,ide,jds,jde &
,ims,ime,jms,jme &
,its,ite,jts,jte &
,lnsh &
,my_domain_id &
,mype

integer(kind=kint),intent(out) :: &
 lnsad &
,nboco

real(kind=kfpt),intent(in) :: &
 codamp &    ! divergence damping coefficient
,dt &        ! fundamental dynamics timestep (sec)
,dlmd &      ! grid increment, delta lambda, degrees
,dphd &      ! grid increment, delta phi, degrees
,pt &        ! Pressure at top of domain (Pa)
,sbd &       ! degrees from center of domain to southern boundary
,smag2 &     ! Smagorinsky coefficient for 2nd order diffusion 
,tlm0d &
,tph0d &
,wbd &       ! degrees from center of domain to western boundary
,wcor

real(kind=kfpt),intent(out) :: &
 ddmpv &
,dyh &
,dyv &
,ef4t &
,glat_sw &   ! geographic latitude (radians) of domain's SW corner
,glon_sw &   ! geographic longitude (radians) of domain's SW corner (positive east)
,rdyh &
,rdyv &
,tboco
 
real(kind=kfpt),dimension(jds:jde),intent(out) :: &
 curv &
,dare &
,ddmpu &
,ddv &
,dxh &
,dxv &
,fad &
,fah &
,fcp &
,fdiv &
,rare &
,rddv &
,rdxh &
,rdxv &
,wpdar
 
real(kind=kfpt),dimension(ims:ime,jms:jme),intent(out) :: &
 f &
,glat &
,glon &
,vlat &
,vlon &
,hdacx &
,hdacy &
,hdacvx &
,hdacvy
 
logical(kind=klog),intent(in) :: &
 global      ! global forecast if true
 
logical(kind=klog),intent(out) :: &
 adv_standard &              ! is task in standard advec region?
,adv_upstream &              ! is task in upstream advec region?
,e_bdy &                     ! is task on domain's eastern boundary?
,n_bdy &                     ! is task on domain's northern boundary?
,s_bdy &                     ! is task on domain's southern boundary?
,w_bdy                       ! is task on domain's western boundary?
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      if(global) then
!-----------------------------------------------------------------------
!---global branch-------------------------------------------------------
!-----------------------------------------------------------------------
!
        lnsbc=lnsh
        bofac=0.
!
        lnsad=0
!-----------------------------------------------------------------------
      else
!-----------------------------------------------------------------------
!---regional branch-----------------------------------------------------
!-----------------------------------------------------------------------
! 
        lnsbc=lnsh
        bofac=4.
!
        lnsad=lnsh+2
!-----------------------------------------------------------------------
      endif
!-----------------------------------------------------------------------
      adv_upstream=.false.
!!!   if(jts<jds+1+lnsad.or.jte>jde-1-lnsad.or. &
!!!      its<ids+1+lnsad.or.ite>ide-1-lnsad)then
      if(jts<jds+1+lnsad.or.jte>=jde-1-lnsad.or. &
         its<ids+1+lnsad.or.ite>=ide-1-lnsad)then
        adv_upstream=.true.
      endif
!
      adv_standard=.false.
      if(jte>=jds+1+lnsad.and.jts<=jde-1-lnsad.and. &
         ite>=ids+1+lnsad.and.its<=ide-1-lnsad)then
        adv_standard=.true.
      endif
!-----------------------------------------------------------------------
      if(.not.global.and.my_domain_id==1)then
        ihrbc=0
        write(infile,'(a,i4.4)')'boco.',ihrbc
        nbc=18
        open(unit=nbc,file=infile,status='old',form='unformatted')
        read (nbc) runbc,idatbc,ihrstbc,tboco
!	write(0,*) 'runbc: ', runbc
!	write(0,*) 'idatbc: ', idatbc
!	write(0,*) 'ihrstbc: ', ihrstbc
!	write(0,*) 'tboco: ', tboco
        rewind nbc
        close(unit=nbc)
        if(mype==0)then
          write(0,*)'*** Read tboco in CONSTS from ',infile
        endif
        nboco=nint(tboco/dt)
      endif
!
!-----------------------------------------------------------------------
!---derived vertical grid constants-------------------------------------
!-----------------------------------------------------------------------
      ef4t=0.5*dt/cp
!
!-----------------------------------------------------------------------
!-------------grid related arrays---------------------------------------
!-----------------------------------------------------------------------
!
      call grid_consts(global &
      ,dt,smag2,codamp,wcor &
      ,tph0d,tlm0d &
      ,sbd,wbd &
      ,dphd,dlmd &
      ,dxh,rdxh &
      ,dxv,rdxv &
      ,dyh,rdyh &
      ,dyv,rdyv &
      ,ddv,rddv &
      ,ddmpu,ddmpv &
      ,wpdar &
      ,fcp,fdiv &
      ,curv,f &
      ,fad,fah &
      ,dare,rare &
      ,glat,glon &
      ,glat_sw,glon_sw &
      ,vlat,vlon &
      ,hdacx,hdacy &
      ,hdacvx,hdacvy &
      ,e_bdy,n_bdy,s_bdy,w_bdy &
      ,its,ite,jts,jte &
      ,ims,ime,jms,jme &
      ,ids,ide,jds,jde )
!
!-----------------------------------------------------------------------
!-------------look-up tables--------------------------------------------
!-----------------------------------------------------------------------
!
      call exptbl
!
!-----------------------------------------------------------------------
      pl=max(pt,1000.)
      call tablep
      call tablet
!-----------------------------------------------------------------------
      call psitbl
!-----------------------------------------------------------------------
!
                        end subroutine consts       
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
                        subroutine grid_consts &
      ( global &
      ,dt,smag2,codamp,wcor &
      ,tph0d,tlm0d &
      ,sbd,wbd &
      ,dphd,dlmd &
      ,dxh,rdxh &
      ,dxv,rdxv &
      ,dyh,rdyh &
      ,dyv,rdyv &
      ,ddv,rddv &
      ,ddmpu,ddmpv &
      ,wpdar &
      ,fcp,fdiv &
      ,curv,f &
      ,fad,fah &
      ,dare,rare &
      ,glat,glon &
      ,glat_sw,glon_sw &
      ,vlat,vlon &
      ,hdacx,hdacy &
      ,hdacvx,hdacvy &
      ,e_bdy,n_bdy,s_bdy,w_bdy &
      ,its,ite,jts,jte &
      ,ims,ime,jms,jme &
      ,ids,ide,jds,jde )
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
integer(kind=kint),intent(in) :: &
 its &   ! starting integration index in i
,ite &   ! ending integration index in i
,ims &   ! starting memory index in i
,ime &   ! ending memory index in i
,ids &   ! starting domain index in i
,ide &   ! ending domain index in i
,jts &   ! starting integration index in j
,jte &   ! ending integration index in j
,jms &   ! starting memory index in j
,jme &   ! ending memory index in j
,jds &   ! starting domain index in j
,jde     ! ending domain index in j

real(kind=kfpt),intent(in) :: &
 codamp &    ! divergence damping coefficient
,dt &        ! fundamental dynamics timestep (sec)
,sbd &       ! degrees from center of domain to southern boundary
,smag2 &     ! Smagorinsky coefficient for 2nd order diffusion 
,tlm0d &
,tph0d &
,wbd &       ! degrees from center of domain to western boundary
,dlmd &      ! grid increment, delta lambda, degrees
,dphd &      ! grid increment, delta phi, degrees
,wcor

real(kind=kfpt),intent(out) :: &
 ddmpv &
,dyh &
,dyv &
,glat_sw &   ! geographic latitude (radians) of domain's SW corner
,glon_sw &   ! geographic longitude (radians) of domain's SW corner (positive east)
,rdyh &
,rdyv
 
real(kind=kfpt),dimension(jds:jde),intent(out) :: &
 curv &
,dare &
,ddmpu &
,ddv &
,dxh &
,dxv &
,fad &
,fah &
,fcp &
,fdiv &
,rare &
,rddv &
,rdxh &
,rdxv &
,wpdar
 
real(kind=kfpt),dimension(ims:ime,jms:jme),intent(out) :: &
 f &
,glat &
,glon &
,vlat &
,vlon &
,hdacx &
,hdacy &
,hdacvx &
,hdacvy
 
logical(kind=klog),intent(in) :: &
 global      ! global forecast if true      

logical(kind=klog),intent(inout) :: &
 e_bdy &                     ! is task on domain's eastern boundary?
,n_bdy &                     ! is task on domain's northern boundary?
,s_bdy &                     ! is task on domain's southern boundary?
,w_bdy                       ! is task on domain's western boundary?
       
!
!-----------------------------------------------------------------------
!--local variables------------------------------------------------------
!-----------------------------------------------------------------------
!      

integer(kind=kint):: &
 i &                         ! index in x direction
,i_hi &                      ! max i loop limit (cannot be > ide)
,i_lo &                      ! min i loop limit (cannot be < ids)
,j &                         ! index in y direction
,j_hi &                      ! max j loop limit (cannot be > jde)
,j_lo                        ! min j loop limit (cannot be < jds)
  
real(kind=kfpt),save :: &
 eps=1.e-5

real(kind=kfpt):: &
 acdt &                      ! diffusion coefficient parameter
,alm &                       ! lambda
,anum &                      ! numerator
,aph &                       ! phi
,arg &                       ! 
,ave &                       ! average
,cddamp &                    ! divergence damping factor
,coac &                      ! nonlinear diffusion coefficient
,ctlm &                      ! temporary
,ctph &                      ! temporary
,ctph0 &                     ! cos(tph0)
,denom &                     ! denominator
,dlm &                       ! grid size, delta lambda, radians
,dph &                       ! grid size, delta phi, degrees
,fpole &                     ! factor to modify polar area
,rdlm &                      ! 1 / delta lambda
,rdph &                      ! 1 / delta phi
,relm &                      ! temporary
,sb &                        ! radians from center to southern boundary
,sph &                       ! temporary
,stlm &                      ! temporary
,stph &                      ! temporary
,stph0 &                     ! sin(tph0)
,tlm &                       ! rotated lambda
,tlm_base &                  ! temporary
,tph &                       ! rotated phi
,tph_base &                  ! temporary
,tpv &                       ! rotated phi, v points
,tph0 &                      ! radians grid rotated in phi direction
,wb                          ! radians from center to western boundary

real(kind=kfpt),dimension(jds:jde):: &
 hdacxj &
,hdacyj &
,hdacvxj &
,hdacvyj     

!-----------------------------------------------------------------------
!***  Because subdomains that lie along the global domain boundary
!***  may have haloes that extend beyond the global limits, create
!***  limits here that keep loops from reaching beyond those
!***  global limits.
!-----------------------------------------------------------------------
!
      i_lo=max(ims,ids)
      i_hi=min(ime,ide)
      j_lo=max(jms,jds)
      j_hi=min(jme,jde)
!
!-----------------------------------------------------------------------
!---to be read from a namelist in the future----------------------------
!-----------------------------------------------------------------------
!
      coac=smag2*smag2 !second order
!
!-----------------------------------------------------------------------
!---derived physical constants------------------------------------------
!-----------------------------------------------------------------------
!
      acdt=coac*dt
      cddamp=codamp*dt
!
!-----------------------------------------------------------------------
!***  Each task identifies whether or not it is adjacent to a boundary
!***  of the full domain.
!-----------------------------------------------------------------------
!
      w_bdy=(its==ids)  ! This task is on the western boundary
      e_bdy=(ite==ide)  ! This task is on the eastern boundary
      s_bdy=(jts==jds)  ! This task is on the southern boundary
      n_bdy=(jte==jde)  ! This task is on the northern boundary
!
!-----------------------------------------------------------------------
!--------------derived geometrical constants----------------------------
!-----------------------------------------------------------------------
!
      tph0=tph0d*dtr
      wb=wbd*dtr
      sb=sbd*dtr
      dlm=dlmd*dtr
      dph=dphd*dtr
      rdlm=1./dlm
      rdph=1./dph
!
      stph0=sin(tph0)
      ctph0=cos(tph0)
!-----------------------------------------------------------------------
!---derived horizontal grid constants-----------------------------------
!-----------------------------------------------------------------------
      dyh=a*dph
      dyv=a*dph
      rdyh=1./dyh
      rdyv=1./dyv
!-----------------------------------------------------------------------
      do j=jds,jde
        dxh(j)=0.
        rdxh(j)=0.
        dare(j)=0.
        rare(j)=0.
        wpdar(j)=0.
        fah(j)=0.
        fcp(j)=0.
        fdiv(j)=0.
        dxv(j)=0.
        rdxv(j)=0.
        curv(j)=0.
        fad(j)=0.
        ddmpu(j)=0.
        hdacxj(j)=0.
        hdacyj(j)=0.
        hdacvxj(j)=0.
        hdacvyj(j)=0.
      enddo
!-----------------------------------------------------------------------
!
      global_regional_setup: if(global) then
!
!-----------------------------------------------------------------------
!---global branch-------------------------------------------------------
!-----------------------------------------------------------------------
        tph=sb-dph
        fpole=4.0                                                                      
!-----------------------------------------------------------------------
!----south pole---------------------------------------------------------
!-----------------------------------------------------------------------
        tph=tph+dph
        tpv=tph+dph*0.5                                               
        dxh(jds+1)=0.
        dxv(jds+1)=a*dlm*cos(tpv)
        ddmpv=cddamp*dyv/(2.*dyv)                                                         
!
        rdxh(jds+1)=0.
        rdxv(jds+1)=1./dxv(jds+1)
        dare(jds+1)=dxv(jds+1)*dyv*0.5*fpole !for ghost area
        rare(jds+1)=1./dare(jds+1)
        wpdar(jds+1)=-wcor*(dyh)**2 &
                     /(dt*dxv(jds+1)*dyv*0.5*fpole)/100000.00
        curv(jds+1)=tan(tpv)/a
        fad(jds+1)=-dt/(3.*dxv(jds+1)*dyv*2.*2.)
        fah(jds+1)=-dt/(3.*dxv(jds+1)*dyv*0.5*fpole)
        fcp(jds+1)= dt/(3.*dxv(jds+1)*dyv*0.5*cp*fpole)
        fdiv(jds+1)=2./(3.*dxv(jds+1)*dyv*0.5*fpole)
!
        hdacxj(jds+1)=0.
        hdacyj(jds+1)=acdt*dyh**2 &
                  /(4.*dxv(jds+1)*dyv*0.5*fpole)
        hdacvxj(jds+1)=acdt*dyv**2 &
                  /(4.*dxv(jds+1)*dyv)
        hdacvyj(jds+1)=acdt*dyv**2 &
                  /(4.*dxv(jds+1)*dyv)
!
!        ddmpu(jds+1)=cddamp*dxv(jds+1)/(2.*dxv(jds+1))
        ddmpu(jds+1)=cddamp*dyv/(2.*dxv(jds+1))
!-----------------------------------------------------------------------
!-------------between the poles----------------------------------------- 
!-----------------------------------------------------------------------
        do j=jds+2,jde-2
          tph=sb+(j-jds-1)*dph
          tpv=tph+dph*0.5
          dxh(j)=a*dlm*cos(tph)
          dxv(j)=a*dlm*cos(tpv)
          rdxh(j)=1./dxh(j)
          rdxv(j)=1./dxv(j)
          dare(j)=dxh(j)*dyh
          rare(j)=1./dare(j)
          wpdar(j)=-wcor*(dyh)**2 &
                  /(dt*dxh(j)*dyh)/100000.00
          curv(j)=tan(tpv)/a
          fad(j)=-dt/(3.*dxv(j)*dyv*2.*2.)
          fah(j)=-dt/(3.*dxh(j)*dyh)
          fcp(j)= dt/(3.*dxh(j)*dyh*cp)
          fdiv(j)=2./(3.*dxh(j)*dyh)
!
          hdacxj(j)= acdt*dyh*max(dxh(j),dyh)/(4.*dxh(j)*dyh)
          hdacyj(j)= acdt*dyh*max(dxh(j),dyh)/(4.*dxh(j)*dyh)
          hdacvxj(j)=acdt*dyv*max(dxv(j),dyv)/(4.*dxv(j)*dyv)
          hdacvyj(j)=acdt*dyv*max(dxv(j),dyv)/(4.*dxv(j)*dyv)
!
!          ddmpu(j)=cddamp*dxv(j)/(2.*dxv(j))
          ddmpu(j)=cddamp*dyv/(2.*dxv(j))
        enddo
!-----------------------------------------------------------------------
!-------------ghost line beyond the south pole---------------------------- 
!-----------------------------------------------------------------------
        dxh(jds)=dxh(jds+2)
        dxv(jds)=dxv(jds+1)
        rdxh(jds)=rdxh(jds+2)
        rdxv(jds)=rdxv(jds+1)
        dare(jds)=dare(jds+2)
        rare(jds)=rare(jds+2)
        wpdar(jds)=wpdar(jds+2)
        curv(jds)=curv(jds+1)
        fad(jds)=fad(jds+1)
        fah(jds)=fah(jds+2)
        fcp(jds)=fcp(jds+2)
        fdiv(jds)=fdiv(jds+2)
        hdacxj(jds)=hdacxj(jds+2)
        hdacyj(jds)=hdacyj(jds+2)
        hdacvxj(jds)=hdacvxj(jds+1)
        hdacvyj(jds)=hdacvyj(jds+1)
        ddmpu(jds)=ddmpu(jds+1)
!-----------------------------------------------------------------------
!-------------north pole------------------------------------------------
!-----------------------------------------------------------------------
        tph=tph+dph
        tpv=tph+dph*0.5
        dxh(jde-1)=0.
        rdxh(jde-1)=0.
        dare(jde-1)=dxv(jde-2)*dyv*0.5*fpole
        rare(jde-1)=1./dare(jde-1)
        wpdar(jde-1)=-wcor*(dyh)**2 &
                     /(dt*dxv(jde-2)*dyv*0.5*fpole)/100000.00
        fah(jde-1)=-dt/(3.*dxv(jde-2)*dyv*0.5*fpole)
        fcp(jde-1)= dt/(3.*dxv(jde-2)*dyv*0.5*fpole*cp)
        fdiv(jde-1)=2./(3.*dxv(jde-2)*dyv*0.5*fpole)
        hdacxj(jde-1)=0.
        hdacyj(jde-1)=acdt*dyh**2 &
                    /(4.*dxv(jde-2)*dyv*0.5*fpole)
!-----------------------------------------------------------------------
!-------------ghost line beyond north pole------------------------------
!-----------------------------------------------------------------------
        dxh(jde)=dxh(jde-2)
        dxv(jde-1)=dxv(jde-2)
        dxv(jde)=dxv(jde-2)
        rdxh(jde)=rdxh(jde-2)
        rdxv(jde-1)=rdxv(jde-2)
        rdxv(jde)=rdxv(jde-2)
        dare(jde)=dare(jde-2)
        rare(jde)=rare(jde-2)
        wpdar(jde)=wpdar(jde-2)
        curv(jde-1)=curv(jde-2)
        curv(jde)=curv(jde-2)
        fad(jde-1)=fad(jde-2)
        fad(jde)=fad(jde-2)
        fah(jde)=fah(jde-2)
        fcp(jde)=fcp(jde-2)
        fdiv(jde)=fdiv(jde-2)
        hdacxj(jde)=hdacxj(jde-2)
        hdacyj(jde)=hdacyj(jde-2)
        hdacvxj(jde-1)=hdacvxj(jde-2)
        hdacvyj(jde-1)=hdacvyj(jde-2)
        hdacvxj(jde)=hdacvxj(jde-2)
        hdacvyj(jde)=hdacvyj(jde-2)
        ddmpu(jde-1)=ddmpu(jde-2)
        ddmpu(jde)=ddmpu(jde-2)
!-----------------------------------------------------------------------
!-------------averaging over height latitudes for accuracy-------------- 
!-----------------------------------------------------------------------
        do j=jds,jde/2
          ave=(dxh(j)+dxh(jde+1-j))*0.5
          dxh(j)=ave
          dxh(jde+1-j)=ave
          ave=(rdxh(j)+rdxh(jde+1-j))*0.5
          rdxh(j)=ave
          rdxh(jde+1-j)=ave
          ave=(dare(j)+dare(jde+1-j))*0.5
          dare(j)=ave
          dare(jde+1-j)=ave
          ave=(rare(j)+rare(jde+1-j))*0.5
          rare(j)=ave
          rare(jde+1-j)=ave
          ave=(wpdar(j)+wpdar(jde+1-j))*0.5
          wpdar(j)=ave
          wpdar(jde+1-j)=ave
          ave=(fah(j)+fah(jde+1-j))*0.5
          fah(j)=ave
          fah(jde+1-j)=ave
          ave=(fcp(j)+fcp(jde+1-j))*0.5
          fcp(j)=ave
          fcp(jde+1-j)=ave
          ave=(fdiv(j)+fdiv(jde+1-j))*0.5
          fdiv(j)=ave
          fdiv(jde+1-j)=ave
          ave=(hdacxj(j)+hdacxj(jde+1-j))*0.5
          hdacxj(j)=ave
          hdacxj(jde+1-j)=ave
          ave=(hdacyj(j)+hdacyj(jde+1-j))*0.5
          hdacyj(j)=ave
          hdacyj(jde+1-j)=ave
        enddo
!-----------------------------------------------------------------------
!-------------averaging over wind latitudes for accuracy--------------- 
!-----------------------------------------------------------------------
        do j=jds,(jde-1)/2
          ave=(dxv(j)+dxv(jde-j))*0.5
          dxv(j)=ave
          dxv(jde-j)=ave
          ave=(rdxv(j)+rdxv(jde-j))*0.5
          rdxv(j)=ave
          rdxv(jde-j)=ave
          ave=(fad(j)+fad(jde-j))*0.5
          fad(j)=ave
          fad(jde-j)=ave
          ave=(hdacvxj(j)+hdacvxj(jde-j))*0.5
          hdacvxj(j)=ave
          hdacvxj(jde-j)=ave
          ave=(hdacvyj(j)+hdacvyj(jde-j))*0.5
          hdacvyj(j)=ave
          hdacvyj(jde-j)=ave
          ave=(ddmpu(j)+ddmpu(jde-j))*0.5
          ddmpu(j)=ave
          ddmpu(jde-j)=ave
        enddo
!-----------------------------------------------------------------------
!-------------diagonal distances at v points----------------------------
!-----------------------------------------------------------------------
        if(s_bdy)then
          ddv(jds+1)=dyv
          ddv(jds)=ddv(jds+1)
          rddv(jds+1)=1./ddv(jds+1)
          rddv(jds)=rddv(jds+1)
        endif
!
        do j=max(jms,jds+2),min(jme,jde-3)
          ddv(j)=sqrt(dxv(j)**2+dyv**2)
          rddv(j)=1./ddv(j)
        enddo
!
        if(n_bdy)then
          ddv(jde-2)=dyv
          ddv(jde-1)=ddv(jde-2)
          rddv(jde-2)=1./ddv(jde-2)
          rddv(jde-1)=rddv(jde-2)
        endif
!-----------------------------------------------------------------------
!-------------defining the diffusion coefficient inside the domain------
!-----------------------------------------------------------------------
!       do j=jms,jme
!         do i=ims,ime
        do j=j_lo,j_hi
          do i=i_lo,i_hi
            hdacx(i,j)=hdacxj(j)
            hdacy(i,j)=hdacyj(j)
            hdacvx(i,j)=hdacvxj(j)
            hdacvy(i,j)=hdacvyj(j)
          enddo
        enddo
!-----------------------------------------------------------------------
!-------------coriolis parameter in tll system--------------------------
!-----------------------------------------------------------------------
        tph_base=sb-dph-dph*.5
!       do j=jms,jme
        do j=j_lo,j_hi
          tph=tph_base+(j-jds+1)*dph
          stph=sin(tph)
          ctph=cos(tph)
!
          tlm_base=wb-dlm-dlm*0.5
!         do i=ims,ime
          do i=i_lo,i_hi
            tlm=tlm_base+(i-ids+1)*dlm
            f(i,j)=twom*(ctph0*stph+stph0*ctph*cos(tlm))
          enddo
        enddo
!-----------------------------------------------------------------------
!--------------latitudes and longitudes of h points in radians----------     
!-----------------------------------------------------------------------
!       tph_base=sb-dph-dph
!       do j=j_lo,j_hi
!         tlm_base=wb-dlm-dlm
!         tph=tph_base+(j-jds+1)*dph
!         stph=sin(tph)
!         ctph=cos(tph)                                                                       
!         do i=i_lo,i_hi
!           tlm=tlm_base+(i-ids+1)*dlm
!           stlm=sin(tlm)
!           ctlm=cos(tlm)
!           sph=ctph0*stph+stph0*ctph*ctlm
!           aph=asin(sph)
!           glat(i,j)=aph
!           anum=ctph*stlm
!           denom=(ctlm*ctph-stph0*sph)/ctph0
!           relm=atan2(anum,denom)
!           alm=relm+tlm0d*dtr
!           if(alm> pi) alm=alm-pi-pi
!           if(alm< -pi) alm=alm+pi+pi
!           glon(i,j)=alm
!         enddo
!       enddo
        tph_base=sb-dph-dph
        do j=j_lo,j_hi
          tlm_base=wb-dlm-dlm
          aph=tph_base+(j-jds+1)*dph
          do i=i_lo,i_hi
            alm=tlm_base+(i-ids+1)*dlm
            if(alm> pi)then
              if(i<=ite-2)then
                alm=alm-pi-pi
              else
                alm=pi+(i-ite+1)*dlm
              endif
            endif
            if(alm<-pi)then
              if(i>=its+2)then
                alm=alm+pi+pi
              else
                alm=-pi+(i-its-1)*dlm
              endif
            endif
            glat(i,j)=aph
            glon(i,j)=alm
          enddo
        enddo
!
!***  Repeat the preceding loop for the offset V points
!
        tph_base=sb-dph-0.5*dph
        do j=j_lo,j_hi
          tlm_base=wb-dlm-0.5*dlm
          aph=tph_base+(j-jds+1)*dph
          do i=i_lo,i_hi
            alm=tlm_base+(i-ids+1)*dlm
            if(alm> pi) alm=alm-pi-pi
            if(alm<-pi) alm=alm+pi+pi
            vlat(i,j)=aph
            vlon(i,j)=alm
          enddo
        enddo
!
!-----------------------------------------------------------------------
!***  Save the geographic lat/lon (radians) of this domain's SW corner.
!-----------------------------------------------------------------------
!
        arg=sin(sb)*ctph0+cos(sb)*stph0*cos(wb)
        arg=min(arg,1.)
        arg=max(arg,-1.)
        glat_sw=asin(arg)
!
        arg=(cos(sb)*cos(wb))/(cos(glat_sw)*ctph0)                      &
            -tan(glat_sw)*tan(tph0d*dtr)
        arg=min(arg,1.)
        arg=max(arg,-1.)
        glon_sw=tlm0d*dtr+sign(1.,wb)*acos(arg)
!
!-----------------------------------------------------------------------
      else !regional
!-----------------------------------------------------------------------
!
!---regional branch-----------------------------------------------------
!
!-----------------------------------------------------------------------
!-------------between the poles----------------------------------------- 
!-----------------------------------------------------------------------
!
        ddmpv=cddamp*dyv/(2.*dyv)
!
        do j=jds,jde
          tph=sb+(j-jds)*dph
          tpv=tph+dph*0.5
          dxh(j)=a*dlm*cos(tph)
          dxv(j)=a*dlm*cos(tpv)
          rdxh(j)=1./dxh(j)
          rdxv(j)=1./dxv(j)
          dare(j)=dxh(j)*dyh
          rare(j)=1./dare(j)
          wpdar(j)=-wcor*(dyh)**2 &
                  /(dt*dxh(j)*dyh)/100000.00
          curv(j)=tan(tpv)/a
          fad(j)=-dt/(3.*dxv(j)*dyv*2.*2.)
          fah(j)=-dt/(3.*dxh(j)*dyh)
          fcp(j)= dt/(3.*dxh(j)*dyh*cp)
          fdiv(j)=2./(3.*dxh(j)*dyh)
!
          hdacxj(j)= acdt*dyh*max(dxh(j),dyh)/(4.*dxh(j)*dyh)
          hdacyj(j)= acdt*dyh*max(dxh(j),dyh)/(4.*dxh(j)*dyh)
          hdacvxj(j)=acdt*dyv*max(dxv(j),dyv)/(4.*dxv(j)*dyv)
          hdacvyj(j)=acdt*dyv*max(dxv(j),dyv)/(4.*dxv(j)*dyv)
!
!          ddmpu(j)=cddamp*dxv(j)/(2.*dxv(j))
          ddmpu(j)=cddamp*dyv/(2.*dxv(j))
        enddo
!-----------------------------------------------------------------------
!-------------diagonal distances at v points----------------------------
!-----------------------------------------------------------------------
        do j=jds,jde
          ddv(j)=sqrt(dxv(j)**2+dyv**2)
          rddv(j)=1./ddv(j)
        enddo
!-----------------------------------------------------------------------
!---defining the diffusion coefficient inside the domain----------------
!-----------------------------------------------------------------------
!       do j=jms,jme
!         do i=ims,ime
        do j=j_lo,j_hi
          do i=i_lo,i_hi
            hdacx(i,j)=hdacxj(j)
            hdacy(i,j)=hdacyj(j)
            hdacvx(i,j)=hdacvxj(j)
            hdacvy(i,j)=hdacvyj(j)
          enddo
        enddo
!-----------------------------------------------------------------------
!---enhancing diffusion along boundaries--------------------------------
!-----------------------------------------------------------------------
        if(lnsbc>=2)then
!
          if(s_bdy)then
            do j=jts+1,jts-1+lnsbc
              do i=max(its,ids+1),min(ite,ide-1)
                hdacx(i,j)=hdacxj(j)*bofac
                hdacy(i,j)=hdacyj(j)*bofac
              enddo
            enddo
          endif
!
          if(n_bdy)then
            do j=jte+1-lnsbc,jte-1
              do i=max(its,ids+1),min(ite,ide-1)
                hdacx(i,j)=hdacxj(j)*bofac
                hdacy(i,j)=hdacyj(j)*bofac
              enddo
            enddo
          endif
!
        endif
!
        if(w_bdy)then
          do j=max(jts,jds+lnsbc),min(jte,jde-lnsbc)
            do i=its+1,its-1+lnsbc
              hdacx(i,j)=hdacxj(j)*bofac
              hdacy(i,j)=hdacyj(j)*bofac
            enddo
          enddo
        endif
!
        if(e_bdy)then
          do j=max(jts,jds+lnsbc),min(jte,jde-lnsbc)
            do i=ite+1-lnsbc,ite-1
              hdacx(i,j)=hdacxj(j)*bofac
              hdacy(i,j)=hdacyj(j)*bofac
            enddo
          enddo
        endif
!
!-----------------------------------------------------------------------
        if(s_bdy)then
          do j=jts+1,jts-1+lnsbc
            do i=max(its,ids+1),min(ite,ide-2)
              hdacvx(i,j)=hdacvxj(j)*bofac
              hdacvy(i,j)=hdacvyj(j)*bofac
            enddo
          enddo
        endif
!
        if(n_bdy)then
          do j=jte-lnsbc,jte-2
            do i=max(its,ids+1),min(ite,ide-2)
              hdacvx(i,j)=hdacvxj(j)*bofac
              hdacvy(i,j)=hdacvyj(j)*bofac
            enddo
          enddo
        endif
!
        if(w_bdy)then
          do j=max(jts,jds+lnsbc),min(jte,jde-1-lnsbc)
            do i=its+1,its-1+lnsbc
              hdacvx(i,j)=hdacvxj(j)*bofac
              hdacvy(i,j)=hdacvyj(j)*bofac
            enddo
          enddo
        endif
!
        if(e_bdy)then
          do j=max(jts,jds+lnsbc),min(jte,jde-1-lnsbc)
            do i=ite-lnsbc,ite-2
              hdacvx(i,j)=hdacvxj(j)*bofac
              hdacvy(i,j)=hdacvyj(j)*bofac
            enddo
          enddo
        endif
!
!-----------------------------------------------------------------------
!-------------coriolis parameter in tll system--------------------------
!-----------------------------------------------------------------------
      
        tph_base=sb-dph*.5
!
        do j=jts,jte
          tph=tph_base+(j-jds+1)*dph
          stph=sin(tph)
          ctph=cos(tph)
!
          tlm_base=wb-dlm*0.5
          do i=its,ite
            tlm=tlm_base+(i-ids+1)*dlm
            f(i,j)=twom*(ctph0*stph+stph0*ctph*cos(tlm))
          enddo
        enddo
!-----------------------------------------------------------------------
!--------------latitudes and longitudes of h points in radians----------     
!-----------------------------------------------------------------------
        tph_base=sb-dph
!	write(0,*) 'sb,dph: ', sb, dph
!	write(0,*) 'wb,dlm: ', wb, dlm
!
        do j=j_lo,j_hi
          tph=tph_base+(j-jds+1)*dph
          stph=sin(tph)
          ctph=cos(tph)                                                                       
!
          tlm_base=wb-dlm
          do i=i_lo,i_hi
            tlm=tlm_base+(i-ids+1)*dlm
            stlm=sin(tlm)
            ctlm=cos(tlm)
            sph=ctph0*stph+stph0*ctph*ctlm
            aph=asin(sph)
            glat(i,j)=aph
            anum=ctph*stlm
            denom=(ctlm*ctph-stph0*sph)/ctph0
            relm=atan2(anum,denom)
            alm=relm+tlm0d*dtr
            if(alm>  pi) alm=alm-pi-pi
            if(alm< -pi) alm=alm+pi+pi
            glon(i,j)=alm
          enddo
        enddo
!
!***  Repeat preceding loop for V lat/lon offset
!
        tph_base=sb-0.5*dph
!
        do j=jts,jte
          tph=tph_base+(j-jds+1)*dph
          stph=sin(tph)
          ctph=cos(tph)
!
          tlm_base=wb-0.5*dlm
          do i=its,ite
            tlm=tlm_base+(i-ids+1)*dlm
            stlm=sin(tlm)
            ctlm=cos(tlm)
            sph=ctph0*stph+stph0*ctph*ctlm
            aph=asin(sph)
            vlat(i,j)=aph
            anum=ctph*stlm
            denom=(ctlm*ctph-stph0*sph)/ctph0
            relm=atan2(anum,denom)
            alm=relm+tlm0d*dtr
            if(alm>  pi) alm=alm-pi-pi
            if(alm< -pi) alm=alm+pi+pi
            vlon(i,j)=alm
          enddo
        enddo
!
!-----------------------------------------------------------------------
!***  Save the geographic lat/lon (radians) of this domain's SW corner.
!-----------------------------------------------------------------------
!
        arg=sin(sb)*ctph0+cos(sb)*stph0*cos(wb)
        arg=min(arg,1.)
        arg=max(arg,-1.)
        glat_sw=asin(arg)
!
        arg=(cos(sb)*cos(wb))/(cos(glat_sw)*ctph0)                      &
            -tan(glat_sw)*tan(tph0d*dtr)
        arg=min(arg,1.)
        arg=max(arg,-1.)
        glon_sw=tlm0d*dtr+sign(1.,wb)*acos(arg)
!
!-----------------------------------------------------------------------
!
      endif global_regional_setup

                        end subroutine grid_consts

!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
                        subroutine boundary_init &
      (its,ite,jts,jte,lm &
      ,ims,ime,jms,jme &
      ,ids,ide,jds,jde &
      ,lnsh,lnsv &
      ,nvars_2d_h,nvars_3d_h,nvars_4d_h &
      ,nvars_2d_v,nvars_3d_v &
      ,bnd_vars_h &
      ,bnd_vars_v &
       )
!
!-----------------------------------------------------------------------
!***  Initialize boundary variable arrays for nested domains.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
 
!------------------------
!***  Argument variables
!------------------------
 
integer(kind=kint),intent(in) :: &
 ids &
,ide &
,ims &
,ime &
,its &
,ite &
,jds &
,jde &
,jms &
,jme &
,jts &
,jte &
,lm &
,lnsh &
,lnsv

integer(kind=kint),intent(in) :: &
 nvars_2d_h &
,nvars_3d_h &
,nvars_4d_h &
,nvars_2d_v &
,nvars_3d_v

type(bc_h_all),intent(inout) :: bnd_vars_h

type(bc_v_all),intent(inout) :: bnd_vars_v
 
!---------------------
!***  Local variables
!---------------------
 
integer(kind=kint) :: &
 i &
,j &
,k &
,l &
,l1 &
,l2 &
,n &
,nv
 
real(kind=kfpt),dimension(:,:),pointer :: var2d
real(kind=kfpt),dimension(:,:,:),pointer :: bnd2d,var3d
real(kind=kfpt),dimension(:,:,:,:),pointer :: bnd3d,var4d
real(kind=kfpt),dimension(:,:,:,:,:),pointer :: bnd4d
 
logical(kind=klog) :: &
 e_bdy &
,n_bdy &
,s_bdy &
,w_bdy
 
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      w_bdy=(its==ids)  ! This task is on the western boundary
      e_bdy=(ite==ide)  ! This task is on the eastern boundary
      s_bdy=(jts==jds)  ! This task is on the southern boundary
      n_bdy=(jte==jde)  ! This task is on the northern boundary
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  South
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!
      if(s_bdy)then
!
!-----------------------------------------------------------------------
!
!------------------------
!***  2-d h variables
!------------------------
!
        if(nvars_2d_h>0)then 
!
          do nv=1,nvars_2d_h         
!
            bnd2d=>bnd_vars_h%var_2d(nv)%south
            var2d=>bnd_vars_h%var_2d(nv)%full_var
!
            n=0
            do j=1,lnsh
              n=n+1
              do i=ims,ime
                bnd2d(i,j,1)=var2d(i,jds-1+n)
              enddo
            enddo
!
          enddo
!
        endif
!
!------------------------
!***  3-d h variables
!------------------------
!
        if(nvars_3d_h>0)then 
!
          do nv=1,nvars_3d_h         
!
            bnd3d=>bnd_vars_h%var_3d(nv)%south
            var3d=>bnd_vars_h%var_3d(nv)%full_var
!
            do k=1,lm
              n=0
              do j=1,lnsh
                n=n+1
                do i=ims,ime
                  bnd3d(i,j,k,1)=var3d(i,jds-1+n,k)
                enddo
              enddo
            enddo
!
          enddo
!
        endif
!
!------------------------
!***  4-d h variables
!------------------------
!
        if(nvars_4d_h>0)then
!
          do nv=1,nvars_4d_h
!
            bnd4d=>bnd_vars_h%var_4d(nv)%south
            var4d=>bnd_vars_h%var_4d(nv)%full_var
!
            l1=lbound(bnd_vars_h%var_4d(nv)%full_var,4)
            l2=ubound(bnd_vars_h%var_4d(nv)%full_var,4)
            do l=l1,l2
              do k=1,lm
                n=0
                do j=1,lnsh
                  n=n+1
                  do i=ims,ime
                    bnd4d(i,j,k,l,1)=var4d(i,jds-1+n,k,l)
                  enddo
                enddo
              enddo
            enddo
!
          enddo
!
        endif
!
!------------------------
!***  2-d v variables
!------------------------
!
        if(nvars_2d_v>0)then 
!
          do nv=1,nvars_2d_v         
!
            bnd2d=>bnd_vars_v%var_2d(nv)%south
            var2d=>bnd_vars_v%var_2d(nv)%full_var
!
            n=0
            do j=1,lnsv
              n=n+1
              do i=ims,ime
                bnd2d(i,j,1)=var2d(i,jds-1+n)
              enddo
            enddo
!
          enddo
!
        endif
!
!------------------------
!***  3-d v variables
!------------------------
!
        if(nvars_3d_v>0)then 
!
          do nv=1,nvars_3d_v         
!
            bnd3d=>bnd_vars_v%var_3d(nv)%south
            var3d=>bnd_vars_v%var_3d(nv)%full_var
!
            do k=1,lm
              n=0
              do j=1,lnsv
                n=n+1
                do i=ims,ime
                  bnd3d(i,j,k,1)=var3d(i,jds-1+n,k)
                enddo
              enddo
            enddo
!
          enddo
!
        endif
!
      endif
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  North
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      if(n_bdy)then
!
!------------------------
!***  2-d h variables
!------------------------
!
        if(nvars_2d_h>0)then 
!
          do nv=1,nvars_2d_h         
!
            bnd2d=>bnd_vars_h%var_2d(nv)%north
            var2d=>bnd_vars_h%var_2d(nv)%full_var
!
            n=0
            do j=1,lnsh
              n=n+1
              do i=ims,ime
                bnd2d(i,j,1)=var2d(i,jde-lnsh+n)
              enddo
            enddo
!
          enddo
!
        endif
!
!------------------------
!***  3-d h variables
!------------------------
!
        if(nvars_3d_h>0)then 
!
          do nv=1,nvars_3d_h         
!
            bnd3d=>bnd_vars_h%var_3d(nv)%north
            var3d=>bnd_vars_h%var_3d(nv)%full_var
!
            do k=1,lm
              n=0
              do j=1,lnsh
                n=n+1
                do i=ims,ime
                  bnd3d(i,j,k,1)=var3d(i,jde-lnsh+n,k)
                enddo
              enddo
            enddo
!
          enddo
!
        endif
!
!------------------------
!***  4-d h variables
!------------------------
!
        if(nvars_4d_h>0)then
!
          do nv=1,nvars_4d_h
!
            bnd4d=>bnd_vars_h%var_4d(nv)%north
            var4d=>bnd_vars_h%var_4d(nv)%full_var
!
            l1=lbound(bnd_vars_h%var_4d(nv)%full_var,4)
            l2=ubound(bnd_vars_h%var_4d(nv)%full_var,4)
            do l=l1,l2
              do k=1,lm
                n=0
                do j=1,lnsh
                  n=n+1
                  do i=ims,ime
                    bnd4d(i,j,k,l,1)=var4d(i,jde-lnsh+n,k,l)
                  enddo
                enddo
              enddo
            enddo
!
          enddo
!
        endif
!
!------------------------
!***  2-d v variables
!------------------------
!
        if(nvars_2d_v>0)then 
!
          do nv=1,nvars_2d_v         
!
            bnd2d=>bnd_vars_v%var_2d(nv)%north
            var2d=>bnd_vars_v%var_2d(nv)%full_var
!
            n=0
            do j=1,lnsv
              n=n+1
              do i=ims,ime
                bnd2d(i,j,1)=var2d(i,jde-1-lnsv+n)
              enddo
            enddo
!
          enddo
!
        endif
!
!------------------------
!***  3-d v variables
!------------------------
!
        if(nvars_3d_v>0)then 
!
          do nv=1,nvars_3d_v         
!
            bnd3d=>bnd_vars_v%var_3d(nv)%north
            var3d=>bnd_vars_v%var_3d(nv)%full_var
!
            do k=1,lm
              n=0
              do j=1,lnsv
                n=n+1
                do i=ims,ime
                  bnd3d(i,j,k,1)=var3d(i,jde-1-lnsv+n,k)
                enddo
              enddo
            enddo
!
          enddo
!
        endif
!
      endif
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  West
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      if(w_bdy)then
!
!------------------------
!***  2-d h variables
!------------------------
!
        if(nvars_2d_h>0)then 
!
          do nv=1,nvars_2d_h         
!
            bnd2d=>bnd_vars_h%var_2d(nv)%west
            var2d=>bnd_vars_h%var_2d(nv)%full_var
!
            do j=jms,jme
              n=0
              do i=1,lnsh
                n=n+1
                bnd2d(i,j,1)=var2d(ids-1+n,j)
              enddo
            enddo
!
          enddo
!
        endif
!
!------------------------
!***  3-d h variables
!------------------------
!
        if(nvars_3d_h>0)then 
!
          do nv=1,nvars_3d_h         
!
            bnd3d=>bnd_vars_h%var_3d(nv)%west
            var3d=>bnd_vars_h%var_3d(nv)%full_var
!
            do k=1,lm
              do j=jms,jme
                n=0
                do i=1,lnsh
                  n=n+1
                  bnd3d(i,j,k,1)=var3d(ids-1+n,j,k)
                enddo
              enddo
            enddo
!
          enddo
!
        endif
!
!------------------------
!***  4-d h variables
!------------------------
!
        if(nvars_4d_h>0)then
!
          do nv=1,nvars_4d_h
!
            bnd4d=>bnd_vars_h%var_4d(nv)%west
            var4d=>bnd_vars_h%var_4d(nv)%full_var
!
            l1=lbound(bnd_vars_h%var_4d(nv)%full_var,4)
            l2=ubound(bnd_vars_h%var_4d(nv)%full_var,4)
            do l=l1,l2
              do k=1,lm
                do j=jms,jme
                  n=0
                  do i=1,lnsh
                    n=n+1
                    bnd4d(i,j,k,l,1)=var4d(ids-1+n,j,k,l)
                  enddo
                enddo
              enddo
            enddo
!
          enddo
!
        endif
!
!------------------------
!***  2-d v variables
!------------------------
!
        if(nvars_2d_v>0)then 
!
          do nv=1,nvars_2d_v         
!
            bnd2d=>bnd_vars_v%var_2d(nv)%west
            var2d=>bnd_vars_v%var_2d(nv)%full_var
!
            do j=jms,jme
              n=0
              do i=1,lnsv
                n=n+1
                bnd2d(i,j,1)=var2d(ids-1+n,j)
              enddo
            enddo
!
          enddo
!
        endif
!
!------------------------
!***  3-d v variables
!------------------------
!
        if(nvars_3d_v>0)then 
!
          do nv=1,nvars_3d_v         
!
            bnd3d=>bnd_vars_v%var_3d(nv)%west
            var3d=>bnd_vars_v%var_3d(nv)%full_var
!
            do k=1,lm
              do j=jms,jme
                n=0
                do i=1,lnsv
                  n=n+1
                  bnd3d(i,j,k,1)=var3d(ids-1+n,j,k)
                enddo
              enddo
            enddo
!
          enddo
!
        endif
!
      endif
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  East
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      if(e_bdy)then
!
!------------------------
!***  2-d h variables
!------------------------
!
        if(nvars_2d_h>0)then 
!
          do nv=1,nvars_2d_h         
!
            bnd2d=>bnd_vars_h%var_2d(nv)%east
            var2d=>bnd_vars_h%var_2d(nv)%full_var
!
            do j=jms,jme
              n=0
              do i=1,lnsh
                n=n+1
                bnd2d(i,j,1)=var2d(ide-lnsh+n,j)
              enddo
            enddo
!
          enddo
!
        endif
!
!------------------------
!***  3-d h variables
!------------------------
!
        if(nvars_3d_h>0)then 
!
          do nv=1,nvars_3d_h         
!
            bnd3d=>bnd_vars_h%var_3d(nv)%east
            var3d=>bnd_vars_h%var_3d(nv)%full_var
!
            do k=1,lm
              do j=jms,jme
                n=0
                do i=1,lnsh
                  n=n+1
                  bnd3d(i,j,k,1)=var3d(ide-lnsh+n,j,k)
                enddo
              enddo
            enddo
!
          enddo
!
        endif
!
!------------------------
!***  4-d h variables
!------------------------
!
        if(nvars_4d_h>0)then
!
          do nv=1,nvars_4d_h
!
            bnd4d=>bnd_vars_h%var_4d(nv)%east
            var4d=>bnd_vars_h%var_4d(nv)%full_var
!
            l1=lbound(bnd_vars_h%var_4d(nv)%full_var,4)
            l2=ubound(bnd_vars_h%var_4d(nv)%full_var,4)
            do l=l1,l2
              do k=1,lm
                do j=jms,jme
                  n=0
                  do i=1,lnsh
                    n=n+1
                    bnd4d(i,j,k,l,1)=var4d(ids-1+n,j,k,l)
                  enddo
                enddo
              enddo
            enddo
!
          enddo
!
        endif
!
!------------------------
!***  2-d v variables
!------------------------
!
        if(nvars_2d_v>0)then 
!
          do nv=1,nvars_2d_v         
!
            bnd2d=>bnd_vars_v%var_2d(nv)%east
            var2d=>bnd_vars_v%var_2d(nv)%full_var
!
            do j=jms,jme
              n=0
              do i=1,lnsv
                n=n+1
                bnd2d(i,j,1)=var2d(ide-1-lnsv+n,j)
              enddo
            enddo
!
          enddo
!
        endif
!
!------------------------
!***  3-d v variables
!------------------------
!
        if(nvars_3d_v>0)then 
!
          do nv=1,nvars_3d_v         
!
            bnd3d=>bnd_vars_v%var_3d(nv)%east
            var3d=>bnd_vars_v%var_3d(nv)%full_var
!
            do k=1,lm
              do j=jms,jme
                n=0
                do i=1,lnsv
                  n=n+1
                  bnd3d(i,j,k,1)=var3d(ide-1-lnsv+n,j,k)
                enddo
              enddo
            enddo
!
          enddo
!
        endif
!
      endif
!
!-----------------------------------------------------------------------
!
      end subroutine boundary_init
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
                        subroutine exptbl
!     ******************************************************************
!     *                                                                *
!     *               exponential function table                       *
!     *               responsible person: z.janjic                     *
!     *                                                                *
!     ******************************************************************
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!--local variables------------------------------------------------------
!-----------------------------------------------------------------------
integer(kind=kint):: &
 k                           ! index

real(kind=kfpt):: &
 x &                         ! argument
,xrng                        ! argument range
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
      xmax= 30.
      xmin=-30.
!
      kexm2=kexm-2
      xrng=xmax-xmin
!
      dex=xrng/(kexm-1)
      rdex=1./dex
!--------------function definition loop---------------------------------
      x=xmin-dex
      do k=1,kexm
        x=x+dex
        expf(k)=exp(x)
      enddo
!-----------------------------------------------------------------------
!
                        end subroutine exptbl
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
                        function zjexp(x)
!     ******************************************************************
!     *                                                                *
!     *               exponential function table                       *
!     *               responsible person: z.janjic                     *
!     *                                                                *
!     ******************************************************************
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
real(kind=kfpt):: &
 zjexp

!-----------------------------------------------------------------------
!--local variables------------------------------------------------------
!-----------------------------------------------------------------------
integer(kind=kint):: &
 k                           ! index

real(kind=kfpt):: &
 ak &                        ! position in table
,x                           ! argument
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
      ak=(x-xmin)*rdex
      k=int(ak)
      k=max(k,0)
      k=min(k,kexm2)
!
      zjexp=(expf(k+2)-expf(k+1))*(ak-real(k))+expf(k+1)
!-----------------------------------------------------------------------
!
                        end function zjexp
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
                        subroutine tablep
!     ******************************************************************
!     *                                                                *
!     *    generates the table for finding pressure from               *
!     *    saturation specific humidity and potential temperature      *
!     *                                                                *
!     ******************************************************************
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
real(kind=kfpt),parameter:: &
 eps=1.e-10
!-----------------------------------------------------------------------
!--local variables------------------------------------------------------
!-----------------------------------------------------------------------
integer(kind=kint):: &
 kth &                       ! index
,kp                          ! index

real(kind=kfpt):: &
 ape &                       ! exner function
,dth &                       ! potential temperature step
,dp &                        ! pressure step
,dqs &                       ! saturation specific humidity step
,p &                         ! pressure
,qs0k &                      ! base value for saturation humidity
,sqsk &                      ! saturation spec. humidity range
,th                          ! potential temperature

real(kind=kfpt),dimension(1:itb):: &
 app &                       ! temporary
,aqp &                       ! temporary
,pnew &                      ! new pressures
,pold &                      ! old pressure
,qsnew &                     ! new saturation spec. humidity
,qsold &                     ! old saturation spec. humidity
,y2p                         ! temporary
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
      dth=(thh-thl)/real(jtb-1)
      dp=(ph-pl)/real(itb-1)
      rdth=1./dth
!-----------------------------------------------------------------------
      th=thl-dth
      do kth=1,jtb
        th=th+dth
        p=pl-dp
        do kp=1,itb
          p=p+dp
          ape=(100000./p)**cappa
          qsold(kp)=pq0/p*exp(a2*(th-a3*ape)/(th-a4*ape))
          pold(kp)=p
        enddo
!
        qs0k=qsold(1)
        sqsk=qsold(itb)-qsold(1)
        qsold(1)=0.
        qsold(itb)=1.
!
        do kp=2,itb-1
          qsold(kp)=(qsold(kp)-qs0k)/sqsk
!wwwwwwwwwwwwww fix due to 32 bit precision limitation wwwwwwwwwwwwwwwww
          if((qsold(kp)-qsold(kp-1)).lt.eps) then
            qsold(kp)=qsold(kp-1)+eps
          endif
!mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
        enddo
!
        qs0(kth)=qs0k
        sqs(kth)=sqsk
!
        qsnew(1)=0.
        qsnew(itb)=1.
        dqs=1./real(itb-1)
        rdq=1./dqs
!
        do kp=2,itb-1
          qsnew(kp)=qsnew(kp-1)+dqs
        enddo
!
        y2p(1)=0.
        y2p(itb)=0.
!
        call spline(jtb,itb,qsold,pold,y2p,itb,qsnew,pnew,app,aqp)
!
        do kp=1,itb
          ptbl(kp,kth)=pnew(kp)
        enddo
!-----------------------------------------------------------------------
      enddo
!-----------------------------------------------------------------------
!
                        end subroutine tablep
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
                        subroutine tablet
!     ******************************************************************
!     *                                                                *
!     *    generates the table for finding temperature from            *
!     *    pressure and equivalent potential temperature               *
!     *                                                                *
!     ******************************************************************
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
real(kind=kfpt),parameter:: &
 eps=1.e-10
!-----------------------------------------------------------------------
!--local variables------------------------------------------------------
!-----------------------------------------------------------------------
integer(kind=kint):: &
 kth &                       ! index
,kp                          ! index

real(kind=kfpt):: &
 ape &                       ! exner function
,dth &                       ! potential temperature step
,dp &                        ! pressure step
,dthe &                      ! equivalent pot. temperature step
,p &                         ! pressure
,qs &                        ! saturation specific humidity
,the0k &                     ! base value for equivalent pot. temperature
,sthek &                     ! equivalent pot. temperature range
,th                          ! potential temperature

real(kind=kfpt),dimension(1:jtb):: &
 apt &                       ! temporary
,aqt &                       ! temporary
,tnew &                      ! new temperature
,told &                      ! old temperature
,thenew &                    ! new equivalent potential temperature
,theold &                    ! old equivalent potential temperature
,y2t                         ! temporary
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
      dth=(thh-thl)/real(jtb-1)
      dp=(ph-pl)/real(itb-1)
      rdp=1./dp
!-----------------------------------------------------------------------
      p=pl-dp
      do kp=1,itb
        p=p+dp
        th=thl-dth
        do kth=1,jtb
          th=th+dth
          ape=(100000./p)**cappa
          qs=pq0/p*exp(a2*(th-a3*ape)/(th-a4*ape))
          told(kth)=th/ape
          theold(kth)=th*exp(eliwv*qs/(cp*told(kth)))
        enddo
!
        the0k=theold(1)
        sthek=theold(jtb)-theold(1)
        theold(1)=0.
        theold(jtb)=1.
!
        do kth=2,jtb-1
          theold(kth)=(theold(kth)-the0k)/sthek
!wwwwwwwwwwwwww fix due to 32 bit precision limitation wwwwwwwwwwwwwwwww
          if((theold(kth)-theold(kth-1)).lt.eps) then
            theold(kth)=theold(kth-1)+eps
          endif
!mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
        enddo
!
        the0(kp)=the0k
        sthe(kp)=sthek
!
        thenew(1)=0.
        thenew(jtb)=1.
        dthe=1./real(jtb-1)
        rdthe=1./dthe
!
        do kth=2,jtb-1
          thenew(kth)=thenew(kth-1)+dthe
        enddo
!
        y2t(1)=0.
        y2t(jtb)=0.
!
        call spline(jtb,jtb,theold,told,y2t,jtb,thenew,tnew,apt,aqt)
!
        do kth=1,jtb
          ttbl(kth,kp)=tnew(kth)
        enddo
!-----------------------------------------------------------------------
      enddo
!-----------------------------------------------------------------------
!
                        end subroutine tablet
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
      subroutine spline(jtb,nold,xold,yold,y2,nnew,xnew,ynew,p,q)           
!     ******************************************************************
!     *                                                                *
!     *  this is a one-dimensional cubic spline fitting routine        *
!     *  programed for a small scalar machine.                         *
!     *                                                                *
!     *  programer: z. janjic, yugoslav fed. hydromet. inst., beograd  *
!     *                                                                *
!     *                                                                *
!     *                                                                *
!     *  nold - number of given values of the function.  must be ge 3. *
!     *  xold - locations of the points at which the values of the     *
!     *         function are given.  must be in ascending order.       *
!     *  yold - the given values of the function at the points xold.   *
!     *  y2   - the second derivatives at the points xold.  if natural *
!     *         spline is fitted y2(1)=0. and y2(nold)=0. must be      *
!     *         specified.                                             *
!     *  nnew - number of values of the function to be calculated.     *
!     *  xnew - locations of the points at which the values of the     *
!     *         function are calculated.  xnew(k) must be ge xold(1)   *
!     *         and le xold(nold).                                     *
!     *  ynew - the values of the function to be calculated.           *
!     *  p, q - auxiliary vectors of the length nold-2.                *
!     *                                                                *
!     ******************************************************************
!-----------------------------------------------------------------------
!
      integer,intent(in) :: jtb,nnew,nold
!
      real,dimension(jtb),intent(in) :: xnew,xold,yold
      real,dimension(jtb),intent(inout) :: p,q,y2
      real,dimension(jtb),intent(out) :: ynew
!
      integer :: k,k1,k2,kold,noldm1
!
      real :: ak,bk,ck,den,dx,dxc,dxl,dxr,dydxl,dydxr                   &
             ,rdx,rtdxc,x,xk,xsq,y2k,y2kp1
!
!-----------------------------------------------------------------------
      noldm1=nold-1                                                     
!                                                                       
      dxl=xold(2)-xold(1)                                               
      dxr=xold(3)-xold(2)                                               
      dydxl=(yold(2)-yold(1))/dxl                                       
      dydxr=(yold(3)-yold(2))/dxr                                       
      rtdxc=.5/(dxl+dxr)                                                
!                                                                       
      p(1)= rtdxc*(6.*(dydxr-dydxl)-dxl*y2(1))                          
      q(1)=-rtdxc*dxr                                                   
!                                                                       
      if(nold.eq.3) go to 700                                           
!-----------------------------------------------------------------------
      k=3                                                               
!                                                                       
 100  dxl=dxr                                                           
      dydxl=dydxr                                                       
      dxr=xold(k+1)-xold(k)                                             
      dydxr=(yold(k+1)-yold(k))/dxr                                     
      dxc=dxl+dxr                                                       
      den=1./(dxl*q(k-2)+dxc+dxc)                                       
!                                                                       
      p(k-1)= den*(6.*(dydxr-dydxl)-dxl*p(k-2))                         
      q(k-1)=-den*dxr                                                   
!                                                                       
      k=k+1                                                             
      if(k.lt.nold) go to 100                                           
!-----------------------------------------------------------------------
 700  k=noldm1                                                          
!                                                                       
 200  y2(k)=p(k-1)+q(k-1)*y2(k+1)                                       
!                                                                       
      k=k-1                                                             
      if(k.gt.1) go to 200                                              
!-----------------------------------------------------------------------
      k1=1                                                              
!                                                                       
 300  xk=xnew(k1)                                                       
!                                                                       
      do 400 k2=2,nold                                                  
      if(xold(k2).le.xk) go to 400                                      
      kold=k2-1                                                         
      go to 450                                                         
 400  continue                                                          
      ynew(k1)=yold(nold)                                               
      go to 600                                                         
!                                                                       
 450  if(k1.eq.1)   go to 500                                           
      if(k.eq.kold) go to 550                                           
!                                                                       
 500  k=kold                                                            
!                                                                       
      y2k=y2(k)                                                         
      y2kp1=y2(k+1)                                                     
      dx=xold(k+1)-xold(k)                                              
      rdx=1./dx                                                         
!                                                                       
      ak=.1666667*rdx*(y2kp1-y2k)                                       
      bk=.5*y2k                                                         
      ck=rdx*(yold(k+1)-yold(k))-.1666667*dx*(y2kp1+y2k+y2k)            
!                                                                       
 550  x=xk-xold(k)                                                      
      xsq=x*x                                                           
!                                                                       
      ynew(k1)=ak*xsq*x+bk*xsq+ck*x+yold(k)                             
!                                                                       
 600  k1=k1+1                                                           
      if(k1.le.nnew) go to 300                                          
!-----------------------------------------------------------------------
!
                        endsubroutine spline    
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
      subroutine psitbl
!     ******************************************************************
!     *                                                                *
!     *               surface layer integral functions                 *
!     *               responsible person: z.janjic                     *
!     *                                                                *
!     ******************************************************************
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
real(kind=kfpt),parameter:: &
 eps=0.000001
!--local variables------------------------------------------------------
integer(kind=kint):: &
 k                           ! index
 
real(kind=kfpt):: &
 x &                         ! temporary
,zeta1 &                     ! z/L, sea
,zeta2 &                     ! z/L, land
,zrng1 &                     ! z/L range, sea
,zrng2                       ! z/L range, land
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
      kztm2=kztm-2
!
      fh01=1.
      fh02=1.
!
      ztmin1=-5.0
      ztmax1= 1.0
!
      ztmin2=-5.0
      ztmax2= 1.0
!
      zrng1=ztmax1-ztmin1
      zrng2=ztmax2-ztmin2
!
      dzeta1=zrng1/(kztm-1)
      dzeta2=zrng2/(kztm-1)
!--------------function definition loop---------------------------------
      zeta1=ztmin1
      zeta2=ztmin2
      do k=1,kztm
!--------------unstable range-------------------------------------------
        if(zeta1.lt.0.)then
!--------------paulson 1970 functions-----------------------------------
          x=sqrt(sqrt(1.-16.*zeta1))
          psim1(k)=-2.*log((x+1.)/2.)-log((x*x+1.)/2.)+2.*atan(x)-pihf
          psih1(k)=-2.*log((x*x+1.)/2.)
!--------------stable range---------------------------------------------
        else
!--------------holtslag and de bruin 1988-------------------------------
          psim1(k)=0.7*zeta1+0.75*zeta1*(6.-0.35*zeta1)*exp(-0.35*zeta1)
          psih1(k)=0.7*zeta1+0.75*zeta1*(6.-0.35*zeta1)*exp(-0.35*zeta1)
!-----------------------------------------------------------------------
        endif
!--------------unstable range-------------------------------------------
        if(zeta2.lt.0.)then
!--------------paulson 1970 functions-----------------------------------
          x=sqrt(sqrt(1.-16.*zeta2))
          psim2(k)=-2.*log((x+1.)/2.)-log((x*x+1.)/2.)+2.*atan(x)-pihf
          psih2(k)=-2.*log((x*x+1.)/2.)
!--------------stable range---------------------------------------------
        else
!--------------holtslag and de bruin 1988-------------------------------
          psim2(k)=0.7*zeta2+0.75*zeta2*(6.-0.35*zeta2)*exp(-0.35*zeta2)
          psih2(k)=0.7*zeta2+0.75*zeta2*(6.-0.35*zeta2)*exp(-0.35*zeta2)
!-----------------------------------------------------------------------
        endif
!-----------------------------------------------------------------------
        if(k.eq.kztm)then
          ztmax1=zeta1
          ztmax2=zeta2
        endif
!
        zeta1=zeta1+dzeta1
        zeta2=zeta2+dzeta2
!-----------------------------------------------------------------------
      enddo
!-----------------------------------------------------------------------
      ztmax1=ztmax1-eps
      ztmax2=ztmax2-eps
!-----------------------------------------------------------------------
!
                        endsubroutine psitbl
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
      FUNCTION TIMEF()
!
!-----------------------------------------------------------------------
!
      REAL*8 TIMEF
      INTEGER(kind=KINT) :: IC,IR

      TIMEF=MPI_Wtime()
!
!-----------------------------------------------------------------------
!
      END FUNCTION TIMEF
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
      subroutine NMMB_Finalize
        integer irc
        CALL MPI_Finalize(irc)
        stop 911
!
      end subroutine NMMB_Finalize
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
                       end module module_control
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
