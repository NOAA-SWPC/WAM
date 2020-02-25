
!
      subroutine idea_visc_conduct(im,ix,levs,grav,prsi,prsl, rdelp,
     & up,dudt, dtp, vurho, murho, ahs_i, cp)
!-----------------------------------------------------------------------
!
! solve for  "temp, UV-wind" tendencies caused by  
!       molecular + eddy viscosity and thermal conductivity
!       all O-O2-N2 in "cm3"
!       to make this transform in GFS you need to "mmr" => 1/cm3
!
!vay-2015 take-out (levs+1)-needless extension for local arrays
!-----------------------------------------------------------------------

      use machine, only : kind_phys
      use idea_composition
      implicit none
     
!
! Argument
      integer, intent(in) :: im    ! number of data points in up,dudt(first dim)
      integer, intent(in) :: ix    ! max data points in fields
      integer, intent(in) :: levs  ! number of pressure levels
      real,    intent(in) :: dtp   ! time step in second

      real, intent(in)    :: prsi(ix,levs+1) ! interface pressure in KPa
      real, intent(in)    :: prsl(ix,levs)   ! layer pressure in KPa
 
      real, intent(in)     :: grav(ix,levs)   ! (m/s2) => to KpA  ....9.81e-3
      real, intent(in)     :: cp(ix, levs)
      real, intent(in)     :: up(ix,levs,3)                               ! input      u v T         at dt=0
      real, intent(out)    :: dudt(ix,levs,3)                             ! tendency   du/dt dv/dt  dT/dt    

      real, intent(in)     :: ahs_i(im, levs+1)
      real, intent(in)     :: rdelp(im,levs)  ! 1/dp layer thickness in KPa
      real, intent(in), dimension(im, levs+1) ::  vurho, murho        !

!
! Local variables.... Except "vurho" & "murho" all other VARS is defined fro, [1:levs]
!                     ac(1) = 0;     cc(levs) =0; 
!                     ec & dc ...alpha*Y(k+1) + beta(k) = Y(k) ... make only sense at
!                     mid-points where Y-[u,v,T] are defined in WAM
!       
      real ::  ac(levs),cc(levs),ec_i(levs),dc_i(levs)
      real ::  ynew(levs)
      real ::  coef_i(levs,2)
      real ::  hs_i(levs+1)
      real ::  partb_i(levs),  parta(levs,2)
      real ::  hold1, dtp1, hold2
      real :: sc_grav          !1.e-3*gravit
      integer j, k,i,kk,kk1
      integer, parameter :: Sw_dif_solver  = 0
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! OLD-WAM-1D-diff solver
!
!          k=levs             !,1,-1
!            hold1=1./(1.+ac(k)+cc(k)-cc(k)*ec_i(k+1))
!             hold1 = 1./(b[levs])
!            ec_i(k)=ac(k)*hold1  = ac/(1+ac+cc)
!            dc_i(k)=(cc(k)*0+ Y(levs))*hold1 = Y(levs)/(1+ac+cc)
!
!            Ynew(levs) = [ac/(1+ac+cc)]*Ynew(levs-1)+ Yold(levs)/(1+ac+cc)
!
! incorrect "DIFF-solver" at the Top LID
!             b*Ynew(levs) -a*Ynew(levs-1) = Yold(levs)
!===============================================================
! NEW-WAM-1D-diff solver
!
! set boundary conditions at level =1    (surface)
!                            level=levs  (top lid)
!
! std-BC  surface: "Dirichle" - same values if PBL is not active
!         top lid: " Neiman"  - zero vertical grafients or fluxes
! Dirichle
!         Y(1) = Ec(1)*Y(2)+Dc(1)
!         Ec(1) =0
!         Dc(1) = Y(1)
! Neiman
!         Y(N-1) = Ec(N-1)*Y(N) + Dc(N-1)
!         Ec(N-1) = 1.
!         Dc(N-1) = 0.
! Algorithm:
!           "upward-forward"    computaions of EC(2-levs) & DC(2-levs)
!           "downward-backward" computaions of Ynew & tendencies (Ynew-Yold)/dt
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      partb_i(1)=0.
      ac(1)=0.
      cc(levs)=0.

      dtp1=1./dtp
      if (Sw_dif_solver == 1) then

      ec_i(1)= 0. ! Ynew(1) = Yold(1) in the limit of the "weak" eddy+ molecular diffusion => dY/dt ~0. near the surface 
                  ! Zero vertical flux at the surface w/o PBL-diffusion is ?-able
!     dc_i(1)  = Yold
!
! for each longitude
!
      do i=1,im
 
!================================================================
!  solve tridiagonal problem with dynamical diff-n coeff-ts 
!                                 weighted by rho-density
!========================= can be done in the single height loop
!
! d(rho*U) = d([rho*Visc]*dU/dz)
!
!   expressions in code: dPK = p(k)- p(k+1)
!                        dPkm= p(k)- p(k-1)
!
!   A(K, 2:levs)   =dtime* 1./dpK*[rhoVu_K] *Pi(K)  /Hs(k)  /dPkm 
!   C(K, 1:levs-1) =dtime* 1./dpK*[rhoVu_kP]*Pi(K+1)/Hs(k+1)/dPk 
!
!   d(rho*cp*T) = d([rho*Cond]*dT/dz)
!
!    Factorization
!
!=================================================================
!    -A(k)*Y(k-1) + (1.+A(k)+C(k))*Y(k) -C(k)*Y(k+1) = F(K)=up(k)
!
!    E(k) = A(k)/[ B(k) -C(k)*E(K+1)], k= levs => 1 with step =-1
!
!    D(k) = [C(k)*D(K+1) + F(k)]/[B(k) -C(k)*E(K+1)]
!    Y(1) = D(1)
!    Y(K) = D(K) + E(K)*Y(K-1) k=2, levs, 1
!
! Verification of DISCRETE EQ-ns
!      Y[k-1].. Ac => 
!========================================================
       do k=1, levs
          hs_i(k)= ahs_i(i,k)
          coef_i(k,1)=vurho(i,k)
          coef_i(k,2)=murho(i,k)
!
!                                                 !parta(k)=dtp*grav(i,k)*.001/(prsi(i,k)-prsi(i,k+1))    ! delp - peredaet 1.e-3
          parta(k,1)=dtp*grav(i,k)*rdelp(i,k)     ! prsi(k) -prsi(k+1)
          parta(k,2)=parta(k,1)                   ! divided by Cp in "K-heat"
        enddo
        do kk=1,3
           kk1=kk/3+1                             ! if kk=1,2 kk1 =1  else kk=3 kk1=2
          do k=2,levs
          partb_i(k)=coef_i(k,kk1)*prsi(i,k)/hs_i(k)/
     &      (prsl(i,k-1)-prsl(i,k))
!
!  More accurate expressions with scale-height Hp-fix for Keddy
!        partb_i(k)=Vu_coef(k,kk1)*prsi(i,k)*prsi(i,k)/(prsl(i,k-1)-prsl(i,k)))/Hp/Hp
!
             ac(k)=parta(k,kk1)*partb_i(k)
          enddo
          do k=1,levs-1
            cc(k)=parta(k,kk1)*partb_i(k+1)
          enddo
!          do k=levs,1,-1                        Y(1)*(1+cc(1)) - cc(1)*Y(2)  =  Yold(1)
           dc_i(1) =  up(i,1,kk) / (1. +cc(1) )
           ec_i(1) =  cc(1)/(1. +cc(1) )
!
!           dc_i(1) =  up(i,1,kk)
!           ec_i(1) =  0.
!
          do k=2,levs-1
            hold1=1./(1.+ac(k)+cc(k)-ac(k)*ec_i(k-1))
            ec_i(k)=cc(k)*hold1 
            dc_i(k)=(ac(k)*dc_i(k-1)+up(i,k,kk))*hold1 
          enddo
!
! UBC
!
            k= levs
            hold1=1./(1.+ac(k)-ac(k)*ec_i(k-1))
            dc_i(k)=(ac(k)*dc_i(k-1)+up(i,k,kk))*hold1      
            ynew(k) = dc_i(k)
            dudt(i,K,kk)=(ynew(K)-up(i,K,kk))*dtp1
!
! backward factorization ...compute dc_i.....from lev 2 => top_lid
!
          do k=levs-1,1,-1
            ynew(k)=dc_i(k)+ec_i(k)*ynew(k+1)
            dudt(i,k,kk)=(ynew(k)-up(i,k,kk))*dtp1
          enddo
        enddo         !kk=1,2 (u,v)
!
!  kinetic energy visc. loss add to the temperature tendency due to energy conservation 
!
!!       do k=1,levs
!!        dudt(i,k,3)=dudt(i,k,3)
!!    &   -(up(i,k,1)*dudt(i,k,1)+up(i,k,2)*dudt(i,k,2))/cp(i,k)        
!       enddo
!
      enddo !i-horiz loop
!
      ELSE
!                           ! WAM-old if (Sw_dif_solver =/= 1) then
!
! set boundary
      partb_i(1)=0.
!      partb_i(levs+1)=0.
!      ec_i(levs+1)=0.
!      dc_i(levs+1)=0.

      ac(1)=0.
      cc(levs)=0.
      dtp1=1./dtp
!
! for each longitude
      do i=1, im
        do k=1,levs 
          hs_i(k)= ahs_i(i,k)          
          coef_i(k,1)=vurho(i,k)
          coef_i(k,2)=murho(i,k)
!          vurho(i,k) = coef_i(k,1)
!          murho(i,k) = coef_i(k,2)*cp1(k)
        enddo
 
! solve tridiagonal problem
        do k=1,levs
          parta(k,1)=dtp*grav(i,k)/(prsi(i,k)-prsi(i,k+1))   ! *0.001 due kPa
          parta(k,2)=parta(k,1)
        enddo
        do kk=1,3
          kk1=kk/3+1
          do k=2,levs
            partb_i(k)=coef_i(k,kk1)*prsi(i,k)/                         
     &       (hs_i(k)*(prsl(i,k-1)-prsl(i,k)))
            ac(k)=parta(k,kk1)*partb_i(k)
          enddo
          do k=1,levs-1
            cc(k)=parta(k,kk1)*partb_i(k+1)
          enddo
          do k=levs,1,-1
            hold1=1./(1.+ac(k)+cc(k)-cc(k)*ec_i(k+1))
            ec_i(k)=ac(k)*hold1 
            dc_i(k)=(cc(k)*dc_i(k+1)+up(i,k,kk))*hold1 
          enddo
          dudt(i,1,kk)=(dc_i(1)-up(i,1,kk))*dtp1
! recompute dc_i
          do k=2,levs
            dc_i(k)=dc_i(k)+ec_i(k)*dc_i(k-1)
            dudt(i,k,kk)=(dc_i(k)-up(i,k,kk))*dtp1
          enddo
        enddo  !kk
!        do k=1,levs
!          dt6dt(i,k,5)=dudt(i,k,3)
!        enddo
! u v changes add to temperature tendency due to energy conservation 
!       do k=1,levs
!         dudt(i,k,3)=dudt(i,k,3)-cp1(k)*(up(i,k,1)*dudt(i,k,1)         
!    &    +up(i,k,2)*dudt(i,k,2))
!         dt6dt(i,k,6)= -1.*cp1(k)*(up(i,k,1)*dudt(i,k,1)               
!    &    +up(i,k,2)*dudt(i,k,2))
!       enddo
      enddo !i
      ENDIF 
      return
      end subroutine idea_visc_conduct
!
      subroutine idea_vdiff_tracers(im,ix,levs,grav,prsi,prsl,
     & prslk, exner, dtp, adr, adt, ntrac, ktemp, ktrac, ahs_i, rdelp)
!-----------------------------------------------------------------------
!
! solve eddy diffusion of major WAM-tracers + Potential temp-re  caused by  
!       "eddy" diffusion and thermal conductivity
!      
!       to make this transform in GFS you need to "mmr" => 1/cm3
!-----------------------------------------------------------------------
 
      use machine, only : kind_phys
      use idea_composition
      implicit none
      
!
! Argument
!
      integer, intent(in) :: ntrac ! number of tracers for Keddy in vertical
      integer, intent(in) :: im    ! number of data points in up,dudt(first dim)
      integer, intent(in) :: ix    ! max data points in fields
      integer, intent(in) :: levs  ! number of pressure levels
      real,    intent(in) :: dtp   ! time step in second
!
! updated variables tracers and potential temperature (kinetic temp-re)
!
      real, intent(inout)    :: adr(ix,levs, ntrac) ! tracers in "mmr" solver of equations in "vmr" 
      real, intent(inout)    :: adt(ix,levs)        ! kin-c temp-re K, solver for potential temp-re 
      real, intent(in)    :: ktemp(im,levs+1)      ! eddy heat conductivity, rho*m2/s
      real, intent(in)    :: ktrac(im,levs+1)      ! eddy tracer diffusion,  rho*m2/s
      real, intent(in)    :: ahs_i(im,levs+1)      ! eddy tracer diffusion,  rho*m2/s
      real, intent(in)    :: rdelp(im,levs+1)   
      real, intent(in)    :: exner(ix,levs)        ! PT = exner*T
      real, intent(in)    :: prslk(ix,levs)        ! T = prslk*PT

      real, intent(in)    :: prsi(ix,levs+1) ! interface pressure in KPa
      real, intent(in)    :: prsl(ix,levs)   ! layer pressure in KPa
      real, intent(in)    :: grav(ix,levs)   ! (m/s2)
 
! locals
      real                :: adtpt(levs), atrac(levs)
! 
! Local variables
      
      real ac(levs),cc(levs),ec_i(levs+1),dc_i(levs+1)
      real acT(levs),ccT(levs),ecT_i(levs+1),dcT_i(levs+1)

      real t_i(levs+1),hs_i(levs+1)
      real partb(levs+1),parta(levs),hold1(levs)
      real partbT(levs+1),partaT(levs)
      real Thold1, dtp1, hold2, sc_grav 
      real ratpdt
      integer j, k,i,kk,kk1
!
! set boundary
      partb(1)=0.
      partb(levs+1)=0.
      ec_i(levs+1)=0.
      dc_i(levs+1)=0.
      ac(1)=0.
      cc(levs)=0.

      partbT(1)=0.
      partbT(levs+1)=0.
      ecT_i(levs+1)=0.
      dcT_i(levs+1)=0.
      acT(1)=0.
      ccT(levs)=0.


      dtp1=1./dtp
!
! for each longitude
!
      do i=1,im
! prepare coeffs of the  tridiagonal problem
        do k=1,levs
          sc_grav = grav(i,k)               !*.001
          parta(k)=dtp*sc_grav*rdelp(i,k)   ! delp - peredat' agam 1.e-3/factor
          partat(k)=parta(k)                ! *cp1(k)      
          hs_i(k)= ahs_i(i,k)
!
          adtpt(k) = adt(i,k)*exner(i,k)  
!
        enddo
!
! compute ac and cc for PT and tracers .......put in "adr(....4:5)" - updated O & O2
! Realistic influence of Ktr_eddy => H2O, O3, O, O2  , no CLW
  
          do k=2,levs
            ratpdt = prsi(i,k)/(prsl(i,k-1)-prsl(i,k))/hs_i(k)                     
            partb(k) =ktrac(i,k)*ratpdt
            partbt(k)=ktemp(i,k)*ratpdt
!
            ac(k) =parta(k)*partb(k)
            act(k)=partaT(k)*partbT(k)
          enddo
          do k=1,levs-1
            cc(k) =parta(k)*partb(k+1)
            ccT(k)=partaT(k)*partbT(k+1)
          enddo
!
          do k=levs,1,-1
            hold1(k)=1./(1.+ac(k)+cc(k)-cc(k)*ec_i(k+1))           
            ec_i(k)=ac(k)*hold1(k)
          enddo
!
        do j=1, ntrac             !==============ntrac
           do k=levs,1,-1        
            dc_i(k)=(cc(k)*dc_i(k+1)+adr(i,k,j))*hold1(k) 
          enddo
          if (dc_i(1).lt.0) dc_i(1) =0.0
          adr(i,1,j)=dc_i(1)
! backward step for  dc_i >= 0. (non-negative tracers)
          do k=2,levs
            dc_i(k)=dc_i(k)+ec_i(k)*dc_i(k-1)
            if (dc_i(k).lt.0) dc_i(k) =0.0
            adr(i,k,j) = dc_i(k)
          enddo
!
        enddo                      !==============ntrac
!
!  put eddy-conduction on PT .....Temp-re > 90K
!
          do k=levs,1,-1
            thold1=1./(1.+act(k)+ccT(k)-ccT(k)*ecT_i(k+1))
            ecT_i(k)=acT(k)*Thold1 
!
            dcT_i(k)=(ccT(k)*dcT_i(k+1)+adtPT(k))*Thold1 
          enddo
!
! Temp-re > 90K
!
           adtPT(1)=dcT_i(1)      !no changes at the surface due GW-eddies
!
!  backward steps
!
          do k=2,levs
            dcT_i(k)=dcT_i(k)+ecT_i(k)*dcT_i(k-1)
          enddo
!  back to kinetic temp-re
        do k=1,levs
           adt(i,k)=dcT_i(k)*prslk(i,k)
        enddo
!
      enddo         !i
      return
      end subroutine idea_vdiff_tracers
!
! WSH-PBL TRIDIAGs
!
      subroutine gwunif_wam_viscond(im,ix,levs,
     &     grav,prsi,prsl, adt, o_n,o2_n, n2_n, cp,
     &     vurho, murho, ahs_i, amol_i, rho_i)
!     call gwunif_wam_viscond(im,ix,levs,grav,prsi,prsl,adt,
!     &   o_n,o2_n, n2_n, cp, vurho, murho, ahs_i, amol_i, rho_i)
!    
!
! compute vurho, murho  - coefficients of dynamical "molecular" viscosity/conduction
!                         temp-re & comp-n dependnet in the thermosphere
!
!con_rgas   =8.314472 ! molar gas constant  (J/mol/K) ..
!            WRONG value for P= R/mu*(dens*temp)
!
      use physcons,  rgas=>con_rgas, amo2=>con_amo2
      use machine, only : kind_phys
      use idea_composition, only : amo, amn2

      IMPLICIT NONE
! Argument
!      integer, intent(in) :: me              ! PE-id number
      integer, intent(in) :: im              ! number of columns actual at ix(me)=im
      integer, intent(in) :: ix              ! max-reserved  points/columns
      integer, intent(in) :: levs            ! number of pressure levels

   
      real, intent(in) ::   PRSI(IX,levs+1)   ! interface pressure in Pa
      real, intent(in) ::   PRSL(IX,levs)    ! 1/pressure thickness
      real, intent(in)    :: grav(ix,levs)   ! (m/s2)
      real, intent(in)    :: o_n(ix,levs)    ! number density (/m3) of O
      real, intent(in)    :: o2_n(ix,levs)   ! number density (/m3) of O2
      real, intent(in)    :: n2_n(ix,levs)   ! number density (/m3) of N2
!
      real, intent(in)    :: adt(ix,levs)    ! temperature, K
      real, intent(in)    ::  cp(ix,levs)       
!
      real, intent(out), dimension(ix, levs+1)::  vurho, murho
      real, intent(out), dimension(ix, levs+1)::  ahs_i, amol_i, rho_i
!locals
! define some constants
      real (kind=kind_phys), parameter:: muo =3.9e-7    ! viscosity coefficient of O (kg/m/s) 
      real (kind=kind_phys), parameter:: muo2=4.03e-7   ! viscosity coefficient of O2 (kg/m/s) 
      real (kind=kind_phys), parameter:: mun2=3.43e-7   ! viscosity coefficient of N2 (kg/m/s) 
!                                     
      real (kind=kind_phys), parameter:: lao =75.9e-5   ! thermal conductivity coeff of O (W/m/K)
      real (kind=kind_phys), parameter:: lao2=56.e-5    ! thermal conductivity coefficient of O2(W/m/K)
      real (kind=kind_phys), parameter:: lan2=56.e-5    ! thermal conductivity coefficient of N2(W/m/K)
!                                                      
      real (kind=kind_phys), parameter:: cpo =2.5       !specific heats of o   5/2
      real (kind=kind_phys), parameter:: cpo2=3.5       !specific heats of o2  7/2
      real (kind=kind_phys), parameter:: cpn2=3.5       !specific heats of n2  7/2
!
      real, dimension(levs+1)   ::  o_ni, o2_ni, n2_ni, t_i, cp1, ma_i
!
      real, dimension(levs+1)   ::  mu_i               ! mumb-dens weighted viscosity kg/m/s
      real, dimension(levs+1)   ::  la_i               ! mumb-dens weighted comduction J/m/s/K
      integer ::   i, k
      real    :: hold1, hold2, tdep
      real    :: R8314,  Rmumol                        ! Adequate R-gas constants
      R8314 = 1000.*RGAS
!
      do i=1,im
! get compositions at interface pressure levels
        o_ni(1)=o_n(i,1)
        o2_ni(1)=o2_n(i,1)
        n2_ni(1)=n2_n(i,1)
!
        do k=2,levs
          o_ni(k)=(o_n(i,k-1)+o_n(i,k))*.5
          o2_ni(k)=(o2_n(i,k-1)+o2_n(i,k))*.5
          n2_ni(k)=(n2_n(i,k-1)+n2_n(i,k))*.5
        enddo
! calculate coefficients of mu, lambda, 1./cp, at interface pressure
        do k=1,levs
          hold1=1./(o_ni(k)+o2_ni(k)+n2_ni(k))
!
          mu_i(k)=(o_ni(k)*muo+o2_ni(k)*muo2+n2_ni(k)*mun2)*hold1
          la_i(k)=(o_ni(k)*lao+o2_ni(k)*lao2+n2_ni(k)*lan2)*hold1
!
          hold2=o_ni(k)*amo+o2_ni(k)*amo2+n2_ni(k)*amn2
          amol_i(i,k)=hold2*hold1*grav(i,k)          
          ma_i(k)=hold2*hold1  
          
        enddo
         amol_i(i,levs+1)=  amol_i(i,levs) 
! at layer
        do k=1,levs-1
!          cp1(k)=.5/(cp(i,k)+cp(i,k+1))
         cp1(k)=1./cp(i,k)
        enddo
          cp1(levs) = 1./cp(i,levs)
          cp1(levs+1) = cp1(levs)
! calculate temp in interface pressure levels
! calculate scale height
        t_i(1)=adt(i,1)
        t_i(levs+1)=adt(i,levs)
        do k=2,levs
          t_i(k)=(adt(i,k-1)+adt(i,k))*.5
        enddo
! now use t_i**0.69
! calculate viscosity put in vurho(i,k)
! calculate thermal conductivity put in murho(i,k)
        do k=1,levs 
          tdep=t_i(k)**(0.69)
          Rmumol = R8314/amol_i(i,k)  
!          ahs_i(i,k) = Rmumol*t_i(k)   ! (RT/(Mu*Grav))
!identical
          ahs_i(i,k)=1000.*rgas*t_i(k)/(ma_i(k)*grav(i,k))
!
          vurho(i,k) = mu_i(k)*tdep
          murho(i,k) = la_i(k)*tdep*cp1(k)
          rho_i(i,k) = prsi(i,k)/ahs_i(i,k)/grav(i,k)
        enddo

          k = levs+1
          ahs_i(i,k) = ahs_i(i,k-1)        !Rmumol*t_i(k)   ! (RT/(Mu*Grav))
          vurho(i,k) = vurho(i,k-1)
          murho(i,k) = murho(i,k-1)
          rho_i(i,k) = prsi(i,k)/ahs_i(i,k-1)/grav(i,k-1)

      ENDDO ! im
!
      end  subroutine gwunif_wam_viscond
!
      subroutine tridiag_uv(l,n,cl,cm,cu,r1,r2,au,a1,a2)
!
! Tridiag Form:  cl*Y(k-1) +cm*Y(k) +cu*Y(k+1) = [r1, r2]
!                fk = 1/denom , denom = cm(i,k)-cl(i,k)*au(i,k-1)
!                cm = 1.-cl-cm
!                au(k) = cu(k)/denom
!         [a1, a2](k) = ([r1, r2] -cl(k)*[a1,a2](k-1))/denom
!
      use machine     , only : kind_phys
      implicit none
      integer             k,n,l,i
      real(kind=kind_phys) fk

      real(kind=kind_phys) cl(l,2:n),cm(l,n),cu(l,n-1),r1(l,n),r2(l,n),
     &          au(l,n-1),a1(l,n),a2(l,n)
!-----------------------------------------------------------------------
      do i=1,l
        fk      = 1./cm(i,1)
        au(i,1) = fk*cu(i,1)
        a1(i,1) = fk*r1(i,1)
        a2(i,1) = fk*r2(i,1)
      enddo
      do k=2,n-1
        do i=1,l
          fk      = 1./(cm(i,k)-cl(i,k)*au(i,k-1))
          au(i,k) = fk*cu(i,k)
          a1(i,k) = fk*(r1(i,k)-cl(i,k)*a1(i,k-1))
          a2(i,k) = fk*(r2(i,k)-cl(i,k)*a2(i,k-1))
        enddo
      enddo
!
! UBC
!
      do i=1,l
        fk      = 1./(cm(i,n)-cl(i,n)*au(i,n-1))
        a1(i,n) = fk*(r1(i,n)-cl(i,n)*a1(i,n-1))
        a2(i,n) = fk*(r2(i,n)-cl(i,n)*a2(i,n-1))
      enddo
!
! from top to the surface update "U-V"-values
! for system with -a*Y+bY-c*Y = r
! Ynew(k) = F(k) + E(k)*Ynew(k+1)
!         F(k)  = Denom * ( RHS(k) + a(k)*F(k-1))
!         E(k)  = C(k)*Denom
!         Denom = 1./(cm - a*E(k-1))
! to "save" memory Ynew === F on the backward progon
!
      do k=n-1,1,-1
        do i=1,l
          a1(i,k) = a1(i,k)-au(i,k)*a1(i,k+1)
          a2(i,k) = a2(i,k)-au(i,k)*a2(i,k+1)
        enddo
      enddo
!-----------------------------------------------------------------------
      return
      end subroutine tridiag_uv
!
!
!-----------------------------------------------------------------------
      subroutine tridiag_ptntrac(l,n,nt,cl,cm,cu,r1,r2,au,a1,a2)
!
! for tracers + PT
!
      use machine     , only : kind_phys
      implicit none
      integer             is,k,kk,n,nt,l,i
      real(kind=kind_phys) fk(l)
!
      real(kind=kind_phys) cl(l,2:n), cm(l,n), cu(l,n-1),
     &                     r1(l,n),   r2(l,n*nt),
     &                     au(l,n-1), a1(l,n), a2(l,n*nt),
     &                     fkk(l,2:n-1)
!-----------------------------------------------------------------------
      do i=1,l
        fk(i)   = 1./cm(i,1)
        au(i,1) = fk(i)*cu(i,1)
        a1(i,1) = fk(i)*r1(i,1)
      enddo
      do k = 1, nt
        is = (k-1) * n
        do i = 1, l
          a2(i,1+is) = fk(i) * r2(i,1+is)
        enddo
      enddo
      do k=2,n-1
        do i=1,l
          fkk(i,k) = 1./(cm(i,k)-cl(i,k)*au(i,k-1))
          au(i,k)  = fkk(i,k)*cu(i,k)
          a1(i,k)  = fkk(i,k)*(r1(i,k)-cl(i,k)*a1(i,k-1))
        enddo
      enddo
      do kk = 1, nt
        is = (kk-1) * n
        do k=2,n-1
          do i=1,l
            a2(i,k+is) = fkk(i,k)*(r2(i,k+is)-cl(i,k)*a2(i,k+is-1))
          enddo
        enddo
      enddo
      do i=1,l
        fk(i)   = 1./(cm(i,n)-cl(i,n)*au(i,n-1))
        a1(i,n) = fk(i)*(r1(i,n)-cl(i,n)*a1(i,n-1))
      enddo
      do k = 1, nt
        is = (k-1) * n
        do i = 1, l
          a2(i,n+is) = fk(i)*(r2(i,n+is)-cl(i,n)*a2(i,n+is-1))
        enddo
      enddo
      do k=n-1,1,-1
        do i=1,l
          a1(i,k) = a1(i,k) - au(i,k)*a1(i,k+1)
        enddo
      enddo
      do kk = 1, nt
        is = (kk-1) * n
        do k=n-1,1,-1
          do i=1,l
            a2(i,k+is) = a2(i,k+is) - au(i,k)*a2(i,k+is+1)
          enddo
        enddo
      enddo
!-----------------------------------------------------------------------
      return
      end  subroutine tridiag_ptntrac  
