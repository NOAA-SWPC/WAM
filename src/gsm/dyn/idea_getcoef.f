      subroutine idea_getcoef(levs,ntrac,q,plyr,mur,lam,d12)

! calculate global avg viscosity, conductivity, and diffusion coeffs
! Apr 06 2012    Henry Juang, initial implement
! Dec    2020    VAY, correction of global coefficients and moluular mass if 
!                WAM-IC is not balanced (unrealistic H2O)
!                we apply calculations only for "major" species O+O2+N2=1
!                to avoid "bad" values of O3 and H2O above 80 km
! to do for testing  "Cp*T" vs "T" Eul-code perfoemance
! gfs_dyn_tracer_const.f  define # of tracers ri(0:max_num_tracer=50???) & cpi(0:max_num_tracer)
!                       it would be better to put WAM-values w/o reading them from "exglobal"
!
!  computation of am(levs) is here  for initialisation of global AMOL and "mur,lam,d12"
!                 it can be "bad" for the WAM spin-up with different SOLAR activity 
!
!  am(levs): idea_composition, only: am=>amgm define inn-t radiation ans spectral disispation
!
      use gfs_dyn_physcons,  amo2 => con_amo2,
     &               avgd => con_avgd, p0 => con_p0,
     &               amh2o =>con_amw, amo3 =>con_amo3
!hmhj gfs_dyn_coordinate_def                            ! cpi-definition
      use gfs_dyn_tracer_const                          ! cpi-definition
      use idea_composition, only: am=>amgm
      implicit none

! subroutine params

      integer, intent(in) :: levs   ! number of pressure model layers
      integer, intent(in) :: ntrac   ! number of tracers
      real,    intent(in) :: q(levs,0:ntrac) ! temp*cp,h2o,o3,cld,o,o2 tracer
      real,    intent(in) :: plyr(levs) ! global mean presure at model layer
      real,    intent(out):: mur(levs) ! mu/rho (m2/s)
      real,    intent(out):: lam(levs) ! lambda/rho/cp (m2/s)
      real,    intent(out):: d12(levs) ! d12/n

! local params

      real amo,amn2,muo,muo2,mun2,lao,lao2,lan2
      parameter (amo=15.9994, amn2=28.013)               !g/mol
      parameter (muo=3.9e-7, muo2=4.03e-7, mun2=3.43e-7) !kg/m/s
      parameter (lao=75.9e-5, lao2=56.e-5, lan2=56.e-5)  !kg/m/s         
      real, parameter:: bz=1.3806505e-23                 ! Boltzmann constant 
      real, parameter:: s12=0.774,a12=9.69e18            ! O-O2 diffusion params
!vay-2021      
      real, parameter:: dked_min =0.1                  ! minval for eddy-horizontal additional stab-n for HR-runs
      real, parameter:: amu_dry= 28.9647               !g/mol    
! local variables

      real hold1,la,rho,o_n,o2_n,cp,t69,mu,n,n2_n,q_n2,tem,rn
      real rrn, rbz, amol, nair, nsum
      integer k
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!      print *,'in idea_getcoef,cp=',cpi(1:5)

        rrn=1.e-3/(avgd)
	rbz = 1./bz
        do k=1,levs
         q_n2=1.- q(k,4)-q(k,5)	
	 cp= cpi(4)*q(k,4)+cpi(5)*q(k,5)+ cpi(0)*q_n2	
	 tem= q(k,0)/cp
	 amol=1./(q(k,4)/amo+q(k,5)/amo2+q_n2/amn2)        ! g/mol
         nair=rbz*plyr(k)/tem                              ! kg/m3
         rn=avgd*bz
         rho=amol*nair *rrn                                ! kg/m3
         o_n=q(k,4)*amol/amo                               ! 1/m3    
         o2_n=q(k,5)*amol/amo2         
         n2_n=q_n2*amol/amn2   
	 nsum =1./(o_n + o2_n + n2_n)      
!
! kinematic diffusion
!
         d12(k)=a12* tem**(s12)/nair  + dked_min                     !d12 
!
! kinematic viscosity and heat conductivity see idea_dissipation.f
!
         mu=(o_n*muo+o2_n*muo2+n2_n*mun2)*nsum
         la=(o_n*lao+o2_n*lao2+n2_n*lan2)*nsum
	 
         t69=tem**(0.69)/rho

         mur(k)=mu*t69    + dked_min
         lam(k)=la*t69/cp + dked_min
	 amol = min( amol, amu_dry)                  ! to avoid amol > amu_dry 
	 am(k)= max( amol, amo)  	 
	enddo 

        return
!============================ OLD code ==============      
      do k=1,levs
!         am(k) = 28.9647                   !g/mol

! get global mean pressure

!VAF-2021 check that plyr in Pa    plyr=(ak5(k)+ak5(k+1)+p0*.5*(bk5(k)+bk5(k+1)))

! get n2 kg/kg

       if( ntrac.eq.5 ) then
         q_n2=1.-q(k,1)-q(k,2)-q(k,4)-q(k,5)
       else
         q_n2=1.
         do n=1,ntrac
           if(cpi(n).ne.0.0) q_n2=q_n2-q(k,n)
         enddo
       endif

! get cp, recover temp from enthalpy

       if( ntrac.eq.5 ) then
         cp=cpi(1)*q(k,1)+cpi(2)*q(k,2)+cpi(4)*q(k,4)+cpi(5)*q(k,5)+
     &        cpi(0)*q_n2
       else
         cp=cpi(0)*q_n2
         do n=1,ntrac
           if(cpi(n).ne.0.0) cp=cp+cpi(n)*q(k,n)
         enddo
       endif
       if(abs(cp)<1.0e-10) then
        print *,'in idea_getcoef,k=',k,'q(k1)=',q(k,1),q(k,2),
     &   q(k,3),q(k,4),q(k,5),'q_n2=',q_n2
       endif
         tem=q(k,0)/cp    ! K

! get molecular mass, number densities

       if( ntrac.eq.5 ) then
         am(k)=1./(q(k,4)/amo+q(k,5)/amo2+q_n2/amn2+q(k,1)/amh2o+
     &        q(k,2)/amo3)                      ! g/mol
         n=plyr(k)/tem                          ! 1/m3/bz
         rn=avgd*bz
         rho=1e-3*am(k)*n/rn                     ! kg/m3
         o_n=q(k,4)*am(k)/amo            
         o2_n=q(k,5)*am(k)/amo2         
         n2_n=q_n2*am(k)/amn2          
         hold1=tem/plyr(k)
!! get d12

         d12(k)=(a12*bz)*tem**(s12)*hold1      !d12 

! calculate mass nu=mu/rho,lamda/(rho*cp)

         mu=o_n*muo+o2_n*muo2+n2_n*mun2
         la=o_n*lao+o2_n*lao2+n2_n*lan2

! now use tem**0.69

         t69=tem**(0.69)

         mur(k)=mu*t69/rho
         lam(k)=la*t69/rho/cp

!      print*,'ntrac=5,k=',k,plyr(k),mur(k),lam(k),d12(k),
!     &  'q(k,0)=',q(k,0),'q(k,1)=',q(k,1),'q(k,2)=',q(k,2),
!     &  'q(k,3)=',q(k,3),'q(k,4)=',q(k,4),'q(k,5)=',q(k,5),
!     &  'cp=',cp
!
!------------------------------------------------------
       else  ! assume ntrac=3
!       
! ABZAATCH of Jun Wang N2 +H2O !!!!  VAY-2021 
!       
         am(k)=1./(q_n2/amn2+q(k,1)/amh2o+q(k,2)/amo3)
         o_n=0.0
         o2_n=0.0
!jw
         n=plyr(k)/tem                          ! 1/m3/bz
         rn=bz*avgd                          ! 1/m3
         rho=1e-3*am(k)*n/rn                   ! kg/m3
         n2_n=q_n2*am(k)/amn2              !
         hold1=tem/plyr(k)

!       print *,'ntrac=3,k am rho n2_n ',k,am(k),rho,n2_n,n,rn,tem,      &
!     &  'q=',q_n2,q(k,1:ntrac),cp,'cpi=',cpi(0:ntrac),q(k,0)

! get d12

         d12(k)=(a12*bz)*tem**(s12)*hold1   

! calculate mass nu=mu/rho,lamda/(rho*cp)

         mu=n2_n*mun2
         la=n2_n*lan2

! now use tem**0.69

         t69=tem**(0.69)
         
         mur(k)=mu*t69/rho
         lam(k)=la*t69/rho/cp

       endif

      enddo
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

      end subroutine idea_getcoef

