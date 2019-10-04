!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        module module_dynamics_routines
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
use mpi
use module_control,only:klog,kint,kfpt,kdbl,timef
 
use module_my_domain_specs

use module_clocktimes,only : timers
 
use module_dm_parallel,only : looplimits

use module_exchange,only: halo_exch
use module_fltbnds,only: polehn,polewn,swaphn,swapwn
use module_constants

private
  
public :: adv1,adv2 &
,cdwdt,cdzdt,ddamp,dht &
,hdiff &
,mono,pdtsdt,pgforce &
,updates,updatet,updateuv &
,vsound,vtoa

integer(kind=kint),save :: &
 jstart &
,jstop 

real(kind=kdbl) :: &
 btim
 
!-----------------
#ifdef ENABLE_SMP
!-----------------
integer(kind=kint) :: &
 nth &
,omp_get_num_threads &
,omp_get_thread_num &
,tid
 
external omp_get_num_threads,omp_get_thread_num
!------
#endif
!------
 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
                        subroutine pgforce &
(first,global,restart &
,lm &
,dt,ntimestep,rdyv &
,dsg2,pdsg1 &
,rdxv,wpdar &
,fis,pd &
,t,q,cw &
,pint &
,rtop &
!---temporary arguments-------------------------------------------------
,div &
,pcne,pcnw &
,pcx,pcy &
,tcu,tcv)
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
implicit none
!-----------------------------------------------------------------------
integer(kind=kint),parameter:: &
 jw=01                       ! rows w/o correction next to poles

real(kind=kfpt),parameter:: &
 cfc=1.533 &                 ! adams-bashforth positioning in time
,bfc=1.-cfc                  ! adams bashforth positioning in time
!-----------------------------------------------------------------------
logical(kind=klog),intent(in):: &
 first &                     ! first pass
,global &                    ! global domain
,restart                     ! restart case

integer(kind=kint),intent(in):: &
 lm &                        ! total # of levels
,ntimestep                   ! the current timestep

real(kind=kfpt),intent(in):: &
 dt &                        ! dynamics time step
,rdyv                        ! 1/deltay

real(kind=kfpt),dimension(1:lm),intent(in):: &
 dsg2 &                      ! delta sigmas
,pdsg1                       ! delta pressures

real(kind=kfpt),dimension(jds:jde),intent(in):: &
 rdxv &                      ! 1/deltax
,wpdar                       ! divergence correction weight

real(kind=kfpt),dimension(ims:ime,jms:jme),intent(in):: &
 fis &                       ! surface geopotential
,pd                          ! sigma range pressure difference

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(in):: &
 t &                         ! temperature
,q &                         ! specific humidity
,cw                          ! condensate

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm+1),intent(in):: &
 pint                        ! pressure at interfaces

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(out):: &
 rtop                        ! RT/p

!-----------------------------------------------------------------------
!---temporary arguments-------------------------------------------------
!-----------------------------------------------------------------------
real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(out) :: &
 div &                       ! horizontal mass divergence
,pcne &                      ! second term of pgf, ne direction
,pcnw &                      ! second term of pgf, nw direction
,pcx &                       ! second term of pgf, x direction
,pcy                         ! second term of pgf, y direction

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(inout) :: &
 tcu &                       ! time change of u
,tcv                         ! time change of v
!-----------------------------------------------------------------------
!--local variables------------------------------------------------------
!-----------------------------------------------------------------------
integer(kind=kint):: &
 i &                         ! index in x direction
,icl &                       !
,ich &                       !
,ip &                        !
,j &                         ! index in y direction
,jcl &                       ! lower bound for no divergence correction
,jch &                       ! upper bound for no divergence correction
,jp &                        !
,l                           ! index in p direction

real(kind=kfpt):: &
 apd &                       ! hydrostatic pressure difference at the point
,apelp &                     ! pressure at the point
,dfip &                      ! delta phi
,dfdp &                      ! dfi/dp
,fiup &                      ! geopotential at the upper interface
,ppne &                      ! first term of pgf, ne direction
,ppnw &                      ! first term of pgf, nw direction
,ppx &                       ! first term of pgf, x direction
,ppy &                       ! first term of pgf, y direction
,rdu &                       !
,rdv &                       !
,rpdp &                      !
,wprp                        ! divergence modification weight at the point

real(kind=kfpt),dimension(its_h2:ite_h2,jts_h2:jte_h2):: &
 apel &                      ! scratch, pressure in the middle of the layer
,dfi &                       ! scratch, delta phi
,filo &                      ! scratch, geopotential at lower interface
,fim                         ! scratch, geopotential in the middle of the layer

real(kind=kfpt),dimension(its_h2:ite_h2,jts_h2:jte_h2,1:lm):: &
 apel_3d &                   ! scratch, 3d copy of pressure in the middle of the layer
,fim_3d  &                   ! scratch, 3d copy of geopotential in the middle of the layer
,dfi_3d                      ! scratch, 3d copy of delta phi

real(kind=kfpt),dimension(its_b1:ite_h2,jts_b1:jte_h2):: &
 pgne &                      ! scratch, pgf, ne direction
,pgnw                        ! scratch, pgf, nw direction

real(kind=kfpt),dimension(its_b1:ite_h2,jts:jte_h2):: &
 pgx                         ! scratch, pgf, x direction

real(kind=kfpt),dimension(its:ite_h2,jts_b1:jte_h2):: &
 pgy                         ! scratch, pgf, y direction
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------

      do j=jts_h2,jte_h2
        do i=its_h2,ite_h2
          filo(i,j)=fis(i,j)
        enddo
      enddo
!-----------------------------------------------------------------------
!---vertical grand loop-------------------------------------------------
!-----------------------------------------------------------------------
!-----------------
#ifdef ENABLE_SMP
!-----------------
!.......................................................................
!$omp parallel private(apelp,dfdp,dfip,fiup,i,j,jstart,jstop,l,nth,tid)
!.......................................................................
      nth = omp_get_num_threads()
      tid = omp_get_thread_num()
      call looplimits(tid, nth, jts_h2, jte_h2, jstart, jstop)
!-----------------
#else
!-----------------
      jstart = jts_h2
      jstop = jte_h2
!-----------------
#endif
!-----------------
!-----------------------------------------------------------------------
      do l=lm,1,-1
!-----------------------------------------------------------------------
        do j=jstart,jstop
          do i=its_h2,ite_h2
            apelp=(pint(i,j,l)+pint(i,j,l+1))*0.5
            apel_3d(i,j,l)=apelp
            dfdp=(q(i,j,l)*0.608+(1.-cw(i,j,l)))*t(i,j,l)*r/apelp
            dfip=dfdp*(dsg2(l)*pd(i,j)+pdsg1(l))
            rtop(i,j,l)=dfdp
            fiup=filo(i,j)+dfip
            dfi_3d(i,j,l)=dfip*0.5
            fim_3d(i,j,l)=(filo(i,j)+fiup)*0.5
            filo(i,j)=fiup
          enddo
        enddo
      enddo
!-----------------
#ifdef ENABLE_SMP
!-----------------
!.......................................................................
!$omp end parallel
!.......................................................................
!-----------------
#endif
!-----------------
!.......................................................................
!$omp parallel do &
!$omp private (apd,apel,dfi,fim,i,j,jch,jcl,l,pgne,pgnw,pgx,pgy, &
!$omp          ppne,ppnw,ppx,ppy,rdu,rdv,rpdp,wprp,icl,ich,ip,jp)
!.......................................................................
!-----------------------------------------------------------------------
!---vertical grand loop-------------------------------------------------
!-----------------------------------------------------------------------
      vertical_loop: do l=lm,1,-1
        do j=jts_h1,jte_h2
          do i=its_h1,ite_h2
           apel(i,j)= apel_3d(i,j,l)
           fim(i,j) = fim_3d(i,j,l)
           dfi(i,j) = dfi_3d(i,j,l)
          enddo
        enddo
!-----------------------------------------------------------------------
!---pressure gradient force components, behind h points-----------------
!-----------------------------------------------------------------------
        do j=jts,jte_h2
          do i=its_b1,ite_h2
            ppx=(fim(i,j)-fim(i-1,j)) &
               *(pint(i  ,j ,l+1)+pint(i-1,j  ,l+1) &
                -pint(i  ,j ,l  )-pint(i-1,j  ,l  ))*0.5
            pcx(i,j,l)=(dfi (i-1,j  )+dfi (i  ,j  )) &
                      *(apel(i  ,j  )-apel(i-1,j  ))
            pgx(i,j)=ppx+pcx(i,j,l)
          enddo
        enddo
!
        do j=jts_b1,jte_h2
          do i=its,ite_h2
            ppy=(fim(i,j)-fim(i,j-1)) &
               *(pint(i ,j ,l+1)+pint(i  ,j-1,l+1) &
                -pint(i ,j ,l  )-pint(i  ,j-1,l  ))*0.5
            pcy(i,j,l)=(dfi (i  ,j-1)+dfi (i  ,j  )) &
                      *(apel(i  ,j  )-apel(i  ,j-1))
            pgy(i,j)=ppy+pcy(i,j,l)
          enddo
        enddo
!
        do j=jts_b1,jte_h2
          do i=its_b1,ite_h2
            ppne=(fim(i,j)-fim(i-1,j-1)) &
                *(pint(i-1,j-1,l+1)+pint(i  ,j  ,l+1) &
                 -pint(i-1,j-1,l  )-pint(i  ,j  ,l  ))*0.5
            ppnw=(fim(i-1,j)-fim(i,j-1)) &
                *(pint(i-1,j  ,l+1)+pint(i  ,j-1,l+1) &
                 -pint(i-1,j  ,l  )-pint(i  ,j-1,l  ))*0.5
            pcne(i,j,l)=(dfi (i-1,j-1)+dfi (i  ,j  )) &
                       *(apel(i  ,j  )-apel(i-1,j-1))
            pcnw(i,j,l)=(dfi (i  ,j-1)+dfi (i-1,j  )) &
                       *(apel(i-1,j  )-apel(i  ,j-1))
            pgne(i,j)=ppne+pcne(i,j,l)
            pgnw(i,j)=ppnw+pcnw(i,j,l)
          enddo
        enddo
!-----------------------------------------------------------------------
!---divergence correction, janjic 1974 mwr, 1979 beitrage---------------
!-----------------------------------------------------------------------
        if(global) then
!-----------------------------------------------------------------------      
          jcl=jds
          jch=jds-1+jw
!
          if(jts<=jch)then
            do j=max(jts,jcl),min(jte,jch)
              do i=its,ite
                div(i,j,l)=0.
              enddo
            enddo
          endif
!
          jcl=jde-jw+1
          jch=jde
!
          if(jte>=jcl)then
            do j=max(jts,jcl),min(jte,jch)
              do i=its,ite
                div(i,j,l)=0.
              enddo
            enddo
          endif
!
          jcl=jds+jw
          if(jts<=jcl.and.jte>=jcl)then
            wprp=wpdar(jcl)
            do i=its_b1,ite_b1
              div(i,jcl,l)=(  pgx(i+1 ,jcl)-pgx(i,jcl)   &
                            + pgy(i   ,jcl+1)            &
                            -(pgne(i+1,jcl+1)            &
                            + pgnw(i  ,jcl+1))*0.5)*wprp
            enddo
          endif
!
          jch=jde-jw
          if(jts<=jch.and.jte>=jch)then
            wprp=wpdar(jch)
            do i=its_b1,ite_b1
              div(i,jch,l)=( pgx(i+1,jch)-pgx(i,jch)     &
                            -pgy(i  ,jch)                &
                          -(-pgne(i  ,jch)               &
                            -pgnw(i+1,jch))*0.5)*wprp
            enddo
          endif
!
          jcl=jds+jw+1
          jch=jde-jw-1
          if(jte>=jcl.and.jts<=jch)then
            do j=max(jts,jcl),min(jte,jch)
              wprp=wpdar(j)
              do i=its_b1,ite_b1
                div(i,j,l)=(pgx(i+1,j)-pgx(i,j)          &
                           +pgy(i,j+1)-pgy(i,j)          &
                          -(pgne(i+1,j+1)-pgne(i,j)      &
                           +pgnw(i,j+1)-pgnw(i+1,j))*0.5)*wprp
              enddo
            enddo
          endif
!-----------------------------------------------------------------------
        else ! regional
!-----------------------------------------------------------------------
          icl=ids+2
          ich=ide-2
!
          jcl=jds+2
          jch=jde-2
!
          if(s_bdy) then
            jp=jds+1
            wprp=wpdar(jp)
            do i=max(its,icl),min(ite,ich)
              div(i,jp,l)=(pgx(i+1,jp)-pgx(i,jp)         &
                          +pgy(i,jp+1)                   &
                         -(pgne(i+1,jp+1)                &
                          +pgnw(i,jp+1))*0.5)*wprp
            enddo
          endif
!
          if(n_bdy) then
            jp=jde-1
            wprp=wpdar(jp)
            do i=max(its,icl),min(ite,ich)
              div(i,jp,l)=(pgx(i+1,jp)-pgx(i,jp)         &
                          -pgy(i,jp)                     &
                        -(-pgne(i,jp)                    &
                          -pgnw(i+1,jp))*0.5)*wprp
            enddo
          endif
!
          if(w_bdy) then
            ip=ids+1
            do j=max(jts,jcl),min(jte,jch)
              wprp=wpdar(j)
              div(ip,j,l)=(pgx(ip+1,j)                   &
                          +pgy(ip,j+1)-pgy(ip,j)         &
                         -(pgne(ip+1,j+1)                &
                          -pgnw(ip+1,j))*0.5)*wprp
            enddo
          endif
!
          if(e_bdy) then
            ip=ide-1
            do j=max(jts,jcl),min(jte,jch)
              wprp=wpdar(j)
              div(ip,j,l)=(-pgx(ip,j)                    &
                           +pgy(ip,j+1)-pgy(ip,j)        &
                         -(-pgne(ip,j)                   &
                           +pgnw(ip,j+1))*0.5)*wprp
            enddo
          endif
!
          if(s_bdy.and.w_bdy) then
            ip=ids+1
            jp=jds+1
            wprp=wpdar(jp)
            div(ip,jp,l)=(pgx(ip+1,jp)                   &
                         +pgy(ip,jp+1)                   &
                        -(pgne(ip+1,jp+1))*0.5)*wprp
          endif
!
          if(s_bdy.and.e_bdy) then
            ip=ide-1
            jp=jds+1
            wprp=wpdar(jp)
            div(ip,jp,l)=(-pgx(ip,jp)                    &
                          +pgy(ip,jp+1)                  &
                         -(pgnw(ip,jp+1))*0.5)*wprp
          endif
!
          if(n_bdy.and.w_bdy) then
            ip=ids+1
            jp=jde-1
            wprp=wpdar(jp)
            div(ip,jp,l)=(pgx(ip+1,jp)                   &
                         -pgy(ip,jp)                     &
                       -(-pgnw(ip+1,jp))*0.5)*wprp
          endif
!
          if(n_bdy.and.e_bdy) then
            ip=ide-1
            jp=jde-1
            wprp=wpdar(jp)
            div(ip,jp,l)=(-pgx(ip,jp)                    &
                          -pgy(ip,jp)                    &
                        -(-pgne(ip,jp))*0.5)*wprp
          endif
!
          do j=jts_b2,jte_b2
            wprp=wpdar(j)
            do i=its_b2,ite_b2
              div(i,j,l)=(pgx(i+1,j)-pgx(i,j)            &
                         +pgy(i,j+1)-pgy(i,j)            &
                        -(pgne(i+1,j+1)-pgne(i,j)        &
                         +pgnw(i,j+1)-pgnw(i+1,j))*0.5)*wprp
            enddo
          enddo
!-----------------------------------------------------------------------
        endif ! global/regional branching
!-----------------------------------------------------------------------
!---first pass switch---------------------------------------------------
!-----------------------------------------------------------------------
        if(.not.first.or.restart) then
!-----------------------------------------------------------------------
!---updating u and v due to pgf force, end of time-step for u and v-----
!-----------------------------------------------------------------------
          rdv=rdyv*dt
          do j=jts_b1,jte_b2
            rdu=rdxv(j)*dt
            do i=its_b1,ite_b2
              apd=(pd(i,j)+pd(i+1,j)+pd(i,j+1)+pd(i+1,j+1))*0.25
              rpdp=0.3333333333/(dsg2(l)*apd+pdsg1(l))
!
              tcu(i,j,l)=-((pgx(i+1,j)+pgx(i+1,j+1)) &
                          +(pgne(i+1,j+1)-pgnw(i+1,j+1))*0.5)*rdu*rpdp &
                         +tcu(i,j,l)
              tcv(i,j,l)=-((pgy(i,j+1)+pgy(i+1,j+1)) &
                          +(pgne(i+1,j+1)+pgnw(i+1,j+1))*0.5)*rdv*rpdp &
                         +tcv(i,j,l)
            enddo
          enddo
!-----------------------------------------------------------------------
        else
!-----------------------------------------------------------------------
          do j=jts_b1,jte_b2
            do i=its_b1,ite_b2
              tcu(i,j,l)=0.
              tcv(i,j,l)=0.
            enddo
          enddo
!-----------------------------------------------------------------------
        endif
!-----------------------------------------------------------------------
!
      enddo vertical_loop
!.......................................................................
!$omp end parallel do
!.......................................................................
!
!-----------------------------------------------------------------------
!
                        endsubroutine pgforce
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
                        subroutine dht &
(global &
,lm &
,dyv &
,dsg2,pdsg1 &
,dxv,fcp,fdiv &
,pd,pdo &
,u,v &
,omgalf &
!---temporary arguments-------------------------------------------------
,pcne,pcnw,pcx,pcy,pfne,pfnw,pfx,pfy,div,tdiv)
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
implicit none
!-----------------------------------------------------------------------
real(kind=kfpt),parameter:: &
 cfc=1.533 &                 ! adams-bashforth positioning in time
,bfc=1.-cfc                  ! adams bashforth positioning in time
!-----------------------------------------------------------------------
logical(kind=klog),intent(in):: &
 global

integer(kind=kint),intent(in):: &
 lm                          ! total # of levels
 
real(kind=kfpt),intent(in):: &
 dyv                         ! deltay, v point
 
real(kind=kfpt),dimension(1:lm),intent(in):: &
 dsg2 &                      ! delta sigmas
,pdsg1                       ! delta pressures
 
real(kind=kfpt),dimension(jds:jde),intent(in):: &
 dxv &                       ! deltax, v point
,fcp &                       !
,fdiv                        !
 
real(kind=kfpt),dimension(ims:ime,jms:jme),intent(in):: &
 pd &                        ! sigma range pressure difference
,pdo                         ! sigma range pressure difference
 
real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(in):: &
 u &                         ! u wind component
,v                           ! v wind component
 
real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(out):: &
 omgalf                      ! omega-alfa (horizontal)
!-----------------------------------------------------------------------
!---temporary arguments-------------------------------------------------
!-----------------------------------------------------------------------
real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(in) :: &
 pcne &                      ! second term of pgf, ne direction
,pcnw &                      ! second term of pgf, nw direction
,pcx &                       ! second term of pgf, x direction
,pcy                         ! second term of pgf, y direction
 
real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(out) :: &
 pfne &                      ! mass flux, ne direction
,pfnw &                      ! mass flux, nw direction
,pfx &                       ! mass flux, x direction
,pfy &                       ! mass flux, y direction
,tdiv                        ! integrated horizontal mass divergence

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(inout) :: &
 div                         ! horizontal mass divergence
!-----------------------------------------------------------------------
!--local variables------------------------------------------------------
!-----------------------------------------------------------------------
integer(kind=kint):: &
 i &                         ! index in x direction
,j &                         ! index in y direction
,l                           ! index in p direction
 
real(kind=kfpt):: &
 dyp1 &                      !
,dyp &                       !
,dxp1 &                      !
,dxp &                       !
,fcpp &                      !
,fdp &                       !
,pdp &                       ! hydrostatic pressure difference at the point
,pdxp &                      ! hydrostatic pressure at the point
,pdyp &                      ! hydrostatic pressure at the point
,pdnep &                     ! hydrostatic pressure at the point
,pdnwp &                     ! hydrostatic pressure at the point
,udy &                       !
,vdx &                       !
,udy1 &                      !
,vdx1                        !
 
real(kind=kfpt),dimension(its_b1:ite_h2,jts_b1:jte_h2):: &
 pdne &                      ! hydrostatic pressure difference at the point
,pdnw &                      ! hydrostatic pressure difference at the point
,pdx &                       ! hydrostatic pressure difference at the point
,pdy &                       ! hydrostatic pressure difference at the point
,tne &                       ! temperature flux, ne direction
,tnw &                       ! temperature flux, nw direction
,tx &                        ! temperature flux, x direction
,ty                          ! temperature flux, y direction
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      dyp1=dyv
      dyp=dyp1*0.5
!
!-----------------
#ifdef ENABLE_SMP
!-----------------
!.......................................................................
!$omp parallel &
!$omp private(i,j,jstart,jstop,nth,tid)
!.......................................................................
      nth = omp_get_num_threads()
      tid = omp_get_thread_num()
      call looplimits(tid, nth,jts_b1,jte_h2,jstart,jstop)
!-----------------
#else
!-----------------
      jstart = jts_b1
      jstop  = jte_h2
!-----------------
#endif 
!-----------------
!
      do j=jstart,jstop
        do i=its_b1,ite_h2
          pdx (i,j)=((pd (i-1,j)+pd (i,j))*cfc &
                    +(pdo(i-1,j)+pdo(i,j))*bfc)*0.5
        enddo
      enddo
!
      do j=jstart,jstop
        do i=its_b1,ite_h2
          pdy (i,j)=((pd (i,j-1)+pd (i,j))*cfc &
                    +(pdo(i,j-1)+pdo(i,j))*bfc)*0.5
        enddo
      enddo
!
      do j=jstart,jstop
        do i=its_b1,ite_h2
          pdne(i,j)=((pd (i-1,j-1)+pd (i,j))*cfc &
                    +(pdo(i-1,j-1)+pdo(i,j))*bfc)*0.5
          pdnw(i,j)=((pd (i,j-1)+pd (i-1,j))*cfc &
                    +(pdo(i,j-1)+pdo(i-1,j))*bfc)*0.5
        enddo
      enddo
!-----------------
#ifdef ENABLE_SMP
!-----------------
!.......................................................................
!$omp end parallel
!.......................................................................
!-----------------
#endif
!-----------------
!-----------------------------------------------------------------------
!---vertical grand loop-------------------------------------------------
!-----------------------------------------------------------------------
!
!.......................................................................
!$omp parallel do                                                       &
!$omp private (dxp,dxp1,fcpp,fdp,i,j,l,pdnep,pdnwp,pdp,pdxp,pdyp,       &
!$omp          tne,tnw,tx,ty,udy,udy1,vdx,vdx1)
!.......................................................................
!-----------------------------------------------------------------------
      vertical_loop: do l=1,lm
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!---mass and pressure fluxes, on h points-------------------------------
!-----------------------------------------------------------------------
        pdp=pdsg1(l)
        do j=jts_b1,jte_h2
          dxp1=dxv(j-1)
          dxp=dxp1*0.5
          do i=its_b1,ite_h2
!-----------------------------------------------------------------------
            pdxp=dsg2(l)*pdx(i,j)+pdp
            pdyp=dsg2(l)*pdy(i,j)+pdp
            pdnep=dsg2(l)*pdne(i,j)+pdp
            pdnwp=dsg2(l)*pdnw(i,j)+pdp
!
            udy=(u(i-1,j-1,l)+u(i-1,j,l))*dyp
            vdx=(v(i-1,j-1,l)+v(i,j-1,l))*dxp
!
            pfx(i,j,l)=udy*pdxp
            pfy(i,j,l)=vdx*pdyp
!
            udy1=u(i-1,j-1,l)*dyp1
            vdx1=v(i-1,j-1,l)*dxp1
!
            pfne(i,j,l)=( udy1+vdx1)*pdnep
            pfnw(i,j,l)=(-udy1+vdx1)*pdnwp
!
            tx(i,j)=pcx(i,j,l)*udy
            ty(i,j)=pcy(i,j,l)*vdx
!
            tne(i,j)=pcne(i,j,l)*( udy1+vdx1)
            tnw(i,j)=pcnw(i,j,l)*(-udy1+vdx1)
          enddo
        enddo
!-----------------------------------------------------------------------
!---divergence and hor. pressure advection in t eq.---------------------
!-----------------------------------------------------------------------
        do j=jts_b1,jte_b1_h1
          fdp=fdiv(j)
          fcpp=fcp(j)
          do i=its_b1,ite_b1_h1
            div(i,j,l)=((pfx (i+1,j  ,l)-pfx (i  ,j  ,l) &
                        +pfy (i  ,j+1,l)-pfy (i  ,j  ,l)) &
                       +(pfne(i+1,j+1,l)-pfne(i  ,j  ,l) &
                        +pfnw(i  ,j+1,l)-pfnw(i+1,j  ,l))*0.25)*fdp  &
                      +div(i,j,l)
            tdiv(i,j,l)=div(i,j,l)
            omgalf(i,j,l)=((tx (i  ,j  )+tx (i+1,j  ) &
                           +ty (i  ,j  )+ty (i  ,j+1)) &
                          +(tne(i+1,j+1)+tne(i  ,j  ) &
                           +tnw(i  ,j+1)+tnw(i+1,j  ))*0.25) &
                         *fcpp/(dsg2(l)*pd(i,j)+pdsg1(l))
          enddo
        enddo
!-----------------------------------------------------------------------
!---zero divergence along regional domain boundaries--------------------
!-----------------------------------------------------------------------
        if(.not.global) then
          if(s_bdy)then
            do i=ims,ime
              div (i,jds,l)=0.
              tdiv(i,jds,l)=0.
            enddo
          endif
!
          if(n_bdy)then
            do i=ims,ime
              div (i,jde,l)=0.
              tdiv(i,jde,l)=0.
            enddo
          endif
!
          if(w_bdy)then
            do j=jms,jme
              div (ids,j,l)=0.
              tdiv(ids,j,l)=0.
            enddo
          endif
!
          if(e_bdy)then
            do j=jms,jme
              div (ide,j,l)=0.
              tdiv(ide,j,l)=0.
            enddo
          endif
       endif
!-----------------------------------------------------------------------
      enddo vertical_loop
!-----------------------------------------------------------------------
!.......................................................................
!$omp end parallel do
!.......................................................................
!-----------------
#ifdef ENABLE_SMP
!-----------------
!.......................................................................
!$omp parallel &
!$omp private(i,j,jstart,jstop,nth,tid)
!.......................................................................
      nth = omp_get_num_threads()
      tid = omp_get_thread_num()
      call looplimits(tid, nth,jts_b1,jte_b1,jstart,jstop)
!-----------------
#else
!-----------------
      jstart = jts_b1
      jstop  = jte_b1
!-----------------
#endif 
!-----------------
      do l=2,lm
        do j=jstart,jstop
          do i=its_b1,ite_b1
            tdiv(i,j,l)=tdiv(i,j,l-1)+tdiv(i,j,l)
          enddo
        enddo
      enddo
!-----------------
#ifdef ENABLE_SMP
!-----------------
!.......................................................................
!$omp end parallel
!.......................................................................
!-----------------
#endif
!-----------------
     
!-----------------------------------------------------------------------
!
                        endsubroutine dht
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
                        subroutine ddamp &
(lm &
,ddmpv,pdtop &
,dsg2,pdsg1 &
,sg1,sg2 &
,ddmpu &
,freerun &
,pd,pdo &
,u,v &
!---temporary arguments-------------------------------------------------
,div)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
implicit none
!-----------------------------------------------------------------------
real(kind=kfpt),parameter:: &
 cfc=1.533 &                 ! adams-bashforth positioning in time
,bfc=1.-cfc                  ! adams bashforth positioning in time
!-----------------------------------------------------------------------
logical(kind=klog),intent(in):: &
 freerun

integer(kind=kint),intent(in):: &
 lm                          ! total # of levels

real(kind=kfpt),intent(in):: &
 ddmpv &                     ! divergence damping, v component
,pdtop                       ! pressure coordinate depth

real(kind=kfpt),dimension(1:lm),intent(in):: &
 dsg2 &                      ! delta sigmas
,pdsg1                       ! delta pressures

real(kind=kfpt),dimension(1:lm+1),intent(in):: &
 sg1 &                       !
,sg2                         !

real(kind=kfpt),dimension(jds:jde),intent(in):: &
 ddmpu                       ! divergence damping, u direction

real(kind=kfpt),dimension(ims:ime,jms:jme),intent(in):: &
 pd &                        ! sigma range pressure difference
,pdo                         ! sigma range pressure difference

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(inout):: &
 u &                         ! u wind component
,v                           ! v wind component

!---temporary arguments-------------------------------------------------
real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm):: &
 div                         ! horizontal mass divergence

!--local variables------------------------------------------------------
logical(kind=klog),parameter:: &
 extmod=.true.

real(kind=kfpt),parameter:: &
 assimfc=5.0

integer(kind=kint):: &
 i &                         ! index in x direction
,j &                         ! index in y direction
,l                           ! index in p direction

real(kind=kfpt):: &
 dfac &                      ! fcim enhancement factor at top
,dpb &                       !
,dpl &                       !
,fcim &                      ! relative weight of internal mode damping
,fcxm &                      ! blow up factor for external mode damping
,fint &                      !
,rdpd &                      !
,dud &                       !
,dvd                         !

real(kind=kfpt),dimension(ims:ime,jms:jme):: &
 apd &                       !
,dive &                      !
,rddu &                      !
,rddv                        !
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
      dfac=9.
      dpb=1500.
      fint=(1.-dfac)/dpb
!-----------------------------------------------------------------------
!
      dvd=ddmpv
      do j=jts_b1,jte_b2
        do i=its_b1,ite_b2
          apd(i,j)=((pd (i  ,j  )+pd (i+1,j  ) &
                    +pd (i  ,j+1)+pd (i+1,j+1))*cfc &
                   +(pdo(i  ,j  )+pdo(i+1,j  ) &
                    +pdo(i  ,j+1)+pdo(i+1,j+1))*bfc)*0.25
        enddo
      enddo
!
!-----------------------------------------------------------------------
!---external mode-------------------------------------------------------
!-----------------------------------------------------------------------
!
      if(extmod) then
        if(freerun) then
          fcxm=1.0
        else
          fcxm=assimfc
        endif
!
!-----------------
#ifdef ENABLE_SMP
!-----------------
!.......................................................................
!$omp parallel private(i,j,l,jstart,jstop,nth,tid)
!.......................................................................
        nth = omp_get_num_threads()
        tid = omp_get_thread_num()
        call looplimits(tid, nth,jts_b1,jte_b1_h2,jstart,jstop)
!-----------------
#else
!-----------------
        jstart = jts_b1
        jstop  = jte_b1_h2
!-----------------
#endif
!-----------------
        do j=jstart,jstop
          do i=its_b1,ite_b1_h2
            dive(i,j)=div(i,j,1)
          enddo
        enddo
        do l=2,lm
        do j=jstart,jstop
          do i=its_b1,ite_b1_h2
              dive(i,j)=div(i,j,l)+dive(i,j)
            enddo
          enddo
        enddo
!-----------------
#ifdef ENABLE_SMP
!-----------------
!.......................................................................
!$omp end parallel
!.......................................................................
!-----------------
#endif
!-----------------
!
        do j=jts_b1,jte_b2
          dud=ddmpu(j)
          do i=its_b1,ite_b2
            rdpd=fcxm/(sg1(lm+1)*pdtop+sg2(lm+1)*apd(i,j))
            rddu(i,j)=(dive(i+1,j)+dive(i+1,j+1) &
                      -dive(i  ,j)-dive(i  ,j+1))*dud*rdpd
            rddv(i,j)=(dive(i,j+1)+dive(i+1,j+1) &
                      -dive(i,j  )-dive(i+1,j  ))*dvd*rdpd
          enddo
        enddo
!
        jstart = jts_b1
        jstop = jte_b2
!.......................................................................
!$omp parallel do private(dpl,dud,fcim,i,j,l,rdpd)
!.......................................................................
        do l=1,lm
          if(freerun) then
            dpl=sg1(l+1)*pdtop+sg2(l+1)*10000.
            if(dpl.lt.dpb) then
              fcim=fint*dpl+dfac
            else
              fcim=1.
            endif
          else
            fcim=assimfc
          endif
!
          do j=jstart,jstop
            dud=ddmpu(j)
            do i=its_b1,ite_b2
              rdpd=fcim/(dsg2(l)*apd(i,j)+pdsg1(l))
              u(i,j,l)=((div(i+1,j,l)+div(i+1,j+1,l) &
                        -div(i  ,j,l)-div(i  ,j+1,l))*dud)*rdpd &
                      +rddu(i,j)+u(i,j,l)
              v(i,j,l)=((div(i,j+1,l)+div(i+1,j+1,l) &
                        -div(i,j  ,l)-div(i+1,j  ,l))*dvd)*rdpd &
                      +rddv(i,j)+v(i,j,l)
            enddo
          enddo
        enddo      ! end of vertical loop
!.......................................................................
!$omp end parallel do
!.......................................................................

!-----------------------------------------------------------------------
!
      else
!
!-----------------------------------------------------------------------
!---divergence damping--------------------------------------------------
!-----------------------------------------------------------------------
!.......................................................................
!$omp parallel do private(dpl,dud,dvd,fcim,i,j,l,rdpd)
!.......................................................................
        do l=1,lm
          if(freerun) then
            dpl=sg1(l+1)*pdtop+sg2(l+1)*10000.
            if(dpl.lt.dpb) then
              fcim=fint*dpl+dfac
            else
              fcim=1.
            endif
          else
            fcim=assimfc
          endif
!
          dvd=ddmpv
          do j=jts_b1,jte_b2
            dud=ddmpu(j)
            do i=its_b1,ite_b2
              rdpd=fcim/(dsg2(l)*apd(i,j)+pdsg1(l))
              u(i,j,l)=(div(i+1,j,l)+div(i+1,j+1,l) &
                       -div(i  ,j,l)-div(i  ,j+1,l)) &
                      *dud*rdpd+u(i,j,l)
              v(i,j,l)=(div(i,j+1,l)+div(i+1,j+1,l) &
                       -div(i,j  ,l)-div(i+1,j  ,l)) &
                      *dvd*rdpd+v(i,j,l)
            enddo
          enddo
        enddo      ! end of vertical loop
!.......................................................................
!$omp end parallel do
!.......................................................................
!-----------------------------------------------------------------------
!
      endif
!
!-----------------------------------------------------------------------
!
                        endsubroutine ddamp
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
                        subroutine pdtsdt &
(lm &
,dt &
,sg2 &
,pd &
,pdo,psdt &
,psgdt &
!---temporary arguments-------------------------------------------------
,div,tdiv)
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
implicit none
!-----------------------------------------------------------------------
integer(kind=kint),intent(in):: &
 lm                          ! total # of levels

real(kind=kfpt),intent(in):: &
 dt                          ! time step

real(kind=kfpt),dimension(1:lm+1),intent(in):: &
 sg2                        !

real(kind=kfpt),dimension(ims:ime,jms:jme),intent(inout):: &
 pd                          ! sigma range pressure difference

real(kind=kfpt),dimension(ims:ime,jms:jme),intent(out):: &
 pdo &                       ! sigma range pressure difference
,psdt                        ! hydrostatic pressure tendency

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm-1),intent(out):: &
 psgdt                       ! vertical mass flux
!-----------------------------------------------------------------------
!---temporary arguments-------------------------------------------------
!-----------------------------------------------------------------------
real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm):: &
 div &                       ! horizontal mass divergence
,tdiv                        ! integrated unfiltered mass divergence
!-----------------------------------------------------------------------
!--local variables------------------------------------------------------
!-----------------------------------------------------------------------
integer(kind=kint):: &
 i &                         ! index in x direction
,j &                         ! index in y direction
,l &                         ! index in p direction
,loc_npts &                  ! local point counts for diag
,glb_npts &                  ! global point counts for diag
,iret &                    
,rc                         


real(kind=kfpt) :: &
 task_change &               ! task sum of abs(surface pressure change)
,global_change               ! domain total sum
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!-----------------
#ifdef ENABLE_SMP
!-----------------
!.......................................................................
!$omp parallel &
!$omp private(i,j,jstart,jstop,l,nth,tid)
!.......................................................................
      nth = omp_get_num_threads()
      tid = omp_get_thread_num()
      call looplimits(tid, nth,jts_b1,jte_b1,jstart,jstop)
!-----------------
#else
!-----------------
      jstart = jts_b1
      jstop  = jte_b1
!-----------------
#endif
!-----------------
      do l=2,lm
        do j=jstart,jstop
          do i=its_b1,ite_b1
            div(i,j,l)=div(i,j,l-1)+div(i,j,l)
          enddo
        enddo
      enddo
!-----------------
#ifdef ENABLE_SMP
!-----------------
!.......................................................................
!$omp end parallel
!.......................................................................
!-----------------
#endif
!-----------------

!-----------------------------------------------------------------------
      do j=jts_h2,jte_h2
        do i=its_h2,ite_h2
          psdt(i,j)=0.
          pdo(i,j)=pd(i,j)
        enddo
      enddo
!
!-----------------
#ifdef ENABLE_SMP
!-----------------
!.......................................................................
!$omp parallel &
!$omp private(i,j,jstart,jstop,l,nth,tid)
!.......................................................................
      nth = omp_get_num_threads()
      tid = omp_get_thread_num()
      call looplimits(tid, nth,jts_b1,jte_b1,jstart,jstop)
!-----------------
#else
!-----------------
      jstart = jts_b1
      jstop  = jte_b1
!-----------------
#endif
!-----------------

!      loc_npts=0
!      task_change=0.

      do j=jstart,jstop
        do i=its_b1,ite_b1
          psdt(i,j)=-div(i,j,lm)
          pd(i,j)=psdt(i,j)*dt+pd(i,j)
        enddo
      enddo

!      do j=jstart,jstop
!        do i=its_b1,ite_b1
!          if (abs(psdt(i,j)*dt) .ge. 250. .or. psdt(i,j) .ne. psdt(i,j)) then
!            write(0,*) 'big PD change...I,J, change: ', I,J, psdt(i,j)*dt
!	    write(0,*) 'previous PD was: ', PD(I,J)-psdt(i,j)*dt
!          endif
!
!          if (abs(psdt(i,j)*dt) .ge. 5000.) then
!            write(0,*) 'huge PD change...I,J, change: ', I,J, psdt(i,j)*dt
!            call ESMF_Finalize(endflag=ESMF_END_ABORT)
!          endif
!
!          loc_npts=loc_npts+1
!          task_change=task_change+abs(PSDT(I,J)) ! * (108./abs(dt)) )  ! .01*10800/dt (hPa/3 h)
!
!        enddo
!      enddo
!
!      call mpi_reduce(task_change, global_change, 1, mpi_real, mpi_sum,0, &
!                      mpi_comm_comp, iret)
!
!      call mpi_reduce(loc_npts, glb_npts, 1, mpi_integer, mpi_sum,0, &
!                      mpi_comm_comp, iret)
!
!      if (mype == 0) then
!        if (dt .gt. 0) then
!          write(0,*) 'FWD avg global change (Pa/timestep): ', GLOBAL_CHANGE/ GLB_NPTS
!        else
!          write(0,*) 'BCKWD avg global change (Pa/timestep): ', GLOBAL_CHANGE/ GLB_NPTS
!        endif
!      endif

!-----------------------------------------------------------------------
!---boundary conditions-------------------------------------------------
!-----------------------------------------------------------------------
      do l=1,lm-1
        do j=jstart,jstop
          do i=its_b1,ite_b1
            psgdt(i,j,l)=-(-tdiv(i,j,lm)*sg2(l+1)+tdiv(i,j,l))
          enddo
        enddo
      enddo
!-----------------
#ifdef ENABLE_SMP
!-----------------
!.......................................................................
!$omp end parallel
!.......................................................................
!-----------------
#endif
!-----------------
!-----------------------------------------------------------------------
!
                        endsubroutine pdtsdt
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
                        subroutine adv1 &
(global,secadv &
,lm,lnsad,inpes,jnpes &
,dt,dyv,rdyh,rdyv &
,dsg2,pdsg1 &
,curv,dxv,fad,fah,rdxh,rdxv &
,f,pd,pdo &
,omgalf,psgdt &
,t,u,v &
,tp,up,vp &
!---temporary arguments-------------------------------------------------
,pfne,pfnw,pfx,pfy,tct,tcu,tcv)
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
implicit none
!-----------------------------------------------------------------------
real(kind=kfpt),parameter:: &
 cfc=1.533 &                 ! adams-bashforth positioning in time
,bfc=1.-cfc &                ! adams bashforth positioning in time
,epscm=2.e-6 &               ! a floor value (not used)
,pfc=1.+4./6. &              ! 4th order momentum advection
,sfc=-1./6. &                ! 4th order momentum advection
,w1=0.9 &                    ! crank-nicholson uncentering
!,w1=0.0 &                    ! crank-nicholson uncentering
,w2=2.-w1                    ! crank-nicholson uncentering

logical(kind=klog),intent(in):: &
 global &                    ! global or regional
,secadv                      ! second order momentum advection

integer(kind=kint),intent(in):: &
 lm &                        ! total # of levels
,lnsad &                     ! # of boundary lines w. upstream advection
,inpes &                     ! domain decomposition parameter
,jnpes                       ! domain decomposition parameter

real(kind=kfpt),intent(in):: &
 dt &                        ! dynamics time step
,dyv &                       ! deltay
,rdyh &                      ! 1/deltay
,rdyv                        ! 1/deltay

real(kind=kfpt),dimension(1:lm),intent(in):: &
 dsg2 &                      ! delta sigmas
,pdsg1                       ! delta pressures

real(kind=kfpt),dimension(jds:jde),intent(in):: &
 curv &                      ! curvature
,dxv &                       ! dxv
,fad &                       ! grid factor
,fah &                       ! grid factor
,rdxh &                      ! 1/deltax
,rdxv                        ! 1/deltax

real(kind=kfpt),dimension(ims:ime,jms:jme),intent(in):: &
 f &                         ! coriolis parameter
,pd &                        ! sigma range pressure difference
,pdo                         ! sigma range pressure difference

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm-1),intent(in):: &
 psgdt                       ! vertical mass flux

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(in):: &
 omgalf &                    !
,t &                         ! temperature
,u &                         ! u wind component
,v                           ! v wind component

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(inout):: &
 tp &                        ! old temperature
,up &                        ! old u
,vp                          ! old v
!-----------------------------------------------------------------------
!---temporary arguments-------------------------------------------------
!-----------------------------------------------------------------------
real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(in) :: &
 pfne &                      ! mass flux, ne direction
,pfnw &                      ! mass flux, nw direction
,pfx &                       ! mass flux, x direction
,pfy                         ! mass flux, y direction

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(out) :: &
 tct &                       ! time change of temperature
,tcu &                       ! time change of u
,tcv                         ! time change of v
!-----------------------------------------------------------------------
!--local variables------------------------------------------------------
!-----------------------------------------------------------------------
integer(kind=kint):: &
 i &                         ! index in x direction
,iap &                       ! offset in x direction
,ibeg &                      ! starting i in some horiz advec loops
,iend &                      ! ending i in some horiz advec loops
,j &                         ! index in y direction
,jap &                       ! offset in y direction
,jbeg &                      ! starting j in some horiz advec loops
,jend &                      ! ending j in some horiz advec loops
,l                           ! index in p direction

real(kind=kfpt):: &
 cf &                        ! temporary
,cmt &                       ! temporary
,cmw &                       ! temporary
,crv &                       ! curvature
,dtq &                       ! dt/4
,dts &                       ! dt/16
,dux1 &                      ! momentum advection component
,duy1 &                      ! momentum advection component
,dvx1 &                      ! momentum advection component
,dvy1 &                      ! momentum advection component
,dxody &                     !
,dyodx &                     !
,emhp &                      ! scaling in x direction
,enhp &                      ! scaling in y direction
,emvp &                      ! scaling in x direction
,envp &                      ! scaling in y direction
,fadp &                      ! grid factor at the point
,fahp &                      ! temporary grid factor
,fdpp &                      ! temporary grid factor
,fp &                        ! coriolis parameter factor
,fpp &                       ! coriolis with curvature
,pp &                        ! scaled trajectory, x direction
,qq &                        ! scaled trajectory, y direction
,rdp &                       ! 1/deltap
,rdxp &                      !
,rdyp &                      !
,vvlo &                      ! vertical velocity, lower interface
,vvup &                      ! vertical velocity, upper interface
,pvvup                       ! vertical mass flux, upper interface

real(kind=kfpt),dimension(its:ite,jts:jte):: &
 pdop &                      ! hydrostatic pressure difference at v points
,pvvlo                       ! vertical mass flux, lower interface

real(kind=kfpt),dimension(its_b1:ite_b1_h1,jts_b1:jte_b1_h1):: &
 pfnex1 &                    ! average mass flux for momentum advection
,pfney1 &                    ! average mass flux for momentum advection
,pfnwx1 &                    ! average mass flux for momentum advection
,pfnwy1 &                    ! average mass flux for momentum advection
,pfxx1 &                     ! average mass flux for momentum advection
,pfxy1 &                     ! average mass flux for momentum advection
,pfyx1 &                     ! average mass flux for momentum advection
,pfyy1 &                     ! average mass flux for momentum advection
,ufnex1 &                    ! average mass flux for momentum advection
,ufney1 &                    ! average mass flux for momentum advection
,ufnwx1 &                    ! average mass flux for momentum advection
,ufnwy1 &                    ! average mass flux for momentum advection
,ufxx1 &                     ! average mass flux for momentum advection
,ufxy1 &                     ! average mass flux for momentum advection
,ufyx1 &                     ! average mass flux for momentum advection
,ufyy1 &                     ! average mass flux for momentum advection
,vfnex1 &                    ! average mass flux for momentum advection
,vfney1 &                    ! average mass flux for momentum advection
,vfnwx1 &                    ! average mass flux for momentum advection
,vfnwy1 &                    ! average mass flux for momentum advection
,vfxx1 &                     ! average mass flux for momentum advection
,vfxy1 &                     ! average mass flux for momentum advection
,vfyx1 &                     ! average mass flux for momentum advection
,vfyy1                       ! average mass flux for momentum advection

real(kind=kfpt),dimension(its_b1:ite_h1,jts_b1:jte_h1):: &
 tne &                       ! temperature flux, ne direction
,tnw &                       ! temperature flux, nw direction
,tx &                        ! temperature flux, x direction
,ty                          ! temperature flux, y direction

real(kind=kfpt),dimension(its_h1:ite_h1,jts_h1:jte_h1):: &
 t1                          ! extrapolated temperature between time levels

real(kind=kfpt),dimension(ims:ime,jms:jme):: &
 u2d &                       ! 4th order diagonal u between time levels
,v2d &                       ! 4th order diagonal v between time levels

!real(kind=kfpt),dimension(its_h1:ite_h1,jts_h1:jte_h1):: &
,u1d &                       ! extrapolated diagonal u between time levels
,v1d                         ! extrapolated diagonal v between time levels

real(kind=kfpt),dimension(its_b1:ite_b1,jts_b1:jte_b1,1:lm):: &
 crt &                       ! vertical advection temporary
,rcmt &                      ! vertical advection temporary
,rstt                        ! vertical advection temporary

real(kind=kfpt),dimension(its_b1:ite_b2,jts_b1:jte_b2,1:lm):: &
 crw &                       ! vertical advection temporary
,rcmw &                      ! vertical advection temporary
,rstu &                      ! vertical advection temporary
,rstv                        ! vertical advection temporary
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
      do j=jts,jte
        do i=its,ite
          pdop(i,j)=(pd(i,j)+pdo(i,j))*0.5
        enddo
      enddo
!-----------------------------------------------------------------------
!---crank-nicholson vertical advection----------------------------------
!-----------------------------------------------------------------------
!
      dtq=dt*0.25
!
!-----------------
#ifdef ENABLE_SMP
!-----------------
!.......................................................................
!$omp parallel private(cf,cmt,i,j,jstart,jstop,l,nth,pvvup,rdp,tid,     &
!$omp                  vvlo,vvup)
!.......................................................................
      nth = omp_get_num_threads()
      tid = omp_get_thread_num()
      call looplimits(tid, nth, jts_b1, jte_b1, jstart, jstop)
!-----------------
#else
!-----------------
      jstart = jts_b1
      jstop  = jte_b1
!-----------------
#endif
!-----------------
       do j=jstart, jstop
        do i=its_b1,ite_b1
          pvvlo(i,j)=psgdt(i,j,1)*dtq
          vvlo=pvvlo(i,j)/(dsg2(1)*pdop(i,j)+pdsg1(1))
          cmt=-vvlo*w2+1.
!          if(abs(cmt).lt.epscm) cmt=epscm
          rcmt(i,j,1)=1./cmt
          crt(i,j,1)=vvlo*w2
          rstt(i,j,1)=(-vvlo*w1)*(t(i,j,2)-t(i,j,1))+t(i,j,1)
        enddo
      enddo
!
      do l=2,lm-1
       do j=jstart, jstop
          do i=its_b1,ite_b1
            rdp=1./(dsg2(l)*pdop(i,j)+pdsg1(l))
            pvvup=pvvlo(i,j)
            pvvlo(i,j)=psgdt(i,j,l)*dtq
!
            vvup=pvvup*rdp
            vvlo=pvvlo(i,j)*rdp
!
            cf=-vvup*w2*rcmt(i,j,l-1)
            cmt=-crt(i,j,l-1)*cf+((vvup-vvlo)*w2+1.)
!            if(abs(cmt).lt.epscm) cmt=epscm
            rcmt(i,j,l)=1./cmt
            crt(i,j,l)=vvlo*w2
            rstt(i,j,l)=-rstt(i,j,l-1)*cf+t(i,j,l) &
                        -((t(i,j,l)-t(i,j,l-1))*vvup &
                         +(t(i,j,l+1)-t(i,j,l))*vvlo)*w1
          enddo
        enddo
      enddo
!
       do j=jstart, jstop
        do i=its_b1,ite_b1
          pvvup=pvvlo(i,j)
          vvup=pvvup/(dsg2(lm)*pdop(i,j)+pdsg1(lm))
!
          cf=-vvup*w2*rcmt(i,j,lm-1)
          cmt=-crt(i,j,lm-1)*cf+(vvup*w2+1.)
!          if(abs(cmt).lt.epscm) cmt=epscm
          rcmt(i,j,lm)=1./cmt
          crt(i,j,lm)=0.
          rstt(i,j,lm)=-rstt(i,j,lm-1)*cf+t(i,j,lm) &
     &                 -(t(i,j,lm)-t(i,j,lm-1))*vvup*w1
!
          tct(i,j,lm)=rstt(i,j,lm)*rcmt(i,j,lm)-t(i,j,lm)
       enddo
      enddo
!
      do l=lm-1,1,-1
       do j=jstart, jstop
          do i=its_b1,ite_b1
            tct(i,j,l)=(-crt(i,j,l)*(t(i,j,l+1)+tct(i,j,l+1)) &
                        +rstt(i,j,l)) &
                      *rcmt(i,j,l)-t(i,j,l)
          enddo
        enddo
      enddo
!-----------------
#ifdef ENABLE_SMP
!-----------------
!.......................................................................
!$omp end parallel
!.......................................................................
!-----------------
#endif
!-----------------

!.......................................................................
!$omp parallel do &
!$omp private (emhp,enhp,fahp,i,iap,ibeg,iend,j,jap,jbeg,jend,l,        &
!$omp          pp,qq,t1,tne,tnw,tx,ty)
!.......................................................................
!-----------------------------------------------------------------------
      vertical_loop1: do l=1,lm
!-----------------------------------------------------------------------
        do j=jts_h1,jte_h1
          do i=its_h1,ite_h1
            t1(i,j)=t(i,j,l)*cfc+tp(i,j,l)*bfc
            tp(i,j,l)=t(i,j,l)
          enddo
        enddo
!-----------------------------------------------------------------------
!---temperature fluxes, on h points-------------------------------------
!-----------------------------------------------------------------------
        do j=jts_b1,jte_h1
          do i=its_b1,ite_h1
            tx(i,j)=(t1(i,j)-t1(i-1,j))*pfx(i,j,l)
            ty(i,j)=(t1(i,j)-t1(i,j-1))*pfy(i,j,l)
            tne(i,j)=(t1(i,j)-t1(i-1,j-1))*pfne(i,j,l)
            tnw(i,j)=(t1(i-1,j)-t1(i,j-1))*pfnw(i,j,l)
          enddo
        enddo
!-----------------------------------------------------------------------
!---advection of temperature--------------------------------------------
!-----------------------------------------------------------------------
        if(adv_standard)then
          ibeg=max(its,ids+1+lnsad)
          iend=min(ite,ide-1-lnsad)
          jbeg=max(jts,jds+1+lnsad)
          jend=min(jte,jde-1-lnsad)
!
          do j=jbeg,jend
            fahp=fah(j)
            do i=ibeg,iend
              tct(i,j,l)=(((tx (i  ,j  )+tx (i+1,j  ) &
                           +ty (i  ,j  )+ty (i  ,j+1)) &
                          +(tne(i+1,j+1)+tne(i  ,j  ) &
                           +tnw(i  ,j+1)+tnw(i+1,j  ))*0.25)*fahp) &
                        /(dsg2(l)*pdop(i,j)+pdsg1(l)) &
                        +omgalf(i,j,l)+tct(i,j,l)
            enddo
          enddo
        endif
!-----------------------------------------------------------------------
!---regional branch-----------------------------------------------------
!-----------------------------------------------------------------------
        if(.not.global.and.adv_upstream) then
!-----------------------------------------------------------------------
          enhp=-dt*rdyh*0.25
!
!***  Upstream advection along southern rows
!
          do j=jts_b1,min(jte,jds+lnsad)
            emhp=-dt*rdxh(j)*0.25
            do i=its_b1,ite_b1
              pp=(u(i-1,j-1,l)+u(i,j-1,l)+u(i-1,j,l)+u(i,j,l))*emhp
              qq=(v(i-1,j-1,l)+v(i,j-1,l)+v(i-1,j,l)+v(i,j,l))*enhp
!
              if(pp.le.0.) then
                iap=-1
                pp=-pp
              else
                iap=1
              endif
!
              if(qq.le.0.) then
                jap=-1
                qq=-qq
              else
                jap=1
              endif
!
              tct(i,j,l)=(t(i+iap,j,l)-t(i,j,l))*pp &
                        +(t(i,j+jap,l)-t(i,j,l))*qq &
                        +(t(i,j,l)-t(i+iap,j,l) &
                         -t(i,j+jap,l)+t(i+iap,j+jap,l))*pp*qq &
                        +omgalf(i,j,l)+tct(i,j,l)
            enddo
          enddo
!
!***  Upstream advection along northern rows
!
          do j=max(jts,jde-lnsad),jte_b1
            emhp=-dt*rdxh(j)*0.25
            do i=its_b1,ite_b1
              pp=(u(i-1,j-1,l)+u(i,j-1,l)+u(i-1,j,l)+u(i,j,l))*emhp
              qq=(v(i-1,j-1,l)+v(i,j-1,l)+v(i-1,j,l)+v(i,j,l))*enhp
!
              if(pp.le.0.) then
                iap=-1
                pp=-pp
              else
                iap=1
              endif
!
              if(qq.le.0.) then
                jap=-1
                qq=-qq
              else
                jap=1
              endif
!
              tct(i,j,l)=(t(i+iap,j,l)-t(i,j,l))*pp &
                        +(t(i,j+jap,l)-t(i,j,l))*qq &
                        +(t(i,j,l)-t(i+iap,j,l) &
                         -t(i,j+jap,l)+t(i+iap,j+jap,l))*pp*qq &
                        +omgalf(i,j,l)+tct(i,j,l)
            enddo
          enddo
!
!***  Upstream advection along western rows
!
          do j=max(jts,jds+1+lnsad),min(jte,jde-1-lnsad)
            emhp=-dt*rdxh(j)*0.25
            do i=its_b1,min(ite,ids+lnsad)
              pp=(u(i-1,j-1,l)+u(i,j-1,l)+u(i-1,j,l)+u(i,j,l))*emhp
              qq=(v(i-1,j-1,l)+v(i,j-1,l)+v(i-1,j,l)+v(i,j,l))*enhp
!
              if(pp.le.0.) then
                iap=-1
                pp=-pp
              else
                iap=1
              endif
!
              if(qq.le.0.) then
                jap=-1
                qq=-qq
              else
                jap=1
              endif
!
              tct(i,j,l)=(t(i+iap,j,l)-t(i,j,l))*pp &
                        +(t(i,j+jap,l)-t(i,j,l))*qq &
                        +(t(i,j,l)-t(i+iap,j,l) &
                         -t(i,j+jap,l)+t(i+iap,j+jap,l))*pp*qq &
                        +omgalf(i,j,l)+tct(i,j,l)
            enddo
          enddo
!
!***  Upstream advection along eastern rows
!
          do j=max(jts,jds+1+lnsad),min(jte,jde-1-lnsad)
            emhp=-dt*rdxh(j)*0.25
            do i=max(its,ide-lnsad),ite_b1
              pp=(u(i-1,j-1,l)+u(i,j-1,l)+u(i-1,j,l)+u(i,j,l))*emhp
              qq=(v(i-1,j-1,l)+v(i,j-1,l)+v(i-1,j,l)+v(i,j,l))*enhp
!
              if(pp.le.0.) then
                iap=-1
                pp=-pp
              else
                iap=1
              endif
!
              if(qq.le.0.) then
                jap=-1
                qq=-qq
              else
                jap=1
              endif
!
              tct(i,j,l)=(t(i+iap,j,l)-t(i,j,l))*pp &
                        +(t(i,j+jap,l)-t(i,j,l))*qq &
                        +(t(i,j,l)-t(i+iap,j,l) &
                         -t(i,j+jap,l)+t(i+iap,j+jap,l))*pp*qq &
                        +omgalf(i,j,l)+tct(i,j,l)
            enddo
          enddo
!-----------------------------------------------------------------------
        endif ! regional lateral boundaries
!-----------------------------------------------------------------------
!
      enddo vertical_loop1
!.......................................................................
!$omp end parallel do
!.......................................................................
!
!-----------------------------------------------------------------------
!---crank-nicholson vertical advection----------------------------------
!-----------------------------------------------------------------------
!
      dts=dt*(0.25*0.25)
!
!-----------------
#ifdef ENABLE_SMP
!-----------------
!.......................................................................
!$omp parallel private(cf,cmw,i,j,jstart,jstop,l,nth,pvvup,rdp,tid,     &
!$omp                  vvlo,vvup)
!.......................................................................
       nth = omp_get_num_threads()
       tid = omp_get_thread_num()
       call looplimits(tid, nth, jts_b1, jte_b2, jstart, jstop)
!-----------------
#else
!-----------------
      jstart = jts_b1
      jstop  = jte_b2
!-----------------
#endif
!-----------------
      do j=jstart, jstop
        do i=its_b1,ite_b2
          pdop(i,j)=(pd (i,j)+pd (i+1,j)+pd (i,j+1)+pd (i+1,j+1) &
                    +pdo(i,j)+pdo(i+1,j)+pdo(i,j+1)+pdo(i+1,j+1))*0.125
          pvvlo(i,j)=(psgdt(i,j  ,1)+psgdt(i+1,j  ,1) &
                     +psgdt(i,j+1,1)+psgdt(i+1,j+1,1))*dts
          vvlo=pvvlo(i,j)/(dsg2(1)*pdop(i,j)+pdsg1(1))
          cmw=-vvlo*w2+1.
!          if(abs(cmw).lt.epscm) cmw=epscm
          rcmw(i,j,1)=1./cmw
          crw(i,j,1)=vvlo*w2
          rstu(i,j,1)=(-vvlo*w1)*(u(i,j,2)-u(i,j,1))+u(i,j,1)
          rstv(i,j,1)=(-vvlo*w1)*(v(i,j,2)-v(i,j,1))+v(i,j,1)
        enddo
      enddo
!
      do l=2,lm-1
        do j=jstart,jstop
          do i=its_b1,ite_b2
            rdp=1./(dsg2(l)*pdop(i,j)+pdsg1(l))
            pvvup=pvvlo(i,j)
            pvvlo(i,j)=(psgdt(i,j,l)+psgdt(i+1,j,l) &
                       +psgdt(i,j+1,l)+psgdt(i+1,j+1,l))*dts
            vvup=pvvup*rdp
            vvlo=pvvlo(i,j)*rdp
            cf=-vvup*w2*rcmw(i,j,l-1)
            cmw=-crw(i,j,l-1)*cf+((vvup-vvlo)*w2+1.)
!            if(abs(cmw).lt.epscm) cmw=epscm
            rcmw(i,j,l)=1./cmw
            crw(i,j,l)=vvlo*w2
            rstu(i,j,l)=-rstu(i,j,l-1)*cf+u(i,j,l) &
                        -((u(i,j,l)-u(i,j,l-1))*vvup &
                         +(u(i,j,l+1)-u(i,j,l))*vvlo)*w1
            rstv(i,j,l)=-rstv(i,j,l-1)*cf+v(i,j,l) &
                        -((v(i,j,l)-v(i,j,l-1))*vvup &
                         +(v(i,j,l+1)-v(i,j,l))*vvlo)*w1
          enddo
        enddo
      enddo
!
      do j=jstart,jstop
        do i=its_b1,ite_b2
          pvvup=pvvlo(i,j)
          vvup=pvvup/(dsg2(lm)*pdop(i,j)+pdsg1(lm))
          cf=-vvup*w2*rcmw(i,j,lm-1)
          cmw=-crw(i,j,lm-1)*cf+(vvup*w2+1.)
!          if(abs(cmw).lt.epscm) cmw=epscm
          rcmw(i,j,lm)=1./cmw
          crw(i,j,lm)=vvlo
          rstu(i,j,lm)=-rstu(i,j,lm-1)*cf+u(i,j,lm) &
                       -(u(i,j,lm)-u(i,j,lm-1))*vvup*w1
          rstv(i,j,lm)=-rstv(i,j,lm-1)*cf+v(i,j,lm) &
                       -(v(i,j,lm)-v(i,j,lm-1))*vvup*w1
          tcu(i,j,lm)=rstu(i,j,lm)*rcmw(i,j,lm)-u(i,j,lm)
          tcv(i,j,lm)=rstv(i,j,lm)*rcmw(i,j,lm)-v(i,j,lm)
        enddo
      enddo
!
      do l=lm-1,1,-1
        do j=jstart,jstop
          do i=its_b1,ite_b2
            tcu(i,j,l)=(-crw(i,j,l)*(u(i,j,l+1)+tcu(i,j,l+1)) &
                        +rstu(i,j,l)) &
                      *rcmw(i,j,l)-u(i,j,l)
            tcv(i,j,l)=(-crw(i,j,l)*(v(i,j,l+1)+tcv(i,j,l+1)) &
                        +rstv(i,j,l)) &
                      *rcmw(i,j,l)-v(i,j,l)
          enddo
        enddo
      enddo
!-----------------
#ifdef ENABLE_SMP
!-----------------
!.......................................................................
!$omp end parallel
!.......................................................................
!-----------------
#endif
!-----------------

!-----------------------------------------------------------------------
!---grand vertical loop-------------------------------------------------
!-----------------------------------------------------------------------
!
!.......................................................................
!$omp parallel do &
!$omp private (crv,dux1,duy1,dvx1,dvy1,dxody,dyodx,emvp,envp,           &
!$omp          fadp,fdpp,fp,fpp,i,iap,ibeg,iend,j,jap,jbeg,jend,l,      &
!$omp          pfnex1,pfney1,pfnwx1,pfnwy1,pfxx1,pfxy1,pfyx1,pfyy1,pp,  &
!$omp          qq,u1d,u2d,ufnex1,ufney1,ufnwx1,ufxx1,ufyx1,ufnwy1,      &
!$omp          ufxy1,ufyy1,v1d,v2d,vfnex1,vfney1,vfnwx1,vfnwy1,         &
!$omp          vfxx1,vfxy1,vfyx1,vfyy1)
!.......................................................................
!-----------------------------------------------------------------------
!
      vertical_loop2: do l=1,lm
!
!-----------------------------------------------------------------------
!
        do j=jts_h2,jte_h2
          do i=its_h2,ite_h2
            u1d(i,j)=u(i,j,l)*cfc+up(i,j,l)*bfc
            v1d(i,j)=v(i,j,l)*cfc+vp(i,j,l)*bfc
!
            up(i,j,l)=u(i,j,l)
            vp(i,j,l)=v(i,j,l)
          enddo
        enddo
!-----------------------------------------------------------------------
!---coriolis force tendency---------------------------------------------
!-----------------------------------------------------------------------
        do j=jts_b1,jte_b2
          crv=curv(j)*dt
          do i=its_b1,ite_b2
            fp=f(i,j)*dt
            fpp=u1d(i,j)*crv+fp
!
            tcu(i,j,l)= fpp*v1d(i,j)+tcu(i,j,l)
            tcv(i,j,l)=-fpp*u1d(i,j)+tcv(i,j,l)
          enddo
        enddo
!-----------------------------------------------------------------------
!---diagonally averaged fluxes on v points------------------------------
!-----------------------------------------------------------------------
        do j=jts_b1,jte_b1_h1
          do i=its_b1,ite_b1_h1
            pfxx1 (i,j)=pfx (i  ,j  ,l)+pfx (i+1,j+1,l)
            pfyx1 (i,j)=pfy (i  ,j  ,l)+pfy (i+1,j+1,l)
            pfnex1(i,j)=pfne(i  ,j  ,l)+pfne(i+1,j+1,l)
            pfnwx1(i,j)=pfnw(i  ,j  ,l)+pfnw(i+1,j+1,l)
!
            pfxy1 (i,j)=pfx (i+1,j  ,l)+pfx (i  ,j+1,l)
            pfyy1 (i,j)=pfy (i+1,j  ,l)+pfy (i  ,j+1,l)
            pfney1(i,j)=pfne(i+1,j  ,l)+pfne(i  ,j+1,l)
            pfnwy1(i,j)=pfnw(i+1,j  ,l)+pfnw(i  ,j+1,l)
          enddo
        enddo
!-----------------------------------------------------------------------
!---4th order momentum advection----------------------------------------
!-----------------------------------------------------------------------
        if(.not.secadv) then
          do j=jts_b1_h1,jte_b1_h1
            do i=its_b1_h1,ite_b1_h1
              u2d(i,j)=(u1d(i,j-1)+u1d(i-1,j) &
                       +u1d(i+1,j)+u1d(i,j+1))*sfc &
                      +(u1d(i,j)*pfc)
              v2d(i,j)=(v1d(i,j-1)+v1d(i-1,j) &
                       +v1d(i+1,j)+v1d(i,j+1))*sfc &
                      +(v1d(i,j)*pfc)
            enddo
          enddo
          if(global) then
            btim=timef()
            call swapwn(u2d,ims,ime,jms,jme,1,inpes)
            call swapwn(v2d,ims,ime,jms,jme,1,inpes)
            timers(my_domain_id)%swapwn_tim=timers(my_domain_id)%swapwn_tim+(timef()-btim)
!
            btim=timef()
            call polewn(u2d,v2d,ims,ime,jms,jme,1,inpes,jnpes)
            timers(my_domain_id)%polewn_tim=timers(my_domain_id)%polewn_tim+(timef()-btim)
          else
            if(s_bdy)then
              do i=ims,ime
                u2d(i,jds)=u1d(i,jds)
                v2d(i,jds)=v1d(i,jds)
              enddo
            endif
            if(n_bdy)then
              do i=ims,ime
                u2d(i,jde-1)=u1d(i,jde-1)
                v2d(i,jde-1)=v1d(i,jde-1)
                u2d(i,jde)=u1d(i,jde-1)
                v2d(i,jde)=v1d(i,jde-1)
              enddo
            endif
            if(w_bdy)then
              do j=jms,jme
                u2d(ids,j)=u1d(ids,j)
                v2d(ids,j)=v1d(ids,j)
              enddo
            endif
            if(e_bdy)then
              do j=jms,jme
                u2d(ide-1,j)=u1d(ide-1,j)
                v2d(ide-1,j)=v1d(ide-1,j)
                u2d(ide,j)=u1d(ide-1,j)
                v2d(ide,j)=v1d(ide-1,j)
              enddo
            endif
          endif
!
          do j=jts_h1,jte_h1
            do i=its_h1,ite_h1
              u1d(i,j)=u2d(i,j)
              v1d(i,j)=v2d(i,j)
            enddo
          enddo
        endif
!-----------------------------------------------------------------------
!---horizontal fluxes of momentum components on v points----------------
!-----------------------------------------------------------------------
        if(global) then
          if(s_bdy) then
            do i=its_b1,ite_b1_h1
              pfyx1 (i,jts_b1)=0.
              pfyy1 (i,jts_b1)=0.
              pfnex1(i,jts_b1)=0.
              pfnwx1(i,jts_b1)=0.
              pfney1(i,jts_b1)=0.
              pfnwy1(i,jts_b1)=0.
            enddo
          endif
!
          if(n_bdy) then
            do i=its_b1,ite_b1_h1
              pfyx1 (i,jte_b1_h1)=0.
              pfyy1 (i,jte_b1_h1)=0.
              pfnex1(i,jte_b1_h1)=0.
              pfnwx1(i,jte_b1_h1)=0.
              pfney1(i,jte_b1_h1)=0.
              pfnwy1(i,jte_b1_h1)=0.
            enddo
          endif
        endif
!-----------------------------------------------------------------------
        do j=jts_b1,jte_b1_h1
          do i=its_b1,ite_b1_h1
            ufxx1 (i,j)=(u1d(i  ,j  )-u1d(i-1,j  ))*pfxx1 (i,j)
            ufyx1 (i,j)=(u1d(i  ,j  )-u1d(i  ,j-1))*pfyx1 (i,j)
            ufnex1(i,j)=(u1d(i  ,j  )-u1d(i-1,j-1))*pfnex1(i,j)
            ufnwx1(i,j)=(u1d(i-1,j  )-u1d(i  ,j-1))*pfnwx1(i,j)
!
            ufxy1 (i,j)=(u1d(i  ,j  )-u1d(i-1,j  ))*pfxy1 (i,j)
            ufyy1 (i,j)=(u1d(i  ,j  )-u1d(i  ,j-1))*pfyy1 (i,j)
            ufney1(i,j)=(u1d(i  ,j  )-u1d(i-1,j-1))*pfney1(i,j)
            ufnwy1(i,j)=(u1d(i-1,j  )-u1d(i  ,j-1))*pfnwy1(i,j)
!
            vfxx1 (i,j)=(v1d(i  ,j  )-v1d(i-1,j  ))*pfxx1 (i,j)
            vfyx1 (i,j)=(v1d(i  ,j  )-v1d(i  ,j-1))*pfyx1 (i,j)
            vfnex1(i,j)=(v1d(i  ,j  )-v1d(i-1,j-1))*pfnex1(i,j)
            vfnwx1(i,j)=(v1d(i-1,j  )-v1d(i  ,j-1))*pfnwx1(i,j)
!
            vfxy1 (i,j)=(v1d(i  ,j  )-v1d(i-1,j  ))*pfxy1 (i,j)
            vfyy1 (i,j)=(v1d(i  ,j  )-v1d(i  ,j-1))*pfyy1 (i,j)
            vfney1(i,j)=(v1d(i  ,j  )-v1d(i-1,j-1))*pfney1(i,j)
            vfnwy1(i,j)=(v1d(i-1,j  )-v1d(i  ,j-1))*pfnwy1(i,j)
          enddo
        enddo
!-----------------------------------------------------------------------
!---advection of u1d and v1d--------------------------------------------
!-----------------------------------------------------------------------
        if(adv_standard) then
          ibeg=max(its,ids+1+lnsad)
          iend=min(ite,ide-2-lnsad)
          jbeg=max(jts,jds+1+lnsad)
          jend=min(jte,jde-2-lnsad)
!
          do j=jbeg,jend
            dxody=dxv(j)*rdyv
            dyodx=dyv   *rdxv(j)
            fadp=fad(j)
            do i=ibeg,iend
!
              fdpp=fadp/(dsg2(l)*pdop(i,j)+pdsg1(l))
!
              dux1=(ufnex1(i+1,j+1)+ufnex1(i  ,j  ) &
                   +ufnwx1(i  ,j+1)+ufnwx1(i+1,j  ))*0.25 &
                  +(ufxx1 (i  ,j  )+ufxx1 (i+1,j  ) &
                   +ufyx1 (i  ,j  )+ufyx1 (i  ,j+1))
!
              dvx1=(vfnex1(i+1,j+1)+vfnex1(i  ,j  ) &
                   +vfnwx1(i  ,j+1)+vfnwx1(i+1,j  ))*0.25 &
                  +(vfxx1 (i  ,j  )+vfxx1 (i+1,j  ) &
                   +vfyx1 (i  ,j  )+vfyx1 (i  ,j+1))
!
              duy1=(ufney1(i+1,j+1)+ufney1(i  ,j  ) &
                   +ufnwy1(i  ,j+1)+ufnwy1(i+1,j  ))*0.25 &
                  +(ufxy1 (i  ,j  )+ufxy1 (i+1,j  ) &
                   +ufyy1 (i  ,j  )+ufyy1 (i  ,j+1))
!
              dvy1=(vfney1(i+1,j+1)+vfney1(i  ,j  ) &
                   +vfnwy1(i  ,j+1)+vfnwy1(i+1,j  ))*0.25 &
                  +(vfxy1 (i  ,j  )+vfxy1 (i+1,j  ) &
                   +vfyy1 (i  ,j  )+vfyy1 (i  ,j+1))
!
              tcu(i,j,l)=((dvx1-dvy1)*dxody+(dux1+duy1))*fdpp &
                        +tcu(i,j,l)
              tcv(i,j,l)=((dux1-duy1)*dyodx+(dvx1+dvy1))*fdpp &
                        +tcv(i,j,l)
            enddo
          enddo
        endif
!-----------------------------------------------------------------------
!---regional branch-----------------------------------------------------
!-----------------------------------------------------------------------
!
        if(.not.global.and.adv_upstream) then
!
!-----------------------------------------------------------------------
!---upstream advection of momentum along lateral boundaries-------------
!-----------------------------------------------------------------------
!
          envp=-dt*rdyv
!
!***  Upstream advection along southern rows.
!
          jbeg=max(jts,jds+1)
          jend=min(jte,jds+lnsad)
!
          do j=jbeg,jend
            emvp=-dt*rdxv(j)
            do i=its_b1,ite_b2
              pp=u(i,j,l)*emvp
              qq=v(i,j,l)*envp
!
              if(pp.le.0.) then
                iap=-1
                pp=-pp
              else
                iap=1
              endif
!
              if(qq.le.0.) then
                jap=-1
                qq=-qq
              else
                jap=1
              endif
!
              tcu(i,j,l)=(u(i+iap,j,l)-u(i,j,l))*pp &
                        +(u(i,j+jap,l)-u(i,j,l))*qq &
                        +(u(i,j,l)-u(i+iap,j,l) &
                         -u(i,j+jap,l)+u(i+iap,j+jap,l))*pp*qq &
                        +tcu(i,j,l)
              tcv(i,j,l)=(v(i+iap,j,l)-v(i,j,l))*pp &
                        +(v(i,j+jap,l)-v(i,j,l))*qq &
                        +(v(i,j,l)-v(i+iap,j,l) &
                         -v(i,j+jap,l)+v(i+iap,j+jap,l))*pp*qq &
                        +tcv(i,j,l)
            enddo
          enddo
!
!***  Upstream advection along northern rows.
!
          do j=max(jts,jde-1-lnsad),jte_b2
            emvp=-dt*rdxv(j)
            do i=its_b1,ite_b2
              pp=u(i,j,l)*emvp
              qq=v(i,j,l)*envp
!
              if(pp.le.0.) then
                iap=-1
                pp=-pp
              else
                iap=1
              endif
!
              if(qq.le.0.) then
                jap=-1
                qq=-qq
              else
                jap=1
              endif
!
              tcu(i,j,l)=(u(i+iap,j,l)-u(i,j,l))*pp &
                        +(u(i,j+jap,l)-u(i,j,l))*qq &
                        +(u(i,j,l)-u(i+iap,j,l) &
                         -u(i,j+jap,l)+u(i+iap,j+jap,l))*pp*qq &
                        +tcu(i,j,l)
              tcv(i,j,l)=(v(i+iap,j,l)-v(i,j,l))*pp &
                        +(v(i,j+jap,l)-v(i,j,l))*qq &
                        +(v(i,j,l)-v(i+iap,j,l) &
                         -v(i,j+jap,l)+v(i+iap,j+jap,l))*pp*qq &
                        +tcv(i,j,l)
            enddo
          enddo
!
!***  Upstream advection along western rows.
!
          do j=max(jts,jds+1+lnsad),min(jte,jde-2-lnsad)
            emvp=-dt*rdxv(j)
            do i=its_b1,ids+lnsad
              pp=u(i,j,l)*emvp
              qq=v(i,j,l)*envp
!
              if(pp.le.0.) then
                iap=-1
                pp=-pp
              else
                iap=1
              endif
!
              if(qq.le.0.) then
                jap=-1
                qq=-qq
              else
                jap=1
              endif
!
              tcu(i,j,l)=(u(i+iap,j,l)-u(i,j,l))*pp &
                        +(u(i,j+jap,l)-u(i,j,l))*qq &
                        +(u(i,j,l)-u(i+iap,j,l) &
                         -u(i,j+jap,l)+u(i+iap,j+jap,l))*pp*qq &
                        +tcu(i,j,l)
              tcv(i,j,l)=(v(i+iap,j,l)-v(i,j,l))*pp &
                        +(v(i,j+jap,l)-v(i,j,l))*qq &
                        +(v(i,j,l)-v(i+iap,j,l) &
                         -v(i,j+jap,l)+v(i+iap,j+jap,l))*pp*qq &
                        +tcv(i,j,l)
            enddo
          enddo
!
!***  Upstream advection along eastern rows.
!
          do j=max(jts,jds+1+lnsad),min(jte,jde-2-lnsad)
            emvp=-dt*rdxv(j)
            do i=max(its,ide-1-lnsad),ite_b2
              pp=u(i,j,l)*emvp
              qq=v(i,j,l)*envp
!
              if(pp.le.0.) then
                iap=-1
                pp=-pp
              else
                iap=1
              endif
!
              if(qq.le.0.) then
                jap=-1
                qq=-qq
              else
                jap=1
              endif
!
              tcu(i,j,l)=(u(i+iap,j,l)-u(i,j,l))*pp &
                        +(u(i,j+jap,l)-u(i,j,l))*qq &
                        +(u(i,j,l)-u(i+iap,j,l) &
                         -u(i,j+jap,l)+u(i+iap,j+jap,l))*pp*qq &
                        +tcu(i,j,l)
              tcv(i,j,l)=(v(i+iap,j,l)-v(i,j,l))*pp &
                        +(v(i,j+jap,l)-v(i,j,l))*qq &
                        +(v(i,j,l)-v(i+iap,j,l) &
                         -v(i,j+jap,l)+v(i+iap,j+jap,l))*pp*qq &
                        +tcv(i,j,l)
            enddo
          enddo
!-----------------------------------------------------------------------
        endif ! regional lateral boundaries
!-----------------------------------------------------------------------
!
      enddo vertical_loop2
!.......................................................................
!$omp end parallel do
!.......................................................................
!
!-----------------------------------------------------------------------
!
                        endsubroutine adv1
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
                        subroutine vtoa &
(lm &
,dt,ef4t,pt &
,sg2 &
,psdt &
,dwdt,rtop &
,omgalf &
,pint &
!---temporary arguments-------------------------------------------------
,tdiv,tct)

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
implicit none
!-----------------------------------------------------------------------
integer(kind=kint),intent(in):: &
 lm                          ! total # of levels

real(kind=kfpt),intent(in):: &
 dt &                        ! dynamics time step
,ef4t &                      ! vertical grid parameter
,pt                          ! pressure at the top of the model's atmosphere

real(kind=kfpt),dimension(1:lm+1),intent(in):: &
 sg2                         ! delta sigmas

real(kind=kfpt),dimension(ims:ime,jms:jme),intent(in):: &
 psdt                        ! surface pressure tendency

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(in):: &
 dwdt &                      ! nonhydrostatic correction factor
,rtop                        ! RT/p

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(inout):: &
 omgalf                      !

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm+1),intent(inout):: &
 pint                        ! pressure at interfaces
!-----------------------------------------------------------------------
!---temporary arguments-------------------------------------------------
!-----------------------------------------------------------------------
real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(in):: &
 tdiv                        ! integrated horizontal mass divergence

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(inout):: &
 tct                         ! time change of temperature
!-----------------------------------------------------------------------
!--local variables------------------------------------------------------
!-----------------------------------------------------------------------
integer(kind=kint):: &
 i &                         ! index in x direction
,j &                         ! index in y direction
,l                           ! index in p direction

real(kind=kfpt):: &
 dwdtp &                     ! nonhydrostatic correction factor at the point
,toa &                       ! omega-alpha temporary
,tpmp                        ! pressure temporary at the point

real(kind=kfpt),dimension(its:ite,jts:jte):: &
 tpm                         ! pressure temporary
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------
#ifdef ENABLE_SMP
!-----------------
!.......................................................................
!$omp parallel private(dwdtp,i,j,jstart,jstop,l,nth,tid,toa,tpmp)
!.......................................................................
      nth = omp_get_num_threads()
      tid = omp_get_thread_num()
      call looplimits(tid, nth, jts_b1,jte_b1,jstart,jstop)
!-----------------
#else
!-----------------
      jstart = jts_b1
      jstop  = jte_b1
!-----------------
#endif 
!-----------------


!do l=1,lm
!  do j=jstart,jstop
!    do i=its_b1,ite_b1
!      tct(i,j,l)=0.
!    enddo
!  enddo
!enddo


      do j=jstart,jstop
        do i=its_b1,ite_b1
          pint(i,j,1)=pt
          tpm(i,j)=pt+pint(i,j,2)
!
          dwdtp=dwdt(i,j,1)
!
          tpmp=pint(i,j,2)+pint(i,j,3)
          toa=-tdiv(i,j,1)*rtop(i,j,1)*dwdtp*ef4t
!
          omgalf(i,j,1)=omgalf(i,j,1)+toa
          tct(i,j,1)=tct(i,j,1)+toa
!
          pint(i,j,2)=psdt(i,j)*(sg2(1)+sg2(2))*dwdtp*dt &
                     +tpm(i,j)-pint(i,j,1)
!
          tpm(i,j)=tpmp
        enddo
      enddo
!
      do l=2,lm-1
        do j=jstart,jstop
          do i=its_b1,ite_b1
            dwdtp=dwdt(i,j,l)
!
            tpmp=pint(i,j,l+1)+pint(i,j,l+2)
            toa=-(tdiv(i,j,l-1)+tdiv(i,j,l))*rtop(i,j,l)*dwdtp*ef4t
!
            omgalf(i,j,l)=omgalf(i,j,l)+toa 
            tct(i,j,l)=tct(i,j,l)+toa
!
            pint(i,j,l+1)=psdt(i,j)*(sg2(l)+sg2(l+1))*dwdtp*dt &
                         +tpm(i,j)-pint(i,j,l)
!
            tpm(i,j)=tpmp
          enddo
        enddo
      enddo
!
      do j=jstart,jstop
        do i=its_b1,ite_b1
          dwdtp=dwdt(i,j,lm)
!
          toa=-(tdiv(i,j,lm-1)+tdiv(i,j,lm))*rtop(i,j,lm)*dwdtp*ef4t
!
          omgalf(i,j,lm)=omgalf(i,j,lm)+toa
          tct(i,j,lm)=tct(i,j,lm)+toa
!
          pint(i,j,lm+1)=psdt(i,j)*(sg2(lm)+sg2(lm+1))*dwdtp*dt &
                        +tpm(i,j)-pint(i,j,lm)
        enddo
      enddo
!-----------------
#ifdef ENABLE_SMP
!-----------------
!.......................................................................
!$omp end parallel
!.......................................................................
!-----------------
#endif
!-----------------

!-----------------------------------------------------------------------
!
                        endsubroutine vtoa
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
                        subroutine updates &
(lm,ntracers,kss,kse,s &
!---temporary arguments-------------------------------------------------
,tcs)
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
implicit none
!-----------------------------------------------------------------------
integer(kind=kint),intent(in):: &
 lm &                        ! total # of levels
,ntracers &                  ! total # of tracers
,kss &                       ! starting index of tracers to be updated
,kse                         ! ending index of tracers to be updated

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm,1:ntracers),intent(inout):: &
 s                           ! tracer
!-----------------------------------------------------------------------
!---temporary arguments-------------------------------------------------
!-----------------------------------------------------------------------
real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm,1:ntracers):: &
 tcs                         ! tracer time change
!-----------------------------------------------------------------------
!--local variables------------------------------------------------------
!-----------------------------------------------------------------------
integer(kind=kint):: &
 i &                         ! index in x direction
,j &                         ! index in y direction
,l &                         ! index in p direction
,k                           ! tracers index
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!.......................................................................
!$omp parallel do private(i,j,l,k)
!.......................................................................
      do k=kss,kse
      do l=1,lm
        do j=jts_b1,jte_b1
          do i=its_b1,ite_b1
            s(i,j,l,k)=s(i,j,l,k)+tcs(i,j,l,k)
!
            tcs(i,j,l,k)=0.
          enddo
        enddo
      enddo
      enddo
!.......................................................................
!$omp end parallel do
!.......................................................................
!-----------------------------------------------------------------------
!
                        endsubroutine updates
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
                        subroutine updatet &
(lm,t &
!---temporary arguments-------------------------------------------------
,tct)
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
implicit none
!-----------------------------------------------------------------------
integer(kind=kint),intent(in):: &
 lm                          ! total # of levels

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(inout):: &
 t                           ! temperature
!-----------------------------------------------------------------------
!---temporary arguments-------------------------------------------------
!-----------------------------------------------------------------------
real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm):: &
 tct                         ! temperature time change
!-----------------------------------------------------------------------
!--local variables------------------------------------------------------
!-----------------------------------------------------------------------
integer(kind=kint):: &
 i &                         ! index in x direction
,j &                         ! index in y direction
,l                           ! index in p direction
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!.......................................................................
!$omp parallel do private(i,j,l)
!.......................................................................
      do l=1,lm
        do j=jts_b1,jte_b1
          do i=its_b1,ite_b1
            t(i,j,l)=t(i,j,l)+tct(i,j,l)
!
            tct(i,j,l)=0.
          enddo
        enddo
      enddo
!.......................................................................
!$omp end parallel do
!.......................................................................
!-----------------------------------------------------------------------
!
                        endsubroutine updatet
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++                            
!-----------------------------------------------------------------------
                        subroutine updateuv &
(lm,u,v &
!---temporary arguments-------------------------------------------------
,tcu,tcv)
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
implicit none
!-----------------------------------------------------------------------
real(kind=kfpt),parameter:: &
 cflfc=1./(140.*140.)        ! cfl limit
           
integer(kind=kint),intent(in):: &
 lm                          ! total # of levels
           
real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(inout):: &
 u &                         ! u
,v                           ! v
!-----------------------------------------------------------------------
!---temporary arguments-------------------------------------------------
!-----------------------------------------------------------------------
real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm):: &
 tcu &                       ! u time change
,tcv                         ! v time change
!-----------------------------------------------------------------------
!--local variables------------------------------------------------------
!-----------------------------------------------------------------------
integer(kind=kint):: &
 i &                         ! index in x direction
,j &                         ! index in y direction
,l                           ! index in p direction
           
real(kind=kfpt):: &
 cflc &                      !
,rcflc                       !
           
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!.......................................................................
!$omp parallel do private (cflc,i,j,l,rcflc)
!.......................................................................
      do l=1,lm
        do j=jts_b1,jte_b2
          do i=its_b1,ite_b2
            u(i,j,l)=u(i,j,l)+tcu(i,j,l)
            v(i,j,l)=v(i,j,l)+tcv(i,j,l)
!
            tcu(i,j,l)=0.
            tcv(i,j,l)=0.
!
            cflc=(u(i,j,l)*u(i,j,l)+v(i,j,l)*v(i,j,l))*cflfc
            if(cflc.gt.1.) then
              rcflc=sqrt(1./cflc)
              u(i,j,l)=u(i,j,l)*rcflc
              v(i,j,l)=v(i,j,l)*rcflc
            endif
          enddo
        enddo
      enddo
!.......................................................................
!$omp end parallel do
!.......................................................................
!-----------------------------------------------------------------------
!
                        endsubroutine updateuv
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++                            
!-----------------------------------------------------------------------
                        subroutine hdiff &
(global,hydro &
,inpes,jnpes,lm,lpt2 &
,dyh,rdyh &
,epsq2 &
,dxv,rare,rdxh &
,sice,sm &
,hdacx,hdacy,hdacvx,hdacvy &
,w,z &
,cw,q,q2,t,u,v,def3d)
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
implicit none
!-----------------------------------------------------------------------
real(kind=kfpt),parameter:: &
 defvfc=2. &                 ! vertical deformation weight
,scq2=50. &                  ! 2tke weighting factor
,epsq=1.e-20 &               ! floor value for specific humidity
,slopec=.05                  ! critical slope
           
!-----------------------------------------------------------------------
logical(kind=klog),intent(in):: &
 global &                    ! global or regional
,hydro                       ! logical switch for nonhydrostatic dynamics
           
integer(kind=kint),intent(in):: &
 inpes &                     ! w-e # of subdomains
,jnpes &                     ! n-s # of subdomains
,lm &                        ! total # of levels
,lpt2                        ! # of levels in the pressure range

real(kind=kfpt),intent(in):: &
 dyh &                       !
,rdyh                        ! 1/deltay
 
real(kind=kfpt),dimension(1:lm),intent(in):: &
 epsq2                       ! floor value of 2tke

real(kind=kfpt),dimension(jds:jde),intent(in):: &
 dxv &                       !
,rare &                      !
,rdxh                        ! 1/deltax
           
real(kind=kfpt),dimension(ims:ime,jms:jme),intent(in):: &
 hdacx &                     ! exchange coefficient for mass points
,hdacy &                     ! exchange coefficient for mass points
,hdacvx &                    ! exchange coefficient for velocity points
,hdacvy &                    ! exchange coefficient for velocity points
,sice &                      ! sea-ice mask
,sm                          ! sea mask
           
real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(in):: &
 w &                         ! w wind component
,z                           ! height at mass points
           
real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(inout):: &
 cw &                        ! condensate
,q &                         ! specific humidity
,q2 &                        ! 2tke
,t &                         ! temperature
,u &                         ! u wind component
,v                           ! v wind component

real(kind=kfpt),dimension(ims:ime,jms:jme,lm),intent(out) :: &
 def3d          
           
!-----------------------------------------------------------------------
!--local variables------------------------------------------------------
!-----------------------------------------------------------------------
logical(kind=klog):: &
 cilinx &                    ! coast/ice line
,ciliny                      ! coast/ice line
 
integer(kind=kint):: &
 i &                         ! index in x direction
,j &                         ! index in y direction
,l                           ! index in p direction
           
real(kind=kfpt):: &
 def1 &                      ! component of deformation
,def2 &                      ! component of deformation
,def3 &                      ! component of deformation
,def4 &                      ! component of deformation
,defc &                      ! deformation floor
,defm &                      ! deformation cap
,defp &                      ! deformation at the point
,defs &                      ! component of deformation
,deft &                      ! component of deformation
,defvp &
,defhp &
,facdif &                    ! diffusion factor
,hkfx &                      ! def with slope factor
,hkfy &                      ! def with slope factor
,q2trm &                     !
,slopx &                     ! x slope
,slopy                       ! y slope

real(kind=kfpt),dimension(ims:ime,jms:jme):: &
 cdif &                      ! condensate 2nd order diffusion
,cx &                        ! condensate difference, x direction
,cy &                        ! condensate difference, y direction
,def &                       ! deformation (in halo exchange)
,fmlx &                      ! x slope mask
,fmly &                      ! y slope mask
,hkx &                       ! deformation sum, x direction
,hky &                       ! deformation sum, y direction
,qdif &                      ! specific humidity 2nd order diffusion
,qx &                        ! specific humidity difference, x direction
,qy &                        ! specific humidity difference, y direction
,q2dif &                     ! 2tke 2nd order diffusion
,q2x &                       ! 2tke difference, x direction
,q2y &                       ! 2tke difference, y direction
,tdif &                      ! temperature 2nd order diffusion
,tx &                        ! temperature difference, x direction
,ty &                        ! temperature difference, y direction
,udif &                      ! u wind component 2nd order diffusion
,ux &                        ! u wind component difference, x direction
,uy &                        ! u wind component difference, y direction
,vdif &                      ! v wind component 2nd order diffusion
,vx &                        ! v wind component difference, x direction
,vy                          ! v wind component difference, y direction

!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
      if(global) then
        defm=1.6875e-4 ! 1.35e-3/2./4. ! deformation cap, /4. for smag2=0.4
        facdif=1.0
      else
        defm=9999.
        facdif=1.0   !bsf: was 4.0 in 10/18/2011 NAM implementation
      endif
!-----------------------------------------------------------------------
!
      do l=1,lm
       if(s_bdy)then
         do i=ims,ime
           def3d(i,jds,l)=0.
         enddo
       endif
!
       if(n_bdy)then
         do i=ims,ime
           def3d(i,jde,l)=0.
         enddo
       endif
!
       if(w_bdy)then
         do j=jms,jme
           def3d(ids,j,l)=0.
         enddo
       endif
!
       if(e_bdy)then
         do j=jms,jme
           def3d(ide,j,l)=0.
         enddo
       endif
     enddo
!
!-----------------------------------------------------------------------
!---grand vertical loop-------------------------------------------------
!-----------------------------------------------------------------------
!
!.......................................................................
!$omp parallel do &
!$omp private(def1,def2,def3,def4,defc,defp,defs,deft,i,j,l,q2trm)
!.......................................................................
!-----------------------------------------------------------------------
!
      vertical_loop_1: do l=1,lm
!
!-----------------------------------------------------------------------
        if(global)then
          if(l.gt.lpt2/3) then
            defc=0.
          else
            defc=defm*0.01
          endif
        else ! regional
          if (l.le.3) then
            defc=0.015*(3-l+1.)/3.
          else
            defc=0.
          endif
        endif
!-----------------------------------------------------------------------
!        do j=jts_h1,jte_h2
!          do i=its_h1,ite_h2
!            q2(i,j,l)=max(q2(i,j,l),epsq2)
!          enddo
!        enddo
!-----------------------------------------------------------------------
        do j=jts_b1_h1,jte_b1_h2
          do i=its_b1_h1,ite_b1_h2
            deft=((u(i  ,j-1,l)+u(i  ,j  ,l) &
                  -u(i-1,j-1,l)-u(i-1,j  ,l))*dyh &
                 -(v(i-1,j  ,l)+v(i  ,j  ,l))*dxv(j  ) &
                 +(v(i-1,j-1,l)+v(i  ,j-1,l))*dxv(j-1))*rare(j)
            defs=((u(i-1,j  ,l)+u(i  ,j  ,l))*dxv(j  ) &
                 -(u(i-1,j-1,l)+u(i  ,j-1,l))*dxv(j-1) &
                 +(v(i  ,j-1,l)+v(i  ,j  ,l) &
                  -v(i-1,j-1,l)-v(i-1,j  ,l))*dyh     )*rare(j)
!            defz=(-(u(i-1,j  ,l)+u(i  ,j  ,l))*dxv(j  ) &
!                  +(u(i-1,j-1,l)+u(i  ,j-1,l))*dxv(j-1) &
!                  +(v(i  ,j-1,l)+v(i  ,j  ,l) &
!                   -v(i-1,j-1,l)-v(i-1,j  ,l))*dyh     )*rare(j)  *10.
!            if(defz.gt.0.) defz=0.
!
            if(.not.hydro) then
              def1=(w(i,j,l)-w(i,j-1,l))*rdyh
              def2=(w(i,j,l)-w(i-1,j,l))*rdxh(j)
              def3=(w(i+1,j,l)-w(i,j,l))*rdxh(j)
              def4=(w(i,j+1,l)-w(i,j,l))*rdyh
            else
              def1=0.
              def2=0.
              def3=0.
              def4=0.
            endif

            if(q2(i,j,l).gt.epsq2(L)) then
              q2trm=scq2*q2(i,j,l)*rare(j)
            else
              q2trm=0.
            endif


            defhp=deft*deft+defs*defs+q2trm
            defvp=def1*def1+def2*def2+def3*def3+def4*def4 
!
            defp=sqrt(defhp*2.+defvp*defvfc)
            defp=max (defp,defc)
            defp=min (defp,defm)
!
            def3d(i,j,l)=defp
          enddo
        enddo
!-----------------------------------------------------------------------
!
      enddo vertical_loop_1
!.......................................................................
!$omp end parallel do
!.......................................................................
!
!-----------------------------------------------------------------------
!
      if(global) then
        btim=timef()
        call swaphn(def3d,ims,ime,jms,jme,lm,inpes)
        timers(my_domain_id)%swaphn_tim=timers(my_domain_id)%swaphn_tim+(timef()-btim)
!
        btim=timef()
        call polehn(def3d,ims,ime,jms,jme,lm,inpes,jnpes)
        timers(my_domain_id)%polehn_tim=timers(my_domain_id)%polehn_tim+(timef()-btim)
      endif
!
!-----------------------------------------------------------------------
!.......................................................................
!$omp parallel do &
!$omp private(cdif,cilinx,ciliny,cx,cy,fmlx,fmly, &
!$omp         hkfx,hkfy,hkx,hky,i,j,l, &
!$omp         q2dif,q2x,q2y,qdif,qx,qy,slopx,slopy, &
!$omp         tdif,tx,ty,udif,ux,uy,vdif,vx,vy)
!.......................................................................
!-----------------------------------------------------------------------
      vertical_loop_3: do l=1,lm
!-----------------------------------------------------------------------
        if(l.gt.lpt2.and..not.hydro) then
          do j=jts_b1,jte_h2
            do i=its_b1,ite_h2
              cilinx=sice(i-1,j).ne.sice(i,j) &
                 .or.sm  (i-1,j).ne.sm  (i,j)
              ciliny=sice(i,j-1).ne.sice(i,j) &
                 .or.sm  (i,j-1).ne.sm  (i,j)
              slopx=abs((z(i,j,l)-z(i-1,j,l))*rdxh(j))
              slopy=abs((z(i,j,l)-z(i,j-1,l))*rdyh   )
!
              if(slopx.le.slopec.or.cilinx) then
                fmlx(i,j)=1.
              else
                fmlx(i,j)=0.
              endif
!
              if(slopy.le.slopec.or.ciliny) then
                fmly(i,j)=1.
              else
                fmly(i,j)=0.
              endif
!
            enddo
          enddo
        else
          do j=jts_b1,jte_h2
            do i=its_b1,ite_h2
              fmlx(i,j)=1.
              fmly(i,j)=1.
            enddo
          enddo
        endif
!-----------------------------------------------------------------------
!---contributions behind mass points------------------------------------
!-----------------------------------------------------------------------
        do j=jts_b1,jte_h2
          do i=its_b1,ite_h2
            hkx(i,j)=(def3d(i-1,j,l)+def3d(i,j,l))
            hky(i,j)=(def3d(i,j-1,l)+def3d(i,j,l))
            hkfx=hkx(i,j)*fmlx(i,j)
            hkfy=hky(i,j)*fmly(i,j)
!
            tx (i,j)=(t (i,j,l)-t (i-1,j,l))*hkfx
            qx (i,j)=(q (i,j,l)-q (i-1,j,l))*hkfx
            cx (i,j)=(cw(i,j,l)-cw(i-1,j,l))*hkfx
            q2x(i,j)=(q2(i,j,l)-q2(i-1,j,l))*hkfx
!
            ty (i,j)=(t (i,j,l)-t (i,j-1,l))*hkfy
            qy (i,j)=(q (i,j,l)-q (i,j-1,l))*hkfy
            cy (i,j)=(cw(i,j,l)-cw(i,j-1,l))*hkfy
            q2y(i,j)=(q2(i,j,l)-q2(i,j-1,l))*hkfy
          enddo
        enddo
!-----------------------------------------------------------------------
!---u,v, contributions, behind v points---------------------------------
!-----------------------------------------------------------------------
        do j=jts_b1,jte_b1_h1
          do i=its_b1,ite_b1_h1
            ux(i,j)=(u(i,j,l)-u(i-1,j,l))*hky(i,j+1)
            vx(i,j)=(v(i,j,l)-v(i-1,j,l))*hky(i,j+1)
            uy(i,j)=(u(i,j,l)-u(i,j-1,l))*hkx(i+1,j)
            vy(i,j)=(v(i,j,l)-v(i,j-1,l))*hkx(i+1,j)
          enddo
        enddo
!-----------------------------------------------------------------------
        do j=jts_b1,jte_b1
          do i=its_b1,ite_b1
            tdif (i,j)=(tx (i+1,j)-tx (i,j))*hdacx(i,j) &
                      +(ty (i,j+1)-ty (i,j))*hdacy(i,j)
            qdif (i,j)=(qx (i+1,j)-qx (i,j))*hdacx(i,j) &
                      +(qy (i,j+1)-qy (i,j))*hdacy(i,j)
            cdif (i,j)=(cx (i+1,j)-cx (i,j))*hdacx(i,j) &
                      +(cy (i,j+1)-cy (i,j))*hdacy(i,j)
            q2dif(i,j)=(q2x(i+1,j)-q2x(i,j))*hdacx(i,j) &
                      +(q2y(i,j+1)-q2y(i,j))*hdacy(i,j)
          enddo
        enddo
!-----------------------------------------------------------------------
        do j=jts_b1,jte_b2
          do i=its_b1,ite_b2
            udif(i,j)=(ux(i+1,j)-ux(i,j))*hdacvx(i,j) &
                     +(uy(i,j+1)-uy(i,j))*hdacvy(i,j)
            vdif(i,j)=(vx(i+1,j)-vx(i,j))*hdacvx(i,j) &
                     +(vy(i,j+1)-vy(i,j))*hdacvy(i,j)
          enddo
        enddo
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-------------2-nd order diffusion--------------------------------------
!-----------------------------------------------------------------------
          do j=jts_b1,jte_b1
            do i=its_b1,ite_b1
              t (i,j,l)=t (i,j,l)+tdif (i,j)
!-- Enhanced diffusion for Q & CWM if facdif>1    !bsf:
              q (i,j,l)=q (i,j,l)+facdif*qdif (i,j)
              cw(i,j,l)=cw(i,j,l)+facdif*cdif (i,j)
              q2(i,j,l)=q2(i,j,l)+q2dif(i,j)
            enddo
          enddo
          do j=jts_b1,jte_b2
            do i=its_b1,ite_b2
              u(i,j,l)=u(i,j,l)+udif(i,j)
              v(i,j,l)=v(i,j,l)+vdif(i,j)
            enddo
          enddo
!-----------------------------------------------------------------------
!
        enddo vertical_loop_3
!.......................................................................
!$omp end parallel do
!.......................................................................
!
!-----------------------------------------------------------------------
!
                        endsubroutine hdiff
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
                        subroutine cdzdt &
(global,hydro &
,lm &
,dt &
,dsg2,pdsg1 &
,fah &
,fis,pd,pdo &
,psgdt &
,cw,q,rtop,t &
,pint &
,dwdt,pdwdt,w,baro &
,z &
!---temporary arguments-------------------------------------------------
,pfne,pfnw,pfx,pfy)
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
implicit none
!-----------------------------------------------------------------------
real(kind=kfpt),parameter:: &
 cfc=1.533 &                 ! adams-bashforth positioning in time
,bfc=1.-cfc                  ! adams bashforth positioning in time
!-----------------------------------------------------------------------
logical(kind=klog),intent(in):: &
 global &                    ! global or regional
,hydro                       ! hydrostatic or nonhydrostatic
           
integer(kind=kint),intent(in):: &
 lm                          ! total # of levels
           
real(kind=kfpt),intent(in):: &
 dt                          ! dynamics time step
           
real(kind=kfpt),dimension(1:lm),intent(in):: &
 dsg2 &                      ! delta sigmas
,pdsg1                       ! delta pressures
           
real(kind=kfpt),dimension(jds:jde),intent(in):: &
 fah                         !
           
real(kind=kfpt),dimension(ims:ime,jms:jme),intent(in):: &
 fis &                       ! surface geopotential
,pd &                        ! sigma range pressure difference
,pdo                         ! sigma range pressure difference
      
real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm-1),intent(in):: &
 psgdt                       ! vertical mass flux
           
real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(in):: &
 cw &                        ! condensate
,q &                         ! specific humidity
,rtop &                      ! rt/p
,t                           ! temperature
           
real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm+1),intent(in):: &
 pint                        ! pressure at interfaces
           
real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(inout):: &
 dwdt &                      ! nonhydrostatic correction factor
,pdwdt &                     ! previous nonhydrostatic correction factor
,w                           ! w wind component
           
real(kind=kfpt),dimension(ims:ime,jms:jme),intent(out):: &
 baro

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(out):: &
 z                           ! heights in the middle of the layers
!-----------------------------------------------------------------------
!---temporary arguments-------------------------------------------------
!-----------------------------------------------------------------------
real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(in) :: &
 pfne &                      ! mass flux, ne direction
,pfnw &                      ! mass flux, nw direction
,pfx &                       ! mass flux, x direction
,pfy                         ! mass flux, y direction
!-----------------------------------------------------------------------
!--local variables------------------------------------------------------
!-----------------------------------------------------------------------
integer(kind=kint):: &
 i &                         ! index in x direction
,j &                         ! index in y direction
,l                           ! index in p direction
           
real(kind=kfpt):: &
 dpup &                      !
,dw &                        ! w wind component increment
,dz &                        ! layer thickness
,fahp &                      !
,rg &                        !
,rdt &                       !
,trog &                      !
,wup &                       ! w wind component at upper interface
,zup                         ! height of upper interface
           
real(kind=kfpt),dimension(its_b1:ite_h2,jts_b1:jte_h2):: &
 zne &                       ! height flux, ne direction
,znw &                       ! height flux, nw direction
,zx &                        ! height flux, x direction
,zy                          ! height flux, y direction
           
real(kind=kfpt),dimension(its_b1:ite_b1,jts_b1:jte_b1):: &
 tta &                       ! advection through upper interface
,ttb                         ! advection through lower interface
           
real(kind=kfpt),dimension(its_h1:ite_h2,jts_h1:jte_h2):: &
 wlo &                       ! w wind component at lower interface
,zlo                         ! height at lower interface
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      rg=1./g
      rdt=1./dt
      trog=2.*r/g
!
!-----------------------------------------------------------------------
!-----------------
#ifdef ENABLE_SMP
!-----------------
!.......................................................................
!$omp parallel private(dw,dz,i,j,l,jstart,jstop,nth,tid,wup,zup)
        nth = omp_get_num_threads()
        tid = omp_get_thread_num()
        call looplimits(tid, nth, jts_h1,jte_h2,jstart,jstop)
!.......................................................................
!-----------------
#else
!-----------------
      jstart = jts_h1
      jstop  = jte_h2
!-----------------
#endif 
!-----------------
      do j=jstart,jstop
        do i=its_h1,ite_h2
          wlo(i,j)=0.
          zlo(i,j)=fis(i,j)*rg
        enddo
      enddo
!-----------------------------------------------------------------------
!---nonhydrostatic equation---------------------------------------------
!-----------------------------------------------------------------------
      do l=lm,1,-1
        do j=jstart,jstop
          do i=its_h1,ite_h2
            pdwdt(i,j,l)=dwdt(i,j,l)
            dwdt(i,j,l)=w(i,j,l)
!
            dz=(q(i,j,l)*0.608+(1.-cw(i,j,l)))*t(i,j,l)*trog &
              *(dsg2(l)*pd(i,j)+pdsg1(l))/(pint(i,j,l)+pint(i,j,l+1))
            dw=(dz-rtop(i,j,l)*(dsg2(l)*pdo(i,j)+pdsg1(l))*rg)*rdt
!
            zup=zlo(i,j)+dz
            wup=wlo(i,j)+dw
            z(i,j,l)=dz*0.5+zlo(i,j)
            w(i,j,l)=dw*0.5+wlo(i,j)
!
            zlo(i,j)=zup
            wlo(i,j)=wup
          enddo
        enddo
      enddo
!-----------------
#ifdef ENABLE_SMP
!-----------------
!.......................................................................
!$omp end parallel
!.......................................................................
!-----------------
#endif
!-----------------
!-----------------------------------------------------------------------
!
      if(hydro) return
!
!-----------------------------------------------------------------------
!-----------------
#ifdef ENABLE_SMP
!-----------------
!.......................................................................
!$omp parallel private(i,j,l,jstart,jstop,nth,tid)
!.......................................................................
      nth = omp_get_num_threads()
      tid = omp_get_thread_num()
      call looplimits(tid, nth, jts_b1,jte_b1,jstart,jstop)
!-----------------
#else
!-----------------
      jstart = jts_b1
      jstop  = jte_b1
!-----------------
#endif 
!-----------------
!
      do j=jstart,jstop
        do i=its_b1,ite_b1
          ttb(i,j)=0.
        enddo
      enddo
!
      do l=1,lm-1
        do j=jstart,jstop
          do i=its_b1,ite_b1
            tta(i,j)=(z(i,j,l+1)-z(i,j,l))*psgdt(i,j,l)*0.5
            w(i,j,l)=(tta(i,j)+ttb(i,j)) &
                    /(dsg2(l)*pdo(i,j)+pdsg1(l)) &
                    +w(i,j,l)
            ttb(i,j)=tta(i,j)
          enddo
        enddo
      enddo
!
      do j=jstart,jstop
        do i=its_b1,ite_b1
          w(i,j,lm)=ttb(i,j)/(dsg2(lm)*pdo(i,j)+pdsg1(lm))+w(i,j,lm)
        enddo
      enddo
!
!-----------------------------------------------------------------------
!---grand horizontal loop-----------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------
#ifdef ENABLE_SMP
!-----------------
!.......................................................................
!$omp end parallel
!.......................................................................
!$omp parallel do private(dpup,fahp,i,j,l,zne,znw,zx,zy)
!.......................................................................
!-----------------
#endif
!-----------------
!-----------------------------------------------------------------------
!
      vertical_loop: do l=1,lm
!
!-----------------------------------------------------------------------
!
        dpup=pdsg1(l)
!
!-----------------------------------------------------------------------
!-------------mass fluxes, on h points----------------------------------
!-----------------------------------------------------------------------
        do j=jts_b1,jte_h1
          do i=its_b1,ite_h1
            zx(i,j)=(z(i,j,l)-z(i-1,j,l))*pfx(i,j,l)
            zy(i,j)=(z(i,j,l)-z(i,j-1,l))*pfy(i,j,l)
            zne(i,j)=(z(i,j,l)-z(i-1,j-1,l))*pfne(i,j,l)
            znw(i,j)=(z(i-1,j,l)-z(i,j-1,l))*pfnw(i,j,l)
          enddo
        enddo
!-----------------------------------------------------------------------
!---advection of height-------------------------------------------------
!-----------------------------------------------------------------------
        do j=jts_b1,jte_b1
          fahp=-fah(j)/dt
          do i=its_b1,ite_b1
            w(i,j,l)=((zx(i,j)+zx(i+1,j)+zy(i,j)+zy(i,j+1)) &
                    +(zne(i+1,j+1)+zne(i,j) &
                     +znw(i,j+1)+znw(i+1,j))*0.25)*fahp &
                    /(dsg2(l)*pdo(i,j)+pdsg1(l)) &
                    +w(i,j,l)
          enddo
        enddo
!-----------------------------------------------------------------------
!---setting w on h points to 0. along boundaries------------------------
!-----------------------------------------------------------------------
        if(.not.global) then
!-----------------------------------------------------------------------
          if(s_bdy)then
            do j=jts,jts+1
              do i=its,ite
                w(i,j,l)=0.
              enddo
            enddo
          endif
!
          if(n_bdy)then
            do j=jte-1,jte
              do i=its,ite
                w(i,j,l)=0.
              enddo
            enddo
          endif
!
          if(w_bdy)then
            do j=jts,jte
              do i=its,its+1
                w(i,j,l)=0.
              enddo
            enddo
          endif
!
          if(e_bdy)then
            do j=jts,jte
              do i=ite-1,ite
                w(i,j,l)=0.
              enddo
            enddo
          endif
!-----------------------------------------------------------------------
        endif
!-----------------------------------------------------------------------
      enddo vertical_loop
!------------------
#ifdef ENABLE_SMP
!------------------
!.......................................................................
!$omp end parallel do
!.......................................................................
!------------------
#endif
!------------------
!-----------------------------------------------------------------------
!---taking external mode out--------------------------------------------
!-----------------------------------------------------------------------
      do j=jts_b1,jte_b1
        do i=its_b1,ite_b1
!         ttb(i,j)=0.
          baro(i,j)=0.
        enddo
      enddo
!
      do l=1,lm
        do j=jts_b1,jte_b1
          do i=its_b1,ite_b1
!           ttb(i,j)=(dsg2(l)*pdo(i,j)+pdsg1(l))*w(i,j,l)+ttb(i,j)
            baro(i,j)=(dsg2(l)*pdo(i,j)+pdsg1(l))*w(i,j,l)         &
                      +baro(i,j)
          enddo
        enddo
      enddo
!
      do j=jts_b1,jte_b1
        do i=its_b1,ite_b1
!         ttb(i,j)=ttb(i,j)/pdo(i,j)
          baro(i,j)=baro(i,j)/pdo(i,j)
        enddo
      enddo
!
      do l=1,lm
        do j=jts_b1,jte_b1
          do i=its_b1,ite_b1
!           w(i,j,l)=w(i,j,l)-ttb(i,j)
            w(i,j,l)=w(i,j,l)-baro(i,j)
          enddo
        enddo
      enddo
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
                        endsubroutine cdzdt
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
                        subroutine cdwdt &
(global,hydro,restart &
,inpes,jnpes,lm,ntsd &
,dt,g &
,dsg2,pdsg1,psgml1 &
,fah &
,hdacx,hdacy &
,pd,pdo &
,psgdt &
,dwdt,pdwdt,w &
,pint &
!---temporary arguments-------------------------------------------------
,def,pfx,pfy,pfne,pfnw)
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
implicit none
!-----------------------------------------------------------------------
integer(kind=kint),parameter:: &
 nsmud=00 &                  ! number of smoothing iterations
,lnsnh=02                    ! # of rows with smoothing along boundaries

real(kind=kfpt),parameter:: &
 epsfc=9.80 &                ! limiter value
,epsn=-epsfc &               ! floor value for vertical acceleration
,epsp=epsfc &                ! upper limit for vertical acceleration
,fwhy=1. &                   ! dwdt control factor
,slpd=1500. &                ! dwdt control layer depth (Pa)
,epsvw=9999. &               ! limit on horizontal advection of w
,wa=0.125 &                  ! weighting factor
,wb=0.5 &                    ! weighting factor
!,wad=0.125 &                 ! lateral smoothing weight 
!,wad=0.0625 &                ! lateral smoothing weight
!,wad=0.075 &                 ! lateral smoothing weight
!,wad=0.050 &                 ! lateral smoothing weight
,wad=0. &                    ! lateral smoothing weight
!,wp=0.075 &                  ! time smoothing weight
,wp=0.                       ! time smoothing weight
!-----------------------------------------------------------------------

integer(kind=kint):: &
 lsltp                       ! # of layers within dwdt contol range

logical(kind=klog),intent(in):: &
 global &                    ! global or regional
,hydro  &                    ! hydrostatic or nonhydrostatic
,restart                     ! restart

integer(kind=kint),intent(in):: &
 lm &                        ! total # of levels
,ntsd &                      ! time step
,inpes &                     ! tasks in x direction
,jnpes                       ! tasks in y direction

real(kind=kfpt),intent(in):: &
 dt &                        ! dynamics time step
,g                           ! gravity

real(kind=kfpt),dimension(1:lm),intent(in):: &
 dsg2 &                      ! delta sigmas
,pdsg1 &                     ! delta pressures
,psgml1                      ! midlayer pressure part

real(kind=kfpt),dimension(jds:jde),intent(in):: &
 fah                         ! delta sigmas

real(kind=kfpt),dimension(ims:ime,jms:jme),intent(in):: &
 pd &                        ! sigma range pressure difference
,hdacx &                     ! exchange coefficient mass point
,hdacy &                     ! exchange coefficient mass point
,pdo                         ! old sigma range pressure difference

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm-1),intent(in):: &
 psgdt                       ! vertical mass flux

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(inout):: &
 dwdt &                      ! nonhydrostatic correction factor
,pdwdt &                     ! previous nonhydrostatic correction factor
,w                           ! w wind component

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm+1),intent(inout):: &
 pint                        ! pressure
!-----------------------------------------------------------------------
!---temporary arguments-------------------------------------------------
!-----------------------------------------------------------------------
real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(in) :: &
 def &                       ! deformation
,pfx &                       ! mass flux, x direction
,pfy &                       ! mass flux, y direction
,pfne &                      ! mass flux, ne direction
,pfnw                        ! mass flux, nw direction
!-----------------------------------------------------------------------
!--local variables------------------------------------------------------
!-----------------------------------------------------------------------
logical(kind=klog) :: diffw      ! turn horizontal diffusion of w on/off

integer(kind=kint):: &
 i &                         ! index in x direction
,imn &                       !
,imx &                       !
,j &                         ! index in y direction
,jmn &                       !
,jmx &                       !
,l &                         ! index in p direction
,lmn &                       !
,lmx &                       !
,kn &                        ! counter
,kp &                        ! counter
,ks                          ! smoothing counter

real(kind=kfpt):: &
 advec &                     !
,arg &                       !
,dwdtmn &                    ! minimum value of dwdt
,dwdtmx &                    ! maximum value of dwdt
,dwdtp &                     ! nonhydrostatic correction factor at the point
,fahp &                      ! grid factor
,rdt &                       ! 1/dt
,rg &                        ! 1/g
,sltp                        ! uppermost midlayer pressure

real(kind=kfpt),dimension(1:lm):: &
 why                         ! dwdt control weight

real(kind=kfpt),dimension(its_b1:ite_b1,jts_b1:jte_b1):: &
 tta                         ! advection through upper interface
real(kind=kfpt),dimension(its:ite,jts:jte):: &
 ttb                         ! advection through lower interface
real(kind=kfpt),dimension(its_b1:ite_h1,jts_b1:jte_h1):: &
 wne &                       ! height flux, ne direction
,wnw &                       ! height flux, nw direction
,wx &                        ! height flux, x direction
,wy                          ! height flux, y direction

real(kind=kfpt),dimension(ims:ime,jms:jme):: &
 ww                          ! temporary for lateral smoothing
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
!--- cannot diffuse w when in the backward step of digital filtering 
!
      if (dt .gt. 0) then
        diffw=.true.
      else
        diffw=.false.
      endif
!
!-----------------------------------------------------------------------
      if(hydro.or.(.not.hydro.and..not.restart.and.ntsd.lt.2)) then
!-----------------------------------------------------------------------
        do l=1,lm
          do j=jts,jte
            do i=its,ite
              dwdt(i,j,l)=1.
              pdwdt(i,j,l)=1.
              pint(i,j,l+1)=dsg2(l)*pd(i,j)+pdsg1(l)+pint(i,j,l)
            enddo
          enddo
        enddo
!-----------------------------------------------------------------------
        return
!-----------------------------------------------------------------------
      endif
!-----------------------------------------------------------------------
      if(.not.global) then 
!---smoothing of w on h points along boundaries-------------------------
        if(nsmud.gt.0) then
!-----------------------------------------------------------------------
          do ks=1,nsmud
!-----------------------------------------------------------------------
            do l=1,lm
!-----------------------------------------------------------------------
              if(s_bdy)then
                do j=jts+1,lnsnh
                  do i=its_b1,ite_b1
                    ww(i,j)=(w(i,j-1,l)+w(i-1,j,l) &
                            +w(i+1,j,l)+w(i,j+1,l))*wa &
                            +w(i,j,l)*wb       
                  enddo
                enddo
              endif
!
              if(n_bdy)then
                do j=jte-lnsnh+1,jte-1
                  do i=its_b1,ite_b1
                    ww(i,j)=(w(i,j-1,l)+w(i-1,j,l) &
                            +w(i+1,j,l)+w(i,j+1,l))*wa &
                            +w(i,j,l)*wb       
                  enddo
                enddo
              endif
!
              if(w_bdy)then
                do j=max(jts,jds+lnsnh),min(jte,jde-lnsnh)
                  do i=its+1,its-1+lnsnh
                    ww(i,j)=(w(i,j-1,l)+w(i-1,j,l) &
                            +w(i+1,j,l)+w(i,j+1,l))*wa &
                            +w(i,j,l)*wb       
                  enddo
                enddo
              endif
!
              if(e_bdy)then
                do j=max(jts,jds+lnsnh),min(jte,jde-lnsnh)
                  do i=ite-lnsnh+1,ite-1
                    ww(i,j)=(w(i,j-1,l)+w(i-1,j,l) &
                            +w(i+1,j,l)+w(i,j+1,l))*wa &
                            +w(i,j,l)*wb       
                  enddo
                enddo
              endif
!
              if(s_bdy)then
                do j=its+1,its-1+lnsnh
                  do i=its_b1,ite_b1
                    w(i,j,l)=ww(i,j)
                  enddo
                enddo
              endif
!
              if(n_bdy)then
                do j=jte-lnsnh+1,jte-1
                  do i=its_b1,ite_b1
                    w(i,j,l)=ww(i,j)
                  enddo
                enddo
              endif
!
              if(w_bdy)then
                do j=max(jts,jds+lnsnh),min(jte,jde-lnsnh)
                  do i=its+1,lnsnh
                    w(i,j,l)=ww(i,j)
                  enddo
                enddo
              endif
!
              if(e_bdy)then
                do j=max(jts,jds+lnsnh),min(jte,jde-lnsnh)
                  do i=ite-lnsnh+1,ite-1
                    w(i,j,l)=ww(i,j)
                  enddo
                enddo
              endif
!-----------------------------------------------------------------------
            enddo
!-----------------------------------------------------------------------
          enddo
!-----------------------------------------------------------------------
        endif ! end of smoothing of w along lateral boundaries
!-----------------------------------------------------------------------
!
      endif
!
!-----------------------------------------------------------------------
!---local derivative of w on h points-----------------------------------
!-----------------------------------------------------------------------
!
      rdt=1./dt
!
!.......................................................................
!$omp parallel do private(i,j,l)
!.......................................................................
      do l=1,lm
        do j=jts,jte
          do i=its,ite
            dwdt(i,j,l)=(w(i,j,l)-dwdt(i,j,l))*rdt
          enddo
        enddo
      enddo
!.......................................................................
!$omp end parallel do 
!.......................................................................

!-----------------------------------------------------------------------
!---lateral diffusion of w----------------------------------------------
!-----------------------------------------------------------------------
!
      if(diffw) then
!
!-----------------------------------------------------------------------
!
        do l=1,lm
!
!-----------------------------------------------------------------------
!---w fluxes, on h points-----------------------------------------------
!-----------------------------------------------------------------------
!
          do j=jts_b1,jte_h1
            do i=its_b1,ite_h1
              wx(i,j)=(w(i,j,l)-w(i-1,j,l))*(def(i-1,j,l)+def(i,j,l))
              wy(i,j)=(w(i,j,l)-w(i,j-1,l))*(def(i,j-1,l)+def(i,j,l))
            enddo
          enddo
!
!-----------------------------------------------------------------------
!---diffusion of w------------------------------------------------------
!-----------------------------------------------------------------------
!
          do j=jts_b1,jte_b1
            do i=its_b1,ite_b1
              dwdt(i,j,l)=-((wx(i+1,j)-wx(i,j))*hdacx(i,j) &
                           +(wy(i,j+1)-wy(i,j))*hdacy(i,j))*rdt &
                         +dwdt(i,j,l)
            enddo
          enddo
!
!-----------------------------------------------------------------------
!
        enddo
!
!-----------------------------------------------------------------------
!
      endif
!-----------------------------------------------------------------------
!---vertical advection of w---------------------------------------------
!-----------------------------------------------------------------------
      do j=jts,jte
        do i=its,ite
          ttb(i,j)=0.
        enddo
      enddo
!
      do l=1,lm-1
        do j=jts_b1,jte_b1
          do i=its_b1,ite_b1
            tta(i,j)=(w(i,j,l+1)-w(i,j,l))*psgdt(i,j,l)*0.5
            dwdt(i,j,l)=(tta(i,j)+ttb(i,j)) &
                       /(dsg2(l)*pdo(i,j)+pdsg1(l)) &
                       +dwdt(i,j,l)
            ttb(i,j)=tta(i,j)
          enddo
        enddo
      enddo
!
      do j=jts,jte
        do i=its,ite
          dwdt(i,j,lm)=ttb(i,j)/(dsg2(lm)*pdo(i,j)+pdsg1(lm)) &
                      +dwdt(i,j,lm)
        enddo
      enddo

      kn=0
      kp=0
      imn=0
      jmn=0
      lmn=0
      imx=0
      jmx=0
      lmx=0
      dwdtmx=0.
      dwdtmn=0.

!-----------------------------------------------------------------------
!
!.......................................................................
!$omp parallel do private (fahp,i,j,l,wne,wnw,wx,wy)
!.......................................................................
!-----------------------------------------------------------------------
      do l=1,lm
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!---w fluxes, on h points-----------------------------------------------
!-----------------------------------------------------------------------
        do j=jts_b1,jte_h1
          do i=its_b1,ite_h1
            wx(i,j)=(w(i,j,l)-w(i-1,j,l))*pfx(i,j,l)
            wy(i,j)=(w(i,j,l)-w(i,j-1,l))*pfy(i,j,l)
            wne(i,j)=(w(i,j,l)-w(i-1,j-1,l))*pfne(i,j,l)
            wnw(i,j)=(w(i-1,j,l)-w(i,j-1,l))*pfnw(i,j,l)
          enddo
        enddo
!-----------------------------------------------------------------------
!---advection of w------------------------------------------------------
!-----------------------------------------------------------------------
        do j=jms,jme
          do i=ims,ime
            ww(i,j)=0.
          enddo
        enddo

        do j=jts_b1,jte_b1
          fahp=-fah(j)/dt
          do i=its_b1,ite_b1
            advec=((wx(i,j)+wx(i+1,j)+wy(i,j)+wy(i,j+1)) &
                  +(wne(i+1,j+1)+wne(i,j) &
                   +wnw(i,j+1)+wnw(i+1,j))*0.25)*fahp &
                 /(dsg2(l)*pdo(i,j)+pdsg1(l))

            dwdtp=advec

            if(dwdtp.gt.dwdtmx) then
              dwdtmx=dwdtp
              imx=i
              jmx=j
              lmx=l
            endif
            if(dwdtp.lt.dwdtmn) then
              dwdtmn=dwdtp
              imn=i
              jmn=j
              lmn=l
            endif
            if(dwdtp.gt. epsvw) then
              kp=kp+1
            endif
            if(dwdtp.lt.-epsvw) then
              kn=kn+1
            endif

!
            if(advec.gt. epsvw) advec= epsvw
            if(advec.lt.-epsvw) advec=-epsvw
!
            dwdt(i,j,l)=advec &
                       +dwdt(i,j,l)

          enddo
        enddo
!-----------------------------------------------------------------------
!
      enddo
!.......................................................................
!$omp end parallel do 
!.......................................................................
 1300 format(' **** advecmx=',f9.5,' kp=',i6,' imx=',i4,' jmx=',i4,' lmx=',i2)
 1400 format(' **** advecmn=',f9.5,' kn=',i6,' imn=',i4,' jmn=',i4,' lmn=',i2)
        do l=1,lm
          why(l)=-99.
        enddo
!
        lsltp=0
        sltp=(psgml1(1)+psgml1(2))*0.5
        do l=1,lm-1
          arg=((psgml1(l)+psgml1(l+1))*0.5-sltp)/slpd
          if(arg.gt.1.) exit
          why(l)=1.-fwhy*cos(arg*pi*0.5)**2
          lsltp=l
        enddo
!
      do l=1,lsltp
        do j=jts_b1,jte_b1
          do i=its_b1,ite_b1
            dwdt(i,j,l)=dwdt(i,j,l)*why(l)
          enddo
        enddo
      enddo
!-----------------------------------------------------------------------
!---taking external mode out--------------------------------------------
!-----------------------------------------------------------------------
      do j=jts_b1,jte_b1
        do i=its_b1,ite_b1
          ttb(i,j)=0.
        enddo
      enddo
!
      do l=1,lm
        do j=jts_b1,jte_b1
          do i=its_b1,ite_b1
            ttb(i,j)=(dsg2(l)*pdo(i,j)+pdsg1(l))*dwdt(i,j,l)+ttb(i,j)
          enddo
        enddo
      enddo
!
      do j=jts_b1,jte_b1
        do i=its_b1,ite_b1
          ttb(i,j)=ttb(i,j)/pdo(i,j)
        enddo
      enddo
!
      do l=1,lm
        do j=jts_b1,jte_b1
          do i=its_b1,ite_b1
            dwdt(i,j,l)=dwdt(i,j,l)-ttb(i,j)
          enddo
        enddo
      enddo
!-----------------------------------------------------------------------
      if(global) then
        btim=timef()
        call swaphn(dwdt,ims,ime,jms,jme,lm,inpes)
        timers(my_domain_id)%swaphn_tim=timers(my_domain_id)%swaphn_tim+(timef()-btim)
!
        btim=timef()
        call polehn(dwdt,ims,ime,jms,jme,lm,inpes,jnpes)
        timers(my_domain_id)%polehn_tim=timers(my_domain_id)%polehn_tim+(timef()-btim)
      endif
!
      btim=timef()
      call halo_exch(dwdt,lm,1,1)
      timers(my_domain_id)%exch_dyn=timers(my_domain_id)%exch_dyn+(timef()-btim)
!
!-----------------------------------------------------------------------
!------------spatial filtering of dwdt----------------------------------
!-----------------------------------------------------------------------
!
      if(wad.gt.0.) then
!
!.......................................................................
!$omp parallel do private (i,j,l,ww)
!.......................................................................
        do l=1,lm
          do j=jts_b1,jte_b1
            do i=its_b1,ite_b1
              ww(i,j)=(dwdt(i,j-1,l)+dwdt(i-1,j,l) &
                      +dwdt(i+1,j,l)+dwdt(i,j+1,l)-dwdt(i,j,l)*4.0)
            enddo
          enddo
!
          do j=jts_b1,jte_b1
            do i=its_b1,ite_b1
              dwdt(i,j,l)=ww(i,j)*wad+dwdt(i,j,l)
            enddo
          enddo
        enddo
!.......................................................................
!$omp end parallel do 
!.......................................................................
      endif
!
!-----------------------------------------------------------------------
!
      if(global) then
        btim=timef()
        call swaphn(dwdt,ims,ime,jms,jme,lm,inpes)
        timers(my_domain_id)%swaphn_tim=timers(my_domain_id)%swaphn_tim+(timef()-btim)
!
        btim=timef()
        call polehn(dwdt,ims,ime,jms,jme,lm,inpes,jnpes)
        timers(my_domain_id)%polehn_tim=timers(my_domain_id)%polehn_tim+(timef()-btim)
      endif
!-----------------------------------------------------------------------
      rg=1./g
!
      kn=0
      kp=0
      imn=0
      jmn=0
      lmn=0
      imx=0
      jmx=0
      lmx=0
      dwdtmx=0.
      dwdtmn=0.
!
      do l=1,lm
        do j=jts,jte
          do i=its,ite
            dwdtp=dwdt(i,j,l)
            if(dwdtp.gt.dwdtmx) then
              dwdtmx=dwdtp
              imx=i
              jmx=j
              lmx=l
            endif
            if(dwdtp.lt.dwdtmn) then
              dwdtmn=dwdtp
              imn=i
              jmn=j
              lmn=l
            endif
            if(dwdtp.lt.epsn) then
              dwdtp=epsn
              kn=kn+1
            endif
            if(dwdtp.gt.epsp) then
              dwdtp=epsp
              kp=kp+1
            endif
            dwdt(i,j,l)=(dwdtp*rg+1.)*(1.-wp)+pdwdt(i,j,l)*wp
          enddo
        enddo
      enddo
 1100 format(' dwdtmx=',f9.5,' kp=',i6,' imx=',i4,' jmx=',i4,' lmx=',i2)
 1200 format(' dwdtmn=',f9.5,' kn=',i6,' imn=',i4,' jmn=',i4,' lmn=',i2)
!-----------------------------------------------------------------------
!
      if(.not.global) then 
!
!-----------------------------------------------------------------------
!---setting dwdt on h points to 1. along boundaries---------------------
!-----------------------------------------------------------------------
        do l=1,lm
          if(s_bdy)then
            do j=jts,jts+1
              do i=its,ite
                dwdt(i,j,l)=1.
              enddo
            enddo
          endif
!
          if(n_bdy)then
            do j=jte-1,jte
              do i=its,ite
                dwdt(i,j,l)=1.
              enddo
            enddo
          endif
!
          if(w_bdy)then
            do j=jts,jte
              do i=its,its+1
                dwdt(i,j,l)=1.
              enddo
            enddo
          endif
!
          if(e_bdy)then
            do j=jts,jte
              do i=ite-1,ite
                dwdt(i,j,l)=1.
              enddo
            enddo
          endif
        enddo
!-----------------------------------------------------------------------
      endif ! regional 
!-----------------------------------------------------------------------
!
                        endsubroutine cdwdt
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
                        subroutine vsound &
(global,hydro,restart &
,lm,ntsd &
,cp,dt,pt &
,dsg2,pdsg1 &
,pd &
,cw,q,rtop &
,dwdt,t,w,w_tot,baro &
,pint)
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
implicit none
!-----------------------------------------------------------------------
integer(kind=kint),parameter:: &
 nsmud=0                     ! number of smoothing iterations

real(kind=kfpt),parameter:: &
 wght=0.35                   ! first guess weight
!-----------------------------------------------------------------------
logical(kind=klog),intent(in):: &
 global &                    ! global or regional
,hydro &                     ! hydrostatic or nonhydrostatic
,restart                     ! restart case

integer(kind=kint),intent(in):: &
 lm &                        ! total # of levels
,ntsd                        ! time step

real(kind=kfpt),intent(in):: &
 cp &                        ! cp
,dt &                        ! dynamics time step
,pt                          ! pressure at the top of the model's atmosphere

real(kind=kfpt),dimension(1:lm),intent(in):: &
 dsg2 &                      ! delta sigmas
,pdsg1                       ! delta pressures

real(kind=kfpt),dimension(ims:ime,jms:jme),intent(in):: &
 baro 

real(kind=kfpt),dimension(ims:ime,jms:jme),intent(in):: &
 pd                          ! sigma range pressure difference

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(in):: &
 cw &                        ! condensate
,q &                         ! specific humidity
,rtop                        ! rt/p

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(inout):: &
 dwdt &                      ! nonhydrostatic correction factor
,t &                         ! previous nonhydrostatic correction factor
,w                           ! w wind component

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(out):: &
 w_tot                       ! total w wind component for output

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm+1),intent(inout):: &
 pint
!-----------------------------------------------------------------------
!--local variables------------------------------------------------------
!-----------------------------------------------------------------------
integer(kind=kint):: &
 i &                         ! index in x direction
,j &                         ! index in y direction
,l                           ! index in p direction

real(kind=kfpt):: &
 cappa &                     ! R/cp
,cofl &                      !
,dp &                        !
,delp &                      !
,dppl &                      !
,dpstr &                     !
,dptl &                      !
,fcc &                       ! 
,ffc &                       ! 
,gdt &                       ! g*dt
,gdt2 &                      ! gdt**2
,pp1 &                       !
,rcph &                      ! 
,rdpdn &                     !
,rdpup &                     !
,tfc &                       !
,tmp                         !

real(kind=kfpt),dimension(its_b1:ite_b1,jts_b1:jte_b1):: &
 dptu                        !

real(kind=kfpt),dimension(its_b1:ite_b1,jts_b1:jte_b1,1:lm):: &
 b1 &                        !
,b2 &                        !
,b3 &                        !
,c0 &                        !
,rdpp                        !

real(kind=kfpt),dimension(its_b1:ite_b1,jts_b1:jte_b1,1:lm+1):: &
 chi &                       !
,coff &                      !
,dfrh &                      !
,pnp1 &                      !
,pone &                      !
,pstr                        !
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      if(hydro.or.ntsd.lt.2) return
!
!-----------------------------------------------------------------------
      cappa=r/cp
      gdt=g*dt
      gdt2=gdt*gdt
      ffc=-r*0.25/gdt2
      rcph=0.5/cp
!-----------------------------------------------------------------------
!-----------------
#ifdef ENABLE_SMP
!-----------------
!.......................................................................
!$omp parallel &
!$omp private(i,j,jstart,jstop,nth,tid)
!.......................................................................
         nth = omp_get_num_threads()
         tid = omp_get_thread_num()
         call looplimits(tid, nth, jts_b1, jte_b1, jstart, jstop)
!-----------------
#else
!-----------------
      jstart = jts_b1
      jstop  = jte_b1
!-----------------
#endif
!-----------------
!
      do j=jstart,jstop
        do i=its_b1,ite_b1
          pone(i,j,1)=pt
          pstr(i,j,1)=pt
          pnp1(i,j,1)=pt
        enddo
      enddo
!
!-----------------
#ifdef ENABLE_SMP
!-----------------
!.......................................................................
!$omp end parallel
!.......................................................................
!-----------------
#endif
!-----------------
!
      jstart = jts_b1
      jstop  = jte_b1
!
      do l=2,lm+1
        do j=jstart,jstop
          do i=its_b1,ite_b1
            dppl=dsg2(l-1)*pd(i,j)+pdsg1(l-1)
            rdpp(i,j,l-1)=1./dppl
            dpstr=dwdt(i,j,l-1)*dppl
            pstr(i,j,l)=pstr(i,j,l-1)+dpstr
            pp1=pnp1(i,j,l-1)+dpstr
            pone(i,j,l)=pint(i,j,l)
            dp=(pp1-pone(i,j,l))*wght
            pnp1(i,j,l)=pone(i,j,l)+dp
            tfc=q(i,j,l-1)*0.608-cw(i,j,l-1)+1.
            fcc=(1.-cappa*tfc)*tfc*ffc
            cofl=t(i,j,l-1)*dppl*fcc &
                /((pnp1(i,j,l-1)+pnp1(i,j,l))*0.5)**2
            coff(i,j,l-1)=cofl
            dfrh(i,j,l)=(pstr(i,j,l-1)+pstr(i,j,l) &
                        -pone(i,j,l-1)-pone(i,j,l))*cofl
          enddo
        enddo
      enddo
!-----------------------------------------------------------------------
!-----------------
#ifdef ENABLE_SMP
!-----------------
!.......................................................................
!$omp parallel &
!$omp private(delp,dptl,i,j,jstart,jstop,l,nth,rdpdn,rdpup,tid,tmp)
!.......................................................................
         nth = omp_get_num_threads()
         tid = omp_get_thread_num()
         call looplimits(tid, nth, jts_b1, jte_b1, jstart, jstop)
!-----------------
#else
!-----------------
      jstart = jts_b1
      jstop  = jte_b1
!-----------------
#endif
!-----------------
! 
      do l=2,lm
        do j=jstart,jstop
          do i=its_b1,ite_b1
            rdpdn=rdpp(i,j,l)
            rdpup=rdpp(i,j,l-1)
            b1(i,j,l)=coff(i,j,l-1)+rdpup
            b2(i,j,l)=coff(i,j,l-1)+coff(i,j,l)-rdpup-rdpdn
            b3(i,j,l)=coff(i,j,l)+rdpdn
            c0(i,j,l)=-dfrh(i,j,l)-dfrh(i,j,l+1)
          enddo
        enddo
      enddo
!-----------------------------------------------------------------------
       do j=jstart,jstop
        do i=its_b1,ite_b1
          b2(i,j,lm)=b2(i,j,lm)+b3(i,j,lm)
        enddo
      enddo
!
      do l=3,lm
        do j=jstart,jstop
          do i=its_b1,ite_b1
            tmp=-b1(i,j,l)/b2(i,j,l-1)
            b2(i,j,l)=b3(i,j,l-1)*tmp+b2(i,j,l)
            c0(i,j,l)=c0(i,j,l-1)*tmp+c0(i,j,l)
          enddo
        enddo
      enddo
!-----------------------------------------------------------------------
       do j=jstart,jstop
        do i=its_b1,ite_b1
          chi(i,j,1)=0.
          chi(i,j,lm)=c0(i,j,lm)/b2(i,j,lm)
          chi(i,j,lm+1)=chi(i,j,lm)
        enddo
      enddo
!
      do l=lm-1,2,-1
        do j=jstart,jstop
          do i=its_b1,ite_b1
            chi(i,j,l)=(-b3(i,j,l)*chi(i,j,l+1)+c0(i,j,l))/b2(i,j,l)
          enddo
        enddo
      enddo
!-----------------------------------------------------------------------
      do l=1,lm+1
       do j=jstart,jstop
           do i=its_b1,ite_b1
            pnp1(i,j,l)=chi(i,j,l)+pstr(i,j,l)
            pint(i,j,l)=pnp1(i,j,l)
          enddo
        enddo
      enddo
!-----------------------------------------------------------------------
      do j=jstart,jstop
        do i=its_b1,ite_b1
          dptu(i,j)=0.
        enddo
      enddo
!
      do l=1,lm
        do j=jstart,jstop
          do i=its_b1,ite_b1
            dptl=pnp1(i,j,l+1)-pone(i,j,l+1)
            t(i,j,l)=(dptu(i,j)+dptl)*rtop(i,j,l)*rcph+t(i,j,l)
            delp=(pnp1(i,j,l+1)-pnp1(i,j,l))*rdpp(i,j,l)
            w(i,j,l)=(delp-dwdt(i,j,l))*gdt+w(i,j,l)
            w_tot(i,j,l)=w(i,j,l)+baro(i,j)
            dwdt(i,j,l)=delp
            dptu(i,j)=dptl
          enddo
        enddo
      enddo
!-----------------
#ifdef ENABLE_SMP
!-----------------
!.......................................................................
!$omp end parallel
!.......................................................................
!-----------------
#endif
!-----------------
!-----------------------------------------------------------------------
!
      if(.not.global) then
!
!-----------------------------------------------------------------------
!---setting dwdt on h points to 1. along boundaries---------------------
!-----------------------------------------------------------------------
        if(s_bdy)then
          do l=1,lm
            do j=jts,jts+1
              do i=its,ite
                dwdt(i,j,l)=1.
              enddo
            enddo
          enddo
        endif
!
        if(n_bdy)then
          do l=1,lm
            do j=jte-1,jte
              do i=its,ite
                dwdt(i,j,l)=1.
              enddo
            enddo
          enddo
        endif
!
        if(w_bdy)then
          do l=1,lm
            do j=jts,jte
              do i=its,its+1
                dwdt(i,j,l)=1.
              enddo
            enddo
          enddo
        endif
!
        if(e_bdy)then
          do l=1,lm
            do j=jts,jte
              do i=ite-1,ite
                dwdt(i,j,l)=1.
              enddo
            enddo
          enddo
        endif
!-----------------------------------------------------------------------
      endif ! regional
!-----------------------------------------------------------------------
!
                        endsubroutine vsound
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        subroutine adv2 &
(global &
,idtadt,kss,kse,lm,lnsad &
,dt,rdyh &
,dsg2,pdsg1 &
,epsq2 &
,fah,rdxh &
,pd,pdo &
,psgdt &
,up,vp &
,indx_q2 &
,s,sp &
!---temporary arguments-------------------------------------------------
,pfne,pfnw,pfx,pfy,s1,tcs)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

implicit none

!-----------------------------------------------------------------------
real(kind=kfpt),parameter:: &
 cfc=1.533 &                 ! adams-bashforth positioning in time
,bfc=1.-cfc &                ! adams bashforth positioning in time
,epsq=1.e-20 &               ! floor value for specific humidity
,pfc=1.+4./6. &              ! 4th order momentum advection
,sfc=-1./6. &                ! 4th order momentum advection
,epscm=2.e-6 &               ! a floor value (not used)
,w1=0.9 &                    ! crank-nicholson uncentering
!,w1=1.0 &                    ! crank-nicholson uncentering
!,w1=0.80 &                   ! crank-nicholson uncentering
,w2=2.-w1                    ! crank-nicholson uncentering

logical(kind=klog),intent(in):: &
 global

integer(kind=kint),intent(in):: &
 idtadt &                    !
,kse &                       ! terminal species index
,kss &                       ! initial species index
,lm &                        ! total # of levels
,lnsad &                     !
,indx_q2                     ! location of q2 in tracer arrays

real(kind=kfpt),intent(in):: &
 dt &                        ! dynamics time step
,rdyh                        !

real(kind=kfpt),dimension(1:lm),intent(in):: &
 dsg2 &                      ! delta sigmas
,pdsg1 &                     ! delta pressures
,epsq2                       ! floor value of 2tke

real(kind=kfpt),dimension(jds:jde),intent(in):: &
 fah &                       ! grid factor
,rdxh                        !

real(kind=kfpt),dimension(ims:ime,jms:jme),intent(in):: &
 pd &                        ! sigma range pressure difference
,pdo                         ! sigma range pressure difference

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm-1),intent(in):: &
 psgdt                       ! vertical mass flux

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(in):: &
 up &                        !
,vp                          !

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm,1:kse),intent(inout):: &
 s &                         ! tracers
,sp                          ! s at previous time level

!---temporary arguments-------------------------------------------------
real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm):: &
 pfne &                      ! mass flux, ne direction
,pfnw &                      ! mass flux, nw direction
,pfx &                       ! mass flux, x direction
,pfy                         ! mass flux, y direction

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm,1:kse),intent(inout):: &
 s1 &                        ! intermediate value of sqrt(s)
,tcs                         ! timechange of s

!--local variables------------------------------------------------------
real(kind=kfpt),parameter:: &
 wb=1.0 &                    ! weighting factor
,wa=(1.0-wb)*0.25            ! weighting factor

integer(kind=kint):: &
 i &                         !
,iap &                       !
,ibeg &                      !
,iend &                      !
,j &                         !
,jap &                       !
,jbeg &                      !
,jend &                      !
,ks &                        !
,l                           !

real(kind=kfpt):: &
 cf &                        ! temporary
,cms &                       ! temporary
,dtq &                       ! dt/4
,emhp &                      !
,enhp &                      !
,fahp &                      ! temporary grid factor
,pp &                        !
,qq &                        !
,rdp &                       ! 1/deltap
,vvlo &                      ! vertical velocity, lower interface
,vvup &                      ! vertical velocity, upper interface
,pvvup                       ! vertical mass flux, upper interface

real(kind=kfpt),dimension(ims:ime,jms:jme):: &
 pdop &                      ! hydrostatic pressure difference at v points
,pdops &                     ! smoothed hydrostatic pressure difference at h points
,pvvlo &                     ! vertical mass flux, lower interface
,ss1 &                       ! extrapolated species between time levels
,ssne &                      ! flux, ne direction
,ssnw &                      ! flux, nw direction
,ssx &                       ! flux, x direction
,ssy                         ! flux, y direction

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm):: &
 crs &                       ! vertical advection temporary
,rcms                        ! vertical advection temporary

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm,1:kse):: &
 rsts                        ! vertical advection temporary
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------
#ifdef ENABLE_SMP
!-----------------
!.......................................................................
!$omp parallel private(i,j,jstart,jstop,ks,l,nth,tid)
!.......................................................................
      nth = omp_get_num_threads()
      tid = omp_get_thread_num()
      call looplimits(tid, nth, jts_h1, jte_h1, jstart, jstop)
!-----------------
#else
!-----------------
      jstart = jts_h1
      jstop  = jte_h1
!-----------------
#endif
!-----------------
      do ks=kss,kse
        do l=1,lm
          do j=jstart,jstop
            do i=its_h1,ite_h1
              s(i,j,l,ks)=max(s(i,j,l,ks),epsq)
            enddo
          enddo
        enddo
      enddo
!
!***  Interpolate q2 (2*TKE) from interfaces to midlayers
!
q2_check: if (kss<=indx_q2 .and. indx_q2<=kse) then
      do l=lm,2,-1
        do j=jstart,jstop
          do i=its_h1,ite_h1
            s(i,j,L,indx_q2)=max((s(i,j,L,indx_q2)+s(i,j,L-1,indx_q2))*0.5  &
                                ,(epsq2(L)+epsq2(L-1))*0.5)
          enddo
        enddo
      enddo
      do j=jstart,jstop
        do i=its_h1,ite_h1
          s(i,j,1,indx_q2)=max((s(i,j,1,indx_q2)+epsq2(1))*0.5,epsq2(1))
        enddo
      enddo
    endif  q2_check
!-----------------------------------------------------------------------
      do ks=kss,kse ! loop by species
!-----------------------------------------------------------------------
        do l=1,lm
          do j=jstart,jstop
            do i=its_h1,ite_h1
              s1(i,j,l,ks)=sqrt(s(i,j,l,ks))
            enddo
          enddo
        enddo
!-----------------------------------------------------------------------
!
      enddo ! end of the loop by species
!
!-----------------
#ifdef ENABLE_SMP
!-----------------
!.......................................................................
!$omp end parallel
!.......................................................................
!-----------------
#endif
!-----------------
!-----------------------------------------------------------------------
      do j=jts_h1,jte_h1
        do i=its_h1,ite_h1
          pdop(i,j)=(pd(i,j)+pdo(i,j))*0.5
        enddo
      enddo
!-----------------
#ifdef ENABLE_SMP
!-----------------
!.......................................................................
!$omp parallel private(cf,cms,dtq,i,j,jstart,jstop,l,nth,pvvup,rdp, &
!$omp                  tid,vvlo,vvup)
!.......................................................................
      nth = omp_get_num_threads()
      tid = omp_get_thread_num()
      call looplimits(tid, nth, jts_b1, jte_b1, jstart, jstop)
!-----------------
#else
!-----------------
      jstart = jts_b1
      jstop  = jte_b1
!-----------------
#endif
!-----------------
!-----------------------------------------------------------------------
!---crank-nicholson vertical advection----------------------------------
!-----------------------------------------------------------------------
      dtq=dt*0.25*idtadt
      do j=jstart,jstop
        do i=its_b1,ite_b1
          pdops(i,j)=(pdop(i,j-1)+pdop(i-1,j) &
                    +pdop(i+1,j)+pdop(i,j+1))*wa &
                    +pdop(i,j)*wb
        enddo
      enddo
      do j=jstart,jstop
        do i=its_b1,ite_b1
          pvvlo(i,j)=((psgdt(i,j-1,1)+psgdt(i-1,j,1) &
                      +psgdt(i+1,j,1)+psgdt(i,j+1,1))*wa + &
                      psgdt(i,j,1)*wb)*dtq
          vvlo=pvvlo(i,j)/(dsg2(1)*pdops(i,j)+pdsg1(1))
          cms=-vvlo*w2+1.
          rcms(i,j,1)=1./cms
          crs(i,j,1)=vvlo*w2
!
          do ks=kss,kse
            rsts(i,j,1,ks)=(-vvlo*w1) &
                          *(s1(i,j,2,ks)-s1(i,j,1,ks)) &
                          +s1(i,j,1,ks)
          enddo
        enddo
      enddo
      do l=2,lm-1
        do j=jstart,jstop
          do i=its_b1,ite_b1
            rdp=1./(dsg2(l)*pdops(i,j)+pdsg1(l))
            pvvup=pvvlo(i,j)
            pvvlo(i,j)=((psgdt(i,j-1,l)+psgdt(i-1,j,l) &
                        +psgdt(i+1,j,l)+psgdt(i,j+1,l))*wa + &
                        psgdt(i,j,l)*wb)*dtq
!
            vvup=pvvup*rdp
            vvlo=pvvlo(i,j)*rdp
!
            cf=-vvup*w2*rcms(i,j,l-1)
            cms=-crs(i,j,l-1)*cf+((vvup-vvlo)*w2+1.)
            rcms(i,j,l)=1./cms
            crs(i,j,l)=vvlo*w2
!
            do ks=kss,kse
              rsts(i,j,l,ks)=-rsts(i,j,l-1,ks)*cf+s1(i,j,l,ks) &
                             -(s1(i,j,l  ,ks)-s1(i,j,l-1,ks))*vvup*w1 &
                             -(s1(i,j,l+1,ks)-s1(i,j,l  ,ks))*vvlo*w1
            enddo
          enddo
        enddo
      enddo
      do j=jstart,jstop
        do i=its_b1,ite_b1
          pvvup=pvvlo(i,j)
          vvup=pvvup/(dsg2(lm)*pdops(i,j)+pdsg1(lm))
!
          cf=-vvup*w2*rcms(i,j,lm-1)
          cms=-crs(i,j,lm-1)*cf+(vvup*w2+1.)
          rcms(i,j,lm)=1./cms
          crs(i,j,lm)=0.
!
          do ks=kss,kse
            rsts(i,j,lm,ks)=-rsts(i,j,lm-1,ks)*cf+s1(i,j,lm,ks) &
                           -(s1(i,j,lm,ks)-s1(i,j,lm-1,ks))*vvup*w1
!
            tcs(i,j,lm,ks)=rsts(i,j,lm,ks)*rcms(i,j,lm)-s1(i,j,lm,ks)
          enddo
        enddo
      enddo
      do ks=kss,kse
        do l=lm-1,1,-1
          do j=jstart,jstop
            do i=its_b1,ite_b1
              tcs(i,j,l,ks)=(-crs(i,j,l)*(s1(i,j,l+1,ks)+tcs(i,j,l+1,ks)) &
                             +rsts(i,j,l,ks)) &
                           *rcms(i,j,l)-s1(i,j,l,ks)
            enddo
          enddo
        enddo
      enddo
!-----------------
#ifdef ENABLE_SMP
!-----------------
!.......................................................................
!$omp end parallel
!.......................................................................
!-----------------
#endif
!-----------------
!-----------------------------------------------------------------------
!.......................................................................
!$omp parallel do &
!$omp private(fahp,i,ibeg,iend,j,jbeg,jend,ks,l,ss1,ssne,ssnw,ssx,ssy)
!.......................................................................
!-----------------------------------------------------------------------
!
      do ks=kss,kse ! loop by species
!
!-----------------------------------------------------------------------
        do l=1,lm
          do j=jts_h1,jte_h1
            do i=its_h1,ite_h1
              ss1(i,j)=s1(i,j,l,ks)*cfc+sp(i,j,l,ks)*bfc
              sp(i,j,l,ks)=s1(i,j,l,ks)
            enddo
          enddo
!---temperature fluxes, on h points-------------------------------------
          do j=jts_b1,jte_h1
            do i=its_b1,ite_h1
              ssx(i,j)=(ss1(i,j)-ss1(i-1,j))*pfx(i,j,l)
              ssy(i,j)=(ss1(i,j)-ss1(i,j-1))*pfy(i,j,l)
!
              ssne(i,j)=(ss1(i,j)-ss1(i-1,j-1))*pfne(i,j,l)
              ssnw(i,j)=(ss1(i-1,j)-ss1(i,j-1))*pfnw(i,j,l)
            enddo
          enddo
!---advection of species------------------------------------------------
          if(adv_standard)then
            ibeg=max(its,ids+1+lnsad)
            iend=min(ite,ide-1-lnsad)
            jbeg=max(jts,jds+1+lnsad)
            jend=min(jte,jde-1-lnsad)
!
            do j=jbeg,jend
              fahp=fah(j)*idtadt
              do i=ibeg,iend
                tcs(i,j,l,ks)=(((ssx (i  ,j  )+ssx (i+1,j  ) &
                                +ssy (i  ,j  )+ssy (i  ,j+1)) &
                               +(ssne(i+1,j+1)+ssne(i  ,j  ) &
                                +ssnw(i  ,j+1)+ssnw(i+1,j  ))*0.25) &
                               *fahp) &
                             /(dsg2(l)*pdop(i,j)+pdsg1(l)) &
                             +tcs(i,j,l,ks)
              enddo
            enddo
          endif
        enddo
!-----------------------------------------------------------------------
!
      enddo ! end of the loop by the species
!.......................................................................
!$omp end parallel do
!.......................................................................
!
!-----------------------------------------------------------------------
!---regional branch-----------------------------------------------------
!-----------------------------------------------------------------------
      if(.not.global.and.adv_upstream) then
!-----------------------------------------------------------------------
        enhp=-dt*rdyh*0.25*idtadt
!-----------------------------------------------------------------------
        do l=1,lm
!-----------------------------------------------------------------------
!
!***  Upstream advection along southern rows
!
          do j=jts_b1,min(jte,jds+lnsad)
            emhp=-dt*rdxh(j)*0.25*idtadt
            do i=its_b1,ite_b1
              pp=(up(i-1,j-1,l)+up(i  ,j-1,l) &
                 +up(i-1,j  ,l)+up(i  ,j  ,l))*emhp
              qq=(vp(i-1,j-1,l)+vp(i  ,j-1,l) &
                 +vp(i-1,j  ,l)+vp(i  ,j  ,l))*enhp
!
              if(pp.le.0.) then
                iap=-1
                pp=-pp
              else
                iap=1
              endif
!
              if(qq.le.0.) then
                jap=-1
                qq=-qq
              else
                jap=1
              endif
!
              do ks=kss,kse
                tcs(i,j,l,ks)=(s1(i+iap,j,l,ks)-s1(i,j,l,ks))*pp &
                             +(s1(i,j+jap,l,ks)-s1(i,j,l,ks))*qq &
                             +(s1(i,j,l,ks)-s1(i+iap,j,l,ks) &
                              -s1(i,j+jap,l,ks)+s1(i+iap,j+jap,l,ks)) &
                             *pp*qq &
                             +tcs(i,j,l,ks)
              enddo
            enddo
          enddo
!
!***  Upstream advection along northern rows
!
          do j=max(jts,jde-lnsad),jte_b1
            emhp=-dt*rdxh(j)*0.25*idtadt
            do i=its_b1,ite_b1
              pp=(up(i-1,j-1,l)+up(i  ,j-1,l) &
                 +up(i-1,j  ,l)+up(i  ,j  ,l))*emhp
              qq=(vp(i-1,j-1,l)+vp(i  ,j-1,l) &
                 +vp(i-1,j  ,l)+vp(i  ,j  ,l))*enhp
!
              if(pp.le.0.) then
                iap=-1
                pp=-pp
              else
                iap=1
              endif
!
              if(qq.le.0.) then
                jap=-1
                qq=-qq
              else
                jap=1
              endif
!
              do ks=kss,kse
                tcs(i,j,l,ks)=(s1(i+iap,j,l,ks)-s1(i,j,l,ks))*pp &
                             +(s1(i,j+jap,l,ks)-s1(i,j,l,ks))*qq &
                             +(s1(i,j,l,ks)-s1(i+iap,j,l,ks) &
                              -s1(i,j+jap,l,ks)+s1(i+iap,j+jap,l,ks)) &
                             *pp*qq &
                             +tcs(i,j,l,ks)
              enddo
            enddo
          enddo
!
!***  Upstream advection along western rows
!
          do j=max(jts,jds+1+lnsad),min(jte,jde-1-lnsad)
            emhp=-dt*rdxh(j)*0.25*idtadt
            do i=its_b1,min(ite,ids+lnsad)
              pp=(up(i-1,j-1,l)+up(i  ,j-1,l) &
                 +up(i-1,j  ,l)+up(i  ,j  ,l))*emhp
              qq=(vp(i-1,j-1,l)+vp(i  ,j-1,l) &
                 +vp(i-1,j  ,l)+vp(i  ,j  ,l))*enhp
!
              if(pp.le.0.) then
                iap=-1
                pp=-pp
              else
                iap=1
              endif
!
              if(qq.le.0.) then
                jap=-1
                qq=-qq
              else
                jap=1
              endif
!
              do ks=kss,kse
                tcs(i,j,l,ks)=(s1(i+iap,j,l,ks)-s1(i,j,l,ks))*pp &
                             +(s1(i,j+jap,l,ks)-s1(i,j,l,ks))*qq &
                             +(s1(i,j,l,ks)-s1(i+iap,j,l,ks) &
                              -s1(i,j+jap,l,ks)+s1(i+iap,j+jap,l,ks)) &
                             *pp*qq &
                             +tcs(i,j,l,ks)
              enddo
            enddo
          enddo
!
!***  Upstream advection along eastern rows
!
          do j=max(jts,jds+1+lnsad),min(jte,jde-1-lnsad)
            emhp=-dt*rdxh(j)*0.25*idtadt
            do i=max(its,ide-lnsad),ite_b1
              pp=(up(i-1,j-1,l)+up(i  ,j-1,l) &
                 +up(i-1,j  ,l)+up(i  ,j  ,l))*emhp
              qq=(vp(i-1,j-1,l)+vp(i  ,j-1,l) &
                 +vp(i-1,j  ,l)+vp(i  ,j  ,l))*enhp
!
              if(pp.le.0.) then
                iap=-1
                pp=-pp
              else
                iap=1
              endif
!
              if(qq.le.0.) then
                jap=-1
                qq=-qq
              else
                jap=1
              endif
!
              do ks=kss,kse
                tcs(i,j,l,ks)=(s1(i+iap,j,l,ks)-s1(i,j,l,ks))*pp &
                             +(s1(i,j+jap,l,ks)-s1(i,j,l,ks))*qq &
                             +(s1(i,j,l,ks)-s1(i+iap,j,l,ks) &
                              -s1(i,j+jap,l,ks)+s1(i+iap,j+jap,l,ks)) &
                             *pp*qq &
                             +tcs(i,j,l,ks)
              enddo
            enddo
          enddo
!-----------------------------------------------------------------------
        enddo
!-----------------------------------------------------------------------
      endif ! regional lateral boundaries
!-----------------------------------------------------------------------
!
                        endsubroutine adv2
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
                        subroutine mono &
(idtadt,kss,kse,lm &
,dsg2,pdsg1 &
,epsq2 &
,dare &
,pd &
,indx_q2 &
,s &
,inpes,jnpes &
,use_allreduce &
!---temporary arguments-------------------------------------------------
,s1,tcs)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

implicit none

!-----------------------------------------------------------------------
real(kind=kfpt),parameter:: &
 epsq=1.e-20                 ! floor value for specific humidity

integer(kind=kint),intent(in):: &
 idtadt &                    !
,inpes &                     ! number of tasks in x direction
,jnpes &                     ! number of tasks in y direction
,kse &                       ! terminal species index
,kss &                       ! initial species index
,lm &                        ! total # of levels
,indx_q2                     ! location of q2 in 4-d tracer arrays

real(kind=kfpt),dimension(1:lm),intent(in):: &
 dsg2 &                      ! delta sigmas
,pdsg1 &                     ! delta pressures
,epsq2                       ! floor value of 2tke

real(kind=kfpt),dimension(jds:jde),intent(in):: &
 dare                        ! grid box area

real(kind=kfpt),dimension(ims:ime,jms:jme),intent(in):: &
 pd                          ! sigma range pressure difference

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm,1:kse),intent(inout):: &
 s                           ! s at previous time level

logical(kind=klog) :: &
 use_allreduce

!---temporary arguments-------------------------------------------------
real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm,1:kse),intent(inout):: &
 s1 &                        ! intermediate value of s
,tcs                         ! timechange of s
!--local variables------------------------------------------------------
integer(kind=kint):: &
 i &                         !
,ierr &                      !
,irecv &                     !
,j &                         !
,ks &                        !
,l &                         !
,lngth &                     !
,loc_its_b1 &                !
,loc_ite_b1 &                !
,loc_jts_b1 &                !
,loc_jte_b1 &                !
,loc_len &                   !
,loc_lngth &                 !
,n &                         !
,nn &                        !
,pe                          !

integer,dimension(mpi_status_size) :: jstat

real(kind=kfpt):: &
 s1p &                       !
,smax &                      ! local maximum
,smin &                      ! local minimum
,smaxh &                     ! horizontal local maximum
,sminh &                     ! horizontal local minimum
,smaxv &                     ! vertical local maximum
,sminv &                     ! vertical local minimum
,sn &                        !
,steep                       !

real(kind=kdbl):: &
 dsp &                       !
,rfacs &                     !
,sfacs &                     !
,sumns &                     !
,sumps                       !

real(kind=kdbl),dimension(ide*jde*kse):: &
 s1_glob                     !

real(kind=kfpt),dimension((ite_b1-its_b1+1)*(jte_b1-jts_b1+1)*kse):: &
 s1_loc                      !

real(kind=kfpt),dimension(:), allocatable :: &
 s1_pe_loc                   !

real(kind=kdbl),dimension(1:2*kse):: &
 gsums &                     ! sum of neg/pos changes all global fields
,xsums                       ! sum of neg/pos changes all global fields

real(kind=kfpt),dimension(its_b1:ite_b1,jts_b1:jte_b1,1:kse):: &
 s1l_sum

real(kind=kfpt),dimension(its_b1:ite_b1,jts_b1:jte_b1,1:lm):: &
 dvol &                      ! grid box volume
,rdvol                       ! 1./grid box volume
!-----------------------------------------------------------------------
integer(kind=kint) :: &
 istat

logical(kind=klog) :: &
 opened

logical(kind=klog),save :: &
 sum_file_is_open=.false.

character(10) :: &
 fstatus
!-----------------------------------------------------------------------
real(kind=kdbl),save :: sumdo3=0.
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      steep=1.-0.040*idtadt
!.......................................................................
!$omp parallel 
!$omp do private(i,j,l)
!.......................................................................
      do l=1,lm
        do j=jts_b1,jte_b1
          do i=its_b1,ite_b1
            dvol (i,j,l)=(dsg2(l)*pd(i,j)+pdsg1(l))*dare(j)
            rdvol(i,j,l)=1./dvol(i,j,l)
          enddo
        enddo
      enddo
!.......................................................................
!$omp end do
!$omp end parallel
!.......................................................................
!
!-----------------------------------------------------------------------
!---monotonization------------------------------------------------------
!-----------------------------------------------------------------------
!
      if(use_allreduce)then
!
!-----------------------------------------------------------------------
!.......................................................................
!$omp parallel 
!$omp do private (dsp,i,j,ks,l,s1p,smax,smaxh,smaxv,smin,sminh,sminv,sn)
!.......................................................................
!-----------------------------------------------------------------------
      do ks=kss,kse ! loop by species
          gsums(2*ks-1)=0.
          gsums(2*ks  )=0.
          xsums(2*ks-1)=0.
          xsums(2*ks  )=0.
!
        do j=jts_b1,jte_b1
          do i=its_b1,ite_b1
            s1l_sum(i,j,ks)=0.
          enddo
        enddo
!-----------------------------------------------------------------------
        do l=1,lm
          do j=jts_b1,jte_b1
            do i=its_b1,ite_b1
              s1p=(s1(i,j,l,ks)+tcs(i,j,l,ks))**2
              tcs(i,j,l,ks)=s1p-s(i,j,l,ks)
!
              sminh=min(s(i-1,j-1,l,ks) &
                       ,s(i  ,j-1,l,ks) &
                       ,s(i+1,j-1,l,ks) &
                       ,s(i-1,j  ,l,ks) &
                       ,s(i  ,j  ,l,ks) &
                       ,s(i+1,j  ,l,ks) &
                       ,s(i-1,j+1,l,ks) &
                       ,s(i  ,j+1,l,ks) &
                       ,s(i+1,j+1,l,ks))
              smaxh=max(s(i-1,j-1,l,ks) &
                       ,s(i  ,j-1,l,ks) &
                       ,s(i+1,j-1,l,ks) &
                       ,s(i-1,j  ,l,ks) &
                       ,s(i  ,j  ,l,ks) &
                       ,s(i+1,j  ,l,ks) &
                       ,s(i-1,j+1,l,ks) &
                       ,s(i  ,j+1,l,ks) &
                       ,s(i+1,j+1,l,ks))
!
              if(l.gt.1.and.l.lt.lm) then
                sminv=min(s(i,j,l-1,ks),s(i,j,l  ,ks),s(i,j,l+1,ks))
                smaxv=max(s(i,j,l-1,ks),s(i,j,l  ,ks),s(i,j,l+1,ks))
              elseif(l.eq.1) then
                sminv=min(s(i,j,l  ,ks),s(i,j,l+1,ks))
                smaxv=max(s(i,j,l  ,ks),s(i,j,l+1,ks))
              elseif(l.eq.lm) then
                sminv=min(s(i,j,l-1,ks),s(i,j,l  ,ks))
                smaxv=max(s(i,j,l-1,ks),s(i,j,l  ,ks))
              endif
!
              smin=min(sminh,sminv)
              smax=max(smaxh,smaxv)
!
              sn=s1p
              if(sn.gt.steep*smax) sn=smax
              if(sn.lt.     smin) sn=smin
!
              dsp=(sn-s1p)*dvol(i,j,l)
!
              if(dsp.gt.0.) then
                xsums(2*ks-1)=xsums(2*ks-1)+dsp
              else
                xsums(2*ks  )=xsums(2*ks  )+dsp
              endif
!
            enddo
          enddo
        enddo
      enddo ! end of the loop by species
!.......................................................................
!$omp end do
!$omp end parallel
!.......................................................................
!-----------------------------------------------------------------------
!
      else ! use send/recv
!
!-----------------------------------------------------------------------
!.......................................................................
!$omp parallel 
!$omp do private (dsp,i,j,ks,l,s1p,smax,smaxh,smaxv,smin,sminh,sminv,sn)
!.......................................................................
!-----------------------------------------------------------------------
      do ks=kss,kse ! loop by species
          gsums(2*ks-1)=0.
          gsums(2*ks  )=0.
!
        do j=jts_b1,jte_b1
          do i=its_b1,ite_b1
            s1l_sum(i,j,ks)=0.
          enddo
        enddo
!-----------------------------------------------------------------------
        do l=1,lm
          do j=jts_b1,jte_b1
            do i=its_b1,ite_b1
              s1p=(s1(i,j,l,ks)+tcs(i,j,l,ks))**2
              tcs(i,j,l,ks)=s1p-s(i,j,l,ks)
!
              sminh=min(s(i-1,j-1,l,ks) &
                       ,s(i  ,j-1,l,ks) &
                       ,s(i+1,j-1,l,ks) &
                       ,s(i-1,j  ,l,ks) &
                       ,s(i  ,j  ,l,ks) &
                       ,s(i+1,j  ,l,ks) &
                       ,s(i-1,j+1,l,ks) &
                       ,s(i  ,j+1,l,ks) &
                       ,s(i+1,j+1,l,ks))
              smaxh=max(s(i-1,j-1,l,ks) &
                       ,s(i  ,j-1,l,ks) &
                       ,s(i+1,j-1,l,ks) &
                       ,s(i-1,j  ,l,ks) &
                       ,s(i  ,j  ,l,ks) &
                       ,s(i+1,j  ,l,ks) &
                       ,s(i-1,j+1,l,ks) &
                       ,s(i  ,j+1,l,ks) &
                       ,s(i+1,j+1,l,ks))
!
              if(l.gt.1.and.l.lt.lm) then
                sminv=min(s(i,j,l-1,ks),s(i,j,l  ,ks),s(i,j,l+1,ks))
                smaxv=max(s(i,j,l-1,ks),s(i,j,l  ,ks),s(i,j,l+1,ks))
              elseif(l.eq.1) then
                sminv=min(s(i,j,l  ,ks),s(i,j,l+1,ks))
                smaxv=max(s(i,j,l  ,ks),s(i,j,l+1,ks))
              elseif(l.eq.lm) then
                sminv=min(s(i,j,l-1,ks),s(i,j,l  ,ks))
                smaxv=max(s(i,j,l-1,ks),s(i,j,l  ,ks))
              endif
!
              smin=min(sminh,sminv)
              smax=max(smaxh,smaxv)
!
              sn=s1p
              if(sn.gt.steep*smax) sn=smax
              if(sn.lt.     smin) sn=smin
!
              dsp=(sn-s1p)*dvol(i,j,l)
!
              s1(i,j,l,ks)=dsp
              s1l_sum(i,j,ks) = s1l_sum(i,j,ks) + dsp
!
            enddo
          enddo
        enddo
      enddo ! end of the loop by species
!.......................................................................
!$omp end do
!$omp end parallel
!.......................................................................
!-----------------------------------------------------------------------
!
      endif
!
!-----------------------------------------------------------------------
!***  Global reductions
!-----------------------------------------------------------------------
!
      global_reduce: if(use_allreduce)then
!
!-----------------------------------------------------------------------
!***  Skip computing the global reduction if they are to be read in
!***  from another run to check bit reproducibility.
!-----------------------------------------------------------------------
        lngth=2*kse
        call mpi_allreduce(xsums,gsums,lngth &
                          ,mpi_double_precision &
                          ,mpi_sum,mpi_comm_comp,irecv)
!-----------------------------------------------------------------------
!
      else  global_reduce
!
!-----------------------------------------------------------------------
        if (mype==0) then

          do ks=kss,kse
            do j=jts_b1,jte_b1
              do i=its_b1,ite_b1
                n=i+(j-1)*ide+(ks-1)*ide*jde
                s1_glob(n) = s1l_sum(i,j,ks)
              enddo
            enddo
          enddo

          do pe=1,inpes*jnpes-1
            loc_its_b1 = max(local_istart(pe),ids+1)
            loc_ite_b1 = min(local_iend  (pe),ide-1)
            loc_jts_b1 = max(local_jstart(pe),jds+1)
            loc_jte_b1 = min(local_jend  (pe),jde-1)
            loc_len    = (loc_ite_b1-loc_its_b1+1)*    &
                         (loc_jte_b1-loc_jts_b1+1)*kse
            allocate(s1_pe_loc(1:loc_len))
            call mpi_recv(s1_pe_loc(1:loc_len),loc_len &
                         ,mpi_real,pe,pe,mpi_comm_comp,jstat,ierr)
            nn=0
            do ks=kss,kse
              do j=loc_jts_b1,loc_jte_b1
                do i=loc_its_b1,loc_ite_b1
                  nn=nn+1
                  n=i+(j-1)*ide+(ks-1)*ide*jde
                  s1_glob(n) = s1_pe_loc(nn)
                enddo
              enddo
            enddo
            deallocate(s1_pe_loc)
          end do

          do ks=kss,kse
            do j=jds+1,jde-1
              do i=ids+1,ide-1
                n=i+(j-1)*ide+(ks-1)*ide*jde
                if(s1_glob(n).gt.0.0d0) then
                  gsums(2*ks-1) = gsums(2*ks-1) + s1_glob(n)
                else
                  gsums(2*ks  ) = gsums(2*ks  ) + s1_glob(n)
                endif
              enddo
            enddo
          enddo

        else
          loc_lngth=(ite_b1-its_b1+1)*(jte_b1-jts_b1+1)*kse

          n=0
          do ks=kss,kse
            do j=jts_b1,jte_b1
              do i=its_b1,ite_b1
                n=n+1
                s1_loc(n) = s1l_sum(i,j,ks)
              enddo
            enddo
          enddo

          call mpi_send(s1_loc(1:loc_lngth),loc_lngth &
                       ,mpi_real, 0, mype, mpi_comm_comp, ierr)

        endif

        call mpi_bcast(gsums,kse*2 &
                      ,mpi_double_precision,0,mpi_comm_comp,ierr)

!----------------------------------------------------------------------
!
      endif  global_reduce
!
!----------------------------------------------------------------------
!.......................................................................
!$omp parallel do private (dsp,i,j,ks,l,rfacs,sfacs,sumns,sumps)
!.......................................................................
      do ks=kss,kse
        sumps=gsums(2*ks-1)
        sumns=gsums(2*ks  )
!
        if(sumps*(-sumns).gt.1.) then
          sfacs=-sumns/sumps
          rfacs=1./sfacs
        else
          sfacs=0.
          rfacs=0.
        endif
!
        do l=1,lm
          do j=jts_b1,jte_b1
            do i=its_b1,ite_b1
              dsp=s1(i,j,l,ks)*rdvol(i,j,l)
              if(sfacs.lt.1.) then
                if(dsp.gt.0.) dsp=dsp*sfacs
              else
                if(dsp.lt.0.) dsp=dsp*rfacs
              endif
              tcs(i,j,l,ks)=tcs(i,j,l,ks)+dsp
            enddo
          enddo
        enddo
!-----------------------------------------------------------------------
!
      enddo ! end of the loop by species
!.......................................................................
!$omp end parallel do 
!.......................................................................
!
!-----------------------------------------------------------------------
!
!***  Interpolate q2 tendencies and q2 itself from midlayers back to 
!     interfaces
!
q2_check: if (kss<=indx_q2 .and. indx_q2<=kse) then
      do l=1,lm
        do j=jts_b1,jte_b1
          do i=its_b1,ite_b1
            tcs(i,j,l,indx_q2)=(dsg2(l)*pd(i,j)+pdsg1(l))*tcs(i,j,l,indx_q2)
          enddo
        enddo
      enddo
      do l=1,lm-1
        do j=jts_b1,jte_b1
          do i=its_b1,ite_b1
            tcs(i,j,l,indx_q2)=(tcs(i,j,l,indx_q2)+tcs(i,j,l+1,indx_q2)) &
                        /((dsg2(l  )*pd(i,j)+pdsg1(l  )) &
                         +(dsg2(l+1)*pd(i,j)+pdsg1(l+1)))
            s(i,j,l,indx_q2)=max(0.5*(s(i,j,l,indx_q2)+s(i,j,l+1,indx_q2)),epsq2(l))
          enddo
        enddo
      enddo
      do j=jts_b1,jte_b1
        do i=its_b1,ite_b1
          tcs(i,j,lm,indx_q2)=0.
          s(i,j,lm,indx_q2)=epsq2(lm)
        enddo
      enddo
    endif  q2_check
!
!-----------------------------------------------------------------------
!
                        endsubroutine mono
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
!----------------------------------------------------------------------
!######################################################################
!----------------------------------------------------------------------
!
                        end module module_dynamics_routines
!
!-----------------------------------------------------------------------
