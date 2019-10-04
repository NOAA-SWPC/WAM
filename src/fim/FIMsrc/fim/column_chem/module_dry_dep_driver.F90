!WRF:MODEL_LAYER:CHEMICS
!
MODULE module_dry_dep_driver
  IMPLICIT NONE

CONTAINS

    subroutine dry_dep_driver(ktau,dtstep,                &
               moist,p8w,alt,                                             &
               chem,rho_phy,dz8w,exch_h,hfx,                              &
               ivgtyp,tsk,pbl,ust,znt,z,z_at_w,                           &
               xland,dep_vel_o3,g,                                        &
               e_co,kemit,numgas,                                         &
               num_chem,num_moist,                                        &
               ids,ide, jds,jde, kds,kde,                                 &
               ims,ime, jms,jme, kms,kme,                                 &
               its,ite, jts,jte, kts,kte                                  )
!----------------------------------------------------------------------
! USE module_model_constants
! USE module_configure
! USE module_state_description
! USE module_dep_simple
  USE module_initial_chem_namelists,only:p_o3,p_dust_1,p_vash_1,p_vash_4,p_vash_10,p_dms
  USE module_vertmx_wrf
  USE module_chemvars,only:epsilc
  USE module_initial_chem_namelists,only: chem_opt,drydep_opt,wesely,     &
            chem_tracer,gocart_simple,GOCARTRACM_KPP,RADM2SORG,RADM2SORG, &
            RADM2SORG_AQ,RADM2SORG_KPP,RACMSORG_AQ,RACMSORG_KPP, &
            CBMZ_MOSAIC_4BIN, CBMZ_MOSAIC_8BIN, CBMZ_MOSAIC_4BIN_AQ,      &
            CBMZ_MOSAIC_8BIN_AQ
! USE module_data_sorgam
! USE module_aerosols_sorgam
! USE module_gocart_settling
  USE module_gocart_drydep,only: gocart_drydep_driver
! USE module_mosaic_drydep, only:  mosaic_drydep_driver
! USE module_mixactivate_wrappers, only: mosaic_mixactivate, sorgam_mixactivate
  IMPLICIT NONE

   INTEGER,      INTENT(IN   ) :: numgas,                       &
               num_chem,num_moist,                                        &
                                  ids,ide, jds,jde, kds,kde,    &
                                  ims,ime, jms,jme, kms,kme,    &
                                  its,ite, jts,jte, kts,kte
   INTEGER,      INTENT(IN   ) ::                               &
                                  ktau
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_moist ),        &
         INTENT(IN ) ::                                   moist
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),         &
         INTENT(INOUT ) ::                                 chem

   INTEGER,      INTENT(IN   ) :: kemit
   REAL, DIMENSION( ims:ime, kms:kemit, jms:jme ),            &
         INTENT(IN ) ::                                                    &
          e_co




   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme )         ,    &
          INTENT(IN   ) ::                                      &
                                                        alt,    &
                                                      dz8w,     &
                                              p8w,z_at_w ,  &
                                              exch_h,rho_phy,z
   INTEGER,DIMENSION( ims:ime , jms:jme )                  ,    &
          INTENT(IN   ) ::                                      &
                                                     ivgtyp
   REAL,  DIMENSION( ims:ime , jms:jme )                   ,    &
          INTENT(INOUT) ::                                      &
                                                     tsk,       &
                                                     pbl,       &
                                                     ust,       &
                                                     hfx,       &
                                                   xland,znt
   REAL,  DIMENSION( ims:ime , jms:jme )                   ,    &
          INTENT(OUT) ::                                      &
                                                     dep_vel_o3

      REAL,      INTENT(IN   ) ::                               &
                             dtstep,g

!--- deposition and emissions stuff
! .. Parameters ..
! ..
! .. Local Scalars ..
      REAL ::  clwchem,  dvfog, dvpart,  &
        rad, rhchem, ta, ustar, z1,zntt

      INTEGER :: iland, iprt, iseason, jce, jcs,  &
                 n, nr, ipr, jpr, nvr,   &
                 idrydep_onoff

      LOGICAL :: highnh3, rainflag, vegflag, wetflag
!     CHARACTER (4) :: luse_typ,mminlu_loc
! ..
! .. Local Arrays ..
      REAL :: p(kts:kte)
   REAL, DIMENSION( its:ite, jts:jte, num_chem ) ::   ddvel

!  REAL,  DIMENSION( ims:ime , kms:kme , jms:jme ) :: dryrho_phy
   REAL,  DIMENSION( kms:kme ) :: dryrho_1d

! turbulent transport
      real :: pblst(kts:kte),ekmfull(kts:kte+1),zzfull(kts:kte+1),zz(kts:kte)
      integer :: ii,jj,kk,i,j,k,nv
!
! necessary for aerosols (module dependent)
!
!  REAL, DIMENSION( its:ite, jts:jte ) ::   aer_res

! ..
! .. Intrinsic Functions ..
      INTRINSIC max, min

!
! compute dry deposition velocities = ddvel
!
! 28-jun-2005 rce - initialize ddvel=0; call aerosol drydep routine
!           only when drydep_opt == WESELY
!       the wesely_driver routine computes aer_res, and currently
!	    you cannot compute aerosol drydep without it !!
! 08-jul-2005 rce - pass idrydep_onoff to mixactivate routines
!
!  write(6,*)'call dry dep driver'
   dep_vel_o3(:,:)=0.
   ddvel(:,:,:) = 0.0
   idrydep_onoff = 0

!  drydep_select: SELECT CASE(drydep_opt)

!    CASE ( WESELY )
!
! drydep_opt == WESELY means 
!     wesely for gases 
!     other (appropriate) routine for aerosols
!
!      CALL wrf_debug(15,'DOING DRY DEP VELOCITIES WITH WESELY METHOD')

!      IF( chem_opt /= CHEM_TRACER  .and. chem_opt /= GOCART_SIMPLE ) THEN
!         call wesely_driver(id,ktau,dtstep,                                 &
!              config_flags,                                              &
!              gmt,julday,t_phy,moist,p8w,t8w,                            &
!              p_phy,chem,rho_phy,dz8w,ddvel,aer_res,                     &
!              ivgtyp,tsk,gsw,vegfra,pbl,rmol,ust,znt,xlat,xlong,z,z_at_w,&
!              numgas,                                                    &
!              ids,ide, jds,jde, kds,kde,                                 &
!              ims,ime, jms,jme, kms,kme,                                 &
!              its,ite, jts,jte, kts,kte                                  )
       IF (( chem_opt == GOCART_SIMPLE ) .or.            &
              ( chem_opt == GOCARTRACM_KPP)  .or.            &
              ( chem_opt == 316)  .or.            &
              ( chem_opt == 317)  .or.            &
!             ( chem_opt == 502)  .or.            &
                (chem_opt == 304          )) then
!
! this does aerosol species (dust,seas, bc,oc) for gocart only
!,
         call gocart_drydep_driver(0,                     &
               moist,p8w,chem,rho_phy,dz8w,ddvel,xland,hfx,    &
               ivgtyp,tsk,pbl,ust,znt,g,                       &
               num_moist,num_chem,                             &
               ids,ide, jds,jde, kds,kde,                      &
               ims,ime, jms,jme, kms,kme,                      &
               its,ite, jts,jte, kts,kte                       )
!         ddvel(:,:,p_vash_1:num_chem) = 0.
          ddvel(:,:,p_dms) = 0.
       ELSE if (chem_opt == 501 ) then
! for caesium .1cm/s
!
          ddvel(:,:,:)=.001
       ELSE
          !Set dry deposition velocity to zero when using the
          !chemistry tracer mode.
          ddvel(:,:,:) = 0.
!         write(6,*)'no dry deposition '
       END IF

       idrydep_onoff = 1





!   This will be called later from subgrd_transport_driver.F !!!!!!!!
!
!
      do 100 j=jts,jte
      do 100 i=its,ite
      if(p_dust_1.gt.1)dep_vel_o3(i,j)=ddvel(i,j,p_dust_1)
      pblst=0.
!
!
!-- start with vertical mixing
!
      do k=kts,kte+1
         zzfull(k)=z_at_w(i,k,j)-z_at_w(i,kts,j)
      enddo
      do k=kts,kte
         ekmfull(k)=max(1.e-6,exch_h(i,k,j))
      enddo
      ekmfull(kts)=0.
      ekmfull(kte+1)=0.

!!$! UNCOMMENT THIS AND FINE TUNE LEVELS TO YOUR DOMAIN IF YOU WANT TO
!!$! FORCE MIXING TO A CERTAIN DEPTH:
!!$!
!!$! --- Mix the emissions up several layers
!!$!     if e_co > 0., the grid cell should not be over water
!!$!     if e_co > 200, the grid cell should be over a large urban region
!!$!
!       if (e_co(i,kts,j) .gt. 0) then
!          ekmfull(kts:kts+10) = max(ekmfull(kts:kts+10),1.)
!       endif
!       if (e_co(i,kts,j) .gt. 200) then
!          ekmfull(kts:kte/2) = max(ekmfull(kts:kte/2),2.)
!       endif
!
!
      do k=kts,kte
         zz(k)=z(i,k,j)-z_at_w(i,kts,j)
      enddo
!
!   vertical mixing routine (including deposition)
!   need to be careful here with that dumm tracer in spot 1
!   do not need lho,lho2
!   (03-may-2006 rce - calc dryrho_1d and pass it to vertmx)
!
!     if(j.eq.681)write(6,*)ddvel(1,681,1:num_chem)
!     if(p_o3.gt.1)dep_vel_o3(i,j)=ddvel(i,j,p_o3)
      do nv=1,num_chem-0
         do k=kts,kte
            pblst(k)=max(epsilc,chem(i,k,j,nv))
            dryrho_1d(k) = 1./alt(i,k,j)
!           if(j.eq.681.and.nv.eq.10)then
!             write(6,*)k,chem(i,k,j,nv),exch_h(i,k,j),ddvel(i,j,nv)
!             write(6,*)dryrho_1d(k),zz(k),zzfull(k)
!           endif
         enddo

         mix_select: SELECT CASE(chem_opt)
         CASE (RADM2SORG_AQ, RACMSORG_AQ, CBMZ_MOSAIC_4BIN_AQ, CBMZ_MOSAIC_8BIN_AQ)
!           if(.not.is_aerosol(nv))then ! mix gases not aerosol
               call vertmx(dtstep,pblst,ekmfull,dryrho_1d, &
                           zzfull,zz,ddvel(i,j,nv),kts,kte)

!           endif

         CASE DEFAULT
               call vertmx(dtstep,pblst,ekmfull,dryrho_1d, &
                           zzfull,zz,ddvel(i,j,nv),kts,kte)

         END SELECT mix_select

         do k=kts,kte-1
!           if(j.eq.681.and.nv.eq.10)then
!             write(6,*)pblst(k)
!           endif
            chem(i,k,j,nv)=max(epsilc,pblst(k))
!           if(j.eq.681.and.nv.eq.10)then
!             write(6,*)dtstep,pblst(k),chem(i,k,j,nv)
!           endif
!           if(j.eq.75.and.nv.eq.16)write(6,*)zzfull(k),chem(i,k,j,nv),alt(i,k,j)
         enddo
      enddo
100   continue
!
!  vertical mixing and activation of aerosol
!
!  where( alt(its:ite,kts:kte,jts:jte) /= 0. )  !get dry density to conserve mass in mixactivate, wig, 24-apr-2006
!     dryrho_phy(its:ite,kts:kte,jts:jte) = 1./alt(its:ite,kts:kte,jts:jte)
!  elsewhere
!     dryrho_phy(its:ite,kts:kte,jts:jte) = 0.
!  end where
!  dryrho_phy(its:ite,kte+1,jts:jte) = 0. !wig: testing, should never need this

!  mixactivate_select: SELECT CASE(config_flags%chem_opt)

!  CASE (RADM2SORG_AQ, RACMSORG_AQ)
!     call sorgam_mixactivate (                        &
!	id, ktau, dtstep, config_flags, idrydep_onoff,   &
!	dryrho_phy, t_phy, w, cldfra, cldfra_old, &
!	ddvel, z, dz8w, p8w, t8w, exch_h,         &
!	moist(ims,kms,jms,P_QV), moist(ims,kms,jms,P_QC), moist(ims,kms,jms,P_QI), &
!       scalar(ims,kms,jms,P_QNDROP), f_qc, f_qi, chem, &
!       ccn1, ccn2, ccn3, ccn4, ccn5, ccn6, nsource,       &
!	ids,ide, jds,jde, kds,kde,                        &
!	ims,ime, jms,jme, kms,kme,                        &
!	its,ite, jts,jte, kts,kte                         )
!  CASE (CBMZ_MOSAIC_4BIN_AQ, CBMZ_MOSAIC_8BIN_AQ)
!      CALL wrf_debug(15,'call mixactive for mosaic aerosol')
!     call mosaic_mixactivate (                        &
!	id, ktau, dtstep, config_flags, idrydep_onoff,   &
!	dryrho_phy, t_phy, w, cldfra, cldfra_old, &
!	ddvel, z, dz8w, p8w, t8w, exch_h,         &
!	moist(ims,kms,jms,P_QV), moist(ims,kms,jms,P_QC), moist(ims,kms,jms,P_QI), &
!       scalar(ims,kms,jms,P_QNDROP), f_qc, f_qi, chem,   &
!       ccn1, ccn2, ccn3, ccn4, ccn5, ccn6, nsource,      &
!	ids,ide, jds,jde, kds,kde,                        &
!	ims,ime, jms,jme, kms,kme,                        &
!	its,ite, jts,jte, kts,kte                         )
!  CASE DEFAULT
!  END SELECT mixactivate_select
!  settling_select: SELECT CASE(config_flags%chem_opt)
!  CASE (GOCART_SIMPLE,GOCARTRACM_KPP)
!      CALL wrf_debug(15,'call gocart settling routine')
!        call gocart_settling_driver(dtstep,config_flags,t_phy,moist,  &
!        chem,rho_phy,dz8w,p8w,p_phy,         &
!        dx,g, &
!        ids,ide, jds,jde, kds,kde,                                        &
!        ims,ime, jms,jme, kms,kme,                                        &
!        its,ite, jts,jte, kts,kte                                         )
!  CASE DEFAULT
!      CALL wrf_debug(15,'no settling routine')
!  END SELECT settling_select

!      CALL wrf_debug(15,'end of dry_dep_driver')

END SUBROUTINE dry_dep_driver

END MODULE module_dry_dep_driver
