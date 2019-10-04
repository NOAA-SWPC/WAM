MODULE module_wetdep_ls
  USE module_initial_chem_namelists
CONTAINS
subroutine wetdep_ls(dt,var,rain,moist,rho,var_rmv,num_moist, &
         num_chem,numgas,p_qc,dz8w,vvel,chem_opt,             &
         ids,ide, jds,jde, kds,kde,                                        &
         ims,ime, jms,jme, kms,kme,                                        &
         its,ite, jts,jte, kts,kte                                         )
  IMPLICIT NONE

   INTEGER,      INTENT(IN   ) :: num_chem,numgas,num_moist,p_qc,                          &
                                  chem_opt,ids,ide, jds,jde, kds,kde,               &
                                  ims,ime, jms,jme, kms,kme,               &
                                  its,ite, jts,jte, kts,kte
   real, INTENT(IN ) :: dt
    REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_moist ),                &
         INTENT(IN ) ::                                   moist
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme ),                        &
          INTENT(IN   ) :: rho,dz8w,vvel        
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme ,1:num_chem),                        &
          INTENT(INOUT) :: var        
   REAL,  DIMENSION( jms:jme ),                                  &
          INTENT(IN   ) :: rain
   REAL,  DIMENSION( ims:ime ,  jms:jme,num_chem ),                                  &
          INTENT(INOUT   ) :: var_rmv
   REAL,  DIMENSION( its:ite ,  jts:jte ) :: var_sum
   REAL,  DIMENSION( its:ite ,  kts:kte, jts:jte ) :: var_rmvl
   REAL,  DIMENSION( its:ite ,  jts:jte ) :: frc,var_sum_clw,rain_clw     
    real :: dvar,factor,clsum,alpha,rho_water
   integer :: nv,i,j,k,km,kb,kbeg
    rho_water = 1000.
    var_rmv (:,:,:)=0.
!   write(6,*) 'in wetdepls, p_qc = ',p_qc
!   nv=p_bc1
    do nv=1,num_chem
    alpha = .5    ! scavenging factor
    if(chem_opt >= 300 .and. chem_opt < 500)then
       if(nv.le. numgas .and. nv.ne.p_sulf)cycle
       if(nv.eq.p_bc1 .or. nv.eq.p_oc1 .or. nv.eq.p_dms)cycle
       if(nv.eq.p_sulf .or. nv.eq.p_seas_1 .or. nv.eq.p_seas_2 .or. &
          nv.eq.p_seas_3 .or. nv.eq.p_seas_4)alpha=1.
       if(nv.eq.p_bc2 .or. nv.eq.p_oc2)alpha=0.8
    endif
    do i=its,ite
    do j=jts,jte
     var_sum_clw(i,j)=0.
     var_sum(i,j)=0.
     var_rmvl(i,:,j)=0.
     frc(i,j)=0.
     rain_clw(i,j)=0.
     if(rain(j).gt.1.e-6)then
! convert rain back to rate
!
        rain_clw(i,j)=rain(j)/dt
! total cloud water
!
        do k=1,kte-1
           dvar=max(0.,moist(i,k,j,p_qc)*rho(i,k,j)*vvel(i,k,j)*dz8w(i,k,j))
           var_sum_clw(i,j)=var_sum_clw(i,j)+dvar
           var_sum(i,j)=var_sum(i,j)+var(i,k,j,nv)*rho(i,k,j)
        enddo
        if(var_sum(i,j).gt.1.e-8 .and. var_sum_clw(i,j).gt.1.e-6 ) then
!        assuming that frc is onstant, it is my conversion factor 
!       (just like in convec. parameterization
           frc(i,j)=rain_clw(i,j)/var_sum_clw(i,j)
!    print *,'frc ', frc(i,j),var_sum_clw(i,j),var_sum(i,j)
           frc(i,j)=max(1.e-6,min(frc(i,j),.5))
        endif
     endif
    enddo
    enddo
!
! get rid of it
!
    do i=its,ite
    do j=jts,jte
     if(rain(j).gt.1.e-6 .and. var_sum(i,j).gt.1.e-8 .and. var_sum_clw(i,j).gt.1.e-6)then
       do k=kts,kte-2
        if(var(i,k,j,nv).gt.1.e-10 .and. moist(i,k,j,p_qc).gt.0.)then
        factor = max(0.,frc(i,j)*rho(i,k,j)*dz8w(i,k,j)*vvel(i,k,j))
!       print *,'var before ',k,km,var(i,k,j,nv),factor
!       dvar=.05*alpha*factor/(1+factor)*var(i,k,j,nv)
        dvar=max(0.,.05*alpha*factor/(1+factor)*var(i,k,j,nv))
        dvar=min(dvar,var(i,k,j,nv))
        var_rmvl(i,k,j)=dvar
        if((var(i,k,j,nv)-dvar).lt.1.e-12)then
           dvar=var(i,k,j,nv)-1.e-12
           var(i,k,j,nv)=var(i,k,j,nv)-dvar
        else
           var(i,k,j,nv)=var(i,k,j,nv)-dvar
        endif
        var_rmv(i,j,nv)=var_rmv(i,j,nv)+var_rmvl(i,k,j)
!       print *,'var after ',km,var(i,k,j,nv),dvar
        endif
       enddo
!      var_rmv(i,j)=var_rmv(i,j)+var_rmvl(i,j)
    endif
    enddo
    enddo
    enddo
END SUBROUTINE WETDEP_LS
END MODULE module_wetdep_ls


