!*********************************************************************
module module_wrfphys_alloc
!	This module allocates variables used in wrfphys, used chem_alloc2 from
!       T Henderson as example...
!  	GG   April 2009
!*********************************************************************
contains




subroutine wrfphys_alloc

  USE module_wrf_control, only: ims,ime,jms,jme,kms,kme,  &
                                num_moist,num_soil_layers,num_scalar
  ! yes, use ALL of it
  use module_wrfphysvars

  implicit none

  ALLOCATE( moist( ims:ime, kms:kme, jms:jme, num_moist ) )
  moist=0.
  ALLOCATE( scalar( ims:ime, kms:kme, jms:jme, num_scalar ) )
  scalar=0.
  ALLOCATE( tsk( ims:ime, jms:jme) )
  ALLOCATE( rri( ims:ime , kms:kme , jms:jme ) ) 
  ALLOCATE( t_phy( ims:ime , kms:kme , jms:jme ) )
  ALLOCATE( th_phy( ims:ime , kms:kme , jms:jme ) )
  ALLOCATE( p_phy( ims:ime , kms:kme , jms:jme ) )
  ALLOCATE( pi_phy( ims:ime , kms:kme , jms:jme ) )
  ALLOCATE( dz8w( ims:ime , kms:kme , jms:jme ) ) 
  ALLOCATE( t8w( ims:ime , kms:kme , jms:jme ) )
  ALLOCATE( p8w( ims:ime , kms:kme , jms:jme ) )
  ALLOCATE( z_at_w ( ims:ime , kms:kme , jms:jme ) )
  ALLOCATE( zmid ( ims:ime , kms:kme , jms:jme ) )
  ALLOCATE( u_phy( ims:ime , kms:kme , jms:jme ) )
  ALLOCATE( v_phy( ims:ime , kms:kme , jms:jme ) )
  ALLOCATE( vvel( ims:ime , kms:kme , jms:jme ) )
  ALLOCATE( rho_phy( ims:ime , kms:kme , jms:jme ) )
  ALLOCATE( exch_h( ims:ime , kms:kme , jms:jme ) )
  ALLOCATE( cldfra( ims:ime , kms:kme , jms:jme ) )
  ALLOCATE( rqvcuten( ims:ime , kms:kme , jms:jme ) )
  ALLOCATE( rqvblten( ims:ime , kms:kme , jms:jme ) )
  ALLOCATE( rqvften( ims:ime , kms:kme , jms:jme ) )
  ALLOCATE( rthcuten( ims:ime , kms:kme , jms:jme ) )
  ALLOCATE( rthblten( ims:ime , kms:kme , jms:jme ) )
  ALLOCATE( rthraten( ims:ime , kms:kme , jms:jme ) )
  ALLOCATE( rthften( ims:ime , kms:kme , jms:jme ) )
  ALLOCATE( rqccuten( ims:ime , kms:kme , jms:jme ) )
  ALLOCATE( rqrcuten( ims:ime , kms:kme , jms:jme ) )
  ALLOCATE( rqscuten( ims:ime , kms:kme , jms:jme ) )
  ALLOCATE( rqicuten( ims:ime , kms:kme , jms:jme ) )
  ALLOCATE( rqgcuten( ims:ime , kms:kme , jms:jme ) )
  
  ALLOCATE( ivgtyp( ims:ime , jms:jme ) )
  ALLOCATE( isltyp( ims:ime , jms:jme ) )
  ALLOCATE( u10( ims:ime , jms:jme ) )
  ALLOCATE( v10( ims:ime , jms:jme ) )
  ALLOCATE( gsw( ims:ime , jms:jme ) )
  ALLOCATE( sr( ims:ime , jms:jme ) )
  ALLOCATE( pbl( ims:ime , jms:jme ) )
  ALLOCATE( hfx( ims:ime , jms:jme ) )
  ALLOCATE( vegfra( ims:ime , jms:jme ) )
  ALLOCATE( rmol( ims:ime , jms:jme ) )
  ALLOCATE( ust( ims:ime , jms:jme ) )
  ALLOCATE( xland( ims:ime , jms:jme ) )
  ALLOCATE( xlat( ims:ime , jms:jme ) )
  ALLOCATE( xlong( ims:ime , jms:jme ) )
  ALLOCATE( znt( ims:ime , jms:jme ) )
  ALLOCATE( ht( ims:ime , jms:jme ) )
! for convective schemes
  ALLOCATE(  rainc( ims:ime , jms:jme ))
  ALLOCATE(  apr_gr( ims:ime , jms:jme ))
  ALLOCATE(  apr_w( ims:ime , jms:jme ))
  ALLOCATE(  apr_mc( ims:ime , jms:jme ))
  ALLOCATE(  apr_as( ims:ime , jms:jme ))
  ALLOCATE(  apr_st( ims:ime , jms:jme ))
  ALLOCATE(  apr_capma( ims:ime , jms:jme ))
  ALLOCATE(  apr_capme( ims:ime , jms:jme ))
  ALLOCATE(  apr_capmi( ims:ime , jms:jme ))
  ALLOCATE(  mass_flux( ims:ime , jms:jme ))
  ALLOCATE(  cugd_tten( ims:ime , kms:kme , jms:jme ))
  ALLOCATE(  cugd_ttens( ims:ime , kms:kme , jms:jme ))
  ALLOCATE(  cugd_qvten( ims:ime , kms:kme , jms:jme ))
  ALLOCATE(  cugd_qcten( ims:ime , kms:kme , jms:jme ))
  ALLOCATE(  cugd_qvtens( ims:ime , kms:kme , jms:jme ))
  ALLOCATE(  gd_cloud( ims:ime , kms:kme , jms:jme ))
  ALLOCATE(  gd_cloud2( ims:ime , kms:kme , jms:jme ))
  ALLOCATE(  raincv( ims:ime , jms:jme ))

  ALLOCATE ( rainnc( ims:ime , jms:jme ))
  ALLOCATE ( rainncv( ims:ime , jms:jme ))
  ALLOCATE ( snownc( ims:ime , jms:jme ))
  ALLOCATE ( snowncv( ims:ime , jms:jme ))
  ALLOCATE ( graupelnc( ims:ime , jms:jme ))
  ALLOCATE ( graupelncv( ims:ime , jms:jme ))

  ALLOCATE( dxy( ims:ime , jms:jme ) )
 ALLOCATE( smois( ims:ime, num_soil_layers, jms:jme ) )
return
end subroutine  wrfphys_alloc

end module module_wrfphys_alloc
