!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read required dynamics fields from the restart file. SMS will modify the code to do the
! appropriate single-task reading of the restart file, and subsequent scattering of the
! data to other MPI tasks.
!
! CRITICAL: If you modify this file, you MUST also modify write_restart_dyn.F90 in the
! same way. Otherwise restart will be broken.
!
! read_restart_dyn and write_restart_dyn belong in a module, but SMS doesn't like multiple
! subroutines in a file.
!
! readarr32 assumes 32-bit data.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine read_restart_dyn (unitno)
  use module_sfc_variables, only: rn2d, rc2d, rn2d0, rc2d0, rg2d0, flxswavg2d, flxlwavg2d, flxlwtoa2d
  use module_variables,     only: curr_write_time, nf, of, vof, adbash1, adbash2, adbash3, psrf, &
                                  ptdcy, pw2d, u_tdcy, v_tdcy, dp_tdcy, dpl_tdcy, tr3d, trdp, &
                                  trc_tdcy, trl_tdcy, us3d, vs3d, ws3d, mp3d, tk3d, dp3d, rh3d, &
                                  vor, dpinit, pr3d, ex3d, ph3d, sdot, massfx, cumufx
  use module_constants,     only: dpsig, thetac, lat, lon, nprox, proxs, area, cs, sn, sidevec_c, &
                                  sidevec_e, sideln, rsideln, rprox_ln, area, rarea, corio, &
                                  deg_lat, deg_lon
  use module_control,       only: dt, nabl, ntra, ntrb, npp, nvl, nvlp1
  use module_globsum,       only: qmstr, qmstrc, qmstrn, qdtr_set

  implicit none

  integer, intent(in) :: unitno ! unit number to read from

  integer :: n, t ! indices

!SMS$SERIAL BEGIN    
  read (unitno, err=90) curr_write_time, nf, of, vof, adbash1, adbash2, adbash3, thetac, dpsig
  read (unitno, err=90) lat, lon, nprox, proxs, area, cs, sn, psrf, ptdcy
!SMS$SERIAL END
!SMS$SERIAL BEGIN
  read (unitno, err=90) qmstr, qmstrc, qmstrn, qdtr_set
!SMS$SERIAL END
!SMS$SERIAL BEGIN    
  read (unitno, err=90) sidevec_c, sidevec_e, sideln, rprox_ln
!SMS$SERIAL END
!SMS$SERIAL BEGIN
  read (unitno, err=90) rn2d0, rc2d0, rg2d0, flxswavg2d, flxlwavg2d, flxlwtoa2d
!SMS$SERIAL END
!SMS$SERIAL BEGIN    
  read (unitno, err=90) rarea, rsideln, corio, deg_lat, deg_lon, rn2d, pw2d, rc2d
!SMS$SERIAL END

! These arrays are dimensioned (nvl,nip[,other dimensions]):
  do n=1,nabl
    call readarr32 (u_tdcy(:,:,n), nvl, unitno)
    call readarr32 (v_tdcy(:,:,n), nvl, unitno)
    call readarr32 (dp_tdcy(:,:,n), nvl, unitno)
    call readarr32 (dpl_tdcy(:,:,n), nvl, unitno)
  end do

  do t=1,ntra+ntrb
    call readarr32 (tr3d(:,:,t), nvl, unitno)
    call readarr32 (trdp(:,:,t), nvl, unitno)
    do n=1,nabl
      call readarr32 (trc_tdcy(:,:,n,t), nvl, unitno)
      call readarr32 (trl_tdcy(:,:,n,t), nvl, unitno)
    end do
  end do

  call readarr32 (us3d, nvl, unitno)
  call readarr32 (vs3d, nvl, unitno)
  call readarr32 (ws3d, nvl, unitno)
  call readarr32 (mp3d, nvl, unitno)
  call readarr32 (tk3d, nvl, unitno)
  call readarr32 (dp3d, nvl, unitno)
  call readarr32 (rh3d, nvl, unitno)
  call readarr32 (vor, nvl, unitno)
  call readarr32 (dpinit, nvl, unitno)

! These arrays are dimensioned (nvlp1,nip):
  call readarr32 (pr3d, nvlp1, unitno)
  call readarr32 (ex3d, nvlp1, unitno)
  call readarr32 (ph3d, nvlp1, unitno)
  call readarr32 (sdot, nvlp1, unitno)

! These arrays are dimensioned (nvl,npp,nip[,other dimensions]):
! Simplest coding folds nvl*npp into a single dimension. Will need to rewrite massfx
! and cumufx to a transpose if root process can't hold npp 3-d fields in memory
  do n=1,nabl
    call readarr32 (massfx(:,:,:,n), nvl*npp, unitno)
  end do
  call readarr32 (cumufx(:,:,:), nvl*npp, unitno)

  write(6,*) 'read_restart_dyn: successfully read dynamics fields from restart file'

  return

90 write(6,*)'read_restart_dyn: Error reading from unit ', unitno, '. Stopping'
  stop
end subroutine read_restart_dyn
