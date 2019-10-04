  program newname
  use module_control  ,only: nip,glvl,curve,control
!SMS$ignore begin
  USE slint, ONLY: bilinear_init, bl_int
!SMS$ignore end
  implicit none
  integer, parameter :: imax=288,jmax=181,nspecies=25,imax2=360, &
                        klev=55,iklev=56

  real(4) f1(imax*jmax),f2(imax*jmax),f3(imax*jmax),f4(imax*jmax)
  real(4) tmp_h2o2(imax,jmax,klev),tmp_oh(imax,jmax,klev), &
          tmp_no3(imax,jmax,klev),gocart_lev(iklev)
  real(4) emissions(imax*jmax,nspecies),minv,maxv
  real(4) p_gocart(iklev)
  integer k,nv,nv_g,itime

  integer ipn,nvp,ios
  real(4), allocatable :: ps_i(:)   ! ps_i(nip)
  real(4), allocatable :: oh(:,:)   ! oh(klev,nip)
  CHARACTER (LEN=7) :: ename(nspecies)
  CHARACTER (LEN=20) :: dname
  CHARACTER(len=80)   :: grid_file1,grid_file2,g3
  DATA  ename/'e_so2','e_no','e_ald','e_hcho','e_ora2','e_nh3','e_hc3','e_hc5','e_hc8',   &
               'e_eth','e_co','e_ol2','e_olt','e_oli','e_tol','e_xyl','e_ket','e_csl',     &
               'e_iso','e_pm_25','e_pm_10','e_oc','e_bc','e_dms','e_sulf'/

! read FIM namelists
  call control(.true.)

! allocate arrays
  ALLOCATE(ps_i(nip))
  ALLOCATE(oh(klev,nip))

! set up interpolation from chem grid to icos
  grid_file1 = "chemltln.dat"
  grid_file2 = "glvl.dat"
  write(6,*) 'newname:  generating interpolation weights from ',TRIM(grid_file1),' to ',TRIM(grid_file2),' ...'
  OPEN(66,file=grid_file2,status='old',form='unformatted',iostat=ios)
  if (ios /= 0) then
    write(6,*) 'ERROR:  newname failed to open file ',TRIM(grid_file1)
    stop
  endif
  call TestGlvlHeader (66,     grid_file2,'newname',glvl)
  call TestCurveHeader(66,     grid_file2,'newname',curve)
  CALL bilinear_init(grid_file1, imax*jmax, 66, nip)
  close(66)
!
! first do dust erosion map
!
  write(6,*) 'newname:  interpolating erod_binary ...'
  open(unit=21,file='erod_binary',form='unformatted',status='old',iostat=ios)
  if (ios /= 0) then
    write(6,*) 'ERROR:  failed to open file erod_binary'
    stop
  endif
  read(21)f1,f2,f3
  write(6,*)maxval(f1)
  write(6,*)minval(f1)

  CALL bl_int (f1, ps_i)  ! unit in pascal
  g3 = "erod1.dat"
  write(6,*)g3
  open(unit=23,file=g3,form='unformatted',iostat=ios)
  if (ios /= 0) then
    write(6,*) 'ERROR:  failed to open file ',TRIM(g3)
    stop
  endif
  write(23)ps_i
  write(6,*)maxval(ps_i)
  write(6,*)minval(ps_i)


  write(6,*)maxval(f2)
  write(6,*)minval(f2)
  CALL bl_int (f2, ps_i)  ! unit in pascal
  g3 = "erod2.dat"
  write(6,*)g3
  open(unit=24,file=g3,form='unformatted',iostat=ios)
  if (ios /= 0) then
    write(6,*) 'ERROR:  failed to open file ',TRIM(g3)
    stop
  endif
  write(24)ps_i
  write(6,*)maxval(ps_i)
  write(6,*)minval(ps_i)


  write(6,*)maxval(f3)
  write(6,*)minval(f3)
  CALL bl_int (f3, ps_i)  ! unit in pascal
  g3 = "erod3.dat"
  write(6,*)g3
  open(unit=25,file=g3,form='unformatted',iostat=ios)
  if (ios /= 0) then
    write(6,*) 'ERROR:  failed to open file ',TRIM(g3)
    stop
  endif
  write(25)ps_i

  close (21)
  close (23)
  close (24)
  close (25)
!
!  dms emissions
!
  write(6,*) 'newname:  interpolating dm0_binary ...'
  open(unit=22,file='dm0_binary',form='unformatted',status='old',iostat=ios)
  if (ios /= 0) then
    write(6,*) 'ERROR:  failed to open file dm0_binary'
    stop
  endif
  read(22)f1
  close (22)
  CALL bl_int (f1, ps_i)  ! unit in pascal
  g3 = "dm0.dat"
  write(6,*)g3
  open(unit=26,file=g3,form='unformatted',iostat=ios)
  if (ios /= 0) then
    write(6,*) 'ERROR:  failed to open file ',TRIM(g3)
    stop
  endif
  write(26)ps_i
  close (26)
  write(6,*)'dms emissions ',maxval(f1)
  write(6,*)minval(f1)
  write(6,*)maxval(ps_i)
  write(6,*)minval(ps_i)
!
!  gocart background fields
!
  write(6,*) 'newname:  interpolating gocart_backgd_littlee ...'
  open(unit=29,file='gocart_backgd_littlee',convert='little_endian',form='unformatted',status='old',iostat=ios)
  if (ios /= 0) then
    write(6,*) 'ERROR:  failed to open file gocart_backgd_littlee'
    stop
  endif
   write(6,*) 'opened file gocart_backgd '
   read(29)p_gocart
   write(6,*)p_gocart
   read(29)tmp_oh
   write(6,*)'read oh'
   read(29)tmp_h2o2
   write(6,*)'read h2o2'
   read(29)tmp_no3
   write(6,*)'read no3'
   close(29)
!
! loop over levels
!
  do nv=1,klev
  print *,'read level ',nv
  minv=minval(tmp_oh(:,:,nv))
  maxv=maxval(tmp_oh(:,:,nv))
  print *,'minv,maxv = ',minv,maxv
  CALL bl_int (tmp_oh(:,:,nv), ps_i)  ! unit in pascal
  write(6,*)'after interpolation ',minval(ps_i),maxval(ps_i)
  oh(nv,:)=ps_i(:)
  enddo
!
! end level loop, now write
!
  g3 = "oh.dat"
  open(unit=26,file=g3,form='unformatted',iostat=ios)
  if (ios /= 0) then
    write(6,*) 'ERROR:  failed to open file ',TRIM(g3)
    stop
  endif
  write(26)p_gocart
  write(26)oh
  close (26)
!
! next variable
!
!
! loop over levels
!
  do nv=1,klev
  print *,'read level ',nv
  minv=minval(tmp_h2o2(:,:,nv))
  maxv=maxval(tmp_h2o2(:,:,nv))
  print *,'minv,maxv = ',minv,maxv
  CALL bl_int (tmp_h2o2(:,:,nv), ps_i)  ! unit in pascal
  write(6,*)'after interpolation ',minval(ps_i),maxval(ps_i)
  oh(nv,:)=ps_i(:)
  enddo
!
! end level loop, now write
!
  g3 = "h2o2.dat"
  open(unit=26,file=g3,form='unformatted',iostat=ios)
  if (ios /= 0) then
    write(6,*) 'ERROR:  failed to open file ',TRIM(g3)
    stop
  endif
  write(26)p_gocart
  write(26)oh
  close (26)
!
! next variable
!
!
! loop over levels
!
  do nv=1,klev
  print *,'read level ',nv
  minv=minval(tmp_no3(:,:,nv))
  maxv=maxval(tmp_no3(:,:,nv))
  print *,'minv,maxv = ',minv,maxv
  CALL bl_int (tmp_no3(:,:,nv), ps_i)  ! unit in pascal
  write(6,*)'after interpolation ',minval(ps_i),maxval(ps_i)
  oh(nv,:)=ps_i(:)
  enddo
!
! end level loop, now write
!
  g3 = "no3.dat"
  open(unit=26,file=g3,form='unformatted',iostat=ios)
  if (ios /= 0) then
    write(6,*) 'ERROR:  failed to open file ',TRIM(g3)
    stop
  endif
  write(26)p_gocart
  write(26)oh
  close (26)
!
!  anhropogenic emissions
!
  write(6,*) 'newname:  interpolating anthro_binary ...'
  open(unit=29,file='anthro_binary',convert='big_endian',form='unformatted',status='old',iostat=ios)
  if (ios /= 0) then
    write(6,*) 'ERROR:  failed to open file anthro_binary'
    stop
  endif
   write(6,*) 'opened file anthro'
  read(29)nv_g
  print *,nv_g
  read(29)dname
  print *,dname
  read(29)itime
  print *,itime
  do nv=1,nspecies
  read(29)f4 !emissions
  minv=minval(f4)
  maxv=maxval(f4)
  write(6,*)'read max,min = ',maxv,minv
! f1(1:imax*jmax)=emissions(1:imax*jmax,nv)
  CALL bl_int (f4, ps_i)  ! unit in pascal
  g3 = TRIM(ename(nv)) // ".dat"
  write(6,*)'species ',nv,'file = ',TRIM(g3),' emissions ',TRIM(ename(nv))
  open(unit=26,file=g3,form='unformatted',iostat=ios)
  if (ios /= 0) then
    write(6,*) 'ERROR:  failed to open file ',TRIM(g3)
    stop
  endif
  write(26)ps_i
  close (26)
! write(6,*)ename(nv), ' emissions ',maxval(f1)
! write(6,*)minval(f1)
  write(6,*)maxval(ps_i)
  write(6,*)minval(ps_i)
  enddo
  close (29)
!
!   horizontal interp of chemistry background for gocart
!

! deallocate arrays
  DEALLOCATE(ps_i)

  end program
