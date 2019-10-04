subroutine readenkfanal(nvp,bkgFile,anlFile,bkgFileSig,us3d,vs3d,dp3d,mp3d,pr3d,ex3d,ph3d,tr3d)

! read data produced by EnKF analysis 

  use module_control,only: glvl,nvl,nvlp1,nip,ntra,ntrb,curve,NumCacheBLocksPerPE,&
                           PrintIpnDiag,PrintDiagProgVars,PrintDiagNoise,PrintDiags,ptop,pure_sig
  use module_constants,only:sigak,sigbk
  implicit none

  integer          ,intent(IN)     :: nvp
  CHARACTER(len=80),intent(IN)     :: anlFile,bkgFile,bkgFileSig

!SMS$DISTRIBUTE(dh,NIP) BEGIN
  real,intent(out) :: us3d(nvl  ,nip)	 ! zonal wind (m/s), layer
  real,intent(out) :: vs3d(nvl  ,nip)	 ! meridional wind (m/s), layer
  real,intent(out) :: dp3d(nvl  ,nip)	 ! del p between coord levels (pascals)
  real,intent(out) :: mp3d(nvl  ,nip)	 ! Montgomery Potential (m^2/s^2)
  real,intent(out) :: pr3d(nvlp1,nip)	 ! pressure (pascal)
  real,intent(out) :: ex3d(nvlp1,nip)	 ! exner function
  real,intent(out) :: ph3d(nvlp1,nip)	 ! geopotential (=gz), m^2/s^2
  real,intent(out) :: tr3d(nvl  ,nip,ntra+ntrb)! 1=pot.temp, 2=water vapor, 3=cloud condensate, 4=ozone
  real(4)          :: hs_i(      nip)    ! geopotential (g*zs) at surface
  real(4)          :: ps_i(      nip)    ! surface pressure in pascals
  real(4)          :: t_i (nvp  ,nip)    ! temperature in Kelvins
  real(4)          :: qv_i(nvp  ,nip)    ! specific humidity
  real(4)          :: qc_i(nvp  ,nip)    ! cloud condensate
  real(4)          :: u_i (nvp  ,nip)    ! zonal velocity
  real(4)          :: v_i (nvp  ,nip)    ! meridional velocity
  real(4)          :: o3_i(nvp  ,nip)    ! ozone mixing ratio
  real(4)          :: p_i (nvp+1,nip)    ! interface pressure
  real(4)          :: z_i (nvp+1,nip)    ! geopotential height 
  integer :: k
!SMS$DISTRIBUTE END

! don't use 82, since it's assumed to be big endian
  integer,parameter :: luanl=83 

  allocate(sigak(nvlp1),sigbk(nvlp1))
  sigak(1:65) =                                         		&
  (/  0.000000,     0.000000,     0.575000,     5.741000,    21.516001, &
     55.712002,   116.899002,   214.014999,   356.222992,   552.719971, &
    812.489014,  1143.988037,  1554.788940,  2051.149902,  2637.552979, &
   3316.217041,  4086.614014,  4945.028809,  5884.206055,  6893.117188, &
   7956.908203,  9057.050781, 10171.711914, 11276.347656, 12344.490234, &
  13348.670898, 14261.434570, 15056.341797, 15708.892578, 16197.315430, &
  16503.144531, 16611.603516, 16511.736328, 16197.966797, 15683.489258, &
  14993.074219, 14154.316406, 13197.065430, 12152.936523, 11054.852539, &
   9936.614258,  8832.537109,  7777.149902,  6804.874023,  5937.049805, &
   5167.145996,  4485.493164,  3883.052002,  3351.459961,  2883.038086, &
   2470.788086,  2108.365967,  1790.051025,  1510.711060,  1265.751953, &
   1051.079956,   863.057983,   698.456970,   554.424011,   428.433990, &
    318.265991,   221.957993,   137.789993,    64.247002,     0.000000 /)

  sigbk(1:65) =                                         		&
(/ 1.000000000,  0.994671166,  0.988626599,  0.981742263,  0.973867595, &
   0.964827597,  0.954434097,  0.942491055,  0.928797305,  0.913151026, &
   0.895354986,  0.875223577,  0.852590680,  0.827318847,  0.799309731, &
   0.768514693,  0.734945238,  0.698682904,  0.659887016,  0.618799627, &
   0.575746655,  0.531134844,  0.485443324,  0.439210802,  0.393018246, &
   0.347468495,  0.303164124,  0.260685444,  0.220570192,  0.183296233, &
   0.149268776,  0.118812189,  0.092166908,  0.069474578,  0.050646842, &
   0.035441618,  0.023555880,  0.014637120,  0.008294020,  0.004106710, &
   0.001635910,  0.000431060,  0.000036970,  0.000000000,  0.000000000, &
   0.000000000,  0.000000000,  0.000000000,  0.000000000,  0.000000000, &
   0.000000000,  0.000000000,  0.000000000,  0.000000000,  0.000000000, &
   0.000000000,  0.000000000,  0.000000000,  0.000000000,  0.000000000, &
   0.000000000,  0.000000000,  0.000000000,  0.000000000,  0.000000000 /)

  sigak(65)=ptop
!SMS$SERIAL (<hs_i,ps_i,z_i,p_i,t_i,qv_i,u_i,v_i,o3_i,qc_i,OUT> : default=ignore)  BEGIN
  if ( pure_sig) then
  call readgriddata_sig(luanl,        anlFile,           hs_i,ps_i,p_i,z_i,t_i,u_i,v_i,qv_i,o3_i,qc_i,nvp,nip) 
  else
  call readgriddata(    luanl,bkgFile,anlFile,bkgFileSig,hs_i,ps_i,p_i,z_i,t_i,u_i,v_i,qv_i,o3_i,qc_i,nvp,nip)
  endif
  
  print*,'returned from readgriddata'
  print*,hs_i(100),ps_i(100),nvp,nip
  DO k=1,nvp,6
     print*,p_i(k,100),z_i(k,100),t_i(k,100)
  ENDDO
!SMS$SERIAL END
! Now do vertical interpolation.
  call fimini(nvp,hs_i,ps_i,z_i,p_i,t_i,qv_i,u_i,v_i,o3_i,qc_i,	&
              us3d,vs3d,dp3d,mp3d,pr3d,ex3d,ph3d,tr3d)
  print*,'returned from fimini'
  DO k=1,nvp,6
     print*,pr3d(k,100),ex3d(k,100),tr3d(k,100,1)
  ENDDO


  return 
 
end subroutine readenkfanal

 subroutine readgriddata(iunit,filename,filename1,filename2,topo,psg,pr3d,ph3d,tempg,u3d,v3d,qg1,qg2,qg3,nlevs,npts)
 
  use module_constants,only: rd,cp,grvity,p1000,qvmin,sigak,sigbk
  use module_control  ,only: ptop,PrintDiags,PrintIpnDiag
  use physcons        ,only: con_fvirt
  implicit none
  integer,                       intent(in)  :: nlevs,npts,iunit
  character,                     intent(in)  :: filename*80,filename1*80,filename2*80
  real, dimension(npts),         intent(out) :: topo,psg
  real, dimension(nlevs,npts),   intent(out) :: tempg,u3d,v3d,qg1,qg2,qg3
  real, dimension(nlevs+1,npts), intent(out) :: pr3d,ph3d
  
! locals  
  real, dimension(npts) :: lons,lats,tmp1,tmp2
  real, dimension(nlevs,npts) :: pslg
  real, dimension(nlevs,npts) :: inc,tmp3d,sig
  real, dimension(nlevs+1,npts) :: pr3da,ex3d
  real :: exn,ptopin
  integer i,k,nlevsin,ntracin,nptsin,ierr,iunit1,iunit2,ii
  real, dimension(nlevs,npts)     :: uswrk,vswrk,thwrk,qvwrk,o3wrk,qcwrk
  real, dimension(nlevs+1,npts)     :: exwrk
  real :: theta_lyrs(nlevs)
  real :: targ_in(nlevs,npts)
  real :: th_in (nlevs,npts)

!SMS$ignore begin
  iunit1=44
  iunit2=45
  open(iunit,file=trim(filename),form="unformatted")
  open(iunit1,file=trim(filename1),form="unformatted")
  open(iunit2,file=trim(filename2),form="unformatted")
  read(iunit,iostat=ierr) nptsin,nlevsin,ntracin,ptopin
  read(iunit1,iostat=ierr) nptsin,nlevsin,ntracin,ptopin
  read(iunit2,iostat=ierr) nptsin,nlevsin,ntracin,ptopin
  if (npts .ne. nptsin .or. nlevs .ne. nlevsin .or. ntracin .ne. 3) then
    print *,'error reading input file - npts,nlevs,ntrac !=',&
             npts,nlevs,3,nptsin,nlevsin,ntracin
    stop
  end if
  print *,'ptop = ',ptop
  ! read lons, lats on model grid (radians) - not used here.
  read(iunit) lons
  read(iunit1) 
  read(iunit2) 
  read(iunit) lats
  read(iunit1) 
  read(iunit2) 
  ! read surface orography.
  read(iunit) topo
  read(iunit1) 
  read(iunit2) 
  print *,'min/max topo',minval(topo),maxval(topo)
  ! read pressure (hPa) on model layer midpoints.
  do k=1,nlevs
     read(iunit)  pslg(k,:)
     read(iunit1) 
     read(iunit2) 
  enddo
  print *,'min/max pslg',minval(pslg),maxval(pslg)
  ! read pressure (hPa) on model layer interaces (including
  ! surface pressure (k=1) but not model top (k=nlevs+1)).  
  ! Model top pressure is assumed constant = ptop.
  do k=1,nlevs
     read(iunit)  pr3d(k,:)
     read(iunit1) pr3da(k,:)
     read(iunit2)
  enddo
! calculate sigma values for pr3d,  update surface pressure from analysis, then re-calculate pr3d for rest of interfaces
  ii=100
  DO k=2,nlevs
      sig(k,:)=pr3d(k,:)/pr3d(1,:)
  ENDDO
  print*,'pr3d  at ',ii,'=',pr3d(:,ii)
  print*,'pr3da at ',ii,'=',pr3da(:,ii)
  pr3d(1,:)=pr3da(1,:)
  DO k=2,nlevs
      pr3d(k,:)=sig(k,:)*pr3d(1,:)
  ENDDO
  print*,'pr3d  at ',ii,'=',pr3d(:,ii)
!  surface pressure does not need to be interpolated
  pr3d  = 1.0e2*pr3d
  pr3da = 1.0e2*pr3da
  pr3d(nlevs+1,:)=ptop
  pr3da(nlevs+1,:)=ptopin
  !reset pr3da to be on the gfs coordimate
  do k=2,nlevs+1
      pr3da(k,:) =sigak(k)+sigbk(k)*pr3da(1,:)
  enddo
  ! re-diagnose layer pressures from interface pressures
  ! using FIM's algorithm (from output.F90).
  ex3d = cp*(pr3d/p1000)**(rd/cp)
  do k=1,nlevs
     do i=1,npts
        if (pr3d(k,i).gt.pr3d(k+1,i)+0.1) then
           exn = (ex3d(k  ,i)*pr3d(k  ,i)                     &
                   -ex3d(k+1,i)*pr3d(k+1,i))/              &
                   ((cp+rd)*(pr3d(k,i)-pr3d(k+1,i)))
        else
           exn = .5*(ex3d(k,i)+ex3d(k+1,i))/cp
        end if
        pslg(k,i)=p1000*(exn)**(cp/rd) ! layer pressure
     enddo
  enddo
  ! put pr3da back to mb for interpolation
  print *,'min/max pr3d',minval(pr3d),maxval(pr3d)
  print *,'min/max psf',minval(pr3d(1,:)),maxval(pr3d(1,:))
  DO k=1,nlevs
      print*,pr3d(k,ii),pr3da(k,ii),pslg(k,ii) 
  ENDDO
  ! read virtual temperature.
  do k=1,nlevs
      read(iunit) tempg(k,:)
      read(iunit1) tmp1
      read(iunit2) tmp2
      inc(k,:)=tmp1-tmp2
      print*,'temp inc =',k,inc(k,ii),tmp1(ii),tmp2(ii),tempg(k,ii)
  enddo
  call vlint2coor(npts, nlevs, nlevs+1, pr3da, inc,tmp3d, pslg, nlevs)
  tempg(:,:)   =     tmp3d(:,:) + tempg(:,:)          ! virt.pot.temperature
  do k=1,nlevs
      print*,'temp inc =',k,inc(k,ii),tmp3d(k,ii),tempg(k,ii)
  enddo
  print *,'min/max tempg',minval(tempg),maxval(tempg)
  ! read u and v winds.
  do k=1,nlevs
     read(iunit) u3d(k,:)
     read(iunit1) tmp1
     read(iunit2) tmp2
     inc(k,:)=tmp1-tmp2
  enddo
  call vlint2coor(npts, nlevs, nlevs+1, pr3da, inc, tmp3d, pslg, nlevs)
  u3d(:,:)    =     tmp3d(:,:) + u3d(:,:)
  do k=1,nlevs
     read(iunit) v3d(k,:)
     read(iunit1) tmp1
     read(iunit2) tmp2
     inc(k,:)=tmp1-tmp2
  enddo
  call vlint2coor(npts, nlevs, nlevs+1, pr3da, inc, tmp3d, pslg, nlevs)
  v3d(:,:)    =     tmp3d(:,:) + v3d(:,:)
  ! read "tracers" (vapor, ozone, cloud condensate)
  do k=1,nlevs
     read(iunit) qg1(k,:)
     read(iunit1) tmp1
     read(iunit2) tmp2
     inc(k,:)=tmp1-tmp2
  enddo
  call vlint2coor(npts, nlevs, nlevs+1, pr3da, inc,tmp3d, pslg, nlevs)
  qg1(:,:)   = max(tmp3d(:,:) + qg1(:,:),qvmin) ! water vapor
  do k=1,nlevs
     read(iunit) qg2(k,:)
     read(iunit1) tmp1
     read(iunit2) tmp2
     inc(k,:)=tmp1-tmp2
  enddo
  call vlint2coor(npts, nlevs, nlevs+1, pr3da, inc,tmp3d, pslg, nlevs)
  qg2(:,:)   = max(tmp3d(:,:) + qg2(:,:),0.)    ! liquid water/condensate
  do k=1,nlevs
     read(iunit) qg3(k,:)
     read(iunit1) tmp1
     read(iunit2) tmp2
     inc(k,:)=tmp1-tmp2
  enddo
  call vlint2coor(npts, nlevs, nlevs+1, pr3da, inc,tmp3d, pslg, nlevs)
  qg3(:,:)   = max(tmp3d(:,:) + qg3(:,:),0.)    ! ozone
  close(iunit) 
  close(iunit1) 
  close(iunit2) 

!SMS$ignore end
  
  topo=topo*grvity
  psg = pr3d(1,:)
  call temptoz(npts,nlevs,rd,cp,p1000,pr3d,pslg,topo,tempg,ph3d)
! convert virtual t to sensible t.
  tempg=tempg/(1.+con_fvirt*qg1(:,:))

 end subroutine readgriddata

 subroutine readgriddata_sig(iunit,filename,topo,psg,pr3d,ph3d,tempg,u3d,v3d,qg1,qg2,qg3,nlevs,npts)
 
  use module_constants,only: rd,cp,grvity,p1000,qvmin,sigak,sigbk
  use module_control  ,only: ptop,PrintDiags,PrintIpnDiag
  use physcons        ,only: con_fvirt
  implicit none
  integer,                       intent(in)  :: nlevs,npts,iunit
  character,                     intent(in)  :: filename*80
  real, dimension(npts),         intent(out) :: topo,psg
  real, dimension(nlevs,npts),   intent(out) :: tempg,u3d,v3d,qg1,qg2,qg3
  real, dimension(nlevs+1,npts), intent(out) :: pr3d,ph3d
  
! locals  
  real, dimension(npts) :: lons,lats,tmp1,tmp2
  real, dimension(nlevs,npts) :: pslg
  real, dimension(nlevs,npts) :: inc,tmp3d,sig
  real, dimension(nlevs+1,npts) :: ex3d
  real :: exn,ptopin
  integer i,k,nlevsin,ntracin,nptsin,ierr,ii
  real, dimension(nlevs,npts)     :: uswrk,vswrk,thwrk,qvwrk,o3wrk,qcwrk
  real, dimension(nlevs+1,npts)     :: exwrk
  real :: theta_lyrs(nlevs)
  real :: targ_in(nlevs,npts)
  real :: th_in (nlevs,npts)

!SMS$ignore begin
  open(iunit,file=trim(filename),form="unformatted")
  read(iunit,iostat=ierr) nptsin,nlevsin,ntracin,ptopin
  if (npts .ne. nptsin .or. nlevs .ne. nlevsin .or. ntracin .ne. 3) then
    print *,'error reading input file - npts,nlevs,ntrac !=',&
             npts,nlevs,3,nptsin,nlevsin,ntracin
    stop
  end if
  print *,'ptop = ',ptop
  ! read lons, lats on model grid (radians) - not used here.
  read(iunit) lons
  read(iunit) lats
  ! read surface orography.
  read(iunit) topo
  print *,'min/max topo',minval(topo),maxval(topo)
  ! read pressure (hPa) on model layer midpoints.
  do k=1,nlevs
     read(iunit)  pslg(k,:)
  enddo
  print *,'min/max pslg',minval(pslg),maxval(pslg)
  ! read pressure (hPa) on model layer interaces (including
  ! surface pressure (k=1) but not model top (k=nlevs+1)).  
  ! Model top pressure is assumed constant = ptop.
  do k=1,nlevs
     read(iunit)  pr3d(k,:)
  enddo
  ii=100
  print*,'before pr3d  at ',ii,'=',pr3d(:,ii) 
  pr3d(1,:)  = 1.0e2*pr3d(1,:)
  !reset pr3d to be on the gfs coordimate
  do k=2,nlevs
      pr3d(k,:) =sigak(k)+sigbk(k)*pr3d(1,:)
  enddo
  pr3d(nlevs+1,:)=ptop
  print*,'after pr3d  at ',ii,'=',pr3d(:,ii)
  ! re-diagnose layer pressures from interface pressures
  ! using FIM's algorithm (from output.F90).
  ex3d = cp*(pr3d/p1000)**(rd/cp)
  do k=1,nlevs
     do i=1,npts
        if (pr3d(k,i).gt.pr3d(k+1,i)+0.1) then
           exn = (ex3d(k  ,i)*pr3d(k  ,i)                     &
                   -ex3d(k+1,i)*pr3d(k+1,i))/              &
                   ((cp+rd)*(pr3d(k,i)-pr3d(k+1,i)))
        else
           exn = .5*(ex3d(k,i)+ex3d(k+1,i))/cp
        end if
        pslg(k,i)=p1000*(exn)**(cp/rd) ! layer pressure
     enddo
  enddo 
  print*,'after pslg  at ',ii,'=',pslg(:,ii)

  ! read virtual temperature.
  do k=1,nlevs
      read(iunit) tempg(k,:)
  enddo
  print *,'min/max tempg',minval(tempg),maxval(tempg)
  ! read u and v winds.
  do k=1,nlevs
     read(iunit) u3d(k,:)
  enddo
  do k=1,nlevs
     read(iunit) v3d(k,:)
  enddo
  ! read "tracers" (vapor, ozone, cloud condensate)
  do k=1,nlevs
     read(iunit) qg1(k,:)
  enddo
  do k=1,nlevs
     read(iunit) qg2(k,:)
  enddo
  do k=1,nlevs
     read(iunit) qg3(k,:)
  enddo
  close(iunit) 

!SMS$ignore end
  
  topo=topo*grvity
  psg = pr3d(1,:)
  call temptoz(npts,nlevs,rd,cp,p1000,pr3d,pslg,topo,tempg,ph3d)
! convert virtual t to sensible t.
  tempg=tempg/(1.+con_fvirt*qg1(:,:))

 end subroutine readgriddata_sig

!=============================================================================
! Level and layer interpolation, target coordinate: passed in (v_coor) 
!
!
! Figure 1.                          | Figure 2.
! Level variables at interface:      | Layer variables: 
!                                    |
! -------------- int. level nvl + 1  | ---------------- int. level nvl + 1 
!       :                            |         :
!       :                            |         :
! -------------- int. level k + 1    | ---------------- int. level k + 1
!                                    | //////////////// variable at layer k
! -------------- int. level k        | ---------------- int. level k
!       :                            |         :
!       :                            |         :
! -------------- int. level 1        | ---------------- int. level 1
!                                    |
! 
!
! N. Wang, Feb. 2008
!=============================================================================

      SUBROUTINE vlint2coor(nip, nvl, nvlp1, data_pr, data_var, data, v_coor, nvc)
        IMPLICIT NONE

        INTEGER, intent(in) ::  nip, nvl, nvlp1, nvc
        REAL, intent(in) ::  data_pr(nvlp1,nip),  data_var(nvl,nip), v_coor(nvl,nip)
        REAL, intent(out) :: data(nvc,nip)

        REAL pi_dn, pi_up, pi_co, dn_val, up_val
        INTEGER k,i, l
        DO i = 1, nip   
            k = 1
            DO l = 1, nvc 
              IF (v_coor(l,i) >= data_pr(1,i)) THEN
                data(l,i) = data_var(1,i)
                CYCLE
              END IF
              IF (v_coor(l,i) <= data_pr(nvlp1,i)) THEN
                data(l,i) = data_var(nvl,i)
                CYCLE
              END IF
              DO WHILE (v_coor(l,i) < data_pr(k+1,i)) 
                IF (k == nvlp1 - 1) THEN
                  EXIT 
                ELSE
                  k = k + 1
                ENDIF
              END DO ! k and k+1 are the current indexes for interpolation  
              IF (nvl == nvlp1) THEN ! level variables, see fig. 1.
                pi_dn = (data_pr(k,i) / 1000.00)**0.286
                pi_up = (data_pr(k+1,i) / 1000.00)**0.286 
                pi_co = (v_coor(l,i) / 1000.00)**0.286
                dn_val = data_var(k,i)
                up_val = data_var(k+1,i)
              ELSE  ! layer variables, see fig. 2.
                pi_dn = (data_pr(k,i) / 1000.00)**0.286
                pi_up = (data_pr(k+1,i) / 1000.00)**0.286 
                pi_co = (v_coor(l,i) / 1000.00)**0.286
                IF (pi_co > (pi_dn + pi_up) / 2.0) THEN ! lower half of the layer
                  IF (k == 1) THEN
                    pi_dn = (data_pr(k,i) / 1000.00)**0.286
                    pi_up = ((data_pr(k,i) / 1000.00)**0.286 + (data_pr(k+1,i) / 1000.00)**0.286)/ 2.0
                    pi_co = (v_coor(l,i) / 1000.00)**0.286
                    dn_val = data_var(k,i)
                    up_val = dn_val 
                  ELSE
                    pi_dn = ((data_pr(k,i) / 1000.00)**0.286 + (data_pr(k-1,i) / 1000.00)**0.286)/ 2.0 
                    pi_up = ((data_pr(k,i) / 1000.00)**0.286 + (data_pr(k+1,i) / 1000.00)**0.286)/ 2.0 
                    pi_co = (v_coor(l,i) / 1000.00)**0.286
                    dn_val = data_var(k-1,i)
                    up_val = data_var(k,i)
                  ENDIF
                ELSE ! upper half of the layer
                  IF (k == nvl) THEN
                    pi_dn = ((data_pr(k,i) / 1000.00)**0.286 + (data_pr(k+1,i) / 1000.00)**0.286)/ 2.0 
                    pi_up = (data_pr(k+1,i) / 1000.00)**0.286
                    pi_co = (v_coor(l,i) / 1000.00)**0.286
                    dn_val = data_var(k,i)
                    up_val = dn_val 
                  ELSE
                    pi_dn = ((data_pr(k,i) / 1000.00)**0.286 + (data_pr(k+1,i) / 1000.00)**0.286)/ 2.0
                    pi_up = ((data_pr(k+1,i) / 1000.00)**0.286 + (data_pr(k+2,i) / 1000.00)**0.286)/ 2.0 
                    pi_co = (v_coor(l,i) / 1000.00)**0.286
                    dn_val = data_var(k,i)
                    up_val = data_var(k+1,i)
                  ENDIF

                ENDIF
              ENDIF
              data(l,i) = up_val +  (pi_co - pi_up) * &
                (dn_val - up_val) / (pi_dn - pi_up) 
            END DO
           IF (i.EQ.100) THEN
                print*,'in vlint'
                print*,'data_pr=',data_pr(:,i)
                print*,'data_var=',data_var(:,i)
                print*,'data=',data(:,i)
                print*,'v_coor=',v_coor(:,i)
           ENDIF
        END DO
      END SUBROUTINE vlint2coor
 subroutine temptoz(npts,nlevs,rgas,cp,p1000,pint,pl,zs,tv,z)
  implicit none
  integer, intent(in) :: npts,nlevs
  real, dimension(nlevs, npts) :: thetav,pil
  real, dimension(nlevs+1,npts) :: pii
  real, intent(in), dimension(nlevs,npts) :: tv,pl
  real, intent(in), dimension(nlevs+1,npts) :: pint
  real, intent(out), dimension(nlevs+1,npts) :: z
  real, intent(in), dimension(npts) :: zs
  real, intent(in) :: rgas,cp,p1000
  integer j,k

  pii = cp*(pint/p1000)**(rgas/cp)
  pil = cp*(pl/p1000)**(rgas/cp)
  thetav = cp*tv/pil
  do j=1,npts
     z(1,j) = zs(j)
     do k=1,nlevs
        z(k+1,j) = z(k,j)  - thetav(k,j) * (pii(k+1,j)-pii(k,j))
     end do
  end do

 end subroutine temptoz
