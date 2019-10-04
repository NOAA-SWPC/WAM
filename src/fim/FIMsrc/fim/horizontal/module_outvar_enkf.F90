module module_outvar_enkf

contains

subroutine outvar_enkf(time,pr3d,ex3d,us3d,vs3d,tr3d,ph3d)
  use module_constants       ,only: lat,lon,rd,cp,p1000,grvity
  use module_control         ,only: &
  nip,nvl,nvlp1,ntra,FixedGridOrder,dt,ptop,glvl,curve
  use module_savesfc, only : savesfc
  implicit none
  ! External variable declarations:
!SMS$DISTRIBUTE (dh,nip) BEGIN
  real   ,intent(IN)    :: us3d(nvl  ,nip),vs3d(nvl,nip)
  real   ,intent(IN)    :: pr3d(nvlp1,nip)
  real   ,intent(IN)    :: ph3d(nvlp1,nip), ex3d(nvlp1,nip)
  real   ,intent(IN)    :: tr3d(nvl,nip,ntra)
!SMS$DISTRIBUTE END
  real                  :: tv3d(nvl  ,nip),td3d(nvl,nip) ! local vars
  character(len=6)::timestr
  ! local vars
  real exn
  integer ivl,ipn, time, idx,lunout

  ! compute virtual temperature, layer pressure.
!SMS$SERIAL (<ex3d,pr3d,tr3d,IN>,<td3d,tv3d,OUT> : default=ignore) BEGIN
  do ipn=1,nip
  do ivl=1,nvl
      if (pr3d(ivl,ipn).gt.pr3d(ivl+1,ipn)+0.1) then
        exn = (ex3d(ivl  ,ipn)*pr3d(ivl  ,ipn)                     &
                   -ex3d(ivl+1,ipn)*pr3d(ivl+1,ipn))/                   &
                   ((cp+rd)*(pr3d(ivl,ipn)-pr3d(ivl+1,ipn)))
      else
        exn = .5*(ex3d(ivl,ipn)+ex3d(ivl+1,ipn))/cp
      end if
      !exn = .5*(ex3d(ivl,ipn)+ex3d(ivl+1,ipn))/cp
      td3d(ivl,ipn)=p1000*(exn)**(cp/rd) ! layer pressure
      tv3d(ivl,ipn) = tr3d(ivl,ipn,1) * exn                    ! virtual temp
      !tv3d(ivl,ipn)=tv3d(ivl,ipn)/ (1.+0.6078*tr3d(ivl,ipn,1)) ! temperature
  enddo
  enddo
!SMS$SERIAL END

  ! write only fields for EnKF DA cycle to a single file.
  lunout = 77

!SMS$SERIAL BEGIN
  print *,'writing out EnKFIO file'
  ! save model state.
  write (timestr,'(i6.6)') time
  open (lunout,file='fim_out_'//timestr,form="unformatted")
  write(lunout) nip,nvl,3,0.01*ptop
  write(lunout) lon
  write(lunout) lat
  ! orography.
  write(lunout) ph3d(1,:)/grvity
  ! pressure (hPa) on model layer midpoints.
  do ivl=1,nvl
    write(lunout) 0.01*td3d(ivl,:)
  enddo
  ! pressure (hPa) on model layer interaces (including
  ! surface pressure (k=1) but not model top (k=nlevs+1)).  
  ! Model top pressure is assumed constant = ptop.
  ! pressure (hPa) on model layer midpoints.
  do ivl=1,nvl
    write(lunout) 0.01*pr3d(ivl,:)
  enddo
  ! virtual temperature.
  do ivl=1,nvl
    write(lunout) tv3d(ivl,:)
  enddo
  ! u and v winds.
  do ivl=1,nvl
    write(lunout) us3d(ivl,:)
  enddo
  do ivl=1,nvl
    write(lunout) vs3d(ivl,:)
  enddo
  ! "tracers" (vapor, ozone, cloud condensate)
  ! note: order is different in GFS/FIM (ozone comes before clw in GFS)
  do ivl=1,nvl
    write(lunout) tr3d(ivl,:,2)
  enddo
  do ivl=1,nvl
    write(lunout) tr3d(ivl,:,4)
  enddo
  do ivl=1,nvl
    write(lunout) tr3d(ivl,:,3)
  enddo
  close(lunout)
!SMS$SERIAL END
  ! save surface data.
  call savesfc('fimsfc_out_'//timestr)
  return
end subroutine outvar_enkf

end module module_outvar_enkf
