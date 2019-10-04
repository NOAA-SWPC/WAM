!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        module module_fltbnds
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use mpi

use module_kinds

use module_control,only : klog,kint,kfpt 
 
use module_my_domain_specs

use module_derived_types,only: bc_h_all,bc_v_all,filt_4d

use module_dm_parallel,only : gather_layers,scatter_layers

use module_exchange, only: halo_exch

public :: presmud

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
real(kind=kfpt),parameter:: &
 pi=3.141592653589793238462643383279502884197169399375105820 &
,pih=pi/2. &
,tpi=pi+pi &
,dtr=pi/180.
!
integer(kind=kint) :: &
 ipe_start_north &
,ipe_start_south &
,ipe_end_north &
,ipe_end_south &
,jh_start_fft_north &
,jh_end_fft_north &
,jh_start_fft_south &
,jh_end_fft_south &
,jv_start_fft_north &
,jv_end_fft_north &
,jv_start_fft_south &
,jv_end_fft_south &
,lm_fft &   ! Max number of model layers per task for FFTs
,msize_dummy_fft 
!
integer(kind=kint),dimension(mpi_status_size),private :: &
 istatw
!
integer(kind=kint),allocatable,dimension(:) :: &
 k1_fft &
,k2_fft &
,my_jrow_start_h &
,my_jrow_start_v &
,my_jrow_end_h &
,my_jrow_end_v
!
real(kind=kfpt),allocatable,dimension(:,:,:) :: &
 hn &
,un &
,vn &
,wn
!
logical(kind=klog) :: &
 fft_north &
,fft_south 
!
logical(kind=klog),allocatable,dimension(:) :: &
 my_domain_has_fft_lats_h &
,my_domain_has_fft_lats_v
!
!-----------------------------------------------------------------------
      contains
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
                        subroutine prefft &
(dlmd,dphd,sbd,lm &
,khfilt,kvfilt &
,hfilt,vfilt &
,wfftrh,nfftrh,wfftrw,nfftrw &
,inpes,jnpes,mype)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------

real(kind=kfpt),parameter :: &
 cxnc=0.0 &
,rwind=1./2.6 &
,cfilt=1.

integer(kind=kint),intent(in) :: &
 inpes &
,jnpes &
,mype &
,lm

integer(kind=kint),dimension(jds:jde),intent(out) :: &
 khfilt &
,kvfilt

integer(kind=kint),dimension(1:15),intent(out) :: &
 nfftrh &
,nfftrw 

real(kind=kfpt),intent(in) :: &
 dlmd &
,dphd &
,sbd
 
real(kind=kfpt),dimension(ids:ide,jds:jde),intent(out) :: &
 hfilt &
,vfilt

real(kind=kfpt),dimension(1:2*(ide-3)),intent(out) :: &
 wfftrh &
,wfftrw

!-----------------------------------------------------------------------
!***  Local Variables
!-----------------------------------------------------------------------

logical(kind=klog) :: &
 first

integer(kind=kint) :: &
 icycle &
,imax &
,istat &
,j &
,jmax &
,jmax_north &
,jmax_south &
,k &
,k2 &
,kount_layers &
,kount_pes &
,ks &
,l_remain &
,lyr_frac_north &
,lyr_frac_south &
,n &
,n_extra &
,n_factor &
,n_group1 &
,n_group2 &
,n_remain &
,n_remainder_h_group1 &
,n_remainder_h_group2 &
,n_remainder_v_group1 &
,n_remainder_v_group2 &
,nnew &
,npe &
,npe_next &
,npes &
,npes_north &
,npes_south &
,nrow_x &
,nrows_fft_north_h &
,nrows_fft_south_h &
,nrows_fft_north_v &
,nrows_fft_south_v &
,nrows_group1_h &
,nrows_group2_h &
,nrows_group1_v &
,nrows_group2_v &
,nsmud

real(kind=kfpt) :: &
 cpf &
,cph &
,cx &
,cxn &
,dlm &
,dph &
,flt &
,rcph &   
,rcycle &
,sb &
,sub_j &
,sxl &
,tph &
,x &
,xl

!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
      icycle=ide-3
      rcycle=1./icycle
!                                                      
      dlm=dlmd*dtr                                                      
      dph=dphd*dtr
      sb=sbd*dtr
!
      cpf=dph/dlm
!-----------------------------------------------------------------------
!***  Maximum number of layers used per task for FFTs.                 
!***  Since the FFT regions cap the poles, there are 2*LM layers
!***  in play.
!-----------------------------------------------------------------------
!
      npes=inpes*jnpes     
!
!-----------------------------------------------------------------------
!-------------filters on h and v rows-----------------------------------
!-----------------------------------------------------------------------
!
      jh_start_fft_south=-10
      jh_end_fft_south=-10
      jh_start_fft_north=-10
      jh_end_fft_north=-10
!
      sub_j=real(jds+1)
!
!-----------------------------------------------------------------------
!
      h_filters: do j=jds+2,jde-2
!
!-----------------------------------------------------------------------
        tph=sb+(j-sub_j)*dph
        cph=cos(tph)
        rcph=(cph/cpf)
        ks=icycle
        nsmud=0
!-----------------------------------------------------------------------
!***  Perform Fourier filtering where cos(phi) < dph/dlm .
!-----------------------------------------------------------------------
!
        if(rcph.lt.cfilt) then
!
          if(tph<0.)then
            if(jh_start_fft_south<-5)jh_start_fft_south=j
            jh_end_fft_south=j
          else
            if(jh_start_fft_north<-5)jh_start_fft_north=j
            jh_end_fft_north=j
          endif
!
          khfilt(j)=icycle
          hfilt(1,j)=rcycle
          first=.true.
!
          do k=2,icycle-1,2
            x=k*dlm*0.25
            xl=min(x/rcph*cfilt,pih)
            sxl=sin(xl)+sin(2.*xl)*0.5*rwind
!
            if(rcph/sxl.gt.cfilt) then
              flt=rcycle
              ks=k
            else
!--filter definition----------------------------------------------------
                cx=cos(x)
                if(abs(cx).lt.1.e-7) first=.false.
                if(first) then
                  nnew=(log(sxl/(xl*(1.+rwind)))/log(cx)+0.5)
                  if(nnew.gt.nsmud) nsmud=nnew
                  if(mod(nsmud,2).gt.0) nsmud=nsmud+1
                  if(xl.eq.pih) first=.false.
                endif
!-----------------------------------------------------------------------
              cxn=cx**nsmud
!
              if(cxn.lt.cxnc) then
                khfilt(j)=k-1
!
                do k2=k+1,icycle+1
                  hfilt(k2,j)=0.
                enddo
!
                cycle h_filters
              endif
!
              flt=min(cxn*rcycle,rcycle)
            endif
!
            hfilt(k  ,j)=flt
            hfilt(k+1,j)=flt
!
          enddo
!
         
          x=icycle*dlm*0.25
          xl=min(x/rcph*cfilt,pih)
          sxl=sin(xl)+sin(2.*xl)*0.5*rwind
!
          if(rcph/sxl.gt.cfilt) then
            flt=rcycle
            ks=k
          else
!--filter definition----------------------------------------------------
                cx=cos(x)
                if(abs(cx).lt.1.e-7) first=.false.
                if(first) then
                  nnew=(log(sxl/(xl*(1.+rwind)))/log(cx)+0.5)
                  if(nnew.gt.nsmud) nsmud=nnew
                  if(mod(nsmud,2).gt.0) nsmud=nsmud+1
                  if(xl.eq.pih) first=.false.
                endif
!-----------------------------------------------------------------------
              cxn=cx**nsmud
!
            if(cxn.gt.cxnc) then
              khfilt(j)=icycle
            else
              cxn=0.
            endif
!
            flt=min(cxn*rcycle,rcycle)
          endif
!
          hfilt(icycle,j)=flt
        else
          khfilt(j)=icycle+1
          do k=1,icycle
            hfilt(k  ,j)=rcycle
          enddo
        endif
!-----------------------------------------------------------------------
!
      enddo h_filters
!
!-----------------------------------------------------------------------
!
      if(ks.gt.khfilt(j)) ks=0
!
      khfilt(jds)=khfilt(jds+2)
      khfilt(jds+1)=1
      khfilt(jde-1)=1
      khfilt(jde  )=khfilt(jde-2)
!
      do k=1,icycle
        hfilt(k,jds)=hfilt(k,jds+2)
        hfilt(k,jds+1)=0.
        hfilt(k,jde-1)=0.
        hfilt(k,jde  )=hfilt(k,jde-2)
      enddo
!
      hfilt(1,jds+1)=1.
      hfilt(1,jde-1)=1.
!-----------------------------------------------------------------------
!
      jv_start_fft_south=-10
      jv_end_fft_south=-10
      jv_start_fft_north=-10
      jv_end_fft_north=-10
!
      sub_j=jds+0.5
!
!-----------------------------------------------------------------------
!
      v_filters: do j=jds+1,jde-2
!
!-----------------------------------------------------------------------
        tph=sb+(j-sub_j)*dph
        cph=cos(tph)
        rcph=(cph/cpf)
        ks=icycle
        nsmud=0
!-----------------------------------------------------------------------
!***  Perform Fourier filtering where cos(phi) < dph/dlm .
!-----------------------------------------------------------------------
        if(rcph.lt.cfilt) then
!
          if(tph<0.)then
            if(jv_start_fft_south<-5)jv_start_fft_south=j
            jv_end_fft_south=j
          else
            if(jv_start_fft_north<-5)jv_start_fft_north=j
            jv_end_fft_north=j
          endif
!
          kvfilt(j)=icycle
          vfilt(1,j)=rcycle
          first=.true.
!
          do k=2,icycle-1,2
            x=k*dlm*0.25
            xl=min(x/rcph*cfilt,pih)
            sxl=sin(xl)+sin(2.*xl)*0.5*rwind
!
            if(rcph/sxl.gt.cfilt) then
              flt=rcycle
              ks=k
            else
!--filter definition----------------------------------------------------
                cx=cos(x)
                if(abs(cx).lt.1.e-7) first=.false.
                if(first) then
                  nnew=(log(sxl/(xl*(1.+rwind)))/log(cx)+0.5)
                  if(nnew.gt.nsmud) nsmud=nnew
                  if(mod(nsmud,2).gt.0) nsmud=nsmud+1
                  if(xl.eq.pih) first=.false.
                endif
!-----------------------------------------------------------------------
              cxn=cx**nsmud
!
              if(cxn.lt.cxnc) then
                kvfilt(j)=k-1
!
                do k2=k+1,icycle+1
                  vfilt(k2,j)=0.
                enddo
!
                cycle v_filters
              endif
!
              flt=min(cxn*rcycle,rcycle)
            endif
!
            vfilt(k  ,j)=flt
            vfilt(k+1,j)=flt
!
          enddo
!
          x=icycle*dlm*0.25
          xl=min(x/rcph*cfilt,pih)
          sxl=sin(xl)+sin(2.*xl)*0.5*rwind
!
          if(rcph/sxl.gt.cfilt) then
            flt=rcycle
            ks=k
          else
!--filter definition----------------------------------------------------
                cx=cos(x)
                if(abs(cx).lt.1.e-7) first=.false.
                if(first) then
                  nnew=(log(sxl/(xl*(1.+rwind)))/log(cx)+0.5)
                  if(nnew.gt.nsmud) nsmud=nnew
                  if(mod(nsmud,2).gt.0) nsmud=nsmud+1
                  if(xl.eq.pih) first=.false.
                endif
!-----------------------------------------------------------------------
              cxn=cx**nsmud
!
            if(cxn.gt.cxnc) then
              kvfilt(j)=icycle
            else
              cxn=0.
            endif
!
            flt=min(cxn*rcycle,rcycle)
          endif
!
          vfilt(icycle,j)=flt
!
        else
          kvfilt(j)=icycle+1
          do k=1,icycle
            vfilt(k  ,j)=rcycle
          enddo
        endif
!-----------------------------------------------------------------------
!
      enddo v_filters
!
!-----------------------------------------------------------------------
!
      if(ks.gt.kvfilt(j)) ks=0
!
      kvfilt(jds)=kvfilt(jds+1)
      kvfilt(jde-1)=kvfilt(jde-2)
!
      do k=1,icycle
        vfilt(k,jds)=vfilt(k,jds+1)
        vfilt(k,jde-1)=vfilt(k,jde-2)
      enddo
!
      call rffti(icycle,wfftrh,nfftrh)
      call rffti(icycle,wfftrw,nfftrw)
!
!-----------------------------------------------------------------------
!***  Much preparation is needed for setting up the sharing of
!***  FFT computations by all compute tasks.  If only a single
!***  compute task has been designated for the domain then
!***  no prep is needed.
!-----------------------------------------------------------------------
!
      if(inpes==1.and.jnpes==1)then
        return
      endif
!
!-----------------------------------------------------------------------
!
      if(jh_start_fft_south==jds+2)jh_start_fft_south=jds+1
      if(jh_end_fft_north==jde-2)jh_end_fft_north=jde-1
!
!-----------------------------------------------------------------------
!***  Identify tasks as handling Northern or Southern Hemipshere
!***  model layers for FFTs.
!-----------------------------------------------------------------------
!
      fft_south=.false.
      fft_north=.false.
      ipe_start_south=0
      ipe_end_north=npes-1
!
      if(jts<jde/2)then
        fft_south=.true.
      else
        fft_north=.true.
      endif
!
      south: do n=0,npes-1
        if(local_jstart(n)<jde/2)then
          ipe_end_south=n
        else
          exit south
        endif
      enddo south
!
      north: do n=npes-1,0,-1
        if(local_jstart(n)<jde/2)then
          ipe_start_north=n+1
          exit north
        endif
      enddo north
!
!-----------------------------------------------------------------------
!***  Does this task's subdomain contain any FFT latitudes?
!-----------------------------------------------------------------------
!
      allocate(my_domain_has_fft_lats_h(0:npes-1),stat=istat)
      allocate(my_domain_has_fft_lats_v(0:npes-1),stat=istat)
!
      do n=0,npes-1
        my_domain_has_fft_lats_h(n)=.false.
        my_domain_has_fft_lats_v(n)=.false.
      enddo
!
      do n=0,npes-1
!
        if(fft_south.and.local_jstart(n)<=jh_end_fft_south.or. &
           fft_north.and.local_jend(n)>=jh_start_fft_north)then
          my_domain_has_fft_lats_h(n)=.true.
        endif
!
        if(fft_south.and.local_jstart(n)<=jv_end_fft_south.or. &
           fft_north.and.local_jend(n)>=jv_start_fft_north)then
          my_domain_has_fft_lats_v(n)=.true.
        endif
!
      enddo
!
!-----------------------------------------------------------------------
!***  Assign layers to each task.
!***  k1_fft and k2_fft are the first and last model layers in
!***  a task's group of layers over which it will apply FFTs.
!***  Groups of model layers will be assigned from top down
!***  in the southern or northern hemisphere and then divided
!***  if there are more than 2*LM MPI tasks being used.
!***  When there are "remainder" layers, give them one at a time
!***  to each task in the row until they are used up.
!
!***  If there are more than LM MPI tasks in a hemisphere then
!***  the layers themselves begin to be divided up to continue
!***  to ensure that all tasks will receive some of the FFT work.
!-----------------------------------------------------------------------
!
      allocate(k1_fft(0:npes-1),stat=istat)
      allocate(k2_fft(0:npes-1),stat=istat)
      allocate(my_jrow_start_h(0:npes-1),stat=istat)
      allocate(my_jrow_end_h(0:npes-1)  ,stat=istat)
      allocate(my_jrow_start_v(0:npes-1),stat=istat)
      allocate(my_jrow_end_v(0:npes-1)  ,stat=istat)
!
!-------------------------
!***  Southern Hemisphere
!-------------------------
!
      npes_south=ipe_end_south-ipe_start_south+1
!
!----------------------------------------------------
!***  The number of tasks in the Southern Hemisphere
!***  does not exceed the number of model layers.
!----------------------------------------------------
!
      limits_south: if(npes_south<=lm)then                  
!
        lyr_frac_south=lm/npes_south
        l_remain=lm-npes_south*lyr_frac_south
!
        k2=0
        do npe=0,ipe_end_south
          k1_fft(npe)=k2+1
          k2=k1_fft(npe)+lyr_frac_south-1
          if(l_remain>0)then
            k2=k2+1
            l_remain=l_remain-1
          endif
          k2_fft(npe)=k2
!
          my_jrow_start_h(npe)=jh_start_fft_south
          my_jrow_end_h(npe)=jh_end_fft_south
          my_jrow_start_v(npe)=jv_start_fft_south
          my_jrow_end_v(npe)=jv_end_fft_south
        enddo
!
!----------------------------------------------------
!***  If there are more tasks in the hemisphere
!***  than there are model layers then divide the
!***  layers into n_factor pieces for tasks in 
!***  n_group1 and divide the remaining layers into
!***  n_factor+1 pieces for tasks in n_group2.
!----------------------------------------------------
!
      else                                   
        lyr_frac_south=0
        n_factor=npes_south/lm                      
        n_remain=npes_south-n_factor*lm
        n_group1=n_factor*(lm-n_remain)                                   !<-- This many tasks get layers divided into n_factor pieces
        n_group2=npes_south-n_group1                                      !<-- This many tasks get layers divided into n_factor+1 pieces
!
!----------------------------------------------------
!***  Divide layers of FFTs into n_factor pieces
!***  for tasks in n_group1.
!***  Divide remaining layers into n_factor+1 pieces
!***  for tasks in n_group2.
!----------------------------------------------------
!
        nrows_fft_south_h=jh_end_fft_south-jh_start_fft_south+1
!
        nrows_group1_h=nrows_fft_south_h/n_factor                          !<-- Each task in group 1 handles this many H lat rows 
        n_remainder_h_group1=nrows_fft_south_h-nrows_group1_h*n_factor     !    or one additional row to take care of remainders.
!
        nrows_group2_h=nrows_fft_south_h/(n_factor+1)                      !<-- Each task in group 2 handles this many H lat rows
        n_remainder_h_group2=nrows_fft_south_h-nrows_group2_h*(n_factor+1) !    or one additional row to take care of remainders.
!
        nrows_fft_south_v=jv_end_fft_south-jv_start_fft_south+1
!
        nrows_group1_v=nrows_fft_south_v/n_factor                          !<-- Each task in group 1 handles this many V lat rows
        n_remainder_v_group1=nrows_fft_south_v-nrows_group1_v*n_factor     !    or one additional row to take care of remainders.
!
        nrows_group2_v=nrows_fft_south_v/(n_factor+1)                      !<-- Each task in group 2 handles this many V lat rows
        n_remainder_v_group2=nrows_fft_south_v-nrows_group2_v*(n_factor+1) !    or one additional row to take care of remainders.
!
!---------------------------
!*** Tasks in group 1 for H
!---------------------------
!
        kount_pes=0
        kount_layers=1
        nrow_x=jh_start_fft_south
        n_extra=n_remainder_h_group1
!
        do npe=0,ipe_end_south
          my_jrow_start_h(npe)=nrow_x
          my_jrow_end_h(npe)=min(nrow_x+nrows_group1_h-1,jh_end_fft_south)
!
          if(n_extra>0)then                                               !<-- Use up remainder H lat rows.
            my_jrow_end_h(npe)=my_jrow_end_h(npe)+1
            n_extra=n_extra-1
          endif
!
          k1_fft(npe)=kount_layers
          k2_fft(npe)=kount_layers
          kount_pes=kount_pes+1
!
          if(kount_pes==n_group1)then 
            npe_next=npe+1
            kount_layers=kount_layers+1
            exit                                                          !<-- Now move on to group 2 tasks
          endif
!
          if(my_jrow_end_h(npe)==jh_end_fft_south)then                    !<-- Ready to move down to next model layer.
            nrow_x=jh_start_fft_south
            kount_layers=kount_layers+1
            n_extra=n_remainder_h_group1
          else                                                            !<-- Still divvying up this model layer.
            nrow_x=my_jrow_end_h(npe)+1
          endif
!
        enddo
!
!---------------------------
!*** Tasks in group 2 for H
!---------------------------
!
        if(npe_next<=ipe_end_south)then
          kount_pes=0
          nrow_x=jh_start_fft_south
          n_extra=n_remainder_h_group2
!
          do npe=npe_next,ipe_end_south
            my_jrow_start_h(npe)=nrow_x
            my_jrow_end_h(npe)=min(nrow_x+nrows_group2_h-1,jh_end_fft_south)
!
            if(n_extra>0)then                                             !<-- Use up remainder H lat rows.
              my_jrow_end_h(npe)=my_jrow_end_h(npe)+1
              n_extra=n_extra-1
            endif
!
            k1_fft(npe)=kount_layers
            k2_fft(npe)=kount_layers
            kount_pes=kount_pes+1
            if(kount_pes==n_group2)then                                   !<-- All Southern Hemisphere tasks are assigned to FFT's
              exit                                                        !    for H points.
            endif
            if(my_jrow_end_h(npe)==jh_end_fft_south)then                  !<-- Ready to move down to next model layer.
              nrow_x=jh_start_fft_south
              kount_layers=kount_layers+1
              n_extra=n_remainder_h_group2
            else                                                          !<-- Still divvying up this model layer.
              nrow_x=my_jrow_end_h(npe)+1
            endif
          enddo
        endif
!
!---------------------------
!*** Tasks in group 1 for V
!---------------------------
!
        kount_pes=0
        kount_layers=1
        nrow_x=jv_start_fft_south
        n_extra=n_remainder_v_group1
!
        do npe=0,ipe_end_south
          my_jrow_start_v(npe)=nrow_x
          my_jrow_end_v(npe)=min(nrow_x+nrows_group1_v-1,jv_end_fft_south)
!
          if(n_extra>0)then                                               !<-- Use up remainder V lat rows.
            my_jrow_end_v(npe)=my_jrow_end_v(npe)+1
            n_extra=n_extra-1
          endif
!
          k1_fft(npe)=kount_layers
          k2_fft(npe)=kount_layers
          kount_pes=kount_pes+1
!
          if(kount_pes==n_group1)then
            npe_next=npe+1
            kount_layers=kount_layers+1
            exit                                                          !<-- Move on to group 2 tasks.
          endif
!
          if(my_jrow_end_v(npe)==jv_end_fft_south)then                    !<-- Ready to move down to next model layer.
            nrow_x=jv_start_fft_south
            kount_layers=kount_layers+1
            n_extra=n_remainder_v_group1
          else                                                            !<-- Still divvying up this model layer.
            nrow_x=my_jrow_end_v(npe)+1
          endif
!
        enddo
!
!---------------------------
!*** Tasks in group 2 for V
!---------------------------
!
        if(npe_next<=ipe_end_south)then
          kount_pes=0
          nrow_x=jv_start_fft_south
          n_extra=n_remainder_v_group2
!
          do npe=npe_next,ipe_end_south
            my_jrow_start_v(npe)=nrow_x
            my_jrow_end_v(npe)=min(nrow_x+nrows_group2_v-1,jv_end_fft_south)
!
            if(n_extra>0)then                                             !<-- Use up remainder V lat rows.
              my_jrow_end_v(npe)=my_jrow_end_v(npe)+1
              n_extra=n_extra-1
            endif
!
            k1_fft(npe)=kount_layers
            k2_fft(npe)=kount_layers
            kount_pes=kount_pes+1
!
            if(kount_pes==n_group2)then                                   !<-- All Southern Hemisphere tasks are assigned to FFTs
              exit                                                        !    for V points.
            endif
!
            if(my_jrow_end_v(npe)==jv_end_fft_south)then                  !<-- Ready to move down to next model layer.
              nrow_x=jv_start_fft_south
              kount_layers=kount_layers+1
              n_extra=n_remainder_v_group2
            else                                                          !<-- Still divvying up this model layer.
              nrow_x=my_jrow_end_v(npe)+1
            endif
!
          enddo
        endif
!
      endif limits_south
!
!-------------------------
!***  Northern Hemisphere
!-------------------------
!
      npes_north=ipe_end_north-ipe_start_north+1
!
!----------------------------------------------------
!***  The number of tasks in the Northern Hemisphere
!***  does not exceed the number of model layers.
!----------------------------------------------------
!
      limits_north: if(npes_north<=lm)then
!
        lyr_frac_north=lm/npes_north
        l_remain=lm-npes_north*lyr_frac_north
!
        k2=0
        do npe=ipe_start_north,ipe_end_north
          k1_fft(npe)=k2+1
          k2=k1_fft(npe)+lyr_frac_north-1
          if(l_remain>0)then
            k2=k2+1
            l_remain=l_remain-1
          endif
          k2_fft(npe)=k2
!
          my_jrow_start_h(npe)=jh_start_fft_north
          my_jrow_end_h(npe)=jh_end_fft_north
          my_jrow_start_v(npe)=jv_start_fft_north
          my_jrow_end_v(npe)=jv_end_fft_north
        enddo
!
!----------------------------------------------------
!***  If there are more tasks in the hemisphere
!***  than there are model layers then divide the
!***  layers into n_factor pieces for tasks in
!***  n_group1 and divide the remaining layers into
!***  n_factor+1 pieces for tasks in n_group2.
!----------------------------------------------------
!
      else
        lyr_frac_north=0
        n_factor=npes_north/lm
        n_remain=npes_north-n_factor*lm
        n_group1=n_factor*(lm-n_remain)                                   !<-- This many tasks get layers divided into n_factor pieces
        n_group2=npes_north-n_group1                                      !<-- This many tasks get layers divided into n_factor+1 pieces
!
!----------------------------------------------------
!***  Divide layers of FFTs into n_factor pieces
!***  for tasks in n_group1.
!***  Divide remaining layers into n_factor+1 pieces
!***  for tasks in n_group2.
!----------------------------------------------------
!
        nrows_fft_north_h=jh_end_fft_north-jh_start_fft_north+1
!
        nrows_group1_h=nrows_fft_north_h/n_factor                          !<-- Each task in group 1 handles this many H lat rows
        n_remainder_h_group1=nrows_fft_north_h-nrows_group1_h*n_factor     !    or one additional row to take care of remainders.
!
        nrows_group2_h=nrows_fft_north_h/(n_factor+1)                      !<-- Each task in group 2 handles this many H lat rows
        n_remainder_h_group2=nrows_fft_north_h-nrows_group2_h*(n_factor+1) !    or one additional row to take care of remainders.
!
        nrows_fft_north_v=jv_end_fft_north-jv_start_fft_north+1
!
        nrows_group1_v=nrows_fft_north_v/n_factor                          !<-- Each task in group 1 handles this many V lat rows
        n_remainder_v_group1=nrows_fft_north_v-nrows_group1_v*n_factor     !    or one additional row to take care of remainders.
!
        nrows_group2_v=nrows_fft_north_v/(n_factor+1)                      !<-- Each task in group 2 handles this many V lat rows
        n_remainder_v_group2=nrows_fft_north_v-nrows_group2_v*(n_factor+1) !    or one additional row to take care of remainders
!
!---------------------------
!*** Tasks in group 1 for H
!---------------------------
!
        kount_pes=0
        kount_layers=1
        nrow_x=jh_start_fft_north
        n_extra=n_remainder_h_group1
!
        do npe=ipe_start_north,ipe_end_north
          my_jrow_start_h(npe)=nrow_x
          my_jrow_end_h(npe)=min(nrow_x+nrows_group1_h-1,jh_end_fft_north)
!
          if(n_extra>0)then                                               !<-- Use up remainder H lat rows.
            my_jrow_end_h(npe)=my_jrow_end_h(npe)+1
            n_extra=n_extra-1
          endif
!
          k1_fft(npe)=kount_layers
          k2_fft(npe)=kount_layers
          kount_pes=kount_pes+1
!
          if(kount_pes==n_group1)then 
            npe_next=npe+1
            kount_layers=kount_layers+1
            exit                                                          !<-- Now move on to group 2 tasks.
          endif
!
          if(my_jrow_end_h(npe)==jh_end_fft_north)then                    !<-- Ready to move down to next model layer.
            nrow_x=jh_start_fft_north
            kount_layers=kount_layers+1
            n_extra=n_remainder_h_group1
          else
            nrow_x=my_jrow_end_h(npe)+1                                   !<-- Still divvying up this model layer.
          endif
!
        enddo
!
!---------------------------
!*** Tasks in group 2 for H
!---------------------------
!
        if(npe_next<=ipe_end_north)then
          kount_pes=0
          nrow_x=jh_start_fft_north
          n_extra=n_remainder_h_group2
!
          do npe=npe_next,ipe_end_north
            my_jrow_start_h(npe)=nrow_x
            my_jrow_end_h(npe)=min(nrow_x+nrows_group2_h-1,jh_end_fft_north)
!
            if(n_extra>0)then                                             !<-- Use up remainder H lat rows.
              my_jrow_end_h(npe)=my_jrow_end_h(npe)+1
              n_extra=n_extra-1
            endif
!
            k1_fft(npe)=kount_layers
            k2_fft(npe)=kount_layers
            kount_pes=kount_pes+1
!
            if(kount_pes==n_group2)then                                   !<-- All Northern Hemisphere tasks are assigned to FFT's
              exit                                                        !    for H points.
            endif
!
            if(my_jrow_end_h(npe)==jh_end_fft_north)then
              nrow_x=jh_start_fft_north
              kount_layers=kount_layers+1                                 !<-- Ready to move down to next model layer.
              n_extra=n_remainder_h_group2
            else
              nrow_x=my_jrow_end_h(npe)+1                                 !<-- Still divvying up rows in this model layer.                
            endif
!
          enddo
        endif
!
!---------------------------
!*** Tasks in group 1 for V
!---------------------------
!
        kount_pes=0
        kount_layers=1
        nrow_x=jv_start_fft_north
        n_extra=n_remainder_v_group1
!
        do npe=ipe_start_north,ipe_end_north
          my_jrow_start_v(npe)=nrow_x
          my_jrow_end_v(npe)=min(nrow_x+nrows_group1_v-1,jv_end_fft_north)
!
          if(n_extra>0)then                                               !<-- Use up remainder H lat rows.
            my_jrow_end_v(npe)=my_jrow_end_v(npe)+1                       
            n_extra=n_extra-1
          endif
!
          k1_fft(npe)=kount_layers
          k2_fft(npe)=kount_layers
          kount_pes=kount_pes+1
!
          if(kount_pes==n_group1)then
            npe_next=npe+1
            kount_layers=kount_layers+1
            exit                                                          !<-- Now move on to group 2 tasks.
          endif
!
          if(my_jrow_end_v(npe)==jv_end_fft_north)then
            nrow_x=jv_start_fft_north
            kount_layers=kount_layers+1                                   !<-- Ready to move down to next model layer.
            n_extra=n_remainder_v_group1
          else
            nrow_x=my_jrow_end_v(npe)+1                                   !<-- Still divvying up rows in this model layer.
          endif
!
        enddo
!
!---------------------------
!*** Tasks in group 2 for V
!---------------------------
!
        if(npe_next<=ipe_end_north)then
          kount_pes=0
          nrow_x=jv_start_fft_north
          n_extra=n_remainder_v_group2
!
          do npe=npe_next,ipe_end_north
            my_jrow_start_v(npe)=nrow_x
            my_jrow_end_v(npe)=min(nrow_x+nrows_group2_v-1,jv_end_fft_north)
!
            if(n_extra>0)then                                             !<-- Use up remainder H lat rows.
              my_jrow_end_v(npe)=my_jrow_end_v(npe)+1
              n_extra=n_extra-1
            endif
!
            k1_fft(npe)=kount_layers
            k2_fft(npe)=kount_layers
            kount_pes=kount_pes+1
!
            if(kount_pes==n_group2)then                                   !<-- All Northern Hemisphere tasks are assigned to FFTs
              exit                                                        !    for V points.
            endif
!
            if(my_jrow_end_v(npe)==jv_end_fft_north)then
              nrow_x=jv_start_fft_north
              kount_layers=kount_layers+1                                 !<-- Ready to move down to next model layer.
              n_extra=n_remainder_v_group2
            else
              nrow_x=my_jrow_end_v(npe)+1                                 !<-- Still divvying up rows in this model layer.
            endif
!
          enddo
        endif
!
      endif limits_north
!
!-----------------------------------------------------------------------
!***  Allocate the working arrays for the FFTs depending on which
!***  hemisphere this task is in.
!-----------------------------------------------------------------------
!
      lm_fft=max(lyr_frac_north,lyr_frac_south)+1
!
      allocate(hn(ids:ide,my_jrow_start_h(mype):my_jrow_end_h(mype),1:lm_fft) &
              ,stat=istat)
      allocate(wn(ids:ide,my_jrow_start_v(mype):my_jrow_end_v(mype),1:lm_fft) &
              ,stat=istat)
      allocate(un(ids:ide,my_jrow_start_v(mype):my_jrow_end_v(mype),1:lm_fft) &
              ,stat=istat)
      allocate(vn(ids:ide,my_jrow_start_v(mype):my_jrow_end_v(mype),1:lm_fft) &
              ,stat=istat)
!
!-----------------------------------------------------------------------
!***  Maximum size of dummy space for FFTs
!-----------------------------------------------------------------------
!
      imax=(ide-ids+1)/inpes+1
      jmax_south=max(jh_end_fft_south-jh_start_fft_south &
                    ,jv_end_fft_south-jv_start_fft_south)+1
      jmax_north=max(jh_end_fft_north-jh_start_fft_north &
                    ,jv_end_fft_north-jv_start_fft_north)+1
      jmax=max(jmax_south,jmax_north)
      msize_dummy_fft=imax*jmax*lm_fft
!
!-----------------------------------------------------------------------
!
                        endsubroutine prefft
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
                        subroutine fftfhn &
(km &
,khfilt &
,hfilt &
,field_h &
,wfftrh,nfftrh &
,npes,mype,mpi_comm_comp)

!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------

integer(kind=kint),intent(in) :: &
 km &
,mpi_comm_comp & ! This domain's fcst task intracommunicator
,mype &          ! Rank of this task in the fcst intracommunicator
,npes            ! Number of compute tasks

integer(kind=kint),dimension(1:15),intent(in):: &
 nfftrh

integer(kind=kint),dimension(jds:jde),intent(in):: &
 khfilt

real(kind=kfpt),dimension(ids:ide,jds:jde),intent(in):: &
 hfilt

real(kind=kfpt),dimension(ims:ime,jms:jme,1:km),intent(inout):: &
 field_h

real(kind=kfpt),dimension(1:2*(ide-3)),intent(in) :: &
 wfftrh  

!-----------------------------------------------------------------------
!***  Local Variables
!-----------------------------------------------------------------------

integer(kind=kint) :: &
 i &
,icycle &
,ipe_end &
,ipe_start &
,j &
,jend &
,jstart &
,jh_end_fft &
,jh_start_fft &
,k1 &
,k2 &
,l &
,n &
,nend

real(kind=kfpt) :: &
 an &
,as &
,rcycle

!real(kind=kfpt),dimension(ids:ide,jts:jte,1:lm_fft) :: &
!real(kind=kfpt),dimension(ids:ide,jh_start_fft:jh_end_fft,1:lm_fft) :: &
! hn

real(kind=kfpt),dimension(1:ide-3):: &
 buff

integer :: ierr,ixx,jxx,kxx
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
      icycle=ide-3
      rcycle=1./icycle
!-----------------------------------------------------------------------
!***  First take care of the special case in which only a single
!***  MPI task has been designated for the entire domain.
!-----------------------------------------------------------------------
!
      if(npes==1)then
!
!-----------------------------------------------------------------------
!
        as=0.
        an=0.
!
        do l=1,km
!
          do i=ids+1,ide-2
            as=field_h(i,jds+1,l)+as
            an=field_h(i,jde-1,l)+an
          enddo
!
          as=as*rcycle
          an=an*rcycle
!
          do i=ids,ide
            field_h(i,jds+1,l)=as
            field_h(i,jde-1,l)=an
          enddo
!
          do j=jds+2,jde-2
!
            if(khfilt(j)<=icycle) then
!
              do i=ids+1,ide-2
                buff(i-1)=field_h(i,j,l)
              enddo
!
              call rfftf(icycle,buff,wfftrh,nfftrh)
!
              do i=1,khfilt(j)-1
                buff(i)=buff(i)*hfilt(i,j)
              enddo
!
              do i=khfilt(j),icycle
                buff(i)=0.
              enddo
!
              call rfftb(icycle,buff,wfftrh,nfftrh)
!
              do i=ids+1,ide-2
                field_h(i,j,l)=buff(i-1)
              enddo
!
              field_h(ide-1,j,l)=buff(1)
!
            endif
!
          enddo
!
        enddo
!
        return
!
!-----------------------------------------------------------------------
!
      endif
!
!-----------------------------------------------------------------------
!***  The remainder of this subroutine is relevant when more than 
!***  a single task has been specified.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  k1 and k2 are starting/ending vertical indices of model layers
!***  that this task will use in applying FFTs.
!-----------------------------------------------------------------------
!
      k1=k1_fft(mype)
      k2=k2_fft(mype)
!
!-----------------------------------------------------------------------
!***  Select hemisphere-dependent variables.
!-----------------------------------------------------------------------
!
      if(fft_south)then
        jh_start_fft=jh_start_fft_south
        jh_end_fft=jh_end_fft_south
        ipe_start=ipe_start_south
        ipe_end=ipe_end_south
      elseif(fft_north)then
        jh_start_fft=jh_start_fft_north
        jh_end_fft=jh_end_fft_north
        ipe_start=ipe_start_north
        ipe_end=ipe_end_north
      endif
!
!-----------------------------------------------------------------------
!***  Gather the model layer data from full latitude circles
!***  onto appropriate tasks for the FFTs.
!-----------------------------------------------------------------------
!
      call gather_layers(field_h,km,npes,msize_dummy_fft &
                        ,lm_fft,k1_fft,k2_fft &
                        ,local_istart,local_iend &
                        ,local_jstart,local_jend &
                        ,my_jrow_start_h(mype),my_jrow_end_h(mype) &
                        ,my_jrow_start_h,my_jrow_end_h &
                        ,ipe_start,ipe_end &
                        ,my_domain_has_fft_lats_h &
                        ,mype,mpi_comm_comp &
                        ,its,ite,ims,ime,ids,ide,jts,jte,jms,jme &
                        ,hn)
!
!-----------------------------------------------------------------------
!
      nend=k2-k1+1
      kloop1: do n=1,nend
!
!-----------------------------------------------------------------------
        if(fft_south)then
!-----------------------------------------------------------------------
          as=0.
!
          if(lbound(hn,2)==jds+1)then
            do i=ids+1,ide-2
              as=hn(i,jds+1,n)+as
            enddo
!
            as=as*rcycle
!
            do i=ids,ide
              hn(i,jds+1,n)=as
            enddo
          endif
!
!-----------------------------------------------------------------------
        elseif(fft_north)then
!-----------------------------------------------------------------------
          an=0.
!
          if(ubound(hn,2)==jde-1)then
            do i=ids+1,ide-2
              an=hn(i,jde-1,n)+an
            enddo
!
            an=an*rcycle
!
            do i=ids,ide
              hn(i,jde-1,n)=an
            enddo
          endif
!
!-----------------------------------------------------------------------
        endif   
!-----------------------------------------------------------------------
!
      enddo kloop1
!
!-----------------------------------------------------------------------
!***  jstart and jend are the starting/ending rows on which this task
!***  will apply FFTs
!-----------------------------------------------------------------------
!
      jstart=max(my_jrow_start_h(mype),jds+2)
      jend=min(my_jrow_end_h(mype),jde-2)
!
!-----------------------------------------------------------------------
!
      kloop2: do n=1,nend
!
!-----------------------------------------------------------------------
!
        do j=jstart,jend
!
!-----------------------------------------------------------------------
          if(khfilt(j)<=icycle) then
!-----------------------------------------------------------------------
            do i=ids+1,ide-2
              buff(i-1)=hn(i,j,n)
            enddo
!
            call rfftf(icycle,buff,wfftrh,nfftrh)
!
            do i=1,khfilt(j)-1
              buff(i)=buff(i)*hfilt(i,j)
            enddo
!
            do i=khfilt(j),icycle
              buff(i)=0.
            enddo
!
            call rfftb(icycle,buff,wfftrh,nfftrh)
!
            do i=ids+1,ide-2
              hn(i,j,n)=buff(i-1)
            enddo
!
            hn(ide-1,j,n)=buff(1)
!-----------------------------------------------------------------------
          endif
!-----------------------------------------------------------------------
        enddo
!-----------------------------------------------------------------------
!
      enddo kloop2
!
!-----------------------------------------------------------------------
!***  Now scatter the model layer data from full latitude circles
!***  back to the appropriate tasks.
!-----------------------------------------------------------------------
!
      call scatter_layers(hn,km,npes,msize_dummy_fft &
                         ,lm_fft,k1_fft,k2_fft &
                         ,local_istart,local_iend &
                         ,local_jstart,local_jend &
                         ,my_jrow_start_h(mype),my_jrow_end_h(mype) &
                         ,my_jrow_start_h,my_jrow_end_h &
                         ,ipe_start,ipe_end &
                         ,my_domain_has_fft_lats_h &
                         ,mype,mpi_comm_comp &
                         ,its,ite,ims,ime,ids,ide,jts,jte,jms,jme &
                         ,field_h)
!
!-----------------------------------------------------------------------
!
                        endsubroutine fftfhn
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
                        subroutine fftfwn &
(km &
,kvfilt &
,vfilt &
,field_w &
,wfftrw,nfftrw &
,npes)
!
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------

integer(kind=kint),intent(in) :: &
 km &
,npes      ! Number of compute tasks

integer(kind=kint),dimension(1:15),intent(in):: &
 nfftrw

integer(kind=kint),dimension(jds:jde),intent(in):: &
 kvfilt

real(kind=kfpt),dimension(ids:ide,jds:jde),intent(in):: &
 vfilt

real(kind=kfpt),dimension(ims:ime,jms:jme,1:km),intent(inout):: &
 field_w

!-----------------------------------------------------------------------
!***  Local Variables
!-----------------------------------------------------------------------

integer(kind=kint) :: &
 i &
,ipe_start &
,ipe_end &
,icycle &
,iloc_mype &
,j &
,jend &
,jstart &
,jv_end_fft &
,jv_start_fft &
,k1 &
,k2 &
,l &
,n &
,nend

real(kind=kfpt) :: &
 rcycle

!real(kind=kfpt),dimension(ids:ide,jts:jte,1:lm_fft) :: &
!real(kind=kfpt),dimension(ids:ide,jv_start_fft:jv_end_fft,1:lm_fft) :: &
! wn

real(kind=kfpt),dimension(1:ide-3):: &
 buff

real(kind=kfpt),dimension(1:2*(ide-3)):: &
 wfftrw

!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
      icycle=ide-3
      rcycle=1./icycle
!-----------------------------------------------------------------------
!***  k1 and k2 are starting/ending vertical indices of model layers
!***  that this task will use in applying FFTs.
!-----------------------------------------------------------------------
!
      k1=k1_fft(mype)
      k2=k2_fft(mype)
!
!-----------------------------------------------------------------------
!***  Select hemisphere-dependent variables.
!-----------------------------------------------------------------------
!
      if(fft_south)then
        jv_start_fft=jv_start_fft_south
        jv_end_fft=jv_end_fft_south
        ipe_start=ipe_start_south
        ipe_end=ipe_end_south
      elseif(fft_north)then
        jv_start_fft=jv_start_fft_north
        jv_end_fft=jv_end_fft_north
        ipe_start=ipe_start_north
        ipe_end=ipe_end_north
      endif
!
!-----------------------------------------------------------------------
!***  Gather the model layer data from full latitude circles
!***  onto appropriate tasks for the FFTs.
!-----------------------------------------------------------------------
!
      call gather_layers(field_w,km,npes,msize_dummy_fft &
                        ,lm_fft,k1_fft,k2_fft &
                        ,local_istart,local_iend &
                        ,local_jstart,local_jend &
                        ,my_jrow_start_v(mype),my_jrow_end_v(mype) &
                        ,my_jrow_start_v,my_jrow_end_v &
                        ,ipe_start,ipe_end &
                        ,my_domain_has_fft_lats_v &
                        ,mype,mpi_comm_comp &
                        ,its,ite,ims,ime,ids,ide,jts,jte,jms,jme &
                        ,wn)
!
!-----------------------------------------------------------------------
!***  jstart and jend are the starting/ending rows on which this task
!***  will apply FFTs
!-----------------------------------------------------------------------
!
      jstart=max(my_jrow_start_v(mype),jds+1)
      jend=min(my_jrow_end_v(mype),jde-2)
!
!-----------------------------------------------------------------------
!
      nend=k2-k1+1
      kloop1: do n=1,nend
!
        do j=jstart,jend
!
!-----------------------------------------------------------------------
          if(kvfilt(j)<=icycle) then
!-----------------------------------------------------------------------
            do i=ids+1,ide-2
              buff(i-1)=wn(i,j,n)
            enddo
!
            call rfftf(icycle,buff,wfftrw,nfftrw)
!
            do i=1,kvfilt(j)-1
              buff(i)=buff(i)*vfilt(i,j)
            enddo
            do i=kvfilt(j),icycle
              buff(i)=0.
            enddo
!
            call rfftb(icycle,buff,wfftrw,nfftrw)
!
            do i=ids+1,ide-2
              wn(i,j,n)=buff(i-1)
            enddo
!
            wn(ide-1,j,n)=buff(1)
!-----------------------------------------------------------------------
          endif
!-----------------------------------------------------------------------
        enddo
!-----------------------------------------------------------------------
      enddo kloop1
!-----------------------------------------------------------------------
!***  Now scatter the model layer data from full latitude circles
!***  back to the appropriate tasks.
!-----------------------------------------------------------------------
!
      call scatter_layers(wn,km,npes,msize_dummy_fft &
                         ,lm_fft,k1_fft,k2_fft &
                         ,local_istart,local_iend &
                         ,local_jstart,local_jend &
                         ,my_jrow_start_v(mype),my_jrow_end_v(mype) &
                         ,my_jrow_start_v,my_jrow_end_v &
                         ,ipe_start,ipe_end &
                         ,my_domain_has_fft_lats_v &
                         ,mype,mpi_comm_comp &
                         ,its,ite,ims,ime,ids,ide,jts,jte,jms,jme &
                         ,field_w)
!
!-----------------------------------------------------------------------
!
                        endsubroutine fftfwn
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
                        subroutine fftfuvn &
(km &
,kvfilt &
,vfilt &
,u,v &
,wfftrw,nfftrw &
,npes,mype,mpi_comm_comp)
!
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------

integer(kind=kint),intent(in) :: &
 km &
,mpi_comm_comp & ! This domain's fcst task intracommunicator
,mype &          ! Rank of this task in the fcst intracommunicator
,npes            ! Number of compute tasks
                                                                                                                                              
integer(kind=kint),dimension(1:15),intent(in):: &
 nfftrw
                                                                                                                                              
integer(kind=kint),dimension(jds:jde),intent(in):: &
 kvfilt
                                    
real(kind=kfpt),dimension(ids:ide,jds:jde),intent(in):: &
 vfilt
                                    
real(kind=kfpt),dimension(ims:ime,jms:jme,1:km),intent(inout):: &
 u &
,v
                                    
!-----------------------------------------------------------------------
!***  Local Variables
!-----------------------------------------------------------------------
                                    
integer(kind=kint) :: &
 i &
,ipe_start &
,ipe_end &
,icycle &
,iloc_mype &
,j &
,jend &
,jstart &
,jv_end_fft &
,jv_start_fft &
,k1 &
,k2 &
,l &
,n
                                    
real(kind=kfpt) :: &
 anu &
,anv &
,asu &
,asv &
,rcycle
     
real(kind=kfpt),dimension(1:ide-3):: &
 buffu &
,buffv
       
real(kind=kfpt),dimension(1:2*(ide-3)):: &
 wfftrw
       
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
      icycle=ide-3
      rcycle=1./icycle
!-----------------------------------------------------------------------
!***  First take care of the special case in which only a single
!***  MPI task has been designated for the entire domain.
!-----------------------------------------------------------------------
!
      if(npes==1)then
!
!-----------------------------------------------------------------------
!
        do l=1,km
!
          do j=jds+1,jde-2
!
            if(kvfilt(j)<=icycle) then
!
              do i=ids+1,ide-2
                buffu(i-1)=u(i,j,l)
                buffv(i-1)=v(i,j,l)
              enddo
!
              call rfftf(icycle,buffu,wfftrw,nfftrw)
              call rfftf(icycle,buffv,wfftrw,nfftrw)
!
              do i=1,kvfilt(j)-1
                buffu(i)=buffu(i)*vfilt(i,j)
                buffv(i)=buffv(i)*vfilt(i,j)
              enddo
              do i=kvfilt(j),icycle
                buffu(i)=0.
                buffv(i)=0.
              enddo
!
              call rfftb(icycle,buffu,wfftrw,nfftrw)
              call rfftb(icycle,buffv,wfftrw,nfftrw)
!
              do i=ids+1,ide-2
                u(i,j,l)=buffu(i-1)
                v(i,j,l)=buffv(i-1)
              enddo
!
              u(ide-1,j,l)=buffu(1)
              v(ide-1,j,l)=buffv(1)
!
            endif
!
          enddo
!
        enddo
!
        return
!
!-----------------------------------------------------------------------
!
      endif
!
!-----------------------------------------------------------------------
!***  The remainder of this routine is relevant for the use of 
!***  multiple MPI tasks on the domain.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  k1 and k2 are starting/ending vertical indices of model layers
!***  that this task will use in applying FFTs.
!-----------------------------------------------------------------------
!
      k1=k1_fft(mype)
      k2=k2_fft(mype)
!
!-----------------------------------------------------------------------
!***  Select hemisphere-dependent variables.
!-----------------------------------------------------------------------
!
      if(fft_south)then
        jv_start_fft=jv_start_fft_south
        jv_end_fft=jv_end_fft_south
        ipe_start=ipe_start_south
        ipe_end=ipe_end_south
      elseif(fft_north)then
        jv_start_fft=jv_start_fft_north
        jv_end_fft=jv_end_fft_north
        ipe_start=ipe_start_north
        ipe_end=ipe_end_north
      endif
!
!-----------------------------------------------------------------------
!***  Gather the model layer data from full latitude circles
!***  onto appropriate tasks for the FFTs.
!-----------------------------------------------------------------------
!
      call gather_layers(u,km,npes,msize_dummy_fft &
                        ,lm_fft,k1_fft,k2_fft &
                        ,local_istart,local_iend &
                        ,local_jstart,local_jend &
                        ,my_jrow_start_v(mype),my_jrow_end_v(mype) &
                        ,my_jrow_start_v,my_jrow_end_v &
                        ,ipe_start,ipe_end &
                        ,my_domain_has_fft_lats_v &
                        ,mype,mpi_comm_comp &
                        ,its,ite,ims,ime,ids,ide,jts,jte,jms,jme &
                        ,un)
      call gather_layers(v,km,npes,msize_dummy_fft &
                        ,lm_fft,k1_fft,k2_fft &
                        ,local_istart,local_iend &
                        ,local_jstart,local_jend &
                        ,my_jrow_start_v(mype),my_jrow_end_v(mype) &
                        ,my_jrow_start_v,my_jrow_end_v &
                        ,ipe_start,ipe_end &
                        ,my_domain_has_fft_lats_v &
                        ,mype,mpi_comm_comp &
                        ,its,ite,ims,ime,ids,ide,jts,jte,jms,jme &
                        ,vn)
!
!-----------------------------------------------------------------------
!
!      n=0
!      kloop1: do l=k1,k2
!
!-----------------------------------------------------------------------
!        asu=0.
!        asv=0.
!        anu=0.
!        anv=0.
!        n=n+1
!-----------------------------------------------------------------------
!        if(fft_south)then
!-----------------------------------------------------------------------
!          do i=ids+1,ide-2
!            asu=un(i,jds+1,n)+asu
!            asv=vn(i,jds+1,n)+asv
!          enddo
!
!          asu=asu*rcycle
!          asv=asv*rcycle
!
!          do i=ids,ide
!            un(i,jds+1,n)=un(i,jds+1,n)-asu
!            vn(i,jds+1,n)=vn(i,jds+1,n)-asv
!          enddo
!-----------------------------------------------------------------------
!        elseif(fft_north)then
!-----------------------------------------------------------------------
!          do i=ids+1,ide-2
!            anu=un(i,jde-2,n)+anu
!            anv=vn(i,jde-2,n)+anv
!          enddo
!
!          anu=anu*rcycle
!          anv=anv*rcycle
!
!          do i=ids,ide
!            un(i,jde-2,n)=un(i,jde-2,n)-anu
!            vn(i,jde-2,n)=vn(i,jde-2,n)-anv
!          enddo
!-----------------------------------------------------------------------
!        endif
!-----------------------------------------------------------------------
!
!      enddo kloop1
!
!-----------------------------------------------------------------------
!***  jstart and jend are the starting/ending rows on which this task
!***  will apply FFTs
!-----------------------------------------------------------------------
!
      jstart=max(my_jrow_start_v(mype),jds+1)
      jend=min(my_jrow_end_v(mype),jde-2)
!
!-----------------------------------------------------------------------
!
      n=0
      kloop: do l=k1,k2
!
        n=n+1
!
        do j=jstart,jend
!-----------------------------------------------------------------------
          if(kvfilt(j).le.icycle) then
!-----------------------------------------------------------------------
            do i=ids+1,ide-2
              buffu(i-1)=un(i,j,n)
              buffv(i-1)=vn(i,j,n)
            enddo
!
            call rfftf(icycle,buffu,wfftrw,nfftrw)
            call rfftf(icycle,buffv,wfftrw,nfftrw)
!
            do i=1,kvfilt(j)-1
              buffu(i)=buffu(i)*vfilt(i,j)
              buffv(i)=buffv(i)*vfilt(i,j)
            enddo
            do i=kvfilt(j),icycle
              buffu(i)=0.
              buffv(i)=0.
            enddo
!
            call rfftb(icycle,buffu,wfftrw,nfftrw)
            call rfftb(icycle,buffv,wfftrw,nfftrw)
!
            do i=ids+1,ide-2
              un(i,j,n)=buffu(i-1)
              vn(i,j,n)=buffv(i-1)
            enddo
!
            un(ide-1,j,n)=buffu(1)
            vn(ide-1,j,n)=buffv(1)
!-----------------------------------------------------------------------
          endif
!-----------------------------------------------------------------------
        enddo
!-----------------------------------------------------------------------
      enddo kloop
!-----------------------------------------------------------------------
!***  Now scatter the model layer data from full latitude circles
!***  back to the appropriate tasks.
!-----------------------------------------------------------------------
!
      call scatter_layers(un,km,npes,msize_dummy_fft &
                         ,lm_fft,k1_fft,k2_fft &
                         ,local_istart,local_iend &
                         ,local_jstart,local_jend &
                         ,my_jrow_start_v(mype),my_jrow_end_v(mype) &
                         ,my_jrow_start_v,my_jrow_end_v &
                         ,ipe_start,ipe_end &
                         ,my_domain_has_fft_lats_v &
                         ,mype,mpi_comm_comp &
                         ,its,ite,ims,ime,ids,ide,jts,jte,jms,jme &
                         ,u)
      call scatter_layers(vn,km,npes,msize_dummy_fft &
                         ,lm_fft,k1_fft,k2_fft &
                         ,local_istart,local_iend &
                         ,local_jstart,local_jend &
                         ,my_jrow_start_v(mype),my_jrow_end_v(mype) &
                         ,my_jrow_start_v,my_jrow_end_v &
                         ,ipe_start,ipe_end &
                         ,my_domain_has_fft_lats_v &
                         ,mype,mpi_comm_comp &
                         ,its,ite,ims,ime,ids,ide,jts,jte,jms,jme &
                         ,v)
!
!-----------------------------------------------------------------------
!
                        endsubroutine fftfuvn
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!    subroutine rffti(n,wsave,nsave)                                                 
!                                                                              
!     SUBROUTINE RFFTI INITIALIZES THE ARRAY WSAVE WHICH IS USED IN             
!     BOTH RFFTF AND RFFTB. THE PRIME FACTORIZATION OF N TOGETHER WITH          
!     A TABULATION OF THE TRIGONOMETRIC FUNCTIONS ARE COMPUTED AND              
!     STORED IN WSAVE.                                                          
!                                                                               
!     INPUT PARAMETER                                                           
!                                                                               
!     N       THE LENGTH OF THE SEQUENCE TO BE TRANSFORMED.                     
!                                                                               
!     OUTPUT PARAMETER                                                          
!                                                                               
!     WSAVE   A WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 2*N+15.           
!             THE SAME WORK ARRAY CAN BE USED FOR BOTH RFFTF AND RFFTB          
!             AS LONG AS N REMAINS UNCHANGED. DIFFERENT WSAVE ARRAYS            
!             ARE REQUIRED FOR DIFFERENT VALUES OF N. THE CONTENTS OF           
!             WSAVE MUST NOT BE CHANGED BETWEEN CALLS OF RFFTF OR RFFTB.        
!                                                                               
      subroutine rffti (n,wsave,nsave)                                                
      dimension wsave(*),nsave(*)                                                  
!                                                                               
      if (n .eq. 1) return                                                      
      call rffti1 (n,wsave(n+1),nsave)                                   
      return                                                                    
                        endsubroutine rffti                                                                       
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE RFFTI1 (N,WA,IFAC)                                             
      DIMENSION       WA(*)      ,IFAC(*)    ,NTRYH(4)                          
      DATA NTRYH(1),NTRYH(2),NTRYH(3),NTRYH(4)/4,2,3,5/                         
      NL = N                                                                    
      NF = 0                                                                    
      J = 0                                                                     
  101 J = J+1                                                                   
      if(j.le.4) then
        ntry=ntryh(j)
      else
        ntry=ntry+2
      endif                                                             
  104 NQ = NL/NTRY                                                              
      NR = NL-NTRY*NQ                                                           
      if(nr.ne.0) go to 101
      nf=nf+1                                                       
      IFAC(NF+2) = NTRY                                                         
      NL = NQ                                                                   
      IF (NTRY .NE. 2) GO TO 107                                                
      IF (NF .EQ. 1) GO TO 107                                                  
      DO 106 I=2,NF                                                             
         IB = NF-I+2                                                            
         IFAC(IB+2) = IFAC(IB+1)                                                
  106 CONTINUE                                                                  
      IFAC(3) = 2                                                               
  107 IF (NL .NE. 1) GO TO 104                                                  
      IFAC(1) = N                                                               
      IFAC(2) = NF                                                              
!      TPI = 2.0*3.141592653589793238462643383279502884197169399375105820                                                     
      ARGH = TPI/FLOAT(N)                                                       
      IS = 0                                                                    
      NFM1 = NF-1                                                               
      L1 = 1                                                                    
      IF (NFM1 .EQ. 0) RETURN                                                   
      DO 110 K1=1,NFM1                                                          
         IP = IFAC(K1+2)                                                        
         LD = 0                                                                 
         L2 = L1*IP                                                             
         IDO = N/L2                                                             
         IPM = IP-1                                                             
         DO 109 J=1,IPM                                                         
            LD = LD+L1                                                          
            I = IS                                                              
            ARGLD = FLOAT(LD)*ARGH                                              
            FI = 0.                                                             
            DO 108 II=3,IDO,2                                                   
               I = I+2                                                          
               FI = FI+1.                                                       
               ARG = FI*ARGLD                                                   
               WA(I-1) = COS(ARG)                                               
               WA(I) = SIN(ARG)                                                 
  108       CONTINUE                                                            
            IS = IS+IDO                                                         
  109    CONTINUE                                                               
         L1 = L2                                                                
  110 CONTINUE                                                                  
      RETURN                                                                    
      endsubroutine                                                                       
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     subroutine rfftf(n,r,wsave,nsave)                                               
!                                                                               
!     SUBROUTINE RFFTF COMPUTES THE FOURIER COEFFICIENTS OF A REAL              
!     PERODI! SEQUENCE (FOURIER ANALYSIS). THE TRANSFORM IS DEFINED             
!     BELOW AT OUTPUT PARAMETER R.                                              
!                                                                               
!     INPUT PARAMETERS                                                          
!                                                                               
!     N       THE LENGTH OF THE ARRAY R TO BE TRANSFORMED.  THE METHOD          
!             IS MOST EFFICIENT WHEN N IS A PRODUCT OF SMALL PRIMES.            
!             N MAY CHANGE SO LONG AS DIFFERENT WORK ARRAYS ARE PROVIDED        
!                                                                               
!     R       A REAL ARRAY OF LENGTH N WHICH CONTAINS THE SEQUENCE              
!             TO BE TRANSFORMED                                                 
!                                                                               
!     WSAVE   A WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 2*N+15.           
!             IN THE PROGRAM THAT CALLS RFFTF. THE WSAVE ARRAY MUST BE          
!             INITIALIZED BY CALLING SUBROUTINE RFFTI(N,WSAVE) AND A            
!             DIFFERENT WSAVE ARRAY MUST BE USED FOR EACH DIFFERENT             
!             VALUE OF N. THIS INITIALIZATION DOES NOT HAVE TO BE               
!             REPEATED SO LONG AS N REMAINS UNCHANGED THUS SUBSEQUENT           
!             TRANSFORMS CAN BE OBTAINED FASTER THAN THE FIRST.                 
!             THE SAME WSAVE ARRAY CAN BE USED BY RFFTF AND RFFTB.              
!                                                                               
!                                                                               
!     OUTPUT PARAMETERS                                                         
!                                                                               
!     R       R(1) = THE SUM FROM I=1 TO I=N OF R(I)                            
!                                                                               
!             IF N IS EVEN SET L =N/2   , IF N IS ODD SET L = (N+1)/2           
!                                                                               
!               THEN FOR K = 2,...,L                                            
!                                                                               
!                  R(2*K-2) = THE SUM FROM I = 1 TO I = N OF                    
!                                                                               
!                       R(I)*COS((K-1)*(I-1)*2*PI/N)                            
!                                                                               
!                  R(2*K-1) = THE SUM FROM I = 1 TO I = N OF                    
!                                                                               
!                      -R(I)*SIN((K-1)*(I-1)*2*PI/N)                            
!                                                                               
!             IF N IS EVEN                                                      
!                                                                               
!                  R(N) = THE SUM FROM I = 1 TO I = N OF                        
!                                                                               
!                       (-1)**(I-1)*R(I)                                        
!                                                                               
!      *****  NOTE                                                              
!                  THIS TRANSFORM IS UNNORMALIZED SINCE A CALL OF RFFTF         
!                  FOLLOWED BY A CALL OF RFFTB WILL MULTIPLY THE INPUT          
!                  SEQUENCE BY N.                                               
!                                                                               
!     WSAVE   CONTAINS RESULTS WHICH MUST NOT BE DESTROYED BETWEEN              
!             CALLS OF RFFTF OR RFFTB.                                          
!                                                                               
      subroutine rfftf (n,r,wsave,nsave)                                              
      dimension r(*),wsave(*),nsave(*)                                      
!                                                                               
      if (n .eq. 1) return                                                      
      call rfftf1 (n,r,wsave,wsave(n+1),nsave)                           
      return                                                                    
                        endsubroutine rfftf                                                                     
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE RFFTF1 (N,C,CH,WA,IFAC)                                        
      DIMENSION       CH(*)      ,C(*)       ,WA(*)      ,IFAC(*)               
      NF = IFAC(2)                                                              
      NA = 1                                                                    
      L2 = N                                                                    
      IW = N                                                                    
      DO 111 K1=1,NF                                                            
         KH = NF-K1                                                             
         IP = IFAC(KH+3)                                                        
         L1 = L2/IP                                                             
         IDO = N/L2                                                             
         IDL1 = IDO*L1                                                          
         IW = IW-(IP-1)*IDO                                                     
         NA = 1-NA                                                              
         IF (IP .NE. 4) GO TO 102                                               
         IX2 = IW+IDO                                                           
         IX3 = IX2+IDO                                                          
         IF (NA .NE. 0) GO TO 101                                               
         CALL RADF4 (IDO,L1,C,CH,WA(IW),WA(IX2),WA(IX3))                        
         GO TO 110                                                              
  101    CALL RADF4 (IDO,L1,CH,C,WA(IW),WA(IX2),WA(IX3))                        
         GO TO 110                                                              
  102    IF (IP .NE. 2) GO TO 104                                               
         IF (NA .NE. 0) GO TO 103                                               
         CALL RADF2 (IDO,L1,C,CH,WA(IW))                                        
         GO TO 110                                                              
  103    CALL RADF2 (IDO,L1,CH,C,WA(IW))                                        
         GO TO 110                                                              
  104    IF (IP .NE. 3) GO TO 106                                               
         IX2 = IW+IDO                                                           
         IF (NA .NE. 0) GO TO 105                                               
         CALL RADF3 (IDO,L1,C,CH,WA(IW),WA(IX2))                                
         GO TO 110                                                              
  105    CALL RADF3 (IDO,L1,CH,C,WA(IW),WA(IX2))                                
         GO TO 110                                                              
  106    IF (IP .NE. 5) GO TO 108                                               
         IX2 = IW+IDO                                                           
         IX3 = IX2+IDO                                                          
         IX4 = IX3+IDO                                                          
         IF (NA .NE. 0) GO TO 107                                               
         CALL RADF5 (IDO,L1,C,CH,WA(IW),WA(IX2),WA(IX3),WA(IX4))                
         GO TO 110                                                              
  107    CALL RADF5 (IDO,L1,CH,C,WA(IW),WA(IX2),WA(IX3),WA(IX4))                
         GO TO 110                                                              
  108    IF (IDO .EQ. 1) NA = 1-NA                                              
         IF (NA .NE. 0) GO TO 109                                               
         CALL RADFG (IDO,IP,L1,IDL1,C,C,C,CH,CH,WA(IW))                         
         NA = 1                                                                 
         GO TO 110                                                              
  109    CALL RADFG (IDO,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW))                        
         NA = 0                                                                 
  110    L2 = L1                                                                
  111 CONTINUE                                                                  
      IF (NA .EQ. 1) RETURN                                                     
      DO 112 I=1,N                                                              
         C(I) = CH(I)                                                           
  112 CONTINUE                                                                  
      RETURN                                                                    
      endsubroutine
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE RADF2 (IDO,L1,CC,CH,WA1)                                       
      DIMENSION       CH(IDO,2,L1)           ,CC(IDO,L1,2)           , &        
                      WA1(*)                                                    
      DO 101 K=1,L1                                                             
         CH(1,1,K) = CC(1,K,1)+CC(1,K,2)                                        
         CH(IDO,2,K) = CC(1,K,1)-CC(1,K,2)                                      
  101 CONTINUE                                                                  
      IF (IDO-2) 107,105,102                                                    
  102 IDP2 = IDO+2                                                              
      DO 104 K=1,L1                                                             
         DO 103 I=3,IDO,2                                                       
            IC = IDP2-I                                                         
            TR2 = WA1(I-2)*CC(I-1,K,2)+WA1(I-1)*CC(I,K,2)                       
            TI2 = WA1(I-2)*CC(I,K,2)-WA1(I-1)*CC(I-1,K,2)                       
            CH(I,1,K) = CC(I,K,1)+TI2                                           
            CH(IC,2,K) = TI2-CC(I,K,1)                                          
            CH(I-1,1,K) = CC(I-1,K,1)+TR2                                       
            CH(IC-1,2,K) = CC(I-1,K,1)-TR2                                      
  103    CONTINUE                                                               
  104 CONTINUE                                                                  
      IF (MOD(IDO,2) .EQ. 1) RETURN                                             
  105 DO 106 K=1,L1                                                             
         CH(1,2,K) = -CC(IDO,K,2)                                               
         CH(IDO,1,K) = CC(IDO,K,1)                                              
  106 CONTINUE                                                                  
  107 RETURN                                                                    
      endsubroutine                                                                      
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE RADF3 (IDO,L1,CC,CH,WA1,WA2)                                   
      DIMENSION       CH(IDO,3,L1)           ,CC(IDO,L1,3)           , &          
                      WA1(*)     ,WA2(*)                                        
      DATA TAUR,TAUI /-.5,.866025403784439/                                     
      DO 101 K=1,L1                                                             
         CR2 = CC(1,K,2)+CC(1,K,3)                                              
         CH(1,1,K) = CC(1,K,1)+CR2                                              
         CH(1,3,K) = TAUI*(CC(1,K,3)-CC(1,K,2))                                 
         CH(IDO,2,K) = CC(1,K,1)+TAUR*CR2                                       
  101 CONTINUE                                                                  
      IF (IDO .EQ. 1) RETURN                                                    
      IDP2 = IDO+2                                                              
      DO 103 K=1,L1                                                             
         DO 102 I=3,IDO,2                                                       
            IC = IDP2-I                                                         
            DR2 = WA1(I-2)*CC(I-1,K,2)+WA1(I-1)*CC(I,K,2)                       
            DI2 = WA1(I-2)*CC(I,K,2)-WA1(I-1)*CC(I-1,K,2)                       
            DR3 = WA2(I-2)*CC(I-1,K,3)+WA2(I-1)*CC(I,K,3)                       
            DI3 = WA2(I-2)*CC(I,K,3)-WA2(I-1)*CC(I-1,K,3)                       
            CR2 = DR2+DR3                                                       
            CI2 = DI2+DI3                                                       
            CH(I-1,1,K) = CC(I-1,K,1)+CR2                                       
            CH(I,1,K) = CC(I,K,1)+CI2                                           
            TR2 = CC(I-1,K,1)+TAUR*CR2                                          
            TI2 = CC(I,K,1)+TAUR*CI2                                            
            TR3 = TAUI*(DI2-DI3)                                                
            TI3 = TAUI*(DR3-DR2)                                                
            CH(I-1,3,K) = TR2+TR3                                               
            CH(IC-1,2,K) = TR2-TR3                                              
            CH(I,3,K) = TI2+TI3                                                 
            CH(IC,2,K) = TI3-TI2                                                
  102    CONTINUE                                                               
  103 CONTINUE                                                                  
      RETURN                                                                    
      endsubroutine                                                                       
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE RADF4 (IDO,L1,CC,CH,WA1,WA2,WA3)                               
      DIMENSION       CC(IDO,L1,4)           ,CH(IDO,4,L1)           , &         
                      WA1(*)     ,WA2(*)     ,WA3(*)                            
      DATA HSQT2 /.7071067811865475/                                            
      DO 101 K=1,L1                                                             
         TR1 = CC(1,K,2)+CC(1,K,4)                                              
         TR2 = CC(1,K,1)+CC(1,K,3)                                              
         CH(1,1,K) = TR1+TR2                                                    
         CH(IDO,4,K) = TR2-TR1                                                  
         CH(IDO,2,K) = CC(1,K,1)-CC(1,K,3)                                      
         CH(1,3,K) = CC(1,K,4)-CC(1,K,2)                                        
  101 CONTINUE                                                                  
      IF (IDO-2) 107,105,102                                                    
  102 IDP2 = IDO+2                                                              
      DO 104 K=1,L1                                                             
         DO 103 I=3,IDO,2                                                       
            IC = IDP2-I                                                         
            CR2 = WA1(I-2)*CC(I-1,K,2)+WA1(I-1)*CC(I,K,2)                       
            CI2 = WA1(I-2)*CC(I,K,2)-WA1(I-1)*CC(I-1,K,2)                       
            CR3 = WA2(I-2)*CC(I-1,K,3)+WA2(I-1)*CC(I,K,3)                       
            CI3 = WA2(I-2)*CC(I,K,3)-WA2(I-1)*CC(I-1,K,3)                       
            CR4 = WA3(I-2)*CC(I-1,K,4)+WA3(I-1)*CC(I,K,4)                       
            CI4 = WA3(I-2)*CC(I,K,4)-WA3(I-1)*CC(I-1,K,4)                       
            TR1 = CR2+CR4                                                       
            TR4 = CR4-CR2                                                       
            TI1 = CI2+CI4                                                       
            TI4 = CI2-CI4                                                       
            TI2 = CC(I,K,1)+CI3                                                 
            TI3 = CC(I,K,1)-CI3                                                 
            TR2 = CC(I-1,K,1)+CR3                                               
            TR3 = CC(I-1,K,1)-CR3                                               
            CH(I-1,1,K) = TR1+TR2                                               
            CH(IC-1,4,K) = TR2-TR1                                              
            CH(I,1,K) = TI1+TI2                                                 
            CH(IC,4,K) = TI1-TI2                                                
            CH(I-1,3,K) = TI4+TR3                                               
            CH(IC-1,2,K) = TR3-TI4                                              
            CH(I,3,K) = TR4+TI3                                                 
            CH(IC,2,K) = TR4-TI3                                                
  103    CONTINUE                                                               
  104 CONTINUE                                                                  
      IF (MOD(IDO,2) .EQ. 1) RETURN                                             
  105 CONTINUE                                                                  
      DO 106 K=1,L1                                                             
         TI1 = -HSQT2*(CC(IDO,K,2)+CC(IDO,K,4))                                 
         TR1 = HSQT2*(CC(IDO,K,2)-CC(IDO,K,4))                                  
         CH(IDO,1,K) = TR1+CC(IDO,K,1)                                          
         CH(IDO,3,K) = CC(IDO,K,1)-TR1                                          
         CH(1,2,K) = TI1-CC(IDO,K,3)                                            
         CH(1,4,K) = TI1+CC(IDO,K,3)                                            
  106 CONTINUE                                                                  
  107 RETURN                                                                    
      endsubroutine                                                                       
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE RADF5 (IDO,L1,CC,CH,WA1,WA2,WA3,WA4)                           
      DIMENSION       CC(IDO,L1,5)           ,CH(IDO,5,L1)           , &          
                      WA1(*)     ,WA2(*)     ,WA3(*)     ,WA4(*)                
      DATA TR11,TI11,TR12,TI12 /.309016994374947,.951056516295154, &             
      -.809016994374947,.587785252292473/                                       
      DO 101 K=1,L1                                                             
         CR2 = CC(1,K,5)+CC(1,K,2)                                              
         CI5 = CC(1,K,5)-CC(1,K,2)                                              
         CR3 = CC(1,K,4)+CC(1,K,3)                                              
         CI4 = CC(1,K,4)-CC(1,K,3)                                              
         CH(1,1,K) = CC(1,K,1)+CR2+CR3                                          
         CH(IDO,2,K) = CC(1,K,1)+TR11*CR2+TR12*CR3                              
         CH(1,3,K) = TI11*CI5+TI12*CI4                                          
         CH(IDO,4,K) = CC(1,K,1)+TR12*CR2+TR11*CR3                              
         CH(1,5,K) = TI12*CI5-TI11*CI4                                          
  101 CONTINUE                                                                  
      IF (IDO .EQ. 1) RETURN                                                    
      IDP2 = IDO+2                                                              
      DO 103 K=1,L1                                                             
         DO 102 I=3,IDO,2                                                       
            IC = IDP2-I                                                         
            DR2 = WA1(I-2)*CC(I-1,K,2)+WA1(I-1)*CC(I,K,2)                       
            DI2 = WA1(I-2)*CC(I,K,2)-WA1(I-1)*CC(I-1,K,2)                       
            DR3 = WA2(I-2)*CC(I-1,K,3)+WA2(I-1)*CC(I,K,3)                       
            DI3 = WA2(I-2)*CC(I,K,3)-WA2(I-1)*CC(I-1,K,3)                       
            DR4 = WA3(I-2)*CC(I-1,K,4)+WA3(I-1)*CC(I,K,4)                       
            DI4 = WA3(I-2)*CC(I,K,4)-WA3(I-1)*CC(I-1,K,4)                       
            DR5 = WA4(I-2)*CC(I-1,K,5)+WA4(I-1)*CC(I,K,5)                       
            DI5 = WA4(I-2)*CC(I,K,5)-WA4(I-1)*CC(I-1,K,5)                       
            CR2 = DR2+DR5                                                       
            CI5 = DR5-DR2                                                       
            CR5 = DI2-DI5                                                       
            CI2 = DI2+DI5                                                       
            CR3 = DR3+DR4                                                       
            CI4 = DR4-DR3                                                       
            CR4 = DI3-DI4                                                       
            CI3 = DI3+DI4                                                       
            CH(I-1,1,K) = CC(I-1,K,1)+CR2+CR3                                   
            CH(I,1,K) = CC(I,K,1)+CI2+CI3                                       
            TR2 = CC(I-1,K,1)+TR11*CR2+TR12*CR3                                 
            TI2 = CC(I,K,1)+TR11*CI2+TR12*CI3                                   
            TR3 = CC(I-1,K,1)+TR12*CR2+TR11*CR3                                 
            TI3 = CC(I,K,1)+TR12*CI2+TR11*CI3                                   
            TR5 = TI11*CR5+TI12*CR4                                             
            TI5 = TI11*CI5+TI12*CI4                                             
            TR4 = TI12*CR5-TI11*CR4                                             
            TI4 = TI12*CI5-TI11*CI4                                             
            CH(I-1,3,K) = TR2+TR5                                               
            CH(IC-1,2,K) = TR2-TR5                                              
            CH(I,3,K) = TI2+TI5                                                 
            CH(IC,2,K) = TI5-TI2                                                
            CH(I-1,5,K) = TR3+TR4                                               
            CH(IC-1,4,K) = TR3-TR4                                              
            CH(I,5,K) = TI3+TI4                                                 
            CH(IC,4,K) = TI4-TI3                                                
  102    CONTINUE                                                               
  103 CONTINUE                                                                  
      RETURN                                                                    
      endsubroutine                                                                       
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE RADFG (IDO,IP,L1,IDL1,CC,C1,C2,CH,CH2,WA)                      
      DIMENSION       CH(IDO,L1,IP)          ,CC(IDO,IP,L1)          , &          
                      C1(IDO,L1,IP)          ,C2(IDL1,IP), &                    
                      CH2(IDL1,IP)           ,WA(*)                             
!      TPI = 2.0*3.141592653589793238462643383279502884197169399375105820                                                     
      ARG = TPI/FLOAT(IP)                                                       
      DCP = COS(ARG)                                                            
      DSP = SIN(ARG)                                                            
      IPPH = (IP+1)/2                                                           
      IPP2 = IP+2                                                               
      IDP2 = IDO+2                                                              
      NBD = (IDO-1)/2                                                           
      IF (IDO .EQ. 1) GO TO 119                                                 
      DO 101 IK=1,IDL1                                                          
         CH2(IK,1) = C2(IK,1)                                                   
  101 CONTINUE                                                                  
      DO 103 J=2,IP                                                             
         DO 102 K=1,L1                                                          
            CH(1,K,J) = C1(1,K,J)                                               
  102    CONTINUE                                                               
  103 CONTINUE                                                                  
      IF (NBD .GT. L1) GO TO 107                                                
      IS = -IDO                                                                 
      DO 106 J=2,IP                                                             
         IS = IS+IDO                                                            
         IDIJ = IS                                                              
         DO 105 I=3,IDO,2                                                       
            IDIJ = IDIJ+2                                                       
            DO 104 K=1,L1                                                       
               CH(I-1,K,J) = WA(IDIJ-1)*C1(I-1,K,J)+WA(IDIJ)*C1(I,K,J)          
               CH(I,K,J) = WA(IDIJ-1)*C1(I,K,J)-WA(IDIJ)*C1(I-1,K,J)            
  104       CONTINUE                                                            
  105    CONTINUE                                                               
  106 CONTINUE                                                                  
      GO TO 111                                                                 
  107 IS = -IDO                                                                 
      DO 110 J=2,IP                                                             
         IS = IS+IDO                                                            
         DO 109 K=1,L1                                                          
            IDIJ = IS                                                           
            DO 108 I=3,IDO,2                                                    
               IDIJ = IDIJ+2                                                    
               CH(I-1,K,J) = WA(IDIJ-1)*C1(I-1,K,J)+WA(IDIJ)*C1(I,K,J)          
               CH(I,K,J) = WA(IDIJ-1)*C1(I,K,J)-WA(IDIJ)*C1(I-1,K,J)            
  108       CONTINUE                                                            
  109    CONTINUE                                                               
  110 CONTINUE                                                                  
  111 IF (NBD .LT. L1) GO TO 115                                                
      DO 114 J=2,IPPH                                                           
         JC = IPP2-J                                                            
         DO 113 K=1,L1                                                          
            DO 112 I=3,IDO,2                                                    
               C1(I-1,K,J) = CH(I-1,K,J)+CH(I-1,K,JC)                           
               C1(I-1,K,JC) = CH(I,K,J)-CH(I,K,JC)                              
               C1(I,K,J) = CH(I,K,J)+CH(I,K,JC)                                 
               C1(I,K,JC) = CH(I-1,K,JC)-CH(I-1,K,J)                            
  112       CONTINUE                                                            
  113    CONTINUE                                                               
  114 CONTINUE                                                                  
      GO TO 121                                                                 
  115 DO 118 J=2,IPPH                                                           
         JC = IPP2-J                                                            
         DO 117 I=3,IDO,2                                                       
            DO 116 K=1,L1                                                       
               C1(I-1,K,J) = CH(I-1,K,J)+CH(I-1,K,JC)                           
               C1(I-1,K,JC) = CH(I,K,J)-CH(I,K,JC)                              
               C1(I,K,J) = CH(I,K,J)+CH(I,K,JC)                                 
               C1(I,K,JC) = CH(I-1,K,JC)-CH(I-1,K,J)                            
  116       CONTINUE                                                            
  117    CONTINUE                                                               
  118 CONTINUE                                                                  
      GO TO 121                                                                 
  119 DO 120 IK=1,IDL1                                                          
         C2(IK,1) = CH2(IK,1)                                                   
  120 CONTINUE                                                                  
  121 DO 123 J=2,IPPH                                                           
         JC = IPP2-J                                                            
         DO 122 K=1,L1                                                          
            C1(1,K,J) = CH(1,K,J)+CH(1,K,JC)                                    
            C1(1,K,JC) = CH(1,K,JC)-CH(1,K,J)                                   
  122    CONTINUE                                                               
  123 CONTINUE                                                                  
!                                                                               
      AR1 = 1.                                                                  
      AI1 = 0.                                                                  
      DO 127 L=2,IPPH                                                           
         LC = IPP2-L                                                            
         AR1H = DCP*AR1-DSP*AI1                                                 
         AI1 = DCP*AI1+DSP*AR1                                                  
         AR1 = AR1H                                                             
         DO 124 IK=1,IDL1                                                       
            CH2(IK,L) = C2(IK,1)+AR1*C2(IK,2)                                   
            CH2(IK,LC) = AI1*C2(IK,IP)                                          
  124    CONTINUE                                                               
         DC2 = AR1                                                              
         DS2 = AI1                                                              
         AR2 = AR1                                                              
         AI2 = AI1                                                              
         DO 126 J=3,IPPH                                                        
            JC = IPP2-J                                                         
            AR2H = DC2*AR2-DS2*AI2                                              
            AI2 = DC2*AI2+DS2*AR2                                               
            AR2 = AR2H                                                          
            DO 125 IK=1,IDL1                                                    
               CH2(IK,L) = CH2(IK,L)+AR2*C2(IK,J)                               
               CH2(IK,LC) = CH2(IK,LC)+AI2*C2(IK,JC)                            
  125       CONTINUE                                                            
  126    CONTINUE                                                               
  127 CONTINUE                                                                  
      DO 129 J=2,IPPH                                                           
         DO 128 IK=1,IDL1                                                       
            CH2(IK,1) = CH2(IK,1)+C2(IK,J)                                      
  128    CONTINUE                                                               
  129 CONTINUE                                                                  
!                                                                               
      IF (IDO .LT. L1) GO TO 132                                                
      DO 131 K=1,L1                                                             
         DO 130 I=1,IDO                                                         
            CC(I,1,K) = CH(I,K,1)                                               
  130    CONTINUE                                                               
  131 CONTINUE                                                                  
      GO TO 135                                                                 
  132 DO 134 I=1,IDO                                                            
         DO 133 K=1,L1                                                          
            CC(I,1,K) = CH(I,K,1)                                               
  133    CONTINUE                                                               
  134 CONTINUE                                                                  
  135 DO 137 J=2,IPPH                                                           
         JC = IPP2-J                                                            
         J2 = J+J                                                               
         DO 136 K=1,L1                                                          
            CC(IDO,J2-2,K) = CH(1,K,J)                                          
            CC(1,J2-1,K) = CH(1,K,JC)                                           
  136    CONTINUE                                                               
  137 CONTINUE                                                                  
      IF (IDO .EQ. 1) RETURN                                                    
      IF (NBD .LT. L1) GO TO 141                                                
      DO 140 J=2,IPPH                                                           
         JC = IPP2-J                                                            
         J2 = J+J                                                               
         DO 139 K=1,L1                                                          
            DO 138 I=3,IDO,2                                                    
               IC = IDP2-I                                                      
               CC(I-1,J2-1,K) = CH(I-1,K,J)+CH(I-1,K,JC)                        
               CC(IC-1,J2-2,K) = CH(I-1,K,J)-CH(I-1,K,JC)                       
               CC(I,J2-1,K) = CH(I,K,J)+CH(I,K,JC)                              
               CC(IC,J2-2,K) = CH(I,K,JC)-CH(I,K,J)                             
  138       CONTINUE                                                            
  139    CONTINUE                                                               
  140 CONTINUE                                                                  
      RETURN                                                                    
  141 DO 144 J=2,IPPH                                                           
         JC = IPP2-J                                                            
         J2 = J+J                                                               
         DO 143 I=3,IDO,2                                                       
            IC = IDP2-I                                                         
            DO 142 K=1,L1                                                       
               CC(I-1,J2-1,K) = CH(I-1,K,J)+CH(I-1,K,JC)                        
               CC(IC-1,J2-2,K) = CH(I-1,K,J)-CH(I-1,K,JC)                       
               CC(I,J2-1,K) = CH(I,K,J)+CH(I,K,JC)                              
               CC(IC,J2-2,K) = CH(I,K,JC)-CH(I,K,J)                             
  142       CONTINUE                                                            
  143    CONTINUE                                                               
  144 CONTINUE                                                                  
      RETURN                                                                    
      endsubroutine                                                                                                                                              
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!    subroutine rfftb(n,r,wsave,nsave)                                               
!                                                                              
!    SUBROUTINE RFFTB COMPUTES THE REAL PERODIC SEQUENCE FROM ITS              
!    FOURIER COEFFICIENTS (FOURIER SYNTHESIS). THE TRANSFORM IS DEFINED        
!    BELOW AT OUTPUT PARAMETER R.                                              
!                                                                              
!    INPUT PARAMETERS                                                          
!                                                                              
!    N       THE LENGTH OF THE ARRAY R TO BE TRANSFORMED.  THE METHOD          
!            IS MOST EFFICIENT WHEN N IS A PRODUCT OF SMALL PRIMES.            
!            N MAY CHANGE SO LONG AS DIFFERENT WORK ARRAYS ARE PROVIDED        
!                                                                              
!    R       A REAL ARRAY OF LENGTH N WHICH CONTAINS THE SEQUENCE              
!            TO BE TRANSFORMED                                                 
!                                                                              
!    WSAVE   A WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 2*N+15.           
!            IN THE PROGRAM THAT CALLS RFFTB. THE WSAVE ARRAY MUST BE          
!            INITIALIZED BY CALLING SUBROUTINE RFFTI(N,WSAVE) AND A            
!            DIFFERENT WSAVE ARRAY MUST BE USED FOR EACH DIFFERENT             
!            VALUE OF N. THIS INITIALIZATION DOES NOT HAVE TO BE               
!            REPEATED SO LONG AS N REMAINS UNCHANGED THUS SUBSEQUENT           
!            TRANSFORMS CAN BE OBTAINED FASTER THAN THE FIRST.                 
!            THE SAME WSAVE ARRAY CAN BE USED BY RFFTF AND RFFTB.              
!                                                                              
!                                                                              
!    OUTPUT PARAMETERS                                                         
!                                                                              
!    R       FOR N EVEN AND FOR I = 1,...,N                                    
!                                                                              
!                 R(I) = R(1)+(-1)**(I-1)*R(N)                                 
!                                                                              
!                      PLUS THE SUM FROM K=2 TO K=N/2 OF                       
!                                                                              
!                       2.*R(2*K-2)*COS((K-1)*(I-1)*2*PI/N)                    
!                                                                              
!                      -2.*R(2*K-1)*SIN((K-1)*(I-1)*2*PI/N)                    
!                                                                              
!            FOR N ODD AND FOR I = 1,...,N                                     
!                                                                              
!                 R(I) = R(1) PLUS THE SUM FROM K=2 TO K=(N+1)/2 OF            
!                                                                              
!                      2.*R(2*K-2)*COS((K-1)*(I-1)*2*PI/N)                     
!                                                                              
!                     -2.*R(2*K-1)*SIN((K-1)*(I-1)*2*PI/N)                     
!                                                                              
!     *****  NOTE                                                              
!                 THIS TRANSFORM IS UNNORMALIZED SINCE A CALL OF RFFTF         
!                 FOLLOWED BY A CALL OF RFFTB WILL MULTIPLY THE INPUT          
!                 SEQUENCE BY N.                                               
!                                                                              
!    WSAVE   CONTAINS RESULTS WHICH MUST NOT BE DESTROYED BETWEEN              
!            CALLS OF RFFTB OR RFFTF.                                          
!                                                                              
!                                                                              
      subroutine rfftb (n,r,wsave,nsave)               
      dimension r(*),wsave(*),nsave(*)                                      
!                                                                             
      if (n .eq. 1) return                                                      
      call rfftb1 (n,r,wsave,wsave(n+1),nsave)                           
      return                                                                    
                        endsubroutine rfftb                                         
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE RFFTB1 (N,C,CH,WA,IFAC)                                        
      DIMENSION       CH(*)      ,C(*)       ,WA(*)      ,IFAC(*)               
      NF = IFAC(2)                                                              
      NA = 0                                                                    
      L1 = 1                                                                    
      IW = 1                                                                    
      DO 116 K1=1,NF                                                            
         IP = IFAC(K1+2)                                                        
         L2 = IP*L1                                                             
         IDO = N/L2                                                             
         IDL1 = IDO*L1                                                          
         IF (IP .NE. 4) GO TO 103                                               
         IX2 = IW+IDO                                                           
         IX3 = IX2+IDO                                                          
         IF (NA .NE. 0) GO TO 101                                               
         CALL RADB4 (IDO,L1,C,CH,WA(IW),WA(IX2),WA(IX3))                        
         GO TO 102                                                              
  101    CALL RADB4 (IDO,L1,CH,C,WA(IW),WA(IX2),WA(IX3))                        
  102    NA = 1-NA                                                              
         GO TO 115                                                              
  103    IF (IP .NE. 2) GO TO 106                                               
         IF (NA .NE. 0) GO TO 104                                               
         CALL RADB2 (IDO,L1,C,CH,WA(IW))                                        
         GO TO 105                                                              
  104    CALL RADB2 (IDO,L1,CH,C,WA(IW))                                        
  105    NA = 1-NA                                                              
         GO TO 115                                                              
  106    IF (IP .NE. 3) GO TO 109                                               
         IX2 = IW+IDO                                                           
         IF (NA .NE. 0) GO TO 107                                               
         CALL RADB3 (IDO,L1,C,CH,WA(IW),WA(IX2))                                
         GO TO 108                                                              
  107    CALL RADB3 (IDO,L1,CH,C,WA(IW),WA(IX2))                                
  108    NA = 1-NA                                                              
         GO TO 115                                                              
  109    IF (IP .NE. 5) GO TO 112                                               
         IX2 = IW+IDO                                                           
         IX3 = IX2+IDO                                                          
         IX4 = IX3+IDO                                                          
         IF (NA .NE. 0) GO TO 110                                               
         CALL RADB5 (IDO,L1,C,CH,WA(IW),WA(IX2),WA(IX3),WA(IX4))                
         GO TO 111                                                              
  110    CALL RADB5 (IDO,L1,CH,C,WA(IW),WA(IX2),WA(IX3),WA(IX4))                
  111    NA = 1-NA                                                              
         GO TO 115                                                              
  112    IF (NA .NE. 0) GO TO 113                                               
         CALL RADBG (IDO,IP,L1,IDL1,C,C,C,CH,CH,WA(IW))                         
         GO TO 114                                                              
  113    CALL RADBG (IDO,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW))                        
  114    IF (IDO .EQ. 1) NA = 1-NA                                              
  115    L1 = L2                                                                
         IW = IW+(IP-1)*IDO                                                     
  116 CONTINUE                                                                  
      IF (NA .EQ. 0) RETURN                                                     
      DO 117 I=1,N                                                              
         C(I) = CH(I)                                                           
  117 CONTINUE                                                                  
      RETURN                                                                    
      endsubroutine                                                                       
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE RADB2 (IDO,L1,CC,CH,WA1)                                       
      DIMENSION       CC(IDO,2,L1)           ,CH(IDO,L1,2)           , &          
                      WA1(*)                                                    
      DO 101 K=1,L1                                                             
         CH(1,K,1) = CC(1,1,K)+CC(IDO,2,K)                                      
         CH(1,K,2) = CC(1,1,K)-CC(IDO,2,K)                                      
  101 CONTINUE                                                                  
      IF (IDO-2) 107,105,102                                                    
  102 IDP2 = IDO+2                                                              
      DO 104 K=1,L1                                                             
         DO 103 I=3,IDO,2                                                       
            IC = IDP2-I                                                         
            CH(I-1,K,1) = CC(I-1,1,K)+CC(IC-1,2,K)                              
            TR2 = CC(I-1,1,K)-CC(IC-1,2,K)                                      
            CH(I,K,1) = CC(I,1,K)-CC(IC,2,K)                                    
            TI2 = CC(I,1,K)+CC(IC,2,K)                                          
            CH(I-1,K,2) = WA1(I-2)*TR2-WA1(I-1)*TI2                             
            CH(I,K,2) = WA1(I-2)*TI2+WA1(I-1)*TR2                               
  103    CONTINUE                                                               
  104 CONTINUE                                                                  
      IF (MOD(IDO,2) .EQ. 1) RETURN                                             
  105 DO 106 K=1,L1                                                             
         CH(IDO,K,1) = CC(IDO,1,K)+CC(IDO,1,K)                                  
         CH(IDO,K,2) = -(CC(1,2,K)+CC(1,2,K))                                   
  106 CONTINUE                                                                  
  107 RETURN                                                                    
      endsubroutine                                                                       
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE RADB3 (IDO,L1,CC,CH,WA1,WA2)                                   
      DIMENSION       CC(IDO,3,L1)           ,CH(IDO,L1,3)           , &          
                      WA1(*)     ,WA2(*)                                        
      DATA TAUR,TAUI /-.5,.866025403784439/                                     
      DO 101 K=1,L1                                                             
         TR2 = CC(IDO,2,K)+CC(IDO,2,K)                                          
         CR2 = CC(1,1,K)+TAUR*TR2                                               
         CH(1,K,1) = CC(1,1,K)+TR2                                              
         CI3 = TAUI*(CC(1,3,K)+CC(1,3,K))                                       
         CH(1,K,2) = CR2-CI3                                                    
         CH(1,K,3) = CR2+CI3                                                    
  101 CONTINUE                                                                  
      IF (IDO .EQ. 1) RETURN                                                    
      IDP2 = IDO+2                                                              
      DO 103 K=1,L1                                                             
         DO 102 I=3,IDO,2                                                       
            IC = IDP2-I                                                         
            TR2 = CC(I-1,3,K)+CC(IC-1,2,K)                                      
            CR2 = CC(I-1,1,K)+TAUR*TR2                                          
            CH(I-1,K,1) = CC(I-1,1,K)+TR2                                       
            TI2 = CC(I,3,K)-CC(IC,2,K)                                          
            CI2 = CC(I,1,K)+TAUR*TI2                                            
            CH(I,K,1) = CC(I,1,K)+TI2                                           
            CR3 = TAUI*(CC(I-1,3,K)-CC(IC-1,2,K))                               
            CI3 = TAUI*(CC(I,3,K)+CC(IC,2,K))                                   
            DR2 = CR2-CI3                                                       
            DR3 = CR2+CI3                                                       
            DI2 = CI2+CR3                                                       
            DI3 = CI2-CR3                                                       
            CH(I-1,K,2) = WA1(I-2)*DR2-WA1(I-1)*DI2                             
            CH(I,K,2) = WA1(I-2)*DI2+WA1(I-1)*DR2                               
            CH(I-1,K,3) = WA2(I-2)*DR3-WA2(I-1)*DI3                             
            CH(I,K,3) = WA2(I-2)*DI3+WA2(I-1)*DR3                               
  102    CONTINUE                                                               
  103 CONTINUE                                                                  
      RETURN                                                                    
      endsubroutine                                                                       
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE RADB4 (IDO,L1,CC,CH,WA1,WA2,WA3)                               
      DIMENSION       CC(IDO,4,L1)           ,CH(IDO,L1,4)           , &          
                      WA1(*)     ,WA2(*)     ,WA3(*)                            
      DATA SQRT2 /1.414213562373095/                                            
      DO 101 K=1,L1                                                             
         TR1 = CC(1,1,K)-CC(IDO,4,K)                                            
         TR2 = CC(1,1,K)+CC(IDO,4,K)                                            
         TR3 = CC(IDO,2,K)+CC(IDO,2,K)                                          
         TR4 = CC(1,3,K)+CC(1,3,K)                                              
         CH(1,K,1) = TR2+TR3                                                    
         CH(1,K,2) = TR1-TR4                                                    
         CH(1,K,3) = TR2-TR3                                                    
         CH(1,K,4) = TR1+TR4                                                    
  101 CONTINUE                                                                  
      IF (IDO-2) 107,105,102                                                    
  102 IDP2 = IDO+2                                                              
      DO 104 K=1,L1                                                             
         DO 103 I=3,IDO,2                                                       
            IC = IDP2-I                                                         
            TI1 = CC(I,1,K)+CC(IC,4,K)                                          
            TI2 = CC(I,1,K)-CC(IC,4,K)                                          
            TI3 = CC(I,3,K)-CC(IC,2,K)                                          
            TR4 = CC(I,3,K)+CC(IC,2,K)                                          
            TR1 = CC(I-1,1,K)-CC(IC-1,4,K)                                      
            TR2 = CC(I-1,1,K)+CC(IC-1,4,K)                                      
            TI4 = CC(I-1,3,K)-CC(IC-1,2,K)                                      
            TR3 = CC(I-1,3,K)+CC(IC-1,2,K)                                      
            CH(I-1,K,1) = TR2+TR3                                               
            CR3 = TR2-TR3                                                       
            CH(I,K,1) = TI2+TI3                                                 
            CI3 = TI2-TI3                                                       
            CR2 = TR1-TR4                                                       
            CR4 = TR1+TR4                                                       
            CI2 = TI1+TI4                                                       
            CI4 = TI1-TI4                                                       
            CH(I-1,K,2) = WA1(I-2)*CR2-WA1(I-1)*CI2                             
            CH(I,K,2) = WA1(I-2)*CI2+WA1(I-1)*CR2                               
            CH(I-1,K,3) = WA2(I-2)*CR3-WA2(I-1)*CI3                             
            CH(I,K,3) = WA2(I-2)*CI3+WA2(I-1)*CR3                               
            CH(I-1,K,4) = WA3(I-2)*CR4-WA3(I-1)*CI4                             
            CH(I,K,4) = WA3(I-2)*CI4+WA3(I-1)*CR4                               
  103    CONTINUE                                                               
  104 CONTINUE                                                                  
      IF (MOD(IDO,2) .EQ. 1) RETURN                                             
  105 CONTINUE                                                                  
      DO 106 K=1,L1                                                             
         TI1 = CC(1,2,K)+CC(1,4,K)                                              
         TI2 = CC(1,4,K)-CC(1,2,K)                                              
         TR1 = CC(IDO,1,K)-CC(IDO,3,K)                                          
         TR2 = CC(IDO,1,K)+CC(IDO,3,K)                                          
         CH(IDO,K,1) = TR2+TR2                                                  
         CH(IDO,K,2) = SQRT2*(TR1-TI1)                                          
         CH(IDO,K,3) = TI2+TI2                                                  
         CH(IDO,K,4) = -SQRT2*(TR1+TI1)                                         
  106 CONTINUE                                                                  
  107 RETURN                                                                    
      endsubroutine                                                                       
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE RADB5 (IDO,L1,CC,CH,WA1,WA2,WA3,WA4)                           
      DIMENSION       CC(IDO,5,L1)           ,CH(IDO,L1,5)           , &          
                      WA1(*)     ,WA2(*)     ,WA3(*)     ,WA4(*)                
      DATA TR11,TI11,TR12,TI12 /.309016994374947,.951056516295154, &             
      -.809016994374947,.587785252292473/                                       
      DO 101 K=1,L1                                                             
         TI5 = CC(1,3,K)+CC(1,3,K)                                              
         TI4 = CC(1,5,K)+CC(1,5,K)                                              
         TR2 = CC(IDO,2,K)+CC(IDO,2,K)                                          
         TR3 = CC(IDO,4,K)+CC(IDO,4,K)                                          
         CH(1,K,1) = CC(1,1,K)+TR2+TR3                                          
         CR2 = CC(1,1,K)+TR11*TR2+TR12*TR3                                      
         CR3 = CC(1,1,K)+TR12*TR2+TR11*TR3                                      
         CI5 = TI11*TI5+TI12*TI4                                                
         CI4 = TI12*TI5-TI11*TI4                                                
         CH(1,K,2) = CR2-CI5                                                    
         CH(1,K,3) = CR3-CI4                                                    
         CH(1,K,4) = CR3+CI4                                                    
         CH(1,K,5) = CR2+CI5                                                    
  101 CONTINUE                                                                  
      IF (IDO .EQ. 1) RETURN                                                    
      IDP2 = IDO+2                                                              
      DO 103 K=1,L1                                                             
         DO 102 I=3,IDO,2                                                       
            IC = IDP2-I                                                         
            TI5 = CC(I,3,K)+CC(IC,2,K)                                          
            TI2 = CC(I,3,K)-CC(IC,2,K)                                          
            TI4 = CC(I,5,K)+CC(IC,4,K)                                          
            TI3 = CC(I,5,K)-CC(IC,4,K)                                          
            TR5 = CC(I-1,3,K)-CC(IC-1,2,K)                                      
            TR2 = CC(I-1,3,K)+CC(IC-1,2,K)                                      
            TR4 = CC(I-1,5,K)-CC(IC-1,4,K)                                      
            TR3 = CC(I-1,5,K)+CC(IC-1,4,K)                                      
            CH(I-1,K,1) = CC(I-1,1,K)+TR2+TR3                                   
            CH(I,K,1) = CC(I,1,K)+TI2+TI3                                       
            CR2 = CC(I-1,1,K)+TR11*TR2+TR12*TR3                                 
            CI2 = CC(I,1,K)+TR11*TI2+TR12*TI3                                   
            CR3 = CC(I-1,1,K)+TR12*TR2+TR11*TR3                                 
            CI3 = CC(I,1,K)+TR12*TI2+TR11*TI3                                   
            CR5 = TI11*TR5+TI12*TR4                                             
            CI5 = TI11*TI5+TI12*TI4                                             
            CR4 = TI12*TR5-TI11*TR4                                             
            CI4 = TI12*TI5-TI11*TI4                                             
            DR3 = CR3-CI4                                                       
            DR4 = CR3+CI4                                                       
            DI3 = CI3+CR4                                                       
            DI4 = CI3-CR4                                                       
            DR5 = CR2+CI5                                                       
            DR2 = CR2-CI5                                                       
            DI5 = CI2-CR5                                                       
            DI2 = CI2+CR5                                                       
            CH(I-1,K,2) = WA1(I-2)*DR2-WA1(I-1)*DI2                             
            CH(I,K,2) = WA1(I-2)*DI2+WA1(I-1)*DR2                               
            CH(I-1,K,3) = WA2(I-2)*DR3-WA2(I-1)*DI3                             
            CH(I,K,3) = WA2(I-2)*DI3+WA2(I-1)*DR3                               
            CH(I-1,K,4) = WA3(I-2)*DR4-WA3(I-1)*DI4                             
            CH(I,K,4) = WA3(I-2)*DI4+WA3(I-1)*DR4                               
            CH(I-1,K,5) = WA4(I-2)*DR5-WA4(I-1)*DI5                             
            CH(I,K,5) = WA4(I-2)*DI5+WA4(I-1)*DR5                               
  102    CONTINUE                                                               
  103 CONTINUE                                                                  
      RETURN                                                                    
      endsubroutine                                                                     
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE RADBG (IDO,IP,L1,IDL1,CC,C1,C2,CH,CH2,WA)                      
      DIMENSION       CH(IDO,L1,IP)          ,CC(IDO,IP,L1)          , &        
                      C1(IDO,L1,IP)          ,C2(IDL1,IP),  &                    
                      CH2(IDL1,IP)           ,WA(*)                             
!      TPI = 2.0*3.141592653589793238462643383279502884197169399375105820                                                     
      ARG = TPI/FLOAT(IP)                                                       
      DCP = COS(ARG)                                                            
      DSP = SIN(ARG)                                                            
      IDP2 = IDO+2                                                              
      NBD = (IDO-1)/2                                                           
      IPP2 = IP+2                                                               
      IPPH = (IP+1)/2                                                           
      IF (IDO .LT. L1) GO TO 103                                                
      DO 102 K=1,L1                                                             
         DO 101 I=1,IDO                                                         
            CH(I,K,1) = CC(I,1,K)                                               
  101    CONTINUE                                                               
  102 CONTINUE                                                                  
      GO TO 106                                                                 
  103 DO 105 I=1,IDO                                                            
         DO 104 K=1,L1                                                          
            CH(I,K,1) = CC(I,1,K)                                               
  104    CONTINUE                                                               
  105 CONTINUE                                                                  
  106 DO 108 J=2,IPPH                                                           
         JC = IPP2-J                                                            
         J2 = J+J                                                               
         DO 107 K=1,L1                                                          
            CH(1,K,J) = CC(IDO,J2-2,K)+CC(IDO,J2-2,K)                           
            CH(1,K,JC) = CC(1,J2-1,K)+CC(1,J2-1,K)                              
  107    CONTINUE                                                               
  108 CONTINUE                                                                  
      IF (IDO .EQ. 1) GO TO 116                                                 
      IF (NBD .LT. L1) GO TO 112                                                
      DO 111 J=2,IPPH                                                           
         JC = IPP2-J                                                            
         DO 110 K=1,L1                                                          
            DO 109 I=3,IDO,2                                                    
               IC = IDP2-I                                                      
               CH(I-1,K,J) = CC(I-1,2*J-1,K)+CC(IC-1,2*J-2,K)                   
               CH(I-1,K,JC) = CC(I-1,2*J-1,K)-CC(IC-1,2*J-2,K)                  
               CH(I,K,J) = CC(I,2*J-1,K)-CC(IC,2*J-2,K)                         
               CH(I,K,JC) = CC(I,2*J-1,K)+CC(IC,2*J-2,K)                        
  109       CONTINUE                                                            
  110    CONTINUE                                                               
  111 CONTINUE                                                                  
      GO TO 116                                                                 
  112 DO 115 J=2,IPPH                                                           
         JC = IPP2-J                                                            
         DO 114 I=3,IDO,2                                                       
            IC = IDP2-I                                                         
            DO 113 K=1,L1                                                       
               CH(I-1,K,J) = CC(I-1,2*J-1,K)+CC(IC-1,2*J-2,K)                   
               CH(I-1,K,JC) = CC(I-1,2*J-1,K)-CC(IC-1,2*J-2,K)                  
               CH(I,K,J) = CC(I,2*J-1,K)-CC(IC,2*J-2,K)                         
               CH(I,K,JC) = CC(I,2*J-1,K)+CC(IC,2*J-2,K)                        
  113       CONTINUE                                                            
  114    CONTINUE                                                               
  115 CONTINUE                                                                  
  116 AR1 = 1.                                                                  
      AI1 = 0.                                                                  
      DO 120 L=2,IPPH                                                           
         LC = IPP2-L                                                            
         AR1H = DCP*AR1-DSP*AI1                                                 
         AI1 = DCP*AI1+DSP*AR1                                                  
         AR1 = AR1H                                                             
         DO 117 IK=1,IDL1                                                       
            C2(IK,L) = CH2(IK,1)+AR1*CH2(IK,2)                                  
            C2(IK,LC) = AI1*CH2(IK,IP)                                          
  117    CONTINUE                                                               
         DC2 = AR1                                                              
         DS2 = AI1                                                              
         AR2 = AR1                                                              
         AI2 = AI1                                                              
         DO 119 J=3,IPPH                                                        
            JC = IPP2-J                                                         
            AR2H = DC2*AR2-DS2*AI2                                              
            AI2 = DC2*AI2+DS2*AR2                                               
            AR2 = AR2H                                                          
            DO 118 IK=1,IDL1                                                    
               C2(IK,L) = C2(IK,L)+AR2*CH2(IK,J)                                
               C2(IK,LC) = C2(IK,LC)+AI2*CH2(IK,JC)                             
  118       CONTINUE                                                            
  119    CONTINUE                                                               
  120 CONTINUE                                                                  
      DO 122 J=2,IPPH                                                           
         DO 121 IK=1,IDL1                                                       
            CH2(IK,1) = CH2(IK,1)+CH2(IK,J)                                     
  121    CONTINUE                                                               
  122 CONTINUE                                                                  
      DO 124 J=2,IPPH                                                           
         JC = IPP2-J                                                            
         DO 123 K=1,L1                                                          
            CH(1,K,J) = C1(1,K,J)-C1(1,K,JC)                                    
            CH(1,K,JC) = C1(1,K,J)+C1(1,K,JC)                                   
  123    CONTINUE                                                               
  124 CONTINUE                                                                  
      IF (IDO .EQ. 1) GO TO 132                                                 
      IF (NBD .LT. L1) GO TO 128                                                
      DO 127 J=2,IPPH                                                           
         JC = IPP2-J                                                            
         DO 126 K=1,L1                                                          
            DO 125 I=3,IDO,2                                                    
               CH(I-1,K,J) = C1(I-1,K,J)-C1(I,K,JC)                             
               CH(I-1,K,JC) = C1(I-1,K,J)+C1(I,K,JC)                            
               CH(I,K,J) = C1(I,K,J)+C1(I-1,K,JC)                               
               CH(I,K,JC) = C1(I,K,J)-C1(I-1,K,JC)                              
  125       CONTINUE                                                            
  126    CONTINUE                                                               
  127 CONTINUE                                                                  
      GO TO 132                                                                 
  128 DO 131 J=2,IPPH                                                           
         JC = IPP2-J                                                            
         DO 130 I=3,IDO,2                                                       
            DO 129 K=1,L1                                                       
               CH(I-1,K,J) = C1(I-1,K,J)-C1(I,K,JC)                             
               CH(I-1,K,JC) = C1(I-1,K,J)+C1(I,K,JC)                            
               CH(I,K,J) = C1(I,K,J)+C1(I-1,K,JC)                               
               CH(I,K,JC) = C1(I,K,J)-C1(I-1,K,JC)                              
  129       CONTINUE                                                            
  130    CONTINUE                                                               
  131 CONTINUE                                                                  
  132 CONTINUE                                                                  
      IF (IDO .EQ. 1) RETURN                                                    
      DO 133 IK=1,IDL1                                                          
         C2(IK,1) = CH2(IK,1)                                                   
  133 CONTINUE                                                                  
      DO 135 J=2,IP                                                             
         DO 134 K=1,L1                                                          
            C1(1,K,J) = CH(1,K,J)                                               
  134    CONTINUE                                                               
  135 CONTINUE                                                                  
      IF (NBD .GT. L1) GO TO 139                                                
      IS = -IDO                                                                 
      DO 138 J=2,IP                                                             
         IS = IS+IDO                                                            
         IDIJ = IS                                                              
         DO 137 I=3,IDO,2                                                       
            IDIJ = IDIJ+2                                                       
            DO 136 K=1,L1                                                       
               C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)-WA(IDIJ)*CH(I,K,J)          
               C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)+WA(IDIJ)*CH(I-1,K,J)            
  136       CONTINUE                                                            
  137    CONTINUE                                                               
  138 CONTINUE                                                                  
      GO TO 143                                                                 
  139 IS = -IDO                                                                 
      DO 142 J=2,IP                                                             
         IS = IS+IDO                                                            
         DO 141 K=1,L1                                                          
            IDIJ = IS                                                           
            DO 140 I=3,IDO,2                                                    
               IDIJ = IDIJ+2                                                    
               C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)-WA(IDIJ)*CH(I,K,J)          
               C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)+WA(IDIJ)*CH(I-1,K,J)            
  140       CONTINUE                                                            
  141    CONTINUE                                                               
  142 CONTINUE                                                                  
  143 RETURN                                                                    
      endsubroutine                                                                       
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
                        subroutine presmud &
(dlmd,dphd,sbd &
,nhsmud)
!-----------------------------------------------------------------------
!
      implicit none
!-----------------------------------------------------------------------

integer(kind=kint),dimension(jms:jme),intent(out) :: &
 nhsmud
!
real(kind=kfpt),intent(in) :: &
 dlmd &
,dphd &
,sbd
!
!-----------------------------------------------------------------------
!*** Local variables
!-----------------------------------------------------------------------
integer(kind=kint) :: &
 j
real(kind=kfpt) :: &
 dlm &
,dph &
,sb &
,tph
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
      dlm=dlmd*dtr
      dph=dphd*dtr
      sb=sbd*dtr
!-----------------------------------------------------------------------
      do j=jts_b2,jte_b2
        tph=sb+(j-jds-1)*dph
        nhsmud(j)=.99*dph/(dlm*cos(tph))
      enddo
!-----------------------------------------------------------------------
!
                        endsubroutine presmud
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
                        subroutine poavhn &
(i_start,i_end,j_start,j_end,km,hn,inpes,jnpes &
,use_allreduce)
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------

integer(kind=kint),intent(in) :: &
 i_start &
,i_end &
,j_start &
,j_end &
,inpes &
,jnpes &
,km 
!
real(kind=kfpt),dimension(i_start:i_end,j_start:j_end,km),intent(inout):: &
 hn
!
logical(kind=klog),intent(in) :: &
 use_allreduce
!
!-----------------------------------------------------------------------
!***  Local Variables
!-----------------------------------------------------------------------
!
integer(kind=kint) :: &
 i &
,ierr &
,irecv &
,l &
,loc_its_b1 &
,loc_ite_b2 &
,loc_min &
,loc_max &
,n &
,num_pes &
,pe
!
integer,dimension(mpi_status_size) :: jstat
!
real(kind=kfpt) :: &
 rcycle
!
real(kind=kfpt),dimension(ids+1:(ide-3)*km+1) :: &
 hn_glob
!
real(kind=kfpt),dimension((ite_b2-its_b1+1)*km):: &
 hn_loc
!
real(kind=kfpt),dimension(1:km) :: &
 an &
,an_g &
,as &
,as_g
 
integer(kind=kint) :: istat
logical(kind=klog) :: opened
logical(kind=klog),save :: sum_file_is_open=.false.
character(10) :: fstatus
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      num_pes=inpes*jnpes
!
!-----------------------------------------------------------------------
!
      rcycle=1./(ide-3)
      do l=1,km
        as(l)=0.
        an(l)=0.
      enddo
!
!-----------------------------------------------------------------------
!***  Sum the values along the row north of the southern boundary.
!-----------------------------------------------------------------------
!
      if(use_allreduce)then
!
!-----------------------------------------------------------------------
!***  Generate the global sum of the previous sums from each task
!***  using mpi_allreduce
!-----------------------------------------------------------------------
!
        if(s_bdy)then
          do l=1,km
            do i=its_b1,ite_b2
              as(l)=hn(i,jts+1,l)+as(l)
            enddo
          enddo
        endif
!
        call mpi_allreduce(as,as_g,km,mpi_real,mpi_sum,mpi_comm_comp &
                          ,irecv)
!
      else
!
!-----------------------------------------------------------------------
!***  Generate the global sum of the previous sums from each task
!***  using mpi_recv & mpi_send
!-----------------------------------------------------------------------
!
        as_g(1:km)=0

        if (mype==0) then

          do i=its_b1,ite_b2
            do l=1,km
              n=its_b1+l-1+km*(i-its_b1)
              hn_glob(n) = hn(i,jts+1,l)
            enddo
          enddo

          do pe=1,inpes-1
            loc_its_b1 = max(local_istart(pe),ids+1)
            loc_ite_b2 = min(local_iend(pe),ide-2)
            loc_min = (loc_its_b1-its_b1)*km+its_b1
            loc_max = loc_min + (loc_ite_b2-loc_its_b1+1)*km-1
            call mpi_recv(hn_glob(loc_min:loc_max),(loc_max-loc_min+1) &
                         ,mpi_real,pe,pe,mpi_comm_comp,jstat,ierr)
          end do

          do i=ids+1,ide-2
            do l=1,km
              n=ids+1+l-1+km*(i-ids-1)
              as_g(l) = hn_glob(n) + as_g(l)
            enddo
          end do

        else

          if(s_bdy)then

            n=0
            do i=its_b1,ite_b2
              do l=1,km
                n=n+1
                hn_loc(n)=hn(i,jts+1,l)
              enddo
            enddo

            call mpi_send(hn_loc(1:(ite_b2-its_b1+1)*km),(ite_b2-its_b1+1)*km &
                         ,mpi_real, 0, mype, mpi_comm_comp, ierr)

          endif
        endif

        call mpi_bcast (as_g,km,mpi_real,0,mpi_comm_comp,ierr)
!
      endif
!
!-----------------------------------------------------------------------
!***  Reset the array values in that same row to the global sum.
!-----------------------------------------------------------------------
!
      if(s_bdy)then
        do l=1,km
          as(l)=as_g(l)*rcycle
          do i=its,ite
            hn(i,jts+1,l)=as(l)
          enddo
        enddo
      endif
!
!-----------------------------------------------------------------------
!***  Now sum the values along the row south of the northern boundary.
!-----------------------------------------------------------------------
!
      if(use_allreduce)then
!
!-----------------------------------------------------------------------
!***  Generate the global sum of the previous sums from each task
!***  using mpi_allreduce
!-----------------------------------------------------------------------
!
        if(n_bdy)then
          do l=1,km
            do i=its_b1,ite_b2
              an(l)=hn(i,jte-1,l)+an(l)
            enddo
          enddo
        endif
!
        call mpi_allreduce(an,an_g,km,mpi_real,mpi_sum,mpi_comm_comp &
                          ,irecv)
!
      else
!
!-----------------------------------------------------------------------
!***  Generate the global sum of the previous sums from each task
!***  using mpi_recv & mpi_send
!-----------------------------------------------------------------------

        an_g(1:km)=0

        if (mype==num_pes-inpes) then

          do i=its_b1,ite_b2
            do l=1,km
              n=its_b1+l-1+km*(i-its_b1)
              hn_glob(n) = hn(i,jte-1,l)
            enddo
          enddo

          do pe=num_pes-inpes+1,num_pes-1
            loc_its_b1 = max(local_istart(pe),ids+1)
            loc_ite_b2 = min(local_iend(pe),ide-2)
            loc_min = (loc_its_b1-its_b1)*km+its_b1
            loc_max = loc_min + (loc_ite_b2-loc_its_b1+1)*km-1
            call mpi_recv(hn_glob(loc_min:loc_max),(loc_max-loc_min+1) &
                         ,mpi_real,pe,pe,mpi_comm_comp,jstat,ierr)
          end do

          do i=ids+1,ide-2
            do l=1,km
              n=ids+1+l-1+km*(i-ids-1)
              an_g(l) = hn_glob(n) + an_g(l)
            enddo
          end do

        else

          if(n_bdy)then

            n=0
            do i=its_b1,ite_b2
              do l=1,km
                n=n+1
                hn_loc(n)=hn(i,jte-1,l)
              enddo
            enddo

            call mpi_send(hn_loc(1:(ite_b2-its_b1+1)*km),(ite_b2-its_b1+1)*km &
                         ,mpi_real, num_pes-inpes, mype, mpi_comm_comp, ierr)

          endif
        endif

        call mpi_bcast (an_g,km,mpi_real,num_pes-inpes,mpi_comm_comp,ierr)


      endif
!
!-----------------------------------------------------------------------
!***  Reset the array values in that same row to the global sum.
!-----------------------------------------------------------------------
!
      if(n_bdy)then
        do l=1,km
          an(l)=an_g(l)*rcycle
          do i=its,ite
            hn(i,jte-1,l)=an(l)
          enddo
        enddo
      endif
!
!-----------------------------------------------------------------------
!
                        endsubroutine poavhn
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
                        subroutine swaphn &
(hn,i_start,i_end,j_start,j_end,km,inpes)
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------

integer(kind=kint),intent(in) :: &
 i_start &
,i_end &
,j_start &
,j_end &
,km &
,inpes
!
real(kind=kfpt),dimension(i_start:i_end,j_start:j_end,km),intent(inout):: &
 hn
!
!-----------------------------------------------------------------------
!***  Local Variables
!-----------------------------------------------------------------------
integer(kind=kint) :: &
 irecv &
,isend &
,j &
,l &
,length &
,ntask
!
integer(kind=kint),dimension(mpi_status_size) :: &
 jstat
!
real(kind=kfpt) :: &
 ave 

real(kind=kfpt),dimension(2,jts_h1:jte_h1,km) :: &
 eastx &
,westx
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Each task along the eastern and western "boundaries" needs to
!***  exchange with its counterpart on the opposite side.
!
!***  What is being said here is that i=1 (west) and i=ite-2 (east)
!***  are coincident, i=3 (west) and i=ite (east) are coincident,
!***  and since i=2 (west) and i=ite-1 (east) are coincident and
!***  in the true integration domain then they will equal their 
!***  mean.
!-----------------------------------------------------------------------
!
      length=2*(jte_h1-jts_h1+1)*km
!
      if(e_bdy)then
        do l=1,km
        do j=jts_h1,jte_h1
          eastx(1,j,l)=hn(ite-2,j,l)
          eastx(2,j,l)=hn(ite-1,j,l)
        enddo
        enddo
!
        ntask=mype-inpes+1
        call mpi_send(eastx,length,mpi_real,ntask,mype &
                     ,mpi_comm_comp,isend)
!
        call mpi_recv(westx,length,mpi_real,ntask,ntask &
                     ,mpi_comm_comp,jstat,irecv)
!
        do l=1,km
        do j=jts_h1,jte_h1
          ave=(westx(1,j,l)+hn(ite-1,j,l))*0.5
          hn(ite-1,j,l)=ave
          hn(ite,j,l)=westx(2,j,l) 
        enddo
        enddo
!-----------------------------------------------------------------------
      elseif(w_bdy)then
!
        ntask=mype+inpes-1  
        call mpi_recv(eastx,length,mpi_real,ntask,ntask &
                     ,mpi_comm_comp,jstat,irecv)
!
        do l=1,km
        do j=jts_h1,jte_h1
          westx(1,j,l)=hn(its+1,j,l)
          westx(2,j,l)=hn(its+2,j,l)
        enddo
        enddo
!
        call mpi_send(westx,length,mpi_real,ntask,mype &
                     ,mpi_comm_comp,isend)
!
        do l=1,km
        do j=jts_h1,jte_h1
          hn(its,j,l)=eastx(1,j,l)
          ave=(hn(its+1,j,l)+eastx(2,j,l))*0.5
          hn(its+1,j,l)=ave
        enddo
        enddo
!
!-----------------------------------------------------------------------
      endif
!-----------------------------------------------------------------------
!
                        endsubroutine swaphn
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
                        subroutine polehn &
(hn,i_start,i_end,j_start,j_end,km,inpes,jnpes)
!
!-----------------------------------------------------------------------
!***  Create polar mass point arrays holding all values around the
!***  southernmost and northernmost latitude circles in the
!***  true integration regions.
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------

integer(kind=kint),intent(in) :: &
 i_end &
,i_start &
,inpes &
,j_end &
,j_start &
,jnpes &
,km
!
real(kind=kfpt),dimension(ims:ime,jms:jme,km),intent(inout):: &
 hn
!
!-----------------------------------------------------------------------
!***  Local Variables
!-----------------------------------------------------------------------
!
integer(kind=kint) :: &
 i &
,l
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!***  First do the tasks at the south pole.
!-----------------------------------------------------------------------
!
      south_tasks: if(s_bdy)then
!
!------------------------------------------------------------
!
        do l=1,km
          do i=i_start,i_end
            hn(i,1,l)=hn(i,3,l)
          enddo
        enddo
!
!-----------------------------------------------------------------
!
      endif south_tasks
!
!-----------------------------------------------------------------
!***  Carry out the same procedure for the tasks at the north pole
!***  using three latitude circles.
!-----------------------------------------------------------------
!
      north_tasks: if(n_bdy)then
!
!--------------------------------------------------------------
!
        do l=1,km
          do i=i_start,i_end
            hn(i,jte,l)=hn(i,jte-2,l)
          enddo
        enddo
!
!-----------------------------------------------------------------------
!
      endif north_tasks
!
!-----------------------------------------------------------------------
!
                        end subroutine polehn
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
                        subroutine swapwn &
(wn,i_start,i_end,j_start,j_end,km,inpes)
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------

integer(kind=kint),intent(in) :: &
 i_start &
,i_end &
,j_start &
,j_end &
,km &
,inpes
!
real(kind=kfpt),dimension(i_start:i_end,j_start:j_end,km),intent(inout):: &
 wn
!
!-----------------------------------------------------------------------
!***  Local Variables
!-----------------------------------------------------------------------
integer(kind=kint) :: &
 irecv &
,isend &
,j &
,l &
,length_e &
,length_w &
,ntask
!
integer(kind=kint),dimension(mpi_status_size) :: &
 jstat
!
real(kind=kfpt),dimension(jts_h1:jte_h1,km) :: &
 eastx 
!
real(kind=kfpt),dimension(2,jts_h1:jte_h1,km) :: &
 westx
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Each task along the eastern and western "boundaries" needs to
!***  exchange with its counterpart on the opposite side.
!
!***  What is being said here is that i=1 (west) and i=ite-2 (east)
!***  are coincident, i=2 (west) and i=ite-1 (east) are coincident,
!***  and i=3 (west) and i=ite (east) are coincident.
!-----------------------------------------------------------------------
!
      length_e=(jte_h1-jts_h1+1)*km
      length_w=2*length_e
!
      if(e_bdy)then
        do l=1,km
        do j=jts_h1,jte_h1
          eastx(j,l)=wn(ite-2,j,l)
        enddo
        enddo
!
        ntask=mype-inpes+1
        call mpi_send(eastx,length_e,mpi_real,ntask,mype &
                     ,mpi_comm_comp,isend)
!
        call mpi_recv(westx,length_w,mpi_real,ntask,ntask &
                     ,mpi_comm_comp,jstat,irecv)
!
        do l=1,km
        do j=jts_h1,jte_h1
          wn(ite-1,j,l)=westx(1,j,l)
          wn(ite,j,l)=westx(2,j,l)
        enddo
        enddo
!-----------------------------------------------------------------------
      elseif(w_bdy)then
!
        ntask=mype+inpes-1
        call mpi_recv(eastx,length_e,mpi_real,ntask,ntask &
                     ,mpi_comm_comp,jstat,irecv)
!
        do l=1,km
        do j=jts_h1,jte_h1
          wn(its,j,l)=eastx(j,l)
        enddo
        enddo
!
        do l=1,km
        do j=jts_h1,jte_h1
          westx(1,j,l)=wn(its+1,j,l)
          westx(2,j,l)=wn(its+2,j,l)
        enddo
        enddo
!
        call mpi_send(westx,length_w,mpi_real,ntask,mype &
                     ,mpi_comm_comp,isend)
!
!-----------------------------------------------------------------------
      endif
!-----------------------------------------------------------------------
!
                        endsubroutine swapwn
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
                        subroutine polewn &
(u,v,i_start,i_end,j_start,j_end,km,inpes,jnpes)
!
!-----------------------------------------------------------------------
!***  Create polar wind arrays holding all values around the
!***  southernmost and northernmost latitude circles in the
!***  true integration regions.
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------

integer(kind=kint),intent(in) :: &
 i_end &
,i_start &
,inpes &
,j_end &
,j_start &
,jnpes &
,km
!
real(kind=kfpt),dimension(ims:ime,jms:jme,km),intent(inout):: &
 u &
,v
!
!-----------------------------------------------------------------------
!***  Local Variables
!-----------------------------------------------------------------------
!
integer(kind=kint) :: &
 i &
,l
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!***  First do the tasks at the south pole.
!-----------------------------------------------------------------------
!
      south_tasks: if(s_bdy)then
!
!--------------------------------------------------------------
!
        do l=1,km
          do i=i_start,i_end
            u(i,1,l)=-u(i,2,l)
            v(i,1,l)=-v(i,2,l)
          enddo
        enddo
!
!-----------------------------------------------------------------
!
      endif south_tasks
!
!-----------------------------------------------------------------
!***  Carry out the same procedure for the tasks at the north pole
!***  using three latitude circles.
!-----------------------------------------------------------------
!
      north_tasks: if(n_bdy)then
!
!--------------------------------------------------------------------
!
        do l=1,km
          do i=i_start,i_end
            u(i,jte-1,l)=-u(i,jte-2,l)
            v(i,jte-1,l)=-v(i,jte-2,l)
            u(i,jte  ,l)=-u(i,jte-2,l)
            v(i,jte  ,l)=-v(i,jte-2,l)
          enddo
        enddo
!
!-----------------------------------------------------------------------
!
      endif north_tasks
!
!-----------------------------------------------------------------------
!
                        end subroutine polewn
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
                        subroutine read_bc &
(lm,lnsh,lnsv,ntsd,dt &
,runbc,idatbc,ihrstbc,tboco &
,nvars_bc_2d_h,nvars_bc_3d_h,nvars_bc_4d_h &
,nvars_bc_2d_v,nvars_bc_3d_v &
,bnd_vars_h &
,bnd_vars_v &
,n_bc_3d_h )
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  This routine reads in boundary update data for the single domain
!***  in the forecast or else for the uppermost parent if there are
!***  nests.  The data is in external files generated during the
!***  pre-processing.
!
!***  This routine is used for nest boundaries as well if digital
!***  filtering is being used.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
implicit none
!
!------------------------
!***  Argument variables
!------------------------

integer(kind=kint),intent(in):: &
 lm  &                       ! total # of levels
,lnsh &                      ! # of boundary h lines for bc in reg. setup
,lnsv &                      ! # of boundary v lines for bc in reg. setup
,ntsd                        ! current timestep

integer(kind=kint),intent(in) :: &
 nvars_bc_2d_h &             ! # of 2-d h-pt boundary variables
,nvars_bc_3d_h &             ! # of 3-d h-pt boundary variables
,nvars_bc_4d_h &             ! # of 4-d h-pt boundary variables
,nvars_bc_2d_v &             ! # of 2-d v-pt boundary variables
,nvars_bc_3d_v               ! # of 3-d v-pt boundary variables

integer(kind=kint),intent(out):: &
 ihrstbc                     ! boundary conditions starting time

integer(kind=kint),dimension(1:3),intent(in):: &
 n_bc_3d_h                   ! order in which to store domain #1's BC vbls

integer(kind=kint),dimension(1:3),intent(out):: &
 idatbc                      ! date of boundary data, day, month, year

real(kind=kfpt),intent(in):: &
 dt                          ! dynamics time step

real(kind=kfpt),intent(out):: &
 tboco                       ! boundary conditions interval, hours

type(bc_h_all),intent(inout) :: &
 bnd_vars_h                  ! boundary variables on h points

type(bc_v_all),intent(inout) :: &
 bnd_vars_v                  ! boundary variables on v points

logical(kind=klog),intent(out) :: runbc

!---------------------
!***  Local variables
!---------------------

integer(kind=kint):: &
 i &                         ! index in x direction
,i_hi &                      ! max i loop limit (cannot be > ide)
,i_lo &                      ! min i loop limit (cannot be < ids)
,ihr &                       ! current hour
,ihrbc &                     ! current hour for bc's
,istat &                     ! return status
,iunit &                     ! file unit number
,j_hi &                      ! max j loop limit (cannot be > jde)
,j_lo &                      ! min j loop limit (cannot be < jds)
,l &
,lbnd &
,n1 &                        ! dimension 1 of working arrays
,n2 &                        ! dimension 2 of working arrays
,n3 &                        ! dimension 3 of working arrays
,n4 &                        ! dimension 4 of working arrays
,n5 &                        ! dimension 5 of working arrays
,nl &
,nv &
,nvx &
,ubnd

real(kind=kfpt),dimension(1:lnsh,jds:jde,1:2) :: global_2d_h_east &
                                                ,global_2d_h_west

real(kind=kfpt),dimension(ids:ide,1:lnsh,1:2) :: global_2d_h_north &
                                                ,global_2d_h_south
 
real(kind=kfpt),dimension(1:lnsh,jds:jde,1:lm,1:2) :: global_3d_h_east &
                                                     ,global_3d_h_west

real(kind=kfpt),dimension(ids:ide,1:lnsh,1:lm,1:2) :: global_3d_h_north &
                                                     ,global_3d_h_south

real(kind=kfpt),dimension(1:lnsv,jds:jde,1:2) :: global_2d_v_east &
                                                ,global_2d_v_west

real(kind=kfpt),dimension(ids:ide,1:lnsv,1:2) :: global_2d_v_north &
                                                ,global_2d_v_south
 
real(kind=kfpt),dimension(1:lnsv,jds:jde,1:lm,1:2) :: global_3d_v_east &
                                                     ,global_3d_v_west

real(kind=kfpt),dimension(ids:ide,1:lnsv,1:lm,1:2) :: global_3d_v_north &
                                                     ,global_3d_v_south

character(64) :: infile
logical(kind=klog) :: opened
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Because subdomains that lie along the global domain boundary
!***  may have haloes that extend beyond the global limits, create
!***  limits here that keep loops from reaching beyond those
!***  global limits.
!-----------------------------------------------------------------------
!
      i_lo=max(ims,ids)
      i_hi=min(ime,ide)
      j_lo=max(jms,jds)
      j_hi=min(jme,jde)
!
!-----------------------------------------------------------------------
!***  Read in boundary variables and arrays. 
!-----------------------------------------------------------------------
!
      ihr=nint(ntsd*abs(dt)/3600.)

      ihrbc=ihr
      write(infile,'(a,i4.4)')'boco.',ihrbc
!
      select_unit: do l=51,59
        inquire(l,opened=opened)
        if(.not.opened)then
          iunit=l
          exit select_unit
        endif
      enddo select_unit
!
      open(unit=iunit,file=infile,status='OLD' &
          ,form='UNFORMATTED',iostat=istat)

!     if (mype_share == 0) write(0,*) 'reading from boco file: ', trim(infile)
!
!     write(6,*) 'DEBUG-GT, reading from boco file: ', trim(infile), 'unit number ', iunit
!     write(6,*) 'DEBUG-GT, domain size is ide,jde,lm = ', ide, jde, lm
!     write(6,*) 'DEBUG-GT, reading logical, 4 INTS, 1 FLOAT | 21 bytes'
      rewind(iunit)
      read(iunit)runbc,idatbc,ihrstbc,tboco
!     write(6,*) 'DEBUG-GT, found: runbc,idatbc,ihrstbc,tboco', runbc,idatbc,ihrstbc,tboco

!     write(6,*) 'DEBUG-GT, assortment of dimensions in subsequent read statements;  i_lo, i_hi, j_lo, j_hi, lm, lnsh, lnsv:', i_lo, i_hi, j_lo, j_hi, lm, lnsh, lnsv
!
!-----------------------------------------------------------------------
!***  If there are 2-d h-pt boundary variables then loop through them
!***  reading them from the input file.
!-----------------------------------------------------------------------
!
      if(nvars_bc_2d_h>0)then
!
        do nv=1,nvars_bc_2d_h
!     write(6,*) 'DEBUG-GT, reading PD.1, 2 arrays with N elements to read: ', (ide-ids+1)*lnsh*2, (ide-ids+1)*lnsh*2*4, ' bytes'
!     write(6,*) 'DEBUG-GT, reading PD.2, 2 arrays with M elements to read: ', (jde-jds+1)*lnsh*2, (jde-jds+1)*lnsh*2*4, ' bytes'
      read(iunit)global_2d_h_south,global_2d_h_north,global_2d_h_west,global_2d_h_east
!     write(6,*) 'DEBUG-GT, read arrays pdbs_g,pdbn_g,pdbw_g,pdbe_g'
!
          do n3=1,2
            do n2=1,lnsh
            do n1=i_lo,i_hi
              bnd_vars_h%var_2d(nv)%south(n1,n2,n3)=global_2d_h_south(n1,n2,n3)
              bnd_vars_h%var_2d(nv)%north(n1,n2,n3)=global_2d_h_north(n1,n2,n3)
            enddo
            enddo
!
            do n2=j_lo,j_hi
            do n1=1,lnsh
              bnd_vars_h%var_2d(nv)%west(n1,n2,n3)=global_2d_h_west(n1,n2,n3)
              bnd_vars_h%var_2d(nv)%east(n1,n2,n3)=global_2d_h_east(n1,n2,n3)
            enddo
            enddo
          enddo
!
        enddo
!
      endif
!
!-----------------------------------------------------------------------
!***  If there are 3-d h-pt boundary variables then loop through them
!***  reading them from the input file.
!***  IMPORTANT:  See the note in BUILD_BC_BUNDLE regarding the 
!***              handling of order in which these 3-D H-pt BC
!***              variables are stored.  If there are nests then
!***              the order must agree with these variables'
!***              order in the nests.txt file.
!-----------------------------------------------------------------------
!
      if(nvars_bc_3d_h>0)then
!
        do nv=1,nvars_bc_3d_h
!
          nvx=n_bc_3d_h(nv)                                                !<-- Order of domain #1's T,Q,CW storage in bndry object
!
          read(iunit)global_3d_h_south,global_3d_h_north,global_3d_h_west,global_3d_h_east
!
          do n4=1,2
            do n3=1,lm
              do n2=1,lnsh
              do n1=i_lo,i_hi
                bnd_vars_h%var_3d(nvx)%south(n1,n2,n3,n4)=global_3d_h_south(n1,n2,n3,n4)
                bnd_vars_h%var_3d(nvx)%north(n1,n2,n3,n4)=global_3d_h_north(n1,n2,n3,n4)
              enddo
              enddo
!
              do n2=j_lo,j_hi
              do n1=1,lnsh
                bnd_vars_h%var_3d(nvx)%west(n1,n2,n3,n4)=global_3d_h_west(n1,n2,n3,n4)
                bnd_vars_h%var_3d(nvx)%east(n1,n2,n3,n4)=global_3d_h_east(n1,n2,n3,n4)
              enddo
              enddo
            enddo
          enddo
!
        enddo
!
      endif
!
!-----------------------------------------------------------------------
!***  If there are 4-d h-pt boundary variables then loop through them
!***  reading them from the input file.
!-----------------------------------------------------------------------
!
      if(nvars_bc_4d_h>0)then
!
        do nv=1,nvars_bc_4d_h
!
          lbnd=lbound(bnd_vars_h%var_4d(nv)%full_var,4)
          ubnd=ubound(bnd_vars_h%var_4d(nv)%full_var,4)
          do nl=lbnd,ubnd
            read(iunit)global_3d_h_south,global_3d_h_north,global_3d_h_west,global_3d_h_east
!
            do n4=1,2
              do n3=1,lm
                do n2=1,lnsh
                do n1=i_lo,i_hi
                  bnd_vars_h%var_4d(nv)%south(n1,n2,n3,n4,nl)=global_3d_h_south(n1,n2,n3,n4)
                  bnd_vars_h%var_4d(nv)%north(n1,n2,n3,n4,nl)=global_3d_h_north(n1,n2,n3,n4)
                enddo
                enddo
!
                do n2=j_lo,j_hi
                do n1=1,lnsh
                  bnd_vars_h%var_4d(nv)%west(n1,n2,n3,n4,nl)=global_3d_h_west(n1,n2,n3,n4)
                  bnd_vars_h%var_4d(nv)%east(n1,n2,n3,n4,nl)=global_3d_h_east(n1,n2,n3,n4)
                enddo
                enddo
              enddo
            enddo
!
          enddo
!
        enddo
!
      endif
!
!-----------------------------------------------------------------------
!***  If there are 2-d v-pt boundary variables then loop through them
!***  reading them from the input file.
!-----------------------------------------------------------------------
!
      if(nvars_bc_2d_v>0)then
!
        do nv=1,nvars_bc_2d_v
!
          read(iunit)global_2d_v_south,global_2d_v_north,global_2d_v_west,global_2d_v_east
!
          do n3=1,2
            do n2=1,lnsv
            do n1=i_lo,i_hi
              bnd_vars_v%var_2d(nv)%south(n1,n2,n3)=global_2d_v_south(n1,n2,n3)
              bnd_vars_v%var_2d(nv)%north(n1,n2,n3)=global_2d_v_north(n1,n2,n3)
            enddo
            enddo
!
            do n2=j_lo,j_hi
            do n1=1,lnsv
              bnd_vars_v%var_2d(nv)%west(n1,n2,n3)=global_2d_v_west(n1,n2,n3)
              bnd_vars_v%var_2d(nv)%east(n1,n2,n3)=global_2d_v_east(n1,n2,n3)
            enddo
            enddo
          enddo
!
        enddo
!
      endif
!
!-----------------------------------------------------------------------
!***  If there are 3-d v-pt boundary variables then loop through them
!***  reading them from the input file.
!-----------------------------------------------------------------------
!
      if(nvars_bc_3d_v>0)then
!
        do nv=1,nvars_bc_3d_v
!
          read(iunit)global_3d_v_south,global_3d_v_north,global_3d_v_west,global_3d_v_east
!
          do n4=1,2
          do n3=1,lm
            do n2=1,lnsv
            do n1=i_lo,i_hi
              bnd_vars_v%var_3d(nv)%south(n1,n2,n3,n4)=global_3d_v_south(n1,n2,n3,n4)
              bnd_vars_v%var_3d(nv)%north(n1,n2,n3,n4)=global_3d_v_north(n1,n2,n3,n4)
            enddo
            enddo
!
            do n2=j_lo,j_hi
            do n1=1,lnsv
              bnd_vars_v%var_3d(nv)%west(n1,n2,n3,n4)=global_3d_v_west(n1,n2,n3,n4)
              bnd_vars_v%var_3d(nv)%east(n1,n2,n3,n4)=global_3d_v_east(n1,n2,n3,n4)
            enddo
            enddo
          enddo
          enddo
!
        enddo
!
      endif
!
!-----------------------------------------------------------------------
!
      close(iunit)
!
!-----------------------------------------------------------------------

                        end subroutine read_bc

!-----------------------------------------------------------------------
                        subroutine write_bc &
(lm,lnsh,lnsv,ntsd,dt &
,runbc,tboco &
,nvars_bc_2d_h &
,nvars_bc_3d_h &
,nvars_bc_4d_h &
,nvars_bc_2d_v &
,nvars_bc_3d_v &
,bnd_vars_h &
,bnd_vars_v &
,recomp_tend)
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
implicit none
!
!-----------------------------------------------------------------------
!
!include 'kind.inc'
!
!-----------------------------------------------------------------------
integer(kind=kint),intent(in):: &
 lm  &                       ! total # of levels
,lnsh &                      ! # of boundary h lines for bc in reg. setup
,lnsv &                      ! # of boundary v lines for bc in reg. setup
,ntsd                        ! current timestep

integer(kind=kint),intent(in) :: &
 nvars_bc_2d_h &
,nvars_bc_3d_h &
,nvars_bc_4d_h &
,nvars_bc_2d_v &
,nvars_bc_3d_v

real(kind=kfpt),intent(in):: &
 dt                          ! dynamics time step

real(kind=kfpt),intent(in):: &
 tboco                       ! boundary conditions interval

logical, intent(in) :: recomp_tend

logical(kind=klog),intent(in) :: runbc

type(bc_h_all),intent(inout) :: &
 bnd_vars_h

type(bc_v_all),intent(inout) :: &
 bnd_vars_v

!-----------------------------------------------------------------------
!---local variables-----------------------------------------------------
!-----------------------------------------------------------------------
integer(kind=kint):: &
 i &                         ! index in x direction
,i_hi &                      ! max i loop limit (cannot be > ide)
,i_lo &                      ! min i loop limit (cannot be < ids)
,ihr &                       ! current hour
,ihrbc &                     ! current hour for bc's
,istat &                     ! return status
,iunit &                     ! file unit number
,j_hi &                      ! max j loop limit (cannot be > jde)
,j_lo &                      ! min j loop limit (cannot be < jds)
,l &
,lbnd &
,n1 &                        ! dimension 1 of working arrays
,n2 &                        ! dimension 2 of working arrays
,n3 &                        ! dimension 3 of working arrays
,n4 &                        ! dimension 4 of working arrays
,n5 &                        ! dimension 5 of working arrays
,nv &
,ubnd

real(kind=kfpt) :: &
 r_tboco

real(kind=kfpt),dimension(:,:,:,:),allocatable :: &
 targ_2d_h_e &
,targ_2d_h_n &
,targ_2d_h_s &
,targ_2d_h_w

real(kind=kfpt),dimension(:,:,:,:,:),allocatable :: &
 targ_3d_h_e &
,targ_3d_h_n &
,targ_3d_h_s &
,targ_3d_h_w

real(kind=kfpt),dimension(:,:,:,:),allocatable :: &
 targ_2d_v_e &
,targ_2d_v_n &
,targ_2d_v_s &
,targ_2d_v_w

real(kind=kfpt),dimension(:,:,:,:,:),allocatable :: &
 targ_3d_v_e &
,targ_3d_v_n &
,targ_3d_v_s &
,targ_3d_v_w

type(filt_4d),dimension(:),allocatable :: &
 targ_4d_h_e &
,targ_4d_h_n &
,targ_4d_h_s &
,targ_4d_h_w

!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      r_tboco=1./tboco
!
!-----------------------------------------------------------------------
!***  Because subdomains that lie along the global domain boundary
!***  may have haloes that extend beyond the global limits, create
!***  limits here that keep loops from reaching beyond those
!***  global limits.
!-----------------------------------------------------------------------
!
      i_lo=max(ims,ids)
      i_hi=min(ime,ide)
      j_lo=max(jms,jds)
      j_hi=min(jme,jde)
!
!     if (mype_share == 0) write(0,*) 'inside write_bc with recomp_tend: ', recomp_tend
      if (recomp_tend) then

!-----------------------------------------------------------------------
!***  Is this task subdomain on the full domain's north side?
!-----------------------------------------------------------------------

        IF (n_bdy) THEN

          if(nvars_bc_2d_h>0)then
            allocate(targ_2d_h_n(ims:ime,1:lnsh,1,1:nvars_bc_2d_h))
            do nv=1,nvars_bc_2d_h
              do n3=1,1
              do n2=1,lnsh
              do n1=i_lo,i_hi
                targ_2d_h_n(n1,n2,n3,nv)=bnd_vars_h%var_2d(nv)%north(n1,n2,1) &
                                        +tboco*bnd_vars_h%var_2d(nv)%north(n1,n2,2)
              enddo
              enddo
              enddo
            enddo
          endif

          if(nvars_bc_3d_h>0)then
            allocate(targ_3d_h_n(ims:ime,1:lnsh,1:lm,1,1:nvars_bc_3d_h))
            do nv=1,nvars_bc_3d_h
              do n4=1,1
              do n3=1,lm
              do n2=1,lnsh
              do n1=i_lo,i_hi
                targ_3d_h_n(n1,n2,n3,n4,nv)=bnd_vars_h%var_3d(nv)%north(n1,n2,n3,1) &
                                           +tboco*bnd_vars_h%var_3d(nv)%north(n1,n2,n3,2)
              enddo
              enddo
              enddo
              enddo
            enddo
          endif

          if(nvars_bc_4d_h>0)then
            allocate(targ_4d_h_n(1:nvars_bc_4d_h))
            do nv=1,nvars_bc_4d_h
              lbnd=lbound(bnd_vars_h%var_4d(nv)%full_var,4)
              ubnd=ubound(bnd_vars_h%var_4d(nv)%full_var,4)
              allocate(targ_4d_h_n(nv)%base(ims:ime,1:lnsh,1:lm,1,lbnd:ubnd))
              do n5=lbnd,ubnd
              do n4=1,1
              do n3=1,lm
              do n2=1,lnsh
              do n1=i_lo,i_hi
                targ_4d_h_n(nv)%base(n1,n2,n3,n4,n5)=bnd_vars_h%var_4d(nv)%north(n1,n2,n3,1,n5) &
                                                    +tboco*bnd_vars_h%var_4d(nv)%north(n1,n2,n3,2,n5)
              enddo
              enddo
              enddo
              enddo
              enddo
            enddo
          endif

          if(nvars_bc_2d_v>0)then
            allocate(targ_2d_v_n(ims:ime,1:lnsv,1,nvars_bc_2d_v))
            do nv=1,nvars_bc_2d_v
              do n3=1,1
              do n2=1,lnsv
              do n1=i_lo,i_hi
                targ_2d_v_n(n1,n2,n3,nv)=bnd_vars_v%var_2d(nv)%north(n1,n2,1) &
                                        +tboco*bnd_vars_v%var_2d(nv)%north(n1,n2,2)
              enddo
              enddo
              enddo
            enddo
          endif

          if(nvars_bc_3d_v>0)then
            allocate(targ_3d_v_n(ims:ime,1:lnsv,1:lm,1,nvars_bc_3d_v))
            do nv=1,nvars_bc_3d_v
              do n4=1,1
              do n3=1,lm
              do n2=1,lnsv
              do n1=i_lo,i_hi
                targ_3d_v_n(n1,n2,n3,n4,nv)=bnd_vars_v%var_3d(nv)%north(n1,n2,n3,1) &
                                           +tboco*bnd_vars_v%var_3d(nv)%north(n1,n2,n3,2)
              enddo
              enddo
              enddo
              enddo
            enddo
          endif

	ENDIF ! N_BDY

!-----------------------------------------------------------------------
!***  Is this task subdomain on the full domain's south side?
!-----------------------------------------------------------------------

        IF (s_bdy) THEN

          if(nvars_bc_2d_h>0)then
            allocate(targ_2d_h_s(ims:ime,1:lnsh,1,1:nvars_bc_2d_h))
            do nv=1,nvars_bc_2d_h
              do n3=1,1
              do n2=1,lnsh
              do n1=i_lo,i_hi
                targ_2d_h_s(n1,n2,n3,nv)=bnd_vars_h%var_2d(nv)%south(n1,n2,1) &
                                        +tboco*bnd_vars_h%var_2d(nv)%south(n1,n2,2)
              enddo
              enddo
              enddo
            enddo
          endif

          if(nvars_bc_3d_h>0)then
            allocate(targ_3d_h_s(ims:ime,1:lnsh,1:lm,1,1:nvars_bc_3d_h))
            do nv=1,nvars_bc_3d_h
              do n4=1,1
              do n3=1,lm
              do n2=1,lnsh
              do n1=i_lo,i_hi
                targ_3d_h_s(n1,n2,n3,n4,nv)=bnd_vars_h%var_3d(nv)%south(n1,n2,n3,1) &
                                           +tboco*bnd_vars_h%var_3d(nv)%south(n1,n2,n3,2)
              enddo
              enddo
              enddo
              enddo
            enddo
          endif

          if(nvars_bc_4d_h>0)then
            allocate(targ_4d_h_s(1:nvars_bc_4d_h))
            do nv=1,nvars_bc_4d_h
              lbnd=lbound(bnd_vars_h%var_4d(nv)%full_var,4)
              ubnd=ubound(bnd_vars_h%var_4d(nv)%full_var,4)
              allocate(targ_4d_h_s(nv)%base(ims:ime,1:lnsh,1:lm,1,lbnd:ubnd))
              do n5=lbnd,ubnd
              do n4=1,1
              do n3=1,lm
              do n2=1,lnsh
              do n1=i_lo,i_hi
                targ_4d_h_s(nv)%base(n1,n2,n3,n4,n5)=bnd_vars_h%var_4d(nv)%south(n1,n2,n3,1,n5) &
                                                    +tboco*bnd_vars_h%var_4d(nv)%south(n1,n2,n3,2,n5)
              enddo
              enddo
              enddo
              enddo
              enddo
            enddo
          endif

          if(nvars_bc_2d_v>0)then
            allocate(targ_2d_v_s(ims:ime,1:lnsv,1,nvars_bc_2d_v))
            do nv=1,nvars_bc_2d_v
              do n3=1,1
              do n2=1,lnsv
              do n1=i_lo,i_hi
                targ_2d_v_s(n1,n2,n3,nv)=bnd_vars_v%var_2d(nv)%south(n1,n2,1) &
                                        +tboco*bnd_vars_v%var_2d(nv)%south(n1,n2,2)
              enddo
              enddo
              enddo
            enddo
          endif

          if(nvars_bc_3d_v>0)then
            allocate(targ_3d_v_s(ims:ime,1:lnsv,1:lm,1,nvars_bc_3d_v))
            do nv=1,nvars_bc_3d_v
              do n4=1,1
              do n3=1,lm
              do n2=1,lnsv
              do n1=i_lo,i_hi
                targ_3d_v_s(n1,n2,n3,n4,nv)=bnd_vars_v%var_3d(nv)%south(n1,n2,n3,1) &
                                           +tboco*bnd_vars_v%var_3d(nv)%south(n1,n2,n3,2)
              enddo
              enddo
              enddo
              enddo
            enddo
          endif

	ENDIF ! S_BDY

!-----------------------------------------------------------------------
!***  Is this task subdomain on the full domain's west side?
!-----------------------------------------------------------------------

        IF (w_bdy) THEN

          if(nvars_bc_2d_h>0)then
            allocate(targ_2d_h_w(1:lnsh,jms:jme,1,1:nvars_bc_2d_h))
            do nv=1,nvars_bc_2d_h
              do n3=1,1
              do n2=j_lo,j_hi
              do n1=1,lnsh
                targ_2d_h_w(n1,n2,n3,nv)=bnd_vars_h%var_2d(nv)%west(n1,n2,1) &
                                        +tboco*bnd_vars_h%var_2d(nv)%west(n1,n2,2)
              enddo
              enddo
              enddo
            enddo
          endif

          if(nvars_bc_3d_h>0)then
            allocate(targ_3d_h_w(1:lnsh,jms:jme,1:lm,1,1:nvars_bc_3d_h))
            do nv=1,nvars_bc_3d_h
              do n4=1,1
              do n3=1,lm
              do n2=j_lo,j_hi
              do n1=1,lnsh
                targ_3d_h_w(n1,n2,n3,n4,nv)=bnd_vars_h%var_3d(nv)%west(n1,n2,n3,1) &
                                           +tboco*bnd_vars_h%var_3d(nv)%west(n1,n2,n3,2)
              enddo
              enddo
              enddo
              enddo
            enddo
          endif

          if(nvars_bc_4d_h>0)then
            allocate(targ_4d_h_w(1:nvars_bc_4d_h))
            do nv=1,nvars_bc_4d_h
              lbnd=lbound(bnd_vars_h%var_4d(nv)%full_var,4)
              ubnd=ubound(bnd_vars_h%var_4d(nv)%full_var,4)
              allocate(targ_4d_h_w(nv)%base(1:lnsh,jms:jme,1:lm,1,lbnd:ubnd))
              do n5=lbnd,ubnd
              do n4=1,1
              do n3=1,lm
              do n2=j_lo,j_hi
              do n1=1,lnsh
                targ_4d_h_w(nv)%base(n1,n2,n3,n4,n5)=bnd_vars_h%var_4d(nv)%west(n1,n2,n3,1,n5) &
                                                    +tboco*bnd_vars_h%var_4d(nv)%west(n1,n2,n3,2,n5)
              enddo
              enddo
              enddo
              enddo
              enddo
            enddo
          endif

          if(nvars_bc_2d_v>0)then
            allocate(targ_2d_v_w(1:lnsv,jms:jme,1,nvars_bc_2d_v))
            do nv=1,nvars_bc_2d_v
              do n3=1,1
              do n2=j_lo,j_hi
              do n1=1,lnsv
                targ_2d_v_w(n1,n2,n3,nv)=bnd_vars_v%var_2d(nv)%west(n1,n2,1) &
                                        +tboco*bnd_vars_v%var_2d(nv)%west(n1,n2,2)
              enddo
              enddo
              enddo
            enddo
          endif

          if(nvars_bc_3d_v>0)then
            allocate(targ_3d_v_w(1:lnsv,jms:jme,1:lm,1,nvars_bc_3d_v))
            do nv=1,nvars_bc_3d_v
              do n4=1,1
              do n3=1,lm
              do n2=j_lo,j_hi
              do n1=1,lnsv
                targ_3d_v_w(n1,n2,n3,n4,nv)=bnd_vars_v%var_3d(nv)%west(n1,n2,n3,1) &
                                           +tboco*bnd_vars_v%var_3d(nv)%west(n1,n2,n3,2)
              enddo
              enddo
              enddo
              enddo
            enddo
          endif

	ENDIF ! W_BDY

!-----------------------------------------------------------------------
!***  Is this task subdomain on the full domain's east side?
!-----------------------------------------------------------------------

        IF (e_bdy) THEN

          if(nvars_bc_2d_h>0)then
            allocate(targ_2d_h_e(1:lnsh,jms:jme,1,1:nvars_bc_2d_h))
            do nv=1,nvars_bc_2d_h
              do n3=1,1
              do n2=j_lo,j_hi
              do n1=1,lnsh
                targ_2d_h_e(n1,n2,n3,nv)=bnd_vars_h%var_2d(nv)%east(n1,n2,1) &
                                        +tboco*bnd_vars_h%var_2d(nv)%east(n1,n2,2)
              enddo
              enddo
              enddo
            enddo
          endif

          if(nvars_bc_3d_h>0)then
            allocate(targ_3d_h_e(1:lnsh,jms:jme,1:lm,1,1:nvars_bc_3d_h))
            do nv=1,nvars_bc_3d_h
              do n4=1,1
              do n3=1,lm
              do n2=j_lo,j_hi
              do n1=1,lnsh
                targ_3d_h_e(n1,n2,n3,n4,nv)=bnd_vars_h%var_3d(nv)%east(n1,n2,n3,1) &
                                           +tboco*bnd_vars_h%var_3d(nv)%east(n1,n2,n3,2)
              enddo
              enddo
              enddo
              enddo
            enddo
          endif

          if(nvars_bc_4d_h>0)then
            allocate(targ_4d_h_e(1:nvars_bc_4d_h))
            do nv=1,nvars_bc_4d_h
              lbnd=lbound(bnd_vars_h%var_4d(nv)%full_var,4)
              ubnd=ubound(bnd_vars_h%var_4d(nv)%full_var,4)
              allocate(targ_4d_h_e(nv)%base(1:lnsh,jms:jme,1:lm,1,lbnd:ubnd))
              do n5=lbnd,ubnd
              do n4=1,1
              do n3=1,lm
              do n2=j_lo,j_hi
              do n1=1,lnsh
                targ_4d_h_e(nv)%base(n1,n2,n3,n4,n5)=bnd_vars_h%var_4d(nv)%east(n1,n2,n3,1,n5) &
                                                    +tboco*bnd_vars_h%var_4d(nv)%east(n1,n2,n3,2,n5)
              enddo
              enddo
              enddo
              enddo
              enddo
            enddo
          endif

          if(nvars_bc_2d_v>0)then
            allocate(targ_2d_v_e(1:lnsv,jms:jme,1,nvars_bc_2d_v))
            do nv=1,nvars_bc_2d_v
              do n3=1,1
              do n2=j_lo,j_hi
              do n1=1,lnsv
                targ_2d_v_e(n1,n2,n3,nv)=bnd_vars_v%var_2d(nv)%east(n1,n2,1) &
                                        +tboco*bnd_vars_v%var_2d(nv)%east(n1,n2,2)
              enddo
              enddo
              enddo
            enddo
          endif

          if(nvars_bc_3d_v>0)then
            allocate(targ_3d_v_e(1:lnsv,jms:jme,1:lm,1,nvars_bc_3d_v))
            do nv=1,nvars_bc_3d_v
              do n4=1,1
              do n3=1,lm
              do n2=j_lo,j_hi
              do n1=1,lnsv
                targ_3d_v_e(n1,n2,n3,n4,nv)=bnd_vars_v%var_3d(nv)%east(n1,n2,n3,1) &
                                           +tboco*bnd_vars_v%var_3d(nv)%east(n1,n2,n3,2)
              enddo
              enddo
              enddo
              enddo
            enddo
          endif

	ENDIF ! E_BDY

      endif  ! recomp_tend

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Define the first portion of the boundary arrays from the full domain fields
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IF (n_bdy) THEN

        if(nvars_bc_2d_h>0)then
          do nv=1,nvars_bc_2d_h
            do n3=1,1
            do n2=1,lnsh
            do n1=i_lo,i_hi
              bnd_vars_h%var_2d(nv)%north(n1,n2,n3)=               &
                bnd_vars_h%var_2d(nv)%full_var(n1,n2+j_hi-lnsh)
            enddo
            enddo
            enddo
          enddo
        endif

        if(nvars_bc_3d_h>0)then
          do nv=1,nvars_bc_3d_h
            do n4=1,1
            do n3=1,lm
            do n2=1,lnsh
            do n1=i_lo,i_hi
              bnd_vars_h%var_3d(nv)%north(n1,n2,n3,n4)=            &
                bnd_vars_h%var_3d(nv)%full_var(n1,n2+j_hi-lnsh,n3)
            enddo
            enddo
            enddo
            enddo
          enddo
        endif

        if(nvars_bc_4d_h>0)then
          do nv=1,nvars_bc_4d_h
            lbnd=lbound(bnd_vars_h%var_4d(nv)%full_var,4)
            ubnd=ubound(bnd_vars_h%var_4d(nv)%full_var,4)
            do n5=lbnd,ubnd
            do n4=1,1
            do n3=1,lm
            do n2=1,lnsh
            do n1=i_lo,i_hi
              bnd_vars_h%var_4d(nv)%north(n1,n2,n3,n4,n5)=          &
                bnd_vars_h%var_4d(nv)%full_var(n1,n2+j_hi-lnsh,n3,n5)
            enddo
            enddo
            enddo
            enddo
            enddo
          enddo
        endif

        if(nvars_bc_2d_v>0)then
          do nv=1,nvars_bc_2d_v
            do n3=1,1
            do n2=1,lnsv
            do n1=i_lo,i_hi
              bnd_vars_v%var_2d(nv)%north(n1,n2,n3)=               &
                bnd_vars_v%var_2d(nv)%full_var(n1,n2+j_hi-lnsv-1)
            enddo
            enddo
            enddo
          enddo
        endif

        if(nvars_bc_3d_v>0)then
          do nv=1,nvars_bc_3d_v
            do n4=1,1
            do n3=1,lm
            do n2=1,lnsv
            do n1=i_lo,i_hi
              bnd_vars_v%var_3d(nv)%north(n1,n2,n3,n4)=            &
                bnd_vars_v%var_3d(nv)%full_var(n1,n2+j_hi-lnsv,n3)
            enddo
            enddo
            enddo
            enddo
          enddo
        endif

      ENDIF ! N_BDY

      IF (s_bdy) THEN

        if(nvars_bc_2d_h>0)then
          do nv=1,nvars_bc_2d_h
            do n3=1,1
            do n2=1,lnsh
            do n1=i_lo,i_hi
              bnd_vars_h%var_2d(nv)%south(n1,n2,n3)=               &
                bnd_vars_h%var_2d(nv)%full_var(n1,j_lo+n2-1)
            enddo
            enddo
            enddo
          enddo
        endif

        if(nvars_bc_3d_h>0)then
          do nv=1,nvars_bc_3d_h
            do n4=1,1
            do n3=1,lm
            do n2=1,lnsh
            do n1=i_lo,i_hi
              bnd_vars_h%var_3d(nv)%south(n1,n2,n3,n4)=            &
                bnd_vars_h%var_3d(nv)%full_var(n1,j_lo+n2-1,n3)
            enddo
            enddo
            enddo
            enddo
          enddo
        endif

        if(nvars_bc_4d_h>0)then
          do nv=1,nvars_bc_4d_h
            lbnd=lbound(bnd_vars_h%var_4d(nv)%full_var,4)
            ubnd=ubound(bnd_vars_h%var_4d(nv)%full_var,4)
            do n5=lbnd,ubnd
            do n4=1,1
            do n3=1,lm
            do n2=1,lnsh
            do n1=i_lo,i_hi
              bnd_vars_h%var_4d(nv)%south(n1,n2,n3,n4,n5)=          &
                bnd_vars_h%var_4d(nv)%full_var(n1,j_lo+n2-1,n3,n5)
            enddo
            enddo
            enddo
            enddo
            enddo
          enddo
        endif

        if(nvars_bc_2d_v>0)then
          do nv=1,nvars_bc_2d_v
            do n3=1,1
            do n2=1,lnsv
            do n1=i_lo,i_hi
              bnd_vars_v%var_2d(nv)%south(n1,n2,n3)=               &
                bnd_vars_v%var_2d(nv)%full_var(n1,j_lo+n2-1)
            enddo
            enddo
            enddo
          enddo
        endif

        if(nvars_bc_3d_v>0)then
          do nv=1,nvars_bc_3d_v
            do n4=1,1
            do n3=1,lm
            do n2=1,lnsv
            do n1=i_lo,i_hi
              bnd_vars_v%var_3d(nv)%south(n1,n2,n3,n4)=            &
                bnd_vars_v%var_3d(nv)%full_var(n1,j_lo+n2-1,n3)
            enddo
            enddo
            enddo
            enddo
          enddo
        endif

      ENDIF ! S_BDY

      IF (w_bdy) THEN

        if(nvars_bc_2d_h>0)then
          do nv=1,nvars_bc_2d_h
            do n3=1,1
            do n2=j_lo,j_hi
            do n1=1,lnsh
              bnd_vars_h%var_2d(nv)%west(n1,n2,n3)=                &
                bnd_vars_h%var_2d(nv)%full_var(i_lo+n1-1,n2)
            enddo
            enddo
            enddo
          enddo
        endif

        if(nvars_bc_3d_h>0)then
          do nv=1,nvars_bc_3d_h
            do n4=1,1
            do n3=1,lm
            do n2=j_lo,j_hi
            do n1=1,lnsh
              bnd_vars_h%var_3d(nv)%west(n1,n2,n3,n4)=             &
                bnd_vars_h%var_3d(nv)%full_var(i_lo+n1-1,n2,n3)
            enddo
            enddo
            enddo
            enddo
          enddo
        endif

        if(nvars_bc_4d_h>0)then
          do nv=1,nvars_bc_4d_h
            lbnd=lbound(bnd_vars_h%var_4d(nv)%full_var,4)
            ubnd=ubound(bnd_vars_h%var_4d(nv)%full_var,4)
            do n5=lbnd,ubnd
            do n4=1,1
            do n3=1,lm
            do n2=j_lo,j_hi
            do n1=1,lnsh
              bnd_vars_h%var_4d(nv)%west(n1,n2,n3,n4,n5)=           &
                bnd_vars_h%var_4d(nv)%full_var(i_lo+n1-1,n2,n3,n5)
            enddo
            enddo
            enddo
            enddo
            enddo
          enddo
        endif

        if(nvars_bc_2d_v>0)then
          do nv=1,nvars_bc_2d_v
            do n3=1,1
            do n2=j_lo,j_hi
            do n1=1,lnsv
              bnd_vars_v%var_2d(nv)%west(n1,n2,n3)=                &
                bnd_vars_v%var_2d(nv)%full_var(i_lo+n1-1,n2)
            enddo
            enddo
            enddo
          enddo
        endif

        if(nvars_bc_3d_v>0)then
          do nv=1,nvars_bc_3d_v
            do n4=1,1
            do n3=1,lm
            do n2=j_lo,j_hi
            do n1=1,lnsv
              bnd_vars_v%var_3d(nv)%west(n1,n2,n3,n4)=             &
                bnd_vars_v%var_3d(nv)%full_var(i_lo+n1-1,n2,n3)
            enddo
            enddo
            enddo
            enddo
          enddo
        endif

      ENDIF ! W_BDY

      IF (e_bdy) THEN

        if(nvars_bc_2d_h>0)then
          do nv=1,nvars_bc_2d_h
            do n3=1,1
            do n2=j_lo,j_hi
            do n1=1,lnsh
              bnd_vars_h%var_2d(nv)%east(n1,n2,n3)=                &
                bnd_vars_h%var_2d(nv)%full_var(n1+i_hi-lnsh,n2)
            enddo
            enddo
            enddo
          enddo
        endif

        if(nvars_bc_3d_h>0)then
          do nv=1,nvars_bc_3d_h
            do n4=1,1
            do n3=1,lm
            do n2=j_lo,j_hi
            do n1=1,lnsh
              bnd_vars_h%var_3d(nv)%east(n1,n2,n3,n4)=             &
                bnd_vars_h%var_3d(nv)%full_var(n1+i_hi-lnsh,n2,n3)
            enddo
            enddo
            enddo
            enddo
          enddo
        endif

        if(nvars_bc_4d_h>0)then
          do nv=1,nvars_bc_4d_h
            lbnd=lbound(bnd_vars_h%var_4d(nv)%full_var,4)
            ubnd=ubound(bnd_vars_h%var_4d(nv)%full_var,4)
            do n5=lbnd,ubnd
            do n4=1,1
            do n3=1,lm
            do n2=j_lo,j_hi
            do n1=1,lnsh
              bnd_vars_h%var_4d(nv)%east(n1,n2,n3,n4,n5)=           &
                bnd_vars_h%var_4d(nv)%full_var(n1+i_hi-lnsh,n2,n3,n5)
            enddo
            enddo
            enddo
            enddo
            enddo
          enddo
        endif

        if(nvars_bc_2d_v>0)then
          do nv=1,nvars_bc_2d_v
            do n3=1,1
            do n2=j_lo,j_hi
            do n1=1,lnsv
              bnd_vars_v%var_2d(nv)%east(n1,n2,n3)=                &
                bnd_vars_v%var_2d(nv)%full_var(n1+i_hi-lnsh,n2)
            enddo
            enddo
            enddo
          enddo
        endif

        if(nvars_bc_3d_v>0)then
          do nv=1,nvars_bc_3d_v
            do n4=1,1
            do n3=1,lm
            do n2=j_lo,j_hi
            do n1=1,lnsv
              bnd_vars_v%var_3d(nv)%east(n1,n2,n3,n4)=             &
                bnd_vars_v%var_3d(nv)%full_var(n1+i_hi-lnsv-1,n2,n3)
            enddo
            enddo
            enddo
            enddo
          enddo
        endif

      ENDIF ! E_BDY

      if (recomp_tend) then
    
        IF (n_bdy) THEN

          if(nvars_bc_2d_h>0)then
            do nv=1,nvars_bc_2d_h
              do n2=1,lnsh
              do n1=i_lo,i_hi
                bnd_vars_h%var_2d(nv)%north(n1,n2,2)=   & 
                  (targ_2d_h_n(n1,n2,1,nv)-bnd_vars_h%var_2d(nv)%north(n1,n2,1))/tboco
!                 (targ_2d_h_n(n1,n2,1,nv)-bnd_vars_h%var_2d(nv)%north(n1,n2,1))*r_tboco
              enddo
              enddo
            enddo
            deallocate(targ_2d_h_n)
          endif

          if(nvars_bc_3d_h>0)then
            do nv=1,nvars_bc_3d_h
              do n3=1,lm
              do n2=1,lnsh
              do n1=i_lo,i_hi
                bnd_vars_h%var_3d(nv)%north(n1,n2,n3,2)=   & 
                  (targ_3d_h_n(n1,n2,n3,1,nv)-bnd_vars_h%var_3d(nv)%north(n1,n2,n3,1))/tboco
!                 (targ_3d_h_n(n1,n2,n3,1,nv)-bnd_vars_h%var_3d(nv)%north(n1,n2,n3,1))*r_tboco
              enddo
              enddo
              enddo
            enddo
            deallocate(targ_3d_h_n)
          endif

          if(nvars_bc_4d_h>0)then
            do nv=1,nvars_bc_4d_h
              lbnd=lbound(bnd_vars_h%var_4d(nv)%full_var,4)
              ubnd=ubound(bnd_vars_h%var_4d(nv)%full_var,4)
              do n4=lbnd,ubnd
              do n3=1,lm
              do n2=1,lnsh
              do n1=i_lo,i_hi
                bnd_vars_h%var_4d(nv)%north(n1,n2,n3,2,n4)=   & 
                  (targ_4d_h_n(nv)%base(n1,n2,n3,1,nv)-bnd_vars_h%var_4d(nv)%north(n1,n2,n3,1,n4))/tboco
!                 (targ_4d_h_n(nv)%base(n1,n2,n3,1,nv)-bnd_vars_h%var_4d(nv)%north(n1,n2,n3,1,n4))*r_tboco
              enddo
              enddo
              enddo
              enddo
              deallocate(targ_4d_h_n(nv)%base)
            enddo
            deallocate(targ_4d_h_n)
          endif

          if(nvars_bc_2d_v>0)then
            do nv=1,nvars_bc_2d_v
              do n2=1,lnsv
              do n1=i_lo,i_hi
                bnd_vars_v%var_2d(nv)%north(n1,n2,2)=   & 
                  (targ_2d_v_n(n1,n2,1,nv)-bnd_vars_v%var_2d(nv)%north(n1,n2,1))/tboco
!                 (targ_2d_v_n(n1,n2,1,nv)-bnd_vars_v%var_2d(nv)%north(n1,n2,1))*r_tboco
              enddo
              enddo
            enddo
            deallocate(targ_2d_v_n)
          endif

          if(nvars_bc_3d_v>0)then
            do nv=1,nvars_bc_3d_v
              do n3=1,lm
              do n2=1,lnsv
              do n1=i_lo,i_hi
                bnd_vars_v%var_3d(nv)%north(n1,n2,n3,2)=   & 
                  (targ_3d_v_n(n1,n2,n3,1,nv)-bnd_vars_v%var_3d(nv)%north(n1,n2,n3,1))/tboco
!                 (targ_3d_v_n(n1,n2,n3,1,nv)-bnd_vars_v%var_3d(nv)%north(n1,n2,n3,1))*r_tboco
              enddo
              enddo
              enddo
            enddo
            deallocate(targ_3d_v_n)
          endif

	ENDIF ! N_BDY

        IF (s_bdy) THEN

          if(nvars_bc_2d_h>0)then
            do nv=1,nvars_bc_2d_h
              do n2=1,lnsh
              do n1=i_lo,i_hi
                bnd_vars_h%var_2d(nv)%south(n1,n2,2)=   &
                  (targ_2d_h_s(n1,n2,1,nv)-bnd_vars_h%var_2d(nv)%south(n1,n2,1))/tboco
!                 (targ_2d_h_s(n1,n2,1,nv)-bnd_vars_h%var_2d(nv)%south(n1,n2,1))*r_tboco
              enddo
              enddo
            enddo
            deallocate(targ_2d_h_s)
          endif

          if(nvars_bc_3d_h>0)then
            do nv=1,nvars_bc_3d_h
              do n3=1,lm
              do n2=1,lnsh
              do n1=i_lo,i_hi
                bnd_vars_h%var_3d(nv)%south(n1,n2,n3,2)=   &
                  (targ_3d_h_s(n1,n2,n3,1,nv)-bnd_vars_h%var_3d(nv)%south(n1,n2,n3,1))/tboco
!                 (targ_3d_h_s(n1,n2,n3,1,nv)-bnd_vars_h%var_3d(nv)%south(n1,n2,n3,1))*r_tboco
              enddo
              enddo
              enddo
            enddo
            deallocate(targ_3d_h_s)
          endif

          if(nvars_bc_4d_h>0)then
            do nv=1,nvars_bc_4d_h
              lbnd=lbound(bnd_vars_h%var_4d(nv)%full_var,4)
              ubnd=ubound(bnd_vars_h%var_4d(nv)%full_var,4)
              do n4=lbnd,ubnd
              do n3=1,lm
              do n2=1,lnsh
              do n1=i_lo,i_hi
                bnd_vars_h%var_4d(nv)%south(n1,n2,n3,2,n4)=   &
                  (targ_4d_h_s(nv)%base(n1,n2,n3,1,nv)-bnd_vars_h%var_4d(nv)%south(n1,n2,n3,1,n4))/tboco
!                 (targ_4d_h_s(nv)%base(n1,n2,n3,1,nv)-bnd_vars_h%var_4d(nv)%south(n1,n2,n3,1,n4))*r_tboco
              enddo
              enddo
              enddo
              enddo
              deallocate(targ_4d_h_s(nv)%base)
            enddo
            deallocate(targ_4d_h_s)
          endif

          if(nvars_bc_2d_v>0)then
            do nv=1,nvars_bc_2d_v
              do n2=1,lnsv
              do n1=i_lo,i_hi
                bnd_vars_v%var_2d(nv)%south(n1,n2,2)=   &
                  (targ_2d_v_s(n1,n2,1,nv)-bnd_vars_v%var_2d(nv)%south(n1,n2,1))/tboco
!                 (targ_2d_v_s(n1,n2,1,nv)-bnd_vars_v%var_2d(nv)%south(n1,n2,1))*r_tboco
              enddo
              enddo
            enddo
            deallocate(targ_2d_v_s)
          endif

          if(nvars_bc_3d_v>0)then
            do nv=1,nvars_bc_3d_v
              do n3=1,lm
              do n2=1,lnsv
              do n1=i_lo,i_hi
                bnd_vars_v%var_3d(nv)%south(n1,n2,n3,2)=   &
                  (targ_3d_v_s(n1,n2,n3,1,nv)-bnd_vars_v%var_3d(nv)%south(n1,n2,n3,1))/tboco
!                 (targ_3d_v_s(n1,n2,n3,1,nv)-bnd_vars_v%var_3d(nv)%south(n1,n2,n3,1))*r_tboco
              enddo
              enddo
              enddo
            enddo
            deallocate(targ_3d_v_s)
          endif

	ENDIF ! S_BDY

        IF (w_bdy) THEN

          if(nvars_bc_2d_h>0)then
            do nv=1,nvars_bc_2d_h
              do n2=j_lo,j_hi
              do n1=1,lnsh
                bnd_vars_h%var_2d(nv)%west(n1,n2,2)=   &
                  (targ_2d_h_w(n1,n2,1,nv)-bnd_vars_h%var_2d(nv)%west(n1,n2,1))/tboco
!                 (targ_2d_h_w(n1,n2,1,nv)-bnd_vars_h%var_2d(nv)%west(n1,n2,1))*r_tboco
              enddo
              enddo
            enddo
            deallocate(targ_2d_h_w)
          endif

          if(nvars_bc_3d_h>0)then
            do nv=1,nvars_bc_3d_h
              do n3=1,lm
              do n2=j_lo,j_hi
              do n1=1,lnsh
                bnd_vars_h%var_3d(nv)%west(n1,n2,n3,2)=   &
                  (targ_3d_h_w(n1,n2,n3,1,nv)-bnd_vars_h%var_3d(nv)%west(n1,n2,n3,1))/tboco
!                 (targ_3d_h_w(n1,n2,n3,1,nv)-bnd_vars_h%var_3d(nv)%west(n1,n2,n3,1))*r_tboco
              enddo
              enddo
              enddo
            enddo
            deallocate(targ_3d_h_w)
          endif

          if(nvars_bc_4d_h>0)then
            do nv=1,nvars_bc_4d_h
              lbnd=lbound(bnd_vars_h%var_4d(nv)%full_var,4)
              ubnd=ubound(bnd_vars_h%var_4d(nv)%full_var,4)
              do n4=lbnd,ubnd
              do n3=1,lm
              do n2=j_lo,j_hi
              do n1=1,lnsh
                bnd_vars_h%var_4d(nv)%west(n1,n2,n3,2,n4)=   &
                  (targ_4d_h_w(nv)%base(n1,n2,n3,1,nv)-bnd_vars_h%var_4d(nv)%west(n1,n2,n3,1,n4))/tboco
!                 (targ_4d_h_w(nv)%base(n1,n2,n3,1,nv)-bnd_vars_h%var_4d(nv)%west(n1,n2,n3,1,n4))*r_tboco
              enddo
              enddo
              enddo
              enddo
              deallocate(targ_4d_h_w(nv)%base)
            enddo
            deallocate(targ_4d_h_w)
          endif

          if(nvars_bc_2d_v>0)then
            do nv=1,nvars_bc_2d_v
              do n2=j_lo,j_hi
              do n1=1,lnsv
                bnd_vars_v%var_2d(nv)%west(n1,n2,2)=   &
                  (targ_2d_v_w(n1,n2,1,nv)-bnd_vars_v%var_2d(nv)%west(n1,n2,1))/tboco
!                 (targ_2d_v_w(n1,n2,1,nv)-bnd_vars_v%var_2d(nv)%west(n1,n2,1))*r_tboco
              enddo
              enddo
            enddo
            deallocate(targ_2d_v_w)
          endif

          if(nvars_bc_3d_v>0)then
            do nv=1,nvars_bc_3d_v
              do n3=1,lm
              do n2=j_lo,j_hi
              do n1=1,lnsv
                bnd_vars_v%var_3d(nv)%west(n1,n2,n3,2)=   &
                  (targ_3d_v_w(n1,n2,n3,1,nv)-bnd_vars_v%var_3d(nv)%west(n1,n2,n3,1))/tboco
!                 (targ_3d_v_w(n1,n2,n3,1,nv)-bnd_vars_v%var_3d(nv)%west(n1,n2,n3,1))*r_tboco
              enddo
              enddo
              enddo
            enddo
            deallocate(targ_3d_v_w)
          endif

	ENDIF ! W_BDY

        IF (e_bdy) THEN

          if(nvars_bc_2d_h>0)then
            do nv=1,nvars_bc_2d_h
              do n2=j_lo,j_hi
              do n1=1,lnsh
                bnd_vars_h%var_2d(nv)%east(n1,n2,2)=   &
                  (targ_2d_h_e(n1,n2,1,nv)-bnd_vars_h%var_2d(nv)%east(n1,n2,1))/tboco
!                 (targ_2d_h_e(n1,n2,1,nv)-bnd_vars_h%var_2d(nv)%east(n1,n2,1))*r_tboco
              enddo
              enddo
            enddo
            deallocate(targ_2d_h_e)
          endif

          if(nvars_bc_3d_h>0)then
            do nv=1,nvars_bc_3d_h
              do n3=1,lm
              do n2=j_lo,j_hi
              do n1=1,lnsh
                bnd_vars_h%var_3d(nv)%east(n1,n2,n3,2)=   &
                  (targ_3d_h_e(n1,n2,n3,1,nv)-bnd_vars_h%var_3d(nv)%east(n1,n2,n3,1))/tboco
!                 (targ_3d_h_e(n1,n2,n3,1,nv)-bnd_vars_h%var_3d(nv)%east(n1,n2,n3,1))*r_tboco
              enddo
              enddo
              enddo
            enddo
            deallocate(targ_3d_h_e)
          endif

          if(nvars_bc_4d_h>0)then
            do nv=1,nvars_bc_4d_h
              lbnd=lbound(bnd_vars_h%var_4d(nv)%full_var,4)
              ubnd=ubound(bnd_vars_h%var_4d(nv)%full_var,4)
              do n4=lbnd,ubnd
              do n3=1,lm
              do n2=j_lo,j_hi
              do n1=1,lnsh
                bnd_vars_h%var_4d(nv)%east(n1,n2,n3,2,n4)=   &
                  (targ_4d_h_e(nv)%base(n1,n2,n3,1,nv)-bnd_vars_h%var_4d(nv)%east(n1,n2,n3,1,n4))/tboco
!                 (targ_4d_h_e(nv)%base(n1,n2,n3,1,nv)-bnd_vars_h%var_4d(nv)%east(n1,n2,n3,1,n4))*r_tboco
              enddo
              enddo
              enddo
              enddo
              deallocate(targ_4d_h_e(nv)%base)
            enddo
            deallocate(targ_4d_h_e)
          endif

          if(nvars_bc_2d_v>0)then
            do nv=1,nvars_bc_2d_v
              do n2=j_lo,j_hi
              do n1=1,lnsv
                bnd_vars_v%var_2d(nv)%east(n1,n2,2)=   &
                  (targ_2d_v_e(n1,n2,1,nv)-bnd_vars_v%var_2d(nv)%east(n1,n2,1))/tboco
!                 (targ_2d_v_e(n1,n2,1,nv)-bnd_vars_v%var_2d(nv)%east(n1,n2,1))*r_tboco
              enddo
              enddo
            enddo
            deallocate(targ_2d_v_e)
          endif

          if(nvars_bc_3d_v>0)then
            do nv=1,nvars_bc_3d_v
              do n3=1,lm
              do n2=j_lo,j_hi
              do n1=1,lnsv
                bnd_vars_v%var_3d(nv)%east(n1,n2,n3,2)=   &
                  (targ_3d_v_e(n1,n2,n3,1,nv)-bnd_vars_v%var_3d(nv)%east(n1,n2,n3,1))/tboco
!                 (targ_3d_v_e(n1,n2,n3,1,nv)-bnd_vars_v%var_3d(nv)%east(n1,n2,n3,1))*r_tboco
              enddo
              enddo
              enddo
            enddo
            deallocate(targ_3d_v_e)
          endif

	ENDIF ! E_BDY

      endif ! recomp_tend
!
!-----------------------------------------------------------------------
                        end subroutine write_bc
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
                        subroutine bocoh &
(lm,lnsh &
,dt,pt &
,pd,dsg2,pdsg1 &
,nvars_bc_2d_h,nvars_bc_3d_h,nvars_bc_4d_h &
,lbnd_4d,ubnd_4d &
,bnd_vars_h &
,pint)
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
implicit none
!-----------------------------------------------------------------------

real(kind=kfpt),parameter:: &
 wa=0.5 &                    ! weighting factor
,w1=wa*0.25 &                ! weighting factor
,w2=1.-wa                    ! weighting factor
!-----------------------------------------------------------------------
integer(kind=kint),intent(in):: &
 lm &                        ! total # of levels
,lnsh                        ! blending area width, h points

integer(kind=kint),intent(in) :: &
 nvars_bc_2d_h &             ! # of 2-d h-pt boundary variables
,nvars_bc_3d_h &             ! # of 3-d h-pt boundary variables
,nvars_bc_4d_h               ! # of 4-d h-pt boundary variables

integer(kind=kint),dimension(*),intent(in) :: &
 lbnd_4d &                   ! lower index of list of 4-D h-pt boundary variables
,ubnd_4d                     ! upper index of list of 4-D h-pt boundary variables

real(kind=kfpt),intent(in):: &
 dt &                        ! dynamics time step
,pt                          ! pressure at the top

real(kind=kfpt),dimension(1:lm),intent(in):: &
 dsg2 &                      ! delta sigmas
,pdsg1                       ! delta pressures

real(kind=kfpt),dimension(ims:ime,jms:jme),intent(in):: &
 pd

type(bc_h_all),intent(inout) :: &
 bnd_vars_h                  ! boundary variables on h points

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm+1),intent(inout):: &
 pint                        ! pressure at interfaces

!-----------------------------------------------------------------------
!---local variables-----------------------------------------------------
!-----------------------------------------------------------------------

integer(kind=kint):: &
 i &                         ! index in x direction
,i1 &
,i2 &
,ib &                        ! index in x direction, boundary zone
,ihe &                       ! ending index in x direction, boundaries
,ihs &                       ! starting index in x direction, boundaries
,j &                         ! index in y direction
,j1 &
,j2 &
,jb &                        ! index in y direction, boundary zone
,jhe &                       ! ending index in x direction, boundaries
,jhs &                       ! starting index in x direction, boundaries
,k &                         ! boundary line counter
,ks &                        ! smoothing counter
,l &                         ! index in p direction
,lbnd &
,lines &                     ! boundary smoothing area
,nl &
,nsmud &                     ! number of smoothing passes
,nv &
,ubnd

real(kind=kfpt),dimension(1:lnsh):: &
 wh(lnsh) &                  ! blending weighting function, temperature
,wq(lnsh)                    ! blending weighting function, moisture

real(kind=kfpt),dimension(:,:),allocatable :: &
 hbc_2d

real(kind=kfpt),dimension(:,:,:),pointer :: temp_3d

!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
      lines=lnsh
!     nsmud=lines-1
      nsmud=0
!-----------------------------------------------------------------------
      wh(1)=1.
      wq(1)=1.
      do k=2,lnsh
         wh(k)=1.-(0.9/real(lnsh-1))*(k-1)
         wq(k)=1.-(0.9/real(lnsh-1))*(k-1)
      enddo
!-----------------------------------------------------------------------
!***  Update values of the boundary working objects at H points.
!-----------------------------------------------------------------------
!
!--------------------------------------
!***  Southern and northern boundaries
!--------------------------------------
!
      if(s_bdy)then
!
        if(nvars_bc_2d_h>0)then
          do nv=1,nvars_bc_2d_h
            do jb=1,lnsh
            do ib=its_h2,ite_h2
              bnd_vars_h%var_2d(nv)%south(ib,jb,1)=bnd_vars_h%var_2d(nv)%south(ib,jb,1)   &
                                                  +bnd_vars_h%var_2d(nv)%south(ib,jb,2)*dt
            enddo
            enddo
          enddo
        endif
!
        if(nvars_bc_3d_h>0)then
          do nv=1,nvars_bc_3d_h
            do l=1,lm
            do jb=1,lnsh
            do ib=its_h2,ite_h2
              bnd_vars_h%var_3d(nv)%south(ib,jb,l,1)=bnd_vars_h%var_3d(nv)%south(ib,jb,l,1)   &
                                                    +bnd_vars_h%var_3d(nv)%south(ib,jb,l,2)*dt
            enddo
            enddo
            enddo
          enddo
        endif
!
        if(nvars_bc_4d_h>0)then
          do nv=1,nvars_bc_4d_h
            lbnd=lbnd_4d(nv)
            ubnd=ubnd_4d(nv)
            do nl=lbnd,ubnd
              do l=1,lm
              do jb=1,lnsh
              do ib=its_h2,ite_h2
                bnd_vars_h%var_4d(nv)%south(ib,jb,l,1,nl)=bnd_vars_h%var_4d(nv)%south(ib,jb,l,1,nl)   &
                                                         +bnd_vars_h%var_4d(nv)%south(ib,jb,l,2,nl)*dt
              enddo
              enddo
              enddo
            enddo
          enddo
        endif
!
      endif
!
      if(n_bdy)then
!
        if(nvars_bc_2d_h>0)then
          do nv=1,nvars_bc_2d_h
            do jb=1,lnsh
            do ib=its_h2,ite_h2
              bnd_vars_h%var_2d(nv)%north(ib,jb,1)=bnd_vars_h%var_2d(nv)%north(ib,jb,1)   &
                                                  +bnd_vars_h%var_2d(nv)%north(ib,jb,2)*dt
            enddo
            enddo
          enddo
        endif
!
        if(nvars_bc_3d_h>0)then
          do nv=1,nvars_bc_3d_h
            do l=1,lm
            do jb=1,lnsh
            do ib=its_h2,ite_h2
              bnd_vars_h%var_3d(nv)%north(ib,jb,l,1)=bnd_vars_h%var_3d(nv)%north(ib,jb,l,1)   &
                                                    +bnd_vars_h%var_3d(nv)%north(ib,jb,l,2)*dt
            enddo
            enddo
            enddo
          enddo
        endif
!
        if(nvars_bc_4d_h>0)then
          do nv=1,nvars_bc_4d_h
            lbnd=lbnd_4d(nv)
            ubnd=ubnd_4d(nv)
            do nl=lbnd,ubnd
              do l=1,lm
              do jb=1,lnsh
              do ib=its_h2,ite_h2
                bnd_vars_h%var_4d(nv)%north(ib,jb,l,1,nl)=bnd_vars_h%var_4d(nv)%north(ib,jb,l,1,nl)   &
                                                         +bnd_vars_h%var_4d(nv)%north(ib,jb,l,2,nl)*dt
              enddo
              enddo
              enddo
            enddo
          enddo
        endif
!
      endif
!
!------------------------------------
!***  Western and eastern boundaries
!------------------------------------
!
      if(w_bdy)then
!
        if(nvars_bc_2d_h>0)then
          do nv=1,nvars_bc_2d_h
            do jb=jts_h2,jte_h2
            do ib=1,lnsh
              bnd_vars_h%var_2d(nv)%west(ib,jb,1)=bnd_vars_h%var_2d(nv)%west(ib,jb,1)   &
                                                 +bnd_vars_h%var_2d(nv)%west(ib,jb,2)*dt
            enddo
            enddo
          enddo
        endif
!
        if(nvars_bc_3d_h>0)then
          do nv=1,nvars_bc_3d_h
            do l=1,lm
            do jb=jts_h2,jte_h2
            do ib=1,lnsh
              bnd_vars_h%var_3d(nv)%west(ib,jb,l,1)=bnd_vars_h%var_3d(nv)%west(ib,jb,l,1)   &
                                                   +bnd_vars_h%var_3d(nv)%west(ib,jb,l,2)*dt
            enddo
            enddo
            enddo
          enddo
        endif
!
        if(nvars_bc_4d_h>0)then
          do nv=1,nvars_bc_4d_h
            lbnd=lbnd_4d(nv)
            ubnd=ubnd_4d(nv)
            do nl=lbnd,ubnd
              do l=1,lm
              do jb=jts_h2,jte_h2
              do ib=1,lnsh
                bnd_vars_h%var_4d(nv)%west(ib,jb,l,1,nl)=bnd_vars_h%var_4d(nv)%west(ib,jb,l,1,nl)   &
                                                        +bnd_vars_h%var_4d(nv)%west(ib,jb,l,2,nl)*dt
              enddo
              enddo
              enddo
            enddo
          enddo
        endif
!
      endif
!
      if(e_bdy)then
!
        if(nvars_bc_2d_h>0)then
          do nv=1,nvars_bc_2d_h
            do jb=jts_h2,jte_h2
            do ib=1,lnsh
              bnd_vars_h%var_2d(nv)%east(ib,jb,1)=bnd_vars_h%var_2d(nv)%east(ib,jb,1)   &
                                                 +bnd_vars_h%var_2d(nv)%east(ib,jb,2)*dt
            enddo
            enddo
          enddo
        endif
!
        if(nvars_bc_3d_h>0)then
          do nv=1,nvars_bc_3d_h
            do l=1,lm
            do jb=jts_h2,jte_h2
            do ib=1,lnsh
              bnd_vars_h%var_3d(nv)%east(ib,jb,l,1)=bnd_vars_h%var_3d(nv)%east(ib,jb,l,1)   &
                                                   +bnd_vars_h%var_3d(nv)%east(ib,jb,l,2)*dt
            enddo
            enddo
            enddo
          enddo
        endif
!
        if(nvars_bc_4d_h>0)then
          do nv=1,nvars_bc_4d_h
            lbnd=lbnd_4d(nv)
            ubnd=ubnd_4d(nv)
            do nl=lbnd,ubnd
              do l=1,lm
              do jb=jts_h2,jte_h2
              do ib=1,lnsh
                bnd_vars_h%var_4d(nv)%east(ib,jb,l,1,nl)=bnd_vars_h%var_4d(nv)%east(ib,jb,l,1,nl)   &
                                                        +bnd_vars_h%var_4d(nv)%east(ib,jb,l,2,nl)*dt
              enddo
              enddo
              enddo
            enddo
          enddo
        endif
!
      endif
!
!-----------------------------------------------------------------------
!***  Now update the actual prognostic variables using the values 
!***  just computed (coming from the parent) blended with the values
!***  already generated in the blending region by the child.
!-----------------------------------------------------------------------
!
!-----------------------
!***  Southern boundary
!-----------------------
!
      if(s_bdy)then
!
        if(nvars_bc_2d_h>0)then
          do nv=1,nvars_bc_2d_h
            do j=1,lnsh
              jb=j
              ihs=jb
              ihe=ide+1-jb
              do i=max(its_h2,ihs),min(ite_h2,ihe)
                ib=i
                bnd_vars_h%var_2d(nv)%full_var(i,j)=bnd_vars_h%var_2d(nv)%south(ib,jb,1)*wh(jb)    &
                                                   +bnd_vars_h%var_2d(nv)%full_var(i,j)*(1.-wh(jb))
              enddo
            enddo
          enddo
        endif
!
        if(nvars_bc_3d_h>0)then
          do nv=1,nvars_bc_3d_h
            do l=1,lm
              do j=1,lnsh
                jb=j
                ihs=jb
                ihe=ide+1-jb
                do i=max(its_h2,ihs),min(ite_h2,ihe)
                  ib=i
                  bnd_vars_h%var_3d(nv)%full_var(i,j,l)=bnd_vars_h%var_3d(nv)%south(ib,jb,l,1)*wh(jb)    &
                                                       +bnd_vars_h%var_3d(nv)%full_var(i,j,l)*(1.-wh(jb))
                enddo
              enddo
            enddo
          enddo
        endif
!
        if(nvars_bc_4d_h>0)then
          do nv=1,nvars_bc_4d_h
            lbnd=lbnd_4d(nv)
            ubnd=ubnd_4d(nv)
            do nl=lbnd,ubnd
              do l=1,lm
                do j=1,lnsh
                  jb=j
                  ihs=jb
                  ihe=ide+1-jb
                  do i=max(its_h2,ihs),min(ite_h2,ihe)
                    ib=i
                    bnd_vars_h%var_4d(nv)%full_var(i,j,l,nl)=bnd_vars_h%var_4d(nv)%south(ib,jb,l,1,nl)*wh(jb)    &
                                                            +bnd_vars_h%var_4d(nv)%full_var(i,j,l,nl)*(1.-wh(jb))
                  enddo
                enddo
              enddo
            enddo
          enddo
        endif
!
!--------------------
!***  The pint array
!--------------------
!
        do j=1,lnsh
          jb=j
          ihs=jb
          ihe=ide+1-jb
          do i=max(its_h2,ihs),min(ite_h2,ihe)
            ib=i
            pint(i,j,1)=pt
          enddo
        enddo
!
        do l=1,lm
          do j=1,lnsh
            jb=j
            ihs=jb
            ihe=ide+1-jb
            do i=max(its_h2,ihs),min(ite_h2,ihe)
              ib=i
              pint(i,j,l+1)=(pint(i,j,l) &
                            +(dsg2(l)*pd(i,j)+pdsg1(l)))*wh(jb) &
                           +pint(i,j,l+1)*(1.-wh(jb))
            enddo
          enddo
        enddo
!
      endif
!
!-----------------------
!***  Northern boundary
!-----------------------
!
      if(n_bdy)then
!
        if(nvars_bc_2d_h>0)then
          do nv=1,nvars_bc_2d_h
            do j=jde-lnsh+1,jde
              jb=j-jde+lnsh
              ihs=1-jb+lnsh
              ihe=ide+jb-lnsh
              do i=max(its_h2,ihs),min(ite_h2,ihe)
                ib=i
                bnd_vars_h%var_2d(nv)%full_var(i,j)=bnd_vars_h%var_2d(nv)%north(ib,jb,1)*wh(lnsh+1-jb)    &
                                                   +bnd_vars_h%var_2d(nv)%full_var(i,j)*(1.-wh(lnsh+1-jb))
              enddo
            enddo
          enddo
        endif
!
        if(nvars_bc_3d_h>0)then
          do nv=1,nvars_bc_3d_h
            do l=1,lm
              do j=jde-lnsh+1,jde
                jb=j-jde+lnsh
                ihs=1-jb+lnsh
                ihe=ide+jb-lnsh
                do i=max(its_h2,ihs),min(ite_h2,ihe)
                  ib=i
                  bnd_vars_h%var_3d(nv)%full_var(i,j,l)=bnd_vars_h%var_3d(nv)%north(ib,jb,l,1)*wh(lnsh+1-jb)    &
                                                       +bnd_vars_h%var_3d(nv)%full_var(i,j,l)*(1.-wh(lnsh+1-jb))
                enddo
              enddo
            enddo
          enddo
        endif
!
        if(nvars_bc_4d_h>0)then
          do nv=1,nvars_bc_4d_h
            lbnd=lbnd_4d(nv)
            ubnd=ubnd_4d(nv)
            do nl=lbnd,ubnd
              do l=1,lm
                do j=jde-lnsh+1,jde
                  jb=j-jde+lnsh
                  ihs=1-jb+lnsh
                  ihe=ide+jb-lnsh
                  do i=max(its_h2,ihs),min(ite_h2,ihe)
                    ib=i
                    bnd_vars_h%var_4d(nv)%full_var(i,j,l,nl)=bnd_vars_h%var_4d(nv)%north(ib,jb,l,1,nl)*wh(lnsh+1-jb)    &
                                                            +bnd_vars_h%var_4d(nv)%full_var(i,j,l,nl)*(1.-wh(lnsh+1-jb))
                  enddo
                enddo
              enddo
            enddo
          enddo
        endif
!
!--------------------
!***  The pint array
!--------------------
!
        do j=jde-lnsh+1,jde
          jb=j-jde+lnsh
          ihs=1-jb+lnsh
          ihe=ide+jb-lnsh
          do i=max(its_h2,ihs),min(ite_h2,ihe)
            ib=i
            pint(i,j,1)=pt
          enddo
        enddo
!
        do l=1,lm
          do j=jde-lnsh+1,jde
            jb=j-jde+lnsh
            ihs=1-jb+lnsh
            ihe=ide+jb-lnsh
            do i=max(its_h2,ihs),min(ite_h2,ihe)
              ib=i
              pint(i,j,l+1)=(pint(i,j,l) &
                            +(dsg2(l)*pd(i,j)+pdsg1(l)))*wh(lnsh+1-jb) &
                            +pint(i,j,l+1)*(1.-wh(lnsh+1-jb))
            enddo
          enddo
        enddo
!
      endif
!
!----------------------
!***  Western boundary
!----------------------
!
      if(w_bdy)then
!
        if(nvars_bc_2d_h>0)then
          do nv=1,nvars_bc_2d_h
            do i=1,lnsh
              ib=i
              jhs=1+ib
              jhe=jde-ib
              do j=max(jts_h2,jhs),min(jte_h2,jhe)
                jb=j
                bnd_vars_h%var_2d(nv)%full_var(i,j)=bnd_vars_h%var_2d(nv)%west(ib,jb,1)*wh(ib)  &
                                                   +bnd_vars_h%var_2d(nv)%full_var(i,j)*(1.-wh(ib))
              enddo
            enddo
          enddo
        endif
!
        if(nvars_bc_3d_h>0)then
          do nv=1,nvars_bc_3d_h
            do l=1,lm
              do i=1,lnsh
                ib=i
                jhs=1+ib
                jhe=jde-ib
                do j=max(jts_h2,jhs),min(jte_h2,jhe)
                  jb=j
                  bnd_vars_h%var_3d(nv)%full_var(i,j,l)=bnd_vars_h%var_3d(nv)%west(ib,jb,l,1)*wh(ib)  &
                                                       +bnd_vars_h%var_3d(nv)%full_var(i,j,l)*(1.-wh(ib))
                enddo
              enddo
            enddo
          enddo
        endif
!
        if(nvars_bc_4d_h>0)then
          do nv=1,nvars_bc_4d_h
            lbnd=lbnd_4d(nv)
            ubnd=ubnd_4d(nv)
            do nl=lbnd,ubnd
              do l=1,lm
                do i=1,lnsh
                  ib=i
                  jhs=1+ib
                  jhe=jde-ib
                  do j=max(jts_h2,jhs),min(jte_h2,jhe)
                    jb=j
                    bnd_vars_h%var_4d(nv)%full_var(i,j,l,nl)=bnd_vars_h%var_4d(nv)%west(ib,jb,l,1,nl)*wh(ib)  &
                                                            +bnd_vars_h%var_4d(nv)%full_var(i,j,l,nl)*(1.-wh(ib))
                  enddo
                enddo
              enddo
            enddo
          enddo
        endif
!
!--------------------
!***  The pint array
!--------------------
!
        do i=1,lnsh
          ib=i
          jhs=1+ib
          jhe=jde-ib
          do j=max(jts_h2,jhs),min(jte_h2,jhe)
            jb=j
            pint(i,j,1)=pt
          enddo
        enddo
!
        do l=1,lm
          do i=1,lnsh
            ib=i
            jhs=1+ib
            jhe=jde-ib
            do j=max(jts_h2,jhs),min(jte_h2,jhe)
              jb=j
              pint(i,j,l+1)=(pint(i,j,l) &
                            +(dsg2(l)*pd(i,j)+pdsg1(l)))*wh(ib) &
                            +pint(i,j,l+1)*(1.-wh(ib))
            enddo
          enddo
        enddo
!
      endif
!
!----------------------
!***  Eastern boundary
!----------------------
!
      if(e_bdy)then
!
        if(nvars_bc_2d_h>0)then
          do nv=1,nvars_bc_2d_h
            do i=ide+1-lnsh,ide
              ib=i-ide+lnsh
              jhs=2+lnsh-ib
              jhe=jde-lnsh-1+ib
              do j=max(jts_h2,jhs),min(jte_h2,jhe)
                jb=j
                bnd_vars_h%var_2d(nv)%full_var(i,j)=bnd_vars_h%var_2d(nv)%east(ib,jb,1)*wh(lnsh+1-ib)  &
                                                   +bnd_vars_h%var_2d(nv)%full_var(i,j)*(1.-wh(lnsh+1-ib))
              enddo
            enddo
          enddo
        endif
!
        if(nvars_bc_3d_h>0)then
          do nv=1,nvars_bc_3d_h
            do l=1,lm
              do i=ide+1-lnsh,ide
                ib=i-ide+lnsh
                jhs=2+lnsh-ib
                jhe=jde-lnsh-1+ib
                do j=max(jts_h2,jhs),min(jte_h2,jhe)
                  jb=j
                  bnd_vars_h%var_3d(nv)%full_var(i,j,l)=bnd_vars_h%var_3d(nv)%east(ib,jb,l,1)*wh(lnsh+1-ib)  &
                                                       +bnd_vars_h%var_3d(nv)%full_var(i,j,l)*(1.-wh(lnsh+1-ib))
                enddo
              enddo
            enddo
          enddo
        endif
!
        if(nvars_bc_4d_h>0)then
          do nv=1,nvars_bc_4d_h
            lbnd=lbnd_4d(nv)
            ubnd=ubnd_4d(nv)
            do nl=lbnd,ubnd
              do l=1,lm
                do i=ide+1-lnsh,ide
                  ib=i-ide+lnsh
                  jhs=2+lnsh-ib
                  jhe=jde-lnsh-1+ib
                  do j=max(jts_h2,jhs),min(jte_h2,jhe)
                    jb=j
                    bnd_vars_h%var_4d(nv)%full_var(i,j,l,nl)=bnd_vars_h%var_4d(nv)%east(ib,jb,l,1,nl)*wh(lnsh+1-ib)  &
                                                            +bnd_vars_h%var_4d(nv)%full_var(i,j,l,nl)*(1.-wh(lnsh+1-ib))
                  enddo
                enddo
              enddo
            enddo
          enddo
        endif
!
!--------------------
!***  The pint array
!--------------------
!
        do i=ide+1-lnsh,ide
          ib=i-ide+lnsh
          jhs=2+lnsh-ib
          jhe=jde-lnsh-1+ib
          do j=max(jts_h2,jhs),min(jte_h2,jhe)
            jb=j
            pint(i,j,1)=pt
          enddo
        enddo
!
        do l=1,lm
          do i=ide+1-lnsh,ide
            ib=i-ide+lnsh
            jhs=2+lnsh-ib
            jhe=jde-lnsh-1+ib
            do j=max(jts_h2,jhs),min(jte_h2,jhe)
              jb=j
              pint(i,j,l+1)=(pint(i,j,l) &
                            +(dsg2(l)*pd(i,j)+pdsg1(l)))*wh(lnsh+1-ib) &
                            +pint(i,j,l+1)*(1.-wh(lnsh+1-ib))
            enddo
          enddo
        enddo
!
      endif
!
!-----------------------------------------------------------------------
!***  The four corner points.
!-----------------------------------------------------------------------
!
!-----------------------------
      if(s_bdy.and.w_bdy)then
!-----------------------------
!
        if(nvars_bc_2d_h>0)then
          do nv=1,nvars_bc_2d_h
            bnd_vars_h%var_2d(nv)%full_var(its+1,jts+1)=                &
              (bnd_vars_h%var_2d(nv)%full_var(its+1,jts)                &
              +bnd_vars_h%var_2d(nv)%full_var(its,jts+1)                &
              +bnd_vars_h%var_2d(nv)%full_var(its+2,jts+1)              &
              +bnd_vars_h%var_2d(nv)%full_var(its+1,jts+2))*0.25
          enddo
        endif
!
        if(nvars_bc_3d_h>0)then
          do nv=1,nvars_bc_3d_h
            do l=1,lm
              bnd_vars_h%var_3d(nv)%full_var(its+1,jts+1,l)=            &
                (bnd_vars_h%var_3d(nv)%full_var(its+1,jts,l)            &
                +bnd_vars_h%var_3d(nv)%full_var(its,jts+1,l)            &
                +bnd_vars_h%var_3d(nv)%full_var(its+2,jts+1,l)          &
                +bnd_vars_h%var_3d(nv)%full_var(its+1,jts+2,l))*0.25
            enddo
          enddo
        endif
!
        if(nvars_bc_4d_h>0)then
          do nv=1,nvars_bc_4d_h
            lbnd=lbnd_4d(nv)
            ubnd=ubnd_4d(nv)
            do nl=lbnd,ubnd
              do l=1,lm
                bnd_vars_h%var_4d(nv)%full_var(its+1,jts+1,l,nl)=       &
                  (bnd_vars_h%var_4d(nv)%full_var(its+1,jts,l,nl)       &
                  +bnd_vars_h%var_4d(nv)%full_var(its,jts+1,l,nl)       &
                  +bnd_vars_h%var_4d(nv)%full_var(its+2,jts+1,l,nl)     &
                  +bnd_vars_h%var_4d(nv)%full_var(its+1,jts+2,l,nl))*0.25
              enddo
            enddo
          enddo
        endif
!
      endif
!
!-----------------------------
      if(s_bdy.and.e_bdy)then
!-----------------------------
!
        if(nvars_bc_2d_h>0)then
          do nv=1,nvars_bc_2d_h
            bnd_vars_h%var_2d(nv)%full_var(ite-1,jts+1)=                &
              (bnd_vars_h%var_2d(nv)%full_var(ite-1,jts)                &
              +bnd_vars_h%var_2d(nv)%full_var(ite-2,jts+1)              &
              +bnd_vars_h%var_2d(nv)%full_var(ite,jts+1)                &
              +bnd_vars_h%var_2d(nv)%full_var(ite-1,jts+2))*0.25
          enddo
        endif
!
        if(nvars_bc_3d_h>0)then
          do nv=1,nvars_bc_3d_h
            do l=1,lm
              bnd_vars_h%var_3d(nv)%full_var(ite-1,jts+1,l)=            &
                (bnd_vars_h%var_3d(nv)%full_var(ite-1,jts,l)            &
                +bnd_vars_h%var_3d(nv)%full_var(ite-2,jts+1,l)          &
                +bnd_vars_h%var_3d(nv)%full_var(ite,jts+1,l)            &
                +bnd_vars_h%var_3d(nv)%full_var(ite-1,jts+2,l))*0.25
            enddo
          enddo
        endif
!
        if(nvars_bc_4d_h>0)then
          do nv=1,nvars_bc_4d_h
            lbnd=lbnd_4d(nv)
            ubnd=ubnd_4d(nv)
            do nl=lbnd,ubnd
              do l=1,lm
                bnd_vars_h%var_4d(nv)%full_var(ite-1,jts+1,l,nl)=       &
                  (bnd_vars_h%var_4d(nv)%full_var(ite-1,jts,l,nl)       &
                  +bnd_vars_h%var_4d(nv)%full_var(ite-2,jts+1,l,nl)     &
                  +bnd_vars_h%var_4d(nv)%full_var(ite,jts+1,l,nl)       &
                  +bnd_vars_h%var_4d(nv)%full_var(ite-1,jts+2,l,nl))*0.25
              enddo
            enddo
          enddo
        endif
!
      endif
!
!-----------------------------
      if(n_bdy.and.w_bdy)then
!-----------------------------
!
        if(nvars_bc_2d_h>0)then
          do nv=1,nvars_bc_2d_h
            bnd_vars_h%var_2d(nv)%full_var(its+1,jte-1)=                &
              (bnd_vars_h%var_2d(nv)%full_var(its+1,jte-2)              &
              +bnd_vars_h%var_2d(nv)%full_var(its,jte-1)                &
              +bnd_vars_h%var_2d(nv)%full_var(its+2,jte-1)              &
              +bnd_vars_h%var_2d(nv)%full_var(its+1,jte))*0.25
          enddo
        endif
!
        if(nvars_bc_3d_h>0)then
          do nv=1,nvars_bc_3d_h
            do l=1,lm
              bnd_vars_h%var_3d(nv)%full_var(its+1,jte-1,l)=            &
                (bnd_vars_h%var_3d(nv)%full_var(its+1,jte-2,l)          &
                +bnd_vars_h%var_3d(nv)%full_var(its,jte-1,l)            &
                +bnd_vars_h%var_3d(nv)%full_var(its+2,jte-1,l)          &
                +bnd_vars_h%var_3d(nv)%full_var(its+1,jte,l))*0.25
            enddo
          enddo
        endif
!
        if(nvars_bc_4d_h>0)then
          do nv=1,nvars_bc_4d_h
            lbnd=lbnd_4d(nv)
            ubnd=ubnd_4d(nv)
            do nl=lbnd,ubnd
              do l=1,lm
                bnd_vars_h%var_4d(nv)%full_var(its+1,jte-1,l,nl)=      &
                  (bnd_vars_h%var_4d(nv)%full_var(its+1,jte-2,l,nl)    & 
                  +bnd_vars_h%var_4d(nv)%full_var(its,jte-1,l,nl)      & 
                  +bnd_vars_h%var_4d(nv)%full_var(its+2,jte-1,l,nl)    &
                  +bnd_vars_h%var_4d(nv)%full_var(its+1,jte,l,nl))*0.25
              enddo
            enddo
          enddo
        endif
!
      endif
!
!-----------------------------
      if(n_bdy.and.e_bdy)then
!-----------------------------
!
        if(nvars_bc_2d_h>0)then
          do nv=1,nvars_bc_2d_h
            bnd_vars_h%var_2d(nv)%full_var(ite-1,jte-1)=                &
              (bnd_vars_h%var_2d(nv)%full_var(ite-1,jte-2)              &
              +bnd_vars_h%var_2d(nv)%full_var(ite-2,jte-1)              &
              +bnd_vars_h%var_2d(nv)%full_var(ite,jte-1)                &
              +bnd_vars_h%var_2d(nv)%full_var(ite-1,jte))*0.25
          enddo
        endif
!
        if(nvars_bc_3d_h>0)then
          do nv=1,nvars_bc_3d_h
            do l=1,lm
              bnd_vars_h%var_3d(nv)%full_var(ite-1,jte-1,l)=            &
                (bnd_vars_h%var_3d(nv)%full_var(ite-1,jte-2,l)          &
                +bnd_vars_h%var_3d(nv)%full_var(ite-2,jte-1,l)          &
                +bnd_vars_h%var_3d(nv)%full_var(ite,jte-1,l)            &
                +bnd_vars_h%var_3d(nv)%full_var(ite-1,jte,l))*0.25
            enddo
          enddo
        endif
!
        if(nvars_bc_4d_h>0)then
          do nv=1,nvars_bc_4d_h
            lbnd=lbnd_4d(nv)
            ubnd=ubnd_4d(nv)
            do nl=lbnd,ubnd
              do l=1,lm
                bnd_vars_h%var_4d(nv)%full_var(ite-1,jte-1,l,nl)=       &
                  (bnd_vars_h%var_4d(nv)%full_var(ite-1,jte-2,l,nl)     &
                  +bnd_vars_h%var_4d(nv)%full_var(ite-2,jte-1,l,nl)     &
                  +bnd_vars_h%var_4d(nv)%full_var(ite,jte-1,l,nl)       &
                  +bnd_vars_h%var_4d(nv)%full_var(ite-1,jte,l,nl))*0.25
              enddo
            enddo
          enddo
        endif
!
      endif
!
!-----------------------------------------------------------------------
!
      if(s_bdy)then
        do j=jts,jts+1
          do i=its_h2,ite_h2
            pint(i,j,1)=pt
          enddo
        enddo
!
        do l=1,lm
          do j=jts,jts+1
            do i=its_h2,ite_h2
              pint(i,j,l+1)=pint(i,j,l)+(dsg2(l)*pd(i,j)+pdsg1(l))
            enddo
          enddo
        enddo
      endif
!
      if(n_bdy)then
        do j=jte-1,jte
          do i=its_h2,ite_h2
            pint(i,j,1)=pt
          enddo
        enddo
!
        do l=1,lm
          do j=jte-1,jte
            do i=its_h2,ite_h2
              pint(i,j,l+1)=pint(i,j,l)+(dsg2(l)*pd(i,j)+pdsg1(l))
            enddo
          enddo
        enddo
      endif
!
      if(w_bdy)then
        do j=jts_h2,jte_h2
          do i=its,its+1
            pint(i,j,1)=pt
          enddo
        enddo
!
        do l=1,lm
          do j=jts_h2,jte_h2
            do i=its,its+1
              pint(i,j,l+1)=pint(i,j,l)+(dsg2(l)*pd(i,j)+pdsg1(l))
            enddo
          enddo
        enddo
      endif
!
      if(e_bdy)then
        do j=jts_h2,jte_h2
          do i=ite-1,ite
            pint(i,j,1)=pt
          enddo
        enddo
!
        do l=1,lm
          do j=jts_h2,jte_h2
            do i=ite-1,ite
              pint(i,j,l+1)=pint(i,j,l)+(dsg2(l)*pd(i,j)+pdsg1(l))
            enddo
          enddo
        enddo
      endif
!
!-----------------------------------------------------------------------
!***  The following is executed if spatial smoothing is invoked.
!-----------------------------------------------------------------------
!
      if(nsmud>=1)then
!
        smooth: do ks=1,nsmud
!
!-----------------------------------------------------------------------
!
!-----------
!***  South
!-----------
!
          if(s_bdy)then
            i1=max(its_h2,ids+1)
            i2=min(ite_h2,ide-1)
            j1=jts+1
            j2=jts-1+lines
            allocate(hbc_2d(i1:i2,j1:j2))
!
            if(nvars_bc_2d_h>0)then
              do nv=1,nvars_bc_2d_h
                do j=j1,j2
                do i=i1,i2
                  hbc_2d(i,j)=(bnd_vars_h%var_2d(nv)%full_var(i,j-1) &
                              +bnd_vars_h%var_2d(nv)%full_var(i-1,j) &
                              +bnd_vars_h%var_2d(nv)%full_var(i+1,j) &
                              +bnd_vars_h%var_2d(nv)%full_var(i,j+1))*w1 &
                              +w2*bnd_vars_h%var_2d(nv)%full_var(i,j)
                enddo
                enddo
!
                do j=j1,j2
                do i=i1,i2
                  bnd_vars_h%var_2d(nv)%full_var(i,j)=hbc_2d(i,j)
                enddo
                enddo
              enddo
            endif
!
            if(nvars_bc_3d_h>0)then
              do nv=1,nvars_bc_3d_h
                do l=1,lm
                  do j=j1,j2
                  do i=i1,i2
                    hbc_2d(i,j)=(bnd_vars_h%var_3d(nv)%full_var(i,j-1,l) &
                                +bnd_vars_h%var_3d(nv)%full_var(i-1,j,l) &
                                +bnd_vars_h%var_3d(nv)%full_var(i+1,j,l) &
                                +bnd_vars_h%var_3d(nv)%full_var(i,j+1,l))*w1 &
                                +w2*bnd_vars_h%var_3d(nv)%full_var(i,j,l)
                  enddo
                  enddo
                  do j=j1,j2
                  do i=i1,i2
                    bnd_vars_h%var_3d(nv)%full_var(i,j,l)=hbc_2d(i,j)
                  enddo
                  enddo
                enddo
              enddo
            endif
!
            if(nvars_bc_4d_h>0)then
              do nv=1,nvars_bc_4d_h
                lbnd=lbnd_4d(nv)
                ubnd=ubnd_4d(nv)
                do nl=lbnd,ubnd
                  do l=1,lm
                    do j=j1,j2
                    do i=i1,i2
                      hbc_2d(i,j)=(bnd_vars_h%var_4d(nv)%full_var(i,j-1,l,nl) &
                                  +bnd_vars_h%var_4d(nv)%full_var(i-1,j,l,nl) &
                                  +bnd_vars_h%var_4d(nv)%full_var(i+1,j,l,nl) &
                                  +bnd_vars_h%var_4d(nv)%full_var(i,j+1,l,nl))*w1 &
                                  +w2*bnd_vars_h%var_4d(nv)%full_var(i,j,l,nl)
                    enddo
                    enddo
                    do j=j1,j2
                    do i=i1,i2
                      bnd_vars_h%var_4d(nv)%full_var(i,j,l,nl)=hbc_2d(i,j)
                    enddo
                    enddo
                  enddo
                enddo
              enddo
            endif
!
            do l=1,lm
              do j=j1,j2
              do i=i1,i2
                hbc_2d(i,j)=(pint(i,j-1,l+1)+pint(i-1,j,l+1) &
                            +pint(i+1,j,l+1)+pint(i,j+1,l+1))*w1 &
                            +w2*pint(i,j,l+1)       
              enddo
              enddo
              do j=j1,j2
              do i=i1,i2
                pint(i,j,l+1)=hbc_2d(i,j)
              enddo
              enddo
            enddo
!
            deallocate(hbc_2d)
          endif
!
!-----------
!***  North
!-----------
!
          if(n_bdy)then
            i1=max(its_h2,ids+1)
            i2=min(ite_h2,ide-1)
            j1=jte-lines+1
            j2=jte-1
            allocate(hbc_2d(i1:i2,j1:j2))
!
            if(nvars_bc_2d_h>0)then
              do nv=1,nvars_bc_2d_h
                do j=j1,j2
                do i=i1,i2
                  hbc_2d(i,j)=(bnd_vars_h%var_2d(nv)%full_var(i,j-1) &
                              +bnd_vars_h%var_2d(nv)%full_var(i-1,j) &
                              +bnd_vars_h%var_2d(nv)%full_var(i+1,j) &
                              +bnd_vars_h%var_2d(nv)%full_var(i,j+1))*w1 &
                              +w2*bnd_vars_h%var_2d(nv)%full_var(i,j)
                enddo
                enddo
!
                do j=j1,j2
                do i=i1,i2
                  bnd_vars_h%var_2d(nv)%full_var(i,j)=hbc_2d(i,j)
                enddo
                enddo
              enddo
            endif
!
            if(nvars_bc_3d_h>0)then
              do nv=1,nvars_bc_3d_h
                do l=1,lm
                  do j=j1,j2
                  do i=i1,i2
                    hbc_2d(i,j)=(bnd_vars_h%var_3d(nv)%full_var(i,j-1,l) &
                                +bnd_vars_h%var_3d(nv)%full_var(i-1,j,l) &
                                +bnd_vars_h%var_3d(nv)%full_var(i+1,j,l) &
                                +bnd_vars_h%var_3d(nv)%full_var(i,j+1,l))*w1 &
                                +w2*bnd_vars_h%var_3d(nv)%full_var(i,j,l)
                  enddo
                  enddo
                  do j=j1,j2
                  do i=i1,i2
                    bnd_vars_h%var_3d(nv)%full_var(i,j,l)=hbc_2d(i,j)
                  enddo
                  enddo
                enddo
              enddo
            endif
!
            if(nvars_bc_4d_h>0)then
              do nv=1,nvars_bc_4d_h
                lbnd=lbnd_4d(nv)
                ubnd=ubnd_4d(nv)
                do nl=lbnd,ubnd
                  do l=1,lm
                    do j=j1,j2
                    do i=i1,i2
                      hbc_2d(i,j)=(bnd_vars_h%var_4d(nv)%full_var(i,j-1,l,nl) &
                                  +bnd_vars_h%var_4d(nv)%full_var(i-1,j,l,nl) &
                                  +bnd_vars_h%var_4d(nv)%full_var(i+1,j,l,nl) &
                                  +bnd_vars_h%var_4d(nv)%full_var(i,j+1,l,nl))*w1 &
                                  +w2*bnd_vars_h%var_4d(nv)%full_var(i,j,l,nl)
                    enddo
                    enddo
                    do j=j1,j2
                    do i=i1,i2
                      bnd_vars_h%var_4d(nv)%full_var(i,j,l,nl)=hbc_2d(i,j)
                    enddo
                    enddo
                  enddo
                enddo
              enddo
            endif
!
            do l=1,lm
              do j=j1,j2
              do i=i1,i2
                hbc_2d(i,j)=(pint(i,j-1,l+1)+pint(i-1,j,l+1) &
                            +pint(i+1,j,l+1)+pint(i,j+1,l+1))*w1 &
                            +w2*pint(i,j,l+1)       
              enddo
              enddo
              do j=j1,j2
              do i=i1,i2
                pint(i,j,l+1)=hbc_2d(i,j)
              enddo
              enddo
            enddo
!
            deallocate(hbc_2d)
          endif
!
!----------
!***  West
!----------
!
          if(w_bdy)then
            i1=its+1
            i2=its-1+lines
            j1=max(jts_h2,jds+lines)
            j2=min(jte_h2,jde-lines)
            allocate(hbc_2d(i1:i2,j1:j2))
!
            if(nvars_bc_2d_h>0)then
              do nv=1,nvars_bc_2d_h
                do j=j1,j2
                do i=i1,i2
                  hbc_2d(i,j)=(bnd_vars_h%var_2d(nv)%full_var(i,j-1) &
                              +bnd_vars_h%var_2d(nv)%full_var(i-1,j) &
                              +bnd_vars_h%var_2d(nv)%full_var(i+1,j) &
                              +bnd_vars_h%var_2d(nv)%full_var(i,j+1))*w1 &
                              +w2*bnd_vars_h%var_2d(nv)%full_var(i,j)
                enddo
                enddo
!
                do j=j1,j2
                do i=i1,i2
                  bnd_vars_h%var_2d(nv)%full_var(i,j)=hbc_2d(i,j)
                enddo
                enddo
              enddo
            endif
!
            if(nvars_bc_3d_h>0)then
              do nv=1,nvars_bc_3d_h
                do l=1,lm
                  do j=j1,j2
                  do i=i1,i2
                    hbc_2d(i,j)=(bnd_vars_h%var_3d(nv)%full_var(i,j-1,l) &
                                +bnd_vars_h%var_3d(nv)%full_var(i-1,j,l) &
                                +bnd_vars_h%var_3d(nv)%full_var(i+1,j,l) &
                                +bnd_vars_h%var_3d(nv)%full_var(i,j+1,l))*w1 &
                                +w2*bnd_vars_h%var_3d(nv)%full_var(i,j,l)
                  enddo
                  enddo
                  do j=j1,j2
                  do i=i1,i2
                    bnd_vars_h%var_3d(nv)%full_var(i,j,l)=hbc_2d(i,j)
                  enddo
                  enddo
                enddo
              enddo
            endif
!
            if(nvars_bc_4d_h>0)then
              do nv=1,nvars_bc_4d_h
                lbnd=lbnd_4d(nv)
                ubnd=ubnd_4d(nv)
                do nl=lbnd,ubnd
                  do l=1,lm
                    do j=j1,j2
                    do i=i1,i2
                      hbc_2d(i,j)=(bnd_vars_h%var_4d(nv)%full_var(i,j-1,l,nl) &
                                  +bnd_vars_h%var_4d(nv)%full_var(i-1,j,l,nl) &
                                  +bnd_vars_h%var_4d(nv)%full_var(i+1,j,l,nl) &
                                  +bnd_vars_h%var_4d(nv)%full_var(i,j+1,l,nl))*w1 &
                                  +w2*bnd_vars_h%var_4d(nv)%full_var(i,j,l,nl)
                    enddo
                    enddo
                    do j=j1,j2
                    do i=i1,i2
                      bnd_vars_h%var_4d(nv)%full_var(i,j,l,nl)=hbc_2d(i,j)
                    enddo
                    enddo
                  enddo
                enddo
              enddo
            endif
!
            do l=1,lm
              do j=j1,j2
              do i=i1,i2
                hbc_2d(i,j)=(pint(i,j-1,l+1)+pint(i-1,j,l+1) &
                            +pint(i+1,j,l+1)+pint(i,j+1,l+1))*w1 &
                            +w2*pint(i,j,l+1)       
              enddo
              enddo
              do j=j1,j2
              do i=i1,i2
                pint(i,j,l+1)=hbc_2d(i,j)
              enddo
              enddo
            enddo
!
            deallocate(hbc_2d)
          endif
!
!----------
!***  East
!----------
!
          if(e_bdy)then
            i1=ite-lines+1
            i2=ite-1
            j1=max(jts_h2,jds+lines)
            j2=min(jte_h2,jde-lines)
            allocate(hbc_2d(i1:i2,j1:j2))
!
            if(nvars_bc_2d_h>0)then
              do nv=1,nvars_bc_2d_h
                do j=j1,j2
                do i=i1,i2
                  hbc_2d(i,j)=(bnd_vars_h%var_2d(nv)%full_var(i,j-1) &
                              +bnd_vars_h%var_2d(nv)%full_var(i-1,j) &
                              +bnd_vars_h%var_2d(nv)%full_var(i+1,j) &
                              +bnd_vars_h%var_2d(nv)%full_var(i,j+1))*w1 &
                              +w2*bnd_vars_h%var_2d(nv)%full_var(i,j)
                enddo
                enddo
!
                do j=j1,j2
                do i=i1,i2
                  bnd_vars_h%var_2d(nv)%full_var(i,j)=hbc_2d(i,j)
                enddo
                enddo
              enddo
            endif
!
            if(nvars_bc_3d_h>0)then
              do nv=1,nvars_bc_3d_h
                do l=1,lm
                  do j=j1,j2
                  do i=i1,i2
                    hbc_2d(i,j)=(bnd_vars_h%var_3d(nv)%full_var(i,j-1,l) &
                                +bnd_vars_h%var_3d(nv)%full_var(i-1,j,l) &
                                +bnd_vars_h%var_3d(nv)%full_var(i+1,j,l) &
                                +bnd_vars_h%var_3d(nv)%full_var(i,j+1,l))*w1 &
                                +w2*bnd_vars_h%var_3d(nv)%full_var(i,j,l)
                  enddo
                  enddo
                  do j=j1,j2
                  do i=i1,i2
                    bnd_vars_h%var_3d(nv)%full_var(i,j,l)=hbc_2d(i,j)
                  enddo
                  enddo
                enddo
              enddo
            endif
!
            if(nvars_bc_4d_h>0)then
              do nv=1,nvars_bc_4d_h
                lbnd=lbnd_4d(nv)
                ubnd=ubnd_4d(nv)
                do nl=lbnd,ubnd
                  do l=1,lm
                    do j=j1,j2
                    do i=i1,i2
                      hbc_2d(i,j)=(bnd_vars_h%var_4d(nv)%full_var(i,j-1,l,nl) &
                                  +bnd_vars_h%var_4d(nv)%full_var(i-1,j,l,nl) &
                                  +bnd_vars_h%var_4d(nv)%full_var(i+1,j,l,nl) &
                                  +bnd_vars_h%var_4d(nv)%full_var(i,j+1,l,nl))*w1 &
                                  +w2*bnd_vars_h%var_4d(nv)%full_var(i,j,l,nl)
                    enddo
                    enddo
                    do j=j1,j2
                    do i=i1,i2
                      bnd_vars_h%var_4d(nv)%full_var(i,j,l,nl)=hbc_2d(i,j)
                    enddo
                    enddo
                  enddo
                enddo
              enddo
            endif
!
            do l=1,lm
              do j=j1,j2
              do i=i1,i2
                hbc_2d(i,j)=(pint(i,j-1,l+1)+pint(i-1,j,l+1) &
                            +pint(i+1,j,l+1)+pint(i,j+1,l+1))*w1 &
                            +w2*pint(i,j,l+1)       
              enddo
              enddo
              do j=j1,j2
              do i=i1,i2
                pint(i,j,l+1)=hbc_2d(i,j)
              enddo
              enddo
            enddo
!
            deallocate(hbc_2d)
          endif
!
!-----------------------------------------------------------------------
!
          if(nvars_bc_2d_h>0)then
            do nv=1,nvars_bc_2d_h
              call halo_exch(bnd_vars_h%var_2d(nv)%full_var,1,1,1)
            enddo
          endif
!
          if(nvars_bc_3d_h>0)then
            do nv=1,nvars_bc_3d_h
              call halo_exch(bnd_vars_h%var_3d(nv)%full_var,lm,1,1)
            enddo
          endif
!
          if(nvars_bc_4d_h>0)then
            do nv=1,nvars_bc_4d_h
              lbnd=lbnd_4d(nv)
              ubnd=ubnd_4d(nv)
              do nl=lbnd,ubnd
                temp_3d=>bnd_vars_h%var_4d(nv)%full_var(:,:,:,nl)     
                call halo_exch(temp_3d,lm,1,1)
              enddo
            enddo
          endif
!
          call halo_exch(pint,lm+1,1,1)
!
!-----------------------------------------------------------------------
!
        enddo smooth
!
      endif
!
!-----------------------------------------------------------------------
!
                        endsubroutine bocoh
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
                        subroutine bocov &
(lm,lnsv &
,dt &
,nvars_bc_2d_v,nvars_bc_3d_v &
,bnd_vars_v &
   )
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
implicit none
!-----------------------------------------------------------------------

real(kind=kfpt),parameter:: &
 wa=0.5 &                    ! weighting factor
,w1=wa*0.25 &                ! weighting factor
,w2=1.-wa                    ! weighting factor

!-----------------------------------------------------------------------

integer(kind=kint),intent(in):: &
 lm &                        ! total # of levels
,lnsv                        ! blending area width, v points

integer(kind=kint),intent(in):: &
 nvars_bc_2d_v &             ! # of 2-d v-pt boundary variables
,nvars_bc_3d_v               ! # of 3-d v-pt boundary variables

real(kind=kfpt),intent(in):: &
 dt                          ! dynamics time step

type(bc_v_all),intent(inout) :: &
 bnd_vars_v                  ! boundary variables on v points

!-----------------------------------------------------------------------
!---local variables-----------------------------------------------------
!-----------------------------------------------------------------------
integer(kind=kint):: &
 i &                         ! index in x direction
,i1 &
,i2 &
,ib &                        ! index in x direction, boundary zone
,ive &                       ! ending index in x direction, boundaries
,ivs &                       ! starting index in x direction, boundaries
,j &                         ! index in y direction
,j1 &
,j2 &
,jb &                        ! index in y direction, boundary zone
,jve &                       ! ending index in x direction, boundaries
,jvs &                       ! starting index in x direction, boundaries
,k &                         ! boundary line counter
,ks &                        ! smoothing counter
,l &                         ! index in p direction
,lines &                     ! boundary smoothing area
,nsmud &                     ! number of smoothing passes
,nv

real(kind=kfpt),dimension(1:lnsv):: &
 wv(lnsv)                    ! blending weighting function, wind

real(kind=kfpt),dimension(:,:),allocatable :: &
 vbc_2d
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
      lines=lnsv
!     nsmud=lines-1
      nsmud=0
!-----------------------------------------------------------------------
      wv(1)=1.
      do k=2,lnsv
        wv(k)=1.-(0.9/real(lnsv-1))*(k-1)
      enddo
!-----------------------------------------------------------------------
!***  Update values of the boundary working objects at V points.
!-----------------------------------------------------------------------
!
!--------------------------------------
!***  Southern and northern boundaries
!--------------------------------------
!
      if(s_bdy)then
!
        if(nvars_bc_2d_v>0)then
          do nv=1,nvars_bc_2d_v
            do jb=1,lnsv
            do ib=its_h2,min(ite_h2,ide-1)
              bnd_vars_v%var_2d(nv)%south(ib,jb,1)=bnd_vars_v%var_2d(nv)%south(ib,jb,1)   &
                                                  +bnd_vars_v%var_2d(nv)%south(ib,jb,2)*dt
            enddo
            enddo
          enddo
        endif
!
        if(nvars_bc_3d_v>0)then
          do nv=1,nvars_bc_3d_v
            do l=1,lm
            do jb=1,lnsv
            do ib=its_h2,min(ite_h2,ide-1)
              bnd_vars_v%var_3d(nv)%south(ib,jb,l,1)=bnd_vars_v%var_3d(nv)%south(ib,jb,l,1)   &
                                                    +bnd_vars_v%var_3d(nv)%south(ib,jb,l,2)*dt
            enddo
            enddo
            enddo
          enddo
        endif
!
      endif
!
      if(n_bdy)then
!
        if(nvars_bc_2d_v>0)then
          do nv=1,nvars_bc_2d_v
            do jb=1,lnsv
            do ib=its_h2,min(ite_h2,ide-1)
              bnd_vars_v%var_2d(nv)%north(ib,jb,1)=bnd_vars_v%var_2d(nv)%north(ib,jb,1)   &
                                                  +bnd_vars_v%var_2d(nv)%north(ib,jb,2)*dt
            enddo
            enddo
          enddo
        endif
!
        if(nvars_bc_3d_v>0)then
          do nv=1,nvars_bc_3d_v
            do l=1,lm
            do jb=1,lnsv
            do ib=its_h2,min(ite_h2,ide-1)
              bnd_vars_v%var_3d(nv)%north(ib,jb,l,1)=bnd_vars_v%var_3d(nv)%north(ib,jb,l,1)   &
                                                    +bnd_vars_v%var_3d(nv)%north(ib,jb,l,2)*dt
            enddo
            enddo
            enddo
          enddo
        endif
!
      endif
!
!------------------------------------
!***  Western and eastern boundaries
!------------------------------------
!
      if(w_bdy)then
!
        if(nvars_bc_2d_v>0)then
          do nv=1,nvars_bc_2d_v
            do jb=jts_h2,min(jte_h2,jde-1)
            do ib=1,lnsv
              bnd_vars_v%var_2d(nv)%west(ib,jb,1)=bnd_vars_v%var_2d(nv)%west(ib,jb,1)   &
                                                 +bnd_vars_v%var_2d(nv)%west(ib,jb,2)*dt
            enddo
            enddo
          enddo
        endif
!
        if(nvars_bc_3d_v>0)then
          do nv=1,nvars_bc_3d_v
            do l=1,lm
            do jb=jts_h2,min(jte_h2,jde-1)
            do ib=1,lnsv
              bnd_vars_v%var_3d(nv)%west(ib,jb,l,1)=bnd_vars_v%var_3d(nv)%west(ib,jb,l,1)   &
                                                   +bnd_vars_v%var_3d(nv)%west(ib,jb,l,2)*dt
            enddo
            enddo
            enddo
          enddo
        endif
!
      endif
!
      if(e_bdy)then
!
        if(nvars_bc_2d_v>0)then
          do nv=1,nvars_bc_2d_v
            do jb=jts_h1,min(jte_h2,jde-1)
            do ib=1,lnsv
              bnd_vars_v%var_2d(nv)%east(ib,jb,1)=bnd_vars_v%var_2d(nv)%east(ib,jb,1)   &
                                                 +bnd_vars_v%var_2d(nv)%east(ib,jb,2)*dt
            enddo
            enddo
          enddo
        endif
!
        if(nvars_bc_3d_v>0)then
          do nv=1,nvars_bc_3d_v
            do l=1,lm
            do jb=jts_h1,min(jte_h2,jde-1)
            do ib=1,lnsv
              bnd_vars_v%var_3d(nv)%east(ib,jb,l,1)=bnd_vars_v%var_3d(nv)%east(ib,jb,l,1)   &
                                                   +bnd_vars_v%var_3d(nv)%east(ib,jb,l,2)*dt
            enddo
            enddo
            enddo
          enddo
        endif
!
      endif
!
!-----------------------------------------------------------------------
!***  Now update the actual prognostic variables using the values
!***  just computed (coming from the parent) blended with the values
!***  already generated in the blending region by the child.
!-----------------------------------------------------------------------
!
!-----------------------
!***  Southern boundary
!-----------------------
!
      if(s_bdy)then
!
        if(nvars_bc_2d_v>0)then
          do nv=1,nvars_bc_2d_v
            do j=jts,jts-1+lnsv
              jb=j
              ivs=max(its_h1,jb)
              ive=min(ite_h1,ide-jb)
              do i=ivs,ive
                ib=i
                bnd_vars_v%var_2d(nv)%full_var(i,j)=bnd_vars_v%var_2d(nv)%south(ib,jb,1)*wv(jb)    &
                                                   +bnd_vars_v%var_2d(nv)%full_var(i,j)*(1.-wv(jb))
              enddo
            enddo
          enddo
        endif
!
        if(nvars_bc_3d_v>0)then
          do nv=1,nvars_bc_3d_v
            do l=1,lm
              do j=jts,jts-1+lnsv
                jb=j
                ivs=max(its_h1,jb)
                ive=min(ite_h1,ide-jb)
                do i=ivs,ive
                  ib=i
                  bnd_vars_v%var_3d(nv)%full_var(i,j,l)=bnd_vars_v%var_3d(nv)%south(ib,jb,l,1)*wv(jb)    &
                                                       +bnd_vars_v%var_3d(nv)%full_var(i,j,l)*(1.-wv(jb))
                enddo
              enddo
            enddo
          enddo
        endif
!
      endif
!
!-----------------------
!***  Northern boundary
!-----------------------
!
      if(n_bdy)then
!
        if(nvars_bc_2d_v>0)then
          do nv=1,nvars_bc_2d_v
            do j=jte-lnsv,jte-1
              jb=j-jte+lnsv+1
              ivs=max(its_h1,lnsv-jb+1)
              ive=min(ite_h1,ide+jb-lnsv-1)
              do i=ivs,ive
                ib=i
                bnd_vars_v%var_2d(nv)%full_var(i,j)=bnd_vars_v%var_2d(nv)%north(ib,jb,1)*wv(lnsv+1-jb)    &
                                                   +bnd_vars_v%var_2d(nv)%full_var(i,j)*(1.-wv(lnsv+1-jb))
              enddo
            enddo
          enddo
        endif
!
        if(nvars_bc_3d_v>0)then
          do nv=1,nvars_bc_3d_v
            do l=1,lm
              do j=jte-lnsv,jte-1
                jb=j-jte+lnsv+1
                ivs=max(its_h1,lnsv-jb+1)
                ive=min(ite_h1,ide+jb-lnsv-1)
                do i=max(its_h2,ivs),min(ite_h2,ive)
                  ib=i
                  bnd_vars_v%var_3d(nv)%full_var(i,j,l)=bnd_vars_v%var_3d(nv)%north(ib,jb,l,1)*wv(lnsv+1-jb)    &
                                                       +bnd_vars_v%var_3d(nv)%full_var(i,j,l)*(1.-wv(lnsv+1-jb))
                enddo
              enddo
            enddo
          enddo
        endif
!
      endif
!
!
!----------------------
!***  Western boundary
!----------------------
!
      if(w_bdy)then
!
        if(nvars_bc_2d_v>0)then
          do nv=1,nvars_bc_2d_v
            do i=its,its-1+lnsv
              ib=i
              jvs=max(jts_h1,1+ib)
              jve=min(jte_h1,jde-1-ib)
              do j=jvs,jve
                jb=j
                bnd_vars_v%var_2d(nv)%full_var(i,j)=bnd_vars_v%var_2d(nv)%west(ib,jb,1)*wv(ib)  &
                                                   +bnd_vars_v%var_2d(nv)%full_var(i,j)*(1.-wv(ib))
              enddo
            enddo
          enddo
        endif
!
        if(nvars_bc_3d_v>0)then
          do nv=1,nvars_bc_3d_v
            do l=1,lm
              do i=its,its-1+lnsv
                ib=i
                jvs=max(jts_h1,1+ib)
                jve=min(jte_h1,jde-1-ib)
                do j=jvs,jve
                  jb=j
                  bnd_vars_v%var_3d(nv)%full_var(i,j,l)=bnd_vars_v%var_3d(nv)%west(ib,jb,l,1)*wv(ib)  &
                                                       +bnd_vars_v%var_3d(nv)%full_var(i,j,l)*(1.-wv(ib))
                enddo
              enddo
            enddo
          enddo
        endif
!
      endif
!
!----------------------
!***  Eastern boundary
!----------------------
!
      if(e_bdy)then
!
        if(nvars_bc_2d_v>0)then
          do nv=1,nvars_bc_2d_v
            do i=ite-lnsv,ite-1
              ib=i-ide+lnsv+1
              jvs=max(jts_h1,lnsv-ib+2)
              jve=min(jte_h1,jde-lnsv+ib-2)
              do j=jvs,jve
                jb=j
                bnd_vars_v%var_2d(nv)%full_var(i,j)=bnd_vars_v%var_2d(nv)%east(ib,jb,1)*wv(lnsv+1-ib)  &
                                                   +bnd_vars_v%var_2d(nv)%full_var(i,j)*(1.-wv(lnsv+1-ib))
              enddo
            enddo
          enddo
        endif
!
        if(nvars_bc_3d_v>0)then
          do nv=1,nvars_bc_3d_v
            do l=1,lm
              do i=ite-lnsv,ite-1
                ib=i-ide+lnsv+1
                jvs=max(jts_h1,lnsv-ib+2)
                jve=min(jte_h1,jde-lnsv+ib-2)
                do j=jvs,jve
                  jb=j
                  bnd_vars_v%var_3d(nv)%full_var(i,j,l)=bnd_vars_v%var_3d(nv)%east(ib,jb,l,1)*wv(lnsv+1-ib)  &
                                                       +bnd_vars_v%var_3d(nv)%full_var(i,j,l)*(1.-wv(lnsv+1-ib))
                enddo
              enddo
            enddo
          enddo
        endif
!
      endif
!
!-----------------------------------------------------------------------
!***  The following is executed if spatial smoothing is invoked.
!-----------------------------------------------------------------------
!
      if(nsmud>=1)then
!
        smooth: do ks=1,nsmud
!
!-----------------------------------------------------------------------
!
!-----------
!***  South
!-----------
!
          if(s_bdy)then
            i1=max(its_h1,ids+1)
            i2=min(ite_h1,ide-2)
            j1=jts+1
            j2=jts-1+lines
            allocate(vbc_2d(i1:i2,j1:j2))
!
            if(nvars_bc_2d_v>0)then
              do nv=1,nvars_bc_2d_v
                do j=j1,j2
                do i=i1,i2
                  vbc_2d(i,j)=(bnd_vars_v%var_2d(nv)%full_var(i,j-1) &
                              +bnd_vars_v%var_2d(nv)%full_var(i-1,j) &
                              +bnd_vars_v%var_2d(nv)%full_var(i+1,j) &
                              +bnd_vars_v%var_2d(nv)%full_var(i,j+1))*w1 &
                              +w2*bnd_vars_v%var_2d(nv)%full_var(i,j)
                enddo
                enddo
!
                do j=j1,j2
                do i=i1,i2
                  bnd_vars_v%var_2d(nv)%full_var(i,j)=vbc_2d(i,j)
                enddo
                enddo
              enddo
            endif
!
            if(nvars_bc_3d_v>0)then
              do nv=1,nvars_bc_3d_v
                do l=1,lm
                  do j=j1,j2
                  do i=i1,i2
                    vbc_2d(i,j)=(bnd_vars_v%var_3d(nv)%full_var(i,j-1,l) &
                                +bnd_vars_v%var_3d(nv)%full_var(i-1,j,l) &
                                +bnd_vars_v%var_3d(nv)%full_var(i+1,j,l) &
                                +bnd_vars_v%var_3d(nv)%full_var(i,j+1,l))*w1 &
                                +w2*bnd_vars_v%var_3d(nv)%full_var(i,j,l)
                  enddo
                  enddo
                  do j=j1,j2
                  do i=i1,i2
                    bnd_vars_v%var_3d(nv)%full_var(i,j,l)=vbc_2d(i,j)
                  enddo
                  enddo
                enddo
              enddo
            endif
!
            deallocate(vbc_2d)
          endif
!
!-----------
!***  North
!-----------
!
          if(n_bdy)then
            i1=max(its_h1,ids+1)
            i2=min(ite_h1,ide-2)
            j1=jte-lines
            j2=jte-2
            allocate(vbc_2d(i1:i2,j1:j2))
!
            if(nvars_bc_2d_v>0)then
              do nv=1,nvars_bc_2d_v
                do j=j1,j2
                do i=i1,i2
                  vbc_2d(i,j)=(bnd_vars_v%var_2d(nv)%full_var(i,j-1) &
                              +bnd_vars_v%var_2d(nv)%full_var(i-1,j) &
                              +bnd_vars_v%var_2d(nv)%full_var(i+1,j) &
                              +bnd_vars_v%var_2d(nv)%full_var(i,j+1))*w1 &
                              +w2*bnd_vars_v%var_2d(nv)%full_var(i,j)
                enddo
                enddo
!
                do j=j1,j2
                do i=i1,i2
                  bnd_vars_v%var_2d(nv)%full_var(i,j)=vbc_2d(i,j)
                enddo
                enddo
              enddo
            endif
!
            if(nvars_bc_3d_v>0)then
              do nv=1,nvars_bc_3d_v
                do l=1,lm
                  do j=j1,j2
                  do i=i1,i2
                    vbc_2d(i,j)=(bnd_vars_v%var_3d(nv)%full_var(i,j-1,l) &
                                +bnd_vars_v%var_3d(nv)%full_var(i-1,j,l) &
                                +bnd_vars_v%var_3d(nv)%full_var(i+1,j,l) &
                                +bnd_vars_v%var_3d(nv)%full_var(i,j+1,l))*w1 &
                                +w2*bnd_vars_v%var_3d(nv)%full_var(i,j,l)
                  enddo
                  enddo
                  do j=j1,j2
                  do i=i1,i2
                    bnd_vars_v%var_3d(nv)%full_var(i,j,l)=vbc_2d(i,j)
                  enddo
                  enddo
                enddo
              enddo
            endif
!
            deallocate(vbc_2d)
          endif
!
!----------
!***  West
!----------
!
          if(w_bdy)then
            i1=ite-lines
            i2=ite-2
            j1=max(jts_h1,jds+lines)
            j2=min(jte_h1,jde-lines-1)
            allocate(vbc_2d(i1:i2,j1:j2))
!
            if(nvars_bc_2d_v>0)then
              do nv=1,nvars_bc_2d_v
                do j=j1,j2
                do i=i1,i2
                  vbc_2d(i,j)=(bnd_vars_v%var_2d(nv)%full_var(i,j-1) &
                              +bnd_vars_v%var_2d(nv)%full_var(i-1,j) &
                              +bnd_vars_v%var_2d(nv)%full_var(i+1,j) &
                              +bnd_vars_v%var_2d(nv)%full_var(i,j+1))*w1 &
                              +w2*bnd_vars_v%var_2d(nv)%full_var(i,j)
                enddo
                enddo
!
                do j=j1,j2
                do i=i1,i2
                  bnd_vars_v%var_2d(nv)%full_var(i,j)=vbc_2d(i,j)
                enddo
                enddo
              enddo
            endif
!
            if(nvars_bc_3d_v>0)then
              do nv=1,nvars_bc_3d_v
                do l=1,lm
                  do j=j1,j2
                  do i=i1,i2
                    vbc_2d(i,j)=(bnd_vars_v%var_3d(nv)%full_var(i,j-1,l) &
                                +bnd_vars_v%var_3d(nv)%full_var(i-1,j,l) &
                                +bnd_vars_v%var_3d(nv)%full_var(i+1,j,l) &
                                +bnd_vars_v%var_3d(nv)%full_var(i,j+1,l))*w1 &
                                +w2*bnd_vars_v%var_3d(nv)%full_var(i,j,l)
                  enddo
                  enddo
!
                  do j=j1,j2
                  do i=i1,i2
                    bnd_vars_v%var_3d(nv)%full_var(i,j,l)=vbc_2d(i,j)
                  enddo
                  enddo
                enddo
              enddo
            endif
!
            deallocate(vbc_2d)
          endif
!
!----------
!***  East
!----------
!
          if(e_bdy)then
            i1=ite-lines
            i2=ite-2
            j1=max(jts_h1,jds+lines)
            j2=min(jte_h1,jde-lines-1)
            allocate(vbc_2d(i1:i2,j1:j2))
!
            if(nvars_bc_2d_v>0)then
              do nv=1,nvars_bc_2d_v
                do j=j1,j2
                do i=i1,i2
                  vbc_2d(i,j)=(bnd_vars_v%var_2d(nv)%full_var(i,j-1) &
                              +bnd_vars_v%var_2d(nv)%full_var(i-1,j) &
                              +bnd_vars_v%var_2d(nv)%full_var(i+1,j) &
                              +bnd_vars_v%var_2d(nv)%full_var(i,j+1))*w1 &
                              +w2*bnd_vars_v%var_2d(nv)%full_var(i,j)
                enddo
                enddo
!
                do j=j1,j2
                do i=i1,i2
                  bnd_vars_v%var_2d(nv)%full_var(i,j)=vbc_2d(i,j)
                enddo
                enddo
              enddo
            endif
!
            if(nvars_bc_3d_v>0)then
              do nv=1,nvars_bc_3d_v
                do l=1,lm
                  do j=j1,j2
                  do i=i1,i2
                    vbc_2d(i,j)=(bnd_vars_v%var_3d(nv)%full_var(i,j-1,l) &
                                +bnd_vars_v%var_3d(nv)%full_var(i-1,j,l) &
                                +bnd_vars_v%var_3d(nv)%full_var(i+1,j,l) &
                                +bnd_vars_v%var_3d(nv)%full_var(i,j+1,l))*w1 &
                                +w2*bnd_vars_v%var_3d(nv)%full_var(i,j,l)
                  enddo
                  enddo
!
                  do j=j1,j2
                  do i=i1,i2
                    bnd_vars_v%var_3d(nv)%full_var(i,j,l)=vbc_2d(i,j)
                  enddo
                  enddo
                enddo
              enddo
            endif
!
            deallocate(vbc_2d)
          endif
!
!-----------------------------------------------------------------------
!
          if(nvars_bc_2d_v>0)then
            do nv=1,nvars_bc_2d_v
              call halo_exch(bnd_vars_v%var_2d(nv)%full_var,1,1,1)
            enddo
          endif
!
          if(nvars_bc_3d_v>0)then
            do nv=1,nvars_bc_3d_v
              call halo_exch(bnd_vars_v%var_3d(nv)%full_var,lm,1,1)
            enddo
          endif
!
!-----------------------------------------------------------------------
!
        enddo smooth
!
!-----------------------------------------------------------------------
      endif
!-----------------------------------------------------------------------
!
                        endsubroutine bocov
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
                        endmodule module_fltbnds
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

