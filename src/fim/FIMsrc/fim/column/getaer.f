      subroutine getaer(nf,kprfg,idxcg,cmixg,denng
     &,                 nxc, ndn, imxae, jmxae, IMON)
!****************************************************************
!   read in opac global aerosol data set (1998) for sw radiation
!
!   input variables:
!     nf        - input file unit nuber
!
!   output variables:
!     kprfg     - aerosol profile type index
!     idxcg     - aerosol components types indices
!     cmixg     - aerosol components mixing ratioes
!     denng     - first two layers aerosol number densities
!
!****************************************************************
      use machine
      implicit none
! - input variable:
      integer nf, nxc, ndn, imxae, jmxae, imon
! - output variables:
      integer idxcg(nxc,imxae,jmxae), kprfg(imxae,jmxae)
      real (kind=kind_io8) cmixg(nxc,imxae,jmxae),denng(ndn,imxae,jmxae)
! - local variables:
      integer idxc(nxc), kprf, i, j, k, nc
      real (kind=kind_io8)    cmix(nxc), denn
      real (kind=kind_io8)    temp
      character cline*80, ctyp*3, aerosol_file*40
!
      write(aerosol_file,101) imon
  101 format('aeropac3a.m',i2.2)
!
      open (unit=nf, file=aerosol_file, status='OLD', form='FORMATTED')
!
      read(nf,10) cline
  10  format(a80)
!
      do j=1,jmxae
        do i=1,imxae
          read(nf,20) (idxc(k),cmix(k),k=1,nxc),kprf,denn,nc,ctyp
  20      format(5(i2,e11.4),i2,f8.2,i3,1x,a3)
!
          kprfg(i,j)     = kprf
          denng(1,i,j)   = denn       ! num density of 1st layer
          if (kprf .ge. 6) then
            denng(2,i,j) = cmix(nxc)  ! num density of 2dn layer
          else
            denng(2,i,j) = 0.0
          end if
!
          temp = 1.0
          do k=1,nxc-1
            idxcg(k,i,j) = idxc(k)    ! component index
            cmixg(k,i,j) = cmix(k)    ! component mixing ratio
            temp         = temp - cmix(k)
          end do
          idxcg(nxc,i,j) = idxc(nxc)
          cmixg(nxc,i,j) = temp       ! to make sure all add to 1.
        end do
      end do
      close(nf)
!
      return
      end
