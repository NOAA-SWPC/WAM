      subroutine GET_PRS_gc(im,ix,levs,rkap,rd,nu,
     &                   t,q,prsi,prki,prsl,prkl,phii,phil,del)
!
! hmhj : this is modified hybrid by finite difference from henry
!
      USE MACHINE , ONLY : kind_phys
      implicit none
!
      integer im, ix, levs
      real(kind=kind_phys) rkap, rd, nu
      real(kind=kind_phys) prsi(ix,levs+1), prki(ix,levs+1)
     &,                    phii(ix,levs+1), phil(ix,levs)
     &,                    prsl(ix,levs),   prkl(ix,levs)
     *,                    del(ix,levs),    T(ix,levs),   q(ix,levs)
      real(kind=kind_phys) tem, dphi
      integer i, k
cc
cc--------------------------------------------------------------------
cc
      real(kind=kind_phys) cons_0               !constant
cc

!     print *,' enter get_prs_v_gc_fd '

      cons_0          =         0.d0            !constant
cc
cc--------------------------------------------------------------------
cc
!                                    Pressure is in centibars!!!!
      do k=1,levs
        do i=1,im
          del(i,k) = PRSI(i,k) - PRSI(i,k+1)
          prsl(i,k) =(PRSI(i,k) + PRSI(i,k+1))*0.5
          prkl(i,k) =(prsl(i,k)*0.01) ** rkap
        enddo
      enddo
      do k=1,levs+1
      do i=1,im
        prki(i,k) = (prsi(i,k)*0.01) ** rkap
      enddo
      enddo
      do i=1,im
        phii(i,1)   = 0.0           ! Ignoring topography height here
      enddo
      DO k=1,levs
        do i=1,im
          TEM         = rd * T(i,k) * (1.0 + NU * max(Q(i,k),cons_0))     !constant
          DPHI        = (PRSI(i,k) - PRSI(i,k+1)) * TEM
     &                 /(PRSI(i,k) + PRSI(i,k+1))
          phil(i,k)   = phii(i,k) + DPHI
          phii(i,k+1) = phil(i,k) + DPHI
        ENDDO
      ENDDO
!
      return
      end
      subroutine GET_PHI_gc(im,ix,levs,rd,nu,t,q,prsi,phii,phil)
!
      USE MACHINE , ONLY : kind_phys
      implicit none
!
      integer im, ix, levs
      real(kind=kind_phys) rd, nu
      real(kind=kind_phys) prsi(ix,levs+1)
     &,                    phii(ix,levs+1), phil(ix,levs)
     *,                    T(ix,levs),      q(ix,levs)
      real(kind=kind_phys) tem, dphi
      integer i, k
cc
cc--------------------------------------------------------------------
cc
      real(kind=kind_phys) cons_0               !constant
cc
      cons_0          =         0.d0            !constant
cc
cc--------------------------------------------------------------------
cc
      do i=1,im
        phii(i,1)   = 0.0           ! Ignoring topography height here
      enddo
      DO k=1,levs
        do i=1,im
          TEM         = RD * T(i,k) * (1.0 + NU * max(Q(i,k),cons_0))
          DPHI        = (PRSI(i,k) - PRSI(i,k+1)) * TEM
     &                 /(PRSI(i,k) + PRSI(i,k+1))
          phil(i,k)   = phii(i,k) + DPHI
          phii(i,k+1) = phil(i,k) + DPHI
        ENDDO
      ENDDO
!

!     print *,' end of get_prs_v_fd '

      return
      end
