      module idea_das_utils

 
      contains
!
      subroutine compute_dthick_obs(nz, pmb, zkm, dpmb, dzkm)
       implicit none
       integer :: nz
       real, dimension(nz) :: pmb, zkm,  dpmb, dzkm
       integer :: k
       real :: p1
      if(pmb(1).ge.pmb(2)) then
       do k=2, nz-1
        dpmb(k) = .5*(pmb(k-1)-pmb(k+1))
        dzkm(k) = .5*(zkm(k+1)-zkm(k-1))
       enddo
       P1 =pmb(1)  !(1.+exp(dzkm(2)/7.))
       dzkm(1)  = dzkm(2)
       dpmb(1) = (P1-pmb(2))
       dzkm(nz)  = dzkm(nz-1)
       dpmb(nz) = .5*pmb(nz-1)

       else

       do k=2, nz-1
        dpmb(k) = .5*(pmb(k+1)-pmb(k-1))
        dzkm(k) = .5*(zkm(k-1)-zkm(k+1))
      enddo
! reverse-GCM style
         P1 =pmb(nz)  !(1.+exp(dzkm(2)/7.))
         dzkm(nz)  = dzkm(nz-1)
         dpmb(nz) = .5*pmb(nz-1)
         dzkm(1)  = dzkm(2)
         dpmb(1) = .5*pmb(1)

      endif
!print *
!do k=1, nz
!print *, pmb(k), dpmb(k), zkm(k), dzkm(k), k
!enddo
!print *
      end subroutine compute_dthick_obs
!
      Subroutine good_data_1d(Y, nd, ng, Igd, m1, m2)
       implicit none
       integer :: nd, ng, k
       real :: m1, m2
       real :: Y(nd)
       integer :: Igd(nd)
       Igd(1:nd)=-99.
       ng=0
       do k=1, nd
       if (Y(k).ge.m1.and.Y(k).le.m2) then 
        ng=ng+1
        Igd(ng)=k
       endif
      enddo
      end subroutine good_data_1d
!
      SUBROUTINE WHERE_1D(XL, X, XR, N, Ip, Nobs, ix1, ix2)
!
!  For Monotonic XL=>XR
!  XL < X < XR + indices Ip
!  X(ip) = X-interval
!

       Implicit NONE
       integer :: N
       integer :: Ip(N)
       real :: X(N)
       real :: XL, XR                          ! BIN
       integer :: Nobs, ix1, ix2
!
       integer i, j, ik, Ng
       ng = 0
       ip(1:N)=-99
       do ik=1, n
!
       if(x(ik).ge.xL.and.x(ik).le.XR) then
        ng = ng+1
        ip(ng) = ik
       endif
       enddo
      Nobs =Ng
      ix2 = maxval(Ip(1:Nobs))
      ix1 = minval(Ip(1:Nobs))
!
!     print *, ix1, ix2
!     pause '  ix1-ix2 '
!
       end subroutine WHERE_1D
!
       SUBROUTINE WHERE_NUMB( X, Nx, aX1, aX2, N12, ix1, ix2)
         Integer :: Nx
         real    :: X( Nx)          ! montonic array 
         real:: aX1, aX2      ! bounds ax1 < ax2
         Integer :: N12, ix1, ix2 ! # and coordinates in X(nx) 
         Integer :: i
         real :: X1, X2, dx, xmax
         integer, parameter :: iulog=6
!        print *, Nx, aX1, aX2
          N12 =0
            ix1 =1

            x1 = ax1
            x2 = ax2
   
            ix2 = size(X)
            dx = X(2) -X(1)
            xmax = maxval(X)

            if((x2-x1).lt.dx) then

               if (x1.gt.dx) x1 = x1-dx
               if (x2.gt.xmax-dx) x2 = x2-dx
!      write(iulog,*) ' X1:X2 interval needs special treatment ', X1, X2, dX
            endif
            N12 =0

         do i=2, Nx-1

           if (X(i).ge.X1.and.X(i).le.X2)  N12 = N12+1

           if (X(i).ge.X1.and.X(i-1).lt.X1)  ix1 = i
           if (X(i).le.X2.and.X(i+1).gt.X2)  ix2 = i

         enddo
!   Ix1 = Ix1+1
!   Ix2 = Ix2-1
!
            If ((Ix2-Ix1+1).ne.N12) then
              write(iulog,*)  aX1, aX2,' aX1-aX2 '
              write(iulog,*)  X1, X2,  ' X1-X2 .... Ix2-Ix1+1 =/=N12'
              write(iulog,*) Ix2, Ix1, ' ix2-ix1 ', N12,  ' n12 '
            endif
   
            if (Ix1.gt.Nx.or.Ix2.lt.1) then
!
            write(iulog,*) ' error in WHERE_NUMB '
            write(iulog,*)   X1, X2, ' X-coord' , minval(X), maxval(X)
            write(iulog,*)   ' N12 ', N12
            write(iulog,*)    'ix1', ix1, ' ix2 ', ix2
            write(iulog,*)  ' ++++++++++++++++++++++++++++ '
            write(iulog,*)  X
            write(iulog,*)    ' ++++++++++++++++++++++++++++ '
!            call endrun(' das_utils.F90 WHERE_NUMB  .... errors')
            endif
         end Subroutine WHERE_NUMB
!


!
         end module idea_das_utils 
!
