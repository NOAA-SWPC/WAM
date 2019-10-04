module module_digifilt
implicit none
contains

subroutine digifilt_wts(wts,nwts)
! compute wts for digital filter.
! the type of window is controlled by the
! wts_type namelist parameter.
! if wts_type=1, use Lanczos window.
! if wts_type=2, use Hamming window.
! if wts_type=3, use Dolph window.
! the half-width of the filter window is controlled by the
! tfiltwin namelist parameter (units are seconds).
! on input, wts is an unallocated allocatable array.
! on return, wts(nwts) contains the filter weights
! (normalized so the sum=1).
use module_constants, only : pi
use module_control, only : &
dt,numphr,PhysicsInterval,RadiationInterval,tfiltwin,wts_type

implicit none
real, intent(out), allocatable, dimension(:) :: wts
real, allocatable, dimension(:) :: wttmp
integer, intent(out) :: nwts
integer k,kk,CallPhysics,mm
real sx,hk,sumwts

nwts = numphr*tfiltwin/3600 ! truncated to nearest model time step.
mm = nwts+1 ! index of filter midpoint.
nwts = 2*nwts+1 ! total number of weights

CallPhysics   = max(1,numphr*PhysicsInterval/3600)
! check to see that middle of filter window is on physics time step.
!SMS$SERIAL BEGIN
print *,nwts,'digital filter weights'
print *,'middle of filter window at t = ',mm*dt/3600.,' hours'
if (mod(mm,CallPhysics) .ne. 0) then
 print *,'warning: middle of digital filter window not on a physics time step'
end if
!SMS$SERIAL END

! setup up digital filter weights
allocate(wts(nwts))
allocate(wttmp(0:mm-1))
if (wts_type .eq. 1) then
!SMS$SERIAL BEGIN
   print *,'Lanczos window for digital filter'
!SMS$SERIAL END
   call lanczos(mm-1,wttmp)
else if (wts_type .eq. 2) then
!SMS$SERIAL BEGIN
   print *,'Hamming window for digital filter'
!SMS$SERIAL END
   call hamming(mm-1,wttmp)
else if (wts_type .eq. 3) then
!SMS$SERIAL BEGIN
   print *,'Dolph window for digital filter'
!SMS$SERIAL END
   call dolphwin(mm-1,wttmp)
else
!SMS$SERIAL BEGIN
   print *,'Fatal error: illegal wts_type in digifilt_wts'
!SMS$SERIAL END
   stop
end if
wttmp(mm-1)=0.

sumwts = 0.
do k=1,nwts
   kk = k-mm
   sx = pi*real(kk)/real(mm-1)
   if (kk .ne. 0) then
      hk = sin(sx)/(pi*real(kk)) ! hk --> 1./(mm-1) as sx --> 0
   else
      hk = 1./real(mm-1)
   end if
   if (k .le. mm) then
      wts(k) = wttmp(mm-k)*hk
   else
      wts(k) = wts(2*mm-k)
   end if
   sumwts = sumwts + wts(k)
end do
deallocate(wttmp)
wts = wts/sumwts ! normalize weights so sum=1

!SMS$SERIAL BEGIN
do k=1,nwts
   print *,'filter wt',k,wts(k)
end do
!SMS$SERIAL END

end subroutine digifilt_wts

   SUBROUTINE dolphwin(m, window)

!     calculation of dolph-chebyshev window or, for short,
!     dolph window, using the expression in the reference:
!
!     antoniou, andreas, 1993: digital filters: analysis,
!     design and applications. mcgraw-hill, inc., 689pp.
!
!     the dolph window is optimal in the following sense:
!     for a given main-lobe width, the stop-band attenuation
!     is minimal; for a given stop-band level, the main-lobe
!     width is minimal.
!
!     it is possible to specify either the ripple-ratio r
!     or the stop-band edge thetas.

      IMPLICIT NONE

      ! Arguments
      INTEGER, INTENT(IN)                  ::  m
      REAL, DIMENSION(0:M), INTENT(OUT)    ::  window

      ! local data
      REAL, DIMENSION(0:2*M)               :: t
      REAL, DIMENSION(0:M)                 :: w, time
      REAL    :: pi, thetas, x0, term1, term2, rr, r, db, sum, arg
      INTEGER :: n, nm1, nt, i

      PI = 4*ATAN(1.D0)
      THETAS = 2*PI/M

      N = 2*M+1
      NM1 = N-1
      X0 = 1/COS(THETAS/2)

      TERM1 = (X0 + SQRT(X0**2-1))**(FLOAT(N-1))
      TERM2 = (X0 - SQRT(X0**2-1))**(FLOAT(N-1))
      RR = 0.5*(TERM1+TERM2)
      R = 1/RR
      DB = 20*LOG10(R)
      WRITE(*,'(1X,''DOLPH: M,N='',2I8)')M,N
      WRITE(*,'(1X,''DOLPH: THETAS (STOP-BAND EDGE)='',F10.3)')THETAS
      WRITE(*,'(1X,''DOLPH: R,DB='',2F10.3)')R, DB

      DO NT=0,M
        SUM = RR
        DO I=1,M
          ARG = X0*COS(I*PI/N)
          CALL CHEBY(T,NM1,ARG)
          TERM1 = T(NM1)
          TERM2 = COS(2*NT*PI*I/N)
          SUM = SUM + 2*TERM1*TERM2
        ENDDO
        W(NT) = SUM/N
        TIME(NT) = NT
      ENDDO

!     fill up the array for return
      DO NT=0,M
        WINDOW(NT) = W(NT)
      ENDDO

      RETURN

   END SUBROUTINE dolphwin

   SUBROUTINE cheby(t, n, x)

!     calculate all chebyshev polynomials up to order n
!     for the argument value x.

!     reference: numerical recipes, page 184, recurrence
!         t_n(x) = 2xt_{n-1}(x) - t_{n-2}(x) ,  n>=2.

      IMPLICIT NONE

      ! Arguments
      INTEGER, INTENT(IN)  :: n
      REAL, INTENT(IN)     :: x
      REAL, DIMENSION(0:N) :: t
 
      integer  :: nn

      T(0) = 1
      T(1) = X
      IF(N.LT.2) RETURN
      DO NN=2,N
         T(NN) = 2*X*T(NN-1) - T(NN-2)
      ENDDO

      RETURN

   END SUBROUTINE cheby

   SUBROUTINE LANCZOS(NSTEPS,WW)

   ! define (genaralised) lanczos window function.

      implicit none 

      integer,  parameter                      :: nmax = 1000
      integer,  intent(in)                     :: nsteps
      real   ,  dimension(0:nmax), intent(out) :: ww
      integer  ::  n
      real     :: power, pi, w

      ! (for the usual lanczos window, power = 1 )
      POWER = 1

      PI=4*ATAN(1.)
      DO N=0,NSTEPS
         IF ( N .EQ. 0 ) THEN
            W = 1.0
         ELSE
            W = SIN(N*PI/(NSTEPS+1)) / ( N*PI/(NSTEPS+1))
         ENDIF
         WW(N) = W**POWER
      ENDDO

      RETURN

   END SUBROUTINE lanczos


   SUBROUTINE HAMMING(NSTEPS,WW)

   ! define (genaralised) hamming window function.

      implicit none

      integer, intent(in)           :: nsteps
      real, dimension(0:nsteps)    :: ww
      integer   ::   n
      real      :: alpha, pi, w

      ! (for the usual hamming window, alpha=0.54,
      !      for the hann window, alpha=0.50).
      ALPHA=0.54

      PI=4*ATAN(1.)
      DO N=0,NSTEPS
         IF ( N .EQ. 0 ) THEN
            W = 1.0
         ELSE
            W = ALPHA + (1-ALPHA)*COS(N*PI/(NSTEPS))
         ENDIF
         WW(N) = W
      ENDDO

      RETURN

   END SUBROUTINE hamming

end module module_digifilt
