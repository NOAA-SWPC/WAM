      subroutine mymaxmin(a,im,ix,kx,ch)
c    intro of implicit none > >
	implicit none
	integer :: i
	integer :: k
c
      integer ix,kx
      real a(ix,kx)
      character*(*) ch
      integer im
      real fmin,fmax
      do k=1,kx
        fmin = a(1,k)
        fmax = a(1,k)
        do i=1,im
          fmin = min(fmin,a(i,k))
          fmax = max(fmax,a(i,k))
        enddo
        write(0,*)' max=',fmax,' min=',fmin,' at k=',k,' for ',ch
      enddo
      return
      end subroutine mymaxmin
