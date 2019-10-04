! wgt = 0.2 for now
subroutine smooth(a,mx,my,mz,wgt)
!
! --- 3-dimensional smoothing routine
!
real a(mx,my,mz)
real, allocatable, dimension(:,:) :: a2d
real, allocatable, dimension(:,:,:) :: b
!
allocate(a2d(mx, my))
allocate(b(mx, my, mz))
!
do 6 k=1,mz
!
! --- Smooth in  i  direction
do 4 i=1,mx
ia=2*max0( 1,i-1)-(i-1)
ib=2*min0(mx,i+1)-(i+1)
do 4 j=1,my
4 a2d(i,j)=(1.-wgt-wgt)*a(i,j,k)+wgt*(a(ia,j,k)+a(ib,j,k))
!
! --- Smooth in  j  direction
do 5 j=1,my
ja=2*max0( 1,j-1)-(j-1)
jb=2*min0(my,j+1)-(j+1)
do 5 i=1,mx
5 b(i,j,k)=(1.-wgt-wgt)*a2d(i,j)+wgt*(a2d(i,ja)+a2d(i,jb))
6 continue

a = b

deallocate(a2d,b)
!
return
end

