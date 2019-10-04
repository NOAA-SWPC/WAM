program readGLVL
  implicit none
  integer,parameter :: npp=6
  integer,parameter :: nd=2	! number of directions (x,y)
  integer,parameter :: glvl=6
  integer,parameter :: nip=10*(2**glvl)**2+2
  real lat(nip),lon(nip)		! lat and lon in radians
  integer nprox      (nip)      ! Holds number of proximity points
  integer proxs      (npp,nip)  ! Holds index of proximity sides
  integer prox       (npp,nip)  ! Holds index of proximity points
  real area(nip)	      	   ! the area of cell polygon (m**2)
  real cs(4,npp,nip),sn(4,npp,nip)
  real sidevec_c(nd,npp,nip) ! side vectors projected from center
  real sidevec_e(nd,npp,nip) ! side vectors projected from edge
  real sideln   (npp,nip)	   ! the length of side vectors (m)
  real rprox_ln (npp,nip)    ! reciprocal of distance cell cent to prox pts
  integer isn,ipn

  open(unit=28,file="glvl.dat",form="unformatted")
  read(28)lat,lon,nprox,proxs,prox,area,cs,sn, &
          sidevec_c,sidevec_e,sideln,rprox_ln
  close(28)
  write(76,100) (ipn,(prox(isn,ipn),isn=1,npp),ipn=1,nip)
100 format(7i10)

end program readGLVL
