PROGRAM GridInfoGen

!==========================================================
!  This program reads in the icosahedren grid (grid point 
!  sequencial number, lat and lon values), and generates 
!  the necessary information arrays associated with the grid
!  (Vorinoi cells mainly) to allow efficient numerical 
!  integration of the model (finite volume model).  
!
!  Ning Wang, Sep. 2006, Partially adapted from Yuanfu Xie's 
!             GridGen.F90 program.            
!==========================================================

      USE datastru      ,only: perm, inv_perm
      use module_control,only: control,glvl,curve,NumCacheBlocksPerPE
      use read_queue_namelist,only: GetNprocs
      IMPLICIT NONE
 
      INTEGER :: n,i,j,k,l, jmax, nb(6), seq, diam_sz, edge_ln, gp_type
      INTEGER :: diam_idx, pg_idx, r_i, sgp_type
      INTEGER :: grid_point_type
      REAL*8  :: d, d2r, dist, distance
      REAL*8  :: lat, lon
      REAL*8  :: cos_theta,e,pi2, topgrid(12,2),eps
      REAL*8  :: v(3,6)
      LOGICAL :: double_precision, statistics

      ! FIM grid info arrays:
      REAL*8,  ALLOCATABLE :: icos_grid(:,:)     ! grid location (ll)
      REAL,    ALLOCATABLE :: icos_grid4(:,:)    ! grid location (ll), single precision 
      INTEGER, ALLOCATABLE :: icos_nprox(:)      ! number of neighbors(5 or 6)
      INTEGER, ALLOCATABLE :: icos_prox(:,:)     ! neighbor seq nidex
      REAL*8,  ALLOCATABLE :: icos_edge(:,:,:,:) ! end points of edges (ll)
      REAL,    ALLOCATABLE :: icos_edge4(:,:,:,:)! end points of edges (ll), single precision
      REAL*8,  ALLOCATABLE :: xyz(:,:,:)

      ! temp arrays to help move the grid point around
      REAL*8,  ALLOCATABLE :: icos_grid_tmp(:,:)     ! grid location (ll)
      INTEGER, ALLOCATABLE :: icos_nprox_tmp(:)      ! number of neighbors(5 or 6)
      INTEGER, ALLOCATABLE :: icos_prox_tmp(:,:)     ! neighbor seq nidex
      REAL*8,  ALLOCATABLE :: icos_edge_tmp(:,:,:,:) ! end points of edges (ll)
      REAL*8               :: v1(3), v2(3)
      INTEGER              :: isn,i2,i3
      INTEGER              :: nprocs,PE
      INTEGER              :: HaloSize
      CHARACTER(len=80)    :: DecompInfoFile
      integer              :: glvlh
      integer,allocatable  :: Rstart(:),Rend(:),RegionSize(:)
      LOGICAL              :: write_human_readable

      call control ! To read the namelist

      d2r = 4.0*ATAN(1.0)/180.0
      pi2 = 8.0*ATAN(1.0)
      eps = 1.0e-4

      double_precision = .false.
      statistics = .false.

      ! Read in 12 gridpoints with 5 neighbors:
      OPEN(unit=11,file='top_grid',status='old')
      DO i=1,12
        READ(11,*) topgrid(i,1),topgrid(i,2)
        topgrid(i,1:2) = topgrid(i,1:2)*d2r ! To radians
      ENDDO
      CLOSE(11)

      OPEN(10,file='icos_grid_level.dat',status='old')
      READ(10,*) n

      ALLOCATE(icos_grid(2,n))	
      ALLOCATE(icos_prox(6,n))
      ALLOCATE(icos_nprox(n))
      ALLOCATE(icos_edge(6,2,2,n))	
      ALLOCATE(xyz(n,0:6,3))	
      icos_edge = 0.0

      ! Default 6 neighbors:
      icos_nprox = 6

      DO i = 1, n
        READ(10, *) seq, lat, lon
        IF (lon < 0) THEN
          lon = lon + 360.0
        ENDIF
        icos_grid(1,seq) = lat*d2r	! To radians
        icos_grid(2,seq) = lon*d2r	! To radians
      END DO
      CLOSE(10)

      PRINT*,'Number of gridpoints: ',n

      diam_sz = (n - 2) / 10
      edge_ln = SQRT(REAL(diam_sz)) + 0.5  

      ! Assign the neighbors to each grid point

!            s *
!          /\    /\
!       1 /0 \4 /  \
!  s ___ /____\/____\ ___ s
!        \ 0  /\    /\
!        2\  /0 \5 /  \
!    s ___ \/____\/____\ ___ s
!           \ 0  /\    /
!           3\  /6 \  /
!             \/    \/
!               s *

      DO i = 1, n
       CALL special_points(i, nb, diam_sz, edge_ln, sgp_type)
       IF (sgp_type == 1) THEN
         icos_nprox(i) = 5
       END IF 
       IF (sgp_type == 0) THEN
         gp_type = grid_point_type(i, diam_sz, edge_ln)
         IF (gp_type == 0) THEN  ! internal gridpoint
           nb(1) = i - 1
           nb(2) = i - edge_ln
           nb(3) = i - edge_ln + 1
           nb(4) = i + 1
           nb(5) = i + edge_ln
           nb(6) = i + edge_ln - 1
         ELSE IF (gp_type == 1) THEN ! see the fig above 
           nb(1) = i - 1
           pg_idx = (i - 1) / (2 * diam_sz)
           IF (pg_idx >= 1) THEN
             pg_idx = pg_idx - 1
           ELSE
             pg_idx = 4
           END IF
           r_i = MOD((i - 1), 2 * diam_sz) 
           nb(2) = pg_idx * 2 * diam_sz + (r_i - 1) * edge_ln + 2 
           nb(3) = nb(2) + edge_ln
           nb(4) = i + 1
           nb(5) = i + edge_ln
           nb(6) = i + edge_ln - 1
         ELSE IF (gp_type == 2) THEN ! see the fig above 
           nb(1) = i - edge_ln
           pg_idx = (i - 1) / (2 * diam_sz)
           IF (pg_idx >= 1) THEN
             pg_idx = pg_idx - 1
           ELSE
             pg_idx = 4
           END IF
           r_i = MOD((i - 1), 2 * diam_sz) 
           nb(2) = pg_idx * 2 * diam_sz + diam_sz - 2 * edge_ln + r_i + 2  
           nb(3) = nb(2) + edge_ln
           nb(4) = i + edge_ln  
           nb(5) = nb(4) - 1
           nb(6) = i - 1
         ELSE IF (gp_type == 3) THEN ! see the fig above
           nb(1) = i - edge_ln
           pg_idx = (i - 1) / (2 * diam_sz)
           IF (pg_idx >= 1) THEN
             pg_idx = pg_idx - 1
           ELSE
             pg_idx = 4
           END IF
           r_i = MOD((i - 1), 2 * diam_sz) 
           nb(2) = (pg_idx + 1) * 2 * diam_sz - edge_ln + (r_i - (diam_sz + 2 * edge_ln)) / edge_ln + 2     
           nb(3) = nb(2) + 1
           nb(4) = i + edge_ln
           nb(5) = nb(4) - 1 
           nb(6) = i - 1
         ELSE IF (gp_type == 4) THEN ! see the fig above
           nb(1) = i - edge_ln
           nb(2) = nb(1) + 1
           nb(3) = i + 1
           nb(4) = i + edge_ln
           pg_idx = (i - 1) / (2 * diam_sz)
           IF (pg_idx < 4) THEN
             pg_idx = pg_idx + 1
           ELSE
             pg_idx = 0
           END IF
           r_i = MOD((i - 1), 2 * diam_sz) 
           nb(5) = pg_idx * 2 * diam_sz + 2 + ((r_i - 1) / edge_ln)   
           nb(6) = nb(5) - 1
         ELSE IF (gp_type == 5) THEN ! see the fig above
           nb(1) = i - edge_ln
           nb(2) = nb(1) + 1
           nb(3) = nb(2) + edge_ln
           nb(4) = i + edge_ln
           pg_idx = (i - 1) / (2 * diam_sz)
           IF (pg_idx < 4) THEN
             pg_idx = pg_idx + 1
           ELSE
             pg_idx = 0
           END IF
           r_i = MOD((i - 1), 2 * diam_sz) 
           nb(5) = pg_idx * 2 * diam_sz + r_i - diam_sz + 2 * edge_ln   
           nb(6) = nb(5) - edge_ln
         ELSE IF (gp_type == 6) THEN ! see the fig above
           nb(1) = i - 1
           nb(2) = i - edge_ln
           nb(3) = nb(2) + 1
           nb(4) = i + 1
           pg_idx = (i - 1) / (2 * diam_sz)
           IF (pg_idx < 4) THEN
             pg_idx = pg_idx + 1
           ELSE
             pg_idx = 0
           END IF
           r_i = MOD((i - 1), 2 * diam_sz) 
           nb(5) = pg_idx * 2 * diam_sz  + diam_sz + (r_i - 2 * diam_sz + edge_ln + 1) * edge_ln + 1     
           nb(6) = nb(5) - edge_ln
         END IF
       ENDIF 

       xyz(i, 0, 1) = COS(icos_grid(1,i))*COS(icos_grid(2,i))
       xyz(i, 0, 2) = COS(icos_grid(1,i))*SIN(icos_grid(2,i))
       xyz(i, 0, 3) = SIN(icos_grid(1,i))
       IF (nb(6) == -1) THEN
         jmax = 5
       ELSE
         jmax = 6
       ENDIF
       DO j = 1, jmax
         xyz(i, j, 1) = COS(icos_grid(1,nb(j)))*COS(icos_grid(2,nb(j)))
         xyz(i, j, 2) = COS(icos_grid(1,nb(j)))*SIN(icos_grid(2,nb(j)))
         xyz(i, j, 3) = SIN(icos_grid(1,nb(j)))
       END DO
      
       icos_prox(1:6,i) = nb
        

!       DO j = 1, icos_nprox(i) 
!         distance = dist(icos_grid(1, i), icos_grid(2, i), &
!         icos_grid(1, icos_prox(j,i)) , icos_grid(2, icos_prox(j,i)) )
!         if (i  < 30) THEN
!         if (i  == 0) THEN
!           print*, 'i = ', i, 'j = ', icos_prox(j,i) 
!           print*, icos_grid(1, i), icos_grid(2, i),  icos_grid(1, icos_prox(j,i)), icos_grid(2, icos_prox(j,i)) 
!           print*, 'distance = ', distance  
!           print*, ' '
!         ENDIF
!       END DO


       ! Voronoi corners:
       DO j=1,6
         IF (nb(j) .GT. 0) THEN
           ! Do not do the extra point for a 5-neighbor point:
           k = j+1
           IF (k .GT. 6) k = k-6
           IF (nb(k) .LT. 0) k = k+1
           IF (k .GT. 6) k = k-6

           v1 = xyz(i,j,1:3)-xyz(i,0,1:3)
           v2 = xyz(i,k,1:3)-xyz(i,0,1:3)
           CALL cross_product(v(1,j),v1,v2) 
           d = SQRT(DOT_PRODUCT(v(1:3,j),v(1:3,j)))
           ! Normalize the Voronoi points to the unit sphere:
           v(1:3,j) = v(1:3,j)/d
         ENDIF
       ENDDO

        ! Edges connecting these Voronoi points:
       DO j=1,6
         IF (nb(j) .GT. 0) THEN

         ! Do not do the extra point for a 5-neighbor point:
           k = j+1
           IF (k .GT. 6) k = k-6
           IF (nb(k) .LT. 0) k = k+1
           IF (k .GT. 6) k = k-6

        ! Convert to lat/lon:
           icos_edge(j,1,1,i) = ASIN(v(3,j))
           IF (ABS(icos_edge(j,1,1,i)) .GT. pi2/4.0-eps) THEN
             icos_edge(j,1,2,i) = 0.0
           ELSE 
             IF (v(1,j)/COS(icos_edge(j,1,1,i)) .GE. 1.0) THEN
               icos_edge(j,1,2,i) = 0.0
             ELSE IF (v(1,j)/COS(icos_edge(j,1,1,i)) .LE. -1.0) THEN
               icos_edge(j,1,2,i) = pi2/2.0
             ELSE
               icos_edge(j,1,2,i) = ACOS(v(1,j)/COS(icos_edge(j,1,1,i)))
               IF (v(2,j) .LT. 0.0) then
                 icos_edge(j,1,2,i) = pi2-icos_edge(j,1,2,i)
               ENDIF
             ENDIF
           ENDIF

           icos_edge(j,2,1,i) = ASIN(v(3,k))
           IF (ABS(icos_edge(j,2,1,i)) .GT. pi2/4.0-eps) THEN
             icos_edge(j,2,2,i) = 0.0
           ELSE 
             IF (v(1,k)/COS(icos_edge(j,2,1,i)) .GE. 1.0) THEN
               icos_edge(j,2,2,i) = 0.0
             ELSE IF (v(1,k)/COS(icos_edge(j,2,1,i)) .LE. -1.0) THEN
               icos_edge(j,2,2,i) = pi2/2.0
             ELSE
               icos_edge(j,2,2,i)=ACOS(v(1,k)/COS(icos_edge(j,2,1,i)))
               IF (v(2,k) .LT. 0.0) then
                 icos_edge(j,2,2,i) = pi2-icos_edge(j,2,2,i)
               ENDIF
             ENDIF
           ENDIF

         ENDIF
       ENDDO
      ENDDO

      ! Output lat and lon in IJ cordinates for post
      open(28,file="latlonIJ.dat", form="unformatted")
      call WriteGlvlHeader(28,glvl)
      ALLOCATE(icos_grid4(2,n))
      icos_grid4 = icos_grid
      write(28) icos_grid4(1,:),icos_grid4(2,:)
      close(28)
      DEALLOCATE(icos_grid4)
      print*, 'done saving latlonIJ.dat'

      ! move the grid points around to follow the 2D filling curve order
      call GetNprocs(nprocs)
      allocate(Rstart(nprocs),Rend(nprocs))
      CALL mk_perm  (n,curve,nprocs,NumCacheBlocksPerPE,Rstart,Rend) 

      ALLOCATE(icos_grid_tmp(2,n))	
      DO i = 1, n
        icos_grid_tmp(1:2, i) = icos_grid(1:2,perm(i))
      END DO
      DO i = 1, n
        icos_grid(1:2, i) = icos_grid_tmp(1:2,i)
      END DO

      DEALLOCATE(icos_grid_tmp) 

      ALLOCATE(icos_prox_tmp(6,n))
      ALLOCATE(icos_nprox_tmp(n))
      DO i = 1, n
        icos_prox_tmp(1:6, i) = icos_prox(1:6,perm(i))
        icos_nprox_tmp(i) = icos_nprox(perm(i))
        DO j = 1, icos_nprox_tmp(i)
          icos_prox_tmp(j, i) = inv_perm(icos_prox_tmp(j, i))
!         icos_prox_tmp(j, i) =     perm(icos_prox_tmp(j, i))
        END DO
      END DO
      DO i = 1, n
        icos_prox(1:6, i) = icos_prox_tmp(1:6,i)
        icos_nprox(i) = icos_nprox_tmp(i)
      END DO
      DEALLOCATE(icos_prox_tmp)
      DEALLOCATE(icos_nprox_tmp)

      call UnstructuredHaloCalc(6,n,nprocs,icos_prox,Rstart,Rend,HaloSize) !Get halo size
      allocate(RegionSize(nprocs))
      do PE=1,nprocs
        RegionSize(PE) = Rend(PE)-Rstart(PE)+1
      enddo
      deallocate(Rstart,Rend)
      print*,'HaloSize,minval(RegionSize)',HaloSize,minval(RegionSize)
      if(HaloSize > minval(RegionSize)) then
        !Halo size greater than interior size causes problems in GET_PPP_TRANS_INDX
        !and perhaps other places because it was not allowed by SMS and therefore not tested.
        print*,'Error in GridGen: Halo size > interior size', &
                HaloSize,minval(RegionSize)
        print*,'Reduce the number of processor and/or increase ', &
               'the G level and try again'
        stop
      endif

      ALLOCATE(icos_edge_tmp(6,2,2,n))	
      DO i = 1, n
        icos_edge_tmp(1:6, 1:2, 1:2, i) = icos_edge(1:6,1:2, 1:2, perm(i))
      END DO
      DO i = 1, n
        icos_edge(1:6, 1:2, 1:2, i) = icos_edge_tmp(1:6,1:2, 1:2, i)
      END DO
      DEALLOCATE(icos_edge_tmp)

      IF (statistics) THEN
       CALL minmax_values(n, icos_grid, icos_prox, icos_nprox, icos_edge)
      END IF
      DEALLOCATE(xyz)

      ! Write out the grid information:

      ! halo size is used for SMS parallel runs.  
      write(DecompInfoFile,"('DecompInfo_',i0,'.dat')") nprocs
      open (10,file=TRIM(DecompInfoFile))
      write(10,*) HaloSize
      write(10,*) RegionSize
      close(10)
      deallocate(RegionSize)

      write_human_readable = .FALSE.
      IF (write_human_readable) THEN
      !TBH:  write a more easily human-readable file
      OPEN (40,file="icos_grid_info_human.txt",form='formatted')
      WRITE(40,*) 'number of grid cells = ',n
      WRITE(40,*) 'glvl = ',glvl
      WRITE(40,*) 'curve = ',curve
      IF (.not.double_precision) THEN
        ALLOCATE(icos_grid4(2,n))
        ALLOCATE(icos_edge4(6,2,2,n))
        icos_grid4 = icos_grid
        icos_edge4 = icos_edge
      ENDIF
      do i=1,n
        WRITE(40,*) 'Grid cell ',i
        WRITE(40,*) ' Number of sides = ',icos_nprox(i)
        do isn = 1,icos_nprox(i)
          WRITE(40,*) ' Neighbor across edge #',isn,' = ',icos_prox(isn,i)
        enddo
        IF (double_precision) THEN
          WRITE(40,*) ' Latitude  of cell center (radians) = ',icos_grid(1,i)
          WRITE(40,*) ' Longitude of cell center (radians) = ',icos_grid(2,i)
          do isn = 1,icos_nprox(i)
            WRITE(40,*) ' Edge #',isn,':'
            WRITE(40,*) '  Lat-lon of 1st edge endpoint = ', &
                        icos_edge(isn,1,1,i),',',icos_edge(isn,1,2,i)
            WRITE(40,*) '  Lat-lon of 2nd edge endpoint = ', &
                        icos_edge(isn,2,1,i),',',icos_edge(isn,2,2,i)
          enddo
        ELSE
          WRITE(40,*) ' Latitude  of cell center (radians) = ',icos_grid4(1,i)
          WRITE(40,*) ' Longitude of cell center (radians) = ',icos_grid4(2,i)
          do isn = 1,icos_nprox(i)
            WRITE(40,*) ' Edge #',isn,':'
            WRITE(40,*) '  Lat-lon of 1st edge endpoint = ', &
                        icos_edge4(isn,1,1,i),',',icos_edge4(isn,1,2,i)
            WRITE(40,*) '  Lat-lon of 2nd edge endpoint = ', &
                        icos_edge4(isn,2,1,i),',',icos_edge4(isn,2,2,i)
          enddo
        ENDIF
      enddo
      IF (.not.double_precision) THEN
        DEALLOCATE(icos_grid4)
        DEALLOCATE(icos_edge4)
      ENDIF
      CLOSE (40)
      ENDIF   ! write_human_readable

      OPEN (10,file="icos_grid_info_level.dat",form='unformatted')
      call WriteGlvlHeader (10,glvl )
      call WriteCurveHeader(10,curve)

      IF (double_precision) THEN
        WRITE(10) icos_grid(1,1:n),icos_grid(2,1:n)
        do isn = 1,size(icos_prox,1)
          WRITE(10) icos_prox(isn,1:n)
        enddo
        WRITE(10) icos_nprox(1:n)
        do i3 = 1,size(icos_edge,3)
          do i2 = 1,size(icos_edge,2)
            do isn = 1,size(icos_edge,1)
              WRITE(10) icos_edge(isn,i2,i3,1:n)
            enddo
          enddo
        enddo
        do i=1,n
          if(perm(inv_perm(i)) /= i) then
            print*,'inv_perm error in GridInfoGen',i,inv_perm(i),perm(inv_perm(i))
            stop
          endif
        enddo
        WRITE(10) inv_perm
        DEALLOCATE(icos_grid)
        DEALLOCATE(icos_prox)
        DEALLOCATE(icos_nprox)
        DEALLOCATE(icos_edge)
      ELSE 
        ALLOCATE(icos_grid4(2,n))	
        icos_grid4 = icos_grid
        DEALLOCATE(icos_grid)
        ALLOCATE(icos_edge4(6,2,2,n))	
        icos_edge4 = icos_edge
        DEALLOCATE(icos_edge)
        WRITE(10) icos_grid4(1,1:n),icos_grid4(2,1:n)
        do isn = 1,size(icos_prox,1)
          WRITE(10) icos_prox(isn,1:n)
        enddo
        WRITE(10) icos_nprox(1:n)
        do i3 = 1,size(icos_edge4,3)
          do i2 = 1,size(icos_edge4,2)
            do isn = 1,size(icos_edge4,1)
              WRITE(10) icos_edge4(isn,i2,i3,1:n)
            enddo
          enddo
        enddo
        do i=1,n
          if(perm(inv_perm(i)) /= i) then
            print*,'inv_perm error in GridInfoGen',i,inv_perm(i),perm(inv_perm(i))
            stop
          endif
        enddo
        WRITE(10) inv_perm
        DEALLOCATE(icos_grid4)
        DEALLOCATE(icos_prox)
        DEALLOCATE(icos_nprox)
        DEALLOCATE(icos_edge4)
      END IF
     
      CLOSE(10)

END PROGRAM GridInfoGen

! functions and subroutines 

INTEGER FUNCTION grid_point_type(seq, diam_sz, edge_ln)
     IMPLICIT NONE
     INTEGER seq, diam_sz, edge_ln
     INTEGER i

     i = MOD((seq - 1), 2 * diam_sz)

     IF (i < edge_ln ) THEN
       grid_point_type = 1
     ELSE IF (MOD(i, edge_ln) == 0 .AND. i <= diam_sz) THEN
       grid_point_type = 2
     ELSE IF (MOD(i, edge_ln) == 0 .AND. i < 2 * diam_sz) THEN
       grid_point_type = 3
     ELSE IF (MOD(i, edge_ln) == 1 .AND. i < diam_sz) THEN
       grid_point_type = 4
     ELSE IF (MOD(i, edge_ln) == 1 .AND. i < 2 * diam_sz) THEN
       grid_point_type = 5
     ELSE IF (i >= 2 * diam_sz - edge_ln + 1) THEN
       grid_point_type = 6
     ELSE 
       grid_point_type = 0
     END IF 

END FUNCTION grid_point_type


SUBROUTINE special_points(i, nb, diam_sz, edge_ln, type)
       INTEGER i, diam_sz, edge_ln, type
       INTEGER nb(6)

       INTEGER r_i, pg_idx

       type = 0
       r_i = MOD((i - 1), 2 * diam_sz)
       IF (i == 1) THEN ! North pole
         nb(1) = 2
         DO k = 2, 5
           nb(k) = nb(k-1) + 2 * diam_sz
         END DO
         nb(6) = -1
         type = 1
       ELSE IF (i == 10 * diam_sz + 2) THEN ! South pole
         nb(1) = 10 * diam_sz + 1
         DO k = 2, 5
           nb(k) = nb(k-1) - 2 * diam_sz
         END DO
         nb(6) = -1
         type = 1
       ELSE IF (r_i == 1) THEN ! top end point of the parallel graph
         nb(1) = 1    !North pole
         IF ((i - 2 * diam_sz) > 0) THEN 
           nb(2) = i - 2 * diam_sz
         ELSE 
           nb(2) = i + 8 * diam_sz 
         END IF
         nb(3) =  nb(2) + edge_ln 
         nb(4) = i + 1
         nb(5) = i + edge_ln 
         nb(6) = MOD(i + 2 * diam_sz, 10 * diam_sz)
         type = 2
       ELSE IF (r_i == 0) THEN ! bottom end point of the parallel graph
         nb(1) = i - edge_ln 
         IF ((i - 4 * diam_sz) > 0) THEN 
           nb(2) = i - 2 * diam_sz - 1
         ELSE 
           nb(2) = i + 8 * diam_sz - 1
         END IF
         nb(3) = nb(2) + 1 
         nb(4) = 10 * diam_sz + 2    ! South pole
         IF ((i + 2 * diam_sz ) <= 10 * diam_sz + 1) THEN 
           nb(5) = i + 2 * diam_sz  
         ELSE
           nb(5) = i + 2 * diam_sz - 10 * diam_sz 
         END IF
         nb(6) = i - 1
         type = 2
       ELSE IF (r_i == edge_ln) THEN ! 5 neighbors grid point in the first diamond
         nb(1) = i - 1
         pg_idx = i / (2 * diam_sz)
         if (pg_idx >= 1) THEN 
           pg_idx = pg_idx - 1
         ELSE 
           pg_idx = 4 
         END IF
         nb(2) = pg_idx * 2 * diam_sz  + diam_sz - edge_ln + 2
         nb(3) = nb(2) + edge_ln
         nb(4) = i + edge_ln
         nb(5) = nb(4) - 1     
         nb(6) = -1 
         type = 1
       ELSE IF (r_i == diam_sz + edge_ln) THEN ! 5 neighbors grid point in the second diamond
         nb(1) = i - 1
         nb(2) = i - edge_ln
         pg_idx = i / (2 * diam_sz)
         if (pg_idx >= 1) THEN 
           pg_idx = pg_idx - 1
         ELSE 
           pg_idx = 4 
         END IF
         nb(3) = pg_idx * 2 * diam_sz  + 2 * diam_sz - edge_ln + 2
         nb(4) = i + edge_ln
         nb(5) = nb(4) - 1     
         nb(6) = -1 
         type = 1
       ELSE IF (r_i == 2 * diam_sz - edge_ln + 1) THEN   
         nb(1) = i - edge_ln
         nb(2) = nb(1) + 1
         nb(3) = nb(2) + edge_ln
         pg_idx = (i - 1) / (2 * diam_sz)
         IF (pg_idx < 4) THEN
           pg_idx = pg_idx + 1
         ELSE
           pg_idx = 0
         END IF
         r_i = MOD((i - 1), 2 * diam_sz)
         nb(4) = pg_idx * 2 * diam_sz + r_i - diam_sz + 3 * edge_ln
         nb(5) = nb(4) - edge_ln 
         nb(6) = nb(5) - edge_ln
         type = 2
       END IF
END SUBROUTINE special_points


REAL*8 FUNCTION cos_theta(v,w)

IMPLICIT NONE

REAL, INTENT(IN) :: v(3),w(3)

 cos_theta =      DOT_PRODUCT(v,w) / &
             SQRT(DOT_PRODUCT(v,v))/ &
             SQRT(DOT_PRODUCT(w,w))

END FUNCTION cos_theta

SUBROUTINE cross_product(p,v,w)

IMPLICIT NONE

REAL*8, INTENT(IN) :: v(3),w(3)
REAL*8, INTENT(OUT) :: p(3)

p(1) = v(2)*w(3)-v(3)*w(2)
p(2) = v(3)*w(1)-v(1)*w(3)
p(3) = v(1)*w(2)-v(2)*w(1)

END SUBROUTINE cross_product

SUBROUTINE minmax_values(m, icos_grid, icos_prox, icos_nprox, icos_edge)

      IMPLICIT NONE

      REAL*8    :: dist
      INTEGER i, j, k
      INTEGER :: m
      INTEGER :: icos_prox(6,m),icos_nprox(m)
      REAL*8    :: icos_grid(2,m)		! grid location (ll)
      REAL*8    :: icos_edge(6,2,2,m)		! end points of edges (ll)

      REAL*8 max_value, min_value, value(6), ratio, dist_corner(6)
      REAL*8 border_value, max_value_wc, min_value_wc, max_ratio_wc, max_ratio_wt, edge_sum
      REAL*8 avg_grid_dist, avg_grid_dist_sq
      REAL*8 a, b, c, s, area
      REAL*8 lat, lon, lat1, lon1, lat2, lon2, lat3, lon3, min_apr, avg_apr
      INTEGER ::  idx1, idx2, next

      avg_grid_dist = 0.0
      avg_grid_dist_sq = 0.0
      max_value = 0.0
      min_value = 1000000.0
      max_value_wc = 0.0
      min_value_wc = 1000000.0
      max_ratio_wc = 0.0
      max_ratio_wt = 0.0
      border_value = dist(icos_grid(1, 1), icos_grid(2, 1), icos_grid(1, 2), icos_grid(2, 2))
      DO i = 1, m
       DO j = 1, icos_nprox(i)
         value(j) = dist(icos_grid(1, i), icos_grid(2, i), icos_grid(1, icos_prox(j, i)), icos_grid(2, icos_prox(j, i))) 
!         PRINT*, i, icos_prox(j, i), value(j) / border_value
         avg_grid_dist = avg_grid_dist + value(j)
         avg_grid_dist_sq = avg_grid_dist_sq + value(j) * value(j)
         IF (value(j) < min_value) THEN
           min_value = value(j)
         END IF
         IF (value(j) > max_value) THEN
           max_value = value(j)
         END IF
         IF (value(j) < min_value_wc) THEN
           min_value_wc = value(j)
         END IF
         IF (value(j) > max_value_wc) THEN
           max_value_wc = value(j)  
         END IF
       END DO
       DO j = 1, icos_nprox(i)
         IF (j == icos_nprox(i)) THEN
           ratio = max(value(j), value(1)) / min(value(j), value(1))
           idx2 = 1
         ELSE
           ratio = max(value(j), value(j+1)) / min(value(j), value(j+1))
           idx2 = j + 1
         END IF
         IF (max_ratio_wt < ratio) THEN
           max_ratio_wt = ratio 
           idx1 = i
         END IF
       END DO
       IF (max_ratio_wc < max_value_wc / min_value_wc) THEN
         max_ratio_wc = max_value_wc / min_value_wc
       END IF
       max_value_wc = 0.0
       min_value_wc = 1000000.0
      END DO
      PRINT*,'min max distance between neighboring grid points = ', min_value, max_value, &
              max_value / min_value 
      PRINT*,'min max distance ratio within a triangle = ', max_ratio_wt, 'at ', idx1  
      PRINT*,'min max distance ratio within a cell = ', max_ratio_wc  
      PRINT*,'average distance between neighboring grid points = ', avg_grid_dist / (6 * m)
      PRINT*,'average distance square between neighboring grid points =', &
              avg_grid_dist_sq / (6 * m)
     
      max_value = 0.0
      min_value = 1000000.0
      min_apr = 1000000.0
      avg_apr = 0.0
      DO i = 1, m
       DO j = 1, icos_nprox(i)
         value(j) = dist(icos_edge(j, 1, 1, i), icos_edge(j, 1, 2, i), &
           icos_edge(j, 2, 1, i), icos_edge(j, 2, 2, i)) 
         dist_corner(j) = dist(icos_edge(j, 1, 1, i), icos_edge(j, 1, 2, i), &
           icos_grid(1, i), icos_grid(2, i)) 
           IF (value(j) < min_value) THEN
             min_value = value(j)
             lat = icos_grid(1, i)
             lon = icos_grid(2, i)
             lat1 = icos_grid(1, icos_prox(j, i))
             lon1 = icos_grid(2, icos_prox(j, i))
             IF (j < icos_nprox(i)) THEN
               lat2 = icos_grid(1, icos_prox(j+1, i))
               lon2 = icos_grid(2, icos_prox(j+1, i))
             ELSE
               lat2 = icos_grid(1, icos_prox(1, i))
               lon2 = icos_grid(2, icos_prox(1, i))
             END IF
             IF (j < icos_nprox(i) - 1) THEN
               lat3 = icos_grid(1, icos_prox(j+2, i))
               lon3 = icos_grid(2, icos_prox(j+2, i))
             ELSE
               lat3 = icos_grid(1, icos_prox(2, i))
               lon3 = icos_grid(2, icos_prox(2, i))
             END IF
           END IF
           IF (value(j) > max_value) THEN
             max_value = value(j)
           END IF
       END DO
       edge_sum = 0
       DO j = 1, icos_nprox(i)
         edge_sum = edge_sum + value(j)  
         IF (j == icos_nprox(i)) THEN
           next = 1
         ELSE
           next = j+1
         END IF 
         a = value(j) / 6371.220
         b = dist_corner(j) / 6371.220
         c = dist_corner(next) / 6371.220
         s = (a + b + c) / 2.0
         area = tan(s) * tan(s-a) * tan (s-b) * tan(s-c)
         area = atan(area)
         area = area * 4.0
         avg_apr = avg_apr + area / edge_sum
         IF (min_apr > area / edge_sum) THEN
           min_apr = area / edge_sum
         END IF 
       END DO
      END DO
      avg_apr = avg_apr / m
      PRINT*, 'min max distance of the Voronoi cell edges = ', min_value, max_value  
      PRINT*, 'distances of three adjacent neighbors, which yields minimum Voronoi cell edge', &
              dist(lat, lon, lat1, lon1),dist(lat, lon, lat2, lon2),dist(lat, lon, lat3, lon3)
      PRINT*, 'average area / cell perimeter ratio = ', avg_apr, &
              'minmum area / cell perimeter ratio = ', min_apr 

END SUBROUTINE minmax_values


REAL*8 FUNCTION dist(lat1, lon1, lat2, lon2)

IMPLICIT NONE

REAL*8, INTENT(IN) :: lat1, lon1, lat2, lon2
REAL*8             :: x

x = COS(lat1) * COS(lat2) * COS(lon1 - lon2) + SIN(lat1) * SIN(lat2)
if(abs(x) >= 1.0D0) then
  dist = 0.0D0
else
  dist = 6371.220D0 * ACOS(x)
endif

END FUNCTION dist

