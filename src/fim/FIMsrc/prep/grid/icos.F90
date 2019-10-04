!----------------------------------------------------
!  This program generates an icosahedron grid using
!  triangle_bisect subroutines. 
!
!  Ning Wang, Aug. 2006, partially adpated from 
!  Yuanfu Xie's icos.F90 files. The main feature of
!  this new version is that the grid points created
!  are in an order on the globe.
!   
!----------------------------------------------------

   PROGRAM Icosahedral_grid

     USE DataStru 
     USE read_queue_namelist,only: ReturnGLVL, ReturnSubdivNum

     IMPLICIT NONE

     ! Grid spec variables
     INTEGER           :: glvl ! The grid level
     INTEGER           :: SubdivNum(20)  ! subdivision specs. 

     ! Position of plot title:
     CHARACTER(len=80) :: datadir
     REAL*8 latlon_d(2), theta, lambda 
     INTEGER i, j, seqnum, nip, buf_sz

     CHARACTER(len=128) :: icos_grid_file  
     CHARACTER(len=1) :: gls

     ! Initial Icosahedron always has 12 grid points 
     TYPE(GridPointWnb) :: top_grid(12)
     TYPE(Triangle) :: top_triangle(20)

! define and read in the name list
     CALL ReturnGLVL(glvl)
     CALL ReturnSubdivNum(SubdivNum) 

! Set up initial 12 icosahedron points:
     CALL set_top_gridpoints(top_grid)

! Rotate the top grid to a desired orientation. It consists of two intrinsic 
! rotations. First rotate theta degrees about the z axis, then lambda degrees
! about the vector v, (-sin(theta), cos(theta), 0).
     theta = 10.0
     lambda = 0.0
     CALL rotate(theta, lambda, top_grid, 12)

! Create the top triangles, with the top_grid values
     CALL set_top_triangle(top_triangle, top_grid) 

     ml = glvl
     sl = 1
     i = 1
     DO WHILE (SubdivNum(i) /= 0 .AND. i <= ml) 
       sl = sl * SubdivNum(i)
       i = i + 1
     ENDDO 

     IF (i /= ml + 1) THEN
       PRINT*, "Namelist variable 'SubdivNum' specification is incomplete."
       STOP
     ENDIF
     SubdivNum(i) = 2

     buf_sz = 20 * (sl + 2) * (sl + 1) / 2 
     ALLOCATE(gp(buf_sz))
     
! num of icos grid points
     nip = 10 * sl * sl + 2

! compute icos grid starting from the top grid
     CALL comp_grid(top_triangle, SubdivNum, 10)

! post process the grid foir each triangle, including remove overlapped grid points, combine triangles
! to form diamond, etc. 
     CALL proc_grid(top_triangle)

! save the results
     icos_grid_file = "icos_grid_level.dat"
     OPEN(10,file=icos_grid_file)
     WRITE(10,*) nip
     ! North pole and South pole
     latlon_d(1) = top_grid(1)%latlon(1)
     latlon_d(2) = top_grid(1)%latlon(2) 
     WRITE(10,*) 1, latlon_d(1),latlon_d(2)
     latlon_d(1) = top_grid(7)%latlon(1)
     latlon_d(2) = top_grid(7)%latlon(2) 
     WRITE(10,*) nip, latlon_d(1),latlon_d(2)
     DO j = 1, 20 
       DO i= top_triangle(j)%start_idx, top_triangle(j)%end_idx
         latlon_d = gp(i)%latlon
         seqnum = gp(i)%seqnum + 1
         WRITE(10,*) seqnum, latlon_d(1),latlon_d(2)
       END DO
     END DO
     CLOSE(10)

     DEALLOCATE(gp)
     PRINT*,'End of the grid point generation.'

   END PROGRAM Icosahedral_grid


SUBROUTINE comp_grid(top_triangle, SubdivNum, type)
     USE DataStru

     IMPLICIT NONE
 
     TYPE(Triangle) :: top_triangle(20)
     INTEGER type, SubdivNum(20)

     INTEGER i
     TYPE(GridPoint) :: ph
     
     init_a_sn = 0
     init_b_sn = sl
     init_c_sn = (sl + 2) * (sl + 1) / 2 - 1
! for each of the 20 triangles, recursively bisect them, and
! order them long the line that is parallel to the side (a, b)
     DO i = 1,20
       PRINT*, i 
       offset = (i - 1) * (init_c_sn + 1) + 1
       top_triangle(i)%vertex(1)%seqnum = init_a_sn 
       top_triangle(i)%vertex(2)%seqnum = init_b_sn 
       top_triangle(i)%vertex(3)%seqnum = init_c_sn
       top_triangle(i)%start_idx = offset + init_a_sn
       top_triangle(i)%end_idx = offset + init_c_sn
       IF (type == 0) THEN
         CALL bisect_triangle(top_triangle(i)%vertex(1), top_triangle(i)%vertex(2), top_triangle(i)%vertex(3), 0)
       ELSE IF (type == 2) THEN
         CALL bisect_triangle_new2(top_triangle(i)%vertex(1), top_triangle(i)%vertex(2), &
                                   top_triangle(i)%vertex(3), ph, ph, ph, 0, 0)
       ELSE IF (type == 10) THEN
         IF (SubdivNum(1) == 2) THEN
           CALL bisect_triangle_nasr(top_triangle(i)%vertex(1), top_triangle(i)%vertex(2), top_triangle(i)%vertex(3), 0, SubdivNum)
         ENDIF
         IF (SubdivNum(1) == 3) THEN
           CALL trisect_triangle_nasr(top_triangle(i)%vertex(1), top_triangle(i)%vertex(2), top_triangle(i)%vertex(3), 0, SubdivNum)
         ENDIF
       ENDIF
     END DO

END SUBROUTINE comp_grid


SUBROUTINE proc_grid(top_triangle)
     USE DataStru
     IMPLICIT NONE

     TYPE(Triangle) :: top_triangle(20)

     INTEGER i, j, j_offset

! for each triangle remove the grid points along the edge that are overlapped with adjacent triangles  
     DO i=1,5
       CALL remove_edge1(top_triangle(i)%start_idx, top_triangle(i)%end_idx)
     END DO
     DO i=11,15
       CALL remove_edge1(top_triangle(i)%start_idx, top_triangle(i)%end_idx)
     END DO
     DO i=6,10
       CALL remove_edge2(top_triangle(i)%start_idx, top_triangle(i)%end_idx)
     END DO
     DO i=16,20
       CALL remove_edge2(top_triangle(i)%start_idx, top_triangle(i)%end_idx)
     END DO

! combine them together to create diamond
     DO i=1,5
       CALL combine(top_triangle(i), top_triangle(i + 15))      
     END DO

     DO i=6,9
       CALL combine(top_triangle(i + 6), top_triangle(i))      
     END DO
     CALL combine(top_triangle(11), top_triangle(10))

! add offset to each diamond
     offset = sl * sl
     DO j = 1, 5 
       j_offset = (j - 1) * offset * 2 + 1 
       DO i= top_triangle(j)%start_idx, top_triangle(j)%end_idx
         gp(i)%seqnum = gp(i)%seqnum + j_offset
       END DO
     END DO 
     DO j = 16, 20 
       j_offset = (j - 16) * offset * 2 + 1 
       DO i= top_triangle(j)%start_idx, top_triangle(j)%end_idx
         gp(i)%seqnum = gp(i)%seqnum + j_offset
       END DO
     END DO 
     DO j = 12, 15 
       j_offset = (j - 12) * offset * 2 + offset + 1 
       DO i= top_triangle(j)%start_idx, top_triangle(j)%end_idx
         gp(i)%seqnum = gp(i)%seqnum + j_offset
       END DO
     END DO 
     j_offset = 4 * offset * 2 + offset + 1 
     DO i= top_triangle(11)%start_idx, top_triangle(11)%end_idx
       gp(i)%seqnum = gp(i)%seqnum + j_offset
     END DO
     DO j = 6, 10 
       j_offset = (j - 6) * offset * 2 + offset + 1 
       DO i= top_triangle(j)%start_idx, top_triangle(j)%end_idx
         gp(i)%seqnum = gp(i)%seqnum + j_offset
       END DO
     END DO 

END SUBROUTINE proc_grid 


