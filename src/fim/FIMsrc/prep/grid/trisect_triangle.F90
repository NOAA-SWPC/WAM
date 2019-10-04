!=======================================================================
!  This set of subroutines trisect the given triangle and recursively 
!  calls a appropriate multisection routines to generate the icos grid.
!
!  HISTORY: 
!  Jun. 2009:  Original version, Ning Wang. 
!              Created to allow combined multi-section refinements.  
!=======================================================================

! Non-automatic self recursive trisection

RECURSIVE SUBROUTINE trisect_triangle_nasr(a, b, c, level, msec) 
     USE DataStru
     IMPLICIT NONE
                                                                                                                                                              
     TYPE(GridPoint) :: a, b, c
     INTEGER :: level, side, msec(20)

     TYPE(GridPoint) :: a_b_1, a_b_2, a_c_1, a_c_2, b_c_1, b_c_2, o 
     TYPE(GridPoint) :: t1, t2, t3 
     INTEGER :: clevel

     clevel = level + 1
    
! store the vertex grid point in an array according to their seq num.
     gp(a%seqnum + offset)%latlon = a%latlon(:)
     gp(b%seqnum + offset)%latlon = b%latlon(:)
     gp(c%seqnum + offset)%latlon = c%latlon(:)
     gp(a%seqnum + offset)%seqnum = a%seqnum
     gp(b%seqnum + offset)%seqnum = b%seqnum
     gp(c%seqnum + offset)%seqnum = c%seqnum

! check to see if recursive call bottoms out.
     IF (clevel > ml) THEN
       RETURN
     END IF 

!                   a
!                   /\
!                  /  \
!                 /    \
!                /      \
!               /        \
!        a_b_1 /__________\ a_c_1 
!             / \        / \
!            /   \      /   \
!           /     \    /     \
!          /       \o /       \
!   a_b_2 /_________\/_________\ a_c_2
!        / \        /\        / \
!       /   \      /  \      /   \
!      /     \    /    \    /     \
!     /       \  /      \  /       \
!   b/_________\/________\/_________\ c
!            b_c_1     b_c_2
!

     CALL compute_index_tri(a%seqnum, b%seqnum, c%seqnum, &
          a_b_1%seqnum, a_b_2%seqnum, a_c_1%seqnum, a_c_2%seqnum,&
          b_c_1%seqnum,  b_c_2%seqnum, o%seqnum, clevel, 1) 
!print*, a%seqnum , b%seqnum, c%seqnum,a_b_1%seqnum, a_b_2%seqnum, a_c_1%seqnum, a_c_2%seqnum, b_c_1%seqnum,  b_c_2%seqnum, o%seqnum, clevel

     CALL trisect(a%latlon, b%latlon, a_b_1%latlon, a_b_2%latlon)
     CALL trisect(a%latlon, c%latlon, a_c_1%latlon, a_c_2%latlon)
     CALL trisect(b%latlon, c%latlon, b_c_1%latlon, b_c_2%latlon)

     CALL middle(a_b_1%latlon, b_c_2%latlon, t1%latlon)
     CALL middle(a_c_1%latlon, b_c_1%latlon, t2%latlon)
     CALL middle(a_b_2%latlon, a_c_2%latlon, t3%latlon)
     o%latlon = (t1%latlon + t2%latlon + t3%latlon) / 3.0

     IF (msec(clevel + 1) == 2) THEN
       CALL bisect_triangle_nasr(a, a_b_1, a_c_1, clevel, msec)
       CALL bisect_triangle_nasr(a_b_1, a_b_2, o, clevel, msec)
       CALL bisect_r_triangle_nasr(a_b_1, a_c_1, o, clevel, msec)
       CALL bisect_triangle_nasr(a_c_1, o, a_c_2, clevel, msec)
       CALL bisect_triangle_nasr(a_b_2, b, b_c_1, clevel, msec)
       CALL bisect_r_triangle_nasr(a_b_2, o, b_c_1, clevel, msec)
       CALL bisect_triangle_nasr(o, b_c_1, b_c_2, clevel, msec)
       CALL bisect_r_triangle_nasr(o, a_c_2, b_c_2, clevel, msec)
       CALL bisect_triangle_nasr(a_c_2, b_c_2, c, clevel, msec)
     ENDIF
     IF (msec(clevel + 1) == 3) THEN
       CALL trisect_triangle_nasr(a, a_b_1, a_c_1, clevel, msec)
       CALL trisect_triangle_nasr(a_b_1, a_b_2, o, clevel, msec)
       CALL trisect_r_triangle_nasr(a_b_1, a_c_1, o, clevel, msec)
       CALL trisect_triangle_nasr(a_c_1, o, a_c_2, clevel, msec)
       CALL trisect_triangle_nasr(a_b_2, b, b_c_1, clevel, msec)
       CALL trisect_r_triangle_nasr(a_b_2, o, b_c_1, clevel, msec)
       CALL trisect_triangle_nasr(o, b_c_1, b_c_2, clevel, msec)
       CALL trisect_r_triangle_nasr(o, a_c_2, b_c_2, clevel, msec)
       CALL trisect_triangle_nasr(a_c_2, b_c_2, c, clevel, msec)
     ENDIF

END SUBROUTINE trisect_triangle_nasr


RECURSIVE SUBROUTINE trisect_r_triangle_nasr(a, b, c, level, msec)

     USE DataStru
     IMPLICIT NONE

     TYPE(GridPoint) :: a, b, c
     INTEGER :: level, msec(20)

     TYPE(GridPoint) :: a_b_1, a_b_2, a_c_1, a_c_2, b_c_1, b_c_2, o 
     TYPE(GridPoint) :: t1, t2, t3 
     INTEGER :: clevel

     gp(a%seqnum + offset)%latlon = a%latlon(:)
     gp(b%seqnum + offset)%latlon = b%latlon(:)
     gp(c%seqnum + offset)%latlon = c%latlon(:)

     gp(a%seqnum + offset)%seqnum = a%seqnum
     gp(b%seqnum + offset)%seqnum = b%seqnum
     gp(c%seqnum + offset)%seqnum = c%seqnum

     clevel = level + 1

! check to see if recursive call bottoms out.
     IF (clevel > ml) THEN
       RETURN
     END IF 
!            a_b_1    a_b_2
!  a ______________________________ b
!    \        /\        /\        /
!     \      /  \      /  \      /
!      \    /    \    /    \    /
!       \  /      \o /      \  /
!  a_c_1 \/________\/________\/ b_c_1
!         \        /\        /
!          \      /  \      /
!           \    /    \    /
!            \  /      \  /
!       a_c_2 \/________\/ b_c_2
!              \        /
!               \      /
!                \    /
!                 \  /
!                  \/
!                  c

     CALL compute_index_tri(a%seqnum, b%seqnum, c%seqnum, & 
          a_b_1%seqnum, a_b_2%seqnum, a_c_1%seqnum,a_c_2%seqnum, &
          b_c_1%seqnum, b_c_2%seqnum, o%seqnum,  clevel, 0) 

     CALL trisect(a%latlon, b%latlon, a_b_1%latlon, a_b_2%latlon)
     CALL trisect(a%latlon, c%latlon, a_c_1%latlon, a_c_2%latlon)
     CALL trisect(b%latlon, c%latlon, b_c_1%latlon, b_c_2%latlon)

     CALL middle(a_b_1%latlon, b_c_2%latlon, t1%latlon)
     CALL middle(a_c_1%latlon, b_c_1%latlon, t2%latlon)
     CALL middle(a_b_2%latlon, a_c_2%latlon, t3%latlon)
     o%latlon = (t1%latlon + t2%latlon + t3%latlon) / 3.0
     IF (msec(clevel + 1) == 2) THEN
       CALL bisect_r_triangle_nasr(a, a_b_1, a_c_1, clevel, msec)
       CALL bisect_triangle_nasr(a_b_1, a_c_1, o, clevel, msec)
       CALL bisect_r_triangle_nasr(a_b_1, a_b_2, o, clevel, msec)
       CALL bisect_triangle_nasr(a_b_2, o, b_c_1, clevel, msec)
       CALL bisect_r_triangle_nasr(a_b_2, b, b_c_1, clevel, msec)
       CALL bisect_r_triangle_nasr(a_c_1, o, a_c_2, clevel, msec)
       CALL bisect_triangle_nasr(o, a_c_2, b_c_2, clevel, msec)
       CALL bisect_r_triangle_nasr(o, b_c_1, b_c_2, clevel, msec)
       CALL bisect_r_triangle_nasr(a_c_2, b_c_2, c, clevel, msec)
     ENDIF
     IF (msec(clevel + 1) == 3) THEN
       CALL trisect_r_triangle_nasr(a, a_b_1, a_c_1, clevel, msec)
       CALL trisect_triangle_nasr(a_b_1, a_c_1, o, clevel, msec)
       CALL trisect_r_triangle_nasr(a_b_1, a_b_2, o, clevel, msec)
       CALL trisect_triangle_nasr(a_b_2, o, b_c_1, clevel, msec)
       CALL trisect_r_triangle_nasr(a_b_2, b, b_c_1, clevel, msec)
       CALL trisect_r_triangle_nasr(a_c_1, o, a_c_2, clevel, msec)
       CALL trisect_triangle_nasr(o, a_c_2, b_c_2, clevel, msec)
       CALL trisect_r_triangle_nasr(o, b_c_1, b_c_2, clevel, msec)
       CALL trisect_r_triangle_nasr(a_c_2, b_c_2, c, clevel, msec)
     ENDIF

END SUBROUTINE trisect_r_triangle_nasr


SUBROUTINE compute_index_tri(a_sn, b_sn, c_sn, a_b_1_sn, a_b_2_sn, &
           a_c_1_sn, a_c_2_sn, b_c_1_sn, b_c_2_sn, o_sn, level, delta_triangle)   

     USE DataStru

     IMPLICIT NONE

     INTEGER :: a_sn, b_sn, c_sn, a_b_1_sn, a_b_2_sn,  a_c_1_sn, a_c_2_sn  
     INTEGER :: b_c_1_sn, b_c_2_sn, o_sn, level, delta_triangle
     INTEGER :: init_stride, nrtt, nrt2t, stride1_1, stride1_2, stride2_1, stride2_2
     REAL :: row

     init_stride = (sl + 1) - row(a_sn)
     nrtt = (row(c_sn) - row(a_sn)) / 3
     nrt2t = 2 * nrtt 
     stride1_1 = (2 * init_stride - (nrtt -1)) * nrtt / 2 
     stride1_2 = (2 * init_stride - (nrt2t -1)) * nrt2t / 2 
     stride2_1 = (2 * (init_stride - 1) - (nrtt -1)) * nrtt / 2 
     stride2_2 = (2 * (init_stride - 1) - (nrt2t -1)) * nrt2t / 2 
     a_c_1_sn = a_sn + stride1_1
     a_c_2_sn = a_sn + stride1_2

     IF (delta_triangle == 1) THEN 
       a_b_1_sn = a_sn + (b_sn - a_sn) / 3
       a_b_2_sn = a_sn + 2 * (b_sn - a_sn) / 3
       b_c_1_sn = b_sn + stride2_1
       b_c_2_sn = b_sn + stride2_2
       o_sn = (a_c_1_sn + b_c_1_sn) / 2
     ELSE
       b_c_1_sn = b_sn + (c_sn - b_sn) / 3
       b_c_2_sn = b_sn + 2 * (c_sn - b_sn) / 3
       a_b_1_sn = a_sn + stride2_1
       a_b_2_sn = a_sn + stride2_2
       o_sn = (a_b_2_sn + a_c_2_sn) / 2
     END IF

END SUBROUTINE compute_index_tri

