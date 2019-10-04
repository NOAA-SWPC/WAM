!==================================================================
!  This set of subroutines recursively bisect the given triangle 
!  and generate the icos grid.
!
!  In addition to create the latlon coordinates for the grid points, 
!  the subroutines also order the grid points according to their 
!  geographic locations. There are two purposes to do this: 
!
!  a). to increase the cache hit rate and to facilitate the further 
!      work on ordring the grid points using various 2D surface 
!      filling curves; 
!  b). to eliminate the need to search for the six (five) nearest
!      neighbors.
!  
!
!  HISTORY: 
!  Aug. 2006:  Original version,  Ning Wang.
!  Oct. 2007:  Added a new scheme for recursive refinement. N.W.
!  Jun. 2009:  added new bisect_triangle and bisect_r_triangle
!              subroutines to allow combined bisection and 
!              tri-sectiuon refinements. N.W. 
!==================================================================
RECURSIVE SUBROUTINE bisect_triangle(a, b, c, level)

     USE DataStru
     IMPLICIT NONE

     TYPE(GridPoint) :: a, b, c
     INTEGER :: level

     TYPE(GridPoint) :: a_b, a_c, b_c 
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

!              a
!              /\
!             /  \
!            /    \
!           /  s1  \
!          /        \
!     a_b /__________\ a_c 
!        / \        / \
!       /   \  s3  /   \
!      / s2  \    / s4  \
!     /       \  /       \
!   b/_________\/_________\c
!            b_c 

     CALL compute_index(a%seqnum, b%seqnum, c%seqnum, a_b%seqnum, a_c%seqnum, b_c%seqnum, clevel, 1) 

     CALL middle(a%latlon, b%latlon, a_b%latlon)
     CALL middle(a%latlon, c%latlon, a_c%latlon)
     CALL middle(b%latlon, c%latlon, b_c%latlon)

     CALL bisect_triangle(a, a_b, a_c, clevel)
     CALL bisect_triangle(a_b, b, b_c, clevel)
     CALL bisect_r_triangle(a_b, a_c, b_c, clevel)
     CALL bisect_triangle(a_c, b_c, c, clevel)

END SUBROUTINE bisect_triangle

RECURSIVE SUBROUTINE bisect_r_triangle(a, b, c, level)

     USE DataStru
     IMPLICIT NONE

     TYPE(GridPoint) :: a, b, c
     INTEGER :: level

     TYPE(GridPoint) :: a_b, a_c, b_c 
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
!            a_b
!  a ____________________ b
!    \        /\        /
!     \ s1   /  \  s3  /
!      \    / s2 \    /
!       \  /      \  /
!    a_c \/________\/ b_c
!         \        /
!          \  s4  /
!           \    /
!            \  /
!             \/
!             c

     CALL compute_index(a%seqnum, b%seqnum, c%seqnum, a_b%seqnum, a_c%seqnum, b_c%seqnum, clevel, 0) 

     CALL middle(a%latlon, b%latlon, a_b%latlon)
     CALL middle(a%latlon, c%latlon, a_c%latlon)
     CALL middle(b%latlon, c%latlon, b_c%latlon)

     CALL bisect_r_triangle(a, a_b, a_c, clevel)
     CALL bisect_triangle(a_b, a_c, b_c, clevel)
     CALL bisect_r_triangle(a_b, b, b_c, clevel)
     CALL bisect_r_triangle(a_c, b_c, c, clevel)

END SUBROUTINE bisect_r_triangle


! Non-automatic self recursive bisection

RECURSIVE SUBROUTINE bisect_triangle_nasr(a, b, c, level, msec) 
     USE DataStru
     IMPLICIT NONE
                                                                                                                                                              
     TYPE(GridPoint) :: a, b, c
     INTEGER :: level, side, msec(20)

     TYPE(GridPoint) :: a_b, a_c, b_c 
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

!              a
!              /\
!             /  \
!            /    \
!           /  s1  \
!          /        \
!     a_b /__________\ a_c 
!        / \        / \
!       /   \  s3  /   \
!      / s2  \    / s4  \
!     /       \  /       \
!   b/_________\/_________\c
!            b_c 

     CALL compute_index(a%seqnum, b%seqnum, c%seqnum, a_b%seqnum, a_c%seqnum, b_c%seqnum, clevel, 1) 

     CALL middle(a%latlon, b%latlon, a_b%latlon)
     CALL middle(a%latlon, c%latlon, a_c%latlon)
     CALL middle(b%latlon, c%latlon, b_c%latlon)

    IF (msec(clevel + 1) == 2) THEN
       CALL bisect_triangle_nasr(a, a_b, a_c, clevel, msec)
       CALL bisect_triangle_nasr(a_b, b, b_c, clevel, msec)
       CALL bisect_r_triangle_nasr(a_b, a_c, b_c, clevel, msec)
       CALL bisect_triangle_nasr(a_c, b_c, c, clevel, msec)
    ENDIF
    IF (msec(clevel + 1) == 3) THEN
       CALL trisect_triangle_nasr(a, a_b, a_c, clevel, msec)
       CALL trisect_triangle_nasr(a_b, b, b_c, clevel, msec)
       CALL trisect_r_triangle_nasr(a_b, a_c, b_c, clevel, msec)
       CALL trisect_triangle_nasr(a_c, b_c, c, clevel, msec)
    ENDIF

END SUBROUTINE bisect_triangle_nasr


RECURSIVE SUBROUTINE bisect_r_triangle_nasr(a, b, c, level, msec)

     USE DataStru
     IMPLICIT NONE

     TYPE(GridPoint) :: a, b, c
     INTEGER :: level, msec(20)

     TYPE(GridPoint) :: a_b, a_c, b_c 
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
!            a_b
!  a ____________________ b
!    \        /\        /
!     \ s1   /  \  s3  /
!      \    / s2 \    /
!       \  /      \  /
!    a_c \/________\/ b_c
!         \        /
!          \  s4  /
!           \    /
!            \  /
!             \/
!             c

     CALL compute_index(a%seqnum, b%seqnum, c%seqnum, a_b%seqnum, a_c%seqnum, b_c%seqnum, clevel, 0) 

     CALL middle(a%latlon, b%latlon, a_b%latlon)
     CALL middle(a%latlon, c%latlon, a_c%latlon)
     CALL middle(b%latlon, c%latlon, b_c%latlon)

     IF (msec(clevel + 1) == 2) THEN
       CALL bisect_r_triangle_nasr(a, a_b, a_c, clevel, msec)
       CALL bisect_triangle_nasr(a_b, a_c, b_c, clevel, msec)
       CALL bisect_r_triangle_nasr(a_b, b, b_c, clevel, msec)
       CALL bisect_r_triangle_nasr(a_c, b_c, c, clevel, msec)
     ENDIF
     IF (msec(clevel + 1) == 3) THEN
       CALL trisect_r_triangle_nasr(a, a_b, a_c, clevel, msec)
       CALL trisect_triangle_nasr(a_b, a_c, b_c, clevel, msec)
       CALL trisect_r_triangle_nasr(a_b, b, b_c, clevel, msec)
       CALL trisect_r_triangle_nasr(a_c, b_c, c, clevel, msec)
     ENDIF

END SUBROUTINE bisect_r_triangle_nasr


RECURSIVE SUBROUTINE bisect_triangle_new2(a, b, c, a1, a2, a3, level, side)

     USE DataStru
     IMPLICIT NONE

     TYPE(GridPoint) :: a, b, c, a1, a2, a3
     INTEGER :: level, side

     TYPE(GridPoint) :: a_b, a_c, b_c, a_a_b, a_b_b, a_a_c, a_c_c, b_b_c, b_c_c 
     TYPE(GridPoint) :: ph, anc1, anc2, anc3, t1, t2, t3
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

     CALL compute_index(a%seqnum, b%seqnum, c%seqnum, a_b%seqnum, a_c%seqnum, b_c%seqnum, clevel, 1) 

! compute the bisect points a_b, a_c, b_c
     CALL middle(a%latlon, b%latlon, a_b%latlon)
     CALL middle(a%latlon, c%latlon, a_c%latlon)
     CALL middle(b%latlon, c%latlon, b_c%latlon)
     IF (side == 1) THEN
       CALL oneThird(b_c%latlon, a1%latlon, t1%latlon)
       CALL oneThird(a_c%latlon, a2%latlon, t2%latlon)
       CALL middle(a%latlon, b%latlon, t3%latlon)
       CALL average(t1%latlon, t2%latlon, t3%latlon, a_b%latlon)
       a3 = a_b ! save the result and return it to the caller
     ELSE IF (side == 2) THEN
       CALL oneThird(b_c%latlon, a1%latlon, t1%latlon)
       CALL oneThird(a_b%latlon, a2%latlon, t2%latlon)
       CALL middle(a%latlon, c%latlon, t3%latlon)
       CALL average(t1%latlon, t2%latlon, t3%latlon, a_c%latlon)
       a3 = a_c ! save the result and return it to the caller
     ELSE IF (side == 3) THEN
       CALL oneThird(a_c%latlon, a1%latlon, t1%latlon)
       CALL oneThird(a_b%latlon, a2%latlon, t2%latlon)
       CALL middle(b%latlon, c%latlon, t3%latlon)
       CALL average(t1%latlon, t2%latlon, t3%latlon, b_c%latlon)
       a3 = b_c ! save the result and return it to the caller
     ELSE IF (side == 4) THEN
       ! don't need to compute, it is passed in
       a_b = a1
       a_c = a2
       b_c = a3
     ENDIF

! compute some anchor points for the next level
     CALL middle(a%latlon, a_b%latlon, a_a_b%latlon)
     CALL middle(a_b%latlon, b%latlon, a_b_b%latlon)
     CALL middle(a%latlon, a_c%latlon, a_a_c%latlon)
     CALL middle(a_c%latlon, c%latlon, a_c_c%latlon)
     CALL middle(b%latlon, b_c%latlon, b_b_c%latlon)
     CALL middle(b_c%latlon, c%latlon, b_c_c%latlon)

! recursive calls to create next level grid
     CALL bisect_triangle_new2(a, a_b, a_c, b_b_c, b_c_c, anc1, clevel, 3)
     CALL bisect_triangle_new2(a_b, b, b_c, a_a_c, a_c_c, anc2, clevel, 2)
     CALL bisect_triangle_new2(a_c, b_c, c, a_a_b, a_b_b, anc3, clevel, 1)
     CALL bisect_r_triangle_new2(a_b, a_c, b_c, anc1, anc2, anc3, clevel, 4)

END SUBROUTINE bisect_triangle_new2

RECURSIVE SUBROUTINE bisect_r_triangle_new2(a, b, c, a1, a2, a3, level, side)

     USE DataStru
     IMPLICIT NONE

     TYPE(GridPoint) :: a, b, c, a1, a2, a3
     INTEGER :: level, side

     TYPE(GridPoint) :: a_b, a_c, b_c, a_a_b, a_b_b, a_a_c, a_c_c, b_b_c, b_c_c 
     TYPE(GridPoint) :: ph, anc1, anc2, anc3, t1, t2, t3
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

     CALL compute_index(a%seqnum, b%seqnum, c%seqnum, a_b%seqnum, a_c%seqnum, b_c%seqnum, clevel, 0) 

     CALL middle(a%latlon, b%latlon, a_b%latlon)
     CALL middle(a%latlon, c%latlon, a_c%latlon)
     CALL middle(b%latlon, c%latlon, b_c%latlon)
     IF (side == 1) THEN
       CALL oneThird(b_c%latlon, a1%latlon, t1%latlon)
       CALL oneThird(a_c%latlon, a2%latlon, t2%latlon)
       !CALL middle(t1%latlon, t2%latlon, a_b%latlon)
       CALL middle(a%latlon, b%latlon, t3%latlon)
       CALL average(t1%latlon, t2%latlon, t3%latlon, a_b%latlon)
!       CALL gcc(b_c%latlon, a1%latlon, a_c%latlon, a2%latlon, a%latlon, b%latlon, a_b%latlon)
       a3 = a_b
     ELSE IF (side == 2) THEN
       CALL oneThird(b_c%latlon, a1%latlon, t1%latlon)
       CALL oneThird(a_b%latlon, a2%latlon, t2%latlon)
       !CALL middle(t1%latlon, t2%latlon, a_c%latlon)
       CALL middle(a%latlon, c%latlon, t3%latlon)
       CALL average(t1%latlon, t2%latlon, t3%latlon, a_c%latlon)
!       CALL gcc(b_c%latlon, a1%latlon, a_b%latlon, a2%latlon, a%latlon, c%latlon, a_c%latlon)
       a3 = a_c ! save the result and return it to the caller
     ELSE IF (side == 3) THEN
       CALL oneThird(a_c%latlon, a1%latlon, t1%latlon)
       CALL oneThird(a_b%latlon, a2%latlon, t2%latlon)
       !CALL middle(t1%latlon, t2%latlon, b_c%latlon)
       CALL middle(b%latlon, c%latlon, t3%latlon)
       CALL average(t1%latlon, t2%latlon, t3%latlon, b_c%latlon)
!       CALL gcc(a_c%latlon, a1%latlon, a_b%latlon, a2%latlon, b%latlon, c%latlon, b_c%latlon)
       a3 = b_c ! save the result and return it to the caller
     ELSE IF (side == 4) THEN
       ! don't need to compute, it is passed in
       a_b = a1
       a_c = a2
       b_c = a3
     ENDIF

! compute some anchor points for the next level
     CALL middle(a%latlon, a_b%latlon, a_a_b%latlon)
     CALL middle(a_b%latlon, b%latlon, a_b_b%latlon)
     CALL middle(a%latlon, a_c%latlon, a_a_c%latlon)
     CALL middle(a_c%latlon, c%latlon, a_c_c%latlon)
     CALL middle(b%latlon, b_c%latlon, b_b_c%latlon)
     CALL middle(b_c%latlon, c%latlon, b_c_c%latlon)

! recursive calls to create next level grid
     CALL bisect_r_triangle_new2(a, a_b, a_c, b_b_c, b_c_c, anc1, clevel, 3)
     CALL bisect_r_triangle_new2(a_b, b, b_c, a_a_c, a_c_c, anc2, clevel, 2)
     CALL bisect_r_triangle_new2(a_c, b_c, c, a_a_b, a_b_b, anc3, clevel, 1)
     CALL bisect_triangle_new2(a_b, a_c, b_c, anc1, anc2, anc3, clevel, 4)

END SUBROUTINE bisect_r_triangle_new2


RECURSIVE SUBROUTINE bisect_triangle_new(a, b, c, a1, a2, a3, level, side)

     USE DataStru
     IMPLICIT NONE

     TYPE(GridPoint) :: a, b, c, a1, a2, a3
     INTEGER :: level, side

     TYPE(GridPoint) :: a_b, a_c, b_c
     TYPE(GridPoint) :: ph
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

     CALL compute_index(a%seqnum, b%seqnum, c%seqnum, a_b%seqnum, a_c%seqnum, b_c%seqnum, clevel, 1) 

     CALL middle(a%latlon, b%latlon, a_b%latlon)
     CALL middle(a%latlon, c%latlon, a_c%latlon)
     CALL middle(b%latlon, c%latlon, b_c%latlon)
     IF (side == 1) THEN
       CALL middle(a1%latlon, a2%latlon, a_b%latlon)
     ELSE IF (side == 2) THEN
       CALL middle(a1%latlon, a2%latlon, a_c%latlon)
     ELSE IF (side == 3) THEN
       CALL middle(a1%latlon, a2%latlon, b_c%latlon)
     ELSE IF (side == 4) THEN
       CALL middle(a1%latlon, c%latlon, a_b%latlon)
       CALL middle(b%latlon, a2%latlon, a_c%latlon)
       CALL middle(a%latlon, a3%latlon, b_c%latlon)
     ENDIF

     CALL bisect_triangle_new(a, a_b, a_c, a, b_c, ph, clevel, 3)
     CALL bisect_triangle_new(a_b, b, b_c, b, a_c, ph, clevel, 2)
     CALL bisect_r_triangle_new(a_b, a_c, b_c, a, b, c, clevel, 4)
     CALL bisect_triangle_new(a_c, b_c, c, c, a_b, ph, clevel, 1)

END SUBROUTINE bisect_triangle_new

RECURSIVE SUBROUTINE bisect_r_triangle_new(a, b, c, a1, a2, a3, level, side)

     USE DataStru
     IMPLICIT NONE

     TYPE(GridPoint) :: a, b, c, a1, a2, a3
     INTEGER :: level, side

     TYPE(GridPoint) :: a_b, a_c, b_c 
     TYPE(GridPoint) :: ph
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

     CALL compute_index(a%seqnum, b%seqnum, c%seqnum, a_b%seqnum, a_c%seqnum, b_c%seqnum, clevel, 0) 

     CALL middle(a%latlon, b%latlon, a_b%latlon)
     CALL middle(a%latlon, c%latlon, a_c%latlon)
     CALL middle(b%latlon, c%latlon, b_c%latlon)
     IF (side == 1) THEN
       CALL middle(a1%latlon, a2%latlon, a_b%latlon)
     ELSE IF (side == 2) THEN
       CALL middle(a1%latlon, a2%latlon, a_c%latlon)
     ELSE IF (side == 3) THEN
       CALL middle(a1%latlon, a2%latlon, b_c%latlon)
     ELSE IF (side == 4) THEN
       CALL middle(a1%latlon, c%latlon, a_b%latlon)
       CALL middle(b%latlon, a2%latlon, a_c%latlon)
       CALL middle(a%latlon, a3%latlon, b_c%latlon)
     ENDIF

     CALL bisect_r_triangle_new(a, a_b, a_c, a, b_c, ph, clevel, 3)
     CALL bisect_triangle_new(a_b, a_c, b_c, a, b, c, clevel, 4)
     CALL bisect_r_triangle_new(a_b, b, b_c, b, a_c, ph, clevel, 2)
     CALL bisect_r_triangle_new(a_c, b_c, c, c, a_b, ph, clevel, 1)

END SUBROUTINE bisect_r_triangle_new

SUBROUTINE compute_index(a_sn, b_sn, c_sn, a_b_sn, a_c_sn, b_c_sn, level, delta_triangle)   

     USE DataStru
     IMPLICIT NONE
     INTEGER :: a_sn, b_sn, c_sn, a_b_sn, a_c_sn, b_c_sn, level, delta_triangle
     INTEGER :: init_stride, nrtm, stride1, stride2
     REAL :: row

     init_stride = (sl + 1) - row(a_sn)
     nrtm = (row(c_sn) - row(a_sn)) / 2
     stride1 = (2 * init_stride - (nrtm -1)) * nrtm / 2 
     stride2 = (2 * (init_stride - 1) - (nrtm -1)) * nrtm / 2 
     a_c_sn = a_sn + stride1

     IF (delta_triangle == 1) THEN 
       a_b_sn = (a_sn + b_sn) / 2
       b_c_sn = b_sn + stride2
     ELSE
       b_c_sn = (b_sn + c_sn) / 2
       a_b_sn = a_sn + stride2
     END IF

END SUBROUTINE compute_index

REAL FUNCTION row(x)

     USE DataStru

     IMPLICIT NONE

     INTEGER :: x, s, n

     s = init_c_sn - x
     n = INT((sqrt(REAL(1 + 8 * s)) - 1) / 2)
     row = sl - n

END FUNCTION row

SUBROUTINE remove_edge1(start_idx, end_idx)

     USE DataStru
     IMPLICIT NONE
     
     INTEGER start_idx, end_idx
     INTEGER stride, idx, idx2
     
! mark the grid points that need to be removed
     stride = sl + 1
     idx = start_idx
     DO WHILE (idx <= end_idx .AND. stride >= 0) 
       gp(idx)%seqnum = -1
       idx = idx + stride
       stride = stride - 1
     END DO

! skip marked grid points, and ... 
     stride = 0
     idx2 = start_idx
     DO idx  = start_idx, end_idx
       IF (gp(idx)%seqnum == -1) THEN 
         stride = stride + 1
         CYCLE
       END IF
       gp(idx2)%latlon = gp(idx)%latlon 
       gp(idx2)%seqnum = gp(idx)%seqnum - stride 
       idx2 = idx2 + 1
     END DO
  
     end_idx = end_idx - stride

END SUBROUTINE remove_edge1


SUBROUTINE remove_edge2(start_idx, end_idx)

     USE DataStru
     IMPLICIT NONE
     
     INTEGER start_idx, end_idx
     INTEGER stride, idx, idx2
  
! mark the grid points that need to be removed
     DO idx = start_idx, start_idx + sl 
       gp(idx)%seqnum = -1
     END DO

     stride = sl - 1
     idx = start_idx + sl * 2 
     DO WHILE (idx <= end_idx .AND. stride >= 0) 
       gp(idx)%seqnum = -1
       idx = idx + stride
       stride = stride - 1
     END DO

! skip marked grid points, and ... 
     stride = 0
     idx2 = start_idx
     DO idx  = start_idx, end_idx
       IF (gp(idx)%seqnum == -1) THEN 
         stride = stride + 1
         CYCLE
       END IF
       gp(idx2)%latlon = gp(idx)%latlon(:) 
       gp(idx2)%seqnum = gp(idx)%seqnum - stride 
       idx2 = idx2 + 1
     END DO
  
     end_idx = end_idx - stride

END SUBROUTINE remove_edge2

SUBROUTINE combine(triangle1, triangle2)

     USE DataStru
     IMPLICIT NONE

     TYPE(Triangle) :: triangle1, triangle2
     INTEGER idx, idx2, stride, gap, delta

     DO idx = triangle2%start_idx, triangle2%end_idx
       gp(idx)%seqnum = gp(triangle2%end_idx)%seqnum - gp(idx)%seqnum
     END DO

! for triangle one ...
     stride = sl - 3
     gap = 1
     delta = 1
     
     idx = triangle1%start_idx + (2 * sl - 1)
     DO WHILE (idx <= triangle1%end_idx)     
       DO idx2 = idx, idx + stride 
         gp(idx2)%seqnum = gp(idx2)%seqnum + gap
       END DO
       stride = stride - 1
       delta = delta + 1
       gap = gap + delta
       idx  = idx2
     END DO 

! for triangle two ...
     
     stride = 0
     gap = 2 * sl - 1
     delta = sl - 1
     
     idx = triangle2%end_idx
     DO WHILE (idx >= triangle2%start_idx)     
       DO idx2 = idx, idx - stride, -1 
         gp(idx2)%seqnum = gp(idx2)%seqnum + gap
       END DO
       stride = stride + 1
       delta = delta - 1
       gap = gap + delta
       idx  = idx2
     END DO 


END SUBROUTINE combine

