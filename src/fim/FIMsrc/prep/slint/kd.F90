!==========================================================================
! Module for the k-d tree
!
!  The algorithm is based on Bentley's optimal k-d tree, plus 
!  some modifications of mine. The priority queue was omitted
!  since the number of nearest neighbors to be searched is 
!  always relatively small (3-11). 
!
!  Ning Wang, Nov. 2006:  Initial implementation.
!
!  Ning Wang, Jul. 2011:  Created a new version of the k-d tree, which
!    includes following modifications and enhancements.
!    a). wrapped kd_datastru.F90 and kd.F90 into one module;
!    b). simplified interface to the recursive subroutine SearchRec;
!    c). added a new capability to the k-d tree, such that the search
!        space can be bounded with a pair of hyper-planes for each
!        query (search).
!==========================================================================
MODULE kd
     TYPE Node
       LOGICAL :: bucket   ! true if it is a bucket node
!       INTEGER :: npoints  ! number of points of the node 
       INTEGER :: discrim  ! discriminant dimension  
       REAL :: cutval      ! the cut value for the space
       INTEGER :: lopt     ! lower point  
       INTEGER :: hipt     ! high point  
       TYPE(Node), POINTER :: loson  
       TYPE(Node), POINTER :: hison 
     END TYPE Node

     INTEGER n, d, gc, bottom_level, dir, num_k
     REAL, ALLOCATABLE :: hp1(:), hp2(:), tc1(:), tc2(:), qry(:)

     TYPE(Node), TARGET, ALLOCATABLE :: nodes(:)
     INTEGER, ALLOCATABLE :: perm(:)
     INTEGER, TARGET, ALLOCATABLE :: nni(:)
     REAL, ALLOCATABLE :: nnd(:)
     REAL, ALLOCATABLE :: q2cut(:)
     REAL, POINTER:: ppoints(:,:)
     TYPE(Node), POINTER :: root

     REAL :: curDistSq, partsum
CONTAINS

SUBROUTINE BuildTree(num, dim, points)

    IMPLICIT NONE

    INTEGER :: num, dim
    REAL, TARGET :: points(dim, num)
    INTEGER :: i
        
    IF (ALLOCATED(nodes)) THEN
      DEALLOCATE(nodes)
    END IF
    IF (ALLOCATED(perm)) THEN
      DEALLOCATE(perm)
    END IF
    IF (ALLOCATED(q2cut).AND.ALLOCATED(hp1).AND.ALLOCATED(hp2)) THEN
      DEALLOCATE(q2cut, hp1, hp2)
    END IF
    IF (ALLOCATED(qry).AND.ALLOCATED(tc1).AND.ALLOCATED(tc2)) THEN
      DEALLOCATE(qry, tc1, tc2)
    END IF

    ALLOCATE(nodes(2 * num))
    ALLOCATE(perm(num))
    ALLOCATE(q2cut(dim), hp1(dim), hp2(dim))
    ALLOCATE(qry(dim), tc1(dim), tc2(dim))

    DO i = 1, num
      perm(i) = i
    END DO

    n = num
    d = dim
    num_k = 1
    dir = 0
    gc = 1

    ppoints => points
    bottom_level = log(REAL(n)) / log(2.0) + 1
    root => BuildTreeRec(points, 1, n, 1)

END SUBROUTINE BuildTree

SUBROUTINE DeleteTree()

    DEALLOCATE(nodes)
    DEALLOCATE(perm)
    DEALLOCATE(q2cut)
    DEALLOCATE(hp1, hp2, tc1, tc2, qry)
    DEALLOCATE(ppoints)

END SUBROUTINE DeleteTree

RECURSIVE FUNCTION BuildTreeRec(points, l, u, level) RESULT(theNode)

    IMPLICIT NONE

    REAL points(d, n)
    INTEGER :: l, u, level
    TYPE(Node), POINTER :: theNode

    INTEGER :: m

    theNode => NewNode() ! get a new node

    IF (level >= bottom_level) THEN
      theNode%bucket = .true.
      theNode%lopt = l
      theNode%hipt = u
    ELSE 
      theNode%bucket = .false.
      m = (l + u ) / 2 
      theNode%discrim = dir2cut(l, u, points)
      CALL partition(points,l, u, m, theNode%discrim) 
      theNode%cutval = points(theNode%discrim, perm(m))
      theNode%loson => BuildTreeRec(points, l, m, level + 1)
      theNode%hison => BuildTreeRec(points, m+1, u, level + 1) 
    END IF 

END FUNCTION BuildTreeRec

SUBROUTINE Set_k(k)
    IMPLICIT NONE

    INTEGER::k

    num_k = k

    IF (ALLOCATED(nni)) THEN
      DEALLOCATE(nni)
    END IF
    IF (ALLOCATED(nnd)) THEN
      DEALLOCATE(nnd)
    END IF

    ALLOCATE(nni(k))
    ALLOCATE(nnd(k))

END SUBROUTINE Set_k

SUBROUTINE Set_qry(q)
    IMPLICIT NONE

    REAL::q(d)

    qry = q

END SUBROUTINE Set_qry

SUBROUTINE Set_hps(hp_1, hp_2)
    IMPLICIT NONE

    REAL::hp_1(d), hp_2(d)

    hp1 = hp_1
    hp2 = hp_2

END SUBROUTINE Set_hps

SUBROUTINE Set_tcs(hp_1, hp_2)
    IMPLICIT NONE

    REAL::hp_1(d), hp_2(d)
    INTEGER i

    DO i = 1, d
      tc1(i) =  REAL(sgn(hp_1(i)))
      tc2(i) =  REAL(sgn(hp_2(i)))
    END DO

END SUBROUTINE Set_tcs

FUNCTION sgn(x)
    IMPLICIT NONE
    REAL x
    INTEGER sgn
    IF (x .GE. 0.0) THEN
        sgn = 1
    ELSE
        sgn = -1
    END IF
END FUNCTION sgn

SUBROUTINE Search(query) 
  
    IMPLICIT NONE

    REAL :: query(d)

    INTEGER :: i

    partsum = 0
    DO i = 1, d
      q2cut(i) = 0.0
    END DO
    
    DO i = 1, num_k
      nnd(i) = 10.0 * 10.0
      nni(i) = 0
    END DO
    curDistSq = 10.0 * 10.0 

    CALL SearchRec(query, root)

END SUBROUTINE Search
 
RECURSIVE SUBROUTINE SearchRec(query, tnode) 

    IMPLICIT NONE

    REAL :: query(d), cur_q2c, cur_ps, cur_tc1, cur_tc2
    TYPE(Node) :: tnode

    INTEGER :: i
    REAL :: distSq, d2o
    IF (tnode%bucket) THEN
      DO i = tnode%lopt, tnode%hipt
        IF (in_ss(ppoints(1:3, perm(i)), ppoints(1:3, perm(i)))) THEN 
        distSq = inp(query, ppoints(1:3, perm(i))) 
        IF (distSq < curDistSq) THEN
          curDistSq = distSq
          CALL insert(perm(i))
        END IF
        END IF
      END DO   
    ELSE 
      d2o = query(tnode%discrim) - tnode%cutval
      IF (d2o < 0.0) THEN
        cur_ps = partsum
        cur_q2c = q2cut(tnode%discrim)
        cur_tc1 = tc1(tnode%discrim)
        cur_tc2 = tc2(tnode%discrim)
        CALL SearchRec(query, tnode%loson)
        partsum = cur_ps + d2o * d2o - q2cut(tnode%discrim) 
        q2cut(tnode%discrim) = d2o * d2o 
        tc1(tnode%discrim) = max(tnode%cutval, tc1(tnode%discrim))
        tc2(tnode%discrim) = max(tnode%cutval, tc2(tnode%discrim))
        IF (partsum < curDistSq .AND. in_ss(tc1, tc2) ) THEN
          CALL SearchRec(query, tnode%hison) 
        END IF
        q2cut(tnode%discrim) = cur_q2c
        tc1(tnode%discrim) = cur_tc1
        tc2(tnode%discrim) = cur_tc2
      ELSE
        cur_ps = partsum
        cur_q2c = q2cut(tnode%discrim)
        cur_tc1 = tc1(tnode%discrim)
        cur_tc2 = tc2(tnode%discrim)
        CALL SearchRec(query, tnode%hison)
        partsum = cur_ps + d2o * d2o - q2cut(tnode%discrim) 
        q2cut(tnode%discrim) = d2o * d2o 
        tc1(tnode%discrim) = min(tnode%cutval, tc1(tnode%discrim))
        tc2(tnode%discrim) = min(tnode%cutval, tc2(tnode%discrim))
        IF (partsum < curDistSq .AND. in_ss(tc1, tc2)) THEN
          CALL SearchRec(query, tnode%loson) 
        END IF
        q2cut(tnode%discrim) = cur_q2c
        tc1(tnode%discrim) = cur_tc1
        tc2(tnode%discrim) = cur_tc2
      END IF
    END IF

 END SUBROUTINE SearchRec

 FUNCTION in_ss(tc1, tc2)
     IMPLICIT NONE

     LOGICAL in_ss
     REAL tc1(d), tc2(d)

     REAL ip, tc(d)
     INTEGER i

     tc = tc1 - qry
     ip = 0.0
     DO i = 1, d
       ip = ip + tc(i) * hp1(i)
     END DO
     IF (ip < 0) THEN
       in_ss = .false.
       RETURN
     ENDIF

     tc = tc2 - qry
     ip = 0.0
     DO i = 1, d
       ip = ip + tc(i) * hp2(i)
     END DO
     IF (ip < 0) THEN
       in_ss = .false.
       RETURN
     ENDIF

     in_ss = .true.
 END FUNCTION in_ss

! Subroutines and functions that are helps creatation and searching k-d tree.
! get a new node
FUNCTION NewNode() 
    IMPLICIT NONE

    TYPE(Node), POINTER :: NewNode

    NewNode => nodes(gc)
    gc = gc + 1
   
END FUNCTION NewNode

! partition the points along dir 'discrim' into lower and upper parts 
SUBROUTINE partition(points, l, u, m, discrim)
    IMPLICIT NONE

    REAL :: points(d, n)
    INTEGER :: l, u, m, discrim

    REAL :: v
    INTEGER :: i, j, t, r, lo

    r = u
    lo = l

    DO WHILE ( r > lo)
      v = points(discrim, perm(r))
      i = lo
      j = r - 1
      DO WHILE (.true.)
        DO WHILE (points(discrim,perm(i)) < v) 
          i = i + 1
        END DO
        DO WHILE (points(discrim, perm(j)) >= v .AND. j > lo)
          j = j - 1
        END DO
        IF (i >= j) EXIT
        t = perm(i)
        perm(i) = perm(j)
        perm(j) = t
      END DO
      t = perm(i)
      perm(i) = perm(r)
      perm(r) = t
      IF (i >= m) r = i - 1;
      IF (i <= m) lo = i + 1
    END DO
      
END SUBROUTINE partition

! function returns the direction to divide
FUNCTION dir2cut(l, u, points)
    IMPLICIT NONE

    REAL :: points(d, n)
    INTEGER :: l, u

    INTEGER :: dir2cut

    dir = dir + 1
    IF (dir > d) THEN
      dir = 1
    END IF
    dir2cut = dir
   
END FUNCTION dir2cut

! function to compute the inner product of p1-p2
FUNCTION inp(p1, p2) 
    IMPLICIT NONE
        
    REAL :: p1(d), p2(d)
    REAL :: inp

    REAL sum, dif
    INTEGER i

    sum = 0
    DO i = 1, d
      dif = p1(i) - p2(i)
      sum = sum + dif * dif
    END DO 

    inp = sum

END FUNCTION inp

! subroutine to insert the current nn  
SUBROUTINE insert(pt_idx)
    IMPLICIT NONE

    INTEGER :: pt_idx
 
    INTEGER :: i, j

    DO i = 1, num_k	
      IF (curDistSq < nnd(i)) THEN
        DO j = num_k, i + 1, -1
          nni(j) = nni(j - 1)
          nnd(j) = nnd(j - 1)
        END DO
        nni(i) = pt_idx
        nnd(i) = curDistSq
        EXIT            
      END IF
    END DO
      
    curDistSq = nnd(num_k)
          
END SUBROUTINE insert
          

! function to return the index array
SUBROUTINE result()
    IMPLICIT NONE
    INTEGER :: i
      
    DO i = 1, num_k
      PRINT*, nni(i)
      PRINT*, ppoints(:, nni(i)), nnd(i)
    END DO 

END SUBROUTINE result

SUBROUTINE init_kd_tree(llpoints, n, k)
    IMPLICIT NONE
    integer,intent(in) :: n,k
    REAL   ,intent(in) :: llpoints(n,2) 
    REAL, ALLOCATABLE, SAVE :: points(:,:) 
    INTEGER :: i, seq, dim

    dim = 3
    IF (ALLOCATED(points)) THEN
      DEALLOCATE(points)
    ENDIF
    ALLOCATE(points(dim,n))
    CALL lls2xyzs(llpoints, points, n)
    CALL BuildTree(n, dim, points)
    CALL Set_k(k)
END SUBROUTINE init_kd_tree

SUBROUTINE close_kd_tree()

    CALL DeleteTree()

END SUBROUTINE close_kd_tree
         

SUBROUTINE knn_search(ll, nn, min_dist, hp1, hp2)
    IMPLICIT NONE
    
    REAL ll(2) 
    INTEGER nn(num_k)
    REAL hp1(3), hp2(3), min_dist

    REAL :: q(3)
    INTEGER i
    
    CALL ll2xyz(ll, q)
    CALL Set_qry(q)
    CALL Set_hps(hp1, hp2)
    CALL Set_tcs(hp1, hp2)
    CALL Search (q) 
    DO i = 1, num_k
      nn(i) = nni(i)
    END DO 
    min_dist = nnd(1)

END SUBROUTINE knn_search

SUBROUTINE knn_search_e(p, nn, nbs, min_dist)
    IMPLICIT NONE

    REAL, INTENT(IN) :: p(3)
    INTEGER, INTENT(OUT) :: nn(num_k)
    REAL, INTENT(OUT) :: nbs(3,3), min_dist

    REAL q(3)
    INTEGER i

    q = p
    CALL Search (q)
    DO i = 1, num_k
      nn(i) = nni(i)
      nbs(1:3,i) = ppoints(1:3, nni(i))
    END DO
    min_dist = nnd(1)

END SUBROUTINE knn_search_e

SUBROUTINE ll2xyz(p, e)
    IMPLICIT NONE
    REAL p(2)
    REAL e(3)

    e(1) = cos(p(1)) * cos(p(2))
    e(2) = cos(p(1)) * sin(p(2))
    e(3) = sin(p(1))

END SUBROUTINE ll2xyz

SUBROUTINE lls2xyzs(llpts, xyzpts, n)
    IMPLICIT NONE

    INTEGER :: n
    REAL :: llpts(n, 2), xyzpts(3, n)  
   
    INTEGER :: i

    DO i = 1, n
      xyzpts(1,i) = cos(llpts(i,1)) * cos(llpts(i,2))
      xyzpts(2,i) = cos(llpts(i,1)) * sin(llpts(i,2))
      xyzpts(3,i) = sin(llpts(i,1))
    END DO 

END SUBROUTINE lls2xyzs

END MODULE kd
