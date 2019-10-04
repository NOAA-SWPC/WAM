    PROGRAM ints_test

    REAL*8 p1(2), p2(2), p3(2), p4(2), p(2)

    REAL*8 d2r, r2d

    d2r = 4.0 * atan(1.0) / 180.0
    r2d = 1.0 / d2r

    p1(1) = 10.0 * d2r
    p1(2) = 30.0 * d2r
    p2(1) = 10.0 * d2r
    p2(2) = 60.0 * d2r

    p3(1) = 30.0 * d2r
    p3(2) = 40.0 * d2r
    p4(1) = 12.00 * d2r
    p4(2) = 60.0 * d2r


    CALL intersection(p1, p2, p3, p4, p)

    PRINT*, p1 * r2d, p2 * r2d
    PRINT*, p3 * r2d, p4 * r2d
    PRINT*, p * r2d

    END PROGRAM ints_test
