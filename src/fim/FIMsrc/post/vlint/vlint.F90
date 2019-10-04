!=============================================================================
! subroutines to perform vertical interpolations
!
!
! Ning Wang, June 2007 
!
! 
!=============================================================================
! layer interpolation, target coordinate: 1100mb -- 0mb with -10mb increment 
! the interpolated data will have constant layers and transition layers 
      SUBROUTINE v_interp(mx, my, nvl, data_xyz_pr, data_xyz_var, nvgp, data_xyz, nodata, mode)
        IMPLICIT NONE
        INTEGER mx, my, nvl, nvgp, mode
        REAL data_xyz_pr(mx, my,nvl + 1),  data_xyz_var(mx, my, nvl), data_xyz(mx, my, nvgp)
        REAL nodata;
      
        REAL max_pr, min_pr, intv, extent, up_pr, dn_pr, mid_pr, last_up_val
        REAL l_dn_pr, l_up_pr, l_pr
        INTEGER i, j, k, l, up_idx, dn_idx, last_up_idx, up_b_idx

        IF (mode == 1) THEN
          extent = 0.85
        ELSE
          extent = 0.0
        ENDIF
        max_pr = 110000 !1100 mb
        min_pr = 0 !0 mb

        intv = (max_pr - min_pr) / (nvgp - 1)
 
        DO i = 1, mx   
          DO j = 1, my   
            last_up_val = data_xyz_var(i, j, 1)  
            last_up_idx = nvgp - (data_xyz_pr(i, j, 1) - min_pr) / intv  
            IF (last_up_idx >= 1) THEN
              data_xyz(i, j, last_up_idx) = last_up_val
            ENDIF

            DO k = 1, nvl  
              ! pressures for the upper and lower side of the current layer
              up_pr = data_xyz_pr(i, j, k + 1)
              dn_pr = data_xyz_pr(i, j, k)
              mid_pr = (up_pr + dn_pr) / 2.0 
              up_pr = mid_pr - extent * 0.5 * (dn_pr - up_pr) 
              dn_pr = mid_pr + extent * 0.5 * (dn_pr - up_pr) 

              ! fill in the constant layer value (middle)
              up_idx = nvgp - INT((up_pr - min_pr) / intv)
              dn_idx = nvgp - INT((dn_pr - min_pr) / intv) 
              DO l = dn_idx, up_idx
                data_xyz(i, j, l) = data_xyz_var(i, j, k)
              END DO 

              ! linearly interpolate the transition layer (down, between two layers)
              l_dn_pr = log(REAL((nvgp - last_up_idx) * intv + min_pr))
              l_up_pr = log(REAL((nvgp - dn_idx + 1) * intv + min_pr)) 
              DO l = last_up_idx + 1, dn_idx - 1
                l_pr = log(REAL((nvgp - l) * intv + min_pr))
                data_xyz(i, j, l) = last_up_val + (l_pr - l_dn_pr) *  &
                 (data_xyz(i, j, dn_idx) - last_up_val) / & 
                 (l_up_pr - l_dn_pr)   
              END DO 
              
              ! if it is the lowest layer, fill in the no value for the data points below
              IF (k == 1) THEN 
                DO l = 1, last_up_idx
                  data_xyz(i, j, l) = nodata
                END DO

              ! if it is the highest layer, interpolate to the interface 
              ELSE IF (k == nvl) THEN
                up_b_idx = nvgp - (data_xyz_pr(i, j, nvl + 1) - min_pr) / intv
                l_dn_pr = log(REAL((nvgp - up_idx) * intv + min_pr))
                l_up_pr = log(REAL((nvgp - up_b_idx) * intv + min_pr)) 
                DO l = up_idx + 1, up_b_idx
                  l_pr = log(REAL((nvgp - l) * intv + min_pr))
                  data_xyz(i, j, l) = data_xyz(i, j, up_idx) + (l_pr - l_dn_pr) *  &
                   (data_xyz_var(i, j, nvl) - data_xyz(i, j, up_idx)) / & 
                    (l_up_pr - l_dn_pr)   
                END DO 
                ! then fill in the nodata value for the data points above
                DO l = up_b_idx + 1, nvgp
                  data_xyz(i, j, l) = nodata
                END DO
              ENDIF 
              last_up_idx = up_idx
              last_up_val = data_xyz(i, j, up_idx) 
            END DO
          END DO
        END DO

      END SUBROUTINE

! level interpolation, target coordinate: 1100mb -- 0mb with -10mb increment 
      SUBROUTINE v_interp_lvlvar(mx, my, nvl, data_xyz_pr, data_xyz_var, nvgp, data_xyz, nodata, mode)
        IMPLICIT NONE
        INTEGER mx, my, nvl, nvgp, mode 	
        REAL data_xyz_pr(mx, my,nvl + 1),  data_xyz_var(mx, my, nvl), data_xyz(mx, my, nvgp)
        REAL nodata;
      
        REAL max_pr, min_pr, intv, extent, up_pr, dn_pr, mid_pr, last_up_val
        REAL l_dn_pr, l_up_pr, l_pr, pi_dn, pi_up, pi_pr
        INTEGER i, j, k, l, up_idx, dn_idx, last_up_idx, up_b_idx

        max_pr = 110000.00 !1100 mb
        min_pr = 0.00 !0 mb

        intv = (max_pr - min_pr) / (nvgp - 1)
 
        DO i = 1, mx   
          DO j = 1, my   
            DO k = 1, nvl - 1 
              ! pressures for the upper and lower side of the current layer
              up_pr = data_xyz_pr(i, j, k + 1)
              dn_pr = data_xyz_pr(i, j, k)
              up_idx = nvgp - (up_pr - min_pr) / intv
              dn_idx = nvgp - (dn_pr - min_pr) / intv 
              ! linearly interpolate the transition layer
              pi_dn = (dn_pr / 100000.00)**0.286
              pi_up = (up_pr / 100000.00)**0.286 
              DO l = dn_idx, up_idx
                l_pr = (nvgp - l) * intv + min_pr
                pi_pr = (l_pr / 100000.00)**0.286
                data_xyz(i, j, l) = data_xyz_var(i, j, k) + (pi_pr - pi_dn) *  &
                 (data_xyz_var(i, j, k + 1) - data_xyz_var(i, j, k)) / & 
                 (pi_up - pi_dn)   
              END DO 
              
              ! if it is the lowest layer, fill in the no value for the data points below
              IF (k == 1) THEN 
                DO l = 1, dn_idx - 1
                  data_xyz(i, j, l) = nodata
                END DO
              ! then fill in the nodata value for the data points above
              ELSE IF (k == nvl) THEN
                DO l = up_idx + 1, nvgp
                  data_xyz(i, j, l) = nodata
                END DO
              ENDIF 
            END DO
          END DO
        END DO

      END SUBROUTINE

!=============================================================================
! Level and layer interpolation, target coordinate: passed in (v_coor) 
!
!
! Figure 1.                          | Figure 2.
! Level variables at interface:      | Layer variables: 
!                                    |
! -------------- int. level nvl + 1  | ---------------- int. level nvl + 1 
!       :                            |         :
!       :                            |         :
! -------------- int. level k + 1    | ---------------- int. level k + 1
!                                    | //////////////// variable at layer k
! -------------- int. level k        | ---------------- int. level k
!       :                            |         :
!       :                            |         :
! -------------- int. level 1        | ---------------- int. level 1
!                                    |
! 
!
! N. Wang, Feb. 2008
!=============================================================================

      SUBROUTINE vlint2coor(mx, my, nvl, nvlp1, data_xyz_pr, data_xyz_var, data_xyz, v_coor, nvc)
        IMPLICIT NONE
        INTEGER mx, my, nvl, nvlp1, nvc
        REAL data_xyz_pr(mx, my, nvlp1),  data_xyz_var(mx, my, nvl), data_xyz(mx, my, nvc), v_coor(nvc)
        REAL v_coor_pa(nvc), pi_dn, pi_up, pi_co, dn_val, up_val
        INTEGER i, j, k, l

        v_coor_pa = v_coor * 100.0
        DO i = 1, mx   
          DO j = 1, my   
            k = 1
            DO l = 1, nvc 
              IF (v_coor_pa(l) >= data_xyz_pr(i, j, 1)) THEN
                data_xyz(i, j, l) = data_xyz_var(i, j, 1)
                CYCLE
              END IF
              IF (v_coor_pa(l) <= data_xyz_pr(i, j, nvlp1)) THEN
                data_xyz(i, j, l) = data_xyz_var(i, j, nvl)
                CYCLE
              END IF
              DO WHILE (v_coor_pa(l) < data_xyz_pr(i, j, k + 1)) 
                IF (k == nvlp1 - 1) THEN
                  EXIT 
                ELSE
                  k = k + 1
                ENDIF
              END DO ! k and k+1 are the current indexes for interpolation  
              IF (nvl == nvlp1) THEN ! level variables, see fig. 1.
                pi_dn = (data_xyz_pr(i, j, k) / 100000.00)**0.286
                pi_up = (data_xyz_pr(i, j, k + 1) / 100000.00)**0.286 
                pi_co = (v_coor_pa(l) / 100000.00)**0.286
                dn_val = data_xyz_var(i, j, k)
                up_val = data_xyz_var(i, j, k + 1)
              ELSE  ! layer variables, see fig. 2.
                pi_dn = (data_xyz_pr(i, j, k) / 100000.00)**0.286
                pi_up = (data_xyz_pr(i, j, k + 1) / 100000.00)**0.286 
                pi_co = (v_coor_pa(l) / 100000.00)**0.286
                IF (pi_co > (pi_dn + pi_up) / 2.0) THEN ! lower half of the layer
                  IF (k == 1) THEN
                    pi_dn = (data_xyz_pr(i, j, k) / 100000.00)**0.286
                    pi_up = ((data_xyz_pr(i, j, k) / 100000.00)**0.286 + (data_xyz_pr(i, j, k + 1) / 100000.00)**0.286)/ 2.0
                    pi_co = (v_coor_pa(l) / 100000.00)**0.286
                    dn_val = data_xyz_var(i, j, k)
                    up_val = dn_val 
                  ELSE
                    pi_dn = ((data_xyz_pr(i, j, k) / 100000.00)**0.286 + (data_xyz_pr(i, j, k - 1) / 100000.00)**0.286)/ 2.0 
                    pi_up = ((data_xyz_pr(i, j, k) / 100000.00)**0.286 + (data_xyz_pr(i, j, k + 1) / 100000.00)**0.286)/ 2.0 
                    pi_co = (v_coor_pa(l) / 100000.00)**0.286
                    dn_val = data_xyz_var(i, j, k - 1)
                    up_val = data_xyz_var(i, j, k)
                  ENDIF
                ELSE ! upper half of the layer
                  IF (k == nvl) THEN
                    pi_dn = ((data_xyz_pr(i, j, k) / 100000.00)**0.286 + (data_xyz_pr(i, j, k + 1) / 100000.00)**0.286)/ 2.0 
                    pi_up = (data_xyz_pr(i, j, k + 1) / 100000.00)**0.286
                    pi_co = (v_coor_pa(l) / 100000.00)**0.286
                    dn_val = data_xyz_var(i, j, k)
                    up_val = dn_val 
                  ELSE
                    pi_dn = ((data_xyz_pr(i, j, k) / 100000.00)**0.286 + (data_xyz_pr(i, j, k + 1) / 100000.00)**0.286)/ 2.0
                    pi_up = ((data_xyz_pr(i, j, k + 1) / 100000.00)**0.286 + (data_xyz_pr(i, j, k + 2) / 100000.00)**0.286)/ 2.0 
                    pi_co = (v_coor_pa(l) / 100000.00)**0.286
                    dn_val = data_xyz_var(i, j, k)
                    up_val = data_xyz_var(i, j, k + 1)
                  ENDIF

                ENDIF
              ENDIF
              data_xyz(i, j, l) = up_val +  (pi_co - pi_up) * &
                (dn_val - up_val) / (pi_dn - pi_up) 
            END DO
          END DO
        END DO
      END SUBROUTINE

