MODULE grib_datastru

        INTEGER tbl_sz, grib_lun
        PARAMETER(tbl_sz=300) !maximum table size
        CHARACTER*80 varnames(tbl_sz)
        CHARACTER*10 varabv(tbl_sz)
        INTEGER parm(tbl_sz), ztype(tbl_sz), iz1(tbl_sz),iz2(tbl_sz)
        INTEGER itime_range(tbl_sz), dscal(tbl_sz)
        INTEGER nvars_in_tbl

END MODULE grib_datastru
