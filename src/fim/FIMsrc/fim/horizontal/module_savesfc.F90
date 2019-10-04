module module_savesfc
contains
subroutine savesfc(filename)
 use module_control     ,only: &
 nts,nvl,nip,glvl,curve
 use module_sfc_variables
 implicit none
 character(len=17) filename
 integer lunout,idx,ipn
 ! save sfc variables at center of filter window.
 print*,'in savesfc',filename
!SMS$SERIAL BEGIN
 lunout = 77
 open (lunout,file=filename,form="unformatted",status="replace")
 call WriteGlvlHeader (lunout,glvl )
 call WriteCurveHeader(lunout,curve)
 do idx = 1,size(st3d,1)
   write(lunout) st3d(idx,:)
 enddo
 do idx = 1,size(sm3d,1)
   write(lunout) sm3d(idx,:)
 enddo
 do idx = 1,size(slc3d,1)
   write(lunout) slc3d(idx,:)
 enddo
 write(lunout) ts2d
 write(lunout) sheleg2d
 write(lunout) tg32d
 write(lunout) zorl2d
 write(lunout) cv2d
 write(lunout) cvb2d
 write(lunout) cvt2d
 write(lunout) alvsf2d
 write(lunout) alvwf2d
 write(lunout) alnsf2d
 write(lunout) alnwf2d
 write(lunout) slmsk2d
 write(lunout) vfrac2d
 write(lunout) canopy2d
 write(lunout) f10m2d
 write(lunout) t2m2d
 write(lunout) q2m2d
 write(lunout) vtype2d
 write(lunout) stype2d
 write(lunout) facsf2d
 write(lunout) facwf2d
 write(lunout) uustar2d
 write(lunout) ffmm2d
 write(lunout) ffhh2d
 write(lunout) hice2d
 write(lunout) fice2d
 write(lunout) tprcp2d
 write(lunout) srflag2d
 write(lunout) snwdph2d
 write(lunout) slc2d
 write(lunout) shdmin2d
 write(lunout) shdmax2d
 write(lunout) slope2d
 write(lunout) snoalb2d
 do idx = 1,size(hprm2d,1)
   write(lunout) hprm2d(idx,:)
 enddo
 close(lunout)
!SMS$SERIAL END
end subroutine savesfc
end module module_savesfc
