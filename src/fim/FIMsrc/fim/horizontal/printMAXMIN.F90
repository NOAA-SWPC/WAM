module module_printMAXMIN
contains
subroutine printMAXMIN(its,nvl,nip,ntr,tr,dp3d)
implicit none
integer,intent(IN) :: its, nvl,nip,ntr
!SMS$DISTRIBUTE (dh,nip) BEGIN
real   ,intent(IN) :: tr  (nvl,nip,ntr)
real   ,intent(IN) :: dp3d(nvl,nip)
!SMS$DISTRIBUTE END
real               :: maxqv3d,minqv3d,aveqv3d,maxdp3d,mindp3d,avedp3d

maxqv3d = maxval(tr(:,1:nip,2))
minqv3d = minval(tr(:,1:nip,2))
aveqv3d = sum   (tr(:,1:nip,2))
maxdp3d = maxval(dp3d(:,1:nip))
mindp3d = minval(dp3d(:,1:nip))
avedp3d = sum   (dp3d(:,1:nip))
!SMS$REDUCE(maxqv3d,maxdp3d,max)
!SMS$REDUCE(minqv3d,mindp3d,min)
!SMS$REDUCE(aveqv3d,avedp3d,sum)
aveqv3d = aveqv3d/(nvl*nip)
avedp3d = avedp3d/(nvl*nip)
print"('MAXMIN -    ITS, qv(ivl=2)- max/ave/min,      dp(all levels)- max/ave/min')"
print"('MAXMIN',i10,1p6e10.3)",its,maxqv3d,aveqv3d,minqv3d,maxdp3d,avedp3d,mindp3d !All should be > 0
return
end subroutine printMAXMIN
end module module_printMAXMIN
