   subroutine lin2stp(xold,yold,kold,xnew,ynew,knew,vrbos,ipn)
!
! --- convert piecewise linear curve (xold,yold) into stairstep curve by
! --- integrating -yold- over consecutive -xnew- intervals
!
   implicit none
   integer,intent(IN) :: kold,knew
   integer,intent(IN) :: ipn		!  current location in horiz.grid
   real,intent(IN)    :: xold(kold),yold(kold),xnew(knew+1)
   real,intent(OUT)   :: ynew(knew)
   logical,intent(IN) :: vrbos		!  if true, print diagnostics
   integer k,ko
   real*8 colin,clout,colmx,yinteg,xlo,xhi,xa,xb,ya,yb,wgt
   real*8 xxold(kold),yyold(kold),xxnew(knew+1),yynew(knew)
   real,parameter     :: acurcy=1.e-9
!
   if (vrbos)								&
    write (*,101) ipn,'  lin2stp -- old profile:     x           y',	&
     (k,xold(k),yold(k),k=1,kold)
 101  format (i8,a/(i30,f12.1,es13.4))
!
   if (xold(1).lt.xold(kold)) then
    xxold= xold
    xxnew= xnew
   else
    xxold=-xold
    xxnew=-xnew
   end if
   yyold=yold
!
! --- column integrals (colin/clout) are computed for diagnostic purposes only
   if (xold(1).ne.xnew(1))						&
    print *,ipn,'  lin2stp warning - bottom xold and xnew differ',	&
     xold(1),xnew(1)
   if (xold(kold).ne.xnew(knew+1))					&
    print *,ipn,'  lin2stp warning - top xold and xnew differ',		&
     xold(kold),xnew(knew+1)
   colin=0.
   clout=0.
   colmx=0.
   do 3 k=1,kold-1
   colmx=max(colmx,abs(yyold(k)))
 3 colin=colin+.5*(yyold(k)+yyold(k+1))*(xxold(k+1)-xxold(k))
!
   do 4 k=1,knew
   xlo=xxnew(k  )
   xhi=xxnew(k+1)
   yynew(k)=yyold(1)
   if (xhi.gt.xlo) then
    yinteg=0.
    do ko=1,kold-1
     xa=max(xlo,min(xhi,xxold(ko  )))
     xb=max(xlo,min(xhi,xxold(ko+1)))
     if (xb.gt.xa) then
      wgt=(xa-xxold(ko))/(xxold(ko+1)-xxold(ko))
      ya=yyold(ko+1)*wgt+yyold(ko)*(1.-wgt)
      wgt=(xb-xxold(ko))/(xxold(ko+1)-xxold(ko))
      yb=yyold(ko+1)*wgt+yyold(ko)*(1.-wgt)
!     if (vrbos) print '(2i3,a,4f7.1,2es11.4)',k,ko,			&
!      ' xlo,xhi,xa,xb,ya,yb:',xlo,xhi,xa,xb,ya,yb
      yinteg=yinteg+.5*(ya+yb)*(xb-xa)
     end if
    end do
    yynew(k)=yinteg/(xhi-xlo)
!   if (vrbos) print '(i3,a,es12.4)',k,' ynew:',yynew(k)
    clout=clout+yinteg
   end if
 4 continue
   ynew=yynew
!
   if (abs(clout-colin).gt.acurcy*colmx*xold(kold/2))			&
    write (*,100) ipn,' lin2stp - column intgl.error',			&
      colin,clout,(clout-colin)/colin
 100  format (i8,a,2es14.6,es9.1)
!
   if (vrbos)								&
    write (*,101) ipn,'  lin2stp -- new profile:     x           y',	&
     (k,xnew(k),ynew(k),k=1,knew),knew+1,xnew(knew+1)
   return
   end subroutine lin2stp
