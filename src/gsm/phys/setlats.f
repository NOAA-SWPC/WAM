      subroutine setlats_r(lats_nodes_r,lats_nodes_ext,global_lats_r,
     &                     global_lats_ext,iprint,lonsperlar)
!
      use resol_def, ONLY: latr, jintmx, nypt, lonr
      use layout1,   ONLY: nodes
!jw      use mpi_def,   ONLY: icolor, liope
      implicit none
!
      integer              lats_nodes_r(nodes)
!
      integer              global_lats_r(latr)
      integer              lats_nodes_ext(nodes)
      integer        global_lats_ext(latr+2*jintmx+2*nypt*(nodes-1))
!
      integer              iprint,opt,ifin,nodesio
!
      integer              lonsperlar(latr)
!
      integer              ijk,jcount,jpt,lat,lats_sum,node,i,j
      integer              ILATPE,ngrptg,ngrptl,ipe,irest,idp
!
      integer,allocatable :: lats_hold(:,:)
!
!     print *,' enter setlats_r ',latr,nodes
      write(0,*)' lonsperlar in setlats_r=',lonsperlar

      allocate ( lats_hold(latr,nodes) )
!
!     iprint=1
      iprint=0
      OPT=1
      lats_nodes_r=0
!jw      if (liope) then
!jw         if (icolor.eq.2) then
!jw           nodesio=1
!jw         else
!jw           nodesio=nodes
!jw         endif
!jw      else
         nodesio=nodes
!jw      endif
!
      ngrptg=0
      do lat=1,latr
         if (opt.eq.1) then
           ifin=lonsperlar(lat)
         elseif(opt.eq.2) then
           ifin=lonr
         endif
         do i=1,ifin
            ngrptg=ngrptg+1
         enddo
      enddo
!
!     distribution of the grid
      ILATPE=ngrptg/nodesio
      ngrptl=0
      ipe=0
      irest=0
      idp=1
      do lat=1,latr
         if (opt.eq.1) then
           ifin=lonsperlar(lat)
         elseif(opt.eq.2) then
           ifin=lonr
         endif
          ngrptl=ngrptl+ifin
       if (ngrptl*nodesio.le.ngrptg+irest) then
           lats_nodes_r(ipe+1) = lats_nodes_r(ipe+1) + 1
           lats_hold(idp,ipe+1)=lat
           idp=idp+1
       else
          ipe=ipe+1
          if (ipe.le.nodesio) lats_hold(1,ipe+1) = lat
          idp=2
          irest=irest+ngrptg-(ngrptl-ifin)*nodesio
          ngrptl=ifin
          lats_nodes_r(ipe+1) = lats_nodes_r(ipe+1) + 1
       endif
      enddo
!!
      write(0,*)' lonsperlar in setlats_r2=',lonsperlar
!     do node=1,nodesio
!       if (nodesio.eq.1) then
!         lats_nodes_ext(node)=lats_nodes_r(node)+2*jintmx
!       else
!         if (node.eq.1.or.node.eq.nodesio) then
!          lats_nodes_ext(node)=lats_nodes_r(node)+jintmx+nypt
!         else
!          lats_nodes_ext(node)=lats_nodes_r(node)+2*nypt
!         endif
!       endif
!     enddo
!     write(0,*)' lonsperlar in setlats_r3=',lonsperlar
!........................................................
!
      jpt=0
      do node=1,nodesio
         if ( lats_nodes_r(node) .gt. 0 ) then
            do jcount=1,lats_nodes_r(node)
               global_lats_r(jpt+jcount) = lats_hold(jcount,node)
            enddo
         endif
         jpt=jpt+lats_nodes_r(node)
      enddo
!
!     jpt=0
!     do node=1,nodesio
!      if (nodesio.eq.1) then
!       do i=1,jintmx
!         global_lats_ext(i)=global_lats_r(1)
!       enddo
!       do i=1,jintmx
!         global_lats_ext(jintmx+latr+i)=global_lats_r(latr)
!       enddo
!       do i=1,latr
!         global_lats_ext(i+jintmx)=global_lats_r(i)
!       enddo
!      else
!       do jcount=1,lats_nodes_r(node)
!        global_lats_ext(jpt+jintmx+jcount+2*nypt*(node-1))=
!    &               global_lats_r(jpt+jcount)
!       enddo
!     write(0,*)' lonsperlar in setlats_r4=',lonsperlar
!       if (node.eq.1) then
!        do i=1,jintmx
!          global_lats_ext(i)=global_lats_r(1)
!        enddo
!        do i=1,nypt
!          global_lats_ext(jintmx+lats_nodes_r(node)+i)=
!    &               global_lats_r(lats_nodes_r(node))+i
!        enddo
!       elseif (node.eq.nodesio) then
!        do i=1,jintmx
!          global_lats_ext(latr+jintmx+2*nypt*(nodesio-1)+i)=
!    &                    global_lats_r(latr)
!        enddo
!        do i=nypt,1,-1
!          global_lats_ext(jpt+jintmx+2*nypt*(node-1)-i+1)=
!    &                    global_lats_r(jpt)-i+1
!        enddo
!       else
!        do i=nypt,1,-1
!          global_lats_ext(jpt+jintmx+2*nypt*(node-1)-i+1)=
!    &                    global_lats_r(jpt)-i+1
!        enddo
!        do i=1,nypt
!        global_lats_ext(jpt+jintmx+2*nypt*(node-1)+
!    &                    lats_nodes_r(node)+i)=
!    &              global_lats_r(jpt+lats_nodes_r(node))+i
!        enddo
!       endif
!     write(0,*)' lonsperlar in setlats_r5=',lonsperlar
!      endif
!       jpt=jpt+lats_nodes_r(node)
!     enddo
!
      write(0,*)' lonsperlar in setlats_r6=',lonsperlar
      if ( iprint .ne. 1 ) return
!
      jpt=0
      do node=1,nodesio
         if ( lats_nodes_r(node) .gt. 0 ) then
            print 600
            lats_sum=0
            do jcount=1,lats_nodes_r(node)
               lats_sum=lats_sum + lonsperlar(global_lats_r(jpt+jcount))
               print 700, node-1,
     &                    node,    lats_nodes_r(node),
     &                    jpt+jcount, global_lats_r(jpt+jcount),
     &                     lonsperlar(global_lats_r(jpt+jcount)),
     &                    lats_sum
            enddo
         endif
         jpt=jpt+lats_nodes_r(node)
      enddo
!
      print 600
!
  600 format ( ' ' )
!
  700 format (  'setlats  me=', i4,
     x          '  lats_nodes_r(',  i4, ' )=', i4,
     x          '  global_lats_r(', i4, ' )=', i4,
     x          '  lonsperlar=', i5,
     x          '  lats_sum=',   i6 )
!
      deallocate ( lats_hold )
!
      return
      end
