      subroutine setlats_r_ext(lats_nodes_r,lats_nodes_ext,
     .                   global_lats_r,
     &                   global_lats_ext,iprint,lonsperlat)
cc
      use resol_def, ONLY: latr, jintmx, nypt
      use layout1,   ONLY: nodes
!jw      use mpi_def,   ONLY: icolor, liope
      implicit none
cc
      integer              lats_nodes_r(nodes)
cc
      integer              global_lats_r(latr)
      integer              lats_nodes_ext(nodes)
      integer        global_lats_ext(latr+2*jintmx+2*nypt*(nodes-1))
cc
      integer              iprint,opt,ifin,nodesio
cc
      integer              lonsperlat(latr)
cc
      integer              ijk,jcount,jpt,lat,lats_sum,node,i,j
      integer              ILATPE,ngrptg,ngrptl,ipe,irest,idp
cc
cc
      OPT=1
!jw      if (liope) then
!jw         if (icolor.eq.2) then
!jw           nodesio=1
!jw         else
!jw           nodesio=nodes
!jw         endif
!jw      else
         nodesio=nodes
!jw      endif
!!
      do node=1,nodesio
        if (nodesio.eq.1) then
          lats_nodes_ext(node)=lats_nodes_r(node)+2*jintmx
        else
          if (node.eq.1.or.node.eq.nodesio) then
           lats_nodes_ext(node)=lats_nodes_r(node)+jintmx+nypt
          else
           lats_nodes_ext(node)=lats_nodes_r(node)+2*nypt
          endif
        endif
      enddo
cc........................................................
cc
cc
      jpt=0
      do node=1,nodesio
       if (nodesio.eq.1) then
        do i=1,jintmx
          global_lats_ext(i)=global_lats_r(1)
        enddo
        do i=1,jintmx
          global_lats_ext(jintmx+latr+i)=global_lats_r(latr)
        enddo
        do i=1,latr
          global_lats_ext(i+jintmx)=global_lats_r(i)
        enddo
       else
        do jcount=1,lats_nodes_r(node)
         global_lats_ext(jpt+jintmx+jcount+2*nypt*(node-1))=
     &               global_lats_r(jpt+jcount)
        enddo
        if (node.eq.1) then
         do i=1,jintmx
           global_lats_ext(i)=global_lats_r(1)
         enddo
         do i=1,nypt
           global_lats_ext(jintmx+lats_nodes_r(node)+i)=
     &               global_lats_r(lats_nodes_r(node))+i
         enddo
        elseif (node.eq.nodesio) then
         do i=1,jintmx
           global_lats_ext(latr+jintmx+2*nypt*(nodesio-1)+i)=
     &                    global_lats_r(latr)
         enddo
         do i=nypt,1,-1
           global_lats_ext(jpt+jintmx+2*nypt*(node-1)-i+1)=
     &                    global_lats_r(jpt)-i+1
         enddo
        else
         do i=nypt,1,-1
           global_lats_ext(jpt+jintmx+2*nypt*(node-1)-i+1)=
     &                    global_lats_r(jpt)-i+1
         enddo
         do i=1,nypt
         global_lats_ext(jpt+jintmx+2*nypt*(node-1)+
     &                    lats_nodes_r(node)+i)=
     &              global_lats_r(jpt+lats_nodes_r(node))+i
         enddo
        endif
       endif
        jpt=jpt+lats_nodes_r(node)
      enddo
cc
cc
 
      if ( iprint .ne. 1 ) return
cc
      jpt=0
      do node=1,nodesio
         if ( lats_nodes_r(node) .gt. 0 ) then
            print 600
            lats_sum=0
            do jcount=1,lats_nodes_r(node)
               lats_sum=lats_sum + lonsperlat(global_lats_r(jpt+jcount))
               print 701, node-1,
     x                    node,    lats_nodes_r(node),
     x                    jpt+jcount, global_lats_r(jpt+jcount)
!selax                     lonsperlat(global_lats_r(jpt+jcount)),
!selax                    lats_sum
            enddo
         endif
         jpt=jpt+lats_nodes_r(node)
      enddo
cc
      print 600
cc
  600 format ( ' ' )
cc
  701 format (  'setlats  me=', i4,
     x          '  lats_nodes_r(',  i4, ' )=', i4,
     x          '  global_lats_r(', i4, ' )=', i4)
  700 format (  'setlats  me=', i4,
     x          '  lats_nodes_r(',  i4, ' )=', i4,
     x          '  global_lats_r(', i4, ' )=', i4,
     x          '  lonsperlat=', i5,
     x          '  lats_sum=',   i6 )
cc
      return
      end
c
      subroutine setlats_r_ext_shuff(lats_nodes_r,lats_nodes_ext,
     .                   global_lats_r,
     &                   global_lats_ext,iprint,lonsperlat)
cc
      use resol_def, ONLY: latr, jintmx, nypt
      use layout1,   ONLY: nodes
!jw      use mpi_def,   ONLY: icolor, liope
      implicit none
cc
      integer              lats_nodes_r(nodes)
cc
      integer              global_lats_r(latr)
      integer              lats_nodes_ext(nodes)
      integer        global_lats_ext(latr+2*jintmx+2*nypt*(nodes-1))
cc
      integer              iprint,opt,ifin,nodesio
cc
      integer              lonsperlat(latr)
cc
      integer              ijk,jcount,jpt,lat,lats_sum,node,i,j
      integer              ILATPE,ngrptg,ngrptl,ipe,irest,idp
cc
cc
      OPT=1
!jw      if (liope) then
!jw         if (icolor.eq.2) then
!jw           nodesio=1
!jw         else
!jw           nodesio=nodes
!jw         endif
!jw      else
         nodesio=nodes
!jw      endif
!!
      do node=1,nodesio
        if (nodesio.eq.1) then
          lats_nodes_ext(node)=lats_nodes_r(node)+2*jintmx
        else
          if (node.eq.1.or.node.eq.nodesio) then
           lats_nodes_ext(node)=lats_nodes_r(node)+jintmx+nypt
          else
           lats_nodes_ext(node)=lats_nodes_r(node)+2*nypt
          endif
        endif
      enddo
cc........................................................
cc
cc
c$$$      jpt=0
c$$$      do node=1,nodesio
c$$$       if (nodesio.eq.1) then
c$$$        do i=1,jintmx
c$$$          global_lats_ext(i)=global_lats_r(1)
c$$$        enddo
c$$$        do i=1,jintmx
c$$$          global_lats_ext(jintmx+latr+i)=global_lats_r(latr)
c$$$        enddo
c$$$        do i=1,latr
c$$$          global_lats_ext(i+jintmx)=global_lats_r(i)
c$$$        enddo
c$$$       else
c$$$        do jcount=1,lats_nodes_r(node)
c$$$         global_lats_ext(jpt+jintmx+jcount+2*nypt*(node-1))=
c$$$     &               global_lats_r(jpt+jcount)
c$$$        enddo
c$$$        if (node.eq.1) then
c$$$         do i=1,jintmx
c$$$           global_lats_ext(i)=global_lats_r(1)
c$$$         enddo
c$$$         do i=1,nypt
c$$$           global_lats_ext(jintmx+lats_nodes_r(node)+i)=
c$$$     &               global_lats_r(lats_nodes_r(node))+i
c$$$         enddo
c$$$        elseif (node.eq.nodesio) then
c$$$         do i=1,jintmx
c$$$           global_lats_ext(latr+jintmx+2*nypt*(nodesio-1)+i)=
c$$$     &                    global_lats_r(latr)
c$$$         enddo
c$$$         do i=nypt,1,-1
c$$$           global_lats_ext(jpt+jintmx+2*nypt*(node-1)-i+1)=
c$$$     &                    global_lats_r(jpt)-i+1
c$$$         enddo
c$$$        else
c$$$         do i=nypt,1,-1
c$$$           global_lats_ext(jpt+jintmx+2*nypt*(node-1)-i+1)=
c$$$     &                    global_lats_r(jpt)-i+1
c$$$         enddo
c$$$         do i=1,nypt
c$$$         global_lats_ext(jpt+jintmx+2*nypt*(node-1)+
c$$$     &                    lats_nodes_r(node)+i)=
c$$$     &              global_lats_r(jpt+lats_nodes_r(node))+i
c$$$         enddo
c$$$        endif
c$$$       endif
c$$$        jpt=jpt+lats_nodes_r(node)
c$$$      enddo
cc
cc
 
      if ( iprint .ne. 1 ) return
cc
      jpt=0
      do node=1,nodesio
         if ( lats_nodes_r(node) .gt. 0 ) then
            print 600
            lats_sum=0
            do jcount=1,lats_nodes_r(node)
               lats_sum=lats_sum + lonsperlat(global_lats_r(jpt+jcount))
               print 701, node-1,
     x                    node,    lats_nodes_r(node),
     x                    jpt+jcount, global_lats_r(jpt+jcount)
!selax                     lonsperlat(global_lats_r(jpt+jcount)),
!selax                    lats_sum
            enddo
         endif
         jpt=jpt+lats_nodes_r(node)
      enddo
cc
      print 600
cc
  600 format ( ' ' )
cc
  701 format (  'setlats  me=', i4,
     x          '  lats_nodes_r(',  i4, ' )=', i4,
     x          '  global_lats_r(', i4, ' )=', i4)
  700 format (  'setlats  me=', i4,
     x          '  lats_nodes_r(',  i4, ' )=', i4,
     x          '  global_lats_r(', i4, ' )=', i4,
     x          '  lonsperlat=', i5,
     x          '  lats_sum=',   i6 )
cc
      return
      end
