module MODULE_RELAX4E
  implicit none
contains
  subroutine relax4e(work,mask,relax_coeff,nrelax,  &
                     IDS,IDE,JDS,JDE, &
                     IMS,IME,JMS,JME, &
                     ITS,ITE,JTS,JTE)
    use MODULE_EXCHANGE, only: HALO_EXCH
    implicit none
    real, intent(in) :: relax_coeff
    integer, intent(in) :: nrelax
    real, intent(inout) :: work(ims:ime,jms:jme)
    integer, intent(in) :: mask(ims:ime,jms:jme)
    integer,intent(in) :: IDS,IDE,JDS,JDE, IMS,IME,JMS,JME, ITS,ITE,JTS,JTE

    real :: next(its:ite,jts:jte)
    real :: r,r1
    integer :: i,j,irelax,iter,i1,j1,i2,j2

    ! Caller should change domain dimension for V grid data, but does
    ! not have to change tile dimension.  That means we need to
    ! calculate new boundaries i=i1:i2, j=j1:j2 based on max/min of
    ! tile and domain dimensions.
    j1=max(jts,jds+1)
    j2=min(jte,jde-1)
    i1=max(its,ids+1)
    i2=min(ite,ide-1)

    ! Aliases to simplify the expressions below, we'll use "r" instead
    ! of relax_coeff, and r1 instead of 1-r:
    r=relax_coeff
    r1=1.0-r

    !$OMP PARALLEL DO PRIVATE(i,j)
    initj: do j=jts,jte
       initi: do i=its,ite
          next(i,j)=work(i,j)
       end do initi
    end do initj
    !$OMP END PARALLEL DO

    relaxloop: do irelax=1,nrelax
       call HALO_EXCH(work,1,1,1)

       !$omp parallel do private(i,j)
       bigj: do j=j1,j2
          bigi: do i=i1,i2
             if(mask(i,j)/=0) &
                next(i,j) = r1*work(i,j) + r*(&
                     work(i-1,j) + work(i,j-1) + work(i+1,j) + work(i,j+1) )/4
          enddo bigi
       enddo bigj
       ! Handle boundary points next.
       ! SOUTH:
       if(jts<=jds) then
          j=1
          !$omp parallel do private(i)
          do i=i1,i2
             if(mask(i,j)/=0) &
                next(i,j) = r1 * work(i,j) + r * &
                     (work(i-1,j) + work(i+1,j) + work(i,j+1))/3
          enddo
       endif
       ! NORTH:
       if(jte>=jde) then
          j=jde
          !$omp parallel do private(i)
          do i=i1,i2
             if(mask(i,j)/=0) &
                next(i,j) = r1 * work(i,j) + r * &
                     (work(i-1,j) + work(i+1,j) + work(i,j-1))/3
          enddo
       endif
       ! WEST:
       if(its<=ids) then
          i=1
          !$omp parallel do private(j)
          do j=j1,j2
             if(mask(i,j)/=0) &
                next(i,j) = r1 * work(i,j) + r * &
                     (work(i+1,j) + work(i,j-1) + work(i,j+1))/3
          enddo
       endif
       ! EAST:
       if(ite>=ide) then
          i=ide
          !$omp parallel do private(j)
          do j=j1,j2
             if(mask(i,j)/=0) &
                next(i,j) = r1 * work(i,j) + r * &
                     (work(i-1,j) + work(i,j-1) + work(i,j+1))/3
          enddo
       endif

       ! Finally, handle corner points:
       ! SOUTHWEST:
       if(its<=ids .and. jts<=jds) then
          if(mask(ids,jds)/=0) &
             next(ids,jds) = r1 * work(ids,jds) + r * &
                  (work(ids+1,jds) + work(ids,jds+1))/2
       endif
       ! SOUTHEAST:
       if(ite>=ide .and. jts<=jds) then
          if(mask(ide,jds)/=0) &
             next(ide,jds) = r1 * work(ide,jds) + r * &
                  (work(ide-1,jds) + work(ide,jds+1))/2
       endif
       ! NORTHWEST:
       if(its<=ids .and. jte>=jde) then
          if(mask(ids,jde)/=0) &
             next(ids,jde) = r1 * work(ids,jde) + r * &
                  (work(ids+1,jde) + work(ids,jde-1))/2
       endif
       ! NORTHEAST:
       if(ite>=ide .and. jte>=jde) then
          if(mask(ide,jde)/=0) &
             next(ide,jde) = r1 * work(ide,jde) + r * &
                  (work(ide-1,jde) + work(ide,jde-1))/2
       endif

       !$OMP PARALLEL DO PRIVATE(i,j)
       backj: do j=jts,jte
          backi: do i=its,ite
             work(i,j)=next(i,j)
          end do backi
       end do backj
       !$OMP END PARALLEL DO
    enddo relaxloop
  end subroutine relax4e
end module MODULE_RELAX4E
