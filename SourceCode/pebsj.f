      subroutine pebsj(nr1,nr2,dk)
      use pealloc
      implicit none
      integer*4 nr1,nr2
      real*8 dk
c
      integer*4 i,ir,ik,ierr
      real*8 k,x
      real*8 bessj0,bessj1
c
      nbsj=nbsjmax
      do ik=1,nbsjmax
        k=dble(ik)*dk
        x=k*r0
        disk(ik)=(2.d0*bessj1(x)-x*bessj0(x))*(2.d0/x)**3
        if(dmax1(x,k*r(nr1)).ge.xbsjmax.and.ik.ge.nbsjmin)then
          nbsj=ik
          goto 100
        endif
      enddo
100   continue
c
c     re-allocate bsj
c
      deallocate(bsj)
      allocate(bsj(nbsj,0:1,nr),stat=ierr)
      if(ierr.ne.0)stop 'Error in pebsj: bsj not allocated!'
c
      do ik=1,nbsj
        k=dble(ik)*dk
        do ir=nr1,nr2
          x=k*r(ir)
          bsj(ik,0,ir)=bessj0(x)
          bsj(ik,1,ir)=bessj1(x)
        enddo
      enddo
c
      return
      end