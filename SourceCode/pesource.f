      subroutine pesource(cs)
      use pealloc
      implicit none
      complex*16 cs
c
      integer*4 i,n
      real*8 qbiot
      complex*16 spv
c
c     constants
c
      real*8 pi,pi2
      data pi,pi2/3.14159265358979d0,6.28318530717959d0/
c
      do i=1,6
        sfct(i)=(0.d0,0.d0)
      enddo
      if(istype.eq.0)then
        if(disksource)then
          sfct(2)=dcmplx(1.d0/pi2,0.d0)
        else
          sfct(2)=dcmplx(1.d0/pi2/slength,0.d0)
        endif
        return
      else if(istype.eq.1)then
        n=nno(lstop)
        sfct(6)=dcmplx(alf(n)/pi2/slength,0.d0)
        return
      else if(istype.eq.2)then
        if(disksource)then
          sfct(6)=dcmplx(-1.d0/pi2,0.d0)
        else
          sfct(6)=dcmplx(-1.d0/pi2/slength,0.d0)
        endif
        return
      else if(istype.eq.3)then
        sfct(6)=dcmplx(-1.d0/pi2/slength,0.d0)
        call pepvol(cs,spv)
        sfct(6)=dcmplx(-1.d0/pi2/slength,0.d0)/spv
      else
        stop ' Error in pesource: wrong source geometry selected!'
      endif
c
      return
      end