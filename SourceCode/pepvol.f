      subroutine pepvol(cs,spv)
      use pealloc
      implicit none
      complex*16 cs,spv
c
      integer*4 i,j,ik,n,nk
      real*8 k,dk,x,sdisk,bsj0,bsj1
      complex*16 cps,spv0
      complex*16 dmai(6,6)
c
      real*8 bessj0,bessj1
c
c     constants
c
      real*8 pi,pi2,eps
      integer*4 dnk
c
      data dnk/128/
      data pi,pi2,eps/3.14159265358979d0,6.28318530717959d0,1.0d-03/
c
c     integrate pore pressure over the source volume
c
      dk=pi2/sradius/dble(dnk)
      n=nno(lstop)
c
      ik=0
      nk=dnk
      spv0=(0.d0,0.d0)
      spv=(0.d0,0.d0)
20    continue
      do i=1,nk
        ik=ik+1
        k=dble(ik)*dk
        x=k*sradius
        bsj0=bessj0(x)
        bsj1=bessj1(x)
        sdisk=(2.d0*bsj1-x*bsj0)*(2.d0/x)**3
        call pekern(cs,k)
        if(disksource)then
          cps=ypsv(5,lstop)
        else
          call pedifmai(dmai,cs,k,n)
          ypsv(2,lstop)=ypsv(2,lstop)
     &                 -dcmplx(alf(n),0.d0)*ypsv(5,lstop)
          ypsv(6,lstop)=ypsv(6,lstop)
     &                 +dcmplx(alf(n),0.d0)*cs*ypsv(1,lstop)
c
          ypsv(2,lsbtm)=ypsv(2,lsbtm)
     &                 -dcmplx(alf(n),0.d0)*ypsv(5,lsbtm)
          ypsv(6,lsbtm)=ypsv(6,lsbtm)
     &                 +dcmplx(alf(n),0.d0)*cs*ypsv(1,lsbtm)
c
          cps=(0.d0,0.d0)
          do j=1,6
            cps=cps+dmai(5,j)*((ypsv(j,lsbtm)-ypsv(j,lstop))
     &                        /dcmplx(slength,0.d0)
     &         +sfct(j))
          enddo
        endif
        spv=spv+cps*dcmplx(sdisk*bsj1,0.d0)
      enddo
      if(cdabs(spv-spv0).gt.eps*cdabs(spv).and.ik.lt.nbsjmax)then
        spv0=spv
        nk=2*nk
        goto 20
      endif
      spv=spv*dcmplx(pi2*slength*sradius*dk,0.d0)
c
      return
      end