      subroutine pekern(cs,k)
      use pealloc
      implicit none
c
c     calculation of response function in Laplace domain
c     k: wave number
c     cs: complex Laplace variable
c     uz,ur = displacement
c     p = pore pressure
c     ezz, ess, ezr = normal, surface, shear strain
c     dpz = pore pressure vertical gradient
c
      real*8 k
      complex*16 cs
c
      integer*4 i,n,ly,izr
      real*8 beti
      complex*16 ck,cla,cmu,chi,cla0,cmu0,chi0
c
      real*8 eps
      data eps/1.0d-03/
c
      ck=dcmplx(k,0.d0)
c
      do n=1,n0
        beti=(1.d0+alf(n)*qa(n)/(la(n)+2.d0*mu(n)))/dm(n)
        ka(n)=cdsqrt(ck*ck+cs*dcmplx(beti,0.d0))
        smalls(n)=cdabs(ka(n)-ck)/k.lt.eps
        if(smalls(n))then
          ck3(n)=ck
        else
          ck3(n)=ka(n)
        endif
      enddo
c
      call pepsv(cs,k)
c
      do izr=1,nzr
        ly=lzr(izr)
        n=nno(ly)
        cla=dcmplx(la(n),0.d0)
        cmu=dcmplx(mu(n),0.d0)
        chi=dcmplx(dm(n)*alf(n)/qa(n),0.d0)
c
        cuz(izr)=ypsv(1,ly)
        cur(izr)=ypsv(3,ly)
        cp(izr)=ypsv(5,ly)
c
        cezz(izr)=(ypsv(2,ly)+cla*ck*ypsv(3,ly))/(cla+(2.d0,0.d0)*cmu)
        cess(izr)=-ck*ypsv(3,ly)
        cezr(izr)=(0.5d0,0.d0)*ypsv(4,ly)/cmu
        cdpz(izr)=-ypsv(6,ly)/chi
      enddo
c
      return
      end