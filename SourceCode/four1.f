      subroutine four1(dat,nn,isign)
      implicit none
      integer*4 nn,isign
      real*8 dat(2*nn)
c
c     fast fourier transform (fft)
c     convention: f(t)=\int F(f)e^{-i2\pi ft} df 
c     replace data by its discrete fourier transform, if isign is input
c     as 1; or replace data by nn times its inverse discrete fourier
c     transform, if isign is input as -1. data is a double complex array of
c     length nn or, equivalently, a double precision array of length 2*nn.
c     nn must be an integer power of 2 (this is not checked!)
c
c     note for convention: f(t)=\int F(f)e^{i2\pi f t} df, t-domain ==>
c     f-domain, if isign=-1, and f-domain ==> t-domain, if isign=1.
c
      integer*4 i,j,n,m,mmax,istep
      real*8 tempr,tempi
      real*8 wr,wi,wpr,wpi,wtemp,theta
      n=2*nn
      j=1
      do i=1,n,2
        if(j.gt.i)then
          tempr=dat(j)
          tempi=dat(j+1)
          dat(j)=dat(i)
          dat(j+1)=dat(i+1)
          dat(i)=tempr
          dat(i+1)=tempi
        endif
        m=n/2
1       if((m.ge.2).and.(j.gt.m))then
          j=j-m
          m=m/2
        goto 1
        endif
        j=j+m
      enddo
      mmax=2
2     if(n.gt.mmax)then
        istep=2*mmax
        theta=6.28318530717959d0/dble(isign*mmax)
        wpr=-2.d0*dsin(0.5d0*theta)**2
        wpi=dsin(theta)
        wr=1.d0
        wi=0.d0
        do m=1,mmax,2
          do i=m,n,istep
              j=i+mmax
              tempr=wr*dat(j)-wi*dat(j+1)
              tempi=wr*dat(j+1)+wi*dat(j)
              dat(j)=dat(i)-tempr
              dat(j+1)=dat(i+1)-tempi
              dat(i)=dat(i)+tempr
              dat(i+1)=dat(i+1)+tempi
          enddo
          wtemp=wr
          wr=wr*wpr-wi*wpi+wr
          wi=wi*wpr+wtemp*wpi+wi
        enddo
        mmax=istep
      goto 2
      endif
      return
      end
