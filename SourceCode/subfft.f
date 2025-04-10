      subroutine subfft(fgrn,tgrn,nfg,nfmax,dfg,dtg,alpha,
     &                  ts,sinj,nts,obs,nt,dt,fgrn0,intdif)
      implicit none
c
      integer*4 nfg,nfmax,nts,nt,intdif
      real*8 dfg,dtg,dt,alpha,delta
      real*8 ts(nts),sinj(nts)
      complex*16 fgrn0
      real*8 tgrn(2*nfmax)
      real*8 obs(4*nfmax)
      complex*16 fgrn(nfmax)
c
c     local memory
c
      integer*4 i,it,itga,itgb,lf,jf,ntg,its,jts
      real*8 a,b,t,pi2,coseis
c
      real*8 eps
      data eps/1.0d-02/
c
      pi2=8.d0*datan(1.d0)
      ntg=2*nfg
c
      coseis=dreal(fgrn0)
c
      do lf=1,nfg
        fgrn(lf)=fgrn(lf)-dcmplx(coseis,0.d0)
      enddo
c
      do lf=1,nfg
        obs(2*lf-1)=dreal(fgrn(lf))
        obs(2*lf  )=dimag(fgrn(lf))
      enddo
c
      jf=1
      do lf=2*nfg,nfg+2,-1
        jf=jf+1
        obs(2*lf-1)= dreal(fgrn(jf))
        obs(2*lf  )=-dimag(fgrn(jf))
      enddo
      obs(2*nfg+1)=0.d0
      obs(2*nfg+2)=0.d0
c
c	convention for Fourier transform:
c	f(t)=\int F(f) exp(i2\pi f t) df
c
      call four1(obs(1),2*nfg,+1)
c
c     calculate response to the Heaviside source history
c
      tgrn(1)=0.d0
      do it=2,ntg
        tgrn(it)=tgrn(it-1)+obs(2*it-1)
     &          *dtg*dfg*dexp(alpha*dble(it-1)*dtg)
      enddo
c
c     convolution with source time function
c
      do it=1,2*((1+nt)/2)
        t=dble(it-1)*dt
        obs(it)=0.d0
        do its=1,nts
          if(ts(its).le.t)then
            b=dmod((t-ts(its))/dtg,1.d0)
            a=1.d0-b
            itga=min0(1+idint((t-ts(its))/dtg),ntg)
            itgb=min0(2+idint((t-ts(its))/dtg),ntg)
            obs(it)=obs(it)+sinj(its)*(coseis+tgrn(itga)*a+tgrn(itgb)*b)
          endif
        enddo
      enddo
c
      if(intdif.eq.1)then
        do it=2,2*((1+nt)/2)-1
          obs(it)=(obs(it+1)-obs(it))/dt
        enddo
        obs(1)=2.d0*obs(2)-obs(3)
        obs(nt)=obs(nt-1)
      else if(intdif.eq.-1)then
        obs(1)=0.d0
        do it=2,2*((1+nt)/2)
          obs(it)=obs(it-1)+obs(it)*dt
        enddo
      endif
c
      do it=1,(1+nt)/2
        fgrn(it)=dcmplx(obs(2*it-1),obs(2*it))
      enddo
c
      return
      end