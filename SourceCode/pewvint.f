      subroutine pewvint(accuracy,nr1,nr2,lf,cs,
     &                   nkc,nk,dk,converge)
      use pealloc
      implicit none
c
      integer*4 nr1,nr2,lf,nkc,nk
      real*8 accuracy,dk
      complex*16 cs
      logical*2 converge
c
c     disp: 1=uz, 2=ur,
c           3=ezz, 4=err, 5=ett, 6=ezr, 7=-dur/dz
c           8=p, 9=vz, 10=vr
c
      integer*4 i,id,ir,ik,izr,nnk
      real*8 k,k0,fac,kdk,kfdk
      complex*16 ck,cdisk,cdk,c2dk,cdk2,ckf,ckf2
      complex*16 cfac
      logical*2 finish
c
c     constants
c
      real*8 pi,pi2
      complex*16 c2
      data pi,pi2/3.14159265358979d0,6.28318530717959d0/
      data c2/(2.d0,0.d0)/
c
c     initialization
c
      do izr=1,nzr
        do ir=nr1,nr2
          grnuz(lf,ir,izr)=(0.d0,0.d0)
          grnur(lf,ir,izr)=(0.d0,0.d0)
          grnezz(lf,ir,izr)=(0.d0,0.d0)
          grness(lf,ir,izr)=(0.d0,0.d0)
          grnerr(lf,ir,izr)=(0.d0,0.d0)
          grnett(lf,ir,izr)=(0.d0,0.d0)
          grnezr(lf,ir,izr)=(0.d0,0.d0)
          grntlt(lf,ir,izr)=(0.d0,0.d0)
          grnpp(lf,ir,izr)=(0.d0,0.d0)
          grndvz(lf,ir,izr)=(0.d0,0.d0)
          grndvr(lf,ir,izr)=(0.d0,0.d0)
        enddo
      enddo
c
      k0=dmax1(dble(nd)*100.d0*dk,10.d0*pi2/(sradius+r(nr2)))
c
      cm2(1)=(0.d0,0.d0)
      cm2(2)=(1.d0,0.d0)
      cm2(3)=(0.d0,0.d0)
      cm2(4)=(0.d0,0.d0)
      cm2(5)=(0.d0,0.d0)
      cm2(6)=(1.d0,0.d0)
      cm2(7)=(1.d0,0.d0)
      cm2(8)=(0.d0,0.d0)
      cm2(9)=(0.d0,0.d0)
      cm2(10)=(1.d0,0.d0)
      do izr=1,nzr
        do id=0,nd
          do i=1,10
            dym(i,id,izr)=(0.d0,0.d0)
            dy0(i,id,izr)=(0.d0,0.d0)
            dyp(i,id,izr)=(0.d0,0.d0)
          enddo
        enddo
        do ir=nr1,nr2
          wpreuz(ir,izr)=(0.d0,0.d0)
          wpreur(ir,izr)=(0.d0,0.d0)
          wpreezz(ir,izr)=(0.d0,0.d0)
          wpreerr(ir,izr)=(0.d0,0.d0)
          wpreett(ir,izr)=(0.d0,0.d0)
          wpreess(ir,izr)=(0.d0,0.d0)
          wpreezr(ir,izr)=(0.d0,0.d0)
          wpretlt(ir,izr)=(0.d0,0.d0)
          wprepp(ir,izr)=(0.d0,0.d0)
          wpredvz(ir,izr)=(0.d0,0.d0)
          wpredvr(ir,izr)=(0.d0,0.d0)
c
          uz0(ir,izr)=(0.d0,0.d0)
          ur0(ir,izr)=(0.d0,0.d0)
          ezz0(ir,izr)=(0.d0,0.d0)
          ess0(ir,izr)=(0.d0,0.d0)
          err0(ir,izr)=(0.d0,0.d0)
          ett0(ir,izr)=(0.d0,0.d0)
          ezr0(ir,izr)=(0.d0,0.d0)
          tlt0(ir,izr)=(0.d0,0.d0)
          pp0(ir,izr)=(0.d0,0.d0)
          dvz0(ir,izr)=(0.d0,0.d0)
          dvr0(ir,izr)=(0.d0,0.d0)
        enddo
      enddo
c
      nnk=nkc
      cdk=dcmplx(dk,0.d0)
      c2dk=dcmplx(2.d0*dk,0.d0)
      cdk2=dcmplx(dk*dk,0.d0)
      do ik=1,nbsj
        k=dble(ik)*dk
        kdk=k*dk
        ck=dcmplx(k,0.d0)
        cdisk=dcmplx(disk(ik),0.d0)
c
        call pekern(cs,k)
c
        do izr=1,nzr
c
c         frequency-wavenumber spectra:
c
c         cyobs(1) = uz
c         cyobs(2) = ur
c         cyobs(3) = ezz
c         cyobs(4) = ess
c         cyobs(5) = ett (will be derived later)
c         cyobs(6) = ezr
c         cyobs(7) = duz/dr
c         cyobs(8) = p
c         cyobs(9) = dp/dz
c         cyobs(10) = dp/dr
c
          cyobs(1,izr)=cuz(izr)*cdisk
          cyobs(2,izr)=cur(izr)*cdisk
          cyobs(3,izr)=cezz(izr)*cdisk
          cyobs(4,izr)=cess(izr)*cdisk
          cyobs(5,izr)=(0.d0,0.d0)
          cyobs(6,izr)=cezr(izr)*cdisk
          cyobs(7,izr)=ck*cuz(izr)*cdisk
          cyobs(8,izr)=cp(izr)*cdisk
          cyobs(9,izr)=cdpz(izr)*cdisk
          cyobs(10,izr)=ck*cp(izr)*cdisk
        enddo
c
        do izr=1,nzr
          do i=1,10
            dym(i,0,izr)=dy0(i,0,izr)
            dy0(i,0,izr)=dyp(i,0,izr)
          enddo
        enddo
c
        if(ik.le.nd)then
          do izr=1,nzr
            do i=1,10
              dyp(i,0,izr)=(0.d0,0.d0)
            enddo
          enddo
        else
          cfac=dcmplx(dexp(-((k-dble(nd)*dk)/k0)**2),0.d0)
          do izr=1,nzr
            do i=1,10
              dyp(i,0,izr)=cyobs(i,izr)*((1.d0,0.d0)-cfac)
              cyobs(i,izr)=cyobs(i,izr)*cfac
            enddo
          enddo
        endif
c
        do ir=nr1,nr2
          cbs(0)=dcmplx( bsj(ik,0,ir)*kdk,0.d0)
          cbs(1)=dcmplx(-bsj(ik,1,ir)*kdk,0.d0)
c
          do izr=1,nzr
            uz0(ir,izr)=uz0(ir,izr)+cyobs(1,izr)*cbs(0)
            ur0(ir,izr)=ur0(ir,izr)+cyobs(2,izr)*cbs(1)
            ezz0(ir,izr)=ezz0(ir,izr)+cyobs(3,izr)*cbs(0)
            ess0(ir,izr)=ess0(ir,izr)+cyobs(4,izr)*cbs(0)
            ezr0(ir,izr)=ezr0(ir,izr)+cyobs(6,izr)*cbs(1)
            tlt0(ir,izr)=tlt0(ir,izr)+cyobs(7,izr)*cbs(1)
            pp0(ir,izr)=pp0(ir,izr)+cyobs(8,izr)*cbs(0)
            dvz0(ir,izr)=dvz0(ir,izr)+cyobs(9,izr)*cbs(0)
            dvr0(ir,izr)=dvr0(ir,izr)+cyobs(10,izr)*cbs(1)
          enddo
        enddo
c
        do id=1,min0(ik-1,nd)
          ckf=dcmplx(k-dble(id)*dk,0.d0)
          ckf2=ckf*ckf
          do izr=1,nzr
            do i=1,10
              dym(i,id,izr)=dy0(i,id,izr)
              dy0(i,id,izr)=dyp(i,id,izr)
              dyp(i,id,izr)=dy0(i,id-1,izr)*(zrs2(izr)+cm2(i)/ckf2)
     &         -(dyp(i,id-1,izr)-dym(i,id-1,izr))/c2dk/ckf
     &         -(dyp(i,id-1,izr)+dym(i,id-1,izr)
     &          -c2*dy0(i,id-1,izr))/cdk2
            enddo
          enddo
        enddo
c
        if(ik.gt.nd.and.ik.lt.nbsj-nd)then
          kfdk=(k-dble(nd)*dk)*dk
          do ir=nr1,nr2
            cbs(0)=dcmplx( bsj(ik-nd,0,ir)*kfdk,0.d0)
            cbs(1)=dcmplx(-bsj(ik-nd,1,ir)*kfdk,0.d0)
            do izr=1,nzr
c
              grnuz(lf,ir,izr)=grnuz(lf,ir,izr)
     &                        +dyp(1,nd,izr)*cbs(0)
              grnur(lf,ir,izr)=grnur(lf,ir,izr)
     &                        +dyp(2,nd,izr)*cbs(1)
              grnezz(lf,ir,izr)=grnezz(lf,ir,izr)
     &                         +dyp(3,nd,izr)*cbs(0)
              grness(lf,ir,izr)=grness(lf,ir,izr)
     &                         +dyp(4,nd,izr)*cbs(0)
              grnezr(lf,ir,izr)=grnezr(lf,ir,izr)
     &                         +dyp(6,nd,izr)*cbs(1)
              grntlt(lf,ir,izr)=grntlt(lf,ir,izr)
     &                         +dyp(7,nd,izr)*cbs(1)
              grnpp(lf,ir,izr)=grnpp(lf,ir,izr)
     &                        +dyp(8,nd,izr)*cbs(0)
              grndvz(lf,ir,izr)=grndvz(lf,ir,izr)
     &                         +dyp(9,nd,izr)*cbs(0)
              grndvr(lf,ir,izr)=grndvr(lf,ir,izr)
     &                         +dyp(10,nd,izr)*cbs(1)
            enddo
          enddo
        endif
c
        if(ik.eq.nnk)then
          nnk=2*nnk
          finish=.true.
c
          do izr=1,nzr
            do ir=nr1,nr2
              finish=finish.and.
     &         dsqrt(cdabs(grnuz(lf,ir,izr)-wpreuz(ir,izr))**2
     &              +cdabs(grnur(lf,ir,izr)-wpreur(ir,izr))**2)
     &        .le.accuracy*
     &         dsqrt(cdabs(grnuz(lf,ir,izr))**2
     &              +cdabs(grnur(lf,ir,izr))**2)
c
              finish=finish.and.
     &         dsqrt(cdabs(grnezz(lf,ir,izr)-wpreezz(ir,izr))**2
     &              +cdabs(grness(lf,ir,izr)-wpreess(ir,izr))**2
     &              +cdabs(grnezr(lf,ir,izr)-wpreezr(ir,izr))**2
     &              +cdabs(grntlt(lf,ir,izr)-wpretlt(ir,izr))**2)
     &        .le.accuracy*
     &         dsqrt(cdabs(grnezz(lf,ir,izr))**2
     &              +cdabs(grness(lf,ir,izr))**2
     &              +cdabs(grnezr(lf,ir,izr))**2
     &              +cdabs(grntlt(lf,ir,izr))**2)
c
              if(lzr(izr).ne.1.or.isurfcon.ne.1)then
                finish=finish.and.
     &          cdabs(grnpp(lf,ir,izr)-wprepp(ir,izr)).le.
     &          accuracy*cdabs(grnpp(lf,ir,izr))
              endif
c
              finish=finish.and.
     &         dsqrt(cdabs(grndvz(lf,ir,izr)-wpredvz(ir,izr))**2
     &              +cdabs(grndvr(lf,ir,izr)-wpredvr(ir,izr))**2)
     &        .le.accuracy*
     &         dsqrt(cdabs(grndvz(lf,ir,izr))**2
     &              +cdabs(grndvr(lf,ir,izr))**2)
            enddo
          enddo
c
          if(finish)goto 200
c
          do izr=1,nzr
            do ir=nr1,nr2
              wpreuz(ir,izr)=grnuz(lf,ir,izr)
              wpreur(ir,izr)=grnur(lf,ir,izr)
              wpreezz(ir,izr)=grnezz(lf,ir,izr)
              wpreess(ir,izr)=grness(lf,ir,izr)
              wpreezr(ir,izr)=grnezr(lf,ir,izr)
              wpretlt(ir,izr)=grntlt(lf,ir,izr)
              wprepp(ir,izr)=grnpp(lf,ir,izr)
              wpredvz(ir,izr)=grndvz(lf,ir,izr)
              wpredvr(ir,izr)=grndvr(lf,ir,izr)
            enddo
          enddo
        endif
      enddo
      ik=ik-1
200   nk=ik
      converge=nk.le.nbsj
c
      do ir=nr1,nr2
        do izr=1,nzr
          cfac=dcmplx(1.d0/(zrs2(izr)+r(ir)**2)**nd,0.d0)
          grnuz(lf,ir,izr)=uz0(ir,izr)+grnuz(lf,ir,izr)*cfac
          grnur(lf,ir,izr)=ur0(ir,izr)+grnur(lf,ir,izr)*cfac
          grnezz(lf,ir,izr)=ezz0(ir,izr)+grnezz(lf,ir,izr)*cfac
          grness(lf,ir,izr)=ess0(ir,izr)+grness(lf,ir,izr)*cfac
          grnezr(lf,ir,izr)=ezr0(ir,izr)+grnezr(lf,ir,izr)*cfac
          grntlt(lf,ir,izr)=tlt0(ir,izr)+grntlt(lf,ir,izr)*cfac
          grnpp(lf,ir,izr)=pp0(ir,izr)+grnpp(lf,ir,izr)*cfac
          grndvz(lf,ir,izr)=dvz0(ir,izr)+grndvz(lf,ir,izr)*cfac
          grndvr(lf,ir,izr)=dvr0(ir,izr)+grndvr(lf,ir,izr)*cfac
        enddo
      enddo
c
      do izr=1,nzr
        do ir=nr1,nr2
          if(r(ir).le.0.d0)then
c
c           ett = ess/2 
c
            grnett(lf,ir,izr)=(0.5d0,0.d0)*grness(lf,ir,izr)
c 
c           err = ess/2
c
            grnerr(lf,ir,izr)=(0.5d0,0.d0)*grness(lf,ir,izr)
          else
c
c           ett = ur/r 
c
            grnett(lf,ir,izr)=grnur(lf,ir,izr)/dcmplx(r(ir),0.d0)
c 
c           err = ess - ett
c
            grnerr(lf,ir,izr)=grness(lf,ir,izr)-grnett(lf,ir,izr)
          endif
c
c         Tilt = -dur/dz
c
          grntlt(lf,ir,izr)=grntlt(lf,ir,izr)
     &                     -(2.d0,0.d0)*grnezr(lf,ir,izr)
        enddo
      enddo
c
c     end of total wavenumber integral
c
      write(*,'(i6,a,E12.5,a,i7)')lf,'.',dimag(cs)/pi2,
     &     'Hz: wavenumbers used = ',nk
      return
      end