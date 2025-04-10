      subroutine pespectra(accuracy,nr1,nr2,ierr)
      use pealloc
      implicit none
c
      integer*4 nr1,nr2,ierr
      real*8 accuracy
c
      integer*4 i,ir,izr,l,n,na,izs,nzs,lf
      integer*4 it,ntg,nfg,nkc,nk,chi
      real*8 f,dfg,dtg,fcut,dk,dk0,am
      real*8 zs,dzs,alpha,dismax,dismin,delta
      complex*16 cs
      logical*2 converge,again
c
      real*8 pi2,eps
      data pi2,eps/6.28318530717959d0,1.0d-03/
c
      ierr=0
      print *,' ==============================================='
      write(*,'(i3,a,F12.4,a,f12.4,a)')ipr,'. sub-profile: ',
     &        r(nr1),' -> ',r(nr2),' m'
c
c     determine wavenumber sampling rate
c
      dismin=rsdis(nr1,1)
      dismax=rsdis(nr2,1)
      do izr=2,nzr
        dismin=dmin1(dismin,rsdis(nr1,izr))
        dismax=dmax1(dismax,rsdis(nr2,izr))
      enddo
c
      write(*,'(a)')'  ==============================================='
      write(*,'(a)')'  Finding required sampling rate and range'
c
      dfg=0.5d0/twindow
      alpha=-dfg*dlog(0.1d0)
c
      do izr=1,nzr
        do ir=nr1,nr2
          spreuz(ir,izr)=(0.d0,0.d0)
          spreur(ir,izr)=(0.d0,0.d0)
          spreezz(ir,izr)=(0.d0,0.d0)
          spreerr(ir,izr)=(0.d0,0.d0)
          spreett(ir,izr)=(0.d0,0.d0)
          spreezr(ir,izr)=(0.d0,0.d0)
          spretlt(ir,izr)=(0.d0,0.d0)
          sprepp(ir,izr)=(0.d0,0.d0)
          spredvz(ir,izr)=(0.d0,0.d0)
          spredvr(ir,izr)=(0.d0,0.d0)
        enddo
      enddo
c
      lf=1
      cs=dcmplx(pi2*dfg,0.d0)
c
      dk0=0.5d0*pi2/(r0+0.5d0*(dismin+dismax))
      nkc=128
      dk=dk0
c
      na=0
      call pesource(cs)
100   na=na+1
      call pebsj(nr1,nr2,dk)
      call pewvint(accuracy,nr1,nr2,lf,cs,nkc,nk,dk,converge)
      again=.false.
c
      do izr=1,nzr
        do ir=nr1,nr2
c
c         check convergence of displacement
c
          again=again.or.
     &        dsqrt(cdabs(grnuz(lf,ir,izr)-spreuz(ir,izr))**2
     &             +cdabs(grnur(lf,ir,izr)-spreur(ir,izr))**2)
     &        .gt.accuracy*dsqrt(cdabs(grnuz(lf,ir,izr))**2
     &             +cdabs(grnur(lf,ir,izr))**2)
c
c         check convergence of strain
c
          again=again.or.
     &        cdabs(grnezz(lf,ir,izr)-spreezz(ir,izr)).gt.
     &        accuracy*cdabs(grnezz(lf,ir,izr))
          again=again.or.
     &        cdabs(grnerr(lf,ir,izr)-spreerr(ir,izr)).gt.
     &        accuracy*cdabs(grnerr(lf,ir,izr))
          again=again.or.
     &      dsqrt(cdabs(grnezz(lf,ir,izr)-spreezz(ir,izr))**2
     &           +cdabs(grnerr(lf,ir,izr)-spreerr(ir,izr))**2
     &           +cdabs(grnett(lf,ir,izr)-spreett(ir,izr))**2
     &           +cdabs(grnezr(lf,ir,izr)-spreezr(ir,izr))**2
     &           +cdabs(grntlt(lf,ir,izr)-spretlt(ir,izr))**2)
     &      .gt.accuracy*
     &      dsqrt(cdabs(grnezz(lf,ir,izr))**2
     &           +cdabs(grnerr(lf,ir,izr))**2
     &           +cdabs(grnett(lf,ir,izr))**2
     &           +cdabs(grnezr(lf,ir,izr))**2
     &           +cdabs(grntlt(lf,ir,izr))**2)
c
c         check convergence of pressure
c
          if(lzr(izr).ne.1.or.isurfcon.ne.1)then
            again=again.or.
     &          cdabs(grnpp(lf,ir,izr)-sprepp(ir,izr)).gt.
     &          accuracy*cdabs(grnpp(lf,ir,izr))
          endif
c
c         check convergence of flow
c
          again=again.or.
     &        dsqrt(cdabs(grndvz(lf,ir,izr)-spredvz(ir,izr))**2
     &             +cdabs(grndvr(lf,ir,izr)-spredvr(ir,izr))**2)
     &   .gt.accuracy*
     &        dsqrt(cdabs(grndvz(lf,ir,izr))**2
     &             +cdabs(grndvr(lf,ir,izr))**2)
c
          spreuz(ir,izr)=grnuz(lf,ir,izr)
          spreur(ir,izr)=grnur(lf,ir,izr)
          spreezz(ir,izr)=grnezz(lf,ir,izr)
          spreerr(ir,izr)=grnerr(lf,ir,izr)
          spreett(ir,izr)=grnett(lf,ir,izr)
          spreezr(ir,izr)=grnezr(lf,ir,izr)
          spretlt(ir,izr)=grntlt(lf,ir,izr)
          sprepp(ir,izr)=grnpp(lf,ir,izr)
          spredvz(ir,izr)=grndvz(lf,ir,izr)
          spredvr(ir,izr)=grndvr(lf,ir,izr)
        enddo
      enddo
c
      if(again.and.nbsj.lt.nbsjmax)then
        dk=0.5d0*dk
	  goto 100
      endif
c
c     wavenumber sampling rate dk found
c
      call pebsj(nr1,nr2,dk)
c
c     determine cut-off frequency
c
      write(*,'(a)')'  ==============================================='
      write(*,'(a)')'  Finding necessary cut-off frequency'
c
      do izr=1,nzr
        do ir=nr1,nr2
          maxu(ir,izr)=dsqrt(cdabs(grnuz(lf,ir,izr))**2
     &                      +cdabs(grnur(lf,ir,izr))**2)
          maxe(ir,izr)=dsqrt(cdabs(grnezz(lf,ir,izr))**2
     &                      +cdabs(grnerr(lf,ir,izr))**2
     &                      +cdabs(grnett(lf,ir,izr))**2
     &                      +cdabs(grnezr(lf,ir,izr))**2
     &                      +cdabs(grntlt(lf,ir,izr))**2)
          maxp(ir,izr)=cdabs(grnpp(lf,ir,izr))
          maxv(ir,izr)=dsqrt(cdabs(grndvz(lf,ir,izr))**2
     &                      +cdabs(grndvr(lf,ir,izr))**2)
          spreuz(ir,izr)=grnuz(lf,ir,izr)
          spreur(ir,izr)=grnur(lf,ir,izr)
          spreezz(ir,izr)=grnezz(lf,ir,izr)
          spreerr(ir,izr)=grnerr(lf,ir,izr)
          spreett(ir,izr)=grnett(lf,ir,izr)
          spreezr(ir,izr)=grnezr(lf,ir,izr)
          spretlt(ir,izr)=grntlt(lf,ir,izr)
          spredvz(ir,izr)=grndvz(lf,ir,izr)
          spredvr(ir,izr)=grndvr(lf,ir,izr)
        enddo
      enddo
c
      nfg=4
300   fcut=dble(nfg)*dfg
      cs=dcmplx(pi2*fcut,0.d0)
      call pesource(cs)
      call pewvint(accuracy,nr1,nr2,lf,cs,nkc,nk,dk,converge)
      again=.false.
c
      do izr=1,nzr
        do ir=nr1,nr2
          am=dsqrt(cdabs(grnuz(lf,ir,izr))**2
     &            +cdabs(grnur(lf,ir,izr))**2)
          maxu(ir,izr)=dmax1(eps*maxu(ir,izr),am)
c
          am=dsqrt(cdabs(grnuz(lf,ir,izr)-spreuz(ir,izr))**2
     &            +cdabs(grnur(lf,ir,izr)-spreur(ir,izr))**2)
          again=again.or.am.gt.accuracy*maxu(ir,izr)
c
          am=dsqrt(cdabs(grnezz(lf,ir,izr))**2
     &            +cdabs(grnerr(lf,ir,izr))**2
     &            +cdabs(grnett(lf,ir,izr))**2
     &            +cdabs(grnezr(lf,ir,izr))**2
     &            +cdabs(grntlt(lf,ir,izr))**2)
          maxe(ir,izr)=dmax1(eps*maxe(ir,izr),am)
          am=dsqrt(cdabs(grnezz(lf,ir,izr)-spreezz(ir,izr))**2
     &            +cdabs(grnerr(lf,ir,izr)-spreerr(ir,izr))**2
     &            +cdabs(grnett(lf,ir,izr)-spreett(ir,izr))**2
     &            +cdabs(grnezr(lf,ir,izr)-spreezr(ir,izr))**2
     &            +cdabs(grntlt(lf,ir,izr)-spretlt(ir,izr))**2)
          again=again.or.am.gt.accuracy*maxe(ir,izr)
c
          if(lzr(izr).ne.1.or.isurfcon.ne.1)then
            am=cdabs(grnpp(lf,ir,izr))
            maxp(ir,izr)=dmax1(eps*maxp(ir,izr),am)
            am=cdabs(grnpp(lf,ir,izr)-sprepp(ir,izr))
            again=again.or.am.gt.accuracy*maxp(ir,izr)
          endif
c
          am=dsqrt(cdabs(grndvz(lf,ir,izr))**2
     &            +cdabs(grndvr(lf,ir,izr))**2)
          maxv(ir,izr)=dmax1(eps*maxv(ir,izr),am)
          am=dsqrt(cdabs(grndvz(lf,ir,izr)-spredvz(ir,izr))**2
     &            +cdabs(grndvr(lf,ir,izr)-spredvr(ir,izr))**2)
          again=again.or.am.gt.accuracy*maxv(ir,izr)
c
          spreuz(ir,izr)=grnuz(lf,ir,izr)
          spreur(ir,izr)=grnur(lf,ir,izr)
          spreezz(ir,izr)=grnezz(lf,ir,izr)
          spreerr(ir,izr)=grnerr(lf,ir,izr)
          spreett(ir,izr)=grnett(lf,ir,izr)
          spreezr(ir,izr)=grnezr(lf,ir,izr)
          spretlt(ir,izr)=grntlt(lf,ir,izr)
          sprepp(ir,izr)=grnpp(lf,ir,izr)
          spredvz(ir,izr)=grndvz(lf,ir,izr)
          spredvr(ir,izr)=grndvr(lf,ir,izr)
        enddo
      enddo
c
      if(again.and.nfg.lt.nfmin*nfmax.or.nfg.lt.nfmin)then
        lf=lf+1
        nfg=2*nfg
        goto 300
      endif
c
      if(nfg.gt.nfmax)nfg=nfmax
c
      fcut=dble(nfg-1)*dfg
      write(*,'(a)')'  ==============================================='
      write(*,'(i6,a,E12.5,a)')nfg,' frequencies used,'
     &                             //' cut-off: ',fcut,'Hz'
      write(*,'(a)')'  ==============================================='
c
      do lf=1,nfg
        f=dble(lf-1)*dfg
        cs=dcmplx(alpha,pi2*f)
        call pesource(cs)
        call pewvint(accuracy,nr1,nr2,lf,cs,nkc,nk,dk,converge)
      enddo
c
c     convert to time domain by FFT
c
      dtg=1.d0/(dble(2*nfg)*dfg)
      ntg=2*nfg
c
      call pesinterp(dtg)
c
      do izr=1,nzr
        n=nno(lzr(izr))
        chi=dm(n)*alf(n)/qa(n)
        do ir=nr1,nr2
          call subfft(grnuz(1,ir,izr),tgrn(1),nfg,nfmax,dfg,dtg,alpha,
     &           ts(1),sinj(1),nts,obs,nt,dt,spreuz(ir,izr),intdif)
          call subfft(grnur(1,ir,izr),tgrn(1),nfg,nfmax,dfg,dtg,alpha,
     &           ts(1),sinj(1),nts,obs,nt,dt,spreur(ir,izr),intdif)
          call subfft(grnezz(1,ir,izr),tgrn(1),nfg,nfmax,dfg,dtg,alpha,
     &           ts(1),sinj(1),nts,obs,nt,dt,spreezz(ir,izr),intdif)
          call subfft(grnerr(1,ir,izr),tgrn(1),nfg,nfmax,dfg,dtg,alpha,
     &           ts(1),sinj(1),nts,obs,nt,dt,spreerr(ir,izr),intdif)
          call subfft(grnett(1,ir,izr),tgrn(1),nfg,nfmax,dfg,dtg,alpha,
     &           ts(1),sinj(1),nts,obs,nt,dt,spreett(ir,izr),intdif)
          call subfft(grnezr(1,ir,izr),tgrn(1),nfg,nfmax,dfg,dtg,alpha,
     &           ts(1),sinj(1),nts,obs,nt,dt,spreezr(ir,izr),intdif)
          call subfft(grntlt(1,ir,izr),tgrn(1),nfg,nfmax,dfg,dtg,alpha,
     &           ts(1),sinj(1),nts,obs,nt,dt,spretlt(ir,izr),intdif)
          call subfft(grnpp(1,ir,izr),tgrn(1),nfg,nfmax,dfg,dtg,alpha,
     &           ts(1),sinj(1),nts,obs,nt,dt,sprepp(ir,izr),intdif)
          call subfft(grndvz(1,ir,izr),tgrn(1),nfg,nfmax,dfg,dtg,alpha,
     &           ts(1),sinj(1),nts,obs,nt,dt,spredvz(ir,izr),intdif)
          call subfft(grndvr(1,ir,izr),tgrn(1),nfg,nfmax,dfg,dtg,alpha,
     &           ts(1),sinj(1),nts,obs,nt,dt,spredvr(ir,izr),intdif)
          do it=1,(1+nt)/2
            grndvz(it,ir,izr)=-dcmplx(chi,0.d0)*grndvz(it,ir,izr)
            grndvr(it,ir,izr)=-dcmplx(chi,0.d0)*grndvr(it,ir,izr)
          enddo
        enddo
      enddo
c
      return
      end