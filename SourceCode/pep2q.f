      subroutine pep2q(ierr)
      use pealloc
      implicit none
      integer*4 ierr
c
c     calculate equivalent source time function
c
      integer*4 i,n,it,lf,nfg,ntg
      real*8 f,dfg,dtg,fcut,diffus,alpha,spvmax
      complex*16 cs,spv,spv0,fp2q0
c
      real*8 pi,pi2,eps
      data pi,pi2,eps/3.14159265358979d0,6.28318530717959d0,1.0d-03/
c
      n=nno(lstop)
      diffus=dm(n)*(la(n)+2.d0*mu(n))/(la(n)+2.d0*mu(n)+alf(n)*qa(n))
c
      if(istype.eq.1)then
        nfg=nfmax
        dfg=0.25d0/(twindow)
        dtg=1.d0/(dble(2*nfg)*dfg)
        ntg=2*nfg
c
        allocate(fp2q(2*nfg),stat=ierr)
        if(ierr.ne.0)stop ' Error in pep2q: fpv2q not allocated!'
        allocate(tp2q(4*nfg),stat=ierr)
        if(ierr.ne.0)stop ' Error in pep2q: tpv2q not allocated!'
        allocate(tp2qout(4*nfg),stat=ierr)
        if(ierr.ne.0)stop ' Error in pep2q: tpv2qout not allocated!'
c
        alpha=-dfg*dlog(0.1d0)
c
        do lf=1,nfg
          f=dble(lf-1)*dfg
          cs=dcmplx(alpha,pi2*f)
          fp2q(lf)=dcmplx(-alf(n),0.d0)
        enddo
      else
        do i=1,6
          sfct(i)=(0.d0,0.d0)
        enddo
        sfct(6)=dcmplx(-1.d0/pi2/slength,0.d0)
c
        dfg=dmax1(0.25d0/twindow,eps*eps*diffus/(slength+sradius)**2)
        alpha=-dfg*dlog(0.1d0)
c
        nfg=2
        spv0=(0.d0,0.d0)
        spvmax=0.d0
100     fcut=dble(nfg-1)*dfg
        cs=dcmplx(pi2*fcut,0.d0)
        call pepvol(cs,spv)
        spv=(1.d0,0.d0)/spv
        spvmax=dmax1(spvmax,cdabs(spv))
        if(cdabs(spv-spv0).gt.eps*spvmax.and.nfg.lt.4*nfmax)then
          spv0=spv
          nfg=2*nfg
          goto 100
        endif
        fp2q0=spv
c
        dtg=1.d0/(dble(2*nfg)*dfg)
        ntg=2*nfg
c
        allocate(fp2q(2*nfg),stat=ierr)
        if(ierr.ne.0)stop ' Error in pep2q: fpv2q not allocated!'
        allocate(tp2q(4*nfg),stat=ierr)
        if(ierr.ne.0)stop ' Error in pep2q: tpv2q not allocated!'
        allocate(tp2qout(4*nfg),stat=ierr)
        if(ierr.ne.0)stop ' Error in pep2q: tpv2qout not allocated!'
c
        do lf=1,nfg
          f=dble(lf-1)*dfg
          cs=dcmplx(alpha,pi2*f)
          call pepvol(cs,spv)
          if(istype.eq.2)then
            fp2q(lf)=spv
          else if(istype.eq.3)then
            fp2q(lf)=(1.d0,0.d0)/spv
          else
            fp2q(lf)=(0.d0,0.d0)
          endif
        enddo
      endif
c
      call pesinterp(dtg)
c
      call subfft(fp2q,tp2q,nfg,nfg,dfg,dtg,alpha,
     &            ts,sinj,nts,tp2qout,nt,dt,fp2q0,intdif)
c
      open(30,file=filep2q,status='unknown')
      if(istype.eq.2)then
        write(30,'(a)')'     Time[s]      Injection_moment[Nm]'
      else
        write(30,'(a)')'     Time[s]      Injection_Rate[m^3/s]'
      endif
      do it=1,nt
        write(30,'(3E16.8)')dble(it-1)*dt,tp2qout(it)
      enddo
      close(30)
c
      return
      end