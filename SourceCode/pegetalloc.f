      subroutine pegetalloc(ierr)
      use pealloc
      implicit none
c
      integer*4 ierr
c
      allocate(maup(6,6,nzmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in pegetalloc: maup not allocated!'
      allocate(maiup(6,6,nzmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in pegetalloc: maiup not allocated!'
      allocate(malw(6,6,nzmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in pegetalloc: malw not allocated!'
      allocate(mailw(6,6,nzmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in pegetalloc: mailw not allocated!'
c
      allocate(grnuz(nfmax,nr,nzr),stat=ierr)
      if(ierr.ne.0)stop ' Error in pegetalloc: grnuz not allocated!'
      allocate(grnur(nfmax,nr,nzr),stat=ierr)
      if(ierr.ne.0)stop ' Error in pegetalloc: grnur not allocated!'
      allocate(grnezz(nfmax,nr,nzr),stat=ierr)
      if(ierr.ne.0)stop ' Error in pegetalloc: grnezz not allocated!'
      allocate(grness(nfmax,nr,nzr),stat=ierr)
      if(ierr.ne.0)stop ' Error in pegetalloc: grnss not allocated!'
      allocate(grnerr(nfmax,nr,nzr),stat=ierr)
      if(ierr.ne.0)stop ' Error in pegetalloc: grnerr not allocated!'
      allocate(grnett(nfmax,nr,nzr),stat=ierr)
      if(ierr.ne.0)stop ' Error in pegetalloc: grnett not allocated!'
      allocate(grnezr(nfmax,nr,nzr),stat=ierr)
      if(ierr.ne.0)stop ' Error in pegetalloc: grnzr not allocated!'
      allocate(grnpp(nfmax,nr,nzr),stat=ierr)
      if(ierr.ne.0)stop ' Error in pegetalloc: grnpp not allocated!'
      allocate(grntlt(nfmax,nr,nzr),stat=ierr)
      if(ierr.ne.0)stop ' Error in pegetalloc: grntlt not allocated!'
      allocate(grndvz(nfmax,nr,nzr),stat=ierr)
      if(ierr.ne.0)stop ' Error in pegetalloc: grndvz not allocated!'
      allocate(grndvr(nfmax,nr,nzr),stat=ierr)
      if(ierr.ne.0)stop ' Error in pegetalloc: grnvr not allocated!'
      allocate(tgrn(2*nfmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in pegetalloc: tgrn not allocated!'
      allocate(obs(4*nfmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in pegetalloc: obs not allocated!'
      allocate(fgrn(2*nfmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in pegetalloc: fgrn not allocated!'
c
      allocate(cuz(nzr),stat=ierr)
      if(ierr.ne.0)stop ' Error in pegetalloc: cuz not allocated!'
      allocate(cur(nzr),stat=ierr)
      if(ierr.ne.0)stop ' Error in pegetalloc: cur not allocated!'
      allocate(cp(nzr),stat=ierr)
      if(ierr.ne.0)stop ' Error in pegetalloc: cp not allocated!'
      allocate(cezz(nzr),stat=ierr)
      if(ierr.ne.0)stop ' Error in pegetalloc: cezz not allocated!'
      allocate(cess(nzr),stat=ierr)
      if(ierr.ne.0)stop ' Error in pegetalloc: cess not allocated!'
      allocate(cezr(nzr),stat=ierr)
      if(ierr.ne.0)stop ' Error in pegetalloc: cezr not allocated!'
      allocate(cdpz(nzr),stat=ierr)
      if(ierr.ne.0)stop ' Error in pegetalloc: cdpz not allocated!'
c
      allocate(spreuz(nr,nzr),stat=ierr)
      if(ierr.ne.0)stop ' Error in pegetalloc: spreuz not allocated!'
      allocate(spreur(nr,nzr),stat=ierr)
      if(ierr.ne.0)stop ' Error in pegetalloc: spreur not allocated!'
      allocate(spreezz(nr,nzr),stat=ierr)
      if(ierr.ne.0)stop ' Error in pegetalloc: spreezz not allocated!'
      allocate(spreess(nr,nzr),stat=ierr)
      if(ierr.ne.0)stop ' Error in pegetalloc: spreess not allocated!'
      allocate(spreerr(nr,nzr),stat=ierr)
      if(ierr.ne.0)stop ' Error in pegetalloc: spreerr not allocated!'
      allocate(spreett(nr,nzr),stat=ierr)
      if(ierr.ne.0)stop ' Error in pegetalloc: spreett not allocated!'
      allocate(spreezr(nr,nzr),stat=ierr)
      if(ierr.ne.0)stop ' Error in pegetalloc: spreezr not allocated!'
      allocate(spretlt(nr,nzr),stat=ierr)
      if(ierr.ne.0)stop ' Error in pegetalloc: spretlt not allocated!'
      allocate(sprepp(nr,nzr),stat=ierr)
      if(ierr.ne.0)stop ' Error in pegetalloc: sprepp not allocated!'
      allocate(spredvz(nr,nzr),stat=ierr)
      if(ierr.ne.0)stop ' Error in pegetalloc: spredvz not allocated!'
      allocate(spredvr(nr,nzr),stat=ierr)
      if(ierr.ne.0)stop ' Error in pegetalloc: spredvr not allocated!'
c
      allocate(wpreuz(nr,nzr),stat=ierr)
      if(ierr.ne.0)stop ' Error in pegetalloc: wpreuz not allocated!'
      allocate(wpreur(nr,nzr),stat=ierr)
      if(ierr.ne.0)stop ' Error in pegetalloc: wpreur not allocated!'
      allocate(wpreezz(nr,nzr),stat=ierr)
      if(ierr.ne.0)stop ' Error in pegetalloc: wpreezz not allocated!'
      allocate(wpreess(nr,nzr),stat=ierr)
      if(ierr.ne.0)stop ' Error in pegetalloc: wpreess not allocated!'
      allocate(wpreerr(nr,nzr),stat=ierr)
      if(ierr.ne.0)stop ' Error in pegetalloc: wpreerr not allocated!'
      allocate(wpreett(nr,nzr),stat=ierr)
      if(ierr.ne.0)stop ' Error in pegetalloc: wpreett not allocated!'
      allocate(wpreezr(nr,nzr),stat=ierr)
      if(ierr.ne.0)stop ' Error in pegetalloc: wpreezr not allocated!'
      allocate(wpretlt(nr,nzr),stat=ierr)
      if(ierr.ne.0)stop ' Error in pegetalloc: wpretlt not allocated!'
      allocate(wprepp(nr,nzr),stat=ierr)
      if(ierr.ne.0)stop ' Error in pegetalloc: wprepp not allocated!'
      allocate(wpredvz(nr,nzr),stat=ierr)
      if(ierr.ne.0)stop ' Error in pegetalloc: wpredvz not allocated!'
      allocate(wpredvr(nr,nzr),stat=ierr)
      if(ierr.ne.0)stop ' Error in pegetalloc: wpredvr not allocated!'
c
      allocate(uz0(nr,nzr),stat=ierr)
      if(ierr.ne.0)stop ' Error in pegetalloc: uz0 not allocated!'
      allocate(ur0(nr,nzr),stat=ierr)
      if(ierr.ne.0)stop ' Error in pegetalloc: ur0 not allocated!'
      allocate(ezz0(nr,nzr),stat=ierr)
      if(ierr.ne.0)stop ' Error in pegetalloc: ezz0 not allocated!'
      allocate(ess0(nr,nzr),stat=ierr)
      if(ierr.ne.0)stop ' Error in pegetalloc: ss0 not allocated!'
      allocate(err0(nr,nzr),stat=ierr)
      if(ierr.ne.0)stop ' Error in pegetalloc: err0 not allocated!'
      allocate(ett0(nr,nzr),stat=ierr)
      if(ierr.ne.0)stop ' Error in pegetalloc: ett0 not allocated!'
      allocate(ezr0(nr,nzr),stat=ierr)
      if(ierr.ne.0)stop ' Error in pegetalloc: ezr0 not allocated!'
      allocate(tlt0(nr,nzr),stat=ierr)
      if(ierr.ne.0)stop ' Error in pegetalloc: tlt0 not allocated!'
      allocate(pp0(nr,nzr),stat=ierr)
      if(ierr.ne.0)stop ' Error in pegetalloc: pp0 not allocated!'
      allocate(dvz0(nr,nzr),stat=ierr)
      if(ierr.ne.0)stop ' Error in pegetalloc: dvz0 not allocated!'
      allocate(dvr0(nr,nzr),stat=ierr)
      if(ierr.ne.0)stop ' Error in pegetalloc: dvr0 not allocated!'
c
      allocate(cyobs(10,nzr),stat=ierr)
      if(ierr.ne.0)stop ' Error in pegetalloc: cy not allocated!'
      allocate(dyp(10,0:ndmax,nzr),stat=ierr)
      if(ierr.ne.0)stop ' Error in pegetalloc: dyp not allocated!'
      allocate(dy0(10,0:ndmax,nzr),stat=ierr)
      if(ierr.ne.0)stop ' Error in pegetalloc: dy0 not allocated!'
      allocate(dym(10,0:ndmax,nzr),stat=ierr)
      if(ierr.ne.0)stop ' Error in pegetalloc: dym not allocated!'
c
      allocate(ypsv(6,nzmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in pegetalloc: y not allocated!'
      allocate(y0(6,nzmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in pegetalloc: y0 not allocated!'
      allocate(ym(6,nzmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in pegetalloc: ym not allocated!'
      allocate(yup(6,3,nzmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in pegetalloc: yup not allocated!'
      allocate(ylw(6,3,nzmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in pegetalloc: ylw not allocated!'
      allocate(ysup(6,3,nzmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in pegetalloc: ysup not allocated!'
      allocate(yslw(6,3,nzmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in pegetalloc: yslw not allocated!'
      allocate(yuplw(6,3,nzmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in pegetalloc: yuplw not allocated!'
c
      allocate(maxu(nr,nzr),stat=ierr)
      if(ierr.ne.0)stop ' Error in pegetalloc: maxu not allocated!'
      allocate(maxe(nr,nzr),stat=ierr)
      if(ierr.ne.0)stop ' Error in pegetalloc: maxe not allocated!'
      allocate(maxp(nr,nzr),stat=ierr)
      if(ierr.ne.0)stop ' Error in pegetalloc: maxp not allocated!'
      allocate(maxv(nr,nzr),stat=ierr)
      if(ierr.ne.0)stop ' Error in pegetalloc: maxv not allocated!'
c
      return
      end