      subroutine pelayer(ierr)
      use pealloc
      implicit none
      integer ierr
c
      integer i,l,n,li,lp0
      real*8 zswap
c
      real*8,allocatable:: z(:),z0(:)
c
      real*8 pi
      data pi/3.14159265358979d0/
c
      nzmax=n0+nzr+2
      allocate(z(nzmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in pelayer: z not allocated!'
      allocate(z0(nzmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in pelayer: z0 not allocated!'
      allocate(hp(nzmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in pelayer: hp not allocated!'
      allocate(nno(nzmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in pelayer: nno not allocated!'
      allocate(zrs2(nzmax),stat=ierr)
      if(ierr.ne.0)stop 'Error in pelayer: zrs2 not allocated!'
c
      lp0=1
      z0(lp0)=0.d0
      do n=1,n0-1
        lp0=lp0+1
        if(lp0.gt.nzmax)then
          stop ' Error in pelayer: nzmax defined too small!'
        endif
        z0(lp0)=z0(lp0-1)+h(n)
      enddo
      do i=1,nzr
        lp0=lp0+1
        if(lp0.gt.nzmax)then
          stop ' Error in pelayer: nzmax defined too small!'
        endif
        z0(lp0)=zr(i)
      enddo
      lp0=lp0+1
      if(lp0.gt.nzmax)then
        stop ' Error in pelayer: nzmax defined too small!'
      endif
      z0(lp0)=zstop
      if(zsbtm.gt.zstop)then
        lp0=lp0+1
        if(lp0.gt.nzmax)then
          stop ' Error in pelayer: nzmax defined too small!'
        endif
        z0(lp0)=zsbtm
      endif
c
c     sort the z0-profile
c
      do l=1,lp0-1
        do li=l+1,lp0
          if(z0(li).lt.z0(l))then
            zswap=z0(l)
            z0(l)=z0(li)
            z0(li)=zswap
          endif
        enddo
      enddo
c
c     delete duplicates
c
      lp=1
      z(lp)=0.d0
      do l=2,lp0
        if(z0(l).gt.z(lp))then
          hp(lp)=z0(l)-z(lp)
          lp=lp+1
          z(lp)=z0(l)
        endif
      enddo
      hp(lp)=0.d0
c
c     determine lstop,lsbtm,lzr_1 to nzr
c
      do l=1,lp
        if(l.eq.1.or.l.eq.lp)then
          if(z(l).eq.zstop)lstop=l
          if(z(l).eq.zsbtm)lsbtm=l
          do i=1,nzr
            if(z(l).eq.zr(i))lzr(i)=l
          enddo
        else
          if(dabs(z(l)-zstop).le.0.25*dmin1(hp(l-1),hp(l)))lstop=l
          if(dabs(z(l)-zsbtm).le.0.25*dmin1(hp(l-1),hp(l)))lsbtm=l
          do i=1,nzr
            if(dabs(z(l)-zr(i)).le.0.25*dmin1(hp(l-1),hp(l)))lzr(i)=l
          enddo
        endif
      enddo
c
c     determine layer no of each depth
c
      li=1
      zswap=h(1)
      nno(1)=1
      do l=2,lp
        if(z(l).ge.zswap.and.li.lt.n0)then
          li=li+1
          zswap=zswap+h(li)
        endif
        nno(l)=li
      enddo
c
      deallocate(z,z0)
c
      return
      end