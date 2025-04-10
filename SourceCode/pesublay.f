      subroutine pesublay(ierr)
      use pealloc
      implicit none
c
      integer*4 ierr
c
c     work space
c
      integer*4 i,l,n
      real*8 dh,dla,dmu,dalf,dqa,ddm,z,dz
c
      integer*4,allocatable:: i0(:)
c
      allocate(i0(l0),stat=ierr)
      if(ierr.ne.0)stop ' Error in pegsublay: i0 not allocated!'
c
      n0=0
      do l=1,l0-1
	  dz=z2(l)-z1(l)
        dla=2.d0*dabs(la2(l)-la1(l))/(la2(l)+la1(l))
        dmu=2.d0*dabs(mu2(l)-mu1(l))/(mu2(l)+mu1(l))
        dalf=2.d0*dabs(alf2(l)-alf1(l))/(alf2(l)+alf1(l))
        dqa=2.d0*dabs(qa2(l)-qa1(l))/(qa2(l)+qa1(l))
        ddm=2.d0*dabs(dm2(l)-dm1(l))/(dm2(l)+dm1(l))
        i0(l)=idnint(dmax1(1.d0,dla/reslm,dmu/reslm,dqa/reslm,
     &                         dalf/resla,ddm/resld))
        n0=n0+i0(l)
      enddo
c
      n0=n0+1
c
c     re-allocate h,la,mu,alf,qa,dm
c
      deallocate(h,la,mu,alf,qa,dm)
c
      allocate(h(n0),stat=ierr)
      if(ierr.ne.0)stop ' Error in pesublay: h not allocated!'
      allocate(la(n0),stat=ierr)
      if(ierr.ne.0)stop ' Error in pesublay: la not allocated!'
      allocate(mu(n0),stat=ierr)
      if(ierr.ne.0)stop ' Error in pesublay: mu not allocated!'
      allocate(alf(n0),stat=ierr)
      if(ierr.ne.0)stop ' Error in pesublay: alf not allocated!'
      allocate(qa(n0),stat=ierr)
      if(ierr.ne.0)stop ' Error in pesublay: qa not allocated!'
      allocate(dm(n0),stat=ierr)
      if(ierr.ne.0)stop ' Error in pesublay: dm not allocated!'
c
      allocate(ka(n0),stat=ierr)
      if(ierr.ne.0)stop ' Error in pesublay: ka not allocated!'
      allocate(ck3(n0),stat=ierr)
      if(ierr.ne.0)stop ' Error in pesublay: ck3 not allocated!'
      allocate(smalls(n0),stat=ierr)
      if(ierr.ne.0)stop ' Error in pesublay: smalls not allocated!'
c
      n=0
      do l=1,l0-1
        dz=z2(l)-z1(l)
        dla=(la2(l)-la1(l))/dz
        dmu=(mu2(l)-mu1(l))/dz
        dalf=(alf2(l)-alf1(l))/dz
        dqa=(qa2(l)-qa1(l))/dz
        ddm=(dm2(l)-dm1(l))/dz
        dh=dz/dble(i0(l))
        do i=1,i0(l)
          n=n+1
          h(n)=dh
          z=(dble(i)-0.5d0)*dh
          la(n)=la1(l)+dla*z
          mu(n)=mu1(l)+dmu*z
          alf(n)=alf1(l)+dalf*z
          qa(n)=qa1(l)+dqa*z
          dm(n)=dm1(l)+ddm*z
        enddo
      enddo
c
c     last layer is halfspace
c
      n=n+1
      h(n)=0.d0
      la(n)=la1(l0)
      mu(n)=mu1(l0)
      alf(n)=alf1(l0)
      qa(n)=qa1(l0)
      dm(n)=dm1(l0)
c
      write(*,'(7a)')'  no','    thick(m)    ','  la(Pa)    ',
     &    '  mu(Pa)    ','   alpha ','     Qa(Pa)',
     &    '   dm(m^2/s)  '
      do i=1,n0
        write(*,1001)i,h(i),la(i),mu(i),
     &               alf(i),qa(i),dm(i)
      enddo
1001  format(i4,6E12.4)
c
      deallocate(i0)
c
      ierr=0
      return
      end