	subroutine calcfs(lat0,lon0,latr,lonr,mu,nu,
     &                  ezz,err,ett,ezr,ert,etz,
     &                  pp,friction,strike,dip,rake,
     &                  snn,see,sdd,sne,sed,sdn,
     &                  cfsrec,cfsrecopt,sigrec,cfsopt,sigopt,
     &                  stkopt1,dipopt1,rakopt1,stkopt2,dipopt2,rakopt2,
     &                  rakerecopt)
	implicit none
c
c     Input parameters:
c          lat0,lon0[deg] = source latitude and longitude
c          latr,lonr[deg] = receiver latitude and longitude
c          mu[Pa],nu = receiver site shear modulus and Poisson ratio
c          ezz,err,ett,ezr,ert,etz = induced strain components at receiver site
c          pp[Pa] = induced pore pressure change at receiver site
c          friction = fault friction coefficient
c          strike,dip,rake[deg] = receiver fault mechanism
c
      real*8 lat0,lon0,latr,lonr,mu,nu,
     &       ezz,err,ett,ezr,ert,etz,
     &       pp,friction,strike,dip,rake
c
c     Output paramters:
c          snn,see,sdd,sne,sed,sdn[Pa] = induced stress tensor in local Cartesian coordinate system
c                                        (n/e/d = north/east/down)
c          cfsrec[Pa] = induced Coulomb stress change on receiver fault
c          cfsrecopt[Pa] = induced optimal Coulomb stress change on receiver fault
c          sigrec[Pa] = induced normal stress change on receiver fault
c          cfsopt[Pa] = induced Coulomb stress change on optimal-oriented fault at receiver site
c          sigopt[Pa] = induced normal stress change on optimal-oriented fault at receiver site
c          tkopt1,dipopt1,rakopt1[deg] = strike, dip and rake of 1. optimal-oriented fault
c          tkopt2,dipopt2,rakopt2[deg] = strike, dip and rake of 2. optimal-oriented fault
c          rakrecopt[deg] = optimal rake on receiver fault
c
      real*8 snn,see,sdd,sne,sed,sdn,
     &       cfsrec,sigrec,cfsrecopt,cfsopt,sigopt,
     s       stkopt1,dipopt1,rakopt1,stkopt2,dipopt2,rakopt2,rakerecopt
c
      integer*4 key
      real*8 x,y,lam,azi,backazi
	real*8 scar(3,3),scyl(3,3)
c
	double precision PI,PI2,DEG2RAD,REARTH
	data PI,PI2/3.14159265358979d0,6.28318530717959d0/
	data DEG2RAD,REARTH/1.745329252d-02,6.371d+06/
c
      call disazi(1.d0,latr,lonr,lat0,lon0,x,y)
      azi=datan2(y,x)/DEG2RAD
      backazi=azi+180.d0
c
      lam=mu*nu/(0.5d0-nu)
c
c     Cylindrical coordinates: r,t,z
c
      scyl(1,1)=lam*(err+ett+ezz)+2.d0*mu*err
      scyl(1,2)=2.d0*mu*ert
      scyl(1,3)=2.d0*mu*ezr
      scyl(2,1)=2.d0*mu*ert
      scyl(2,2)=lam*(err+ett+ezz)+2.d0*mu*ett
      scyl(2,3)=2.d0*mu*etz
      scyl(3,1)=2.d0*mu*ezr
      scyl(3,2)=2.d0*mu*etz
      scyl(3,3)=lam*(err+ett+ezz)+2.d0*mu*ezz
c
      call scyl2car(backazi,scyl,scar)
c
      snn=scar(1,1)
      see=scar(2,2)
      sdd=scar(3,3)
      sne=scar(1,2)
      sed=scar(2,3)
      sdn=scar(3,1)
c
c     Coulomb stress change on receiver fault
c
      call cfsmas2d(scar(1,1),scar(2,2),scar(3,3),
     &              scar(1,2),scar(2,3),scar(3,1),
     &              pp,friction,cfsrec,cfsrecopt,sigrec,
     &              strike,dip,rake,rakerecopt)
c
c     Coulomb stress change on optimal oriented fault
c
      key=1
      call cfs3dopt(scar(1,1),scar(2,2),scar(3,3),
     &              scar(1,2),scar(2,3),scar(3,1),
     &              pp,friction,key,
     &              strike,dip,rake,
     &              cfsopt,sigopt,
     &              stkopt1,dipopt1,rakopt1,stkopt2,dipopt2,rakopt2)
	return
      end