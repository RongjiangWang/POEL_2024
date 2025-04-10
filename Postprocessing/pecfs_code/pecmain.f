	program pecmain
	implicit none
c
      integer*4 inp,unit
      real*8 lat0,lon0,latr,lonr,mu,nu
      real*8 friction,ezz,err,ett,ezr,ert,etz
	real*8 pp,cfsrec,sigrec,cfsopt,sigopt
      real*8 strike,dip,rake,snn,see,sdd,sne,sed,sdn
      real*8 stkopt1,dipopt1,rakopt1,stkopt2,dipopt2,rakopt2
      real*8 cfsrecopt,rakerecopt
      character*180 inpfile
c
      print *,'*******************************************************'
      print *,'*Injection-induced Coulomb-Failure-Stress (CFS) change*'
      print *,'*******************************************************'
      print *,' '
      write(*,'(a,$)')' Select manual input (0) or via input file (1): '
      read(*,*)inp
      if(inp.eq.0)then
        print *,'=================== Input ============================'
        print *,' '
        write(*,'(a,$)')' Location of injection site (lat,lon)[deg]: '
        read(*,*)lat0,lon0
        write(*,'(a,$)')' ... and receiver fault (lat,lon)[deg]: '
        read(*,*)latr,lonr
        write(*,'(a,$)')' Receiver site shear modulus (mu)[Pa]: '
        read(*,*)mu
        write(*,'(a,$)')' ... and drained Poisson ratio (nu): '
        read(*,*)nu
        write(*,'(a,$)')' Fault mechanism (strike,dip,rake)[deg]: '
        read(*,*)strike,dip,rake
        write(*,'(a,$)')' ... and friction coefficient (friction): '
        read(*,*)friction
        write(*,'(a,$)')' Injection-induced pore pressure (Pp)[Pa]: '
        read(*,*)pp
        write(*,'(a)')' ... and strain tensor in cylindrical'
        write(*,'(a)')' coordinate (z/r/t = down/radial/azimuth)'
        write(*,'(a,$)')' (Ezz,Err,Ett,Ezr,Ert,Etz)[Pa]: '
        read(*,*)ezz,err,ett,ezr,ert,etz
      else
        write(*,'(a,$)')' Input file name: '
        read(*,'(a)')inpfile
        unit=10
        open(unit,file=inpfile,status='old')
c
c       read location of injection site (lat,lon)[deg]
c
        call skipdoc(unit)
        read(unit,*)lat0,lon0
c
c       read location of receiver fault (lat,lon)[deg]
c
        call skipdoc(unit)
        read(unit,*)latr,lonr
c
c       read receiver-site shear modulus (mu)[Pa]
c
        call skipdoc(unit)
        read(unit,*)mu
c
c       read receiver-site drained Poisson ratio (nu)
c
        call skipdoc(unit)
        read(unit,*)nu
c
c       read receiver fault mechanism (strike,dip,rake)[deg]
c
        call skipdoc(unit)
        read(unit,*)strike,dip,rake
c
c       read friction coefficient (friction)
c
        call skipdoc(unit)
        read(unit,*)friction
c
c       read injection-induced pore pressure (Pp)[Pa]
c
        call skipdoc(unit)
        read(unit,*)pp
c
c       read strain tensor in cylindrical coordinate system
c       (Ezz,Err,Ett,Ezr,Ert,Etz)[Pa]
c
        call skipdoc(unit)
        read(unit,*)ezz,err,ett,ezr,ert,etz
c
        close(unit)
      endif
c
      call calcfs(lat0,lon0,latr,lonr,mu,nu,
     &            ezz,err,ett,ezr,ert,etz,
     &            pp,friction,strike,dip,rake,
     &            snn,see,sdd,sne,sed,sdn,
     &            cfsrec,cfsrecopt,sigrec,cfsopt,sigopt,
     &            stkopt1,dipopt1,rakopt1,stkopt2,dipopt2,rakopt2,
     &            rakerecopt)
c
      print *,' '
      print *,'=================== Output ============================'
      print *,' '
c
      write(*,'(a)')' Induced (drained) stress tensor at receiver site '
      write(*,'(a)')'         (n/e/d = north/east/down) '
      write(*,'(a)')'         Snn         See         Sdd'//
     &              '         Sne         Sed         Sdn[Pa]'
      write(*,'(6E12.4)')snn,see,sdd,sne,sed,sdn
c
      print *,' '
      write(*,'(a)')' Coulomb & normal stress change'//
     &              ' along given & optimal rake on receiver fault'
      write(*,'(a)')'     CFS_Rec     SIG_Rec CFS_Rec_Opt[Pa]'
     &            //' RAKE_Rec_Opt[deg]'
      write(*,'(3E12.4,f17.4)')cfsrec,sigrec,cfsrecopt,rakerecopt
c
      print *,' '
      write(*,'(a)')' Coulomb & normal stress change'//
     &              ' on optimal-oriented fault in 3D space'
      write(*,'(a)')'  CFS_3D_Opt  SIG_3D_Opt[Pa]'
      write(*,'(2E12.4)')cfsopt,sigopt
c
      print *,' '
      write(*,'(a)')' Strike, dip & rake of'//
     &              ' optimal-oriented fault in 3D space'
      write(*,'(a)')'          Strike         Dip        Rake[deg]'
      write(*,'(a,3f12.4)')' 1. ',stkopt1,dipopt1,rakopt1
      write(*,'(a,3f12.4)')' 2. ',stkopt2,dipopt2,rakopt2
c
	stop
      end