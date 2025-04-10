      program poel2024
      use pealloc
      implicit none
c
c     work space
c
      integer*4 i,ir,izr,nr1,nr2,ierr
      real*8 accuracy,d1,d2,dbeta
c
      integer*4 runtime
      integer*4 time
c
c     read input file file
c
      print *,'######################################################'
      print *,'#                                                    #'
      print *,'#                Welcome to Program                  #'
      print *,'#                                                    #'
      print *,'#                                                    #'
      print *,'#         PPPP     OOO     EEEEE    L                #'
      print *,'#         P   P   O   O    E        L                #'
      print *,'#         PPPP    O   O    EEEE     L                #'
      print *,'#         P       O   O    E        L                #'
      print *,'#         P        OOO     EEEEE    LLLLLL           #'
      print *,'#                                                    #'
      print *,'#                                                    #'
      print *,'#                   Version-2024                     #'
      print *,'#                   ============                     #'
      print *,'#             Last modified: June 2024               #'
      print *,'#                                                    #'
      print *,'#                   * * * * *                        #'
      print *,'#         A semi-analytical code for simulating      #'
      print *,'#    fully coupled deformation-diffusion processes   #'
      print *,'#          in a layered poroelastic half-space       #'
      print *,'#          induced by an injection (pump) test       #'
      print *,'#          or an initial pore pressure anomaly       #'
      print *,'#                   * * * * *                        #'
      print *,'#                                                    #'
      print *,'#                   Copyright                        #'
      print *,'#                      by                            #'
      print *,'#                 Rongjiang Wang                     #'
      print *,'#              (wang@gfz-potsdam.de)                 #'
      print *,'#           GeoForschungsZentrum Potsdam             #'
      print *,'#                                                    #'
      print *,'######################################################'
      print *,'                                                      '
      write(*,'(a,$)')' Please type the file name of input data: '
      read(*,'(a)')inputfile
      runtime=time()
c
c     get input data and construct layered model
c
      call pegetinp(accuracy)
c
      write(*,'(a)')'  Wavenumber integration ...'
c
      d1=sradius+dmax1(rsdismin,r(1))
      d2=r(nr)
      i=1+idint(dlog(d2/d1)/dlog(5.d0))
      dbeta=dexp(dlog(d2/d1)/dble(i))
      r0=dmax1(sradius,0.1d0*rsdismin)
c
      ninterp=0
c
      ipr=1
      nr1=1
100   nr2=nr1
      do ir=nr1+1,nr
        if(r(ir).le.dmax1(r0,dbeta*r(nr1),dbeta*rsdismin))then
          nr2=ir
        endif
      enddo
c
      nd=min0(max0(0,
     &   idnint(dlog10(r(nr2)/dsqrt(r0**2+zrsmin**2)))),ndmax)
c
      if(ipr.eq.1)then
        call pegetalloc(ierr)
      endif
c
      call pespectra(accuracy,nr1,nr2,ierr)
c
      if(nr2.lt.nr)then
        ipr=ipr+1
        nr1=nr2+1
        goto 100
      endif
c
c     output time series and snapshots
c
      write(*,'(a)')'  Output ...'
      call peoutput(ierr)
c
      if(istype.ne.0)then
        call pep2q(ierr)
      endif
c
      runtime=time()-runtime
c
      write(*,'(a)')'################################################'
      write(*,'(a)')'#                                              #'
      write(*,'(a)')'#      End of computations with POEL2024       #'
      write(*,'(a)')'#                                              #'
      write(*,'(a,i10,a)')'#        Run time: ',runtime,
     &                                           ' sec              #'
      write(*,'(a)')'################################################'
      stop
      end
