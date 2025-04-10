      module pealloc
c
c     GLOBAL CONSTANTS
c     ================
c     nfmin: min. no of frequency samples for Green's functions (power of 2)
c     nfmax: max. no of frequency samples  for Green's functions (power of 2)
c     ndmax: max. order of differential transform
c     nbsjmin/max: min/max. number of Bessel function samples
c
      integer*4 nfmin,ndmax,nbsjmin,nbsjmax
      parameter(nfmin=512,ndmax=2,nbsjmin=512,nbsjmax=16384)
c
      integer*4 nbsj
      real*8 xbsjmax
      parameter(xbsjmax=1.0d+03)
c
      integer*4 nfmax
c
c     DISCRETISATION ACCURACY FOR LAYERS WITH CONSTANT GRADIENT
c     =========================================================
c     reslm: for moduls
c     resla: for alpha
c     resld: for diffusivity
c
      real*8 reslm,resla,resld
      parameter(reslm=0.05d0,resla=0.05d0,resld=0.250d0)
c
c     original model parameters
c
      integer*4 l0,lmax
      real*8,allocatable:: z1(:),z2(:),la1(:),la2(:),mu1(:),mu2(:),
     &                      alf1(:),alf2(:),qa1(:),qa2(:),dm1(:),dm2(:)
c       
c     model parameter:
c     n0: number of homogeneous layers
c
      integer*4 n0
      real*8,allocatable:: h(:),la(:),mu(:),alf(:),qa(:),dm(:)
c
      complex*16,allocatable:: ka(:),ck3(:),lpf(:)
      logical*2,allocatable:: smalls(:)
c
      integer*4 lp,nzmax
      integer*4,allocatable:: nno(:)
      real*8,allocatable:: hp(:)
c
c     zrec: receiver depth
c     lzrec: sublayer no of receiver
c
      integer*4 nzr
      integer*4,allocatable:: lzr(:)
      real*8,allocatable:: zr(:)
c
      integer*4 nr,nd
      real*8,allocatable:: r(:)
c
      integer*4 nt
      real*8 twindow,dt
c
c     source parameters
c
      integer*4 istype,intdif,lstop,lsbtm
      real*8 zstop,zsbtm,sradius,slength,r0
      complex*16 sfct(6),cbs(0:1),cm2(10)
      logical*2 disksource
c
c     table of J_n(x), n = 0, 1
c
      real*8,allocatable:: disk(:),bsj(:,:,:)
c
      real*8 rsdismax,rsdismin,zrsmin
c
      real*8,allocatable:: zrs2(:)
      real*8,allocatable:: rsdis(:,:)
c
      integer*4 isurfcon
c
c     source functions
c
      integer*4 ipr,nts0,nts
      real*8,allocatable:: ts0(:),sinj0(:)
      real*8,allocatable:: ts(:),sinj(:),tp2q(:),tp2qout(:)
      complex*16,allocatable:: fp2q(:)
c
c     psv layer matrics
c
      complex*16,allocatable:: maup(:,:,:),maiup(:,:,:)
      complex*16,allocatable:: malw(:,:,:),mailw(:,:,:)
c
c     title text
c
      character*12 txttime
      character*4 txtuz,txtur,txtpp,txtezz,txterr,
     &        txtett,txtezr,txttlt,txtdvz,txtdvr
      character*6,allocatable:: txtdep(:),txtdis(:)
c
c     memories for time series
c
      integer*4 seluz,selur,selpp,selezz,selerr,
     &        selett,selezr,seltlt,seldvz,seldvr,selp2q
      character*80 inputfile,fileuz,fileur,filepp,
     &        fileezz,fileerr,fileett,fileezr,filetlt,
     &        filedvz,filedvr,filep2q
      complex*16,allocatable:: grnuz(:,:,:),grnur(:,:,:),grnpp(:,:,:),
     &                       grnezz(:,:,:),grness(:,:,:),grnerr(:,:,:),
     &                       grnett(:,:,:),grnezr(:,:,:),grntlt(:,:,:),
     &                       grndvz(:,:,:),grndvr(:,:,:)
c
c     memories for snapshots
c
      integer*4 nsn,ninterp
      real*8,allocatable:: timesn(:)
      character*80,allocatable:: filesn(:)
c
      real*8,allocatable:: tgrn(:),obs(:)
      real*8,allocatable:: maxu(:,:),maxe(:,:)
      real*8,allocatable:: maxp(:,:),maxv(:,:)
      complex*16,allocatable:: fgrn(:)
      complex*16,allocatable:: spreuz(:,:),spreur(:,:)
      complex*16,allocatable:: spreezz(:,:),spreerr(:,:)
      complex*16,allocatable:: spreett(:,:),spreezr(:,:)
      complex*16,allocatable:: spreess(:,:)
      complex*16,allocatable:: spretlt(:,:),sprepp(:,:)
      complex*16,allocatable:: spredvz(:,:),spredvr(:,:)
      complex*16,allocatable:: wpreuz(:,:),wpreur(:,:)
      complex*16,allocatable:: wpreezz(:,:),wpreerr(:,:)
      complex*16,allocatable:: wpreett(:,:),wpreezr(:,:)
      complex*16,allocatable:: wpreess(:,:)
      complex*16,allocatable:: wpretlt(:,:),wprepp(:,:)
      complex*16,allocatable:: wpredvz(:,:),wpredvr(:,:)
      complex*16,allocatable:: uz0(:,:),ur0(:,:)
      complex*16,allocatable:: ezz0(:,:),err0(:,:)
      complex*16,allocatable:: ett0(:,:),ezr0(:,:)
      complex*16,allocatable:: ess0(:,:)
      complex*16,allocatable:: tlt0(:,:),pp0(:,:)
      complex*16,allocatable:: dvz0(:,:),dvr0(:,:)
      complex*16,allocatable:: cyobs(:,:)
      complex*16,allocatable:: dyp(:,:,:)
      complex*16,allocatable:: dy0(:,:,:)
      complex*16,allocatable:: dym(:,:,:)
      complex*16,allocatable:: ypsv(:,:),y0(:,:),ym(:,:)
      complex*16,allocatable:: ysup(:,:,:),yslw(:,:,:),
     &                         yup(:,:,:),ylw(:,:,:),yuplw(:,:,:)
      complex*16,allocatable:: cuz(:),cur(:),cp(:),cezz(:),
     &                         cess(:),cezr(:),cdpz(:)
      end module