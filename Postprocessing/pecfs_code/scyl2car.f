	subroutine scyl2car(backazi,scyl,scar)
	implicit none
c
c     Input parameters:
c          backazi = back azimuth of receiver site
c          scyl = tensor in cylindrical system
c          srr=scyl(1,1)
c          stt=scyl(2,2)
c          szz=scyl(3,3)
c          srt=scyl(1,2)
c          stz=scyl(2,3)
c          szr=scyl(3,1)
c
      real*8 backazi
      real*8 scyl(3,3)
c
c     Output paramters:
c          scar = tensor in cartesian system
c          sxx=scar(1,1)
c          syy=scar(2,2)
c          szz=scar(3,3)
c          sxy=scar(1,2)
c          syz=scar(2,3)
c          szx=scar(3,1)
c
      real*8 scar(3,3)
c
      integer*4 i,j,k
	real*8 cs,ss
      real*8 exyz(3,3),swap(3,3)
c
	double precision DEG2RAD
	data DEG2RAD/1.745329252d-02/
c
c     Direction vector of cartesian axes in cylindrical coordinates: ex,ey,ez
c     (x/y/z = north/east/down)
c
      cs=dcos(backazi*DEG2RAD)
      ss=dsin(backazi*DEG2RAD)
c
c     ex_vector
c
      exyz(1,1)=cs
      exyz(2,1)=-ss
      exyz(3,1)=0.d0
c
c     ey_vector
c
      exyz(1,2)=ss
      exyz(2,2)=cs
      exyz(3,2)=0.d0
c
c     ez_vector
c
      exyz(1,3)=0.d0
      exyz(2,3)=0.d0
      exyz(3,3)=1.d0
c
c     Catesian stress vector
c
      do i=1,3
        do j=1,3
          swap(i,j)=0.d0
          do k=1,3
            swap(i,j)=swap(i,j)+scyl(i,k)*exyz(k,j)
          enddo
        enddo
      enddo
c
      do i=1,3
        do j=1,3
          scar(i,j)=0.d0
          do k=1,3
            scar(i,j)=scar(i,j)+exyz(i,k)*swap(k,j)
          enddo
        enddo
      enddo
c
	return
      end