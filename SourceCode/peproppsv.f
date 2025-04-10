      subroutine peproppsv(l1,l2,cs,k)
      use pealloc
      implicit none
c
c     propagation of p-sv vectors
c
      integer*4 l1,l2
      real*8 k
      complex*16 cs
c
c     work space
c
      integer*4 i,j,l,ly,n
      real*8 h0
      complex*16 cdet,wave(3)
      complex*16 c0(6,3),c1(6,3),y1(6,3),orth(3,3)
c
      if(l1.eq.l2)then
        return
      else if(l1.lt.l2)then
        do l=l1+1,l2
          h0=hp(l-1)
          n=nno(l-1)
          wave(1)=dcmplx(dexp(-k*h0),0.d0)
          wave(2)=wave(1)
          wave(3)=cdexp(-ck3(n)*dcmplx(h0,0.d0))
c
          call caxcb(maiup(1,1,l-1),yuplw(1,1,l-1),6,6,3,c0)
c
c         orthonormalization of the p-sv modes
c
          cdet=c0(1,1)*c0(3,2)*c0(5,3)
     &        +c0(3,1)*c0(5,2)*c0(1,3)
     &        +c0(5,1)*c0(1,2)*c0(3,3)
     &        -c0(5,1)*c0(3,2)*c0(1,3)
     &        -c0(3,1)*c0(1,2)*c0(5,3)
     &        -c0(1,1)*c0(5,2)*c0(3,3)
          orth(1,1)=(c0(3,2)*c0(5,3)-c0(3,3)*c0(5,2))/cdet
          orth(2,1)=(c0(3,3)*c0(5,1)-c0(3,1)*c0(5,3))/cdet
          orth(3,1)=(c0(3,1)*c0(5,2)-c0(3,2)*c0(5,1))/cdet
          orth(1,2)=(c0(1,3)*c0(5,2)-c0(1,2)*c0(5,3))/cdet
          orth(2,2)=(c0(1,1)*c0(5,3)-c0(1,3)*c0(5,1))/cdet
          orth(3,2)=(c0(1,2)*c0(5,1)-c0(1,1)*c0(5,2))/cdet
          orth(1,3)=(c0(1,2)*c0(3,3)-c0(1,3)*c0(3,2))/cdet
          orth(2,3)=(c0(1,3)*c0(3,1)-c0(1,1)*c0(3,3))/cdet
          orth(3,3)=(c0(1,1)*c0(3,2)-c0(1,2)*c0(3,1))/cdet
          call caxcb(c0,orth,6,3,3,c1)
          do ly=l1,l-1
c
c           orthonormalization of the receiver vectors
c
            call caxcb(yuplw(1,1,ly),orth,6,3,3,y1)
            call cmemcpy(y1,yuplw(1,1,ly),18)
            do j=1,3
              do i=1,6
                yuplw(i,j,ly)=yuplw(i,j,ly)*wave(j)
              enddo
            enddo
          enddo
c
          c1(1,1)=(1.d0,0.d0)
          c1(2,1)=c1(2,1)*wave(1)*wave(1)
          c1(3,1)=(0.d0,0.d0)
          c1(4,1)=c1(4,1)*wave(1)*wave(2)
          c1(5,1)=(0.d0,0.d0)
          c1(6,1)=c1(6,1)*wave(1)*wave(3)
c
          c1(1,2)=(0.d0,0.d0)
          c1(2,2)=c1(2,2)*wave(2)*wave(1)
          c1(3,2)=(1.d0,0.d0)
          c1(4,2)=c1(4,2)*wave(2)*wave(2)
          c1(5,2)=(0.d0,0.d0)
          c1(6,2)=c1(6,2)*wave(2)*wave(3)
c
          c1(1,3)=(0.d0,0.d0)
          c1(2,3)=c1(2,3)*wave(3)*wave(1)
          c1(3,3)=(0.d0,0.d0)
          c1(4,3)=c1(4,3)*wave(3)*wave(2)
          c1(5,3)=(1.d0,0.d0)
          c1(6,3)=c1(6,3)*wave(3)*wave(3)
c
          call caxcb(maup(1,1,l-1),c1,6,6,3,yuplw(1,1,l))
        enddo
      else
        do l=l1-1,l2,-1
          h0=hp(l)
          n=nno(l)
c
c         determination of propagation matrix
c
          wave(1)=dcmplx(dexp(-k*h0),0.d0)
          wave(2)=wave(1)
          wave(3)=cdexp(-ck3(n)*dcmplx(h0,0.d0))
c
          call caxcb(mailw(1,1,l),yuplw(1,1,l+1),6,6,3,c0)
c
c         orthonormalization of the p-sv modes
c
          cdet=c0(2,1)*c0(4,2)*c0(6,3)
     &        +c0(4,1)*c0(6,2)*c0(2,3)
     &        +c0(6,1)*c0(2,2)*c0(4,3)
     &        -c0(6,1)*c0(4,2)*c0(2,3)
     &        -c0(4,1)*c0(2,2)*c0(6,3)
     &        -c0(2,1)*c0(6,2)*c0(4,3)
          orth(1,1)=(c0(4,2)*c0(6,3)-c0(4,3)*c0(6,2))/cdet
          orth(2,1)=(c0(4,3)*c0(6,1)-c0(4,1)*c0(6,3))/cdet
          orth(3,1)=(c0(4,1)*c0(6,2)-c0(4,2)*c0(6,1))/cdet
          orth(1,2)=(c0(2,3)*c0(6,2)-c0(2,2)*c0(6,3))/cdet
          orth(2,2)=(c0(2,1)*c0(6,3)-c0(2,3)*c0(6,1))/cdet
          orth(3,2)=(c0(2,2)*c0(6,1)-c0(2,1)*c0(6,2))/cdet
          orth(1,3)=(c0(2,2)*c0(4,3)-c0(2,3)*c0(4,2))/cdet
          orth(2,3)=(c0(2,3)*c0(4,1)-c0(2,1)*c0(4,3))/cdet
          orth(3,3)=(c0(2,1)*c0(4,2)-c0(2,2)*c0(4,1))/cdet
          call caxcb(c0,orth,6,3,3,c1)
          do ly=l1,l+1,-1
c
c           orthonormalization of the receiver vectors
c
            call caxcb(yuplw(1,1,ly),orth,6,3,3,y1)
            call cmemcpy(y1,yuplw(1,1,ly),18)
            do j=1,3
              do i=1,6
                yuplw(i,j,ly)=yuplw(i,j,ly)*wave(j)
              enddo
            enddo
          enddo
c
          c1(1,1)=c1(1,1)*wave(1)*wave(1)
          c1(2,1)=(1.d0,0.d0)
          c1(3,1)=c1(3,1)*wave(1)*wave(2)
          c1(4,1)=(0.d0,0.d0)
          c1(5,1)=c1(5,1)*wave(1)*wave(3)
          c1(6,1)=(0.d0,0.d0)
c
          c1(1,2)=c1(1,2)*wave(2)*wave(1)
          c1(2,2)=(0.d0,0.d0)
          c1(3,2)=c1(3,2)*wave(2)*wave(2)
          c1(4,2)=(1.d0,0.d0)
          c1(5,2)=c1(5,2)*wave(2)*wave(3)
          c1(6,2)=(0.d0,0.d0)
c
          c1(1,3)=c1(1,3)*wave(3)*wave(1)
          c1(2,3)=(0.d0,0.d0)
          c1(3,3)=c1(3,3)*wave(3)*wave(2)
          c1(4,3)=(0.d0,0.d0)
          c1(5,3)=c1(5,3)*wave(3)*wave(3)
          c1(6,3)=(1.d0,0.d0)
c
          call caxcb(malw(1,1,l),c1,6,6,3,yuplw(1,1,l))
        enddo
      endif
      return
      end