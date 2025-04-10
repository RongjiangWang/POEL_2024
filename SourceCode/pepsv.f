      subroutine pepsv(cs,k)
      use pealloc
      implicit none
c
c     calculation of response to p-sv source
c     ypsv(6): solution vector (complex)
c     k: wave number
c     cs: complex Laplace variable
c
      real*8 k
      complex*16 cs
c
c     work space
c
      integer*4 i,j,l,ls,n,key
      complex*16 ck
      complex*16 ma0(6,6),dmai(6,6),c0(6,3)
      complex*16 coef(6,6),b(6),coefl(12,12),bl(12)
c
      ck=dcmplx(k,0.d0)
c
c===============================================================================
c
      do l=1,lsbtm-1
        call pematrix(maup(1,1,l),cs,k,hp(l),nno(l))
        call pematinv(maiup(1,1,l),cs,k,nno(l))
      enddo
      do l=lstop,lp-1
        call pematrix(malw(1,1,l),cs,k,-hp(l),nno(l))
        call pematinv(mailw(1,1,l),cs,k,nno(l))
      enddo
c
c     matrix propagation from surface to source
c
      do j=1,3
        do i=1,6
          c0(i,j)=(0.d0,0.d0)
          yup(i,j,1)=(0.d0,0.d0)
        enddo
      enddo
c
      if(isurfcon.eq.0)then
c
c       without free surface: full-space theory
c
        c0(1,1)=(1.d0,0.d0)
        c0(3,2)=(1.d0,0.d0)
        c0(5,3)=(1.d0,0.d0)
        n=nno(1)
        call pematrix(ma0,cs,k,0.d0,n)
        call caxcb(ma0,c0,6,6,3,yup(1,1,1))
      else if(isurfcon.eq.1)then
c
c       free surface: p = 0
c
        yup(1,1,1)=(1.d0,0.d0)
        yup(3,2,1)=(1.d0,0.d0)
        yup(6,3,1)=(1.d0,0.d0)
      else if(isurfcon.eq.2)then
c
c       confined surface surface: dp/dz = 0
c
        n=nno(1)
        yup(1,1,1)=(1.d0,0.d0)
        yup(6,1,1)=dcmplx(alf(n),0.d0)*cs*yup(1,1,1)
        yup(3,2,1)=(1.d0,0.d0)
        yup(5,3,1)=(1.d0,0.d0)
      else
        print *,'error in pepsv: undefined surface conditions!'
        return
      endif
c
      do i=1,6
        do j=1,3
          yuplw(i,j,1)=yup(i,j,1)
        enddo
      enddo
      call peproppsv(1,lstop,cs,k)
      do l=1,lstop
        do i=1,6
          do j=1,3
            yup(i,j,l)=yuplw(i,j,l)
          enddo
        enddo
      enddo
c
c===============================================================================
c
c     matrix propagation from half-space to source
c
      do j=1,3
        do i=1,6
          c0(i,j)=(0.d0,0.d0)
        enddo
      enddo
c
c     coefficient vectors in the half-space
c
      n=nno(lp)
c
      c0(2,1)=(1.d0,0.d0)
      c0(4,2)=(1.d0,0.d0)
      c0(6,3)=(1.d0,0.d0)
      call pematrix(ma0,cs,k,0.d0,n)
      call caxcb(ma0,c0,6,6,3,ylw(1,1,lp))
c
      do i=1,6
        do j=1,3
          yuplw(i,j,lp)=ylw(i,j,lp)
        enddo
      enddo
      call peproppsv(lp,lsbtm,cs,k)
      do l=lp,lsbtm,-1
        do i=1,6
          do j=1,3
            ylw(i,j,l)=yuplw(i,j,l)
          enddo
        enddo
      enddo
c
c===============================================================================
c
c     conditions on the source surface
c
c
c     source function
c
      if(disksource)then
        ls=lstop
        do i=1,6
          b(i)=sfct(i)
        enddo
c
        do i=1,6
          do j=1,3
            coef(i,j)=yup(i,j,ls)
            coef(i,j+3)=-ylw(i,j,ls)
          enddo
        enddo
        key=0
        call cdsvd500(coef,b,6,1,1.d-30,key)
        if(key.eq.0)then
          print *,' Warning in pepsv: anormal exit from cdsvd500!'
          return
        endif
        do l=1,ls-1
          do i=1,6
            ypsv(i,l)=(0.d0,0.d0)
            do j=1,3
              ypsv(i,l)=ypsv(i,l)+b(j)*yup(i,j,l)
            enddo
          enddo
        enddo
        do l=ls+1,lp
          do i=1,6
            ypsv(i,l)=(0.d0,0.d0)
            do j=1,3
              ypsv(i,l)=ypsv(i,l)+b(j+3)*ylw(i,j,l)
            enddo
          enddo
        enddo
        do i=1,6
          ypsv(i,ls)=(0.d0,0.d0)
          do j=1,3
            ypsv(i,ls)=ypsv(i,ls)+(0.5d0,0.d0)
     &            *(b(j)*yup(i,j,ls)+b(j+3)*ylw(i,j,ls))
          enddo
        enddo
      else
c
c       vertical line source
c
        do j=1,3
          do i=1,6
            ysup(i,j,lstop)=yup(i,j,lstop)
            yuplw(i,j,lstop)=yup(i,j,lstop)
          enddo
        enddo
        call peproppsv(lstop,lsbtm,cs,k)
        do l=lstop,lsbtm
          do i=1,6
            do j=1,3
              ysup(i,j,l)=yuplw(i,j,l)
            enddo
          enddo
        enddo
c
        do j=1,3
          do i=1,6
            yslw(i,j,lsbtm)=ylw(i,j,lsbtm)
            yuplw(i,j,lsbtm)=ylw(i,j,lsbtm)
          enddo
        enddo
        call peproppsv(lsbtm,lstop,cs,k)
        do l=lsbtm,lstop,-1
          do i=1,6
            do j=1,3
              yslw(i,j,l)=yuplw(i,j,l)
            enddo
          enddo
        enddo
c
        call pedifmai(dmai,cs,k,nno(lstop))
        call caxcb(dmai,sfct,6,6,1,b)
c
        do i=1,6
c
          bl(i  )=b(i)
          bl(i+6)=b(i)
c
          do j=1,3
            coefl(i,j  )= yup(i,j,lstop)
            coefl(i,j+3)=-ysup(i,j,lstop)
            coefl(i,j+6)=-yslw(i,j,lstop)
            coefl(i,j+9)= (0.d0,0.d0)
          enddo
c
          do j=1,3
            coefl(i+6,j  )= (0.d0,0.d0)
            coefl(i+6,j+3)=-ysup(i,j,lsbtm)
            coefl(i+6,j+6)=-yslw(i,j,lsbtm)
            coefl(i+6,j+9)= ylw(i,j,lsbtm)
          enddo
        enddo
        key=0
        call cdsvd500(coefl,bl,12,1,1.d-30,key)
        if(key.eq.0)then
          print *,' Warning in pepsv: anormal exit from cdsvd500!'
          return
        endif
        do l=1,lstop
          do i=1,6
            ypsv(i,l)=(0.d0,0.d0)
            do j=1,3
              ypsv(i,l)=ypsv(i,l)+bl(j)*yup(i,j,l)
            enddo
          enddo
        enddo
        do l=lsbtm,lp
          do i=1,6
            ypsv(i,l)=(0.d0,0.d0)
            do j=1,3
              ypsv(i,l)=ypsv(i,l)+bl(j+9)*ylw(i,j,l)
            enddo
          enddo
        enddo
        do l=lstop+1,lsbtm-1
          do i=1,6
            ypsv(i,l)=b(i)
            do j=1,3
              ypsv(i,l)=ypsv(i,l)+bl(j+3)*ysup(i,j,l)
     &                           +bl(j+6)*yslw(i,j,l)
            enddo
          enddo
        enddo
      endif
c
c     convert poroelastic normal tress to matrix stress
c     and y6 = -xi*dp/dz + alpha*s*y1 to y6 =-xi*dp/dz
c
      do l=1,lp
        n=nno(l)
        ypsv(2,l)=ypsv(2,l)+dcmplx(alf(n),0.d0)*ypsv(5,l)
        ypsv(6,l)=ypsv(6,l)-dcmplx(alf(n),0.d0)*cs*ypsv(1,l)
      enddo
c
      return
      end