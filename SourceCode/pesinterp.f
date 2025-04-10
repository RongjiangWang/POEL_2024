      subroutine pesinterp(dtg)
      use pealloc
      implicit none
      real*8 dtg
c
      integer*4 its,jts,its0,dnts,ierr
      real*8 ddtg,dsinj,tswin
c
      nts=0
      tswin=0.d0
      do its0=2,nts0
        if(sinj0(its0).eq.sinj0(its0-1))then
c
c         ignored if no change in injection
c
        else if(ts0(its0).eq.ts0(its0-1))then		
          nts=nts+1
          tswin=ts0(its0)
        else if(ts0(its0).gt.ts0(its0-1))then
          dnts=1+idint((ts0(its0)-ts0(its0-1))/dmin1(dt,dtg))
          ddtg=(ts0(its0)-ts0(its0-1))/dble(dnts)
          do its=1,dnts
            nts=nts+1
            tswin=tswin+dble(its)*ddtg
          enddo
        else
          stop ' Error in sinterp: wrong injection series!'
        endif
      enddo
c
      if(ninterp.gt.0)then
        deallocate(ts,sinj)
      endif
c
      allocate(ts(nts),stat=ierr)
      if(ierr.ne.0)stop 'Error in pesinterp: ts not allocated!'
      allocate(sinj(nts),stat=ierr)
      if(ierr.ne.0)stop 'Error in pesinterp: sinj not allocated!'
      ninterp=ninterp+1
c
      its=0
c
      do its0=2,nts0
        if(sinj0(its0).eq.sinj0(its0-1))then
c
c         ignored if no change in injection
c
        else if(ts0(its0).eq.ts0(its0-1))then		
          its=its+1
          if(its.gt.nts)then
            stop ' Unkonwn error #1 in sinterp!'
          endif
          ts(its)=ts0(its0)
          sinj(its)=sinj0(its0)-sinj0(its0-1)
        else if(ts0(its0).gt.ts0(its0-1))then
          dnts=1+idint((ts0(its0)-ts0(its0-1))/dmin1(dt,dtg))
          ddtg=(ts0(its0)-ts0(its0-1))/dble(dnts)
          dsinj=(sinj0(its0)-sinj0(its0-1))/dble(dnts)
          do jts=1,dnts
            its=its+1
            if(its.gt.nts)then
              stop ' Unkonwn error #2 in sinterp!'
            endif
            ts(its)=ts0(its0-1)+dble(jts)*ddtg
            sinj(its)=dsinj
          enddo
        else
          stop ' Error in sinterp: wrong injection series!'
        endif
      enddo
c
      if(its.lt.nts)then
        stop ' Unkonwn error #3 in sinterp!'
      endif
c
	return
      end