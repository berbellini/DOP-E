c-----
c       CHANGES
c
c       02 10 2000 - added a sumav1 and a sumag1 
c           for goodness of fit information
c           This will be the sum of | residual |
c           this will be an indicator of goodness of fit
c
c       we will compute average SUM | res |
c-----
        subroutine jramat(id,nfilt,inarg,wref,invcsl,nf1,invdep,
     1      lstinv,sigr)
c-----
c     a general purpose reformatting file
c-----
c
c     list partial derivatives or computed dispersion values
c
      parameter(LER=0,LIN=5,LOT=6)
        integer NL, NL2, NLAY
        parameter(NL=200,NLAY=200,NL2=NL+NL)

c-----
c       LIN - unit for FORTRAN read from terminal
c       LOT - unit for FORTRAN write to terminal
c       LER - unit for FORTRAN error output to terminal
c       NL  - number of layers in model
c       NL2 - number of columns in model (first NL2/2 are
c           - velocity parameters, second NL2/2 are Q values)
c-----
        common/ctrl/numa,d(NL),a(NL),b(NL),r(NL),rat(NL),dd(NL2),x(NL2),
     $      h(NL2),u(NL2),c(NL2),v(NL2,NL2),qbinv(NL),qainv(NL),
     2      wc(NL2)
        logical wc

        real gaussalp
        integer nlow,nhgh,mm,mm2
        real sum2o, sum2p, sum2r,redv

        common/isomod/dl(NLAY),va(NLAY),vb(NLAY),rho(NLAY),
     1      qa(NLAY),qb(NLAY),etap(NLAY),etas(NLAY), 
     2      frefp(NLAY), frefs(NLAY)
        common/depref/refdep
        common/modtit/title
        character*80 title

        character kstnm*8
        integer nzyear, nzjday, nzhour, nzmin, nzmon, nzday

        real invwgt
c-----
c       COMMON FOR WEIGHTING
c-----
        common/datvar/numr, sumrr, sumrd, rnorm,
     1                nums, sumss, sumsd, snorm

c-----
c       get earth and q model to set up causal partials
c---- 
        call getmod(2,'tmpsrfi.17',mmax,title,iunit,iiso,iflsph,
     1      idimen,icnvel,ierr,.false.)
        m = mmax
c-----
c       get inversion control
c-----
        call ddcur(m,m2)
c-----
c       open prediction/partial file
c-----
        if(id.eq.2.or.id.eq.4)then
c-----
c           tmpsrfi.01 is for temporary storage 
c-----
            open(2,file='tmpsrfi.01',form='unformatted',
     1          access='sequential',status='unknown')
            rewind 2
c-----
c       set up file for inversion
c           tmpmrgr.9 is the input file to surfinv
c-----
            open(3,file='tmpmrgr.9',form='unformatted',
     1          access='sequential')
            rewind 3
            open(4,file='tmpsrfi.18',form='unformatted',
     1          access='sequential',status='unknown')
            rewind 4
            write(3)m2,nfilt
            write(3)(dd(i),i=1,m2)
            write(3)(wc(i),i=1,m2)
            sumv0 = 0.0
            sumv1 = 0.0
            sumv2 = 0.0
        numr = 0
        sumwt = 0.0
        sumwr = 0.0
        sumwr2 = 0.0
        sumrr = 0.0
        sumrd = 0.0
        sumrw2 = 0.0
        rnorm = 0.0
 1000   continue
        read(4,end=9000)kstnm
        read(4,end=9000)nzyear, nzjday, nzhour, nzmin, nzmon, nzday,rayp
        read(4,end=9000)nlow,nhgh,mm,mm2,iinvdep, 
     1      gaussalp, sum2o, sum2p, sum2r,redv,invwgt
            wt = abs(invwgt)
        do 1100 iii=nlow,nhgh
c-----
c       invdep = 1 invert for velocity only
c-----
            if(iinvdep.eq.1)then
C               read(4)(dd(i),i=1,m2),obs,pre,dy
                read(4)(dd(i),i=1,m),obs,pre,dy
c-----
c       force use only of dcdb
c-----
C               do 1101 ii=1,m2
C                   dd(ii) = dd(2*ii-1)
C 1101          continue
                do 1101 ii=1,m
                    dd(ii) = dd(ii)*wt
 1101           continue
                call zero(dd,m+1,m2)
                write(2)(dd(j),j=1,m2),dy*wt
                numr = numr + 1
                sumwt = sumwt + wt
                sumwr  = sumwr  + dy*wt
                sumwr2 = sumwr2 + dy*dy*wt
                sumrr = sumrr + ( dy*wt)**2
                sumrd = sumrd + (obs*wt)**2
                sumrw2 = sumrw2 + (wt)**2
c-----
c       invdep = 0 invert for layer thickness change only
c-----
            else if(iinvdep.eq.0)then
                read(4)(dd(i),i=1,m),obs,pre,dy
                do 1102 ii=1,m
                    dd(ii) = dd(ii)*wt
 1102           continue
                call zero(dd,m+1,m2)
                write(2)(dd(j),j=1,m2),dy*wt
                numr = numr + 1
                sumwt = sumwt + wt
                sumwr  = sumwr  + dy*wt
                sumwr2 = sumwr2 + dy*dy*wt
                sumrr = sumrr + ( dy*wt)**2
                sumrd = sumrd + (obs*wt)**2
                sumrw2 = sumrw2 + (wt)**2
            endif
 1100   continue
        go to 1000
 9000   continue
c-----
c       estimate standard deviation
c-----
c       if unweighted
c
c            N SUM x^2 - ( SUM X ) ^2
c       sd = ------------------------
c                     N^2
c       
c       for weighted
c
c            SUM w ( SUM w x^2 ) - ( SUM w  X ) ^2
c       sd = ------------------------
c                     (SUM w ) ^ 2
c       note this would give the proper units in terms of the
c       orignial receiver function  
c-----
c       adjust everything so that I get correct variance
c       this is done by scaling sumwt up to numr
c-----
        sdr = ( sumwt * sumwr2 - sumwr * sumwr )/(sumwt * sumwt)
        sdr = abs(sdr)
        write(LOT,'(a,f12.6)')
     1      'Receiver function fit         std err  :',
     1      sdr
c-----
c       now weight by expected std err - not this could eventually
c       be done by time point
c-----
c-----
c       now apply the priors sigr 
c-----
            if(sigr .gt. sdr)then
                sdr = sigr
            endif
c-----
c       in applying the correction to get a unit variance, we
c       actually need to use 1, w dy and ( w dy )^2 for sums
c       since we have already corrected for things but mess to use the
c       floor
c
c       or,
c       the computed sdr is from the entire data set with mixed alphas.
c       I have tried to correct for the alpha effect on anplitude by
c       dividing the partials and residuals by alpha before computing
c       sdr - Now when I compare to the user sigr the question is
c       what does sigr mean - It is the expected sigma on the
c       gauss alpha = 1.0 traces.  Note this is OK since
c       by using many traces the sigma should be on the distribution
c       of the scatter between traces and not the scatter on the mean,
c       If only one stacked function is fit then the sigma should be
c       from the variability of the stack
c-----
        fac = ( numr * sumrr - sumwr*sumwr )/(numr*numr)
        fac = sqrt(abs(fac))
        sdr = sdr/fac
        rewind (2)
 2000   continue
                read(2,end=9100)(dd(j),j=1,m2),dy
                dy = dy / sdr
                do 2100 i=1,m2
                    dd(i) = dd(i) / sdr
 2100           continue
                write(3)(dd(j),j=1,m2),dy
c-----
c       build row norm
c-----
        rowr = 0.0
        do i=1,m2
          rowr = rowr + abs(dd(j))
        enddo
        if(rowr .gt. rnorm)rnorm = rowr
        go to 2000
 9100   continue
c-----
c       close output devices
c-----
        close(1)
        close(2)
        close(3)
        close(4)
        endif
        return
        end

