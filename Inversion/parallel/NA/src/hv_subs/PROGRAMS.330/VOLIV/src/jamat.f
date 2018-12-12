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
        subroutine jamat(id,nfilt,inarg,wref,invcsl,nf1,invdep,
     1      lstinv,pval,sigv,sigr,sigg,
     1              nurftn)
c-----
c       a general purpose reformatting file
c-----
c
c       list partial derivatives or computed dispersion values
c-----
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
        character*50 names
        character ostr*12

        integer nlow,nhgh,mm,mm2,iid
        real sum2o, sum2p, sum2r,redv

        common/isomod/dl(NLAY),va(NLAY),vb(NLAY),rho(NLAY),
     1      qa(NLAY),qb(NLAY),etap(NLAY),etas(NLAY), 
     2      frefp(NLAY), frefs(NLAY)
        common/depref/refdep
        common/modtit/title
        character*80 title

        character kstnm*8
        integer nzyear, nzjday, nzhour, nzmin, nzmon, nzday
c-----
c       COMMON FOR WEIGHTING
c-----
        common/datvar/numr, sumrr, sumrd, rnorm,
     1                nums, sumss, sumsd, snorm
c-----
c       dummy array for reading partial derivatives
c-----
        real ddd(NL2)
c-----
c       The reduction of variance parameter is defined
c       as RedVar = SUM (obs - pre)^2 / SUM obs ^2
c-----
c       Now estimate the Standard deviation
c-----
        write(6,'(a,f10.5,a,i5,a,1x,i5,a)')
     1      'Percent of Signal Power Fit (RFTN)     :',
     1      (1.0 - sumrr/sumrd)*100.0, '% for',nurftn, ' RFTNs',
     2      numr,' points'
        if(id.ne.3)then
        write(6,'(a,f10.5,a,i5,a)')
     1      'Percent of Signal Power Fit (Disp)     :',
     1      (1.0 - sumss/sumsd)*100.0, '% for',nums  , ' SW Obs'
        endif
c-----
c       To create the weighting for the joint functional,
c       note that we have already divided the residual by the
c       std err of observation
c-----
        p = pval
        facs = sqrt(    p  ) / sqrt( real(nums))
        facr = sqrt(1.0 - p) / sqrt( real(numr))
c-----
c       open the output file for srfinv96
c-----
            open(3,file='tmpsrfi.09',form='unformatted',
     1          access='sequential')
            rewind 3
        sumn = 0
        sumx = 0.0
        sumxx = 0.0
c-----
c       get the current value of dd and save it for the inversion
c-----
        call ddcur(m,m2)
c-----
c       read the surface wave information
c-----
        open(2,file='tmpmrgs.9',form='unformatted',access='sequential')
c-----
c	two pass get maximum fow sum for the partials
c-----
        rowsum = 0
        rewind 2
            read(2)m2,nfilt
            read(2)(ddd(i),i=1,m2)
            read(2)(ddd(i),i=1,m2)
 1001   continue
            read(2,end=1998) (ddd(j),j=1,m2),dy
            rows = 0.0
            do 1101 j=1,m2
                rows = rows + abs(ddd(j))
 1101       continue
            if(rows.gt.rowsum)rowsum = rows
        go to 1001
 1998   continue
        write(0,*)'DISP: facs=',facs,' rowsum=',rowsum
        if(rowsum.eq.0.0)rowsum = 1.0
        facs = facs/rowsum
        write(0,*)'    facs=',facs
c-----
c       now set up the inversion
c-----
        rewind 2
            read(2)m2,nfilt
            read(2)(ddd(i),i=1,m2)
            read(2)(ddd(i),i=1,m2)
            write(3)m2,nfilt
            write(3)(dd(i),i=1,m2)
            write(3)(wc(i),i=1,m2)
 1000   continue
            read(2,end=1999) (ddd(j),j=1,m2),dy
            do 1100 j=1,m2
                ddd(j) = ddd(j) * facs
 1100       continue
            dy = dy * facs
            write(3)(ddd(j),j=1,m2),dy

        go to 1000
 1999   continue
c-----
c       read the receiver function information
c-----
        if(id.eq.2 .or. id.eq.4)then
        open(2,file='tmpmrgr.9',form='unformatted',access='sequential')
c-----
c	two pass get maximum fow sum for the partials
c-----
        rowsum = 0
        rewind 2
            read(2)m2,nfilt
            read(2)(ddd(i),i=1,m2)
            read(2)(ddd(i),i=1,m2)
 2001   continue
            read(2,end=2998) (ddd(j),j=1,m2),dy
            rows = 0.0
            do 2101 j=1,m2
                rows = rows + abs(ddd(j))
 2101       continue
            if(rows.gt.rowsum)rowsum = rows
        go to 2001
 2998   continue
        write(0,*)'RFTN: facr=',facr,' rowsum=',rowsum
        if(rowsum.eq.0.0)rowsum = 1.0
        facr = facr/rowsum
        write(0,*)'    facr=',facr
c-----
c-----
c       now set up the inversion
c-----
        rewind 2
c-----
c           do not have to write this since we have just done it
c-----
            read(2)m2,nfilt
            read(2)(ddd(i),i=1,m2)
            read(2)(ddd(i),i=1,m2)
c-----
c       now write out the partials and residuals for the inversion
c-----
 2000   continue
            read(2,end=2999) (ddd(j),j=1,m2),dy
            do 2100 j=1,m2
                ddd(j) = ddd(j) *facr
 2100       continue
            dy = dy * facr
            write(3)(ddd(j),j=1,m2),dy
        go to 2000
 2999   continue
        close(2)
        endif
        close(3)
        return
        end

