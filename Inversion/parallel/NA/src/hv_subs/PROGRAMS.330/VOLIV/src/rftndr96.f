        program rftndr96
c----------------------------------------------------------------------c
c                                                                      c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                 c
c      VOLUME IV                                                       c
c                                                                      c
c      PROGRAM: RFTNDR96                                               c
c                                                                      c
c      COPYRIGHT 2001                                                  c
c      D. R. Russell, R. B. Herrmann                                   c
c      Department of Earth and Atmospheric Sciences                    c
c      Saint Louis University                                          c
c      221 North Grand Boulevard                                       c
c      St. Louis, Missouri 63103                                       c
c      U. S. A.                                                        c
c                                                                      c
c----------------------------------------------------------------------c
c       This program checks the input file
c       robs.d to process a list of receiver functions and outputs the
c       partial derivatives
c
c       NOTE: the program rftnpr96 has already checked the files for
c       their validity
c-----
c       CHANGES
c       23 NOV 2002 - carry on the inversion weight in tmpsrfi.15
c               to tmpsrfi.18
c       23 JUL 2009 - ensured consistent call to FFT routine with
c               complex arguments. We now call zfour
c-----
        implicit none

        integer LOT
        parameter (LOT=6)
        character fname*80

        logical dop, dotwo
        real delay, dt, rayp, gaussalp, invwgt
        integer n, iout
        character mname*80
        character kstnm*8
        integer nzyear, nzjday, nzhour, nzmin, nzmon, nzday

        integer invdep
        real twnmin,twnmax

        integer mnmarg
        logical outbin

        integer kerr, ssytem
        integer ls
        integer lgstr
        integer idtwo
        
c-----
c       call machine dependent initialization
c-----
        call mchdep()
c-----
c       parse command line arguments
c-----
c       If there are any command line arguments, program is run
c       only under command line control, else the program is
c       run under the command file control
c-----
        if(mnmarg() .gt.0)then
c-----
c           set up default earth model file
c           can be overridden by the command line
c-----
            mname = 'tmpsrfi.17'
            call gcmdln(dop,dotwo,delay,dt,rayp,gaussalp,n,
     1          mname,iout)
            twnmin = -abs(delay)
            twnmax = twnmin + (n-1)*dt
c-----
c           force incident P rftn, RFTN, double length FFT
c-----
            dop = .true.
            iout = 1
            dotwo = .true.
            outbin = .true.
c-----
c           clean up
c-----
            kerr = ssytem('rm -f rf*.sac')
            kerr = ssytem('rm -f ZZ*.sac')
            kerr = ssytem('rm -f dz*.sac')
            kerr = ssytem('rm -f PP*.sac')
            kerr = ssytem('rm -f SS*.sac')
            kerr = ssytem('rm -f dz*.sac')
            kerr = ssytem('rm -f da*.sac')
            kerr = ssytem('rm -f db*.sac')
c-----
c           compute partials with respect to layer thickness
c-----
            call doit(' ', n, dt, rayp, gaussalp, delay,1.0,0,
     1          twnmin,twnmax,dop,dotwo,mname,
     2              iout,outbin,kstnm, nzyear, nzjday, 
     3              nzhour, nzmin, nzmon, nzday)
c-----
c           compute partials with respect to layer velocities
c-----
            call doit(' ', n, dt, rayp, gaussalp, delay,1.0,1,
     1          twnmin,twnmax,dop,dotwo,mname,
     2              iout,outbin,kstnm, nzyear, nzjday, 
     3              nzhour, nzmin, nzmon, nzday)
        else
c-----
c       read in control file information concerning the variance? and
c       number of data files and the actual data window to be used
c       and whether we invert for velocity or for layer thickness
c-----
            open(9,file='rftn.tmp',status='unknown',
     1          form='formatted',access='sequential')
            rewind 9
            open(1,file='tmpsrfi.16',status='unknown',
     1          form='formatted',access='sequential')
            rewind 1
            read(1,'(i5,2f20.10,i5)')invdep,twnmin,twnmax,idtwo
            close (1)
        if(invdep.eq.0)then
            write(LOT,2)twnmin,twnmax
        else
            write(LOT,1)twnmin,twnmax
        endif
    1   format('Processing RFTN Partials for Layer Velocity ',
     1      '- Window [',f10.2,',',f10.2,']')
    2   format('Processing RFTN Partials for Layer Thickness',
     1      '- Window [',f10.2,',',f10.2,']')
        write(LOT,3)
    3   format('IncP Delay    Dt     Rayp   Gauss  Npts',
     1      ' 2x   Fit(%) Sta      Weight File')

c-----
c           open binary file for the partial derivatives
c           open tmpsrfi.18
c           output for each time sample
c               DCDA DCDB OBS PRED RES
c           or
c               DCDH OBS PRED RES
c-----
            open(4,file='tmpsrfi.18',status='unknown',
     1          form='unformatted',access='sequential')
c-----
c           force incident P rftn, RFTN, double length FFT
c-----
            dop = .true.
            iout = 1
            if(idtwo.eq.0)then
                dotwo = .false.
            else
                dotwo = .true.
            endif
            mname = 'tmpsrfi.17'
            outbin = .false.
            
c-----
c           open list of SAC RFTN files with file information
c-----

            open(1,file='tmpsrfi.15',access='sequential',
     1          form='formatted',status='unknown')
            rewind 1

 1000       continue
                read(1,'(a)',end=9000)fname
                read(1,'(a)')kstnm
                read(1,*)nzyear, nzjday, nzhour, nzmin, nzmon, nzday
                read(1,*)n, dt, rayp, gaussalp, delay, invwgt
                ls = lgstr(fname)
                call doit(fname, n, dt, rayp, gaussalp, delay,
     1              invwgt,
     1              invdep,twnmin,twnmax,dop,dotwo,mname,
     2              iout,outbin,kstnm, nzyear, nzjday, 
     3              nzhour, nzmin, nzmon, nzday)
                go to 1000
 9000       continue
            close (1)
            rewind 4
            close(4)
            close(9)
c-----
c       end of COMMAND LINE/COMMAND FILE PROCESSING
c-----
        endif
        end

        subroutine doit(fname, n, dt, rayp, gaussalp, delay,invwgt,
     1      invdep,twmin,twmax,dop,dotwo,mname,iout,outbin,
     2              kstnm, nzyear, nzjday, 
     3              nzhour, nzmin, nzmon, nzday)
        implicit none
        character fname*80
        character kstnm*8
        real rayp, delay, dt, gaussalp,twmin, twmax, invwgt
        integer n, invdep, iout
        integer nzyear, nzjday, nzhour, nzmin, nzmon, nzday
        logical dop, dotwo, outbin
        character mname*(*)
c-----
c       fname   Ch* name of SAC binary RFTN file
c       n   I   number of data points
c       dt  R   sample interval
c       rayp    R   ray parameter in sec/km
c       gaussalpR   gaussian filter parameter
c       delay   R   seconds before t=0
c       invwgt  R   weight >= to be applied to this receiver function
c       invdep  I   1 invert for velocity
c               0 invert for layer thickness
c       twnmin  R   window for receiver function fit <= 0
c       twnmax  R   window for receiver function fit > 0
c       dop L:  .true. get P response
c               .false. get response for incident S
c       dotwo   L   .true. use double length FFT
c       mname   Ch* Name of earth model file
c       iout    I   1 receiver function
c               2 vertical response
c               3 radial response
c       outbin  L   .false. output multiplexed partials on tmpsrfi.18
c-----
c       compute the partial derivatives for the current model
c       out put is on tmpsrfi.16 on unit 3
c       unit 1 is used by main routine
c       unit 2 is used by read sac
c-----
c       20 DEC 2001 - create this program based on hspec96
c-----
c       This program creates a SAC file for the 
c       surface receiver function
c
c-----
        integer LER, LOT, LIN
        parameter(LER=0, LIN=5, LOT=6)
c-----
c       earth model information
c-----
        integer NL
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        real d, a, b, rho, qa, qb, etap, etas, frefp, frefs

        integer mmax
        common/modlly/mmax
        integer iunit, iiso, iflsph, idimen, icnvel, ierr
        character title*80

        common/cntrl/ishank,hnkarg,dstcor,dokjar
        integer*4 dstcor
        real*4 hnkarg
        logical ishank,dokjar

        common/depref/refdep
        real refdep

        integer lgstr
        integer lm, lf
        
        integer nlow, nhgh
        real twnmin,twnmax
        real sum2o, sum2p, sum2r,fitv

c-----
c       time series information
c-----
        integer NSAMP, NFREQ
        parameter (NSAMP=8192,NFREQ=4097)
        real x(NSAMP)
        real y(NSAMP)
c-----
c       time series for the partials observed, predicted and residual
c-----
        real obs(NSAMP), pre(NSAMP), res(NSAMP)
        real drda(NL,NSAMP), drdb(NL,NSAMP), drdh(NL,NSAMP)


c-----
c       fref    R*4 - reference frequency for Causal Q
c       dokjar  L   - .false. use Futterman causal Q
c                 .true.  use Kjartansson Causal Q
c-----


c-----
c       
c-----
        real oldh, newh
        real olda, newa
        real oldb, newb
        integer i, j, k, n21, isign, nd2
        real df, dfac, fac
        integer nerr

        logical ext
        common/damp/alpha,ieqex
        real alpha
        integer ieqex

        real data(NSAMP), eata(NSAMP)
        real pdata(NSAMP, NL+1)
        complex Z(NL+1), Zin, Zsrc
        real freq
        integer kl

        dokjar = .false.
        if(twmax.le.0.0)then
            twnmax = (n-1)*dt - abs(delay)
        else
            twnmax = twmax
        endif
        if(twmin.gt.0.0)then
            twnmin = -abs(delay)
        else
            twnmin = twmin
        endif
        
c-----
c       get the current earth model
c-----
        lm = lgstr(mname)
        inquire(file=mname,exist=ext)
        if(.not. ext)then 
            write(LER,*)'Model file does not exist:',mname(1:lm)
            return
        endif
        call getmod(2,mname,mmax,title,iunit,iiso,iflsph,
     1          idimen,icnvel,ierr,.false.)
        if(ierr .lt. 0)return
c-----
c       make sure that we use 1/Q
c-----
        do 3007 i=1,mmax
            if(qa(i).lt.0.0)qa(i) = 0.0
            if(qb(i).lt.0.0)qb(i) = 0.0
            if(qa(i) .gt. 1.0)qa(i) = 1.0/qa(i)
            if(qb(i) .gt. 1.0)qb(i) = 1.0/qb(i)
            if(frefp(i) .le. 0.0)frefp(i) = 1.0
            if(frefs(i) .le. 0.0)frefs(i) = 1.0
 3007   continue
c-----
c       output control parameters if this is a command line run
c-----
        if(outbin)then
            write(LOT,*)'dop          :',dop
            write(LOT,*)'delay        :',delay
            write(LOT,*)'dt           :',dt   
            write(LOT,*)'rayp         :',rayp 
            write(LOT,*)'gaussalp     :',gaussalp 
            write(LOT,*)'n            :',n 
            write(LOT,*)'Double Length:',dotwo
            lm = lgstr(mname)
            write(LOT,*)'Model        :',mname(1:lm)
        endif

        if(dokjar)then
            write(LOT,*)'Kjartansson Constant Q operator used'
        endif
c-----
c       ensure that the number of points is a power of 2
c-----
        call npow2(n)
        if(dotwo)then
            n = 2*n
        endif
c-----
c       generate a zero phase pulse, with a zero at the 
c       Nyquist Frequency
c-----
        delay = abs(delay)
c-----
c       define the alpha parameter for complex frequency
c-----
        alpha = 2.5/(n*dt)
c-----
c       note this will not work for symmetric pulses with 0 lag ???
c       we remove the lag at a latter stage
c-----
        fac = 1.0
        dfac = exp(-alpha*dt)
        do 1000 i=1,n
            if(i.eq.2)then
                data(i) = (0.25/dt)*fac
            else if(i.eq.3)then
                data(i) = (0.50/dt)*fac
            else if(i.eq.4)then
                data(i) = (0.25/dt)*fac
            else
                data(i) = 0.0
            endif
            fac = fac * dfac
 1000   continue
        isign = -1
        df = 1.0/(n*dt)
        nd2 = n/2
c-----
c       eata is a temporary array used by realft
c----
        call realft(data,eata,nd2,isign,dt)
c-----
c       now do all the work
c-----
c-----
c       now process in the frequency domain
c       include the desired time shift, but not that the 
c       pulse is already
c       shifted 2*dt units 
c-----
        n21 = n / 2 + 1
        do 4000 i=1,n21
            freq=(i-1)*df
            fac = - 6.2831853*freq*(delay - 2*dt)
            if(i.eq.1)then
                Zsrc = cmplx(data(1),0.0)
            else if(i.eq.n21)then
                Zsrc = cmplx(data(2),0.0)
            else
                j = 2*i -1
                k = j + 1
                Zsrc = cmplx(data(j),data(k))
            endif
            Zsrc = Zsrc * cmplx(cos(fac),sin(fac))
            Zin = cmplx(1.0,0.0)
            call excit(freq,dop,rayp,Z,Zin,iout,invdep)
            do 4100 kl=1,mmax+1
            Z(kl) = Z(kl) * Zsrc 
            if(i.eq.1)then
                pdata(1,kl) = real(Z(kl))
            else if(i.eq.n21)then
                pdata(2,kl) = real(Z(kl))
            else
                j = 2*i -1
                k = j + 1
                pdata(j,kl) =  real(Z(kl))
                pdata(k,kl) = aimag(Z(kl))
            endif
 4100       continue
 4000   continue
c-----
c       now inverse Fourier Transform everything
c       this is not too bad since we can use Fortran array indexing
c----- 
        do 4200 kl=1,mmax+1
c-----
c           inverse Fourier transform
c-----
            isign = +1
            nd2 = n / 2
            call realft(pdata(1,kl),eata,nd2,isign,dt)
c-----
c           undamp the time series
c-----
            fac = exp(-alpha*delay)
            dfac = exp(alpha*dt) 
            do 425 i = 1,n 
                pdata(i,kl)= pdata(i,kl) * fac
                fac = fac * dfac 
  425       continue 
c-----
c           now Gaussian filter
c-----
            isign = -1
            df = 1.0/(n*dt)
            nd2 = n / 2
            call realft(pdata(1,kl),eata,nd2,isign,dt)
            do 426 i=1,n21
                freq=(i-1)*df
                fac = (6.2831853*freq)/(2.0*gaussalp)
                if(fac.gt.25.0)then
                    fac = 0.0
                else
                    fac = exp( - fac * fac)
                endif
c-----
c           multiplication since we
c           multiply a complex number by a real
c-----
                if(i.eq.1)then
                    pdata(1,kl) = pdata(1,kl) * fac
                else if(i.eq.n21)then
c-----
c               Source pulse has Zero at Nyquist
c-----
                    pdata(2,kl) = pdata(2,kl) * fac
                else
                    j = 2*i -1
                    k = j + 1
                    pdata(j,kl) = fac * pdata(j,kl)
                    pdata(k,kl) = fac * pdata(k,kl)
                endif
  426       continue
c-----
c           inverse Fourier transform
c-----
            isign = +1
            nd2 = n / 2
            call realft(pdata(1,kl),eata,nd2,isign,dt)
 4200   continue
c-----
c       get extreme values of time series for SAC header
c-----
        if(dotwo)then
            n = n / 2
        endif
c-----
c       Now pdata(1,1)      --- pdata(nd2,1) is perturbed 1st layer
c       Now pdata(1,2)      --- pdata(nd2,2) is perturbed 2nd layer
c       Now pdata(1,mmax)   --- pdata(nd2,mmax) is perturbed halfspace
c       Now pdata(1,mmax+1) --- pdata(nd2,mmax+1) is unperturbed RFTN
c-----
        call copyts(pdata(1,mmax+1),x,n)
        if(outbin)then
            call outsac('rf',0,delay,n,dt,gaussalp,rayp,x)
        else
            call brsac(3,NSAMP,fname,obs,nerr)
c-----
c           save the prediction
c-----
            call copyts(x,pre,n)
c-----
c           compute the residual = observed - predicted
c-----
            call copyts(obs,res,n)
            call dodiff(pre,res,n,1.0)
        endif
        if(invdep.eq.0)then
c-----
c           get partials with respect to layer thickness
c-----
            do 8000 i=1,mmax-1
                oldh = d(i)
                newh = 1.01 * d(i)
                call copyts(pdata(1,i),y,n)
C               if(outbin)then
C                   call outsac('ZZ',i,delay,n,dt,
C     1                 gaussalp,rayp,y)
C               endif
                call dodiff(x,y,n,newh-oldh)
c-----
c               reset
c-----
                if(outbin)then
                    call outsac('dz',i,delay,n,dt,
     1                  gaussalp,rayp,y)
                endif
                if(.not. outbin)then
                    call copyit(y,drdh,i,n,NL,NSAMP)
                endif
 8000       continue
            do 8101 i=1,n
                y(i) = 0.0
 8101       continue
                if(.not. outbin)then
                    call copyit(y,drdh,mmax,n,NL,NSAMP)
                endif
        else if(invdep.eq.1)then
c-----
c           get partials with respect to P-velocity and S-velocity
c-----
            do 8100 i=1,mmax
                oldb = b(i)
                newb = 1.01 * b(i)
                call copyts(pdata(1,i),y,n)
C               if(outbin)then
C                   call outsac('SS',i,delay,n,dt,
C     1                 gaussalp,rayp,y)
C               endif
                call dodiff(x,y,n,newb-oldb)
c-----
c               reset
c-----
                if(outbin)then
                    call outsac('db',i,delay,n,dt,
     1                  gaussalp,rayp,y)
                endif
                if(.not. outbin)then
                    call copyit(y,drdb,i,n,NL,NSAMP)
                endif
 8100       continue
C           do 1200 i=1,mmax
C               olda = a(i)
C               newa = 1.01 * a(i)
C               call copyts(pdata(1,i),y,n)
C               if(outbin)then
C                   call outsac('PP',i,delay,n,dt,
C     1                 gaussalp,rayp,y)
C               endif
C               call dodiff(x,y,n,newa-olda)
Cc-----
Cc              reset
Cc-----
C               if(outbin)then
C                   call outsac('da',i,delay,n,dt,
C     1                 gaussalp,rayp,y)
C               endif
C               if(.not. outbin)then
C                   call copyit(y,drda,i,n,NL,NSAMP)
C               endif
C 1200      continue
        endif
        if(.not. outbin)then
c-----
c           output the windowed information in the file tmpsrfi.18
c           First time sample is at time '-delay'
c           twnmin <=0, twnmax > 0
c-----
            nlow = (delay + twnmin)/dt + 1
            nhgh = (delay + twnmax)/dt + 1
            if(nlow .lt. 1)nlow = 1
            if(nhgh .gt. n)nhgh = n
c-----
c           get sum of squares
c-----
            sum2o = 0.0
            sum2p = 0.0
            sum2r = 0.0
            do 4900 i=nlow,nhgh
                sum2o = sum2o + obs(i)*obs(i)
                sum2p = sum2p + pre(i)*pre(i)
                sum2r = sum2r + res(i)*res(i)
 4900       continue
            fitv = 100.0 * ( 1.0 - sum2r/ sum2o)
c-----
c           output processing parameters and results for this run
c-----
            lf = lgstr(fname)
            write(LOT,1)dop,delay,dt,rayp,gaussalp,n,
     1          dotwo,fitv,kstnm,invwgt,fname(1:lf)
            write(9,'(f8.2,f10.2,1x,a)')gaussalp,fitv,fname(1:lf)
    1   format(L2,f8.2,f7.2,f8.3,f8.2,i6,L3,f9.4,1x,a8,1x,f6.2,1x,a)
c-----
c           output to the processing summary file used by rmodls
c-----

            write(4)kstnm
            write(4)nzyear, nzjday, nzhour, nzmin, nzmon, nzday,rayp
            write(4)nlow,nhgh,mmax,2*mmax,invdep, 
     1          gaussalp, sum2o, sum2p, sum2r,fitv,invwgt
            if(invdep.eq.1)then
                do 5000 i=nlow,nhgh
C           write(4)(drda(j,i),drdb(j,i),j=1,mmax),obs(i),pre(i),res(i)
            write(4)(drdb(j,i),j=1,mmax),obs(i),pre(i),res(i)
 5000           continue
            else if(invdep.eq.0)then
                do 5100 i=nlow,nhgh
            write(4)(drdh(j,i),j=1,mmax),obs(i),pre(i),res(i)
 5100           continue
            endif
        endif
        return
        end

        subroutine  copyit(y,arr2d,l,n,NL,NSAMP)
        implicit none
c-----
c       copy the vector y into the 2-D arr
        integer l, n, NL, NSAMP
        real arr2d(NL,NSAMP)
        real y(NSAMP)
        integer i
        do 1000 i=1,n
            arr2d(l,i) = y(i)
 1000   continue
        return
        end

        subroutine dodiff(x,y,n,dv)
c-----
c       compute the derivative of a function
c-----
        implicit none
        integer n
        real x(n), y(n), dv
        integer i
        do 1000 i=1,n
            y(i) = ( y(i) - x(i) )/dv
 1000   continue
        return
        end

        subroutine copyts(x,y,n)
c-----
c       copy array x to array y
c-----
        implicit none
        integer n
        real x(n), y(n)
        integer i
        do 1000 i=1,n
            y(i) =  x(i) 
 1000   continue
        return
        end

        subroutine outsac(cp,ii,delay,n,dt,gaussalp,rayp,x)
        character*2 cp
        integer ii, n
        real delay, dt, gaussalp, rayp
        real x(n)

        real depmin, depmax
        character ofile*9
        
        write(ofile,'(a2,i3.3,a4)')cp,ii,'.sac'
        call scmxmn(x,n,depmax,depmin,depmen,indmax,indmin)
c-----
c       output the time series as a SAC file
c-----
        call newhdr()
c-----
c       set real header value
c-----
            call setfhv('DELTA   ',dt,nerr)
            call setfhv('DEPMIN  ',depmin,nerr)
            call setfhv('DEPMAX  ',depmax,nerr)
            call setfhv('DEPMEN  ',depmen,nerr)
            beg = -delay 
            e = beg + (n-1)*dt
            call setfhv('B       ',beg   ,nerr)
            call setfhv('TIMMAX  ',beg+indmax*dt,nerr)
            call setfhv('TIMMIN  ',beg+indmin*dt,nerr)
            call setfhv('E       ',e     ,nerr)
            o = 0.0
            call setfhv('O       ',o     ,nerr)
c-----
c       set integer header value
c-----
            ksyear = 1970
            kdoy = 1
            kshour = 0
            ksmin = 0
            isec = 0
            jsmsec = 0
            call setnhv('NZYEAR  ',ksyear,nerr)
            call setnhv('NZJDAY  ',kdoy  ,nerr)
            call setnhv('NZHOUR  ',kshour,nerr)
            call setnhv('NZMIN   ',ksmin ,nerr)
            call setnhv('NZSEC   ',isec  ,nerr)
            call setnhv('NZMSEC  ',jsmsec,nerr)
            call setnhv('NPTS    ',n, nerr)
c-----
c       set logical header value
c-----
            call setlhv('LEVEN   ',.true.,nerr)
            call setlhv('LPSPOL  ',.true.,nerr)
            call setlhv('LOVROK  ',.true.,nerr)
            call setlhv('LCALDA  ',.true.,nerr)
c-----
c       set character header value
c-----
            call setkhv('KSTNM    ','RFTN    ' ,nerr)
c-----
c       set enumerated header value
c-----
            call setihv('IFTYPE   ','ITIME   ',nerr)
            call setihv('IZTYPE   ','IB      ',nerr)
c-----
c       set deconvolution specific headers
c-----
        call setfhv('USER0   ',gaussalp,nerr)
        call setfhv('USER4   ',rayp,nerr)
        call setkhv('KUSER0  ','Rftn    ',nerr)
        call setkhv('KUSER1  ','FD_DECON',nerr)
        call setkhv('KEVNM   ','Rftn    ',nerr)
        call setfhv('USER5   ',100.0,nerr)
        call setkhv('KCMPNM  ','Rftn    ',nerr)
            
c-----
c       output the time series
c-----
        call bwsac(3,n,ofile,x)
        return
        end


        subroutine gcmdln(dop,dotwo,delay,dt,rayp,gaussalp,
     1      nsamp,mname,iout)
c-----
c       parse command line arguments
c       requires subroutine mgtarg() and function mnmarg()
c-----
c       dop L   - .true. P-wave incident
c                 .false. SV-wave incident
c       dotwo   L   - .true. computations done with double length,
c                   but only single length is output
c       delay   R   - padding before the beginning of the 
c           receiver function
c       rayp    R   - ray parameter in sec/km
c       gaussalp R  - Gaussian filter parameter
c       nsamp   I   - number of samples (poer of 2)
c       mname   Ch* - earth model name - required
c       dokjar  L   - 
c       iout    I   - output 1 = rftn, 2 = z, 3 = r time series
c-----
        logical dop, dotwo
        real delay, dt, rayp, gaussalp
        integer nsamp, iout
        character mname*(*)

        common/cntrl/ishank,hnkarg,dstcor,dokjar
        logical ishank,dokjar
        integer*4 dstcor
        real*4 hnkarg



        integer*4 mnmarg
        character*50 name

        IOUT = 1
        DOKJAR = .FALSE.
        ISHANK = .FALSE.
        DSTCOR = 0
        HNKARG = 3.0
        dop = .true.
        dotwo = .false.
        delay = 5.0
        dt = 1.0
        rayp = 0.05
        gaussalp = 1.0
        nsamp = 512
C-----
c       no default for the model name - we set this externally
c----
C       mname = ' '
        nmarg=mnmarg()
        i = 0
   11   i = i + 1
            if(i.gt.nmarg)goto 13
            call mgtarg(i,name)
            if(name(1:2).eq.'-P' )then
                dop = .true.
            else if(name(1:2).eq.'-S' )then
                dop = .false.
            else if(name(1:2).eq.'-D' .and. name(1:3).ne.'-DT')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,f10.0)')delay
            else if(name(1:3).eq.'-DT')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,f10.0)')dt
            else if(name(1:2).eq.'-N')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,i10)')nsamp
            else if(name(1:5).eq.'-RAYP')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,f10.0)')rayp
            else if(name(1:3).eq.'-DT')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,f10.0)')dt
            else if(name(1:4).eq.'-ALP')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,f10.0)')gaussalp
            else if(name(1:2).eq.'-M')then
                i = i + 1
                call mgtarg(i,mname)
            else if(name(1:2).eq.'-2')then
                dotwo = .true.
            else if(name(1:2).eq.'-z')then
                iout = 2
            else if(name(1:2).eq.'-r')then
                iout = 3
            else if(name(1:2).eq.'-?')then
                call usage(' ')
            else if(name(1:2).eq.'-h')then
                call usage(' ')
            endif
        goto 11
   13   continue
        return
        end

        subroutine usage(ostr)
c------
c       write out program syntax
c-----
        character ostr*(*)
        integer LER,LOT,LIN
        parameter(LER=0, LIN=5, LOT=6)
        if(ostr.ne. ' ' )then
            lostr = lgstr(ostr)
            write(LER,*)ostr(1:lostr)
        endif
        write(LER,*)'USAGE: ',
     1  'rftndr96 [-P] [-S] [-2] [-r] [-z] -RAYP p',
     1      ' -ALP alpha -DT dt -NSAMP nsamp'
        write(LER,*)
     1  '-P           (default true )    Incident P wave'
        write(LER,*)
     1  '-S           (default false)    Incident S wave'
        write(LER,*)
     1  '-RAYP p      (default 0.05 )    Ray parameter in sec/km'
        write(LER,*)
     1  '-DT dt       (default 1.0  )    Sample interval for synthetic'
        write(LER,*)
     1  '-NSAMP nsamp (default 512  )    Number samples for synthetic'
        write(LER,*)
     1  '-M   model   (default none )    Earth model name'
        write(LER,*)
     1  '-ALP alp     (default 1.0  )    Number samples for synthetic'
        write(LER,*)
     1  '-2           (default false)    Use 2x length internally'
        write(LER,*)
     1  '-r           (default false)    Output radial   time series'
        write(LER,*)
     1  '-z           (default false)    Output vertical time series'
        write(LER,*)
     1  '     H(f) = exp( - (pi freq/alpha)**2) '
        write(LER,*)
     1  '     Filter corner ~ alpha/pi '
        write(LER,*)
     1  '     -2  (default false) use double length FFT to'
        write(LER,*)
     1  '     avoid FFT wrap around in convolution '
        write(LER,*)
     1  '-D delay     (default 5 sec)    output delay sec before t=0'
        write(LER,*)
     1  '-?                   Display this usage message'
        write(LER,*)
     1  '-h                   Display this usage message'
        stop 
        end


        subroutine aten(om,qa,qb,xka,xkb,alpha,a,b,atna,atnb,iwat,
     2      frefp,frefs)
c-----
c       make velocities complex, using Futterman causality operator
c-----
        real*4 qa,qb,alpha,a,b
        complex*16 om,at,atna,atnb,xka,xkb
        real*8 pi, om1p, om1s, oml, fac, pi2
        common/cntrl/ishank,hnkarg,dstcor,dokjar
        integer*4 dstcor
        real*4 hnkarg
        logical ishank,dokjar
        real*8 CDABS
        complex*16 CDLOG
c-----
c       reference frequency is fref hz
c-----
        om1p=6.2831853*frefp
        om1s=6.2831853*frefs
        pi2 = 1.5707963
        pi=3.1415927d+00
        if(dokjar)then
c-----
c       Kjartansson Constant Q, causal Q operator
c       Kjartansson, E. (1979). 
c           Constant Q-wave propagation and attenuation,
c       J. Geophys. Res. 84, 4737-4748.
c-----
            gama = atan(qa)/pi
            gamb = atan(qb)/pi
            if(gama.le.0.0)then
                atna = cmplx(1.0,0.0)
            else
                fac = pi2*gama
                rfac = sin(fac)/cos(fac)
                atna = dcmplx(1.0d+00,0.0d+00)/
     1              (( (om/om1p)**dble(-gama) ) *
     2               dcmplx(1.0d+00,-dble(rfac)))
            endif
            if(b.gt.1.0e-04*a)then
                if(gamb.le.0.0)then
                    atnb = cmplx(1.0,0.0)
                else
                    fac = pi2*gamb
                    rfac = sin(fac)/cos(fac)
                    atnb = dcmplx(1.0d+00,0.0d+00)/
     1              (( (om/om1s)**dble(-gamb) ) *
     2               dcmplx(1.0d+00,-dble(rfac)))
                endif
            endif
        else
c-----
c       Futterman Causal Q
c-----
c           low frequency cutoff is 0.01 hz
c-----
            oml=0.062831853d+00
            atna=dcmplx(1.0d+00,0.0d+00)
            atnb=dcmplx(1.0d+00,0.0d+00)
            if(qa.gt.0.0)then
              at=dcmplx(0.0d+00,0.0d+00)
              if(CDABS(om).gt.oml) at=CDLOG(om/om1p)/pi
              if(CDABS(om).le.oml) then
                fac=dsqrt(oml*oml + dble(alpha*alpha))/oml
              at=CDLOG(dcmplx(dble(oml),-dble(alpha))/(om1p*fac))/pi
              endif
              atna=(1.+dble(qa)*at+dcmplx(0.0d+00,dble(qa/2.)))
            endif
            if(qb.gt.0.0)then
              at=dcmplx(0.0d+00,0.0d+00)
              if(CDABS(om).gt.oml) at=CDLOG(om/om1s)/pi
              if(CDABS(om).le.oml) then
                fac=dsqrt(oml*oml + dble(alpha*alpha))/oml
              at=CDLOG(dcmplx(dble(oml),-dble(alpha))/(om1s*fac))/pi
              endif
               atnb=(1.+dble(qb)*at+dcmplx(0.0d+00,dble(qb/2.)))
            endif
        endif
        xka=om/(dble(a)*atna)
        if(b.le.1.0e-04*a)then
            iwat = 1
            xkb = dcmplx(0.0d+00,0.0d+00)
        else
            iwat = 0
            xkb=om/(dble(b)*atnb)
        endif
        return
        end

        subroutine hska(aa,w,x,y,z,cosp,cosq,wvno,wvno2,gam,gamm1,rho,
     1      iwat,ex,om2)
        implicit none
c-----
c       command line arguments
c-----
        complex*16 aa(4,4),w,x,y,z,cosp,cosq
        complex*16 wvno,wvno2,gam,gamm1
        real*4 rho
        real*8 ex 
        complex*16 om2
        integer iwat
c-----
c       internal variables
c-----
        complex*16 cpq, gcpq, zw2, gzw2, g1w, g1y, gx
        real*8 dfac
        complex*16 zrho
        integer i,j

        zrho = dcmplx(dble(rho),0.0d+00)
        if(iwat.eq.1)then
c-----
c       fluid layer
c-----
            do 100 j=1,4
                do 101 i=1,4
                    aa(i,j) = dcmplx(0.0d+00,0.0d+00)
  101           continue
  100       continue
            if(ex.gt.35.0d+00)then
                dfac = 0.0d+00
            else
                dfac = dexp(-ex)
            endif
            aa(1,1) = dfac
            aa(4,4) = dfac
            aa(2,2) = cosp
            aa(3,3) = cosp
            aa(2,3) = -x/zrho
            aa(3,2) = - zrho*w
        else
c-----
c       elastic layer
c-----
c       W = Sa/ra
c       X = ra*Sa
c       Y = Sb/rb
c       Z = rb*Sb
c-----
            cpq = cosp-cosq
            gcpq = gam*cpq
            zw2 = z/wvno2
            gzw2 = gam*zw2
            g1w = gamm1*w
            g1y = gamm1*y
            gx = gam*x
            aa(1,1)=   gcpq + cosq
                aa(1,3)= - wvno * cpq/(zrho*om2)
            aa(1,2)=   wvno*(-g1w+gzw2)
            aa(1,4)=   (wvno2*w-z)/(zrho*om2)
            aa(2,1)=   (gx - wvno2*g1y)/wvno
            aa(2,2)= - gcpq + cosp
            aa(2,3)=   (-x+wvno2*y)/(zrho*om2)
            aa(2,4)= - aa(1,3)
            aa(3,1)=   zrho*om2*gamm1*gcpq/wvno
            aa(3,2)=   zrho*om2*((-gamm1*g1w)+(gam*gzw2))
            aa(3,3)=   aa(2,2)
            aa(3,4)= - aa(1,2)
            aa(4,1)=   zrho*om2*(((gam*gx)/wvno2) - (gamm1*g1y))
            aa(4,2)= - aa(3,1)
            aa(4,3)= - aa(2,1)
            aa(4,4)=   aa(1,1)
        endif
        return
        end

        subroutine var(p,q,ra,rb,w,x,y,z,cosp,cosq,ex,
     1      exa,exl,yl,zl,cosql,iwat)
c     not modified for negative p,q
c     this assumes that real p and real q have same signs
        common/ovrflw/a0,cpcq,cpy,cpz,cqw,cqx,xy,xz,wy,wz
        complex*16 cpcq,cpy,cpz,cqw,cqx,xy,xz,wy,wz
        complex*16 p,q,ra,rb,w,x,y,z,cosp,cosq
        complex*16 yl,zl,cosql
        complex*16 eqp,eqm,epp,epm,sinp,sinq
        real *8 a0,pr,pi,qr,qi,fac,qmp,ex,exa,exl
c-----
c       form terms such as cos(p), sin(p), cos(q), sin(q)
c       and cos(p)*cos(q)
c
c       Introduce a factorization of exponentials to
c       make a pseudo floating point ssytem
c
c       ex is the exponent in cosp
c       exl is the exponent in cosq for SH
c       exa is the exponent in cosp*cosq
c-----
        real*8 DREAL
      ex=0.0d+00
      exl = 0.0d+00
      a0=0.0d+00
      pr=dreal(p)
      pi=dimag(p)
      epp=dcmplx(dcos(pi),dsin(pi))/2.
      epm=dconjg(epp)
      ex=pr
      fac=0.0
      if(pr.lt.15.) fac=dexp(-2.*pr)
      cosp=epp + fac*epm
      sinp=epp - fac*epm
      w=sinp/ra
      x=ra*sinp
        if(iwat.eq.1)then
c-----
c       fluid layer
c-----
            a0 = 1.0d+00
            exa = ex
            cosq = 1.0d+00
            y = 0.0d+00
            z = 0.0d+00
            cosql = 1.0d+00
            yl = 0.0d+00
            zl = 0.0d+00
            exl = 0.0d+00
        else
c-----
c       elastic layer
c-----
            qr=dreal(q)
            qi=dimag(q)
            eqp=dcmplx(dcos(qi),dsin(qi))/2.
            eqm=dconjg(eqp)
            exl=qr
            fac=0.0d+00
            if(qr.lt.15.) fac=dexp(-2.*qr)
            cosql=eqp + fac*eqm
            sinq=eqp - fac*eqm
            yl=sinq/rb
            zl=rb*sinq
c-----
c       form factors for compound P-SV matrix
c-----
            exa=pr + qr
            cpcq=cosp*cosql
            cpy=cosp*yl
            cpz=cosp*zl
            cqw=cosql*w
            cqx=cosql*x
            xy=x*yl
            xz=x*zl
            wy=w*yl
            wz=w*zl
            fac=0.0d+00
            qmp=qr-pr
            if(qmp.gt.-40.) fac=dexp(qmp)
            cosq=cosql*fac
            y=fac*yl
            z=fac*zl
            fac=0.0d+00
            if(exa.lt.60.) a0=dexp(-exa)
        endif
        return
        end

        subroutine npow2(npts)
c-----
c       Given npts, determine the N=2**m such that N >= npts
c       return the new ntps
c-----
        integer*4 nsamp, npts
        nsamp = npts
        npts = 1
 1000   continue
        if(npts.ge.nsamp)return
        npts = 2*npts
        go to 1000
        end

        subroutine copy(da,aa)
c-----
c       copy da array to aa
c-----
        implicit none
        complex*16 da(4,4), aa(4,4)
        integer i, j
        
        do 1000 j=1,4
            do 1100 i=1,4
                aa(i,j) = da(i,j)
 1100       continue
 1000   continue
        return
        end

        subroutine excit(freq,dop,rayp,Z,Zin,iout,iptrb)
c-----
c       Compute the medium response at the given frequency
c
c       freq    R   - frequency in Hz
c       dop L   - .true.    P - wave incidenc
c                 .false.   S - wave incident
c       rayp    P   - ray parameter in sec/km
c       Z   C   - array of (mmax+1) RFTN's 
c       Zin C   - source pulse spectrum at this frequency
c       iout    I   1 output Ur/Uz
c               2 output Uz
c               3 output Ur
c       iptrb   I   0 get partials for layer thickness
c               1 get partials for S velocity
c               2 get partials for P velocity
c
c-----
        implicit none
c-----
c       command line arguments - define
c-----
        real freq, rayp
        logical dop
        integer iout, iptrb
        complex Z(*)
        complex Zin
c-----
c       model definition
c-----
        integer NL
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
            real d, a, b, rho, qa, qb, etap, etas, frefp, frefs
        common/modlly/mmax
            integer mmax
c-----
c       matrix components in layers and boundaries saved
c-----
        common/damp/alpha,ieqex
            real alpha
            integer ieqex

c-----
c       internal variables
c       to do computations we need 1st two columns of 4x4 matrix
c-----
C       complex*16 hv(4,4,NL)
        complex*16 hv(4,2,NL)
        integer k,n
        complex*16 aa(4,4)
        complex*16 paa(4,4)
        real*8 ex(NL)
        real*8 pex(NL)
        real*8 exa, pexa
        complex*16 om, om2, wvno, wvno2
        real omega
        real*8 exl(NL)
c-----
c       initialize vectors
c       In what follows, k=1,...,mmax-1 refers to the layer k
c                        k=mmax         refers to the halfspace
c                        k=mmax+1       refers to the unperturbed 
c           solution
c       we compute the response for the oroginal and perturbed medium
c       and then compute the partials later by first order one sided
c       differences. Do not worry about correctness of 
c           partial derivative
c       estimate since the estiamte will get better as the iterative
c       RFTN inversion converges to a solution
c-----
        do 1000 k=1,mmax+1
            hv(1,1,k) = Zin
            hv(1,2,k) = cmplx(0.0,0.0)
            hv(2,1,k) = cmplx(0.0,0.0)
            hv(2,2,k) = Zin
            hv(3,1,k) = cmplx(0.0,0.0)
            hv(3,2,k) = cmplx(0.0,0.0)
            hv(4,1,k) = cmplx(0.0,0.0)
            hv(4,2,k) = cmplx(0.0,0.0)

C           Hv(1,3,k) = cmplx(0.0,0.0)
C           Hv(1,4,k) = cmplx(0.0,0.0)
C           Hv(2,3,k) = cmplx(0.0,0.0)
C           Hv(2,4,k) = cmplx(0.0,0.0)
C           Hv(3,3,k) = cmplx(0.0,0.0)
C           Hv(3,4,k) = cmplx(0.0,0.0)
C           Hv(4,3,k) = cmplx(0.0,0.0)
C           Hv(4,4,k) = cmplx(0.0,0.0)
            exl(k) = 0.0d+00
 1000   continue
c-----
c       define complex angular frequency and wavenumber
c-----
        omega = 6.2831853*freq
        om =  dcmplx(dble(omega), dble(-alpha))
        om2 = om * om
        wvno  = dcmplx(dble(rayp),0.0d+00)*om
        wvno2 = wvno * wvno
c-----
c       multiply the Haskell matrices from top down
c       taking special care for the perturbed layer
c-----
        do 2000 n=1,mmax-1
            call gthska (n, ex(n), aa,      om,om2,wvno,wvno2)
            call pgthska(n,pex(n),paa,iptrb,om,om2,wvno,wvno2)
            do 3000 k=1,mmax+1
                if(k.eq.n)then
c-----
c                   propagate using perturbed model in layer
c-----
                    call hvm(k,hv,paa)
                    exl(k) = exl(k) + pex(k)
                else
c-----
c                   propagate using original  model in layer
c-----
                    call hvm(k,hv, aa)
                    exl(k) = exl(k) + ex(k)
                endif
 3000       continue
 2000   continue
c-----
c       multiply the E sub N sup -1
c-----
            call gteni (mmax, aa,om,om2,wvno,wvno2)
            call pgteni(mmax,paa,om,om2,wvno,wvno2,iptrb)
            do 4000 k=1,mmax+1
                if(k.eq.mmax)then
c-----
c                   propagate using perturbed model in halfspace
c-----
                    call hvm(k,hv,paa)
                else
c-----
c                   propagate using original model in layer
c-----
                    call hvm(k,hv, aa)
                endif
 4000       continue
c-----
c       Now form the solutions - note that I do everything in loops
c       eventually to vectorize
c-----
        if(dop)then
            if(iout.eq.1)then
c-----
c           UR/UZ
c-----
            do 6000 k=1,mmax+1
            Z(k) = - dcmplx(0.0d+00, 1.0d+00) * hv(2,2,k) 
     1          / hv(2,1,k)
 6000       continue
            else if(iout.eq.2)then
c-----
c           UZ
c-----
            do 6010 k=1,mmax+1
            Z(k) = -  hv(2,1,k)*exp(-exl(k))/
     1          (hv(1,1,k)*hv(2,2,k) - hv(2,1,k)*hv(1,2,k))
 6010       continue
            else if(iout.eq.3)then
c-----
c           Ur
c-----
            do 6020 k=1,mmax+1
            Z(k) = dcmplx(0.0d+00, 1.0d+00) * hv(2,2,k)*exp(-exl(k))/
     1          (hv(1,1,k)*hv(2,2,k) - hv(2,1,k)*hv(1,2,k))
 6020       continue
            endif
        else
            if(iout.eq.1)then
c-----
c           UR/UZ
c-----
            do 6030 k=1,mmax+1
            Z(k) = - dcmplx(0.0d+00, 1.0d+00) * hv(1,2,k) 
     1          / hv(1,1,k)
 6030       continue
            else if(iout.eq.2)then
c-----
c           UZ
c-----
            do 6040 k=1,mmax+1
            Z(k) = dcmplx(0.0d+00, 1.0d+00) * hv(1,1,k)*exp(-exl(k))/
     1          (hv(1,1,k)*hv(2,2,k) - hv(2,1,k)*hv(1,2,k))
 6040       continue
            else if(iout.eq.3)then
c-----
c           Ur
c-----
            do 6050 k=1,mmax+1
            Z(k) = -                          hv(1,2,k)*exp(-exl(k))/
     1          (hv(1,1,k)*hv(2,2,k) - hv(2,1,k)*hv(1,2,k))
 6050       continue
            endif
        endif
            
        return
        end

        subroutine gthska(m,ex,aa,om,om2,wvno,wvno2)
        implicit none
        integer m
        real*8 ex
        complex*16 aa(4,4)
        integer iptrb
        complex*16 om, om2, wvno, wvno2

        integer NL
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
            real d, a, b, rho, qa, qb, etap, etas, frefp, frefs
        common/modlly/mmax
            integer mmax
        common/damp/alpha,ieqex
            real alpha
            integer ieqex

        complex*16 w,x,y,z,cosp,cosq,gam,gamm1
        complex*16 xka,xkb,ra,rb
        complex*16 atna, atnb,p,q
        complex*16 yl,zl,cosql
        real*8 exa, exb
        integer iwat

        call aten(om,qa(m),qb(m),xka,xkb,
     1      alpha,a(m),b(m),atna,atnb,iwat,
     2      frefp(m),frefs(m))
        gam=dble(b(m))*(wvno/om)
        gam = gam * atnb
        gam = dcmplx(2.0d+00,0.0d+00)*gam*gam
        gamm1 = gam - dcmplx(1.0d+00,0.0d+00)
        ra=CDSQRT(wvno2-xka*xka)
        rb=CDSQRT(wvno2-xkb*xkb)
        p=ra*dble(d(m))
        q=rb*dble(d(m))
        call var(p,q,ra,rb,w,x,y,z,cosp,cosq,
     1          ex,exa,exb,yl,zl,cosql,iwat)
        call hska(aa,w,x,y,z,cosp,cosq,wvno,wvno2,gam,gamm1,rho(m),
     1      iwat,ex,om2)
        return
        end

        subroutine pgthska(n,ex,aa,iptrb,om,om2,wvno,wvno2)
        implicit none
        integer n
        real*8 ex
        complex*16 aa(4,4)
        integer iptrb
        complex*16 om,om2, wvno, wvno2

        integer NL
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
            real d, a, b, rho, qa, qb, etap, etas, frefp, frefs
        common/modlly/mmax
            integer mmax
        real hsav, asav, bsav
c-----
c       perturb the model
c-----
        if(iptrb.eq.0)then
            hsav = d(n)
            d(n) = 1.01 * d(n)
        else if(iptrb.eq.1)then
            bsav = b(n)
            b(n) = 1.01 * b(n)
        else if(iptrb.eq.2)then
            asav = a(n)
            a(n) = 1.01 * a(n)
        endif
c-----
        call gthska(n,ex,aa,om,om2,wvno,wvno2)
c-----
c       un perturb the model
c-----
        if(iptrb.eq.0)then
            d(n) = hsav
        else if(iptrb.eq.1)then
            b(n) = bsav
        else if(iptrb.eq.2)then
            a(n) = asav
        endif
        return
        end

        subroutine gteni(m,g,om,om2,wvno,wvno2)
        implicit none
        integer m
        complex*16 g(4,4)
        complex*16 om, om2, wvno, wvno2
        integer NL
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
            real d, a, b, rho, qa, qb, etap, etas, frefp, frefs
        common/modlly/mmax
            integer mmax
        common/damp/alpha,ieqex
            real alpha
            integer ieqex
c-----
c       get elements of E sub N sup -1 matrix
c-----
        complex*16 xka,xkb,ra,rb,gam,gamm1
        complex*16 atna, atnb
        complex*16 CDSQRT
        integer iwat
c-----
c       set up halfspace conditions
c-----
            call aten(om,qa(m),qb(m),xka,xkb,
     1          alpha,a(m),b(m),atna,atnb,iwat,
     2          frefp(m),frefs(m))
            gam=dble(b(m))*(wvno/om)
            gam = gam * atnb
            gam = dcmplx(2.0d+00,0.0d+00)*gam*gam
            gamm1 = gam - dcmplx(1.0d+00,0.0)
            ra=CDSQRT(wvno2-xka*xka)
            rb=CDSQRT(wvno2-xkb*xkb)
            g(1,1) =   gam/wvno
            g(2,1) = - gamm1/rb
            g(3,1) =   g(1,1)
            g(4,1) =   g(2,1)
            g(1,2) =  - gamm1 / ra
            g(2,2) =   g(1,1)
            g(3,2) =    gamm1 / ra
            g(4,2) = - g(1,2)
            g(1,3) = - 1.0d+00/(rho(m)*om2)
            g(2,3) =   wvno / ( rb * rho(m)*om2)
            g(3,3) =   g(1,3)
            g(4,3) =   g(3,3)
            g(1,4) =   wvno / ( ra * rho(m)*om2)
            g(2,4) =   g(1,3)
            g(3,4) = - wvno / ( ra * rho(m)*om2)
            g(4,4) = - g(2,4)
        return
        end

        subroutine pgteni(n,aa,om,om2,wvno,wvno2,iptrb)
        implicit none
        integer n
        complex*16 aa(4,4)
        complex*16 om, om2, wvno, wvno2
        integer iptrb
        integer NL
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
            real d, a, b, rho, qa, qb, etap, etas, frefp, frefs
        common/modlly/mmax
            integer mmax
        real hsav, asav, bsav
c-----
c       matrix components in layers and boundaries saved
c-----
c       perturb the model
c       note that if iptrb = 0 for layer thickness we do nothing
c-----
        if(iptrb.eq.1)then
            bsav = b(n)
            b(n) = 1.01 * b(n)
        else if(iptrb.eq.2)then
            asav = a(n)
            a(n) = 1.01 * a(n)
        endif
        call gteni(n,aa,om,om2,wvno,wvno2)
c-----
c       un perturb the model
c-----
        if(iptrb.eq.1)then
            b(n) = bsav
        else if(iptrb.eq.2)then
            a(n) = asav
        endif
        return
        end

        subroutine hvm(n,hv,aa)
c-----
c       FORM hv = aa hv
c-----
        implicit none
        integer NL
        parameter(NL=200)
        complex*16 hv(4,2,NL)
C       complex*16 hv(4,4,NL)
        integer n
        complex*16 aa(4,4)

        complex*16 a11, a21, a31, a41
        complex*16 a12, a22, a32, a42
C       COMPLEX*16 A13, A23, A33, A43
C       COMPLEX*16 A14, A24, A34, A44

        a11 = aa(1,1)*hv(1,1,n) + aa(1,2)*hv(2,1,n) 
     1      + aa(1,3)*hv(3,1,n) + aa(1,4)*hv(4,1,n)
        a12 = aa(1,1)*hv(1,2,n) + aa(1,2)*hv(2,2,n) 
     1      + aa(1,3)*hv(3,2,n) + aa(1,4)*hv(4,2,n)
        a21 = aa(2,1)*hv(1,1,n) + aa(2,2)*hv(2,1,n) 
     1      + aa(2,3)*hv(3,1,n) + aa(2,4)*hv(4,1,n)
        a22 = aa(2,1)*hv(1,2,n) + aa(2,2)*hv(2,2,n) 
     1      + aa(2,3)*hv(3,2,n) + aa(2,4)*hv(4,2,n)
        a31 = aa(3,1)*hv(1,1,n) + aa(3,2)*hv(2,1,n) 
     1      + aa(3,3)*hv(3,1,n) + aa(3,4)*hv(4,1,n)
        a32 = aa(3,1)*hv(1,2,n) + aa(3,2)*hv(2,2,n) 
     1      + aa(3,3)*hv(3,2,n) + aa(3,4)*hv(4,2,n)
        a41 = aa(4,1)*hv(1,1,n) + aa(4,2)*hv(2,1,n) 
     1      + aa(4,3)*hv(3,1,n) + aa(4,4)*hv(4,1,n)
        a42 = aa(4,1)*hv(1,2,n) + aa(4,2)*hv(2,2,n) 
     1      + aa(4,3)*hv(3,2,n) + aa(4,4)*hv(4,2,n)

C       A13 = AA(1,1)*HV(1,3,N) + AA(1,2)*HV(2,3,N) 
C     1     + AA(1,3)*HV(3,3,N) + AA(1,4)*HV(4,3,N)
C       A14 = AA(1,1)*HV(1,2,N) + AA(1,2)*HV(4,N,N) 
C     1     + AA(1,3)*HV(3,2,N) + AA(1,4)*HV(4,2,N)
C       A23 = AA(2,1)*HV(1,3,N) + AA(2,2)*HV(2,3,N) 
C     1     + AA(2,3)*HV(3,3,N) + AA(2,4)*HV(4,3,N)
C       A24 = AA(2,1)*HV(1,4,N) + AA(2,2)*HV(2,4,N) 
C     1     + AA(2,3)*HV(3,4,N) + AA(2,4)*HV(4,4,N)
C       A33 = AA(3,1)*HV(1,3,N) + AA(3,2)*HV(2,3,N) 
C     1     + AA(3,3)*HV(3,3,N) + AA(3,4)*HV(4,3,N)
C       A34 = AA(3,1)*HV(1,4,N) + AA(3,2)*HV(2,4,N) 
C     1     + AA(3,3)*HV(3,4,N) + AA(3,4)*HV(4,4,N)
C       A43 = AA(4,1)*HV(1,3,N) + AA(4,2)*HV(2,3,N) 
C     1     + AA(4,3)*HV(3,3,N) + AA(4,4)*HV(4,3,N)
C       A44 = AA(4,1)*HV(1,4,N) + AA(4,2)*HV(2,4,N) 
C     1     + AA(4,3)*HV(3,4,N) + AA(4,4)*HV(4,4,N)
        hv(1,1,n) = a11
        hv(1,2,n) = a12
        hv(2,1,n) = a21
        hv(2,2,n) = a22
        hv(3,1,n) = a31
        hv(3,2,n) = a32
        hv(4,1,n) = a41
        hv(4,2,n) = a42

C       HV(1,3,N) = A13
C       HV(1,4,N) = A14
C       HV(2,3,N) = A23
C       HV(2,4,N) = A24
C       HV(3,3,N) = A33
C       HV(3,4,N) = A34
C       HV(4,3,N) = A43
C       HV(4,4,N) = A44
        return
        end
