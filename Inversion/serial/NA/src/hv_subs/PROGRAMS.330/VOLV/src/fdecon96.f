        program fdecon96
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME V                                                       c
c                                                                     c
c      PROGRAM: FDECON96                                              c
c                                                                     c
c      COPYRIGHT 1998 R. B. Herrmann                                  c
c                                                                     c
c      Department of Earth and Atmospheric Sciences                   c
c      Saint Louis University                                         c
c      221 North Grand Boulevard                                      c
c      St. Louis, Missouri 63103                                      c
c      U. S. A.                                                       c
c                                                                     c
c---------------------------------------------------------------------c
c-----
c       11 SEP 2000 add TP, TSV TSH to FILE96 format
c       14 JAN 2001 add A, C, F, L, N and density information
c       23 JUL 2009 - ensured consistent call to FFT routine with
c               complex arguments. We now call zfour
c-----
c-----
c       variables in common blocks
c
c       iftype  I*4 File type
c               1 - single trace
c               3 - three component
c               16 - Green's function
c               21 - Green's function
c
c       iobsyn  I*4 1 - observed
c               2 - synthetic
c
c       itmfrq  I*4 1 - time series
c               2 - Fourier spectra (not implemented)
c
c       iunit   I*4 1   - counts 
c               2   - cm
c               3   - cm/sec
c               4   - cm/sec/sec
c               5   - m
c               6   - m/sec
c               7   - m/sec/sec
c               8   - microns
c               9   - microns/sec
c               10  - microns/sec/sec
c               13  - 1/sec
c       junit   I*4
c               11  - Pa  (nt/m^2)
c               12  - MPa  (mega Pascal)
c
c       cfilt   C*80    comment on filtering operations
c
c       keyear  I*4 event year
c       kemon   I*4 event mon
c       keday   I*4 event day
c       kehour  I*4 event hour
c       kemin   I*4 event minute
c       esec    R*4 event second
c       evlat   R*4 event latitude
c       evlon   R*4 event longitude
c       evdep   R*4 event depth
c       
c       stname  C*8 station name
c               NOTE: AH(6), SEED(5), CSS(6), SAC(8)
c       
c       stlat   R*4 station latitude
c       stlon   R*4 station longitude
c       stelev  R*4 station elevation
c
c
c       distkm  R*4 epicentral distance in km
c       distdg  R*4 epicentral distance in degrees (111.195 km/degree)
c       evstaz  R*4 event -> station azimuth, degrees east of north
c       stevaz  R*4 station -> event azimuth, degrees east of north
c       
c       cpulse  C*80    pulse description
c       
c       ccomnt  C*80    comment
c
c       jsrc    I*4 Array of indices to indicate if traces are 
c                   present
c               iftype =  1  jsrc(1) indicates if trace is 
c                   present
c               iftype =  3  jsrc(i),i=1,5 indicates if trace 
c                   is present
c                   i = Z, 2=N, 3=E, 4=R, 5=T
c                   (we could carry all five, but the real 
c                   purpose  is to permit no more than 
c                   3 traces, but 2 may all that are 
c                   present in a real data set)
c
c
c-----

        common/ihdr96/ ftype, iftype, iobsyn, itmfrq, iunit, idecon,
     1      keyear, kemon, keday, kehour, kemin,
     2      ksyear, ksmon, ksday, kshour, ksmin,
     3      jsrc, junit
        integer*4 ftype, iftype, iobsyn, itmfrq, iunit, idecon
        integer*4 keyear, kemon, keday, kehour, kemin
        integer*4 ksyear, ksmon, ksday, kshour, ksmin
        integer*4 jsrc(21), junit

        common/rhdr96/ esec, distkm, distdg, evstaz, stevaz,
     1      evlat, evlon, evdep, stlat, stlon, stelev,
     2      TP, TSV, TSH, SA, SC, SF, SL, SN, SR
        real*4 esec, distkm, distdg, evstaz, stevaz
        real*4 evlat, evlon, evdep 
        real*4 stlat, stlon, stelev
        REAL*4 TP, TSV, TSH
        REAL*4 SA, SC, SF, SL, SN, SR

        common/chdr96/ stname, cfilt, cpulse, ccomnt
        character*8 stname*8, cfilt*80, cpulse*80, ccomnt*80    

        character stcomp*8
        real*4 cmpaz, cmpinc, cmpdt
        integer*4 npts
        real*4 ssec
c-----
c       command line arguments
c-----
        character*80 fnum, fden
        logical dowin
c-----
c       internal variables
c-----
        parameter (LER=6, LIN=5, LOT=6)
        parameter (NSAMP=16384)
        real x(NSAMP),y(NSAMP)
        complex zn(NSAMP), zd(NSAMP), zr(NSAMP)
        character svcomp*8
        real*8 sepoch, eepoch, o
c-----
c       call machine dependent initialization
c-----
        call mchdep()
c-----
c       parse command line arguments
c-----
        call gcmdln(fnum,fden,dowin)
        open(1,file=fnum,access='sequential',form='formatted',
     1      status='unknown')
        rewind 1
        open(2,file=fden,access='sequential',form='formatted',
     1      status='unknown')
        rewind 2
c-----
c       determine which Greens functions computed
c-----
        alpha = 0.01
c-----
c       process Green's Functions
c-----
  100   continue
            call rdhd96(1,nerr)
            if(nerr .lt. 0)go to 9999
            call rdhd96(2,nerr)
            if(nerr .lt. 0)go to 9999
                
c-----
c       identify the output stream
c-----
            iunit = 13
            call wrhd96(LOT,nerr)
            if(nerr .lt. 0)go to 9999
            do 200 j=1,21
                if(jsrc(j) .ne. 0)then
                    do 300 i=1,NSAMP
                        y(i)=0.0
  300               continue
                    call rdtr96(1,stcomp,cmpinc,cmpaz,
     1                  cmpdt, npts, ksyear, ksmon, 
     2                  ksday, kshour, ksmin, ssec, 
     3                  y,nerr,NSAMP)
                    if(nerr .lt. 0 )go to 9999
c-----
c                   set time of first sample
c-----
                    if(ksyear.eq.0)ksyear = 1970
                    if(keyear.eq.0)keyear = 1970
                    if(ksmon.eq.0)ksmon = 1
                    if(kemon.eq.0)kemon = 1
                    if(ksday.eq.0)ksday = 1
                    if(keday.eq.0)keday = 1
                    call htoe(ksyear,ksmon,ksday,kshour,
     1                  ksmin,ssec,sepoch)
                    call htoe(keyear,kemon,keday,kehour,
     1                  kemin,esec,eepoch)
                    o = sepoch - eepoch
                    ti = sngl(o)
                    dt = cmpdt
                    nn = npts
                    call npow2(nn)
                    afac = alpha/npts
                    fac = 1.0
                    dfac = exp(-afac)
        call domxmn(y,npts,depmax,depmin)
        write(0,*)'Numerator: ',depmin,depmax
                    do 1101 i=1,nn
                        if(i.le.npts)then
                        zn(i) = cmplx(fac*y(i),0.0)
                        fac = fac * dfac
                        else
                            zn(i) = cmplx(0.0,0.0)
                        endif
 1101               continue
                    call zfour(zn,nn,-1,dt,df)
c-----
c       work on denominator
c-----
                    do 301 i=1,NSAMP
                        y(i)=0.0
  301               continue
                    call rdtr96(2,stcomp,cmpinc,cmpaz,
     1                  cmpdt, npts, ksyear, ksmon, 
     2                  ksday, kshour, ksmin, ssec, 
     3                  y,nerr,NSAMP)
                    if(nerr .lt. 0 )go to 9999
                    afac = alpha/npts
                    fac = 1.0
                    dfac = exp(-afac)
        call domxmn(y,npts,depmax,depmin)
        write(0,*)'Denominator: ',depmin,depmax
                    do 1102 i=1,nn
                        if(i.le.npts)then
                        zd(i) = cmplx(fac*y(i),0.0)
                        fac = fac * dfac
                        else
                            zd(i) = cmplx(0.0,0.0)
                        endif
 1102               continue
                    call zfour(zd,nn,-1,dt,df)
c-----
c       do the decon
c-----
                    n21 = nn/2 + 1
                    zdmax = 0.0
                    do 1103 i=1,n21
                        if(cabs(zd(i)).gt.zdmax)then
                            zdmax = cabs(zd(i))
                        endif
 1103               continue
                    if(zdmax.eq.0.0)zdmax = 1.0
                    water = 0.001 * zdmax
        write(0,*)'zdmax,water:',zdmax,water
c-----
c       decon, but also shift output by about 1/16 number of points
c-----
                    tshift = nn*dt/16
        write(0,*)'nn,n21,tshift:',nn,n21,tshift
                    do 1104 i=1,n21
                        facd = cabs(zd(i))**2 + water**2
                        zr(i) = zn(i)*conjg(zd(i))/facd
c-----
c       time shift
c-----
                    omega=6.2831853*(i-1)*df
                    theta=omega*tshift
                    zr(i) = zr(i) * exp(-alpha*tshift)
     1                  *cmplx(cos(theta),-sin(theta))
c-----
c       cosine taper
c-----
                    if(dowin)then
                    tap = 0.5 * 
     1                  (1.0 + cos(3.1415927*(i-1)/nn))
                        zr(i) = zr(i) * tap
                    endif
                        if(i.gt.1)then
                            zr(nn+2-i)=conjg(zr(i))
                        endif
 1104               continue
                    zr(n21)=cmplx(real(zr(n21)),0.0)
                    call zfour(zr,nn,+1,dt,df)
                    afac = alpha/npts
                    fac = 1.0
                    dfac = exp(+afac)
                    do 1105 i=1,npts
                        y(i) = fac*real(zr(i))
                        fac = fac * dfac
 1105               continue
        call domxmn(y,npts,depmax,depmin)
        write(0,*)'Decon: ',depmin,depmax
                    call zfour(zr,nn,-1,dt,df)
c-----
c                               output the filtered trace
c-----
                                call wrtr96(LOT,stcomp,cmpinc,cmpaz,
     1                                  cmpdt, npts, ksyear, ksmon,
     2                                  ksday, kshour, ksmin, ssec,
     3                                  y,nerr,NSAMP)

                endif

  200       continue
        goto 100
 9999   continue
        end


        subroutine gcmdln(fnum,fden,dowin)
c-----
c       parse command line arguments
c
c       requires subroutine targ() and funtion mnmarg()
c
c-----
c-----
c       fnum    - Ch*   file name for numerator
c       fden    - Ch*   file name for denominator
c       dowin   - L .true. use a cosine taper
c-----
        
        character*20 name
        character fnum*(*), fden*(*)
        logical dowin
        integer mnmarg

        nmarg = mnmarg()
        dowin = .false.
        i = 0
        nfile = 0
   11   i = i + 1
        if(i.gt.nmarg)go to 13
            call mgtarg(i,name)
            if(name(1:3).eq.'-FN'.or.
     1          name(1:3).eq.'-fn')then
                i=i+1
                call mgtarg(i,fnum)
                nfile = nfile + 1
            else if(name(1:3).eq.'-FD'.or.
     1          name(1:3).eq.'-fd')then
                i=i+1
                call mgtarg(i,fden)
                nfile = nfile + 1
            else if(name(1:2).eq.'-W')then
                dowin = .true.
            else if(name(1:2).eq.'-?')then
                call usage()
            else if(name(1:2).eq.'-h')then
                call usage()
            endif
        go to 11
   13   continue
        if(nfile.ne.2)call usage()
        return
        end

        subroutine usage()
        integer*4 LER, LIN, LOT
        parameter(LER=6, LIN=5, LOT=6)
        write(LER,*)'USAGE: fdecon96  ',
     1      ' -FN file_num -FD file_den',
     2      ' [-W] [-?] [-h]'
        write(LER,*)
     1      ' -FN file_num (default none) numerator'
        write(LER,*)
     1      ' -FD file_den (default none) denominator'
        write(LER,*)
     1      ' -W  (default false) cosine taper in freq domain'
        write(LER,*)
     1      ' -?  This help message'
        write(LER,*)
     1      ' -h  This help message'
        stop
        end

        subroutine etoh(epoch,date,str,doy)
c-----
c       convert from epoch time to human time
c
c       epoch   - R*8 time in seconds relative to 0000 1 Jan 1970
c       date    - I*4 Julian date
c       str - C*  printable string
c       diy - I*4 dayin the year
c-----
        real*8 epoch
        integer*4 date
        character str*(*)
        integer*4 diy, doy, hour, minute, year, month
        integer*4 day
        real*4 second
        real*8 seclft

        str=' '
        doy = (epoch/86400.0d+00)
        seclft = dmod(epoch, 86400.0d+00)
        hour = 0
        minute = 0
        second = 0.00
c-----
c       compute hours minutes seconds
c-----
        if(seclft .ne. 0.00d+00)then
c-----
c                   before 1970 subtract and add a day
c-----
            if(seclft .lt. 0.0d+00)then
                doy = doy - 1
                seclft = seclft + 86400.00d+00
            endif
            hour = (seclft/3600.00d+00)
            seclft = dmod(seclft,3600.0d+00)
            minute = seclft/60.0d+00
            second = dmod(seclft,60.0d+00)
        endif

        if(doy .ge. 0)then
            year = 1970
 1000       continue
                diy =  leapdy(year)
                if(doy .lt. diy)go to 2000
                doy = doy - diy
                year = year + 1
            go to 1000
        else
            year = 1969
 1100       continue
                diy =  leapdy(year)
                doy = doy + diy
                if( doy .gt. 0 ) go to 2000
                year = year - 1
            go to 1100
        endif
 2000   continue
        doy = doy + 1
        date = year*1000 + doy
        call mnthdy(year,doy,month,day)
        write(str,110) year,month,day,hour,minute,second
  110   format(i4,i2,i2,i2,i2,f6.3)
c-----
c       guarantee that there are no blanks in the string str
c-----
        do 2100 i=1,17
            if(str(i:i).eq.' ')str(i:i)='0'
 2100   continue
        return
        end

        function leapdy(yr)
        integer*4 yr
        logical t1, t2, t3
        t1 = mod(yr,4).ne.0
        t2 = mod(yr,100).ne.0
        t3 = mod(yr,400).ne.0
        if( .not.t1 .and. t2)then
            isleap = 1
            leapdy = 366
        elseif( .not.t3)then
            isleap = 1
            leapdy = 366
        else
            isleap = 0
            leapdy = 365
        endif
        return
        end

        subroutine mnthdy(year,doy,month,day)
        integer*4 year, doy, month, day
        integer*4 i, dim, leap
        integer*4 dmnth(12)
        data dmnth/31,28,31,30,31,30,31,31,30,31,30,31/
        if(leapdy(year).eq.366)then
            leap = 1
        else
            leap = 0
        endif
        day = doy
        do 100 i=1,12
            month = i
            dim = dmnth(i)
            if(leap.eq.1 .and. i.eq.2)dim = dim + 1
            if(day .le.dim)goto 1000
            day = day - dim 
  100   continue
 1000   continue
        return
        end


        subroutine htoe(year,month,day,hour,minute,second,epoch)
c-----
c       convert calendar date to epoch time since January 1, 1970
c-----
c       year    - I*4   year
c       month   - I*4   month
c       day - I*4   day
c       hour    - I*4   hour
c       minute  - I*4   minute c    second  - I*4   second
c       second  - R*4   seconds
c       epoch   - R*8   time in seconds relative to 00:00 01 Jan 1970
c-----
        integer*4 year, month, day, hour, minute, date, diy
        real*4 second
        real*8 epoch, dtoepo
        integer*4 daymon(12)
        data daymon/0, 31, 59, 90, 120, 151, 181, 212, 243, 273,
     1      304, 334/
        diy = daymon(month) + day
        if(leapdy(year).eq.366 .and. month.gt.2)diy=diy+1
        date = 1000*year + diy
c       write(6,*)'date=',date
        epoch = dtoepo(date) + hour * 3600.0d+00 + 
     1      minute * 60.0d+00 +dble(second)
        return
        end

c-----
c       convert julian date to epoch time
c-----
        function dtoepo(date)
        real*8 dtoepo
        integer*4 date, diy, cnt, days

        cnt = date / 1000
        days = 0
        if (cnt .gt. 1970)then
 1000       continue
            cnt = cnt -1
            if(cnt.lt.1970)go to 2000
                days = days + leapdy(cnt)
            go to 1000
        else if (cnt .lt. 1970)then
 1100       continue
            if(cnt.ge.1970)goto 2000
                days = days - leapdy(cnt)
                cnt = cnt + 1
            go to 1100
        endif
 2000   continue
        diy = (date -1) / 1000
        diy = (date -1 ) -  1000*diy
c       write(6,*)'days=',days,' diy=',diy
        dtoepo = (days + diy) * 86400.0d+00
        return
        end


        subroutine zfour(zarr,nn,isign,dt,df) 
c-----
c     THE input is a complex array
c     which has numbers stored in memory as
c     R1, I1, R2, I2, ..., Rnn, Inn
c     where nn must be a power of 2 R and I are the real and imaginary
c     parts of the complex number
c
c     For isign -1 this is a complex time series
c     For isign +1 this is a complex frequency series with
c        index 1 (in fortran corresponding to f=0
c              2                              f=df
c            nn/2 + 1                         f = 1/2dt = Nyquist
c            nn - 1                           f = -2df
c             nn                              f = -df

c-----
c     the cooley-tookey fast fourier transform in usasi basic fortran
c     transform(j) = sum(datc(i)*w**((i-1)(j-1)), where i and j run
c     from 1 to nn and w = exp(isign*2*pi*sqrt(-1)/nn).  datc is a one-
c     dimensional complex array (i.e., the real and imaginary parts of
c     datc are located immediately adjacent in storage, such as fortran
c     places them) whose length nn is a power of two.  isign
c     is +1 or -1, giving the sign of the transform.  transform values
c     are returned in array datc, replacing the input datc.  the time is
c     proportional to n*log2(n), rather than the usual n**2
c     rms resolution error being bounded by 6*sqrt(i)*log2(nn)*2**(-b),
c     b is the number of bits in the floating point fraction.
c
c     the program computes df from dt, dt from df and checks to see
c     if they are consistent. In addition, the transforms are multiplied
c     by dt or df to make the results dimensionally correct
c
c     This is a slightly modified version of the original Brenner routine
c     The changes were to add physical dimensions to the transform
c     and to make it all complex
c-----
        complex zarr(*) 
        integer nn, isign
        real dt, df

        complex ztemp
c-----
c       ensure that the dt and df are defined and
c       consistent
c-----
        if(dt.eq.0.0) dt = 1./(nn*df) 
        if(df.eq.0.0) df = 1./(nn*dt) 
        if(dt.ne.(nn*df)) df = 1./(nn*dt) 
c-----
c       now begin the transform
c-----
        jj = 1
        do 5 ii=1,nn 
        if(ii .lt. jj) then
              ztemp = zarr(jj)
              zarr(jj) = zarr(ii)
              zarr(ii) = ztemp
        endif
        mm = nn/2
    3   continue
        if(jj.le.mm) then
            go to 55
        else 
              jj = jj - mm
              mm = mm /2
              if(mm.lt.1)then
                  go to 55
              else
                  go to 3
              endif
        endif
   55   continue
        jj = jj + mm
    5   continue
        mmax = 1 
c-----
    6 continue
        if(mmax .lt. nn)then
            go to 7
        else if(mmax .ge. nn)then
            go to 10
        endif
    7   continue
        istep= 2*mmax 
        theta = 6.283185307/(isign*2.0*mmax) 
        sinth=sin(theta/2.) 
        wstpr=-2.*sinth*sinth 
        wstpi=sin(theta) 
        wr=1.0 
        wi=0.0 
        do 9 mm=1,mmax
              do 8 ii=mm,nn,istep
                    jj=ii+mmax
                    ztemp=cmplx(wr,wi)*zarr(jj)
                    zarr(jj) = zarr(ii) - ztemp
                    zarr(ii) = zarr(ii) + ztemp
    8         continue
c-----
c       use trig relations to compute the next sin/cos
c       without actually calling sin() or cos()
c-----
              tempr = wr 
              wr = wr*wstpr-wi*wstpi + wr 
              wi = wi*wstpr+tempr*wstpi + wi 
    9   continue
        mmax = istep 
        go to 6 
c-----
c       transform is done
c-----
   10   continue 
c-----
c     give the arrays the proper physical dimensions
c-----
        if(isign.lt.0)then
c-----
c             time to frequency domain
c-----
              do  ii = 1,nn
                    zarr(ii) = zarr(ii) * dt
              enddo
        else
c-----
c             frequency to time domain
c-----
              do ii = 1,nn
                    zarr(ii) = zarr(ii) * df
              enddo
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
        npts = 2*npts
        if(npts.lt.nsamp)goto 1000
        return
        end

        subroutine chtofp(str,fout)
c------
c       routine to convert string to floating point
c       The E format is accepted as well as free form
c       input
c
c       This assumes that the string is 20 characters long
c-----
        character*(*) str
        real*4 fout
        integer*4 lgstr
        logical hase
        l = lgstr(str)
c------
c       If the string str contains an E or e, then
c       we parse the result using an E format.
c------
        hase = .false.
        do 100 i=1,l
            if(str(i:i) .eq.'e' .or. str(i:i).eq.'E')then
                hase = .true.
            endif
  100   continue
c-----
c       read floating point number
c-----
        if(hase)then
            read(str,'(bn,e20.13)')fout
        else
            read(str,'(bn,f20.0)') fout
        endif
        return
        end

        subroutine domxmn(x,npts,depmax,depmin)
c-----
c       get extremal values of the time series
c-----
        real*4 x(*)
        real*4 depmax,depmin,depmen
        integer*4 npts
        depmax = -1.0e+38
        depmin =  1.0e+38
        sum = 0.0
        do 1000 i=1, npts
            if( x(i) .gt. depmax) depmax = x(i)
            if( x(i) .lt. depmin) depmin = x(i)
            sum = sum + x(i)
 1000   continue
        return
        end
