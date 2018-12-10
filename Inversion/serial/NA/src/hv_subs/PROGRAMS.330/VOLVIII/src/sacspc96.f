        program sacspc96
c----------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME II                                                      c
c                                                                     c
c      PROGRAM: SACSPC96                                              c
c                                                                     c
c      COPYRIGHT 2000                                                 c
c      R. B. Herrmann                                                 c
c      Department of Earth and Atmospheric Sciences                   c
c      Saint Louis University                                         c
c      221 North Grand Boulevard                                      c
c      St. Louis, Missouri 63103                                      c
c      U. S. A.                                                       c
c                                                                     c
c----------------------------------------------------------------------c
c
c      Plot the Fourier Amplitude spectra of a SAC time series
c
c      Real purpose is to overlay with spectral amplitude plots
c      from other programs in this package
c-----
c       CHANGES
c         01 NOV 2007 - changed the color flag from -P kolor to -K kolor for consistency
c            with other programs
c       23 JUL 2009 - ensured consistent call to FFT routine with
c               complex arguments. We now call zfour
c-----
        integer LIN, LOT, LER
        parameter (LIN=5, LER=0, LOT=6)
        integer NPTSER, NPTFRQ
        parameter (NPTSER=8200, NPTFRQ=5000)
        complex z(NPTSER)
        real  tarr(NPTSER)

c-----
c      command line arguments
c-----
        logical dofreq, dobox, doxlog, doylog
        real*4 xmin, xmax, ymin, ymax, x0, y0, xlen, ylen
        integer*4 kolor

        integer*4 n , n21 
        real*4  deg , baz , t0
        character sta*8, comp*8, cdate*12
        integer ierr

        character sacfile*100
c-----
c      machine dependent
c-----
        call mchdep()
c-----
c      get command line arguments
c-----
        sacfile = ' '
        call gcmdln(dofreq, dobox, xmin, xmax, ymin, ymax, 
     1      x0, y0, xlen, ylen, 
     2      kolor, doxlog,doylog,sacfile)
        if(sacfile .eq. ' ')call usage('No SAC file specified')
c-----
C       write(LER,*)'dofreq   :', dofreq
C       write(LER,*)'dobox    :', dobox 
C       write(LER,*)'xmin     :', xmin  
C       write(LER,*)'xmax     :', xmax  
C       write(LER,*)'ymin     :', ymin  
C       write(LER,*)'ymax     :', ymax  
C       write(LER,*)'x0       :', x0    
C       write(LER,*)'y0       :', y0    
C       write(LER,*)'xlen     :', xlen    
C       write(LER,*)'ylen     :', ylen    
C       write(LER,*)'kolor    :', kolor    
C       write(LER,*)'doxlog   :', doxlog    
C       write(LER,*)'doylog   :', doylog    
c-----
c      initialize CALPLOT graphics
c-----
        call pinitf('SACSPC96.PLT')
c-----
c      open the SAC file, here assume that it is binary
c-----
        call gettrc(sacfile,n,n21,npts,dist,deg,az,baz,t0,dt,sta,
     1      comp,cdate,z,tarr,NPTSER,ierr)
c-----
c      make a spectra
c-----
        tmin =  1.0e+38
        tmax = -1.0e+38
        do 1000 i=1,n21
            tarr(i) = cabs(z(i))
            if(tarr(i).gt.tmax)tmax = tarr(i)
            if(tarr(i).lt.tmin)tmin = tarr(i)
 1000   continue
c-----
c      verify the plot limits
c-----
        call pltlim(n21,n,dt,xmin,xmax,ymin,ymax,doxlog,doylog,
     1      dofreq,tmax,tmin)
C       write(LER,*)'tmin     :', tmin  
C       write(LER,*)'tmax     :', tmax  
C       write(LER,*)'xmin     :', xmin  
C       write(LER,*)'xmax     :', xmax  
C       write(LER,*)'ymin     :', ymin  
C       write(LER,*)'ymax     :', ymax  
c-----
c      plot it
c-----
        call pltspc(dofreq, dobox, xmin, xmax, ymin, ymax, 
     1      x0, y0, xlen, ylen, kolor, doxlog,doylog,
     2      n,n21,npts,tarr,dt)
        call pend()
        end

        subroutine gcmdln(dofreq, dobox, xmin, xmax, 
     1      ymin, ymax, 
     1      x0, y0, xlen, ylen, 
     2      kolor, doxlog,doylog,sacfile)
        logical dofreq, dobox, doxlog, doylog
        real*4 xmin, xmax, ymin, ymax
        real*4 x0, y0, xlen, ylen
        integer*4 kolor

        character names*80
        character sacfile*(*)
c-----
c      defaults
c-----
        dofreq = .true.
        dobox = .true.
        xmin = -1.0
        xmax = -1.0
        ymin = -1.0
        ymax = -1.0
        x0 = 2.0
        y0 = 1.0
        xlen = 6.0
        ylen = 6.0
        kolor = 1
        doxlog = .false.
        doylog = .false.

        nmarg = mnmarg()
        i = 0
 1000   continue
            i = i + 1
            if(i.gt.nmarg)go to 2000
            call mgtarg(i,names)
            if(names(1:5).eq.'-FREQ')then
                dofreq = .true. 
            else if(names(1:4).eq.'-PER')then
                dofreq = .false.    
            else if(names(1:5).eq.'-XMIN')then
                i = i + 1
                call mgtarg(i,names)
                read(names,'(bn,g20.7)')xmin
            else if(names(1:5).eq.'-XMAX')then
                i = i + 1
                call mgtarg(i,names)
                read(names,'(bn,g20.7)')xmax
            else if(names(1:5).eq.'-YMIN')then
                i = i + 1
                call mgtarg(i,names)
                read(names,'(bn,g20.7)')ymin
            else if(names(1:5).eq.'-YMAX')then
                i = i + 1
                call mgtarg(i,names)
                read(names,'(bn,g20.7)')ymax
            else if(names(1:3).eq.'-X0')then
                i = i + 1
                call mgtarg(i,names)
                read(names,'(bn,g20.7)')x0
            else if(names(1:3).eq.'-Y0')then
                i = i + 1
                call mgtarg(i,names)
                read(names,'(bn,g20.7)')y0
            else if(names(1:5).eq.'-XLEN')then
                i = i + 1
                call mgtarg(i,names)
                read(names,'(bn,g20.7)')xlen
            else if(names(1:5).eq.'-YLEN')then
                i = i + 1
                call mgtarg(i,names)
                read(names,'(bn,g20.7)')ylen
            else if(names(1:2).eq.'-K' .and.
     1              names(1:4).ne.'-PER')then
                i = i + 1
                call mgtarg(i,names)
                read(names,'(bn,i20)')kolor
            else if(names(1:6).eq.'-NOBOX')then
                dobox = .false.
            else if(names(1:5).eq.'-XLOG')then
                doxlog = .true.
            else if(names(1:5).eq.'-YLOG')then
                doylog = .true.
            else if(names(1:2).eq.'-f')then
                i = i + 1
                call mgtarg(i,sacfile)
            else if(names(1:2).eq.'-?')then
                call usage(' ')
            else if(names(1:2).eq.'-h')then
                call usage(' ')
            endif
        go to 1000
 2000   continue
        return
        end

        subroutine usage(str)
        character str*(*)
        parameter (LER=0,LIN=5,LOT=6)
        write(LER,*)'sacspc96: ',str
        write(LER,*)'USAGE:',
     1  'sacspc96 ',
     1  '[-FREQ -PER]',
     1  '-XMIN xmin -XMAX xmax -YMIN ymin -YMAX ymax ',
     1  '-X0 x0 -Y0 y0 -K kolor -NOBOX -XLOG -YLOG -? -h',
     1  '-f sacfile_name'
        write(LER,*)
     1  '-FREQ      (default true) X-Axis is frequency'
        write(LER,*)
     1  '-PER       (default false)X-Axis is period'
        write(LER,*)
     1  '-XMIN xmin (default 0.0)  minimum value of X-Axis'
        write(LER,*)
     1  '-XMAX xmax (default    )  maximum value of X-Axis'
        write(LER,*)
     1  '-YMIN ymin (default 0.0)  minimum value of Y-Axis'
        write(LER,*)
     1  '-YMAX ymax (default 0.0)  maximum value of Y-Axis'
        write(LER,*)
     1  '-X0 x0     (default 2.0)  lower left corner of plot'
        write(LER,*)
     1  '-Y0 y0     (default 1.0)  bottom left corner of plot'
        write(LER,*)
     1  '-XLEN xlen (default 6.0)  length of X-Axis'
        write(LER,*)
     1  '-YLEN ylen (default 6.0)  length of Y-Axis'
        write(LER,*)
     1  '-P kolor   (default 1  )  color for curves'
        write(LER,*)
     1  '-NOBOX     (default false) do not plot axes'
        write(LER,*)
     1  '-XLOG      (default linear) X axis is logarithmic'
        write(LER,*)
     1  '-YLOG      (default linear) Y axis is logarithmic'
        write(LER,*)
     1  '-f sacfilename  (default none) Binary SAC file name'
        write(LER,*)
     1  '-h         (default false) online help'
        write(LER,*)
     1  '-?         (default false) online help'
        stop 
        end

        subroutine gettrc(xfile,n,n21,npts,dist,deg,az,baz,t0,dt,sta,
     1      comp,cdate,z,tarr,nptser,ierr)
c-----
c      get parameters required for processing
c-----
c      xfile   Ch* - name of data file
c      n   I*4 - array size, power of 2
c      n21 I*4 - n/2 + 1
c      npts    I*4 - original number of points (n >= npts)
c      dist    R*4 - epicentral distance (km)
c      deg R*4 - epicentral distance (deg)
c      az  R*4 - src -> rec azimuth
c      baz R*4 - rec -> src azimuth
c      t0  R*4 - time of first sample after origin time
c      dt  R*4 - sample interval
c      sta Ch*8    - station name
c      comp    Ch*8    - station component
c      cdate   Ch*12   - date string
c      z   C*4 - complex Fourier Transform array of samples
c      tarr    R*4 - temporary array  for trace
c      nptser  I*4 - length of z() and tarr() arrays
c      ierr    I*4 - error code
c-----
        character xfile*(*)
        integer n, n21, npts, nptser, ierr
        real dist, deg, az, baz, t0, dt
        character sta*8, comp*8, cdate*12
        complex z(nptser) 
        real tarr(nptser)
c-----
c      iunit   I*4 Physical units of amplitude
c              -1 counts (default)
c              0 meters
c              1 centimeters
c              2 nanometers
c      idva    I*4 input time series (note output spectra will be 
c              spectra is displacement in cm, e.g., cm-sec of count-sec
c              0   displacment
c              1   velocity
c              2   acceleration
c      unitcor R*4 unit correction factor to get to cm
c-----
        common/units/idva,iunit,unitcor
        integer*4 idva, iunit
        real*4 unitcor

        call getsac(xfile,npts,dist,deg,az,baz,t0,dt,sta,
     1                 comp,cdate,tarr)
        ls = lgstr(xfile)
        if(npts.gt.nptser)npts = nptser
        call npow2(npts,n,n21)
c-----
c      fill up the waveform array, use the opportunity to set
c      the units correctly
c-----
        do 100 i=1,n
            if(i.le.npts)then
                z(i) = cmplx(tarr(i),0.0)
            else
                z(i) = cmplx(0.0,0.0)
            endif
  100   continue
        call zfour(z,n,-1,dt,df)
        return
        end
c
        subroutine getsac(name,npts,dist,deg,az,baz,t0,dt,sta,
     1                 comp,cdate,seis)
c-----
c
c      name    - file name to write
c      n   - number of points in FFT must be power of 2
c      n21 - number of frequencies = n/2 + 1
c      npts    - number of points in original time series
c          - which may have been zero filled to make power of 2
c      dist    - epicentral distance in km
c      deg - epicentral distance in degrees
c      az  - source - receiver azimuth in degrees
c      baz - receiver-source back azimuth
c      t0  - time of first sample after origin
c      dt  - sampling interval
c      sta - C*4 station name string
c      comp    - C*4 component name string
c      cdate   - C*12 date string
c      z   - COMPLEX array of spectra
c-----
        integer MXPTS
        parameter (MXPTS = 64000)
        integer  npts
        real dist, deg, az, baz, t0, dt
        character sta*8, comp*8, cdate*12
        character kstnm*8, kcmpnm*8
        character name*(*)
        real seis(MXPTS)
C
        real evla, evlo, evdp, stla, stlo, stel, beg, origtime
        integer ntimes(6)
        common /shdr/ evla,evlo,evdp,stla,stlo,stel,
     &               beg,ntimes,origtime

        common/sachdr/rhdr,ihdr,chdr
        real*4 rhdr(70)
        integer*4 ihdr(40)
        character*4 chdr(48)
C
        call brsac (1,MXPTS,name,seis,ierr)
C
        call getfhv('AZ      ',az,nerr)
        call getfhv('BAZ     ',baz,nerr)
        call getfhv('DIST    ',dist,nerr)
        call getfhv('GCARC   ',deg,nerr)
        call getfhv('DELTA   ', dt, nerr)
        call getnhv('NPTS    ',npts,nerr)
        call getfhv('B       ', beg, nerr)
        call getfhv('O       ',origtime,nerr)
        call getfhv('EVLA    ',evla,nerr)
        call getfhv('EVLO    ',evlo,nerr)
        call getfhv('STLA    ',stla,nerr)
        call getfhv('STLO    ',stlo,nerr)
        call getkhv('KSTNM   ',kstnm,nerr)
        call getkhv('KCMPNM  ',kcmpnm,nerr)
C
        if(nerr .eq. 0 .and. origtime .ne. -12345)then
            t0 = beg - origtime
        else
            t0 = beg
        end if
C
        tp = tp - origtime
        ts = ts - origtime
C       write(6,*)'t0,tp,ts:',t0,tp,ts
        sta = kstnm(1:8)
        comp = kcmpnm(1:8)
C       write(6,*)'name,npts,dist,deg,az,baz,t0,dt,sta',
C     1     name,npts,dist,deg,az,baz,t0,dt,sta
        cdate = ' '
C
C
        return
        end

        subroutine npow2(nsamp,npts,npts21)
c-----
c      Given nsamp, find npts >= nsamp such that npts is a power of 2
c-----  
        integer*4 nsamp, npts, npts21
        npts = 1
 1000   continue
        npts = 2*npts
        if(npts.lt.nsamp)go to 1000
        npts21 = npts/2 + 1
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

        subroutine pltspc(dofreq, dobox, xmin, xmax, ymin, ymax, 
     1      x0, y0, xlen, ylen, kolor, doxlog,doylog,
     2      n,n21,npts,tarr,dt)
        logical dofreq, dobox, doxlog, doylog
        real xmin, xmax, ymin, ymax, x0, y0, xlen, ylen, dt
        real tarr(n)
        integer n, n21, npts

C       write(6,*)dofreq,dobox,xmin,xmax,ymin,ymax,x0,y0
C       write(6,*)xlen,ylen,kolor,doxlog,doylog
C       write(6,*)n,n21,npts,dt
c-----
        df = 1.0/(n*dt)
c-----
c      put up the axes
c-----
        call newpen(1)
        xlow = x0 
        ylow = y0
        xhgh = x0 + xlen
        yhgh = y0 + ylen
        call gbox(xlow,ylow,xhgh,yhgh)
        if(dobox)then
            if(doxlog)then
            if(dofreq)then
                call dologx(x0,y0,xlen,xmax,xmin,0.10,
     1          .false.,.false.,.true.,14,'Frequency (Hz)')
            else
                call dologx(x0,y0,xlen,xmax,xmin,0.10,
     1          .false.,.false.,.true.,10,'Period (s)')
            endif
            else
            if(dofreq)then
                call dolinx(x0,y0,xlen,xmax,xmin,0.10,
     1          .false.,.false.,.true.,14,'Frequency (Hz)')
            else
                call dolinx(x0,y0,xlen,xmax,xmin,0.10,
     1          .false.,.false.,.true.,10,'Period (s)')
            endif
            endif

            if(doylog)then
                call dology(x0,y0,ylen,ymax,ymin,0.07,
     1          .true.,.true.,.true., 0,' ')
            else
                call doliny(x0,y0,ylen,ymax,ymin,0.07,
     1          .true.,.true.,.true., 0,' ')
            endif
        else
c-----
c      put in corner markers
c-----
            call plot(x0+0.14,y0+0.00,3)
            call plot(x0+0.00,y0+0.00,2)
            call plot(x0+0.00,y0+0.14,2)

            call plot(x0+xlen-0.14,y0+0.00,3)
            call plot(x0+xlen-0.00,y0+0.00,2)
            call plot(x0+xlen-0.00,y0+0.14,2)

            call plot(x0+xlen-0.14,y0+ylen-0.00,3)
            call plot(x0+xlen-0.00,y0+ylen-0.00,2)
            call plot(x0+xlen-0.00,y0+ylen-0.14,2)

            call plot(x0+0.14,y0+ylen-0.00,3)
            call plot(x0+0.00,y0+ylen-0.00,2)
            call plot(x0+0.00,y0+ylen-0.14,2)

        endif
c-----
c      for plotting we invoke clipping
c-----
        call gclip('on' , xlow,ylow,xhgh,yhgh)
c-----
c      plot with the correct color
c-----
        call newpen(kolor)
c-----
c      do the plotting
c-----
        ipen = 3
        do 1000 i=1,n21
            freq = (i-1)*df
c-----
c          UGH!! GO TO's for SAFETY
c-----
            if(i.eq.1 .and. doxlog .or. i.eq.1 
     1          .and. .not. dofreq)go to 1000
            if(dofreq)then
                xval = freq
            else
                xval = 1.0/freq
            endif
            yval = tarr(i)
c-----
c      now map into the plot space
c-----
            if(doxlog)then
                xx = x0 + xlen*alog10(xval/xmin)
     1                                  /alog10(xmax/xmin)  
            else
                xx = x0 + xlen*(xval - xmin)/(xmax - xmin)
            endif
            if(doylog)then
c-----
c              beware of log of 0.0
c-----
                if(yval.lt. 0.001*ymin)then
                    yval = 0.001*ymin
                endif
                yy = y0 + ylen*alog10(yval/ymin)
     1                                  /alog10(ymax/ymin)  
            else
                yy = y0 + ylen*(yval - ymin)/(ymax - ymin)
            endif
            call plot(xx,yy,ipen)
            ipen = 2
 1000   continue
        call plot(xx,yy,3)
c-----
c      reset the plot state
c-----
        call gclip('off', xlow,ylow,xhgh,yhgh)
        call newpen(1)
        return
        end

        subroutine pltlim(n21,n,dt,xmin,xmax,ymin,ymax,doxlog,
     1      doylog,dofreq,tmax,tmin)
c-----
c      if xmax, xmin, ymax, ymin not given from the command line, 
c      generate here - also ensure that xmax > xmin, ymax > ymin
c
c      plot limit is easy - since the plot space is x>-0, y>=0
c-----
        real tmin, tmax
        real dt, xmin, xmax, ymin, ymax
        integer n, n21
        logical doxlog, doylog, dofreq

        df = 1.0/(n*dt)
c-----
c      check the X-axis 
c-----
        if(xmax.lt.0.0)then
            if(dofreq)then
                xmax = n21*df
            else
                xmax = 1.0/df
            endif
        endif
        if(xmin.lt.0.0)then
            if(doxlog)then
c-----
c              force three cycles for log axis
c-----
                xmin = xmax/1000.0
            else
                xmin = 0.0
            endif
        endif
c-----
c      check the Y-axis 
c-----
        if(ymax.lt.0.0)then
            ymax = tmax
        endif
        if(ymin.lt.0.0)then
            if(doylog)then
c-----
c              force three cycles for log axis
c-----
                ymin = ymax/1000.0
            else
                ymin = tmin
            endif
        endif

c-----
c      check order
c-----
        if(xmin.gt.xmax)then
            tmp = xmin
            xmin = xmax
            xmax = tmp
        endif
        if(ymin.gt.ymax)then
            tmp = ymin
            ymin = ymax
            ymax = tmp
        endif
        return
        end
