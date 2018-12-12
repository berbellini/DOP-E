        program pltsac
c---------------------------------------------------------------------c 
c                                                                    c 
c     COMPUTER PROGRAMS IN SEISMOLOGY                                c 
c     VOLUME V                                                       c 
c                                                                    c
c     PROGRAM: PLTSAC                                                c
c                                                                    c
c     COPYRIGHT 2002                                                 c
c     R. B. Herrmann                                                 c
c     Department of Earth and Atmospheric Sciences                   c
c     Saint Louis University                                         c
c     221 North Grand Boulevard                                      c
c     St. Louis, Missouri 63103                                      c
c     U. S. A.                                                       c
c                                                                    c
c---------------------------------------------------------------------c
c      CHANGES
c      15 SEP 2008 - increased the size of the amplitude and time shift
c              strings
c      07 MAY 2010 - added -TMIN tmin and -TMAX tmax flags to make
c              time scale prettier
c      05 SEP 2010 - plot also the percentage of fit given in USER5
c-----
        integer MXPTS
        parameter (MXPTS = 64000)
        real seis(MXPTS)
        character sta*8, comp*8, cdate*12
        integer iplmn, kolor
        integer MXFILE
        parameter (MXFILE=100)
        character*80 fname(MXFILE)
        real tmin, tmax
        logical dotmnmx
        integer nfile
        logical doabs
        logical doover
        integer ipen
        logical doamp
        logical dotmtp, dotmbt
        logical dosctp, doscbt
        logical doovershd
        logical douser9

c-----
c      get plot controls
c-----
        call gcmdln(x0,y0,xlen,ylen,nmin,nmax,MXPTS,ybot,
     1      iplmn,kolor,fname,nfile,MXFILE,doabs,pcy,doover,doamp,
     1      dotmtp, dotmbt,dosctp, doscbt, doovershd,douser9,
     2      tmin,tmax,dotmnmx)
        y0keep = y0
        
        
        call pinitf('PLTSAC.PLT')
        call newpen(1)
c-----
c      if absolute scaling is desired, get individual trace extrema
c-----
        ampmax = 0.0
        nfiletoplt = 0
        if(doabs)then
        if(nfile.gt.0)then
            do 990 i=1,nfile
            call getsac(fname(i),npts,dist,deg,az,baz,t0,dt,sta,
     1                      comp,cdate,seis,tp,ts,user9,user5,
     2              nzyear,nzjday,nzhour,nzmin,nzsec,beg)
            if(npts.gt.nmax)npts=nmax
            if(npts.lt.nmin)go to 990
            call dmxmn(seis,nmin,nmax,ampmx,ampmn)
C       WRITE(0,*)ampmax,ampmx,ampmn
            if(ampmx.gt.ampmax)ampmax = ampmx
            if(abs(ampmn).gt.ampmax)ampmax = abs(ampmn)
            nfiletoplt = nfiletoplt + 1
  990       continue
        endif
        if(ampmax.eq.0.0)ampmax = 1.0e-3
        endif

c-----
c      now plot the traces
c-----
        dy = ylen/nfile
        do 1000 i=1,nfile
            
        call getsac(fname(i),npts,dist,deg,az,baz,t0,dt,sta,
     1                 comp,cdate,seis,tp,ts,user9,user5,
     2          nzyear,nzjday,nzhour,nzmin,nzsec,beg)
        nh = npts
        nl = 1
        if(nh.gt.nmax)nh=nmax
        if(nmin.lt.1)then
            nl = 1
        else
            nl = nmin
        endif
        if(nl.gt.nh)then
            kk = nh
            nh = nl
            nl = kk
        endif
        if(kolor.ge.0)then
                ipen = kolor
        else
            if(nfile.eq.1)then
                ipen = 2
            else
                ipen = 1000 + real(100 * 
     1              real(i-1)/real(nfile-1)) 
                if(ipen.gt.1100)ipen = 1100
                if(ipen.lt.1000)ipen = 1000
            endif
        endif
        call plotit(seis,nl,nh,npts,x0,y0,xlen,ylen,deg,user9,user5,
     2          nzyear,nzjday,nzhour,nzmin,nzsec,fname,
     3          iplmn,ipen,doabs,ampmax,pcy,dy,ybot,doamp,
     4          douser9)
        if(.not.doover)then
            y0 = y0 - dy
        endif
        if(doovershd)then
            kolor = 1
            y0 = y0 - 0.25 *dy
        endif
 1000   continue
        call newpen(1)
        twnmin = beg + (nl -1)*dt
        twnmax = beg + (nh -1)*dt
        if(dotmnmx)then
           twnmin = tmin
           twnmax = tmax
        endif
            ypos = y0
        if(dotmbt)then
            if(doscbt)then
        call dolinx(x0,ypos,xlen,twnmax,twnmin,0.07,
     1      .true.,.false.,.true.,10,'Time (sec)')
            else
        call dolinx(x0,ypos,xlen,twnmax,twnmin,0.07,
     1      .true.,.false.,.false.,1 ,' ')
            endif
C        write(6,*)twnmax,twnmin,y0,ypos
        endif
        if(dotmtp)then
            if(dosctp)then
        call dolinx(x0,y0keep+0.5,xlen,twnmax,twnmin,0.07,
     1      .false.,.true.,.true.,10,'Time (sec)')
            else
        call dolinx(x0,y0keep+0.5,xlen,twnmax,twnmin,0.07,
     1      .false.,.true.,.false.,1 ,' ')
            endif
        endif
        call pend()
        end

        subroutine plotit(x,nl,nh,npts,x0,y0,xlen,ylen,deg,user9,user5,
     2          nzyear,nzjday,nzhour,nzmin,nzsec,fname,
     3          iplmn,kolor,doabs,ampmax,pcy,dy,ybot,doamp,
     4          douser9)
        real*4 x(npts)
        character fname*(*)
        character outstr*14
        logical doabs
        logical doamp
        logical douser9
        real ymxval
        ymax = abs(x(nl))
        do 1001 i=nl,nh
            if(abs(x(i)) .gt.ymax)ymax = abs(x(i))
 1001   continue
        ymxval = ymax
        if(ymax.eq.0.0)ymax = 1.0
        if(doabs)then
                ymax = ampmax
        endif
        if(ymax.lt.ybot)ymax = 1.0e+30
        if(npts.le.0)return
        dx = xlen/(nh-nl+1)
C       write(6,*)ymax,x0,y0,xlen,ylen,pcy,dy,iplmn,nl,nh,npts,dx
        ipen = 3
        xx = x0
        if(iplmn.ne.0)then
            if(iplmn.lt.0) call shdsei(x0,y0,0,1,1)
            if(iplmn.gt.0) call shdsei(x0,y0,0,1,0)
        endif
        call newpen(kolor)
        do 1000 i=nl,nh
            xx = xx + dx
            yy = y0 + pcy*dy*(x(i) )/(ymax )
            call plot(xx,yy,ipen)
            ipen= 2
 1000   continue
        if(iplmn.ne.0)then
            call shdsei(x0,y0,0,0,0)
        endif
            ht = 0.1*dy
            HT = 0.15*DY
        if(doamp)then
            if(ht.lt.0.07)ht = 0.07
            if(doabs)then
            call  number(x0,y0+0.4*dy,ht,ampmax,0.0,2002)
            else
            call  number(x0,y0+0.4*dy,ht,ymxval,0.0,2002)
            endif
            call plot(x0,y0,3)
        endif
        call newpen(1)
C        write(0,*)douser9,user9,fname
        if(douser9 .and. user9 .ne. -12345.)then
            call number(x0+xlen-5.*ht,y0+0.4*dy,ht,user9,0.0,2)
        endif
        if(douser9 .and. user5 .ne. -12345.)then
            call number(x0+xlen-5.*ht,y0+0.2*dy,ht,user5,0.0,-1)
            call symbol(999.,y0+0.2*dy,ht,'%',0.0,1)
        endif
c-----
c      after shading, plot trace in black as outline
c-----
        if(iplmn.ne.0)then
        xx = x0
        call plot(x0,y0,3)
        do 2000 i=nl,nh
            xx = xx + dx
            yy = y0 + pcy*dy*(x(i) )/(ymax )
            call plot(xx,yy,ipen)
            ipen= 2
 2000   continue
        endif
        return
        end


        subroutine getsac(name,npts,dist,deg,az,baz,t0,dt,sta,
     1                 comp,cdate,seis,tp,ts,user9,user5,
     2          nzyear,nzjday,nzhour,nzmin,nzsec,beg)
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
c          nzyear,nzjday,nzhour,nzmin,nzsec)
c-----
        integer MXPTS
        parameter (MXPTS = 64000)
        integer  npts
        real dist, deg, az, baz, t0, dt
        character sta*8, comp*8, cdate*12
        character kstnm*8, kcmpnm*8
        character name*(*)
        real seis(MXPTS)
        integer*4  nzyear,nzjday,nzhour,nzmin,nzsec
*
        real evla, evlo, evdp, stla, stlo, stel, beg, origtime

        common/sachdr/rhdr,ihdr,chdr
        real*4 rhdr(70)
        integer*4 ihdr(40)
        character*4 chdr(48)
*
        
        call brsac (3,MXPTS,name,seis,ierr)
*
        call getfhv('AZ      ',az,nerr)
        call getfhv('BAZ     ',baz,nerr)
        call getfhv('DIST    ',dist,nerr)
        call getfhv('GCARC   ',deg,nerr)
        call getfhv('DELTA   ', dt, nerr)
        call getnhv('NPTS    ',npts,nerr)
        call getfhv('B       ', beg, nerr)
        call getfhv('O       ',origtime,nerr)
        call getfhv('A       ',tp,nerr)
        call getfhv('T0      ',ts,nerr)
        call getfhv('USER5   ',user5,nerr)
        call getfhv('USER9   ',user9,nerr)
        call getfhv('EVLA    ',evla,nerr)
        call getfhv('EVLO    ',evlo,nerr)
        call getfhv('STLA    ',stla,nerr)
        call getfhv('STLO    ',stlo,nerr)
        call getkhv('KSTNM   ',kstnm,nerr)
        call getkhv('KCMPNM  ',kcmpnm,nerr)
        call getnhv('NZYEAR', nzyear,nerr)
        call getnhv('NZJDAY', nzjday,nerr)
        call getnhv('NZHOUR', nzhour,nerr)
        call getnhv('NZMIN', nzmin,nerr)
        call getnhv('NZSEC', nzsec,nerr)
C       write(6,*)'az,baz,dist,deg,dt,npts,b,o,nerr,tp,ts',
C     1     az,baz,dist,deg,dt,npts,beg,origtime,nerr,tp,ts
C       write(6,*)'evla,evlo,stla,stlo,kstnm,kcmpnm',
C     1     evla,evlo,stla,stlo,kstnm,kcmpnm
*
        if(nerr .eq. 0 .and. origtime .ne. -12345)then
            t0 = beg - origtime
        else
            t0 = beg
        end if
        ttp = tp - origtime
*
C       write(6,*)'ttp=',ttp
C       write(6,*)'dist=',dist
        sta = kstnm(1:8)
        comp = kcmpnm(1:8)
C       write(6,*)'name,npts,dist,deg,az,baz,t0,dt,sta',
C     1     name,npts,dist,deg,az,baz,t0,dt,sta
        cdate = ' '
*
*
        return
        end

        subroutine dmxmn(seis,n1,n2,ampmx,ampmn)
        real seis(*)
        integer n1, n2
        ampmx = -1.0e+38
        ampmn =  1.0e+38
        do 1000 i=n1,n2
            if(seis(i).gt.ampmx)ampmx=seis(i)
            if(seis(i).lt.ampmn)ampmn=seis(i)
 1000   continue
        return 
        end

        subroutine gcmdln(x0,y0,xlen,ylen,nmin,nmax,MXPTS,ybot,
     1      iplmn,kolor,fname,nfile,MXFILE,doabs,pcy,doover,doamp,
     1      dotmtp, dotmbt,dosctp, doscbt,doovershd,douser9,
     2      tmin,tmax,dotmnmx)
c-----
c      parse the command line arguments
c-----
c      x0  R*4 - lower left corner of plot frame
c      y0  R*4 - lower left corner of plot frame
c      xlen    R*4 - width  of plot frame on page
c      ylen    R*4 - height of plot frame on page
c      dokolor K   - for multi file plots, use different colors
c                  and annotate the plot
c      nmin    R*4 - minimum time index for plot
c      nmax    R*4 - maximum time index for plot
c      MXPTS   I*4 - maximum numebr of time points permitted 
c      ybot    R*4 - of ymax < ybot, ymax set to ybot
c      iplmn   I*4 - 0 no shading
c                1 shade positive amplitudes
c               -1 shade negative amplitudes
c      fname   C*80    - character array of SAC file names
c      nfile   I*4 - actual number of SAC files
c      MXFILE  I*4 - maximum number of SAC files permitted
c      doabs   L   - .true. scale according to maximum ofr traces
c      pcy R*4 - percent of DY for trace amplitude
c      doover  L   - .true. overlay and do not space in Y direction
c      doamp   L   - .true. annotate plotwith amplitude at far left
c      dotmtp  L   - .true. time scale at top
c      dosctp  L       - .true. time scale with caption at top
c      dotmbt  L       - .true. time scale at bottom
c      doscbt  L       - .true. time scale with caption at bottom
c      doovershd   L   - .true. Use shading when overlaying KOLOR 
c                  then black
c      douser9 L   - .true. plot time shift in USER9 header
c      tmin    R   - minimum value for time scale
c      tmax    R   - maximum value for time scale
c      dotmnmx L   - if true use tmin, tmax pair
c-----
        real*4 x0, y0, xlen, ylen
        logical dokolor
        integer nmin, nmax
        integer MXPTS
        integer nfile
        integer MXFILE
        character*80 fname(MXFILE)
        real tmin, tmax
        logical dotmnmx
        logical doabs
        logical doover
        logical doamp
        logical dotmtp, dotmbt
        logical dosctp, doscbt
        logical doovershd
        logical douser9

        integer mnmarg
        character name*80
        integer i

        x0 = 2.0
        y0 = 7.0
        xlen = 6.0
        ylen = 6.0
        iplmn = 0
        nmin = 1
        nmax = MXPTS
        kolor = 1
        nfile = 0
        dokolor = .true.
        doabs = .false.
        pcy = 0.6
        ybot = 0.0
        doover = .false.
        doamp = .false.
        dotmtp = .false.
        dotmbt = .false.
        dosctp = .false.
        doscbt = .false.
        doovershd = .false.
        douser9 = .false.
c-----
c       totally unreasonable numbers
c-----
        tmin = -1.0e+38
        tmax = -1.0e+38
        dotmnmx = .false.
        
        nmarg = mnmarg()
        i = 0
 1000   continue
            i = i + 1
            if(i.gt.nmarg)go to 2000
            call mgtarg(i,name)
            if(name(1:2).eq.'-K')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(i10)')kolor
            else if(name(1:3).eq.'-X0')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,f10.0)')x0
            else if(name(1:3).eq.'-Y0')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,f10.0)')y0
            else if(name(1:5).eq.'-XLEN')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,f10.0)')xlen
            else if(name(1:5).eq.'-YLEN')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,f10.0)')ylen
            else if(name(1:5).eq.'-YBOT')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,f10.0)')ybot
            else if(name(1:4).eq.'-PCY')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,f10.0)')pcy
            else if(name(1:5).eq.'-NMIN')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(i10)')nmin
            else if(name(1:5).eq.'-NMAX')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(i10)')nmax
            else if(name(1:6).eq.'-DOAMP')then
                doamp = .true.
            else if(name(1:5).eq.'-DOPL')then
                iplmn =  1
            else if(name(1:5).eq.'-DOMI')then
                iplmn = -1
            else if(name(1:2) .eq. '-O'.and.
     1          name(1:3).ne.'-OS')then
                doover = .true.
            else if(name(1:3).eq.'-OS')then
                doover = .true.
                doovershd = .true.
            else if(name(1:2) .eq. '-?')then
                call usage(' ')
            else if(name(1:2) .eq. '-h')then
                call usage(' ')
            else if(name(1:4).eq.'-ABS')then
                doabs = .true.
            else if(name(1:5).eq.'-TSCT')then
                dosctp = .true.
                dotmtp = .true.
            else if(name(1:5).eq.'-TSCB')then
                doscbt = .true.
                dotmbt = .true.
            else if(name(1:5).eq.'-TTMT')then
                dotmtp = .true.
            else if(name(1:5).eq.'-TTMB')then
                dotmbt = .true.
            else if(name(1:6).eq.'-USER9')then
                douser9 = .true.
            else if(name(1:5).eq.'-TMIN')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,f10.0)')tmin
            else if(name(1:5).eq.'-TMAX')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,f10.0)')tmax
            else
                nfile = nfile + 1
                fname(nfile) = name
            endif
        go to 1000
 2000   continue
c-----
c      safety
c----
        if(nmin.gt.nmax)then
            n = nmin
            nmin = nmax
            nmax = n
        endif
        if(nmin.lt.1)nmin = 1
        if(nmax.lt.1 .or.nmax.gt.MXPTS)nmax=MXPTS
        if(nmin.gt.nmax)then
            n = nmin
            nmin = nmax
            nmax = n
        endif
        if(tmin .gt. -1.0e+38 .and. tmax .gt. -1.0e+38)then
            dotmnmx = .true.
            if(tmax .lt. tmin)then
               tmp = tmax
               tmax = tmin
               tmin  = tmp
            endif
        else
            dotmnmx = .false.
        endif

        return
        end

        subroutine usage(emsg)
        implicit none
        character emsg*(*)
        integer LER
        parameter (LER=0)
        integer lgstr
        integer ls
        ls = lgstr(emsg)
        if(emsg .ne. ' ')write(LER,*)emsg(1:ls)
        write(LER,*)'Usage: pltsac -XLEN xlen -YLEN ylen',
     1      ' -X0 x0 -Y0 y0 ',
     2      '-nmin NMIN -NMAX nmax -ABS -PCY pcy -YBOT ybot ',
     4      '-TSCT -TSCB -TTMT -TTMB -USER9 -TMIN tmin -TMAX tmax',
     3      '-KOLOR IPEN -DOAMP  -DOPLUS -DOMINUS -O -OS sac_files'
        write(LER,*)'Plot SAC files'
        write(LER,*)
     1  '-XLEN xlen (default 6.0  ) Length X-axis'
        write(LER,*)
     1  '-YLEN ylen (default 6.0  ) Length Y-axis'
        write(LER,*)
     1  '-NMIN nmin (default 1    )  First point to plot'
        write(LER,*)
     1  '-NMAX nmax (default 64000)  last point to plot'
        write(LER,*)
     1  '-PCY pcy   (default 0.6  ) Fraction of dy for trace amplitude'
        write(LER,*)
     1  '-YBOT ybot (default 0.0) if ymax < ybot trace is zero'
        write(LER,*)
     1  '-ABS       (default false) plot absolute amplitudes'
        write(LER,*)
     1  '-X0 x0     (default  2.0 )  x-position of lower left corner'
        write(LER,*)
     1  '-Y0 y0     (default  7.0 )  y-position of lower left corner'
        write(LER,*)
     1  '-K PEN     (default  1)  Use color for '
        write(LER,*)
     1  '           if kolor < 0 use red->blue progression'
        write(LER,*)
     1  '-DOAMP     (default false)  Annotate with amplitude value '
        write(LER,*)
     1  '-DOPLUS    (default false)  Shade positive values '
        write(LER,*)
     1  '-DOMINUS   (default false)  Shade negative values '
        write(LER,*)
     1  '-O         (default false)  Overlay traces no y space '
        write(LER,*)
     1  '-OS        (default false)  Overlay traces and shade '
        write(LER,*)
     1  'sac_files                   Sac binary trace files'
        write(LER,*)
     1  '-TSCT      (default false)  Put time scale/label at top'
        write(LER,*)
     1  '-TSCB      (default false)  Put time scale/label at bottom'
        write(LER,*)
     1  '-TTMT      (default false)  Put time scale grid at top'
        write(LER,*)
     1  '-TTMB      (default false)  Put time scale grid at bottom'
        write(LER,*)
     1  '-TMIN tmin (default from file) minimum time for time scale'
        write(LER,*)
     1  '-TMAX tmin (default from file) maximum time for time scale'
        write(LER,*)
     1  '-USER9     (default false)  Show wvfgrd96 time shift'
        write(LER,*)
     1  '-?        (default none )  this help message '
        write(LER,*)
     1  '-h        (default none )  this help message '
        stop
        end
        
