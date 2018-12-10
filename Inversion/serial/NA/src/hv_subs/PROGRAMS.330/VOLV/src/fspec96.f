        program fspec96
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME V                                                       c
c                                                                     c
c      PROGRAM: FSPEC96                                               c
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
c       19 NOV 2000 added newest dolinx from grphsubf.f
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
        integer kolor
        logical dofreq, lxislg, lyislg, lxset, lyset
        real xmin, xmax, ymin, ymax
c-----
c       internal variables
c-----
        parameter (LER=6, LIN=5, LOT=6)
        parameter (NSAMP=16384)
        real x(NSAMP),y(NSAMP)
        complex z(NSAMP)
        character svcomp*8
        real*8 sepoch, eepoch, o
c-----
c       graphics control
c-----
        logical ticup, labtop, dopow
        logical ticlft, lablft
        real sizex, sizey

        integer lgstr

        parameter (NUNIT=13)
        character*16 ctunit(NUNIT)
        character*16 csunit(NUNIT)
        character*40 outstr
        data ctunit/ 'COUNTS          ',
     1  'cm              ', 'cm/sec          ', 'cm/sec/sec      ',
     2  'm               ', 'm/sec           ', 'm/sec/sec       ',
     3  'micron          ', 'micron/sec      ', 'micron/sec/sec  ',
     4  'Pa              ', 'MPa             ', '1/sec           '/
        data csunit/ 'COUNTS          ',
     1  'cm-sec          ', 'cm              ', 'cm/sec          ',
     2  'm-sec           ', 'm               ', 'm/sec           ',
     3  'micron-sec      ', 'micron          ', 'micron/sec      ',
     4  'Pa-sec          ', 'MPa-sec         ', 'none            '/

c-----
c       call machine dependent initialization
c-----
        call mchdep()
c-----
c       parse command line arguments
c-----
        call gcmdln(kolor,dofreq,lxislg, lyislg, lxset, lyset,
     1      xminsv,xmaxsv,yminsv,ymaxsv,nx,ny)
c-----
c       graphic initialization
c-----
        call pinitf('FSPEC96.PLT')
c-----
c       read in data, and plot all traces of a given 'station'
c       on the same page
c-----
c-----
c       determine which Greens functions computed
c-----
c-----
c       process Green's Functions
c-----
        ipage = 0
c-----
c       basic graphics control
c-----
        xaxlen = 5.0
        yaxlen = 5.0
        x0 = 2.0
        y0 = 1.0
        yy0 = 6.5
        yyaxlen = 1.0
  100   continue
            call rdhd96(LIN,nerr)
            if(nerr .lt. 0)go to 9999
            if(ipage .eq.0)then
                ipage = 1
            else
                call frame()
            endif
c-----
c       only plot upto iftype traces
c-----
            ndone = 0
            do 200 j=1,16
                if(jsrc(j) .ne. 0)then
                    do 300 i=1,NSAMP
                        y(i)=0.0
  300               continue
                    call rdtr96(LIN,stcomp,cmpinc,cmpaz,
     1                  cmpdt, npts, ksyear, ksmon, 
     2                  ksday, kshour, ksmin, ssec, 
     3                  y,nerr,NSAMP)
                    if(nerr .lt. 0 )go to 9999
                    ndone = ndone + 1
                    if(ndone.gt.iftype)go to 200
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
                    svcomp = stcomp
                    nn = npts
                    call npow2(nn)
                    do 1101 i=1,nn
                        x(i) = ti + (i-1)*dt
                        if(i.le.npts)then
                            z(i) = cmplx(y(i),0.0)
                        else
                            z(i) = cmplx(0.0,0.0)
                        endif
 1101               continue
                    call domxmn(x,npts,xxmax,xxmin)
                    call domxmn(y,npts,yymax,yymin)
c------
c       force scale to be the same, e.g., zero is at yy0+0.5*yyaxlen
c-----
                    if(abs(yymax).gt.abs(yymin))then
                        tmp = abs(yymax)
                        yymax = tmp
                        yymin = - yymax
                    else
                        tmp = abs(yymin)
                        yymax = tmp
                        yymin = - yymax
                    endif
                    sizex = 0.07
                    ticup =.false.
                    labtop = .false.
                    dopow = .true.
                    outstr = ' '
                    llt = lgstr(outstr)
                    call dolinx(x0,yy0-0.10,xaxlen,
     1                  xxmax,xxmin,sizex,
     2                  ticup,labtop,dopow,
     3                  llt,outstr)
                    sizey = 0.07
                    ticlft =.true.
                    lablft = .true.
                    dopow = .true.
                    llt = lgstr(ctunit(iunit))
                    call doliny(x0-0.10,yy0,yyaxlen,
     1                  yymax,yymin,sizey,
     2                  ticlft,lablft,dopow,
     3                  llt,ctunit(iunit))
                    call xyplot(x,y,npts,x0,yy0,xaxlen,
     1                  yyaxlen,kolor,.false.,.false.,
     2                  xxmin,xxmax,yymin,yymax,.false.,
     3                  .true.)

c-----
c       get the spectra and plot it
c-----
                    call zfour(z,nn,-1,dt,df)
                    n21 = nn/2 + 1
                    do 1102 i=1,n21
                        x(i) = (i-1)*df
                        y(i) = cabs(z(i))
 1102               continue
                    xmin = xminsv
                    ymin = yminsv
                    xmax = xmaxsv
                    ymax = ymaxsv
                    if(.not.lxset)then
                        if(lxislg)then
                        call domxmn(x(2),n21 -1 ,xmax,xmin)
                        else
                        call domxmn(x(1),n21,xmax,xmin)
                        endif
                        if(.not.dofreq)then
                            tmp = 1./xmax
                            xmax = 1./xmin
                            xmin = tmp
                        endif
                        if(nx.gt.0)then
                            xmin = xmax/10.0**nx
                        endif
                    endif
                    if(.not.lyset)then
                        call domxmn(y,n21,ymax,ymin)
                        if(ny.gt.0)then
                            ymin = ymax/10.0**ny
                        endif
                    endif
c-----
c       safety
c-----
                    if(lxislg .and. xmin.eq.0.0)then
                        xmin = xmax/10000.
                    endif
                    if(lyislg .and. ymin.eq.0.0)then
                        ymin = ymax/10000.
                    endif
                    sizex = 0.10
                    naxdig = 3
                    if(dofreq)then
                        outstr='Frequency (Hz)'
                    else
                        outstr='Period (sec)'
                    endif
                    ls = lgstr(outstr)
                    if(lxislg)then
                        ticup =.true.
                        labtop = .false.
                        dopow = .true.
                        call dologx(x0,y0,xaxlen,
     1                      xmax,xmin,sizex,
     2                      ticup,labtop,dopow,
     3                      ls,outstr)
                        ticup =.false.
                        labtop = .false.
                        dopow = .false.
                        call dologx(x0,y0+yaxlen,xaxlen,
     1                      xmax,xmin,sizex,
     2                      ticup,labtop,dopow,
     3                      ls,outstr)
                    else
                        ticup =.true.
                        labtop = .false.
                        dopow = .true.
                        call dolinx(x0,y0,xaxlen,
     1                      xmax,xmin,sizex,
     2                      ticup,labtop,dopow,
     3                      ls,outstr)
                        ticup =.false.
                        labtop = .false.
                        dopow = .false.
                        call dolinx(x0,y0+yaxlen,xaxlen,
     1                      xmax,xmin,sizex,
     2                      ticup,labtop,dopow,
     3                      ls,outstr)
                    endif
                    sizey = 0.10
                    naxdig = 3
                    outstr = csunit(iunit)
                    ls = lgstr(outstr)
                    if(lyislg)then
                        ticlft =.false.
                        lablft = .true.
                        dopow = .true.
                        call dology(x0,y0,yaxlen,
     1                      ymax,ymin,sizey,
     2                      ticlft,lablft,dopow,
     3              10+ls,'Spectra ('//outstr(1:ls)//')')
                        ticlft =.true.
                        lablft = .false.
                        dopow = .false.
                        call dology(x0+xaxlen,y0,yaxlen,
     1                      ymax,ymin,sizey,
     2                      ticlft,lablft,dopow,
     3              10+ls,'Spectra ('//outstr(1:ls)//')')
                    else
                        ticlft =.false.
                        lablft = .true.
                        dopow = .true.
                        call doliny(x0,y0,yaxlen,
     1                      ymax,ymin,sizey,
     2                      ticlft,lablft,dopow,
     3              10+ls,'Spectra ('//outstr(1:ls)//')')
                        ticlft =.true.
                        lablft = .false.
                        dopow = .false.
                        call doliny(x0+xaxlen,y0,yaxlen,
     1                      ymax,ymin,sizey,
     2                      ticlft,lablft,dopow,
     3              10+ls,'Spectra ('//outstr(1:ls)//')')
                    endif
                    if(lxislg)then
                    call xyplot(x(2),y(2),n21 -1,x0,y0,xaxlen,
     1                  yaxlen,kolor,lxislg,lyislg,
     2                  xmin,xmax,ymin,ymax,.true.,
     3                  dofreq)
                    else
                    call xyplot(x,y,n21,x0,y0,xaxlen,
     1                  yaxlen,kolor,lxislg,lyislg,
     2                  xmin,xmax,ymin,ymax,.true.,
     3                  dofreq)
                    endif
c-----
c       put on plot labels
c-----
                l = lgstr(ccomnt)
                if(l.gt.1)then
                call symbol(1.0,0.30,0.10,ccomnt(1:l),0.0,l)
                endif
                endif
  200       continue
        goto 100
 9999   continue
        call pend()
        end

        subroutine xyplot(x,y,n,x0,y0,xaxlen,yaxlen,kolor,
     1      lxislg,lyislg,xmin,xmax,ymin,ymax,dobox,dofreq)

c-----
c       x() R   - array of x-values
c       y() R   - array of y-values
c       n   I   - number of points to plot
c       x0,y0   R   - lower left corner of plot
c       xaxlen  R   - x-axis length
c       yaxlen  R   - y-axis length
c       kolor   I   - pen color
c       lxislg  L   - .true. x-axis is logarithmic
c       lyislg  L   - .true. y-axis is logarithmic
c       xmin    R   - minimum value of x-axis
c       xmax    R   - maximum value of x-axis
c       ymin    R   - minimum value of y-axis
c       ymax    R   - maximum value of y-axis
c       dobox   L   - .true. put box all around
c                 .false. no box
c       dofreq  L   - .false. plot y vs 1/x
c-----
        real x(n), y(n)
        integer n
        real x0,y0,xaxlen,yaxlen,xmin,xmax,ymin,ymax
        logical lxislg, lyislg, dobox, dofreq

        xlow = x0
        ylow = y0
        xhgh = x0 + xaxlen
        yhgh = y0 + yaxlen
        if(dobox)call gbox(xlow,ylow,xhgh,yhgh)
        call gclip('ON',xlow,ylow,xhgh,yhgh)
c-----
c       plot the values
c       be very careful for small values
c-----
        call newpen(kolor)
        do 1000 i=1,n
            if(lyislg)then
                if(y(i).le.0.0)then
                    yy = -100000.
                else
                    yy = y0 + yaxlen*alog10(y(i)/ymin)/
     1                  alog10(ymax/ymin)
                endif
            else
                yy = y0 + yaxlen*(y(i) - ymin)/(ymax - ymin)
            endif
            if(lxislg)then
                if(x(i).le.0.0)then
                    xx = -100000.
                else
                    if(dofreq)then
                        xval = x(i)
                    else
                        xval = 1.0/x(i)
                    endif
                    xx = x0 + xaxlen*alog10(xval/xmin)/
     1                  alog10(xmax/xmin)
                endif
            else
                if(dofreq)then
                    xval = x(i)
                else
                    xval = 1.0/x(i)
                endif
                xx = x0 + xaxlen*(xval - xmin)/(xmax - xmin)
            endif
            if(i.eq.1)then
                call plot(xx,yy,3)
            else
                call plot(xx,yy,2)
            endif
 1000   continue
        call gclip('OFF',xlow,ylow,xhgh,yhgh)
        call newpen(1)
        return
        end
        

        subroutine gcmdln(kolor,dofreq,lxislg, lyislg, lxset, lyset,
     1      xmin,xmax,ymin,ymax,nx,ny)
c-----
c       parse command line arguments
c
c       requires subroutine targ() and funtion mnmarg()
c
c-----
c-----
c       kolor   - I*4   CALPLOT pen number default = 1
c       dofreq  - L .true. x-axis is frequency
c       lxislg  - L .true. x-axis is linear
c       lyislg  - L .true. y-axis is linear
c       lxset   - L .true. axis limits are given, else
c                   determined from data
c       lyset   - L .true. axis limits are given, else
c                   determined from data
c       xmin    - R minimum value of x-axis
c       xmax    - R maximum value of x-axis
c       ymin    - R minimum value of y-axis
c       ymax    - R maximum value of y-axis
c       nx  - I for log plot, maximum number of cycles
c       ny  - I for log plot, maximum number of cycles
c-----
        character*20 name
        integer kolor
        logical dofreq, lxislg, lyislg, lxset, lyset
        real xmin, xmax, ymin, ymax
        integer nocx, nocy
        integer mnmarg

        kolor = 1
        lxislg = .false.
        lyislg = .false.
        lxset = .false.
        lyset = .false.
        xmin =  1.0e+38
        xmax = -1.0e+38
        ymin =  1.0e+38
        ymax = -1.0e+38
        defmax = -1.0e+38
        defmin =  1.0e+38
        nocx = -1
        nocy = -1
        nmarg = mnmarg()
        dofreq = .true.
        i = 0
   11   i = i + 1
        if(i.gt.nmarg)go to 13
            call mgtarg(i,name)
            if(name(1:2).eq.'-K')then
                i=i+1
                call mgtarg(i,name)
                read(name,'(bn,i10)')kolor
                if(kolor.lt.1)kolor = 1
            else if(name(1:3).eq.'-NX')then
                i=i+1
                call mgtarg(i,name)
                read(name,'(bn,i10)')nx
            else if(name(1:3).eq.'-NY')then
                i=i+1
                call mgtarg(i,name)
                read(name,'(bn,i10)')ny
            else if(name(1:2).eq.'-?')then
                call usage()
            else if(name(1:2).eq.'-h')then
                call usage()
            else if (name(1:5).eq.'-XLOG')then
                lxislg = .true.
            else if (name(1:5).eq.'-XLIN')then
                lxislg = .false.
            else if (name(1:5).eq.'-YLOG')then
                lyislg = .true.
            else if (name(1:5).eq.'-YLIN')then
                lyislg = .false.
            else if (name(1:5).eq.'-XMIN')then
                i = i + 1
                call mgtarg(i,name)
                call chtofp(name,xmin)
            else if (name(1:5).eq.'-XMAX')then
                i = i + 1
                call mgtarg(i,name)
                call chtofp(name,xmax)
            else if (name(1:5).eq.'-YMIN')then
                i = i + 1
                call mgtarg(i,name)
                call chtofp(name,ymin)
            else if (name(1:5).eq.'-YMAX')then
                i = i + 1
                call mgtarg(i,name)
                call chtofp(name,ymax)
            elseif (name(1:4).eq.'-PER' .or.
     1          name(1:4).eq.'-per')then
                dofreq = .false.
            elseif (name(1:5).eq.'-FREQ' .or.
     1          name(1:5).eq.'-freq')then
                dofreq = .true.
            endif
        go to 11
   13   continue
c-----
c       do quality control
c-----
        if(xmin.ne. defmin .and. xmax.ne.defmax .and. xmin.gt.xmax)then
            tmp = xmin
            xmin = xmax
            xmax = tmp
        endif
        if(ymin.ne. defmin .and. ymax.ne.defmax .and. ymin.gt.ymax)then
            tmp = ymin
            ymin = ymax
            ymax = tmp
        endif
        if(xmin.eq.defmin .or. xmax.eq.defmax)then
            lxset = .false.
        else
            lxset = .true.
        endif
        if(ymin.eq.defmin .or. ymax.eq.defmax)then
            lyset = .false.
        else
            lyset = .true.
        endif
        if(xmin.eq.defmin .or. xmax.eq.defmax)then
            lxset = .false.
        else
            lxset = .true.
        endif
        if(ymin.eq.defmin .or. ymax.eq.defmax)then
            lyset = .false.
        else
            lyset = .true.
        endif
        return
        end

        subroutine usage
        integer*4 LER, LIN, LOT
        parameter(LER=6, LIN=5, LOT=6)
        write(LER,*)'USAGE: fspec96  ',
     1      ' [-K kolor] [-XLOG] [-XLIN] [-YLIN] [-YLOG]',
     2      ' [-XMIN xmin -XMAX xmax] [-YMIN ymin -YMAX ymax]',     
     3      ' [-NX nx] [-NY ny] [-FREQ -PER] [-?] [-h]'
        write(LER,*)
     1      ' -K kolor    (default 1)  pen color'
        write(LER,*)
     1      ' -XLOG     (default .false.) '
        write(LER,*)
     1      ' -XLIN     (default .true.) '
        write(LER,*)
     1      ' -YLOG     (default .false.) '
        write(LER,*)
     1      ' -YLIN     (default .true.) '
        write(LER,*)     
     1      ' -XMIN xmin  (default auto) minimum value of x-axis'
        write(LER,*)     
     1      ' -XMAX xmax  (default auto) maximum value of x-axis'
        write(LER,*)     
     1      ' -YMIN ymin  (default auto) minimum value of y-axis'
        write(LER,*)     
     1      ' -YMAX ymax  (default auto) maximum value of y-axis'
        write(LER,*)     
     1      ' -NX nx  (default auto ) maximum cyles in log x-axis' 
        write(LER,*)     
     1      ' -NY ny  (default auto ) maximum cyles in log y-axis' 
        write(LER,*)
     1      ' -FREQ   (default .true.) x-axis is frequency'
        write(LER,*)
     1      ' -PER   (default .false.) x-axis is period'
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
        real*4 depmax,depmin
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
