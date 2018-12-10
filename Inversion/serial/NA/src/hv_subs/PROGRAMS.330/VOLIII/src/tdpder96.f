      program tdpder96
c----------------------------------------------------------------------c
c                                                                      c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                 c
c      VOLUME III                                                      c
c                                                                      c
c      PROGRAM: TDPDER96                                               c
c                                                                      c
c      COPYRIGHT 1996 ,2010                                            c
c      R. B. Herrmann, Chien-Ying Wang                                 c
c      Department of Earth and Atmospheric Sciences                    c
c      Saint Louis University                                          c
c      221 North Grand Boulevard                                       c
c      St. Louis, Missouri 63103                                       c
c      U. S. A.                                                        c
c                                                                      c
c----------------------------------------------------------------------c
c Revision history:
c       27 DEC 2000 - removed commented out code 
c       22 JUL 2002 - TXT files use MODE = 0 for fundamental
c       16 OCT 2006 - fixed depth scale
c----------------------------------------------------------------------c
        parameter (LIN=5,LOT=6,LER=0)

        parameter (NL=1000)
        common/tmodl/d(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      Rho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL),
     3      frefp(NL), frefs(NL)
        real d, ta, tc, tl, tn, tf, rho, qa, qb, etap,etas,frefp, frefs


        common/ctrl/ tmin,tmax,cmin,cmax

        character ofile*20, dfile*20
        logical ext
c-----
c       command line arguments
c-----
        logical dobox, dotext, doclean
        integer*4 ilorr
        real*4 xmin, xmax, ymin, ymax, x0, y0, xlen, ylen
        integer*4 kolor
c-----
c       machine dependent initialization
c-----
        call mchdep()
c-----
c       get control parameters from command line
c-----
        call gcmdln(ilorr, dobox, xmin, xmax, ymin, ymax, 
     1      x0, y0, xlen, ylen, 
     2      kolor,dotext,doclean)
        if(ilorr.lt.1)call usage('specify -L or -R')
c------
c     get control parameters and open input file.
c------
        pi=2.*3.141592653
        if(ilorr.eq.1)then
            dfile = 'tlegn96.der'
            ofile = 'TLDER.TXT'
        else if(ilorr.eq.2)then
            dfile = 'tregn96.der'
            ofile = 'TRDER.TXT'
        endif
        inquire(file=dfile,exist=ext)
        if(.not.ext)then
            call usage('Eigenfunction file '//dfile//
     1      ' does not exist')
        endif
        open(1,file=dfile,status='unknown',
     1      access='sequential',form='unformatted')
        rewind 1
        

        if(dotext)then
            call outtxt(ofile,ilorr)
        endif

        call outplt( dobox, xmin, xmax, ymin, ymax, 
     1      x0, y0, xlen, ylen, 
     2      kolor,ilorr,doclean)

        close (1)
        end

        subroutine  outplt( dobox, xmin, xmax, ymin, ymax, 
     1      xx0, yy0, xxlen, yylen, 
     2      kolor,ilorr,doclean)
c-----
c       plot the dispersion curves
c
c       ilorr   I*4 1 = Love, 2 = Rayleigh
c       dobox   L   .true. put frame around plot
c       xmin    R*4 minimum value of independent variable
c       xmax    R*4 maximum value of independent variable
c       ymin    R*4 minimum value of   dependent variable
c       ymax    R*4 maximum value of   dependent variable
c       xx0 R*4 
c       yy0 R*4 lower left coordinates of plot box
c       xxlen   R*4 length of X axis in inches
c       yylen   R*4 length of Y axis in inches
c       kolor   I*4 CALPLOT pen value
c       ilorr   I*4 1 = Love
c               2 = Rayleigh
c-----
        integer*4 ilorr
        logical dobox, doxlog, doylog
        real*4 xmin, xmax, ymin, ymax
        real*4 xx0, yy0, xxlen, yylen
        integer*4 kolor
        logical doclean

        
        parameter(LIN=5, LOT=6, LER=0)

        parameter (NL=1000)
        common/tmodl/d(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      Rho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL),
     3      frefp(NL), frefs(NL)
        real d, ta, tc, tl, tn, tf, rho, qa, qb, etap,etas,frefp, frefs


        integer*4 nmax,nper,ifunc,kmode


        real*4 dcdh(NL), dcdav(NL), dcdah(NL), dcdbv(NL),
     1     dcdbh(NL), dcdn(NL), dcdr(NL)
        real*4 ur(NL), tur(NL), uz(NL), tuz(NL)

        real*4 z(NL)


        logical dopowy
        common/ctrl/ tmin,tmax,cmin,cmax

        character mname*80
        real*4 fpar(20)
        integer*4 ipar(20)

        doxlog = .false.
        doylog = .false.
c-----
c           begin reading the output file of sdisp96 
c           and do integrity check
c-----
            rewind 1
            call getmdt(1,nmax,d,ta,tc,tf,
     1           tl,tn,rho,qa,qb,nper,ddpths,depthr,
     1           mname,ipar,fpar)
            call gethed(1,ifunc,kmode,t0,ierr)
            if(ierr.eq.200)then
                write(LOT,*)'End of File reading header'
                stop
            endif
            if(ierr.eq.1001)then
                write(LOT,*)'Error reading header'
                stop
            endif
            if(ifunc.eq.5 .and. ilorr.ne.1)then
                write(LOT,*)'Data file is not for Love waves'
                stop
            endif
            if(ifunc.eq.6 .and. ilorr.ne.2)then
                write(LOT,*)'Data file is not for Rayleigh wave'
                stop
            endif
            
            write(LOT,*) ' '
            write(LOT,*) 'Total period=',nper
            write(LOT,*) ' '
            write(LOT,*)
     *      'Number of modes at the lowest period:',t0,' = ',kmode
c-----
c       print the dispersion values.
c       Since the data files are all modes for each frequency, we will
c       repeatedly rwind the input to make a list of all frequencies
c       for each mode
c-----
            refdep = fpar(1)
            z(1) = -refdep
            do 100 i=2,nmax
                z(i) = z(i-1) + d(i-1)
  100       continue
c-----
c           now start processing from the beginning
c-----
        if(ilorr.eq.1)then
            call pinitf('TLDER.PLT')
        else if(ilorr.eq.2)then
            call pinitf('TRDER.PLT')
        endif
            rewind 1
            call getmdt(1,nmax,d,ta,tc,tf,
     1           tl,tn,rho,qa,qb,nper,dphs,dphr,
     1           mname,ipar,fpar)
c-----
c       set up counters
c-----
c------
c     get the eigenfunction information
c     pick up the phase velocities at different frequencies for
c     a particular mode.
c------
        lorr = ilorr + 4
  200   continue
            call gethed(1,ifunc,nmode,t0,ierr)
            omega = 6.2831853/t0
            if(ierr.eq.1001)go to 2001
            if(ierr.eq.200)then
                write(LOT,*) 
     1      'No ifunc=-1. Data end at perd=',t0
                ifunc=-1
            endif
            if(ifunc.le.0) go to 1000
            if(nmode.le.0) go to 200
            do 300 j=1,nmode
            call getdrt(1,lorr,wvno,u,gammal,
     1      sur,sdur,sd2ur,suz,sduz,sd2uz,sare,wvnsrc,sur0,
     2      rur,rtr,ruz,rtz,rare,wvnrec,rur0,
     3      sumkr,sumgr,sumgv,ierr,nmax,dcdh,
     4      dcdav,dcdah, dcdbv, dcdbh, dcdn,dcdr,
     5      ur,tur,uz,tuz,ipar)


                if(ierr.eq.200)go to 2001
                if(ierr.eq.1001)go to 1001
                f0=1./t0
                c = omega/wvno
c-----
c       systematically plot the data
c-----
                xlen =  xxlen
                ylen = -yylen
                x0 = xx0
                y0 = yy0 + yylen
                dopowy = .true.
                if(ilorr.eq.1)then
                    if(ipar(4).eq.1)then
                       call doit1(z,ur,nmax,5,
     1                  'UT','DEPTH',   
     1                  xlen,ylen,doxlog,doylog,dopowy,
     2                  x0,y0,c,t0,j,dobox,kolor,doclean)
                       call doit1(z,tur,nmax,1,
     1                  'TT','DEPTH',   
     1                  xlen,ylen,doxlog,doylog,dopowy,
     2                  x0,y0,c,t0,j,dobox,kolor,doclean)
                    endif
                    if(ipar(5).eq.1)then
                       call doit1(z,dcdh,nmax,2,
     1                  'DC/DH','DEPTH',    
     1                  xlen,ylen,doxlog,doylog,dopowy,
     2                  x0,y0,c,t0,j,dobox,kolor,doclean)
                    endif
                    if(ipar(7).eq.1)then
                       call doit1(z,dcdbv,nmax,2,
     1                  'DC/DBv','DEPTH',    
     1                  xlen,ylen,doxlog,doylog,dopowy,
     2                  x0,y0,c,t0,j,dobox,kolor,doclean)
                    endif
                    if(ipar(8).eq.1)then
                       call doit1(z,dcdr,nmax,2,
     1                  'DC/DR','DEPTH',    
     1                  xlen,ylen,doxlog,doylog,dopowy,
     2                  x0,y0,c,t0,j,dobox,kolor,doclean)
                    endif
                    if(ipar(11).eq.1)then
                       call doit1(z,dcdbh,nmax,2,
     1                  'DC/DBh','DEPTH',    
     1                  xlen,ylen,doxlog,doylog,dopowy,
     2                  x0,y0,c,t0,j,dobox,kolor,doclean)
                    endif
                else if(ilorr.eq.2)then
                    if(ipar(4).eq.1)then
                       call doit2(z,ur,tuz,nmax,4,
     1                  'UR','DEPTH',   
     1                  xlen,ylen,doxlog,doylog,dopowy,
     2                  x0,y0,c,t0,j,dobox,kolor,
     3                  omega,wvno,doclean)
                       call doit1(z,uz,nmax,1,
     1                  'UZ','DEPTH',   
     1                  xlen,ylen,doxlog,doylog,dopowy,
     2                  x0,y0,c,t0,j,dobox,kolor,doclean)
                       call doit1(z,tur,nmax,1,
     1                  'TR','DEPTH',   
     1                  xlen,ylen,doxlog,doylog,dopowy,
     2                  x0,y0,c,t0,j,dobox,kolor,doclean)
                       call doit1(z,tuz,nmax,1,
     1                  'TZ','DEPTH',   
     1                  xlen,ylen,doxlog,doylog,dopowy,
     2                  x0,y0,c,t0,j,dobox,kolor,doclean)
                    endif
                    if(ipar(5).eq.1)then
                       call doit1(z,dcdh,nmax,2,
     1                  'DC/DH','DEPTH',    
     1                  xlen,ylen,doxlog,doylog,dopowy,
     2                  x0,y0,c,t0,j,dobox,kolor,doclean)
                    endif
                    if(ipar(6).eq.1)then
                       call doit1(z,dcdav,nmax,2,
     1                  'DC/DAv','DEPTH',    
     1                  xlen,ylen,doxlog,doylog,dopowy,
     2                  x0,y0,c,t0,j,dobox,kolor,doclean)
                    endif
                    if(ipar(7).eq.1)then
                       call doit1(z,dcdbv,nmax,2,
     1                  'DC/DBv','DEPTH',    
     1                  xlen,ylen,doxlog,doylog,dopowy,
     2                  x0,y0,c,t0,j,dobox,kolor,doclean)
                    endif
                    if(ipar(8).eq.1)then
                       call doit1(z,dcdr,nmax,2,
     1                  'DC/DR','DEPTH',    
     1                  xlen,ylen,doxlog,doylog,dopowy,
     2                  x0,y0,c,t0,j,dobox,kolor,doclean)
                    endif
                    if(ipar(9).eq.1)then
                       call doit1(z,dcdah,nmax,2,
     1                  'DC/DAh','DEPTH',    
     1                  xlen,ylen,doxlog,doylog,dopowy,
     2                  x0,y0,c,t0,j,dobox,kolor,doclean)
                    endif
                    if(ipar(20).eq.1)then
                       call doit1(z,dcdn,nmax,2,
     1                  'DC/DN','DEPTH',    
     1                  xlen,ylen,doxlog,doylog,dopowy,
     2                  x0,y0,c,t0,j,dobox,kolor,doclean)
                    endif
                    if(ipar(11).eq.1)then
                       call doit1(z,dcdbh,nmax,2,
     1                  'DC/DBh','DEPTH',    
     1                  xlen,ylen,doxlog,doylog,dopowy,
     2                  x0,y0,c,t0,j,dobox,kolor,doclean)
                    endif
                endif
  300       continue
        go to 200
 2001   continue
 1001       continue
        write(LOT,*) ' '
        write(LOT,*)
     *      '(l/r)egn96 data file terminates at  perd=',t0
        write(LOT,*) ' '
 1000   continue
        call pend()
        return
        end


        subroutine plotit(x,y,n,x0,y0,titlex,titley,
     1      c,t0,jmode,kolor,dobox,
     1      sxmin,sxmax,symin,symax,xlen,ylen,
     2      doxlog,doylog,dopowy,doclean)
c-----
c       plot the array
c
c       x   R*4     array of x-values
c       y   R*4     array of y-values
c       n   I*4 number of (x,y) pairs
c       x0  R*4 origin of plot
c       y0  R*4 origin
c       titlex  C*(*)   x-axis title
c       titley  C*(*)   y-axis title
c       c   R*4 phase velocity
c       t0  R*4 period
c       jmode   I*4 mode number
c       xmin    R*4 minimum x value
c       xmax    R*4 maximum x value
c       ymin    R*4 minimum y value
c       ymax    R*4 maximum y value
c       xlen    R*4 length of x-axis
c       ylen    R*4 length of y-ayis
c       doxlog  L   x-axis is logarithmic
c       doylog  L   y-ayis is logarithmic
c       dopowy  L   put in the numbers on y-axis
c       doclean L   .true. put in period, mode, phase velocity
c-----
        real*4 x(n), y(n)
        integer*4 n
        real*4 xmin, xmax, ymin, ymax, xlen, ylen
        logical doxlog, doylog
        logical dobox, dopowy, doclean
        character titlex*(*), titley*(*)
c-----
c       safety on extreme values
c-----
        if(sxmin.eq.sxmax)then
            xmin = sxmin
            xmax = xmin + 1.0
        else
            xmin = sxmin
            xmax = sxmax
        endif
        if(symin.eq.symax)then
            ymin = symin
            ymax = ymin + 1.0
        else
            ymin = symin
            ymax = symax
        endif


c-----
c       To ensure useful logarithmic plots, ensure that there
c       is at least one complete cycle
c-----
        if(doylog)then
            if( ymin .gt. 0.1*ymax)ymin = 0.99*ymax
        endif
        if(doxlog)then
            if( xmin .gt. 0.1*xmax)xmin = 0.99*xmax
        endif

        call start(dobox,x0,y0,xlen,ylen,
     1      xmin,xmax,ymin,ymax,
     1      titlex, titley,
     1      doxlog,doylog,dopowy)
c-----
c       annotate the base of the plot
c-----
        if(.not.doclean)then
            call symbol(x0,y0+ylen-0.25,0.10,'T0=',0.0,3)
            call number(999.0,999.0,0.10,t0,0.0,1003)
            call symbol(999.0,999.0,0.10,' MODE=',0.0,6)
            call number(999.0,999.0,0.10,real(jmode-1),0.0,-1)
            call symbol(999.0,999.0,0.10,' C=',0.0,3)
            call number(999.0,999.0,0.10,c,0.0,1003)
        endif

        call newpen(kolor)
        ipen = 3
        do 100 i=1,n
            xc1 = x(i)
            yc1 = y(i)
            if(doxlog)then
                 xx = x0 + xlen*alog10(xc1/xmin)
     1              /alog10(xmax/xmin)
            else
                 xx = x0 + xlen*(xc1 - xmin)/(xmax - xmin)
            endif
            if(doylog)then
                 yy = y0 + ylen*alog10(yc1/ymin)
     1              /alog10(ymax/ymin)
            else
                 yy = y0 + ylen*(yc1 - ymin)/(ymax - ymin)
            endif
            call plot(xx,yy,ipen)
            ipen = 2
  100   continue
        call newpen(1)
        return 
        end
c
c-----------------------------------------------------------------
c
        subroutine start(dobox,x0,y0,xlen,ylen,
     1      xmn,xmx,ymn,ymx,
     1      titlex, titley,
     2      doxlog,doylog,dopowy)
        parameter (LIN=5,LOT=6,LER=0)

        logical dobox, doxlog, doylog, dopowy
        character titlex*(*), titley*(*)
c-----
c       ensure that the plot limits are in increasing order
c-----
        if(xmn.gt.xmx)then
            xmin = xmx
            xmax = xmn
        else
            xmin = xmn
            xmax = xmx
        endif
        if(ymn.gt.ymx)then
            ymin = ymx
            ymax = ymn
        else
            ymin = ymn
            ymax = ymx
        endif
c-----
c       Plot the frame.
c-----
        lx = lgstr(titlex)
        ly = lgstr(titley)
        if(dobox)then
            call gbox(x0+0.0,y0+0.0,x0+xlen,y0+ylen)
            if(doylog)then
                call dology(x0+0.0 ,y0+0.0,ylen,ymax,ymin,
     1              0.10,.false.,.true.,dopowy,ly,
     2              titley)
                call dology(x0+xlen,y0+0.0,ylen,ymax,ymin,
     1              0.10,.true.,.false.,.false.,ly,
     2              titley)
            else
                call dnliny(x0+0.0 ,y0+ylen,-ylen,-ymin,-ymax,
     1              0.10,.false.,.true.,dopowy,ly,
     2              titley)
                call dnliny(x0+xlen,y0+ylen,-ylen,-ymin,-ymax,
     1              0.10,.true.,.false.,.false.,ly,
     2              titley)
            endif
            if(doxlog)then
                call dologx(x0+0.0,y0+0.0,xlen,xmax,xmin,
     1              0.14,.true.,.false.,.true.,lx,
     2              titlex)
                call dologx(x0+0.0,y0+ylen,xlen,xmax,xmin,
     1              0.14,.false.,.false.,.false.,lx,
     2              titlex)
            else
                call dolinx(x0+0.0 ,y0+0.0,xlen,xmax,xmin,
     1              0.14,.false.,.true.,.true.,lx,
     2              titlex)
                call dolinx(x0+0.0,y0+ylen,xlen,xmax,xmin,
     1              0.14,.true.,.false.,.false.,lx,
     2              titlex)
            endif
        else
c-----
c       put in corner markers
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
        return
        end


        subroutine gcmdln(ilorr, dobox, xmin, xmax, 
     1      ymin, ymax, x0, y0, xlen, ylen, 
     2      kolor,dotext,doclean)
        integer*4 ilorr
        logical dobox, dotext, doclean
        real*4 xmin, xmax
        real*4 x0, y0, xlen, ylen
        integer*4 kolor

        character names*80
c-----
c       defaults
c-----
        ilorr = 0
        dobox = .true.
        xmin = -1.0
        xmax = -1.0
        ymin = -1.0e+38
        ymax = -1.0e+38
        x0 = 2.0
        y0 = 1.0
        xlen = 5.0
        ylen = 5.0
        kolor = 1
        dotext = .false.
        doclean = .false.

        nmarg = mnmarg()
        i = 0
 1000   continue
            i = i + 1
            if(i.gt.nmarg)go to 2000
            call mgtarg(i,names)
            if(names(1:2).eq.'-L')then
                ilorr = 1
            else if(names(1:2).eq.'-R')then
                ilorr = 2
            else if(names(1:5).eq.'-XMIN')then
                i = i + 1
                call mgtarg(i,names)
                read(names,'(bn,f20.0)')xmin
            else if(names(1:5).eq.'-XMAX')then
                i = i + 1
                call mgtarg(i,names)
                read(names,'(bn,f20.0)')xmax
            else if(names(1:5).eq.'-YMIN')then
                i = i + 1
                call mgtarg(i,names)
                read(names,'(bn,f20.0)')ymin
            else if(names(1:5).eq.'-YMAX')then
                i = i + 1
                call mgtarg(i,names)
                read(names,'(bn,f20.0)')ymax
            else if(names(1:3).eq.'-X0')then
                i = i + 1
                call mgtarg(i,names)
                read(names,'(bn,f20.0)')x0
            else if(names(1:3).eq.'-Y0')then
                i = i + 1
                call mgtarg(i,names)
                read(names,'(bn,f20.0)')y0
            else if(names(1:5).eq.'-XLEN')then
                i = i + 1
                call mgtarg(i,names)
                read(names,'(bn,f20.0)')xlen
            else if(names(1:5).eq.'-YLEN')then
                i = i + 1
                call mgtarg(i,names)
                read(names,'(bn,f20.0)')ylen
            else if(names(1:2).eq.'-K')then
                i = i + 1
                call mgtarg(i,names)
                read(names,'(bn,i20)')kolor
            else if(names(1:6).eq.'-NOBOX')then
                dobox = .false.
            else if(names(1:4).eq.'-TXT')then
                dotext = .true.
            else if(names(1:4).eq.'-CLE')then
                doclean = .true.
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
        write(LER,*)'sdpder96: ',str
        write(LER,*)'USAGE:',
     1  'tdpder96 ',
     1  '[-L -R] ',
     1  '-XMIN xmin -XMAX xmax -YMIN ymin -YMAX ymax ',
     1  '-X0 x0 -Y0 y0 -K kolor -NOBOX -TXT -CLEAN -? -h'

        write(LER,*)
     1  '-L                        Love waves'
        write(LER,*)
     1  '-R                        Rayleigh waves'
        write(LER,*)
     1  '     Note one of -L or -R is required'
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
     1  '-K kolor   (default 1  )  color for curves'
        write(LER,*)
     1  '-NOBOX     (default false) do not plot axes'
        write(LER,*)
     1  '-TXT       (default false) output text file'
        write(LER,*)
     1  '-CLEAN     (default false) No period,mode annotation'
        write(LER,*)
     1  '-?         (default false) online help'
        write(LER,*)
     1  '-h         (default false) online help'
        stop 
        end

        subroutine outtxt(ofile, ilorr)
c-----
c       make a listing of all dispersion parameters
c-----
c       ofile   C*(*)   name of output file
c       ilorr   I*4 1 Love
c               2 Rayleigh
c-----
        parameter (LIN=5,LOT=6,LER=0)
        character ofile*(*)

        parameter (NL=1000)
        common/tmodl/d(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      Rho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL),
     3      frefp(NL), frefs(NL)
        real d, ta, tc, tl, tn, tf, rho, qa, qb, etap,etas,frefp, frefs

        
        real*4 dcdh(NL), dcdav(NL), dcdah(NL), dcdbv(NL),
     1     dcdbh(NL), dcdn(NL), dcdr(NL)
        real*4 ur(NL), tur(NL), uz(NL), tuz(NL)

        character mname*80
        integer ipar(20)
        real*4 fpar(20)
c-----
c       FORMAT STATEMENTS
c-----
   14   format( '  LAYER     H(km)       TA          TC          TF',
     1     '          TL          TN        Density',
     2     '     QA(inv)     QB(inv)')

   15   format(i5,9f12.5)
   16   format(' -SURFACE ','- - - - - ','- - - - - ',
     1      '- - - - - ','- - - - - ','- - - - - ','- - -')
   21   format(//17x,'  LOVE WAVE        MODE #',i3/
     1  '        T =',e11.4,' C =   ',e11.4,' U   =',e11.4/
     2  '        AL=',e11.4,' GAMMA=',e11.4,' ZREF=',e11.4/
     3  '    M       UT         TT       ',
     4  'DC/DH      DC/DBv     DC/DBh     DC/DR')
   31    format(//22x,'RAYLEIGH WAVE      MODE #',i3/
     1  '        T =',e11.4,' C =   ',e11.4,' U   =',e11.4/
     2  '        AR=',e11.4,' GAMMA=',e11.4,' ZREF=',e11.4/
     3  '    M       UR         TR        UZ'
     4  ,'         TZ        DC/DH      DC/DAv'
     5  ,'     DC/DAh     DC/DN      DC/DBv     DC/DBh     DC/DR')
c-----
c       open the output file
c-----
            open(2,file=ofile,status='unknown',
     1          access='sequential',form='formatted')
            rewind 2
c-----
c           begin reading the output file of sdisp96 
c           and do integrity check
c-----
            rewind 1
            call getmdt(1,nmax,d,ta,tc,tf,
     1           tl,tn,rho,qa,qb,nper,depths,depthr,
     1           mname,ipar,fpar)
            call gethed(1,ifunc,kmode,t0,ierr)
            if(ierr.eq.200)then
                write(LOT,*)'End of File reading header'
                stop
            endif
            if(ierr.eq.1001)then
                write(LOT,*)'Error reading header'
                stop
            endif
            if(ifunc.eq.5 .and. ilorr.ne.1)then
                write(LOT,*)'Data file is not for Love waves'
                stop
            endif
            if(ifunc.eq.6 .and. ilorr.ne.2)then
                write(LOT,*)'Data file is not for Rayleigh wave'
                stop
            endif
            
            write(2,*) ' '
            write(2,*) 'Model:'
            write(2,14)
            refdep = fpar(1)
            dep = 0.0
            do 100 i=1,nmax
                write(2,15) i,d(i),ta(i),tc(i), 
     1            tf(i), tl(i),tn(i),rho(i),qa(i),qb(i)
                dep = dep + d(i)
                if(dep.eq.refdep)write(2,16)
  100       continue
c-----
c       print the dispersion values.
c       Since the data files are all modes for each frequency, we will
c       repeatedly rwind the input to make a list of all frequencies
c       for each mode
c-----
c-----
c           now start processing from the beginning
c-----
            rewind 1
            call getmdt(1,nmax,d,ta,tc,tf,
     1           tl,tn,rho,qa,qb,nper,depths,depthr,
     1           mname,ipar,fpar)

c-----
c       set up counters
c-----
c------
c     get the eigenfunction information
c     pick up the phase velocities at different frequencies for
c     a particular mode.
c------
        lorr = ilorr + 4
  200   continue
            call gethed(1,ifunc,nmode,t0,ierr)
            omega = 6.2831853/t0
            if(ierr.eq.1001)go to 2001
            if(ierr.eq.200)then
                write(LOT,*) 
     1      'No ifunc=-1. Data end at perd=',t0
                ifunc=-1
            endif
            if(ifunc.le.0) go to 1000
            if(nmode.le.0) go to 200
            do 300 j=1,nmode
            call getdrt(1,lorr,wvno,u,gammal,
     1      sur,sdur,sd2ur,suz,sduz,sd2uz,sare,wvnsrc,sur0,
     2      rur,rtr,ruz,rtz,rare,wvnrec,rur0,
     3      sumkr,sumgr,sumgv,ierr,nmax,dcdh,
     4      dcdav,dcdah, dcdbv, dcdbh, dcdn,dcdr,
     5      ur,tur,uz,tuz,ipar)


                if(ierr.eq.200)go to 2001
                if(ierr.eq.1001)go to 1001
                f0=1./t0
                c = omega/wvno
                if(ilorr.eq.1)then
        write(2,21)j-1,t0,c,u,sare,gammal,refdep
        write(2,1)(i,ur(i),tur(i),dcdh(i),
     1        dcdbv(i),dcdbh(i),dcdr(i),i=1,nmax)
    1   format(i5,1x,6e11.3)
                else if(ilorr.eq.2)then
        write(2,31)j-1,t0,c,u,sare,gammal,refdep
        write(2,2)(i,ur(i),tur(i),uz(i),tuz(i),dcdh(i),
     1     dcdav(i),dcdah(i),dcdn(i),dcdbv(i),dcdbh(i),dcdr(i),i=1,nmax)
    2   format(i5,1x,11e11.3)
                endif
  300       continue
        go to 200
 2001   continue
 1001       continue
        write(LOT,*) ' '
        write(LOT,*)
     *      '(l/r)egn96 data file terminates at  perd=',t0
        write(LOT,*) ' '
 1000   continue
        close(2)
        return
        end

        subroutine limtxy(x,y,n,xmax,xmin,ymax,ymin)
c-----
c       x   R*4 - array of x-values
c       y   R*4 - array of y-values
c       n   I*4 - number of (x,y) pairs
c       xmax    R*4 - maximum value of x
c       xmin    R*4 - minimum value of x
c       ymax    R*4 - maximum value of y
c       ymin    R*4 - minimum value of y
c-----
        real*4 x(n), y(n)
        integer*4 n
        real*4 xmin,xmax, ymin,ymax

        xmin =  1.0e+38
        xmax = -1.0e+38
        ymin =  1.0e+38
        ymax = -1.0e+38
        do 100 i=1,n
            if(abs(x(i)).gt.xmax)xmax = abs(x(i))
            if(x(i).lt.xmin)xmin = x(i)
            if(y(i).gt.ymax)ymax = y(i)
            if(y(i).lt.ymin)ymin = y(i)
  100   continue
c-----
c       now guarantee a default and make it symmetric
c       we foce the limits to the +- so plot is symmetric about
c       x = 0
c-----
        xmax = 1.2*xmax
        xmin = -abs(xmax)
        
        end

        subroutine fillxy(x,y,n,z,ur,nmax,icmd)
c-----
c       fill the plot array with the eigenfunctions
c
c       icmd    I*4 = 1 just plot the points
c               = 2 force constant in a layer
c                   add extra at end
c               = 3 for delta e.g., dc/dh which is partial
c                   for boundary change not layer thickness
c               = 5 for UT if water this is zero at base of layer too
c-----
        real*4 x(n), y(n)
        integer*4 n
        real*4 z(n), ur(n)
        integer*4 nmax
        parameter (NL=1000)
        common/tmodl/d(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      Rho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL),
     3      frefp(NL), frefs(NL)
        real d, ta, tc, tl, tn, tf, rho, qa, qb, etap,etas,frefp, frefs


        logical watsol, solwat
        n = 0
c-----
c       if there is a fluid/solid layer boundary, then add one 
c       additional point at the boundary and generate the 
c       correct ur from the tz
c-----
        watsol = .false.
        solwat = .false.
        do 100 i=1,nmax
            if(i.ne.nmax)then
                if(tl(i).eq.0.0 .and. tl(i+1).ne.0.0)then
                    watsol = .true.
                else
                    watsol = .false.
                endif
            endif
            if(i.ne.1)then
                if(tl(i-1).ne.0.0 .and. tl(i).eq.0.0)then
                    solwat = .true.
                else
                    solwat = .false.
                endif
            endif
            if(icmd.eq.1)then
                n = n + 1
                x(n) = ur(i)
                y(n) = z(i) 
                if(i.ne.nmax)then
                    n = n + 1
                    x(n) = ur(i+1)
                    y(n) = z(i+1) 
                endif
            else if(icmd.eq.2)then
                n = n + 1
                x(n) = ur(i)
                y(n) = z(i) 
                if(i.ne.nmax)then
                    n=n+1
                    x(n) = x(n-1)
                    y(n) = z(i+1)
                endif
            else if(icmd.eq.3)then
                n = n + 1
                x(n) = 0.0
                y(n) = z(i) 
                n = n + 1
                x(n) = ur(i)
                y(n) = z(i)
                n = n + 1
                x(n) = 0.0
                y(n) = z(i) 
            else if(icmd.eq.5)then
                if(solwat)then
                    n = n + 1
                    x(n) = ur(i)
                    y(n) = z(i-1)
                endif
                n = n + 1
                x(n) = ur(i)
                y(n) = z(i) 
                if(watsol)then
                    n = n + 1
                    x(n) = ur(i)
                    y(n) = z(i+1)
                endif
                if(i.ne.nmax)then
                    n = n + 1
                    x(n) = ur(i+1)
                    y(n) = z(i+1) 
                endif
            endif
  100   continue
        return
        end

        subroutine fllxyu(x,y,n,z,ur,tz,nmax,wvno,omega,icmd)
c-----
c       fill the plot array with the eigenfunctions
c
c       icmd    I*4 
c               = 4 special fix for water layer for Ur
c-----
        real*4 x(n), y(n)
        integer*4 n
        real*4 z(n), ur(n), tz(n)
        integer*4 nmax

        parameter (NL=1000)
        common/tmodl/d(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      Rho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL),
     3      frefp(NL), frefs(NL)
        real d, ta, tc, tl, tn, tf, rho, qa, qb, etap,etas,frefp, frefs


        logical watsol, solwat
        n = 0
c-----
c       if there is a fluid/solid layer boundary, then add one 
c       additional point at the boundary and generate the 
c       correct ur from the tz
c-----
        
        do 100 i=1,nmax
c-----
c       watsol .true. This boundary is aa transition from water to solid
c       solwat .true. This boundary is aa transition from solid to water
c-----
            watsol = .false.
            solwat = .false.
            if(i.ne.1)then
                if(tl(i-1).ne.0.0 .and. tl(i).eq.0.0)then
                    solwat = .true.
                else
                    solwat = .false.
                endif
                if(tl(i-1).eq.0.0 .and. tl(i).ne.0.0)then
                    watsol = .true.
                else
                    watsol = .false.
                endif
            endif
            if(icmd.eq.4)then
                if(.not. watsol .and. .not. solwat)then
                    n = n + 1
                    x(n) = ur(i)
                    y(n) = z(i) 
                else if(watsol)then
c-----
c       in a solid layer - the last was a fluid layer
c-----
                    n = n + 1
                    x(n) = -wvno*tz(i)/(omega*omega*rho(i))
                    y(n) = z(i) 
                    n = n + 1
                    x(n) = ur(i)
                    y(n) = z(i) 
                else if(solwat)then
c-----
c       in a water layer - the last was a solid layer
c-----
                    n = n + 1
                    x(n) = ur(i)
                    y(n) = z(i) 
                    n = n + 1
                    x(n) = -wvno*tz(i)/(omega*omega*rho(i))
                    y(n) = z(i) 
                endif
            endif
  100   continue
        return
        end

        subroutine doit1(z,u,nmax,jcmd,strx,stry,
     1      xlen,ylen,doxlog,doylog,dopowy,
     2      x0,y0,c,t0,j,dobox,kolor,doclean)
c-----
c       interface to plotting routines
c-----
        real*4 u(*), z(*)
        character strx*(*), stry*(*)
        logical doxlog,doylog,dopowy,dobox,doclean
        parameter (NL=1000)
        real*4 x(NL+NL+NL), y(NL+NL+NL)
            call fillxy(x,y,n,z,u,nmax,jcmd)
            call limtxy(x,y,n,xmax,xmin,ymax,ymin)
            call plotit(x,y,n,x0,y0,strx,stry,
     1          c,t0,j,kolor,dobox,
     1          xmin,xmax,ymin,ymax,
     2          xlen,ylen,doxlog,doylog,dopowy,doclean)
            call frame()
        return
        end

        subroutine doit2(z,u,t,nmax,jcmd,strx,stry,
     1      xlen,ylen,doxlog,doylog,dopowy,
     2      x0,y0,c,t0,j,dobox,kolor,
     3      omega,wvno,doclean)
c-----
c       interface to plotting routines
c-----
        real*4 u(*), z(*), t(*)
        character strx*(*), stry*(*)
        logical doxlog,doylog,dopowy,dobox,doclean
        parameter (NL=1000)
        real*4 x(NL+NL+NL), y(NL+NL+NL)
            call fllxyu(x,y,n,z,u,t,nmax,wvno,omega,jcmd)
            call limtxy(x,y,n,xmax,xmin,ymax,ymin)
            call plotit(x,y,n,x0,y0,strx,stry,
     1          c,t0,j,kolor,dobox,
     1          xmin,xmax,ymin,ymax,
     2          xlen,ylen,doxlog,doylog,dopowy,doclean)
            call frame()
        return
        end
