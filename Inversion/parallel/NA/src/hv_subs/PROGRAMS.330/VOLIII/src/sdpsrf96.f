c----------------------------------------------------------------------c
c                                                                      c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                 c
c      VOLUME III                                                      c
c                                                                      c
c      PROGRAM: SDPSRF96                                               c
c                                                                      c
c      COPYRIGHT 1996                                                  c
c      R. B. Herrmann, Chien-Ying Wang                                 c
c      Department of Earth and Atmospheric Sciences                    c
c      Saint Louis University                                          c
c      221 North Grand Boulevard                                       c
c      St. Louis, Missouri 63103                                       c
c      U. S. A.                                                        c
c                                                                      c
c----------------------------------------------------------------------c
c Revision history:
c       27 DEC 2000 - removed code commented out
c----------------------------------------------------------------------c
        parameter (LIN=5,LOT=6,LER=6)

        parameter (NL=1000)
        common/modl/d(NL), a(NL), b(NL), rho(NL), qa(NL), qb(NL)

        common/dim/ cp
        parameter(MAXMOD=2000)
        real*8 cp(MAXMOD)
        common/ctrl/ tmin,tmax,cmin,cmax

        character ofile*20, dfile*20, afile*20
        logical ext
c-----
c       command line arguments
c-----
        logical dofreq, dobox, dotext
        logical dotmp, doxlog
        logical doasc

c
c-----
c       machine dependent initialization
c-----
        call mchdep()
c-----
c       get control parameters from command line
c-----
        call gcmdln(ilorr, dofreq, dobox, xmin, xmax, ymin, ymax, 
     1      x0, y0, xlen, ylen, 
     2      kolor, dotext,dotmp,doxlog,doasc)
        if(ilorr.lt.1)call usage('specify -L or -R')
c
c------
c     get control parameters and open input file.
c------
        pi=2.*3.141592653
        
        if(ilorr.eq.1)then
            if(dotmp)then
                dfile = 'tsdisp96.lov'
            else
                dfile = 'sdisp96.lov'
            endif
            ofile = 'SDISPL.TXT'
            afile = 'SDISPL.ASC'
        else if(ilorr.eq.2)then
            if(dotmp)then
                dfile = 'tsdisp96.ray'
            else
                dfile = 'sdisp96.ray'
            endif
            ofile = 'SDISPR.TXT'
            afile = 'SDISPR.ASC'
        endif
        inquire(file=dfile,exist=ext)
        if(.not.ext)then
            call usage('Dispersion file '//dfile//
     1      ' does not exist')
        endif
        open(1,file=dfile,status='unknown',
     1      access='sequential',form='unformatted')

        if(dotext)then
            rewind 1
            call outtxt(ofile,dofreq,cp)
        endif
        if(doasc )then
            rewind 1
            call outasc(afile,dofreq,cp)
        endif

        call outplt(ilorr, dofreq, dobox, xmin, xmax, ymin, ymax, 
     1      x0, y0, xlen, ylen, 
     2      kolor,doxlog)

        close (1)
        end

        subroutine  outplt(ilorr, dofreq, dobox, xmin, xmax, ymin, ymax,
     1      x0, y0, xlen, ylen, 
     2      kolor,doxlog)
c-----
c       plot the dispersion curves
c
c       ilorr   I*4 1 = Love, 2 = Rayleigh
c       dofreq  L   .true. horizontal axis is frequency
c       dobox   L   .true. put frame around plot
c       xmin    R*4 minimum value of independent variable
c       xmax    R*4 maximum value of independent variable
c       ymin    R*4 minimum value of   dependent variable
c       ymax    R*4 maximum value of   dependent variable
c       x0  R*4 
c       y0  R*4 lower left coordinates of plot box
c       xlen    R*4 length of X axis in inches
c       ylen    R*4 length of Y axis in inches
c-----
        integer ilorr, iporg
        logical dofreq, dobox, doxlog
        real xmin, xmax, ymin, ymax
        real x0, y0, xlen, ylen
        integer kolor

        parameter (LIN=5,LOT=6,LER=6)

        parameter (NL=1000)
        common/modl/d(NL), a(NL), b(NL), rho(NL), qa(NL), qb(NL)

        integer*4 nmax,nper,ifunc,kmode

        parameter(MAXMOD=2000)
        common/dim/ cp
        real*8 cp(MAXMOD), t0


        common/ctrl/ tmin,tmax,cmin,cmax


        character mname*80
        integer ipar(10)
        real*4 fpar(10)

        logical doylog
        doylog = .false.
c------
c       plot the dispersion values.
c------
        if(dofreq)then
            ipg = 2
        else
            ipg = 1
        endif
c------
c       plot the frame and find the mode range to be plotted.
c------
        call limt(dofreq,modmax)
        if(xmin .lt.0 .or. xmax .lt.0)then
            xmin = tmin
            xmax = tmax
        endif
        if(ymin .lt.0 .or. ymax .lt.0)then
            ymin = cmin
            ymax = cmax
        endif
c-----
c       To ensure useful logatirhmic plots, ensure that there
c       is at least one complete cycle
c-----
        if(doylog)then
            if( ymin .gt. 0.1*ymax)ymin = 0.099*ymax
        endif
        if(doxlog)then
            if( xmin .gt. 0.1*xmax)xmin = 0.099*xmax
        endif

        iporg = 1
        call start(iporg,ilorr,x0,y0,xlen,ylen,
     1      xmin,xmax,ymin,ymax,
     2      dobox,doxlog,doylog, dofreq)
c-----
c       plot the dispersion curves
c-----
        call newpen(kolor)
        icnt=0
  400   continue
            icnt=icnt+1
            if(icnt.gt.modmax) go to 900
c-----
c       go to beginning of data file
c-----
            rewind 1
            call gtsmdl(1,nmax,d,a,b,rho,qa,qb,nper,
     1      mname,ipar,fpar)
c------
c     Pick up the desired dispersion curve from data file.
c       and plot in groups  of 1024 pairs
c------
c       jplt    output of rclip, plot line segment if jplt > 0
c       (xlow, ylow) and (xup,yup) are the original
c           unclipped segment ends. If clipping is
c           successful, clipped segment is defined by 
c           (xc1,yc1)  (xc2,yc2)
c       iscont > 0 indicates that this segment is connected 
c           to the previous
c-----
        jplt = 0
        xlow = -1
        xup = -1
        ylow = -1
        yup = -1
        iscont = -1
  500   continue
            call gtshed(1,ifunc,kmode,t0,ierr)
            if(ierr.eq.200)then
            write(LOT,*) 'No ifunc=-1. Data end at perd=',t0
                ifunc=-1
            else if(ierr.eq.1001)then
                go to 1001
            endif
            if(ifunc.le.0) then
                go to 400
            endif
            if(kmode.le.0) go to 500
            call gtsval(1,cp,kmode,ierr)    
            if(ierr.eq.200)then
                go to 1001
            else if(ierr.eq.1001)then
                go to 1001
            endif
c-----
c       get the proper transformation
c-----
            if(ipg.eq.2)then
                t = 1.0/t0
            else
                t = t0
            endif
            if(icnt.le.kmode)then
            v = cp(icnt)
            xup  = t
            yup  = v
            
            if(jplt.gt.0)then
                call rclip(xmin,xmax,ymin,ymax,
     2              xc1,yc1,xc2,yc2,
     3              xlow,ylow,xup,yup,iplt)
                if(iplt.gt.0)then
                if(doxlog)then
                     xx = x0 + xlen*alog10(xc1/xmin)
     1                  /alog10(xmax/xmin)
                else
                     xx = x0 + xlen*(xc1 - xmin)/(xmax - xmin)
                endif
                if(doylog)then
                     yy = y0 + ylen*alog10(yc1/ymin)
     1                  /alog10(ymax/ymin)
                else
                     yy = y0 + ylen*(yc1 - ymin)/(ymax - ymin)
                endif
                if(iscont.lt.0)call plot(xx,yy,3)
                if(doxlog)then
                     xx = x0 + xlen*alog10(xc2/xmin)
     1                  /alog10(xmax/xmin)
                else
                     xx = x0 + xlen*(xc2  - xmin)/(xmax - xmin)
                endif
                if(doylog)then
                     yy = y0 + ylen*alog10(yc2/ymin)
     1                  /alog10(ymax/ymin)
                else
                     yy = y0 + ylen*(yc2  - ymin)/(ymax - ymin)
                endif
                call plot(xx,yy,2)
                iscont = 1
                else
                iscont = -1
                
                endif
            endif
            if(icnt.gt.kmode)then
                jplt = 0
            else
                jplt = jplt + 1
            endif
            xlow = xup
            ylow = yup
            endif
        go to 500
  900   continue
        write(LOT,*) 'Last mode plotted=',icnt-1
        go to 1100
 1001   continue
        write(LOT,*) ' '
        write(LOT,*)
     *  'sdisp96 data file terminates at  perd=',t0
        write(LOT,*) ' '
 1100   continue
        call newpen(1)
        call pend()
        end

        subroutine start(iporg,ilorr,x0,y0,xlen,ylen,
     1      xmn,xmx,ymn,ymx,
     2      dobox,doxlog,doylog, dofreq)
        parameter (LIN=5,LOT=6,LER=6)

        logical dobox, doxlog, doylog, dofreq
        character ylabel*12, xlabel*14
c------
c       open the plotter and read in the plot control parameter.
c------
        if(dofreq)then
            xlabel = 'Frequency (Hz)'
        else
            xlabel = 'Period (sec)'
        endif
        if(ilorr.eq.1)then
            if(iporg.eq.1)then
                call pinitf('SDISPL.PLT')
                ylabel = 'C (km/s)'
            endif
        else if(ilorr.eq.2)then
            if(iporg.eq.1)then
                call pinitf('SDISPR.PLT')
                ylabel = 'C (km/s)'
            endif
        endif
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
        lx = lgstr(xlabel)
        ly = lgstr(ylabel)
        if(dobox)then
            call gbox(x0+0.0,y0+0.0,x0+xlen,y0+ylen)
            if(doylog)then
                call dology(x0+0.0 ,y0+0.0,ylen,ymax,ymin,
     1              0.14,.false.,.true.,.true.,ly,
     2              ylabel)
                call dology(x0+xlen,y0+0.0,ylen,ymax,ymin,
     1              0.14,.true.,.false.,.false.,ly,
     2              ylabel)
            else
                call doliny(x0+0.0 ,y0+0.0,ylen,ymax,ymin,
     1              0.14,.false.,.true.,.true.,ly,
     2              ylabel)
                call doliny(x0+xlen,y0+0.0,ylen,ymax,ymin,
     1              0.14,.true.,.false.,.false.,ly,
     2              ylabel)
            endif
            if(doxlog)then
                call dologx(x0+0.0,y0+0.0,xlen,xmax,xmin,
     1              0.14,.true.,.false.,.true.,lx,
     2              xlabel)
                call dologx(x0+0.0,y0+ylen,xlen,xmax,xmin,
     1              0.14,.false.,.false.,.false.,lx,
     2              xlabel)
            else
                call dolinx(x0+0.0 ,y0+0.0,xlen,xmax,xmin,
     1              0.14,.true.,.false.,.true.,lx,
     2              xlabel)
                call dolinx(x0+0.0,y0+ylen,xlen,xmax,xmin,
     1              0.14,.false.,.true.,.false.,lx,
     2              xlabel)
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

        subroutine limt(dofreq,modmax)
c-----
c       find the range of mode numbers to be plotted.
c       and also plotting limits
c-----
        parameter (LIN=5,LOT=6,LER=6)
        logical dofreq

        parameter (NL=1000)
        common/modl/d(NL), a(NL), b(NL), rho(NL), qa(NL), qb(NL)

        integer*4 ifunc,kmode

        parameter(MAXMOD=2000)
        common/dim/  cp
        real*8 cp(MAXMOD)
        real*8 t0
        real*4 tp,vp

        common/ctrl/ tmin,tmax,cmin,cmax

        character mname*80
        integer ipar(10)
        real*4 fpar(10)
c-----
c       go to beginning of data file
c-----
        rewind 1
        call gtsmdl(1,nmax,d,a,b,rho,qa,qb,nper,
     1      mname,ipar,fpar)
        kfirst = 0
c------
c       find tmin and tmax in which data will be plotted.
c------
        tmin=1.0e+10
        tmax=0.0
        cmin = 1.0e+10
        cmax = 0.0
        modmax = 0
  100   continue
            call gtshed(1,ifunc,kmode,t0,ierr)
            if(ierr.eq.200)then
                write(LOT,*) 'No ifunc=-1. Data end at perd=',t0
                ifunc=-1
            else if(ierr.eq.1001)then
                go to 1001
            endif
            if(ifunc.le.0) go to 400
            if(kmode.le.0) go to 100
            call gtsval(1,cp,kmode,ierr)    
            if(ierr.eq.200)then
                go to 1000
            else if(ierr.eq.1001)then
                go to 1001
            endif
            if(kmode .gt. modmax) modmax=kmode
            if(dofreq)then
                tp = 1.0/t0
            else
                tp = t0
            endif
            if(tp.gt.tmax)tmax = tp
            if(tp.lt.tmin)tmin = tp
            do 300 i=1,kmode
                vp=cp(i)
                if(vp.gt.cmax)cmax = vp
                if(vp.lt.cmin)cmin = vp
  300       continue
        go to 100
  400   continue
        return
 1000   continue
        if(kfirst.eq.0) then
            kfirst=1
            write(LOT,*) 'Data end, not normal',
     1          ' termination at perd=',t0
        endif
        return
 1001   continue
        write(LOT,*) ' '
        write(LOT,*)
     *  'surface85 data file terminates at  perd=',t0
        write(LOT,*) ' '
        stop
        end

        subroutine gcmdln(ilorr, dofreq, dobox, xmin, xmax, cmin, cmax, 
     1      x0, y0, xlen, ylen, 
     2      kolor, dotext,dotmp,doxlog,doasc)
c       dotmp   L   - .true. use file tsdisp96.lov
        integer*4 ilorr
        logical dofreq, dobox, dotext,dotmp, doxlog
        logical doasc
        real*4 xmin, xmax, cmin, cmax
        real*4 x0, y0, xlen, ylen
        integer*4 kolor

        character names*80
c-----
c       defaults
c-----
        ilorr = 0
        dofreq = .true.
        dobox = .true.
        xmin = -1.0
        xmax = -1.0
        cmin = -1.0
        cmax = -1.0
        x0 = 2.0
        y0 = 1.0
        xlen = 6.0
        ylen = 6.0
        kolor = 1
        dotext = .false.
        dotmp  = .false.
        doxlog = .false.
        doasc  = .false.

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
            else if(names(1:5).eq.'-FREQ')then
                dofreq = .true. 
            else if(names(1:4).eq.'-PER')then
                dofreq = .false.    
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
                read(names,'(bn,f20.0)')cmin
            else if(names(1:5).eq.'-YMAX')then
                i = i + 1
                call mgtarg(i,names)
                read(names,'(bn,f20.0)')cmax
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
            else if(names(1:4).eq.'-ASC')then
                doasc  = .true.
            else if(names(1:2) .eq. '-T'.and.
     1          names(1:4).ne.'-TXT')then
                    dotmp = .true.
            else if(names(1:5).eq.'-XLOG')then
                doxlog = .true.
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
        parameter (LER=6,LIN=5,LOT=6)
        write(LER,*)'sdpsrf96 : ',str
        write(LER,*)'USAGE:',
     1  'sdpsrf96 ',
     1  '[-L -R]  [-FREQ -PER]',
     1  '-XMIN xmin -XMAX xmax -YMIN ymin -YMAX ymax ',
     1  '-X0 x0 -Y0 y0 -K kolor -NOBOX -TXT -XLOG -ASC -? -h'

        write(LER,*)
     1  '-L                        Love waves'
        write(LER,*)
     1  '-R                        Rayleigh waves'
        write(LER,*)
     1  '     Note one of -L or -R is required'
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
     1  '-K kolor   (default 1  )  color for curves'
        write(LER,*)
     1  '-NOBOX     (default false) do not plot axes'
        write(LER,*)
     1  '-TXT       (default false) output text file'
        write(LER,*)
     1  '-ASC       (default false) output ASCII column file'
        write(LER,*)
     1  '-T          (default false) use tdisp96.lov or tdisp96.ray'
        write(LER,*)
     1  '-XLOG      (default linear) X axis is logarithmic'
        write(LER,*)
     1  '-?         (default false) online help'
        write(LER,*)
     1  '-h         (default false) online help'
        stop 
        end

        
        subroutine outtxt(ofile, dofreq, cp)
c-----
c       make a listing of all dispersion parameters
c-----
c       ofile   C*(*)   name of output file
c       dofreq  L   .true. output frequency
c                   else period
c       cp  R*8 array of phase velocities 
c-----
        parameter (LIN=5,LOT=6,LER=6)
        character ofile*(*)
        logical dofreq

        parameter (NL=1000)
        common/modl/d(NL), a(NL), b(NL), rho(NL), qa(NL), qb(NL)

        parameter(MAXMOD=2000)
        real*8 cp(MAXMOD), t0

        character mname*80
        integer ipar(10)
        real*4 fpar(10)
c-----
c       FORMAT STATEMENTS
c-----

   10   format(///19x,'LOVE WAVE      MODE #',i3//21x,
     1   '    PERIOD         PHASE VEL')
   11   format(///19x,'LOVE WAVE      MODE #',i3//21x,
     1   '      FREQ         PHASE VEL')
   15   format(4f12.5)
   20   format(1x,12x,i5,1x,f15.10,2x,f15.10)
   30   format(///18x,'RAYLEIGH WAVE      MODE #',i3//21x,
     1   '    PERIOD          PHASE VEL')
   31   format(///18x,'RAYLEIGH WAVE      MODE #',i3//21x,
     1   '      FREQ          PHASE VEL')
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
            call gtsmdl(1,nmax,d,a,b,rho,qa,qb,nper,
     1      mname,ipar,fpar)
            call gtshed(1,ifunc,kmode,t0,ierr)
            if(ierr.eq.200)then
                write(LOT,*)'End of File reading header'
                go to 1000
            endif
            if(ierr.eq.1001)then
                write(LOT,*)'Error reading header'
                go to 1000
            endif
            
            write(LOT,*) ' '
            write(LOT,*) 'Total period=',nper
            write(LOT,*) ' '
            write(LOT,*)
     *      'Number of modes at the lowest period:',t0,' = ',kmode
            write(LOT,*) ' '
            write(LOT,*) 'Model:'
            write(LOT,15) (d(i),a(i),b(i),rho(i),i=1,nmax)
c-----
c       print the dispersion values.
c       Since the data files are all modes for each frequency, we will
c       repeatedly rewind the input to make a list of all frequencies
c       for each mode
c-----
            icnt = 0
            ilorr = ifunc
  100       continue 
            icnt = icnt + 1
            if(icnt.gt.kmode)go to 1000
c-----
c           now start processing from the beginning
c-----
                rewind 1
                call gtsmdl(1,nmax,d,a,b,rho,qa,qb,nper,
     1      mname,ipar,fpar)
                if(ilorr.eq.1) then
                    if(dofreq)then
                        write(2,11) icnt-1
                    else
                        write(2,10) icnt-1
                    endif
                else
                    if(dofreq)then
                        write(2,31) icnt-1
                    else
                        write(2,30) icnt-1
                    endif
                endif
c-----
c       set up counters
c-----
            nn=0
c------
c     pick up the phase velocities at different frequencies for
c     a particular mode.
c------
  200       continue
                call gtshed(1,ifunc,nmode,t0,ierr)
                if(ierr.eq.1001)go to 2001
                if(ierr.eq.200)then
                    write(LOT,*) 
     1          'No ifunc=-1. Data end at perd=',t0
                    ifunc=-1
                endif
                if(ifunc.le.0) go to 100
                if(nmode.le.0) go to 200
                call gtsval(1,cp,nmode,ierr)
                if(ierr.eq.200)go to 2001
                if(ierr.eq.2001)go to 1001
                nn = nn + 1
                t00=t0
                if(dofreq) then
                    t00=1./t0
                else
                    t00=t0
                endif
                if(icnt.le.nmode)then
                    write(2,20) nn,t00,cp(icnt)
                endif
            go to 200
 2001   continue
 1001       continue
        write(LOT,*) ' '
        write(LOT,*)
     *      'sdisp96 data file terminates at  perd=',t0
        write(LOT,*) ' '
 1000   continue
        close(2)
        return
        end

        
        subroutine outasc(afile, dofreq, cp)
c-----
c       make a listing of all dispersion parameters
c-----
c       afile   C*(*)   name of output file
c       dofreq  L   .true. output frequency
c                   else period
c       cp  R*8 array of phase velocities 
c-----
        parameter (LIN=5,LOT=6,LER=6)
        character afile*(*)
        logical dofreq

        parameter (NL=1000)
        common/modl/d(NL), a(NL), b(NL), rho(NL), qa(NL), qb(NL)

        parameter(MAXMOD=2000)
        real*8 cp(MAXMOD), t0

        character mname*80
        integer ipar(10)
        real*4 fpar(10)
c-----
c       FORMAT STATEMENTS
c-----

c-----
c       open the output file
c-----
            open(2,file=afile,status='unknown',
     1          access='sequential',form='formatted')
            rewind 2
c-----
c           begin reading the output file of sdisp96 
c           and do integrity check
c-----
            call gtsmdl(1,nmax,d,a,b,rho,qa,qb,nper,
     1      mname,ipar,fpar)
            call gtshed(1,ifunc,kmode,t0,ierr)
            if(ierr.eq.200)then
                write(LOT,*)'End of File reading header'
                go to 1000
            endif
            if(ierr.eq.1001)then
                write(LOT,*)'Error reading header'
                go to 1000
            endif
c-----
c       output column header
c-----
            if(ifunc.eq.1 )then
            write(2,'(a,a)')
     1  'LMODE NFREQ       PERIOD(S)     FREQUENCY(Hz)     ',
     2  '      C(KM/S)'
            else
            write(2,'(a,a)')
     1  'RMODE NFREQ       PERIOD(S)     FREQUENCY(Hz)     ',
     2  '      C(KM/S)'
            endif
            
c-----
c       print the dispersion values.
c       Since the data files are all modes for each frequency, we will
c       repeatedly rwind the input to make a list of all frequencies
c       for each mode
c-----
            icnt = 0
            ilorr = ifunc
  100       continue 
            icnt = icnt + 1
            if(icnt.gt.kmode)go to 1000
c-----
c           now start processing from the beginning
c-----
                rewind 1
                call gtsmdl(1,nmax,d,a,b,rho,qa,qb,nper,
     1      mname,ipar,fpar)
c-----
c       set up counters
c-----
            nn=0
c------
c     pick up the phase velocities at different frequencies for
c     a particular mode.
c------
  200       continue
                call gtshed(1,ifunc,nmode,t0,ierr)
                if(ierr.eq.1001)go to 2001
                if(ierr.eq.200)then
                    write(LOT,*) 
     1          'No ifunc=-1. Data end at perd=',t0
                    ifunc=-1
                endif
                if(ifunc.le.0) go to 100
                if(nmode.le.0) go to 200
                call gtsval(1,cp,nmode,ierr)
                if(ierr.eq.200)go to 2001
                if(ierr.eq.2001)go to 1001
                nn = nn + 1
c-----
c       t0 is period
c       f0 is frequency
c-----
                f0=1./t0
                icntm1 = icnt-1
   10   format(i5,1x,i5,1x,g18.10,1x,g18.10,1x,g18.10) 
                if(icnt.le.nmode)then
                    write(2,10) 
     1                  icntm1,nn,t0,f0,cp(icnt)
                endif
            go to 200
 2001   continue
 1001       continue
        write(LOT,*) ' '
        write(LOT,*)
     *      'sdisp96 data file terminates at  perd=',t0
        write(LOT,*) ' '
 1000   continue
        close(2)
        return
        end
