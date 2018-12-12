c----------------------------------------------------------------------c
c                                                                      c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                 c
c      VOLUME III                                                      c
c                                                                      c
c      PROGRAM: SDPEGN96                                               c
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
c       13 NOV 2001 - ensured that all modes are out for -S option
c       17 NOV 2003 - add option to plot ARE ALE (-A0) or
c               ARE sqrt(c) ALE sqrt(c) for -Ac for SPAC
c       20 JUL 2004 - changed order of declaration of NOBS in outplt
c       28 MAY 2007 - -TXT -M together give debug output
c       13 MAY 2010 - added -NOBLACK to prevent black outline around
c                    observed dispersion for compatibility with sdpdsp96
c       11 MAY 2011 - added -TICONLY to permit tics but no annotation
c                   - added -ONLYKM to label axis velocty but no C or U
c----------------------------------------------------------------------c
        integer LIN, LER, LOT
        parameter (LIN=5,LOT=6,LER=0)

        parameter (NL=1000)
        common/modl/d(NL), a(NL), b(NL), rho(NL), qa(NL), qb(NL)

        common/ctrl/ tmin,tmax,cmin,cmax

        character ofile*20, dfile*20, afile*20, dsfile*20
        logical ext
c-----
c       command line arguments
c-----
        logical dofreq, dobox, dotext, doxlog, doylog, dolat, dosurf
        logical doasc, doerr, doshwmd, doblack, doticonly, dokmonly
        logical verbose
        integer lintyp
        integer*4 ilorr, iporg
        real*4 xmin, xmax, ymin, ymax, x0, y0, xlen, ylen, width
        integer*4 kolor
        character datfil*80

c
c-----
c       machine dependent initialization
c-----
        call mchdep()
c-----
c       get control parameters from command line
c-----
        call gcmdln(ilorr, dofreq, dobox, xmin, xmax, ymin, ymax, 
     1      x0, y0, xlen, ylen, width,
     2      kolor, dotext,iporg,doxlog,dolat,doticonly,dokmonly,
     3      doasc,dosurf,datfil,doerr,doshwmd,lintyp,doblack,
     4      verbose)
        if(ilorr.lt.1)call usage('specify -L or -R')
        if(iporg.lt.1)call usage('specify -C or -U or -G -A0 -Ac')
        if(iporg.eq.3 .or.iporg.eq.4 .or. iporg.eq.5)then
            doylog = .true.
        else
            doylog = .false.
        endif
        if(verbose)then
            WRITE(LOT,*)'ilorr     :', ilorr
            WRITE(LOT,*)'dofreq    :', dofreq
            WRITE(LOT,*)'dobox     :', dobox 
            WRITE(LOT,*)'xmin      :', xmin  
            WRITE(LOT,*)'xmax      :', xmax  
            WRITE(LOT,*)'ymin      :', ymin  
            WRITE(LOT,*)'ymax      :', ymax  
            WRITE(LOT,*)'x0        :', x0    
            WRITE(LOT,*)'y0        :', y0    
            WRITE(LOT,*)'xlen      :', xlen    
            WRITE(LOT,*)'ylen      :', ylen    
            WRITE(LOT,*)'width     :', width    
            WRITE(LOT,*)'kolor     :', kolor    
            WRITE(LOT,*)'iporg     :', iporg     
            WRITE(LOT,*)'doxlog    :', doxlog    
            WRITE(LOT,*)'dolat     :', dolat     
            WRITE(LOT,*)'doticonly :', doticonly     
            WRITE(LOT,*)'dokmonly  :', dokmonly     
            WRITE(LOT,*)'doasc     :', doasc        
            WRITE(LOT,*)'dosurf    :', dosurf        
            WRITE(LOT,*)'doswhmd   :', doshwmd        
            WRITE(LOT,*)'lintyp    :', lintyp        
            WRITE(LOT,*)'doblack   :', doblack        
            WRITE(LOT,*)'datfil    :', datfil         
        endif
c------
c     get control parameters and open input file.
c------
        pi=2.*3.141592653
        if(ilorr.eq.1)then
            if(dolat)then
                dfile = 'slatl96.egn'
            else
                dfile = 'slegn96.egn'
            endif
            ofile = 'SLEGN.TXT'
            afile = 'SLEGN.ASC'
            dsfile = 'SLEGN.dsp'
        else if(ilorr.eq.2)then
            if(dolat)then
                dfile = 'slatr96.egn'
            else
                dfile = 'sregn96.egn'
            endif
            ofile = 'SREGN.TXT'
            afile = 'SREGN.ASC'
            dsfile = 'SREGN.dsp'
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
            call outtxt(ofile,dofreq,ilorr,doshwmd)
        endif

        if(doasc )then
            call outasc(afile,dofreq,ilorr)
        endif

        if(dosurf)then
            call outdsp(dsfile,dofreq,ilorr)
        endif

        call outplt(ilorr, dofreq, dobox, xmin, xmax, ymin, ymax, 
     1      x0, y0, xlen, ylen, width, doblack,
     2      kolor,iporg,doxlog,doylog,datfil,doerr,lintyp,
     3      dokmonly, doticonly)

        close (1)
        end

        subroutine  outplt(ilorr, dofreq, dobox, xmin, xmax, ymin, ymax,
     1      x0, y0, xlen, ylen, width, doblack,
     2      kolor,iporg,doxlog,doylog,datfil,doerr,lintyp,
     3      dokmonly, doticonly)
        integer LIN, LER, LOT
        parameter (LIN=5,LOT=6,LER=0)
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
c       kolor   I*4 CALPLOT pen value
c       iporg   I*4 1 phase velocity
c               2 group velocity
c               3 gamma
c               4 A0
c               5 A0 sqrt(c)
c       datfil  Ch*(*)  User provided dispersion points
c       lintyp  I*4 0 - solid
c                   1 - short dash
c                   2 - long dash
c-----
        integer*4 ilorr
        logical dofreq, dobox, doxlog, doylog, dokmonly, doticonly
        real*4 xmin, xmax, ymin, ymax
        real*4 x0, y0, xlen, ylen
        integer*4 kolor, iporg
        integer lintyp
        character datfil*(*)
        logical doerr
        logical doblack

        integer NL
        parameter (NL=1000)
        common/modl/d(NL), a(NL), b(NL), rho(NL), qa(NL), qb(NL)

        integer*4 nmax,nper,ifunc,kmode
c-----
c       user provided data
c-----
        integer ierr, ndsp, NOBS
        parameter (NOBS=1500000)
        integer*4 jlorr(NOBS), jobs(NOBS),jobsyn(NOBS), jmode(NOBS)  
        real*4  fper(NOBS)   , fobs(NOBS)   , fobserr(NOBS)



        common/ctrl/ tmin,tmax,cmin,cmax

        character mname*80
        integer ipar(10)
        real*4 fpar(10)

        logical inside
        integer ipat
        real patlen

        pi=2.*3.141592653
c------
c       plot the dispersion values.
c------
c------
c       plot the frame and find the mode range to be plotted.
c------
        call limt(dofreq,modmax,iporg)
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

        call start(dobox,iporg,ilorr,x0,y0,xlen,ylen,
     1      xmin,xmax,ymin,ymax,
     1      doxlog,doylog,dofreq,dokmonly,doticonly)
c-----
c       plot the dispersion curves
c-----
        call newpen(kolor)
        if(width.gt.0.0)call gwidth(width)
        icnt=0
  400   continue
            icnt=icnt+1
            if(icnt.gt.modmax) go to 900
c-----
c       go to beginning of data file
c-----
            rewind 1
            call getmdl(1,nmax,d,a,b,rho,qa,qb,
     1          nper,depths,depthr,
     1          mname,ipar,fpar)
c------
c       Pick up the desired dispersion curve from data file.
c------
c       jplt    output of rclip, plot line segment if jplt > 0
c       (xlow, ylow) and (xup,yup) are the original
c           unclipped segment ends. If clipping is
c           successful, clipped segment is defined by 
c             (xc1,yc1)  (xc2,yc2)
c       iscont > 0 indicates that this segment is connected 
c             to the previous
c-----
        jplt = 0
        xlow = -1
        xup = -1
        ylow = -1
        yup = -1
        iscont = -1
  500   continue
            call gethed(1,ifunc,kmode,t0,ierr)
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
            do 300 i=1,kmode
        call getegn(1,ifunc,1,wvno,u,gamma,
     1      sur,sdur,suz,sduz,sare,wvnsrc,sur0,
     2      rur,rdur,ruz,rduz,rare,wvnrec,rur0,
     3      sumkr,sumgr,sumgv,ierr)
            if(ierr.eq.200)then
                go to 1001
            else if(ierr.eq.1001)then
                go to 1001
            endif
            if(u.gt.0.0)then
            if(i.eq.icnt)then
c-----
c       get the proper transformation
c-----
                if(dofreq)then
                    t = 1.0/t0
                else
                    t = t0
                endif
                if(iporg.eq.1)then
                    v = (pi/t0)/wvno
                else if(iporg.eq.2)then 
                    v = u
                else if(iporg.eq.3)then 
                    v = gamma
                else if(iporg.eq.4)then 
                    v = rare
                else if(iporg.eq.5)then 
                    v = rare*sqrt((pi/t0)/wvno)
                endif
                xup  = t
                yup  = v
            
                if(jplt.gt.0)then
                    call rclip(xmin,xmax,ymin,ymax,
     2                  xc1,yc1,xc2,yc2,
     3                  xlow,ylow,xup,yup,iplt)

                if(iplt.gt.0)then
                    call dotran(doxlog,xmin,xmax,xc1,
     1                  x0,xlen,xx,
     1                  doylog,ymin,ymax,yc1,0.0,
     1                  y0,ylen,yy,yvp,yvm)
                    if(iscont.lt.0)call plot(xx,yy,3)
                    call dotran(doxlog,xmin,xmax,xc2,
     1                  x0,xlen,xx,
     1                  doylog,ymin,ymax,yc2,0.0,
     1                  y0,ylen,yy,yvp,yvm)
                    if(lintyp.eq.1)then
                          ipat = 3
                          patlen = 0.05
                          call plotd(xx,yy,ipat,patlen)
                    else if(lintyp.eq.2)then
                          ipat = 7
                          patlen = 0.05
                          call plotd(xx,yy,ipat,patlen)
                    else 
                          call plot(xx,yy,2)
                    endif
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
            endif
  300       continue
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
        call gwidth(0.001)
        call newpen(1)
c-----
c       plot observation points
c-----
        if(datfil .ne. ' ')then
c-----
c       get dispersion values
c-----
            call rddisp(datfil,3,ierr,ndsp,NOBS,
     1          jlorr,jobs,jobsyn,jmode,
     2          fper,fobs,fobserr)
        write(LOT,*)'NDSP:',ndsp
c-----
c       plot the observed dispersion values
c-----
            do 2000 i=1,ndsp
                if(dofreq)then
                    t = 1.0/fper(i)
                else
                    t = fper(i)
                endif
                if(iporg .eq. jobs(i) .and. ilorr.eq.jlorr(i))then
                    v = fobs(i)
                    dv = fobserr(i)
                    if(.not.doerr)dv = 0.0
                    if(inside(t,v,xmin,xmax,ymin,ymax))then
                    call dotran(doxlog,xmin,xmax,t,
     1                  x0,xlen,xx,
     1                  doylog,ymin,ymax,v,dv,
     1                  y0,ylen,yy,yyp,yym)
                    call gsolid(xx,yy,0.03,jmode(i)+1)
                    if(doblack)then
                    call newpen(1)
                    call gpoly(xx,yy,solsiz,jmode(i)+1)
                    call newpen(kolor)
                    endif
                    if(doerr)then
                        call plot(xx,yym,3)
                        call plot(xx,yyp,2)
                    endif
                    endif
                endif
 2000       continue
        endif
        call pend()
        end
c
c-----------------------------------------------------------------
c
        subroutine start(dobox,iporg,ilorr,x0,y0,xlen,ylen,
     1      xmn,xmx,ymn,ymx,
     2      doxlog,doylog, dofreq,dokmonly,doticonly)
        parameter (LIN=5,LOT=6,LER=0)

        logical dobox, doxlog, doylog, dofreq,dokmonly,doticonly
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
                call pinitf('SLEGNC.PLT')
                if(dokmonly)then
                    ylabel = 'km/s'
                else
                    ylabel = 'C (km/s)'
                endif
            else if(iporg.eq.2)then
                call pinitf('SLEGNU.PLT')
                if(dokmonly)then
                    ylabel = 'km/s'
                else
                    ylabel = 'U (km/s)'
                endif
            else if(iporg.eq.3)then
                call pinitf('SLEGNG.PLT')
                ylabel = 'Gamma (1/km)'
            else if(iporg.eq.4)then
                call pinitf('SLEGNA0.PLT')
                ylabel = 'AL'
            else if(iporg.eq.5)then
                call pinitf('SLEGNAC.PLT')
                ylabel = 'AL sqrt(c)'
            endif
        else if(ilorr.eq.2)then
            if(iporg.eq.1)then
                call pinitf('SREGNC.PLT')
                if(dokmonly)then
                    ylabel = 'km/s'
                else
                    ylabel = 'C (km/s)'
                endif
            else if(iporg.eq.2)then
                call pinitf('SREGNU.PLT')
                if(dokmonly)then
                    ylabel = 'km/s'
                else
                    ylabel = 'U (km/s)'
                endif
            else if(iporg.eq.3)then
                call pinitf('SREGNG.PLT')
                ylabel = 'Gamma (1/km)'
            else if(iporg.eq.4)then
                call pinitf('SREGNA0.PLT')
                ylabel = 'AR'
            else if(iporg.eq.5)then
                call pinitf('SREGNAC.PLT')
                ylabel = 'AR sqrt(c)'
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
                if(doticonly)then
                call dology(x0+0.0 ,y0+0.0,ylen,ymax,ymin,
     1              0.14,.false.,.true.,.false.,ly,
     2              ylabel)
                else
                call dology(x0+0.0 ,y0+0.0,ylen,ymax,ymin,
     1              0.14,.false.,.true.,.true.,ly,
     2              ylabel)
                endif
                call dology(x0+xlen,y0+0.0,ylen,ymax,ymin,
     1              0.14,.true.,.false.,.false.,ly,
     2              ylabel)
            else
                if(doticonly)then
                call doliny(x0+0.0 ,y0+0.0,ylen,ymax,ymin,
     1              0.14,.false.,.true.,.false.,ly,
     2              ylabel)
                else
                call doliny(x0+0.0 ,y0+0.0,ylen,ymax,ymin,
     1              0.14,.false.,.true.,.true.,ly,
     2              ylabel)
                endif
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

        subroutine limt(dofreq,modmax,iporg)
c-----
c       find the range of mode numbers to be plotted.
c       and also plotting limits
c-----
        parameter (LIN=5,LOT=6,LER=0)
        logical dofreq

        parameter (NL=1000)
        common/modl/d(NL), a(NL), b(NL), rho(NL), qa(NL), qb(NL)

        integer*4 ifunc,kmode

        real*4 t0
        real*4 tp,vp

        common/ctrl/ tmin,tmax,cmin,cmax

        character mname*80
        integer ipar(10)
        real*4 fpar(10)
c-----
c       go to beginning of data file
c-----
        rewind 1
        call getmdl(1,nmax,d,a,b,rho,qa,qb,
     1          nper,depths,depthr,
     1          mname,ipar,fpar)
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
            call gethed(1,ifunc,kmode,t0,ierr)
            if(modmax.eq.0)modmax = kmode
            if(ierr.eq.200)then
                write(LOT,*) 'No ifunc=-1. Data end at perd=',t0
                ifunc=-1
            else if(ierr.eq.1001)then
                go to 1001
            endif
            continue
            if(ifunc.le.0) go to 400
            if(kmode.le.0) go to 100
            do 300 i=1,kmode

        call getegn(1,ifunc,1,wvno,u,gamma,
     1      sur,sdur,suz,sduz,sare,wvnsrc,sur0,
     2      rur,rdur,ruz,rduz,rare,wvnrec,rur0,
     3      sumkr,sumgr,sumgv,ierr)

            if(ierr.eq.200)then
                go to 1000
            else if(ierr.eq.1001)then
                go to 1001
            endif
            if(u.gt.0.0)then
            if(kmode .gt. modmax) modmax=kmode
            if(dofreq)then
                tp=1.0/t0
            else
                tp=t0
            endif
            if(tp.gt.tmax)tmax = tp
            if(tp.lt.tmin)tmin = tp
                if(iporg.eq.1)then
                    vp = (6.2831853/t0)/wvno
                else if(iporg.eq.2)then
                    vp = u
                else if(iporg.eq.3)then
                    vp = gamma
                else if(iporg.eq.4)then
                    vp = rare
                else if(iporg.eq.5)then
                    vp = rare*sqrt((6.2831853/t0)/wvno)
                endif
            
                if(vp.gt.cmax)cmax = vp
                if(vp.lt.cmin)cmin = vp
            endif
  300       continue
        go to 100
  400   continue
        WRITE(LOT,*)'xmin-xmax:',tmin,tmax
        WRITE(LOT,*)'ymin-ymax:',cmin,cmax
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

        subroutine gcmdln(ilorr, dofreq, dobox, xmin, xmax, 
     1      cmin, cmax, 
     1      x0, y0, xlen, ylen, width,
     2      kolor, dotext,iporg,doxlog,dolat,doticonly,dokmonly,
     3      doasc,dosurf,datfil,doerr,doshwmd,lintyp,doblack,
     4      verbose)
        integer*4 ilorr, iporg
        logical dofreq, dobox, dotext, doxlog, dolat, dosurf
        logical doasc
        logical doerr,doshwmd, doblack, doticonly, dokmonly
        logical verbose
        integer lintyp
        real*4 xmin, xmax, cmin, cmax
        real*4 x0, y0, xlen, ylen, width
        integer*4 kolor
        character datfil*(*)

        character names*80
c-----
c       defaults
c-----
        ilorr = 0
        iporg = 0
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
        width = -1
        kolor = 1
        dotext = .false.
        doxlog = .false.
        dolat = .false.
        doasc = .false.
        dosurf = .false.
        doerr = .false.
        doshwmd = .false.
        datfil = ' '
        lintyp = 0
        doblack = .true.
        doticonly = .false.
        dokmonly = .false.
        verbose = .false.

        nmarg = mnmarg()
        i = 0
 1000   continue
            i = i + 1
            if(i.gt.nmarg)go to 2000
            call mgtarg(i,names)
            if(names(1:2).eq.'-L' .and. 
     1              names(1:4).ne.'-LAT')then
                ilorr = 1
            else if(names(1:2).eq.'-R')then
                ilorr = 2
            else if(names(1:2).eq.'-C')then
                iporg = 1
            else if(names(1:2).eq.'-U')then
                iporg = 2
            else if(names(1:2).eq.'-G')then
                iporg = 3
            else if(names(1:3).eq.'-A0')then
                iporg = 4
            else if(names(1:3).eq.'-Ac')then
                iporg = 5
            else if(names(1:4).eq.'-LAT')then
                dolat = .true.
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
            else if(names(1:2).eq.'-W')then
                i = i + 1
                call mgtarg(i,names)
                read(names,'(bn,f20.0)')width
                width = abs(width)
            else if(names(1:2).eq.'-K')then
                i = i + 1
                call mgtarg(i,names)
                read(names,'(bn,i20)')kolor
            else if(names(1:8).eq.'-TICONLY')then
                doticonly = .true.
            else if(names(1:7).eq.'-ONLYKM')then
                dokmonly = .true.
            else if(names(1:6).eq.'-NOBOX')then
                dobox = .false.
            else if(names(1:4).eq.'-TXT')then
                dotext = .true.
            else if(names(1:4).eq.'-ASC')then
                doasc  = .true.
            else if(names(1:5).eq.'-XLOG')then
                doxlog = .true.
            else if(names(1:3).eq.'-DT')then
                i = i + 1
                call mgtarg(i,names)
                if(names(1:6).eq."short")then
                        lintyp = 1
                else if(names(1:5).eq."long")then
                        lintyp = 2
                else if(names(1:6).eq."solid")then
                        lintyp = 0
                endif
            else if(names(1:3).eq.'-DE')then
                i = i + 1
                call mgtarg(i,names)
                datfil = names
                doerr = .true.
            else if(names(1:2).eq.'-D')then
                i = i + 1
                call mgtarg(i,names)
                datfil = names
            else if(names(1:2).eq.'-S')then
                dosurf = .true.
            else if(names(1:2).eq.'-M')then
                doshwmd = .true.
            else if(names(1:5).eq.'-NOBL')then
                doblack = .false.
            else if(names(1:2).eq.'-V')then
                verbose = .true.
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
        write(LER,*)'sdpegn96: ',str
        write(LER,*)'USAGE:',
     1  'sdpegn96 ',
     1  '[-L -R] [-C -U -G -A0 -Ac ] [-FREQ -PER]',
     1  '-XMIN xmin -XMAX xmax -YMIN ymin -YMAX ymax ',
     1  '-X0 x0 -Y0 y0 -K kolor -NOBOX -TXT -ASC -XLOG ',
     1  '-D dispfile -DE dispfile -LAT -NOBLACK -? -h -DT linetype',
     1  '-TICONLY -ONLYKM'

        write(LER,*)
     1  '-L                        Love waves'
        write(LER,*)
     1  '-R                        Rayleigh waves'
        write(LER,*)
     1  '     Note one of -L or -R is required'
        write(LER,*)
     1  '-C                        Phase velocity'
        write(LER,*)
     1  '-U                Group velocity'
        write(LER,*)
     1  '-G        Anelastic attenuation coefficient'
        write(LER,*)
     1  '-A0               AR or AL amplitude factor'
        write(LER,*)
     1  '-Ac               AR sqrt(c) or AL sqrt(c) amplitude factor'
        write(LER,*)
     1  '     Note one of -C  -U or -G is required'
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
     1  '-W width (default 0.001)  width of dispersion line'
        write(LER,*)
     1  '-K kolor   (default 1  )  color for curves'
        write(LER,*)
     1  '-NOBOX     (default false) do not plot axes'
        write(LER,*)
     1  '-TXT       (default false) output text file'
        write(LER,*)
     1  '-ASC       (default false) output ASCII column file'
        write(LER,*)
     1  '-XLOG      (default linear) X axis is logarithmic'
        write(LER,*)
     1  '-DT linetype (default solid) linetype= solid short long '
        write(LER,*)
     1  '-D dispfile (default ignore) plot SURF96 dispersion values'
        write(LER,*)
     1  '-DE dispfile (default ignore) plot SURF96 dispersion/error'
        write(LER,*)
     1  '-NOBLACK    (default black) do not but black around symbol'
        write(LER,*)
     1  '-S         (default ignore) create SLEGN.dsp or SREGN.dsp'
        write(LER,*)
     1  '-TICONLY   (default false) no numebrs on y-axis only tics'
        write(LER,*)
     1  '-ONLYKM    (default false) km/s not U(km/s) C(km/s)'
        write(LER,*)
     1  '-LAT       (default false) use output of slat2d96'
        write(LER,*)
     1  '-?         (default false) online help'
        write(LER,*)
     1  '-h         (default false) online help'
        stop 
        end

        subroutine outtxt(ofile, dofreq, ilorr,doshwmd)
c-----
c       make a listing of all dispersion parameters
c-----
c       ofile   C*(*)   name of output file
c       dofreq  L   .true. output frequency
c                   else period
c       ilorr   I*4 1 Love
c               2 Rayleigh
c-----
        parameter (LIN=5,LOT=6,LER=0)
        character ofile*(*)
        logical dofreq,doshwmd

        parameter (NL=1000)
        common/modl/d(NL), a(NL), b(NL), rho(NL), qa(NL), qb(NL)

        character mname*80
        integer ipar(10)
        real*4 fpar(10)
c-----
c       FORMAT STATEMENTS
c-----

   15    format(6f12.5)
   20    format(//19x,'  LOVE WAVE        MODE #',i3//
     1      4x,'N',4x,' T ',7x,'  F  ',10x,'C',
     2      12x,'U',11x,'ARE',8x,'GAMMA')
   30    format(//24x,'RAYLEIGH WAVE      MODE #',i3/
     1      4x,'N',4x,' T ',7x,'  F  ',10x,'C',
     2      12x,'U',11x,'ARE',8x,'GAMMA',8x,'ELP')
   40    format(i5,2g11.4,2g13.5,2e13.5,g13.5)
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
            call getmdl(1,nmax,d,a,b,rho,qa,qb,
     1          nper,depths,depthr,
     1          mname,ipar,fpar)
        if(doshwmd)then
        call tputmdl(2,nmax,d,a,b,rho,qa,qb,
     1          nper,depths,depthr,
     1          mname,ipar,fpar)
        endif
            call gethed(1,ifunc,kmode,t0,ierr)
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
            write(LOT,15) (d(i),a(i),b(i),rho(i),qa(i),qb(i),i=1,nmax)
c-----
c       print the dispersion values.
c       Since the data files are all modes for each frequency, we will
c       repeatedly rwind the input to make a list of all frequencies
c       for each mode
c-----
            icnt = 0
            if(ifunc.eq.1 .or. ifunc.eq.3)then
                ilorr = 1
            else if(ifunc.eq.2 .or. ifunc.eq.4)then
                ilorr = 2
            endif
  100       continue 
            icnt = icnt + 1
            if(icnt.gt.kmode)go to 1000
c-----
c           now start processing from the beginning
c-----
                rewind 1
                call getmdl(1,nmax,d,a,b,rho,qa,qb,
     1              nper,depths,depthr,
     1          mname,ipar,fpar)
                icntm1 = icnt - 1
                if(ilorr.eq.1) then
                    if(dofreq)then
                        write(2,20) icntm1
                    else
                        write(2,20) icntm1
                    endif
                else
                    if(dofreq)then
                        write(2,30) icntm1
                    else
                        write(2,30) icntm1
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
                call gethed(1,ifunc,nmode,t0,ierr)
                omega = 6.2831853/t0
                if(ierr.eq.1001)go to 2001
                if(ierr.eq.200)then
                    write(LOT,*) 
     1          'No ifunc=-1. Data end at perd=',t0
                    ifunc=-1
                endif
                if(ifunc.le.0) go to 100
                if(nmode.le.0) go to 200
                do 300 i=1,nmode
                call getegn(1,ifunc,1,wvno,u,gamma,
     1              sur,sdur,suz,sduz,sare,wvnsrc,sur0,
     2              rur,rdur,ruz,rduz,rare,wvnrec,rur0,
     3              sumkr,sumgr,sumgv,ierr)
        if(doshwmd)then
        write(6,*)ifunc,nmode,t0
        call tputegn(6,ifunc,1,wvno,u,gamma,
     1      sur,sdur,suz,sduz,sare,wvnsrc,sur0,
     2      rur,rtr,ruz,rtz,rare,wvnrec,rur0,
     2      sumkr,sumgr,sumgv)
        endif
                if(ierr.eq.200)go to 2001
                if(ierr.eq.1001)go to 1001
                if(u.gt.0.0)then
                if(i.eq.icnt)then
                nn = nn + 1
                f0=1./t0
                c=omega/wvno
                if(ilorr.eq.1)then
                    write(2,40) nn,t0,f0,c,u,rare,gamma
                else if(ilorr.eq.2)then
                    write(2,40) nn,t0,f0,c,u,rare,gamma,rur0
                endif
                endif
                endif
  300           continue
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

        subroutine outasc(ofile, dofreq, ilorr)
c-----
c       make a listing of all dispersion parameters
c-----
c       ofile   C*(*)   name of output file
c       dofreq  L   .true. output frequency
c                   else period
c       ilorr   I*4 1 Love
c               2 Rayleigh
c-----
        parameter (LIN=5,LOT=6,LER=0)
        character ofile*(*)
        logical dofreq

        parameter (NL=1000)
        common/modl/d(NL), a(NL), b(NL), rho(NL), qa(NL), qb(NL)

        character mname*80
        integer ipar(10)
        real*4 fpar(10)
c-----
c       FORMAT STATEMENTS
c-----

   40   format(i5,1x,i5,1x,7g13.5)
   20   format(' LMODE NFREQ    PERIOD(S) FREQUENCY(Hz)  C(KM/S) ',
     1  '     U(KM/S)      ENERGY       GAMMA(1/KM)   ')
   30   format(' RMODE NFREQ    PERIOD(S) FREQUENCY(Hz)  C(KM/S) ',
     1  '     U(KM/S)      ENERGY       GAMMA(1/KM)  ELLIPTICITY')
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
            call getmdl(1,nmax,d,a,b,rho,qa,qb,
     1          nper,depths,depthr,
     1          mname,ipar,fpar)
            call gethed(1,ifunc,kmode,t0,ierr)
            if(ierr.eq.200)then
                write(LOT,*)'End of File reading header'
                go to 1000
            endif
            if(ierr.eq.1001)then
                write(LOT,*)'Error reading header'
                go to 1000
            endif
            
c-----
c       print the dispersion values.
c       Since the data files are all modes for each frequency, we will
c       repeatedly rwind the input to make a list of all frequencies
c       for each mode
c-----
            icnt = 0
            if(ifunc.eq.1 .or. ifunc.eq.3)then
                ilorr = 1
                write(2,20)
            else if(ifunc.eq.2 .or. ifunc.eq.4)then
                ilorr = 2
                write(2,30)
            endif
  100       continue 
            icnt = icnt + 1
            if(icnt.gt.kmode)go to 1000
c-----
c           now start processing from the beginning
c-----
                rewind 1
                call getmdl(1,nmax,d,a,b,rho,qa,qb,
     1          nper,depths,depthr,
     1          mname,ipar,fpar)
c-----
c       set up counters
c-----
            nn=0
c------
c     pick up the phase velocities at different frequencies for
c     a particular mode.
c------
  200       continue
                call gethed(1,ifunc,nmode,t0,ierr)
                omega = 6.2831853/t0
                if(ierr.eq.1001)go to 2001
                if(ierr.eq.200)then
                    write(LOT,*) 
     1          'No ifunc=-1. Data end at perd=',t0
                    ifunc=-1
                endif
                if(ifunc.le.0) go to 100
                if(nmode.le.0) go to 200
                do 300 i=1,nmode
        call getegn(1,ifunc,1,wvno,u,gamma,
     1      sur,sdur,suz,sduz,sare,wvnsrc,sur0,
     2      rur,rdur,ruz,rduz,rare,wvnrec,rur0,
     3      sumkr,sumgr,sumgv,ierr)
                if(ierr.eq.200)go to 2001
                if(ierr.eq.1001)go to 1001
                if(u.gt.0.0)then
                if(i.eq.icnt)then
                nn = nn + 1
                f0=1./t0
                c=omega/wvno
                icntm1 = icnt - 1
                if(ilorr.eq.1)then
                    write(2,40)
     1              icntm1,nn,t0,f0,c,u,rare,gamma
                else if(ilorr.eq.2)then
                    write(2,40)
     1              icntm1,nn,t0,f0,c,u,rare,gamma,rur0
                endif
                endif
                endif
  300           continue
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

        subroutine dotran(doxlog,xmin,xmax,xc1,x0,xlen,xx,
     1      doylog,ymin,ymax,yc1,ydv,y0,ylen,yy,yyp,yym)
c------
c       transform from user space to plot space
c----
c       doxlog  L   .true. X-axis is logarithmic
c       xmin    R*4 minimum value of X-axis
c       xmax    R*4 maximum value of X-axis
c       xc1 R*4 X-coordinate
c       x0  R*4 Physical position of xmin
c       xlen    R*4 Length of X-axis
c       xx  R*4 Physical position of point xc1
c       doylog  L   .true. Y-axis is logarithmic
c       ymin    R*4 minimum value of Y-axis
c       ymax    R*4 maximum value of Y-axis
c       yc1 R*4 Y-coordinate
c       ydv R*4 error in y-coordinate
c       y0  R*4 Physical position of xmin
c       ylen    R*4 Length of Y-axis
c       yy  R*4 Physical position of point yc1
c       yyp R*4 Physical position of point yc1 + ydv
c       yym R*4 Physical position of point yc1 - ydv
c-----
        logical doxlog, doylog
        real*4 xmin, xmax, xc1, yc1, x0, y0, xlen, ylen, xx, yy
                if(doylog)then
                     yy = y0 + ylen*alog10(yc1/ymin)
     1                  /alog10(ymax/ymin)
                     yyp = y0 + ylen*alog10((yc1+ydv)/ymin)
     1                  /alog10(ymax/ymin)
                     yym = y0 + ylen*alog10((yc1-ydv)/ymin)
     1                  /alog10(ymax/ymin)
                else
                     yy = y0 + ylen*(yc1 - ymin)/(ymax - ymin)
                     yyp = y0 + ylen*(yc1+ydv - ymin)/(ymax - ymin)
                     yym = y0 + ylen*(yc1-ydv - ymin)/(ymax - ymin)
                endif
                if(doxlog)then
                     xx = x0 + xlen*alog10(xc1/xmin)
     1                  /alog10(xmax/xmin)
                else
                     xx = x0 + xlen*(xc1 - xmin)/(xmax - xmin)
                endif
        return
        end

        function inside(x,y,xmin,xmax,ymin,ymax)
        logical inside
        real*4 x, y, xmin, xmax, ymin, ymax
        if(x.ge.xmin .and. x.le.xmax .and. 
     1          y.ge.ymin .and. y.le.ymax)then
            inside = .true.
        else
            inside = .false.
        endif
        return
        end

        subroutine outdsp(ofile, dofreq, ilorr)
c-----
c       make a listing of all dispersion parameters in a SURF96 format
c-----
c       ofile   C*(*)   name of output file
c       dofreq  L   .true. output frequency
c                   else period
c       ilorr   I*4 1 Love
c               2 Rayleigh
c-----
        parameter (LIN=5,LOT=6,LER=0)
        character ofile*(*)
        logical dofreq

        parameter (NL=1000)
        common/modl/d(NL), a(NL), b(NL), rho(NL), qa(NL), qb(NL)

        character mname*80
        integer ipar(10)
        real*4 fpar(10)
c-----
c       FORMAT STATEMENTS
c-----
   40   format('SURF96 L C T ',i5, 3g14.5 )
   41   format('SURF96 L U T ',i5, 3g14.5 )
   42   format('SURF96 L G T ',i5, 3g14.5 )
   43   format('SURF96 R C T ',i5, 3g14.5 )
   44   format('SURF96 R U T ',i5, 3g14.5 )
   45   format('SURF96 R G T ',i5, 3g14.5 )
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
            call getmdl(1,nmax,d,a,b,rho,qa,qb,
     1          nper,depths,depthr,
     1          mname,ipar,fpar)
            call gethed(1,ifunc,kmode,t0,ierr)
            if(ierr.eq.200)then
                write(LOT,*)'End of File reading header'
                go to 1000
            endif
            if(ierr.eq.1001)then
                write(LOT,*)'Error reading header'
                go to 1000
            endif
            
c-----
c       print the dispersion values.
c       Since the data files are all modes for each frequency, we will
c       repeatedly rwind the input to make a list of all frequencies
c       for each mode
c-----
            icnt = 0
  100       continue 
            icnt = icnt + 1
            if(icnt.gt.kmode)go to 1000
c-----
c           now start processing from the beginning
c-----
                rewind 1
                call getmdl(1,nmax,d,a,b,rho,qa,qb,nper,depths,depthr,
     1              mname,ipar,fpar)
c-----
c       set up counters
c-----
            nn=0
c------
c     pick up the phase velocities at different frequencies for
c     a particular mode.
c------
  200       continue
                call gethed(1,ifunc,nmode,t0,ierr)
                omega = 6.2831853/t0
                if(ierr.eq.1001)go to 2001
                if(ierr.eq.200)then
                    write(LOT,*) 
     1          'No ifunc=-1. Data end at perd=',t0
                    ifunc=-1
                endif
                if(ifunc.le.0) go to 100
                if(nmode.le.0) go to 200
                do 300 i=1,nmode
                call getegn(1,ifunc,1,wvno,u,gamma,
     1              sur,sdur,suz,sduz,sare,wvnsrc,sur0,
     2              rur,rdur,ruz,rduz,rare,wvnrec,rur0,
     3              sumkr,sumgr,sumgv,ierr)
                if(ierr.eq.200)go to 2001
                if(ierr.eq.1001)go to 1001
                if(i.eq.icnt)then
                nn = nn + 1
                f0=1./t0
                c=omega/wvno
                icntm1 = icnt - 1
c-----
c       perhaps define negative error if not defined
c       then beware of plotting program
c-----
                errc = 0.001
                erru = 0.001
                errg = 1.0e-7
                if(u.gt.0.0)then
                if(ifunc.eq.1 .or. ifunc.eq.3)then
                    nnmod = icnt - 1
                    write(2,40)nnmod, t0,c,errc
                    write(2,41)nnmod, t0,u,errc
                    write(2,42)nnmod, t0,gamma,errg
                else if(ifunc.eq.2 .or. ifunc.eq.4)then
                    nnmod = i - 1
                    write(2,43)nnmod, t0,c,errc
                    write(2,44)nnmod, t0,u,erru
                    write(2,45)nnmod, t0,gamma,errg
                endif
                endif
                endif
  300           continue
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

        subroutine tputmdl(lun,mmax,d,a,b,rho,qa,qb,nper,dphs,dphr,
     1      mname,ipar,fpar)
        parameter (NL=1000)
        real *4 d(NL), a(NL), b(NL), rho(NL), qa(NL), qb(NL)
        character mname*80
        integer ipar(10)
        real*4 fpar(10)
            write(lun,*)mname
            write(lun,*)ipar
            write(lun,*)fpar
            write(lun,*)mmax
            write(lun,1)(d(i),a(i),b(i),rho(i),
     1          qa(i),qb(i),i=1,mmax)
            write(lun,*)nper,dphs,dphr
    1   format(4f12.4, 2g11.4)
        return
        end

        subroutine tputegn(lun,lorr,intext,wvno,u,gamma,
     1      sur,sdur,suz,sduz,sare,wvnsrc,sur0,
     2      rur,rtr,ruz,rtz,rare,wvnrec,rur0,
     2      sumkr,sumgr,sumgv)
c-----
c       lun I*4 - Logical unit number
c       lorr    I*4 - 1 = Love, 2 = Rayleigh
c       intext  I*4 - 1 write (lr)eigen85 files
c               - 0 write internal files for slat2d96
c       wvno    R*4 - wavenumber for path
c       u   R*4 - group velocity for path
c       gamma   R*4 - anelastic atenuation coefficient for path
c
c       sur R*4 - UR for Rayleigh or UT for Love at source depth
c       sdur    R*4 - d UR/dz for Rayleigh or d UT/dz 
c                   for Love at source depth
c       suz R*4 - UZ for Rayleigh at source depth
c       sduz    R*4 - d UZ/dz for Rayleigh at source depth
c       sare    R*4 - energy integral at source
c       wvnsrc  R*4 - wavenumber at source
c       sur0    R*4 - free surface ellipticity ur/uz at source
c
c       wvnrec  R*4 - wavenumber at receiver
c       rur R*4 - UR for Rayleigh or UT for Love at receiver depth
c       rtr R*4 - d UR/dz for Rayleigh or d UT/dz 
c                   for Love at receiver depth
c       ruz R*4 - UZ for Rayleigh at receiver depth
c       rtz R*4 - d UZ/dz for Rayleigh at receiver depth
c       rare    R*4 - energy integral at receiver
c       rur0    R*4 - ellipticity at receiver
c
c       sumkr   R*4 - sum of kr for slat2d96
c       sumgv   R*4 - sum of r/u for slat2d96
c       sumgr   R*4 - sum of gamma*r for slat2d96
c-----
        if(intext .ne. 0)then
            if(lorr.eq.1) then
        write(LUN,*)'wvno,u,gamma/sut,sdut,rut,rdut/sale'
                write(lun,*)wvno,u,gamma
                write(lun,*)sur,sdur,rur,rtr
                write(lun,*)sare
            else if(lorr.eq.3) then
                write(lun,*)wvno,u,gamma
                write(lun,*)sur,sdur,rur,rtr
                write(lun,*)wvnsrc,wvnrec,sare,rare
            else if(lorr.eq.2) then
        write(LUN,*)'wvno,u,gamma/sur,sdur,suz,sduz/',
     1  'rur,rtr,ruz,rtz,sur0,rur0,sare'
                write(lun,*)wvno,u,gamma
                write(lun,*)sur,sdur,suz,sduz
                write(lun,*)rur,rtr,ruz,rtz
                write(lun,*)sur0, rur0,sare
            else if(lorr.eq.4) then
                write(lun,*)wvno,u,gamma
                write(lun,*)sur,sdur,suz,sduz
                write(lun,*)rur,rtr,ruz,rtz
                write(lun,*)sur0, rur0,sare,rare
                write(lun,*)wvnsrc,wvnrec
            endif
        else
            if(lorr.eq.1) then
                write(lun,*)wvno,are,u,gamma,ur,dur
CTMP     1          ,arer,wvnor,
CTMP     1              sumkr, sumgr,sumgv,
CTMP     2              wvnsrc,wvnrec
            else if(lorr.eq.2) then
                write(lun,*)wvno,are,u,gamma,ur,dur,uz,duz
CTMP     1          ,arer,wvnor,sumkr, sumgr,sumgv, ur0r,
CTMP     2          wvnsrc,wvnrec.sur0,rur0
            endif
        endif
        return
        end

