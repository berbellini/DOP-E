c----------------------------------------------------------------------c
c                                                                      c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                 c
c      VOLUME III                                                      c
c                                                                      c
c      PROGRAM: SDPSPC96                                               c
c                                                                      c
c      COPYRIGHT 2002                                                  c
c      R. B. Herrmann                                                  c
c      Department of Earth and Atmospheric Sciences                    c
c      Saint Louis University                                          c
c      221 North Grand Boulevard                                       c
c      St. Louis, Missouri 63103                                       c
c      U. S. A.                                                        c
c                                                                      c
c----------------------------------------------------------------------c
c Revision history:
c       21 JUL 2002 - programs created
c       18 SEP 2003 - modified to permit external YMAX
c           this will permit shell script to plot Rayleigh
c           and Love to the same scale
c       09 JAN 2004 - corrected getobs so that if we want to plot BRD
c               we do not also plot BRDG (Hyun-Jae Yoo hjyoo5@snu.ac.kr)
c       09 JAN 2005 - defined LOT in subroutine gperam
c       11 AUG 2007 - cleaned up code so that it works correctly when the
c                     moment tensor is specified. In addition added a new
c                     source -CRACK which requires -DIP dip -STK dip 
c                     direction MW/M0 to define a expanding/closing crack
c----------------------------------------------------------------------c
        implicit none
        integer LER, LIN, LOT
        parameter (LIN=5,LOT=6,LER=0)

        real dip, rake, strike,xmom
        real xmt(3,3), fx, fy, fz, hs, rnorm
        logical verby
        integer isds, kmode, ilorr
        character dfile*80
        character path*180
        character sta*8, comp*8
        logical dofreq, dobox, doxlog, doylog
        real xmin, xmax, ymin, ymax, x0, y0, xlen, ylen 
        integer ndat
        integer ierr

c-----
c       internal variables
c-----
        character efile*80
        logical ext
        integer ls, lp
        integer lgstr

        integer i
        integer NOBS
        parameter (NOBS=1000)
        real perarr(NOBS), amparr(NOBS), az, dist
        integer modarr(NOBS)
        integer kolor

c-----
c       the number of theoretical can be very large - however
c       in practice this will never be large since the periods
c       must be the periods of from multiple filter analysis
c-----
        integer NTHEO
        parameter (NTHEO=8192)
        real tparr(NTHEO), taarr(NTHEO)
        integer ntdat

        real ampmin, ampmax
        real permin, permax
        character ostr*80
c-----
c       machine dependent initializationc
c-----
        call mchdep()
c-----
c       parse command line parameters
c-----
        call gcmdln(dip,rake,strike,isds,verby,xmom,ilorr,
     1      xmt,fx,fy,fz,hs,rnorm,dfile,kmode,
     1      sta,comp,dofreq, dobox, xmin, xmax, ymin, ymax,
     1          x0, y0, xlen, ylen,
     2          kolor, doxlog,doylog,path)
C       WRITE(6,*)'dip,rake,strike,isds,verby,xmom,ilorr:',
C     1     dip,rake,strike,isds,verby,xmom,ilorr
C       WRITE(6,*)'fx,fy,fz,hs,r,x0,y0,kmode:',
C     1     fx,fy,fz,hs,r,x0,y0,kmode
        WRITE(6,*)'sta,comp,dofreq, dobox, xmin, xmax, ymin, ymax:',
     1      sta,comp,dofreq, dobox, xmin, xmax, ymin, ymax
        WRITE(6,*)'kolor, doxlog,doylog,path:',
     1      kolor, doxlog,doylog,path
        WRITE(6,*)'x0, y0, xlen, ylen:',
     1      x0, y0, xlen, ylen
c-----
c       decide whether there are observed data. 
c       If there are observed data,
c       the observed data define the period,  mode and data type, 
c       else we
c       use the given period and kmode and choice of ilorr. 
c-----
        if(dfile.ne.' ')then
        call getobs(dfile, kmode,ilorr, modarr, perarr, amparr, 
     1      az, dist, ndat, sta, comp, ierr )
C           WRITE(6,*)ndat
C           WRITE(6,*)(i,az,amparr(i),i=1,ndat)
        else
            ndat = 0
        endif
        if(rnorm.lt.0.0)then
            if(ndat.eq.0)then
            call usage('No observed data and -DIST dist missing')
            else
                rnorm = dist
            endif
        endif
            
c-----
c       go through the observed data to get the amplitude bounds
c-----
        ampmax = 0.0
        ampmin = 1.0e+38
        permax = 0.0
        permin = 1.0e+38
        do 1000 i=1,ndat
            if(amparr(i).gt.ampmax)then
                ampmax = amparr(i)
            endif
            if(amparr(i).lt.ampmin)then
                ampmin = amparr(i)
            endif
            if(perarr(i).gt.permax)then
                permax = perarr(i)
            endif
            if(perarr(i).lt.permin)then
                permin = perarr(i)
            endif
c------
c           normalize for surface wave spreading to rnorm
c-----
            amparr(i) = amparr(i)*sqrt(dist/rnorm)
C       WRITE(6,*)i,perarr(i),amparr(i)
 1000   continue
C       WRITE(6,*)'permin,ampmin,permax,ampmax:',
c       permin,ampmin,permax,ampmax
c-----
c       go through the eigenfunction file to create a list of
c       predicted amplitudes for this mechanizm 
c       at the reference distance - gamma is included here
c       noting always the extremal amplitudes
c-----
        lp =lgstr(path)
        if(ilorr.eq.1)then
            efile = path(1:lp)//'slegn96.der'
        else if(ilorr.eq.2)then
            efile = path(1:lp)//'sregn96.der'
        endif
        inquire(file=efile,exist=ext)
        ls = lgstr(efile)
        if(.not.ext)then
            call usage('Eigenfunction file '//efile(1:ls)//
     1      ' does not exist')
        endif
        open(1,file=efile,status='unknown',
     1      access='sequential',form='unformatted')
        rewind 1
        call gperam(1, kmode, ilorr, hs, ntdat,tparr, taarr, rnorm,
     1      isds, dip, rake, strike, xmt, xmom, fx,fy,fz,az,dist)
c-----
        close(1)
c       continue searching for the extrema
c-----
C       WRITE(6,*)'ntdat:',ntdat
        do 1100 i=1,ntdat
            if(taarr(i).gt.ampmax)then
                ampmax = taarr(i)
            endif
            if(taarr(i).lt.ampmin)then
                ampmin = taarr(i)
            endif
            if(tparr(i).gt.permax)then
                permax = tparr(i)
            endif
            if(tparr(i).lt.permin)then
                permin = tparr(i)
            endif
            WRITE(6,*)i,tparr(i),taarr(i)
 1100   continue
        WRITE(6,*)'permin,ampmin,permax,ampmax:',
     1      permin,ampmin,permax,ampmax
c-----
c       initialize the graphics
c-----
        if(ilorr.eq.1)then
            call pinitf('SSPCL.PLT')
        else if(ilorr.eq.2)then
            call pinitf('SSPCR.PLT')
        endif
c-----
c       do the plot
c-----
        xmin = 1
        xmax = 1000.0
        if(ymax .lt. 0.0)then
c-----
c           base the plot on the actual data
c-----
            ymax = 1.5*ampmax
        else
c-----
c           base the plot on the externally given maximum
c-----
            ymax = 1.5*ymax
        endif
        ymin = ymax / 1000.0
        call pltspc(dofreq, dobox, xmin, xmax, ymin, ymax, 
     1      x0, y0, xlen, ylen, kolor, doxlog,doylog,
     2      tparr, taarr, ntdat, .true., .true.)
        call pltspc(dofreq, dobox, xmin, xmax, ymin, ymax, 
     1      x0, y0, xlen, ylen, kolor, doxlog,doylog,
     2      perarr, amparr, ndat, .false., .false.)
c-----
c       annotate
c-----
        if(ilorr.eq.1)then
            call gleft(1.0,0.1,0.1,'Love', 0.0)
        else
            call gleft(1.0,0.1,0.1,'Rayleigh', 0.0)
        endif
        write(ostr,1)STA,COMP,dist,az,rnorm
    1   format(a8,1x,a8,1x,'Dist=',f10.2,' Az=',f5.0,' Rnorm=',f10.2)
        call gleft(2.0,0.1,0.10,ostr,0.0)
        
c-----
c       terminate the plot
c-----
        call pend()
        end
        


        subroutine gperam(lun, kmode, ilorr,hs,ntdat,tparr,taarr,rnorm,
     1      isds, dip, rake, strike, xmt, xmom, fx,fy,fz,az,dist)
c-----
c       build up and array of predicted spectral amplitudes at the
c       reference distance
c-----
c       lun I   logical unit for IO
c       kmode   I   target mode number 1 = fundamental internally
c       ilorr   I   1 - Love, 2 - Rayleigh
c       hs  R   Desired depth in km
c       ntdat   I   Number of amplitude period points computed
c       tparr   R   Array of theoretical periods
c       taarr   R   Array of theoretical amplitudes
c               [ Note this is a reference distance and with gamma ]
c       rnorm   R   Normalization distance for geometrical spreading
c       isds    I   0 strike, dip slip, 1 moment tensor, 2 explosion
c       xmt I   Moment tensor elements
c       fx,fy,fz R  Point force elements
c       az  R   station azimuth
c       dist    R   actual distance of station
c-----
        integer LER, LIN, LOT
        parameter (LIN=5,LOT=6,LER=0)
        integer lun
        integer kmode
        integer ilorr
        real hs
        integer ntdat
        real tparr(ntdat), taarr(ntdat)
        real rnorm
        integer isds
        real dip, rake, strike,xmom
        real xmt(3,3), fx, fy, fz
        real az
        real dist
        
        real fur, fdur, fuz, fduz, gamma

        integer NGRN
        parameter (NGRN=21)
        complex xx(NGRN)
        integer ierr
        real f1, f2, f3, v1, v2, v3
c-----
c       getmdl
c-----
        parameter (NL=1000)
        common/modl/d(NL), a(NL), b(NL), rho(NL), qa(NL), qb(NL)
        real d, a, b, rho, qa, qb

        integer*4 mmax,nper
        real depths, depthr

        character mname*80
        real*4 fpar(10)
        integer*4 ipar(10)
c-----
c       gethed
c-----
        integer ifunc, nmode
        real t0
        real*4 z(NL)
c-----
c       getder
c-----
        real f0, c, omega, wvno, u
        real sur,sdur,sd2ur,suz,sduz,sd2uz,sare,wvnsrc,sur0
        real rur,rtr,ruz,rtz,rare,wvnrec,rur0
        real sumkr,sumgr,sumgv
        real*4 dcdh(NL), dcda(NL), dcdb(NL), dcdr(NL)
        real*4 ur(NL), tur(NL), uz(NL), tuz(NL)
c-----
c       internal variables
c-----
        integer lorr
c-----
c       If strike dip or rake given, compute the double couple
c       orientation and normal
c-----
        if(isds.eq.0)then
            call trans(dip,strike,rake,f1,f2,f3,v1,v2,v3)
C       WRITE(0,*)'f1,f2,f3,v1,v2,v3:',
C     1     f1,f2,f3,v1,v2,v3
        endif
c-----
c       begin reading the output file of sdisp96 
c       and do integrity check
c-----
        rewind lun
        call getmdl(lun,mmax,d,a,b,rho,qa,qb,nper,depths,depthr,
     1      mname,ipar,fpar)
        call gethed(lun,ifunc,nmode,t0,ierr)
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
c-----
c       print the dispersion values.
c       Since the data files are all modes for each frequency, we will
c       repeatedly rewind the input to make a list of all frequencies
c       for each mode
c-----
            refdep = fpar(1)
            z(1) = -refdep
            do 100 i=2,mmax
                z(i) = z(i-1) + d(i-1)
  100       continue
c-----
c       now find the desired period
c-----
        rewind lun

        call getmdl(1,mmax,d,a,b,rho,qa,qb,nper,
     1      depths,depthr,mname,ipar,fpar)
c------
c       get the eigenfunction information
c       pick up the phase velocities at different frequencies for
c       a particular mode.
c------
        lorr = ilorr + 4
        ntdat = 0
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
            call getder(1,lorr,wvno,u,gamma,
     1      sur,sdur,sd2ur,suz,sduz,sd2uz,sare,wvnsrc,sur0,
     2      rur,rtr,ruz,rtz,rare,wvnrec,rur0,
     3      sumkr,sumgr,sumgv,ierr,mmax,dcdh,dcda,dcdb,dcdr,
     4      ur,tur,uz,tuz,ipar)

            if(ierr.eq.200)go to 2001
            if(ierr.eq.1001)go to 1001
            f0=1./t0
            c = omega/wvno
            if(j.eq.kmode)then
c-----
c               compute the eigenfunctions at the desired depth
c-----  
                call dinter(mmax,lorr,hs,wvno,z,ur,tur,uz,tuz,
     1              fur,fdur,fuz,fduz)
                if(lorr.eq.1)then
                    fuz = 0.0
                    fduz = 0.0
                endif
                do 1111 i=1,NGRN
                    xx(i) = cmplx(0.0, 0.0)
 1111           continue                           
C       WRITE(6,*)'ilorr,t0,rnorm,wvno,sare,rare,fur,fdur,fuz,fduz:',
C     1     ilorr,t0,rnorm,wvno,sare,rare,fur,fdur,fuz,fduz
                if(ilorr.eq.1)then
C       WRITE(6,*)'CALL EXCITL:',ipar
                call excitl(t0,rnorm,xx,ipar,wvno,sare,rare,
     1              fur,fdur,rur)
                else if(ilorr.eq.2)then
C       WRITE(6,*)'CALL EXCITR:',ipar
                call excitr(t0,rnorm,xx,ipar,wvno,sare,rare,
     1              fur,fdur,fuz,fduz, rur, ruz, rtz)
                endif
C           WRITE(6,*)'xx:',xx
                call    makamp(az,isds,f1,f2,f3,v1,v2,v3,ipar,
     1              xmt, xmom, fx, fy, fz, 
     2              amp, ilorr, xx, t0)
                amp = amp * exp(-gamma*dist)
C       WRITE(6,*)'az,isds,xmt,xmom,fx,fy,fz,ipar,ilorr,xx:',
C     1     az,isds,xmt,xmom,fx,fy,fz,ipar,ilorr,xx
C       WRITE(6,*)ntdat, fur,fdur,fuz,fduz,fur0,fuz0,
c       ftz0,gamma,f0,t0,dist,amp,rnorm
                ntdat = ntdat + 1
                tparr(ntdat) = t0
                taarr(ntdat) = amp
            endif
  300       continue
        go to 200
 2001   continue
 1001   continue
 1000   continue
        return
        end

        subroutine gperdr(lun, period, kmode, ilorr, hs, ipar, gamma,
     1      rare, sare, wvno, fur, fdur, fuz, fduz, success,
     2      rur, ruz, rtz)
c-----
c       get the eigenfunction information from a
c       sregn96 -DER or slegn96 -DER data file
c
c       Get the necessary eigenfunction information
c       at period, mode, depth for ilorr wave
c       return success =.true. is successful
c-----
        integer lun, kmode, ilorr
        real period, hs, gamma
        logical success
c-----
c       initialize return
c-----
        success = .false.
c-----
        return
        end

        subroutine dinter(mmax,lorr,hs,wvno,z,ur,tur,uz,tuz,
     1              fur,fdur,fuz,fduz)
        implicit none
        integer NL
        parameter (NL=1000)
        common/modl/d(NL), a(NL), b(NL), rho(NL), qa(NL), qb(NL)
        real d, a, b, rho, qa, qb
        real*4 ur(NL), tur(NL), uz(NL), tuz(NL)
        real*4 z(NL)
        integer mmax
        integer lorr
        real hs, wvno
        real fur,fdur,fuz,fduz
        real ftur, ftuz
        real zl, zh, p
        real xmu, xlam
        integer i
        do 1000 i=1,mmax-1
            zl = z(i)
            zh = z(i+1)
c-----
c       note use of =< and <  precludes divide by zero
c       
c       We interpolate the stresses since they are continuous
c       then derive the vertical derivative of the
c       eigenfunction
c-----
            if(hs.ge.zl.and.hs.lt.zh)then
                p = (zh-hs)/(zh-zl)
                fur = p*ur(i) + (1.0-p)*ur(i+1)
                ftur = p*tur(i) + (1.0-p)*tur(i+1)
                xmu = rho(i)*b(i)*b(i)
                xlam = rho(i)*(a(i)*a(i) - 2.0*b(i)*b(i))
                if(lorr.eq.6)then
                    fuz = p*uz(i) + (1.0-p)*uz(i+1)
                    ftuz = p*tuz(i) + (1.0-p)*tuz(i+1)
                    fduz = 
     1              (ftuz + wvno*xlam*fur)/(xlam+2.0*xmu)
                    if(b(i).gt.0.0)then
                        fdur = -wvno*fuz + ftur/xmu
                    else
                        fdur = wvno * fuz
                    endif

                else
                    fdur = ftur / xmu
                endif
                
            endif
 1000   continue
        return
        end

        subroutine gcmdln(dip,rake,strike,isds,verby,xmom,ilorr,
     1      x,fx,fy,fz,hs,rnorm,dfile,kmode,
     1      sta,comp,dofreq, dobox, xmin, xmax, ymin, ymax,
     1          x0, y0, xlen, ylen,
     2          kolor, doxlog,doylog,path)
c-----
c       parse command line arguments
c
c       requires subroutine mgtarg() and function mnmarg()
c
c-----
c       dip R*4 - dip of fault
c       rake    R*4 - rake angle of motion on fault
c       strike  R*4 - strike of fault
c       isds    I*4 - indicator of couple source description
c                 -1 none given
c                  0 strike, dip, rake
c                  1 moment tensor
c                  2 explosion
c                    crack - requires dip and strike and Mw/Mo
c                     actually this is the moment tensor solution
c       xmom    R*4 - seismic moment of strike, dip, rake form
c                 1.0 == 1.0 dyne-cm for km,gm,km/s system
c                 true, output radial and transverse motion
c       ilorr   I*4 - 1 Love wave
c                 2 Rayleigh wave
c       x(3,3)  R*4 - moment tensor components (units of dyne-cm)
c       fx  R*4 - point force components
c       fy  R*4   1.0 == 1.0  dyne for km,gm,km/s system
c       fz  R*4
c       hs  R*4 - source depth
c       rnorm   R*4 - epicentral distance
c       dfile   Ch* name of observed data file
c       x0  R*4 -
c       y0  R*4 - position of center of figure - radius is 1.0
c       kmode   I*4 - desired mode   - dfile has precedence
c       sta Ch*8    - station name required
c       comp    Ch*8    - component name required
c       path    Ch*180  - path to eigenfunction file, default is ./
c-----
        implicit none
        integer LER, LIN, LOT
        parameter (LER=0, LIN=5, LOT=6)
        real dip, rake, strike,xmom
        real x(3,3), fx, fy, fz, hs, rnorm
        logical verby
        integer isds, kmode, ilorr, kolor
        character dfile*(*), path*(*)
        character sta*(*), comp*(*)
        logical dofreq, dobox, doxlog, doylog
        real xmin, xmax, ymin, ymax, x0, y0, xlen, ylen 

        character*25 names
        integer*4 mnmarg
        integer nmarg
        integer i, j
        logical docrack
c-----
c       initialize variables
c-----
        dip = 0.0
        rake = 0.0
        strike = 0.0
        isds=-1
        xmom=1.0
        fx = 0.0
        fy = 0.0
        fz = 0.0
        hs = 0.0
        dofreq = .false.
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
        doxlog = .true.
        doylog = .true.     
c create GReens from spulse96
        rnorm = -12345.0
        kmode = 0.0
        ilorr = 2
        dfile = ' '
        sta = ' '
        comp= ' '
        docrack = .false.
        path = './'
        do 120 i=1,3
            do 121 j=1,3
                x(j,i) = 0.0
  121       continue
  120   continue
c-----
c       process command line arguments
c-----
        nmarg=mnmarg()
        if(nmarg.le.0)then
            call usage(' ')
        endif
        i = 0
   11   i = i + 1
        if(i.gt.nmarg)go to 13
            call mgtarg(i,names)
            if(names(1:2).eq.'-D'.and.names(1:5).ne.'-DIST'
     1          .and.names(1:3).ne.'-DM')then
                i=i+1
                call mgtarg(i,names)
                call chtofp(names,dip)
                isds = 0
            else if(names(1:5).eq.'-DIST')then
                i=i+1
                call mgtarg(i,names)
                call chtofp(names,rnorm)
            else if(names(1:2).eq.'-S'.and.
     1          names(1:4).ne.'-STA')then
                i=i+1
                call mgtarg(i,names)
                call chtofp(names,strike)
                isds = 0
            else if(names(1:4).eq.'-STA')then
                i=i+1
                call mgtarg(i,sta)
            else if(names(1:5).eq.'-COMP')then
                i=i+1
                call mgtarg(i,comp)
            else if(names(1:3).eq.'-RA')then
                i=i+1
                call mgtarg(i,names)
                call chtofp(names,rake)
                isds = 0
            else if(names(1:3).eq.'-M0' .or. names(1:3).eq.'-MO')then
                i=i+1
                call mgtarg(i,names)
                call chtofp(names,xmom)
            else if(names(1:3).eq.'-MW' .or. names(1:3).eq.'-Mw')then
                i=i+1
                call mgtarg(i,names)
                call chtofp(names,xmom)
                xmom = 10.**(1.5*xmom + 16.05)
            else if(names(1:2).eq.'-E')then
                isds = 2
            else if(names(1:2).eq.'-L')then
                ilorr = 1 
            else if(names(1:2).eq.'-R')then
                ilorr = 2 
            else if(names(1:3).eq.'-HS')then
                i=i+1
                call mgtarg(i,names)
                call chtofp(names,hs)
            else if(names(1:3).eq.'-DI')then
                i=i+1
                call mgtarg(i,names)
                call chtofp(names,rnorm)
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
                read(names,'(bn,g20.0)')ymin
            else if(names(1:5).eq.'-YMAX')then
                i = i + 1
                call mgtarg(i,names)
                read(names,'(bn,g20.0)')ymax
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
            else if(names(1:5).eq.'-XLOG')then
                doxlog = .true.
            else if(names(1:5).eq.'-YLOG')then
                doylog = .true.
            else if(names(1:3).eq.'-M'.and. 
     1          (names(1:3).ne.'-M0'.or.names(1:3).ne.'-MO').and.
     2          (names(1:3).ne.'-MW'.or.names(1:3).ne.'-Mw'))then
                i=i+1
                call mgtarg(i,names)
                read(names,'(bn,i10)')kmode
                kmode = kmode + 1
            else if(names(1:2).eq.'-O')then
                i=i+1
                call mgtarg(i,names)
                dfile = names
            else if(names(1:3).eq.'-xx' .or. names(1:3).eq.'-XX')then
                i=i+1
                call mgtarg(i,names)
                call chtofp(names,x(1,1))
                isds = 1
            else if(names(1:3).eq.'-yy' .or. names(1:3).eq.'-YY')then
                i=i+1
                call mgtarg(i,names)
                call chtofp(names,x(2,2))
                isds = 1
            else if(names(1:3).eq.'-zz' .or. names(1:3).eq.'-ZZ')then
                i=i+1
                call mgtarg(i,names)
                call chtofp(names,x(3,3))
                isds = 1
            else if(names(1:3).eq.'-xy' .or. names(1:3).eq.'-XY')then
                i=i+1
                call mgtarg(i,names)
                call chtofp(names,x(1,2))
                x(2,1) = x(1,2)
                isds = 1
            else if(names(1:3).eq.'-xz' .or. names(1:3).eq.'-XZ')then
                i=i+1
                call mgtarg(i,names)
                call chtofp(names,x(1,3))
                x(3,1) = x(1,3)
                isds = 1
            else if(names(1:3).eq.'-yz' .or. names(1:3).eq.'-YZ')then
                i=i+1
                call mgtarg(i,names)
                call chtofp(names,x(2,3))
                x(3,2) = x(2,3)
                isds = 1
            else if(names(1:3).eq.'-fx' .or. names(1:3).eq.'-FX')then
                i=i+1
                call mgtarg(i,names)
                call chtofp(names,fx)
            else if(names(1:3).eq.'-fy' .or. names(1:3).eq.'-FY')then
                i=i+1
                call mgtarg(i,names)
                call chtofp(names,fy)
            else if(names(1:3).eq.'-fz' .or. names(1:3).eq.'-FZ')then
                i=i+1
                call mgtarg(i,names)
                call chtofp(names,fz)
            else if(names(1:3).eq.'-CR' .or. names(1:3).eq.'-cr')then
                docrack = .true.
            else if(names(1:4).eq.'-PAT')then
                i=i+1
                call mgtarg(i,path)
            else if(names(1:2).eq.'-V')then
                verby = .true.
            else if(names(1:2).eq.'-?')then
                call usage(' ')
            else if(names(1:2).eq.'-h')then
                call usage(' ')
            endif
        goto 11
   13   continue
c-----
c       special handling for the crack source. convert it to a moment
c       tensor
c       Nakano, M., and H. Kumagi (2005). Waveform inversion of volcano-seismic 
c            signals assuming possible source geometries,
c            Geophys. Res. Letters 32, L12302, doi:10.1029/2005GL022666.
c       the conversion to moment tensor must be deferred until later in the code
c       since the lambda and mu are required at the source depth
c-----
        if(docrack)then
              isds = 1
              call makecrack(x,strike,dip,rake,xmom)
        endif
c-----
c       convert everything to the units required for the synthetics.
c       The synthetics are generated using KM for distance and
c       layer thickness, KM/SEC for velocity, and GM/CC for density.
c       The conversions below cause accoount for the disjoint between
c       these mixed units and the CM CM/SEC and GM/CC required
c-----

        fx     = fx     / 1.0E+15
        fy     = fy     / 1.0E+15
        fz     = fz     / 1.0E+15
        if(isds.eq.1)then
                x(1,1) = x(1,1) / 1.0E+20
                x(1,2) = x(1,2) / 1.0E+20
                x(1,3) = x(1,3) / 1.0E+20
                x(2,3) = x(2,3) / 1.0E+20
                x(3,3) = x(3,3) / 1.0E+20
                x(2,2) = x(2,2) / 1.0E+20
                x(3,2) = x(2,3)
                x(3,1) = x(1,3)
                x(2,1) = x(1,2)
        else if(isds.ne.1)then
                 xmom   = abs(xmom)   / 1.0E+20
        endif
        return
        end
        
        subroutine usage(ostr)
        implicit none
        character ostr*(*)
        integer LER, LIN, LOT
        parameter (LER=0, LIN=5, LOT=6)
        integer ls
        integer lgstr
        ls = lgstr(ostr)
        if(ostr.ne.' ')write(LER,*)ostr(1:ls)
        write(LER,*)'sdpspc96 -DIP Dip -STK Stk -RAKE Rake -M0 Mom -E'
        write(LER,*)'  -MW mw -DIST dist -HS hs -L -R -STA stanam '
        write(LER,*)'  -COMP cmpnam -PATH path -CRACK'
        write(LER,*)'  -XX Mxx -YY Myy -ZZ Mzz -XY -Mxy -XZ Mxz'
        write(LER,*)'  -YZ Myz  -fx Fx -fy Fy -fz Fz'
        write(LER,*)'  -M mode -O obs',
     1  '-FREQ -PER',
     1  '-XMIN xmin -XMAX xmax -YMIN ymin -YMAX ymax ',
     1  '-X0 x0 -Y0 y0 -K kolor -NOBOX -XLOG -YLOG -? -h'
        write(LER,*)
     1  ' '
        write(LER,*)
     1  ' -DIP Dip               dip of fault plane'
        write(LER,*)
     1  ' -STK Strike            strike of fault plane'
        write(LER,*)
     1  ' -RAKE Rake              slip angle on fault plane'
        write(LER,*)
     1  ' -M0 Moment (def=1.0) Seismic moment in units of dyne-cm'
        write(LER,*)
     1  ' -MW mw               Moment Magnitude  '
        write(LER,*)
     1  ' -E                   Explosion'
        write(LER,*)
     1  ' -DIST dist           Normalization distance (km)'
        write(LER,*)
     1  ' -HS hs               Source depth km'
        write(LER,*)
     1  ' -fx FX -fy Fy -fZ fz  Point force amplitudes ',
     2  ' (N,E,down) in  dynes'
        write(LER,*)
     1  ' -XX Mxx -YY Myy -ZZ Mzz  Moment tensor elements in units of'
        write(LER,*)
     1  ' -XY Mxy -XZ Mxz -YZ Myz    dyne-cm'
        write(LER,*)
     1  '-FREQ      (default false) X-Axis is frequency'
        write(LER,*)
     1  '-PER       (default true)X-Axis is period'
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
     1  '-XLOG      (default linear) X axis is logarithmic'
        write(LER,*)
     1  '-YLOG      (default linear) Y axis is logarithmic'
        write(LER,*)
     1  ' -O observed_data     File with observations single',
     2  ' period and mode data'
        write(LER,*)
     1  ' -PER period  (default 20.0 sec) desired period'
        write(LER,*)
     1  ' -M mode      (default  0) desired mode (0-Fund)'
        write(LER,*)
     1  ' -L           (default  Rayl) Plot Love wave radiation'
        write(LER,*)
     1  ' -R           (default  ) Plot Rayleigh wave radiation'
        write(LER,*)
     1  ' -STA  stanam (Required   ) station name'
        write(LER,*)
     1  ' -COMP cmpnam (Required   ) component name'
        write(LER,*)
     1  ' -PATH path   (default ./ ) path to slegn96.der sregn96.der'
        write(LER,*)
     1  ' -CRACK       (default no ) use crack model: Dip is dip of ',
     2  '                     crack, Strike is dip direction of crack',
     3  '                     rake and Mw or M0, where M0 > 0 and'
        write(LER,*)
     1  '               Rake > 0 for expanding crack',
     2  '                    < 0 for closing crack. '
        write(LER,*)
     1  '               M0= sgn(rake) mu DELTA Volume'
        write(LER,*)
     1  '               Poisson ratio = 0.25'
        write(LER,*)
     1  ' -?                   This online help'
        write(LER,*)
     1  ' -h                   This online help'
        stop
        end

        subroutine chtofp(str,fout)
c------
c       routine to convert string to floating point
c       The E format is accepted as well as free form
c       input
c
c       This assumes that the string is 20 characters long
c-----
        implicit none
        character*(*) str
        real*4 fout
        integer*4 lgstr
        integer i, l
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

        subroutine trans(dip,stk,rake,f1,f2,f3,v1,v2,v3)
            degrad=0.01745329
            sins=sin(stk*degrad)
            coss=cos(stk*degrad)
            sind=sin(dip*degrad)
            cosd=cos(dip*degrad)
            sinf=sin(rake*degrad)
            cosf=cos(rake*degrad)
            a11=cosf*coss + sinf*cosd*sins
            a12=cosf*sins - sinf*cosd*coss
            a13= - sinf*sind
            a21= -sins*sind
            a22= coss*sind
            a23= - cosd
            f1=a11
            f2=a12
            f3=a13
            v1=a21
            v2=a22
            v3=a23
        return
        end

        subroutine makamp(az, isds, f1, f2, f3, v1, v2, v3, ipar,
     1      xmt, xmom, forcex, forcey, forcez, 
     2      amp, lorr, xx, period)
        implicit none
        real az, xmom, forcex, forcey, forcez, amp, period
        real xmt(3,3)
        real f1, f2, f3, v1, v2, v3
        integer isds, lorr, ipar(10)

        real fr(4), fz(4), ft(4), fzh, fth
        real degrad
        real cosa, cos2a, sina, sin2a
        integer i
        integer NGRN
        parameter (NGRN=21)
        complex xx(NGRN)
        complex zzamp, zramp, ztamp, zpamp
c-----
c       generate the predicted amplitude for this azimuth from
c       the basic Green's functions
c-----
        degrad = 3.1415927/180. 0
        cosa = cos(az*degrad)
        sina = sin(az*degrad)
        cos2a = cos(2.0*az*degrad)
        sin2a = sin(2.0*az*degrad)
        if(isds.eq.0)then
            fz(1) = f3 * v3
            fz(2) = (f1*v3+f3*v1)*cosa + (f2*v3+f3*v2)*sina
            fz(3) = (f1*v1-f2*v2)*cos2a + (f1*v2+f2*v1)*sin2a
            fz(4) = 0.0
            fr(1) = fz(1)
            fr(2) = fz(2)
            fr(3) = fz(3)
            fr(4) = 0.0
            ft(1) = 0.0
            ft(2) = (f1*v3+f3*v1)*sina - (f2*v3+f3*v2)*cosa
            ft(3) = (f1*v1-f2*v2)*sin2a - (f1*v2+f2*v1)*cos2a
            ft(4) = 0.0
            do 2001 i=1,4
                fz(i) = fz(i) * xmom
                fr(i) = fr(i) * xmom
                ft(i) = ft(i) * xmom
 2001       continue
        else if(isds.eq.1)then
c-----
c       MOMENT TENSOR SPECIFIED
c-----
            fz(1) = -(xmt(1,1)+xmt(2,2))/6.0 + xmt(3,3)/3.0
            fz(2) = xmt(1,3)*cosa + xmt(2,3)*sina
            fz(3) = 0.5*(xmt(1,1)-xmt(2,2))*cos2a + xmt(1,2)*sin2a
            fz(4) = (xmt(1,1)+xmt(2,2)+xmt(3,3))/3.0
            fr(1) = fz(1)
            fr(2) = fz(2)
            fr(3) = fz(3)
            fr(4) = fz(4)
            ft(1) = 0.0
            ft(2) = -xmt(2,3)*cosa + xmt(1,3)*sina
            ft(3) = 0.5*(xmt(1,1)-xmt(2,2))*sin2a - xmt(1,2)*cos2a
            ft(4) = 0.0
        else if(isds.eq.2)then
c-----
c       EXPLOSION SPECIFIED
c-----
            fz(1) = 0.0
            fz(2) = 0.0
            fz(3) = 0.0
            fz(4) = xmom
            fr(1) = 0.0
            fr(2) = 0.0
            fr(3) = 0.0
            fr(4) = xmom
            ft(1) = 0.0
            ft(2) = 0.0
            ft(3) = 0.0
            ft(4) = 0.0
        endif
        fzh = (forcex*cosa + forcey*sina)
        fth = (forcex*sina - forcey*cosa)
c-----
c       now compute the Green's functions using the code from
c       spulse96
c-----
        zzamp = 0.0
        zramp = 0.0
        ztamp = 0.0
        do 9200 i=1,NGRN
            if(i.eq.1)then
                zzamp=zzamp + fz(1)*xx(i)
            else if(i.eq.2)then
                zramp=zramp + fr(1)*xx(i)
            else if(i.eq.3)then
                zzamp=zzamp + fz(2)*xx(i)
            else if(i.eq.4)then
                zramp=zramp + fr(2)*xx(i)
            else if(i.eq.5)then
                ztamp=ztamp + ft(2)*xx(i)
            else if(i.eq.6)then
                zzamp=zzamp + fz(3)*xx(i)
            else if(i.eq.7)then
                zramp=zramp + fr(3)*xx(i)
            else if(i.eq.8)then
                ztamp=ztamp + ft(3)*xx(i)
            else if(i.eq.9)then
                zzamp=zzamp + fz(4)*xx(i)
            else if(i.eq.10)then
                zramp=zramp + fr(4)*xx(i)
            else if(i.eq.11)then
                zzamp=zzamp + forcez*xx(i)
            else if(i.eq.12)then
                zramp=zramp + forcez*xx(i)
            else if(i.eq.13)then
                zzamp=zzamp + fzh*xx(i)
            else if(i.eq.14)then
                zramp=zramp + fzh*xx(i)
            else if(i.eq.15)then
                ztamp=ztamp + fth*xx(i)
            else if(i.eq.16)then
                zpamp=zpamp + fz(4)*xx(i)
            else if(i.eq.17)then
                zpamp=zpamp + fz(1)*xx(i)
            else if(i.eq.18)then
                zpamp=zpamp + fz(2)*xx(i)
            else if(i.eq.19)then
                zpamp=zpamp + fz(3)*xx(i)
            else if(i.eq.20)then
                zpamp=zpamp + forcez*xx(i)
            else if(i.eq.21)then
                zpamp=zpamp + fzh*xx(i)
            endif
 9200   continue
        if(lorr.eq.1)then
            amp = cabs(ztamp)
        else
            amp = cabs(zzamp)
        endif
c-----
c       put in the STEP source S(omega) = 1 / ( i omega )
c-----
        amp = amp / ( 6.2831853/period)
        return
        end

        subroutine excitr(per,rr,xx,ipar,wvno,ares,arer,
     1      ur,dur,uz,duz, ur0, uz0, tz0)
c-----
c
c     This generates the Z and R components of
c     seismogram for
c      1: 45-deg dip-slip source  2: strike-slip source
c      3: dip-slip source         4: explosion source.
c
c-----
        implicit none
        real per, rr, wvno, ares,arer,ur,dur,uz,duz
        real ur0, uz0, tz0
        integer ipar(10)
        integer NGRN
        parameter (NGRN=21)
        complex xx(NGRN)
        complex cvz(6), cvr(6), cvp(6)

        real d1,v1,w1,u1
        real p1
        real t1,t2,ct1,ct2,st1,st2,xmom
        real pi,pi4,pi34
        real dk(6), dkk(6), dkp(6)
        integer j

        pi=3.141592653589793
c------
        xmom=1.0/sqrt(2.00*pi)
        pi4=pi/4.0
        pi34=3.0*pi/4.0
        do 91 j=1,6
            cvz(j)  = cmplx(0.0,0.0)
            cvr(j)  = cmplx(0.0,0.0)
            cvp(j)  = cmplx(0.0,0.0)
   91   continue


        v1 = sqrt(ares*arer)/sqrt(wvno*rr)
c------
c       DD component.
c------
        d1    =  duz+0.5*wvno*ur
        d1    =  d1*v1
        u1    =  d1*ur0
        w1    =  d1*uz0
        p1    = -d1*tz0
        dk(1) =  w1
        dkk(1)=  u1
        dkp(1)=  p1
c------
c       SS component.
c------
        d1    =  wvno*ur
        d1    =  d1*v1
        u1    =  d1*ur0
        w1    =  d1*uz0
        p1    = -d1*tz0
        dk(2) =  w1
        dkk(2)=  u1
        dkp(2)=  p1
c------
c       DS component.
c------
        d1    =  wvno*uz + dur
        d1    =  d1*v1
        u1    =  d1*ur0
        w1    =  d1*uz0
        p1    = -d1*tz0
        dk(3) = w1
        dkk(3)= u1
        dkp(3)= p1
c------
c       EXPLOSION.
c------
        d1    =  duz - wvno*ur
        d1    =  d1*v1
        u1    =  d1*ur0
        w1    =  d1*uz0
        p1    = -d1*tz0
c-----
        dk(4) = w1
        dkk(4)= u1
        dkp(4)= p1
c------
c       VF component.
c------
        d1    =  uz
        d1    =  d1*v1
        u1    =  d1*ur0
        w1    =  d1*uz0
        p1    = -d1*tz0
        dk(5) =  w1
        dkk(5)=  u1
        dkp(5)=  p1
c------
c       HF component.
c------
        d1    =  ur
        d1    =  d1*v1
        u1    =  d1*ur0
        w1    =  d1*uz0
        p1    = -d1*tz0
        dk(6) =  w1
        dkk(6)=  u1
        dkp(6)=  p1
c------
c       1/sqrt(2 pi).  exp(-i pi/4).  exp(-i 3pi/4).
c       A far-field approximation of Hankel function.
c------
c-----
c       introduce time shift and also phase shift from eigenfunction
c-----
        t1 =-pi4
        t2 =-pi34
        ct1=cos(t1)*xmom
        st1=sin(t1)*xmom
        ct2=cos(t2)*xmom
        st2=sin(t2)*xmom
c-----
c       DD - requires real eigenfunction
c-----
        cvz(1) = cvz(1) - dk(1) *cmplx(ct1,st1)
        cvr(1) = cvr(1) + dkk(1)*cmplx(ct2,st2)
        cvp(1) = cvp(1) - dkp(1)*cmplx(ct1,st1)
c-----
c       DS - requires -i  eigenfunction
c-----
        cvz(3) = cvz(3) + dk(3) *cmplx(st1,-ct1)
        cvr(3) = cvr(3) - dkk(3)*cmplx(st2,-ct2)
        cvp(3) = cvp(3) + dkp(3) *cmplx(st1,-ct1)
c-----
c       SS - requires real eigenfunction
c-----
        cvz(2) = cvz(2) - dk(2) *cmplx(ct1,st1)
        cvr(2) = cvr(2) + dkk(2)*cmplx(ct2,st2)
        cvp(2) = cvp(2) - dkp(2)*cmplx(ct1,st1)
c-----
c       EX - requires real eigenfunction
c
c       Note vr[1,7] = vz[1,7] for a fluid , but vz[1,7] will always
c       be correct for Tz
c-----
        cvz(4) = cvz(4) - dk(4)*cmplx(ct1,st1)
        cvr(4) = cvr(4) + dkk(4)*cmplx(ct2,st2)
        cvp(4) = cvp(4) - dkp(4)*cmplx(ct1,st1)
c-----
c       VF - requires real eigenfunction
c-----
        cvz(5) = cvz(5) - dk(5) *cmplx(ct1,st1)
        cvr(5) = cvr(5) + dkk(5)*cmplx(ct2,st2)
        cvp(5) = cvp(5) - dkp(5)*cmplx(ct1,st1)
c-----
c       HF - requires -i eigenfunction
c-----
        cvz(6) = cvz(6) + dk(6) *cmplx(st1,-ct1)
        cvr(6) = cvr(6) - dkk(6)*cmplx(st2,-ct2)
        cvp(6) = cvp(6) + dkp(6)*cmplx(st2,-ct2)
c------
c       The standard output is impulse response type. i.e.,
c       ground velocity, no matter what source type is specified.
c------
        do 601 j=1,6
            t1 = abs(cvz(j))
            if(t1.le.1.0d-20)cvz(j) = cmplx(0.0d+00,0.0d+00)
            t2 = abs(cvr(j))
            if(t2.le.1.0d-20)cvr(j) = cmplx(0.0d+00,0.0d+00)
            t2 = abs(cvp(j))
            if(t2.le.1.0d-20)cvp(j) = cmplx(0.0d+00,0.0d+00)
  601   continue
        if(ipar(2).eq.0)then
c-----
c       source in solid
c-----
            xx(1) = 2.0*cvz(1)
            xx(2) = 2.0*cvr(1)
            xx(3) =     cvz(3)
            xx(4) =     cvr(3)
            xx(6) =    -cvz(2)
            xx(7) =    -cvr(2)
            xx(9) =     cvz(4)
            xx(10) =     cvr(4)
            xx(11) =     cvz(5)
            xx(12) =     cvr(5)
            xx(13) =     cvz(6)
            xx(14) =     cvr(6)
            if(ipar(3).eq.1)then
c-----
c       receiver in fluid
c-----
                xx(16) =     cvp(4)
                xx(17) = 2.0*cvp(1)
                xx(18) =     cvp(3)
                xx(19) =    -cvp(2)
                xx(20) =     cvp(5)
                xx(21) =     cvp(6)
            else
c-----
c       receiver in solid
c-----
                xx(16) = 0.0
                xx(17) = 0.0
                xx(18) = 0.0
                xx(19) = 0.0
                xx(20) = 0.0
                xx(21) = 0.0
            endif
        else 
c-----
c       source in fluid
c-----
            xx(1) = 0.0
            xx(2) = 0.0
            xx(3) = 0.0
            xx(4) = 0.0
            xx(6) = 0.0
            xx(7) = 0.0
            xx(9) =     cvz(4)
            xx(10)=     cvr(4)
            xx(11)= 0.0
            xx(12)= 0.0
            xx(13)= 0.0
            xx(14)= 0.0
            if(ipar(3).eq.1)then
c-----
c       receiver in fluid
c-----
                xx(16) =    cvp(4)
                xx(17) =    0.0
                xx(18) =    0.0
                xx(19) =    0.0
                xx(20) =    0.0
                xx(21) =    0.0
            else
c-----
c       receiver in solid
c-----
                xx(16) = 0.0
                xx(17) = 0.0
                xx(18) = 0.0
                xx(19) = 0.0
                xx(20) = 0.0
                xx(21) = 0.0
            endif
        endif
            
      return
      end
c
c  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
        subroutine excitl(per,rr,xx,ipar,wvno,ales,aler,
     1      ut,dut,ut0)
c-----
c     This generates the T component of seismograms
c     for  1: dip-slip source     2: strike-slip source.
c---
        real per, rr, wvno, ales, aler, ut, dut , ut0
        integer ipar(10)
        integer NGRN
        parameter (NGRN=21)
        complex xx(NGRN)

        real*8 v1,w1,t1,ct1,st1,pi,pi4
        real*8 dk(3)
        real xmom
        integer  j

        complex cvt(3)

        pi  = 3.141592653589793
        pi4 = pi / 4.0
c
        xx(5) = cmplx(0.0,0.0)
        xx(8) = cmplx(0.0,0.0)
        xx(15) = cmplx(0.0,0.0)
C       WRITE(6,*)'ipar(2),ipar(3):',ipar(2),ipar(3)
C       WRITE(6,*)'per, rr, wvno, ales, aler, ut, dut , ut0:',
C     1     per, rr, wvno, ales, aler, ut, dut , ut0 
c-----
c       if source or receiver is in fluid, 
c       there is no far-field SH motion
c-----
        if(ipar(2).eq.1 .or. ipar(3).eq.1)return
        xmom=1.0/sqrt(2.0*pi)
        do 90 j=1,3
            cvt(j)=cmplx(0.0d+00,0.0d+00)
   90   continue
c-----
c       process this distance
c-----
        v1 = sqrt(ales*aler)/sqrt(wvno*rr)
c------
c       DS component.
c------
            
            w1 = dut*v1
            w1 = ut0*w1
            dk(1) = w1
c------
c       SS component.
c------
            w1 = wvno*ut
            w1 = w1*v1
            w1 = ut0*w1
            dk(2) = w1
c-----
c       HF Component
c-----
            w1 = ut*v1
            w1=  ut0*w1
            dk(3) = w1
c-----
c       Introduce time shift and also phase shift from eigenfunctions
c-----
            t1=pi4
            ct1=cos(t1)
            st1=sin(t1)
c-----
c       TDS - requires -i in front of eigenfunction
c-----
            cvt(1) = cvt(1) + dk(1)*xmom*cmplx(st1,-ct1)
c-----
c       TSS - eigenfunction excitation is pure real
c-----
            cvt(2) = cvt(2) + dk(2)*xmom*cmplx(ct1, st1)
c-----
c       THF - requires -i in front of eigenfunction
c-----
            cvt(3) = cvt(3) + dk(3)*xmom*cmplx(st1,-ct1)
c-----
c       5 - TDS
c       8 - TSS
c       15 - THF
c-----
        xx(5) = -cvt(1)
        xx(8) = -cvt(2)
        xx(15)= -cvt(3)
        return
        end

        subroutine getobs(dfile, kmode, ilorr, modarr, 
     1      perarr, amparr, az, dist, ndat, 
     1      sta, comp, ierr )
c-----
c       read observed data
c----
        implicit none
        character dfile*80
        integer ilorr, ndat, ierr, kmode
        character sta*(*), comp*(*)
        integer NOBS
        parameter (NOBS=1000)
        real perarr(NOBS), amparr(NOBS), az, dist
        integer modarr(NOBS)

        integer ls, lss
        integer lgstr
        logical ext
        character instr*170
        integer mode, lsep, lnobl, lnobl1, j, lorr
        real per, u, du, r, a, amp
        character ic*1
        character tsta*8, tcomp*8
        integer lcomp, lsta, ltsta, ltcomp
c-----
c       first determine if the observation file exists
c-----
        inquire(file=dfile,exist=ext)
        ls = lgstr(dfile)

        if(.not.ext)then
            call usage('Observed Dispersion-Spectra file '
     1      //dfile(1:ls)//' does not exist')
        endif
        open(2,file=dfile,status='unknown',
     1      access='sequential',form='formatted')
        rewind 2
        lcomp = lgstr(comp)
        lsta  = lgstr(sta)
        ierr = 0
        ndat = 0
 1000   continue
        read(2,'(a)',end=2000)instr
        ls = lgstr(instr)
c-----
c       find the pattern MFT96 - if it is not there, then
c       read another line 
c-----
        lsep = index(instr,'MFT96')
        if(lsep.gt.0)then
c-----
c       move pointer to the last character
c-----
            lsep = lsep + 5 -1
            call getblnk(instr,lsep,ls,lnobl)
            ic = instr(lnobl:lnobl)
            if(ic(1:1).eq.'L')then
                lorr = 1
            else if(ic(1:1).eq.'R')then
                lorr = 2
            endif
c-----
c           get the C U A identifier
c-----
            lsep = lnobl
            call getblnk(instr,lsep,ls,lnobl)
c-----
c           now find the mode
c-----
            lsep = lnobl
            call getblnk(instr,lsep,ls,lnobl)
            do 10 j=lnobl,ls
                if(instr(j:j+6).eq.'COMMENT')then
                    lss = j -1
                endif
   10       continue
c-----
c           get the numerical values
c-----
            read(instr(lnobl:lss),*)mode,per,u,du,r,a,amp
c-----
c           get station and component starting with the end of comment
c-----
            call getblnk(instr,lss+9,ls,lnobl)
            tsta = ' '
            do 11 j=lnobl,ls
                if(instr(j:j).ne.' ')then
                    tsta(j-lnobl+1:j-lnobl+1)= instr(j:j)
                else
                    lnobl1 = j
                    go to 12
                endif
   11       continue
   12       continue
            call getblnk(instr,lnobl1,ls,lnobl)
            tcomp = ' '
            do 13 j=lnobl,ls
                if(instr(j:j).ne.' ')then
                    tcomp(j-lnobl+1:j-lnobl+1)= instr(j:j)
                else
                    lnobl1 = j
                    go to 14
                endif
   13       continue
   14       continue
            
            ltcomp = lgstr(tcomp)
            ltsta  = lgstr(tsta)
c-----
c           internally use mode = 1 for fundamental
c-----
c-----
c           only use the data if it is for the correct
c           period and the correct mode
c-----
C     1     ilorr,lorr,mode,kmode,per,u,du,r,a,amp
            mode = mode + 1
            if(sta(1:lsta).eq.tsta(1:ltsta) .and.
     1          comp(1:lcomp).eq.tcomp(1:ltcomp).and.
     1          mode.eq.kmode .and. lorr.eq.ilorr)then
                ndat = ndat + 1
                modarr(ndat) = mode
                perarr(ndat) = per
                az  = a
                dist = r
                amparr(ndat) = amp
            endif
        endif
        go to 1000
 2000   continue
        close(2)
        return
        end

        subroutine getblnk(instr,lsep,ls,lnobl)
c-----
c       determine first non-blank character
c
c       instr   Ch* Character string to be parsed
c       lsep    I*4 index of last non blank character
c       ls  I*4 length of input string
c       lnobl   I*4 index of first non blank character
c-----
        implicit none
        character instr*(*)
        integer lsep,ls,lnobl
        integer igotit, i
        character tab*1
        tab=char(9)
        lnobl = lsep+1
        igotit = 0
        do 1000 i=lsep+1,ls
            if(igotit.eq.0)then
            if(instr(i:i).ne.' ' .and. instr(i:i).ne.tab)then
                lnobl = i
                igotit = 1
            endif
            endif
 1000   continue
        return
        end

c-----
c       plot it
c-----
        subroutine pltspc(dofreq, dobox, xmin, xmax, ymin, ymax, 
     1      x0, y0, xlen, ylen, kolor, doxlog,doylog,
     2      parr, aarr, npts, doaxes, docont)
        logical dofreq, dobox, doxlog, doylog
        real xmin, xmax, ymin, ymax, x0, y0, xlen, ylen
        real parr(npts), aarr(npts)
        integer  npts
        logical doaxes
        logical docont

c-----
c       put up the axes
c-----
        call newpen(1)
        xlow = x0 
        ylow = y0
        xhgh = x0 + xlen
        yhgh = y0 + ylen
        if(doaxes)then
        if(dobox)then
        if(doxlog)then
            if(dofreq)then
                call dologx(x0,y0,xlen,xmax,xmin,0.10,
     1          .false.,.false.,.true.,14,'Frequency (Hz)')
            else
                call dologx(x0,y0,xlen,xmax,xmin,0.10,
     1          .false.,.false.,.true.,10,'Period (s)')
            endif
                call dologx(x0,y0+ylen,xlen,xmax,xmin,0.10,
     1          .true.,.true.,.false., 0,' ')
        else
            if(dofreq)then
                call dolinx(x0,y0,xlen,xmax,xmin,0.10,
     1          .false.,.false.,.true.,14,'Frequency (Hz)')
            else
                call dolinx(x0,y0,xlen,xmax,xmin,0.10,
     1          .false.,.false.,.true.,10,'Period (s)')
            endif
                call dolinx(x0,y0,xlen,xmax,xmin,0.10,
     1          .true.,.true.,.false.,0,' ')
        endif

            if(doylog)then
                call dology(x0,y0,ylen,ymax,ymin,0.07,
     1          .true.,.true.,.true.,12,'Amp (cm-sec)')
                call dology(x0+xlen,y0,ylen,ymax,ymin,0.07,
     1          .false.,.false.,.false.,0,' ')
            else
                call doliny(x0,y0,ylen,ymax,ymin,0.07,
     1          .true.,.true.,.true.,12,'Amp (cm-sec)')
                call doliny(x0,y0,ylen,ymax,ymin,0.07,
     1          .false.,.false.,.false., 0,' ')
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
        endif
c-----
c       for plotting we invoke clipping
c-----
        call gclip('on' , xlow,ylow,xhgh,yhgh)
c-----
c       plot with the correct color
c-----
        call newpen(kolor)
c-----
c       do the plotting
c-----
        ipen = 3
        do 1000 i=1,npts
            if(dofreq)then
                xval = 1.0/parr(i)
            else
                xval = parr(i)
            endif
            yval = aarr(i)
c-----
c       now map into the plot space
c-----
            if(doxlog)then
                xx = x0 + xlen*alog10(xval/xmin)
     1                                  /alog10(xmax/xmin)  
            else
                xx = x0 + xlen*(xval - xmin)/(xmax - xmin)
            endif
            if(doylog)then
c-----
c               beware of log of 0.0
c-----
                if(yval.lt. 0.001*ymin)then
                    yval = 0.001*ymin
                endif
                yy = y0 + ylen*alog10(yval/ymin)
     1                                  /alog10(ymax/ymin)  
            else
                yy = y0 + ylen*(yval - ymin)/(ymax - ymin)
            endif
            if(docont)then
                call plot(xx,yy,ipen)
                ipen = 2
            else
                call newpen(2)
                call fillit('CI',0.032,xx,yy)
                call newpen(1)
                call curvit('CI',0.032,xx,yy) 
            endif
 1000   continue
        call plot(xx,yy,3)
c-----
c       reset the plot state
c-----
        call gclip('off', xlow,ylow,xhgh,yhgh)
        call newpen(1)
        return
        end

        subroutine  makecrack(x,strike,dip,rake,xmom)
        real strike, dip, rake, x(3,3), xmom
c-----
c       dip     R*4      - dip of crack with respect to the horizontal
c       strike  R*4      - strike of crack plane - the plane dips
c                          down to the right
c       rake    R*4      - > 0 opening crack, moment > 0
c                          < 0 closing crack, moment < 0
c        x      R*4      - moment tensor
c-----
c       NOTE THIS ASSUMES THAT lambda = mu, isotropic medium with
c           Poisson ratio = 0.25
c-----
c-----
c       special handling for the crack source. convert it to a moment tensor
c       Nakano, M., and H. Kumagi (2005). Waveform inversion of volcano-seismic 
c            signals assuming possible source geometries,
c            Geophys. Res. Letters 32, L12302, doi:10.1029/2005GL022666.
c       the conversion to moment tensor must be deferred until later in the code
c       since the lambda and mu are required at the source depth
c
c       Note they define a theta as the angle that the plane normal makes with
c       respect to the vertical. By definition, theta = dip
c-----
       real degrad, cs,ss, cd, sd, xm, lamdmu
       real ct, st, s2t, s2s
       degrad = 3.1415927/180.0
       cs = cos(degrad*strike)
       ss = sin(degrad*strike)
       ct = cos(degrad*dip   )
       st = sin(degrad*dip   )

       s2s = 2.*ss*cs
       s2t = 2.*st*ct
c-----
c      opening crack
c-----
       if(rake.gt.0.0)then
              xm = xmom
c-----
c      closing crack
c-----
       else if(rake .le. 0.0)then
              xm = - xmom
       endif
c-----
c      lambda/mu
c-----
       lamdmu = 1.0
       x(1,1) = xm * ( lamdmu + 2.* st*st*cs*cs)
       x(1,2) = xm * ( st*st*s2s )
       x(1,3) = xm * ( s2t*cs )
       x(2,2) = xm * ( lamdmu + 2.*st*st*ss*ss )
       x(2,3) = xm * ( s2t * ss)
       x(3,3) = xm * ( lamdmu + 2.*ct*ct)
       x(2,1) = x(1,2)
       x(3,1) = x(1,3)
       x(3,2) = x(2,3)
       return 
       end
