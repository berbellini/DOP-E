c----------------------------------------------------------------------c
c                                                                      c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                 c
c      VOLUME III                                                      c
c                                                                      c
c      PROGRAM: SDPRAD96                                               c
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
c       29 JAN 2002 - programs created
c       note I could also use the other eigenfunction filed used by 
c       spulse96 - perhaps key on the ilorr to define file type
c       21 JUL 2002 - added distance sieve through the -DMIN dmin and
c           -DMAX dmax flags
c       28 JUL 2002 - added -A flag, also added text output
c       01 AUG 2002 - add code for computing updated gammav value and 
c           and outputting group velocity in SURF96 format
c       08 JAN 2005 - defined LOT  in subroutine gtprdr
c       11 AUG 2007 - cleaned up code so that it works correctly when the
c                     moment tensor is specified. In addition added a new
c                     source -CRACK which requires -DIP dip -STK dip
c                     direction MW/M0 to define a expanding/closing crack
c       03 DEC 2007 - length of efile changed to be consistent
c----------------------------------------------------------------------c
        implicit none
        integer LER, LIN, LOT
        parameter (LIN=5,LOT=6,LER=0)

        real dip, rake, strike,xmom, x0, y0
        real xmt(3,3), fx, fy, fz, hs, r, period
        real dmin, dmax
        logical verby, doauto
        integer isds, kmode, ilorr
        character dfile*80
        character path*190


c-----
c       internal variables
c-----
        character efile*180
        logical ext
        integer ls, lp
        integer lgstr

        real f1, f2, f3, v1, v2, v3
        integer i
        integer  ierr, ndat
        integer NOBS
        parameter (NOBS=2000)
        real dstarr(NOBS), amparr(NOBS), azarr(NOBS)
        real uarr(NOBS)
        integer modarr(NOBS)
        character*8 sta(NOBS), comp(NOBS)

        integer NX
        parameter(NX=100)
        real per(NX)
        integer nper
        real pmin, pmax
c-----
c       machine dependent initializationc
c-----
        call mchdep()
c-----
c       parse command line parameters
c-----
        call gcmdln(dip,rake,strike,isds,verby,xmom,ilorr,
     1      xmt,fx,fy,fz,hs,r,dfile,x0,y0,period,kmode,dmin,dmax,
     2      doauto,path)
        write(LOT,*)'Moments are in multiples of 1+E20 dyne-cm'
        if(isds.eq.0)then
            call trans(dip,strike,rake,f1,f2,f3,v1,v2,v3)
            write(LOT,*)'Earthquake: mom=',xmom
            write(LOT,*)'Dip =',dip
            write(LOT,*)'Stk =',strike
            write(LOT,*)'Rake=',rake
        else if(isds.eq.1)then
            write(LOT,*)'Mxx =',xmt(1,1)
            write(LOT,*)'Mxy =',xmt(1,2)
            write(LOT,*)'Mxz =',xmt(1,3)
            write(LOT,*)'Myy =',xmt(2,2)
            write(LOT,*)'Myz =',xmt(2,3)
            write(LOT,*)'Mzz =',xmt(3,3)
        else if(isds.eq.2)then
            write(LOT,*)'Explosion: mom=',xmom
        endif
c-----
c       decide whether there are observed data. If there are 
c       observed data,
c       the observed data define the period,  mode and data type, 
c       else we
c       use the given period and kmode and choice of ilorr. 
c-----
        if(dfile.ne.' ')then
        call getobs(dfile, period, kmode,ilorr, modarr, dstarr, amparr, 
     1      uarr,azarr, ndat ,dmin, dmax, pmin, pmax, sta, comp, ierr )
        else
            ndat = 0
        endif
c-----
c       open the output text file
c-----
        if(ilorr.eq.1)then
            open(4,file='SRADL.TXT',status='unknown',
     1          access='sequential',form='formatted')
        else if(ilorr.eq.2)then
            open(4,file='SRADR.TXT',status='unknown',
     1          access='sequential',form='formatted')
        endif
        rewind 4
c-----
c       open the eigenfunction / derivative file
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
c-----
c       initialize graphics
c-----
        if(ilorr.eq.1)then
            call pinitf('SRADL.PLT')
        else if(ilorr.eq.2)then
            call pinitf('SRADR.PLT')
        endif
c-----
c       if automatic plot all radiation patterns for this particular
c       wavetype and mode, else just do one
c-----
        if(doauto.and.dfile.ne.' ')then
c-----
c           get the list of master periods
c-----
            call getper(per,NX,pmin,pmax,nper,pmin,pmax) 
c-----
c           cycle through the periods
c-----
            do 1234 i=0,nper-1
                if(mod(i,3).eq.0)then
                    x0 = 2.0
                else if(mod(i,3).eq.1)then
                    x0 = 5.0
                else if(mod(i,3).eq.2)then
                    x0 = 8.0
                endif
                if(mod(i,6).eq.0.and.i.ne.0)then
                    call frame()
                endif
                if(mod(i,6).ge.0.and.mod(i,6).le.2)then
                    y0 = 6.0
                else if(mod(i,6).ge.3.and.mod(i,6).le.5)then
                    y0 = 2.0
                endif
                period = per(i+1)
c-----
c               get the observed 
c-----
                if(dfile.ne.' ')then
                call getobs(dfile, period, kmode,ilorr, modarr, 
     1              dstarr, amparr, uarr,azarr, ndat ,dmin, 
     1              dmax, pmin, pmax, sta, comp, ierr )
                else
                    ndat = 0
                endif
                call doplt(period, kmode, ilorr, hs, dfile, 
     1          ndat, modarr, dstarr, amparr, azarr, uarr,
     2          r, x0, y0, sta, comp,
     3          xmom, xmt, fx, fy, fz, isds,
     4          f1, f2, f3, v1, v2, v3)
 1234       continue
c-----
c       otherwise a single period/mode plot - we already have the
c       data arrays
c-----
        else
            call doplt(period, kmode, ilorr, hs, dfile, 
     1          ndat, modarr, dstarr, amparr, azarr, uarr,
     2          r, x0, y0, sta, comp,
     3          xmom, xmt, fx, fy, fz, isds,
     4          f1, f2, f3, v1, v2, v3)
        endif
        call pend()
        close (1)
        close (4)
        end

                 
        subroutine doplt(period, kmode, ilorr, hs, dfile, ndat,
     1      modarr, dstarr, amparr, azarr, uarr,
     2      r, x0, y0, sta, comp, xmom, xmt, fx, fy, fz, isds,
     3      f1, f2, f3, v1, v2, v3)
        real period, hs
        integer kmode, ilorr, ndat, isds
        character dfile*(*)
        integer NOBS
        parameter (NOBS=2000)
        real dstarr(NOBS), amparr(NOBS), azarr(NOBS)
        real uarr(NOBS)
        integer modarr(NOBS)
        character*8 sta(NOBS), comp(NOBS)
        real r, x0, y0
        real xmt(3,3), fx, fy, fz, xmom
        real f1, f2, f3, v1, v2, v3

        logical success
        real rare, sare, fur, fdur, fuz, fduz, wvno, gammav
        real fur0, fuz0, ftz0
        integer ipar(10)

        integer i 
        integer NGRN
        parameter (NGRN=21)
        complex xx(NGRN)

        real tampar(NOBS)

        rewind 1
c-----
c       go through the eigenfunction derivative file to find
c       the desired period
c-----
        call gperdr(1, period, kmode, ilorr, hs, ipar, gammav,
     1      rare, sare, wvno, fur, fdur, fuz, fduz, success,
     2      fur0, fuz0, ftz0)

            do 1000 i=1,NGRN
                xx(i) = cmplx(0.0, 0.0)
 1000       continue
c-----
c           if observed data, plot them
c-----
            if(dfile.ne.' ')then
                if(ndat.gt.0)then
c-----
c                   normalize geometrical spreading
c                   correct for gammav
c-----
                    do 2000 i=1,ndat
                    tampar(i)=exp(gammav*dstarr(i))*
     1                  sqrt(dstarr(i)/r)*
     2                  amparr(i)
 2000               continue
                endif
            endif
            if(ilorr.eq.1)then
                call excitl(period,r,xx,ipar,wvno,sare,rare,
     1              fur,fdur,fur0)
                    fuz = 0.0
                    fduz = 0.0
                    fuz0 = 0.0
                    ftz0 = 0.0
            else if(ilorr.eq.2)then
                call excitr(period,r,xx,ipar,wvno,sare,rare,
     1              fur,fdur,fuz,fduz, fur0, fuz0, ftz0)
            endif
            call plotth(x0,y0, isds, f1, f2, f3, v1, v2, v3,
     1          xmt, xmom, fx,fy,fz,obmax,r,ilorr,xx,period,
     2          ndat, tampar, amparr, azarr, dstarr, modarr,
     3          sta, comp,uarr,gammav,kmode)
        return
        end


        subroutine plotth(x0,y0, isds, f1, f2, f3, v1, v2, v3,
     1      xmt, xmom, fx, fy, fz, obmax, r,ilorr, xxx, period,
     2      ndat, tampar, amparr, azarr, dstarr, modarr,
     3          sta, comp,uarr,gammav,kmode)
        implicit none
        real x0, y0, f1, f2, f3, v1, v2, v3, period
        real xmt(3,3), xmom, fx, fy, fz, obmax, r
        integer isds, ilorr
        integer ndat
        real amparr(ndat), azarr(ndat), tampar(ndat), dstarr(ndat)
        real uarr(ndat)
        integer modarr(ndat)
        character*8 sta(ndat), comp(ndat)
        real gammav
        integer kmode

        integer NGRN
        parameter (NGRN=21)
        complex xxx(NGRN)
        real az, amp
        real xx, yy, degrad, thmax, valmax
        integer i
        integer NANG
        parameter (NANG=361)
        real xt(NANG), yt(NANG)
        real z, ft, xxxx
        integer mm
        character ic*2
        integer jpen

        real sum, sumx , sumy, sumxx, sumyy, sumxy
        real suml, sumlx , sumly, sumlxx, sumlyy, sumlxy
        real sumu, sumuu
        real a, da, b, db, u, du
        real ampl

        degrad =  0.017453293

        obmax = 0.0
        do 100 i=1,ndat
            if(tampar(i) .gt.obmax)obmax = tampar(i)
  100   continue
c-----
c        create synthetic from fmech96 formulation
c-----
        thmax = 0.0
        do 1000 i=1,NANG
            az = i -1
            call makamp(az, isds, f1, f2, f3, v1, v2, v3,
     1          xmt, xmom, fx, fy, fz, 
     2          amp, ilorr, xxx, period)
            if(amp.gt.thmax)then
                thmax = amp
            endif
            xt(i) = amp * sin(az*degrad)
            yt(i) = amp * cos(az*degrad)
 1000   continue
        if(thmax .gt. obmax)then
            valmax = thmax
        else
            valmax = obmax
        endif
c-----
c       now rescale for the plot
c-----
        do 2000 i=1,NANG
            xx = xt(i) / valmax + x0
            yy = yt(i) / valmax + y0
            if(i.eq.1)then
                call plot(xx,yy,3)
            else
                call plot(xx,yy,2)
            endif
 2000   continue
c-----
c       plot the observations
c-----
        if(ndat.gt.0)then
            write(4,3)
        endif
    3   format(/'STA     ', 'COMP  ',  'Period(s)   ','L/R',
     1  ' MOD','   U(km/s)',' DIST(km)', '   AZ ', '  Amp (obs)',
     2  ' RDIST(km)','  Amp(corr)', '       Pred')

        sum = 0.0
        sumx = 0.0
        sumxx = 0.0
        sumy = 0.0
        sumyy = 0.0
        sumxy = 0.0
        suml = 0.0
        sumlx = 0.0
        sumlxx = 0.0
        sumly = 0.0
        sumlyy = 0.0
        sumlxy = 0.0
        sumu = 0.0
        sumuu = 0.0
        do 200 i=1,ndat
            call makamp(azarr(i), isds, f1, f2, f3, v1, v2, v3,
     1          xmt, xmom, fx, fy, fz, 
     2          amp, ilorr, xxx, period)

            sum = sum + 1
            sumx = sumx + amp
            sumxx = sumxx + amp*amp
            sumy = sumy + tampar(i)
            sumyy= sumyy + tampar(i)*tampar(i)
            sumxy = sumxy + amp*tampar(i)
c-----
c       output text information using mode=0 for Fundamental, 
c       thus modarr(i)-1
c       output observed and predicted and also goodness of fit flagging
c       those amplitudes  exceeding a factor of 3 or 2 of the predicted.
c        The order is
c       important here
c-----
            ic = '  '
            jpen = 2
            if(tampar(i).gt.9.0*amp)then
                ic = ' 9'
                jpen = 2
            else if(tampar(i).gt.8.*amp.and.tampar(i).le.9.*amp)then
                ic = ' 8'
                jpen = 4
            else if(tampar(i).gt.7.*amp.and.tampar(i).le.8.*amp)then
                ic = ' 7'
                jpen = 4
            else if(tampar(i).gt.6.*amp.and.tampar(i).le.7.*amp)then
                ic = ' 6'
                jpen = 4
            else if(tampar(i).gt.5.*amp.and.tampar(i).le.6.*amp)then
                ic = ' 5'
                jpen = 4
            else if(tampar(i).gt.4.*amp.and.tampar(i).le.5.*amp)then
                ic = ' 4'
                jpen = 4
            else if(tampar(i).gt.3.*amp.and.tampar(i).le.4.*amp)then
                ic = ' 3'
                jpen = 4
            else if(tampar(i).gt.2.*amp.and.tampar(i).le.3.*amp)then
                ic = ' 2'
                jpen = 3
            else if(3.*tampar(i).lt.amp)then
                ic = '-3'
                jpen = 4
            else if(2.*tampar(i).lt.amp.and.3.*tampar(i).ge.amp)then
                ic = '-2'
                jpen = 3
            else
                ic = ' '
                jpen = 2
            endif
            if(ilorr.eq.1)then
            write(4,1)sta(i),comp(i),period,modarr(i)-1,uarr(i),
     1          dstarr(i),azarr(i),amparr(i),r,tampar(i),amp,ic
            else if(ilorr.eq.2)then
            write(4,2)sta(i),comp(i),period,modarr(i)-1,uarr(i),
     1          dstarr(i),azarr(i),amparr(i),r,tampar(i),amp,ic
            endif
            xx = tampar(i) * sin(azarr(i)*degrad) / valmax + x0
            yy = tampar(i) * cos(azarr(i)*degrad) / valmax + y0
            call newpen(jpen)
            call fillit('CI',0.032,xx,yy)
            call newpen(1)
            call curvit('CI',0.032,xx,yy)
c-----
c           recompute the gamma value using only the better ratios
c-----
            if(amp.ge.tampar(i)/3.0.and.amp.le.3.0*tampar(i))then
                suml   = suml + 1.0
                sumlx  = sumlx + dstarr(i)
                sumlxx = sumlxx+ dstarr(i)*dstarr(i)
                ampl   = alog(tampar(i)/amp)
                sumly  = sumly + ampl
                sumlyy = sumlyy+ ampl*ampl
                sumlxy = sumlxy+ ampl*dstarr(i)
c-----
c       build up information for dispersion
c-----
                sumu  = sumu  + uarr(i)
                sumuu = sumuu + uarr(i)*uarr(i)
            endif
            
  200   continue
        write(4,*)' MOM_RAT(',period,')',sumxy/sumxx ,
     1     ' R=', sumxy/sqrt(sumxx*sumyy)
        write(4,*)' GAMMA(',period,')  ASSUMED:',gammav
        if(suml.gt.2)then
        call linreg(a,da,b,db,r,suml,sumlx,sumlxx,sumly,sumlyy,sumlxy)
c-----
c       estimate mean grounp velocity and the sigma of the mean
c-----
        u  = sumu/suml
        du = sqrt(abs(sumuu - suml*u*u)/(suml*(suml-1)))
c-----
c       note that the model is exp ( - gamma r)  so slope is b = -gamma
c-----
        write(4,*)' GAMMA(',period,') COMPUTED:',gammav-b,db
        if(ilorr.eq.1)then
              write(4,11)kmode-1,period,gammav-b,db
              write(4,13)kmode-1,period,u, du
        else if(ilorr.eq.2)then
              write(4,12)kmode-1,period,gammav-b,db
              write(4,14)kmode-1,period,u,du
        endif
   11   format('SURF96 L G X',1x,i3,1x,3g11.4)
   12   format('SURF96 R G X',1x,i3,1x,3g11.4)
   13   format('SURF96 L U X',1x,i3,1x,3g11.4)
   14   format('SURF96 R U X',1x,i3,1x,3g11.4)
        
        endif
c-----
c       use the 1P format to make exponential more scientific 
c       but use 0P to reset the 
c       out of r - FORTRAN
c-----
    1   format(2a8,g10.3,' L ',i3,f10.3,f10.2,f6.0,1Pg11.3,
     1      0Pf10.2,1P2g11.3,1x,a2)
    2   format(2a8,g10.3,' R ',i3,f10.3,f10.2,f6.0,1Pg11.3,
     1      0Pf10.2,1P2g11.3,1x,a2)
        call newpen(1)
c-----
c       annotate the figure
c       o place a horizontal scale indicating amplitude
c       o place a North symbol
c       o write the Period
c-----
c       obtain a nice value for the scale
c-----
        z = alog10(valmax)
        mm = z
        ft = z - mm
        ft = 10.**ft
        xxxx = 1./ft
        if(xxxx.gt.1.) mm = mm - 1
        if(xxxx.gt.1.) xxxx = xxxx / 10.
        ft = mm

c-----
c       plot the symbol N at the top oc the figure
c-----
        call symbol(x0-0.028,y0+1.15,0.10,'N',0.0,1)
        call plot(x0,y0+1.10,3)
        call plot(x0,y0+0.97,2)
c-----
c       plot a horizontal scale
c-----
        call plot(x0+xxxx,y0-1.20,3)
        call plot(x0+xxxx,y0-1.30,2)
        call plot(x0+xxxx,y0-1.25,2)
        call plot(x0-xxxx,y0-1.25,2)
        call plot(x0-xxxx,y0-1.20,3)
        call plot(x0-xxxx,y0-1.30,2)
c-----
c       annotate the scale
c-----
        call symbol(x0-0.35,y0-1.15,0.10,'2 x ',0.0,4)
        call symbol(x0+0.05,y0-1.15,0.10,'10',0.0,2)
        call number(x0+0.26,y0-1.07,0.07,ft,0.0,-1)
        call symbol(x0-0.4,y0-1.50,0.10,'T = ',0.0,4)
        call number(x0,y0-1.50,0.10,period,0.0,1)
        return
        end

        subroutine gperdr(lun, period, kmode, ilorr, hs, ipar, gammav,
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
        integer LER, LIN, LOT
        parameter (LIN=5,LOT=6,LER=0)
        integer lun, kmode, ilorr
        real period, hs, gammav
        logical success
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
        integer ifunc, nmode, ierr
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
c       initialize return
c-----
        success = .false.
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
            call getder(1,lorr,wvno,u,gammav,
     1      sur,sdur,sd2ur,suz,sduz,sd2uz,sare,wvnsrc,sur0,
     2      rur,rtr,ruz,rtz,rare,wvnrec,rur0,
     3      sumkr,sumgr,sumgv,ierr,mmax,dcdh,dcda,dcdb,dcdr,
     4      ur,tur,uz,tuz,ipar)

                if(ierr.eq.200)go to 2001
                if(ierr.eq.1001)go to 1001
                f0=1./t0
                c = omega/wvno
            if(j.eq.kmode. and. abs(period-t0).le. 0.001*period)then
c-----
c               compute the eigenfunctions at the desired depth
c-----  
                call dinter(mmax,lorr,hs,wvno,z,ur,tur,uz,tuz,
     1              fur,fdur,fuz,fduz)
                success = .true.
                if(lorr.eq.1)then
                    fuz = 0.0
                    fduz = 0.0
                endif
                return
            endif
  300       continue
        go to 200
 2001   continue
 1001   continue
 1000   continue
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
     1      x,fx,fy,fz,hs,r,dfile,x0,y0,period,kmode,dmin,dmax,
     2      doauto,path)
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
c       r   R*4 - epicentral distance
c       dfile   Ch* name of observed data file
c       x0  R*4 -
c       y0  R*4 - position of center of figure - radius is 1.0
c       period  R*4 - desired period 
c       kmode   I*4 - desired mode   
c       dmin    R*4 - minimum distance for plot
c       dmax    R*4     - maximum distance for plot
c       doauto  L   - .false. only do one period/mode
c                 .true.       do all data for the L or R type
c                              and automatically adjust the positions
c       path    Ch*180  - path to eigenfunction file, default is ./
c-----
        implicit none
        integer LER, LIN, LOT
        parameter (LER=0, LIN=5, LOT=6)
        real dip, rake, strike,xmom, x0, y0
        real x(3,3), fx, fy, fz, hs, r, period
        real dmin, dmax
        logical verby
        integer isds, kmode, ilorr
        character dfile*(*), path*(*)
        logical doauto

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
        x0 = 1.5
        y0 = 1.75
c create GReens from spulse96
        r = 1000.0
        kmode = 0
        period = 20.0
        ilorr = 2
        dfile = ' '
        dmin  = 0.0
        dmax = 100000.0
        doauto = .false.
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
            else if(names(1:2).eq.'-S')then
                i=i+1
                call mgtarg(i,names)
                call chtofp(names,strike)
                isds = 0
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
            else if(names(1:2).eq.'-R'.and.names(1:3).ne.'-RA')then
                ilorr = 2 
            else if(names(1:3).eq.'-HS')then
                i=i+1
                call mgtarg(i,names)
                call chtofp(names,hs)
            else if(names(1:3).eq.'-DI')then
                i=i+1
                call mgtarg(i,names)
                call chtofp(names,r)
            else if(names(1:5).eq.'-DMIN')then
                i=i+1
                call mgtarg(i,names)
                call chtofp(names,DMIN)
            else if(names(1:5).eq.'-DMAX')then
                i=i+1
                call mgtarg(i,names)
                call chtofp(names,DMAX)
            else if(names(1:3).eq.'-X0')then
                i=i+1
                call mgtarg(i,names)
                call chtofp(names,x0)
            else if(names(1:3).eq.'-Y0')then
                i=i+1
                call mgtarg(i,names)
                call chtofp(names,y0)
            else if(names(1:4).eq.'-PER')then
                i=i+1
                call mgtarg(i,names)
                call chtofp(names,period)
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
            else if(names(1:2).eq.'-V')then
                verby = .true.
            else if(names(1:2).eq.'-A')then
                doauto = .true.
            else if(names(1:3).eq.'-CR' .or. names(1:3).eq.'-cr')then
                docrack = .true.
            else if(names(1:4).eq.'-PAT')then
                i=i+1
                call mgtarg(i,path)
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
        WRITE(6,*) 'DOCRACK'
        WRITE(6,*)' xmom  :',xmom  
        WRITE(6,*)' strike:',strike
        WRITE(6,*)' dip   :',dip   
        WRITE(6,*)' rake  :',rake   
        WRITE(6,*)'MT:',x
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
        write(LER,*)'sdprad96 -DIP Dip -STK Stk -RAKE Rake -M0 Mom -E'
        write(LER,*)'  -MW mw -DIST dist -HS hs -L -R -DMIN dmin'
        write(LER,*)'  -DMAX dmax -XX Mxx -YY Myy -ZZ Mzz ',
     1      '-XY -Mxy -XZ Mxz'
        write(LER,*)'  -YZ Myz  -fx Fx -fy Fy -fz Fz -X0 x0 -Y0 y0'
        write(LER,*)'  -PER period -M mode -O obs -PATH path -CRACK'
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
     1  ' -X0 x0     (default=1.5 ) x-coordinate of center of plot'
        write(LER,*)
     1  ' -Y0 y0     (default=1.75) y-coordinate of center of plot'
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
     1  ' -DMIN dmin   (default 0 km     ) minimum for distance sieve'
        write(LER,*)
     1  ' -DMAX dmax   (default 100000 km) maximum for distance sieve'
        write(LER,*)
     1  ' -A           (default false) Plot all periods in one plot'
        write(LER,*)
     1  ' -CRACK       (default no ) use crack model: Dip is  dip of ',
     2  '                     crack, Strike s dip direction of crack',
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

        subroutine makamp(az, isds, f1, f2, f3, v1, v2, v3,
     1      xmt, xmom, forcex, forcey, forcez, 
     2      amp, lorr, xx, period)
        implicit none
        real az, xmom, forcex, forcey, forcez, amp, period
        real xmt(3,3)
        real f1, f2, f3, v1, v2, v3
        integer isds, lorr

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
        integer j,i

        pi=3.141592653589793
        do 1000 i=1,NGRN
            xx(i) = cmplx(0.0,0.0)
 1000   continue
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
        integer  j,i

        complex cvt(3)

        pi  = 3.141592653589793
        pi4 = pi / 4.0
c
        do 1000 i=1,NGRN
            xx(i) = cmplx(0.0,0.0)
 1000   continue
c-----
c       if source or receiver is in fluid, there is no 
c       far-field SH motion
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

        subroutine getobs(dfile, period, kmode, ilorr, modarr, dstarr,
     1      amparr, uarr,azarr, ndat, dmin, dmax, permin, permax, 
     2      sta, comp, ierr )
c-----
c       read observed data
c----
        implicit none
        character dfile*80
        real period
        integer ilorr, ndat, ierr, kmode
        real dmin, dmax
        integer NOBS
        parameter (NOBS=2000)
        real dstarr(NOBS), amparr(NOBS), azarr(NOBS)
        real uarr(NOBS)
        integer modarr(NOBS)
        character*8 sta(NOBS), COMP(NOBS)
        real permin, permax

        integer ls, lss
        integer lgstr
        logical ext
        character instr*170
        integer mode, lsep, lnobl, lnobl1, j, lorr
        real per, u, du, r, a, amp
        character ic*1
        character tsta*8, tcomp*8
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
        ierr = 0
        ndat = 0
        permin = 1.0e+37
        permax = 0.0
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
            
c-----
c           internally use mode = 1 for fundamental
c-----
c-----
c           only use the data if it is for the correct
c           period and the correct mode
c-----
C     1     ilorr,lorr,mode,kmode,per,u,du,r,a,amp
            mode = mode + 1
            if(abs(period - per).lt. 0.001*period .and.
     1          mode.eq.kmode .and. lorr.eq.ilorr
     2          .and. r.ge.dmin .and. r.le.dmax)then
                ndat = ndat + 1
                modarr(ndat) = mode
                period = per
                sta(ndat) = tsta
                comp(ndat) = tcomp
                dstarr(ndat) = r
                azarr(ndat)  = a
                amparr(ndat) = amp
                uarr(ndat) = u
            endif
            if(mode.eq.kmode .and. lorr.eq.ilorr
     2          .and. r.ge.dmin .and. r.le.dmax)then
                if(per.gt.permax)permax = per
                if(per.lt.permin)permin = per
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

        subroutine getper(wn,NX,pmin,pmax,nper,permin,permax) 
c-----
c       automatically create periods from the [pmin,pmax] limits
c       wn()    R*4 array of periods
c       NX  I   dimension of array
c       pmin    R*4 minimum period
c       pmax    R*4 maximum period
c       nper    I*4 number of periods generated
c       permin  R*4 - minimum period to be used in trace
c       permax  R*4 - maximum period to be used in trace
c----
c       this routine is based on the one in sacmft96 which 
c       generates the
c       dispersion values
c-----
        implicit none
        integer NX, nper
        real wn(NX)
        real pmin, pmax
        real permin, permax

        real ymxlog, ymmin, p, tenpow
        integer nocy, iy, ii, jj
        integer NP
        parameter (NP=40)
        real pfac(NP)
        data pfac/1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 
     1      2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9,
     2      3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8,
     3      5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5/

        nper = 0
c-----
c       determine starting power
c-----
        ymxlog = alog10(pmax)   
        ymmin  = alog10(pmin)   
        nocy = ymxlog - ymmin + 1 
        iy = ymmin
        if(ymmin .lt. 0)iy = iy - 1
        do 100 ii=iy,iy+nocy+2
            tenpow = 10.0**ii
            do 200 jj=1,NP
                p = pfac(jj)*tenpow
                if(p .ge. pmin .and. p .le. pmax
     1              .and. p.ge.permin .and. p.le.permax)then
                    nper = nper + 1
                    wn(nper) = p
                    if(nper.eq.NX)go to 1000
                endif
 200        continue
 100    continue
 1000   continue
        return
        end

        subroutine linreg(a,da,b,db,r,summ,sumx,sumxx,sumy,sumyy,sumxy)
        implicit none
        real a,da,b,db,r,summ,sumx,sumxx,sumy,sumyy,sumxy
        real a11, a12, a21, a22, b1, b2, det
        real sum2
        a11  = summ
        a12  = sumx
        a21  = sumx
        a22  = sumxx
        b1   = sumy
        b2   = sumxy
        det  = a11*a22 - a21*a12
        a    = ( a22*b1 - a12*b2)/det
        b    = (-a21*b1 + a11*b2)/det
        sum2 = sumyy - a*sumy - b*sumxy
        da   = sum2 * sumxx/(det* (summ-2))
        db   = sum2 * summ/(det * (summ-2))
        da   = sqrt(abs(da))
        db   = sqrt(abs(db))
        return
        end

        subroutine  makecrack(x,strike,dip,rake,xmom)
        implicit none
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
