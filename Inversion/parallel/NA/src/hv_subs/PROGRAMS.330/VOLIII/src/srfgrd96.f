        program srfgrd96
c----------------------------------------------------------------------c
c                                                                      c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                 c
c      VOLUME III                                                      c
c                                                                      c
c      PROGRAM: SRFGRD96                                               c
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
c       09 JAN 2005 - defined LOT in subroutine  gperdr
c       11 AUG 2007 - cleaned up code so that it works correctly when the
c                     moment tensor is specified. In addition added a new
c                     source -CRACK which requires -DIP dip -STK stk MW/M0 to
c                     define a expanding/closing crack - note the stk is the
c                     not the strike of the crack, but the dip direction
c                     which is 90+strike
c       24 OCT 2011 - add a -TWT flag that weights the fit by
c                     multiplying by sqrt(period) . The prupose is go 
c                     increase the influence of long periods more

c----------------------------------------------------------------------c
        implicit none
        integer LER, LIN, LOT
        parameter (LIN=5,LOT=6,LER=0)

c-----
c       parameters for grid search
c-----
        real dipmn, dipmx, ddip
        integer ndip
        real rakmn, rakmx, drak
        integer nrak
        real stkmn, stkmx, dstk
        integer nstk
        real hmn  , hmx  , dh  
        integer nh
        integer jdip, jrak, jstk, jh
        integer jnorm
        real permin, permax
c-----
c       minimum good ness of fit for output
c-----
        real fitmin
c-----
c       distance range sieve
c-----
        real dmin, dmax 
        integer isds
        logical dowt
        character dfile*80
        character pathdr*180

c-----
c       internal variables
c-----
        real h, dip, stk, rak
        integer ls, lp
        integer lgstr

        real f1, f2, f3, v1, v2, v3
        real xmt(3,3), xmom
        real fx, fy, fz
        integer i
        integer NGRN
        parameter (NGRN=21)
        complex xx(NGRN)
        integer ierr
        character ofile*14
c-----
c       define observation arrays
c-----
c       perarr  Period
c       amparr  Observed amplitudes
c       azarr   Observed azimuth
c       prearr  Predicted amplitudes
c       ilrarr  Indicates L ro R type (note assume UZ UT)
c       modarr  Predicted mode number
c       ptrarr  Pointer to theoretical
c-----
        integer NOBS
        parameter (NOBS=50000)

c-----
c       Love
c-----
        real perarl(NOBS), amparl(NOBS), azarl(NOBS), rl(NOBS)
        integer modarl(NOBS)
        integer ndatl
        real rarel(NOBS), sarel(NOBS), furl(NOBS), fdurl(NOBS)
        real fuzl(NOBS), fduzl(NOBS), gammal(NOBS), wvnol(NOBS)
        real fur0l(NOBS), fuz0l(NOBS), ftz0l(NOBS)
        logical gotitl(NOBS) 
        real prampl(NOBS)
        integer iparl(10)
        real xl(NOBS), yl(NOBS)
        complex xxl(NGRN,NOBS)
c-----
c       Rayleigh - vertical
c-----
        real perarr(NOBS), amparr(NOBS), azarr(NOBS), rr(NOBS)
        integer modarr(NOBS)
        integer ndatr
        real rarer(NOBS), sarer(NOBS), furr(NOBS), fdurr(NOBS)
        real fuzr(NOBS), fduzr(NOBS), gammar(NOBS), wvnor(NOBS)
        real fur0r(NOBS), fuz0r(NOBS), ftz0r(NOBS)
        logical gotitr(NOBS) 
        real prampr(NOBS)
        integer iparr(10)
        real xr(NOBS), yr(NOBS)
        complex xxr(NGRN,NOBS)
        
        common/bstsol/dipsv, raksv, stksv, rrsv, rlsv, xmwrsv, xmwlsv,
     1      bestsv
        real dipsv, raksv, stksv, rrsv, rlsv, xmwrsv, xmwlsv
        real bestsv

c-----
c       machine dependent initialization
c-----
        call mchdep()
c-----
c       define a unit moment
c-----
        xmom = 1.0
        fx = 0.0
        fy = 0.0
        fz = 0.0
c-----
c       parse command line parameters
c-----
        call gcmdln(dipmn,dipmx,ddip,rakmn,rakmx,drak,
     1      stkmn,stkmx,dstk,hmn,hmx,dh,dfile,fitmin,
     2      dmin,dmax,pathdr,jnorm,permin,permax,isds,dowt)
        nh = ( hmx - hmn)/dh + 1
        ndip = ( dipmx - dipmn)/ddip + 1
        nrak = ( rakmx - rakmn)/drak + 1
        nstk = ( stkmx - stkmn)/dstk + 1
        WRITE(LER,*)'nh  ,hmn  ,hmx:',nh,hmn,hmx
        WRITE(LER,*)'ndip,dipmn,dipmx:',ndip,dipmn,dipmx
        WRITE(LER,*)'nrak,rakmn,rakmx:',nrak,rakmn,rakmx
        WRITE(LER,*)'nstk,stkmn,stkmx:',nstk,stkmn,stkmx
        WRITE(LER,*)'dfile           :',dfile
        WRITE(LER,*)'fitmin          :',fitmin
        WRITE(LER,*)'dmin            :',dmin
        WRITE(LER,*)'dmax            :',dmax
        lp =lgstr(pathdr)
        WRITE(LER,*)'pathdr          :',pathdr(1:lp)
        WRITE(LER,*)'search permin   :',permin
        WRITE(LER,*)'search permax   :',permax
        WRITE(LER,*)'jnorm           :',jnorm 
        WRITE(LER,*)'Period weighting:',dowt
c-----
c       open control files for plotting routines
c-----
        open(3,file='fmdfit.dat',access='sequential',
     1          form='formatted',status='unknown')
        rewind 3

        open(4,file='FMFIT.sh',access='sequential',
     1          form='formatted',status='unknown')
        rewind 4
        write(4,'(a)')'#!/bin/sh'
        write(4,'(a)')' '
        write(4,'(a)')'#####'
        write(4,'(a)')'#     Output of srgrid96'
        write(4,'(a)')'#     For use with fmmfit'
        write(4,'(a)')'#     For use with fmdfit and fmdfit'
        write(4,'(a)')'#####'



c-----
c       open eigenfunction files
c-----
C       NOTE CHECK FOR EXISTENCE
C-----
c-----
c       decide whether there are observed data. If there are 
c          observed data,
c       the observed data define the period,  mode and data 
c          type, else we
c       use the given period and kmode and choice of ilorr. 
c-----
        if(dfile.ne.' ')then
        call getobs(dfile, modarr, perarr, amparr, azarr,rr,
     1      modarl, perarl, amparl, azarl,rl,
     1      ndatl, ndatr, dmin, dmax, ierr ,permin, permax)
        else
            ndatl = 0
            ndatr = 0
        endif
        WRITE(0,*)'NDATL,NDATR:',NDATL,NDATR
c-----
c       now match up each observed period mode type with a particular
c       observation - since the periods will be repeated we will
c       establish a keyed array
c-----
        lp =lgstr(pathdr)
        open(1,file=pathdr(1:lp)//'slegn96.der',status='unknown',
     1      access='sequential',form='unformatted')
        rewind 1
        open(2,file=pathdr(1:lp)//'sregn96.der',status='unknown',
     1      access='sequential',form='unformatted')
        rewind 2
        WRITE(0,*)'Using ',pathdr(1:lp)//'sregn96.der'
        WRITE(0,*)'Using ',pathdr(1:lp)//'slegn96.der'

        write(4,7)hmn,hmx
    7   format('fmdfit -HMN',f5.0,' -HMX',f5.0, ' < fmdfit.dat')

        do 1000  jh=1,nh
            h = hmn + (jh-1)*dh
            bestsv = 0.0
            write(ofile,9)jh
    9       format('fmfit',i3.3,'.dat')
            open(8,file=ofile,status='unknown',
     1          access='sequential',form='formatted')
            rewind 8
            ls = lgstr(ofile)
            write(4,6)stkmn,stkmx,dipmn,dipmx,jh,ofile(1:ls)
    6   format('fmmfit -SMN',f6.0,' -SMX',f6.0, ' -DMN', f5.0, ' -DMX',
     1      f5.0, ' -ID ',i3.3, ' < ', a)
        
c-----
c           get the eigenfunction values for this
c           depth - this is brutal since I will do this
c           many times for the same period
c-----
            do 2000 i=1,ndatl
                call gperdr(1, perarl(i), modarl(i), 1, h, 
     1              iparl, gammal(i), rarel(i), sarel(i), 
     2              wvnol(i), furl(i), fdurl(i), fuzl(i), 
     3              fduzl(i), gotitl(i),
     4              fur0l(i), fuz0l(i), ftz0l(i))
                if(gotitl(i))then
                call excitl(perarl(i),1000.0,xx,iparl,
     1              wvnol(i),sarel(i),rarel(i),
     2              furl(i),fdurl(i),fur0l(i))
                call cptoxxx(xx,xxl,i,.true.)
                endif
 2000       continue
            do 2010 i=1,ndatr
                call gperdr(2, perarr(i), modarr(i), 2, h, 
     1              iparr, gammar(i), rarer(i), sarer(i), 
     2              wvnor(i), furr(i), fdurr(i), fuzr(i), 
     3              fduzr(i), gotitr(i),
     4              fur0r(i), fuz0r(i), ftz0r(i))
                if(gotitr(i))then
                call excitr(perarr(i),1000.0,xx,iparr,
     1              wvnor(i),sarer(i),rarer(i),
     1              furr(i),fdurr(i),fuzr(i),fduzr(i), 
     1              fur0r(i), fuz0r(i), ftz0r(i))
                call cptoxxx(xx,xxr,i,.true.)
                endif
 2010       continue

        do 1100  jdip=1,ndip
            dip = dipmn + (jdip-1)*ddip
        do 1200  jstk=1,nstk
            stk = stkmn + (jstk-1)*dstk
        do 1300  jrak=1,nrak
            rak = rakmn + (jrak-1)*drak
            if(isds.eq.0)then
                 call trans(dip,stk,rak,f1,f2,f3,v1,v2,v3)
            else if(isds.eq.3)then
                 call makecrack(xmt,stk,dip,rak,xmom)
            endif
c-----
c           form the array of predictions
c-----
            do 3000 i=1,ndatl
                if(gotitl(i))then
                call cptoxxx(xx,xxl,i,.false.)
                call makamp(azarl(i), isds, f1, f2, f3, v1, v2, v3,
     1              xmt, xmom, fx, fy, fz, 
     2              prampl(i), 1, xx, perarl(i))
c-----
c       correct all observed for gamma and normalize for 
c          geometrical spreading
c       bad aspect gives more emphasis to bad observations 
c          at large distance
c-----
C               xl(i) = prampl(i)
C               yl(i) = amparl(i)*sqrt(rl(i)/1000.0)*
c          exp(gammal(i)*rl(i))
c-----
c       correct all amplitude for geometrical spreading
c       propagate theoretical with gamma
c-----
                xl(i) = prampl(i)*exp(-gammal(i)*rl(i))
                yl(i) = amparl(i)*sqrt(rl(i)/1000.0)
                endif
 3000       continue
            do 3010 i=1,ndatr
                if(gotitr(i))then
                call cptoxxx(xx,xxr,i,.false.)
                call makamp(azarr(i), isds, f1, f2, f3, v1, v2, v3,
     1              xmt, xmom, fx, fy, fz, 
     2              prampr(i), 2, xx, perarr(i))
                if(jnorm.eq.1)then
c-----
c       correct all observed for gamma and normalize 
c          for geometrical spreading
c       bad aspect gives more emphasis to bad observations 
c          at large distance
c-----
                    xr(i) = prampr(i)
                    yr(i) = amparr(i)*sqrt(rr(i)/1000.0)
     1                  *exp(gammar(i)*rr(i))   
c-----
c       correct all amplitude for geometrical spreading
c       propagate theoretical with gamma
c-----
                else if(jnorm.eq.2)then
                    xr(i) = prampr(i)*exp(-gammar(i)*rr(i))
                    yr(i) = amparr(i)*sqrt(rr(i)/1000.0)
                endif
                
                endif
 3010       continue
c-----
c-----
c       compute the observed predicted vectors for all stations
c-----
            call good(h,dip,rak,stk,ndatl,yl,xl,gotitl,
     1          ndatr,yr,xr,gotitr,fitmin,perarl,perarr,dowt)
 1300   continue
 1200   continue
 1100   continue
c-----
c       output the best solution for this depth
c-----
        write(3,8)h,stksv,dipsv,raksv,rrsv,rlsv,xmwrsv,xmwlsv,bestsv
    8   format('SRFGRD96 ',f6.1,3f6.0,2f6.3,2f7.2,f7.4) 
            close(8)
 1000   continue
        close(1)
        close(2)
        close(3)
        close(4)
        end
        
        subroutine good(h,dip,rak,stk,
     1      ndatl,amparl,prampl,gotitl,ndatr,amparr,prampr,gotitr,
     2      fitmin,perarl,perarr,dowt)
c-----
c       compute the goodness of fit for this source specification
c-----
c       h   R   - source depth
c       dip R   - source dip
c       rake    R   - source rake
c       stk R   - source strike
c       ndatl   I   - Number of Love observations
c       amparl  R   - array of observed Love wave spectral amplitudes
c       prampl  R   - array of predicted Love wave spectral amplitudes
c       gotitl  L   - array to indicatge valid Love obs-pre pair
c       ndatr   I   - Number of Rayleigh observations
c       amparr  R   - array of observed Rayleigh wave 
c          spectral amplitudes
c       prampr  R   - array of predicted Rayleigh wave 
c          spectral amplitudes
c       gotitr  L   - array to indicatge valid Rayleigh obs-pre pair
c       fitmin  R   - threshold for output of results, output only
c                 if best > fitmin
c       perarl  R   - array of Love wave periods
c       perarr  R   - array of Rayleigh wave periods
c       dowt    L   - .true. weight by sqrt(period)
c-----
        implicit none
        real h, dip, rak, stk
        integer idip, irak, istk
        integer ndatl
        real amparl(ndatl), prampl(ndatl), perarl(ndatl)
        logical gotitl(ndatl)
        integer ndatr
        real amparr(ndatr), prampr(ndatr), perarr(ndatr)
        logical gotitr(ndatr)
        logical dowt
        real fitmin
        real rl, rvarl, scll
        real rr, rvarr, sclr
        real sfac, best
        real xmwl, xmwr
        
        common/bstsol/dipsv, raksv, stksv, rrsv, rlsv, xmwrsv, xmwlsv,
     1      bestsv
        real dipsv, raksv, stksv, rrsv, rlsv, xmwrsv, xmwlsv
        real bestsv
c-----
c       get fit to Love
c-----
        call gfit(ndatl,prampl,amparl,gotitl,rl,rvarl,scll,perarl,dowt)
c-----
c       get fit to Rayleigh
c-----
        call gfit(ndatr,prampr,amparr,gotitr,rr,rvarr,sclr,perarl,dowt)
c-----
c       determine ratio of Moment(Love)/Moment(Rayl) such that this
c       ratio is always < 1
c-----
        if(scll .gt. sclr)then
            sfac = sclr/scll
        else
            sfac = scll/sclr
        endif
c-----
c       define best fit
c-----
        best = rl*rr*sfac
c-----
c       compute moment magnitudes from the definition
c       log Mo = 1.5 Mw + 16.05  (Hanks and Kanamori)
c       Note that the synthetic Green's functions are for a moment of
c       1.0E+20 dyne-cm
c-----
        xmwl = (20.0+ alog10(scll) - 16.05)/1.5
        xmwr = (20.0+ alog10(sclr) - 16.05)/1.5
c-----
c       output
c-----
C       write(8,1)h,dip,rak,stk,rr,rl,xmwr,xmwl,best
C       write(8,*)'            ',sclr,scll,sfac
        if(best.gt.fitmin)then
        idip = dip
        irak = rak
        istk = stk
        write(8,1)h,istk,idip,irak,rr,rl,xmwr,xmwl,best
    1   format('SRFGRD96 ',f6.1,3i6  ,2f6.3,2f7.2,f7.4) 
        endif
        if(best.gt.bestsv)then
            bestsv = best
            dipsv = dip
            raksv = rak
            stksv = stk
            rrsv = rr
            rlsv = rl
            xmwrsv = xmwr
            xmwlsv = xmwl
        endif
        return
        end

        subroutine gfit(n,x,y,gotit,r,rvar,scl,per,dowt)
c-----
c       solve the problem y = a x
c       n       I   - number of observations
c       x       R   - array of x values
c       y       R   - array of y values
c       gotit   L   - array to indicate if (x,y) pair is valid
c       r       R   - vector dot product  x dot y / sqrt( |x| |y| )
c       rvar    R   - reduction of variance = r^2
c       scl     R   - scale factor a
c       per     R   - array of periods
c       dowt    L   - .true. weight by sqrt(period)
c-----
        implicit none
        integer n
        real x(n), y(n), per(n)
        logical gotit(n), dowt
        real r, rvar, scl

        integer i
        real sum, sumx, sumxx, sumy, sumyy, sumxy
        real wt
        sum  = 0.0
        sumx  = 0.0
        sumy  = 0.0
        sumxx = 0.0
        sumxy = 0.0
        sumyy = 0.0
        do 1000 i=1,n
            if(gotit(i))then
                if(dowt)then
                    wt = sqrt(abs(per(i)))
                else
                    wt = 1.0
                endif
                sum = sum + 1
                sumx  = sumx  + wt*x(i)
                sumy  = sumy  + wt*y(i)
                sumxx = sumxx + wt*x(i)*x(i)
                sumyy = sumyy + wt*y(i)*y(i)
                sumxy = sumxy + wt*x(i)*y(i)
            endif
 1000   continue

        scl = sumxy/sumxx
        r = sumxy/sqrt(sumxx*sumyy)
        rvar = 1.0 - (sumyy + scl*scl*sumxx - 2.0*scl*sumxy)/sumyy
c-----
c       Note that mathematically rvar = r * r
c-----
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
        integer LER, LIN, LOT
        parameter (LIN=5,LOT=6,LER=0)
        integer lun, kmode, ilorr
        real period, hs, gamma
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

        call getmdl(lun,mmax,d,a,b,rho,qa,qb,nper,
     1      depths,depthr,mname,ipar,fpar)
c------
c       get the eigenfunction information
c       pick up the phase velocities at different frequencies for
c       a particular mode.
c------
        lorr = ilorr + 4
  200   continue
            call gethed(lun,ifunc,nmode,t0,ierr)
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
            call getder(lun,lorr,wvno,u,gamma,
     1      sur,sdur,sd2ur,suz,sduz,sd2uz,sare,wvnsrc,sur0,
     2      rur,rtr,ruz,rtz,rare,wvnrec,rur0,
     3      sumkr,sumgr,sumgv,ierr,mmax,dcdh,dcda,dcdb,dcdr,
     4      ur,tur,uz,tuz,ipar)

                if(ierr.eq.200)go to 2001
                if(ierr.eq.1001)go to 1001
                f0=1./t0
                c = omega/wvno
            if(j.eq.kmode. and. abs(period-t0).le. 0.01*period)then
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

        subroutine gcmdln(dipmn,dipmx,ddip,rakmn,rakmx,drak,
     1      stkmn,stkmx,dstk,hmn,hmx,dh,dfile,fitmin,
     2      dmin,dmax,pathdr,jnorm,permin,permax,isds,dowt)
c-----
c       parse command line arguments
c
c       requires subroutine mgtarg() and function mnmarg()
c
c-----
c       dipmn,dipmx,ddip    R*4 - dip range search parameters
c       rakmn,rakmx,drak    R*4 - rake range search parameters
c       stkmn,stkmx,dstk    R*4 - strike range search parameters
c       hmn,hmx,dh      R*4 - depth range search parameters
c       dfile           Ch* - name of observed data file
c       fitmin          R*4 - lower bound of fit for output
c                         this makes the output list smaller
c       dmin,dmax       R*4 - distance range sieve
c       pathdr  Ch*180  - path to eigenfunction file, default is ./
c       jnorm           I   - normalization technique
c                         Given O(obs) and T(theo) corrected for
c                         sqrt(dist) to reference distance and T
c                         has no anelastic attenuation effect
c                         1 =  Obs*exp(gamma r) vs T
c                         2 =  Obs              vs T exp(-gamma r)
c                         jnorm 1 is used for radiation pattern plot
c                         jnorm 2 weights down short period and distant 
c                             perhaps better for region where 
c                                             gamma is poorly known
c       permin          R   - minimum period to use in data set
c       permax          R   - maximum period to use in data set
c       isds    I*4 - indicator of couple source description
c                 -1 none given
c                  0 strike, dip, rake
c                  1 moment tensor
c                  2 explosion
c                  3 crack - requires dip and strike and Mw/Mo
c                     actually this is the moment tensor solution
c       dowt            L - .true. weight the observations by sqrt(period)
c                           .false. (default)
c-----
        implicit none
        integer LER, LIN, LOT
        parameter (LER=0, LIN=5, LOT=6)
c-----
c       parameters for grid search
c-----
        real dipmn, dipmx, ddip
        real rakmn, rakmx, drak
        real stkmn, stkmx, dstk
        real hmn  , hmx  , dh  
        real dmin, dmax 
        real permin, permax
        real fitmin
        integer isds
        logical dowt
        character dfile*(*), pathdr*(*)
        integer jnorm
        real tmp

        character*25 names
        integer*4 mnmarg
        integer nmarg
        integer i
c-----
c       initialize variables
c-----
        isds = 0
        dipmn = 30  
        dipmx = 90
        ddip  = 15
        rakmn = -180
        rakmx = 180
        drak  = 15
        stkmn = 0
        stkmx = 350
        dstk  = 10
        hmn   = 1
        hmx   = 30
        dh    = 2
        dfile = ' '
        dmin = 0.0
        dmax = 100000.0
        fitmin = 0.5
        pathdr = './'
        permin = 0.0
        permax = 2000.0
        jnorm = 2
        dowt = .false.
c-----
c       process command line arguments
c-----
        nmarg=mnmarg()
        if(nmarg.le.0)then
            call usage('Requires observation file: -O obs')
        endif
        i = 0
   11   i = i + 1
        if(i.gt.nmarg)go to 13
            call mgtarg(i,names)
            if(names(1:4).eq.'-DMN')then
                i=i+1
                call mgtarg(i,names)
                call chtofp(names,dipmn)
            else if(names(1:4).eq.'-DMX')then
                i=i+1
                call mgtarg(i,names)
                call chtofp(names,dipmx)
            else if(names(1:3).eq.'-DD')then
                i=i+1
                call mgtarg(i,names)
                call chtofp(names,ddip)
            else if(names(1:4).eq.'-RMN')then
                i=i+1
                call mgtarg(i,names)
                call chtofp(names,rakmn)
            else if(names(1:4).eq.'-RMX')then
                i=i+1
                call mgtarg(i,names)
                call chtofp(names,rakmx)
            else if(names(1:3).eq.'-DR')then
                i=i+1
                call mgtarg(i,names)
                call chtofp(names,drak)
            else if(names(1:4).eq.'-SMN')then
                i=i+1
                call mgtarg(i,names)
                call chtofp(names,stkmn)
            else if(names(1:4).eq.'-SMX')then
                i=i+1
                call mgtarg(i,names)
                call chtofp(names,stkmx)
            else if(names(1:3).eq.'-DS')then
                i=i+1
                call mgtarg(i,names)
                call chtofp(names,dstk)
            else if(names(1:4).eq.'-HMN')then
                i=i+1
                call mgtarg(i,names)
                call chtofp(names,hmn)
            else if(names(1:4).eq.'-HMX')then
                i=i+1
                call mgtarg(i,names)
                call chtofp(names,hmx)
            else if(names(1:3).eq.'-DH')then
                i=i+1
                call mgtarg(i,names)
                call chtofp(names,dh)
            else if(names(1:5).eq.'-DMIN')then
                i=i+1
                call mgtarg(i,names)
                call chtofp(names,DMIN)
            else if(names(1:5).eq.'-DMAX')then
                i=i+1
                call mgtarg(i,names)
                call chtofp(names,DMAX)
            else if(names(1:5).eq.'-FMIN')then
                i=i+1
                call mgtarg(i,names)
                call chtofp(names,fitmin)
            else if(names(1:2).eq.'-O')then
                i=i+1
                call mgtarg(i,names)
                dfile = names
            else if(names(1:4).eq.'-PAT')then
                i=i+1
                call mgtarg(i,pathdr)
            else if(names(1:4).eq.'-PMN')then
                i=i+1
                call mgtarg(i,names)
                call chtofp(names,permin)
            else if(names(1:4).eq.'-PMX')then
                i=i+1
                call mgtarg(i,names)
                call chtofp(names,permax)
            else if(names(1:3).eq.'-N1')then
                jnorm = 1
            else if(names(1:3).eq.'-N2')then
                jnorm = 2
            else if(names(1:3).eq.'-CR' .or. names(1:3).eq.'-cr')then
                isds = 3
            else if(names(1:4).eq.'-WT')then
                dowt = .true.
            else if(names(1:2).eq.'-?')then
                call usage(' ')
            else if(names(1:2).eq.'-h')then
                call usage(' ')
            endif
        goto 11
   13   continue
c-----
c       tests
c----

        if(permin.gt.permax)then
            tmp = permin
            permin = permax
            permax = tmp
        endif
c-----
c       special control for CRACK source
c-----
        if(isds.eq.3)then
            rakmn = 90
            rakmx = 90
            drak = 10
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
        if(ostr.ne.' ') write(LER,*)ostr(1:ls)
        write(LER,*)'srgrid96 -DMN dipmn -DMX dipmx -DD ddip'
        write(LER,*)'         -SMN stkmn -SMX stkmx -DS dstk'
        write(LER,*)'         -RMN rakmn -RMX rakmx -DR drak'
        write(LER,*)'         -HMN hmn   -HMX hmx   -DH dh'
        write(LER,*)'         -DMIN dist_min -DMAX dist_max '
        write(LER,*)'         -FMIN fit_min  -O obs -PATH path'
        write(LER,*)'         -PMN permn -PMX permx '
        write(LER,*)'         -N1 -N2 -CRACK -WT'
        write(LER,*)
     1  ' -DMN dipmn  (default   30)  Minimum dip'
        write(LER,*)
     1  ' -DMX dipmx  (default   90)  Maximum dip'
        write(LER,*)
     1  ' -DD  ddip   (default   15)  Dip increment' 
        write(LER,*)
     1  ' -RMN dipmn  (default -180)  Minimum rake'
        write(LER,*)
     1  ' -RMX dipmx  (default  180)  Maximum rake'
        write(LER,*)
     1  ' -DR  ddip   (default   15)  Rake increment' 
        write(LER,*)
     1  ' -SMN stkmn  (default    0)  Minimum strike'
        write(LER,*)
     1  ' -SMX stkmx  (default  350)  Maximum strike'
        write(LER,*)
     1  ' -DS  dstk   (default   10)  Strike increment' 
        write(LER,*)
     1  ' -HMN dpmin  (default    1)  Minimum strike'
        write(LER,*)
     1  ' -HMX dpmax  (default   30)  Maximum strike'
        write(LER,*)
     1  ' -DH  ddepth (default    2)  Strike increment' 
        write(LER,*)
     1  ' -DMIN dmin  (default 0 km     )',
     2  '  minimum for distance sieve'
        write(LER,*)
     1  ' -DMAX dmax  (default 100000 km)',
     2  '  maximum for distance sieve'
        write(LER,*)
     1  ' -FMAX fitmin (default 0.5) Goodness of fit output',
     2  ' threshold'
        write(LER,*)
     1  ' -O observed_data     File with observations '
        write(LER,*)
     1  ' -PATH path   (default ./ ) path to slegn96.der sregn96.der'
        write(LER,*)
     1  ' -PMN permn  (default  0)  Minimum period to use'
        write(LER,*)
     1  ' -PMX permx  (default  2000)  Maximum period to use'
        write(LER,*)
     1  ' -N1                        Obs*exp(gamma r) vs Theo'
        write(LER,*)
     1  ' -N2          (default    ) Obs vs Theo*exp(-gamma r)'
        write(LER,*)
     1  ' -CRACK       (default no ) use crack model: DIP is dip of ',
     2  '                     crack,  STK is  dip direction of crack '
        write(LER,*)
     1  ' -WT          (default none) weight by sqrt(period)  '
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
        else if(isds.eq.1 .or. isds.eq.3)then
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
c       if source or receiver is in fluid, 
c          there is no far-field SH motion
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

        subroutine getobs(dfile, modarr, perarr, amparr, azarr,rr,
     1      modarl, perarl, amparl, azarl,rl,
     1      ndatl, ndatr, dmin, dmax, ierr ,permin, permax)
c-----
c       read observed data
c----
        implicit none
        character dfile*80
        integer ndatl, ndatr
        integer ierr
        integer NOBS
        parameter (NOBS=50000)
        real perarr(NOBS), amparr(NOBS), azarr(NOBS),rr(NOBS)
        integer modarr(NOBS)
        real perarl(NOBS), amparl(NOBS), azarl(NOBS),rl(NOBS)
        integer modarl(NOBS)
        real dmin, dmax
        real permin, permax

        integer ls
        integer lgstr
        logical ext
        character instr*170
        integer mode, lsep, lnobl, j, lorr
        real per, u, du, r, a, amp
        character ic*1
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
        ndatr = 0
        ndatl = 0
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
                    ls = j -1
                endif
C get station and component
   10       continue
            read(instr(lnobl:ls),*)mode,per,u,du,r,a,amp
c-----
c           internally use mode = 1 for fundamental
c-----
c-----
c           only use the data if it is for the correct
c           period and the correct mode
c-----
C     1     ilorr,lorr,mode,kmode,per,u,du,r,a,amp
            mode = mode + 1
            if(r.ge.dmin.and.r.le.dmax .and.per.ge.permin
     1          .and. per.le.permax)then
                if(lorr.eq.1)then
                    ndatl = ndatl + 1
                    modarl(ndatl) = mode
                    perarl(ndatl) = per
                    azarl(ndatl)  = a
                    amparl(ndatl) = amp
                    rl(ndatl)     = r
                else if(lorr.eq.2)then
                    ndatr = ndatr + 1
                    modarr(ndatr) = mode
                    perarr(ndatr) = per
                    azarr(ndatr)  = a
                    amparr(ndatr) = amp
                    rr(ndatr)     = r
                endif
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

        subroutine cptoxxx(xx,xxx,i,toxxx)
        integer NGRN
        parameter (NGRN=21)
        complex xx(NGRN)         
        integer NOBS
        parameter (NOBS=50000)   
        complex xxx(NGRN,NOBS)
        integer i, j
        logical toxxx

        if(toxxx)then
            do 1000 j=1,NGRN
                xxx(j,i) = xx(j)
 1000       continue
        else
            do 2000 j=1,NGRN
                xx(j) = xxx(j,i) 
 2000       continue
        endif
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
