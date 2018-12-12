c----------------------------------------------------------------------c
c                                                                      c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                 c
c      VOLUME III                                                      c
c                                                                      c
c      PROGRAM: SCOMB96                                                c
c                                                                      c
c      COPYRIGHT 2010                                                  c
c      R. B. Herrmann
c      Department of Earth and Atmospheric Sciences                    c
c      Saint Louis University                                          c
c      221 North Grand Boulevard                                       c
c      St. Louis, Missouri 63103                                       c
c      U. S. A.                                                        c
c                                                                      c
c----------------------------------------------------------------------c
c      Changes
c      20 MAR 2012 - correct type in call to getegsh getsh
c      09 SEP 2012 - correct last line of bsort and dbsort from n=n-m to nn=n-m
c             Thanks to ruoshan at ustc
c-----
        integer LIN, LOT, LER
        parameter (LIN=5,LOT=6,LER=0)

        parameter (NL=1000)
        common/ctrl/ xlog,ylog,tmin,tmax,cmin,cmax
        character*4 xlog,ylog

        character ofile*20, dfile*20
        logical ext
c-----
c       command line arguments
c-----
        logical dofreq, doredo
c-----
c       machine dependent initialization
c-----
        call mchdep()
c-----
c       get control parameters from command line
c-----
        call gcmdln(ilorr, xmin, xmax, ymin, ymax, factor, 
     1      dofreq,doredo)
        if(ilorr.lt.1)call usage('specify -L or -R')
        if(ilorr.eq.1)then
            dfile = 'sdisp96.lov'
            ofile = 'tsdisp96.lov'
        else if(ilorr.eq.2)then
            dfile = 'sdisp96.ray'
            ofile = 'tsdisp96.ray'
        endif
        inquire(file=dfile,exist=ext)
        if(.not.ext)then
            call usage('Dispersion file '//dfile//
     1      ' does not exist')
        endif
        open(1,file=dfile,status='unknown',
     1      access='sequential',form='unformatted')
        rewind 1
        open(2,file=ofile,status='unknown',
     1      access='sequential',form='unformatted')
        rewind 2
c-----
c       do the processing
c-----
        call resrch(ilorr,xmin,xmax,ymin,ymax,dofreq,factor, doredo)
c-----
c       close open files
c-----
        close (1)
        close (2)
        end

        subroutine resrch(ilorr,xmin,xmax,ymin,ymax,dofreq,
     1      factor, doredo)
        parameter (LIN=5,LOT=6,LER=0)

        integer*4 ilorr
        real*4 xmin, xmax, ymin, ymax, factor
        logical dofreq, doredo

        parameter (NL=1000)
        common/modl/d(NL), a(NL), b(NL), rho(NL), qa(NL), qb(NL)
        common/pari/ mmax,mode
        common/water/iwat(NL)
c-----
c       model characterization
c-----
        common/modspec/allfluid
        logical allfluid

        integer*4 ifunc,kmode

        real*8 cp(2000)
        real*8 t0, f0

        character mname*80
        integer ipar(10)
        real*4 fpar(10)
c-----
c       begin basic sequence  to read the dispersion file
c-----
c-----
c       get the model
c-----
        rewind 1
        call gtsmdl(1,nmax,d,a,b,rho,qa,qb,nper,
     1      mname,ipar,fpar)
        mmax = nmax
c-----
c       check for completely fluid problem. Also reset the
c       upper bound for ymax so that it neve exceeds the 
c       halfspace S velocity
c-----
        allfluid = .true.
        do 1200 i=1,mmax
            if(b(i).le.1.0e-4*a(i))then
                iwat(i) = 1
            else
                iwat(i) = 0
            endif
            if(b(i).gt.0.0)then
                allfluid = .false.
            endif
 1200   continue
        if(ymax.gt.b(mmax))then
                ymax = b(mmax)
        endif
        call mdsrch()
c-----
c       save the model
c-----
        rewind 2
        call ptsmdl(2,nmax,d,a,b,rho,qa,qb,nper,
     1      mname,ipar,fpar)
c-----
c       process the data sequence
c-----
        kfirst = 0
  100   continue
            call gtshed(1,ifunc,kmode,t0,ierr)
            if(ierr.eq.200)then
                write(LOT,*) 'No ifunc=-1. Data end at perd=',t0
                ifunc=-1
            else if(ierr.eq.1001)then
                go to 1001
            endif
            if(ifunc.le.0) go to 400
            f0 = 1.0/t0
            if(dofreq)then
                tp = 1.0/t0
            else
                tp = t0
            endif
            if(kmode.gt.0) then
c-----
c       read dispersion values at period t0
c-----
            call gtsval(1,cp,kmode,ierr)    
                if(ierr.eq.200)then
                    go to 1000
                else if(ierr.eq.1001)then
                    go to 1001
                endif
            else
c       we have no roots lets try to get them all
c-----
c       output dispersion values at period t0
c-----
            endif
            if(tp.ge.dble(xmin) 
     1          .and. tp.le.dble(xmax))then
                call refine(f0,ymin,ymax,kmode,cp,factor,ilorr,
     1              doredo)
            endif
            call ptshed(2,ifunc,kmode,t0)
            if(kmode.gt.0) call ptsval(2,cp,kmode)  

        go to 100
 400    continue
            call ptshed(2,ifunc,kmode,t0)
        return
 1000   continue
        if(kfirst.eq.0) then
            kfirst=1
            write(LOT,*) 'Data end, not normal',
     1      ' termination at perd=',t0
        endif
        return
 1001   continue
        write(LOT,*) ' '
        write(LOT,*)
     *  'surface85 data file read error at  perd=',t0
        write(LOT,*) ' '
        return
        end

        subroutine refine(f0,ymin,ymax,kmode,cp,factor,ilorr,
     1              doredo)
c-----
c       f0  R*8 current frequency
c       ymin    R*4 lower phase velocity limit
c       ymax    R*4 upper phase velocity limit
c       kmode   I*4 number of known dispersion points
c       cp  R*8 array of phase velocity dispersion
c       factor  R*4 factor for making search more sense
c       iolrr   I*4 1 Love, 2 Rayleigh
c       doredo  L   .true. ignore previous zeros in search region
c-----  
        real*8 f0, cp(2000)
        real*4 ymin, ymax, factor
        integer*4 kmode, ilorr
        logical doredo
        

        integer LIN, LOT, LER
        parameter (LIN=5,LOT=6,LER=0)

        real*8 eroot(2000)
        real*8 wavmin, wavmax, dwav, omega
        real*8 dcmax, dcmin, dcsum
        real*8 wvn
        
c-----
c       first get wavenumber limits of search region
c-----
        omega =  6.283185307179586d+00 * f0
        wavmax = omega/dble(ymin)
        wavmin = omega/dble(ymax)
        kmodld = kmode
c-----
c       convert phase velocities to wave number
c-----
        if(kmode.gt.0)then
c-----
c       given the wavmax and wavmin limits of the search region
c       purge these entries if the search is to be redone from
c       scratch -I flag
c-----
        write(LOT,*) ' '
        write(LOT,*) '============================================='
            if(doredo)then
                j = 0
                do 100 i=1,kmode
                    wvn = omega/cp(i)
                    if(wvn .lt.wavmin 
     1                  .or. wvn .gt.wavmax)then
                        j = j + 1
                        eroot(j) = wvn
                    else
                        write(LOT,*)
     1                     'Ignoring root:',omega/wvn
                    endif
  100           continue
                ndrop = kmode - j
                write(LOT,*)'Initial roots=',kmode, ' dropping',
     1          ndrop,' between c=[',ymin, ',' ,ymax, ']'
                kmode = j
            else
                do 101 i=1,kmode
                    eroot(i) = omega/cp(i)
  101           continue
            endif
c-----
c           sort wavenumbers in order of increasing wavenumber
c-----
            n = kmode
            call dbsort(n,eroot,+1)
            lmode=kmode
            if(lmode.gt.1) then
                dcmax=-1.0d+30
                dcmin=1.0d+30
                dcsum=0.0d+00
                do 2600 i=1,lmode-1
                    dc=cp(i+1)-cp(i)
                    if(dcmax.le.dc) dcmax=dc
                    if(dcmin.ge.dc) dcmin=dc
                    dcsum=dcsum+dc
 2600           continue
                dcsum=dcsum/float(lmode-1)
        write(LOT,*) 'dk value  at the lowest period=',
     1      1.0/f0
        write(LOT,*) ' between c(1)=',cp(1),' and'
        write(LOT,*) '         c(',lmode,')=',  cp(lmode),' :'
        write(LOT,*) ' mean=',dcsum
        write(LOT,*) ' max =',dcmax
        write(LOT,*) ' min =',dcmin
                dwav = dcmin/factor
        write(LOT,*) ' dk  =',dwav
        write(LOT,*) ' '
            else
                dwav = min(abs(eroot(1)-wavmin),
     1              abs(eroot(1)-wavmax))/factor
            endif
        else
            dwav = (wavmax - wavmin)/factor
        endif
c-----
c       given the wavmax and wavmin limits, search through the
c       eroot array to determine the entries within these limits
c       using the fact that the array is sorted in increasing order
c-----  
        if(kmode .gt. 0)then
            klw = 0
            kup = 0
            do 2000 i=1,kmode
                if(eroot(i).gt.wavmin .and.
     1              eroot(i).lt.wavmax)then
                    if(klw .eq. 0)klw = i
                    kup = i
                endif
 2000       continue
        else
            kup = 0
            klw = 0
        endif
c-----
c       we are now ready to begin processing
c-----
        if(klw.eq.0 .and. kup.eq.0)then
            call search(wavmin,wavmax,eroot,kmode,
     1          dwav,klw,kup,omega,ilorr)
        else 
            call search(wavmin,eroot(klw)-dwav,eroot,kmode,
     1          dwav,klw,kup,omega,ilorr)
            if(kup.gt.klw)then
            do 3000 i=klw,kup-1
            call search(eroot(i)+dwav,eroot(i+1)-dwav,
     1          eroot,kmode,dwav,klw,kup,omega,ilorr)
 3000       continue
            endif
            call search(eroot(kup)+dwav,wavmax,eroot,kmode,
     1          dwav,klw,kup,omega,ilorr)
        endif
        n = kmode
        do 4000 i=1,n
            cp(i) = omega/eroot(i)
 4000   continue
        write(LOT,*)'OLD = ',kmodld,' NEW = ',kmode ,' at T0=',
     1      1.0/f0
        call dbsort(n,cp,1)
        return
        end

        subroutine dbsort(n,x,isign)
c-----
c       do bubble sort.
c       isign=+1  increase order   
c            =-1  decreasing.
c-----
          real*8 x(1), x0
c
        m=0
        do 50 i=2,n
            ii=i-m-1
            do 40 j=ii,1,-1
                x0=x(j+1)-x(j)
                if(isign.le.0) x0=-x0
                if(x0 .lt. 0.0d+00)then
                    x0=x(j)
                    x(j)=x(j+1)
                    x(j+1)=x0
                else if(x0.eq.0.0d+00)then
                    m=m+1
                    do 30 k=j,n-m
                            x(k)=x(k+1)
   30               continue
                    go to 50
                else if(x0.gt.0.0d+00)then
                    go to 50
                endif
   40       continue
   50   continue
        nn = n - m
        return
        end


        subroutine gcmdln( ilorr, xmin, xmax, cmin, cmax, 
     2      factor, dofreq,doredo)
        integer*4 ilorr
        real*4 xmin, xmax, cmin, cmax, factor
        logical dofreq,doredo   

        character names*80
c-----
c       defaults
c-----
        ilorr = 0
        xmin = -1.0
        xmax = -1.0
        cmin = -1.0
        cmax = -1.0
        dofreq = .true.
        doredo = .false.
        factor = 10

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
            else if(names(1:5).eq.'-CMIN')then
                i = i + 1
                call mgtarg(i,names)
                read(names,'(bn,f20.0)')cmin
            else if(names(1:5).eq.'-CMAX')then
                i = i + 1
                call mgtarg(i,names)
                read(names,'(bn,f20.0)')cmax
            else if(names(1:4).eq.'-FAC')then
                i = i + 1
                call mgtarg(i,names)
                read(names,'(bn,f20.0)')xx
                if(xx.lt.5.0)then
                    factor = 10.0
                else
                    factor = xx
                endif
            else if(names(1:2).eq.'-I')then
                doredo = .true.
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
        write(LER,*)'scomb96 : ',str
        write(LER,*)'USAGE:',
     1  'scomb96 ',
     1  '[-L -R]  [-FREQ -PER]',
     1  '-XMIN xmin -XMAX xmax  -CMIN cmin -CMAX cmax ',
     1  '-FAC factor -TXT -? -h'

        write(LER,*)
     1  '-L                        Love waves'
        write(LER,*)
     1  '-R                        Rayleigh waves'
        write(LER,*)
     1  '     Note one of -L or -R is required'
        write(LER,*)
     1  '-FREQ      (default true) X limits are frequency'
        write(LER,*)
     1  '-PER       (default false)X limits are period'
        write(LER,*)
     1  '-XMIN xmin (default    )  minimum value of X-Axis'
        write(LER,*)
     1  '-XMAX xmax (default    )  maximum value of X-Axis'
        write(LER,*)
     1  '-CMIN cmin (default    )  minimum value of phase velocity'
        write(LER,*)
     1  '-CMAX cmax (default    )  maximum value of phase velocity'
        write(LER,*)
     1  ' (these limits define the search region)'
        write(LER,*)
     1  '-I         (default false) ignore previous zeros in ',
     1          'search region'
        write(LER,*)
     1  '-FAC factor (default 10) search factor'
        write(LER,*)
     1  '-?         (default false) online help'
        write(LER,*)
     1  '-h         (default false) online help'
        stop 
        end

        subroutine search(wavmin,wavmax,eroot,kmode,
     1          dwav,klw,kup,omega,ifunc)
c*****  This is the most important subroutine in this program.
c       search pole on wavenumber-frequency domain.
c       Each pole is followed by jumping method.
c       At first two periods, poles are located very precisely.
c
c       However, there doesn't exist a perfect scheme for 
c       such a complicated problem. The program 'scomb96' 
c       is designed to repair this unavoidable fault.
c
c       [eroot]  store the wavenumbers for angular freq [omega0]
c       [cphase] store the phase velocies for [omega0+domega]
c       and they will be replaced by the i'th mode's 
c       wavenumber (at omega) and i'th mode's phase velocity (at omega0)
c-----
        parameter (LIN=5, LOT=6,NL=1000,NL2=NL+NL)
        implicit double precision (a-h,o-z)
        real*4 vts(NL2,2)
        integer*4 nmx,i,j,k,ii,jj,kk
        real*8 wvm(200)
        real*8 eroot(2000)
        common/modl/ d,a,b,rho,qa,qb
        real*4 d(NL),a(NL),b(NL),rho(NL),qa(NL),qb(NL)
        common/pari/ mmax,mode
        common/pard/ twopi,displ,dispr
        common/vels/ mvts(2),vts
c-----
c       this routine finds the roots of period equation using
c       regular halving method to initiate the first two
c       frequencies, then followed by jumping method.
c-----
        if(wavmin .gt. wavmax)return
C       write(6,*)'Search between [',wavmin, ',' ,wavmax, ']',
C     1     ' with nroot = ',kmode ,' and f=',omega/6.2831853,
C     2     ' klw,kup=',klw,kup
        epi  = 1.d-10
        wvmx = wavmax
        wvmn = wavmin
        wvmn1= wavmin/1.01
        disp = 5.0
        nmx = (wavmax - wavmin)/dwav + 1
        mode = 2000
        if(nmx.gt.2*mode)nmx = 2*mode
c-----
c       nmx evenly devides the real k axis between Kmin and Kmax.
c-----
        ndisp=dabs(disp)
        if(ndisp.eq.0) ndisp=1
        nmx=nmx*disp
        rangeK = wvmx-wvmn 
        dk = rangeK/nmx
        do 80 i=1,mvts(ifunc)
            wvm(i)=omega/dble(vts(i,ifunc))
   80   continue
c-----
c       This is the most important routine
c       We know that the phase velocity dispersion is bounded by
c       the maximum shear wave velocity and the
c       minimum of P wave velocity (water layer, shear wave velocity,
c           or surface Rayleigh wave velocity
c
c       The number of modes increases with frequency in the manner
c       of Love waves, e.g.,
c           omega H sqrt ( 1. - B1*B1/B2*B2 ) = N p
c           (this is just the SH higher mode cutoff for a
c           single layer over a halfspace)
c
c       We also know that the higher modes are very dense at fixed
c       frequency, just above a layer velocity. The array vts
c       is a list of all possible layer velocities in the 
c       phase velocity limits between cmin and bmax.
c
c       We initially estimate the density of roots and bracket
c       zeros in the period equation, with further refinement
c       by interval halving (slow but always works) or linear
c       interpolation. Linear interpolation has problems since the
c       behavior of the period equation is not linear near a root.
c
c       In an extreme case we may have in the region of a root
c                   ^
c                   |
c       DEL         X
c        |
c        |
c        |
c        -----------|---|---|---------
c                       x   x
c
c                   c1  c3  c4
c       linear interpolation between c1 and c4 would not converge fast
c       since the next estimate would be near c4
c
c       linear interpolation between c3 and c1 would also not work.
c       However, if c3 and c4 are not identical, a root estimate
c       may be able to be made by linear interpolation on c3 and c4
c       if the result lies between c1 and c4. If this is the case,
c       we use a modified interval halving.
c
c       The other aspect concerns estimating the roots at the next
c       frequency. Assuming that the phase velocity dispersion
c       always decreases with frequency, 
c-----


c*****  The regular halving method.
c-----
c       At the very first two periods or the places jumping method
c       fails, we find the poles using the regular halving method.
c       nmx is chosen for a 40 km crustal model. for shallower
c       thickness a proportionately smaller nmx can be used
c       search for roots of period equation
c-----
        nroot = 0
        c2 = wvmx
        del2 = dltar(c2,omega,ifunc)
        lyr=1
        jj=1
        do 500 i=1,nmx
            jj=jj-1
            if(jj.ne.0) go to 500
            c10 = wvmx-float(i)*dk
            if(i.eq.nmx) c10=wvmn+0.01*dk
            jj = 1
            kk = 1
            if(c10.gt.wvm(lyr)) go to 90
c-----
c           kk and jj represent a denser searching when phase velocity
c           = S wave velocity. Their values can be changed as kk=3*lyr
c           jj=8.
c-----
            if(vts(lyr,ifunc).ne.0.0) 
     1          ra=vts(lyr+1,ifunc)/vts(lyr,ifunc)
            if(ra.gt.3.0) ra=3.01
            if(ra.lt.1.0) ra=1.01
            nra=ra
            kk = 10.0*ra
            kk = kk*ndisp
            lyr = lyr+1
            jj = 4*nra+ndisp
   90       continue
            dk0 = dk/float(kk)
            do 400 j=1,jj
            do 400 k=1,kk
                if(nroot.eq.mode) go to 510
                jk = kk*(j-1)+k
                c1 = c10+dk-dk0*float(jk)
                if(c1.lt.wvmn) go to 510
                del1 = dltar(c1,omega,ifunc)
                if(dsign(1.0d+00,del1)*
     1              dsign(1.0d+00,del2).ge.0) go to 350
                nroot = nroot + 1
                c4 = c2
                del4 = del2
c-----
c       perform interval halving, result is guaranteed to be within
c       be 2**(-15) abs(c4-c1) of the true zero
c-----
                do 200 ii=1,15
                    c3 = 0.5*(c1+c4)
                    del3 = dltar(c3,omega,ifunc)
                    if(dsign(1.0d+00,del1)*
     1                  dsign(1.0d+00,del3).ge.0)then
                        del1 = del3
                        c1 = c3
                    else 
                        del4 = del3
                        c4 = c3
                    endif
                    if(dabs(c4-c1).lt.epi*c1) go to 250
  200           continue
  250           continue
                c3 = 0.5*(c1+c4)
                kmode = kmode + 1
                eroot(kmode) = c3
        write(LOT,*)'New root:',omega/c3
  350           c2 = c1
                del2 = del1
  400       continue
  500   continue
  510   continue
        return
        end

        subroutine mdsrch()
c-----
c       Examine the velocity model. Order the velocity model in terms
c       of increasing velocity. This is done twice, once for S only
c       and then for S and P combined. The objective is to determine
c       regions where a denser search in phase velocity should
c       be made to avoid missing modes
c-----
        parameter (LIN=5, LOT=6, NL=1000, NL2=NL+NL)
        implicit double precision (a-h,o-z)
        real*4 d(NL),a(NL),b(NL),rho(NL),qa(NL),qb(NL)
        real*4 vth(NL2),vtp(NL2)
        real*4 vts(NL2,2)
        equivalence(vth(1),vts(1,1)),(vtp(1),vts(1,2))
        common/modl/ d,a,b,rho,qa,qb
        common/water/iwat(NL)
        common/pari/ nmax,mode
        common/vels/ mvts(2),vts
c-----
c       store the layer velocity to be used to 
c       initiate denser pole searching.
c-----
            m=0
            do 310 i=1,nmax
                if(b(i).gt.0.0) then
                    m=m+1
                    vth(m)=b(i)
                    iwat(i) = 0
                else
                    iwat(i) = 1
                endif
  310       continue
            call bsort(m,vth,1)
            m=m+1
            vth(m)=1.0e+30
            mvts(1)=m
            m=0
            do 320 i=1,nmax
                if(a(i).gt.0.0) then
                    m=m+1
                    vtp(m)=a(i)
                endif
                if(b(i).gt.0.0) then
                    m=m+1
                    vtp(m)=b(i)
                endif
  320       continue
            call bsort(m,vtp,1)
            m=m+1
            vtp(m)=1.0e+30
            mvts(2)=m
        return
        end

        function dltar(wvno,omega,kk)
c-----
c       control the way to P-SV or SH.
c-----
        implicit double precision (a-h,o-z)
        if(kk.eq.1)then
c-----
c           love wave period equation
c-----
            dltar = dltar1(wvno,omega)
C            WRITE(6,*)'T,C,FL:',6.2831853/omega,omega/wvno,dltar
C            DO I=1,10
C            dltar = dltar1(wvno-0.05*I*wvno,omega)
C            WRITE(6,*)'T,C,FL:',6.2831853/omega,
C     1         omega/(wvno-0.05*I*wvno),dltar
C            ENDDO

        else if(kk.eq.2)then
c-----
c           rayleigh wave period equation
c-----
            dltar = dltar4(wvno,omega)
C        WRITE(6,*)6.2831853/omega, omega/wvno,dltar
        endif
        end

        function dltar1(wvno,omega)
c-----
c       find SH dispersion values.
c-----
        parameter (NL=1000)
        implicit double precision (a-h,o-z)
        real*4 d(NL),a(NL),b(NL),rho(NL),qa(NL),qb(NL)
        common/modl/ d,a,b,rho,qa,qb
        common/pari/ mmax,modew
        common/pard/ twopi,displ,dispr
        common/water/iwat(NL)

        complex*16 esh(2,2), einvsh(2,2)
        double precision hl(2,2)
        logical lshimag
        double precision rsh
        complex*16 e10, e20

        real mu

c-----
c       Haskell-Thompson love wave formulation from halfspace
c       to surface.
c-----
c-----
c       get the eigenfunctions for the halfspace
c-----

        mu = rho(mmax)*b(mmax)*b(mmax)
        call gtegsh(mmax,wvno, omega, rsh, lshimag)
        call gtesh(esh,einvsh,rsh,wvno,mu,lshimag)

        e1=dreal(einvsh(1,1))
        e2=dreal(einvsh(1,2))

        mmm1 = mmax - 1
        do 600 m=mmm1,1,-1
            if(iwat(m).eq.0)then
               call gtegsh(m,wvno, omega, rsh, lshimag)
               mu = rho(m)*b(m)*b(m)
               call gtesh(esh,einvsh,rsh,wvno,mu,lshimag)
               call varsh(d(m),rsh,lshimag,cossh,rsinsh,sinshr,exm)
               call hskl(cossh,rsinsh,sinshr,mu,iwat(m),hl,exm,exe)
                e10=e1*hl(1,1)+e2*hl(2,1)
                e20=e1*hl(1,2)+e2*hl(2,2)
                xnor=dabs(dreal(e10))
                ynor=dabs(dreal(e20))
                if(ynor.gt.xnor) xnor=ynor
                if(xnor.lt.1.d-40) xnor=1.0d+00
                e1=dreal(e10)/xnor
                e2=dreal(e20)/xnor
            endif
  600   continue
        dltar1=e1
        return
        end
c
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
        function dltar4(wvno,omga)
c-----
c       find P-SV dispersion values.
c-----
        parameter (NL=1000)
        implicit double precision (a-h,o-z)
        common/modl/ d,a,b,rho,qa,qb
        real*4 d(NL),a(NL),b(NL),rho(NL),qa(NL),qb(NL)


        real*8 e(5),ee(5),ca(5,5)
        common/pari/ mmax,mode
        integer mmax,mode
        common/water/iwat(NL)
        integer iwat

        complex*16 gbr(2,5)
        COMPLEX*16 zA,zB,zC,zD,zE,zF
        COMPLEX*16 rp, rsv
        complex*16 p,q
        COMPLEX*16 NP, NSV
        COMPLEX*16 X11, X21, X31, X41
        COMPLEX*16 X12, X22, X32, X42
        real*8 cosp , cossv
        real*8 rsinp, rsinsv
        real*8 sinpr, sinsvr
        real*8 pex, svex
c-----
c       model characterization
c-----
        common/modspec/allfluid
        logical allfluid

CRBH        integer i,j
c-----
c       set up starting values for bottom halfspace
c-----
        wvno2=wvno*wvno
        om2 = omga*omga
        call evalg(0,mmax,mmax-1,gbr,1,
     1      wvno,omga,om2,wvno2)

c-----
c       CORRECTED 14 OCT 2006
c-----
            e(1) = dreal(gbr(1,1) )
            e(2) = dreal(gbr(1,2) )
            e(3) = dreal(gbr(1,3) )
            e(4) = dreal(gbr(1,4) )
            e(5) = dreal(gbr(1,5) )
C        WRITE(6,'(i4,5e11.3)')mmax,e
c-----
c-----
c       matrix multiplication from bottom layer upward
c-----
        mmm1 = mmax-1
        do 500 m = mmm1,1,-1
                xka = omga/a(m)
                if(b(m).gt.0.0)then
                  xkb = omga/b(m)
                else
                  xkb = 0.0
                endif
                rp =CDSQRT(dcmplx(wvno2-xka*xka,0.0d+00))
                rsv=CDSQRT(dcmplx(wvno2-xkb*xkb,0.0d+00))
                p = rp  * dble(d(m))
                q = rsv * dble(d(m))
                call varsv(p,q,rp, rsv, 
     1              cosp, cossv, rsinp, rsinsv, 
     1              sinpr, sinsvr, pex,svex,iwat(m))
                call dnka(CA,cosp,rsinp,sinpr,cossv,rsinsv,sinsvr,
     1              Rho(m),b(m),iwat(m),pex,pex+svex,wvno,wvno2,om2)

            nmat = 5
            do 200 i=1,nmat
                cr=0.0d+00
                do 100 j=1,nmat
                    cr=cr+e(j)*ca(j,i)
  100           continue
                ee(i)=cr
  200       continue
            call normc(ee,exa,nmat)
            do 300 i = 1,nmat
                e(i)=ee(i)
  300       continue
C        WRITE(6,'(i4,5e11.3)')m,e
  500   continue
        dltar4 = e(1)
        return
        end


        subroutine dnka(CA,cosp,rsinp,sinpr,cossv,rsinsv,sinsvr,
     1              rho,b,iwat,ex,exa,wvno,wvno2,om2)
        implicit none
        real*8 ca(5,5), wvno, wvno2, om2
        real*4 rho,b
        real*8 cosp,rsinp,sinpr,cossv,rsinsv,sinsvr
        integer iwat

        real rho2
        real*8 gam
        real*8 ex, exa
        real*8 a0
        real*8 cpcq,cpy,cpz,cqw,cqx,xy,xz,wy,wz

        real*8 gam2,gamm1,gamm2,a0c,xz2,wy2,temp
        real*8 dfac
        real*8 cqww2, cqxw2, g1wy2, gxz2, g2wy2, g2xz2
        real*8 gg1, a0cgg1
        integer i, j
c-----
c        A11     A12     A13    -A13     A15     A16
c        A21     A22     A23    -A23     A25     A15
c        A31     A32     A33    1-A33   -A23    -A13
c       -A31    -A32    1-A33    A33     A23     A13
c        A51     A52    -A32     A32     A22     A12
c        A61     A51    -A31     A31     A21     A11
c-----
c       this will be multipled on the left by the G matrix
c
c       [ G11   G12 G13 -G13    G15 G16 ]
c
c-----
c       or on the right by
c
c       [ H11   H21 H31 -H31    H51 H61  ] ^T
c-----
c       the number of multiplications can be reduced from 36 to 25 
c       if we define a new matrices
c       related to the original matrices by
c-----
c         A11     A12     A13         A15     A16
c         A21     A22     A23         A25     A15
c        2 A31   2 A32   2 A33 -1   -2 A23  -2 A13
c         A51     A52    -A32         A22     A12
c         A61     A51    -A31         A21     A11
c-----
c
c       [ G11   G12  G13    G15 G16  ]
c       [ H11   H21 2 H31   H51 H61  ] ^T
c
c-----
c       this means that some of the original definitions of the 
c       Aij elements must be changed for the
c       definition of the modified 5x5 compount A matrix
c
c       old 6x6                 new 5x5
c       A11 = 1 - A33 + A22     1 - (1/2)(new A33 + 1) + new A2)
c       A53 = -A32              A43 = - (1/2) new A32
c       A63 = -A31              A53 = - (1/2) new A31
c-----
c       To recover the needed elements, 
c          we note that the old G14 = -old G14 = new G13
c-----
c-----
        if(iwat.eq.1)then
c-----
c       fluid layer
c-----
            do 100 j=1,5
                do 101 i=1,5
                    ca(i,j) = 0.0d+00
  101           continue
  100       continue
            if(ex.gt.35.0d+00)then
                dfac = 0.0d+00
            else
                dfac = dexp(-ex)
            endif
            ca(3,3) = dfac
            ca(1,1) = cosp
            ca(5,5) = cosp
            ca(1,2) = - rsinp/(rho*om2)
            ca(2,1) = - rho*sinpr*om2
            ca(2,2) = cosp
            ca(4,4) = cosp
            ca(4,5) = ca(1,2)
            ca(5,4) = ca(2,1)
        else
            if ( exa.lt. 60.0)then
                a0=dexp(-exa)
            else
                a0 = 0.0d+00
            endif
            cpcq = cosp * cossv
            cpy  = cosp*sinsvr
            cpz  = cosp*rsinsv
            cqw  = cossv*sinpr
            cqx  = cossv*rsinp
            xy   = rsinp*sinsvr
            xz   = rsinp*rsinsv
            wy   = sinpr*sinsvr  
            wz   = sinpr*rsinsv
c-----
c       elastic layer
c-----
            rho2= rho*rho
            gam = 2.0*b*b*wvno2/om2
            gam2  = gam*gam
            gamm1 = gam-1.
            gamm2 = gamm1*gamm1
            cqww2 = cqw * wvno2
            cqxw2 = cqx / wvno2
            gg1 = gam*gamm1
            a0c  = dcmplx(2.0d+00,0.0d+00)*
     1          (dcmplx(a0,0.0d+00)-cpcq)
            xz2  = xz/wvno2
            gxz2 = gam*xz2
            g2xz2 = gam2 * xz2
            a0cgg1 = a0c*(gam+gamm1)
            wy2  = wy*wvno2
            g2wy2 = gamm2 * wy2
            g1wy2 = gamm1 * wy2

c-----
c       OK by symmetry
c----
            temp = a0c*gg1 + g2xz2 + g2wy2
            ca(3,3) = a0 + temp + temp
            ca(1,1) = cpcq-temp
            ca(1,2) = (-cqx + wvno2*cpy)/(rho*om2)
            temp = dcmplx(0.5d+00,0.0d+00)*a0cgg1 + gxz2 + g1wy2
            ca(1,3) = wvno*temp/(rho*om2)

            ca(1,4) = (-cqww2+cpz)/(rho*om2)
            temp = wvno2*(a0c + wy2) + xz
            ca(1,5) = -temp/(rho2*om2*om2)

            ca(2,1) = (-gamm2*cqw + gam2*cpz/wvno2)*rho*om2
            ca(2,2) = cpcq
            ca(2,3) = (gamm1*cqww2 - gam*cpz)/wvno
            ca(2,4) = -wz
            ca(2,5)=ca(1,4)


            temp =dcmplx(0.5d+00,0.0d+00)*a0cgg1*gg1 
     1          + gam2*gxz2 + gamm2*g1wy2
            ca(3,1) = -dcmplx(2.0d+00,0.0d+00)*temp*rho*om2/wvno
            ca(3,2) = -wvno*(gam*cqxw2 - gamm1*cpy)*
     1          dcmplx(2.0d+00,0.0d+00)

            ca(3,4)=-2.0d+00*ca(2,3)
            ca(3,5)=-2.0d+00*ca(1,3)

            ca(4,1) = (-gam2*cqxw2 + gamm2*cpy)*rho*om2
            ca(4,2) = -xy
            ca(4,3)= -ca(3,2)/2.0d+00
            ca(4,4)=ca(2,2)
            ca(4,5)=ca(1,2)

            temp = gamm2*(a0c*gam2 + g2wy2) + gam2*g2xz2
            ca(5,1) = -rho2*om2*om2*temp/wvno2
            ca(5,2)=ca(4,1)
            ca(5,3)=-ca(3,1)/2.0d+00
            ca(5,4)=ca(2,1)
            ca(5,5)=ca(1,1)
        endif
        return
        end

        subroutine evalg(jbdry,m,m1,gbr,in,
     1      wvno,om,om2,wvno2)
        complex*16 gbr(2,5), gbl(2,2)
        parameter(NL=1000)
        common/modl/ d,a,b,rho,qa,qb
        real*4 d(NL),a(NL),b(NL),rho(NL),qa(NL),qb(NL)
        common/shwave/bsh(NL), rhosh(NL), qbsh(NL)
        real bsh, rhosh, qbsh
        common/modlly/mmax
        common/damp/alpha,ieqex
        real*8 om,xka,xkb,gam,gamm1, wvno, wvno2, om2
        complex*16 ra,rb
        common/modspec/allfluid
        logical allfluid
        complex*16 CDSQRT

        complex*16 e(4,4), einv(4,4)
        
c-----
c       set up halfspace conditions
c-----
        xka = om/a(m)
        if(b(m).gt.0.0)then
          xkb = om/b(m)
          iwat = 0
        else
          xkb = 0.0
          iwat = 1
        endif
        ra=CDSQRT(dcmplx(wvno2-xka*xka,0.0d+00))
        rb=CDSQRT(dcmplx(wvno2-xkb*xkb,0.0d+00))
        gam = dble(b(m))*wvno/om
        gam = 2.0d+00 * (gam * gam)
        gamm1 = gam - dcmplx(1.0d+00,0.0d+00)
C        WRITE(6,*)'WVNO :=',wvno
C        WRITE(6,*)'OM2  :=',om2
C        WRITE(6,*)'RA   :=',ra
C        WRITE(6,*)'RB   :=',rb
C        WRITE(6,*)'gam  :=',gam
C        WRITE(6,*)'gamm1:=',gamm1
c-----
c       set up halfspace boundary conditions
c
c       jbdry   = -1  RIGID
c           =  0  ELASTIC
c           = +1  FREE SURFACE
c
c-----
        if(jbdry.lt.0)then
c-----
c       RIGID - check properties of layer above
c-----
            if(b(m) .gt. 0.0)then
c-----
c               ELASTIC ABOVE - RIGID
c-----
                gbr(in,1) = dcmplx(1.0d+00,0.0d+00)
                gbr(in,2) = dcmplx(0.0d+00,0.0d+00)
                gbr(in,3) = dcmplx(0.0d+00,0.0d+00)
                gbr(in,4) = dcmplx(0.0d+00,0.0d+00)
                gbr(in,5) = dcmplx(0.0d+00,0.0d+00)
                gbl(in,1) = dcmplx(1.0d+00,0.0d+00)
                gbl(in,2) = dcmplx(0.0d+00,0.0d+00)
            else
c-----
c               FLUID ABOVE - RIGID
c-----
                gbr(in,1) = dcmplx(0.0d+00,0.0d+00)
                gbr(in,2) = dcmplx(0.0d+00,0.0d+00)
                gbr(in,3) = dcmplx(0.0d+00,0.0d+00)
                gbr(in,4) = dcmplx(0.0d+00,0.0d+00)
                gbr(in,5) = dcmplx(0.0d+00,0.0d+00)
                if(allfluid)then
                    gbr(in,1) = dcmplx(1.0d+00,0.0d+00)
                else
                    gbr(in,4) = dcmplx(1.0d+00,0.0d+00)
                endif
c-----
c               (pseudo SH)
c-----
                gbl(in,1) = dcmplx(1.0d+00,0.0d+00)
                gbl(in,2) = dcmplx(0.0d+00,0.0d+00)
            endif
        else if(jbdry.eq.0)then
c-----
c       HALFSPACE
c-----
            if(iwat.eq.0)then
c-----
c               ELASTIC HALFSPACE
c-----
c       multiply G of Herrmann 2001 by - rho^2 om^4 k^2 ra rb
c       should have no effect since it is in both numerator and
c       denominator -- however will not give the correct potential
c       coefficient -- so rethink?
c-----
            E(1,1)= wvno
            E(1,2)= rb
            E(1,3)= wvno
            E(1,4)= -rb

            E(2,1)= ra
            E(2,2)= wvno
            E(2,3)= -ra
            E(2,4)= wvno

            E(3,1)= rho(m)*om2*gamm1
            E(3,2)= rho(m)*om2*gam*rb/wvno
            E(3,3)= rho(m)*om2*gamm1
            E(3,4)= -rho(m)*om2*gam*rb/wvno

            E(4,1)= rho(m)*om2*gam*ra/wvno
            E(4,2)= rho(m)*om2*gamm1
            E(4,3)= -rho(m)*om2*gam*ra/wvno
            E(4,4)= rho(m)*om2*gamm1

            EINV(1,1)= 0.5*gam/wvno
            EINV(1,2)= -0.5*gamm1/ra
            EINV(1,3)= -0.5*1.0/(rho(m)*om2)
            EINV(1,4)= 0.5*wvno/(rho(m)*om2*ra)

            EINV(2,1)= -0.5*gamm1/rb
            EINV(2,2)= 0.5*gam/wvno
            EINV(2,3)= 0.5*wvno/(rho(m)*om2*rb)
            EINV(2,4)= -0.5*1.0/(rho(m)*om2)

            EINV(3,1)= 0.5*gam/wvno
            EINV(3,2)=  0.5*gamm1/ra
            EINV(3,3)= -0.5*1.0/(rho(m)*om2)
            EINV(3,4)= -0.5*wvno/(rho(m)*om2*ra)

            EINV(4,1)= 0.5*gamm1/rb
            EINV(4,2)= 0.5*gam/wvno
            EINV(4,3)= -0.5*wvno/(rho(m)*om2*rb)
            EINV(4,4)= -0.5*1.0/(rho(m)*om2)

C          do j=1,4
C              do i = 1,4
C                write(6,*)'E   (',I ,',' , J, ')=',E(i,j)
C             enddo
C          enddo
C          do j=1,4
C              do i = 1,4
C                write(6,*)'EINV(',I ,',' , J, ')=',EINV(i,j)
C             enddo
C          enddo
C          do i=1,4
C              do j = 1,4
C                zsum = dcmplx(0.0d+00,0.0d+00)
C                do k=1,4
C                zsum = zsum + E(i,k)*einv(k,j)
C                enddo
C                write(6,*)'E INV(',I ,',' , J, ')=',ZSUM
C             enddo
C          enddo



                gbr(in,1)=dble(rho(m)*rho(m))*om2*om2*
     1              (-gam*gam*ra*rb+wvno2*gamm1*gamm1)
                gbr(in,2)=-dble(rho(m))*(wvno2*ra)*om2
                gbr(in,3)=-dble(rho(m))*(-gam*ra*rb+wvno2*gamm1)
     1              *om2*wvno
                gbr(in,4)=dble(rho(m))*(wvno2*rb)*om2
                gbr(in,5)=wvno2*(wvno2-ra*rb)
        gbr(in,1) = 0.25*gbr(in,1)/(-rho(m)*rho(m)*om2*om2*wvno2*ra*rb)
        gbr(in,2) = 0.25*gbr(in,2)/(-rho(m)*rho(m)*om2*om2*wvno2*ra*rb)
        gbr(in,3) = 0.25*gbr(in,3)/(-rho(m)*rho(m)*om2*om2*wvno2*ra*rb)
        gbr(in,4) = 0.25*gbr(in,4)/(-rho(m)*rho(m)*om2*om2*wvno2*ra*rb)
        gbr(in,5) = 0.25*gbr(in,5)/(-rho(m)*rho(m)*om2*om2*wvno2*ra*rb)

                gbl(in,1) =  dble(rhosh(m))*(dble(bsh(m)*bsh(m))
     1                  )*rb
                gbl(in,2) =  dcmplx(1.0d+00,0.0d+00)
            else if(iwat.eq.1)then
c-----
c               FLUID HALFSPACE
c-----
                if(allfluid)then
                    gbr(in,1) = dble(0.5) / ra
                    gbr(in,2) = dcmplx(0.5d+00,0.0d+00)/
     1                  (-dble(rho(m))*om2)
                    gbr(in,3) = dcmplx(0.0d+00,0.0d+00)
                    gbr(in,4) = dcmplx(0.0d+00,0.0d+00)
                    gbr(in,5) = dcmplx(0.0d+00,0.0d+00)
                else
                    gbr(in,1) = dcmplx(0.0d+00,0.0d+00)
                    gbr(in,2) = dcmplx(0.0d+00,0.0d+00)
                    gbr(in,3) = dcmplx(0.0d+00,0.0d+00)
                    gbr(in,4) = dble(0.5*rho(m)*om2) / ra
                    gbr(in,5) = dcmplx(-0.5d+00,0.0d+00)
                endif
                gbl(in,1) = dcmplx(0.0d+00,0.0d+00)
                gbl(in,2) = dcmplx(1.0d+00,0.0d+00)
            endif
        else if(jbdry.eq.1)then
c-----
c       FREE - check properties of layer above
c-----
            if(b(m) .gt. 0.0)then
                gbr(in,1) = dcmplx(0.0d+00,0.0d+00)
                gbr(in,2) = dcmplx(0.0d+00,0.0d+00)
                gbr(in,3) = dcmplx(0.0d+00,0.0d+00)
                gbr(in,4) = dcmplx(0.0d+00,0.0d+00)
                gbr(in,5) = dcmplx(1.0d+00,0.0d+00)
                gbl(in,1) = dcmplx(0.0d+00,0.0d+00)
                gbl(in,2) = dcmplx(1.0d+00,0.0d+00)
                
            else
                gbr(in,1) = dcmplx(0.0d+00,0.0d+00)
                gbr(in,2) = dcmplx(0.0d+00,0.0d+00)
                gbr(in,3) = dcmplx(0.0d+00,0.0d+00)
                gbr(in,4) = dcmplx(0.0d+00,0.0d+00)
                gbr(in,5) = dcmplx(0.0d+00,0.0d+00)
                if(allfluid)then
                    gbr(in,2) = dcmplx(1.0d+00,0.0d+00)
                else
                    gbr(in,5) = dcmplx(1.0d+00,0.0d+00)
                endif
                gbl(in,1) = dcmplx(0.0d+00,0.0d+00)
                gbl(in,2) = dcmplx(1.0d+00,0.0d+00)
            endif
        endif
        return
        end

        subroutine varsv(p,q, rp, rsv, 
     1      cosp, cosq, rsinp, rsinq, 
     1      sinpr, sinqr, pex,svex,iwat)
c-----
c       p = rp  * h
c       q = rsv * h
c       rp  vertical wave number for P
c       rsv vertical wave number for SV
c       cosp=cosh(p)  rsinp =rp *sinh(p)  sinpr = sinh(p)/rp
c       cosq=cosh(q)  rsinsv=rsv*sinh(p)  sinpq = sinh(p)/rsv
c-----
        implicit none
        COMPLEX*16 p, q
        COMPLEX*16 rp, rsv
        real*8 cosp, cosq
        real*8 rsinp, rsinq
        real*8 sinpr, sinqr
        REAL *8 pex,svex
        integer iwat

        REAL*8 pr, pi, qr, qi
        COMPLEX*16 epp, epm, eqp, eqm
        COMPLEX*16 sinp, sinq

        REAL*8 PFAC, SVFAC
        
        pex  = 0.0d+00
        svex = 0.0d+00
        pr = dreal(p)
        pi = dimag(p)
        qr = dreal(q)
        qi = dimag(q)
        pex   = pr
        if(iwat.eq.1)then
c-----
c       fluid layer
c-----
            epp = dcmplx(dcos(pi), dsin(pi))/2.0
            epm = dconjg(epp)
            if(pr.lt.15.) then
                pfac=dexp(-2.*pr)
            else
                pfac  = 0.0d+00
            endif
            cosp = epp + pfac*epm
            sinp = epp - pfac*epm
            rsinp = rp *sinp
            sinpr = sinp/rp
            cosq  = 1.0d+00
            rsinq = 0.0d+00
            sinqr = 0.0d+00
        else
c-----
c       elastic layer
c-----
            svex = qr
            epp = dcmplx(dcos(pi), dsin(pi))/2.0
            epm = dconjg(epp)
            eqp = dcmplx(dcos(qi), dsin(qi))/2.0
            eqm = dconjg(eqp)
            if(pr.lt.15.) then
                pfac=dexp(-2.*pr)
            else
                pfac  = 0.0d+00
            endif
            cosp = epp + pfac*epm
            sinp = epp - pfac*epm
            rsinp = rp *sinp
            sinpr = sinp/rp

            if(qr.lt.15.) then
                svfac=dexp(-2.*qr)
            else
                svfac  = 0.0d+00
            endif
            cosq = eqp + svfac*eqm
            sinq = eqp - svfac*eqm
            rsinq = rsv*sinq
            sinqr = sinq/rsv

        endif
        return
        end

        subroutine normc(ee,ex,nmat)
c       This routine is an important step to control over- or
c       underflow.
c       The Haskell or Dunkin vectors are normalized before
c       the layer matrix stacking.
c       Note that some precision will be lost during normalization.
c
        implicit double precision (a-h,o-z)
        dimension ee(5)
        ex = 0.0d+00
        t1 = 0.0d+00
        do 10 i = 1,nmat
          if(dabs(ee(i)).gt.t1) t1 = dabs(ee(i))
   10 continue
        if(t1.lt.1.d-40) t1=1.d+00
        do 20 i =1,nmat
          t2=ee(i)
          t2=t2/t1
          ee(i)=t2
   20 continue
c-----
c       store the normalization factor in exponential form.
c-----
        ex=dlog(t1)
        return
        end

c-----
c       notes:
c       For the surface wave the propagator is always positive
c       However the E and E^-1 matrices can be complex
c-----

        subroutine varsh(h,rsh,lshimag,cossh,rsinsh,sinshr,ex)
        implicit none
        real h
        double precision rsh
        logical lshimag
        double precision cossh, rsinsh, sinshr
        double precision ex
        
        double precision q, fac, sinq
        q = rsh*h
        ex =  0.0
        if(lshimag)then
            if(rsh.gt.0.0)then
                 cossh = dcos(q)
                 sinq = sin(q)
                 sinshr = sinq/rsh
                 rsinsh = - rsh*sinq
            else
                 cossh  = 1.0d+00
                 sinshr = dble(h)
                 rsinsh = 0.0
            endif
        else
            if(q.lt.16.0d+00)then
                fac = dexp(-2.0d+00*q)
            else
                fac = 0.0d+00
            endif
            cossh = (1.0d+00 + fac) * 0.5d+00
            sinq   = (1.0d+00 - fac) * 0.5d+00
            sinshr = sinq/rsh
            rsinsh = sinq*rsh
        endif
        return
        end
        

        subroutine hskl(cossh,rsinsh,sinshr,mu,iwat,hl,exm,exe)
c-----
c       True cosh( rd b )    = exp(exm) cossh
c       True sinh( rd b )/rb = exp(exm) sinshr
c       True rb sinh( rd b ) = exp(exm) rsinsh
c-----
        implicit none
        double precision cossh, rsinsh, sinshr
        double precision hl(2,2),exm,exe
        real mu
        integer iwat

        if(iwat.eq.0)then   
            hl(1,1) = cossh
            hl(2,1) = mu*rsinsh
            hl(1,2) = sinshr/(mu)
            hl(2,2) = cossh
            exe = exe + exm
        else
            hl(1,1) = 1.0
            hl(1,2) = 0.0
            hl(2,1) = 0.0
            hl(2,2) = 1.0
        endif
        return
        end 

        subroutine gtesh(esh,einvsh,rsh,wvno,mu,lshimag)
        complex*16 esh(2,2), einvsh(2,2)
        double precision wvno, rsh
        real mu
        logical lshimag
c-----
c       get E and E^-1 for Quasi SH
c       however the E^-1 returned has is really
c       (2k mu rsh) E^-1 to avoid dividing by zero when rsh = 0
c-----
        esh(1,1) =   wvno
        esh(1,2) =   wvno
        if(lshimag)then
            esh(2,1) =   wvno * mu * dcmplx(0.0d+00,rsh)
            esh(2,2) = - wvno * mu * dcmplx(0.0d+00,rsh)
            einvsh(1,1) =   mu * dcmplx(0.0d+00,rsh)
            einvsh(2,1) =   mu * dcmplx(0.0d+00,rsh)
        else
            esh(2,1) =   wvno * mu * dcmplx(rsh,0.0d+00)
            esh(2,2) = - wvno * mu * dcmplx(rsh,0.0d+00)
            einvsh(1,1) =   mu * dcmplx(rsh,0.0d+00)
            einvsh(2,1) =   mu * dcmplx(rsh,0.0d+00)
        endif
        einvsh(1,2) =   1.0
        einvsh(2,2) =  -1.0
        return
        end

        subroutine gtegsh(m,wvno,omega,rsh,lshimag)
        implicit none

        integer m
        double precision wvno, omega, rsh
        logical lshimag

        integer NL
        parameter (NL=1000)
        common/modl/ d,a,b,rho,qa,qb
        real*4 d(NL),a(NL),b(NL),rho(NL),qa(NL),qb(NL)
        common/water/iwat(NL)
        integer iwat


c-----
c       internal variables
c-----
        double precision wvno2, omega2
        wvno2 = wvno*wvno
        omega2 = omega*omega
        if(iwat(m).eq.1)then
c-----
c           check for water layer
c-----
            rsh = 0.0
            lshimag = .false.
        else
            rsh = wvno2 - omega2/(b(m)*b(m))
            if(rsh .lt. 0.0)then
                rsh = dsqrt(dabs(rsh))
                lshimag = .true.
            else
                rsh = dsqrt(dabs(rsh))
                lshimag = .false.
            endif
        endif
        
        return
        end
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
        subroutine bsort(nn,x,isign)
c       do bubble sort.
c       isign=+1  increase   =-1 decrease.
c
          real*4 x(nn)
c
        m=0
        n = nn
        do 50 i=2,n
          ii=i-m-1
          do 40 j=ii,1,-1
            x0=x(j+1)-x(j)
              if(isign.le.0) x0=-x0
            if(abs(x0).lt.1.e-7) go to 20
C              if(x0) 10,20,50
            if(x0 .lt. 0.0)then
                go to 10
            else if(x0 .eq. 0.0)then
                go to 20
            else
                go to 50
            endif
10        continue
            x0=x(j)
            x(j)=x(j+1)
            x(j+1)=x0
            go to 40
20        continue
            m=m+1
            do 30 k=j,n-m
              x(k)=x(k+1)
30        continue
            go to 50
40      continue
50    continue
        nn=n-m
        return
        end
