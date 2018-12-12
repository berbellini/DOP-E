c----------------------------------------------------------------------c
c                                                                      c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                 c
c      VOLUME III                                                      c
c                                                                      c
c      PROGRAM: SCOMB96                                                c
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
c Changes
c      20 MAR 2012 - changed arguments in gtegsh for consistence of type
c      09 SEP 2012 - correct last line of bsort and dbsort from n=n-m to nn=n-m
c             Thanks to ruoshan at ustc
c-----
        integer LIN, LOT, LER
        parameter (LIN=5,LOT=6,LER=0)

        parameter (NL=1000)
        common/timod/d(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      TRho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL), 
     3      frefp(NL), frefs(NL)
        real*4 d,TA,TC,TN,TL,TF,TRho,qa,
     1       qb,etap,etas,frefp,frefs
c-----
c       model characterization
c-----
        common/modspec/allfluid
        logical allfluid

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
            dfile = 'tdisp96.lov'
            ofile = 'ttdisp96.lov'
        else if(ilorr.eq.2)then
            dfile = 'tdisp96.ray'
            ofile = 'ttdisp96.ray'
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
        common/timod/d(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      TRho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL), 
     3      frefp(NL), frefs(NL)
        real*4 d,TA,TC,TN,TL,TF,TRho,qa,
     1       qb,etap,etas,frefp,frefs
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
        integer ipar(20)
        real*4 fpar(20)
c-----
c       begin basic sequence  to read the dispersion file
c-----
c-----
c       get the model
c-----
        rewind 1
        call gtsmdt(1,mmax,d,ta,tc,tf,tl,tn,trho,qai,qbi,nper,
     2              mname,ipar,fpar)

c-----
c       check for  fluid only problem and also reset
c       ymax so that it never exceeds the halfspace S velocity
c-----
        allfluid = .true.
        do 1200 i=1,mmax
            if(TN(i).le.1.0e-4*TA(i))then
                iwat(i) = 1
            else
                iwat(i) = 0
            endif
            if(TN(i).gt.0.0)then
                allfluid = .false.
            endif
 1200   continue
        call mdsrch()
c-----
c       save the model
c-----
        rewind 2
        call ptsmdt(2,mmax,d,ta,tc,tf,tl,tn,Trho,qa,qb,nper,
     1      mname,ipar,fpar)

c-----
c       process the data sequence
c-----
        kfirst = 0
  100   continue
            call gtshed(1,ifunc,kmode,t0,ierr)
            if(ifunc.eq.1)then
                    ymax = sqrt(TN( mmax)/Trho(mmax) )
            else 
                    if(iwat(mmax).eq.1)then
                            ymax = sqrt( TA(mmax)/Trho(mmax))
                    else
                            ymax = sqrt(TL( mmax)/Trho(mmax) )
                    endif
            endif
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
        common/timod/d(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      TRho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL), 
     3      frefp(NL), frefs(NL)
        real*4 d,TA,TC,TN,TL,TF,TRho,qa,
     1       qb,etap,etas,frefp,frefs
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
        real*4 vth(NL2),vtp(NL2)
        real*4 vts(NL2,2)
        equivalence(vth(1),vts(1,1)),(vtp(1),vts(1,2))
        common/timod/d(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      TRho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL), 
     3      frefp(NL), frefs(NL)
        real*4 d,TA,TC,TN,TL,TF,TRho,qa,
     1       qb,etap,etas,frefp,frefs
        common/water/iwat(NL)
        common/pari/ nmax,mode
        common/vels/ mvts(2),vts
c-----
c       store the layer velocity to be used to 
c       initiate denser pole searching.
c-----
            m=0
            do 310 i=1,mmax
                if(sqrt(TL(i)/Trho(i)).gt.0.000) then
                    m=m+1
                    vth(m)=sqrt(TL(i)/Trho(i))
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
            do 320 i=1,mmax
                if(sqrt(TA(i)/TRho(i)).gt.0.0) then
                    m=m+1
                    vtp(m)=sqrt(TA(i)/TRho(i))
                endif
                if(sqrt(TL(i)/Trho(i)).gt.0.000) then
                    m=m+1
                    vtp(m)=sqrt(TL(i)/Trho(i))
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

        common/timod/d(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      TRho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL), 
     3      frefp(NL), frefs(NL)
        real*4 d,TA,TC,TN,TL,TF,TRho,qa,
     1       qb,etap,etas,frefp,frefs
        common/pari/ mmax,modew
        common/pard/ twopi,displ,dispr
        common/water/iwat(NL)

        complex*16 esh(2,2), einvsh(2,2)
        double precision hl(2,2)
        logical lshimag
        double precision rsh
        complex*16 e10, e20

c-----
c       Haskell-Thompson love wave formulation from halfspace
c       to surface.
c-----

c-----
c       get the eigenfunctions for the halfspace
c-----

        call gtegsh(mmax,wvno, omega, rsh, lshimag)
        call gtesh(esh,einvsh,rsh,wvno,dble(TL(mmax)),lshimag,
     1        .false.,iwat(mmax))

        e1=dreal(einvsh(1,1))
        e2=dreal(einvsh(1,2))
        mmm1 = mmax - 1
        exe = 0.0
        do 600 m=mmm1,1,-1
            if(iwat(m).eq.0)then
               call gtegsh(m,wvno, omega, rsh, lshimag)
               call gtesh(esh,einvsh,rsh,wvno,dble(TL(m)),lshimag,
     1         .false.,iwat(m))
               call varsh(dble(d(m)),rsh,lshimag,
     1            cossh,rsinsh,sinshr,exm)
               call hskl(cossh,rsinsh,sinshr,dble(TL(m)),
     1            iwat(m),hl,exm,exe)
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

c-----
c       notes:
c       For the surface wave the propagator is always positive
c       However the E and E^-1 matrices can be complex
c-----

        subroutine varsh(h,rsh,lshimag,cossh,rsinsh,sinshr,ex)
        implicit none
        double precision h
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
            ex = q
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
        

        subroutine hskl(cossh,rsinsh,sinshr,TL,iwat,hl,exm,exe)
c-----
c       True cosh( rd b )    = exp(exm) cossh
c       True sinh( rd b )/rb = exp(exm) sinshr
c       True rb sinh( rd b ) = exp(exm) rsinsh
c-----
        implicit none
        double precision cossh, rsinsh, sinshr
        double precision hl(2,2),exm,exe
        double precision TL
        integer iwat

        if(iwat.eq.0)then   
            hl(1,1) = cossh
            hl(2,1) = TL*rsinsh
            hl(1,2) = sinshr/(TL)
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

        subroutine gtesh(esh,einvsh,rsh,wvno,L,lshimag,ltrueinv,iwat)
        complex*16 esh(2,2), einvsh(2,2)
        double precision wvno, rsh
        double precision L
        logical lshimag,ltrueinv
        integer iwat
        double complex rb
c-----
c       get E and E^-1 for Quasi SH
c       if ltrueinv == .true. then return E^-1
c       else return the normalized avoid divide by zero
c       (2k L rsh) E^-1 to avoid dividing by zero when rsh = 0
c-----
        if(iwat.eq.1)then
c-----
c           for a water layer use an identity matrix. The
c           results have validity only for the extreme case of
c           a model that is fluid-solid-fluid, e.g., ocean,mantle,
c           outer core.  The propagator matrix technique cannot 
c           handle anything more complicated
c-----
             esh(1,1) = dcmplx(1.0d+00, 0.0d+00)
             esh(1,2) = dcmplx(0.0d+00, 0.0d+00)
             esh(2,1) = dcmplx(0.0d+00, 0.0d+00)
             esh(2,2) = dcmplx(1.0d+00, 0.0d+00)
             einvsh(1,1) = dcmplx(1.0d+00, 0.0d+00)
             einvsh(1,2) = dcmplx(0.0d+00, 0.0d+00)
             einvsh(2,1) = dcmplx(0.0d+00, 0.0d+00)
             einvsh(2,2) = dcmplx(1.0d+00, 0.0d+00)
        else
             if(lshimag)then
                 rb =  dcmplx(0.0d+00,rsh)
             else
                 rb = dcmplx(rsh,0.0d+00)
             endif
             esh(1,1) =   wvno
             esh(1,2) =   wvno
             esh(2,1) =   wvno * L * rb
             esh(2,2) = - wvno * L * rb
             einvsh(1,1) =   L * rb
             einvsh(2,1) =   L * rb
             einvsh(1,2) =   1.0
             einvsh(2,2) =  -1.0
             if(ltrueinv)then
                 einvsh(1,1) = einvsh(1,1)/(2.*wvno*L*rb)
                 einvsh(1,2) = einvsh(1,2)/(2.*wvno*L*rb)
                 einvsh(2,1) = einvsh(2,1)/(2.*wvno*L*rb)
                 einvsh(2,2) = einvsh(2,2)/(2.*wvno*L*rb)
             endif
        endif
        return
        end

        subroutine gtegsh(m,wvno,omega,rsh,lshimag)
        implicit none

        integer m
        double precision wvno, omega, rsh
        logical lshimag

        integer NL
        parameter (NL=1000)
        common/timod/d(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      TRho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL), 
     3      frefp(NL), frefs(NL)
        real*4 d,TA,TC,TN,TL,TF,TRho,qa,
     1       qb,etap,etas,frefp,frefs
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
            rsh = TN(m)*wvno2/(TL(m)) - TRho(m)*omega2/(TL(m))
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

        subroutine lmult(d11,d12,d21,d22,hl,iwat,exel,exb,icomp)
        implicit none
c-----
c       multiply SH matrix by a row vector on left
c-----
        complex*16 d11,d12,d21,d22,hl(2,2),e1,e2
        integer iwat
        real*8 exel, exb
        logical icomp
c-----
c       fluid layer do nothing, just return, 
c           equivalent to multiplying by
c       identity matrix
c-----
        if(iwat.eq.0)then
c-----
c       elastic layer
c-----
            e1=d11
            e2=d12
c-----
c           a11 = cosql
c           a12 = yl
c           a21 = zl
c           a22 = cosql
c-----
            d11=e1*hl(1,1) + e2*hl(2,1)
            d12=e1*hl(1,2) + e2*hl(2,2)
            exel = exel + exb
            if(icomp)then
                e1=d21
                e2=d22
                d21=e1*hl(1,1) + e2*hl(2,1)
                d22=e1*hl(1,2) + e2*hl(2,2)
            endif
        endif
        return
        end

        function dltar4(wvno,omga)
c-----
c       find P-SV dispersion values.
c-----
        parameter (NL=1000)
        implicit double precision (a-h,o-z)
        common/timod/d(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      TRho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL),
     3      frefp(NL), frefs(NL)
        real*4 d,TA,TC,TN,TL,TF,TRho,qa,
     1       qb,etap,etas,frefp,frefs
        common/pari/ mmax,mode
        integer mmax,mode
        common/water/iwat(NL)
        integer iwat
        real*8 e(5),ee(5),ca(5,5)

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

CRBH	integer i,j
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
           call  gettiegn(Za,Zb,Zc,Zd,Ze,Zf,om2,wvno2,rp, rsv,
     1      x11,x21,x31,x41,x12,x22,x32,x42,NP, NSV,m,omga,wvno,
     2      dcmplx(1.0d+00,0.0d+00),dcmplx(1.0d+00,0.0d+00))
                p = rp  * dble(d(m))
                q = rsv * dble(d(m))
                call varsv(p,q,rp, rsv, 
     1              cosp, cossv, rsinp, rsinsv, 
     1              sinpr, sinsvr, pex,svex,iwat(m))
                call dnka(CA,cosp,rsinp,sinpr,cossv,rsinsv,sinsvr,
     1              NP,NSV,
     1              x11,x21,x31,x41,x12,x22,x32,x42,
     1              TRho(m),iwat(m),pex+svex,om2)

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
  500   continue
        dltar4 = e(1)
        return
        end

        subroutine dnka(CA,cosp,rsinp,sinpr,cossv,rsinsv,sinsvr,NP,NSV,
     1      x11,x21,x31,x41,x12,x22,x32,x42,
     1      TRho,iwat,ex,om2)
        implicit none
        real*8 CA(5,5)
        real*8 cosp , cossv 
        real*8 rsinp, rsinsv
        real*8 sinpr, sinsvr
        COMPLEX*16 NP, NSV
        COMPLEX*16 X11, X21, X31, X41
        COMPLEX*16 X12, X22, X32, X42
        real TRho
        integer iwat
        real*8 ex, dfac
        real*8  om2
        complex*16 zrho

        integer i, j
        complex*16 TCA(6,6)
c-----
c       introduce conciseness to reduce operations
c-----
        COMPLEX*16 NPNP
        COMPLEX*16 NPNSV
        COMPLEX*16 NSVNSV
        COMPLEX*16 CPCSV, CPSVR, CPRSV, SPRCSV, SPRSVR
        COMPLEX*16 SPRRSV, RSPCSV,  RSPSVR, RSPRSV
        COMPLEX*16 X11X11, X11X21, X11X31, X11X41
        COMPLEX*16         X21X21, X21X31, X21X41
        COMPLEX*16                 X31X31, X31X41
        COMPLEX*16                         X41X41
        COMPLEX*16 X12X12, X12X22, X12X32, X12X42
        COMPLEX*16         X22X22, X22X32, X22X42
        COMPLEX*16                 X32X32, X32X42
        COMPLEX*16                         X42X42

        COMPLEX*16 FAC01, FAC02, FAC03, FAC04, FAC05
        COMPLEX*16 FAC06, FAC07, FAC08, FAC09, FAC10
        COMPLEX*16 FAC11, FAC12, FAC13, FAC14, FAC15
        COMPLEX*16 FAC16, FAC17, FAC18, FAC19, FAC20
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
c       the number of multiplications can be reduced from 
c            36 to 25 if we define a new matrices
c       related to the original matrices by
c-----
c         A11     A12     A13         A15     A16
c         A21     A22     A23         A25     A15
c       2 A31   2 A32    2 A33 -1   -2 A23  -2 A13
c         A51     A52    -A32         A22     A12
c         A61     A51    -A31         A21     A11
c-----
c
c       [ G11   G12  G13    G15 G16  ]
c       [ H11   H21 2 H31   H51 H61  ] ^T
c
c-----
c       this means that some of the original definitions of the 
c           Aij elements must be changed for the
c       definition of the modified 5x5 compount A matrix
c
c       old 6x6                 new 5x5
c       A11 = 1 - A33 + A22     1 - (1/2)(new A33 + 1) + new A2)
c       A53 = -A32              A43 = - (1/2) new A32
c       A63 = -A31              A53 = - (1/2) new A31
c-----
c       To recover the needed elements, we note that the 
c           old G14 = -old G14 = new G13
c-----
c-----
        zrho = dcmplx(dble(TRho),0.0d+00)
        if(ex.gt.35.0d+00)then
            dfac = 0.0d+00
        else
            dfac = dexp(-ex)
        endif
        if(iwat.eq.1)then
c-----
c       fluid layer
c-----
            do 100 j=1,5
                do 101 i=1,5
                    ca(i,j) = 0.0d+00
  101           continue
  100       continue
            ca(3,3) = dfac 
            ca(1,1) = cosp 
            ca(5,5) = cosp 
            ca(1,2) = -rsinp/(zrho*om2) 
            ca(2,1) = - zrho*sinpr*om2 
            ca(2,2) = cosp 
            ca(4,4) = cosp 
            ca(4,5) = ca(1,2) 
            ca(5,4) = ca(2,1) 
        else
c-----
c           introduce some repeated expressions to 
c           reduce multiplications
c-----
            NPNP   = dfac/(NP*NP)
            NSVNSV = dfac/(NSV * NSV)
            NPNSV  = 1.0d+00/(NP * NSV)
            CPCSV  = cosp *cossv
            CPSVR  = cosp *sinsvr
            CPRSV  = cosp *rsinsv
            SPRCSV = sinpr*cossv
            SPRSVR = sinpr*sinsvr
            SPRRSV = sinpr*rsinsv
            RSPCSV = rsinp*cossv
            RSPSVR = rsinp*sinsvr
            RSPRSV = rsinp*rsinsv

            X11X11 = x11*x11
            X11X21 = x11*x21    
            X11X31 = x11*x31    
            X11X41 = x11*x41    
            X21X21 = x21*x21
            X21X31 = x21*x31
            X21X41 = x21*x41
            X31X31 = x31*x31
            X31X41 = x31*x41
            X41X41 = x41*x41

            X12X12 = x12*x12
            X12X22 = x12*x22
            X12X32 = x12*x32
            X12X42 = x12*x42
            X22X22 = x22*x22
            X22X32 = x22*x32
            X22X42 = x22*x42
            X32X32 = x32*x32
            X32X42 = x32*x42
            X42X42 = x42*x42

            FAC01= X11X11*X22X32*RSPCSV*(NPNSV)
            FAC02=X11X11*X22X42*RSPRSV*(NPNSV)
            FAC03=X11X11*X32X42*RSPCSV*(NPNSV)
            FAC04=X11X21*X12X32*CPSVR*(NPNSV)
            FAC05=X11X21*X12X42*CPCSV*(NPNSV)
            FAC06=X11X31*X12X22*RSPCSV*(NPNSV)
            FAC07=X11X31*X12X32*RSPSVR*(NPNSV)
            FAC08=X11X31*X12X42*RSPCSV*(NPNSV)
            FAC09=X11X41*X12X22*CPCSV*(NPNSV)
            FAC10=X11X41*X12X32*CPSVR*(NPNSV)
            FAC11=X11X41*X22X32*CPCSV*(NPNSV)
            FAC12=X21X31*X12X12*CPSVR*(NPNSV)
            FAC13=X21X31*X12X42*CPCSV*(NPNSV)
            FAC14=X21X31*X42X42*CPRSV*(NPNSV)
            FAC15=X21X41*X12X12*SPRSVR*(NPNSV)
            FAC16=X21X41*X32X42*SPRCSV*(NPNSV)
            FAC17=X31X41*X12X12*CPSVR*(NPNSV)
            FAC18=X31X41*X22X42*CPRSV*(NPNSV)
            FAC19=X31X41*X32X42*CPCSV*(NPNSV)
            FAC20=X41X41*X22X32*SPRCSV*(NPNSV)        

c-----
c       repeated terms
c       X11X11*X22X32*RSPCSV*(NPNSV)
c       X11X11*X22X42*RSPRSV*(NPNSV)
c       X11X11*X32X42*RSPCSV*(NPNSV)
c       X11X21*X12X32*CPSVR*(NPNSV)
c       X11X21*X12X42*CPCSV*(NPNSV)
c       X11X31*X12X22*RSPCSV*(NPNSV)
c       X11X31*X12X32*RSPSVR*(NPNSV)
c       X11X31*X12X42*RSPCSV*(NPNSV)
c       X11X41*X12X22*CPCSV*(NPNSV)
c       X11X41*X12X32*CPSVR*(NPNSV)
c       X11X41*X22X32*CPCSV*(NPNSV)
c       X21X31*X12X12*CPSVR*(NPNSV)
c       X21X31*X12X42*CPCSV*(NPNSV)
c       X21X31*X42X42*CPRSV*(NPNSV)
c       X21X41*X12X12*SPRSVR*(NPNSV)
c       X21X41*X32X42*SPRCSV*(NPNSV)
c       X31X41*X12X12*CPSVR*(NPNSV)
c       X31X41*X22X42*CPRSV*(NPNSV)
c       X31X41*X32X42*CPCSV*(NPNSV)
c       X41X41*X22X32*SPRCSV*(NPNSV)        
c-----

c-----
c       ALSO NOTE THAT NONE OF THE TCA(?,4) or TCA(4,?) ARE REQUIRED
c       SO DO NOT COMPUTE
c-----
            
c-----
c       elastic layer
c-----
c CA11 12 12 = 11 22 - 12 21
c CA12 12 13 = 11 23 - 13 21
c CA13 12 14 = 11 24 - 14 21 = - 11 13 - 14 21
c CA14 12 23 = 12 23 - 22 13
c CA15 12 24 = 12 24 - 22 14
c CA16 12 34 = 13 24 - 23 14
c CA21 13 12 = 11 32 - 12 31
c CA22 13 13 = 11 33 - 13 31
c CA23 13 14 = 11 34 - 14 31
c CA24 13 23 = 12 33 - 13 32 = 12 22 - 13 32
c CA25 13 24  = 12 34 - 14 32
c CA26 13 34  = 13 34 - 14 33
c CA31 14 12 = 11 42 - 12 41 = - 11 31 - 12 41
c CA32 14 13 = 11 43 - 13 41
c CA33 14 14 = 11 44 - 14 41 = 11 11 - 14 41
c CA34 14 23 = 12 43 - 13 42 = - 12 21 + 13 31
c CA35 14 24 = 12 44 - 14 42 = 12 11 + 14 31
c CA36 14 34 = 13 44 - 14 43 = 13 11 + 14 21  = - CA13
c CA41 23 12 = 21 32 - 22 31 
c CA42 23 13 = 21 33 - 23 31 
c CA43 23 14 = 21 34 - 24 31 = - 21 12 + 13 31
c CA44 23 23 = 22 33 - 23 32 
c CA45 23 24 = 22 34 - 24 32 = - 22 12 + 13 32 
c                            = - 33 12 + 13 32 = - CA24
c CA46 23 34 = 23 34 - 24 33 = - 23 12 + 13 22 = - CA14 
c CA51 24 12 = 21 42 - 22 41 = - 21 31 - 22 41
c CA52 24 13 = 21 43 - 23 41 = - 21 21 - 23 41
c CA53 24 14 = 21 44 - 24 41 = - 11 43 + 13 41 = - CD32
c CA54 24 23 = 22 43 - 23 42 = - 22 21 + 23 31 
c                            = - 33 21 + 23 31 = -CD42
c CA55 24 24 = 22 44 - 24 42 = 33 11 - 13 31 = CD22 
c CA56 24 34 = 23 44 - 24 43 = 23 11 - 13 21 = CD12
c CA61 34 12 = 31 42 - 32 41 = - 31 31 - 32 41
c CA62 34 13 = 31 43 - 33 41 = - 21 31 - 22 41
c CA63 34 14 = 31 44 - 34 41 = - 11 43 + 13 41 = - CD32
c CA64 34 23 = 32 43 - 33 42 = - 22 21 + 23 31 
c                            = - 33 21 + 23 31 = -CD42
c CA65 34 24 = 32 44 - 34 42 = 32 11 - 12 31 = CD21 
c CA66 34 34 = 33 44 - 34 43 = 22 11 - 12 21 = CD11

        
            TCA(1,1) = - X11X21*X31X41*(NPNP ) 
     1            - X12X22*X32X42*(NSVNSV) 
     1            - FAC11
     1            - FAC13
     1            + X11X31*X22X42*RSPRSV*(NPNSV)
     1            + X21X41*X12X32*SPRSVR*(NPNSV)
            TCA(1,2) = - X11X41*X22X22*CPRSV*(NPNSV) 
     1            - X21X21*X12X42*SPRCSV*(NPNSV) 
     1            + X11X21*X22X42*CPRSV*(NPNSV)
     1            + X21X41*X12X22*SPRCSV*(NPNSV)
            TCA(1,3) =   X11X11*X21X41*(NPNP )
     1            + X12X12*X22X42*(NSVNSV)
     1            + FAC09
     1            + FAC05
     1            - FAC15
     1            - FAC02
C           TCA(1,4) = - X11X21*X21X31*(NPNP )
C     1           - X12X22*X22X32*(NSVNSV)
C     1           + X11X31*X22X22*RSPRSV*(NPNSV)
C     1           + X21X21*X12X32*SPRSVR*(NPNSV)
C     1           - X21X31*X12X22*CPCSV*(NPNSV)
C     1           - X11X21*X22X32*CPCSV*(NPNSV)
            TCA(1,5) =  - FAC06
     1             - FAC04
     1             + FAC12
     1             + FAC01
            TCA(1,6) = - X11X11*X21X21*(NPNP )
     1            - X12X12*X22X22*(NSVNSV)
     1            - 2.0*X11X21*X12X22*CPCSV*(NPNSV)
     1            + X21X21*X12X12*SPRSVR  *(NPNSV)
     1            + X11X11*X22X22*RSPRSV  *(NPNSV)
            TCA(2,1) = - X11X41*X32X32*CPSVR*(NPNSV)
     1            - X31X31*X12X42*RSPCSV*(NPNSV)
     1            + X11X31*X32X42*RSPCSV*(NPNSV)
     1            + X31X41*X12X32*CPSVR*(NPNSV)
            TCA(2,2) = - FAC11
     1            - FAC13
     1            + X11X21*X32X42*CPCSV*(NPNSV)
     1            + X31X41*X12X22*CPCSV*(NPNSV)
            TCA(2,3) =   FAC10
     1            + FAC08
     1            - FAC03
     1            - FAC17
C           TCA(2,4) = + X11X31*X22X32*RSPCSV*(NPNSV)
C     1           + X21X31*X12X32*CPSVR*(NPNSV)
C     1           - X11X21*X32X32*CPSVR*(NPNSV)
C     1           - X31X31*X12X22*RSPCSV*(NPNSV)
            TCA(2,5) = - 2.0d+00 * FAC07
     1            + X11X11*X32X32*RSPSVR*(NPNSV)
     1            + X31X31*X12X12*RSPSVR*(NPNSV)
            TCA(2,6) = - FAC04
     1            - FAC06
     1            + FAC01
     1            + FAC12
            TCA(3,1) = - X11X31*X41X41*(NPNP )
     1            - X12X32*X42X42*(NSVNSV)
     1            - X11X41*X32X42*CPCSV*(NPNSV)
     1            - X31X41*X12X42*CPCSV*(NPNSV)
     1            + X11X31*X42X42*RSPRSV*(NPNSV)
     1            + X41X41*X12X32*SPRSVR*(NPNSV)
            TCA(3,2) = - X11X41*X22X42*CPRSV*(NPNSV)
     1            - X21X41*X12X42*SPRCSV*(NPNSV)
     1            + X11X21*X42X42*CPRSV*(NPNSV)
     1            + X41X41*X12X22*SPRCSV*(NPNSV)
            TCA(3,3) =   X11X11*X41X41*(NPNP )
     1            + X12X12*X42X42*(NSVNSV)
     1            + 2.0*X11X41*X12X42*CPCSV*(NPNSV)
     1            - X11X11*X42X42*RSPRSV*(NPNSV)
     1            - X41X41*X12X12*SPRSVR*(NPNSV)
C           TCA(3,4) = - X11X21*X31X41*(NPNP )
C     1           - X32X42*X12X22*(NSVNSV)
C     1           + X11X31*X22X42*RSPRSV*(NPNSV)
C     1           + X21X41*X12X32*SPRSVR*(NPNSV)
C     1           - X11X21*X32X42*CPCSV*(NPNSV)
C     1           - X31X41*X12X22*CPCSV*(NPNSV)
            TCA(3,5) = - FAC08
     1            - FAC10
     1            + FAC03
     1            + FAC17
            TCA(3,6) =  - TCA(1,3)
C           TCA(4,1) =   X21X41*X31X31*(NPNP )
C     1           + X22X42*X32X32*(NSVNSV)
C     1           - X21X41*X32X32*SPRSVR*(NPNSV)
C     1           - X31X31*X22X42*RSPRSV*(NPNSV)
C     1           + X21X31*X32X42*CPCSV  *(NPNSV)
C     1           + X31X41*X22X32*CPCSV  *(NPNSV)
C           TCA(4,2) = - X21X41*X22X32*SPRCSV*(NPNSV)
C     1           - X21X31*X22X42*CPRSV*(NPNSV)
C     1           + X31X41*X22X22*CPRSV*(NPNSV)
C     1           + X21X21*X32X42*SPRCSV*(NPNSV)
C           TCA(4,3) = - X31X41*X11X21*(NPNP )
C     1           - X12X22*X32X42*(NSV*NSV)
C     1           - X31X41*X12X22*CPCSV*(NPNSV)
C     1           - X11X21*X32X42*CPCSV*(NPNSV)
C     1           + X11X31*X22X42*RSPRSV*(NPNSV)
C     1           + X21X41*X12X32*SPRSVR*(NPNSV)
C           TCA(4,4) =   X21X31*X21X31*(NPNP )
C     1           + X22X32*X22X32*(NSVNSV)
C     1           + X21X31*X22X32*CPCSV  *(NPNSV)
C     1           + X21X31*X22X32*CPCSV  *(NPNSV)
C     1           - X21X21*X32X32*SPRSVR*(NPNSV)
C     1           - X31X31*X22X22*RSPRSV*(NPNSV)
C           TCA(4,5) =   TCA(2,3)
C           TCA(4,6) = - TCA(1,4)
            TCA(5,1) = - FAC16
     1            - FAC18
     1            + FAC14
     1            + FAC20
            TCA(5,2) = - 2.0 * X21X41*X22X42*SPRRSV*(NPNSV)
     1            + X21X21*X42X42*SPRRSV*(NPNSV)
     1            + X41X41*X22X22*SPRRSV*(NPNSV)
            TCA(5,3) = - TCA(3,2)
C           TCA(5,4) = - TCA(4,2)
            TCA(5,5) =   TCA(2,2)
            TCA(5,6) =   TCA(1,2)
            TCA(6,1) = - X31X41*X31X41*(NPNP )
     1            - X32X42*X32X42*(NSVNSV)
     1            - 2.0d+00 *FAC19
     1            + X31X31*X42X42*RSPRSV *(NPNSV)
     1            + X41X41*X32X32*SPRSVR *(NPNSV)
            TCA(6,2) = - FAC18
     1            - FAC16
     1            + FAC14
     1            + FAC20
            TCA(6,3) = - TCA(3,1)
C           TCA(6,4) =   TCA(3,1)
            TCA(6,5) =   TCA(2,1)
            TCA(6,6) =   TCA(1,1)
c-----
c       for development only - clean up later
c       define the CA(5,5)
c-----
            CA(1,1) = dreal( TCA(1,1) )
            CA(1,2) = dreal( TCA(1,2) )
            CA(1,3) = dreal( TCA(1,3) )
            CA(1,4) = dreal( TCA(1,5) )
            CA(1,5) = dreal( TCA(1,6) )

            CA(2,1) = dreal( TCA(2,1) )
            CA(2,2) = dreal( TCA(2,2) )
            CA(2,3) = dreal( TCA(2,3) )
            CA(2,4) = dreal( TCA(2,5) )
            CA(2,5) = dreal( TCA(1,5) )

            CA(3,1) =   2.0d+00*dreal( TCA(3,1) )
            CA(3,2) =   2.0d+00*dreal( TCA(3,2) )
c-----
c       beware of normalization
c-----
            CA(3,3) =   2.0d+00*dreal( TCA(3,3)) - 1.0d+00*dfac 
            CA(3,4) =  -2.0d+00*dreal( TCA(2,3) )
            CA(3,5) =  -2.0d+00*dreal( TCA(1,3) )

            CA(4,1) =   dreal( TCA(5,1) )
            CA(4,2) =   dreal( TCA(5,2) )
            CA(4,3) =  -dreal( TCA(3,2) )
            CA(4,4) =   dreal( TCA(2,2) )
            CA(4,5) =   dreal( TCA(1,2) )

            CA(5,1) =   dreal( TCA(6,1) )
            CA(5,2) =   dreal( TCA(5,1) )
            CA(5,3) =  -dreal( TCA(3,1) )
            CA(5,4) =   dreal( TCA(2,1) )
            CA(5,5) =   dreal( TCA(1,1) )
        endif
        return
        end

        subroutine evalg(jbdry,m,m1,gbr,in,
     1      wvno,om,om2,wvno2)
        implicit none
        integer jbdry, m, m1, in
        complex*16 gbr(2,5)
        real*8 wvno,om,wvno2,om2
        integer NL
        parameter(NL=1000)
        common/timod/d(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      TRho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL),
     3      frefp(NL), frefs(NL)
        real*4 d,TA,TC,TN,TL,TF,TRho,qa,
     1       qb,etap,etas,frefp,frefs
        common/modlly/mmax
        integer mmax
        common/modspec/allfluid
        logical allfluid

        complex*16 cg(6)
        complex*16 g(4,4)
        COMPLEX*16 zA,zB,zC,zD,zE,zF
        COMPLEX*16 rp, rsv
        COMPLEX*16 NP, NSV
        COMPLEX*16 X11, X21, X31, X41
        COMPLEX*16 X12, X22, X32, X42
        integer iwat

        complex*16 e(4,4)



C        WRITE(6,*)'allfluid:',allfluid
c-----
c       set up halfspace conditions
c-----
            if(TL(m).eq. 0.0 .or.TN(m).eq.0.0)then
                iwat = 1
            else
                iwat = 0
            endif
c-----
c       HALFSPACE
c-----
            call gettiegn(za,zb,zc,zd,ze,zf,om2,wvno2,rp, rsv, 
     1          x11,x21,x31,x41,x12,x22,x32,x42,NP, NSV,m,om,wvno,
     2          dcmplx(1.0d+00,0.0d+00),dcmplx(1.0d+00,0.0d+00))
            if(iwat.eq.0)then
c-----
c               ELASTIC HALFSPACE
C       WRITE(6,*)'X:'
C       WRITE(6,*)x11,x21,x31,x41
C       WRITE(6,*)x12,x22,x32,x42
C        WRITE(6,*)'RP,RSV:',RP,RSV
C        WRITE(6,*)'NP,NSV,m,om,wvno:',NP,NSV,m,om,wvno
C       WRITE(6,*)'om2,wvno2:',om2,wvno2


                G(1,1) =    x41/(2.0*NP  *rp  )
                G(1,2) = -  x31/(2.0*NP       )
                G(1,3) = -  x21/(2.0*NP  *rp  )
                G(1,4) =    x11/(2.0*NP       )
                G(2,1) =    x42/(2.0*NSV      )
                G(2,2) = -  x32/(2.0*NSV *rsv )
                G(2,3) = -  x22/(2.0*NSV      )
                G(2,4) =    x12/(2.0*NSV *rsv )
                G(3,1) =    x41/(2.0*NP  *rp  )
                G(3,2) =    x31/(2.0*NP       )
                G(3,3) = -  x21/(2.0*NP  *rp  )
                G(3,4) = -  x11/(2.0*NP       )
                G(4,1) = -  x42/(2.0*NSV      )
                G(4,2) = -  x32/(2.0*NSV *rsv )
                G(4,3) =    x22/(2.0*NSV      )
                G(4,4) =    x12/(2.0*NSV *rsv )
          E(1,1)=  X11*rp
          E(2,1)=  X21
          E(3,1)=  X31*rp
          E(4,1)=  X41
          E(1,2)=  X12
          E(2,2)=  X22*rsv
          E(3,2)=  X32
          E(4,2)=  X42*rsv
          E(1,3)=  X11*rp
          E(2,3)= -X21
          E(3,3)=  X31*rp
          E(4,3)= -x41
          E(1,4)= -X12
          E(2,4)=  X22*rsv
          E(3,4)= -x32
          E(4,4)=  X42*rsv
C          do j=1,4
C              do i = 1,4
C                   write(6,*)'E   (',I ,',' , J, ')=',E(i,J)
C             enddo
C          enddo
C          do j=1,4
C              do i = 1,4
C                   write(6,*)'EINV(',I ,',' , J, ')=',G(i,J)
C             enddo
C          enddo
C          do i=1,4
C              do j = 1,4
C                zsum = dcmplx(0.0d+00,0.0d+00)
C                do k=1,4
C                zsum = zsum + E(i,k)*g(k,j)
C                enddo
C                write(6,*)'E INV(',I ,',' , J, ')=',ZSUM
C             enddo
C          enddo

c CG(1) = G 12 12 = 11 22 - 12 21
c CG(2) = G 12 13 = 11 23 - 13 21
c CG(3) = G 12 14 = 11 24 - 14 21
c CG(4) = G 12 23 = 12 23 - 13 22
c CG(5) = G 12 24 = 12 24 - 14 22
c CG(6) = G 12 34 = 13 24 - 14 23
                CG(1) = G(1,1)*G(2,2) - G(1,2)*G(2,1)
                CG(2) = G(1,1)*G(2,3) - G(1,3)*G(2,1)
                CG(3) = G(1,1)*G(2,4) - G(1,4)*G(2,1)
                CG(4) = G(1,2)*G(2,3) - G(1,3)*G(2,2)
                CG(5) = G(1,2)*G(2,4) - G(1,4)*G(2,2)
                CG(6) = G(1,3)*G(2,4) - G(1,4)*G(2,3)
                gbr(in,1) = CG(1)
                gbr(in,2) = CG(2)
                gbr(in,3) = CG(3)
                gbr(in,4) = CG(5)
                gbr(in,5) = CG(6)

            else if(iwat.eq.1)then
c-----
c               FLUID HALFSPACE
c-----
                if(allfluid)then
                    gbr(in,1) = dble(TRho(m))*om2
                    gbr(in,2) = -rp
                    gbr(in,3) = dcmplx(0.0d+00,0.0d+00)
                    gbr(in,4) = dcmplx(0.0d+00,0.0d+00)
                    gbr(in,5) = dcmplx(0.0d+00,0.0d+00)
                else
                    gbr(in,1) = dcmplx(0.0d+00,0.0d+00)
                    gbr(in,2) = dcmplx(0.0d+00,0.0d+00)
                    gbr(in,3) = dcmplx(0.0d+00,0.0d+00)
                    gbr(in,4) = -dble(TRho(m))*om2
                    gbr(in,5) = rp
                endif
            endif
        return
        end

        subroutine gettiegn(a,b,c,d,e,f,omega2,wvno2,rp, rsv, 
     1      x11,x21,x31,x41,x12,x22,x32,x42,NP, NSV,m,omg,wvn,
     2      atna,atnb)
        implicit none
        COMPLEX*16 A,B,C,D,E,F
        real*8 wvno2, omega2
        COMPLEX*16 rp, rsv
c-----
c       get the modified eigen vectors x11 and x31 have 
c           lambda1 (rp ) factored out
c               modified eigen vectors x22 and x42 have 
c           lambda2 (rsv) factored out
c-----
        COMPLEX*16 X11, X21, X31, X41
        COMPLEX*16 X12, X22, X32, X42
c-----
c       norms
c-----
        COMPLEX*16 NP, NSV
        integer m
        real*8 omg, wvn
        complex*16 atna, atnb

        COMPLEX*16 xka2, xkb2

        integer NL
        parameter (NL=1000)
        common/timod/td(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      TRho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL),
     3      frefp(NL), frefs(NL)
        real*4 td,TA,TC,TN,TL,TF,TRho,qa,
     1       qb,etap,etas,frefp,frefs
c-----
c       internal variables
c-----
        COMPLEX*16 L2(2)
        COMPLEX*16 bb, cc
        COMPLEX*16 CDSQRT
        COMPLEX*16 SRTBB4AC
        COMPLEX*16 ddef, aabc

        COMPLEX*16 ZFAC
c-----
c       first test to see if a fluid layer - if it is fluid, the
c       eigenfunctions are specially computed and we need only the
c       rp
c-----
        if(TL(m).eq.0.0 .or. TN(m).eq.0.0)then
            rp = cdsqrt(dcmplx(wvno2 -omega2*TRho(m)/(TA(m)), 0.0d+00))
            rsv = dcmplx(0.0d+000, 0.0d+00)
            return
        endif


        a = wvn * TF(m) / (TC(m))
        b = 1.0/(TC(m))
        c = - TRho(m)*omg*omg + wvn*wvn *
     1      (TA(m) -TF(m)*TF(m)/(TC(m)))
        d = - wvn
        e = 1.0/(TL(m))
        f = - TRho(m)*omg*omg

c-----
c       do algebra first to avoid numerical problems
c-----
        ddef = wvn*wvn - TRho(m)*omg*omg/(TL(m))
        aabc = wvn*wvn*TA(m)/TC(m) - TRho(m)*omg*omg/(TC(m))

c-----
c       Do the QUASI P and SV - WE MUST BE CAREFUL HERE CONCERNING
c       BRANCH CUTS OF THE SQUARE ROOT
c-----
c-----
c       The characteristic equation to be solved is
c
c       L^4 - L^2[ 2 ad +ec +fb ] + [ (d^2+ef)(a^2+bc)] = 0
c-----
        bb = 2.0d+00 * a*d + e*c +f*b
        cc = ddef * aabc
c----
c       ensure that the SQRT(bb*bb - 4.0D+00*cc) is in the
c       I and II quadrants
c-----

        SRTBB4AC = CDSQRT(bb*bb - 4.0D+00*cc)
        IF(DIMAG(SRTBB4AC) .lt.0.0D+00)THEN
            SRTBB4AC = - SRTBB4AC
        ENDIF
c-----
c       Determine L^2 with care for roundoff
c-----
        IF(DREAL(BB) .LT.0.0D+00 .AND. DREAL(SRTBB4AC).LT.0.0D+00)THEN
            L2(2) = ( bb - SRTBB4AC) / 2.0d+00
            L2(1) = cc/L2(2)
        ELSE
            L2(1) = ( bb + SRTBB4AC) / 2.0d+00
            L2(2) = cc/L2(1)
        ENDIF
c-----
c       Use the Lambda^2 values to form
c       xka^2 == k^2 - L(1)^2
c       xkb^2 == k^2 - L(2)^2
c       Associate the smallest xka, xkb with the P!
c-----
        xka2 = wvno2 - L2(1)
        xkb2 = wvno2 - L2(2)
        if(cdabs(xkb2) .lt. cdabs(xka2))THEN
                ZFAC = L2(1)
                L2(1) = L2(2)
                L2(2) = ZFAC
        endif
        rp  = CDSQRT(L2(1))
        rsv = CDSQRT(L2(2))

c-----
c       get the norms - note that the true norm will be 
c            2  NP amd 2 L(2) NSV
c       The factorization permits us to use the sin nz/n or n sin nz
c-----
C        NP  = (  L2(1)*(-2*a*b*d + 2*a*a*e + b*c*e - b*b*f)
C     1      + (a*a+b*c)*(2*b*d*d - 2*a*d*e + b*e*f - c*e*e) )
C        NSV = (- L2(2)*(2*b*d*d - 2*a*d*e - c*e*e + b*e*f)
C     1      + (d*d+e*f)*(2*a*b*d - 2*a*a*e + b*b*f - b*c*e) )
c-----
c       get the modified eigen vectors x11 and x31 have 
c           lambda1 (rp ) factored out
c               modified eigen vectors x22 and x42 have 
c           lambda2 (rsv) factored out
c-----
        x11 =              (b*d - a*e)
        x21 =  b*L2(1) - e*(b*c + a*a)
        x31 =    L2(1) -   (a*d + c*e)
        x41 = -a*L2(1) + d*(b*c + a*a)

        x12 = -e*L2(2) + b*(d*d + e*f)
        x22 = ( b*d - a*e)
        x32 = d*L2(2) - a*(d*d + e*f)
        x42 = - ( L2(2) -  a*d - b*f)
c-----
c       TEST
c       Force the eigenfunctions to be as given in 7.4.4
c-----
        zfac = rp / x21
        x11  = x11 *zfac
        x21  = rp
        x31  = x31 *zfac
        x41  = x41 *zfac

        zfac = rsv / x12
        x12  = rsv
        x22  = x22 * zfac
        x32  = x32 * zfac
        x42  = x42 * zfac
c-----
c       REDEFINE HERE USING THE adjusted eigenvectors
c       Note that TRUE NP  = 2 * np  * rp
c       Note that TRUE NSV = 2 * nsv * rsv     , where
c       but also note that since the 11 31 22 42 13 33 24 44 have
c       an eigenvalue factored out the normalization using these
c       elements will not use the rp and rsv
c-----  
        np   = x11*x41 - x21*x31
        nsv  = x12*x42 - x22*x32

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
        subroutine hska(AA,tcosp,trsinp,tsinpr,tcossv,trsinsv,tsinsvr,
     1      NP,NSV,
     1      X11, X21, X31, X41,X12, X22, X32, X42,
     2      TRho,iwat,ex,om2)
c-----
c       Changes
c
c       01 May 2002  - defined cosp = tcosp/NP to 
c             reduce nmber of complex divisions
c-----
        implicit none
        real*8 AA(4,4)
        real*8 tcosp , tcossv 
        real*8 trsinp, trsinsv
        real*8 tsinpr, tsinsvr
        COMPLEX*16 NP, NSV
        COMPLEX*16 X11, X21, X31, X41
        COMPLEX*16 X12, X22, X32, X42
        real TRho
        integer iwat
        real*8 ex, dfac
        real*8  om2
        complex*16 zrho

c-----
c       introduce shorthand to reduce work
c-----
        COMPLEX*16 cosp , cossv 
        COMPLEX*16 rsinp, rsinsv
        COMPLEX*16 sinpr, sinsvr

        integer i, j
        zrho = dcmplx(dble(TRho),0.0d+00)
        if(iwat.eq.1)then
c-----
c       fluid layer
c-----
            do 100 j=1,4
                do 101 i=1,4
                    AA(i,j) = 0.0d+00
  101           continue
  100       continue
            if(ex.gt.35.0d+00)then
                dfac = 0.0d+00
            else
                dfac = dexp(-ex)
            endif
            AA(1,1) = dfac
            AA(4,4) = dfac
            AA(2,2) = tcosp
            AA(3,3) = tcosp
            AA(2,3) = -trsinp/(zrho*om2)
            AA(3,2) = - zrho*om2*tsinpr
        else
c-----
c       elastic layer
c-----
            cosp   = tcosp/NP
            sinpr  = tsinpr/NP
            rsinp  = trsinp/NP
            cossv  = tcossv/NSV
            sinsvr = tsinsvr/NSV
            rsinsv = trsinsv/NSV

            AA(1,1) = dreal(  x11*x41*cosp  + x12*x42*cossv  )
            AA(1,2) = dreal(- x11*x31*rsinp - x12*x32*sinsvr )
            AA(1,3) = dreal(- x11*x21*cosp  - x12*x22*cossv  )
            AA(1,4) = dreal(  x11*x11*rsinp + x12*x12*sinsvr )
            AA(2,1) = dreal(  x21*x41*sinpr + x22*x42*rsinsv )
            AA(2,2) = dreal(- x21*x31*cosp  - x22*x32*cossv  )
            AA(2,3) = dreal(- x21*x21*sinpr - x22*x22*rsinsv )
            AA(3,1) = dreal(  x31*x41*cosp  + x32*x42*cossv  )
            AA(3,2) = dreal(- x31*x31*rsinp - x32*x32*sinsvr )
            AA(4,1) = dreal(  x41*x41*sinpr + x42*x42*rsinsv )
            AA(2,4) = - AA(1,3)
            AA(3,3) =   AA(2,2)
            AA(3,4) = - AA(1,2)
            AA(4,2) = - AA(3,1)
            AA(4,3) = - AA(2,1)
            AA(4,4) =   AA(1,1)
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
