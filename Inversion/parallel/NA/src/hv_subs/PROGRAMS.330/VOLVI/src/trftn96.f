        program trftn96
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME VI                                                      c
c                                                                     c
c      PROGRAM: TRFTN96                                               c
c                                                                     c
c      COPYRIGHT 2001                                                 c
c      R. B. Herrmann                                                 c
c      Department of Earth and Atmospheric Sciences                   c
c      Saint Louis University                                         c
c      221 North Grand Boulevard                                      c
c      St. Louis, Missouri 63013                                      c
c      U. S. A.                                                       c
c                                                                     c
c---------------------------------------------------------------------c
c       CHANGES
c       20 DEC 2001 - create this program based on hspec96
c       07 MAY 2002 - program generalized to work
c               with both isotropic and transverse isotropic media
c       25 APR 2002 - implement TI media with vertical symmetry
c       29 OCT 2002 - Fixed evalg and evalh for fluid layers
c       20 JUN 2003 - changed getegn and section evlmat before call hska
c           STILL MUST CORRECT DC CORRECTION SINCE FOR TI MEDIA
c           RP !-> k, RSV !-> k for large k but rather
c       23 JUL 2009 - ensured consistent call to FFT routine with
c               complex arguments. We now call zfour
c
c-----
c       This program creates a SAC file for the 
c           surface receiver function
c
c-----
        integer LER, LIN, LOT
        parameter(LER=0, LIN=5, LOT=6)
c-----
c       earth model information
c-----
        integer NL
        parameter(NL=200)
        common/timodel/D(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      TRho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL), 
     3      frefp(NL), frefs(NL)
        real D, TA, TC, TF, TL, TN, TRho
        real qa, qb, etap, etas, frefp, frefs
        common/modlly/mmax
        integer mmax
        integer iunit, iiso, iflsph, idim, icnvel, ierr
        character mname*80, title*80
c-----
c       time series information
c-----
        integer NSAMP, NFREQ
        parameter (NSAMP=8192,NFREQ=4097)
        complex zdata(NSAMP)
        real x(NSAMP)

        complex Z
        complex ZFAC

        common/cntrl/ishank,hnkarg,dstcor,dokjar
        integer*4 dstcor
        real*4 hnkarg
        logical ishank,dokjar

        common/depref/refdep

        common/damp/alpha,ieqex


c-----
c       fref    R*4 - reference frequency for Causal Q
c       dokjar  L   - .false. use Futterman causal Q
c                 .true.  use Kjartansson Causal Q
c-----


c-----
c       command line arguments
c-----
        logical dop, dotwo
        real delay, dt, rayp, gaussalp
        integer n, iout

        logical ext

c-----
c       call machine dependent initialization
c-----
        call mchdep()
c-----
c       parse command line arguments
c-----
        call gcmdln(dop,dotwo,delay,dt,rayp,gaussalp,n,mname,iout)
c-----
c       get earth model information
c-----
        if(mname.eq. ' ')call usage('Model not specified')  
c-----
c       get the earth model
c-----
        inquire(file=mname,exist=ext)
        if(.not. ext)call usage('Model file does not exist')
        l = lgstr(mname)
        call getmod(1,mname,mmax,title,iunit,iiso,iflsph,
     1          idimen,icnvel,ierr,.false.)
        if(ierr .lt. 0)stop
c-----
c       make sure that we use 1/Q
c-----
        do 3007 i=1,mmax
            if(qa(i).lt.0.0)qa(i) = 0.0
            if(qb(i).lt.0.0)qb(i) = 0.0
            if(qa(i) .gt. 1.0)qa(i) = 1.0/qa(i)
            if(qb(i) .gt. 1.0)qb(i) = 1.0/qb(i)
            if(frefp(i) .le. 0.0)frefp(i) = 1.0
            if(frefs(i) .le. 0.0)frefs(i) = 1.0
 3007   continue
c-----
c       output control parameters
c-----
        write(LOT,*)'dop          :',dop
        write(LOT,*)'delay        :',delay
        write(LOT,*)'dt           :',dt   
        write(LOT,*)'rayp         :',rayp 
        write(LOT,*)'gaussalp     :',gaussalp 
        write(LOT,*)'n            :',n 
        write(LOT,*)'Double Length:',dotwo
        lm = lgstr(mname)
        write(LOT,*)'Model        :',mname(1:lm)
c-----
c       ensure that the number of points is a power of 2
c-----
        call npow2(n)
        if(dotwo)then
            n = 2*n
        endif

        if(dokjar)then
        write(LOT,*)'Kjartansson Constant Q operator used'
        endif
c-----
c       generate a zero phase pulse, with a zero at the 
c           Nyquist Frequency
c-----
        delay = abs(delay)
c-----
c       define the alpha parameter for complex frequency
c-----
        alpha = 2.5/(n*dt)
c-----
c       note this will not work for symmetric pulses with 0 lag ???
c-----
        fac = 1.0
        dfac = exp(-alpha*dt)
        do 1000 i=1,n
            if(i.eq.2)then
                zdata(i) = cmplx(0.25/dt,0.0)*fac
            else if(i.eq.3)then
                zdata(i) = cmplx(0.50/dt,0.0)*fac
            else if(i.eq.4)then
                zdata(i) = cmplx(0.25/dt,0.0)*fac
            else
                zdata(i) = cmplx(0.0,0.0)
            endif
            fac = fac * dfac
 1000   continue
        call zfour(zdata,n,-1,dt,df)
c-----
c       now process in the frequency domain
c       include the desired time shift, but not that the 
c           pulse is already
c       shifted 2*dt units 
c-----
        n21 = n / 2 + 1
        do 4000 i=1,n21
                freq=(i-1)*df
            call excit(freq,dop,rayp,Z,iout)
            zdata(i) = zdata(i) * Z
            fac = - 6.2831853*freq*(delay - 2*dt)
            zdata(i) = zdata(i) * cmplx(cos(fac),sin(fac))
            if(i.gt.1)then
                zdata(n+2 -i ) = conjg(zdata(i))
            endif
 4000   continue
        call zfour(zdata,n,+1,dt,df)
c-----
c       undamp the time series
c-----
            fac = exp(-alpha*delay)
            dfac = exp(alpha*dt) 
            do 425 i = 1,n 
                zdata(i)= zdata(i) * fac
                fac = fac * dfac 
  425       continue 
c-----
c       now Gaussian filter
c-----
        call zfour(zdata,n,-1,dt,df)
        do 426 i=1,n21
            freq=(i-1)*df
            fac = (6.2831853*freq)/(2.0*gaussalp)
            if(fac.gt.25.0)then
                fac = 0.0
            else
                fac = exp( - fac * fac)
            endif
            zdata(i) = zdata(i) * fac
            if(i.gt.1)then
                zdata(n + 2 - i) = conjg(zdata(i))
            endif
  426   continue
c-----
c       The source pulse should have a zero at the Nyquist frequency
c-----
        zdata(n21) = cmplx(0.0,0.0)
        call zfour(zdata,n,+1,dt,df)
        do 427 i=1,n
            x(i) = real(zdata(i))
  427   continue
c-----
c       get extreme values of time series for SAC header
c-----
        if(dotwo)then
            n = n / 2
        endif
        call scmxmn(x,n,depmax,depmin,depmen,indmax,indmin)
c-----
c       output the time series as a SAC file
c-----
        call newhdr()
c-----
c       set real header value
c-----
            call setfhv('DELTA   ',dt,nerr)
            call setfhv('DEPMIN  ',depmin,nerr)
            call setfhv('DEPMAX  ',depmax,nerr)
            call setfhv('DEPMEN  ',depmen,nerr)
            beg = -delay 
            e = beg + (n-1)*dt
            call setfhv('B       ',beg   ,nerr)
            call setfhv('E       ',e     ,nerr)
            call setfhv('TIMMAX  ',beg+indmax*dt, nerr)
            call setfhv('TIMMIN  ',beg+indmin*dt, nerr)
            o = 0.0
            call setfhv('O       ',o     ,nerr)
c-----
c       set integer header value
c-----
            ksyear = 1970
            kdoy = 1
            kshour = 0
            ksmin = 0
            isec = 0
            jsmsec = 0
            call setnhv('NZYEAR  ',ksyear,nerr)
            call setnhv('NZJDAY  ',kdoy  ,nerr)
            call setnhv('NZHOUR  ',kshour,nerr)
            call setnhv('NZMIN   ',ksmin ,nerr)
            call setnhv('NZSEC   ',isec  ,nerr)
            call setnhv('NZMSEC  ',jsmsec,nerr)
            call setnhv('NPTS    ',n, nerr)
c-----
c       set logical header value
c-----
            call setlhv('LEVEN   ',.true.,nerr)
            call setlhv('LPSPOL  ',.true.,nerr)
            call setlhv('LOVROK  ',.true.,nerr)
            call setlhv('LCALDA  ',.true.,nerr)
c-----
c       set character header value
c-----
            call setkhv('KSTNM    ','RFTN    ' ,nerr)
c-----
c       set enumerated header value
c-----
            call setihv('IFTYPE   ','ITIME   ',nerr)
            call setihv('IZTYPE   ','IB      ',nerr)
c-----
c       set deconvolution specific headers
c-----
        call setfhv('USER0   ',gaussalp,nerr)
        call setfhv('USER4   ',rayp,nerr)
        call setkhv('KUSER0  ','Rftn    ',nerr)
        call setkhv('KUSER1  ','FD_DECON',nerr)
        call setkhv('KEVNM   ','Rftn    ',nerr)
        call setfhv('USER5   ',100.0,nerr)
        if(iout.eq.1)then
            call setkhv('KCMPNM  ','Rftn    ',nerr)
        else if(iout.eq.2)then
            call setkhv('KCMPNM  ','Z    ',nerr)
        else if(iout.eq.2)then
            call setkhv('KCMPNM  ','R    ',nerr)
        endif
            
c-----
c       output the time series
c-----
        call bwsac(3,NSAMP,'hrftn96.sac',x)

        end

        subroutine gcmdln(dop,dotwo,delay,dt,rayp,gaussalp,
     1      nsamp,mname,iout)
c-----
c       parse command line arguments
c       requires subroutine mgtarg() and function mnmarg()
c-----
c       dop L   - .true. P-wave incident
c                 .false. SV-wave incident
c       dotwo   L   - .true. computations done with double length,
c                   but only single length is output
c       delay   R   - padding before the beginning of the 
c                     receiver function
c       rayp    R   - ray parameter in sec/km
c       gaussalp R  - Gaussian filter parameter
c       nsamp   I   - number of samples (poer of 2)
c       mname   Ch* - earth model name - required
c       dokjar  L   - 
c       iout    I   - output 1 = rftn, 2 = z, 3 = r time series
c-----
        logical dop, dotwo
        real delay, dt, rayp, gaussalp
        integer nsamp, iout
        character mname*(*)

        common/cntrl/ishank,hnkarg,dstcor,dokjar
        logical ishank,dokjar
        integer*4 dstcor
        real*4 hnkarg



        integer*4 mnmarg
        character*50 name

        IOUT = 1
        DOKJAR = .FALSE.
        ISHANK = .FALSE.
        DSTCOR = 0
        HNKARG = 3.0
        dop = .true.
        dotwo = .false.
        delay = 5.0
        dt = 1.0
        rayp = 0.05
        gaussalp = 1.0
        nsamp = 512
        mname = ' '
        nmarg=mnmarg()
        i = 0
   11   i = i + 1
            if(i.gt.nmarg)goto 13
            call mgtarg(i,name)
            if(name(1:2).eq.'-P' )then
                dop = .true.
            else if(name(1:2).eq.'-S' )then
                dop = .false.
            else if(name(1:2).eq.'-D' .and. name(1:3).ne.'-DT')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,f10.0)')delay
            else if(name(1:3).eq.'-DT')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,f10.0)')dt
            else if(name(1:2).eq.'-N')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,i10)')nsamp
            else if(name(1:5).eq.'-RAYP')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,f10.0)')rayp
            else if(name(1:3).eq.'-DT')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,f10.0)')dt
            else if(name(1:4).eq.'-ALP')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,f10.0)')gaussalp
            else if(name(1:2).eq.'-M')then
                i = i + 1
                call mgtarg(i,mname)
            else if(name(1:2).eq.'-2')then
                dotwo = .true.
            else if(name(1:2).eq.'-z')then
                iout = 2
            else if(name(1:2).eq.'-r')then
                iout = 3
            else if(name(1:2).eq.'-?')then
                call usage(' ')
            else if(name(1:2).eq.'-h')then
                call usage(' ')
            endif
        goto 11
   13   continue
        return
        end

        subroutine usage(ostr)
c------
c       write out program syntax
c-----
        character ostr*(*)
        integer LER,LOT,LIN
        parameter(LER=0, LIN=5, LOT=6)
        if(ostr.ne. ' ' )then
            lostr = lgstr(ostr)
            write(LER,*)ostr(1:lostr)
        endif
        write(LER,*)'USAGE: ',
     1  'trftn96 [-P] [-S] [-2] [-r] [-z] -RAYP p',
     2      ' -ALP alpha -DT dt -NSAMP nsamp'
        write(LER,*)
     1  '-P           (default true )    Incident P wave'
        write(LER,*)
     1  '-S           (default false)    Incident S wave'
        write(LER,*)
     1  '-RAYP p      (default 0.05 )    Ray parameter in sec/km'
        write(LER,*)
     1  '-DT dt       (default 1.0  )    Sample interval for synthetic'
        write(LER,*)
     1  '-NSAMP nsamp (default 512  )    Number samples for synthetic'
        write(LER,*)
     1  '-M   model   (default none )    Earth model name'
        write(LER,*)
     1  '-ALP alp     (default 1.0  )    Number samples for synthetic'
        write(LER,*)
     1  '     H(f) = exp( - (pi freq/alpha)**2) '
        write(LER,*)
     1  '     Filter corner ~ alpha/pi '
        write(LER,*)
     1  '-2           (default false)    Use 2x length internally'
        write(LER,*)
     1  '-r           (default false)    Output radial   time series'
        write(LER,*)
     1  '-z           (default false)    Output vertical time series'
        write(LER,*)
     1  '     -2  (default false) use double length FFT to'
        write(LER,*)
     1  '     avoid FFT wrap around in convolution '
        write(LER,*)
     1  '-D delay     (default 5 sec)    output delay sec before t=0'
        write(LER,*)
     1  '-?                   Display this usage message'
        write(LER,*)
     1  '-h                   Display this usage message'
        stop 
        end

        subroutine dosph(mmax,v,h,ipsvsh)
c-----
c       Transform spherical earth to flat earth
c
c       Schwab, F. A., and L. Knopoff (1972). 
c           Fast surface wave and free
c       mode computations, in  
c           Methods in Computational Physics, Volume 11,
c       Seismology: Surface Waves and Earth Oscillations,  
c           B. A. Bolt (ed),
c       Academic Press, New York
c
c       Love Wave Equations  44, 45 , 41 pp 112-113
c       Rayleigh Wave Equations 102, 108, 109 pp 142, 144
c
c
c-----
c       mmax    I*4 number of layers
c       v() R*4 array of velocities
c       h() R*4 array of layer thicknesses
c       ipsvsh  I*4     1 - get P time
c                       2 - get SV time
c                       3 - get SH time
c-----
        integer*4 mmax
        real*4 v(mmax), h(mmax)
        integer*4 ipsvsh

        double precision z0,z1,r0,r1,dr,ar,tmp

        if(ipsvsh.eq.3)then
            ifunc = 1
        else
            ifunc = 2
        endif

        ar=6370.0d0
        dr=0.0d0
        r0=ar
        h(mmax)=1.0
        do 10 i=1,mmax
C           dr=dr+dble(d(i))
            dr=dr+dble(h(i))
            r1=ar-dr
            z0=ar*dlog(ar/r0)
            z1=ar*dlog(ar/r1)
C           d(i)=z1-z0
            h(i)=z1-z0
            tmp=(ar/r1-ar/r0)*ar/(z1-z0)
            v(i)=v(i)*tmp
C           a(i)=a(i)*tmp
C           b(i)=b(i)*tmp
            if(ifunc.eq.1)then
                tmp=(r0*(r0/ar)**4-r1*(r1/ar)**4)
     1              /(dble(h(i))*5.0d0)
C     1             /(dble(d(i))*5.0d0)
            else if(ifunc.eq.2)then
                tmp=(r0*(r0/ar)-r1*(r1/ar))
     1              /(dble(h(i))*2.0d0)
C     1             /(dble(d(i))*2.0d0)
            endif
C           rho(i)=rho(i)*tmp
            r0 = r1
   10   continue
        h(mmax)=0.0
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
        if(npts.ge.nsamp)return
        npts = 2*npts
        go to 1000
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

        subroutine aten(om,qa,qb,alpha,atna,atnb,frefp,frefs)
        implicit none
c-----
c       om  C*16    complex frequency
c       
c-----
c       make velocities complex, using Futterman or Kjartansson
c       causality operator
c-----
        complex*16 om, atna, atnb
        real*4 qa, qb, frefp, frefs, alpha

        complex*16 at
        real*8 pi, om1p, om1s, oml, fac, pi2
        common/cntrl/ishank,hnkarg,dstcor,dokjar
        integer*4 dstcor
        real*4 hnkarg
        logical ishank,dokjar
        real gama, gamb, rfac
        real*8 CDABS
        complex*16 CDLOG
c-----
c       reference frequency is fref hz
c-----
        om1p=6.2831853*frefp
        om1s=6.2831853*frefs
        pi2 = 1.5707963
        pi=3.1415927d+00
        if(dokjar)then
c-----
c       Kjartansson Constant Q, causal Q operator
c       Kjartansson, E. (1979). 
c           Constant Q-wave propagation and attenuation,
c       J. Geophys. Res. 84, 4737-4748.
c-----
            gama = atan(qa)/pi
            gamb = atan(qb)/pi
            if(gama.le.0.0)then
                atna = dcmplx(1.0d+00,0.0d+00)
            else
                fac = pi2*gama
                rfac = sin(fac)/cos(fac)
                atna = dcmplx(1.0d+00,0.0d+00)/
     1              (( (om/om1p)**dble(-gama) ) *
     2               dcmplx(1.0d+00,-dble(rfac)))
            endif
            if(gamb.le.0.0)then
                atnb = dcmplx(1.0d+00,0.0d+00)
            else
                fac = pi2*gamb
                rfac = sin(fac)/cos(fac)
                atnb = dcmplx(1.0d+00,0.0d+00)/
     1          (( (om/om1s)**dble(-gamb) ) *
     2           dcmplx(1.0d+00,-dble(rfac)))
            endif
        else
c-----
c       Futterman Causal Q
c-----
c           low frequency cutoff is 0.01 hz
c-----
            oml=0.062831853d+00
            atna=dcmplx(1.0d+00,0.0d+00)
            atnb=dcmplx(1.0d+00,0.0d+00)
            if(qa.gt.0.0)then
              at=dcmplx(0.0d+00,0.0d+00)
              if(CDABS(om).gt.oml) at=CDLOG(om/om1p)/pi
              if(CDABS(om).le.oml) then
                fac=dsqrt(oml*oml + dble(alpha*alpha))/oml
              at=CDLOG(dcmplx(dble(oml),-dble(alpha))/(om1p*fac))/pi
              endif
              atna=(1.+dble(qa)*at+dcmplx(0.0d+00,dble(qa/2.)))
            endif
            if(qb.gt.0.0)then
              at=dcmplx(0.0d+00,0.0d+00)
              if(CDABS(om).gt.oml) at=CDLOG(om/om1s)/pi
              if(CDABS(om).le.oml) then
                fac=dsqrt(oml*oml + dble(alpha*alpha))/oml
              at=CDLOG(dcmplx(dble(oml),-dble(alpha))/(om1s*fac))/pi
              endif
               atnb=(1.+dble(qb)*at+dcmplx(0.0d+00,dble(qb/2.)))
            endif
        endif
        return
        end
        subroutine cmult(e,ca,exa,exe)
        implicit none
        common/lwater/lfluid
        logical lfluid
c-----
c       FORM EC where e(1x5) c(5x5)
        complex*16 e(5)
        complex*16 ca(5,5)
        real*8 exa,exe,eval
        real *8 xnorm
        complex*16 c, ee(5)
        integer IUP, i, j

        if(lfluid)then
            IUP = 2
        else
            IUP = 5
        endif
        do 1350 i=1,IUP
            c = dcmplx(0.0d+00,0.0d+00)
            do 1349 j=1,IUP
                c=c+ e(j) * ca(j,i)
 1349       continue
            ee(i)=c
 1350   continue
        exe = exe + exa
C       call normc(ee,eval,xnorm)
C       do 1351 i=1,IUP
C           e(i) = ee(i)*xnorm
C 1351  continue
C       exe = exe + eval
        return
        end

        subroutine copy5(ca,dar,m,itofrm,dex,exa)
        implicit none
        integer NL
        parameter (NL=200)
        complex*16 ca(5,5)
        complex*16 dar(NL,5,5)
        integer m, itofrm
        real*8 dex(NL)
        real*8 exa

        integer i, j
c-----
c       copy from ca to dar
c-----
        if(itofrm.eq.0)then
            do 100 j=1,5
                do 110 i=1,5
                    dar(m,i,j) = ca(i,j)
  110           continue
                dex(m) = exa
  100       continue
c-----
c       copy from dar to ca
c-----
        else
            do 200 j=1,5
                do 210 i=1,5
                    ca(i,j) = dar(m,i,j)
  210           continue
                exa = dex(m)
  200       continue
        endif
        return
        end

        subroutine copy4(aa,har,m,itofrm,hex,ex)
c-----
c       copy between aa and har arrays for m'th layer
c
c       aa  R   - 4x4 Haskell matrix array
c       har R   - NLx4x4 storage array
c       m   I   - layer index
c       itofrm  I   - 0 copy aa to har, ex to hex
c                 !=0 copy from har to aa, hex to ex
c       hex R   - NL array - storage for exponent
c       ex  R   - exponent value
c-----
        implicit none
        integer NL
        parameter (NL=200)
        complex*16 aa(4,4)
        complex*16 har(NL,4,4)
        integer m, itofrm
        real*8 hex(NL)
        real*8 ex

        integer i, j
c-----
c       copy from aa to har
c-----
        if(itofrm.eq.0)then
            do 100 j=1,4
                do 110 i=1,4
                    har(m,i,j) = aa(i,j)
  110           continue
                hex(m) = ex
  100       continue
c-----
c       copy from har to aa
c-----
        else
            do 200 j=1,4
                do 210 i=1,4
                    aa(i,j) = har(m,i,j)
  210           continue
                ex = hex(m)
  200       continue
        endif
        return
        end

        subroutine copy2(hl,hal,m,itofrm,lex,exb)
        implicit none
        integer NL
        parameter (NL=200)
        complex*16 hl(2,2)
        complex*16 hal(NL,2,2)
        integer m, itofrm
        real*8 lex(NL)
        real*8 exb

        integer i, j
c-----
c       copy from hl to hal
c-----
        if(itofrm.eq.0)then
            do 100 j=1,2
                do 110 i=1,2
                    hal(m,i,j) = hl(i,j)
  110           continue
                lex(m) = exb
  100       continue
c-----
c       copy from hal to hl
c-----
        else
            do 200 j=1,2
                do 210 i=1,2
                    hl(i,j) = hal(m,i,j)
  210           continue
                exb = lex(m)
  200       continue
        endif
        return
        end

        subroutine copy(da,aa)
c-----
c       copy da array to aa
c-----
        implicit none
        complex*16 da(4,4), aa(4,4)
        integer i, j
        
        do 1000 j=1,4
            do 1100 i=1,4
                aa(i,j) = da(i,j)
 1100       continue
 1000   continue
        return
        end

        subroutine dmult(da,aa)
c-----
c       propagator up
c       FORM D = DA
c-----
        implicit none
        complex*16 aa(4,4)
        complex*16 sumd,ea(4,4),da(4,4)
        integer i, j, jj
        do 1360 i=1,4
            do 1361 j=1,4
                sumd = dcmplx(0.0d+00,0.0d+00)
                do 1362 jj=1,4
                    sumd=sumd+da(i,jj) * aa(jj,j)
 1362           continue
                ea(i,j)=sumd
 1361       continue
 1360   continue
        do 1363 j=1,4
            do 1364 i=1,4
                da(i,j)=ea(i,j)
 1364       continue
 1363   continue
        return
        end

        subroutine dnka(CA,cosp,rsinp,sinpr,cossv,rsinsv,sinsvr,NP,NSV,
     1      x11,x21,x31,x41,x12,x22,x32,x42,
     1      rho,iwat,ex,om2)
        implicit none
        complex*16 CA(5,5)
        COMPLEX*16 cosp , cossv 
        COMPLEX*16 rsinp, rsinsv
        COMPLEX*16 sinpr, sinsvr
        COMPLEX*16 NP, NSV
        COMPLEX*16 X11, X21, X31, X41
        COMPLEX*16 X12, X22, X32, X42
        real rho
        integer iwat
        real*8 ex, dfac
        complex*16  om2
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
c       the number of multiplications can be reduced from 36 to 25 
c           if we define a new matrices
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
        zrho = dcmplx(dble(rho),0.0d+00)
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
                    ca(i,j) = dcmplx(0.0d+00, 0.0d+00)
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
            CA(1,1) = TCA(1,1)
            CA(1,2) = TCA(1,2)
            CA(1,3) = TCA(1,3)
            CA(1,4) = TCA(1,5)
            CA(1,5) = TCA(1,6)

            CA(2,1) = TCA(2,1)
            CA(2,2) = TCA(2,2)
            CA(2,3) = TCA(2,3)
            CA(2,4) = TCA(2,5)
            CA(2,5) = TCA(1,5)

            CA(3,1) =   2.0d+00*TCA(3,1)
            CA(3,2) =   2.0d+00*TCA(3,2)
c-----
c       beware of normalization
c-----
            CA(3,3) =   2.0d+00*TCA(3,3) - 1.0d+00*dfac
            CA(3,4) =  -2.0d+00*TCA(2,3)
            CA(3,5) =  -2.0d+00*TCA(1,3)

            CA(4,1) =   TCA(5,1)
            CA(4,2) =   TCA(5,2)
            CA(4,3) =  -TCA(3,2)
            CA(4,4) =   TCA(2,2)
            CA(4,5) =   TCA(1,2)

            CA(5,1) =   TCA(6,1)
            CA(5,2) =   TCA(5,1)
            CA(5,3) =  -TCA(3,1)
            CA(5,4) =   TCA(2,1)
            CA(5,5) =   TCA(1,1)
        endif
        return
        end

        subroutine equlck()
        implicit none
        integer NL
        parameter(NL=200)
        common/timodel/D(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      TRho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL), 
     3      frefp(NL), frefs(NL)
        real D, TA, TC, TF, TL, TN, TRho
        real qa, qb, etap, etas, frefp, frefs
        common/modlly/mmax
        integer mmax
        common/updnsm/equalu(NL), equald(NL)
        logical equalu, equald

        integer m
c-----
c       To avoid repeated computation, check to see if 
c           neighboring layers
c       are identical, once for going up and another for going down
c-----
c       First check top down
c-----
        do 100 m=1,mmax
            if(m.eq.1)then
                equald(m) = .false.
            else if(m.gt.1
     1          .and. TA(m).eq.TA(m-1) 
     1          .and. TC(m).eq.TC(m-1) 
     1          .and. TF(m).eq.TF(m-1) 
     1          .and. TL(m).eq.TL(m-1) 
     1          .and. TN(m).eq.TN(m-1) 
     1          .and.  D(m).eq. D(m-1) 
     4          .and. TRho(m).eq.TRho(m-1)
     5          .and. qa(m).eq.qa(m-1)
     6          .and. qb(m).eq.qb(m-1) )then
                equald(m) = .true.
            else
                equald(m) = .false.
            endif
  100   continue
c-----
c       check bottom up
c-----
        do 200 m=1,mmax
            if(m.eq.mmax)then
                equalu(m) = .false.
            else if(m.lt.mmax
     1          .and. TA(m).eq.TA(m+1) 
     1          .and. TC(m).eq.TC(m+1) 
     1          .and. TF(m).eq.TF(m+1) 
     1          .and. TL(m).eq.TL(m+1) 
     1          .and. TN(m).eq.TN(m+1) 
     1          .and.  D(m).eq. D(m+1) 
     4          .and. Trho(m).eq.Trho(m+1)
     5          .and. qa(m).eq.qa(m+1)
     6          .and. qb(m).eq.qb(m+1) )then
                equalu(m) = .true.
            else
                equalu(m) = .false.
            endif
  200   continue
        return
        end

        subroutine evalg(jbdry,m,m1,gbr,gbl,in,
     1      wvno,om,om2,wvno2)
        implicit none
        integer jbdry, m, m1, in
        complex*16 gbr(2,5), gbl(2,2)
        complex*16 wvno,om,wvno2,om2
        integer NL
        parameter(NL=200)
        common/timodel/D(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      TRho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL), 
     3      frefp(NL), frefs(NL)
        real D, TA, TC, TF, TL, TN, TRho
        real qa, qb, etap, etas, frefp, frefs
        common/modlly/mmax
        integer mmax
        common/damp/alpha,ieqex
        real alpha
        integer ieqex
        common/lwater/lfluid
        logical lfluid
        complex*16 CDSQRT

        complex*16 atna, atnb
        complex*16 cg(6)
        complex*16 g(4,4)
        COMPLEX*16 zA,zB,zC,zD,zE,zF
        COMPLEX*16 rp, rsv, rsh
        COMPLEX*16 NP, NSV
        COMPLEX*16 X11, X21, X31, X41
        COMPLEX*16 X12, X22, X32, X42
        integer iwat


c-----
c       set up halfspace conditions
c-----
            if(TL(m).eq. 0.0 .or.TN(m).eq.0.0)then
                iwat = 1
            else
                iwat = 0
            endif
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
            if(TL(m) .gt. 0.0 .and. TN(m).gt.0.0)then
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
                if(lfluid)then
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
            call aten(om,qa(m),qb(m),alpha,atna,atnb,frefp(m),frefs(m))
            call getegn(za,zb,zc,zd,ze,zf,om2,wvno2,rp, rsv, rsh,
     1          x11,x21,x31,x41,x12,x22,x32,x42,NP, NSV,m,om,wvno,
     2          atna,atnb)
            if(iwat.eq.0)then
c-----
c               ELASTIC HALFSPACE
C       WRITE(0,*)'X:'
C       WRITE(0,*)x11,x21,x31,x41
C       WRITE(0,*)x12,x22,x32,x42
C       WRITE(0,*)'RP,RSV,RSH:',RP,RSV,RSH
C       WRITE(0,*)'NP,NSV,m,om,wvno:',NP,NSV,m,om,wvno
C       WRITE(0,*)'om2,wvno2:',om2,wvno2


                G(1,1) =    x41/(2.0*NP  *rp  )
                G(1,2) = -  x31/(2.0*NP       )
                G(1,3) = -  x21/(2.0*NP  *rp  )
                G(1,4) =    x11/(2.0*NP       )
                G(2,1) =    x42/(2.0*NSV      )
                G(2,2) = -  x32/(2.0*NSV *rsv )
                G(2,3) = -  x22/(2.0*NSV      )
                G(2,4) =    x12/(2.0*NSV *rsv )
                G(3,1) = -  x41/(2.0*NP  *rp  )
                G(3,2) = -  x31/(2.0*NP       )
                G(3,3) =    x21/(2.0*NP  *rp  )
                G(3,4) =    x11/(2.0*NP       )
                G(4,1) =    x42/(2.0*NSV      )
                G(4,2) =    x32/(2.0*NSV *rsv )
                G(4,3) = -  x22/(2.0*NSV      )
                G(4,4) = -  x12/(2.0*NSV *rsv )
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

                gbl(in,1) = atnb*atnb*TL(m)*rsh
                gbl(in,2) = dcmplx(1.0d+00,0.0d+00)
            else if(iwat.eq.1)then
c-----
c               FLUID HALFSPACE
c-----
                if(lfluid)then
                    gbr(in,1) = dble(TRho(m))*om2
                    gbr(in,2) = -rp
                    gbr(in,3) = dcmplx(0.0d+00,0.0d+00)
                    gbr(in,4) = dcmplx(0.0d+00,0.0d+00)
                    gbr(in,5) = dcmplx(0.0d+00,0.0d+00)
                else
                    gbr(in,1) = dcmplx(0.0d+00,0.0d+00)
                    gbr(in,2) = dcmplx(0.0d+00,0.0d+00)
                    gbr(in,3) = dcmplx(0.0d+00,0.0d+00)
                    gbr(in,4) = -dble(Trho(m))*om2
                    gbr(in,5) = rp
                endif
                gbl(in,1) = dcmplx(0.0d+00,0.0d+00)
                gbl(in,2) = dcmplx(1.0d+00,0.0d+00)
            endif
        else if(jbdry.eq.1)then
c-----
c       FREE - check properties of layer above
c-----
            if(TL(m) .gt. 0.0 .and. TN(m).gt.0.0)then
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
                if(lfluid)then
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

        subroutine evalh(jbdry,m,m1,hsr,hsl,in,wvno,om,om2,wvno2)
        implicit none
        integer jbdry,m,m1,in
        complex*16 hsr(2,5), hsl(2,2)
        complex*16 wvno,om,wvno2,om2
        integer NL
        parameter(NL=200)
        common/timodel/D(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      TRho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL), 
     3      frefp(NL), frefs(NL)
        real D, TA, TC, TF, TL, TN, TRho
        real qa, qb, etap, etas, frefp, frefs
        common/modlly/mmax
        integer mmax
        common/damp/alpha,ieqex
        real alpha
        integer ieqex
        complex*16 CDSQRT

        complex*16 atna, atnb
        complex*16 HM(4,4)
        complex*16 CH(6)
        COMPLEX*16 zA,zB,zC,zD,zE,zF
        COMPLEX*16 rp, rsv, rsh
        COMPLEX*16 NP, NSV
        COMPLEX*16 X11, X21, X31, X41
        COMPLEX*16 X12, X22, X32, X42
        integer iwat


c-----
c       set up surface conditions
c-----
        call aten(om,qa(m),qb(m),alpha,atna,atnb,frefp(m),frefs(m))
            if(TL(m).eq. 0.0 .or.TN(m).eq.0.0)then
                iwat = 1
            else
                iwat = 0
            endif
        
c-----
c       do top surface condition
c
c           = -1  RIGID
c           =  0  ELASTIC
c           = +1  FREE SURFACE
c-----
        if(iwat.eq.0)then
            if(jbdry.eq.0)then
c-----
c           ELASTIC ABOVE SOLID
c-----
c-----
c       Now give the H matrix
c-----
        call getegn(za,zb,zc,zd,ze,zf,om2,wvno2,rp, rsv, rsh,
     1  x11,x21,x31,x41,x12,x22,x32,x42,NP, NSV,m,om,wvno,
     2  atna,atnb)
                HM(1,1) =   x11*rp
                HM(1,2) =   x12
                HM(1,3) = - x11*rp
                HM(1,4) =   x12
                HM(2,1) =   x21
                HM(2,2) =   x22*rsv
                HM(2,4) =  -x22*rsv
                HM(2,3) =   x21
                HM(3,1) =   x31*rp
                HM(3,2) =   x32
                HM(3,3) = - x31*rp
                HM(3,4) =   x32
                HM(4,1) =   x41
                HM(4,2) =   x42*rsv
                HM(4,3) =   x41
                HM(4,4) =  -x42*rsv
C Compound matrices
c CH(1) = H 12 12 = 11 22 - 12 21
c CH(2) = H 13 12 = 11 32 - 12 31
c CH(3) = H 14 12 = 11 42 - 12 41
c CH(4) = H 23 12 = 21 32 - 31 22
c CH(5) = H 24 12 = 21 42 - 41 22
c CH(6) = H 34 12 = 31 42 - 41 32
                CH(1) = HM(1,1)*HM(2,2) - HM(1,2)*HM(2,1)
                CH(2) = HM(1,1)*HM(3,2) - HM(1,2)*HM(3,1)
                CH(3) = HM(1,1)*HM(4,2) - HM(1,2)*HM(4,1)
                CH(4) = HM(2,1)*HM(3,2) - HM(3,1)*HM(2,2)
                CH(5) = HM(2,1)*HM(4,2) - HM(4,1)*HM(2,2)
                CH(6) = HM(3,1)*HM(4,2) - HM(4,1)*HM(3,2)
                hsr(in,1) =         CH(1)
                hsr(in,2) =         CH(2)
                hsr(in,3) = 2.0d+00*CH(3)
                hsr(in,4) =         CH(5)
                hsr(in,5) =         CH(6)
                hsl(in,1) = dcmplx(1.0d+00,0.0d+00)
                hsl(in,2) = atnb*atnb*TL(m)*rsh
            else if(jbdry.eq.-1)then
c-----
c           RIGID ABOVE SOLID
c-----
                hsr(in,1) = dcmplx(0.0d+00,0.0d+00)
                hsr(in,2) = dcmplx(0.0d+00,0.0d+00)
                hsr(in,3) = dcmplx(0.0d+00,0.0d+00)
                hsr(in,4) = dcmplx(0.0d+00,0.0d+00)
                hsr(in,5) = dcmplx(1.0d+00,0.0d+00)
                hsl(in,1) = dcmplx(0.0d+00,0.00d+00)
                hsl(in,2) = dcmplx(1.0d+00,0.00d+00)
            else if(jbdry.eq.1)then
c-----
c           FREE ABOVE SOLID
c-----
                hsr(in,1) = dcmplx(1.0d+00,0.0d+00)
                hsr(in,2) = dcmplx(0.0d+00,0.0d+00)
                hsr(in,3) = dcmplx(0.0d+00,0.0d+00)
                hsr(in,4) = dcmplx(0.0d+00,0.0d+00)
                hsr(in,5) = dcmplx(0.0d+00,0.0d+00)
                hsl(in,1) = dcmplx(1.0d+00,0.00d+00)
                hsl(in,2) = dcmplx(0.0d+00,0.00d+00)
            endif
        else if(iwat.gt.0)then
            if(jbdry.eq.0)then
c-----
c           HALFSPACE ABOVE FLUID
c-----
                hsr(in,1) = rp
                hsr(in,2) = -dble(TRho(m))*om2
                hsr(in,3) = dcmplx(0.0d+00,0.0d+00)
                hsr(in,4) = dcmplx(0.0d+00,0.0d+00)
                hsr(in,5) = dcmplx(0.0d+00,0.0d+00)
                hsl(in,1) = dcmplx(1.0d+00,0.00d+00)
                hsl(in,2) = dcmplx(0.0d+00,0.00d+00)
            else if(jbdry.eq.-1)then
c-----
c           RIGID ABOVE FLUID
c-----
                hsr(in,1) = dcmplx(0.0d+00,0.0d+00)
                hsr(in,2) = dcmplx(1.0d+00,0.0d+00)
                hsr(in,3) = dcmplx(0.0d+00,0.0d+00)
                hsr(in,4) = dcmplx(0.0d+00,0.0d+00)
                hsr(in,5) = dcmplx(0.0d+00,0.0d+00)
                hsl(in,1) = dcmplx(1.0d+00,0.00d+00)
                hsl(in,2) = dcmplx(0.0d+00,0.00d+00)
            else if(jbdry.eq.1)then
c-----
c           FREE ABOVE FLUID
c-----
                hsr(in,1) = dcmplx(1.0d+00,0.0d+00)
                hsr(in,2) = dcmplx(0.0d+00,0.0d+00)
                hsr(in,3) = dcmplx(0.0d+00,0.0d+00)
                hsr(in,4) = dcmplx(0.0d+00,0.0d+00)
                hsr(in,5) = dcmplx(0.0d+00,0.0d+00)
                hsl(in,1) = dcmplx(1.0d+00,0.00d+00)
                hsl(in,2) = dcmplx(0.0d+00,0.00d+00)
            endif
        endif
        return
        end

        subroutine evemat(om, om2, wvno, wvno2, g,rp,rsv)
c-----
c                
c       This is E
c                 N
c-----
        implicit none
        integer NL
        complex*16 om,wvno,wvno2,om2
        complex*16 g(4,4)
        parameter(NL=200)
        common/timodel/D(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      TRho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL), 
     3      frefp(NL), frefs(NL)
        real D, TA, TC, TF, TL, TN, TRho
        real qa, qb, etap, etas, frefp, frefs
        common/modlly/mmax
        integer mmax
        common/damp/alpha,ieqex
        real alpha
        integer ieqex
c-----
c       matrix components in layers and boundaries saved
c-----
        complex*16 atna, atnb
        COMPLEX*16 zA,zB,zC,zD,zE,zF
        COMPLEX*16 rp, rsv, rsh
        COMPLEX*16 NP, NSV
        COMPLEX*16 X11, X21, X31, X41
        COMPLEX*16 X12, X22, X32, X42

        call aten(om,qa(mmax),qb(mmax),alpha,atna,atnb,
     1      frefp(mmax),frefs(mmax))
        call getegn(za,zb,zc,zd,ze,zf,om2,wvno2,rp, rsv, rsh,
     1      x11,x21,x31,x41,x12,x22,x32,x42,NP, NSV,mmax,om,wvno,
     2      atna,atnb)
        g(1,1) =   x11*rp
        g(1,2) =   x12
        g(1,3) = - x11*rp
        g(1,4) =   x12
        g(2,1) =   x21
        g(2,2) =   x22*rsv
        g(2,3) =   x21
        g(2,4) = - x22*rsv
        g(3,1) =   x31*rp
        g(3,2) =   x32
        g(3,3) = - x31*rp
        g(3,4) =   x32
        g(4,1) =   x41
        g(4,2) =   x42*rsv
        g(4,3) =   x41
        g(4,4) = - x42*rsv
        return
        end

        subroutine evgmat(om, om2, wvno, wvno2, g)
c-----
c                -1
c       This is E
c                 N
c-----
        implicit none
        integer NL
        complex*16 om,wvno,wvno2,om2
        complex*16 g(4,4)
        parameter(NL=200)
        common/timodel/D(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      TRho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL), 
     3      frefp(NL), frefs(NL)
        real D, TA, TC, TF, TL, TN, TRho
        real qa, qb, etap, etas, frefp, frefs
        common/modlly/mmax
        integer mmax
        common/damp/alpha,ieqex
        real alpha
        integer ieqex
c-----
c       matrix components in layers and boundaries saved
c-----
        complex*16 atna, atnb
        COMPLEX*16 zA,zB,zC,zD,zE,zF
        COMPLEX*16 rp, rsv, rsh
        COMPLEX*16 NP, NSV
        COMPLEX*16 X11, X21, X31, X41
        COMPLEX*16 X12, X22, X32, X42

        call aten(om,qa(mmax),qb(mmax),alpha,atna,atnb,
     1      frefp(mmax),frefs(mmax))
        call getegn(za,zb,zc,zd,ze,zf,om2,wvno2,rp, rsv, rsh,
     1      x11,x21,x31,x41,x12,x22,x32,x42,NP, NSV,mmax,om,wvno,
     2      atna,atnb)
        g(1,1) =    x41/(2.0*NP  *rp  )
        g(1,2) = -  x31/(2.0*NP       )
        g(1,3) = -  x21/(2.0*NP  *rp  )
        g(1,4) =    x11/(2.0*NP       )
        g(2,1) =    x42/(2.0*NSV      )
        g(2,2) = -  x32/(2.0*NSV *rsv )
        g(2,3) = -  x22/(2.0*NSV      )
        g(2,4) =    x12/(2.0*NSV *rsv )
        g(3,1) = -  x41/(2.0*NP  *rp  )
        g(3,2) = -  x31/(2.0*NP       )
        g(3,3) =    x21/(2.0*NP  *rp  )
        g(3,4) =    x11/(2.0*NP       )
        g(4,1) =    x42/(2.0*NSV      )
        g(4,2) =    x32/(2.0*NSV *rsv )
        g(4,3) = -  x22/(2.0*NSV      )
        g(4,4) = -  x12/(2.0*NSV *rsv )
        return
        end

        subroutine evlmat(om,wvno,jbdrys,jbdryh,wvno2,om2)
        implicit none
        complex*16 om, wvno, om2, wvno2
        integer jbdrys,jbdryh
        integer NL
        parameter(NL=200)
        common/timodel/D(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      TRho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL), 
     3      frefp(NL), frefs(NL)
        real D, TA, TC, TF, TL, TN, TRho
        real qa, qb, etap, etas, frefp, frefs
        common/modlly/mmax
        integer mmax
        integer NSOURCE, NRECEIVER, NSR
        parameter(NSOURCE=100,NRECEIVER=100,NSR=100)
        common/source/depths(NSOURCE),lmaxs(NSOURCE),mdpths
        real*4 depths
        integer lmaxs, mdpths
        common/receiv/depthr(NRECEIVER),lmaxr(NRECEIVER),mdpthr
        real*4 depthr
        integer lmaxr, mdpthr
        common/damp/alpha,ieqex
        real alpha
        integer ieqex
        complex*16 p,q,r
        complex*16  ca(5,5), hl(2,2)
        complex*16 aa(4,4)
        complex*16 zone
c-----
c       matrix components in layers and boundaries saved
c-----
        common/hamat/har
        complex*16 har(NL,4,4)
        common/damat/dar
        complex*16 dar(NL,5,5)
        common/hsrfr/hsr
        complex*16 hsr(2,5)
        common/gbrfr/gbr
        complex*16 gbr(2,5)
        common/hlmat/hal
        complex*16 hal(NL,2,2)
        common/hsrfl/hsl
        complex*16 hsl(2,2)
        common/gbrfl/gbl
        complex*16 gbl(2,2)
        common/hexex/hex
        real*8 hex(NL)
        common/hexexw/hexw
        real*8 hexw(NL)
        common/dexex/dex
        real*8 dex(NL)
        common/lexex/lex 
        real*8 lex(NL)
        common/water/iwater(NL),iwats(2),iwatb(2)
        integer iwater, iwats, iwatb
        common/updnsm/equalu(NL), equald(NL)
        logical equalu, equald
        logical compute
        complex*16 CDSQRT

        complex*16 atna,atnb
        COMPLEX*16 zA,zB,zC,zD,zE,zF
        COMPLEX*16 rp, rsv, rsh
        real*8 pex,svex,shex, cpex
        COMPLEX*16 cosp , cossv , cossh
        COMPLEX*16 rsinp, rsinsv, rsinsh
        COMPLEX*16 sinpr, sinsvr, sinshr
        COMPLEX*16 NP, NSV
        COMPLEX*16 X11, X21, X31, X41
        COMPLEX*16 X12, X22, X32, X42


        real*8 dfac
        integer m, iwat
        zone  = dcmplx(1.0d+00,0.0d+00)
c-----
c       evaluate the G matrix 
c           gbr(1,x) is for normal stack
c           gbr(2,x) is for inverted stack
c-----
        call evalg(jbdryh,mmax,mmax-1,gbr,gbl,1,wvno,om,om2,wvno2)
        call evalg(jbdrys,1,   2,     gbr,gbl,2,wvno,om,om2,wvno2)
c-----
c       evaluate the H matrix
c           hsr(1,x) is for normal stack
c           hsr(2,x) is for inverted stack
c-----
        call evalh(jbdrys,1,   2,     hsr,hsl,1,wvno,om,om2,wvno2)
        call evalh(jbdryh,mmax,mmax-1,hsr,hsl,2,wvno,om,om2,wvno2)
c-----

c-----
c       matrix multiplication from bottom layer upward
c-----
        do 1340 m = 1,mmax,1
c-----
c       first check to see if computations already done
c-----
            if(equald(m))then
                compute = .false.
            else
                compute = .true.
            endif
            if(compute)then
            if(TL(m).eq. 0.0 .or.TN(m).eq.0.0)then
                iwat = 1
            else
                iwat = 0
            endif
                call aten(om,qa(m),qb(m),alpha,atna,atnb,
     1          frefp(m),frefs(m))
                call getegn(za,zb,zc,zd,ze,zf,om2,wvno2,rp, rsv, rsh,
     1          x11,x21,x31,x41,x12,x22,x32,x42,NP, NSV,m,om,wvno,
     2          atna,atnb)
c-----
c       now evaluate the various hyperbolic functions
c-----
                p = rp  * dble(d(m))
                q = rsv * dble(d(m))
                r = rsh * dble(d(m))
                call var(p,q,r, rp, rsv, rsh,
     1              cosp, cossv, cossh, rsinp, rsinsv, rsinsh,
     1              sinpr, sinsvr, sinshr,pex,svex,shex,iwat)
                call dnka(CA,cosp,rsinp,sinpr,cossv,rsinsv,sinsvr,
     1              NP,NSV,
     1              x11,x21,x31,x41,x12,x22,x32,x42,
     1              Trho(m),iwat,pex+svex,om2)
c-----
c       For isotropic, Re p > Re sv do adjsut only the
c           SV part of the Haskell matrix
c       For TI this is not always true, so carefully adjust
c-----
            if(pex .gt. svex)then
c-----
c               PEX > SVEX, factor out PEX
c-----
                if((pex-svex).gt. 40.0d+00)then
                    dfac = 0.0d+00
                else
                    dfac = dexp(-(pex-svex))
                endif
                cpex = pex
                call hska(AA,cosp,rsinp,sinpr,
     1              dfac*cossv,dfac*rsinsv,dfac*sinsvr,NP,NSV,
     1              X11, X21, X31, X41,X12, X22, X32, X42, 
     2              Trho(m),iwat,pex,om2)
        else
c-----
c               SVEX > PEX, factor out SVEX
c-----
                if((svex-pex).gt. 40.0d+00)then
                    dfac = 0.0d+00
                else
                    dfac = dexp(-(svex-pex))
                endif
                cpex = svex
                call hska(AA,dfac*cosp,dfac*rsinp,dfac*sinpr,
     1              cossv,rsinsv,sinsvr,NP,NSV,
     1              X11, X21, X31, X41,X12, X22, X32, X42, 
     2              Trho(m),iwat,pex,om2)
        endif
                call hskl(cossh,rsinsh,sinshr,TL(m),
     1              iwat,hl,atnb)
            endif
            iwater(m) = iwat
            call copy5(ca,dar,m,0,dex,pex+svex)
            call copy4(aa,har,m,0,hex,cpex)
            call copy2(hl,hal,m,0,lex,shex)
 1340   continue
        return
        end

        subroutine excit(freq,dop,rayp,Z,iout)
        implicit none
        real freq, rayp
        logical dop
        complex Z
        integer iout

        integer NL
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
            real d, a, b, rho, qa, qb, etap, etas, frefp, frefs
        common/modlly/mmax
            integer mmax
c-----
c       matrix components in layers and boundaries saved
c-----
        common/damp/alpha,ieqex
            real alpha
            integer ieqex
c-----
c       matrix components in layers and boundaries saved
c-----
        common/hamat/har
            complex*16 har(NL,4,4)
        common/damat/dar
            complex*16 dar(NL,5,5)
        common/hsrfr/hsr
            complex*16 hsr(2,5)
        common/gbrfr/gbr
            complex*16 gbr(2,5)
        common/hlmat/hal
            complex*16 hal(NL,2,2)
        common/hsrfl/hsl
            complex*16 hsl(2,2)
        common/gbrfl/gbl
            complex*16 gbl(2,2)
        common/hexex/hex
            real*8 hex(NL)
        common/hexexw/hexw
            real*8 hexw(NL)
        common/dexex/dex
            real*8 dex(NL)
        common/lexex/lex 
            real*8 lex(NL)
        common/water/iwater, iwats, iwatb
            integer iwater(NL), iwats(2), iwatb(2)
        common/updnsm/equalu, equald
            logical equalu(NL), equald(NL)
         logical compute

        complex*16 RA, RB
c-----
c       internal variables
c-----
        complex*16 om, om2, wvno, wvno2
        complex*16 aa(4,4), da(4,4), g(4,4)
        real*8 ex, exe, exa, exl
        complex*16   cd(5),e(5),fr, ca(5,5)
        
        real omega
        integer m, i, j
        integer jbdrys, jbdryh


        omega = 6.2831853*freq
        om =  dcmplx(dble(omega), dble(-alpha))
        om2 = om * om
        wvno  = dcmplx(dble(rayp),0.0d+00)*om
        wvno2 = wvno * wvno

        exe = 0.0d+00
        exa = 0.0d+00
        exl = 0.0d+00
        ex  = 0.0d+00
c-----
c       First evaluate all the Haskell matrices
c       Surface is free  , jbdrys = 1
c       Base is halfspace, jbdryh = 0
c-----
        jbdrys = 1
        jbdryh = 0
c-----
c       Check the model to avoid excessive computation
c-----
        call equlck()
c-----
c       evaluate the matrices
c-----  
        call evlmat(om,wvno,jbdrys,jbdryh,wvno2,om2)
C       call evemat(om, om2, wvno, wvno2, g,RA,RB)
        call evgmat(om, om2, wvno, wvno2, g)
c-----
c       initialize the halfspace
c-----
        do 100 i=1,5
            e(i) = gbr(1,i)
            cd(i) = e(i)
  100   continue
c-----
c       now compute the Haskell Product from top down
c-----
        do 1000 j=1,4
            do 1100 i=1,4
                aa(i,j) = dcmplx(0.0d+00, 0.0d+00)
 1100       continue
            aa(j,j) = dcmplx(1.0d+00,0.0d+00)
 1000   continue
                
        do 1340 m=1,mmax-1
c-----
c           get da(m) matrix
c-----
            call copy4(da,har,m,1,hex,ex)
            exl = exl + ex
c-----
c           form da * aa
c-----
            call dmult(da,aa)
c-----
c           copy da to aa
c-----
            call copy(da,aa)
 1340   continue
        do 1341 m=mmax-1,1,-1
c-----
c       handle the compound matrices
c-----
            call copy5(ca,dar,m,1,dex,exa)
            call cmult(e,ca,exa,exe)
 1341   continue
c-----
c       put in the halfspace condition
c-----
        e(1) = e(1)*hsr(1,1) + e(2)*hsr(1,2) + e(3)*hsr(1,3)
     1      + e(4)*hsr(1,4) + e(5)*hsr(1,5)
            call dmult(g ,aa)
c-----
c       compute the receiver function
c-----
        if(dop)then
            if(iout.eq.1)then
c-----
c           UR/UZ
c-----
            Z = - dcmplx(0.0d+00, 1.0d+00) * g(2,2) / g(2,1) 
            else if(iout.eq.2)then
c-----
c           UZ
c-----
            Z = -  g(2,1)*exp(-exl)/(g(1,1)*g(2,2) - g(2,1)*g(1,2))
            Z = Z /(cmplx(0.0,1.0)*om)
            else if(iout.eq.3)then
c-----
c           Ur
c-----
            Z =   dcmplx(0.0d+00, 1.0d+00) * g(2,2)*exp(-exl)/
     1          (g(1,1)*g(2,2) - g(2,1)*g(1,2))
            Z = Z /(cmplx(0.0,1.0)*om)
            endif
        else
            if(iout.eq.1)then
c-----
c           UR/UZ
c-----
            Z = - dcmplx(0.0d+00, 1.0d+00) * g(1,2) / g(1,1)
            else if(iout.eq.2)then
c-----
c           UZ
c-----
            Z =   dcmplx(0.0d+00, 1.0d+00) * g(1,1)*exp(-exl)/
     1          (g(1,1)*g(2,2) - g(2,1)*g(1,2))
            Z = Z /(cmplx(0.0,1.0)*om)
            else if(iout.eq.3)then
c-----
c           Ur
c-----
            Z = -                           g(1,2)*exp(-exl)/
     1          (g(1,1)*g(2,2) - g(2,1)*g(1,2))
            Z = Z /(cmplx(0.0,1.0)*om)
            endif
        endif
        return
        end

        subroutine getegn(a,b,c,d,e,f,omega2,wvno2,rp, rsv, rsh,
     1      x11,x21,x31,x41,x12,x22,x32,x42,NP, NSV,m,omg,wvn,
     2      atna,atnb)
        implicit none
        COMPLEX*16 A,B,C,D,E,F
        COMPLEX*16 wvno2, omega2
        COMPLEX*16 rsh, rp, rsv
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
        COMPLEX*16 omg, wvn
        complex*16 atna, atnb

        COMPLEX*16 xka2, xkb2

        integer NL
        parameter (NL=200)
        common/timodel/H(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      TRho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL), 
     3      frefp(NL), frefs(NL)
        real H, TA, TC, TF, TL, TN, TRho
        real qa, qb, etap, etas, frefp, frefs
c-----
c       internal variables
c-----
        COMPLEX*16 L2(2)
        COMPLEX*16 bb, cc
        COMPLEX*16 CDSQRT 
        COMPLEX*16 SRTBB4AC
        COMPLEX*16 ddef, aabc

        COMPLEX*16 ZFAC

        real mu
        complex*16 CTF
        complex*16 atna2, atnb2

c-----
c       first test to see if a fluid layer - if it is fluid, the
c       eigenfunctions are specially computed and we need only the
c       rp
c-----
        atna2= atna*atna
        atnb2= atnb*atnb
        if(TL(m).eq.0.0 .or. TN(m).eq.0.0)then
            rp = cdsqrt(wvno2 -omega2*TRho(m)/(TA(m)*atna2))
            rsv = dcmplx(0.0d+000, 0.0d+00)
            rsh = dcmplx(0.0d+000, 0.0d+00)
            return
        endif
        
c-----
c       Do the SH
c-----
        rsh = CDSQRT(TN(m)*atnb2*wvno2/(TL(m)*atnb2) 
     1      - Trho(m)*omega2/(TL(m)*atnb2))
c-----
c       I use the artifice of a mu instead of the direct TC
c       since I wish to associate the attenuation with the
c       P and S directly
c-----
        mu=(TC(m)+TA(m)+6.0*TL(m)
     1      +5.0*TN(m)-2.0*TF(m))/15.0
        CTF=0.5d+00*
     1    ((TC(m)+TA(m))*atna*atna
     2      +(6.0*TL(m)+5.0*TN(m)-15.0*mu)*atnb*atnb)
        a = wvn * CTF / (TC(m)*atna2)
        b = 1.0/(TC(m)*atna2)
        c = - TRho(m)*omg*omg + wvn*wvn *
     1      (TA(m)*atna2 -CTF*CTF/(TC(m)*atna2))
        d = - wvn
        e = 1.0/(TL(m)*atnb2)
        f = - TRho(m)*omg*omg

c-----
c       do algebra first to avoid numerical problems
c-----
        ddef = wvn*wvn - TRho(m)*omg*omg/(TL(m)*atna2)
        aabc = wvn*wvn*TA(m)/TC(m) - TRho(m)*omg*omg/(TC(m)*atna2)

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
c       get the norms - note that the true norm will be 2  NP 
c           and 2 L(2) NSV
c       The factorization permits us to use the sin nz/n or n sin nz
c-----
        NP  = (  L2(1)*(-2*a*b*d + 2*a*a*e + b*c*e - b*b*f)
     1      + (a*a+b*c)*(2*b*d*d - 2*a*d*e + b*e*f - c*e*e) )
        NSV = (- L2(2)*(2*b*d*d - 2*a*d*e - c*e*e + b*e*f)
     1      + (d*d+e*f)*(2*a*b*d - 2*a*a*e + b*b*f - b*c*e) )
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
c       Force the eigenfunctions to be as given in 5.4.4
c-----
        zfac = rp / x21
        x11  = x11 *zfac
        x21  = x21 *zfac
        x31  = x31 *zfac
        x41  = x41 *zfac

        zfac = rsv / x12
        x12  = rsv
        x22  = x22 * zfac
        x32  = x32 * zfac
        x42  = x42 * zfac
        
        np   = x11*x41 - x21*x31
        nsv  = x12*x42 - x22*x32

        return
        end
        subroutine hska(AA,tcosp,trsinp,tsinpr,tcossv,trsinsv,tsinsvr,
     1              NP,NSV,
     1      X11, X21, X31, X41,X12, X22, X32, X42, 
     2      rho,iwat,ex,om2)
c-----
c       Changes
c
c       01 May 2002  - defined cosp = tcosp/NP to 
c           reduce nmber of complex divisions
c-----
        implicit none
        complex*16 AA(4,4)
        COMPLEX*16 tcosp , tcossv 
        COMPLEX*16 trsinp, trsinsv
        COMPLEX*16 tsinpr, tsinsvr
        COMPLEX*16 NP, NSV
        COMPLEX*16 X11, X21, X31, X41
        COMPLEX*16 X12, X22, X32, X42
        real rho
        integer iwat
        real*8 ex, dfac
        complex*16  om2
        complex*16 zrho

c-----
c       introduce shorthand to reduce work
c-----
        COMPLEX*16 cosp , cossv 
        COMPLEX*16 rsinp, rsinsv
        COMPLEX*16 sinpr, sinsvr

        integer i, j
        zrho = dcmplx(dble(rho),0.0d+00)
        if(iwat.eq.1)then
c-----
c       fluid layer
c-----
            do 100 j=1,4
                do 101 i=1,4
                    AA(i,j) = dcmplx(0.0d+00,0.0d+00)
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

            AA(1,1) =   x11*x41*cosp  + x12*x42*cossv 
            AA(1,2) = - x11*x31*rsinp - x12*x32*sinsvr
            AA(1,3) = - x11*x21*cosp  - x12*x22*cossv 
            AA(1,4) =   x11*x11*rsinp + x12*x12*sinsvr
            AA(2,1) =   x21*x41*sinpr + x22*x42*rsinsv
            AA(2,2) = - x21*x31*cosp  - x22*x32*cossv 
            AA(2,3) = - x21*x21*sinpr - x22*x22*rsinsv
            AA(3,1) =   x31*x41*cosp  + x32*x42*cossv 
            AA(3,2) = - x31*x31*rsinp - x32*x32*sinsvr
            AA(4,1) =   x41*x41*sinpr + x42*x42*rsinsv
            AA(2,4) = - AA(1,3)
            AA(3,3) =   AA(2,2)
            AA(3,4) = - AA(1,2)
            AA(4,2) = - AA(3,1)
            AA(4,3) = - AA(2,1)
            AA(4,4) =   AA(1,1)
        endif
        return
        end

        subroutine hskl(cossh,rsinsh,sinshr,TL,iwat,hl,atnb)
        implicit none
        complex*16 cossh, rsinsh, sinshr
        real TL
        integer iwat
        complex*16 hl(2,2)
        complex*16 atnb
        if(iwat.eq.0)then   
            hl(1,1) = cossh
            hl(2,1) = TL*rsinsh*atnb*atnb
            hl(1,2) = sinshr/(TL*atnb*atnb)
            hl(2,2) = cossh
        else
            hl(1,1) = dcmplx(1.0d+00,0.0d+00)
            hl(1,2) = dcmplx(0.0d+00,0.0d+00)
            hl(2,1) = dcmplx(0.0d+00,0.0d+00)
            hl(2,2) = dcmplx(1.0d+00,0.0d+00)
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

        subroutine normc(e,ex,xnorm)
        implicit none

        complex*16 e(*)
        real*8 ex, xnorm

        common/lwater/lfluid
        logical lfluid
        real*8 test,testt,x,y,fac
        real*8 DREAL
        integer i, IUP
        test = 0.0D+00
        testt = 0.0D+00
        if(lfluid)then
            IUP = 2
        else
            IUP = 5
        endif
        do 2 i = 1,IUP
            if(dabs(dreal(e(i))).gt.testt) testt =dabs(dreal(e(i)))
            if(dabs(dimag(e(i))).gt.testt) testt =dabs(dimag(e(i)))
    2   continue
        if(testt.lt.1.0e-30)testt=1.0
        do 1 i =1,IUP
            x=dreal(e(i))/testt
            y=dimag(e(i))/testt
            fac = x*x + y*y
            if(test.lt.fac) test = fac
    1   continue
        test = testt*dsqrt(test)
        if(test.lt.1.0d-30) test=1.0
        xnorm = 1./test
        ex =-dlog(xnorm)
        return
        end

        subroutine rcmult(y,c,exa,exe)
        implicit none
        complex*16 y(5,5)
        complex*16 c(5,5)
        real*8 exa, exe
        common/lwater/lfluid
        logical lfluid
c-----
c       FORM YC where y(5x5) c(5x5) RETURN Y
c-----
        real*8 eval
        real*8 xnorm
        complex*16 ztmp, ee(5,5)
        integer IUP, i, j, k
        if(lfluid)then
            IUP = 2
        else
            IUP = 5
        endif
        do 1350 i=1,IUP
            do 1351 j=1,IUP
                ztmp = dcmplx(0.0d+00,0.0d+00)
                do 1349 k=1,IUP
                    ztmp=ztmp+ y(i,k) * c(k,j)
 1349           continue
                ee(i,j)=ztmp
 1351       continue
 1350   continue
        exe = exe + exa
C       call rnormc(ee,eval,xnorm)
C       do 1353 j=1,IUP
C           do 1352 i=1,IUP
C               y(i,j) = ee(i,j)*xnorm
C 1352      continue
C 1353  continue
C       exe = exe + eval
        return
        end

        subroutine rnormc(e,ex,xnorm)
        implicit none

        common/lwater/lfluid
        logical lfluid
        real*8 ex
        real*8 DREAL
        real *8 test,testt,x,y,fac,xnorm
        complex*16 e(5,5)

        integer i, j, IUP

        test = 0.0D+00
        testt = 0.0D+00
        if(lfluid)then
            IUP = 2
        else
            IUP = 5
        endif
        do 3 j=1,IUP
            do 2 i = 1,IUP
            if(dabs(dreal(e(i,j))).gt.testt) testt =dabs(dreal(e(i,j)))
            if(dabs(dimag(e(i,j))).gt.testt) testt =dabs(dimag(e(i,j)))
    2       continue
    3   continue
        if(testt.lt.1.0e-30)testt=1.0
        do 4 j=1,IUP
            do 1 i =1,IUP
                x=dreal(e(i,j))/testt
                y=dimag(e(i,j))/testt
                fac = x*x + y*y
                if(test.lt.fac) test = fac
    1       continue
    4   continue
        test = testt*dsqrt(test)
        if(test.lt.1.0d-30) test=1.0
        xnorm = 1./test
        ex =-dlog(xnorm)
        return
        end

        subroutine var(p,q,r, rp, rsv, rsh,
     1      cosp, cosq, cosr, rsinp, rsinq, rsinr,
     1      sinpr, sinqr, sinrr,pex,svex,shex,iwat)
c-----
c       p = rp  * h
c       q = rsv * h
c       r = rsh * h
c       rp  vertical wave number for P
c       rsv vertical wave number for SV
c       rsh vertical wave number for SH
c       cosp=cosh(p)  rsinp =rp *sinh(p)  sinpr = sinh(p)/rp
c       cosq=cosh(q)  rsinsv=rsv*sinh(p)  sinpq = sinh(p)/rsv
c       cosr=cosh(r)  rsinsh=rsh*sinh(p)  sinpr = sinh(p)/rsh
c-----
        implicit none
        COMPLEX*16 p, q, r
        COMPLEX*16 rp, rsv, rsh
        COMPLEX*16 cosp, cosq, cosr
        COMPLEX*16 rsinp, rsinq, rsinr
        COMPLEX*16 sinpr, sinqr, sinrr
        REAL *8 pex,svex,shex
        integer iwat

        REAL*8 pr, pi, qr, qi, rr, ri
        COMPLEX*16 epp, epm, eqp, eqm, erp, erm
        COMPLEX*16 sinp, sinq, sinr

        REAL*8 PFAC, SVFAC, SHFAC
        
        pex  = 0.0d+00
        svex = 0.0d+00
        shex = 0.0d+00
        pr = dreal(p)
        pi = dimag(p)
        qr = dreal(q)
        qi = dimag(q)
        rr = dreal(r)
        ri = dimag(r)
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
            cosr  = 1.0d+00
            rsinr = 0.0d+00
            sinrr = 0.0d+00
            shex  = 0.0d+00
        else
c-----
c       elastic layer
c-----
            svex = qr
            shex = rr
            epp = dcmplx(dcos(pi), dsin(pi))/2.0
            epm = dconjg(epp)
            eqp = dcmplx(dcos(qi), dsin(qi))/2.0
            eqm = dconjg(eqp)
            erp = dcmplx(dcos(ri), dsin(ri))/2.0
            erm = dconjg(erp)
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

            if(rr.lt.15.) then
                shfac=dexp(-2.*rr)
            else
                shfac  = 0.0d+00
            endif
            cosr = erp + shfac*erm
            sinr = erp - shfac*erm
            rsinr = rsh*sinr
            sinrr = sinr/rsh
        endif
        return
        end

