      program fbutt96
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME V                                                       c
c                                                                     c
c      PROGRAM: FBUTT96                                               c
c                                                                     c
c      COPYRIGHT 1996                                                 c
c      R. B. Herrmann                                                 c
c      Department of Earth and Atmospheric Sciences                   c
c      Saint Louis University                                         c
c      221 North Grand Boulevard                                      c
c      St. Louis, Missouri 63103                                      c
c      U. S. A.                                                       c
c                                                                     c
c---------------------------------------------------------------------c
c-----
c       11 SEP 2000 add TP, TSV TSH to FILE96 format
c       14 JAN 2001 add A, C, F, L, N and density information
c-----
c-----
c
c     hbutt96 -fl FLOW -fh FHIGH -n NORDER
c           where 1 <= NORDER <= 10
c
c     This program applies an n'th order Butterworth bandpass
c     filter, with corner frequencies FLOW and FHIGH, to
c     a time series in a file96(V) output format 
c-----
c-----
c       variables in common blocks
c
c       iftype  I*4 File type
c               1 - single trace
c               3 - three component
c               16 - Green's function
c               21 - Green's function
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
        real*4 ssec
c-----
c       command line arguments
c-----
        integer*4 norder
        real*4 fl, fh
c-----
c       variables local to the program
c-----
        parameter (LER=0, LIN=5, LOT=6)
        integer NSAMP
        parameter (NSAMP=16384)
        integer*4 npts
        common/xx/x
        real*4 x(NSAMP)
c-----
c       call machine dependent initialization
c-----
        call mchdep()
c-----
c       parse command line arguments
c-----
        call gcmdln(fl,fh,norder)
 1000   continue
            call rdhd96(LIN,nerr)
            if(nerr .lt. 0 )go to 9999
c-----
c           modify the comment string
c-----
            call wrhd96(LOT,nerr)
            ifirst = 0
            dtold = 0.0
            nptold = 0
c-----
c           initialize the pulse now that we know dt
c-----
            do 200 jgrn=1,21
c-----
c           initialize input array  
c-----
                if(jsrc(jgrn).ne.0)then
                    call rdtr96(LIN,stcomp,cmpinc,cmpaz,
     1                  cmpdt, npts, ksyear, ksmon, 
     2                  ksday, kshour, ksmin, ssec, 
     3                  x,nerr,NSAMP)
                    if(nerr .lt. 0 )go to 9999
c-----
c               check the sample interval and whether
c               the source has been defined already
c
c               if number of points, sample interval are
c               the same, do not recompute the pulse
c-----
                    if(ifirst.eq.0 .or. dtold.ne.cmpdt 
     1                  .or. npts.ne.nptold)then
c-----
c                       set up filter coefficients which 
c                       depend upon dt
c-----
                        dt = cmpdt
                        fnyq = 0.5/dt
                        if(fh.lt.0.0 .or. 
     1                      fh.gt.fnyq)fh = fnyq
                        call filcof(fl,fh,dt,norder)
                        nptold = npts
                        dtold = cmpdt
                        ifirst = 1
                    endif
c-----
c                   convolve
c-----
c-----
c                   initialize the time series
c-----
                    call rinit(x(1))
                    call filt(x,npts)
c-----
c                   output the filtered trace
c-----
                    call wrtr96(LOT,stcomp,cmpinc,cmpaz,
     1                  cmpdt, npts, ksyear, ksmon, 
     2                  ksday, kshour, ksmin, ssec, 
     3                  x,nerr,NSAMP)
                endif
  200       continue
        go to 1000
 9999   continue
        close (4)
        end

        subroutine gcmdln(fl,fh,norder)
c-----
c       parse command line arguments and return control
c       parameters
c
c       fl  R*4 - low pass corner of filter
c       fl  R*4 - high pass corner of filter
c       n   I*4 - order of filter
c                   low and high frequency asymptotes
c                   vary as f**n and f**-n, respectively
c
c-----
        character*20 name
        integer*4 nmarg
        norder = 1
        fl = 0.0
        fh = -1.0
        nmarg = mnmarg()
        i = 0
   11   i = i + 1
            if(i.gt.nmarg)goto 13
            call mgtarg(i,name)
            if(name(1:2).eq.'-n' .or. name(1:2).eq.'-N')then
                i=i+1
                call mgtarg(i,name)
                read(name,'(bn,i10)')norder
            else if(name(1:3).eq.'-fl'.or.name(1:3).eq.'-FL')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,f10.0)')fl
            else if(name(1:3).eq.'-fh'.or.name(1:3).eq.'-FH')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,f10.0)')fh
            else if(name(1:2).eq.'-?')then
                call usage(' ')
            else if(name(1:2).eq.'-h')then
                call usage(' ')
            endif
        go to 11
   13   continue
        if(norder.lt.1)norder=1
        if(norder.gt.10)norder=10
        return
        end

        subroutine usage(str)
        parameter (LER=0, LIN=5, LOT=6)
        character str*(*)
        write(LER,*)'fbutt96:',str
        write(LER,*)'USAGE: ',
     1  'fbutt96  [ -fl fl ] [ -fh fh ] [ -n norder ] [-?] [-h]'
        write(LER,*)
     1  ' -fl fl      (default 0.0) low  cut corner frequency'
        write(LER,*)
     1  ' -fh fh      (default Nyquist frequency) ',
     2      'high  cut corner frequency'
        write(LER,*)
     1  ' -n norder   (default 1) order of filter',
     2  ' 0 <= n <=10 '
        write(LER,*)
     1  ' -?           Write this help message'
        write(LER,*)
     1  ' -h           Write this help message'
        stop
        end

        subroutine rinit(x0)
c-----
c       initialize past history
c-----
c       set up recursive filter common block
        common/flt/y(3,10),x(3)
        real*8 y,x
        common/flt1/b(2,10),bfact,npole
        real*8 b,bfact
        real*4 x0
        do 1000 i=1,3
            x(i) = x0
            do 1000 j=1,npole
                y(i,j)=0.0d+00
 1000   continue
        return
        end

        subroutine filt(data,n)
c-----
c       data    R*4 - array to be filtered, same array used for
c                   input and output
c       n   I*4 - number of points in time series to be filtered
c-----
        real*4 data(n)

        common/flt/y,x
        real*8 y(3,10),x(3)
        common/flt1/b(2,10),bfact,npole
        real*8 b,bfact
        real*8 dum
        do 1000 k=1,n
            x(1) = data(k)
            y(1,1) = x(1) - x(3) - b(1,1)*y(2,1) - b(2,1)*y(3,1)
            do 100 i = 2,npole
                y(1,i)  = 
     1            y(1,i-1)-y(3,i-1)-b(1,i)*y(2,i)-b(2,i)*y(3,i)
  100       continue
c-----
c           update history
c-----
            x(3) = x(2)
            x(2) = x(1)
            do 200 i = 1,npole
                y(3,i) = y(2,i)
                y(2,i) = y(1,i)
  200       continue
            dum = bfact * y(1,npole)
            if(dabs(dum).lt.1.0d-36)dum = 0.0d+00
            data(k) = dum
 1000   continue
        return
        end

        subroutine filcof(fl,fh,dt,norder)
        common/flt1/b(2,10),bfact,npole
        real*8 b,bfact
        real*8 wdl,wdh,x,wa,s1,sr,si,xn,xd,ang
        complex*16 s2,az,cz,z1,z2,zc,xc,bfac
        COMPLEX*16 CDSQRT
        pi=3.1415926535d+00
        bfac = dcmplx(1.0d+00,0.0d+00)
        npole = norder
        wdl=pi*fl*dt
        wdh=pi*fh*dt
        xn=dcos(wdl+wdh)
        xd=dcos(wdh-wdl)
        x=xn/xd
        wa = dsin(wdh-wdl)/dcos(wdh-wdl)
        s1=dabs(wa)
        npole2 = (npole + 1 ) /2
        do 1000 i=1,npole2,1
            ang = 0.5*pi *(1.0 + (2.0*i-1.0)/npole )
            sr = -s1 * dcos(ang)
            si = -s1 * dsin(ang)
            s2 = dcmplx(sr,si)
            xc=dcmplx(x,0.0d+00)
            az=(1.0d+00,0.0d+00)+s2
            cz=(1.0d+00,0.0d+00)-s2
            zc=xc*xc-az*cz
            zc=CDSQRT(zc)
            z1=(xc+zc)/az
            z2=(xc-zc)/az
            if(i .ne. (npole+1-i))then
                b(1,i) = - 2.*dreal(z1)
                b(2,i) = CDABS(z1)**2
                b(1,npole+1-i) = -2.*dreal(z2)
                b(2,npole+1-i) = CDABS(z2)**2
                bfac = bfac * s2/(1.+s2)
                bfac = bfac * conjg(s2)/(1.0d+00+conjg(s2))
            else
                b(1,i) = -2.*x/(1.0d+00 + s1)
                b(2,i) = ( 1.0d+00 - s1)/(1.0d+00 + s1)
                bfac = bfac * s1/(1.0d+00 + s1)
            endif
 1000   continue
        do 2000 i=1,npole
 2000   continue
        bfact = dreal(bfac)
        return
        end
