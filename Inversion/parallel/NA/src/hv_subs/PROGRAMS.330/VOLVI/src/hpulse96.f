        program hpulse96
c---------------------------------------------------------------------c
c                                                                     c
c        COMPUTER PROGRAMS IN SEISMOLOGY                              c
c        VOLUME V                                                     c
c                                                                     c
c        PROGRAM: HPULSE96                                            c
c                                                                     c
c        COPYRIGHT 1996 R. B. Herrmann                                c
c                                                                     c
c        Department of Earth and Atmospheric Sciences                 c
c        Saint Louis University                                       c
c        221 North Grand Boulevard                                    c
c        St. Louis, Missouri 63103                                    c
c        U. S. A.                                                     c
c                                                                     c
c---------------------------------------------------------------------c
c       CHANGES
c       11 SEP 2000 - build in P, SV and SH first arrival times
c       14 JAN 2001 - put in A, C, F, L and N constants into
c           trace header
c       05 FEB 2004 - modified to be slightly more tolerant about DT for
c           user supplied pulse - now issues WARNING and not termination
c           mlarocca@ov.ingv.it
c       07 FEB 2005 - add a -Z flag to indicate that the 
c           internal parabolic 
c           or triangular pulses are to be zero phase 
c       12 DEC 2007 - set header values evlat, evlon, stlat, stlon to -12345
c           for compatibility of resulting SAC traces
c       23 MAR 2008 - define iprog = 4 for hudson96.f
c       23 JUL 2009 - ensured consistent call to FFT routine with
c               complex arguments. We now call zfour
c       05 JUL 2012 - changed the IO for the user defined pulse to be
c                   list directed IO instead of a specific format
c-----
c
c       hpulse96   [ -t -o -p -i ] -a alpha
c        -l L [ -D -V -A]  [-F rfile ] [ -m mult] 
c        [ -OD -OV -OA ]
c
c       This convolves the output of hspec96 with a source time function
c       to yield an output time history. 
c-----
c-----
c       variables in common blocks
c
c       iftype  I*4 File type
c               1 - single trace
c               3 - three component
c               16 - Green s function
c               21 - Green s function
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
     2      TP, TSV, TSH,
     3      SA, SC, SF, SL, SN, SR
        real*4 esec, distkm, distdg, evstaz, stevaz
        real*4 evlat, evlon, evdep 
        real*4 stlat, stlon, stelev
        REAL*4 TP, TSV, TSH
        REAL*4 SA, SC, SF, SL, SN, SR

        common/chdr96/ stname, cfilt, cpulse, ccomnt
        character*8 stname*8, cfilt*80, cpulse*80, ccomnt*80    

        character stcomp*8
        real*4 cmpaz, cmpinc, cmpdt
        integer*4 npts
        real*4 ssec
c-----
c       command line arguments
c-----
        integer*4 ntau, ipt, idva, iodva
        real*4 alp, xmult
        character*80 rfile
        logical dozero
c-----
c       variables local to the program
c-----
        parameter (LER=6, LIN=5, LOT=6)
        integer NSAMP, NFREQ
        parameter (NSAMP=16384,NFREQ=8193)
        common/xx/x
        real*4 x(NSAMP)
        common/zzx/zdata
        complex zdata(NSAMP)
        common/zzy/datas
        complex datas(NFREQ)
        character*80 ostr
        integer*4 ksrc(21)
        real*4 ar(21), ai(21)
        character*80 mname
        logical ext
        integer*4 iszrt(21)
        character*8 ost(21)
        data iszrt/1,4, 1,4,5, 1,4,5, 1,4, 1,4, 1,4,5,1,1,1,1,1,1/
        data ost/'ZDD     ', 'RDD     ', 'ZDS     ', 'RDS     ',
     1       'TDS     ', 'ZSS     ', 'RSS     ', 'TSS     ',
     2       'ZEX     ', 'REX     ', 'ZVF     ', 'RVF     ',
     3       'ZHF     ', 'RHF     ', 'THF     ', 'PEX     ',
     4       'PDD     ', 'PDS     ', 'PSS     ', 'PVF     ', 
     5       'PHF     '/ 

c-----
c       call machine dependent initialization
c-----
        call mchdep()
c-----
c       parse command line arguments
c-----
        call gcmdln(ntau,ipt,alp,idva,iodva,xmult,rfile,ostr,dozero)
c-----
c       ensure the existence of the Green s functions generated by
c       hspec96(VI)
c-----
        inquire(file='hspec96.grn',exist=ext)
        if(.not. ext)then
                write(LER,*)'Green s function file ',
     1          'hspec96.grn does not exist'
                go to 9000
        endif
        open(4,file='hspec96.grn',access='sequential',
     1      form='unformatted', status='unknown')
        rewind 4
c-----
c       open scratch temporary file
c-----
        open(unit=3,status='scratch',form='unformatted')
        rewind 3
c-----
c       get header information
c-----
        read(4,end=9999,err=9999)iprog
        read(4,end=9999,err=9999) alpha,fl,fu,dt,n,n1,n2,df,nyq2
        read(4,end=9999,err=9999)mname
C       WRITE(0,*)iprog
C       WRITE(0,*)alpha,fl,fu,dt,n,n1,n2,df,nyq2 
C       WRITE(0,*)MNAME

c-----
c       set up file96(V) header parameters
c-----
c-----
c       modify the comment string
c-----
        if(iprog.eq.1)then
            cpulse = 'hspec96  ; '//ostr
        else if(iprog.eq.2)then
            cpulse = 'hspec96p ; '//ostr
        else if(iprog.eq.3)then
            cpulse = 'hwhole96 ; '//ostr
        else if(iprog.eq.4)then
            cpulse = 'hudson96 ; '//ostr
        endif
        if(idva.eq.0)then
            iunit = 2
        else if(idva.eq.1)then
            iunit = 3
        else if(idva.eq.2)then
            iunit = 4
        endif
        if(iodva.eq.0)then
            iunit = 2
        else if(iodva.eq.1)then
            iunit = 3
        else if(iodva.eq.2)then
            iunit = 4
        endif
        junit = 11
c-----
c       default for filter
c-----
        if(ntau.le.0)then
            tau = dt
        else
            tau = ntau * dt
        endif
        if(ipt.eq.0)then
c-----
c           triangular
c-----
            call pultd(tau,dt,n,l)
            duration = 2.0*ntau*dt
        elseif(ipt.eq.1)then
c-----
c           parabolic
c-----
            call pulpd(tau,dt,n,l)
            duration = 4.0*ntau*dt
        elseif(ipt.eq.2)then
c-----
c           Ohnaka
c-----
            call pulod(alp,dt,n,l)
            duration = 0.0
        elseif(ipt.eq.3)then
c-----
c           Dirac
c-----
            call puldd(dt,n,l)
            duration = 0.0
        else if(ipt.eq.4)then
c-----
c           User defined
c-----
            call pulud(rfile,n,tau,ddt)
                if(ddt.ne.dt)then
                    write(LER,*)'Warning: RFILE dt',ddt,
     1              ' is not same as dt',dt,
     1                  ' specified for synthetics'
                endif
            duration = 0.0
c----
c       TEST FOR RICKER WAVELET
c-----
            duration = tau

        endif
        if(.not.dozero)then
            duration = 0.0
        endif
c-----
c       get complex spectra of damped pulse
c-----
        nyq = n/2 + 1
        fac = 1.0 
        dfac = exp(-alpha*dt) 
        do 200 i = 1,n 
            zdata(i) = cmplx(fac*x(i)*xmult,0.0)
            fac = fac * dfac 
  200   continue
        call zfour(zdata,n,-1,dt,df) 
        do 310 i = 1,nyq 
            freq=(i-1)*df 
            if(freq.gt.fu .or. i.lt.n1 .or. i.gt.n2)then
                zdata(i) = cmplx(0.0,0.0)
            endif
  310   continue 
        do 206 i = 1,nyq 
            datas(i) = zdata(i) 
  206   continue
c-----
c       start reading Green s functions. However, only define the pulses
c       with the first read. The data are in the following order:
c           DIST
c               SOURCE_DEPTH
c                   RECEIVER_DEPTH
c                       FREQ
c                           GRN01 .. GRN21
c-----
 1000   continue
            read(4,err=9999,end=9999)r,t0,depths,depthr,
     1          TP,TSV,TSH, 
     2          SA, SC, SF, SL, SN, SR
C           WRITE(0,*)r,t0,depths,depthr,
C     1         TP,TSV,TSH, 
C     1         SA, SC, SF, SL, SN, SR
            if(r.lt.0.0) go to 9999 
            read(4,end=9999,err=9999)ksrc
C           WRITE(0,*)ksrc
c-----
c           form list of output Green functions
c-----
            do 199 i=1,21
                if(ksrc((i)).ne.0)then
                    jsrc(i) = iszrt(i)
                else
                    jsrc(i) = 0
                endif
  199       continue
            iftype = 21
            iobsyn = 2
            itmfrq = 1
            cfilt = 'None'
            keyear = 0
            kemon = 0
            keday = 0
            kehour = 0
            kemin = 0
            esec = 0.0
            evlat = -12345
            evlon = -12345
            evdep = depths

            stname = 'GRN21'
            stlat  = -12345
            stlon = -12345
            stelev = depthr
            distkm = r
            distdg = r/111.195
            evstaz = 0.0
            stevaz = 180.0
            lmnm = lgstr(mname)
            ccomnt = mname(1:lmnm)
c-----
c           Output header
c-----
            call wrhd96(LOT,nerr)
c-----
c           get Green's functions for this distance, 
c           which will have to be demultiplexed
c-----
            rewind 3
            do 300 i=n1,n2
                do 305 j=1,21
                    ar(j)=0.0
                    ai(j)=0.0
                    if(ksrc(j).ne.0)then
                        read(4,err=9999,end=9999)
     1                      ar(j),ai(j)
                    endif
  305           continue
                write(3)ar,ai
  300       continue
        
c-----
c           now work with unit 3 to reorder 
c           everything by Green's function
c-----
            do 1300 jk=1,21
                if(jsrc(jk).ne.0)then
                    do 1301 i=1,n
                        zdata(i)=cmplx(0.0,0.0)
 1301               continue
                    rewind 3
                    do 1302 i=n1,n2
                        read(3)ar,ai
                        zdata(i)=cmplx(ar(jk),ai(jk))
     1                       * datas(i)
c-----
c                   integrate
c-----
                        if(idva.eq.0.and.jk.lt.16)then
                            omega = 6.2831853*(i-1)
     1                         *df      
                            zdata(i) = zdata(i) / 
     1                         cmplx(alpha,omega)
c-----
c               differentiate
c-----
                        elseif(idva.eq.2.and.jk.lt.16)then
                            omega = 6.2831853*(i-1)
     1                         *df
                            zdata(i) = zdata(i) * 
     1                         cmplx(alpha,omega)
                        endif
c-----
c               pressure field in fluid
c-----
                        if(jk.ge.16.and.iodva.lt.0)then
                            omega = 6.2831853*(i-1)
     1                         *df
                            zdata(i) = -zdata(i) 
                        else if(jk.ge.16.and.iodva.ge.0)then
                            omega = 6.2831853*(i-1)
     1                         *df
                            zdata(i) = - zdata(i) * 
     1                         cmplx(alpha,omega)
                        endif
c-----
c               make a real time series
c-----
                        if(i.ne.1)then
                         zdata(n+2-i)=conjg(zdata(i))
                        endif
 1302               continue
                    zdata(nyq)=cmplx(real(zdata(nyq)),0.0)
                    call zfour(zdata,n,+1,dt,df)
c-----
c               correct for damping 
c-----
                    fac = exp(alpha*t0) 
                    dfac = exp(alpha*dt) 
                    do 425 i = 1,n 
                        zdata(i)= zdata(i) * fac
                        fac = fac * dfac 
                        x(i) = real(zdata(i))
  425               continue 
                    if(iszrt(jk).eq.1)then
c-----
c               vertical
c-----
                        cmpinc = -90.0
                        cmpaz  =   0.0
                    else if(iszrt(jk).eq.4)then
c-----
c               radial
c-----
                        cmpinc = 0.0
                        cmpaz  =   0.0
                    else if(iszrt(jk).eq.5)then
c-----
c               transverse
c-----
                        cmpinc = 0.0
                        cmpaz  =  90.0
                    endif
                npts = n
                stcomp = ost(jk)
                cmpdt = dt
                ksyear = 0
                ksmon = 0
                ksday = 0
                kshour = 0
                ksmin = 0
                ssec = t0 - 0.5*duration
                call wrtr96(LOT,stcomp,cmpinc,cmpaz,
     1              cmpdt, npts, ksyear, ksmon, 
     2              ksday, kshour, ksmin, ssec, 
     3              x,nerr,NSAMP)
                endif
 1300       continue
        go to 1000
 9999   continue
        close (4)
        close (3)
 9000   continue
        end

        subroutine gcmdln(ntau,ipt,alp,idva,iodva,xmult,
     1      rfile,ostr,dozero)
c-----
c       parse command line arguments
c       requires subroutine mgtarg() and function mnmarg()
c-----
c
c       ntau    I*4 - pulse duration factor for parabolic and triangular
c               pulses
c       ipt I*4 - pulse type
c               0 - triangular
c               1 - parabolic
c               2 - Ohnaka
c               3 - Dirac delta function
c               4 = user pulse in file rfile
c       alp R*4 - Ohnaka pulse shape factor, fc= alp / 2 pi
c       idva    I*4 - time history type
c               0 - displacement
c               1 - velocity
c               2 - acceleration
c       iodva   I*4 - force output time history type
c               Useful if internal pulse is not to reprsent a
c               step in fault displacement, explosion pressure,
c               or force
c               0 - displacement
c               1 - velocity
c               2 - acceleration
c       xmult   R*4  - moment scaling factor
c       rfile   C*80 - name of user provided pulse 
c       ostr    C*80 - reproduction of command line
c       dozero  L   .true. for triangular and parabolic pulses, center
c               at zero lag
c-----
        character*(*)  rfile
        integer*4 ntau, ipt, idva
        real*4 xmult, alp
        character*80 ostr
        logical dozero

        integer*4 mnmarg
        character*80 name
c-----
c       initialization
c-----
        ntau = -1
        ipt = -1
        alp = -1.0
        idva = 1
        iodva = -1
        rfile = ' '
        xmult = 1.0
        dozero = .false.
        nmarg=mnmarg()
        i = 0
   11   i = i + 1
            if(i.gt.nmarg)go to 13
            call mgtarg(i,name)
            if(name(1:2).eq.'-l')then
                i=i+1
                call mgtarg(i,name)
                read(name,'(bn,i10)')ntau
            else if(name(1:2).eq.'-a')then
                i=i+1
                call mgtarg(i,name)
                read(name,'(bn,f10.0)')alp
            else if(name(1:2).eq.'-t')then
                ipt = 0
            else if(name(1:2).eq.'-p')then
                ipt = 1
            else if(name(1:2).eq.'-o')then
                ipt = 2
            else if(name(1:2).eq.'-i')then
                ipt = 3
            else if(name(1:2).eq.'-F')then
                ipt = 4
                        i=i+1
                        call mgtarg(i,rfile)
            else if(name(1:2).eq.'-D')then
                idva = 0
            else if(name(1:2).eq.'-V')then
                idva = 1
            else if(name(1:2).eq.'-A' .and. name(1:4).ne.'-ALL')then
                idva = 2
            else if(name(1:3).eq.'-OD')then
                iodva = 0
            else if(name(1:3).eq.'-OV')then
                iodva = 1
            else if(name(1:3).eq.'-OA')then
                iodva = 2
            else if(name(1:2).eq.'-m')then
                i=i+1
                call mgtarg(i,name)
                read(name,'(bn,f20.0)')xmult
            else if(name(1:2).eq.'-Z')then
                dozero = .true.
            else if(name(1:2).eq.'-?')then
                call usage('Help')
            else if(name(1:2).eq.'-h')then
                call usage('Help')
            endif
        go to 11
   13   continue
c-----
c     do some error checking, e.g., we cannot generate velocity
c     for triangular pulse
c-----
        if(ipt.ge.0 .and. ipt .le.1 .and . ntau.le.0)ntau = 1
        if(ipt.eq.2 .and. alp .le.0.0)alp = 1.0
        if(ipt.eq.2 .and. alp.lt.0.0)
     1      call usage('No alpha for Ohnaka pulse')
        if(ipt.lt.0)
     1      call usage('No pulse shape defined')
        if(iodva.lt.0)iodva = idva
        ostr = 'hpulse96'
        do 14 i=1,nmarg
            call getarg(i,name)
            L = lgstr(name)
            lo = lgstr(ostr)
            ostr = ostr(1:lo)//' '//name(1:L)
   14   continue
        return
        end

        subroutine usage(str)
        parameter (LER=6, LIN=5, LOT=6)
        character str*(*)
        write(LER,*)'hpulse96:',str
        write(LER,*)'USAGE: ',
     1  'hpulse96   [ -t -o -p -i ] -a alpha',
     2  ' -l L [ -D -V -A]  [-F rfile ] [ -m mult] ',
     3  ' [ -OD -OV -OA ] [-Z] [-?] [-h]'
        write(LER,*)
     1  ' -t           Triangular pulse of base 2 L dt'
        write(LER,*)
     1  ' -p           Parabolic Pulse of base  4 L dt'
        write(LER,*)
     1  ' -o           Ohnaka pulse with parameter alpha'
        write(LER,*)
     1  ' -i           Dirac Delta function'
        write(LER,*)
     1  ' -l L         Shape duration parameter for triangle ',
     2      'and parabolic pulses' 
        write(LER,*)
     1  ' -a alpha     Shape parameter for Ohnaka pulse'
        write(LER,*)
     1  ' -D           Output is ground displacement'
        write(LER,*)
     1  ' -V           Output is ground velocity (default)'
        write(LER,*)
     1  ' -A           Output is ground acceleration'
        write(LER,*)
     1  ' -F rfile     User supplied pulse'
        write(LER,*)
     1  ' -m mult      Multiplier (default 1.0)'
        write(LER,*)
     1  ' -OD           Output is ground displacement'
        write(LER,*)
     1  ' -OV           Output is ground velocity'
        write(LER,*)
     1  ' -OA           Output is ground acceleration'
        write(LER,*)
     1  ' -Z           (default false) zero phase ',
     2  'triangular/parabolic pulse'
        write(LER,*)
     1  ' -?           Write this help message'
        write(LER,*)
     1  ' -h           Write this help message'
        stop
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

        subroutine pulpd(tau,dt,nt,l)
        integer NSAMP
        parameter (NSAMP=16384)
        common/xx/x
        real*4 x(NSAMP)
c-----
c       far field displacement of parabolic pulse of unit area
c       
c       tau R*4     - duration parameter of pulse
c                   pulse duration is 4*tau seconds
c       dt  R*4 - sample interval in seconds
c       nt  I*4 - number of points in time series
c       l   R*4 - non-zero length of pulse
c-----
        ltest = 0
        TL = tau
        t1 = 0.01*dt
        t2 = t1 + tau
        t3 = t2 + tau
        t4 = t3 + tau
        t5 = t4 + tau
        do 100 i = 1,nt
            t=(i-1)*dt
            z = t - t1
            x(i) = 0.0
            if(t.ge.t1 .and. t.lt.t2)then
                x(i) = 0.5*(z/TL)**2
            else if(t.ge.t2 .and. t.lt.t3)then
                x(i)= -0.5*(z/TL)**2 + 2.*(z/TL) - 1
            else if(t.ge.t3 .and. t.lt.t4)then
                x(i)= -0.5*(z/TL)**2 + 2.*(z/TL) - 1
            else if(t.ge.t4 .and. t.lt.t5)then
                x(i)= 0.5*(z/TL)**2 - 4.*(z/TL) + 8.
            else if(t.ge.t5)then
                ltest = ltest + 1
                if(ltest.eq.1) l = i
            endif
  100   continue
c-----
c       pulse normalized so first integral has area of unity
c-----
        do 200 i = 1,nt
            x(i) = x(i)/(2.*TL)
  200   continue
        return
        end

        subroutine pulod(alp,dt,nt,l)
c-----
c       far field displacement of Ohnaka slip history pulse of unit area
c       Harkrider (1976) Geophys J. 47, p 97.
c       
c       alp R*4     - shape parameter
c                   pulse duration is 4*tau seconds
c       dt  R*4 - sample interval in seconds
c       nt  I*4 - number of points in time series
c       l   R*4 - non-zero length of pulse
c-----
        integer NSAMP
        parameter (NSAMP=16384)
        common/xx/x
        real*4 x(NSAMP)
        ltest = 0
        do 100 i=1,nt
            t=(i-1)*dt
            x(i)=0.0
            arg= alp*t
            al2=alp*alp
            if(arg.le.25.)then
                x(i)= al2*t*exp(-arg)
            else
                ltest = ltest +1
                if(ltest.eq.1)l = i
            endif
  100   continue
        return
        end

        subroutine pultd(tau,dt,nt,l)
c-----
c       far field displacement of triangular pulse of unit area
c       
c       tau R*4     - duration parameter of pulse
c                   pulse duration is 4*tau seconds
c       dt  R*4 - sample interval in seconds
c       nt  I*4 - number of points in time series
c       l   R*4 - non-zero length of pulse
c-----
        integer NSAMP
        parameter (NSAMP=16384)
        common/xx/x
        real*4 x(NSAMP)
c-----
c       note re really are working with a triangular pulse
c-----
        ltest = 0
        dt1=tau
        dt2=0.0
        dt3=tau
        t1=dt1
        t2=t1+dt2
        t3=t2+dt3
        fac = 0.5*dt1 + dt2 + 0.5*dt3
        fac = 1./fac
        v1 = fac/dt1
        v2 = 0.0
        v3 = - fac/dt3
        do 100 i=1,nt
        t = (i-1)*dt
        x(i)=0.0
        if(t.le.t1)then
            delt = t - 0.0
            x(i) = delt*fac/dt1
        else if(t.gt.t1.and.t.le.t2)then
            delt = t - t1
            x(i) = fac
        else if(t.gt.t2.and.t.le.t3)then
            delt = t - t2
            x(i)= fac + v3*delt
        else if(t.gt.t3)then
            x(i)=0.0
            ltest = ltest + 1
            if(ltest.eq.1)l = i
        endif
  100   continue
        return
        end

        subroutine puldd(dt,n,l)
c-----
c       Dirac Impulse
c
c       dt  R*4 - sampling interval in seconds
c       n   I*4 - number of samples
c       l   I*4 - duration of pulse , obviously 1 here
c-----
        integer NSAMP
        parameter (NSAMP=16384)
        common/xx/x
        real*4 x(NSAMP)
c-----
c       Dirac Delta Pulse
c-----
        do 100 i=1,n
            x(i) = 0.0
  100   continue
        l = 1
        x(1) = 1.0/dt
        return
        end

        subroutine pulud(rfile,n,tau,dtt)
c-----
c       User defined pulse - hopefully one with unit area
c
c       rfile   C*(*)   - name of user provided pulse
c       n   I*4 - number of points in time series
c       tau R*4 - duration measure
c-----
        character rfile*(*)
        integer NSAMP
        parameter (NSAMP=16384)
        common/xx/x
        real*4 x(NSAMP)
        do 100 i=1,n
            x(i) = 0.0
  100   continue
        open(1,file=rfile,status='unknown',form='formatted',
     1      access='sequential')
        rewind 1
        read(1,*)np,dtt
        read(1,*)(x(i),i=1,np)
   10   format(4e15.7)
        close(1)
        tau = (np-1)*dtt
        return
        end
