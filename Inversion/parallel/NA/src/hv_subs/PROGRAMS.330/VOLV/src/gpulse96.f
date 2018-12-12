        program gpulse96
c---------------------------------------------------------------------c
c                                                                     c
c        COMPUTER PROGRAMS IN SEISMOLOGY                              c
c        VOLUME V                                                     c
c                                                                     c
c        PROGRAM: GPULSE96                                            c
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
c       NOTE this is a no-operation at this time
c       05 FEB 2004 - modified to be slightly more tolerant about DT for
c           user supplied pulse - now issues WARNING and not termination
c           mlarocca@ov.ingv.it
c       09 JAN 2005 - LER not defined in subroutine  source - c
c                     baker@usgs.gov
c       07 FEB 2005 - add a -Z flag to indicate that the c
c                     internal parabolic 
c           or triangular pulses are to be zero phase 
c-----
c
c       gpulse96  [ -v  ] [ -t -o -p -i ] -a alpha
c        -l L [ -D -V -A]  [-F rfile ] [ -m mult] 
c        [ -OD -OV -OA ]
c
c       This convolves the output of genray96 with a 
c       source time function
c       to yield an output time history. The moment tensor Green s
c       functions from genray96 must be convolved with the
c       second derivative of the source time function, while the
c       first derivatives must be convolved with the first derivative
c       of the Green s function. This is mathematically done through
c       the source pulse
c-----
c-----
c       variables in common blocks
c
c       iftype  I*4 File type
c               1 - single trace
c               3 - three component
c               16 - Green s function
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
        real*4 ssec
c-----
c       command line arguments
c-----
        logical ext
        logical lverby
        integer*4 ntau, ipt, idva, iodva
        real*4 alp, xmult
        character*80 rfile
c-----
c       variables local to the program
c-----
        integer LOT, LIN, LER
        parameter (LER=0, LIN=5, LOT=6)
        integer NSAMP
        parameter (NSAMP=2048)
        integer*4 npts
        real*4 x(NSAMP), y(NSAMP)
        character*80 ostr
        logical dozero
c-----
c       call machine dependent initialization
c-----
        call mchdep()
c-----
c       parse command line arguments
c-----
        call gcmdln(ntau,ipt,lverby,alp,idva,iodva,xmult,rfile,
     1      ostr,dozero)
c-----
c       ensure the existence of the Green s functions generated by
c       genray96(V)
c-----
        inquire(file='genray96.grn',exist=ext)
        if(.not. ext)then
                write(LER,*)'Green s function file ',
     1          'genray96.grn does not exist'
                go to 9000
        endif
        open(4,file='genray96.grn',access='sequential',form='formatted',
     1          status='unknown')
        rewind 4

c-----
c       start reading Green s functions. However, only define the pulses
c       with the first read
c-----
        ifirst = 0
 1000   continue
        call rdhd96(4,nerr)
        if(nerr .lt. 0 )then
            go to 9999
        endif
c-----
c       modify the comment string
c-----
        cpulse = 'genray96 ; '//ostr
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
        call wrhd96(LOT,nerr)
            ifirst = 0
            dtold = 0.0
            nptold = 0
c-----
c       initialize the pulse now that we know dt
c-----
            do 200 jgrn=1,16
c-----
c           initialize input array  
c-----
                if(jsrc(jgrn).ne.0)then
                    call rdtr96(4,stcomp,cmpinc,cmpaz,
     1                  cmpdt, npts, ksyear, ksmon, 
     2                  ksday, kshour, ksmin, ssec, 
     3                  x,nerr,NSAMP)
                    if(nerr .lt. 0 )then
                        go to 9999
                    endif
c-----
c               check the sample interval and whether
c               the source has been defined already
c
c               if number of points, sample interval are
c               the same, do not recompute the pulse
c-----
                    if(ifirst.eq.0 .or. dtold.ne.cmpdt 
     1                  .or. npts.ne.nptold)then
                    call source(ntau,ipt,alp,cmpdt,
     1                  npts,rfile,y,l,dozero,
     2                  duration)
                        nptold = npts
                        dtold = cmpdt
                        ifirst = 1
                    endif
c-----
c                   correct for DC offset
c-----
                    xx=x(1)
                    do 201 i=1,npts
                        if(i.le.npts)then
                            x(i)=x(i)-xx
                        else
                            x(i) = 0.0
                        endif
  201               continue
c-----
c                   convolve
c-----
                    call filt(x,y,npts,cmp
     1                  dt,l,xmult,idva,jgrn)
                    call wrtr96(LOT,stcomp,cmpinc,cmpaz,
     1                  cmpdt, npts, ksyear, ksmon, 
     2                  ksday, kshour, ksmin, 
     3                  ssec - 0.5*duration, 
     3                  x,nerr,NSAMP)
                endif
  200       continue
        go to 1000
 9999   continue
        close (4)
 9000   continue
        end

        subroutine filt(x,y,nt,dt,lsamp,mult,idva,jgrn)
c-----
c       time domain convolution
c-----
c       y   R*4 - input time series (source pulse of lsamp samples)
c                   y = source function
c       x   R*4 - output time series derived from x = x*y
c                   x = greens response
c       nt  I*4 - number of samples
c       dt  R*4 - sampling interval
c       lsamp   I*4 - length of source wavelet in samples
c       mult    R*4 - scalar multiplier 
c       idva    I*4 - time history type
c               0 - displacement
c               1 - velocity
c               2 - acceleration
c       jgrn    I*4 - Green s function to be evaluated
c----- 
        real*4 x(*), y(*) 
        integer NSAMP
        parameter (NSAMP=2048)
        common /srce/ datas(NSAMP), dat(NSAMP)
        real*4 mult
        do 100 i=1,nt
            datas(i) = y(i)
  100   continue
c-----
c       take the derivatives of the source pulse
c-----
        if(idva.eq.1)then
c-----
c           ground velocity
c-----
            call deriv(datas,dat,nt,dt)
            lsamp = lsamp + 1
        else if(idva.eq.2)then
c-----
c           ground acceleration
c-----
            call deriv(datas,dat,nt,dt)
            lsamp = lsamp + 1
            call deriv(datas,dat,nt,dt)
            lsamp = lsamp + 1
        endif
c-----
c       one additional derivative for couple and dipole sources
c-----
        if(jgrn.lt.11 .or. jgrn.gt.15)then
            call deriv(datas,dat,nt,dt)
            lsamp = lsamp + 1
        endif
c-----
c       perform the convolution
c
c       note that y(lsamp+k) = 0 for k >= lsamp
c-----
        do 300 k = 1,nt
            imin = k - lsamp + 1
            if(imin.lt.1) imin = 1
            sum = 0.0
            do 200 i = imin,k
                sum = sum + x(i) * datas(k-i+1)
  200       continue
            dat(k)=sum*dt*mult
  300   continue
        do 400 i=1,nt
            x(i) = dat(i)
  400   continue
        return
        end

c-----
c       PULSE DEFINITIONS
c-----

        subroutine source(ntau,ipt,alp,dt,npts,rfile,x,l,dozero,
     1      duration)
c-----
c       define the source time function
c-----
c       ntau    I*4 - pulse duration factor for parabolic and triangular
c               pulses
c       ipt I*4 - pulse type
c               0 - triangular
c               1 - parabolic
c               2 - Ohnaka
c               3 - Dirac delta function
c               4 = user pulse in file rfile
c       alp R*4 - Ohnaka pulse shape factor, fc= alp / 2 pi
c       dt  R*4 - sample rate for time series
c       npts    I*4 - number of samples in time series
c       rfile   C*80- name of user provided pulse 
c       x() R*4 - time series array
c       l   I*4 - non-zero length of pulse c
c                     (used to speed up convolution)
c       dozero  L   - .true. use zero phase triangular/parabolic pulse
c       duration R*4 - duration of triangular/parabolic pulse
c-----
        integer LOT, LIN, LER
        parameter (LER=0, LIN=5, LOT=6)
        integer*4 ntau, ipt, npts,l
        real*4 alp, dt, duration
        character*80 rfile
        logical dozero
        parameter (NSAMP=2048)
        real*4 x(NSAMP)
c-----
c       initialize
c-----
        do 50 i=1,NSAMP
            x(i)=0.0
   50   continue
        if(ipt.eq.0)then
c-----
c           triangular
c-----
            tau = ntau * dt
            nt = 4*ntau + 1
            call pultd(tau,dt,nt,l,x)
            duration = 4*ntau*dt
        else if(ipt.eq.1)then
c-----
c           parabolic
c-----
            tau = ntau * dt
            nt = 8*ntau + 1
            call pulpd(tau,dt,nt,l,x)
            duration = 8*ntau*dt
        else if(ipt.eq.2)then
c-----
c           Ohnaka
c-----
            if(npts.gt.NSAMP)then
                nt=NSAMP
            else
                nt = npts
            endif
            call  pulod(alp,dt,nt,l,x)
            duration = 0.0
        else if(ipt.eq.3)then
c-----
c           Dirac
c-----
            if(npts.gt.NSAMP)then
                nt=NSAMP
            else
                nt = npts
            endif
            call puldd(dt,nt,l,x)
            duration = 0.0
        else if(ipt.eq.4)then
c-----
c           user defined
c-----
            if(npts.gt.NSAMP)then
                nt=NSAMP
            else
                nt = npts
            endif
            call pulud(rfile,nt,tau,ddt,l,x)
                if(ddt.ne.dt)then
                    write(LER,*)'Warning: RFILE dt',ddt,
     1              ' is not same as dt',dt,
     1                  ' specified for synthetic'
                endif
            ntau = tau/dt
            duration = 0.0
        endif
        if(.not. dozero)then
            duration = 0.0
        endif
        return
        end

        subroutine pulpd(tau,dt,nt,l,x)
c-----
c       unit area far field displacement parabolic pulse
c-----
c       tau R*4 - duration parameter, total duration = 4 tau
c       dt  R*4 - sample rate
c       nt  I*4 - number of points for time series
c       l   I*4 - pulse x(i) = 0 for i >=l
c       x   R*4 - time series array
c-----
        real*4 tau, dt
        integer*4 nt, l
        integer NSAMP
        parameter (NSAMP=2048)
        real*4 x(NSAMP)
        ltest = 0
        tl = tau
        t1 = 0.0
        t2 = t1 + tau
        t3 = t2 + tau
        t4 = t3 + tau
        t5 = t4 + tau
        do 100 i = 1,nt
            y=(i-1)*dt
            z = y - t1
            x(i) = 0.0
            if(y.ge.t1 .and. y.lt.t2)then
                x(i) = 0.5*(z/tl)**2
            else if(y.ge.t2 .and. y.lt.t3)then
                x(i)= -0.5*(z/tl)**2 + 2.*(z/tl) - 1
            else if(y.ge.t3 .and. y.lt.t4)then
                x(i)= -0.5*(z/tl)**2 + 2.*(z/tl) - 1.
            else if(y.ge.t4 .and. y.lt.t5)then
                x(i)= 0.5*(z/tl)**2 - 4.*(z/tl) + 8.
            else
                ltest = ltest + 1
                if(ltest.eq.1) l = i
            endif
  100   continue
c-----
c       pulse normalized so first integral has area of unity
c-----
        do 200 i = 1,nt
            x(i) = x(i)/(2.*tl)
  200   continue
        return
        end

        subroutine pulod(alp,dt,nt,l,x)
c-----
c       unit area far field displacement Ohnaka pulse
c           Harkrider (1976) Geophys J. 47, p 97.
c-----
c       alp R*4 - shape parameter, corner frequency
c               fc = alp / 2 pi
c       dt  R*4 - sample rate
c       nt  I*4 - number of points in time series
c       l   I*4 - pulse x(i) = 0 for i >=l
c       x   R*4 - time series array
c-----
        real*4 alp, dt
        integer*4 nt, l
        integer NSAMP
        parameter (NSAMP=2048)
        real*4 x(NSAMP)
        ltest = 0
        al2=alp*alp
        do 100 i=1,nt
            t=(i-1)*dt
            x(i)=0.0
            arg= alp*t
            if(arg.le.25.0)then
                x(i)= al2*t*exp(-arg)
            else
                ltest = ltest +1
                if(ltest.eq.1)l = i
            endif
  100   continue
        if(ltest.eq.0)l=i
        return
        end

        subroutine pultd(tau,dt,nt,l,x)
c-----
c       unit area far field displacement triangular pulse
c-----
c       tau R*4 - duration parameter, total duration = 2 tau
c       dt  R*4 - sample rate
c       nt  I*4 - number of points for time series
c       l   I*4 - pulse x(i) = 0 for i >=l
c       x   R*4 - time series array
c-----
        real*4 tau, dt
        integer*4 nt, l
        integer NSAMP
        parameter (NSAMP=2048)
        real*4 x(NSAMP)
        ltest = 0
        fac = 1./tau
        t1 = tau
        t2 = tau + tau
        do 100 i=1,nt
            t = (i-1)*dt
            x(i)=0.0
            if(t.le.t1)then
                z = (t - 0.0)/tau
                x(i) = z*fac
            elseif(t.gt.t1.and.t.le.t2)then
                z = (t - t1)/tau
                x(i)= fac*(1.0 - z)
            elseif(t.gt.t2)then
                x(i)=0.0
                ltest = ltest + 1
                if(ltest.eq.1)l = i
            endif
  100   continue
        return
        end

        subroutine puldd(dt,n,l,x)
        integer NSAMP
        parameter (NSAMP=2048)
        real*4 x(NSAMP)
c-----
c       Dirac Delta Pulse
c-----
        do 100 i=1,n
            x(i) = 0.0
  100   continue
        x(1) = 1.0/dt
        l = 2
        return
        end

        subroutine pulud(rfile,n,tau,dtt,l,x)
        character rfile*(*)
        integer NSAMP
        parameter (NSAMP=2048)
        real*4 x(NSAMP)
        do 100 i=1,n
            x(i) = 0.0
  100   continue
        open(3,file=rfile,status='unknown',form='formatted',
     1      access='sequential')
        rewind 3
        read(3,'(i5,f10.3)')np,dtt
        if(np.gt.n)np = n
        if(np.gt.NSAMP)np = NSAMP
        read(3,10)(x(i),i=1,np)
   10   format(4e15.7)
        close(3)
        tau = np*dtt
        l = np
        return
        end

        subroutine gcmdln(ntau,ipt,lverby,alp,idva,iodva,xmult,
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
c       lverby  L   - true output additional information
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
        logical lverby
        character*80 ostr
        logical dozero

        integer*4 mnmarg
        character*80 name
c-----
c       initialization
c-----
        ntau = -1
        ipt = -1
        lverby = .false.
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
            else if(name(1:2).eq.'-v')then
                lverby= .true.
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
        ostr = 'gpulse96'
        do 14 i=1,nmarg
            call getarg(i,name)
            L = lgstr(name)
            lo = lgstr(ostr)
            ostr = ostr(1:lo)//' '//name(1:L)
   14   continue
        return
        end

        subroutine usage(str)
        parameter (LER=0, LIN=5, LOT=6)
        character str*(*)
        write(LER,*)'gpulse96:',str
        write(LER,*)'USAGE: ',
     1  'gpulse96  [ -v  ] [ -t -o -p -i ] -a alpha',
     2  ' -l L [ -D -V -A]  [-F rfile ] [ -m mult] ',
     3  ' [ -OD -OV -OA ] [-Z] [-?] [-h]'
        write(LER,*)
     1  ' -v           Verbose output'
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

        subroutine deriv(y,x,n,dt)
c-----
c       Use centered difference to take derivative of time series
c-----
c       y   R*4 - time series to be differentiated and returned as y
c       x   R*4 - temporary storage
c       n   I*4 - length of time series
c       dt  R*4 - sample interval
c-----
        dimension y(n), x(n)
        b0=1.0/dt
        xm1=y(1)
        xm2=0.
        y1=0.
        do 11 i=1,n
            y2=b0*(y(i)-xm1) 
            xm2=xm1
            xm1=y(i)
            y(i)=y2
   11   continue
        return
        end

