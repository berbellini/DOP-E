      program finteg96
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME V                                                       c
c                                                                     c
c      PROGRAM: FINTEG96                                              c
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
c     finteg96
c
c     This program integrates 
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
        logical dodc
c-----
c       variables local to the program
c-----
        parameter (LER=0, LIN=5, LOT=6)
        integer NSAMP
        parameter (NSAMP=16384)
        integer*4 npts
        real*4 x(NSAMP)
c-----
c       call machine dependent initialization
c-----
        call mchdep()
c-----
c       parse command line arguments
c-----
        call gcmdln(dodc)
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
                    call integ(x,npts,dt,dodc)
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

        subroutine gcmdln(dodc)
c-----
c       parse command line arguments and return control
c       parameters
c
c       dodc    L   - .true. subtrace x(1) before integrating
c
c-----
        logical dodc
        character*20 name
        integer*4 nmarg
        nmarg = mnmarg()
        dodc = .false.
        i = 0
   11   i = i + 1
            if(i.gt.nmarg)goto 13
            call mgtarg(i,name)
            if(name(1:3).eq.'-DC')then
                dodc = .true.
            else if(name(1:2).eq.'-?')then
                call usage(' ')
            else if(name(1:2).eq.'-h')then
                call usage(' ')
            endif
        go to 11
   13   continue
        return
        end

        subroutine usage(str)
        parameter (LER=0, LIN=5, LOT=6)
        character str*(*)
        write(LER,*)'finteg96:',str
        write(LER,*)'USAGE: ',
     1  'finteg96  [-DC] [-?] [-h]'
        write(LER,*)
     1  ' -DC  (default false) subtract first point of time series',
     2  ' before integrating'
        write(LER,*)
     1  ' -?           Write this help message'
        write(LER,*)
     1  ' -h           Write this help message'
        stop
        end

        subroutine integ(x,nt,dt,dodc)
        parameter (LER=0, LIN=5, LOT=6)
c-----
c       Form the integral x(t) dt numerically
c       x   R*4 - array to be integrated
c       nt  I*4 - number of data points in time series
c       dt  R*4 - sampling interval
c       dodc    L*4 - .true. remve x(1) before summing
c-----
        
        real x(nt)
        integer*4 nt
        real*4 dt
        logical dodc
        sum=0.0
        if(dodc)then
            do 100 k = 1,nt
                sum = sum + dt*(x(k) - x(1))
                x(k)=sum
  100       continue
        else
            do 200 k = 1,nt
                sum = sum + dt*x(k)
                x(k)=sum
  200       continue
        endif
        return
        end
