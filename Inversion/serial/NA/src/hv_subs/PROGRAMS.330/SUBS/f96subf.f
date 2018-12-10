c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME V                                                       c
c                                                                     c
c      SUBROUTINES:                                                   c
c                                                                     c
c      COPYRIGHT 1996 R. B. Herrmann                                  c
c                                                                     c
c      Department of Earth and Atmospheric Sciences                   c
c      Saint Louis University                                         c
c      221 North Grand Boulevard                                      c
c      St. Louis, Missouri 63103                                      c
c      U. S. A.                                                       c
c                                                                     c
c---------------------------------------------------------------------c
c       06 SEP 2000 - add TP, TSV and TSH travel times
c       also upgrade filetype to FILE01.02 from FILE01.01
c           however recognize old FILE01.01 by
c           setting, TP, TSV and TSH = -12345.
c       16 DEC 2000 - changed format for cmpdt from f10.4 
c           to g10.4 in wrtr96
c       02 MAY 2002 - corrected the read/write order of 
c           sa,sc,sf,sl,sn,sr
c       15 DEC 2006 - extended the format for evlat, evlon, stlat, stlon
c           so that the non-value -12345.0 can be used in the manner of SAC
c       16 MAR 2007 - extended the format for evdep
c           so that the non-value -12345.0 can be used in the manner of SAC
c-----
c-----
c       This is a set of subroutines to performs file96
c       time series I/O
c
c       The routines are
c
c       rdhd96()    - read  a file96 header
c       wrhd96()    - write a file96 header
c       rdtr96()    - read a file96 trace
c       wrtr96()    - read a file96 trace
c
c       The header contains the following:
c
c       FILE_TYPE1      FILE21 FILE16 FILE03 FILE01
c       FILE_TYPE2      OBSERVED SYNTHETIC
c       FILE_TYPE3      TIME_DOMAIN FREQUENCY_DOMAIN
c       UNITS - TIME SERIES          
c       FILTERING COMMENT
c       ORIGIN_TIME     LATITUDE        LONGITUDE       DEPTH
c       STATION_NAME    
c       LATITUDE LONGITUDE      ELEVATION
c       DISTANCE(km)  DISTANCE(DEG) EVSTAZ STEVAZ
c       PULSE DESCRIPTION
c       COMMENT
c       UNITS - PRESSURE
c       TP TSV TSH  (Version .02 else LINE13)
c       LINE14
c       LINE15
c       LINE16
c       JSRC
c-----
c       For each trace there is the following
c-----

        subroutine rdhd96(LUN,nerr)
c-----
c       subroutine arguments
c
c       read a file96 header from the current read position
c
c       LUN I*4 logical unit for read
c       nerr    I*4  0 successful
c               -1 EOF (end of file)
c               -2 ERR (read error)
c               -3 not FILE01 FILE03 FILE16 or FILE21 format
c               -4 not OBSERVED or SYNTHETIC
c               -5 not correct format
c               
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
c               13  - 1/sec (used for ratio)
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
c       evdep   R*4 event depth (down = positive)
c       
c       stname  C*8 station name
c               NOTE: AH(6), SEED(5), CSS(6), SAC(8)
c       
c       stlat   R*4 station latitude
c       stlon   R*4 station longitude
c       stelev  R*4 station elevation (for Green's functions
c                   this is receiver depth (down = positive)
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
c
c       jsrc    I*4 Array of indices to indicate if traces are 
c                   present
c               iftype =  1  jsrc(1) indicates if trace is 
c                   present
c               iftype =  3  jsrc(i),i=1,5 indicates if trace 
c                   is present
c                   1 = Z, 2=N, 3=E, 4=R, 5=T
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
     2      tp, tsv, tsh,
     3      sa, sc, sf, sl, sn, sr
        real*4 esec, distkm, distdg, evstaz, stevaz
        real*4 evlat, evlon, evdep 
        real*4 stlat, stlon, stelev
        real*4 tp, tsv, tsh
        real*4 sa, sc, sf, sl, sn, sr

        common/chdr96/ stname, cfilt, cpulse, ccomnt
        character*8 stname*8, cfilt*80, cpulse*80, ccomnt*80    

        INTEGER VERSION


c-----
c       variables used internally in the routine
c-----
        character str*80
        parameter (NUNIT=13)
        character*16 cunit(NUNIT)
        data cunit/ 'COUNTS          ',
     1  'CM              ', 'CM/SEC          ', 'CM/SEC/SEC      ',
     2  'M               ', 'M/SEC           ', 'M/SEC/SEC       ',
     3  'MICRON          ', 'MICRON/SEC      ', 'MICRON/SEC/SEC  ',
     4  'Pa              ', 'MPa             ', '1/SEC           '/

c-----
c       get file type, later on validate version number
c`      LINE 1: FILE TYPE
c-----
        read(LUN,'(a)',err=9800,end=9900)str
        if(str(1:6).eq.'FILE01' .or. str(1:6).eq.'file01')then
            iftype = 1
        else if(str(1:6).eq.'FILE03' .or. str(1:6).eq.'file03')then
            iftype = 3
        else if(str(1:6).eq.'FILE16' .or. str(1:6).eq.'file16')then
            iftype = 16
        else if(str(1:6).eq.'FILE21' .or. str(1:6).eq.'file21')then
            iftype = 21
        else
            nerr = -3
            return
        endif
C-----
C       TEST OLD FILE??.01 VERSION  WHICH DID NOT HAVE TP TSV and TSH
C-----
        IF(STR(8:9).EQ.'01')THEN
            VERSION = 1
        ELSE IF(STR(8:9).EQ.'02')THEN
            VERSION = 2
        ELSE
            VERSION = 0
        ENDIF
c-----
c       get type of trace
c       LINE 2: TRACE TYPE
c-----
        read(LUN,'(a)',err=9800,end=9900)str
        if(str(1:3).eq.'OBS' .or. str(1:3).eq.'obs')then
            iobsyn = 1
        else if(str(1:3).eq.'SYN' .or. str(1:3).eq.'syn')then
            iobsyn = 2
        else
            nerr = -4
            return
        endif
c-----
c       get time series or spectra
c       LINE 3: DOMAIN
c-----
        read(LUN,'(a)',err=9800,end=9900)str
        if(str(1:3).eq.'TIM' .or. str(1:3).eq.'tim')then
            itmfrq = 1
        else if(str(1:3).eq.'FRE' .or. str(1:3).eq.'fre')then
            itmfrq = 2
        else
            nerr = -5
            return
        endif
c-----
c       get units
c       LINE 4: UNITS
c           for input parse from longest to shortest
c-----
        read(LUN,'(a)',err=9800,end=9900)str
        if(str(1:3).eq.'COU' .or. str(1:3).eq.'cou')then
            iunit = 1
        else if(str(1:10).eq.'CM/SEC/SEC' 
     1      .or. str(1:10).eq.'cm/sec/sec')then
            iunit = 4
        else if(str(1:6).eq.'CM/SEC' 
     1      .or. str(1:6).eq.'cm/sec')then
            iunit = 3
        else if(str(1:2).eq.'CM' 
     1      .or. str(1:2).eq.'cm')then
            iunit = 2
        else if(str(1:9).eq.'M/SEC/SEC' 
     1      .or. str(1:9).eq.'m/sec/sec')then
            iunit = 7
        else if(str(1:5).eq.'M/SEC' 
     1      .or. str(1:5).eq.'m/sec')then
            iunit = 6
        else if(str(1:14).eq.'MICRON/SEC/SEC' 
     1      .or. str(1:14).eq.'micron/sec/sec')then
            iunit = 10
        else if(str(1:10).eq.'MICRON/SEC' 
     1      .or. str(1:10).eq.'micron/sec')then
            iunit = 9
        else if(str(1:6).eq.'MICRON' 
     1      .or. str(1:6).eq.'micron')then
            iunit = 8
        else if(str(1:1).eq.'M' .and.str(1:2).ne.'MI' 
     1      .or. str(1:1).eq.'m'.and. str(1:2).ne.'mi')then
            iunit = 5
        else if(str(1:5).eq.'1/SEC'.or. str(1:5).eq.'1/sec')then
            iunit = 13
        else
            nerr = -5
            return
        endif
c-----
c       get filtering comment
c       LINE 5: FILTERING HISTORY
c-----
        read(LUN,'(a)',err=9800,end=9900)cfilt
c-----
c       get event information   
c       LINE 6: EVENT INFORMATION
c-----
        read(LUN,*,err=9800,end=9900)
     1      keyear, kemon, keday, kehour, kemin, esec,
     2      evlat, evlon, evdep
c-----
c       get station information 
c       LINE 7: STATION NAME
c-----
        read(LUN,'(a)',end=9800,err=9800)stname
c-----
c       LINE 8: STATION COORDINATES
c-----
        read(LUN,*,end=9800,err=9800)
     1      stlat, stlon, stelev
c-----
c       get distance, azimuth if known
c       LINE 9: DISTANCE
c-----
        read(LUN,*,err=9800,end=9900)
     1      distkm, distdg, evstaz, stevaz
c-----
c       get comment on source pulse
c       LINE 10: SOURCE PULSE COMMENT
c-----
        read(LUN,'(a)',err=9800,end=9900)
     1      cpulse
c-----
c       get comment line
c       LINE 11: COMMENT - may be earth model name
c-----
        read(LUN,'(a)',err=9800,end=9900)
     1      ccomnt
c-----
c       get pressure units
c       LINE 12: pressure units
c-----
        read(LUN,'(a)',err=9800,end=9900)str
        if(str(1:2).eq.'Pa' 
     1      .or. str(1:2).eq.'pa')then
            junit = 11
        else if(str(1:3).eq.'MPa' 
     1      .or. str(1:3).eq.'mpa')then
            junit = 12
        endif
c-----
c       get P, SV, and SH first arrival times
c       LINE 13
c-----
        IF(VERSION.GT.1)THEN
            READ(lun,*,ERR=9800,END=9900)TP, TSV, TSH
        ELSE
            READ(lun,'(A)',ERR=9800,END=9900)STR
            TP = -12345.0
            TSV = -12345.0
            TSH = -12345.0
        ENDIF
        if(version.lt.2)then
            read(LUN,'(a)',err=9800,end=9900)str
            sa = 0.0
            sc = 0.0
            sf = 0.0
            sl = 0.0
            sn = 0.0
            sr = 0.0
        else if(version.eq.2)then
            read(LUN,*)sa,sc,sf,sl,sn,sr
        endif
c-----
c       skip over the expansion fields
c       LINE 15-16:
c-----
        read(LUN,'(a)',err=9800,end=9900)str
        read(LUN,'(a)',err=9800,end=9900)str
c-----
c       get jsrc
c       LINE 17:
c-----
        if(iftype.le.16)then
            read(LUN,*,err=9800,end=9900)(jsrc(i),i=1,16)
            jsrc(17) = 0
            jsrc(18) = 0
            jsrc(19) = 0
            jsrc(20) = 0
            jsrc(21) = 0
        else
            read(LUN,*,err=9800,end=9900)(jsrc(i),i=1,21)
        endif
c-----
c       OK HEADER IS COMPLETELY READ
c-----
        nerr = 0
        return
 9800   continue
        nerr = -2
        return
 9900   nerr = -1
        return
        end

        subroutine wrhd96(LUN,nerr)
c-----
c       subroutine arguments
c
c       write a file96 header at the current write position
c
c       LUN I*4 logical unit for read
c       nerr    I*4 -1 error in permitted variable value
c
c       for description of variables see rdhd96()
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
     1      tp, tsv, tsh,
     3      sa, sc, sf, sl, sn, sr
        real*4 esec, distkm, distdg, evstaz, stevaz
        real*4 evlat, evlon, evdep 
        real*4 stlat, stlon, stelev
        real*4 tp, tsv, tsh
        real*4 sa, sc, sf, sl, sn, sr

        common/chdr96/ stname, cfilt, cpulse, ccomnt
        character*8 stname*8, cfilt*80, cpulse*80, ccomnt*80    
        

c-----
c       variables used internally in the routine
c-----
        character*80 fmt1, fmt2, fmt3, fmt4, fmt5, fmt6, fmt7
        parameter (NUNIT=13)
        character*16 cunit(NUNIT)
        data cunit/ 'COUNTS          ',
     1  'CM              ', 'CM/SEC          ', 'CM/SEC/SEC      ',
     2  'M               ', 'M/SEC           ', 'M/SEC/SEC       ',
     3  'MICRON          ', 'MICRON/SEC      ', 'MICRON/SEC/SEC  ',
     4  'Pa              ', 'MPa             ', '1/SEC           '/

c-----
c       OUTPUT FORMATS
c       ON UNIX, there is no CARRIAGE CONTROL
c       ON MSDOS, there is only carriage control to screen, 
c           e.g., unit 6
c       NOTE on UNIX need '' -> '
c       MS FORTRAN 4.0 need '' -> '
c-----
C       if(LUN.eq.6)then
C           fmt1 = '(1x,a)'
C           fmt2 = '(1x,21i3)'
C           fmt3 ='(1x,i5,1x,i2,1x,i2,1x,i2'//
C     1         ',1x,i2,1x,g12.5,1x,f10.6'//
C     2         ',1x,f11.6,1x,f10.6)'
C           fmt4='(1x,a8)'
C           fmt5='(1x,f10.6,1x,f11.6,1x,f10.3)'
C           fmt6='(1x,f15.4,1x,f15.4,1x,'//
C     1         'f15.4,1x,f15.4)'
C           fmt7='(1x,f10.3, 1x, f10.3, 1x, f10.3,3i5)'
C       else    
            fmt1 = '(     a)'
            fmt2 = '(     21i3)'
            fmt3 ='(i5,1x,i2,1x,i2,1x,i2'//
     1          ',1x,i2,1x,g12.5,1x,f13.6'//
     2          ',1x,f13.6,1x,f13.6)'
            fmt4='(a8)'
            fmt5='(f13.6,1x,f13.6,1x,f10.3)'
            fmt6='(f15.4,1x,f15.4,1x,f15.4,1x,f15.4)'
            fmt7='(f10.3, 1x, f10.3, 1x, f10.3,3i5)'
C       endif
c-----
c       write file type, later on validate version number
c-----
        if(iftype.eq.1)then
            write(LUN, fmt1)'FILE01.02'
        else if(iftype.eq.3)then
            write(LUN, fmt1)'FILE03.02'
        else if(iftype.eq.16)then
            write(LUN, fmt1)'FILE16.02'
        else if(iftype.eq.21)then
            write(LUN, fmt1)'FILE21.02'
        else
            nerr = -1
            return
        endif
c-----
c       write type of trace
c-----
        if(iobsyn.eq.1)then
            write(LUN, fmt1)'OBSERVED'
        else if(iobsyn.eq.2)then
            write(LUN, fmt1)'SYNTHETIC'
        else
            nerr = -1
            return
        endif
c-----
c       write time series for spectra
c-----
        if(itmfrq.eq.1)then
            write(LUN, fmt1)'TIME_DOMAIN'
        else if(itmfrq.eq.2)then
            write(LUN, fmt1)'FREQUENCY_DOMAIN'
        else
            nerr = -1
            return
        endif
c-----
c       write units
c-----
        if(iunit.ge.1 .and. iunit.le.10)then
            write(LUN,fmt1)cunit(iunit)
        else if(iunit.eq.13)then
            write(LUN,fmt1)cunit(iunit)
        else
            nerr = -1
            return
        endif
c-----
c       write filtering comment
c-----
        write(LUN,fmt1)cfilt
c-----
c       write event information 
c-----
        write(LUN, fmt3)
     1      keyear, kemon, keday, kehour, kemin, 
     2      esec, evlat, evlon, evdep
c-----
c       write station information   
c-----
        write(LUN,fmt4)
     2      stname
        write(LUN,fmt5)
     2      stlat, stlon, stelev
c-----
c       write distance, azimuth if known
c-----
        write(LUN,fmt6)
     1      distkm, distdg, evstaz, stevaz
c-----
c       write comment on source pulse
c       LINE 10
c-----
        write(LUN,fmt1)cpulse
c-----
c       write comment line
c       LINE 11
c-----
        write(LUN,fmt1)ccomnt
c-----
c       write units for Pressure
c       LINE 12
c-----
        if(junit .lt.11 .or. junit.gt.12)then
            write(LUN,fmt1)'NONE            '
        else
            write(LUN,fmt1)cunit(junit)
        endif
c-----
c       write P, SV, and SH first arrival times
c       LINE 13
c-----
        write(LUN,'(3g20.10)')tp, tsv, tsh
c-----
c       write A, C, F, L, N and rho at source
c       LINE 14
c-----
        write(LUN,'(6e12.4)')sa,sc,sf,sl,sn,sr
c-----
c       output the expansion fields
c       LINE 15-16:
c-----
        write(LUN,fmt1)'LINE15'
        write(LUN,fmt1)'LINE16'
c-----
c       write jsrc
c-----
        write(LUN,fmt2)jsrc
c-----
c       OK HEADER IS COMPLETELY WRITTEN
c-----
        nerr = 0
        return
        end

        subroutine rdtr96(LUN,stcomp,cmpinc,cmpaz,cmpdt, npts,
     1      ksyear, ksmon, ksday, kshour, ksmin, ssec, 
     2      x, nerr, maxarr)
c-----
c       read next trace from the file
c
c       LUN I*4 logical unit for read
c       stcomp  C*8 station component, e.g., BHZ
c               NOTE: AH(6), SEED(3), CSS(8), SAC(8)
c       cmpaz   R*4 Positive motion in this azimuthal direction
c                   0 = north, 90 = east (SEED)
c       cmpinc  R*4 Positive motion in this direction with vertical
c                   -90 = up, 0 = horiz, 90 = down (SEED)
c       cmpdt   R*4 Component sample interval
c       npts    I*4 Number of samples
c       ksyear  I*4 station year - time of first sample
c       ksmon   I*4 station mon
c       ksday   I*4 station day
c       kshour  I*4 station hour
c       ksmin   I*4 station minute
c       ssec    R*4 station second
c       x   R*4 time series array
c       nerr    I*4  0 successful
c               -1 EOF (end of file)
c               -2 ERR (read error)
c       maxarr  I*4 maximum array dimension for x()
c-----
        integer*4 LUN
        character stcomp*8
        real*4 cmpaz, cmpinc, cmpdt
        integer*4 npts
        integer*4 ksyear, ksmon, ksday, kshour, ksmin
        real*4 ssec
        real*4 x(npts)
        integer*4 nerr, maxarr
c-----
c       FORMAT FOR UNIX INPUT
c-----
    1   format(a8)
C    2   format(g10.3, 1x,g10.3, 1x, g10.3,1x,i5)
C    3   format(i5,1x,i2,1x,i2,1x,i2,1x,i2,1x,
C     1      g10.4)
C    4   format(e16.8,1x,e16.8,1x,e16.8,1x,e16.8)

c-----
c       get component name, orientation incidence, orientation azimuth
c-----
        read(LUN,1,err=9800,end=9900)stcomp
        read(LUN,*,err=9800,end=9900)
     1      cmpinc, cmpaz, cmpdt, npts
        if(npts.gt.maxarr)go to 9900
c-----
c       get station start time and sample interval  
c-----
        read(LUN,*,err=9800,end=9900)
     1      ksyear, ksmon, ksday, kshour, ksmin, ssec
c-----
c       get the time history
c-----
        read(LUN,*,err=9800,end=9900)(x(i),i=1,npts)
        nerr = 0
        return
 9800   continue
            nerr = -2
            return
 9900   continue
            nerr = -1
        return
        end


        subroutine wrtr96(LUN,stcomp,cmpinc,cmpaz,cmpdt, npts,
     1      ksyear, ksmon, ksday, kshour, ksmin, ssec, 
     2      x, nerr, maxarr)
c-----
c       write next trace to the file
c
c       LUN I*4 logical unit for read
c       stcomp  C*8 station component, e.g., BHZ
c               NOTE: AH(6), SEED(3), CSS(8), SAC(8)
c       cmpaz   R*4 Positive motion in this azimuthal direction
c                   0 = north, 90 = east (SEED)
c       cmpinc  R*4 Positive motion in this direction with vertical
c                   -90 = up, 0 = horiz, 90 = down (SEED)
c       cmpdt   R*4 Component sample interval
c       npts    I*4 Number of samples
c       ksyear  I*4 station year - time of first sample
c       ksmon   I*4 station mon
c       ksday   I*4 station day
c       kshour  I*4 station hour
c       ksmin   I*4 station minute
c       ssec    R*4 station second
c       x   R*4 time series array
c       nerr    I*4 not used just here for symmetry
c       maxarr  I*4 maximum array dimension for x()
c-----
c       note the output consists of 4 entries per line. 
c           If the number of
c       points is not a multiple of 4, then zero fill. Also
c       assume nothing about the array dimensions
c-----
        integer*4 LUN
        character stcomp*8
        real*4 cmpaz, cmpinc, cmpdt
        integer*4 npts
        integer*4 ksyear, ksmon, ksday, kshour, ksmin
        real*4 ssec
        real*4 x(npts)
        integer*4 nerr, maxarr
        character*80 fmt1, fmt2, fmt3, fmt4
c-----
c       OUTPUT FORMATS
c       ON UNIX, there is no CARRIAGE CONTROL
c       ON MSDOS, there is only carriage control to screen, e.g., 
c           unit 6
c       NOTE on UNIX need '' -> '
c       MS FORTRAN 4.0 need '' -> '
c-----
C       if(LUN.eq.6)then
C           fmt1='(1x,a8)'
C           fmt2='(1x,f10.3, 1x,f10.3,1x, g12.5,1x,i5)'
C           fmt3='(1x,i5,1x,i2,1x,i2,1x,i2 ,1x,i2,1x,g12.5)'
C           fmt4='(1x,e16.8,1x,e16.8 ,1x,e16.8,1x,e16.8)'
C       else
            fmt1='(a8)'
            fmt2='(f10.3, 1x,f10.3, 1x, g12.5,1x,i5)'
            fmt3='(i5,1x,i2,1x,i2,1x,i2 ,1x,i2,1x,g12.5)'
            fmt4='(e16.8,1x,e16.8,1x,e16.8,1x,e16.8)'
C       endif

c-----
c       write component name, orientation incidence, 
c           orientation azimuth
c-----
        write(LUN,fmt1)
     1      stcomp
        write(LUN,fmt2)
     1      cmpinc, cmpaz, cmpdt, npts
c-----
c       get station start time and sample interval  
c-----
        write(LUN,fmt3)
     1      ksyear, ksmon, ksday, kshour, ksmin, ssec
c-----
c       get the time history
c-----
        num4 = npts/4
        n4 = num4 * 4
        if(n4.gt.0)then
            write(LUN,fmt4)(x(i),i=1,n4)
        endif
        ndo = npts - n4
        zero = 0.0
        if(ndo.eq.1)then
            write(LUN,fmt4)(x(i),i=n4+1,npts),zero,zero,zero
        else if(ndo.eq.2)then
            write(LUN,fmt4)(x(i),i=n4+1,npts),zero,zero
        else if(ndo.eq.3)then
            write(LUN,fmt4)(x(i),i=n4+1,npts),zero
        endif
        return
        end
