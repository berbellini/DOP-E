
        program sactof96
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME V                                                       c
c                                                                     c
c      PROGRAM: SACTOF96                                              c
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
c       Convert SAC to file96(V) format 
c
c       sactof96 
c
c       This reads a command language from the standard input
c       The keywords are
c
c       FILE03  Designates type of file96 file
c       FILE21
c       FILE01
c
c       SYNTHETIC
c       OBSERVED
c
c       TIME_DOMAIN
c       FREQUENCY_DOMAIN    (nothing uses or checks this yet)
c
c       UNITS
c           CM      M   MICRON
c           CM/SEC      M/S MICRON/S
c           CM/SEC/SEC  M/S/S   MICRON/S/S
c           COUNTS
c
c       1 ... 21 for position or Z N E R T O
c       HELP

c-----
c-----
c       variables in common blocks
c       command line arguments
c-----
        logical verby
c-----
c       internal program variables
c-----
        parameter (LIN=5, LOT=6, LER=6)

c-----
c       variables for command language
c-----
        logical isbin
        integer*4 ksrc(21)
        character*80 fname(21)
        character*80 in, inn

        parameter (NUNIT=10)
        character*16 cunit(NUNIT)

        data cunit/ 'COUNTS          ',
     1  'CM              ', 'CM/SEC          ', 'CM/SEC/SEC      ',
     2  'M               ', 'M/SEC           ', 'M/SEC/SEC       ',
     3  'MICRON          ', 'MICRON/SEC      ', 'MICRON/SEC/SEC  '/

        data ksrc/21 * 0/
c-----
c       machine dependent initialization
c-----
        call mchdep()
c-----
c       parse command line arguments
c-----
        call gcmdln(verby)
c-----
c       initialize
c-----
        iftype = -1
        iobsyn = -1
        itmfrq = -1
        iunit = -1
        junit = 11
        isbin = .true.

c-----
c       present limitations on use of file96
c       TIME_DOMAIN 
        itmfrq = 1
c-----
c       parse the parameter stream
c-----
 1000   continue
            if(verby)write(LOT,*)'ENTER COMMAND PARAMETER'
            read(LIN,'(a)')in
            call touppr(in)
            if(in .eq. 'HELP')then
                call help()
            else if(in .eq. 'REVIEW')then
                    call review(
     1              isbin, ksrc, fname,
     2              iftype, iobsyn, itmfrq,iunit)
            else if(in .eq. 'BINARY')then
                isbin = .true.
            else if(in .eq. 'ASCII')then
                isbin = .false.
            else if(in .eq. 'OBSERVED')then
                iobsyn = 1
            else if(in .eq. 'SYNTHETIC')then
                iobsyn = 2
            else if(in .eq. 'TIME_DOMAIN')then
                itmfrq = 1
            else if(in .eq. 'FREQUENCY_DOMAIN')then
                itmfrq = 2
            else if(in .eq. 'PROCESS')then
                call process(
     1              isbin, ksrc, fname,
     2              iftype, iobsyn, itmfrq,iunit)
                go to 1100
            else if(in .eq. 'QUIT')then
                go to 1100
            else if (in .eq. 'UNIT')then
                read(LIN,'(a)')inn
                call touppr(inn)
c-----
c       blank fill to make a character*16
c       then search for pattern
c-----
                linn = lgstr(inn)
                if(linn.lt.16)then
                    do 1200 i=linn+1,16
                        inn(i:i) = ' '
 1200               continue
                endif
                iunit = -1
                do 1300 i=1,NUNIT
                    if(cunit(i)(1:16).eq.inn(1:16))then
                        iunit = i
                    endif
 1300           continue
            else if (in .eq. 'FILE01')then
                iftype = 1
            else if (in .eq. 'FILE03')then
                iftype = 3
            else if (in .eq. 'FILE16')then
                iftype = 16
            else if (in .eq. 'FILE21')then
                iftype = 21
            else if(in.eq.'01')then
                ksrc(1) = 1
                read(LIN,'(a)')fname(1)
            else if(in.eq.'02')then
                ksrc(2) = 4
                read(LIN,'(a)')fname(2)
            else if(in.eq.'03')then
                ksrc(3) = 1
                read(LIN,'(a)')fname(3)
            else if(in.eq.'04')then
                ksrc(4) = 4
                read(LIN,'(a)')fname(4)
            else if(in.eq.'05')then
                ksrc(5) = 5
                read(LIN,'(a)')fname(5)
            else if(in.eq.'06')then
                ksrc(6) = 1
                read(LIN,'(a)')fname(6)
            else if(in.eq.'07')then
                ksrc(7) = 4
                read(LIN,'(a)')fname(7)
            else if(in.eq.'08')then
                ksrc(8) = 5
                read(LIN,'(a)')fname(8)
            else if(in.eq.'09')then
                ksrc(9) = 1
                read(LIN,'(a)')fname(9)
            else if(in.eq.'10')then
                ksrc(10) = 4
                read(LIN,'(a)')fname(10)
            else if(in.eq.'11')then
                ksrc(11) = 1
                read(LIN,'(a)')fname(11)
            else if(in.eq.'12')then
                ksrc(12)= 4
                read(LIN,'(a)')fname(12)
            else if(in.eq.'13')then
                ksrc(13) = 1
                read(LIN,'(a)')fname(13)
            else if(in.eq.'14')then
                ksrc(4) = 4
                read(LIN,'(a)')fname(14)
            else if(in.eq.'15')then
                ksrc(15) = 5
                read(LIN,'(a)')fname(15)
            else if(in.eq.'Z')then
                ksrc(1) = 1
                read(LIN,'(a)')fname(1)
            else if(in.eq.'N')then
                ksrc(2) = 2
                read(LIN,'(a)')fname(2)
            else if(in.eq.'E')then
                ksrc(3) = 3
                read(LIN,'(a)')fname(3)
            else if(in.eq.'R')then
                ksrc(4) = 4
                read(LIN,'(a)')fname(4)
            else if(in.eq.'T')then
                ksrc(5) = 5
                read(LIN,'(a)')fname(5)
            else if(in.eq.'O')then
                ksrc(6) = 6
                read(LIN,'(a)')fname(6)
            endif
        go to 1000
 1100   continue
        end

        subroutine touppr(str)
c-----
c       CONVERT THE ASCII STRING str to UPPER CASE
c
c       NOTE THIS IS NON PORTABLE FORTRAN SINCE IT ASSUMES THE
C       ASCII CHARACTER SET
c-----
c       str C*(*)   - ASCII string
c-----
        character str*(*)
        character c*1
        lstr = lgstr(str)
        do 1000 i=1,lstr
            if(str(i:i) .ge. 'a' .and. str(i:i).le.'z')then
                c = str(i:i)
                str(i:i)=char(ichar('A')-ichar('a')+ichar(c))
            endif
 1000   continue
        return
        end
                
        subroutine help()
        parameter (LIN=5, LOT=6, LER=6)
        character*80 str(42)
        str( 1) =' Control the creation of FILE96(V) files by these'
        str( 2) =' commands:'
        str( 3) =' KEYWORD           MEANING'
        str( 4) =' ______________________________________________'
        str( 5) =' HELP              Print this help screen'
        str( 6) =' REVIEW            Print current input list'
        str( 7) =' BINARY            SAC file is in binary form'
        str( 8) =' ASCII             SAC file is in ASCII form'
        str( 9) =' OBSERVED          Series is observed'
        str(10) =' SYNTHETIC         Series is synthetic'
        str(11) =' TIME_DOMAIN       This is a time series'
        str(12) =' FREQUENCY_DOMAIN  This is a time series'
        str(13) =' PROCESS           Begin conversion to FILE96'
        str(14) =' QUIT              Exit this program'
        str(15) =' FILE01            Create a  1 component FILE96'
        str(16) =' FILE03            Create a  3 component FILE96'
        str(17) =' FILE16            Create a 16 component FILE96'
        str(18) =' UNIT              The next line is one of the'
        str(19) ='                   following left justified keywords:'
        str(20) ='        COUNTS CM CM/SEC CM/SEC/SEC M M/SEC '
        str(21) ='        M/SEC/SEC MICRON MICRON/SEC MICRON/SEC/SEC'
        str(22) =' 01           The following FILE is ZDD'
        str(23) =' 02           The following FILE is RDD'
        str(24) =' 03           The following FILE is ZDS'
        str(25) =' 04           The following FILE is RDS'
        str(26) =' 05           The following FILE is TDS'
        str(27) =' 06           The following FILE is ZSS'
        str(28) =' 07           The following FILE is RSS'
        str(29) =' 08           The following FILE is TSS'
        str(30) =' 09           The following FILE is ZEX'
        str(31) =' 10           The following FILE is REX'
        str(32) =' 11           The following FILE is ZVF'
        str(33) =' 12           The following FILE is RVF'
        str(34) =' 13           The following FILE is ZHF'
        str(35) =' 14           The following FILE is RHF'
        str(36) =' 15           The following FILE is THF'
        str(37) =' 16           The following FILE is PEX'
        str(37) =' Z            The following FILE is Vertical'
        str(38) =' N            The following FILE is North'
        str(39) =' E            The following FILE is East'
        str(40) =' R            The following FILE is Radial'
        str(40) =' T            The following FILE is '
     1          //'Transverse'
        str(41) =' O            The following FILE is Other'
        str(42) =' ______________________________________________'
        do 1000 i=1,21
            lstr = lgstr(str(i))
            write(LOT,*)str(i)(1:lstr)
 1000   continue
        write(LOT,*)'_____________________'
        write(LOT,*)'ENTER for next screen'
        read(LIN,'(1x)')
        do 1100 i=22,42
            lstr = lgstr(str(i))
            write(LOT,*)str(i)(1:lstr)
 1100   continue
        write(LOT,*)'_____________________'
        write(LOT,*)'ENTER to continue'
        read(LIN,'(1x)')
        return
        end

        subroutine review(
     1              isbin, ksrc, fname,
     2              iftype, iobsyn, itmfrq,iunit)
c-----
c       show present state of input
c-----  
        logical isbin
        integer*4 ksrc(21)
        character*80 fname(21)
        integer*4 itype, iobsyn, itmfrq, iunit

        parameter (LIN=5, LOT=6, LER=6)

        character ostr*80
        logical ext

        character*4 cgrn(21), cobs(6), ngrn(21)
        parameter (NUNIT=12)
        character*21 cunit(NUNIT)

        data cunit/ 'COUNTS          ',
     1  'CM              ', 'CM/SEC          ', 'CM/SEC/SEC      ',
     2  'M               ', 'M/SEC           ', 'M/SEC/SEC       ',
     3  'MICRON          ', 'MICRON/SEC      ', 'MICRON/SEC/SEC  ',
     4  'Pa              ', 'MPa             '                    /

        data cgrn/'ZDD ', 'RDD ', 'ZDS ', 'RDS ',
     1        'TDS ', 'ZSS ', 'RSS ', 'TSS ',
     2        'ZEX ', 'REX ', 'ZVF ', 'RVF ',
     3        'ZHF ', 'RHF ', 'THF ', 'PEX ',
     4        'PDD ', 'PDS ', 'PSS ', 'PVF ', 'PHF '/

        data ngrn/'01  ', '02  ', '03  ', '04  ',
     1        '05  ', '06  ', '07  ', '08 ',
     2        '09  ', '10  ', '11  ', '12  ',
     3        '13  ', '14  ', '15  ', '16  ',
     3        '17  ', '18  ', '19  ', '20  ', '21  '/

        data cobs/'Z   ', 'N   ', 'E   ', 'R   ',
     1        'T   ', 'O   '/
c-----  
c       file type
c-----
        if(iftype.eq.1)then
            ostr = 'FILE01'
        else if(iftype .eq. 3)then
            ostr = 'FILE03'
        else if(iftype .eq. 16)then
            ostr = 'FILE16'
        else if(iftype .eq. 21)then
            ostr = 'FILE21'
        else
            ostr = 'Not Defined'
        endif
        lostr = lgstr(ostr)
        write(LOT,*)    '     File Type:      ',ostr(1:lostr)
c-----
c       Observed or synthetic
c-----
        if(iobsyn.eq.1)then
            ostr ='OBSERVED'
        else if(iobsyn.eq.2)then
            ostr ='SYNTHETIC'
        else
            ostr = 'Not Defined'
        endif
        lostr = lgstr(ostr)
        write(LOT,*)    '     Data Source:    ',ostr(1:lostr)
c-----
c       Time or frequency domain
c-----
        if(itmfrq.eq.1)then
            ostr ='TIME_DOMAIN'
        else if(itmfrq.eq.2)then
            ostr ='FREQUENCY_DOMAIN'
        else
            ostr = 'Not Defined'
        endif
        lostr = lgstr(ostr)
        write(LOT,*)    '     Domain:         ',ostr(1:lostr)
c-----
c       units
c-----
        if(iunit.ge.1 .and. iunit.le.10)then
            ostr = cunit(iunit)
        else
            ostr = 'Not Defined'
        endif
        lostr = lgstr(ostr)
        write(LOT,*)    '     Units:          ',ostr(1:lostr)
c-----
c       type of SAC file
c-----
        if(isbin)then
            write(LOT,*)    '     SAC File Type:  ','BINARY'
        else
            write(LOT,*)    '     SAC File Type:  ','ASCII'
        endif
c-----
c       cycle through the files
c-----
        do 1000 i=1,21
            if(ksrc(i).gt.0)then
                inquire(file=fname(i),exist=ext)
                lnm = lgstr(fname(i))
                if(ext)then
                    ostr = fname(i)
                    lostr = lgstr(ostr)
                    ostr = ostr(1:lostr)//' '//
     1                  '-- Exists'
                else
                    ostr = fname(i)
                    lostr = lgstr(ostr)
                    ostr = ostr(1:lostr)//' '//
     1                  '-- Does Not Exist'
                endif
                lostr = lgstr(ostr)
                if(iftype.eq.21)then
        write(LOT,*)ngrn(i)//' SAC Data File:  ',ostr(1:lostr)
                else
        write(LOT,*)cobs(i)//' SAC Data File:  ',ostr(1:lostr)
                endif
            
            endif
 1000   continue
        return
        end

        subroutine process(
     1              isbin, ksrc, fname,
     2              kftype, kobsyn, ktmfrq,kkunit)
c-----
c       isbin   L   - .true. binary SAC file
c                 .false. ASCII SAC file
c       ksrc(21)I*4 - array of indices
c       fname(21) C*(*) - file names
c       jftype  I*4 -  1 = FILE01
c                  3 = FILE03
c                 16 = FILE16
c                 21 = FILE21
c       jobsyn  I*4 -  1 = observed time history
c                  2 = synthetic time history
c       jtmfrq  I*4 -  1 = TIME_DOMAIN
c                  2 = FREQUENCY_DOMAIN
c       junit   I*4 - 1 integer from 1 to 21 indicating units
c-----
        logical isbin
        integer*4 ksrc(21)
        character*80 fname(21)
        integer*4 kftype, kobsyn, ktmfrq, kunit
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
        integer*4 npts
        real*4 ssec
c-----
c       internal program variables
c-----
        parameter (LIN=5, LOT=6, LER=6)
        parameter (NSAMP=8192)
        real*4 x(NSAMP)
        logical ext
        real*8 eepoch, sepoch, repoch
        character*50 str
        integer*4 doy, date
        character ofile*80

c-----
c       open the output file
c-----
        open(4,file='file96.fil',form='formatted',
     1      access='sequential',status='unknown')
        rewind 4
c-----
c       set up FILE96 header
c-----
        iftype = kftype
        iunit  = kkunit
        iobsyn = kobsyn
        itmfrq = ktmfrq
        do 1000 i=1,21
            jsrc(i) = ksrc(i)
 1000   continue
c-----
c       scan through the ksrc list to get the first file.
c       This file will set the header
c-----
        ndone = 0
        do 1100 i=1,21
            if(ksrc(i).gt.0 .and. ndone.le.iftype)then
                inquire(file=fname(i),exist=ext)
                lnm = lgstr(fname(i))
                if(ext)then
                    ofile=fname(i)
                    if(isbin) then
                    call brsac (1,NSAMP,ofile,x,nerr)
                    else
                    call arsac (1,NSAMP,ofile,x,nerr)
                    endif
                    if(nerr.lt.0)go to 1150
                    go to 1200
                endif
                ndone = ndone + 1
            endif
 1100   continue
        write(LER,*)'None of the SAC files exist'
        return
 1150   continue
        ll = lgstr(ofile)
        write(LER,*)'Read error on SAC file:',ofile(1:ll)
        return
 1200   continue
c-----
c       we have a valid SAC header. Get the parameters required
c       for the FILE96 header
c-----
        call getfhv('EVDP    ',evdep ,nerr)
c-----
c       SAC depth is METERS, ours is KM
c-----
        evdep = evdep / 1000.0
        call getfhv('DIST    ',distkm,nerr)
        call getfhv('AZ      ',evstaz,nerr)
        call getfhv('BAZ     ',stevaz,nerr)
        call getfhv('GCARC   ',distdg,nerr)
        call getkhv('KSTNM   ',stname ,nerr)
        cpulse = ' '
        ccomnt = ' '
        call getfhv('STLA    ',stlat ,nerr)
        call getfhv('STLO    ',stlon ,nerr)
        call getfhv('STEL    ',stelev,nerr)
c-----
c       SAC elevation is  METERS, ours is KM
c-----
        stelev = stelev / 1000.0
        call getfhv('EVLA    ',evlat ,nerr)
        call getfhv('EVLO    ',evlon ,nerr)
c-----
c       This is the SAC reference time
c-----
        call getnhv('NZYEAR  ',ksyear,nerr)
        call getnhv('NZJDAY  ',doy   ,nerr)
        call getnhv('NZHOUR  ',kshour,nerr)
        call getnhv('NZMIN   ',ksmin ,nerr)
        call getnhv('NZSEC   ',isec  ,nerr)
        call getnhv('NZMSEC  ',jsmsec,nerr)
        ssec = isec + jsmsec/1000.0
c-----
c       get P, SV and SH times
c-----
        call getfhv('A       ',tp    ,nerr)
        call getfhv('T0      ',tsv   ,nerr)
        call getfhv('T1      ',tsh   ,nerr)
c-----
c       set enumerated header value to define reference time
c-----
        call getihv('IZTYPE  ',str, herr)
        if(str .eq. 'IO      ')then
c-----
c           this is the origin time
c-----
        else if(str .eq. 'IB      ')then
c-----
c           this is the start time of the trace
c-----
        endif
c-----
c       get initial value of independent variable
c-----
        call getfhv('B       ',btime ,nerr)
c-----
c       get origin time with respect to reference
c-----
        call getfhv('O       ',otime ,nerr)
c-----
c       convert to epoch for origin time and start of trace time
c       after converting day of year to month and day, then
c       inverse to get proper year month ...
c-----
        call mnthdy(ksyear,doy,ksmon,ksday)
        call htoe(ksyear,ksmon,ksday,kshour,ksmin,ssec,repoch)
        eepoch = repoch + dble(otime)
        call etoh(eepoch,date,str,doy,keyear,
     1      kemon,keday,kehour,kemin,esec)
        if(iobsyn.eq.2)then
            keyear = keyear - 1970
            kemon = kemon -1
            keday = keday - 1
        endif

        lstr = lgstr(str)
        cfilt  = ' '
        cpulse = ' '
        ccomnt = ' '
c-----
c       we need to get the enumerated flag to indicate
c       what the first time point is and what time is
c-----

        call wrhd96(4,nerr)
        if(nerr .lt. 0)then
            write(LER,*)'Error writing FILE96 header'
            return
        endif
c-----
c       Now systematically read the SAC files and write the
c       trace files. error checking has already been done above
c-----
        do 3100 i=1,21
            if(ksrc(i).gt.0)then
                inquire(file=fname(i),exist=ext)
                lnm = lgstr(fname(i))
                if(ext)then
                    ofile=fname(i)
                    if(isbin) then
                    call brsac (1,NSAMP,ofile,x,nerr)
                    else
                    call arsac (1,NSAMP,ofile,x,nerr)
                    endif
                    if(nerr.lt.0)go to 3150
c-----
c       get the SAC header values from this trace and convert
c       for FILE96 trace header
c-----
                    call getkhv('KCMPNM  ',stcomp ,nerr)
c----
c       convert from SAC to SEED component incidence
c-----
                    call getfhv('CMPINC  ',cmpinc,nerr)
                    cmpinc = cmpinc - 90.0
                    call getfhv('CMPAZ   ',cmpaz ,nerr)
                    call getfhv('DELTA   ',cmpdt,nerr)
                    call getnhv('NPTS    ',npts, nerr)
c-----
c       get the start time of this trace
c-----
                    call getnhv('NZYEAR  ',ksyear,nerr)
                    call getnhv('NZJDAY  ',doy   ,nerr)
                    call getnhv('NZHOUR  ',kshour,nerr)
                    call getnhv('NZMIN   ',ksmin ,nerr)
                    call getnhv('NZSEC   ',isec  ,nerr)
                    call getnhv('NZMSEC  ',jsmsec,nerr)
                    ssec = isec + jsmsec/1000.0
c-----
c       set enumerated header value to define reference time
c-----
                    call getihv('IZTYPE  ',str, herr)
                    if(str .eq. 'IO      ')then
c-----
c           this is the origin time
c-----
                    else if(str .eq. 'IB      ')then
c-----
c           this is the start time of the trace
c-----
                    endif
c-----
c       get initial value of independent variable
c-----
                    call getfhv('B       ',btime ,nerr)
c-----
c       convert to epoch start of trace time
c       after converting day of year to month and day, then
c       inverse to get proper year month ...
c-----
                call mnthdy(ksyear,doy,ksmon,ksday)
                call htoe(ksyear,ksmon,ksday,kshour,
     1              ksmin,ssec,repoch)
                sepoch = repoch + btime
                call etoh(sepoch,date,str,doy,ksyear,
     1              ksmon,ksday,kshour,ksmin,ssec)
                if(iobsyn.eq.2)then
                    ksyear = ksyear - 1970
                    ksmon = ksmon -1
                    ksday = ksday - 1
                endif

c-----
c       now write the trace header and trace
c-----
                call wrtr96(4,stcomp,cmpinc,cmpaz,
     1              cmpdt, npts, ksyear, ksmon, 
     2              ksday, kshour, ksmin, ssec, 
     3              x,nerr,NSAMP)
c-----
                endif
            endif
 3100   continue
        return
 3150   continue
        ll = lgstr(ofile)
        write(LER,*)'Read error on SAC file:',ofile(1:ll)
        return
        close(4)
        return
        end

        subroutine gcmdln(verby)
c-----
c       parse command line arguments
c
c       requires subroutine mgtarg() and function mnmarg()
c
c-----
c       verby   L   - verbose output on standard error
c-----
        character*25 name
        logical verby

        verby = .false.
        nmarg=mnmarg()
        i = 0
   11   i = i + 1
        if(i .gt. nmarg)go to 13
            call mgtarg(i,name)
            if(name(1:2).eq.'-v')then
                verby = .true.
            else if(name(1:2).eq.'-?')then
                call usage()
            else if(name(1:2).eq.'-h')then
                call usage()
            endif
        go to 11
   13   continue
        return
        end
        
        subroutine usage()
        parameter (LER=6, LIN=5, LOT=6)
        write(LER,*)'sactof96 [-v] [-?] [-h]'
        write(LER,*)
     1  ' -v          Verbose output about the conversion'
        write(LER,*)
     1  ' -?          This online help'
        write(LER,*)
     1  ' -h          This online help'
        stop
        end

