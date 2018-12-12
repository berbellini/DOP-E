        program sac2loc
c-----
c      Read names of SAC files from the command line
c      Process each by opening to read the header values,
c      searching for P and S arrival time picks in the A and T0 
c      header values which are set by hitting 'P' and 'S' under ppk
c      
c      The output file is called elocate.dat
c
c      Then the program elocate can be run
c-----
c       CHANGES
c       26 MAY 2010 - modified the output part that handles
c          sac/sac2000 phase naming, e.g.,
c          [EI][PS][UD+- X][01234]
c-----
        integer LID
        parameter(LID=4)

        integer NFILEMX
        parameter (NFILEMX=5000)
        character*132 sfile(NFILEMX)

        integer nsfile

        integer NDATA
        parameter(NDATA=2000000)
        real datarr(NDATA)

        integer nerr
        integer nzyear, nzjday, nzhour, nzmin, nzsec, nzmsec
        real a, t0, stla, stlo, stel
        real tv(9)
        character*8 kstnm, kcmpnm
        character*8 ktv(9)
        integer i, j

        integer arid
        character*8 kpha

        character *8 kv(9), kkv(9)
        character *8 ka, kt0
        data kv/'T1      ','T2      ','T3      ','T4      ',
     1      'T5      ','T6      ','T7      ','T8      ',
     2      'T9      '/
        data kkv/'KT1     ','KT2     ','KT3     ','KT4     ',
     1       'KT5     ','KT6     ','KT7     ','KT8     ',
     2       'KT9     '/

        

c-----
c      read the command line to get the name of the sac files
c-----

        call gcmdln(nsfile,sfile,NFILEMX)
        if(nsfile.eq.0)call usage()
c-----
c      open the output file
c-----
        open(LID,file='elocate.dat',access='sequential',
     1      form='formatted',status='unknown')
        rewind LID
c-----
c      process the input files
c-----
        arid = 0
        do 1000 i=1,nsfile
            call brsac(1,NDATA,sfile(i),datarr,nerr)
            if(nerr.ge.0)then
                call getkhv('KSTNM   ',kstnm   ,nerr)
                call getkhv('KCMPNM  ',kcmpnm  ,nerr)
                call getfhv('STLA    ',stla    ,nerr)
                call getfhv('STLO    ',stlo    ,nerr)
                call getfhv('STEL    ',stel    ,nerr)
                call getfhv('A       ',a       ,nerr)
                call getkhv('KA      ',ka      ,nerr)
                call getfhv('T0      ',t0      ,nerr)
                call getkhv('KT0     ',kt0     ,nerr)
                call getnhv('NZYEAR  ',nzyear  ,nerr)
                call getnhv('NZJDAY  ',nzjday  ,nerr)
                call getnhv('NZHOUR  ',nzhour  ,nerr)
                call getnhv('NZMIN   ',nzmin   ,nerr)
                call getnhv('NZSEC   ',nzsec   ,nerr)
                call getnhv('NZMSEC  ',nzmsec  ,nerr)
                do 1001 j=1,9
                call getfhv( kv(j),tv(j)      ,nerr)
                call getkhv(kkv(j),ktv(j)     ,nerr)
 1001           continue
c-----
c          safety check
c-----      
            if(stla.gt.-12345. .and. stlo.gt.-12345.)then
                if(stel.eq. -12345.)then
                    stel =0
                endif
c-----
c      now convert to CSS time NOTE THIS MUST BE CHECKED AGAINST
c      SAC DEFINITIONS OF REFERENCE TIME 
c-----
                call mnthdy(nzyear,nzjday,nzmon,nzday) 
                if(a .gt. -12345.0)then
                    arid = arid + 1
                    call output(LID,kstnm,nzyear,nzmon,
     1                  nzday,nzhour,nzmin,nzsec,nzmsec,
     1                  ka ,a,stla,stlo,stel,arid)
                endif
                if(t0 .gt. -12345.0)then
                    arid = arid + 1
                    call output(LID,kstnm,nzyear,nzmon,
     1                  nzday,nzhour,nzmin,nzsec,nzmsec,
     1                  kt0,t0,stla,stlo,stel,arid)
                endif
                do 1002 j=1,9
                if(tv(j) .gt. -12345.0)then
                    arid = arid + 1
                    kpha = ktv(j)
                    call output(LID,kstnm,nzyear,nzmon,
     1                  nzday,nzhour,nzmin,nzsec,nzmsec,
     1                  kpha,tv(j),stla,stlo,stel,arid)
                endif
 1002           continue

            endif
            endif
 1000   continue
        close (LID)
        end


        subroutine output(LID,kstnm,nzyear,nzmon,nzday,
     1      nzhour,nzmin,nzsec,nzmsec,
     1      kpha,atime,slat,slon,stelv,arid)
        character kstnm*8, kpha*8
        integer LID
        integer nzyear,nzmon,nzday,nzhour,nzmin,nzsec,nzmsec
        real atime
        real slat,slon,stelv

c------
c      output parameters
c-----
        integer*4 asid, arid, mapd
        character*6 sta
        integer*4 ipwt
        real*8 tp
        real*4 xlat,xlon
        character*2 pfm
        character*8 pph
        character*1 pqual
        character*2 pchan
        real*4 elev
        integer*4 ondate, offdat

        character str*34
        integer*4 date
        integer*4 doy, hour, minute, year, month
        integer*4 day
        real*4 second
c------
c      conversions and defaults
c-----
        sta = kstnm(1:6)
        pchan = 'Z'
        if(kpha(1:1).eq.'e')then
            pqual = 'e'
            ipwt = 2
            call getpfm(kpha,pph,pfm)
        else if(kpha(1:1).eq.'i')then
            pqual = 'i'
            ipwt = 0
            call getpfm(kpha,pph,pfm)
c-----
c      default handling for SAC SAC2000
c-----
        else if(kpha(1:1).eq.'E')then
            pph = kpha(2:2)
            pfm = kpha(3:3)
            pqual = 'e'
            if(kpha(4:4).eq.'0')then
                ipwt = 0
            else if(kpha(4:4).eq.'1')then
                ipwt = 1
            else if(kpha(4:4).eq.'2')then
                ipwt = 2
            else if(kpha(4:4).eq.'3')then
                ipwt = 3
            else if(kpha(4:4).eq.'4')then
                ipwt = 4
            endif
        else if(kpha(1:1).eq.'I')then
            pqual = 'i'
            pph = kpha(2:2)
            pfm = kpha(3:3)
            if(kpha(4:4).eq.'0')then
                ipwt = 0
            else if(kpha(4:4).eq.'1')then
                ipwt = 1
            else if(kpha(4:4).eq.'2')then
                ipwt = 2
            else if(kpha(4:4).eq.'3')then
                ipwt = 3
            else if(kpha(4:4).eq.'4')then
                ipwt = 4
            endif
        else
            ipwt = 2
            pqual = 'e'
            pph = kpha(1:1)
            pfm = 'X'
        endif
            
        asid = arid
        xlat = slat
        xlon = slon
        elev = stelv
        ondate = 0000000
        offdat = 9999999
        
        second = nzsec + 0.001*nzmsec
        call htoe(nzyear,nzmon,nzday,nzhour,nzmin,second,tp)
        tp = tp + atime
        call etoh(tp,date,str,doy,
     1      year,month,day,hour,minute,second)

        write(LID,1)sta, pchan, 
     1  year, month, day,hour,minute,second,
     1  arid, pqual,
     2  pph, ipwt, pfm,asid,xlat, xlon,
     3  elev, ondate, offdat
    1   format(a6,1x,a2,1x,
     1  i4.4,1x,i2.2,1x,i2.2,1x,i2.2,1x,i2.2,1x,f6.3,1x,
     1  i8,1x,a1,1x,a8,1x,i2,1x,a2,1x,i8,
     1  1x,f9.4,1x,f9.4,1x,f9.0,1x,i8,1x,i8)

        return
        end
        
        subroutine getpfm(str,phase,pfm)
        character str*8, phase*8, pfm*2
c------
c      for historical reasons we must handle IPU0
c      for GSAC look for _
c-----
        integer ipos
        phase = ' '
        if(str(1:1).eq.'i' .or. str(1:1).eq.'e')then
            ipos = index(str,'_')
            if(ipos.eq. 0)then
                pfm = 'X '
                phase(1:7) = str(2:8)
            else if(ipos.ge.1.and.ipos.lt.8)then
                pfm = ' '
                pfm(1:1) = str(ipos+1:ipos+1)
                phase(1:ipos-2) = str(2:ipos-1)
            endif
        else
            phase = str
            pfm = 'X '
        endif
        return
        end
            

        subroutine gcmdln(nsfile,sfile,NFILEMX)
c-----
c      parse command line arguments
c----
c      nsfile  I   - number of file names
c      sfile   Ch  - array of file names
c      NFILEMX I   - dimension of file name array
c-----
        integer nsfile
        character*132 sfile(NFILEMX)

        integer mnmarg
        integer i, nmarg
        character*132 name
c-----
c      initialize
c-----
        nsfile = 0
        i = 0
        nmarg = mnmarg()

 1000   continue
            i = i + 1
            if(i.gt.nmarg)go to 2000
            call mgtarg(i,name)
            if(name(1:2).eq.'-h')then
                call usage()
            else if(name(1:2).eq.'-?')then
                call usage()
            else
                nsfile = nsfile + 1
                sfile(nsfile) = name
            endif
        go to 1000
 2000   continue
        return
        end

        subroutine usage()
        integer LER
        parameter (LER=0)
        write(LER,*)'Usage: sac2loc [-h] [-?] sacfiles'
        write(LER,*)
     1  ' sacfiles    - names of sacfiles to process',
     2  ' This can be done as  sac2loc *.sac under UNIX'
        write(LER,*)
     1  '-?         (default none )  this help message '
        write(LER,*)
     1  '-h         (default none )  this help message '
        stop
        end
        
