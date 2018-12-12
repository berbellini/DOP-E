        program fplot96
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME V                                                       c
c                                                                     c
c      PROGRAM: FPLOT96                                               c
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
c-----
c       11 SEP 2000 add TP, TSV TSH to FILE96 format
c       28 SEP 2000 add option -A to annotate plots with
c               theoretical P and S times
c       19 NOV 2000 added newest dolinx from grphsubf.f
c       14 JAN 2001 add A, C, F, L, N and density information
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
c       command line arguments
c-----
        integer kolor
        real tmax
        logical doarriv
c-----
c       internal variables
c-----
        parameter (LER=0, LIN=5, LOT=6)
        parameter (NSAMP=16384)
        dimension x(NSAMP),y(NSAMP)
        character svcomp*8
        real*8 sepoch, eepoch, o

        character*4 SYMTP, SYMTS
c-----
c       call machine dependent initialization
c-----
        call mchdep()
c-----
c       parse command line arguments
c-----
        call gcmdln(kolor,tmax,doarriv)
c-----
c       graphic initialization
c-----
        call pinitf('FPLOT96.PLT')
c-----
c       read in data, and plot all traces of a given 'station'
c       on the same page
c-----
c-----
c       determine which Greens functions computed
c-----
c-----
c       process Green's Functions
c-----
        ipage = 0
  100   continue
            call rdhd96(LIN,nerr)
            if(nerr .lt. 0)go to 9999
            if(ipage .eq.0)then
                ipage = 1
            else
                call frame()
            endif
c-----
c       only plot upto iftype traces
c-----
            ndone = 0
            do 200 j=1,21
                do 300 i=1,NSAMP
                    y(i)=0.0
  300           continue
                if(jsrc(j) .ne. 0)then
                    call rdtr96(LIN,stcomp,cmpinc,cmpaz,
     1                  cmpdt, npts, ksyear, ksmon, 
     2                  ksday, kshour, ksmin, ssec, 
     3                  y,nerr,NSAMP)
                    if(nerr .lt. 0 )go to 9999
                    ndone = ndone + 1
                    if(ndone.gt.iftype)go to 200
c-----
c                   set time of first sample
c-----
                    if(ksyear.eq.0)ksyear = 1970
                    if(keyear.eq.0)keyear = 1970
                    if(ksmon.eq.0)ksmon = 1
                    if(kemon.eq.0)kemon = 1
                    if(ksday.eq.0)ksday = 1
                    if(keday.eq.0)keday = 1
                    call htoe(ksyear,ksmon,ksday,kshour,
     1                  ksmin,ssec,sepoch)
                    call htoe(keyear,kemon,keday,kehour,
     1                  kemin,esec,eepoch)
                    o = sepoch - eepoch
                    ti = sngl(o)
                    dt = cmpdt
                    svcomp = stcomp
                    scalex = 6.0
                    scaley = 6.0/(iftype+1)
                    x0 = 1.0
                    y0 = 1.25 + (iftype - ndone + 1 )
     1                  *scaley
                    do 1101 i=1,npts
                        x(i) = ti+(i-1)*cmpdt
 1101               continue
c-----
c       window the output
c-----
                    if(tmax.gt.0.0 .and.
     1                  tmax.le.(npts-1)*cmpdt)then
                        npout = tmax/cmpdt + 1
                        if(npout.lt.2)npout=2
                    else
                        npout = npts
                    endif
c-----
c       get correct timing information
c-----
                    if(iobsyn.eq.2)then
                        if(stcomp(1:1).eq.'Z')then
                            ttp = TP
                            tts = TSV
                            SYMTP='P   '
                            SYMTS='SV  '
                        else if(stcomp(1:1).eq.'R')then
                            ttp = TP
                            tts = TSV
                            SYMTP='P   '
                            SYMTS='SV  '
                        else if(stcomp(1:1).eq.'P')then
                            ttp = TP
                            tts = TSV
                            SYMTP='P   '
                            SYMTS='SV  '
                        else if(stcomp(1:1).eq.'T')then
                            ttp = -12345.
                            tts = TSH
                            SYMTP='    '
                            SYMTS='SH  '
                        endif
                    else
                        ttp = -12345.0
                        tts = -12345.0
                        SYMTP='    '
                        SYMTS='    '
                    endif
                        
                    call xyprof(x,y,x0,y0,npout,
     1                  scalex,
     1                  scaley,distkm,stelev,
     2                  evdep,stcomp,
     2                                  id,kolor,doarriv,ttp,tts,
     3                  SYMTP,SYMTS)
                endif
            if(ndone.eq.1)then
            xmin = ssec
            xmax = ssec + (npout-1)*cmpdt
            call dolinx(x0,1.10,scalex,xmax,xmin,
     1          0.07,.false.,.true.,.true.,1,' ')
            call symbol(1.0,0.8,0.10,'TMIN= ',0.0,6)
            call number(999.,999.,0.10,x(1),0.0,+3)
            call symbol(3.0,0.8,0.10,'TMAX= ',0.0,+6)
            call number(999.,999.,0.10,x(npts),0.0,+3)
            call symbol(5.0,0.80,0.10,'DT  = ',0.0,+6)
            call number(999.,999.,0.10,cmpdt,0.0,+4)
            call symbol(1.0,0.55,0.10,'DIST= ',0.0,+6)
            call number(999.,999.,0.10,distkm,0.0,+2)
            call symbol(3.0,0.55,0.10,'SDEPTH=',0.0,+7)
            call number(999.,999.,0.10,evdep,0.0,+2)
            call symbol(5.0,0.55,0.10,'RDEPTH=',0.0,+7)
            call number(999.,999.,0.10,stelev,0.0,+2)


            l = lgstr(ccomnt)
            if(l.gt.1)then
                call symbol(1.0,0.30,0.10,ccomnt(1:l),0.0,l)
            endif
            endif
  200       continue
        goto 100
 9999   continue
        call pend()
        end

        subroutine xyprof(x,y,x0,y0,n,scalex,scaley,rr,rz,sz,sym,
     1      id,kolor,doarriv,ttp,tts,SYMTP,SYMTS)
        dimension x(*),y(*)
        character*8 sym
        logical doarriv
        real ttp, tts
        character*4 SYMTP, SYMTS
c-----
c       x(i)    R*4 array of x values
c       y(i)    R*4 array of y values
c       n   I*4 number of (x,y) pairs to plot
c       x0  R*4 absolute x coordinate for trace beginning
c       y0  R*4 absolute y coordinate for trace beginning
c       scalex  R*4 length of plotted x axis in inches
c       scaley  R*4 length of plotted y axis in inches
c       rr  R*4 Receiver Radial distance in km
c       rz  R*4 Receiver position in km
c       sz  R*4 Source position in km
c       sym Ch*8    Component name
c       shwabs  L   .true. Traces shown realtive to first 
c                   which gives absolute amplitude plot
c       id  I*4 trace identification number
c       kolor   I*4 Trace pen color
c       doarriv L   .true. annotate with P, SV and SH times
c       ttp R*4 P travel time if > -12345.0
c       tts R*4 S travel time if > -12345.0
c       SYMTP   Ch*4    P arrival symbol
c       SYMTS   Ch*4    S arrival symbol
c-----
c get the extremes
        xmax = -1.0e+35
        xmin = - xmax
        ymin=xmin
        ymax=xmax
        do 100 i=1,n
        if(x(i).gt.xmax)xmax=x(i)
        if(x(i).lt.xmin)xmin=x(i)
        if(y(i).gt.ymax)ymax=y(i)
        if(y(i).lt.ymin)ymin=y(i)
  100   continue
        yymax=abs(ymax)
        if(abs(ymin).gt.yymax)yymax=abs(ymin)
        call symbol(x0+scalex+0.1,y0+0.10,0.10,sym,0.0,+4)
        call number(x0+scalex+0.1,y0-0.05,0.10,yymax,0.0,1003)
c-----
c-----
c       move the pen to the first point and plot
c-----
        ipen=3
        if(yymax.eq.0.0)return
        call newpen(kolor)
        do 200 i=1,n
            xx=x0 + scalex * (x(i)-xmin)/(xmax-xmin)
            yy=y0 + scaley * y(i)/yymax
            if(i.gt.1)ipen=2
            call plot(xx,yy,ipen)
  200   continue
        call newpen(1)
c-----
c       annotate with the P, S times
c-----
        if(doarriv)then
            if(ttp .gt. -12345.0)then
                xpos=x0 + scalex * (ttp-xmin)/(xmax-xmin)
                ypos0 = y0 + 0.2*scaley
                ypos1 = y0 + 0.5*scaley
                ypos2 = y0 + 0.6*scaley
                yht   = 0.20*scaley
                call plot(xpos,ypos0,3)
                call plot(xpos,ypos1,2)
                call symbol(xpos,ypos2,yht,SYMTP,0.0,2)
            endif
            if(tts .gt. -12345.0)then
                xpos=x0 + scalex * (tts-xmin)/(xmax-xmin)
                ypos0 = y0 + 0.2*scaley
                ypos1 = y0 + 0.5*scaley
                ypos2 = y0 + 0.6*scaley
                yht   = 0.20*scaley
                call plot(xpos,ypos0,3)
                call plot(xpos,ypos0,3)
                call plot(xpos,ypos1,2)
                call symbol(xpos,ypos2,yht,SYMTS,0.0,2)
            endif
        endif
c lift the pen
        call plot(xx,yy,3)
        return
        end

        subroutine gcmdln(kolor,tmax,doarriv)
c-----
c       parse command line arguments
c
c       requires subroutine targ() and funtion mnmarg()
c
c-----
c-----
c       kolor   - I*4   CALPLOT pen number default = 1
c       tmax    - R*4   Maximum time window length (default original)
c       doarriv - L .true. Annotate with arrival times
c               default .false.
c-----
        character*20 name
        integer kolor
        real*4 tmax
        logical doarriv
        integer mnmarg

        kolor = 1
        tmax = -1.0
        doarriv = .false.
        nmarg = mnmarg()
        i = 0
   11   i = i + 1
        if(i.gt.nmarg)go to 13
            call mgtarg(i,name)
            if(name(1:2).eq.'-K')then
                i=i+1
                call mgtarg(i,name)
                read(name,'(bn,i10)')kolor
                if(kolor.lt.1)kolor = 1
            else if(name(1:2).eq.'-T')then
                i=i+1
                call mgtarg(i,name)
                read(name,'(bn,f10.0)')tmax
            else if(name(1:2).eq.'-A')then
                doarriv = .true.
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
        integer*4 LER, LIN, LOT
        parameter(LER=0, LIN=5, LOT=6)
        write(LER,*)'USAGE: fplot96  ',
     1      '  [-K kolor] [-T tmax] [-A] [-?] [-h]'
        write(LER,*)
     1      ' -K kolor    (default 1)  pen color'
        write(LER,*)
     1      ' -T tmax  (default all) maximum length of plot in sec'
        write(LER,*)
     1      ' -A       (default off) annotate with 1st arrival time'
        write(LER,*)
     1      ' -?  This help message'
        write(LER,*)
     1      ' -h  This help message'
        stop
        end

        subroutine etoh(epoch,date,str,doy)
c-----
c       convert from epoch time to human time
c
c       epoch   - R*8 time in seconds relative to 0000 1 Jan 1970
c       date    - I*4 Julian date
c       str - C*  printable string
c       diy - I*4 dayin the year
c-----
        real*8 epoch
        integer*4 date
        character str*(*)
        integer*4 diy, doy, hour, minute, year, month
        integer*4 day
        real*4 second
        real*8 seclft

        str=' '
        doy = (epoch/86400.0d+00)
        seclft = dmod(epoch, 86400.0d+00)
        hour = 0
        minute = 0
        second = 0.00
c-----
c       compute hours minutes seconds
c-----
        if(seclft .ne. 0.00d+00)then
c-----
c                   before 1970 subtract and add a day
c-----
            if(seclft .lt. 0.0d+00)then
                doy = doy - 1
                seclft = seclft + 86400.00d+00
            endif
            hour = (seclft/3600.00d+00)
            seclft = dmod(seclft,3600.0d+00)
            minute = seclft/60.0d+00
            second = dmod(seclft,60.0d+00)
        endif

        if(doy .ge. 0)then
            year = 1970
 1000       continue
                diy =  leapdy(year)
                if(doy .lt. diy)go to 2000
                doy = doy - diy
                year = year + 1
            go to 1000
        else
            year = 1969
 1100       continue
                diy =  leapdy(year)
                doy = doy + diy
                if( doy .gt. 0 ) go to 2000
                year = year - 1
            go to 1100
        endif
 2000   continue
        doy = doy + 1
        date = year*1000 + doy
        call mnthdy(year,doy,month,day)
        write(str,110) year,month,day,hour,minute,second
  110   format(i4,i2,i2,i2,i2,f6.3)
c-----
c       guarantee that there are no blanks in the string str
c-----
        do 2100 i=1,17
            if(str(i:i).eq.' ')str(i:i)='0'
 2100   continue
        return
        end

        function leapdy(yr)
        integer*4 yr
        logical t1, t2, t3
        t1 = mod(yr,4).ne.0
        t2 = mod(yr,100).ne.0
        t3 = mod(yr,400).ne.0
        if( .not.t1 .and. t2)then
            isleap = 1
            leapdy = 366
        elseif( .not.t3)then
            isleap = 1
            leapdy = 366
        else
            isleap = 0
            leapdy = 365
        endif
        return
        end

        subroutine mnthdy(year,doy,month,day)
        integer*4 year, doy, month, day
        integer*4 i, dim, leap
        integer*4 dmnth(12)
        data dmnth/31,28,31,30,31,30,31,31,30,31,30,31/
        if(leapdy(year).eq.366)then
            leap = 1
        else
            leap = 0
        endif
        day = doy
        do 100 i=1,12
            month = i
            dim = dmnth(i)
            if(leap.eq.1 .and. i.eq.2)dim = dim + 1
            if(day .le.dim)goto 1000
            day = day - dim 
  100   continue
 1000   continue
        return
        end


        subroutine htoe(year,month,day,hour,minute,second,epoch)
c-----
c       convert calendar date to epoch time since January 1, 1970
c-----
c       year    - I*4   year
c       month   - I*4   month
c       day - I*4   day
c       hour    - I*4   hour
c       minute  - I*4   minute c    second  - I*4   second
c       second  - R*4   seconds
c       epoch   - R*8   time in seconds relative to 00:00 01 Jan 1970
c-----
        integer*4 year, month, day, hour, minute, date, diy
        real*4 second
        real*8 epoch, dtoepo
        integer*4 daymon(12)
        data daymon/0, 31, 59, 90, 120, 151, 181, 212, 243, 273,
     1      304, 334/
        diy = daymon(month) + day
        if(leapdy(year).eq.366 .and. month.gt.2)diy=diy+1
        date = 1000*year + diy
c       write(6,*)'date=',date
        epoch = dtoepo(date) + hour * 3600.0d+00 + 
     1      minute * 60.0d+00 +dble(second)
        return
        end

c-----
c       convert julian date to epoch time
c-----
        function dtoepo(date)
        real*8 dtoepo
        integer*4 date, diy, cnt, days

        cnt = date / 1000
        days = 0
        if (cnt .gt. 1970)then
 1000       continue
            cnt = cnt -1
            if(cnt.lt.1970)go to 2000
                days = days + leapdy(cnt)
            go to 1000
        else if (cnt .lt. 1970)then
 1100       continue
            if(cnt.ge.1970)goto 2000
                days = days - leapdy(cnt)
                cnt = cnt + 1
            go to 1100
        endif
 2000   continue
        diy = (date -1) / 1000
        diy = (date -1 ) -  1000*diy
c       write(6,*)'days=',days,' diy=',diy
        dtoepo = (days + diy) * 86400.0d+00
        return
        end

