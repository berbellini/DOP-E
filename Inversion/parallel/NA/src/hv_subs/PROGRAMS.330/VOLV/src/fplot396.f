        program fplot396
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME V                                                       c
c                                                                     c
c      PROGRAM: PLOT396                                               c
c                                                                     c
c      COPYRIGHT 1985 R. B. Herrmann                                  c
c                   96                                                c
c      Department of Earth and Atmospheric Sciences                   c
c      Saint Louis University                                         c
c      221 North Grand Boulevard                                      c
c      St. Louis, Missouri 63103                                      c
c      U. S. A.                                                       c
c                                                                     c
c------- -------------------------------------------------------------c
c-----
c       11 SEP 2000 add TP, TSV TSH to FILE96 format
c       14 JAN 2001 add A, C, F, L, N and density information
c-----
c-----
c       command line arguments
c-----
        character*80 cmdfil
c-----
c       internal variables
c-----
        integer LER, LIN, LOT
        parameter (LER=0, LIN=5, LOT=6)
        character*80 fname
        integer*4 first

        logical shwdep, shwfil, shwabc,  dovert, shwtit, dorel, shwtim
        character*80 cmd
        character*80 title
        logical ext
        character tmin*10, tmax*10

        common/scl/yscl(16), tscl
        data first/0/
c-----
c       call machine dependent initialization
c-----
        call mchdep()
c-----
c       parse command line arguments
c-----
        call gcmdln(cmdfil)
        if(cmdfil .ne. ' ')then
            inquire(file=cmdfil,exist=ext)
            if(ext)then
                LUN=1
                open(LUN,file=cmdfil,status='unknown',
     1              form='formatted',access='sequential')
                rewind LUN
            else
                write(LER,*)'Specified command file does not',
     1              ' exist'
                call usage()
            endif
        else
            LUN = LIN
        endif
c-----
c       This programs reads 3-component traces in FILE96 format
c       and plots them side by side. The file names are read
c       from a command list on standard input. The syntax is
c       
c       COMMAND
c           action
c       COMMAND
c           action
c       
c       where the action is optional according to the command.
c       The commands are
c       
c       NEWPAGE - which starts a new page
c       SHOWZRT - which puts Z R T at the top of the columns
c       SHOWZNE - which puts Z N E at the top of the columns
c       NOSHOWZRT   - do not put up Z N T or Z R T
c       NOSHOWZNE
c       SHOWABC - puts (a) (b) ... to the right of each trace
c       SHOWFILE    - shows the file name
c       NEWPEN  - get a new pen
c           kolor   integer representing color
c       NEWY0   - initial vertical position on page, 
c                 requires the following 
c               action line
c           y0  - new vertical position
c       RESETABC    - puts resets the counter at (a)
c       FILE    - which indicates a FILE96 to be plotted and 
c                 requires the 
c               following action line
c           'filename' start_time end_time 
c       VERTICAL - traces are aligned vertically on page
c       HORIZONTAL - traces are aligned horizontally on page
c       DOWN - Move the current position down this amount
c           downy
c       DELY - indicates the change in vertical position
c           dely
c       SHOWTITLE  - use this for a title instead of file name
c           title
c       NOSHOWTITLE - do not show title
c       RELATIVE    - gains
c       NORELATIVE
c       TIME    - show time scale at center, this is only done once
c           tmin tmax
c-----
c       initialize plot stream
c-----
        call pinitf('FPLOT396.PLT')
        x0 = 1.0
        y0 = 7.0

        nabc = 0
        iszrt = 0
        shwdep = .false.
        shwfil = .false.
        shwabc = .false.
        dovert = .false.
        shwtit = .false.
        shwtim = .false.
        title = ' '
        kolor = 1
        yy0 = y0
        dely = 0.75
        dorel = .false.
        tscl = -1.0
        icmd = 0
 9996   continue
            read(LUN,'(a)',end=9997,err=9997)cmd
c-----
c           parse the command
c-----
            if(cmd.eq.'NEWPAGE')then
                y0 = yy0
                first = 0
            else if(cmd.eq.'SHOWZRT')then
                iszrt = 1
            else if(cmd.eq.'SHOWZNE')then
                iszrt = 2
            else if(cmd(1:7).eq.'NOSHOWZ')then
                iszrt = 0
            else if(cmd(1:7).eq.'NOSHOWT')then
                shwtit = .false.
            else if(cmd(1:8).eq.'RELATIVE')then
                dorel = .true.
                icmd = 1
            else if(cmd(1:8).eq.'NORELATI')then
                dorel = .false.
                icmd = 0
            else if(cmd.eq.'SHOWABC')then
                shwabc = .true.
                nabc = 0
            else if(cmd(1:4).eq.'TIME')then
                shwtim = .true.
                read(LUN,*,end=9997,err=9997)tmin, tmax
                if(.not.dovert)then
                    call taxis(dovert,tmin,tmax,x0,y0)
                endif
            else if(cmd(1:6).eq.'SHOWDE')then
                shwdep = .true.
            else if(cmd(1:6).eq.'SHOWFI')then
                shwfil = .true.
            else if(cmd(1:6).eq.'SHOWTI')then
                read(LUN,*,end=9997,err=9997)title
                shwtit = .true.
            else if(cmd.eq.'NEWY0')then
                read(LUN,*,end=9997,err=9997)yy0
                y0 = yy0
            else if(cmd.eq.'NEWPEN')then
                read(LUN,*,end=9997,err=9997)kolor
                if(kolor.lt.0)kolor = 1
            else if(cmd.eq.'DOWN')then
                read(LUN,*,end=9997,err=9997)downy
                y0 = y0 - downy
            else if(cmd.eq.'DELY')then
                read(LUN,*,end=9997,err=9997)dely
            else if(cmd.eq.'RESETABC')then
                nabc = 0
            else if(cmd.eq.'FILE')then
                shwdep = .false.
                shwfil = .false.
                read(LUN,*)fname, ts, te
                call proces(dovert,fname,ts,
     1              te,shwdep,shwfil,dely,x0,y0,
     2              nabc, shwabc, iszrt, kolor,
     3              shwtit, title,yy0,first,dorel,
     4              icmd,shwtim,tmin,tmax)
            else if(cmd(1:4).eq.'VERT')then
                dovert = .true.
            else if(cmd(1:4).eq.'HORI')then
                dovert = .false.
            else
                lc = lgstr(cmd)
                write(LER,*)'parse error command:',cmd(1:lc)
            endif
        go to 9996
 9997   continue
        call pend()
        if(LUN.ne.LIN)close (LUN)
        end

        subroutine taxis(dovert,tmin,tmax,x0,y0)
c-----
c       plot the time axis at the current vertical position
c-----
c       dovert  L   - .true. plot traces vertically, 
c                   else side by side
c       x0  R*4 - 
c       y0  R*4 - trace begins at (x0, y0). Note these may
c                   be changed
c       tmin    C*10
c       tmax    C*10       time values to be annotated if shwtim=.true.
c       BEWARE - THIS ROUTINE KNOWS THINGS THAT ARE INTERNAL TO
c           PROCESS
c-----
        logical dovert
        character*10 tmin, tmax 
        real*4 x0, y0
c
        common/scl/yscl(16), tscl
c-----
c       tscl is the number of inches per second from the 
c           last trace plotted
c-----
c       set the x-axis length
c-----
            if(dovert)then
                xlen = 6.0
                xpos = x0 +xlen + 0.1
                xcen = x0 + 0.5*xlen
            else
                xlen = 2.5
                xpos = x0 + 5.6 + xlen + 0.1
                xcen = x0 + 2.8 + 0.5*xlen
            endif
c-----
c       convert from string time to real time
c-----
        read(tmin,'(bn,f10.0)',end=9999,err=9999)ttmin
        read(tmax,'(bn,f10.0)',end=9999,err=9999)ttmax
        if(ttmax .le. ttmin)return
c-----
c       put in the tic marks
c-----
        if(tscl .le. 0.0)return
            xln = (ttmax - ttmin)*tscl
            xlw = xcen - 0.5 * xln
            xhg = xcen + 0.5 * xln
            call plot(xlw,y0,3)

            call plot(xhg,y0,2)
            call plot(xhg,y0+0.2,2)
            call plot(xhg,y0+0.2,3)

            xx = xhg
            yy = y0 + 0.3
            ht = 0.10
            angle = 0.0
            call center(xx,yy,ht,tmax,nchar,angle)
            call symbol(xx,yy,ht,tmax,angle,nchar)
            
            call plot(xlw,y0,3)
            call plot(xlw,y0+0.2,2)
            call plot(xlw,y0+0.2,3)

            xx = xlw
            yy = y0 + 0.3
            ht = 0.10
            angle = 0.0
            call center(xx,yy,ht,tmin,nchar,angle)
            call symbol(xx,yy,ht,tmin,angle,nchar)
            
 9999   continue
        return
        end

        subroutine center(xx,yy,ht,string,il,angle)
            character string*(*)
            il = lgstr(string)
            t = angle*3.1415927/180.0
            ct = cos(t)
            st = sin(t)
            rl = 0.5*il*ht
            
            xx = xx - rl*ct
            yy = yy - rl*st
        return
        end

        subroutine right(xx,yy,ht,string,il,angle)
            character string*(*)
            il = lgstr(string)
            t = angle*3.1415927/180.0
            ct = cos(t)
            st = sin(t)
            rl = 1.0*il*ht
            
            xx = xx - rl*ct
            yy = yy - rl*st
        return
        end
c-----

        subroutine newpag(iszrt,dovert,x0,y0)
c-----
c       annotate the top of the page
c-----
        integer*4 iszrt, first
        logical dovert
        real*4 x0, y0
            if(.not.dovert)then
                if(iszrt.eq.1)then
                    j=1
                    xx = 1.0 + 2.8*(j-1)
                    call symbol(xx+1.4-0.3*0.25,y0+0.7,
     1                  0.25,'Z',0.0,1)
                    j=2
                    xx = 1.0 + 2.8*(j-1)
                    call symbol(xx+1.4-0.3*0.25,y0+0.7,
     1                  0.25,'R',0.0,1)
                    j=3
                    xx = 1.0 + 2.8*(j-1)
                    call symbol(xx+1.4-0.3*0.25,y0+0.7,
     1                  0.25,'T',0.0,1)
                else if(iszrt.eq.2)then
                    j=1
                    xx = 1.0 + 2.8*(j-1)
                    call symbol(xx+1.4-0.3*0.25,y0+0.7,
     1                  0.25,'Z',0.0,1)
                    j=2
                    xx = 1.0 + 2.8*(j-1)
                    call symbol(xx+1.4-0.3*0.25,y0+0.7,
     1                  0.25,'N',0.0,1)
                    j=3
                    xx = 1.0 + 2.8*(j-1)
                    call symbol(xx+1.4-0.3*0.25,y0+0.7,
     1                  0.25,'E',0.0,1)
                endif
            endif
            first = 1
        return
        end
        

        subroutine proces(dovert,fname,ts,te,shwdep,shwfil,dely,x0,y0,
     1          nabc, shwabc, iszrt, kolor, shwtit, title,yyy0,
     2          first,dorel,icmd,shwtim,tmin,tmax)

        integer LER, LIN, LOT
        parameter (LER=0, LIN=5, LOT=6)
c-----
c       set up the plotting
c
c       dovert  L   - .true. plot traces vertically, 
c                   else side by side
c       fname   C*80    - name of file96(V) file
c       ts  R*4 - start time after origin in seconds
c                   if -999.0 then plot original trace
c       te  R*4 - end time to be plotted
c                   a total of (te - ts) 
c                   seconds will be plotted
c       shwdep  L   - .true. show the depth as an annotation
c       shwfil  L   - .true. show the file name
c       dely    R*4 - move the trace down dely units after plotting
c       x0  R*4 - 
c       y0  R*4 - trace begins at (x0, y0). Note these may
c                   be changed
c       nabc    I*4 - value of (a) (nabc=1), (b) (nabc=2) to be shown
c                   incremented here
c       shwabc  L   - .true. annotate with (a) ...
c       iszrt   I*4 - 0 do not put Z R T or Z N E at top of columns
c                 1 put ZRT
c                 2 put ZNE
c       kolor   I*4 - integer indicating color for trace
c       shwtit  L   - .true. show title, overrides shwfil
c       title   C*80    - title to be written at right of trace
c       yyy0    R*4 - initial vertical page position
c       first   I*4 - 0 = first page and do not call an 
c                   initial frame
c       dorel   L   - .true. all first set of traces in file
c                   set up scaling, others use
c       icmd    I*4 - 0 individual trace scaling
c                 1 set scaling by first of each column
c                 2 use previous scaling
c       shwtim  L   - .true. place time scale at center of trace
c       tmin    R*4
c       tmax    R*4    time values to be annotated if shwtim=.true.
c-----
c       subroutine argument parameters
c-----
        logical dovert, shwdep, shwfil, shwabc, shwtit, dorel, shwtim
        character fname*80, title*80
        real*4 ts, te
        real*4 dely, x0, y0, yyy0
        integer*4 nabc, kolor
        integer*4 first, icmd
        character tmin*10, tmax*10
c-----
c       variables in common blocks
c
c       iftype  I*4 File type
c               1 - single trace
c               3 - three component
c               16 - Green's function
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
c       internal variables
c-----
        integer NSAMP
        parameter (NSAMP=16384)
        real*4  y(NSAMP)
        logical ext
        real*8 sepoch, eepoch, o
        character*10 ostr
        integer*4 once
c-----
c       begin processing
c-----
C       y0 = y0 - dely
        if(first.eq.0)then
            call newpag(iszrt,dovert,x0,y0)
            first = 1
        endif
        if(fname .ne. ' ')then
            inquire(file=fname,exist=ext)
            if(ext)then
                LUN=2
                open(LUN,file=fname,status='unknown',
     1              form='formatted',access='sequential')
                rewind LUN
            else
                lfname = lgstr(fname)
                write(LER,*)'Specified file file does not',
     1              ' exist:', fname(1:lfname)
                return
            endif
        endif
c-----
c       set the x-axis length
c-----
        if(dovert)then
            xlen = 6.0
            xpos = x0 +xlen + 0.1
            xcen = x0 + 0.5*xlen
        else
            xlen = 2.5
            xpos = x0 + 5.6 + xlen + 0.1
            xcen = x0 + 2.8 + 0.5*xlen
        endif
c-----
c       open the FILE96 file
c-----
        once = 0

 1000   continue
            if(dovert .and. once.gt.0)then
                call frame()
                call newpag(iszrt,dovert,x0,y0)
                first = 1
            endif
            call rdhd96(LUN,nerr)
            if(iftype.ne.3)return
            if(nerr.lt.0)go to 9000
            do 2000 j=1,16
                if(jsrc(j).ne.0)then
                    call rdtr96(LUN,stcomp,cmpinc,cmpaz,
     1                  cmpdt, npts, ksyear, ksmon, 
     2                  ksday, kshour, ksmin, ssec, 
     3                  y,nerr,NSAMP)
                    if(nerr .lt. 0 )go to 9000
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
                    jjsrc = jsrc (j)
                    if(jjsrc.eq.1)then
                        ostr = 'Z'
c-----
c                       Vertical
c-----
                        if(iszrt.ne.1 
     1                  .and. iszrt.ne.2)go to 2000
                        if(dovert)then
                            xx0 = x0
                            yy0 = y0 
                        else
                            xx0 = x0
                            yy0 = y0 
                        endif
                    else if(jjsrc.eq.2)then
                        ostr = 'N'
c-----
c                       North
c-----
                        if(iszrt.ne.2)go to 2000
                        if(dovert)then
                            xx0 = x0
                            yy0 = y0 -1.5 
                        else
                            xx0 = x0 + 2.8 
                            yy0 = y0 
                        endif
                    else if(jjsrc.eq.3)then
                        ostr = 'E'
c-----
c                       East
c-----
                        if(iszrt.ne.2)go to 2000
                        if(dovert)then
                            xx0 = x0
                            yy0 = y0 - 3.0 
                        else
                            xx0 = x0 + 5.6
                            yy0 = y0 
                        endif
                    else if(jjsrc.eq.4)then
                        ostr = 'R'
c-----
c                       Radial
c-----
                        if(iszrt.ne.1)go to 2000
                        if(dovert)then
                        else
                            xx0 = x0 + 2.8
                            yy0 = y0
                        endif
                    else if(jjsrc.eq.5)then
                        ostr = 'T'
c-----
c                       Transverse
c-----
                        if(iszrt.ne.1)go to 2000
                        if(dovert)then
                        else
                            xx0 = x0 + 5.6
                            yy0 = y0
                        endif
                    endif
c-----
c                   vertical
c-----
                    call showit(xlen,xx0,yy0,y,npts,dt,
     1                  ti,ts,te,icmd,j,
     2                  fname,shwfil,title,shwtit,kolor)
                    if(dovert .and. iszrt.gt.0)then
                        call symbol(x0-0.3,yy0-0.07,
     1                      0.14,ostr,0.0,1)
                    endif
                endif
 2000       continue
            if(icmd.eq.1)icmd = 2
C           if(shwabc .and. first .eq. 0)then
C           call symbol(x0+1.4-0.3*0.25,y0+0.7,0.25,icom(j)(2:2),0.0,1)
            if(j.eq.3)first = 1
            if(shwtit)then
                ls = lgstr(title)
                call symbol(xpos-ls*0.15,y0+0.50,
     1              0.10,title,0.0,ls)
            endif
            if(shwdep)then
            call number(x0+3.4,y0+0.20,0.10,evdep,0.0,2)
            endif
            lf = lgstr(fname)
            if(shwfil .and. .not. shwtit)then
                call symbol(xpos+0.4,y0-0.20,0.07,
     1                  fname(1:lf),0.0,lf)
            endif
            if(shwabc)then
                nabc = nabc + 1
                ostr= '('//char(ichar('a')+nabc-1)//')'
                call symbol(xpos+0.2,y0,0.10,ostr,0.0,3)
            endif
            y0 = y0 - dely
            if(y0.lt.1.0)then
                y0 = yyy0
                call frame()
                call newpag(iszrt,dovert,x0,y0)
                first = 1
            endif
            once = once + 1
        go to 1000
 9000   continue
        return
        end

        subroutine showit(xlen,x0,y0,y,nt,dt,ti,ts,te,icmd,jj,
     1          fname,shwfil,title,shwtit,kolor)
c-----
c       plot a trace and annotate it
c-----
c       xlen    R*4 - length of trace
c       x0  R*4 - absolute position of beginning of trace
c       y0  R*4 - absolute position of beginning of trace
c       y() R*4 - array to be plotted
c       nt  I*4 - number of points to be plotted
c       ti  R*4 - time of first sample relative to origin time
c       ts  R*4 - start time of plot
c       te  R*4 - end time of plot
c       icmd    I*4 - 1 use current traces for scaling, and save
c                   so that other Z R T to same scale as
c                   first Z R T , respectively
c                 2 use previously determined relative scaling
c                 0 trace by trace scaling
c       fname   C*80    - file name to be annotated
c       shwfil  L   - show file name
c       title   C*80    - title name to be annotated
c       shwfil  L   - show title 
c                   overrides shwfil
c       kolor   I*4 - integer for trace color   
c-----
c       subroutine arguments
c-----
        logical shwfil, shwtit
        real*4 xlen, x0, y0
        integer NSAMP
        parameter (NSAMP=16384)
        real*4 y(NSAMP)
        integer*4 nt
        real *4 ti, ts, te
        character fname*80, title*80
        
c-----

        common/scl/yscl(16), tscl
        ymax = 0.0
c-----
c       plot the time series in the window of ts, te
c       this means that part of the time series may be plotted,
c       or even all of it but with blank spots at front or end
c-----
c-----
c       get pointers to first and last points of trace to be plotted
c-----
        if(te.le.ts)return
        if(te.lt.ti)return
        tend = ti + (nt-1)*dt
        if(ts.gt.tend)return
c-----
c       some of the trace must be able to be plotted
c-----
        call gtpnt(ti,ts,nt,dt,nti)
        call gtpnt(ti,te,nt,dt,nte)

        
        tl = te -ts
        if(tl.le.0.0)tl = 1.0
        xfac = xlen/tl
        tscl = xfac
c-----
c       ts, te is the plot window
c       
c-----
        ymax = 0.0
        do 1000 i=nti,nte
            t = ti + (i-1)*dt
            xx0 = x0 + (t - ts)* xfac
            ay = abs(y(i))
            if(ay .gt. ymax)ymax = ay
 1000   continue
        if(ymax.eq.0.0)then
            aay = 1.0
        else
            aay = ymax
        endif
        ylen = 0.50
        if(icmd.eq.1)then
            yfac = ylen/aay
            yscl(jj) = yfac
        else if(icmd.eq.2)then
            yfac = yscl(jj)
        else if(icmd.eq.0)then
            yfac = ylen/aay
        endif
        ipen = 3
        call newpen(kolor)
        do 2000 i=nti,nte
            t = ti + (i-1)*dt
            xx = x0 + (t - ts)* xfac
            yy = y0 + y(i)*yfac
            call plot(xx,yy,ipen)
            if(ipen.eq.3)then
                ipen = 2
            endif
 2000   continue
        call newpen(1)
            call plot(xx,yy,ipen)
            if(ipen.eq.3)then
                ipen = 2
            endif
        call number(x0+0.76*xlen,y0+0.25,0.10,ymax,0.0,1003)
        return
        end

        subroutine gtpnt(ti,te,nt,dt,nte)
c-----
c       determine array index corresponding to the time te
c-----
c       ti  R*4 - time of first sample
c       te  R*4 - desired time
c       nt  I*4 - number of samples in time series
c       dt  R*4 - sample interval
c       nte I*4 - pointer to the time
c-----
        real*4 ti, te, dt
        integer*4 nt, nte
            nte = (te - ti)/dt + 1
            if(nte.lt.1)then
                nte = 1
            else if(nte.gt.nt)then
                nte = nt
            endif
        return
        end

        subroutine gcmdln(cmdfil)
c-----
c       cmdfil  C*80    - name of command file, else stdin
c-----
        character*80 cmdfil
        character*80 name
        cmdfil = ' '
        nmarg = mnmarg()
        i= 0
   11   continue
            i = i + 1
            if(i.gt.nmarg)go to 13
            call mgtarg(i,name)
            if(name(1:2).eq.'-C')then
                i=i+1
                call mgtarg(i,cmdfil)
            else if(name(1:2).eq.'-?')then
                    call usage()
            else if(name(1:2).eq.'-h')then
                    call usage()
            endif
        goto 11
   13   continue
        return
        end

        subroutine usage()
        integer*4 LER, LIN, LOT
        parameter(LER=0, LIN=5, LOT=6)
        write(LER,*)'USAGE: fplot396 [-C cmdfil] [-?] [-h]'
        write(LER,*)
     1      ' -C cmdfil  command file, else standard input'
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
