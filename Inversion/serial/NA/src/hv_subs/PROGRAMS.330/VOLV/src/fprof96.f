        program fprof96
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME V                                                       c
c                                                                     c
c      PROGRAM: FPROF96                                               c
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
c       13 JUN 2003 add flag -S scale to adjust the heights
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
        logical shwabs
        logical doarriv
        integer kolor
        real tmin, tmax, twin
        real scltrc
c-----
c       internal variables
c-----
        integer LER, LIN, LOT
        parameter (LER=0, LIN=5, LOT=6)
        parameter (NSAMP=16384)
        dimension x(NSAMP),y(NSAMP)
        integer jjsrc(21)
        character svcomp*8
        character*4 SYMTP, SYMTS

        data jjsrc/21*0/
c-----
c       call machine dependent initialization
c-----
        call mchdep()
c-----
c       parse command line arguments
c-----
        call gcmdln(kolor,shwabs,tmin,tmax,twin,doarriv,scltrc)
c-----
c       graphic initialization
c-----
        call pinitf('FPROF96.PLT')
c-----
c       read in all data, place on a temporary file
c       which will be reread to plot each green's function in
c       a distance profile
c-----
        open(1,status='scratch',form='unformatted',
     1      access='sequential')
        rewind 1
        ndist = 0
c-----
c       determine which Greens functions computed
c-----
c-----
c       process Green's Functions
c-----
        jftype =0
  100   continue
                call rdhd96(LIN,nerr)
                if(nerr .lt. 0)go to 9999
                ndone = 0
                if(iftype.gt.jftype)jftype = iftype
                do 101 j=1,iftype
                    if(jsrc(j).ne.0 .and. 
     1              ndone.lt.iftype) then
                    jjsrc(j)=1
                    ndone = ndone + 1
                    endif
  101           continue
                write(1)distkm,stelev,evdep,jsrc,iftype,
     1              tp,tsv,tsh
                ndist = ndist + 1
                ndone = 0
            do 200 j=1,iftype
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
                    if(ndone.gt.iftype)go to 201
                    write(1)stcomp,cmpdt,npts,ssec
c-----
c       SUN compiler bug
c-----
c                   write(1)(y(i),i=1,npts)
                    do 9878 il=1,npts,512
                        iu = (il - 1) + 512
                        if(iu .gt. npts)iu = npts
                        write(1)(y(i),i=il,iu)
 9878               continue
  201           continue
                endif
  200       continue
        goto 100
 9999   continue
c-----
c       process each green's function to make a distance profile
c-----
        ifirst = 1
        do 1000 igrn = 1,jftype
            if(jjsrc(igrn).ne.0)then
            if(ifirst.eq.1)then
                ifirst = 0
            else
                call frame()
            endif
            rewind 1
            do 2000 id = 1, ndist
                read(1,end=2000)distkm,stelev,evdep,jsrc,iftype,
     1              tp,tsv,tsh
                ndone = 0
                do 2001 ksrc=1,iftype
                    if(jsrc(ksrc).ne.0 .and.
     1              ndone .lt.iftype)then

                    ndone = ndone + 1
                    read(1,end=2000)stcomp,cmpdt,npts,ssec
c-----
c       SUN compiler bug
c-----
c                   read(1)(y(i),i=1,npts)
                    do 9879 il=1,npts,512
                        iu = (il - 1) + 512
                        if(iu .gt. npts)iu = npts
                        read(1)(y(i),i=il,iu)
 9879               continue
                    else
                        do 9880 i=1,npts
                            y(i) = 0.0
 9880                   continue
                    endif
                    if(igrn.eq.ksrc.and.jsrc(ksrc).ne.0)then
                        svcomp = stcomp
                        scalex = 6.0
                        scaley = 6.0/(ndist+1)
                        x0 = 1.0
                        y0 = 1.0 + (ndist - id + 1 )
     1                      *scaley
                        do 1101 i=1,npts
                            x(i) = ssec+(i-1)*cmpdt
 1101                   continue
c-----
c       window the output
c-----
                        if(twin.gt.0.0 .and.
     1                      twin.le.(npts-1)*cmpdt)then
                            npout = twin/cmpdt + 1
                            if(npout.lt.2)npout=2
                            npb = 1
                            npe = npout
                        else
                            npout = npts
                            npb = 1
                            npe = npts
                            if(tmin.ge.0.0)then
                                npb=tmin/cmpdt+1
                                if(npb.lt.0)npb=1
                                if(npb.gt.npts)then
                                    nbp=npts-1
                                endif
                            endif
                            if(tmax.ge.0.0)then
                                npe=tmax/cmpdt+1
                                if(npe.lt.0)npb=2
                                if(npe.gt.npts)nb2=npts
                            endif
                            if(npe.gt.npb)then
                                npout = npe - npb + 1
                            endif
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
                        else if(stcomp(1:1).eq.'N')then
                            ttp = TP
                            tts = TSV
                            SYMTP='P   '
                            SYMTS='SV  '
                        else if(stcomp(1:1).eq.'E')then
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
                        call xyprof(x,y,x0,y0,npb,npe,
     1                      scalex,
     1                      scaley,distkm,stelev,
     2                      evdep,stcomp,
     2                      shwabs,pscl,id,kolor,
     3                      doarriv,ttp,tts,
     3                      SYMTP,SYMTS,scltrc)

                endif
 2001           continue
 2000       continue
            xmin = ssec + (npb-1)*cmpdt
            xmax = ssec + (npb-1)*cmpdt + (npout-1)*cmpdt
            call dolinx(x0,0.90,scalex,xmax,xmin,
     1          0.07,.false.,.true.,.true.,1,' ')
            call symbol(1.0,0.50,0.14,'SRC = ',0.0,+6)
            call symbol(999.,999.,0.14,svcomp,0.0,+4)
            call symbol(3.0,0.50,0.14,'DT  = ',0.0,+6)
            call number(999.,999.,0.14,cmpdt,0.0,+6)
            l = lgstr(ccomnt)
            if(l.gt.1)then
                call symbol(1.0,0.30,0.10,ccomnt(1:l),0.0,l)
            endif
            endif
 1000   continue
        close (1)
        call pend()
        end

        subroutine xyprof(x,y,x0,y0,npb,npe,scalex,scaley,rr,rz,sz,sym,
     1      shwabs,pscl,id,kolor,doarriv,ttp,tts,SYMTP,SYMTS,scltrc)
        dimension x(*),y(*)
        character*8 sym
        logical shwabs
        logical doarriv
        character*4 SYMTP, SYMTS
        integer LER, LIN, LOT
        parameter (LER=0, LIN=5, LOT=6)
c-----
c       x(i)    R*4 array of x values
c       y(i)    R*4 array of y values
c       npb     I*4 First point to plot
c       npe     I*4 Last point to plot
c       x0      R*4 absolute x coordinate for trace beginning
c       y0      R*4 absolute y coordinate for trace beginning
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
c       pscl    R*4
c       id  I*4
c       kolor   I*4
c       doarriv L   .true. annotate with P, SV and SH times
c       ttp R*4 P travel time if > -12345.0
c       tts R*4 S travel time if > -12345.0
c       SYMTP   Ch*4    P arrival symbol
c       SYMTS   Ch*4    S arrival symbol
c       scltrc  R   scale plotted trace by this factor
c-----
c get the extremes
        xmax = -1.0e+35
        xmin = - xmax
        ymin=xmin
        ymax=xmax
        do 100 i=npb,npe
        if(x(i).gt.xmax)xmax=x(i)
        if(x(i).lt.xmin)xmin=x(i)
        if(y(i).gt.ymax)ymax=y(i)
        if(y(i).lt.ymin)ymin=y(i)
  100   continue
        yymax=abs(ymax)
        if(abs(ymin).gt.yymax)yymax=abs(ymin)
c------
c       put in maximum amplitude at right of plot
c-----
        call number(x0+scalex+0.1+1.0,y0-0.05,0.10,sz,0.0,4)
        call number(x0+scalex+0.1+1.0,y0-0.20,0.10,rz,0.0,4)
c       call symbol(x0+scalex+0.1,y0+0.10,0.10,sym,0.0,+4)
        call number(x0+scalex+0.1,y0-0.05,0.10,yymax,0.0,1003)
        call number(x0+scalex+0.1,y0-0.20,0.10,rr,0.0,4)
c-----
c-----
c       UNIX OUTPUT NO CARRIAGE CONTROL
c-----
c       write(LER,'(a,3f10.2,e10.3,2f11.3)')sym,rr,rz,sz,yymax,ttp,tts
c-----
c       MICROSOFT OUTPUT  CARRIAGE CONTROL
c-----
        write(LER,'(1x,a,3f10.2,e10.3,2f11.3)')
     1      sym,rr,rz,sz,yymax,ttp,tts
c-----
        if(yymax.eq.0.0)yymax=1.0
        if(shwabs)then
            if(id.eq.1)then
                pscl = yymax
            else
                yymax = pscl
            endif
        endif
c-----
c       move the pen to the first point and plot
c-----
        ipen=3
        call newpen(kolor)
        do 200 i=npb,npe
            xx=x0 + scalex * (x(i)-xmin)/(xmax-xmin)
            yy=y0 + scltrc * scaley * y(i)/yymax
            if(i.gt.npb)ipen=2
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

        subroutine gcmdln(kolor,shwabs,tmin,tmax,twin,doarriv,scltrc)
c-----
c       parse command line arguments
c
c       requires subroutine targ() and function mnmarg()
c
c-----
c-----
c       kolor   - I*4   CALPLOT pen number default = 1
c       shwabs  - L first trace defines scaling for 
c                 each Green's function
c       twin    - R*4   Length of plot window in sec
c       tmin    - R*4   Start time relative to begin in sec 
c       tmax    - R*4   End   time relative to begin in sec
c       doarriv - L .true. Annotate with arrival times
c               default .false.
c       scltrc  - R scale plotted trace by this factor default = 1.0
c-----
        character*20 name
        integer kolor
        logical shwabs
        integer mnmarg
        logical doarriv
        real tmin, tmax, twin
        real scltrc

        kolor = 1
        nmarg = mnmarg()
        shwabs = .false.
        tmin = -1.0
        tmax = -1.0
        twin = -1.0
        doarriv = .false.
        scltrc = 1.0
        i = 0
   11   i = i + 1
        if(i.gt.nmarg)go to 13
            call mgtarg(i,name)
            if(name(1:2).eq.'-K')then
                i=i+1
                call mgtarg(i,name)
                read(name,'(bn,i10)')kolor
                if(kolor.lt.1)kolor = 1
            else if(name(1:2).eq.'-T' .and.
     1          name(1:3).ne.'-TM')then
                i=i+1
                call mgtarg(i,name)
                read(name,'(bn,f10.0)')twin
            else if(name(1:5).eq.'-TMIN')then
                i=i+1
                call mgtarg(i,name)
                read(name,'(bn,f10.0)')tmin
            else if(name(1:5).eq.'-TMAX')then
                i=i+1
                call mgtarg(i,name)
                read(name,'(bn,f10.0)')tmax
            else if(name(1:2).eq.'-S')then
                i=i+1
                call mgtarg(i,name)
                read(name,'(bn,f10.0)')scltrc
            else if(name(1:2).eq.'-R')then
                shwabs = .true.
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
        write(LER,*)'USAGE: fprof96 ',
     1  ' [-R] [-K kolor] [-TMAX tmax] [-TMIN tmin] [-A]',
     2  ' [-S scaletrace] [-?] [-h]'
        write(LER,*)
     1      ' -R         (default false) scale by first trace'
        write(LER,*)
     1      ' -K kolor   (default 1)  pen color'
        write(LER,*)
     1      ' -TMAX tmax (default all) end relative to start in sec'
        write(LER,*)
     1      ' -TMIN tmin (default all) begin relative to start in sec'
        write(LER,*)
     1      ' -S scaletrace (def  1.0) scale plotted trace factor'
        write(LER,*)
     1      ' -A         (default off) annotate with 1st arrival time'
        write(LER,*)
     1      ' -?  This help message'
        write(LER,*)
     1      ' -h  This help message'
        stop
        end

