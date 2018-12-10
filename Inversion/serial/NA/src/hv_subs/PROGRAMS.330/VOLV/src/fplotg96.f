        program fplotg96
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME V                                                       c
c                                                                     c
c      PROGRAM: FPLOTG96                                              c
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
        logical shwabs
c-----
c       internal variables
c-----
        parameter (LER=0, LIN=5, LOT=6)
        parameter (NSAMP=16388)
        common/xx/x
        common/yy/y
        real*4 x(NSAMP),y(NSAMP)
        integer*4 kk(2)
        character*50 name(2)
        integer jjsrc(21,2)
        integer first, numgrn
        real*4 xr(2), hr(2), hs(2)
        first = 0
c-----
c       call machine dependent initialization
c-----
        call mchdep()
c-----
c       parse command line
c-----
        call gcmdln(name,jjplt,shwabs,kolor)
c-----
c       initilize open files
c-----
        if(jjplt.gt.0)then
c-----
c       file name or names from command line
c-----
            jplt = jjplt
            do 102 i=1,jplt
                kk(i) = i
                open(kk(i),file=name(i),status='old',
     1          form='formatted',access='sequential')
                rewind kk(i)
  102       continue
        else
c-----
c       single file input from standard input
c-----
            kk(1) = LIN
            jplt = 1
        endif
c-----
c       initialize plot stream
c-----
        call pinitf('FPLOTG96.PLT')
c-----
c       process Green's Functions
c-----
c-----
c       determine which Greens functions computed
c-----
  100   continue
            do 9997 k=1,jplt
                call rdhd96(kk(k),nerr)
                if(nerr .lt. 0)go to 9999
                ndone = 0
                do 101 j=1,21
                    if(jsrc(j).ne.0 .and. 
     1                  ndone.lt.iftype)then
                            jjsrc(j,k)=1
                        ndone = ndone + 1
                    endif
  101           continue
                xr(k) = distkm
                hs(k) = evdep
                hr(k) = stelev
 9997       continue
            if(first.gt.0)call frame()
c-----
c       process Green's Functions
c-----
            numgrn = 0
            do 200 j=1,21
                do 9996 k=1,jplt
                    if(jjsrc(j,k).ne.0)then
                    call rdtr96(kk(k),stcomp,cmpinc,cmpaz,
     1                  cmpdt, npts, ksyear, ksmon, 
     2                  ksday, kshour, ksmin, ssec, 
     3                  y,nerr,NSAMP)
                    if(nerr .lt. 0 )go to 9999
                    do 103 i=1,npts
                        x(i)=ssec + (i-1)*cmpdt
  103               continue
                    numgrn = numgrn +1
                    call seispl(jplt,k,x,y,npts,xr(k),hr(k),
     1                  j,stcomp,hs(k),cmpdt,shwabs,
     1                  ymxx,kolor,numgrn)
                    endif
 9996           continue
  200       continue
            first = first + 1
        go to 100
 9999   continue
        call pend()
        if(jjplt.gt.1)then
            do 9990 k=1,jjplt
                close (kk(k) )
 9990       continue
        endif
        end

        subroutine seispl(jplt,kk,x,y,n,dist,depthr,jgrn,sym,
     1      depths,dt,shwabs,ymxx,
     2                  kolor,numgrn)
        parameter (LER=0,LIN=5, LOT=6)
        character*8 sym
        dimension x(*),y(*)
        parameter (NSAMP=16388)
        dimension xx(NSAMP),yy(NSAMP)
        logical shwabs
        data xlen/2.37/,ylen/1.5/
c-----
c       set positioning for next Green;'s function
c-----
        ii = (jgrn + 1)/2
        irow = mod(ii,4) 
        if(irow .eq.0)irow = 4
        if(ii .gt. 4)then
            icol2 = 3
        else
            icol2 = 1
        endif
        icol = (1-mod(jgrn,2)) + icol2
        x0 = (icol-1) * xlen + 0.5
        y0 = (4-irow) * ylen + 1.0
c-----
c       set up scaling
c-----
        if(jplt.eq.2)then
            if(kk.eq.1)then
                ypt=0.6*ylen
                ysc = 0.38*ylen
                ysn = 0.9*ylen
            else if(kk.eq.2)then
                ypt=0.4*ylen
                ysc = 0.38*ylen
                ysn = 0.15*ylen
            endif
        else if(jplt.eq.1)then
            ysc= 0.45*ylen
            ysn = 0.25*ylen
            ypt = 0.5*ylen
        endif
c-----
c       get extreme values for plot scaling
c-----
        ymin = 1.0e+38
        xmin = 1.0e+38
        xmax =-1.0e+38
        ymax = -1.0e+38
        do 100 i = 1,n
            if(y(i).gt.ymax) ymax = y(i)
            if(x(i).gt.xmax) xmax = x(i)
            if(x(i).lt.xmin) xmin = x(i)
            if(y(i).lt.ymin) ymin = y(i)
  100   continue
        if(shwabs)then
            ysyn = 0.9*ylen
            if(kk.eq.1)then
                ymxx=ymax
                if(abs(ymin).gt.ymax) ymxx = abs(ymin)
                if(ymax.eq.ymin.or.ymxx.eq.0.0)ymxx=1.0
            endif
        else
            ymxx=ymax
            if(abs(ymin).gt.ymax) ymxx = abs(ymin)
            if(ymax.eq.ymin.or.ymxx.eq.0.0)ymxx=1.0
        endif
c-----
c       output max and min values of component
c-----
        write(LOT,8)kk,sym,ymax,ymin,ymxx
c-----
c       MSDOS
c-----
    8   format(1x,i3,1x,a,' ymax=',e10.3,
     1      ' ymin=',e10.3,' ymxx=',e10.3)
c-----
c       put in bounding box
c-----
        if(kk.eq.1)then
            xx(1)=x0
            xx(2)=xlen+x0
            xx(3)=xlen+x0
            xx(4)=x0
            xx(5)=x0
            yy(1)=y0
            yy(2)=y0
            yy(3)=ylen+y0
            yy(4)=ylen+y0
            yy(5)=y0
            xx(6)=0.0
            xx(7)=1.0
            yy(6)=0.0
            yy(7)=1.0
            call line(xx,yy,5,1,0,0)
        endif
c-----
c       scale values for plot
c-----
        do 200 i=1,n
            xx(i)=xlen*(x(i)-xmin)/(xmax-xmin)
            xx(i)=xx(i)+x0
            yy(i)=ysc*y(i)/ymxx + ypt
            yy(i)=yy(i)+y0
  200   continue
        xx(n+1)=0.0
        xx(n+2)=1.0
        yy(n+1)=0.0
        yy(n+2)=1.0
        if(shwabs.and.kk.eq.2)then
            call newpen(2)
        endif
        if(kolor.ne.1)then
            call newpen(kolor)
        endif
        call line(xx,yy,n,1,0,0)
        if(kolor.ne.1)then
            call newpen(1)
        endif
        if(shwabs)then
            if(kk.eq.1)then
            call number(x0+0.8*xlen-0.5,y0+ysn,0.10,ymxx,0.0,1003)
            else if(kk.eq.2)then
                call newpen(1)
            endif
        else
            call number(x0+0.8*xlen-0.5,y0+ysn,0.10,ymxx,0.0,1003)
        endif
            
        ls = lgstr(sym)
        if(kk.eq.1)call symbol(x0+0.8*xlen,y0+0.7*ylen,0.10,sym,0.0,ls)
        if(numgrn.eq.1)then
            call symbol(0.5,0.30,0.14,'TMIN =',0.0,+6)
            call number(1.4,0.30,0.14,xmin,0.0,+3)
            call symbol(2.5,0.30,0.14,'TMAX =',0.0,+6)
            call number(3.4,0.30,0.14,xmax,0.0,+3)
            call symbol(4.5,0.30,0.14,'DIST =',0.0,+6)
            call number(5.4,0.30,0.14,dist,0.0,+3)
            call symbol(6.5,0.30,0.14,'DEPTHR=',0.0,+7)
            call number(7.5,0.30,0.14,depthr,0.0,+3)
            call symbol(0.5,0.05,0.14,'DT =',0.0,+4)
            call number(1.1,0.05,0.14,dt,0.0,+4)
            call symbol(6.5,0.05,0.14,'DEPTHS=',0.0,+7)
            call number(7.5,0.05,0.14,depths,0.0,+3)
        endif
        return
        end

        subroutine gcmdln(names,jjplt,shwabs,kolor)
c-----
c       names(2)    C*(*)   - file names to be plotted
c       jjplt       I*4 - number of files entered on command line
c       shwabs      L   - .true. use f1 file to define scaling of both
c                       traces
c                       (default .false.)
c       kolor       I*4 - color of trace, default = 1
c-----
        integer*4 jjplt
        character names(2)*(*)
        character*40 name
        logical shwabs

        integer mnmarg
        jjplt = 0
        shwabs = .false.
        kolor = 1
        nmarg = mnmarg()
        i = 0
   11   i = i + 1
        if(i.gt.nmarg)goto 13
            call mgtarg(i,name)
            if(name(1:3).eq.'-f1')then
                i = i + 1
                jjplt = jjplt + 1
                call mgtarg(i,names(jjplt))
            else if(name(1:3).eq.'-f2')then
                i = i + 1
                jjplt = jjplt + 1
                call mgtarg(i,names(jjplt))
            else if(name(1:2).eq.'-R')then
                shwabs = .true.
            else if(name(1:2).eq.'-K')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,i10)')kolor 
                if(kolor.lt.0)kolor = 1
            else if(name(1:2).eq.'-?')then
                call usage()
            else if(name(1:2).eq.'-h')then
                call usage()
            else
                call usage()
            endif
        goto 11
   13   continue
        return
        end

        subroutine usage()
        integer*4 LER, LIN, LOT
        parameter(LER=0, LIN=5, LOT=6)
        write(LER,*)'USAGE: fplotg96 [-f1 fname1] [-f2 fname2]',
     1      '[-R] [-K kolor] [-?] [-h]'
        write(LER,*)
     1      ' -f1 fname1 (default none) file96 file name'
        write(LER,*)
     1      ' -f2 fname2 (default none) file96 file name'
        write(LER,*)
     1      ' -R        (default false) scale by fname1'
        write(LER,*)
     1      ' -K kolor    (default 1)  pen color'
        write(LER,*)
     1      ' -?  This help message'
        write(LER,*)
     1      ' -h  This help message'

        stop
        end

