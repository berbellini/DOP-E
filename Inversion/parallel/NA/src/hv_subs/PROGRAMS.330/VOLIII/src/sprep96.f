        program sprep96
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME V                                                       c
c                                                                     c
c      PROGRAM: SPREP96                                               c
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
c       Changes
c       17 OCT 2002 - Added description of dfile format to usage routine
c       09 SEP 2012 - correct last line of bsort from n=n-m to nn=n-m
c             Thanks to ruoshan at ustc
c-----
c       program to prepare input for sdisp96(III)
c-----
c-----
c       LIN I*4 - logical unit for standard input
c       LOT I*4 - logical unit for standard output
c-----
        integer LIN, LOT, LER
        parameter (LIN=5,LOT=6,LER=0)
c-----
        integer NL
        parameter (NL=1000)
c-----
c       mname   C*80    - name of the model file
c       dfile   C*80    - name of the distance file
c       ext L   - logical variable to see if the model file 
c       ierr    I*4 - 0 model file read in correctly
c-----
c       mmax    I*4 - number of layers in the model, last layer is
c                    halfspace
c-----
        character title*80
        logical ext
        integer*4 mmax, iunit, iiso, iflsph

c-----
c       command line information
c-----
        character mname*80
        character dfile*80
        logical dolove, dorayl
        character hsfile*80, hrfile*80
c-----
c       internal variables
c-----
        parameter(NSOURCE=100,NRECEIVER=100)
        common/source/depths(NSOURCE),lmaxs(NSOURCE),mdpths
        common/receiv/depthr(NRECEIVER),lmaxr(NRECEIVER),mdpthr
        real*4 depths, depthr
        integer lmaxs, lmaxr, mdpths, mdpthr

        parameter(NPERIOD=2049)
        real fval(NPERIOD)
        integer nfval

        
        parameter(NDST=100)
        real*4 r(NDST), tshift(NDST), vred(NDST)
c-----
c       machine dependent initialization
c-----
        call mchdep()
c-----
c       parse the command line
c-----
        call gcmdln(dfile,mname,hsfile,hrfile,
     1      hs,hr,ndec,
     2      n,dt,dolove,dorayl,nmodes,faclov,facray,
     3      fval,nfval)
c-----
c       get the earth model
c-----
        inquire(file=mname,exist=ext)
        if(.not. ext)call usage('Model file does not exist')
        l = lgstr(mname)

        write(LOT,*)'Model name: ',mname(1:l)
        call getmod(1,mname,mmax,title,iunit,iiso,iflsph,
     1      idimen,icnvel,ierr,.true.)
c-----
c       use n1 and n2 to flag is user frequencies have been input
c-----
        if(nfval.le.0)then
            n1 = 1
            n2 = n/2 + 1
        else
            n1 = -1
            n2 = -2
            n = nfval
        endif
c-----
c       IF DT and NPTS are NOT DEFINED
c       OPEN the distance file
c-----
        if(dt.le.0.0 .and. n.le.0)then
c-----
c       open the distance file and attempt
c       to arrive at some common DT and NPTS
c-----  
            inquire(file=dfile,exist=ext)
            if(.not. ext)then
                ldf = lgstr(dfile)
                write(LER,*)'Distance file ', dfile(1:ldf),
     1              ' does not exist'
                call usage('Distance file does not exist')
            endif
            open(1,file=dfile,form='formatted',access='sequential',
     1          status='unknown')
            rewind 1
c-----
c           get minimum DT and maximum NPTS
c-----
                dtmin = 1.0E+30
                n = 0
                ndist = 0
                dstmax = 0.0
                dstmin = 1.0e+30
                call npow2(ndec)
 1000           continue
                    read(1,*,end=1001,err=1001)rr,dt,nn,tt,vr
                    ndist = ndist + 1
                    r(ndist) = rr
                    tshift(ndist) = tt
                    vred(ndist) = vr
                    if(ndec .gt. 1)then
                        nn = nn / ndec
                        dt = dt * ndec
                    endif
                    call npow2(nn)
                    if(nn.gt.n)n = nn
                    if(dt.lt.dtmin .and. dt.gt.0.0)dtmin = dt
                    if(rr.gt.dstmax)dstmax = rr
                    if(rr.lt.dstmin)dstmin = rr
                go to 1000
 1001       continue
            close(1)
            dt = dtmin
            delt = dt
        endif
            
c-----
c       output the control information
c-----
        open(2,file='sdisp96.dat',form='formatted',
     1      access='sequential',status='unknown')
        rewind 2
        write(2,1)dt
        write(2,2)n,n1,n2
        lmnm = lgstr(mname)
        write(2,4)mname(1:lmnm)
        write(2,8)dolove, dorayl
        write(2,1)hs,hr
        write(2,3)nmodes
        write(2,1)faclov,facray
        if(nfval.gt.0)then
            do 1010 i=1,nfval
                write(2,5)fval(i)
 1010       continue
        endif

    1   format(2e12.4)
    2   format(3i5)
    3   format(i5)
    4   format(a)
    5   format(e16.8)
    8   format(2l10)
        close (2)
        end


        subroutine iyesno(iyes)
        logical iyes
c-----
c       LIN I*4 - logical unit for standard input
c       LOT I*4 - logical unit for standard output
c-----
        integer LIN, LOT, LER
        parameter (LIN=5,LOT=6,LER=0)
        character ans*2

 1000   continue
            read(LIN,'(a)')ans
            if(ans(1:1).eq.'y' .or. ans(1:1).eq.'Y')then
                iyes = .true.
                return
            else if(ans(1:1).eq.'n' .or. ans(1:1).eq.'N')then
                iyes = .false.
                return
            endif
            write(LOT,*)'Enter (yn) only'
        go to 1000
        end
        
        subroutine gcmdln(dfile,mname,hsfile,hrfile,
     1      hs,hr,ndec,
     2      npts,dt,dolove, dorayl,nmodes,faclov,facray,
     3      fval,nfval)
c-----
c       parse the command line arguments
c-----
c       dfile   C*80    - name of distance file
c       mname   C*80    - name of model file
c       hsfile  C*80    - name of source depth file
c       hrfile  C*80    - name of receiver depth file
c       hs  R*4 source depth (single one specified
c       hr  R*4 receiver depth (single one specified
c       ndec    I*4 time domain decimation
c       npts    I*4 Number of points in time series
c       dt  R*4 Sampling interval
c       dolove  L   Get Love wave dispersion
c       dorayl  L   Get Rayleigh wave dispersion
c       nmodes  I*4 Maximum number of modes
c       faclov  R*4 factor for controlling root search
c       facray  R*4 factor for controlling root search
c       fval    R*4 array of specific periods for evaluation
c       nfval   I*4 number of user specified periods or frequencies
c-----

c-----
        character mname*80, dfile*80
        character hrfile*80, hsfile*80
        logical dolove, dorayl
        integer nmodes
        real*4 facray, faclov

        parameter(NPERIOD=2049)
        real fval(NPERIOD)
        integer nfval

        real pmin, pmax
        real fmin, fmax

        character name*40
        integer mnmarg

        mname = ' '
        dfile = ' '
        hrfile = ' '
        hsfile = ' '
        hs = 0.0
        hr = 0.0
        ndec = 1
        npts = -1
        dt = -1.0
        dolove = .false.
        dorayl = .false.
        nmodes = 1
        faclov = 5.0
        facray = 5.0
        pmin = -1.0
        pmax = -1.0
        fmin = -1.0
        fmax = -1.0
        nmarg = mnmarg()
        nfval = -1
        i = 0
 1000   continue
            i = i + 1
            if(i.gt.nmarg)go to 2000
                call mgtarg(i,name)
                if(name(1:2).eq.'-M')then
                    i = i + 1
                    call mgtarg(i,mname)
                else if(name(1:2).eq.'-d')then
                    i = i + 1
                    call mgtarg(i,dfile)
                else if(name(1:4).eq.'-FHR')then
                    i = i + 1
                    call mgtarg(i,hrfile)
                else if(name(1:4).eq.'-FHS')then
                    i = i + 1
                    call mgtarg(i,hsfile)
                else if(name(1:3).eq.'-HR')then
                    i = i + 1
                    call mgtarg(i,name)
                    read(name,'(bn,f20.0)')hr
                else if(name(1:3).eq.'-HS')then
                    i = i + 1
                    call mgtarg(i,name)
                    read(name,'(bn,f20.0)')hs
                else if(name(1:5).eq.'-NDEC')then
                    i = i + 1
                    call mgtarg(i,name)
                    read(name,'(bn,i10)')ndec
                else if(name(1:3).eq.'-DT')then
                    i = i + 1
                    call mgtarg(i,name)
                    read(name,'(bn,f20.0)')dt
                else if(name(1:5).eq.'-NPTS')then
                    i = i + 1
                    call mgtarg(i,name)
                    read(name,'(bn,i10)')npts
                else if(name(1:5).eq.'-NMOD')then
                    i = i + 1
                    call mgtarg(i,name)
                    read(name,'(bn,i10)')nmodes
                else if(name(1:2).eq.'-L')then
                    dolove = .true.
                else if(name(1:2).eq.'-R')then
                    dorayl = .true.
                else if(name(1:5).eq.'-FACR')then
                    i = i + 1
                    call mgtarg(i,name)
                    read(name,'(bn,f20.0)')facray
                else if(name(1:5).eq.'-FACL')then
                    i = i + 1
                    call mgtarg(i,name)
                    read(name,'(bn,f20.0)')faclov
                else if(name(1:5).eq.'-FREQ')then
                    i = i + 1
                    call mgtarg(i,name)
                    read(name,'(bn,f20.0)')fval(1)
                    nfval = 1.0
                else if(name(1:4).eq.'-PER')then
                    i = i + 1
                    call mgtarg(i,name)
                    read(name,'(bn,f20.0)')fval(1)
                    fval(1) = 1.0/fval(1)
                    nfval = 1.0
                else if(name(1:5).eq.'-FARR')then
                    i = i + 1
                    call mgtarg(i,name)
                    call getfv(name,nfval,fval,.true.)
                else if(name(1:5).eq.'-PARR')then
                    i = i + 1
                    call mgtarg(i,name)
                    call getfv(name,nfval,fval,.false.)
                else if(name(1:4).eq.'-PER')then
                    i = i + 1
                    call mgtarg(i,name)
                    read(name,'(bn,f20.0)')fval(1)
                    fval(1) = 1.0/fval(1)
                    nfval = 1.0
                else if(name(1:5).eq.'-FMIN')then
                    i = i + 1
                    call mgtarg(i,name)
                    read(name,'(bn,f20.0)')fmin
                else if(name(1:5).eq.'-FMAX')then
                    i = i + 1
                    call mgtarg(i,name)
                    read(name,'(bn,f20.0)')fmax
                else if(name(1:5).eq.'-PMIN')then
                    i = i + 1
                    call mgtarg(i,name)
                    read(name,'(bn,f20.0)')pmin
                else if(name(1:5).eq.'-PMAX')then
                    i = i + 1
                    call mgtarg(i,name)
                    read(name,'(bn,f20.0)')pmax
                else if(name(1:2) .eq. '-?')then
                    call usage(' ')
                else if(name(1:2) .eq. '-h')then
                    call usage(' ')
                endif
        go to 1000
 2000   continue
        if(mname .eq. ' ')call usage(' ')
        if(faclov.lt.1)then faclov = 5.0
        if(facray.lt.1)then facray = 5.0
c-----
c       automatically define frequencies if pmin, pmax > 0.0
c-----
        if(pmin.gt.0.0 .and. pmax.gt.0.0)then
                call getfva(nfval,fval,pmin,pmax,.false.)
        endif
        if(fmin.gt.0.0 .and. fmax.gt.0.0)then
                call getfva(nfval,fval,fmin,fmax,.true.)
        endif
        return
        end

        subroutine getfva(nfval,fval,vmin,vmax,isfreq)
c-----
c       define array of periods or frequencies that correspond to
c       sacmft96 values
c-----
c       nfval   I*4 number of frequency values
c       fval    R*4 array of frequency values
c       vmin    R*4 
c       vmax    R*4 bounds on period or frequency
c               IT IS ASSUMED THAT vmin and vmax > 0 ENTERING
c               INTO THIS ROUTINE
c       isfreq  L   .true. array of frequencies, .false. 
c            array of periods
c-----
        logical isfreq
        parameter(NPERIOD=2049)
        real fval(NPERIOD)
        real vmin, vmax
        integer nfval

        real tmp
        real pmin, pmax
c-----
c       convert to period - we really want frequency 
c            but this way we can
c       reuse the getper routine from sdprad96
c-----
        if(isfreq)then
            pmin = 1.0/vmax
            pmax = 1.0/vmin
        else
            pmin = vmin
            pmax = vmax
        endif
c-----
c       order to make pmin the minimum value
c-----
        if(pmin.gt.pmax)then
            tmp = pmin
            pmin = pmax
            pmax = tmp
        endif
c-----
c       now get array of periods
c-----
        WRITE(0,*)'pmin,pmax:',pmin,pmax
        call getper(fval,NPERIOD,pmin,pmax,nfval,pmin,pmax)
c-----
c       convert period to frequency array, ordered from 
c            high frequency to
c       low frequency - this is trivial since the periods are ordered
c       from short to long period
c-----
        do 1000 i=1,nfval
            fval(i) = 1.0/fval(i)
 1000   continue
        
        return
        end

        subroutine getfv(name,nfval,fval,isfreq)
c-----
c       get array of periods or frequencies from a separate file
c-----
c       name    C*40    name of file
c       nfval   I*4 number of frequency values
c       fval    R*4 array of frequency values
c       isfreq  L   .true. array of frequencies, 
c            .false. array of periods
c-----
        character name*(*)
        logical isfreq
        parameter(NPERIOD=2049)
        real fval(NPERIOD)
        integer nfval
        logical ext

        inquire(file=name,exist=ext)
        if(.not. ext)call usage('Named period/freq file soes not exist')
        open(1,file=name,access='sequential',form='formatted',
     1      status='unknown')
        rewind 1
        nfval = 0
 1000   continue
            read(1,*,end=9999,err=9999)xx
            if(xx.le.0.0)go to 1000
            nfval = nfval + 1
            if(isfreq)then
                fval(nfval) = xx
            else
                fval(nfval) = 1.0/xx
            endif
            go to 1000
 9999   continue
        close (1)
c-----
c       sort frequencies in order of decreasing frequency
c-----
        if(nfval.gt.1)then
            call bsort(nfval,fval,-1)
        endif
        return
        end

c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
        subroutine bsort(nn,x,isign)
c       do bubble sort.
c       isign=+1  increase   =-1 decrease.
c
          real*4 x(nn)
c
        n = nn
        m=0
        do 50 i=2,n
          ii=i-m-1
          do 40 j=ii,1,-1
            x0=x(j+1)-x(j)
              if(isign.le.0) x0=-x0
            if(abs(x0).lt.1.e-7) go to 20
C              if(x0) 10,20,50
            if(x0 .lt. 0.0)then
                go to 10
            else if(x0 .eq. 0.0)then
                go to 20
            else
                go to 50
            endif
10        continue
            x0=x(j)
            x(j)=x(j+1)
            x(j+1)=x0
            go to 40
20        continue
            m=m+1
            do 30 k=j,n-m
              x(k)=x(k+1)
30        continue
            go to 50
40      continue
50    continue
        nn=n-m
        return
        end

        subroutine usage(ostr)
        integer LER
        parameter (LER=0)
        character ostr*(*)
        if(ostr .ne. ' ')then
            lostr = lgstr(ostr)
            write(LER,*)ostr(1:lostr)
        endif
        write(LER,*)'Usage: sprep96 -M model  ',
     1      ' -d dfile -HS hs -HR hr ' ,
     2      ' -DT dt -NPTS npts -NMOD nmodes ',
     3      ' -L -R',
     4      ' -FACR facray -FACL faclov',
     5      ' -FREQ freq -PER period -FARR freq_file',
     6      ' -PARR period_file -FMIN fmin -FMAX fmax',
     7      ' -PMIN pmin -PMAX pmax [-?] [-h]'
        write(LER,*)
     1  '-M model   (default none )  Earth model file'
        write(LER,*)
     1  '-d dfile   (default none )  Name of distance file'
        write(LER,*)
     1  '   dfile contains one of more lines with following entries'
        write(LER,*)
     1  '       DIST(km) DT(sec) NPTS T0(sec) VRED(km/s)',
     2  '           first time point is T0 + DIST/VRED',
     3  '           VRED=0 means infinite velocity though'
        write(LER,*)
     1  '-HS hs      (default 0.0 )  Source depth '
        write(LER,*)
     1  '-HR hr      (default 0.0 )  Receiver depth'
        write(LER,*)
     1  '-DT dt      (default 1.0 )  Sampling interval'
        write(LER,*)
     1  '-NPTS npts  (default 1   )  Number of points'
        write(LER,*)
     1  '-NMOD nmodes (default 1   )  Maximum number of modes'
        write(LER,*)
     1  '-L          (default false) Generate Love Waves'
        write(LER,*)
     1  '-R          (default false) Generate Rayleigh Waves'
        write(LER,*)
     1  '-FACL faclov (default 5)    Factor to insure finer search'
        write(LER,*)
     1  '-FACR facray (default 5)    Factor to insure finer search'
        write(LER,*)
     1  '-FREQ freq                  User specified single frequency'
        write(LER,*)
     1  '-PER  period                User specified single period'
        write(LER,*)
     1  '-FARR freq_file             File of user specified frequencies'
        write(LER,*)
     1  '-PARR period_file           File of user specified periods'
        write(LER,*)
     1  '-FMIN fmin -FMAX fmax       Auto generate within freq bounds'
        write(LER,*)
     1  '-PMIN pmin -PMAX pmax       Auto generate within per bounds'
        write(LER,*)
     1  ' NOTE: one of -DT -NPTS, -FARR, -PARR, -PER or -FREQ, ',
     2  ' or -FMIN fmin & -FMAX fmax, -PMIN pmin & -PMAX pmax required'

        write(LER,*)
     1  '-?        (default none )  this help message '
        write(LER,*)
     1  '-h        (default none )  this help message '
        stop
        end
        

        subroutine npow2(npts)
c-----
c       Given npts, determine the N=2**m such that N >= npts
c       return the new ntps
c-----
        integer*4 nsamp, npts
        nsamp = npts
        npts = 1
 1000   continue
        if(npts.ge.nsamp)return
        npts = 2*npts
        go to 1000
        end

        subroutine getper(wn,NX,pmin,pmax,nper,permin,permax) 
c-----
c       automatically create periods from the [pmin,pmax] limits
c       wn()    R*4 array of periods
c       NX  I   dimension of array
c       pmin    R*4 minimum period
c       pmax    R*4 maximum period
c       nper    I*4 number of periods generated
c       permin  R*4 - minimum period to be used in trace
c       permax  R*4 - maximum period to be used in trace
c----
c       this routine is based on the one in sacmft96 
c            which generates the
c       dispersion values
c-----
        implicit none
        integer NX, nper
        real wn(NX)
        real pmin, pmax
        real permin, permax

        real ymxlog, ymmin, p, tenpow
        integer nocy, iy, ii, jj
        integer NP
        parameter (NP=40)
        real pfac(NP)
        data pfac/1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 
     1      2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9,
     2      3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8,
     3      5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5/

        nper = 0
c-----
c       determine starting power
c-----
        ymxlog = alog10(pmax)   
        ymmin  = alog10(pmin)   
        nocy = ymxlog - ymmin + 1 
        iy = ymmin
        if(ymmin .lt. 0)iy = iy - 1
        do 100 ii=iy,iy+nocy+2
            tenpow = 10.0**ii
            do 200 jj=1,NP
                p = pfac(jj)*tenpow
                if(p .ge. pmin .and. p .le. pmax
     1              .and. p.ge.permin .and. p.le.permax)then
                    nper = nper + 1
                    wn(nper) = p
                    if(nper.eq.NX)go to 1000
                endif
 200        continue
 100    continue
 1000   continue
        return
        end
