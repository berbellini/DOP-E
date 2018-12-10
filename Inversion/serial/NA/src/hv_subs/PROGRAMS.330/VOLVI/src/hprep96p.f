        program hprep96p
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME V                                                       c
c                                                                     c
c      PROGRAM: HPREP96P                                              c
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
c       27 APR 2006 - minor change to output format for hspec96.dat 
c           replace e11.4 by e14.7
c-----
c       program to prepare input for hspec96(V)
c-----
c-----
c       LIN I*4 - logical unit for standard input
c       LOT I*4 - logical unit for standard output
c-----
        integer LIN, LOT, LER
        parameter (LIN=5,LOT=6,LER=0)
c-----
        integer NL
        parameter (NL=200)
c-----
c       mname   C*80    - name of the model file
c       dfile   C*80    - name of the distance file
c       ext L   - logical variable to see if the model file 
c       ierr    I*4 - 0 model file read in correctly
c-----
c       mmax    I*4 - number of layers in the model, last layer is
c                    halfspace
c-----
        character mname*80
        character dfile*80
        character title*80
        character hsfile*80, hrfile*80
        integer*4 jbdry, ieqex
        logical ext
        integer*4 mmax, iunit, iiso, iflsph

c-----
c       command line information
c-----
        integer dstcor
c-----
c       internal variables
c-----
        parameter(NSOURCE=100,NRECEIVER=100)
        common/source/depths(NSOURCE),lmaxs(NSOURCE),mdpths
        common/receiv/depthr(NRECEIVER),lmaxr(NRECEIVER),mdpthr
        real*4 depths, depthr
        integer lmaxs, lmaxr, mdpths, mdpthr

        
        parameter(NDST=500)
        real*4 r(NDST), tshift(NDST), vred(NDST)
c-----
c       parse the command line
c-----
        call gcmdln(dfile,mname,hsfile,hrfile,jbdry,ieqex,
     1      pmin,pmax,dp,pcntrl,xleng,xfac,hs,hr,ndec,dstcor,
     2      alphat)
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
c       open the distance file and attempt
c       to arrive at some common DT and NPTS
c-----  
        inquire(file=dfile,exist=ext)
        if(.not. ext)then
            ldf = lgstr(dfile)
            write(LER,*)'Distance file ', dfile(1:ldf),
     1          ' does not exist'
            call usage('Distance file does not exist')
        endif
        open(1,file=dfile,form='formatted',access='sequential',
     1      status='unknown')
        rewind 1
c-----
c       get minimum DT and maximum NPTS
c-----
            dtmin = 1.0E+30
            n = 0
            ndist = 0
            dstmax = 0.0
            dstmin = 1.0e+30
            call npow2(ndec)
 1000       continue
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
 1001   continue
        close(1)
        dt = dtmin
        delt = dt
        n1 = 1
        n2 = n/2 + 1
c-----
c       open the source depth file
c-----  
        if(hsfile .eq. ' ')then
            nsrc = 1
            depths(nsrc) = hs
        else
            inquire(file=hsfile,exist=ext)
            if(.not. ext)then
                ldf = lgstr(hsfile)
                write(LER,*)'Source depth file ', 
     1              hsfile(1:ldf),
     1          ' does not exist'
                call usage(' ') 
            endif
            open(1,file=hsfile,form='formatted',access='sequential',
     1          status='unknown')
            rewind 1
            nsrc = 0
 1100       continue
                read(1,*,end=1101,err=1101)hs
                nsrc = nsrc + 1
                depths(nsrc) = hs
            go to 1100
 1101       continue
            close (1)
        endif
c-----
c       open the receiver depth file
c-----  
        if(hrfile .eq. ' ')then
            nrec = 1
            depthr(nrec) = hr
        else
            inquire(file=hrfile,exist=ext)
            if(.not. ext)then
                ldf = lgstr(hrfile)
                write(LER,*)'Receiver depth file ', 
     1              hrfile(1:ldf),
     1              ' does not exist'
                call usage(' ') 
            endif
            open(1,file=hrfile,form='formatted',access='sequential',
     1          status='unknown')
            rewind 1
            nrec = 0
 1200       continue
                read(1,*,end=1201,err=1201)hr
                nrec = nrec + 1
                depthr(nrec) = hr
            go to 1200
 1201       continue
            close (1)
        endif
c-----
c       if the wavenumber integration sampling is not specified, 
c           attempt to
c       determine this automatically
c-----
c-----
c       open the output file hspec96p.dat
c-----
        open(2,file='hspec96p.dat',form='formatted',
     1      access='sequential',status='unknown')
        rewind 2
        write(2,2)dstcor
        write(2,1)alphat,delt
        write(2,2)n,n1,n2
        write(2,3)ieqex
        write(2,3)jbdry
        lmnm = lgstr(mname)
        write(2,4)mname(1:lmnm)
        write(2,1)xleng, xfac
        write(2,3)nsrc
        do 2008 i=1, nsrc
            write(2,6)depths(i)
 2008   continue
        write(2,3)nrec
        do 2009 i=1, nrec
            write(2,6)depthr(i)
 2009   continue
        write(2,3)ndist
        do 2010 i=1,ndist
            write(2,6)r(i),tshift(i),vred(i)
 2010   continue
        write(2,7)pmin,pmax,dp,pcntrl
    1   format(2e12.4)
    2   format(3i5)
    3   format(i5)
    4   format(a)
    6   format(5x,e14.7,5x,e14.7,5x,e14.7)
    7   format(e11.4,5x,e11.4,5x,e11.4,5x,e11.4)
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
        

        subroutine gcmdln(dfile,mname,hsfile,hrfile,jbdry,ieqex,
     1      pmin,pmax,dp,pcntrl,xleng,xfac,hs,hr,ndec,dstcor,
     2      alphat)
c-----
c       parse the command line arguments
c-----
c       dfile   C*80    - name of distacne file
c       mname   C*80    - name of model file
c       hsfile  C*80    - name of source depth file
c       jbdry   I*4 - 10*surface +halfspace
c           surface = 0 - elastic   halfspace = 0 - elastic
c                     1 - free                  1 - free
c                     2 - rigid                 2 - rigid
c       ieqex   I*4 = 0 EQ+EX
c                 1 PF+EX
c                 2 ALL
c                 3 EX
c                 4 EQ
c                 5 PF
c                 6 SH
c       pmin    R*4 lower bound of ray parameter
c       pmax    R*4 upper bound of ray parameter
c       dp  R*4 ray parameter increment
c       pcntrl  R*4  <= 0 use modified which does pi/2 
c           on some components
c               >0 do actual P-Tau
c       xleng   R*4 wavenumber integration DELTA
c                   DTLTA k = 6.2831853/ xleng
c       xfac    R*4 upper bound in k space at high frequencies
c       hs  R*4 source depth (single one specified
c       hr  R*4 receiver depth (single one specified
c       ndec    I*4 time domain decimation
c       dstcor  R*4 - 0 first time sample is tshift + r/vred
c                 1 first time sample is tshift + z/vred where
c                   z = abs (source depth - receiver depth)
c                 2 first time sample is tshift + sqrt(z*z + r*r) /vred
c       alphat  R*4 - control for FFT wrap around
c-----

c-----
        character mname*80
        character dfile*80
        character hrfile*80, hsfile*80
        integer*4 jbdry

        integer dstcor

        character name*40
        integer mnmarg

        dstcor = 0
        mname = ' '
        hrfile = ' '
        hsfile = ' '
        jbdry = 10
        nmarg = mnmarg()
        jtop = -1
        jbot = -1
        ieqex = 2
        pmin = -1.0
        pmax   = -1.0
        dp   = -1.0
        pcntrl = -1.0
        xfac = 4.0
        xleng = -1.0
        hs = 0.0
        hr = 0.0
        ndec = 1
        alphat = 2.5
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
                else if(name(1:3).eq.'-TF')then
                    jtop = 1
                else if(name(1:3).eq.'-TR' .and.
     1              name(1:5).ne.'-TRUE')then
                    jtop = 2
                else if(name(1:3).eq.'-TH')then
                    jtop = 0
                else if(name(1:3).eq.'-BF')then
                    jbot = 1
                else if(name(1:3).eq.'-BR')then
                    jbot = 2
                else if(name(1:3).eq.'-BH')then
                    jbot = 0
                else if(name(1:4).eq.'-ALL')then
                    ieqex = 2
                else if(name(1:5).eq.'-EQEX')then
                    ieqex = 0
                else if(name(1:4).eq.'-EXF')then
                    ieqex = 1
                else if(name(1:4).eq.'-ALP')then
                    i = i + 1
                    call mgtarg(i,name)
                    read(name,'(bn,f20.0)')alphat
                else if(name(1:5).eq.'-PMIN')then
                    i = i + 1
                    call mgtarg(i,name)
                    read(name,'(bn,f20.0)')pmin
                else if(name(1:5).eq.'-PMAX')then
                    i = i + 1
                    call mgtarg(i,name)
                    read(name,'(bn,f20.0)')pmax
                else if(name(1:3).eq.'-DP')then
                    i = i + 1
                    call mgtarg(i,name)
                    read(name,'(bn,f20.0)')dp
                else if(name(1:5).eq.'-TRUE')then
                    pcntrl = 1
                else if(name(1:5).eq.'-NDEC')then
                    i = i + 1
                    call mgtarg(i,name)
                    read(name,'(bn,i10)')ndec
                else if(name(2:2) .eq. '?')then
                    call usage(' ')
                endif
        go to 1000
 2000   continue
        if(jtop .ge. 0 .and. jbot .ge. 0)then
            jbdry = 10*jtop + jbot
        else
            jbdry = 10
        endif
        if(mname .eq. ' ')call usage(' ')
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
        write(LER,*)'Usage: hprep96p -M model  ',
     1      ' -d dfile -FHR recdep -FHS srcdep -HS hs -HR hr ' ,
     2      ' [-TF -TR -TH ] [-BF -BR -BH]',
     3      ' [-ALL -EQEX -EXF ]',
     4      ' [ -PMIN pmin -PMAX pmax -DP dp -TRUE ]',
     5      ' [-NDEC ndec] [-ALPHA alpha ] ',
     2       '[-?] [-h]'
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
     1  '-FHS srcdep (overrides -HS )  Name of source depth  file'
        write(LER,*)
     1  '-FHR recdep (overrides -HR )  Name of receiver depth  file'
        write(LER,*)
     1  '-HS hs      (default 0.0 )  Source depth  '
        write(LER,*)
     1  '-HR hr      (default 0.0 )  Receiver depth  '
        write(LER,*)
     1  '-TF        (default true )   top surface is free'
        write(LER,*)
     1  '-TR        (default false)   top surface is rigid'
        write(LER,*)
     1  '-TH        (default false)   top surface is halfspace'
        write(LER,*)
     1  '-BF        (default false)   bottom surface is free'
        write(LER,*)
     1  '-BR        (default false)   bottom surface is rigid'
        write(LER,*)
     1  '-BH        (default true )   bottom surface is halfspace'
        write(LER,*)
     1  '-ALL       (default true )   Compute all Green s functions'
        write(LER,*)
     1  '-EQEX      (default false)   Compute earthquake/explosion',
     2       'Green s functions'
        write(LER,*)
     1  '-EXF       (default false)   Compute explosion/point force',
     2       'Green s functions'
        write(LER,*)
     1  '-PMIN pmin -PMAX pmax -DP dp  (default none)',
     2      ' ray parameter sample space in sec/km'
        write(LER,*)
     1  '-TRUE      (default false) use modified p-tau response'
        write(LER,*)
     1  '-NDEC ndec (default 1) decimate the time series'
        write(LER,*)
     1  '-ALPHA alpha (default 2.5) time domain damping control'
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


