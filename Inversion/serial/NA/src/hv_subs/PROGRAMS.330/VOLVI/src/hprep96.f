        program hprep96
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME VI                                                      c
c                                                                     c
c      PROGRAM: HPREP96                                               c
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
c       02 MAY 2002 - extend to use both 
c           isotropic and transverse isotropic
c           models
c       17 OCT 2002 - Added description of dfile format to usage routine
c       16 FEB 2004 - in subroutine fstarr, 
c           change complex*16 pold, pcur, dp
c           to make dp real*8 - a.bitri@brgm.fr
c       20 JUL 2004 - corrected cast line 1210
c           p = dcmplx(i*dp, 0.0d+00) where dp is already complex
c       09 JAN 2005 - vlmm not defined in subroutine 
c           fstarr - replaced with vlmn
c           baker@usgs.gov
c       27 APR 2006 - minor change to output format for hspec96.dat 
c           replace e11.4 by e14.7
c       04 AUG 2006 - corrected error in first arrival pick that
c               falsely gave the refraction time instead of the
c               direct time because the refraction arrival was
c               unphysical
c       14 MAR 2008 - add a multiplier flag to increase the L used
c               in wavenumber integration, e.g., make dk smaller.
c               -ML lmult  The default is 1.0. Mult be >= 1
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
        integer*4 iunit, iiso, iflsph

        common/depref/refdep
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

        integer mmax
        common/modlly/mmax

        common/earth/radius
        real radius

        parameter(NDST=500)
        real*4 r(NDST), tshift(NDST), vred(NDST)

        real SSA, SSC, SSF, SSL,SSN, SSR
        real RRA, RRC, RRF, RRL,RRN, RRR
c-----
c       machine dependent initialization
c-----
        call mchdep()
c-----
c       initialize
c-----
        radius = 6371.
c-----
c       parse the command line
c-----
        call gcmdln(dfile,mname,hsfile,hrfile,jbdry,ieqex,
     1      cmin,c1,c2,cmax,xleng,xfac,hs,hr,ndec,dstcor,
     2      alphat,flmult)
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
                r(ndist) = abs(rr)
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
        if(hsfile(1:1) .eq. ' ')then
            nsrc = 1
            depths(nsrc) = hs
            call insert(hs+refdep)
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
                call insert(hs+refdep)
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
            call insert(hr+refdep)
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
                call insert(hr+refdep)
            go to 1200
 1201       continue
            close (1)
        endif
c-----
c       if the wavenumber integration sampling is not specified, 
c           attempt to
c       determine this automatically
c-----
        if(xleng.le.0.0)then
c-----
c           first cycle through all source, receiver depths and
c           receiver distances to estimate the maximum time to be
c           synthetized
c----- 
            tmax = 0.0
            zmax = 0.0
            do 1501 i=1,nsrc    
                ds = depths(i)
                do 1502 j=1,nrec
                    dr = depthr(j)
                    dz = abs(ds - dr)
                    if(dz.gt.zmax)zmax = dz
                    do 1503 k=1,ndist
                        rr = r(k)
                        tt = tshift(k)
                        vr = vred(k)
                        call gett0(t0,rr,ds,dr,
     1                      tt,vr,dstcor)
                        tm = t0 + (n-1)*delt
                        if(tm.gt.tmax)tmax = tm
 1503               continue
 1502           continue
 1501       continue
            call estlen(xleng,dstmax,dstmin,tmax,
     2          nsrc,nrec,dstcor,zmax,mname)
        endif
c-----
c       open the output file hspec96.dat
c-----
        open(2,file='hspec96.dat',form='formatted',
     1      access='sequential',status='unknown')
        rewind 2
        write(2,2)dstcor
        write(2,1)alphat,delt
        write(2,2)n,n1,n2
        write(2,3)ieqex
        write(2,3)jbdry
        lmnm = lgstr(mname)
        write(2,4)mname(1:lmnm)
        write(2,1)xleng*flmult, xfac
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
c-----
c       read in phase velocity limits for wavenumber filtering
c       Wavenumber filtering will consist of following
c
c       |*ZERO*|-COSINE TAPER-|*ALL PASS*|-COSINE TAPER-|*ZERO
c       |      |              |          |              |
c            omega          omega      omega          omega
c k =   0    -----          -----      -----          -----   infinity
c            cmax            c1         c2             cmin
c-----
c       If c2 or cmin <= 0, then upper wavenumber limit is infinite
c       If c1 or cmax <= 0, then lower wavenumber limit is zero
c-----
c       before writing out results
c       verify the ordering
c-----
        call chkwvl(cmax,c1,c2,cmin)
        write(2,7)cmax,c1,c2,cmin
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
     1      cmin,c1,c2,cmax,xleng,xfac,hs,hr,ndec,dstcor,
     2      alphat,flmult)
c-----
c       parse the command line arguments
c-----
c       dfile   C*80    - name of distance file
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
c       cmin    R*4 phase velocity filter
c       c1  R*4 phase velocity filter
c       c2  R*4 phase velocity filter
c       cmax    R*4 phase velocity filter
c       xleng   R*4 wavenumber integration DELTA
c                   DTLTA k = 6.2831853/ xleng
c       xfac    R*4 upper bound in k space at high frequencies
c       hs  R*4 source depth (single one specified
c       hr  R*4 receiver depth (single one specified
c       ndec    I*4 time domain decimation
c       dstcor  I*4 - 0 first time sample is tshift + r/vred
c                 1 first time sample is tshift + z/vred where
c                   z = abs (source depth - receiver depth)
c                 2 first time sample is tshift + sqrt(z*z + r*r) /vred
c       alphat  R   - alphas value 2.5 is the default
c       flmult  R   - multiplier of the computed length for wavenumber
c                     integration - default is 1.0
c-----

c-----
        character mname*80, dfile*80
        character hrfile*80, hsfile*80
        integer*4 jbdry

        real cmin,c1,c2,cmax,xleng,xfac,hs,hr
        real alphat,flmult

        integer ndec, dstcor

        character name*40
        integer mnmarg

        dstcor = 0
        mname = ' '
        dfile = ' '
        hrfile = ' '
        hsfile = ' '
        jbdry = 10
        nmarg = mnmarg()
        jtop = -1
        jbot = -1
        ieqex = 2
        cmin = -1.0
        c1   = -1.0
        c2   = -1.0
        cmax = -1.0
        xfac = 4.0
        xleng = -1.0
        hs = 0.0
        hr = 0.0
        ndec = 1
        alphat = 2.5
        flmult = 1.0
        i = 0
 1000   continue
            i = i + 1
            if(i.gt.nmarg)go to 2000
                call mgtarg(i,name)
                if(name(1:2).eq.'-M' .and.
     1            name(1:3).ne.'-ML')then
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
                else if(name(1:3).eq.'-TR')then
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
                else if(name(1:3).eq.'-EQ')then
                    ieqex = 0
                else if(name(1:3).eq.'-EX')then
                    ieqex = 1
                else if(name(1:4).eq.'-ALP')then
                    i = i + 1
                    call mgtarg(i,name)
                    read(name,'(bn,f20.0)')alphat
                else if(name(1:5).eq.'-CMIN')then
                    i = i + 1
                    call mgtarg(i,name)
                    read(name,'(bn,f20.0)')cmin
                else if(name(1:3).eq.'-C1')then
                    i = i + 1
                    call mgtarg(i,name)
                    read(name,'(bn,f20.0)')c1
                else if(name(1:3).eq.'-C2')then
                    i = i + 1
                    call mgtarg(i,name)
                    read(name,'(bn,f20.0)')c2
                else if(name(1:5).eq.'-CMAX')then
                    i = i + 1
                    call mgtarg(i,name)
                    read(name,'(bn,f20.0)')cmax
                else if(name(1:3).eq.'-XL')then
                    i = i + 1
                    call mgtarg(i,name)
                    read(name,'(bn,f20.0)')xleng
                else if(name(1:3).eq.'-ML')then
                    i = i + 1
                    call mgtarg(i,name)
                    read(name,'(bn,f20.0)')flmult
                    if(flmult.lt.1.0)flmult = 1.0
                else if(name(1:3).eq.'-XF')then
                    i = i + 1
                    call mgtarg(i,name)
                    read(name,'(bn,f20.0)')xfac
                else if(name(1:5).eq.'-NDEC')then
                    i = i + 1
                    call mgtarg(i,name)
                    read(name,'(bn,i10)')ndec
                else if(name(1:2).eq.'-R')then
                    dstcor = 2
                else if(name(1:2).eq.'-Z')then
                    dstcor = 1
                else if(name(2:2) .eq. '?')then
                    call usage(' ')
                endif
        go to 1000
 2000   continue
        if(jtop .ge. 0 .and. jbot .ge. 0)then
            jbdry = 10*jtop + jbot
        else
            jbdry= 10
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
        write(LER,*)'Usage: hprep96 -M model  ',
     1      ' -d dfile -FHR recdep -FHS srcdep -HS hs -HR hr ' ,
     2      ' [-TF -TR -TH ] [-BF -BR -BH]',
     3      ' [-ALL -EQEX -EXF ]',
     4      ' [ -CMAX cmax -C1 c1 -C2 c2 -CMIN cmin ]',
     5      ' [ -XL xleng -XF xfac -ML xleng_mul ] [-NDEC ndec] ',
     5      ' [-ALP alp ',
     2       '[-Z] [-R] [-?] [-h]'
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
     1  '-HS hs      (default 0.0 )  Source depth '
        write(LER,*)
     1  '-HR hr      (default 0.0 )  Receiver depth'
        write(LER,*)
     1  '-TF         (default true )   top surface is free'
        write(LER,*)
     1  '-TR         (default false)   top surface is rigid'
        write(LER,*)
     1  '-TH         (default false)   top surface is halfspace'
        write(LER,*)
     1  '-BF         (default false)   bottom surface is free'
        write(LER,*)
     1  '-BR         (default false)   bottom surface is rigid'
        write(LER,*)
     1  '-BH         (default true )   bottom surface is halfspace'
        write(LER,*)
     1  '-ALL        (default true )   Compute all Green s functions'
        write(LER,*)
     1  '-EQEX       (default false)   Compute earthquake/explosion',
     2       'Green s functions'
        write(LER,*)
     1  '-EXF       (default false)   Compute explosion/point force',
     2       'Green s functions'
        write(LER,*)
     1  '-CMAX cmax -C1 c1 -C2 c2 -CMIN cmin (default none)',
     2      ' phase velocity filter band'
        write(LER,*)
     1  '-XL xleng (default automatic determination ) ',
     2      'DELTA k = 6.2831853 / xleng' 
        write(LER,*)
     1  '-ML xleng_mul (default 1.0 ) ',
     2      'increase xleng over defined value ',
     3      ' reduce numerical noise, wrap around'
        write(LER,*)
     1  '-XF xfac  (default 4.0) upper bound in wvno parameter'
        write(LER,*)
     1  '-NDEC ndec (default 1) decimate the time series'
        write(LER,*)
     1  '-ALP alp   (default 2.5) time domain damping control'
        write(LER,*)
     1  '-Z t0=tshift+abs(source depth - receiver depth)/vred'
        write(LER,*)
     1  '-R  t0=tshift+sqrt(z*z + r*r)/vred'
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


        subroutine estlen(xleng,dstmax,dstmin,tmax,nsrc,nrec,
     1          dstcor,zmax,mname)
c-----
c       attempt to estimate the wavenumber sampling parameter xleng
c       from the maxmum distance and distance
c
c       The Bouchon criteria  to avoid periodicity problems in
c       wavenumber of integration fo the type
c
c       INT J n (kr) DELTA k , where DELTA k is 6.2831853 / L
c
c       is
c       (1) L > 2 dstmax
c       and
c       (2) first_arrival_time(L - r )  > (  tmax ) 
c
c       The periodicity in the wavenumber sampling maps into the
c       time domain. A record section will like like
c
c       |    \/   
c      T|sig /\
c       |   /  \ noise
c       |  /    \
c       | /      \
c       |/        \
c       |__________\_____________ r
c
c                  L=2r
c
c       The first criteria arises from sampling of the Bessel function
c       and really applies to the t=0 time value
c       The second tries to avoid the numerical noise arrival
c           propagating with negative moveout
c-----
c       xleng   R*4 - the value to be determined, Bouchon's L
c       dstmax  R*4 - maximum distance in record section
c       dstmin  R*4 - minimum distance in record section
c       tmax    R*4 - absolute travel time of largest time
c                   value computed
c       mmax    I*4 - number of layers in the model
c       nsrc    I*4 - number of source layers
c       nrec    I*4 - number of receiver layers
c       dstcor  I*4 - 0 first time sample is tshift + r/vred
c                 1 first time sample is tshift + z/vred where
c                   z = abs (source depth - receiver depth)
c                 2 first time sample is tshift + sqrt(z*z + r*r) /vred
c       zmax    R*4 - maximum separation between source and receiver
c       mname   Ch*80   - model file name
c
c-----
        integer dstcor


        integer NL
        parameter (NL=200)
        common/timodel/D(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      TRho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL), 
     3      frefp(NL), frefs(NL)
        real D, TA, TC, TF, TL, TN, TRho
        real qa, qb, etap, etas, frefp, frefs
        common/modlly/mmax
        integer mmax
        character title*80 

        integer l, lgstr
        integer iunit,iiso,iflsph,idimen,icnvel
        integer ierr

        parameter(NSOURCE=100,NRECEIVER=100)
        common/source/depths(NSOURCE),lmaxs(NSOURCE),mdpths
        common/receiv/depthr(NRECEIVER),lmaxr(NRECEIVER),mdpthr
        real*4 depths, depthr
        integer lmaxs, lmaxr, mdpths, mdpthr

        common/depref/refdep

        character mname*(*)

c-----
c       this is somewhat of an inverse problem
c
c       For refracted arrivals we have
c
c       t = t0 + x/v, so given t, t0 we can solve for x 
c           and use x = L - r
c           to estimate L
c-----
c       cycle through all source and receiver depths
c       depths  - source depth
c       lmaxs   - layer in which source lies
c       dephs   - height of source above bottom of layer
c       
c       depthr  - receiver depth
c       lmaxr   - layer in which receiver lies
c       dephr   - height of receiver above bottom of layer
c-----
        r = dstmax
        xleng = 2.0*r
c-----
c       we will iteratively attempt to determine xleng.
c       First we need an initial estimate. This will be
c       either 2 r or tmax*vfastest
c-----
        vmax = 0.0
        do 1234 i=1,mmax
            avel = sqrt(TA(i)/TRho(i))
            if(avel .gt. vmax)vmax = avel
            avel = sqrt(TC(i)/TRho(i))
            if(avel .gt. vmax)vmax = avel
            avel = sqrt(TF(i)/TRho(i))
            if(avel .gt. vmax)vmax = avel
 1234   continue
        xl = tmax * vmax
        if(xl .gt. xleng)xleng = xl
 1000   continue
            tfirst = 1.0e+38
            do 4000 is=1,nsrc
                do 4100 ir=1,nrec
        CALL FRSTAR(xleng-R,DEPTHS(is),DEPTHR(ir),MNAME,1,TP,
     1      SSA, SSC, SSF, SSL,SSN, SSR,
     1      RRA, RRC, RRF, RRL,RRN, RRR,
     1      rayp, geom, tstar, .false.)
                    if(TP.lt.tfirst)tfirst = TP
 4100       continue
 4000       continue
            if(tfirst .lt. tmax)then
                xleng = xleng + 0.25*xleng 
            else 
                go to 1001
            endif
        go to 1000
c-----
c       safety
c-----
 1001   continue
        if(xleng .lt. 2.0 * r)xleng = 2.0 * r
        return
        end

        subroutine insert(dph)
        implicit none
        real dph
        integer LER, LIN, LOT
        parameter (LER=0, LIN=5, LOT=6)
        integer NL
        parameter(NL=200)
        common/timodel/D(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      TRho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL), 
     3      frefp(NL), frefs(NL)
        real D, TA, TC, TF, TL, TN, TRho
        real qa, qb, etap, etas, frefp, frefs
        common/modlly/mmax
        integer mmax

        integer m
        real dep, dp, dphh, hsave
        integer ls
c-----
c       Insert a depth point into the model by splitting a layer
c       so that the point appears at the top boundary of the layer
c       dph = depth of interest
c-----
c       determine the layer in which the depth dph lies.
c       if necessary, adjust  layer thickness at the base
c-----
c       Here determine layer corresponding to specific depth dph
c       If the bottom layer is not thick enough, extend it
c
c       dep - depth to bottom of layer
c       dphh    - height of specific depth above bottom of the layer
c-----
        if(dph.le.0)then
            d(1) = d(1) - dph
            return
        else if(dph.ge.0)then
            dep = 0.0 
            dp = 0.0 
            dphh = -1.0
            do 100 m = 1,mmax 
                dp = dp + d(m) 
                dphh = dp - dph 
                if(m.eq.mmax)then
                    if(d(mmax).le.0.0 .or. dphh.lt.0.0)then
                        d(mmax) = (dph - dp)
                    endif
                endif
                dep = dep + d(m) 
                dphh = dep - dph 
                ls = m 
                if(dphh.ge.0.0) go to 101 
  100       continue 
  101       continue 
        endif
c-----
c       In the current model, the depth point is in the ls layer
c       with a distance dphh to the bottom of the layer
c
c       Do not create unnecessary layers, e.g., 
c           at surface and internally
c       However do put in a zero thickness layer 
c           at the base if necessary
c-----
        if(dphh .eq. 0.0 .and. ls.ne.mmax)then
            return
        else
c-----
c           adjust layering
c-----
             do 102 m = mmax,ls,-1
                d(m+1) = d(m)
                TA(m+1) = TA(m)
                TC(m+1) = TC(m)
                TF(m+1) = TF(m)
                TL(m+1) = TL(m)
                TN(m+1) = TN(m)
                TRho(m+1) = TRho(m)
                qa(m+1) = qa(m)
                qb(m+1) = qb(m)
                etap(m+1) = etap(m)
                etas(m+1) = etas(m)
                frefp(m+1) = frefp(m)
                frefs(m+1) = frefs(m)
  102       continue
            hsave = d(ls)
            d(ls) = hsave - dphh
            d(ls+1) = dphh
            ls = ls + 1
            mmax = mmax + 1
            if(d(mmax).lt.0.0)d(mmax)=0.0
        endif
        return
        end

        subroutine gett0(t0,r,depths,depthr,tshift,vred,dstcor)
        real*4 t0, depths, depthr, tshift, vred
        integer*4 dstcor
c-----
c       compute time of first sample of the time series
c-----
            if(dstcor.eq.0)then
                rr = r
            else if(dstcor.eq.1)then
                rr = abs(depthr - depths)
            else if(dstcor.eq.2)then
                rr = sqrt(r*r + (depthr-depths)*(depthr-depths))
            endif
            if(vred.eq.0.0)then
                t0 = tshift 
            else
                t0 = tshift + rr/vred
            endif
        return
        end

        subroutine chkwvl(cmax,c1,c2,cmin)
c-----
c       Wavenumber filtering will consist of following
c
c       |*ZERO*|-COSINE TAPER-|*ALL PASS*|-COSINE TAPER-|*ZERO
c       |      |              |          |              |
c            omega          omega      omega          omega
c k =   0    -----          -----      -----          -----   infinity
c            cmax            c1         c2             cmin
c-----
c       If c2 or cmin <= 0, then upper wavenumber limit is infinite
c       If c1 or cmax <= 0, then lower wavenumber limit is zero
c-----
        real tmp1(4), tmp2(4)
        integer key(4)
        if(cmin .lt. 0.0)then
            c2 = -1.0
        endif
        if(cmax .lt. 0.0)then
            c1 = -1.0
        endif
c-----
c       now do sort if possible
c-----
        if(c1.lt.0.0 .and. c2.lt.0.0 
     1      .and. cmin.lt.0.0 .and. cmax.lt.0.0)then
            return
        else
            if(cmin.lt.0.0)then
                xmx = amax1(c1,cmax)
                xmn = amin1(c1,cmax)
                cmax = xmx
                c1   = xmn
            else if(cmax.lt.0.0)then
                xmx = amax1(c2,cmin)
                xmn = amin1(c2,cmin)
                cmin = xmn
                c2   = xmx
            else
                tmp1(1) = cmax
                tmp1(2) = c1
                tmp1(3) = c2
                tmp1(4) = cmin
                key(1) = 1
                key(2) = 2
                key(3) = 3
                key(4) = 4
                tmp2(1) = cmax
                tmp2(2) = c1
                tmp2(3) = c2
                tmp2(4) = cmin
                call sort(tmp1,key,4)
                cmin = tmp2(key(1))
                c2   = tmp2(key(2))
                c1   = tmp2(key(3))
                cmax = tmp2(key(4))
            endif
        endif
        return
        end
        
       subroutine sort(x,key,n)
c-----
c     Starting with x(1) ,,, x(n)
c     return   the xarray sorted in increasing order
c     also return the pointers key to the initial array. 
c     For example given x = [ 3, 1, 2 ]
c     the returned values are
c                       x = [ 1, 2, 3 ]        
c                     key = [ 2, 3, 1 ]
c-----
c        Reference: http://en.wikipedia.org/wiki/Bubble_sort
c-----
       integer n
       real x(n)
       integer key(n)

       do i=1,n
           key(i) = i
       enddo
       do i = n, 1, -1
           do j = 1 , i -1
               if(x(j) .gt. x(j+1))then
                   tmp = x(j)
                   x(j) = x(j+1)
                   x(j+1) = tmp
                   ktmp = key(j)
                   key(j) = key(j+1)
                   key(j+1) = ktmp
                endif
           enddo
       enddo
       return
       end

        subroutine frstar(r,hs,hr,mname,ipsvsh,time,
     1          SSA, SSC, SSF, SSL, SSN, SSR,
     1          RRA, RRC, RRF, RRL, RRN, RRR,
     1          rayp, geom, tstar, dolock)
c-----
c       r   R   Epicentral distance
c       hs  R   Source depth
c       hr  R   Receiver depth
c       mname   Ch*(*)  Name of model file
c       ipsvsh  I*4 1 - get P time
c               2 - get SV time
c               3 - get SH time
c               4 - get pP time
c               5 - get sP time
c       time    R   First arrival time
c       SSA     R   A at the source
c       SSC     R   C at the source
c       SSF     R   F at the source
c       SSL     R   L at the source
c       SSN     R   N at the source
c       SSR     R - density at the source
c       RRA     R   A at the receiver
c       RRC     R   C at the receiver
c       RRF     R   F at the receiver
c       RRL     R   L at the receiver
c       RRN     R   N at the receiver
c       RRR     R - density at the receiver
c       rayp R   Ray parameter in sec/km
c       geom R   geometrical spreading factor
c       tstar R   geometrical spreading factor
c       dolock L .true. apply locked mode which means to ignore the
c                bottom layer
c-----
        implicit none
        real r, hs, hr, time
        real SSA, SSC, SSF, SSL, SSN, SSR
        real RRA, RRC, RRF, RRL, RRN, RRR
        real rayp, geom, tstar
        logical dolock
        character mname*(*)
        integer ipsvsh
        logical ext
c-----
c-----
c       internal variables
c-----
        real depths, depthr
        real dphs, dphr, dphref
        integer lmaxs, lmaxr, lmaxref

        integer NL
        parameter (NL=200)
        common/timodel/D(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      TRho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL),
     3      frefp(NL), frefs(NL)
        real D, TA, TC, TF, TL, TN, TRho
        real qa, qb, etap, etas, frefp, frefs
        common/modlly/mmax
        integer mmax

        common/depref/refdep
        real refdep

        integer l, lgstr
        integer iunit, iiso, iflsph, idimen, icnvel, ierr
        character title*80 
        
c-----
c       get the earth model
c-----
        inquire(file=mname,exist=ext)
        if(.not. ext)call usage('Model file does not exist')
        l = lgstr(mname)

                call getmod(1,mname,mmax,title,iunit,iiso,iflsph,
     1          idimen,icnvel,ierr,.false.)
        if(ierr .lt. 0)return      
                call tdomod()
c-----
c       insert the source and receiver depths into the model
c       placing the source and receiver on a layer boundary
c-----
        call insert(hs+refdep)
        call insert(hr+refdep)       
        call insert(   refdep)       

c-----
c       get the layer in which the source lies
c-----
        call srclyr(hs+refdep, lmaxs, dphs)
        call srclyr(hr+refdep, lmaxr, dphr)
        call srclyr(   refdep, lmaxref, dphref)

        RRA = TA(lmaxr)
        RRC = TC(lmaxr)
        RRF = TF(lmaxr)
        RRL = TL(lmaxr)
        RRN = TN(lmaxr)
        RRR = TRho(lmaxr)
        SSA = TA(lmaxs)
        SSC = TC(lmaxs)
        SSF = TF(lmaxs)
        SSL = TL(lmaxs)
        SSN = TN(lmaxs)
        SSR = TRho(lmaxs)

c-----
c       compute the travel time
c-----
        call fstarr(r,time,lmaxs, lmaxr, lmaxref,
     1      hs+refdep, hr+refdep, ipsvsh,iflsph, rayp,
     2      tstar, dolock)
        return
        end

        subroutine fstarr(dist,tfirst,lmaxs,lmaxr,lmaxref,
     1      depths,depthr,ipsvsh,iflsph, rayp,
     2      tstar, dolock)
c-----
c       given a distance, the source depth, receiver depth,
c       get time of first arrival of P
c-----
c       dist    R   - distance
c       tfirst  R   - first arrival time
c       mmax    I*4 - number of layers in model
c       lmaxs   I*4 - layer index for source
c       lmaxr   I*4 - layer index for receiver
c       lmaxref I*4 - layer index for reference depth,
c                     used only for pP and sS
c       depths  R   - depth of source
c       depthr  R   - depth of receiver
c       ipsvsh  I*4 1 - get P time
c               2 - get SV time
c               3 - get SH time
c               4 - get pP time
c               5 - get sP time
c       iflsph  I*4 - 0 flat earth model
c                 1 spherical earth model
c       rayp    R   - ray parameter in sec/km
c       geom R   geometrical spreading factor
c       dolock L .true. apply locked mode which means to ignore the
c                bottom layer
c-----
c
c       18 JAN 2008 - everything is straightforward. The addition of
c          the request for pP and sP changes the logic in that
c          the direct arrival is ignored, and that the upgoing refraction 
c          from the source is ignored. We handle this by just setting
c          a very large tfirst before trying to do the modified 
c          downward path refraction to avoid another level of
c          if/then/else/endif
c-----
        real dist, tfirst, depths, depthr
        real rayp
        integer lmaxs, lmaxr, lmaxref, ipsvsh, iflsph
        logical dolock

        integer NL
        parameter (NL=200)
        common/timodel/TD(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      TRho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL),
     3      frefp(NL), frefs(NL)
        real TD, TA, TC, TF, TL, TN, TRho
        real qa, qb, etap, etas, frefp, frefs
        common/modlly/mmmax
        integer mmmax
        common/depref/refdep
        real refdep
        common/shwave/TLsh(NL), TNsh(NL), TRhosh(NL),qbsh(NL)
        real TLsh, TNsh, TRhosh,qbsh

        integer mmax

        real*4  h(NL)

        real*8   pupper
        complex*16 p
        integer lmx, lmn
        integer i, l
        real sumx, sumt, tt
        real time

        COMPLEX*16 A,B,C,D,E,F
        COMPLEX*16 wvno2, omega2, wvn, omg
        COMPLEX*16 rsh, rp, rsv
c-----
c       get the modified eigen vectors x11 and x31 
c           have lambda1 (rp ) factored out
c               modified eigen vectors x22 and x42 
c           have lambda2 (rsv) factored out
c-----
        COMPLEX*16 X11, X21, X31, X41
        COMPLEX*16 X12, X22, X32, X42
        COMPLEX*16 NP, NSV
        real getvel

        COMPLEX*16 dtdp
        complex*16 pold, pcur, dp
        complex*16 dtdpold, dtdpcur
        integer ilast

        logical baseisp, layerisp

        common/earth/radius

c-----
c       initialize
c-----
        omg = dcmplx(1.0d+00, 0.0d+00)
        omega2 = omg *omg
        tstar = 0.0

c-----
c       set up default
c-----
        tfirst = 1.0e+30
c-----
c       special case for locked mode
c-----
        if(dolock)then
            mmax = mmmax -1
        else
            mmax = mmmax
        endif

c-----
c       get specifics about upward and downward distances
c       with a layer. We need his to define ray paths
c       We will also use the fact that the source/receiver are
c       on layer boundaries
c
c       lmn = layer number of shallowest of source/receiver
c       lmx = layer number of deepest    of source/receiver
c-----
        lmn = min(lmaxs,lmaxr)
        lmx = max(lmaxs,lmaxr)

c-----
c       perform spherical -> flat earth model transformation
c-----
        if(iflsph.ne.0)then
            call tdosph()
        endif
c-----
c       now fill in velocity array according to desired first arrival
c       for SH there can be no water layer
c       for SV can be a water layer
c       Also define the Q for the T* analysis. Note we define
c        eventually q = 1/Q based on whether the given Q > or < 1
c-----
        do i=1,mmax
            if(qa(i) .gt. 1.0)then
                qa(i) = 1.0 / qa(i)
            endif
            if(qb(i) .gt. 1.0)then
                qb(i) = 1.0 / qb(i)
            endif
            h(i) = td(i)
        enddo

c-----
c       For the computations we look at four cases
c       1) direct path between source and receiver 
c       2) refracted arrivals       
c          a) path is downward from source and then up to
c             receiver
c          b) path is upward from the source and then down to
c             receiver
c          This recognized the possibility that velocity does
c          not increase uniformly with depth
c-----
                    
c-----
c       direct arrival 
c-----
c       Newton Iteration for direct arrival source/receiver at
c           different depths
c           
c           x = SUM h  tan theta
c                    i          i
c
c           t = SUM h  / v  cos theta
c                    i    i          i
c                                                          2 2
c       where sin theta  = p V  , cos theta  = sqrt ( 1 - p V )
c                      i      i                              i
c       and p is the ray parameter bounded by [0, 1/V  ] where V
c                                                    sr         sr
c       is the wave velocity at the starting point of the ray. 
c       Since the ray must also reach the receiver, we consider
c       that possibility too. The upper bound is MIN ( 1/Vs, 1/Vr)
c       Also we test for a real ray path, between source and receiver
c
c       Because source/receiver at top of a layer boundary, we have
c
c           -----------X----------
c           h(lmn)      \
c           ----------------------
c                      ....
c           ----------------------
c           h(lmx-1)        \
c                            \
c           ------------------X---
c            
c-----
c          reflection occurs when dt/dp = 0, so search for the p value 
c          numerically. The travel time is just
c               t = p r + Sum eta h
c-----
            ps = 1.0/getvel(TA,TL,TN,TRho,lmaxs,ipsvsh)
            pr = 1.0/getvel(TA,TL,TN,TRho,lmaxr,ipsvsh)
            if(ps.lt.pr)then
                pupper = ps
            else
                pupper = pr
            endif
            do 1000 l=lmn,lmx
                vl = getvel(TA,TL,TN,TRho,l,ipsvsh)
                if(vl.eq.0.0)return
                p = dcmplx(1.0/vl, 0.0d+00)
                if(dreal(p).lt.pupper)pupper = p
 1000       continue
            pold = dcmplx(0.0d+00, 0.0d+00)
            dp =  dcmplx(pupper/100.0, 0.0d+00)
            do  i=0,100
                ilast = i
                p = i*dp
                wvn = p * omg
                wvno2 = wvn * wvn
            call gttdtdp(omega2,wvno2,omg,wvn,p,ipsvsh,
     1          dist,lmn,lmx-1,time,dtdp,tstar)
            if(dreal(dtdp).lt. 0.0d+00)then
c----
c                       refine the root
c-----
                        pold = p - dp
                        pcur = p
                        dtdpcur = dtdp
                ilast = i -1
                        go to 2000
                endif
                dtdpold = dtdp
        enddo
c-----
c       assume we always get here
c-----
 2000   continue
c-----
c       use interval halving until I can compute the d2t/dp2!
c       also as a fallback, do not do this if the maximum index about was 100
c-----
        if(ilast.ne.100)then
        do 3000 i=1,10
            p = 0.5*(pold + pcur)
            wvn = p*omg
            wvno2 = wvn*wvn
            call gttdtdp(omega2,wvno2,omg,wvn,p,ipsvsh,
     1          dist,lmn,lmx-1,time,dtdp,tstar)
            if(dsign(1.0d+00,dreal(dtdpcur)).eq.
     1          dsign(1.0d+00,dreal(dtdp)))then
                pcur = p
                dtdpcur = dtdp
            else
                pold = pcur
                dtdpold = dtdp
            endif

 3000       continue
        endif
            tfirst = time
            rayp = dreal(p)
c-----
c       now proceed through the possible refracted arrivals
c       considering first upward rays from the source
c-----  
        if(lmn.gt.1)then
        do 3020 m=1,lmn-1
c-----
c       m is the refracting layer
c
c       get velocity of refractor testing if fluid
c-----
            if(ipsvsh.eq.1 .or. ipsvsh.eq.4 .or.ipsvsh.eq.5)then
                 vel = getvel(TA,TL,TN,TRho,m,1)
                 if(vel.eq.0.0)go to 3040
                 baseisp = .true.
            else  
                 vel = getvel(TA,TL,TN,TRho,m,ipsvsh)
                 if(vel.eq.0.0.and.ipsvsh.eq.2)then
c-----
c                     permit S-P conversion, e.g., SKS
c-----
                      vel = getvel(TA,TL,TN,TRho,m,1)
                      baseisp = .true.
                 else
                      baseisp = .false.
                 endif
            endif
            if(vel.eq.0.0)goto 3040
            p = 1.0/vel
            wvn = p * omg
            wvno2 = wvn * wvn
c-----
c
c           --------------------------------
c           h(1)
c           --------------------------------
c                      ....
c           --------------------------------
c           h(m)
c           ----------------...-------------
c           h(m+1)         /   \
c           --------------------------------
c                         /     \
c                      ....
c           --------------------------------
c           h(lmn-1)              \
c           -----------------------X--------
c               
c           h(lmn)     /    
c           --------------------------------
c                      ....
c           --------------------------------
c           h(lmx-1) /
c           --------X-----------------------
c
c       safety check, velocity at source or receiver must be less than
c       refraction velocity
c-----
        if(ipsvsh.le.3)then
             vlmn = getvel(TA,TL,TN,TRho,lmn,ipsvsh)
             vlmx = getvel(TA,TL,TN,TRho,lmx,ipsvsh)
        else
c-----
c       this is just a subterfuge since we will not use the results
c-----
             vlmn = getvel(TA,TL,TN,TRho,lmn,1)
             vlmx = getvel(TA,TL,TN,TRho,lmx,1)
        endif
        if(vlmx.ge.vel)go to 3020
        if(vlmn.ge.vel)go to 3020
c-----
c       single leg
c-----
            sumx = 0.0
            sumt = 0.0
            ts = 0.0
            do 3021 l=lmn,lmx-1
                if(ipsvsh.eq.1 .or. ipsvsh.eq.4 .or.ipsvsh.eq.5)then
                     vl = getvel(TA,TL,TN,TRho,l,1)
                     layerisp = .true.
                else  
                     vl = getvel(TA,TL,TN,TRho,l,ipsvsh)
                     layerisp = .false.
                     if(vl.eq.0.0.and.ipsvsh.eq.2)then
c-----
c                         permit S-P conversion, e.g., SKS
c-----
                          vl = getvel(TA,TL,TN,TRho,l,1)
                          layerisp = .true.
                     endif
                endif

                if(vl.gt.vel)go to 3020
                if(vl.eq.0.0)go to 3040
                call tgetegn(a,b,c,d,e,f,omega2,wvno2,rp, rsv, rsh,
     1              x11,x21,x31,x41,x12,x22,x32,x42,NP, NSV,l,omg,wvn)
                if(ipsvsh.eq.1.or.ipsvsh.eq.4.or.ipsvsh.eq.5)then
C                    if(dimag(rp).lt. 0.0d+00)then
C                         rp = - rp
C                    endif
                    sumt = sumt + h(l) * rp /dcmplx(0.0d+00,1.0d+00)
                    sumx = sumx + h(l) * 
     1                  p/(rp/dcmplx(0.0d+00, 1.0d+00))
                    ts = ts + qa(l)*h(l)*(p*p - rp*rp)
     1                   /(rp/dcmplx(0.0d+00, 1.0d+00))
                else if(ipsvsh.eq.2)then
C                    if(dimag(rsv).lt. 0.0d+00)then
C                         rsv = - rsv
C                    endif
                     if(layerisp)then
                    sumt = sumt + h(l) * rp /dcmplx(0.0d+00,1.0d+00)
                    sumx = sumx + h(l) * 
     1                  p/(rp/dcmplx(0.0d+00, 1.0d+00))
                    ts = ts + qa(l)*h(l)*(p*p - rp*rp)
     1                   /(rp/dcmplx(0.0d+00, 1.0d+00))
                      else
                    sumt = sumt + h(l) * rsv/dcmplx(0.0d+00,1.0d+00)
                    sumx = sumx + h(l) * 
     1                  p/(rsv/dcmplx(0.0d+00, 1.0d+00))
                    ts = ts + qb(l)*h(l)*(p*p - rsv*rsv)
     1                   /(rsv/dcmplx(0.0d+00, 1.0d+00))
                      endif
                else if(ipsvsh.eq.3)then
C                    if(dimag(rsh).lt. 0.0d+00)then
C                         rsh = - rsh
C                    endif
                    sumt = sumt + h(l) * rsh/dcmplx(0.0d+00,1.0d+00)
                    sumx = sumx + h(l) * 
     1                  p/(rsh/dcmplx(0.0d+00, 1.0d+00))
                    ts = ts + qb(l)*h(l)*(p*p - rsh*rsh)
     1                   /(rsh/dcmplx(0.0d+00, 1.0d+00))
                endif
 3021       continue
            do 3022 l=m+1,lmn-1
                if(ipsvsh.eq.1 .or. ipsvsh.eq.4 .or.ipsvsh.eq.5)then
                     vl = getvel(TA,TL,TN,TRho,l,1)
                     layerisp = .true.
                else  
                     vl = getvel(TA,TL,TN,TRho,l,ipsvsh)
                     layerisp = .false.
                     if(vl.eq.0.0.and.ipsvsh.eq.2)then
c-----
c                         permit S-P conversion, e.g., SKS
c-----
                          vl = getvel(TA,TL,TN,TRho,l,1)
                          layerisp = .true.
                     endif
                endif
                if(vl.gt.vel)go to 3020
                if(vl.eq.0.0)go to 3040
                call tgetegn(a,b,c,d,e,f,omega2,wvno2,rp, rsv, rsh,
     1              x11,x21,x31,x41,x12,x22,x32,x42,NP, NSV,l,omg,wvn)
                if(ipsvsh.eq.1.or.ipsvsh.eq.4.or.ipsvsh.eq.5)then
C                    if(dimag(rp).lt. 0.0d+00)then
C                         rp = - rp
C                    endif
                    sumt = sumt + 2.0*h(l) * rp /dcmplx(0.0d+00,1.0d+00)
                    sumx = sumx + 2.0*h(l) * 
     1                  p/(rp/dcmplx(0.0d+00, 1.0d+00))
                    ts = ts + 2.*qa(l)*h(l)*(p*p - rp*rp)
     1                   /(rp/dcmplx(0.0d+00, 1.0d+00))
                else if(ipsvsh.eq.2)then
C                    if(dimag(rsv).lt. 0.0d+00)then
C                         rsv = - rsv
C                    endif
                     if(layerisp)then
                    sumt = sumt + 2.0*h(l) * rp /dcmplx(0.0d+00,1.0d+00)
                    sumx = sumx + 2.0*h(l) * 
     1                  p/(rp/dcmplx(0.0d+00, 1.0d+00))
                    ts = ts + 2.*qa(l)*h(l)*(p*p - rp*rp)
     1                   /(rp/dcmplx(0.0d+00, 1.0d+00))
                      else
                    sumt = sumt + 2.0*h(l) * rsv/dcmplx(0.0d+00,1.0d+00)
                    sumx = sumx + 2.0*h(l) * 
     1                  p/(rsv/dcmplx(0.0d+00, 1.0d+00))
                    ts = ts + 2.*qb(l)*h(l)*(p*p - rsv*rsv)
     1                   /(rsv/dcmplx(0.0d+00, 1.0d+00))
                      endif
                else if(ipsvsh.eq.3)then
C                    if(dimag(rsh).lt. 0.0d+00)then
C                         rsh = - rsh
C                    endif
                    sumt = sumt + 2.0*h(l) * rsh/dcmplx(0.0d+00,1.0d+00)
                    sumx = sumx + 2.0*h(l) * 
     1                  p/(rsh/dcmplx(0.0d+00, 1.0d+00))
                    ts = ts + 2.*qb(l)*h(l)*(p*p - rsh*rsh)
     1                   /(rsh/dcmplx(0.0d+00, 1.0d+00))
                endif
 3022       continue
            tint = sumt
            tt = tint + dist / vel
            if(baseisp)then
                 ts = ts + qa(m)*(dist-sumx)/vel
            else
                 ts = ts + qb(m)*(dist-sumx)/vel
            endif
            if(tt .lt. tfirst .and. dist.ge.sumx)then
                  tfirst = tt
                  rayp = dreal(p)
                 tstar = ts
            endif
 3020       continue
 3040       continue
        endif
c-----
c       For the special case of the depth phases, ignore previous
c       first arrival times
c-----
        if(ipsvsh.eq.4 .or. ipsvsh.eq.5)then
             tfirst = 1.0e+30
        endif
c-----
c       now proceed through the possible refracted arrivals
c       considering first downward rays from the source
c
c       We start considering the deepest point since we place
c       a source/receiver position just below a layer boundary
c       and thus should consider a horizontal ray
c
c       The refraction is accepted only if the desired distance >
c       first refraction from the source - this puts physics in the problem
c           
c           x = SUM h  tan theta
c                    i          i
c
c           t = SUM h  cos theta / V
c                    i          i   i
c                                                          2 2
c       where sin theta  = p V  , cos theta  = sqrt ( 1 - p V )
c                      i      i                              i
c       For the T* computation we need to follow the path, e.g.,
c       SUM h qi / ( cos theta  / V ) + qi (dist -  SUM h tan theta / V )/V
c            i  i             i    i      i              i         i   i   r
c-----  
        do 2020 m=lmx+1, mmax
c-----
c       m is the refracting layer
c
c       get velocity of refractor testing if fluid
c-----
            if(ipsvsh.eq.1 .or. ipsvsh.eq.4 .or.ipsvsh.eq.5)then
                 vel = getvel(TA,TL,TN,TRho,m,1)
                 if(vel.eq.0.0)go to 2040
                 baseisp = .true.
            else  
                 vel = getvel(TA,TL,TN,TRho,m,ipsvsh)
                 if(vel.eq.0.0.and.ipsvsh.eq.2)then
c-----
c                     permit S-P conversion, e.g., SKS
c-----
                      vel = getvel(TA,TL,TN,TRho,m,1)
                      baseisp = .true.
                 else
                      baseisp = .false.
                 endif
            endif
            if(vel.eq.0.0)go to 2040
            p = 1.0/vel
            wvn = p * omg
            wvno2 = wvn * wvn
c-----
c
c           -----------X--------------------
c           h(lmn)      \
c           --------------------------------
c                      ....
c           --------------------------------
c           h(lmx-1)        \             
c                            \           
c           ------------------X--------X----
c           h(lmx)             \       /
c           --------------------\-----/-----
c                      ....      \   /
c           ----------------------...-------
c           h(m)
c
c-----
c       safety check, velocity at source or receiver must be less than
c       refraction velocity otherwise there will be no real ray
c-----
        if(ipsvsh.le.3)then
             vlmn = getvel(TA,TL,TN,TRho,lmn,ipsvsh)
             vlmx = getvel(TA,TL,TN,TRho,lmx,ipsvsh)
        else
c-----
c            this is just a subterfuge since 
c            we will not use the results for pP sP
c-----
             vlmn = getvel(TA,TL,TN,TRho,lmn,1)
             vlmx = getvel(TA,TL,TN,TRho,lmx,1)
        endif
        if(vlmx.ge.vel)go to 2020
        if(vlmn.ge.vel)go to 2020
c-----
c       single leg
c-----
            sumx = 0.0
            sumt = 0.0
            ts = 0.0
c-----
c       special case for depth phases
c-----
            if(ipsvsh.eq.4)then
c-----
c               pP
c-----
                  do  l=lmaxref,lmaxs - 1
                      vp = getvel(TA,TL,TN,TRho,l,1)
                      if(vp.gt.vel)go to 2020
                      if(vp.eq.0.0)go to 2040
                call tgetegn(a,b,c,d,e,f,omega2,wvno2,rp, rsv, rsh,
     1              x11,x21,x31,x41,x12,x22,x32,x42,NP, NSV,l,omg,wvn)
                    sumt = sumt + 2.*h(l) * rp /dcmplx(0.0d+00, 1.0d+00)
                    sumx = sumx + 2.*h(l) *
     1                  p/(rp/dcmplx(0.0d+00, 1.0d+00))
                    ts = ts + qa(l)*2.*h(l)*(p*p - rp*rp)
     1                   /(rp/dcmplx(0.0d+00, 1.0d+00))
                  enddo
            else if(ipsvsh.eq.5)then
c-----
c               sP
c-----
                  do  l=lmaxref,lmaxs - 1
                      vp = getvel(TA,TL,TN,TRho,l,1)
                      vs = getvel(TA,TL,TN,TRho,l,2)
                      if(vp.gt.vel)go to 2020
                      if(vp.eq.0.0)go to 2040
                      if(vs.gt.vel)go to 2020
                      if(vs.eq.0.0)go to 2040
                call tgetegn(a,b,c,d,e,f,omega2,wvno2,rp, rsv, rsh,
     1              x11,x21,x31,x41,x12,x22,x32,x42,NP, NSV,l,omg,wvn)
                    sumt = sumt + h(l) * rp /dcmplx(0.0d+00, 1.0d+00)
     1                          + h(l) * rsv/dcmplx(0.0d+00, 1.0d+00)
                    sumx = sumx 
     1                  + h(l)*p/(rp/dcmplx(0.0d+00, 1.0d+00))
     1                  + h(l)*p/(rsv/dcmplx(0.0d+00, 1.0d+00))
                    ts = ts + qa(l)*h(l)*(p*p - rp*rp)
     1                   /(rp/dcmplx(0.0d+00, 1.0d+00))
     1                      + qb(l)*h(l)*(p*p - rsv*rsv)
     1                   /(rsv/dcmplx(0.0d+00, 1.0d+00))
                  enddo
            endif
            do 2021 l=lmn,lmx - 1
                if(ipsvsh.eq.1 .or. ipsvsh.eq.4 .or.ipsvsh.eq.5)then
                     vl = getvel(TA,TL,TN,TRho,l,1)
                     layerisp = .true.
                else  
                     vl = getvel(TA,TL,TN,TRho,l,ipsvsh)
                     layerisp = .false.
                     if(vl.eq.0.0.and.ipsvsh.eq.2)then
c-----
c                         permit S-P conversion, e.g., SKS
c-----
                          vl = getvel(TA,TL,TN,TRho,l,1)
                          layerisp = .true.
                     endif
                endif

                if(vl.gt.vel)go to 2020
                if(vl.eq.0.0)go to 2040
                call tgetegn(a,b,c,d,e,f,omega2,wvno2,rp, rsv, rsh,
     1              x11,x21,x31,x41,x12,x22,x32,x42,NP, NSV,l,omg,wvn)
c-----
c      KLUDGE - to fix the case when the imaginary part of the wavenumber is
c            negative - it must be positive
c-----
                if(ipsvsh.eq.1.or.ipsvsh.eq.4.or.ipsvsh.eq.5)then
C                    if(dimag(rp).lt. 0.0d+00)then
C                         rp = - rp
C                    endif
                    sumt = sumt + h(l) * rp /dcmplx(0.0d+00, 1.0d+00)
                    sumx = sumx + h(l) * 
     1                  p/(rp/dcmplx(0.0d+00, 1.0d+00))
                    ts = ts + qa(l)*h(l)*(p*p - rp*rp)
     1                   /(rp/dcmplx(0.0d+00, 1.0d+00))
                else if(ipsvsh.eq.2)then
C                    if(dimag(rsv).lt. 0.0d+00)then
C                         rsv = - rsv
C                    endif
                     if(layerisp)then
                    sumt = sumt + h(l) * rp /dcmplx(0.0d+00,1.0d+00)
                    sumx = sumx + h(l) * 
     1                  p/(rp/dcmplx(0.0d+00, 1.0d+00))
                    ts = ts + qa(l)*h(l)*(p*p - rp*rp)
     1                   /(rp/dcmplx(0.0d+00, 1.0d+00))
                      else
                    sumt = sumt + h(l) * rsv/dcmplx(0.0d+00,1.0d+00)
                    sumx = sumx + h(l) * 
     1                  p/(rsv/dcmplx(0.0d+00, 1.0d+00))
                    ts = ts + qb(l)*h(l)*(p*p - rsv*rsv)
     1                   /(rsv/dcmplx(0.0d+00, 1.0d+00))
                      endif
                else if(ipsvsh.eq.3)then
C                    if(dimag(rsh).lt. 0.0d+00)then
C                         rsh = - rsh
C                    endif
                    sumt = sumt + h(l) * rsh/dcmplx(0.0d+00, 1.0d+00)
                    sumx = sumx + h(l) * 
     1                  p/(rsh/dcmplx(0.0d+00, 1.0d+00))
                    ts = ts + qb(l)*h(l)*(p*p - rsh*rsh)
     1                   /(rsh/dcmplx(0.0d+00, 1.0d+00))
                endif
 2021       continue
c-----
c       double leg
c-----

            do 2022 l=lmx,m-1
                if(ipsvsh.eq.1 .or. ipsvsh.eq.4 .or.ipsvsh.eq.5)then
                     vl = getvel(TA,TL,TN,TRho,l,1)
                     layerisp = .true.
                else  
                     vl = getvel(TA,TL,TN,TRho,l,ipsvsh)
                     layerisp = .false.
                     if(vl.eq.0.0.and.ipsvsh.eq.2)then
c-----
c                         permit S-P conversion, e.g., SKS
c-----
                          vl = getvel(TA,TL,TN,TRho,l,1)
                          layerisp = .true.
                     endif
                endif

                if(vl.gt.vel)go to 2020
                if(vl.eq.0.0)go to 2040
                call tgetegn(a,b,c,d,e,f,omega2,wvno2,rp, rsv, rsh,
     1              x11,x21,x31,x41,x12,x22,x32,x42,NP, NSV,l,omg,wvn)
                if(ipsvsh.eq.1.or.ipsvsh.eq.4.or.ipsvsh.eq.5)then
C                    if(dimag(rp).lt. 0.0d+00)then
C                         rp = - rp
C                    endif
                    sumt = sumt + 2.0*h(l) * rp /dcmplx(0.0d+00,1.0d+00)
                    sumx = sumx + 2.0*h(l) * 
     1                  p/(rp/dcmplx(0.0d+00, 1.0d+00))
                    ts = ts + 2.*qa(l)*h(l)*(p*p - rp*rp)
     1                   /(rp/dcmplx(0.0d+00, 1.0d+00))
                else if(ipsvsh.eq.2)then
C                    if(dimag(rsv).lt. 0.0d+00)then
C                         rsv = - rsv
C                    endif
                     if(layerisp)then
                    sumt = sumt + 2.0*h(l) * rp /dcmplx(0.0d+00,1.0d+00)
                    sumx = sumx + 2.0*h(l) * 
     1                  p/(rp/dcmplx(0.0d+00, 1.0d+00))
                    ts = ts + 2.*qa(l)*h(l)*(p*p - rp*rp)
     1                   /(rp/dcmplx(0.0d+00, 1.0d+00))
                      else
                    sumt = sumt + 2.0*h(l) * rsv/dcmplx(0.0d+00,1.0d+00)
                    sumx = sumx + 2.0*h(l) * 
     1                  p/(rsv/dcmplx(0.0d+00, 1.0d+00))
                    ts = ts + 2.*qb(l)*h(l)*(p*p - rsv*rsv)
     1                   /(rsv/dcmplx(0.0d+00, 1.0d+00))
                      endif
                else if(ipsvsh.eq.3)then
C                    if(dimag(rsh).lt. 0.0d+00)then
C                         rsh = - rsh
C                    endif
                    sumt = sumt + 2.0*h(l) * rsh/dcmplx(0.0d+00,1.0d+00)
                    sumx = sumx + 2.0*h(l) * 
     1                  p/(rsh/dcmplx(0.0d+00, 1.0d+00))
                    ts = ts + 2.*qb(l)*h(l)*(p*p - rsh*rsh)
     1                   /(rsh/dcmplx(0.0d+00, 1.0d+00))
                endif
 2022       continue
            tint = sumt
            tt = tint + dist / vel
            if(baseisp)then
                 ts = ts + qa(m)*(dist-sumx)/vel
            else
                 ts = ts + qb(m)*(dist-sumx)/vel
            endif
            if(tt .lt. tfirst .and. dist.ge.sumx)then
                 tfirst = tt
                 rayp = dreal(p)
                 tstar = ts
            endif

                 vp = getvel(TA,TL,TN,TRho,m,1)
                 vsv = getvel(TA,TL,TN,TRho,m,2)
                 vsh = getvel(TA,TL,TN,TRho,m,3)
 2020       continue
 2040       continue
             if(tfirst .eq. 1.0e+30)then
                tfirst = -12345.
                tstar  = -12345.
                rayp   = -12345.
             endif
        return
        end

        subroutine tdosph()
c-----
c       Transform spherical earth to flat earth
c
c       Schwab, F. A., and L. Knopoff (1972). 
c           Fast surface wave and free
c       mode computations, in  
c           Methods in Computational Physics, Volume 11,
c       Seismology: Surface Waves and Earth Oscillations,  
c           B. A. Bolt (ed),
c       Academic Press, New York
c
c       Love Wave Equations  44, 45 , 41 pp 112-113
c       Rayleigh Wave Equations 102, 108, 109 pp 142, 144
c
c
c       We will treat all as P-SV for the heck of it
c       This requires more work
c-----
c       mmax    I*4 number of layers
c       TA     R   A 
c       TC     R   C 
c       TF     R   F 
c       TL     R   L 
c       TN     R   N 
c                  note  density not required
c       TD     R   layer thickness
c       v() R   array of velocities
c       h() R   array of layer thicknesses
c       ipsvsh  I       1 - get P time
c                       2 - get SV time
c                       3 - get SH time
c       refdep R   Reference depth for the model specification
c
c       Note we need the constants here.  Since the velocities
c       must increase with depth, e.g., vf = vs (a/r)
c       and that density  varies
c       as rhof = rhos (a/r)^-P, [not the TI surface wave code has not yet
c        been written], then using the model that m = rho beta^2, we have
c
c       TA = rho VA^2,
c       TAf = rhof * VAf^2 = rhos (a/r)^-P VAs^2 (a/r)^2
c           = (a/r)^2-P TAs
c-----
        integer NL
        parameter (NL=200)
        common/timodel/TD(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      TRho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL),
     3      frefp(NL), frefs(NL)
        real TD, TA, TC, TF, TL, TN, TRho
        real qa, qb, etap, etas, frefp, frefs
        common/modlly/mmax
        integer mmax
        common/depref/refdep
        real refdep
        common/shwave/TLsh(NL), TNsh(NL), TRhosh(NL), qbsh(NL)
        real TLsh, TNsh, TRhosh, qbsh


        double precision z0,z1,r0,r1,ar,tmp

        common/earth/radius
        real radius

        ar=radius
        r0=ar + refdep
        td(mmax)=1.0
        do 10 i=1,mmax
            r1=r0-dble(td(i))
            z0=ar*dlog(ar/r0)
            z1=ar*dlog(ar/r1)
            td(i)=z1-z0
c-----
c        attempt 7 15 2007 - use standard rule but at mid layer depth as per DGH
c-----
            TMP=(ar+ar)/(r0+r1)
c-----
c                SV
c-----
                 rhosph    = trho(i)
                 trho(i)   = rhosph * tmp**(-2.275)
                 trhosh(i) = rhosph * tmp**(-5)

                 ta(i)=ta(i)*tmp**(-0.2750)
                 tc(i)=tc(i)*tmp**(-0.2750)
                 tf(i)=tf(i)*tmp**(-0.2750)

                 elsph = tl(i)
                 tl(i)  =elsph*tmp**(-0.2750)
                 tlsh(i)=elsph*tmp**(-3.0)
                 ensph = tn(i)

                 tn(i)=ensph*tmp**(-0.2750)
                 tnsh(i)=ensph*tmp**(-3.0)
            r0 = r1
   10   continue
        td(mmax)=0.0
        return
        end

        subroutine tdomod()
c-----
c       just fill in the TRhosh, TLsh, TNsh and qbsh arrays
c-----
        integer NL
        parameter (NL=200)
        common/timodel/TD(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      TRho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL),
     3      frefp(NL), frefs(NL)
        real TD, TA, TC, TF, TL, TN, TRho
        real qa, qb, etap, etas, frefp, frefs
        common/modlly/mmax
        integer mmax
        common/depref/refdep
        real refdep
        common/shwave/TLsh(NL), TNsh(NL), TRhosh(NL), qbsh(NL)
        real TLsh, TNsh, TRhosh, qbsh

        do i=1,mmax
           TLsh(i) = TL(i)
           TNsh(i) = TN(i)
           TRhosh(i) = TRho(i)
           qbsh(i) = qbsh(i)
        enddo
        return
        end

        subroutine srclyr(depth,lmax,dph)
        parameter (LER=0, LIN=5, LOT=6)
        parameter(NL=200)
        common/timodel/D(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      TRho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL), 
     3      frefp(NL), frefs(NL)
        real D, TA, TC, TF, TL, TN, TRho
        real qa, qb, etap, etas, frefp, frefs
        common/modlly/mmax
c-----
c       Find source/receiver boundary. It is assumed that
c       it will lie upon a boundary
c
c       lmax = source layer 
c            = 0 is the free surface 
c       depth = source depth 
c       dph = height of  source above lmax + 1 interface 
c       lmax = 0 is the free surface 
c-----
        if(depth.le.0.0)then
            lmax = 1
            dph = 0.0
        else
            dep = 0.0 
            do 100 m = 2,mmax
                dep = dep + d(m-1) 
                dph = dep - depth 
                lmax = m 
                if(abs(dph).lt. 0.0001*d(m-1) .or.
     1              abs(dph).lt.1.0e-6)go to 101
  100       continue 
  101   continue 
        endif
        return 
        end 

        subroutine getabc(m,omg,wvn,a,b,c,d,e,f)
        implicit none
        integer m
        COMPLEX*16 omg,wvn
        COMPLEX*16 a, b, c, d, e, f
        integer NL
        parameter (NL=200)
        common/timodel/H(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      TRho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL), 
     3      frefp(NL), frefs(NL)
        real H, TA, TC, TF, TL, TN, TRho
        real qa, qb, etap, etas, frefp, frefs
        a = wvn * TF(m) / TC(m)
        b = 1.0/TC(m)
        c = - TRho(m)*omg*omg + wvn*wvn *(TA(m) -TF(m)*TF(m)/TC(m))
        d = - wvn
        e = 1.0/TL(m)
        f = - TRho(m)*omg*omg
        return
        end                                               

        subroutine tgetegn(a,b,c,d,e,f,omega2,wvno2,rp, rsv, rsh,
     1      x11,x21,x31,x41,x12,x22,x32,x42,NP, NSV,m,omg,wvn)
        implicit none
        COMPLEX*16 A,B,C,D,E,F
        COMPLEX*16 wvno2, omega2
        COMPLEX*16 rsh, rp, rsv
c-----
c       get the modified eigen vectors x11 and x31 have 
c           lambda1 (rp ) factored out
c               modified eigen vectors x22 and x42 have 
c           lambda2 (rsv) factored out
c-----
        COMPLEX*16 X11, X21, X31, X41
        COMPLEX*16 X12, X22, X32, X42
c-----
c       norms
c-----
        COMPLEX*16 NP, NSV
        integer m
        COMPLEX*16 omg, wvn
        COMPLEX*16 xka2, xkb2

        integer NL
        parameter (NL=200)
        common/timodel/H(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      TRho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL), 
     3      frefp(NL), frefs(NL)
        real H, TA, TC, TF, TL, TN, TRho
        real qa, qb, etap, etas, frefp, frefs
        common/shwave/TLsh(NL), TNsh(NL), TRhosh(NL), qbsh(NL)
        real TLsh, TNsh, TRhosh, qbsh
c-----
c       internal variables
c-----
        COMPLEX*16 L2(2)
        COMPLEX*16 bb, cc
        COMPLEX*16 CDSQRT 

        COMPLEX*16 ZFAC
c-----
c       first test to see if a fluid layer - if it is fluid, the
c       eigenfunctions are specially computed and we need only the
c       rp
c-----
        if(TL(m).eq.0.0 .or. TN(m).eq.0.0)then
            rp = cdsqrt(wvno2 -omega2*TRho(m)/TA(m))
            rsv = dcmplx(0.0d+000, 0.0d+00)
            rsh = dcmplx(0.0d+000, 0.0d+00)
            return
        endif
        
        call getabc(m,omg,wvn,a,b,c,d,e,f)
c-----
c       Do the SH
c-----
        rsh = CDSQRT(TNsh(m)*wvno2/TLsh(m) - Trhosh(m)*omega2/TLsh(m)) 
        if( dimag(rsh) .lt. 0.0)then
                rsh = - rsh
        endif
c-----
c       Do the P and SV
c-----
c-----
c       The characteristic equation to be solved is
c
c       L^4 + L^2[ -2 ad -ec -fb ] + [ (d^2+ef)(a^2+bc)] = 0
c-----
        bb = -2.0d+00 * a*d - e*c -f*b
        cc = ( d*d + e*f)*(a*a + b*c)
        L2(1) = ( - bb + CDSQRT(bb*bb - 4.000*cc))/2.0d+00
        L2(2) = ( - bb - CDSQRT(bb*bb - 4.000*cc))/2.0d+00

        L2(1) = cc/L2(2)
c-----
c       Use the Lambda^2 values to form
c       xka^2 == k^2 - L(1)^2
c       xkb^2 == k^2 - L(2)^2
c       Associate the smallest xka, xkb with the P!
c-----
        xka2 = wvno2 - L2(1)
        xkb2 = wvno2 - L2(2)
        if(cdabs(xkb2) .lt. cdabs(xka2))THEN
                ZFAC = L2(1)
                L2(1) = L2(2)
                L2(2) = ZFAC
        endif
        rp  = CDSQRT(L2(1))
        rsv = CDSQRT(L2(2))
        if( dimag(rp) .lt. 0.0)then
                rp = - rp
        endif
        if( dimag(rsv) .lt. 0.0)then
                rsv = - rsv
        endif
c-----
c       get the norms - note that the true norm will be 
c           2  NP amd 2 L(2) NSV
c       The factorization permits us to use the sin nz/n or n sin nz
c-----
        NP  = (  L2(1)*(-2*a*b*d + 2*a*a*e + b*c*e - b*b*f)
     1      + (a*a+b*c)*(2*b*d*d - 2*a*d*e + b*e*f - c*e*e) )
        NSV = (- L2(2)*(2*b*d*d - 2*a*d*e - c*e*e + b*e*f)
     1      + (d*d+e*f)*(2*a*b*d - 2*a*a*e + b*b*f - b*c*e) )
c-----
c       get the modified eigen vectors x11 and x31 have 
c           lambda1 (rp ) factored out
c               modified eigen vectors x22 and x42 have 
c           lambda2 (rsv) factored out
c-----
        x11 =              (b*d - a*e)
        x21 =  b*L2(1) - e*(b*c + a*a)
        x31 =    L2(1) -   (a*d + c*e)
        x41 = -a*L2(1) + d*(b*c + a*a)

        x12 = -e*L2(2) + b*(d*d + e*f)
        x22 = ( b*d - a*e)
        x32 = d*L2(2) - a*(d*d + e*f)
        x42 = - ( L2(2) -  a*d - b*f)
c-----
c       TEST
c       Force the eigenfunctions to be as given in 5.4.4
c-----
        zfac = rp / x21
        x11  = x11 *zfac
        x21  = x21 *zfac
        x31  = x31 *zfac
        x41  = x41 *zfac

        zfac = rsv / x12
        x12  = rsv
        x22  = x22 * zfac
        x32  = x32 * zfac
        x42  = x42 * zfac
        
        np   = x11*x41 - x21*x31
        nsv  = x12*x42 - x22*x32

        return
        end

        subroutine gttdtdp(omega2,wvno2,omg,wvn,p,ipsvsh,r,
     1      llow,lhgh,time,dtdp,tstar)
        integer ipsvsh, llow, lhgh
        real r, time, tstar
        COMPLEX*16 p
        COMPLEX*16 A,B,C,D,E,F
        COMPLEX*16 wvno2, omega2, wvn, omg
        COMPLEX*16 rsh, rp, rsv
c-----
c       omega2  C - angular frequency squared
c       wvno2   C - wavenumber squared
c       omg     C - angular frequency
c       wvn     C - wavenumber
c       p       C - ray parameter
c       ipsvsh  I - 1 P, 2 SV, 3 SH, 4 pP, 5 sP
c                  since this is for the direct arrival pP and sP not considered
c       r       C - distance
c       llow    I - layer interface indices
c       lhgh    I - layer interface indices
c       time    R - travel time
c       dtdp    C - This must be zero for the direct arrival
c       tstar   R - attenuation operator
c-----
c       get the modified eigen vectors x11 and x31 have 
c           lambda1 (rp ) factored out
c               modified eigen vectors x22 and x42 have 
c           lambda2 (rsv) factored out
c-----
        COMPLEX*16 X11, X21, X31, X41
        COMPLEX*16 X12, X22, X32, X42
        COMPLEX*16 NP, NSV

        COMPLEX*16 detadp
        COMPLEX*16 dtdp

        integer NL
        parameter (NL=200)
        common/timodel/TD(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      TRho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL), 
     3      frefp(NL), frefs(NL)
        real TD, TA, TC, TF, TL, TN, TRho
        real qa, qb, etap, etas, frefp, frefs
        common/modlly/mmax
        integer mmax

        dtdp = dcmplx(dble(r), 0.0d+00)
        time = p*r
        ts = 0.0
        do 1000 l=llow,lhgh

            call tgetegn(a,b,c,d,e,f,omega2,wvno2,rp, rsv, rsh,
     1          x11,x21,x31,x41,x12,x22,x32,x42,NP, NSV,l,omg,wvn)     
        if(ipsvsh.eq.1 .or. ipsvsh.eq.4 .or. ipsvsh.eq.5)then
C               if(dimag(rp).lt. 0.0d+00)then
C                    rp = - rp
C               endif
        dtdp  = dtdp + 
     1      TD(l)*detadp(p,rp *x11,x21,rp *x31,x41,NP ,l,wvn,omg,rp )
        time = time + rp *TD(l)/dcmplx(0.0d+00, 1.0d+00)
C       write(6,*)'l,deta:',l,TD(l),detadp(p,rp *x11,x21,rp *x31,x41,NP ,l,wvn,omg,rp ),dtdp
                    ts = ts + qa(l)*TD(l)*(p*p - rp*rp)
     1                   /(rp/dcmplx(0.0d+00, 1.0d+00))

        else if(ipsvsh.eq.2)then
C               if(dimag(rsv).lt. 0.0d+00)then
C                    rsv = - rsv
C               endif
        dtdp = dtdp + TD(l)*detadp(p,x12,rsv*x22,x32,
     1      rsv*x42,NSV,l,wvn,omg,rsv)
        time = time + rsv*TD(l)/dcmplx(0.0d+00, 1.0d+00)
                    ts = ts + qb(l)*TD(l)*(p*p - rsv*rsv)
     1                   /(rsv/dcmplx(0.0d+00, 1.0d+00))
        else if(ipsvsh.eq.3)then
C               if(dimag(rsh).lt. 0.0d+00)then
C                    rsh = - rsh
C               endif
        dtdp = dtdp + TD(l)*((omg*wvn*TN(l)/TL(l))/rsh)/
     1      cmplx(0.0d+00, 1.0d+00)
        time = time + rsh*TD(l)/dcmplx(0.0d+00, 1.0d+00)
                    ts = ts + qb(l)*TD(l)*(p*p - rsh*rsh)
     1                   /(rsh/dcmplx(0.0d+00, 1.0d+00))
        endif
 1000   continue
        tstar = ts
        return
        end

        function detadp(p,x1,x2,x3,x4,NORM,m,wvn,omg,nu)
c-----
c       van der Hijden  6.109
c-----
        implicit none
        complex*16 p,x1,x2,x3,x4,NORM,nu
        complex*16 detadp
        integer m
        complex*16 wvn, omg
        integer NL
        parameter(NL=200)
        common/timodel/H(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      TRho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL), 
     3      frefp(NL), frefs(NL)
        real H, TA, TC, TF, TL, TN, TRho
        real qa, qb, etap, etas, frefp, frefs
        complex*16 da, db, dc, dd, de, df
        da =   omg *TF(m)/TC(m)
        db =   dcmplx(0.0d+00, 0.0d+00)
        dc =   2.0d+00*wvn*omg*(TA(m) - TF(m)*TF(m)/TC(m)) 
        dd = - omg
        de =   dcmplx(0.0d+00, 0.0d+00)
        df =   dcmplx(0.0d+00, 0.0d+00) 
        detadp =  x4 * (       dd*x2         + de*x4  )
     1      - x3 * (da*x1        + db*x3          )
     1      - x2 * (       df*x2         - dd*x4  )
     1      + x1 * (dc*x1        - da*x3          )
        detadp = detadp /(2.0d+00 * nu * NORM )
        detadp = detadp/dcmplx(00d+00, 1.0d+00)
        return
        end

        function getvel(TA,TL,TN,TRho,m,ipsvsh)
c-----
c     this determines the horizontally propagating velocity
c     This is useful for the refraction and for determining the
c     limits on a reflected arrival
c-----
            real TA(*), TL(*), TN(*), TRho(*)
            integer m, ipsvsh
            real getvel
        
            if(ipsvsh.eq.1)then
                getvel = sqrt(TA(m)/TRho(m))
            else if(ipsvsh.eq.2)then
                getvel = sqrt(TL(m)/TRho(m))
            else if(ipsvsh.eq.3)then
                getvel = sqrt(TN(m)/TRho(m))
            else
                getvel = 1.0
            endif
        return
        end
