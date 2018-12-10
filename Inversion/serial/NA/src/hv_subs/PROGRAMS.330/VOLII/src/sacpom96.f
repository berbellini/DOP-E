        program sacpom96
c----------------------------------------------------------------------c
c                                                                      c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                 c
c      VOLUME II                                                       c
c                                                                      c
c      PROGRAM: SACPOM96                                               c
c                                                                      c
c      COPYRIGHT 1988, 2001                                            c
c      R. B. Herrmann                                                  c
c      Department of Earth and Atmospheric Sciences                    c
c      Saint Louis University                                          c
c      221 North Grand Boulevard                                       c
c      St. Louis, Missouri 63103                                       c
c      U. S. A.                                                        c
c                                                                      c
c----------------------------------------------------------------------c
c       CHANGES
c       19 JUNE 2001 - 
c           - used the permin,permax in USER1 and USER2 to limit
c             the stack frequencies
c       27 JAN 2004 - add flags to increase length of time series for
c           more frequency resolution, 
c           e.g., -2 -4 make 2 x and 4 x longeR
c       PROBLEM AT PRESENT DOES NOT CHECK TO SEE IF 
c           NUMBER OF POINTS AND DT are SAME
c       22 JUL 2004 - line 106 typo pmaxi instead of pmax fixed to
c           pmin,pmax,ido1248
c       23 JUL 2009 - ensured consistent call to FFT routine with
c               complex arguments. We now call zfour
c-----
        integer LIN, LOT, LER
        parameter (LIN=5, LER=0, LOT=6)
        integer NPTSER, NPTFRQ, MAXTRC
        parameter (NPTSER=132000, NPTFRQ=65000)
        parameter (MAXTRC=500)
        complex*8 x(NPTFRQ)
        real*4 dt,df
        real*4 dis(MAXTRC)
        real*4 t0(MAXTRC)
        real*4 vmin,vmax,fmin,fmax
        integer*4 nray
        integer*4 n , k
        integer*4 iunit
        integer*4 lun
        character ntap*50 
        logical ylin, xlin, verbos, query, errbar
        logical shadon
        integer ido1248
        common/pltcnt/laxper
            logical laxper
        common/discnt/kount, lorr
            integer *4 kount, lorr
        common/gtf/doabs,ampmx
        logical doabs
c-----
c       information about plot on the page
c-----
        common/pltinf/gxl(3),gyl(3),gxh(3),gyh(3),
     1      axl(3),ayl(3),axh(3),ayh(3),ixlnlg(3),iylnlg(3)
            real*4 gxl,gyl,gxh,gyh,
     1              axl,ayl,axh,ayh
            integer*4 ixlnlg,iylnlg
        integer MAXFRQ
        parameter (MAXFRQ=NPTFRQ)
        common /fval/ wn, nper, wnin, iwn, mper
            real*4 wn(MAXFRQ), wnin(MAXFRQ)
            integer*4 iwn(MAXFRQ)
            integer*4 nper, mper

        integer MAXRAY
        parameter (MAXRAY=500)

        character*3 lnlg(2)
        character outstr*256
c-----
c       initialize some variables
c-----
        kount = 0
        lnlg(1)='lin'
        lnlg(2)='log'
c-----
c       machine dependent
c-----
        call mchdep()
c-----
c
c       read program paramters from command line 
c
c-----
        call gcmdln(ntap,nray,vmin,vmax,verbos,
     1      xlin,ylin,iunit,query,errbar,shadon,
     2      laxper,pmin,pmax,ido1248)

        if(nray.gt.MAXRAY)nray=MAXRAY
c-----
        if(ntap.ne.' ')then
            lun = 2
            open(lun,file=ntap,status='old',form='formatted',
     1          access='sequential')
            rewind lun
        else
            lun = LIN
        endif
c-----
c       get spectra from spectra files and store on unit 3
c-----
        open(3,access='sequential',form='unformatted',status='scratch')
        rewind 3
        call getspc(x,lun,t0,dis,k,dt,df,n,npts,pmin,pmax,ido1248)
        if(lun.eq.2)close(2)
c-----
c       verbose output here
c-----
        if(verbos)then
            call verbal(ntap,nray,pmin,pmax,vmin,vmax,
     1      xlin,ylin,n,errbar,laxper,shadon,k)
        endif
c-----
c       in sacmft96 USER1 and USER2 are used for 
c       permin and permax limits
c       which are placed in the header by sacmft96 to 
c       indicate the limits
c       of phase match filtering. just use defaults here
c----- 
c-----
c       determine filter periods
c-----
        call getper(pmin,pmax) 
        call mapper(n,df)
        fmin = 1.0 / pmax
        fmax = 1.0 / pmin
c-----
c       now map periods onto the Fourier frequencies
c-----
c-----
c-----
c       process data
c-----

        call pinitf('POM96.PLT')

c-----
c       open file containing dispersion values
c-----
        open(8,file='pom96.dsp',access='sequential',
     1      form='formatted',status='unknown')
        rewind 8
c-----
c       open scratch files
c-----
        open(2,status='scratch',form='unformatted',access='sequential')
        rewind 2
        open(4,status='scratch',form='unformatted',access='sequential')
        rewind 4
c-----
c       begin the phase velocity stacking
c-----
            call pom(x,dt,df,n,npts,nray,k,dis,
     1          fmin,fmax,vmin,vmax,xlin,ylin,errbar,
     2          ntap,t0,iunit,shadon,NFLOW,NFUPR)
        write(LER,*)'Return from pom'
c-----
c       close open files
c-----
        close (2)
        close (4)
        close (8)
        close (3)
        call pend()
c-----
c       open file containing plot control information
c-----
        write(LER,*)'pom96.ctl'
        open(3,file='pom96.ctl',access='sequential',
     1      form='formatted',status='unknown')
        rewind 3
        write(3,*)kount
        if(laxper)then
            write(3,14)
        else
            write(3,13)
        endif
        do 11 i=2,2
        write(3,10) gxl(i),gyl(i),gxh(i),gyh(i),
     1      axl(i),ayl(i),axh(i),ayh(i),
     2      lnlg(ixlnlg(i)),lnlg(iylnlg(i))
   11   continue
        nper = NFUPR - NFLOW + 1
        write(3,*)nper
        do 1234 i=1,nper
            jj = NFLOW + i - 1
            wnin(i) = 1.0/(jj*df)
 1234   continue
        write(3,'(5g15.7)')(wnin(i),i=1,nper)
   10   format(4f7.3,4g11.3,' ',a3,' ',a3,' POM96.PLT')
   13   format('XAXIS-FREQUENCY')
   14   format('XAXIS-PERIOD')
c-----
c       close the plot parameter file
c-----
        close (3)
c-----
c       open a command file for running sdpegn96 for an overlay
c-----
        write(LER,*)'POM96CMP'
        open(3,file='POM96CMP',access='sequential',
     1      form='formatted',status='unknown')
        rewind 3
        write(3,'(a)')'#!/bin/sh'
        write(3,'(a)')' '
        outstr = ' '
   15   format('sdpegn96 -X0 ',f5.2,' -Y0 ',f5.2,' -XLEN ',
     1   f5.2, ' -YLEN ',f5.2,' -XMIN ',g10.3,
     1   ' -XMAX ',g10.3, ' -YMIN ',f5.2,' -YMAX ',f5.2)
        write(outstr,15)gxl(2),gyl(2),gxh(2)-gxl(2),
     1      gyh(2)-gyl(2),axl(2),axh(2), ayl(2), ayh(2)
        ls = lgstr(outstr)
        if(laxper)then
            outstr(ls+1:ls+7) = ' -PER  '
        else
            outstr(ls+1:ls+7) = ' -FREQ '
        endif
        ls = ls + 7
        if(lorr.eq.1)then
            outstr(ls+1:ls+13) = '-L -C -NOBOX -W 0.01 '
        else if(lorr.eq.2)then
            outstr(ls+1:ls+13) = '-R -C -NOBOX -W 0.01 '
        endif
        ls = ls + 13
        if(xlin)then
            outstr(ls+1:ls+6)='-XLIN '
        else
            outstr(ls+1:ls+6)='-XLOG '
        endif
        ls = ls + 6
        if(ylin)then
            outstr(ls+1:ls+6)='-YLIN '
        else
            outstr(ls+1:ls+6)='-YLOG '
        endif
        ls = ls + 6
        write(3,'(a)')outstr(1:ls)
c-----
c       close the SDPEGN96 invocation file
c-----
        close (3)
c-----
c       UNIX
c-----
C       iret = system('chmod +x POM96CMP')
c-----
c       terminate program
c-----
        end

        subroutine usage()
        integer LIN, LOT, LER
        parameter (LIN=5, LER=0, LOT=6)
           write(LER,*)
     :'Usage: sacpom96 -C sacfilelist',
     :' -XLIN -XLOG -PMIN pmin -PMAX pmax',
     :' -VMIN vmin -VMAX vmax',
     :' -nray nray -V -E  -R -L -S',
     :' -FREQ -PER -A -? -h -2 -4 -8'
        write(LER,*)
     :' -XLIN        (default= false)  linear frequency axis'
        write(LER,*)
     :' -XLOG        (default= true)   logarithmic frequency axis'
        write(LER,*)
     :' -S           (default = false) contour shading'
        write(LER,*)
     :' -R       (default = true ) data are Rayleigh'
        write(LER,*)
     :' -L           (default = false) data are Love'
        write(LER,*)
     :' -C  saclist  (default stdin)   input data file list'
        write(LER,*)
     :' -nray  nray  (default = 20)    number of ray parameters '
        write(LER,*)
     :' -PMIN  pmin  (default=4.00)    minimum frequency for plot'
        write(LER,*)
     :' -PMAX  pmax  (default=60.0)    maximum frequency for plot'
        write(LOT,*)
     1' -FREQ        (default false)   x-axis is frequency'
        write(LOT,*)
     1' -PER         (default true)    x-axis is period'
        write(LER,*)
     :' -2       (default = 1) lengthen time zeries by 2x add zero'
        write(LER,*)
     :' -4       (default = 1) lengthen time zeries by 4x add zero'
        write(LER,*)
     :' -8       (default = 1) lengthen time zeries by 8x add zero'
        write(LOT,*)
     :' -VMIN  vmin  (default 2.0)     minimum velocity for plot'
        write(LER,*)
     :' -VMAX  vmax  (default 5.0)     maximum velocity for plot'
        write(LER,*)
     :' -E           (default = false) plot error bars'
        write(LER,*)
     1' -A           (default false)   plot absolute amplitude contours'
        write(LOT,*)
     1' -?           (default false)   usage'
        write(LOT,*)
     1' -h           (default false)   usage'
        stop
        end

        subroutine gcmdln(ntap,
     1      nray,vmin,vmax,verbos,xlin,ylin,
     2      iunit,query,errbar,shadon,
     2      laxper,pmin,pmax,ido1248)
c-----
c     parse command line arguments and return control
c     parameters
c
c     requires subroutine mgtarg(i,name) to return
c           the i'th argument in the string name
c
c     and the function mnmarg() to return the number
c           of arguments excluding program name
c           The first argument is i = 1
c
c-----
c       get command arguments
c
c       ntap    C*120   - input file name
c       nray    I*4 - number of ray parameters
c       fmin    R*4 - minimum frequency for plot
c       fmax    R*4 - maximum frequency for plot
c       vmin    R*4 - minimum velocity for plot
c       vmax    R*4 - maximum velocity for plot
c       verbos  L   - verbose output on standard error
c       xlin    L   - Frequency linear plot
c       ylin    L   - Velocity linear plot
c       iunit   I*4 - 0 for KM, 1 for FT, 2 for M
c       query   L   - if -? flag is set call help, exit
c       errbar  L   - plot velocity error bars if true
c       shadon  L   - true contour color shading
c       laxper  L   - Period if true, frequency if false
c       pmin    R*4 - minimum period
c       pmax    R*4 - maximum period
c       ido1248 I*4 - add zeros to make time series 1 2 4 x longer
c-----
        character ntap*(*)
        integer*4 nray
        real*4 ci, ce,fmin,fmax,vmin,vmax
        real*4 pmin, pmax
        logical verbos,xlin,ylin,query,errbar
        logical shadon, laxper
        integer ido1248
        character*50 name
        integer*4 mnmarg

        common/gtf/doabs,ampmx
        logical doabs

        common/discnt/kount, lorr
        integer *4 kount, lorr

            ntap = ' '
            ci = -1.0
            ce = -1.0
            nray = 20
            fmin = -1.0
            fmax = -1.0
            vmin = -1.0
            vmax = -1.0
            lw = 500
            verbos = .false.
            xlin  = .false.
            ylin  = .true.
            errbar = .false.
            query = .false.
            shadon = .false.
            iunit = 0
            laxper = .true.
            doabs = .false.
            lorr = 0
            pmin = 4.0
            pmax = 60.0
            ido1248 = 1
            nmarg = mnmarg()
            i = 0
   11       i = i + 1
            if(i.gt.nmarg)goto 13
            call mgtarg(i,name)
            if(name(1:2).eq.'-?' .or.name(1:2).eq.'-h')then
                call usage()
            elseif(name(1:5).eq.'-XLIN')then
                xlin = .true.
            elseif(name(1:5).eq.'-XLOG')then
                xlin = .false.
            else if(name(1:2).eq.'-R')then
                lorr = 2
            else if(name(1:2).eq.'-L')then
                lorr = 1
            elseif(name(1:2).eq.'-E')then
                errbar = .true.
            elseif(name(1:2).eq.'-S')then
                shadon = .true.
            elseif(name(1:5).eq.'-nray')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(i10)')nray
            elseif(name(1:5).eq.'-PMIN')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(f10.0)')pmin
            elseif(name(1:5).eq.'-PMAX')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(f10.0)')pmax
            elseif(name(1:5).eq.'-VMIN')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(f10.0)')vmin
            elseif(name(1:5).eq.'-VMAX')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(f10.0)')vmax
            elseif(name(1:2).eq.'-C')then
                i = i + 1
                call mgtarg(i,ntap)
            elseif(name(1:2).eq.'-V')then
                verbos = .true.
            else if(name(1:5).eq.'-FREQ')then
                laxper = .false.
            else if(name(1:4).eq.'-PER')then
                laxper = .true.
            else if(name(1:2).eq.'-A')then
                doabs = .true.
            else if(name(1:2).eq.'-2')then
                ido1248 = 2
            else if(name(1:2).eq.'-4')then
                ido1248 = 4
            else if(name(1:2).eq.'-8')then
                ido1248 = 8
            endif
            goto 11
   13   continue
c-----
c       test for improper command line
c-----
        if(vmin.lt.0.0 .or.vmax.lt.0.0)then
            vmin = 2.0
            vmax = 5.0
        endif
        if(nray.le.1)nray=20
        if(pmin.lt.0.0 .or. pmax.lt.0.0)then
                   pmin = 4.00
                   pmax = 60.0
        endif
        if(pmin .eq. 0.0)pmin = 0.01 * pmax
        return
        end

        subroutine verbal(ntap,nray,pmin,pmax,vmin,vmax,
     1      xlin,ylin,npts,errbar,laxper,shadon,k)
c-----
c       get command arguments
c
c       ntap    C*120   - input file name
c       nray    I*4 - number of ray parameters
c       pmin    R*4 - minimum period for plot
c       pmax    R*4 - maximum period for plot
c       vmin    R*4 - minimum velocity for plot
c       vmax    R*4 - maximum velocity for plot
c       xlin    L   - Frequency linear plot
c       ylin    L   - Velocity linear plot
c       npts    I*4 - number of points in FFT
c       errbar  L   - plot error bars
c       laxper  L   - x-axis is period else freq
c       shadon  L   - true contour color shading
c       k   I*4 - number of traces
c-----
        integer LIN, LOT, LER
        parameter (LIN=5, LER=0, LOT=6)
        character ntap*(*)
        integer*4 nray, k
        real*4 pmin,pmax,vmin,vmax
        logical xlin,ylin,errbar,laxper,shadon

        common/discnt/kount, lorr
            integer *4 kount, lorr

        ls = lgstr(ntap)

        write(LER,*)' Values read from file command line'
        write(LER,*)' Input  data file                    = ',ntap(1:ls)
        write(LER,*)' Number of ray parameters            = ',nray
        write(LER,*)' Minimum period    for plots         = ',pmin
        write(LER,*)' Maximum period    for plots         = ',pmax
        write(LER,*)' Minimum velocity for plots          = ',vmin
        write(LER,*)' Maximum velocity for plots          = ',vmax
        write(LER,*)' Linear frequency contour plot       = ',xlin
        write(LER,*)' Linear velocity  contour plot       = ',ylin
        write(LER,*)' Error bars plotted                  = ',errbar
        write(LER,*)' X-axis is period                    = ',laxper
        write(LER,*)' FFT NPTS                            =  ',NPTS
        if(lorr.eq.0)then
            write(LER,*)' Unknown wave data'
        else if(lorr.eq.1)then
            write(LER,*)' Rayleigh wave data'
        else if(lorr.eq.2)then
            write(LER,*)' Love wave data'
        endif
        write(LER,*)' Contour shading                     = ',shadon
        write(LER,*)' Number of traces                    = ',k
        return
        end

        subroutine getspc(x,lun,t0,dis,k,dt,df,n,npts,pmin,pmax,ido1248)
c-----
c       read file names from lun, open each, read in spectra
c       and write sequentially to unit 3
c-----
c       x   C*8 complex spectra array
c       lun I*4 logical unit for file names
c       t0  R*4 array of initial sample times
c       dis R*4 array of distances
c       k   I*4 total number of spectra
c       dt  R*4 sample interval in seconds
c       df  R*4 frequency spacing
c       n   I*4 FFT points
c       npts    I*4 Points in time series
c       pmin    R*4 Minimum period for plot, but will be adjusted upward
c               according to any limits in the raw spectra
c       pmax    R*4 Maximum period for plot, but will be adjusted down
c               accroding to maximum period in data set
c       ido1248 I*4 Increase time seeries by adding zeros to make
c                   1 2 or 4x longer
c-----
        integer NPTSER, NPTFRQ
        parameter (NPTSER=132000, NPTFRQ=65000)
        complex x(*)
        integer*4 lun,n
        real*4 t0(*), dis(*),dt,df
        character name*256
      real*4 dist, deg, az, baz, tt0
      integer*4  n21 , npts
      character sta *8, comp *8, cdate *12
        real tarr(NPTSER)
        integer nend
        k = 0
 1000   continue
            read(lun,'(a)',end=2000,err=2000)name
c-----
c           process data
c-----
            call getsac(name,npts,dist,deg,az,baz,tt0,dt,sta,
     1                 comp,cdate,tarr,permin,permax)
            if(npts.gt.nptser)npts = nptser
c-----
c       safety for extending time series
c-----
            if(ido1248.gt.1)then
                if(npts .lt. nptser/2 .and. ido1248.eq.2)then
                    nend = 2*npts
                endif
                if(npts .lt. nptser/4 .and. ido1248.eq.4)then
                    nend = 4*npts
                endif
                if(npts .lt. nptser/8 .and. ido1248.eq.8)then
                    nend = 8*npts
                endif
                do 1002 i=npts+1,nend
                    tarr(i) = tarr(npts)
 1002           continue
                npts = nend 
            endif
c-----
c       end of extending time series
c-----

            if(permin.gt.0.0 .and. permin.gt.pmin)pmin = permin
            if(permax.gt.0.0 .and. permax.lt.pmax)pmax = permax
            call npow2(npts,n,n21)
c-----
c       set up time series for spectra computation --
c       be sure to correct for geometrical spreading
c       by multiplying by the sqrt(dis)
c-----
            do 100 i=1,n
                if(i.le.npts)then
                    x(i) = cmplx(tarr(i),0.0)*
     1                  sqrt(dist)
                else
                    x(i) = cmplx(0.0,0.0)
                endif
  100       continue
            call zfour(x,n,-1,dt,df)
            k = k + 1
            t0(k) = tt0
            dis(k) = dist
        ls = lgstr(name)
        write(6,*)'Reading File',k,tt0,dist,' ',name(1:ls)
        do 9875 il=1,n21,512
            iu = (il - 1) + 512
            if(iu .gt. n21)iu = n21
            write(3)(x(i),i=il,iu)
 9875   continue
        go to 1000
 2000   continue
        return
        end


        subroutine pom(x,dt,df,npts,ndpts,nray,ntrc,dis,
     1      fmin,fmax,vmin,vmax,xlin,ylin,errbar,
     2      ntap,t0,iunit,shadon,NFLOW,NFUPR)
c-----
c       subroutine to perform p-omega stack
c       the procedure is to evaluate
c       
c       SUM cos(phi(i) + omega p r(i))
c
c       This is consistent with the definition of the FFT used
c
c       E.g. G(f) = SUM g(t)exp(-i 2 pi f t) dt
c
c-----
c*****
c       WE DO NOT CONSIDER THE TWO LOWEST FFT FREQUENCIES
c       0.0 dt and 1.0df
c*****
c       SUBROUTINE ARGUMENTS
c
c       x   C*8 Two dimensional array of spectra of trace
c               traces(LL,LR) for trace LR and freq LL
c       dt  R*4 Sample interval in seconds
c       df  R*4 Frequency sampling interval in Hz
c       npts    I*4 Number of points for FFT
c       ndpts   I*4 Number of samples in original time series
c       nray    I*4 number of ray parameters
c       ntrc    I*4 Number of traces
c       dis R*4 Array of source-receiver distances
c       fmin    R*4 Minimum frequency for plot
c       fmax    R*4 maximum frequency for plot
c       vmin    R*4 Minimum velocity  for plot
c       vmax    R*4 Maximum velocity  for plot
c       xlin    L   Linear frequency axis if true
c       ylin    L   Linear velocity  axis if true
c       errbar  L   Plot error bars on dispersion plot
c       ntap    C*120   File name
c       nrec    I*4 Record number
c       ktrc1   I*4 Starting trace number
c       ktrc2   I*4 Ending trace number
c       nsi I*4 Sample interval in milliseconds
c       NSAMP   I*4 Number of samples in original time series
c       t0  R*4 Array of times of first sample
c       iunit   I*4 Unit Flag 0=ft, 1=m
c       laxper  L   X-axis is period if true
c       shadon  L   - true contour color shading
c-----
c
c-----
        integer LIN, LOT, LER
        parameter (LIN=5, LER=0, LOT=6)
        integer NPTSER, NPTFRQ, MAXTRC
        parameter (NPTSER=132000, NPTFRQ=65000)
        parameter (MAXTRC=500)
        integer BLUE, BLACK
        parameter (BLUE=4, BLACK=1)
        complex*8 x(*)
        character ntap*(*)
        complex*8 cfac1
        complex*8 cfac(NPTFRQ)
        complex*8 data(NPTFRQ)
        real*4 tri(NPTFRQ)
        real*4 dis(MAXTRC)
        real*4 t0(MAXTRC)
        real*4 amprel
        logical xlin, ylin, errbar, shadon
        equivalence(data(1),tri(1))

        common/pltcnt/laxper
              logical laxper
        common/gtf/doabs,ampmx
        logical doabs
        
        integer MAXAMP
        parameter (MAXAMP=10)
        common/phase/phvel(MAXAMP,NPTFRQ), phamp(MAXAMP,NPTFRQ),
     1      phfrq(MAXAMP,NPTFRQ), numph(NPTFRQ)
            real phvel, phamp, phfrq
            integer numph

c-----
c       array for contour plot
c       real*4  cont(NX,NY) the array to be plotted
c       real*4  x(NX)   x-axis values for contour
c               must be in increasing order
c       real*4  y(NY)   y-axis values for contour
c               must be in increasing order
c-----
        integer NX, NY
        parameter(NX=100,NY=100)
        real cont(NX,NY)
        real*4 wn(NX)



              
c-----
c       initialize array for ploting phase velocity and
c       spectral amplitudes
c-----
        pvmin = vmin
        pvmax = vmax
        DO 5000 I = 1,NPTFRQ
            numph(I) = 0
            DO 5001 J = 1,MAXAMP
                phvel(J,I) = 0.0
                phamp(J,I) = 0.0
                phfrq(J,I) = 0.0
 5001       continue
 5000       continue
c-----
c       initialize array for phase velocity contours
c-----
        do 6001 i=1,NX
            do 6002 j=1,NY
                cont(i,j) = 0.0
 6002       continue
 6001   continue


c-----
c       get number of positive frequency points
c       check maximum plot frequency so it does not exceed nyquist
c-----
        npts21 = npts/2 + 1
        fnyq = 0.5/dt
c-----
c       SPECIFY LIMITS OF FREQUENCY SAMPLING
c       ALSO DO NOT USE FREQ = 0.0 DF and 1.0 DF
c-----
        tmpflw = 2.0/(ndpts*dt)
        if(FMIN .lt. tmpflw)then
            NFLOW = tmpflw/df  + 1
        else
            NFLOW = FMIN/df  + 1
        endif
        
c-----
c       ensure that we have two cycles per original data length
c-----
        NFUPR = FMAX/df +0.5 + 1
        if(NFLOW.LT.3)NFLOW = 3
        if(NFUPR.gt.npts21)NFUPR = npts21
        if(fmax.gt.fnyq)fmax = fnyq
        wmin = (NFLOW-1)*df
        wmax = (NFUPR-1)*df
c-----
c       do p-omega stack
c-----
        ampmx = 0.0
        dc = (vmax -vmin)/(nray -1)
        DO 400 IR = 1,nray
            rewind 3
            c = vmin + (IR-1)*dc
            p = 1./c
            amprel = -1.0e+38
c-----
c       DO STACK FOR EACH FREQUENCY, BUT ANALYTICALLY DIVIDE OUT
c       THE RESPONSE OF THE FIRST TRACE TO REMOVE ANY SOURCE,
c       DISTANT INDEPENDENT PHASE TERM, e.g., Bessel Function
c       If first trace is of low amplitude at any  frequency,
c       do not perform the correction. This may just be a rare
c       numerical event.
c
c       However, if other traces exhibit a very low amplitude,
c       do not include them in the stack [ if(ca .gt.0.0) ] below
c
c-----
            if(mod(IR,5).eq.0)then
            write(6,*)'Processing phase velocity',IR,
     1          ' of ',nray,' c=',c
            endif
            do 200 LR=1,ntrc
c-----
c       SUN compiler bug
c               read(3)(x(i),i=1,npts21)
c-----
                do 9875 il=1,npts21,512
                    iu = (il - 1) + 512
                    if(iu .gt. npts21)iu = npts21
                    read(3)(x(i),i=il,iu)
 9875           continue
                amprel = 0.0
                r = dis(LR)
                do  100 LL=NFLOW,NFUPR
                    omega = 6.2831853*(LL-1)*df
                    omegap= omega*p
                    ompr=omegap*dis(LR)-omega*t0(LR)
                    ca = cabs(x(LL))
c------
c       this actually corrects for geometrical spreading since we only
c       look at the phase term, we form  H(f) / |H(f)|
c-----
                    if(ca.gt.0.0)then
                        cfac1 = x(LL)/ca
                    else
                        cfac1 = cmplx(0.0,0.0)
                    endif
                    cfac1 = cfac1*cmplx(cos(ompr),sin(ompr))
                    if(LR.eq.1)then
                        data(LL) = cmplx(0.0,0.0)
                        cfac(LL)=conjg(cfac1)
                    endif
                    data(LL)=data(LL)+cfac1*cfac(LL)
  100           continue
  200       continue
c-----
c       get maximum values at this ray parameter
c-----
            do 300 LL=NFLOW,NFUPR
                tri(LL) = cabs(data(LL))
                if(tri(LL).gt.amprel)amprel = tri(LL)
  300       continue
            if(amprel.gt.ampmx)ampmx = amprel
            if(amprel.eq.0.0)amprel = 1.0
c-----
c           here map on to the cont array
c-----
            call mappv(IR,nray,tri,npts21,NFLOW,NFUPR,df,cont,NX,NY,
     1          pvmin,pvmax,c,wn)
  400       continue
            if(ampmx.eq.0.0)ampmx=1.0
c-----
c       plot dispersion curve
c------
            if(nray.gt.NY)then
                nnray = NY
            else
                nnray = nray
            endif   
            if(doabs)then
c-----
c           normalize
c-----
                do 4110 j=1,nnray
                    do 4111 i=1,NX
                        cont(i,j) = cont(i,j)/ampmx
                        cont(i,j) = cont(i,j)**2
 4111               continue
 4110           continue
            else
c-----
c           do relative plot
c-----
C               do 4112 i=1,NX
C                   amprel = 0.0
C                   do 4113 j=1,nnray
C                       if(cont(i,j).gt.amprel)then
C                           amprel = cont(i,j)
C                       endif
C 4113              continue
C                   if(amprel.eq.0.0)amprel = 1.0
C                   do 4114 j=1,nray
C                       cont(i,j) = cont(i,j)/amprel
C 4114              continue
C 4112          continue
                do 4112 j=1,nray
                    amprel = 0.0
                    do 4113 i=1,NX
                        if(cont(i,j).gt.amprel)then
                            amprel = cont(i,j)
                        endif
 4113               continue
                    if(amprel.eq.0.0)amprel = 1.0
                    do 4114 i=1,NX
                        cont(i,j) = cont(i,j)/amprel
 4114               continue
 4112           continue
            endif
c-----
c       define the peak stack values for each frequency
c-----
            call getpv(nray,tri,df,NFLOW,NFUPR,pvmin,pvmax,dis,ntrc)
c-----
c       perform contour plots
c-----
        write(6,*)'Call conplt'
            call conplt(n,wn,df,t0,dt,dist,nray,
     1          wmin,wmax,
     2          pvmin,pvmax,xlin,ylin,shadon,fmin,fmax,
     3          cont,pmin,pmax)
c-----
c       clean up open files
c-----
        close (1)
        call plot(-2.0,-1.0,-3)
        return
        end

        subroutine conplt(n,wn,df,t0,dt,dist,nray,
     1       wmin,wmax,
     2  pvmin,pvmax,xlin,ylin,shadon,fmin,fmax,
     3          cont,pmin,pmax)
c-----
c       nray    R*4 number of actual phase velocities
c       cont    R*4 NX x NY array to be contoured
c       NX  I*4 array dimension
c       NY  I*4 array dimension
c       wn  R*4 array of frequencies
c-----
        integer NPTSER, NPTFRQ
        parameter (NPTSER=132000, NPTFRQ=65000)
        integer BLACK, BLUE
        parameter (BLACK=1, BLUE=4)
              integer*4 n
              real*4 df
              real*4 t0,dt,dist
              logical ylin,xlin
        common/pltcnt/laxper
              logical laxper
        integer MAXAMP
        parameter (MAXAMP=10)
        common/phase/phvel(MAXAMP,NPTFRQ), phamp(MAXAMP,NPTFRQ),
     1      phfrq(MAXAMP,NPTFRQ), numph(NPTFRQ)
            real phvel, phamp, phfrq
            integer numph
        logical shadon
c-----
c       information about plot on the page
c-----
        common/pltinf/gxl(3),gyl(3),gxh(3),gyh(3),
     1      axl(3),ayl(3),axh(3),ayh(3),ixlnlg(3),iylnlg(3)
            real*4 gxl,gyl,gxh,gyh,
     1              axl,ayl,axh,ayh
            integer*4 ixlnlg,iylnlg
c-----
c       
c-----
        integer NX, NY
        parameter(NX=100,NY=100)
        real cont(NX,NY), x(NX), y(NY)
        real*4 wn(NX)
        parameter (NC=11)
        integer icol(NC+1)
        real CN(NC)

        xoff = 2.0
        yoff = 1.0
        call plot(xoff,yoff,-3)
c-----
c       normalize each trace individually
c-----
c-----
c       set up plot parameters
c-----
c       if xlin = .true. linear frequency
c             else      logarithmic
c       if ylin = .true. linear phase velocity
c             else logarithmic
c-----
        xaxlen = 6.0
        yaxlen = 6.0

        if(nray.gt.NY)then
            nny = NY
        else
            nny = nray
        endif
        
        call setxy(x,y,nny,wn,fmin,fmax,wmin,wmax,pvmin,pvmax,
     1      xlin,ylin,xlow,xhgh,ylow,yhgh,xaxlen,yaxlen)
c-----
c
c-----
        call plot(0.0,0.0,3)
        dx = 1.0
        dy = 1.0
c-----
c       define color array
c-----
        vmx = 0.95
        vmn = 0.0
        do 2000 i = 1, NC
            fac= 1.0*vmn  + real(i-1)/(real(NC-1))*(vmx - vmn) 
            CN(I) = fac
 2000   continue
        do 126 I=1,NC+1
                  ICOL(I)=1100 - 100.0*real(I)/real(NC)
                  if(ICOL(I).gt.1100)ICOL(I) = 1100
                  if(ICOL(I).lt.1000)ICOL(I) = 1000
  126   continue                    
c-----
c       DRAW CONTOURS WITH SHADING (very simply no triangles exactly )
c-----
c       note X and Y must be ascending arrays
c       so that if frequency is plotted do it on the mappv
c       MODE = 0 shade, 1 contours, 2 contours plus shade
c-----
        if(shadon) then
            call FARB2D(X,NX,Y,nny,cont,NX,CN,ICOL,NC,0)  
        else
            call FARB2D(X,NX,Y,nny,cont,NX,CN,ICOL,NC,1) 
        endif
        call newpen(1)
        call gbox(0.0,0.0,xaxlen,yaxlen)
C       write(6,*)'ICOL:',ICOL
C       write(6,*)'CN  :',CN
C       write(6,*)'X   :', (x(i),i=1,NX)
C       write(6,*)'Y   :', (y(i),i=1,NY)
c-----
c       PLOT THE  PHASE VELOCITIES OF THE FOUR LARGEST SPECTRAL
c       AMPLITUDES ON TOP OF THE CONTOUR PLOT. WE USE THE
c       ARRAY SET UP BY mxval AND STORED IN common/phase/phval
c-----
        call newpen(BLACK)
        call pltpv(xlin,ylin,xhgh,xlow,yhgh,ylow,xaxlen,yaxlen,laxper)
        call plot(-xoff,-yoff,-3)
        gxl(2) = xoff
        gyl(2) = yoff
        gxh(2) = xoff + xaxlen
        gyh(2) = yoff + yaxlen
        axl(2) = xlow
        axh(2) = xhgh
        ayl(2) = ylow
        ayh(2) = yhgh
        if(xlin)then
            ixlnlg(2) = 1
        else
            ixlnlg(2) = 2
        endif
        if(ylin)then
            iylnlg(2) = 1
        else
            iylnlg(2) = 2
        endif
        return
        end

        subroutine setxy(x,y,nray,twn,tfmn,tfmx,twmn,twmx,pvmin,pvmax,
     1      xlin,ylin,xlow,xhgh,ylow,yhgh,xaxlen,yaxlen)
c-----
c       set x,y absolute positions for cont(NX,NY) grid points
c-----
c       x(NX)   R*4 - x-positions which will be frequency/period
c       y(NY)   R*4 - y-positions which will be velocity
c       nny I*4 
c       twn R*4 - array of frequencies to the plotted
c       fmin    R*4 - minimum frequency for x-axis
c       fmin    R*4 - minimum frequency for x-axis
c       wmin    R*4 - actual minimum frequency for x-axis
c       wmax    R*4 - actual maximum frequency for x-axis
c       xlin    L   - .true.  x-axis is linear
c                 .false. x-axis is logarithmic
c       ylin    L   - .true.  y-axis is linear
c                 .false. y-axis is logarithmic
c       wmin    R*4 - minimum value of plotted x-values
c       wmax    R*4 - maximum value of plotted x-values
c       pvmin   R*4 -
c       pvmax   R*4 -
c       xlin    L   - .true. if X-axis linear
c       ylin    L   - .true. if Y-axis linear
c       xlow    R*4 -
c       xhgh    R*4 -
c       ylow    R*4 -
c       yhgh    R*4 -
c       xaxlen  R*4 - length of X-axis
c       yaxlen  R*4 - length of Y-axis
c-----
              logical ylin,xlin
c-----
c       laxper  L   - .true. x-axis is period
c                 .false. x-axis is frequency
c-----
        common/pltcnt/laxper
              logical laxper
        parameter (MAXRAY=500)
c-----
c       information about plot on the page
c-----
        common/pltinf/gxl(3),gyl(3),gxh(3),gyh(3),
     1      axl(3),ayl(3),axh(3),ayh(3),ixlnlg(3),iylnlg(3)
            real*4 gxl,gyl,gxh,gyh,
     1              axl,ayl,axh,ayh
            integer*4 ixlnlg,iylnlg
c-----
c       
c-----
        integer NX, NY
        parameter(NX=100,NY=100)
        real*4 x(NX), y(NY)
        character titlex*80, titley*80
        real*4 wn(NX)
        real*4 twn(NX)

        if(nray.gt.NY)then
            nny = NY
        else
            nny = nray
        endif
c-----
c       set up plot parameters
c-----
c       define the array of increasing values
c-----
            do 900 i=1,NX
                if(laxper)then
                    wn(NX+1-i) = 1.0/twn(i)
                else
                    wn(i) =     twn(i)
                endif
  900       continue
        titley = 'Phase Velocity (km/s)'
        if(laxper)then
            fmin = 1.0/tfmx
            fmax = 1.0/tfmn
            wmin = 1.0/twmx
            wmax = 1.0/twmn
            titlex = 'Period (sec)'
        else
            fmin = tfmn
            fmax = tfmx
            wmin = twmn
            wmax = twmx
            titlex = 'Frequency (Hz)'
        endif
        if(xlin)then
            if(fmax.lt.wmax)fmax = wmax
            if(fmin .lt.0.0 .or. fmin.gt.wmin)then
                x1 = wmin
            else
                x1 = fmin
            endif
            xlow = x1
            xhgh = wmax
            do 1000 i=1,NX
                x(i) = 0.0 + xaxlen*(wn(i)-wn(1))/(wn(NX)-wn(1))
 1000       continue
        else
            xlow = wmin
            xhgh = wmax
            do 1001 i=1,NX
                x(i) = xaxlen*alog10(wn(i)/wn(1))/
     1              alog10(wn(NX)/wn(1))
 1001       continue
        endif
c-----
c       set up y-axis
c-----
        if(ylin)then
            ylow = pvmin
            yhgh = pvmax
            do 1003 i=1,nny
                y(i) =0.0 +  (i-1)*yaxlen/(nny-1)
 1003       continue
        else
            ylow = pvmin
            yhgh = pvmax
            dvel = (pvmax - pvmin)/(nny -1 )
            do 1004 i=1,nny
                vel =  pvmin + (i-1)*dvel
                y(i) =  yaxlen*alog10(vel/pvmin)/
     1                  alog10(pvmax/pvmin)
 1004       continue
        endif
c-----
c       plot the axes
c-----
        lx = lgstr(titlex)
        ly = lgstr(titley)
        call       xyaxes(xlow,xhgh,xlin,ylow,yhgh,ylin,
     1          xaxlen,yaxlen,titlex,titley,lx,ly)
        return
        end

        subroutine xyaxes(xlow,xhgh,xlin,ylow,yhgh,ylin,
     1          xaxlen,yaxlen,titlex, titley,lx,ly)
c-----
c       general routine to put up x,y axes
c-----
c       xlow    R*4 minimum value 
c       xhgh    R*4 maximum value 
c       xlin    L   .true. linear x-axis else logarithmic
c       ylow    R*4 minimum value 
c       yhgh    R*4 maximum value 
c       ylin    L   .true. linear y-axis else logarithmic
c       xaxlen  R*4 length of x-axis
c       yaxlen  R*4 length of y-axis
c       titlex  Ch  x-axis title string
c       titley  Ch  y-axis title string
c-----
        real*4 xlow, xhgh, ylow, yhgh
        logical xlin, ylin
        character titlex*(*), titley*(*)

c-----
c       plot x-axis
c-----
        if(xlin)then
            call dolinx(0.0 ,0.0       ,xaxlen,xhgh,xlow,
     1                  0.10,.false.,.false.,.true.,lx,titlex)
        else
            call dologx(0.0 ,0.0       ,xaxlen,xhgh,xlow,
     1                  0.10,.false.,.false.,.true.,lx,titlex)
        endif
        if(ylin)then
            call doliny(0.0       ,0.0,yaxlen,yhgh,ylow,
     1                 0.07,.true. ,.true.,.true. ,ly,titley)
        else
            call dology(0.0       ,0.0,yaxlen,yhgh,ylow,
     1                 0.14,.true. ,.true.,.true. ,ly,titley)
        endif
        return
        end

        subroutine mappv(iray,nray,x,npts,NFLOW,NFUPR,df,cont,NX,NY,
     1          pvmin,pvmax,c,wn)
c-----
c       map the phase velocities uniformly onto the cont(i,j) array
c       where i=velocity and j = frequency = we map with 1,1 at
c       highest frequency and highest velocity
c-----
c       iray    I*4 - phase velocity index
c       nray    I*4 - number of unique ray
c       x   R*4 - array of stack amplitudes
c       npts    I*4 - number of frequencies >= 0
c       cont    R*4 - NX,NY array of values to be contoured
c       NX  I*4 - array dimension
c       NY  I*4 - array dimension
c       nsamp   I*4 - number of samples in original time series <=n
c       pvmin   R*4 - minimum velocity for plot
c       pvmax   R*4 - maximum velocity for plot
c       laxper  L   - .true. x-axis is period
c                 .false. x-axis is frequency
c       c   R*4 - phase velocity corresponding to x(freq) stack
c       wn  R*4 - array of mapped frequencies
c-----
        real x(npts)
        real cont(NX,NY)
        real wn(NX)

        integer LIN, LOT, LER
        parameter (LIN=5, LER=0, LOT=6)

        common/pltcnt/laxper
            logical laxper

c-----
c           Map velocity axis
c-----
        if(nray.gt.NY)then
            p = (c - pvmin)/(pvmax - pvmin)
            i = (1-p)*real(1) + (p)*real(NY) +0.49
            if(i.lt.0)i=1
            if(i.gt.NY)i=NY
        else
            i=iray
        endif
c-----
c       j = 1 => NFLOW  with frequency (NFLOW-1)*df
c       j = NX => NFUPR  with frequency (NFUPR-1)*df
c-----
        flow = (NFLOW-1)*df
        fupr = (NFUPR-1)*df
        do 1000 j=1,NX
c-----
c           Map frequency axis
c-----
            p = real (j-1)/real(NX-1)
            freq = (1-p)*flow + p*fupr
            wn(j) = freq
            jj = (1-p)*real(NFLOW) + p*real(NFUPR) + 0.49
            if(jj.lt.NFLOW)jj =  NFLOW
            if(jj.gt.NFUPR)jj =  NFUPR
c-----
c           laxper == .true for PERIOD    cont(NX-j+1,NY-i+1) = x(jj)
c           laxper != .true for FREQUENCY cont(   j  ,NY-i+1) = x(jj)
c-----
            if(laxper)then
C               cont(NX-j+1,NY-i+1) = x(jj)
                cont(NX-j+1,i) = x(jj)
            else
C               cont(j     ,NY-i+1) = x(jj)
                cont(j     ,i) = x(jj)
            endif
 1000   continue
        write(2)c,(x(i),i=NFLOW,NFUPR)
        return
        end

        subroutine pltpv(xlin,ylin,xhgh,xlow,yhgh,ylow,
     1      xaxlen,yaxlen,laxper)
c-----
c       map phase velocity value into plot coordinate
c       note we already know the x coordinate because of
c       the filter frequencies
c-----  
c       xlin    L   .true. x-axis is linear
c       ylin    L   .true. y-axis is linear
c       xlow    R*4 lowest value of x-axis
c       xhgh    R*4 highest value of x-axis
c       ylow    R*4 lowest value of y-axis
c       yhgh    R*4 highest value of y-axis
c       yaxlen  R*4 length of y-axis
c       laxper  L   .true. x-axis is period
c-----
        logical xlin, ylin
        real*4 ylow,yhgh
        real*4 yaxlen
        logical laxper

        integer NPTSER, NPTFRQ, MAXTRC
        parameter (NPTSER=132000, NPTFRQ=65000)
        parameter (MAXTRC=500)
        parameter (MAXAMP=10)
        common/phase/phvel(MAXAMP,NPTFRQ), phamp(MAXAMP,NPTFRQ),
     1      phfrq(MAXAMP,NPTFRQ), numph(NPTFRQ)
            real phvel, phamp, phfrq
            integer numph

        do 1000 j=1,NPTFRQ
        if(numph(j).gt.0)then
            do 1100 i=1,numph(j)
                if(phamp(i,j).gt.0.0)then
                u = phvel(i,j)
                f = phfrq(i,j)
                if(laxper)then
                    f = 1.0/f
                endif
                if(xlin)then
                    xx = 0.0 + xaxlen*(f-xlow)/(xhgh-xlow)
                else
                    xx = 0.0 + xaxlen*alog10(f/xlow)
     1                  /alog10(xhgh/xlow)
                endif
                if(ylin)then
                    yy = 0.0 + yaxlen*(u-ylow)/(yhgh-ylow)
                else
                    yy = 0.0 + yaxlen*alog10(u/ylow)
     1                  /alog10(yhgh/ylow)
                endif
                call gsolid(xx,yy,0.03,i-1)
                endif
 1100       continue
            endif
 1000   continue
        return
        end
        
        subroutine geterr(freq,vel,dvel,sum,ndist,dist)
c-----
c       estimate the error in phase velocity by using stacking number
c
c       we assume a normal probability density function
c       for the ray parameter error
c
c       freq    - frequency of interest
c       vel - velocity determined
c       dvel    - desired phase velocity error
c       sum - stack summation
c       ndist   - number of distances
c       dist    - array of distances
c-----
        real*4 dist(*)
c-----
c       first use Newton Raphson to determine error, starting with
c       a zero error estimate
c-----
        sig = 0.0
        do 100 i=1,10
            call fnctn(f,df,ddf,sum,ndist,dist,sig)
            if(df.ne.0.0)then
                sig = sig - f/df
            else
                sig = sig + sqrt(abs(2.*f/ddf))
            endif
  100   continue
c-----
c       now convert error in DELTA(freq ray-parameter) to
c       error in velocity
c-----
        dvel = sig * vel*vel/freq
        if(dvel.gt.100.0 .or. dvel.le.0.0)dvel=99.9999
        return
        end

        subroutine fnctn(f,df,ddf,sum,ndist,dist,sig)
        real*4 dist(*)
        f = 0.0
        df = 0.0
        ddf = 0.0
        sig2 = sig*sig
        do 100 i=1,ndist
            fac = - 2. * (3.1415927*(dist(i)-dist(1)))**2
            fac1 = exp(fac*sig2)
            f = f + fac1
            df = df + 2.0*fac*sig*fac1
            ddf= ddf+ fac1*(2.0*fac + (2.0*fac*sig)**2)
  100   continue
        f = f - sum
        return
        end

        subroutine srt(a,u,iu,ikey,ampl,vel,ii,k,MAX4)
c-----
c       subroutine to maintain a list of the MAX4
c       largest values of ampl in a(i)
c       also saving the corresponding vel in u(i), the 
c       mapping into the original array in ikey
c
c       This is being done to avoid have two very large
c       arrays in overhead and because we do not need
c       to do a full sort
c-----
        real*4 a(MAX4), u(MAX4), vel, ampl
        integer*4 iu(MAX4),ikey(MAX4)
        integer*4 MAX4, k
        integer*4 MAX41
        integer*4 key(11)
        real*4 tmp(11)
        real*4 tamp(11), tvel(11)
        integer*4 itmp(11)
        MAX41 = MAX4 + 1
        if(k .eq. 1)then
            do 100 i=1,MAX4
                a(i) = 0.0
                u(i) = 0.0
                iu(i) = 0.0
  100       continue
        endif
c-----
c       assume amplitudes arranged in decreasing order
c-----
c-----
c       we now know that the value will replace one of the amplitudes,
c       at least the lowest
c-----
        do 200 i=1,MAX4
            tmp(i) = a(i)
            tamp(i) = a(i)
            tvel(i) = u(i)
            itmp(i) = iu(i)
  200   continue
            tmp(MAX41) = ampl
            tamp(MAX41) = ampl
            tvel(MAX41) = vel
            itmp(MAX41) = ii
            call sort(tmp,key,MAX41)
            do 250 i= 1,MAX4
                kk = key(MAX41 +1 - i)
                a(i) = tamp(kk)
                u(i) = tvel(kk)
                iu(i) = itmp(kk)
                ikey(i) = kk
  250       continue
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
c      Reference: http://en.wikipedia.org/wiki/Bubble_sort
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

        subroutine getsac(name,npts,dist,deg,az,baz,t0,dt,sta,
     1                 comp,cdate,seis,permin,permax)
c-----
c
c       name    - file name to write
c       n   - number of points in FFT must be power of 2
c       n21 - number of frequencies = n/2 + 1
c       npts    - number of points in original time series
c           - which may have been zero filled to make power of 2
c       dist    - epicentral distance in km
c       deg - epicentral distance in degrees
c       az  - source - receiver azimuth in degrees
c       baz - receiver-source back azimuth
c       t0  - time of first sample after origin
c       dt  - sampling interval
c       sta - C*4 station name string
c       comp    - C*4 component name string
c       cdate   - C*12 date string
c       z   - COMPLEX array of spectra
c       permin  R*4 - minimum period to be used in trace
c       permax  R*4 - maximum period to be used in trace
c-----
        integer MXPTS
        parameter (MXPTS = 64000)
        integer  npts
        real dist, deg, az, baz, t0, dt
        character sta*8, comp*8, cdate*12
        character kstnm*8, kcmpnm*8
        character name*(*)
        real seis(MXPTS)
        real permin, permax
*
        real evla, evlo, evdp, stla, stlo, stel, beg, origtime
        integer ntimes(6)
        common /shdr/ evla,evlo,evdp,stla,stlo,stel,
     &               beg,ntimes,origtime

        common/sachdr/rhdr,ihdr,chdr
        real*4 rhdr(70)
        integer*4 ihdr(40)
        character*4 chdr(48)
*
C       write(6,*)name
        call brsac (1,MXPTS,name,seis,ierr)
C       write(6,*)'name,ierr:',name,ierr
*
        call getfhv('AZ      ',az,nerr)
        call getfhv('BAZ     ',baz,nerr)
        call getfhv('DIST    ',dist,nerr)
c-----
c       change 04 OCTOBER 2007
c       for this to work with center spread refraction data, use abs(DIST)
c       since in that case there can be negative distances
c-----
        dist = abs(dist)

        call getfhv('GCARC   ',deg,nerr)
        call getfhv('DELTA   ', dt, nerr)
        call getnhv('NPTS    ',npts,nerr)
        call getfhv('B       ', beg, nerr)
        call getfhv('O       ',origtime,nerr)
        call getfhv('EVLA    ',evla,nerr)
        call getfhv('EVLO    ',evlo,nerr)
        call getfhv('STLA    ',stla,nerr)
        call getfhv('STLO    ',stlo,nerr)
        call getkhv('KSTNM   ',kstnm,nerr)
        call getkhv('KCMPNM  ',kcmpnm,nerr)
        call getfhv('USER1', permin, ierr)
        call getfhv('USER2', permax, ierr)
C       write(6,*)'az,baz,dist,deg,dt,npts,b,o,nerr',
C     1     az,baz,dist,deg,dt,npts,beg,origtime,nerr
C       write(6,*)'evla,evlo,stla,stlo,kstnm,kcmpnm',
C     1     evla,evlo,stla,stlo,kstnm,kcmpnm
*
        if(nerr .eq. 0 .and. origtime .ne. -12345)then
            t0 = beg - origtime
        else
            t0 = beg
        end if
*
        tp = tp - origtime
        ts = ts - origtime
C       write(6,*)'t0,tp,ts:',t0,tp,ts
        sta = kstnm(1:8)
        comp = kcmpnm(1:8)
C       write(6,*)'name,npts,dist,deg,az,baz,t0,dt,sta',
C     1     name,npts,dist,deg,az,baz,t0,dt,sta
        cdate = ' '
*
*
        return
        end

        subroutine npow2(nsamp,npts,npts21)
c-----
c       Given nsamp, find npts >= nsamp such that npts is a power of 2
c-----  
        integer*4 nsamp, npts, npts21
        npts = 1
 1000   continue
        npts = 2*npts
        if(npts.lt.nsamp)go to 1000
        npts21 = npts/2 + 1
        return
        end

        subroutine getper(pmin,pmax) 
c-----
c       automatically create periods from the [pmin,pmax] limits
c       wn()    R*4 array of periods
c       NX  I   dimension of array
c       pmin    R*4 minimum period
c       pmax    R*4 maximum period
c       nper    I*4 number of periods generated
c----
        integer NPTSER, NPTFRQ, MAXTRC
        parameter (NPTSER=132000, NPTFRQ=65000)
        parameter (MAXTRC=500)
        integer MAXFRQ
        parameter (MAXFRQ=NPTFRQ)
        common /fval/ wn, nper, wnin, iwn, mper
            real*4 wn(MAXFRQ), wnin(MAXFRQ)
            integer*4 iwn(MAXFRQ)
            integer*4 nper, mper
        integer key(MAXFRQ) 

        real pmin, pmax
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
                if(p .ge. pmin .and. p .le. pmax)then
                    nper = nper + 1
                    mper = nper
                    wn(nper) = 1.0/p
                          wnin(mper) = 1.0/p
                          key(mper) = mper
                    if(nper.eq.MAXFRQ)go to 1000
                endif
 200        continue
 100    continue
 1000   continue
C       write(6,*)(wnin(i),i=1,mper)
              call sort(wnin,key,mper)
C       write(6,*)(wnin(i),i=1,mper)
        return
        end

        subroutine mapper(n,df)
c-----
c       map the desired periods into nearest Fourier frequency
c       
c       n   I*4 number of point power of 2
c       df  R*4 frequency increment
c-----
        integer NPTSER, NPTFRQ, MAXTRC
        parameter (NPTSER=132000, NPTFRQ=65000)
        parameter (MAXTRC=500)
        integer MAXFRQ
        parameter (MAXFRQ=NPTFRQ)
        common /fval/ wn, nper, wnin, iwn, mper
            real*4 wn(MAXFRQ), wnin(MAXFRQ)
            integer*4 iwn(MAXFRQ)
            integer*4 nper, mper

        fnyq = n*df/2.0
        n21 = n/2 + 1
        do 1000 i=1,mper
            j = wnin(i)/df + 1
            if(j.gt.1 .or. j.le.n21)then
                iwn(i) = j
            else
                iwn(i) = -1
            endif
 1000   continue
c-----
c       now systematically go through the list to 
c       remove duplicate values
c-----
        return
        end

        subroutine getpv(nray,tri,df,NFLOW,NFUPR,pvmin,pvmax,dist,ndist)
c-----
c       nray    I*4 number of ray parameters/phase velocities
c       tri R*4 temporary array
c       df  R*4 frequency increment
c       NFLOW   I*4 minimum search frequency
c       NFUPR   I*4 maximum search frequency
c       pvmin   R*4 - minimum velocity for plot
c       pvmax   R*4 - maximum velocity for plot
c       dist    R*4 - array of distances corresponding to each trace
c       ndist   I*4 - number of traces
c-----
        integer NPTSER, NPTFRQ, MAXTRC
        parameter (NPTSER=132000, NPTFRQ=65000)
        parameter (MAXTRC=500)
        real*4 tri(NPTFRQ), df, pvmin, pvmax
        integer*4 nray, NFLOW, NFUPR
        integer*4 ndist
        real*4 dist(ndist)

        parameter (MAXRAY=500, MAXAMP=10)
        common/phase/phvel(MAXAMP,NPTFRQ), phamp(MAXAMP,NPTFRQ),
     1      phfrq(MAXAMP,NPTFRQ), numph(NPTFRQ)
            real phvel, phamp, phfrq
            integer numph

        integer*4 idvel(MAXAMP,NPTFRQ)

        real*4 x(NPTFRQ), c(NPTFRQ)
        real*4 u(MAXAMP),a(MAXAMP)
        integer*4 iu(MAXAMP), ikey(MAXAMP)
        integer*4 MAX4, kk

        real*4 frq(NPTFRQ)
c-----
c       u   R*4 array of phase velocities
c       a   R*4 array of spectral amplitudes
c               ordered from hightest to smallest
c-----

        common/discnt/kount, lorr
            integer *4 kount, lorr
c-----
        MAX4 = MAXAMP
c-----
c       for simplicity do not buffer
c       instead grind through each frequency separately
c       The objective is to fill and array with stack values 
c       at each frequency
c       and then to search for the 10 largest peak values
c-----
        nfreq = NFUPR - NFLOW + 1
        rewind 4
        WRITE(6,*)'FIRST SORT'
        do 1000 ij=NFLOW,NFUPR
            rewind 2
            freq = (ij-1)*df
            if(mod(ij,10).eq.0)then
                write(6,*)'freq=',freq,ij,nflow,nfupr
            endif
            do 2000 j=1,nray
                read(2)c(j),(tri(ii),ii=NFLOW,NFUPR)
                x(j) = tri(ij)
 2000       continue
c-----
c       now search for the MAXAMP maxima
c-----
CRBH    write(8,*)'O',(i,x(i),i=1,nray)
        call getmax(x,c,nray,a,u,iu,ikey,kk,MAX4)
        if(kk.gt.MAX4)then
            kk = MAX4
        endif
CRBH    DO 5008 i=1,kk
CRBH    if(a(i).gt.0.0 .and.freq.gt.0.0)then
CRBH        call geterr(freq,u(i),dvel,a(i),ndist,dist)
CRBH        if(freq.ge.0.0 )then
CRBH            uper = 1.0/freq
CRBH        else
CRBH            uper = -1.0
CRBH        endif
CRBH            
CRBH        if(lorr.eq.0)then
CRBH            write(8,10)
CRBH     1          uper,u(i),dvel,a(i),i
CRBH        else if(lorr.eq.1)then
CRBH            write(8,11)
CRBH     1          uper,u(i),dvel,a(i),i
CRBH        else if(lorr.eq.2)then
CRBH            write(8,12)
CRBH     1          uper,u(i),dvel,a(i),i
CRBH        endif
CRBH        kount = kount + 1
CRBH    endif
CRBH 5008   continue
c-----
c       now zero out the x array and fill only with the peak values
c-----
        do 2500 i=1,nray
            x(i) = 0.0
 2500   continue
        if(kk.gt.MAX4)then
            kk = MAX4
        endif
CRBH    WRITE(8,*)kk,(a(i),iu(i),i=1,kk)
        do 2510 i=1,kk
            kkk = iu(i )
            x(kkk) = a(i)
 2510   continue
CRBH    write(8,*)'M',freq
CRBH    write(8,'(i5,2f10.3)')(i,x(i),c(i),i=1,nray)
        write(4)freq,(x(i),c(i),i=1,nray)
 1000   continue
c-----
c       now search for the MAX4 maxima at each phase velocity
c-----
        WRITE(6,*)'SECOND SORT'
        do 3000 ij = 1,nray
            rewind 4
            cvel = -1.0
            if(mod(ij,5).eq.0)then
                write(6,*)'Ray',ij, ' of ',nray
            endif
            do 4000 ik=NFLOW,NFUPR
                j = ik - NFLOW + 1
                read(4)frq(j),(tri(i),c(i),i=1,nray)
                x(j) = tri(ij)
                if(cvel .eq. -1.0)cvel = c(ij)
 4000   continue
CRBH    if(abs(cvel - 0.7429).lt. 0.005)then
CRBH    write(8,*)'I', cvel
CRBH    write(8,'(i5,2f10.3)')(i,x(i),frq(i),i=1,nfreq)
CRBH    endif
        call getmax(x,frq,nfreq,a,u,iu,ikey,kk,MAX4)
c-----
c       now carefully map into user units
c-----
            
        if(kk.gt.MAX4)then
            kk = MAX4
        endif
CRBH    if(abs(cvel - 0.7429).lt. 0.005)then
CRBH    write(8,*)cvel,(a(i),u(i),iu(i),ikey(i),i=1,kk)
CRBH    endif
        do 501 i=1,kk
            phvel(i,ij) = cvel
            phamp(i,ij) = a(i)
            phfrq(i,ij) = u(i)
            idvel(i,ij) = iu(i)
  501   continue
        numph(ij) = kk
c-----
c       we now have a tabulation by maximum value and frequency
c       can we now perform another seive to look at values at adjacent
c       frequencies for the same velocity??
c
c       so systematically search according to phase velocity
c       sort by frequency to arrange by frequency
c       then look for neighbors
c-----
        DO 5009 i=1,kk
        if(a(i).gt.0.0 .and.freq.gt.0.0)then
            freq = phfrq(i,ij)
            call geterr(freq,phvel(i,ij),dvel,a(i),ndist,dist)
            if(freq.ge.0.0 )then
                uper = 1.0/freq
            else
                uper = -1.0
            endif
                
            if(lorr.eq.0)then
                write(8,10)
     1          uper,phvel(i,ij),dvel,a(i),i
            else if(lorr.eq.1)then
                write(8,11)
     1          uper,phvel(i,ij),dvel,a(i),i
            else if(lorr.eq.2)then
                write(8,12)
     1          uper,phvel(i,ij),dvel,a(i),i
            endif
            kount = kount + 1
        endif
 5009   continue
   10   format('POM96 A C -1 ',g13.5,2f10.4,f13.4,i5)
   11   format('POM96 L C -1 ',g13.5,2f10.4,f13.4,i5)
   12   format('POM96 R C -1 ',g13.5,2f10.4,f13.4,i5)
 3000   continue
        return
        end

        subroutine getmax(x,c,nray,a,u,iu,ikey,kk,MAX4)
        parameter (MAXAMP=10)

        real*4 x(nray), c(nray)
        real*4 u(MAXAMP),a(MAXAMP)
        integer*4 iu(MAXAMP), ikey(MAXAMP)
        integer*4 MAX4, kk
c-----
c-----
c       search the x array to obtain the MAX4 peak amplitudes
c-----
c-----
c       obtain the MAX4 largest peaks
c-----
        kk = 0
        amp1 = x(1)
        amp2 = x(2)
        vel1 = c(1)
        vel2 = c(2)
        nm1 = nray - 1
        do 500 i=3,nm1
            amp3 = x(i)
            vel3 = c(i)
            if( (amp2-amp1).gt.1.0e-7*amp1 .and.
     1          (amp2-amp3).ge. 1.0e-7*amp3)then
                    amp4 = x(i+1)
                    vel4 = c(i+1)
                if((amp2-amp4).ge.1.0e-7*amp4) then
                    kk = kk + 1
                    ampl = amp2 
                    vel = vel2
                call srt(a,u,iu,ikey,ampl,vel,i-1,kk,MAX4)
                endif
            endif
            amp1 = amp2
            amp2 = amp3
            vel1 = vel2
            vel2 = vel3
  500   continue
        return
        end

        subroutine zfour(zarr,nn,isign,dt,df) 
c-----
c     THE input is a complex array
c     which has numbers stored in memory as
c     R1, I1, R2, I2, ..., Rnn, Inn
c     where nn must be a power of 2 R and I are the real and imaginary
c     parts of the complex number
c
c     For isign -1 this is a complex time series
c     For isign +1 this is a complex frequency series with
c        index 1 (in fortran corresponding to f=0
c              2                              f=df
c            nn/2 + 1                         f = 1/2dt = Nyquist
c            nn - 1                           f = -2df
c             nn                              f = -df

c-----
c     the cooley-tookey fast fourier transform in usasi basic fortran
c     transform(j) = sum(datc(i)*w**((i-1)(j-1)), where i and j run
c     from 1 to nn and w = exp(isign*2*pi*sqrt(-1)/nn).  datc is a one-
c     dimensional complex array (i.e., the real and imaginary parts of
c     datc are located immediately adjacent in storage, such as fortran
c     places them) whose length nn is a power of two.  isign
c     is +1 or -1, giving the sign of the transform.  transform values
c     are returned in array datc, replacing the input datc.  the time is
c     proportional to n*log2(n), rather than the usual n**2
c     rms resolution error being bounded by 6*sqrt(i)*log2(nn)*2**(-b),
c     b is the number of bits in the floating point fraction.
c
c     the program computes df from dt, dt from df and checks to see
c     if they are consistent. In addition, the transforms are multiplied
c     by dt or df to make the results dimensionally correct
c
c     This is a slightly modified version of the original Brenner routine
c     The changes were to add physical dimensions to the transform
c     and to make it all complex
c-----
        complex zarr(*) 
        integer nn, isign
        real dt, df

        complex ztemp
c-----
c       ensure that the dt and df are defined and
c       consistent
c-----
        if(dt.eq.0.0) dt = 1./(nn*df) 
        if(df.eq.0.0) df = 1./(nn*dt) 
        if(dt.ne.(nn*df)) df = 1./(nn*dt) 
c-----
c       now begin the transform
c-----
        jj = 1
        do 5 ii=1,nn 
        if(ii .lt. jj) then
              ztemp = zarr(jj)
              zarr(jj) = zarr(ii)
              zarr(ii) = ztemp
        endif
        mm = nn/2
    3   continue
        if(jj.le.mm) then
            go to 55
        else 
              jj = jj - mm
              mm = mm /2
              if(mm.lt.1)then
                  go to 55
              else
                  go to 3
              endif
        endif
   55   continue
        jj = jj + mm
    5   continue
        mmax = 1 
c-----
    6 continue
        if(mmax .lt. nn)then
            go to 7
        else if(mmax .ge. nn)then
            go to 10
        endif
    7   continue
        istep= 2*mmax 
        theta = 6.283185307/(isign*2.0*mmax) 
        sinth=sin(theta/2.) 
        wstpr=-2.*sinth*sinth 
        wstpi=sin(theta) 
        wr=1.0 
        wi=0.0 
        do 9 mm=1,mmax
              do 8 ii=mm,nn,istep
                    jj=ii+mmax
                    ztemp=cmplx(wr,wi)*zarr(jj)
                    zarr(jj) = zarr(ii) - ztemp
                    zarr(ii) = zarr(ii) + ztemp
    8         continue
c-----
c       use trig relations to compute the next sin/cos
c       without actually calling sin() or cos()
c-----
              tempr = wr 
              wr = wr*wstpr-wi*wstpi + wr 
              wi = wi*wstpr+tempr*wstpi + wi 
    9   continue
        mmax = istep 
        go to 6 
c-----
c       transform is done
c-----
   10   continue 
c-----
c     give the arrays the proper physical dimensions
c-----
        if(isign.lt.0)then
c-----
c             time to frequency domain
c-----
              do  ii = 1,nn
                    zarr(ii) = zarr(ii) * dt
              enddo
        else
c-----
c             frequency to time domain
c-----
              do ii = 1,nn
                    zarr(ii) = zarr(ii) * df
              enddo
        endif
        return
        end
