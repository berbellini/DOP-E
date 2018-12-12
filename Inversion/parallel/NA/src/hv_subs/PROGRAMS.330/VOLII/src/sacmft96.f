        program sacmft96
c----------------------------------------------------------------------c
c                                                                      c
c       COMPUTER PROGRAMS IN SEISMOLOGY                                c
c       VOLUME II                                                      c
c                                                                      c
c       PROGRAM: SACMFT96                                              c
c                                                                      c
c       COPYRIGHT 1999                                                 c
c       R. B. Herrmann                                                 c
c       Department of Earth and Atmospheric Sciences                   c
c       Saint Louis University                                         c
c       221 North Grand Boulevard                                      c  
c       St. Louis, Missouri 63103                                      c
c       U. S. A.                                                       c
c                                                                      c
c----------------------------------------------------------------------c
c
c       CHANGES
c
c       10 OCT 2000 - introduced the use of the USER1 and USER2
c                 fields in the SAC header to indicate the
c                 minimum and maximum permissible periods
c                     This is done so that the phase match filtered
c                     trace is not used to get dispersion for 
c                     periods outside these limits
c       30 JAN 2002 - changed mft96.dsp format to include the
c                     source->receiver azimuth between distance and
c                     spectral amplitude
c       01 APR 2002 - line changed in setxy so that linear period axis
c               works correctly (meijian@iag.usp.br)
c       23 AUG 2003 - changes method of picking extrema since
c               it did not work in some cases
c               now just use first difference
c      
c                     changed name of output dispersion file 
c                     from mft96.dsp to
c               mft96.DSP for easier cp later
c       13 SEP 2003 - added -s flag for log/lin color mapping - also 
c               changed
c               original color -> amplitude mapping
c       04 DEC 2003 - correct pltgv routine so that points do not plot
c               outside the plot region
c       20 MAR 2004 - output format change for write(LOT)
c       24 MAY 2004 - if none of the dates are set default to
c               1970 01 01 00 00 00 000
c       03 JUN 2004 - type in line 805 had xmx not the desired vmx
c       30 JUN 2005 - changed name of output dispersion file 
c                     from mft96.DSP to
c               mft96.disp to avoid problems with caseless file systems,
c               e.g., CYGWIN under windows
c               - attempt to 
c       07 FEB 2007 - changed format for the mft96,disp file to permit
c               -12345.0 for AZ, BAZ and GCARC
c       29 MAY 2009 - fixed roundoff problem in the if (p.ge.pmin etc
c-----
        integer LIN, LOT, LER
        parameter (LIN=5, LER=0, LOT=6)
        integer NPTSER, NPTFRQ
        parameter (NPTSER=132000, NPTFRQ=65000)
c-----
c       information about plot on the page
c-----
        common/pltinf/gxl(3),gyl(3),gxh(3),gyh(3),
     1      axl(3),ayl(3),axh(3),ayh(3),ixlnlg(3),iylnlg(3)
            real*4 gxl,gyl,gxh,gyh,
     1              axl,ayl,axh,ayh
            integer*4 ixlnlg,iylnlg
        common/pltcnt/laxper
            logical laxper
        common/discnt/kount, lorr
            integer *4 kount, lorr
        common/tit/npts,az,dt,dist,n,a0,a1
            integer*4 npts
            real*4 az,dt,dist,a0,a1
        integer MAXFRQ
        parameter (MAXFRQ=100)
        common /fval/ wn, nper, wnin, mper
            real*4 wn(MAXFRQ), wnin(MAXFRQ)
            integer*4 nper, mper
        common/timser/datas,data,datad,data2,datad2
        complex datas(NPTFRQ),data(NPTSER),datad(NPTSER)
        complex data2(NPTSER),datad2(NPTSER)
        real tarr(NPTSER)
        common/gtf/doabs,ampmx,doscl
        logical doabs, doscl
        integer*4 n , n21 
        real*4  deg , baz , t0
        character sta*8, comp*8, cdate*12
        character*256 sacfile
        integer ierr
        logical ylin, xlin, verbos, shadon
        character *3 lnlg(2)
c-----
c       iunit   I*4 Physical units of amplitude
c               -1 counts (default)
c               0 meters
c               1 centimeters
c               2 nanometers
c       idva    I*4 input time series (note output spectra will be 
c               spectra is displacement in cm, e.g., cm-sec of count-sec
c               0   displacment
c               1   velocity
c               2   acceleration
c       unitcor R*4 unit correction factor to get to cm
c-----
        common/units/idva,iunit,unitcor
        integer*4 idva, iunit
        real*4 unitcor

        real*4 pmin, pmax
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
c       parse command line
c-----
        call gcmdln(xlin,ylin,pmin,pmax,a0,a1,laxper,
     1      verbos,vmin,vmax,shadon,lorr,idva,iunit,
     2      unitcor,sacfile)
c-----
c     shadon process data
c-----
        call gettrc(sacfile,n,n21,npts,dist,deg,az,baz,t0,dt,sta,
     1      comp,cdate,datas,tarr,NPTSER,ierr,permin,permax)
C       write(6,*)'sacfile,n,n21,npts,dist,deg,az,baz,t0,dt,sta',
C     1     sacfile,n,n21,npts,dist,deg,az,baz,t0,dt,sta
c
c-----
c       determine filter periods
c-----
        call getper(wn,MAXFRQ,pmin,pmax,mper,permin,permax) 
        nper = mper
        if(nper.eq.0)then
        write(LOT,*)
     1   'The command line said to use periods between',pmin,' and ',
     2   pmax
        write(LOT,*)
     1   'The USER1 and USER2 fields of the sac file  said to use ',
     2   ' periods between',permin,' and ', permax
        write(LOT,*)
     1   'This program fails if there is no overlapping period range'
        stop
        endif
        call pinitf('MFT96.PLT')
        ISTART = 1
        ISTOP = npts
        if(verbos)then
              write(LOT,*)'n      =',n
              write(LOT,*)'dt     =',dt
              write(LOT,*)'dist   =',dist
              write(LOT,*)'deg    =',deg 
              write(LOT,*)'a0     =',a0  
              write(LOT,*)'a1     =',a1  
              write(LOT,*)'vmin   =',vmin   
              write(LOT,*)'vmax   =',vmax   
              write(LOT,*)'pmin   =',pmin   
              write(LOT,*)'pmax   =',pmax   
              write(LOT,*)'xlin   =',xlin   
              write(LOT,*)'ylin   =',ylin   
              write(LOT,*)'t0     =',t0      
              write(LOT,*)'mper   =',mper    
              write(LOT,*)'Shading=',shadon
        endif
c-----
c       open file containing dispersion values
c-----
        open(8,file='mft96.disp',access='sequential',
     1      form='formatted',status='unknown')
        rewind 8
c-----
c       begin the multiple filter analysis
c-----
        call gpv(ISTART,ISTOP,vmin,vmax,pmin,pmax,xlin,ylin,t0,shadon)

        close (8)
c-----
c       open file containing plot control information
c-----
c 249
cXAXIS-PERIOD
c 1.300  1.500  4.300  5.500   3.24 0.100E-08 324. 0.100E-05 log log MFT96.PLT
c 5.100  1.500  9.100  5.500   14.0  2.00      75.0  5.00    lin lin MFT96.PLT
c 9.200  1.500  9.700  5.500 -0.250 0.242E+04 0.250  374.    lin lin MFT96.PLT
c 32
c   75.00000       70.00000       65.00000       60.00000       55.00000
c   50.00000       48.00000       46.00000       44.00000       42.00000
c   40.00000       38.00000       36.00000       34.00000       32.00000
c   30.00000       29.00000       28.00000       27.00000       26.00000
c   25.00000       24.00000       23.00000       22.00000       21.00000
c   20.00000       19.00000       18.00000       17.00000       16.00000
c   15.00000       14.00000                                                   
c-----
c Value     Name    Meaning
c 249       kount   Number of dispersion points in mft96.dsp 
c XAXIS-PERIOD  strong  Key word for plot or XAXIS-FREQUENCY
c three lines to describe the three images on the file MFT96.PLT
c   each line has the physical coordinates of the plot (XL,YL) (XH,YH)
c   the physical values at those corners  
c   (XV,YV) at lower left and (XV,YV) at upper right
c   two strongs to describe the plotted axes - for the first figure, 
c   the spectra plot
c       the plot is log-log while for the dispersion curve (center) 
c   and trace plot(thrid) the
c       axes are linear-linear. Finally the graphics file that 
c   contains these images is
c       MFT96.PLT for all three plots
c  32   nper    number of period values
c this is now followed by 32 periods for the plot
c       
c-----
        open(3,file='mft96.ctl',access='sequential',
     1      form='formatted',status='unknown')
        rewind 3
        write(3,*)kount
        if(laxper)then
            write(3,14)
        else
            write(3,13)
        endif
        do 11 i=1,3
        write(3,10) gxl(i),gyl(i),gxh(i),gyh(i),
     1      axl(i),ayl(i),axh(i),ayh(i),
     2      lnlg(ixlnlg(i)),lnlg(iylnlg(i))
   11   continue
   10   format(4f7.3,4g11.3,' ',a3,' ',a3,' MFT96.PLT')
   13   format('XAXIS-FREQUENCY')
   14   format('XAXIS-PERIOD')
        write(3,*)nper
        do 1234 i=1,nper
            wnin(i) = 1.0/wn(i)
 1234   continue
        write(3,'(5g15.7)')(wnin(i),i=1,nper)
c-----
c       close the plot parameter file
c-----
        close (3)
        if(lorr.gt.0)then
c-----
c       open file for invoking  DPEGN96 control information
c-----
C-----
C       DOS
C-----
C       open(3,file='MFT96CMP.BAT',access='sequential',
C     1     form='formatted',status='unknown')
C       rewind 3
C-----
c-----
c       UNIX
c-----
        open(3,file='MFT96CMP',access='sequential',
     1      form='formatted',status='unknown')
        rewind 3
        write(3,'(a)')'#!/bin/sh'
        write(3,'(a)')' '
        lssf = lgstr(sacfile)
        write(3,17)sacfile(1:lssf), a0
        write(3,'(a)')' '
        outstr = ' '
   17 format('# Data file = ',a/'# alpha = ',f10.3)
   15   format('sdpegn96 -X0 ',f5.2,' -Y0 ',f5.2,' -XLEN ',
     1   f5.2, ' -YLEN ',f5.2,' -XMIN ',g10.3,
     1   ' -XMAX ',g10.3, ' -YMIN ',f5.2,' -YMAX ',f5.2)
        write(outstr,15)gxl(2),gyl(2),gxh(2)-gxl(2),gyh(2)-gyl(2),
     1      axl(2), axh(2), ayl(2), ayh(2)
        ls = lgstr(outstr)
        if(laxper)then
            outstr(ls+1:ls+7) = ' -PER  '
        else
            outstr(ls+1:ls+7) = ' -FREQ '
        endif
        ls = ls + 7
        if(lorr.eq.1)then
            outstr(ls+1:ls+13) = '-L -U -NOBOX '
        else if(lorr.eq.2)then
            outstr(ls+1:ls+13) = '-R -U -NOBOX '
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
        write(3,'(a)')' '
        outstr = ' '
   16   format('sacspc96 -X0 ',f5.2,' -Y0 ',f5.2,' -XLEN ',
     1   f5.2, ' -YLEN ',f5.2,' -XMIN ',g10.3,
     1   ' -XMAX ',g10.3, ' -YMIN ',g10.3,' -YMAX ',g10.3)
        write(outstr,16)gxl(1),gyl(1),
     1      gxh(1)-gxl(1),gyh(1)-gyl(1),axl(1),
     1      axh(1), ayl(1), ayh(1)
        ls = lgstr(outstr)
        if(laxper)then
            outstr(ls+1:ls+7) = ' -PER  '
        else
            outstr(ls+1:ls+7) = ' -FREQ '
        endif
        ls = ls + 7
        outstr(ls+1:ls+8) = ' -NOBOX '
        ls = ls + 8
        if(xlin)then
            outstr(ls+1:ls+6)='-XLIN '
        else
            outstr(ls+1:ls+6)='-XLOG '
        endif
        ls = ls + 6
c-----
c       spectral amplitude is always on log scale
c-----
        outstr(ls+1:ls+6)='-YLOG '
        ls = ls + 6
        outstr(ls+1:ls+3)='-f '
        ls = ls + 3
        outstr(ls+1:ls+lssf)= sacfile(1:lssf)
        ls = ls + lssf
        
        write(3,'(a)')outstr(1:ls)

c-----
c       close the SDPEGN96 invocation file
c-----
        close (3)
c-----
c       UNIX
c-----
C       iret = system('chmod +x MFT96CMP')
        endif
c-----
c       terminate program
c-----
        call pend()
        end

        subroutine gpv(ISTART,ISTOP,pumin,pumax,pmin,pmax,xlin,ylin,t0,
     1      shadon)
        integer LOT
        parameter (LOT=6)
        integer NPTSER, NPTFRQ
        parameter (NPTSER=132000, NPTFRQ=65000)
c-----
c       general subroutine for performing 
c       Gaussian Filter Group Velocity analysis
c
c       REFERENCES:
c             Herrmann, R. B. (1973). Some aspects of band-pass
c                   filtering of surface waves, Bull. Seism. Soc.
c                   Am. 63, 663-671 (spectral amplitudes)
c
c             Gaussian filters introduced by Landismann and Dziewonski
c             modified by others, e.g., filter parameter dependent
c             upon frequency
c
c-----
c       SUBROUTINE ARGUMENTS:
c
c       n      I*4      number of points for FFT, power of 2
c       d      R*4      sampling interval in seconds
c       dist   R*4      distance from source to receiver
c       NSAMP  I*4      number of samples in trace <=n
c       ISTART I*4      process time series starting at this sample
c       ISTOP  I*4      process time series ending at this sample
c       a0     R*4      defines Gaussian filter parameter
c       a1     R*4            alpha = a0+a1*centfreq**2
c       pumin  R*4      minimum group velocity for plots
c       pumax  R*4      maximum group velocity for plots
c       xlin  L        .true. linear frequency scale
c                       .false. logarithmic frequency scale
c       ylin   L        .true. linear group velocity scale
c                       .false logarithmic velocity scale
c       t0     R*4      absolute time of first sample
c                       with respect to source origin time
c-----
c       CONTROL COMMON BLOCKS FROM MAIN PROGRAM
c
c-----
c
        common/timser/datas,data,datad,data2,datad2
        complex datas(NPTFRQ),data(NPTSER),datad(NPTSER)
        complex data2(NPTSER),datad2(NPTSER)
        real*4 xamp(NPTSER), xfrq(NPTSER)
        real*4 xampdt(NPTSER)
        equivalence(data(1),xamp(1))
        equivalence(datad(1),xfrq(1))
        integer*4 ISTART, ISTOP
        integer MAXFRQ, MAXVEL
        parameter (MAXFRQ=100, MAXVEL=10)
        common /fval/ wn, nper, wnin, mper
              real*4 wn(MAXFRQ), wnin(MAXFRQ)
              integer*4 nper, mper
        common/tit/nsamp,az,dt,dist,n,a0,a1
              integer*4 nsamp
              real*4 az,dt,dist,a0,a1
        common/pltcnt/laxper
              logical laxper
        common/group/gpvel(MAXVEL,MAXFRQ), gpamp(MAXVEL,MAXFRQ),
     1      numgp(MAXFRQ)
            real gpvel, gpamp
            integer numgp
        real*4 pumin,pumax
        logical ylin,xlin,shadon
        common/gtf/doabs,ampmx,doscl
        logical doabs,doscl
c-----
c       iunit   I*4 Physical units of amplitude
c               -1 counts (default)
c               0 meters
c               1 centimeters
c               2 nanometers
c       idva    I*4 input time series (note output spectra will be 
c               spectra is displacement in cm, e.g., cm-sec of count-sec
c               0   displacment
c               1   velocity
c               2   acceleration
c       unitcor R*4 unit correction factor to get to cm
c-----
        common/units/idva,iunit,unitcor
        integer*4 idva, iunit
        real*4 unitcor

c-----
c       array for contour plot
c       real*4  cont(NX,NY) the array to be plotted
c       real*4  x(NX)   x-axis values for contour
c               must be in increasing order
c       real*4  y(NY)   y-axis values for contour
c               must be in increasing order
c-----
        integer NX, NY
        parameter(NX=100,NY=50)
        real cont(NX,NY)
              
c-----
c       initialize array for ploting group velocity and
c       spectral amplitudes
c-----
        DO 5000 I = 1,MAXFRQ
            numgp(i) = 0
            DO 5001 J = 1,MAXVEL
                gpvel(J,I) = 0.0
                gpamp(J,I) = 0.0
 5001       continue
 5000       continue
c-----
c       initialize array for group velocity contours
c-----
        do 6001 i=1,NX
            do 6002 j=1,NY
                cont(i,j) = 0.0
 6002       continue
 6001   continue
c-----
c       set up frequencies to be used
c       these are in order of increasing frequency
c       these are either automatically gnerated from an
c       equation or provided in a table
c-----

C       write(6,*)'mper,nper,wnin,dt,n,wmax,wmin,df:',
C     1     mper,nper,wnin,dt,n,wmax,wmin,df
        call mkfreq(mper,nper,wn,dt,n,wmax,wmin,df)
C       write(6,*)'mper,nper,wn,dt,n,wmax,wmin,df:',
C     1     mper,nper,wn,dt,n,wmax,wmin,df
            if(nper.eq.0)then
            write(LOT,*)'GPV: nper = 0 check aa,bb,cc'
            return
        endif
C       write (6,*)'NPER,DF:',nper,df
c-----
c             if(laxper)then
c                   dum = wmax
c                   wmax = 1.0/wmin
c                   wmin = 1.0/dum
c             endif
c-----
c       now we are ready to plot
c-----
c-----
c       narrow band pass filter at each frequency
c       and store filtered time series in the array f(4096,80)
c-----
c       narrow band pass filter using gaussian filter, passing
c       rejecting frequencies for which the filter response is
c       less than exp(-3.1415927)
c-----
        n21 = n/2 + 1
        ampmx = 0.0
        do 999 jump = 1,nper
            alpha = a0 + a1*wn(jump)**2
            call dofilt(datas,data,datad,dt,df,
     1          n,n21,alpha,idva,ampmx,wn(jump),pumin,pumax,
     2          dist,t0,nsamp)
c-----
c           for purpose of bias correction, compute
c           also for 2*alpha
c-----
C           call dofilt(datas,data2,datad2,dt,df,
C     1     n,n21,alpha+alpha,idva,ampmx2,wn(jump),pumin,pumax,
C     2         dist,t0)
c-----
c           here map on to the cont array
c-----
        call mapgv(jump,nper,xamp,xampdt,xfrq,n,cont,NX,NY,nsamp,
     1      t0,dist,az,dt,pumin,pumax,laxper,wn(jump),alpha)
  999   continue
        if(ampmx.eq.0.0)ampmx=1.0
c-----
c       plot dispersion curve
c------
            if(doabs)then
c-----
c           normalize
c-----
                do 4110 j=1,NY
                    do 4111 i=1,NX
                        cont(i,j) = cont(i,j)/ampmx
 4111               continue
 4110           continue
            else
c-----
c           do relative plot
c-----
                do 4112 i=1,NX
                    amprel = 0.0
                    do 4113 j=1,NY
                        if(cont(i,j).gt.amprel)then
                            amprel = cont(i,j)
                        endif
 4113               continue
                    if(amprel.eq.0.0)amprel = 1.0
                    do 4114 j=1,NY
                        cont(i,j) = cont(i,j)/amprel
 4114               continue
 4112           continue
c-----
            endif
c-----
c       define mapping of array coordinates to absolute position
c-----
C       write(6,*)'WN  :',(wn(i),i=1,nper)
            call conplt(n,wn,nper,df,t0,dt,dist,
     1          nsamp,wmin,wmax,
     2          pumin,pumax,xlin,ylin,a0,a1,shadon,fmin,fmax,
     3          cont,pmin,pmax)
c-----
c       plot spectral amplitudes
c-----
        xaxlen = 3.0
        yaxlen = 4.0
        call pltamp(wn,nper,.false.,.false.,xaxlen,yaxlen,laxper,
     1      pumin,pumax)
        return
        end

        subroutine band(z,n,i1,i2,nb,nt)
        parameter(NP=4096,NP2=2049)
        complex z(NP)
        real fparz
c-----
c     apply a parzen window to the complex function z
c
c     z - complex array to be windowed
c     n - number of points in the z array
c     i1,i2 - passband of windowing function in units
c              of array index
c            the original array is unchanged between i1 and i2
c     nb - taper width for parzen window about i1, i2
c         the array is zeroed for indices less than i1 - nb + 1
c                   and for indices greater than    i2 + nb - 1
c     nt = 0 assume spectra input, and only window
c             the positive array frequencies, e.g.,
c             the first half of the array. This is
c             in fact not a true parzen window, but
c             we use parzen half-windows to taper the
c             band edges
c         = 1 true parzen window about the center of the array
c             used here for the pseudo-autocorrelation functions
c-----
c-----
c     define center point, corners, and check for valid 
c     choices of corners
c-----
        if(nt.eq.0)then
              m = n / 2
        else
              m = n
        endif
        if(nb.lt.2)then
              kl = 2
        else
              kl = nb
        endif
        ib1 = i1 - kl + 1
        if(ib1.lt.1)ib1=1
        ib2 = i1
        ib3 = i2
        ib4 = i2 + kl -1
        if(ib4.gt.m)ib4=m
c-----
        do 100 i = 1 , m
              if(i.lt.ib1)then
                    z(i) = cmplx(0.0,0.0)
              else if(i.ge.ib1 .and. i.lt.ib2)then
                    fac = fparz(ib2-i,ib2-ib1)
                    z(i) = fac * z(i)
              else if(i.ge.ib2 .and. i.lt.ib3)then
                    z(i) = 1.0*z(i)
              else if(i.ge.ib3 .and. i.lt.ib4)then
                    fac = fparz(i-ib3,ib4-ib3)
                    z(i) = fac * z(i)
              else
                    z(i) = cmplx(0.0,0.0)
              endif
  100   continue
        return
        end

        function fparz(i,iw)
c-----
c       parzen windowing function
c
c        iw = window halfwidth
c        i  = index within window
c             i = 0 corresponds to the passband = 1
c-----
              fi = i
              fiw = iw
              rat = abs(fi/fiw)
              if(2*i .lt. iw)then
                    fparz = 1.0 -6.*rat*rat*(1.-rat)
              else
                    fparz = 2.*(1.0-rat)**3
              endif
        return
        end

        subroutine srt(a,u,ufr,vel,ampl,frq2,k,MAX4)

c-----
c       subroutine to maintain a list of the MAX4
c       largest values of ampl in a(i)
c       also saving the corresponding vel in v(i)
c
c       This is being done to avoid have two very large
c       arrays in overhead and because we do not need
c       to do a full sort
c-----
        parameter (MAXVEL=10)
        real*4 a(MAXVEL), u(MAXVEL), ufr(MAXVEL), vel, ampl
        integer*4 MAX4, k
        integer*4 MAX41
        integer*4 key(11)
        real*4 tmp(11)
        real*4 tamp(11), tvel(11), tfrq(11)
        MAX41 = MAX4 + 1
        if(k .eq. 1)then
            do 100 i=1,MAX4
                a(i) = 0.0
                u(i) = 0.0
                ufr(i) = 0.0
  100       continue
        endif
        if(k.lt.MAX4)then
              kup = k
        else
              kup = MAX4
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
              tfrq(i) = ufr(i)
  200       continue
              tmp(MAX41) = ampl
              tamp(MAX41) = ampl
              tvel(MAX41) = vel
              tfrq(MAX41) = frq2
              call sort(tmp,key,MAX41)
              do 250 i= 1,MAX4
                    kk = key(MAX41 +1 - i)
                    a(i) = tamp(kk)
                    u(i) = tvel(kk)
                    ufr(i) = tfrq(kk)
  250             continue
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

        subroutine conplt(n,wn,nper,df,t0,dt,dist,
     1       nsamp,wmin,wmax,
     2  pumin,pumax,xlin,ylin,a0,a1,shadon,fmin,fmax,
     3          cont,pmin,pmax)
c-----
c       cont    R*4 NX x NY array to be contoured
c       NX  I*4 array dimension
c       NY  I*4 array dimension
c       nper    I*4 actual number of array elements used
c       wn  R*4 array of frequencies
c-----
        integer NPTSER, NPTFRQ
        parameter (NPTSER=132000, NPTFRQ=65000)
        integer BLACK, BLUE
        parameter (BLACK=1, BLUE=4)
        common/timser/datas,data,datad,data2,datad2
        complex datas(NPTFRQ),data(NPTSER),datad(NPTSER)
        complex data2(NPTSER),datad2(NPTSER)
        integer MAXFRQ, MAXVEL
        parameter (MAXFRQ=100, MAXVEL=10)
        real*4 wn(MAXFRQ)
              integer*4 n
              integer*4 nper
              real*4 df
              real*4 t0,dt,dist
              integer*4 nsamp
              logical ylin,xlin
        common/pltcnt/laxper
              logical laxper
        common/group/gpvel(MAXVEL,MAXFRQ), gpamp(MAXVEL,MAXFRQ),
     1      numgp(MAXFRQ)
            real gpvel, gpamp
            integer numgp
        logical shadon
c-----
c       information about plot on the page
c-----
        common/pltinf/gxl(3),gyl(3),gxh(3),gyh(3),
     1      axl(3),ayl(3),axh(3),ayh(3),ixlnlg(3),iylnlg(3)
            real*4 gxl,gyl,gxh,gyh,
     1              axl,ayl,axh,ayh
            integer*4 ixlnlg,iylnlg
        common/gtf/doabs,ampmx,doscl
        logical doabs, doscl
c-----
c       
c-----
        integer NX, NY
        parameter(NX=100,NY=50)
        real cont(NX,NY), x(NX), y(NY)
        parameter (NC=11)
        integer icol(NC+1)
        real CN(NC)


        xoff = 5.1
        yoff = 1.5
        call plot(xoff,yoff,-3)
c-----
c       normalize each trace individually
c-----
c-----
c       set up plot parameters
c-----
c       if xlin = .true. linear frequency
c             else      logarithmic
c       if ylin = .true. linear group velocity
c             else logarithmic
c-----
        xaxlen = 4.0
        yaxlen = 4.0
        call gbox(0.0,0.0,xaxlen,yaxlen)
        call setxy(x,NX,y,NY,wn,nper,pmin,pmax,wmin,wmax,pumin,pumax,
     1      xlin,ylin,xlow,xhgh,ylow,yhgh,xaxlen,yaxlen)
c-----
c
c-----
        call plot(0.0,0.0,3)
        dx = 1.0
        dy = 1.0
        if(doscl)then
        do 1600 i=1,NX
            do 1601 j=1,NY
                    cont(i,j) = 10.0*cont(i,j)
                    if(cont(i,j).gt.0.0)then
                        cont(i,j) = alog10(cont(i,j)) 
                    else
                        cont(i,j) = -40.0
                    endif
 1601       continue
 1600   continue
        endif
        if(doscl)then
            vmx = 1.0
            vmn = -0.7
        else
            vmx = 0.9
            vmn = 0.0
        endif
        
c-----
c       define color array
c-----
        do 2000 i = 1, NC
            CN(I)= 1.0*vmn  + real(i-1)/(real(NC-1))*(vmx - vmn) 
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
c       so that if frequency is plotted do it on the mapgv
c-----
        if(shadon) then
            call FARB2D(X,nper,Y,NY,cont,NX,CN,ICOL,NC,2)  
        else
            call FARB2D(X,nper,Y,NY,cont,NX,CN,ICOL,NC,1) 
        endif
C       write(6,*)'ICOL:',ICOL
C       write(6,*)'CN  :',CN
C       write(6,*)'X   :', (x(i),i=1,nper)
C       write(6,*)'Y   :', (y(i),i=1,NY)
c-----
c       PLOT THE GROUP VELOCITIES OF THE FOUR LARGEST SPECTRAL
c       AMPLITUDES ON TOP OF THE CONTOUR PLOT. WE USE THE
c       ARRAY SET UP BY mxval AND STORED IN common/group/gpval
c-----
        call newpen(BLACK)
        call pltgv(x,nper,ylin,yhgh,ylow,yaxlen,laxper)
c-----
c       obtain original Time series
c       after Parzen Windowing
c-----
              n21 = n/2 + 1
              do 6003 i=1,n21
                    data(i) = datas(i)
 6003       continue
c-----
c     apply Parzen window to spectra using processing limits as a guide
c-----
              nb = 10
              df = 1.0/(n*dt)
              fl = 6.0*df
              fh = 0.5/dt
              i1 = (fl/df + 1.0)
              i2 = (fh/df + 1.0)
              if(i2.gt.n21)i2=n21
              if(i1.lt.0)i1=0
              nt = 0
              if(i1.gt.0)then
                    call band(data,n,i1,i2,nb,nt)
              endif
              do 6004 i=1,n21
                    if(i.gt.1)then
                          data(n+2-i)=conjg(data(i))
                    endif
 6004       continue
              call four(data,n,+1,dt,df)
              trimax = 0.0
        DO 6012 KK=1,nsamp
              if(abs(real(data(KK))).gt.trimax)then
                    trimax = abs(real(data(KK)))
              endif
 6012       continue
        ipen = 3
        fold  = 0.0
        vold =  0.0
c-----
c plot trace to right on a linear scale
c
c NOTE: 02 APR 2002 the figure is replotted to plot the 
c actual trace within the
c       specified group velocity window
c-----
        x1 = xaxlen + 0.1
        x2 = x1 + 0.5 + 0.1
        y1 = 0.0
        y2 = yaxlen
        call box(x1,x2,y1,y2)
        ipen = 3
        tmax=dist/ylow
        tmin=dist/yhgh
        tend=(nsamp -1)*dt + t0
        nmin=(tmin-t0)/(tend - t0)*(nsamp - 1) +1
        nmax=(tmax-t0)/(tend - t0)*(nsamp - 1) +1
        if(nmin.lt.0)nmin = 1
        if(nmax.gt.nsamp)nmax = nsamp

        DO 6014 KK=nmin,nmax
             xx=xaxlen + 0.0 + 0.25 + 0.1 + 0.25*real(data(KK))/trimax
              yy = yaxlen -(KK-nmin)*yaxlen/(nmax-nmin + 1)
              call plot(xx,yy,ipen)
              ipen = 2
 6014       continue
c-----
c       put maximum amplitude at base of trace
c-----
              call number(x1,y1-0.20,0.07,trimax,0.0,1003)
              call plot(0.0,0.0,3)
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
        gxl(3) = xoff + 0.1 +xaxlen 
        gyl(3) = yoff
        gxh(3) = xoff + 0.1 + xaxlen + 0.5
        gyh(3) = yoff + yaxlen
        axl(3) = -0.25
        axh(3) =  0.25
        ayl(3) = (nmax -1)*dt + t0
        ayh(3) = (nmin -1)*dt + t0
        ixlnlg(3) = 1
        iylnlg(3) = 1
        return
        end

        subroutine four(data,nn,isign,dt,df) 
        integer NPTSER
        parameter (NPTSER=132000)
        dimension data(NPTSER*2) 
        n = 2 * nn 
        if(dt.eq.0.0) dt = 1./(nn*df) 
        if(df.eq.0.0) df = 1./(nn*dt) 
        if(dt.ne.(nn*df)) df = 1./(nn*dt) 
        j = 1 
        do 5 i=1,n,2 
C       if(i-j)1,2,2 
        if(i .lt. j) then
            go to 1
        else
            go to 2
        endif
    1 tempr = data(j) 
        tempi = data(j+1) 
        data(j) = data(i) 
        data(j+1)=data(i+1) 
        data(i) = tempr 
        data(i+1) = tempi 
    2 m = n/2 
C    3 if(j-m) 5,5,4 
    3 continue
        if(j.le.m) then
            go to 5
        else 
            go to 4
        endif
    4 j = j-m 
        m = m/2 
C       if(m-2)5,3,3 
        if(m.lt.2)then
            go to 5
        else
            go to 3
        endif
    5 j=j+m 
        mmax = 2 
C    6 if(mmax-n) 7,10,10 
    6 continue
        if(mmax .lt. n)then
            go to 7
        else if(mmax .ge. n)then
            go to 10
        endif
    7 istep= 2 *mmax 
        theta = 6.283185307/float(isign*mmax) 
        sinth=sin(theta/2.) 
        wstpr=-2.*sinth*sinth 
        wstpi=sin(theta) 
        wr=1.0 
        wi=0.0 
        do 9 m=1,mmax,2 
        do 8 i=m,n,istep 
        j=i+mmax 
        tempr=wr*data(j)-wi*data(j+1) 
        tempi=wr*data(j+1)+wi*data(j) 
        data(j)=data(i)-tempr 
        data(j+1)=data(i+1)-tempi 
        data(i)=data(i)+tempr 
    8 data(i+1) = data(i+1)+tempi 
        tempr = wr 
        wr = wr*wstpr-wi*wstpi + wr 
    9 wi = wi*wstpr+tempr*wstpi + wi 
        mmax = istep 
        go to 6 
   10 continue 
        if(isign.lt.0) go to 1002 
c     frequency to time domain 
        do 1001 iiii = 1,n 
 1001 data(iiii) = data(iiii) * df 
        return 
 1002 continue 
c     time to frequency domain 
        do 1003 iiii = 1,n 
 1003 data(iiii) = data(iiii) * dt 
        return 
        end 

        subroutine box(x1,x2,y1,y2)
              call plot(x1,y1,3)
              call plot(x2,y1,2)
              call plot(x2,y2,2)
              call plot(x1,y2,2)
              call plot(x1,y1,2)
        return
        end

        function fmap(x,df)
              real*4 fmap
              real*4 x
              real*4 df
        integer MAXFRQ
        parameter (MAXFRQ=100)
        common /fval/ wn, nper, wnin, mper
              real*4 wn(MAXFRQ), wnin(MAXFRQ)
              integer*4 nper, mper
c-----
c       this general function defines the mapping between 
c       array index I and the center frequency of the 
c gaussian filter
c
c       changing this routine will change the frequency sampling
c       throughout the program
c-----
              il = x
              if(il.lt.1)il=1
              iu = il+1
              if(iu.gt.mper)then
                    iu = mper
                    il = iu -1
              endif
              fmap = wnin(il) + (x-il)*(wnin(iu)-wnin(il))
        return
        end

        subroutine wnstup(cmdfil)
c-----
c     open the command file to get periods or frequencies
c     for the search rather than automatically generating them
c-----
        integer LER
        parameter (LER=0)
        integer MAXFRQ
        parameter (MAXFRQ=100)
        common /fval/ wn, nper, wnin, mper
              real*4 wn(MAXFRQ), wnin(MAXFRQ)
              integer*4 nper, mper
c-----
c       laxper  .true.  plot velocity vs period
c       laxper  .false. plot velocity vs frequency
c-----
        common/pltcnt/laxper
              logical laxper
        common/discnt/kount, lorr
              integer *4 kount, lorr
        character cmdfil*(*)
        integer*4 key(MAXFRQ)
        logical ext
c-----
c     verify that the data file exists
c-----
        inquire(file=cmdfil,exist=ext)
        if(.not.ext)then
              write(LER,*)'Period File does not exist'
              call usage(0)
        endif
              open(2,file=cmdfil,status='old', form='formatted',
     1            access='sequential')
              rewind 2
              read(2,*)ifrper
c-----
c     ifrper = 0 input is period
c             = ! 0 input is frequency
c-----
              mper = 0
              do 100 i=1,MAXFRQ
                    read(2,*,end=101,err=101)xx
                    if(xx.gt.0.0)then
                          if(ifrper.eq.0)xx=1.0/xx
                          mper = mper + 1
                          wnin(mper) = xx
                          key(mper) = mper
                    endif
  100       continue
  101       continue
              call sort(wnin,key,mper)
C             if(ifrper.eq.0)then
C                   laxper = .true.
C             else
C                   laxper = .false.
C             endif
              close(2,status='keep')
        return
        end
        

        subroutine gcmdln(xlin,ylin,pmin,pmax,a0,a1,laxper,
     1    verbos,vmin,vmax,shadon,
     2      lorr,idva,iunit,
     2      unitcor, sacfile)
c-----
c     parse command line arguments and return control
c     parameters
c
c     requires subroutine mgtarg(i,name) to return
c            the i th argument in the string name
c
c     and the function mnmarg() to return the number
c            of arguments excluding program name
c            The first argument is i = 1
c
c-----
c     xlin      L      - if true linear frequency axis
c     ylin      L      - if true kinear velocity axis
c       pmin    R*4 - minimum period
c       pmax    R*4 - maximum period
c     a0         R*4    - Gaussian filter parameter
c     a1         R*4    - Gaussian filter parameter
c       iunit   I*4 Physical units of amplitude
c               -1 counts (default)
c               0 meters
c               1 centimeters
c               2 nanometers
c       idva    I*4 input time series (note output spectra will be 
c               spectra is displacement in cm, e.g., cm-sec of count-sec
c               0   displacment
c               1   velocity
c               2   acceleration
c       unitcor R*4 unit correction factor to get to cm
c     laxper       L      - Period if true, frequency if false
c     verbos     L      - verbose output
c     vmin       R*4    - minimum velocity for plot
c     vmax       R*4    - maximum velocity for plot
c       shadon  L   - true contour color shading
c     lorr  I*4 - 0 not defined, 1 Love, 2 Rayl
c       sacfile Ch*50   - name of sac file
c-----
        logical xlin, ylin, laxper, verbos
        real*4 a0,a1,vmin,vmax
        integer*4 idva, iunit
        real*4 unitcor
            integer *4 mnmarg
            character*50 name
            character*50 unitstr
            character sacfile*(*)
            logical shadon
            common/gtf/doabs,ampmx,doscl
            logical doabs,doscl

c-----
c     set up defaults for poor usage test
c-----
            xlin = .false.
            ylin = .true.
            laxper = .true.
            verbos = .false.
            shadon = .false.
            doabs = .false.
            doscl = .false.
            vmin = -1.0
            vmax = -1.0
            a0 = 50.27
            a1 = 0.0
            lorr = 0
        pmin = 4.0
        pmax = 60.0
c-----
c       first check the file mft96.cmd
c       if this file exists, read in pagameters for processing
c       then permit override by command line
c-----
C       inquire(file='mft96.cmd',exist=ext)
C       if(ext)then
C           call parcmd()
C       endif
c-----
c       read command line options for override
c-----
        nmarg = mnmarg()
        i = 0
   11   continue
            i = i + 1
            if(i.gt.nmarg)go to 13
            call mgtarg(i,name)
            if(name(1:5).eq.'-XLIN')then
                xlin = .true.
            else if(name(1:5).eq.'-XLOG')then
                xlin =.false.
            else if(name(1:2).eq.'-S')then
                shadon = .true.
            else if(name(1:2).eq.'-s')then
                doscl  = .true.
            else if(name(1:2).eq.'-R')then
                lorr = 2
            else if(name(1:2).eq.'-L')then
                lorr = 1
            else if(name(1:3).eq.'-a0')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,f10.0)')a0
            else if(name(1:3).eq.'-a1')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,f10.0)')a1
            else if(name(1:5).eq.'-PMIN')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,f10.0)')pmin
            else if(name(1:5).eq.'-PMAX')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,f10.0)')pmax
            else if(name(1:2).eq.'-V'.and.name(1:3).ne.'-VM')then
                verbos = .true.
            else if(name(1:5).eq.'-FREQ')then
                laxper = .false.
            else if(name(1:4).eq.'-PER')then
                laxper = .false.
            else if(name(1:5).eq.'-VMIN')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(f10.0)')vmin
            else if(name(1:5).eq.'-VMAX')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(f10.0)')vmax
            else if(name(1:2).eq.'-f')then
                i = i + 1
                call mgtarg(i,sacfile)
            else if(name(1:2).eq.'-A')then
                doabs = .true.
            else if(name(1:2).eq.'-U')then
                i = i + 1
                call mgtarg(i,unitstr)
            else if(name(1:2).eq.'-?'.or.name(1:2).eq.'-h')then
                call usage(0)
            endif
        go to 11
   13   continue
c-----
c       test for improper command line
c-----
        if(vmin.lt.0.0 .or.vmax.lt.0.0)then
                    vmin = 2.0
                    vmax = 5.0
        endif
c-----
c       parse units
c       
c-----
        call getunit(idva,iunit,unitcor,unitstr)
        return
        end

        subroutine usage(ierr)
        integer LOT
        parameter (LOT=6)
        integer ierr
c-----
c       prompt for proper usage
c-----
        if(ierr.eq.0)then
            write(LOT,*)'Usage: sacmft96 -XLIN -XLOG -S -s -R -L ',
     1          '-a0 a0 -a a1 -PMIN pmin -PMAX pmax ',
     2          '-V -FREQ -PER -VMIN vmin -VMAX vmax ',
     3          '-f sacfile -A -U units -? -h'
        write(LOT,*)
     1      '-XLIN      (default false)   Linear X-axis'
        write(LOT,*)
     1      '-XLOG      (default true)    Logarithmic X-axis'
        write(LOT,*)
     1      '-S         (default false)   Color dispersion shading'
        write(LOT,*)
     1      '-s         (default false)   Use linear color scale'
        write(LOT,*)
     1      '-R         (default unknown) Rayleigh Wave'
        write(LOT,*)
     1      '-L         (default unknown) Love Wave'
        write(LOT,*)
     1      '-a0        (default 50.27)   Filter parameter'
        write(LOT,*)
     1      '-a1        (default  0.00)   '
        write(LOT,*)
     1      '-PMIN pmin (default  4.00)   Minimum period '
        write(LOT,*)
     1      '-PMAX pmax (default 60.00)   Maximum period'
        write(LOT,*)
     1      '-V         (default false)   verbose'
        write(LOT,*)
     1      '-FREQ      (default false)   x-axis is frequency'
        write(LOT,*)
     1      '-PER       (default true)    x-axis is period'
        write(LOT,*)
     1      '-VMIN vmin (default 2.0 )    minimum velocity for plot'
        write(LOT,*)
     1      '-VMAX vmax (default 5.0 )    maximum velocity for plot'
        write(LOT,*)
     1      '-f sacfile (default none)    sacfile - required'
        write(LOT,*)
     1      '-A         (default false)   absolute amplitude contours'
        write(LOT,*)
     1      '-U units   (default counts)  cm cm/s cm/s/s m m/s',
     1      ' m/s/s ',
     2      'nm nm/s nm/s/s count'
        write(LOT,*)
     1      '-?         (default false)   usage'
        write(LOT,*)
     1      '-h         (default false)   usage'
        
        
        else if(ierr.eq.-1)then
              write(LOT,*)'sacmft96: unexpected EOF on input'
        else if(ierr.eq.-2)then
              write(LOT,*)'sacmft96: read error on input'
        else if(ierr.eq.1)then
              write(LOT,*)'sacmft96: N exceeds internal limits'
        endif
        stop
        end

        subroutine gettrc(xfile,n,n21,npts,dist,deg,az,baz,t0,dt,sta,
     1      comp,cdate,z,tarr,nptser,ierr,permin,permax)
c-----
c       get parameters required for processing
c-----
c       xfile   Ch* - name of data file
c       n   I*4 - array size, power of 2
c       n21 I*4 - n/2 + 1
c       npts    I*4 - original number of points (n >= npts)
c       dist    R*4 - epicentral distance (km)
c       deg R*4 - epicentral distance (deg)
c       az  R*4 - src -> rec azimuth
c       baz R*4 - rec -> src azimuth
c       t0  R*4 - time of first sample after origin time
c       dt  R*4 - sample interval
c       sta Ch*8    - station name
c       comp    Ch*8    - station component
c       cdate   Ch*12   - date string
c       z   C*4 - complex Fourier Transform array of samples
c       tarr    R*4 - temporary array  for trace
c       nptser  I*4 - length of z() and tarr() arrays
c       ierr    I*4 - error code
c       permin  R*4 - minimum period to be used in trace
c       permax  R*4 - maximum period to be used in trace
c-----
        character xfile*(*)
        integer n, n21, npts, nptser, ierr
        real dist, deg, az, baz, t0, dt
        character sta*8, comp*8, cdate*12
        complex z(nptser) 
        real tarr(nptser)
        real permin, permax
c-----
c       iunit   I*4 Physical units of amplitude
c               -1 counts (default)
c               0 meters
c               1 centimeters
c               2 nanometers
c       idva    I*4 input time series (note output spectra will be 
c               spectra is displacement in cm, e.g., cm-sec of count-sec
c               0   displacment
c               1   velocity
c               2   acceleration
c       unitcor R*4 unit correction factor to get to cm
c-----
        common/units/idva,iunit,unitcor
        integer*4 idva, iunit
        real*4 unitcor

        call getsac(xfile,npts,dist,deg,az,baz,t0,dt,sta,
     1                 comp,cdate,tarr,permin,permax)
        ls = lgstr(xfile)
        if(npts.gt.nptser)npts = nptser
        call npow2(npts,n,n21)
c-----
c       fill up the waveform array, use the opportunity to set
c       the units correctly
c-----
        do 100 i=1,n
            if(i.le.npts)then
                z(i) = cmplx(tarr(i)*unitcor,0.0)
            else
                z(i) = cmplx(0.0,0.0)
            endif
  100   continue
        call four(z,n,-1,dt,df)
C       write(6,*)'xfile,n,n21,npts,dist,deg,az,baz,t0,dt,sta',
C     1     xfile,n,n21,npts,dist,deg,az,baz,t0,dt,sta
        return
        end
c
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
        integer lgstr
*
        write(6,*)name
        call brsac (1,MXPTS,name,seis,ierr)
        ls = lgstr(name)
        write(6,*)'FILE: ',name(1:ls)
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
        if(permin .lt. 0.0 .or. permax.lt.0.0)then
            permin = 0.0
            permax = 100000.0
        endif
        write(6,'(a,a)')'    AZ     BAZ          DIST        DEG',
     1        '      DT      NPTS'
        write(6,'(f9.1,f9.1,f13.5,f10.2,f10.5,i7)')
     1      az,baz,dist,deg,dt,npts
        write(6,'(a,a)')
     1  '     EVLA        EVLO       STLA        STLO',
     2  '   KSTNM   KCMPNM'
        write(6,'(f11.4,f12.4,f11.4,f12.4,1x,a8,1x,a8)')
     1      evla,evlo,stla,stlo,kstnm,kcmpnm
*
        if(nerr .eq. 0 .and. origtime .ne. -12345)then
            t0 = beg - origtime
        else
            t0 = beg
        end if
*
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

        subroutine setxy(x,NX,y,NY,twn,nper,tfmn,tfmx,
     1      twmn,twmx,pumin,pumax,
     1      xlin,ylin,xlow,xhgh,ylow,yhgh,xaxlen,yaxlen)
c-----
c       set x,y absolute positions for cont(NX,NY) grid points
c-----
c       x(NX)   R*4 - x-positions
c       y(NY)   R*4 - y-positions
c       NX  I*4 
c       NY  I*4 
c       wn  R*4 - array of frequencies to the plotted
c       nper    I*4 - number of periond
c       tfmn    R*4 - minimum period in data set
c       tfmx    R*4 - maximum period in data set
c       twmn    R*4 - minimum frequency in data set
c       twmx    R*4 - maximum frequency in data set
c       pumin   R*4 - minimum velocity to be plotted
c       pumax   R*4 - maximum velocity to be plotted
c       xlin    L   - .true.  x-axis is linear
c                 .false. x-axis is logarithmic
c       ylin    L   - .true.  y-axis is linear
c                 .false. y-axis is logarithmic
c       xlow    R*4 - value of lower left corner
c       xhgh    R*4 - value of lower left corner
c       ylow    R*4 - value of upper right corner
c       yhgh    R*4 - value of upper right corner
c-----
c       xaxlen  R*4 - length of x-axis
c       yaxlen  R*4 - length of y-axis
c       laxper  L   - .true. x-axis is period
c                 .false. x-axis is frequency
c       fmin    R*4 - minimum frequency for x-axis
c       fmax    R*4 - maximum frequency for x-axis
c       wmin    R*4 - actual minimum frequency for x-axis
c       wmax    R*4 - actual maximum frequency for x-axis

c       wmin    R*4 - minimum value of plotted x-values
c       wmax    R*4 - maximum value of plotted x-values
c-----
              integer*4 nper
              logical ylin,xlin
        common/pltcnt/laxper
              logical laxper
        parameter (MAXFRQ=100)
        real*4 wn(MAXFRQ)
        real*4 twn(MAXFRQ)
c-----
c       information about plot on the page
c-----
        common/pltinf/gxl(3),gyl(3),gxh(3),gyh(3),
     1      axl(3),ayl(3),axh(3),ayh(3),ixlnlg(3),iylnlg(3)
            real*4 gxl,gyl,gxh,gyh,
     1              axl,ayl,axh,ayh
            integer*4 ixlnlg,iylnlg
c-----
c       information about units
c
c       iunit   I*4 Physical units of amplitude
c               -1 counts (default)
c               0 meters
c               1 centimeters
c               2 nanometers
c       idva    I*4 input time series (note output spectra will be 
c               spectra is displacement in cm, e.g., cm-sec of count-sec
c               0   displacment
c               1   velocity
c               2   acceleration
c       unitcor R*4 unit correction factor to get to cm
c-----
        common/units/idva,iunit,unitcor
        integer*4 idva, iunit
        real*4 unitcor
c-----
c       
c-----
        integer NX, NY
        real*4 x(NX), y(NY)
        character titlex*80, titley*80

c-----
c       set up plot parameters
c-----
c       define the array of increasing values
c-----
            do 900 i=1,nper
                if(laxper)then
                    
                    wn(nper+1-i) = 1.0/twn(i)
                else
                    wn(i) =     twn(i)
                endif
  900       continue
        titley = 'Group Velocity (km/s)'
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
C       WRITE(6,*)'fmin,fmax,wmin,wmax:',fmin,fmax,wmin,wmax
        if(xlin)then
            if(fmax.lt.wmax)fmax = wmax
            if(fmin .lt.0.0 .or. fmin.gt.wmin)then
                x1 = wmin
            else
                x1 = fmin
            endif
            xlow = x1
        XLOW = wmin
            xhgh = wmax
            do 1000 i=1,nper
            x(i) = 0.0 + xaxlen*(wn(i)-wn(1))/(wn(nper)-wn(1))
 1000       continue
        else
            xlow = wmin
            xhgh = wmax
            do 1001 i=1,nper
                x(i) = xaxlen*alog10(wn(i)/wn(1))/
     1              alog10(wn(nper)/wn(1))
 1001       continue
        endif
C       WRITE(6,*)'xlin  :',xlin
C       WRITE(6,*)'laxper:',laxper
C       WRITE(6,*)'x(1),wn(1),x(nper),wn(nper):',x(1),wn(1),x(nper),wn(nper)
C       WRITE(6,*)'x1,wmin,fmin               :',x1,wmin,fmin
C       WRITE(6,*)'tfmn,tfmx,twmn,twmx        :',tfmn,tfmx,twmn,twmx
c-----
c       set up y-axis
c-----
        if(ylin)then
            ylow = pumin
            yhgh = pumax
            do 1003 i=1,NY
                y(i) =0.0 +  (i-1)*yaxlen/(NY-1)
 1003       continue
        else
            ylow = pumin
            yhgh = pumax
            dvel = (pumax - pumin)/(NY -1 )
            do 1004 i=1,NY
                vel =  pumin + (i-1)*dvel
                y(i) =  yaxlen*alog10(vel/pumin)/
     1                  alog10(pumax/pumin)
 1004       continue
        endif
c-----
c       plot the axes
c-----
        lx = lgstr(titlex)
        ly = lgstr(titley)
        call xyaxes(xlow,xhgh,xlin,ylow,yhgh,ylin,nscalx,nscaly,
     1          xaxlen,yaxlen,titlex,titley,lx,ly)
        return
        end

        subroutine xyaxes(xlow,xhgh,xlin,ylow,yhgh,ylin,nscalx,nscaly,
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
c       nscalx  I*4 10** scaling factgor for x-axis
c       nscaly  I*4 10** scaling factgor for y-axis
c       xaxlen  R*4 length of x-axis
c       yaxlen  R*4 length of y-axis
c       titlex  Ch  x-axis title string
c       titley  Ch  y-axis title string
c-----
        real*4 xlow, xhgh, ylow, yhgh
        logical xlin, ylin
        integer*4 nscalx, nscaly
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

        subroutine mkfreq(mper,nper,wn,dt,n,wmax,wmin,df)
        parameter (MAXFRQ=100)
        common /fval/ wo, noer, wnin, moer
            real*4 wo(MAXFRQ), wnin(MAXFRQ)
            integer*4 noer, moer
c-----
c       create the center frequencies for filters
c-----
        integer mper, nper, n
        real dt, wmax, wmin, df
        real wn(MAXFRQ)

              fnyq = 1./(2.*dt)
              wmin = 1.0e+38
              wmax = -1.0e+38
              nper = 0
              df = 1./(n*dt)
              DO 5002 II=1,mper
                    wnn = wnin(ii )
                    if(wnn.le.fnyq .and. wnn.gt. 0.0)then
                          nper=nper + 1
                          wn(nper) = wnn
                          if(wn(nper).gt.wmax)wmax=wn(nper)
                          if(wn(nper).lt.wmin)wmin=wn(nper)
                    endif
 5002       continue
        return
        end

        subroutine mapgv(j,nper,x,dxdt,xfrq,n,cont,NX,NY,nsamp,
     1          t0,dist,az,dt,pumin,pumax,laxper,w0,alpha)
c-----
c       map the group velocities uniformly onto the cont(i,j) array
c       where i=frequency and j = velocity = we map with 1,1 at
c       highest frequency and highest velocity
c-----
c       j   I*4 - frequency index
c       nper    I*4 - number of unique filter frequencies
c       x   R*4 - array of amplitudes
c       xfrq    R*4 - array of instantaneous frequencies
c       n   I*4 - number of points in time series
c       cont    R*4 - NX,NY array of values to be contoured
c       NX  I*4 - array dimension
c       NY  I*4 - array dimension
c       nsamp   I*4 - number of samples in original time series <=n
c       t0  R*4 - time for first sample
c       dist    R*4 - distance
c       az  R*4 - source-receiver azimuth
c       dt  R*4 - sample interval
c       pumin   R*4 - minimum velocity for plot
c       pumax   R*4 - maximum velocity for plot
c       laxper  L   - .true. x-axis is period
c                 .false. x-axis is frequency
c       w0  R*4 - filter frequency
c       alpha   R*4 - Gaussian filter parameter used
c-----
        real x(n), xfrq(n)
        real dxdt(n)
        real cont(NX,NY)
        logical laxper
        integer LIN, LOT, LER
        parameter (LIN=5, LER=0, LOT=6)
        parameter (MAXFRQ=100, MAXVEL=10)
        common/group/gpvel(MAXVEL,MAXFRQ), gpamp(MAXVEL,MAXFRQ),
     1      numgp(MAXFRQ)
            real gpvel, gpamp
            integer numgp
        common/discnt/kount, lorr
            integer *4 kount, lorr
        real evla, evlo, evdp, stla, stlo, stel, beg, origtime
        integer ntimes(6)
        common /shdr/ evla,evlo,evdp,stla,stlo,stel,
     &               beg,ntimes,origtime
c-----
c       u   R*4 array of group velocities
c       a   R*4 array of spectral amplitudes
c               ordered from hightest to smallest
c-----
        real*4 u(MAXVEL),a(MAXVEL),ufr(MAXVEL)
        integer*4 MAX4, kk
        integer *4 nzyear, nzjday, nzhour,nzmin
        character*8 kstnm, kcmpnm
c-----
c       get limit indices for pointers into time series
c-----
c       time = t0 + (i-1)*dt
c
c-----
        do 1001 i=1,nsamp-1
            dxdt (i) = x(i+1) - x(i)
 1001   continue
        MAX4 = MAXVEL
        per = 1.0/w0
        if(laxper)then
            wout = per
        else
            wout = w0
        endif
        dvel = (pumax - pumin)/(NY -1 )
C       WRITE(6,*)j,NX,NY,nsamp,t0,dist,dt,pumin,pumax,dvel,laxper,w0
        do 1000 i=1,NY
            velmx =  pumax - (i-1)*dvel
            ii = ((dist/velmx) - t0)/dt + 1
            if(ii.lt.1)go to 1000
            if(ii.gt.n)go to 1000
C           velmn =  pumin - (i-1)*dvel
C           ij = ((dist/velmn) - t0)/dt + 1
C           if(ij.lt.1)go to 1000
C           if(ij.gt.n)go to 1000
c-----
c               for PERIOD    cont(NX-j+1,NY-i+1) = x(ii)
c               for FREQUENCY cont(   j  ,NY-i+1) = x(ii)
c-----
            if(laxper)then
                cont(nper-j+1,NY-i+1) = x(ii)
            else
                cont(   j  ,NY-i+1) = x(ii)
            endif
 1000   continue
c-----
c       obtain the MAX4 largest peaks for each filter periode
c-----
        kk = 0
c-----
c       use the dxdt to search for extrema
c-----
c-----
c       estimate the starting and ending values of the loop
c       according to the group velocity window
c-----
        tuend = dist/pumin - t0
        iuend = tuend/dt + 1
        if(iuend.gt.nsamp)iuend = nsamp
        if(iuend.lt.1)iuend = 1
        tustrt = dist/pumax - t0
        iustrt = tustrt/dt + 1
        if(iustrt.gt.nsamp)iustrt = nsamp
        if(iustrt.lt.1)iustrt = 1
C       do 400 i=2,nsamp-1
        do 400 i=iustrt+1,iuend-1
            if(sign(1.0,dxdt(i-1)).ne.sign(1.0,dxdt(i)))then
c-----
c       prepare to interpolate
c-----
                p = (dxdt(i-1) - 0.0)/(dxdt(i-1) - dxdt(i))
                time = t0 + (i-2 +p) * dt
                if(time.gt.0.0)then
                    vel = dist/time
                    freq = p*xfrq(i-1) +(1-p)*xfrq(i)
                    ampl = p*x(i-1) + (1-p)*x(i)
                    kk = kk + 1
                    call srt(a,u,ufr,vel,ampl,freq,kk,MAX4)
                endif
            endif
 400    continue
CRBH- commented out 23 AUG 2003
CRBH    amp1 = x(1)
CRBH    amp2 = x(2)
CRBH    frq1 = xfrq(1)
CRBH    frq2 = xfrq(2)
CRBH    nm1 = nsamp - 1
CRBH    do 500 i=3,nm1
CRBH        amp3 = x(i)
CRBH        frq3 = xfrq(i)
CRBH        if( (amp2-amp1).gt.1.0e-7*amp1 .and.
CRBH     1          (amp2-amp3).ge. 1.0e-7*amp3)then
CRBH            amp4 = x(i+1)
CRBH            frq4 = xfrq(i+1)
CRBH        if((amp2-amp4).ge.1.0e-7*amp4) then
CRBH            time = t0 + (i-2) * dt
CRBH            if(time.gt.0.0)then
CRBH                vel = dist/time
CRBH        if(vel.ge.pumin .and .vel.le.pumax)then
CRBH                kk = kk + 1
CRBH                ampl = amp2 
CRBH                freq = frq2
CRBH                call srt(a,u,ufr,vel,ampl,frq2,kk,MAX4)
CRBH        endif
CRBH            endif
CRBH        endif
CRBH        endif
CRBH    amp1 = amp2
CRBH    amp2 = amp3
CRBH    frq1 = frq2
CRBH    frq2 = frq3
CRBH  500   continue
        if(kk.gt.MAX4)then
            kk = MAX4
        endif
        do 501 i=1,kk
            gpvel(i,j) = u(i)
            gpamp(i,j) = a(i)
  501   continue
        numgp(j) = kk
c-----
c       output the values - note should save for plot
c       should have a integer array for how many peaks per frequency
c       2d array for amp(frequency,peak) velocity (frequency,peak\)
c       output period always for 
c-----
        call getnhv('NZYEAR',nzyear,ierr)
        call getnhv('NZJDAY',nzjday,ierr)
        call getnhv('NZHOUR',nzhour,ierr)
        call getnhv('NZMIN ',nzmin ,ierr)
        call getkhv('KSTNM ',kstnm ,ierr)
        call getkhv('KCMPNM',kcmpnm,ierr)
c-----
c       safety check on the header values
c-----
        if(nzyear .lt. 0 .or. nzjday .lt. 0 .or. nzhour .lt. 0
     1      .or. nzmin .lt. 0)then
            nzyear = 1970
            nzjday =    1
            nzhour =    0
            nzmin  =    0
        endif

        DO 5008 i=1,kk
        if(a(i).gt.0.0 .and.u(i).ge.pumin.and.u(i).le.pumax)then
            time = dist/u(i)
            err  = u(i)*(per/time)
            ufr(i) = ufr(i)/6.2831853
            if(ufr(i).ge.0.0 )then
                uper = 1.0/ufr(i)
            else
                uper = -1.0
            endif
            if(evla .eq. -12345.0 .or. evlo.eq. -12345.0 .or.
     1          stla .eq. -12345.0 .or. stlo .eq. -12345.0)then
                    eevla = 0.0
                    sstla = dist/111.195
                    eevlo = 0.0
                    sstlo = 0.0
            else
                eevla = evla
                eevlo = evlo
                sstla = stla
                sstlo = stlo
            endif
                
            if(lorr.eq.0)then
                write(8,10)
     1          per,u(i),err,dist,az,a(i),eevla,eevlo,
     1          sstla,sstlo,i,uper,alpha,
     2          kstnm, kcmpnm, nzyear, nzjday, nzhour, nzmin
            else if(lorr.eq.1)then
                write(8,11)
     1          per,u(i),err,dist,az,a(i),eevla,eevlo,
     1          sstla,sstlo,i,uper,alpha,
     2          kstnm, kcmpnm, nzyear, nzjday, nzhour, nzmin
            else if(lorr.eq.2)then
                write(8,12)
     1          per,u(i),err,dist,az,a(i),eevla,eevlo,
     1          sstla,sstlo,i,uper,alpha,
     2          kstnm, kcmpnm, nzyear, nzjday, nzhour, nzmin
            endif
            kount = kount + 1
        endif
 5008   continue
c-----
c      plot the four largest amplitudes
c-----
        if(lorr.eq.0)then
            ifmt=10
        else if(lorr.eq.1)then
            ifmt=11
        else if(lorr.eq.2)then
            ifmt=12
        endif
   10   format('MFT96 A U -1 ',4g11.5,f9.1,e11.4,4f12.6,
     1   ' 0',i3,g12.4,g12.4,' COMMENT: ',2a8,i5,i4,3i3)
   11   format('MFT96 L U -1 ',4g11.5,f9.1,e11.4,4f12.6,
     1   ' 0',i3,g12.4,g12.4,' COMMENT: ',2a8,i5,i4,3i3)
   12   format('MFT96 R U -1 ',4f11.5,f9.1,e11.4,4f12.6,
     1   ' 0',i3,g12.4,g12.4,' COMMENT: ',2a8,i5,i4,3i3)
        return
        end

        subroutine pltgv(x,nper,ylin,yhgh,ylow,yaxlen,laxper)
c-----
c       map group velocity value into plot coordinate
c       note we already know the x coordinate because of
c       the filter frequencies
c-----  
c       x   R*4 array of abscissa
c       nper    I*4 number of unique filter frequencies
c       ylin    L   .true. y-axis is linear
c       ylow    R*4 lowest value of y-axis
c       yhgh    R*4 highest value of y-axis
c       yaxlen  R*4 length of y-axis
c       laxper  L   .true. x-axis is period
c-----
        parameter(NX=100)
        real x(NX)
        integer*4 nper
        logical ylin
        real*4 ylow,yhgh
        real*4 yaxlen
        logical laxper

        parameter (MAXFRQ=100, MAXVEL=10)
        common/group/gpvel(MAXVEL,MAXFRQ), gpamp(MAXVEL,MAXFRQ),
     1      numgp(MAXFRQ)
            real gpvel, gpamp
            integer numgp

        do 1000 j=1,MAXFRQ
            if(laxper)then
                xx = x(nper+1-j)
            else
                xx = x(j)
            endif
            do 1100 i=1,numgp(j)
                u = gpvel(i,j)
                if(u.ge.ylow .and. u.le.yhgh)then
                if(ylin)then
                    yy = 0.0 + yaxlen*(u-ylow)/(yhgh-ylow)
                else
                    yy = 0.0 + yaxlen*alog10(u/ylow)
     1                  /alog10(yhgh/ylow)
                endif
                call gsolid(xx,yy,0.03,i-1)
                endif
 1100       continue
 1000   continue
        return
        end
        
        subroutine pltamp(wn,nper,xlin,ylin,xaxlen,yaxlen,laxper,
     1      pumin,pumax)
c-----
c       map group velocity value into plot coordinate
c       note we already know the x coordinate because of
c       the filter frequencies
c-----  
c       wn  R*4 array of filter frequencies
c       nper    I*4 number of unique filter frequencies
c       xlin    L   .true. x-axis is linear
c       ylin    L   .true. y-axis is linear
c       xaxlen  R*4 length of x-axis
c       yaxlen  R*4 length of y-axis
c       laxper  L   .true. x-axis is period
c       pumin  R*4      minimum group velocity for plots
c       pumax  R*4      maximum group velocity for plots
c-----
        integer*4 nper
        logical xlin, ylin
        real*4 xaxlen, yaxlen
        logical laxper
        real pumin, pumax

        parameter (MAXFRQ=100, MAXVEL=10)
        real wn(MAXFRQ)

        common/group/gpvel(MAXVEL,MAXFRQ), gpamp(MAXVEL,MAXFRQ),
     1      numgp(MAXFRQ)
            real gpvel, gpamp
            integer numgp
c-----
c       information about plot on the page
c-----
        common/pltinf/gxl(3),gyl(3),gxh(3),gyh(3),
     1      axl(3),ayl(3),axh(3),ayh(3),ixlnlg(3),iylnlg(3)
            real*4 gxl,gyl,gxh,gyh,
     1              axl,ayl,axh,ayh
            integer*4 ixlnlg,iylnlg

        common/units/idva,iunit,unitcor
        integer*4 idva, iunit
        real*4 unitcor

        character titlex*80, titley*80
c-----
c       set up unit string
c-----
        if(iunit.lt.0)then
            titley = 'count-sec'
        else
            titley = 'cm-sec'
        endif
c-----
c       get bounds for plot
c-----
c       first search for maximum amplitude
c-----
        ymax = 0.0
        do 1000 j=1,MAXFRQ
            do 1100 i=1,numgp(j)
                amp = gpamp(i,j)
                vel = gpvel(i,j)
                if(vel.ge.pumin .and. vel.le.pumax)then
                    if(amp.gt.ymax)ymax = amp
                endif
 1100       continue
 1000   continue
c-----
c       safety
c-----
        if(ymax.eq.0.0)ymax = 1.0
c-----
c       get frequency/period limits
c-----
        if(laxper)then
            xlow = 1.0/wn(nper)
            xhgh = 1.0/wn(1)
            titlex = 'Period (sec)'
        else
            xlow = wn(1)
            xhgh = wn(nper)
            titlex = 'Frequency (Hz)'
        endif
        if(xhgh.lt.100.0*xlow)then
            rat = 100.0/(xhgh/xlow)
            rat = sqrt(rat)
            xlow = xlow/rat
            xhgh = xhgh*rat
        endif 
        ymax = alog10(ymax)
        lmax = ymax
        if(ymax.gt.real(lmax))lmax = ymax + 1
        yhgh = 10.0**lmax
        ylow = yhgh/1000.0
c-----
c       start plot
c       positioning for spectral amplitude plot
c-----
        xoff = 1.3 
        yoff = 1.5
        call plot( xoff, yoff,-3)
c-----
c       plot axes
c-----
        lx = lgstr(titlex)
        ly = lgstr(titley)
        call xyaxes(xlow,xhgh,.false.,ylow,yhgh,.false.,nscalx,nscaly,
     1          xaxlen,yaxlen,titlex, titley,lx,ly)
c-----
c       plot observed points
c-----
        do 2000 j=1,MAXFRQ
            if(laxper)then
                xv = 1.0/wn(j)
            else
                xv = wn(j)
            endif
            do 2100 i=1,numgp(j)
                amp = gpamp(i,j)
                if(amp.le.yhgh .and. amp.ge.ylow)then
                    xx = xaxlen*(alog10(xv/xlow)/
     1                  alog10(xhgh/xlow))
                    yy = yaxlen*(alog10(amp/ylow)/
     1                  alog10(yhgh/ylow))
                call gsolid(xx,yy,0.03,i-1)
                endif
 2100       continue
 2000   continue
        call plot(-xoff,-yoff,-3)
c-----
c       describe plot positioning
c-----
        gxl(1) = xoff
        gyl(1) = yoff
        gxh(1) = xoff + xaxlen
        gyh(1) = yoff + yaxlen
        axl(1) = xlow
        axh(1) = xhgh
        ayl(1) = ylow
        ayh(1) = yhgh
        ixlnlg(1) = 2
        iylnlg(1) = 2
        return
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
        integer MAXFRQ
        parameter (MAXFRQ=100)
        common /fval/ wo, noer, wnin, mper
            real*4 wo(MAXFRQ), wnin(MAXFRQ)
            integer*4 noer, mper
        integer NX, nper
        real wn(NX)
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
                if(p .ge. 0.99*pmin .and. p .le. 1.01*pmax
     1          .and. p.ge. 0.99*permin .and. p.le.1.01*permax)then
                    nper = nper + 1
                    mper = nper
                    wn(nper) = 1.0/p
                          wnin(mper) = 1.0/p
                          key(mper) = mper
                    if(nper.eq.NX)go to 1000
                endif
 200        continue
 100    continue
 1000   continue
C       write(6,*)(wnin(i),i=1,mper)
              call sort(wnin,key,mper)
C       write(6,*)(wnin(i),i=1,mper)
        return
        end

        subroutine getunit(idva,iunit,unitcor,unitstr)
c-----
c       parse the unit string to define the units of the time series
c       to match observations to synthetics, we would like to
c       convert the input time series to centimeters
c
c       this parsing is crude, just going from the longest to the
c       shortest posswible unique string
c-----
c       iunit   I*4 Physical units of amplitude
c               -1 counts (default)
c               0 meters
c               1 centimeters
c               2 nanometers
c       idva    I*4 input time series (note output spectra will be 
c               spectra is displacement in cm, e.g., cm-sec of count-sec
c               0   displacment
c               1   velocity
c               2   acceleration
c       unitstr CH*(*)  string with units
c-----
        integer*4 idva, iunit
        character unitstr*(*)
c-----
c       internal variables
c       ioff    I*4 offset into string for units
c       ls  I*4 length of string
c       unitcor R*4 unit correction factor to get to cm
c-----
        ls = lgstr(unitstr)
        if(unitstr(1:1).eq.'m')then
                iunit = 0
                ioff = 2
                unitcor = 100.0
        else if(unitstr(1:2).eq.'cm')then
                iunit = 1
                ioff = 3
                unitcor = 1.0
        else if(unitstr(1:2).eq.'nm')then
                iunit = 2
                ioff = 3
                unitcor = 1.0e-7
        else if(unitstr(1:2).eq.'co')then
                iunit = -1
                ioff = 3
                unitcor = 1.0
        else
                iunit = -1
                ioff = 3
                unitcor = 1.0
        endif
        idva = 0
        if(unitstr(ioff:ls).eq.'/s/s')then
                idva = 2
        else if(unitstr(ioff:ls).eq.'/s')then
                idva = 1
        endif
        return
        end

        subroutine dofilt(datas,data,datad,dt,df,
     1          n,n21,alpha,idva,ampmx,wn,pumin,pumax,
     2          dist,t0,nsamp)
c-----
c       perform narrow band pass filter
c
c       datas   C   Complex array of saved spectrum
c       data    C   Complex array of filtered signal
c       datad   C   Complex array of derivative of filtered signal
c       dt  R   sampling interval, seconds
c       df  R   frequency sampling interval
c       n   I   number of points in time series (power of 2)
c       n21 I   n/2 + 1
c       alpha   R   Gaussian filter parameter
c       idva    I   0 displacment time series
c               1 velocity time series
c               2 acceleration time series
c               we want to have displacement spectrum
c       ampmx   R   Maximum amplitude for this filter frequency
c       wn  R   Filter frequency
c       pumin  R*4      minimum group velocity for plots
c       pumax  R*4      maximum group velocity for plots
c       dist    R   epicentral distance in km
c       t0  R   travel time of first sample
c       nsamp   I number of samples in original time series
c-----
        integer NPTSER, NPTFRQ
        parameter (NPTSER=132000, NPTFRQ=65000)
        complex datas(NPTFRQ),data(NPTSER),datad(NPTSER)
        complex s, s2

        integer nsamp
        real pumin, pumax
c-----
c       define cutoff for filter
c-----
            fac = sqrt(3.1415927/alpha)
            frequp = (1.0+fac)*wn
            freqlw = (1.0-fac)*wn
            if(freqlw .le. 0.0)freqlw = df
            do 1002 i=1,n21
                xi = i - 1
                freq = xi *df
                if(freq.ge.freqlw .and. freq.le.frequp)then
                    fact = -alpha*((freq-wn)/
     1                              (wn))**2
                    filt = exp(fact)
                    s = cmplx(0.0, 6.2831853*freq)
                    s2 = s * s
c-----
c       convert from velocity to displacement, or from
c       acceleration to displacement here
c       since we are bandpass filtgering we never divide by zero??
c-----
                    data(i) = filt*datas(i)
                    if(idva.eq.1)then
                        data(i) = data(i)/s
                    else if(idva.eq.2)then
                        data(i) = data(i)/s2
                    endif
c------
c       set up array to compute time derivative which is
c       necessary for computing the instantaneous period/frequency
c-----
                    datad(i) = data(i) * s
                else
                    data(i) = cmplx(0.0,0.0)
                    datad(i) = cmplx(0.0,0.0)
                endif
c-----
c           zero out negative frequencies
c-----
                if(i.gt.1)then
                    data(n+2-i)=cmplx(0.0,0.0)
                    datad(n+2-i)=cmplx(0.0,0.0)
                endif
 1002       continue
            call four(data,n,+1,dt,df)
            call four(datad,n,+1,dt,df)
c-----
c           save envelope for contour plot
c           wn is the filter frequency
c-----
            afac = sqrt(alpha/3.1415927)/wn
            amx = 0.0
c-----
c       compute the spectral amplitude and instantaneous
c       frequency and store into first half of the
c       complex data and datad arrays
c       Note we want the maximum amplitude within a given group velocity window
c-----
        tuend = dist/pumin - t0
        iuend = tuend/dt + 1
        if(iuend.gt.nsamp)iuend = nsamp
        if(iuend.lt.1)iuend = 1
        tustrt = dist/pumax - t0
        iustrt = tustrt/dt + 1
        if(iustrt.gt.nsamp)iustrt = nsamp
        if(iustrt.lt.1)iuend = 1
            do 5006 KK = 1,n
                xamp = cabs(data(KK))*afac
                if(xamp.gt.0)then
                    xfrq = aimag(datad(KK)/data(KK))
                else
                    xfrq = 1.0e+5
                endif
                if(KK.ge.iustrt .and. KK.le. iuend)then
                    if(xamp.gt.amx)amx=xamp
                endif
                data(KK)=cmplx(xamp,0.0)
                datad(KK)=cmplx(xfrq,0.0)
 5006       continue
c-----
c       Now systematically pack the complex array
c       in a manner that is the same as the equivalence
c-----
            KJ = 1
            do 5007 KK=1,n
                xfrq = real(datad(KK))
                xamp = real(data (KK))
                if(mod(KK,2).eq.1)then
                    data (KJ) = cmplx(xamp,0.0)
                    datad(KJ) = cmplx(xfrq,0.0)
                else
                    data (KJ) = cmplx(real(data (KJ)),xamp)
                    datad(KJ) = cmplx(real(datad(KJ)),xfrq)
                    KJ = KJ + 1
                endif
                
 5007       continue
            if(amx.gt.ampmx)ampmx = amx
            return 
            end
