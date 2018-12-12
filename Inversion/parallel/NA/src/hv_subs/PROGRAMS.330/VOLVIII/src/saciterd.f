        program saciterd
c---------------------------------------------------------------------c
c                                                                    c
c     COMPUTER PROGRAMS IN SEISMOLOGY                                c
c     VOLUME V                                                       c
c                                                                    c
c     PROGRAM: SACITERD                                              c
c                                                                    c
c     COPYRIGHT 1996 1998 CHARELS J AMMON                            c
c                                                                    c
c     Department of Earth and Atmospheric Sciences                   c
c     Saint Louis University                                         c
c     221 North Grand Boulevard                                      c
c     St. Louis, Missouri 63103                                      c
c     U. S. A.                                                       c
c     VERSION 1.04                                                   c
c                                                                    c
c---------------------------------------------------------------------c
c      Changes:
c
c      22 JAN 2002 - add -RAYP rayp to the command line so that
c          the ray parameter is placed in the USER4 SAC variable
c      21 NOV 2002 - modified KEVNM header as per 
c          Mejian An Dept Geophysics, IAG, Sao Paulo, Br   
c      09 JAN 2005 - nbumps undefined in  subroutine tdomain, changed
c          to nshift. tshift undefined 
c          in subroutine tdomain  changed the theshift, baker@usgs.gov
c      23 APR 2007 - repaired the doubling construct at line  line 490
c          which was incorrectly implemented on early 2007 - the code now
c          recreates the iamges of 2006-11-10
c      24 APR 2009 - permit ray parameter in E format
c       23 JUL 2009 - ensured consistent call to FFT routine with
c               complex arguments. We now call zfour
c-----
c
c    Based on the Kikuchi and Kanamori (1982) iterative
c     deconvolution algorithm. The output of this code
c     is the deconvolution of the "denominator" from
c     the "numerator"
c    Kikuchi, M., and H. Kanamori (1982). Inversion of complex body waves, 
c    Bull. Seism. Soc. Am. 72, 491-506.
c
c    Ligorria, J. P., and C. A. Ammon (1999). Iterative deconvolution and 
c    receiver-function estimation,  Bull. Seism. Soc. Am. 89, 1395-1400.
c
c    The final deconvolution is called "decon.out"
c      the observed time series (original convolved with Gaussian)
c       is called "observed"
c      the predicted is call "predicted"
c
c    Header values from the "numerator" are copied into the
c      decon.out file. The gwidth is stored in user0.
c
c-----
c  INTERNAL PARAMETERS
c
c  rhdr    = SAC waveform real header fields (70).
c  ihdr    = SAC waveform integer header fields (40).
c  chdr    = SAC waveform character header fields (48).
c  x(i)    = array containing the waveform data.
c  delta   = waveform sampling interval in seconds.
c  btime   = waveform begin time.
c-----
        common/sachdr/rhdr,ihdr,chdr
        real*4 rhdr(70)
        integer*4 ihdr(40)
        character*8 chdr(24)
c-----
c      command line arguments
c-----
        character*80 fnum, fden
        real gwidth, theshift, tol
        logical lpositive, sacbin, verbose, dotwice
        integer maxbumps
c-----
c      internal variables
c-----
        parameter (LER=6, LIN=5, LOT=6)
        parameter (NSAMP=131072)
        real f(NSAMP),g(NSAMP), p(NSAMP), r(NSAMP)
        common/zval/z1, z2
        complex z1(NSAMP), z2(NSAMP)
c-----
c      call machine dependent initialization
c-----
        call mchdep()
c-----
c      Author's credit
c-----
        write(LOT,*) ' '
        write(LOT,*) 'Program iterdeconfd - Version 1.04, 1997-98'
        write(LOT,*) 'Chuck Ammon, Saint Louis University'
        write(LOT,*) ' '
c-----
c      parse command line arguments
c-----
        call gcmdln(fnum,fden,maxbumps,theshift,tol,lpositive,gwidth,
     1      sacbin,verbose,dotwice,rayp)

        ls = lgstr(fnum)
        write(LOT,*)'Numerator   file:',fnum(1:ls)
        ls = lgstr(fden)
        write(LOT,*)'Denominator file:',fden(1:ls)
        write(LOT,*)'Maxbumps        :',maxbumps
        write(LOT,*)'Time shift      :',theshift
        write(LOT,*)'Tolerance       :',tol
        write(LOT,*)'Postive Amp Only:',lpositive
        write(LOT,*)'Filter parameter:',gwidth
        write(LOT,*)'Verbose         :',verbose
        write(LOT,*)'SAC binary files:',sacbin
c-----
c      do the decon
c-----
        call tdomain(fnum,fden,f,g,p,r,nn,npts,dt,gwidth,theshift,
     1      verbose,maxbumps,lpositive,tol,sacbin,dotwice,rayp)
        end

        subroutine gcmdln(fnum,fden,niter,delay,error,
     1      dopos,alpha,sacbin,verbose,dotwice,rayp)
c-----
c      parse command line arguments
c
c      requires subroutine targ() and funtion mnmarg()
c
c-----
c-----
c      fnum    Ch* - file name for numerator
c      fden    Ch* - file name for denominator
c      niter   I*4 - numebr of iterations
c      delay   R*4 - time shift (seconds before P)
c      error   R*4 - convergence criteria
c      dopos   L   - .true. permit only positive bumps
c      alpha   R   - Gaussian filter parameter
c      sacbin  L   - .true. if SAC binary output, else ascii
c      verbose L   - .true. output many intermediate files
c      dotwice L   - .true. use double length FFT to avoid wrap around
c      rayp    R   - ray parameter in sec/km, default of -12345.
c-----
        
        character*80 name
        character fnum*(*), fden*(*)
        logical dopos, sacbin, verbose, dotwice
        real alpha, error, delay, rayp
        integer niter
        integer mnmarg, i
        

        nmarg = mnmarg()
        i = 0
        nfile = 0
        alpha = 1.0
        niter = 100
        error = 0.001
        dopos = .false.
        sacbin = .true.
        verbose = .false.
        dotwice = .false.
        delay = 5.0
        rayp = -12345.
        fnum = ' '
        fden = ' '
   11   i = i + 1
        if(i.gt.nmarg)go to 13
            call mgtarg(i,name)
            if(name(1:3).eq.'-FN'.or.
     1          name(1:3).eq.'-fn')then
                i=i+1
                call mgtarg(i,fnum)
                nfile = nfile + 1
            else if(name(1:3).eq.'-FD'.or.
     1          name(1:3).eq.'-fd')then
                i=i+1
                call mgtarg(i,fden)
                nfile = nfile + 1
            else if(name(1:4).eq.'-ALP')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,f10.0)')alpha
            else if(name(1:2).eq.'-E')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,f10.0)')error
            else if(name(1:2).eq.'-D')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,f10.0)')delay
            else if(name(1:2).eq.'-N')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,i10)')niter
            else if(name(1:4).eq.'-POS')then
                dopos = .true.
            else if(name(1:2).eq.'-V')then
                verbose = .true.
C           else if(name(1:2).eq.'-B')then
C               sacbin = .true.
C           else if(name(1:2).eq.'-A' .and. name(1:3).ne.'-ALP')then
C               sacbin = .false.
            else if(name(1:2).eq.'-2')then
                dotwice = .true.
            else if(name(1:5).eq.'-RAYP')then
                i = i + 1
                call mgtarg(i,name)
                call chtofp(name,rayp)
            else if(name(1:2).eq.'-?')then
                call usage()
            else if(name(1:2).eq.'-h')then
                call usage()
            endif
        go to 11
   13   continue
        if(nfile.ne.2)call usage()
        if(fnum .eq. ' ' .or. fden .eq. ' ')call usage()
        return
        end

        subroutine usage()
        integer*4 LER, LIN, LOT
        parameter(LER=6, LIN=5, LOT=6)
        write(LER,*)'Usage: saciterd  ',
     1      ' -FN file_num -FD file_den [-E error] [-N niter]',
     2      ' [-POS] [-D delay] [-ALP alpha] [-2]',
     3      ' [-RAYP rayp]  [-?] [-h]'
        write(LER,*)
     1      ' -FN file_num (default none) numerator SAC binary'
        write(LER,*)
     1      ' -FD file_den (default none) denominator SAC binary'
        write(LER,*)
     1      ' -E error (0.001) convergence criteria'
        write(LER,*)
     1      ' -ALP alpha (default 1.0) Gaussian Filter Width'
        write(LER,*)
     1      '      H(f) = exp( - (pi freq/alpha)**2) '
        write(LER,*)
     1      '      Filter corner ~ alpha/pi '
        write(LER,*)
     1      ' -N niter (default 100, 1000 max) Number iterations/bumps'
        write(LER,*)
     1      ' -D delay (5 sec) Begin output delay sec before t=0'
        write(LER,*)
     1      ' -POS (default false) Only permit positive amplitudes'
C       write(LER,*)
C     1     ' -A  (default false) data are SAC ascii'
C       write(LER,*)
C     1     ' -B  (default true) data are SAC binary'
        write(LER,*)
     1      ' -2  (default false) use double length FFT to ',
     2      ' avoid FFT wrap around in convolution'
        write(LER,*)
     1      ' -RAYP rayp (default -12345.0) ray parameter in ',
     2      ' sec/km used by rftn96/joint96'
        write(LER,*)
     1      ' -?  This help message'
        write(LER,*)
     1      ' -h  This help message'
        write(LER,*)'Output files:'
        write(LER,*)'    observed: original numerator convolved ',
     1          'with Gaussian'
        write(LER,*)'   numerator: original numerator convolved ',
     1          'with Gaussian'
        write(LER,*)' denominator: original numerator convolved ',
     1          'with Gaussian'
        write(LER,*)'   decon.out: Receiver function for Gaussian'
        write(LER,*)'   predicted: Receiver function for Gaussian'

        write(LER,*)'SAC header values set'
        write(LER,*)' B     :  delay                          '
        write(LER,*)' USERO :  gwidth        KUSER0:  Rftn    '
        write(LER,*)' USER4 :  rayp (sec/km) '
        write(LER,*)' USER5 :  fit in %      '
        write(LER,*)' KEVNM :  Rftn          KUSER1:  IT_DECON'
        stop
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

        subroutine npow2(npts)
c-----
c      Given npts, determine the N=2**m such that N >= npts
c      return the new ntps
c-----
        integer*4 nsamp, npts
        nsamp = npts
        npts = 1
 1000   continue
        npts = 2*npts
        if(npts.lt.nsamp)goto 1000
        return
        end

        subroutine chtofp(str,fout)
c------
c      routine to convert string to floating point
c      The E format is accepted as well as free form
c      input
c
c      This assumes that the string is 20 characters long
c-----
        character*(*) str
        real*4 fout
        integer*4 lgstr
        logical hase
        l = lgstr(str)
c------
c      If the string str contains an E or e, then
c      we parse the result using an E format.
c------
        hase = .false.
        do 100 i=1,l
            if(str(i:i) .eq.'e' .or. str(i:i).eq.'E')then
                hase = .true.
            endif
  100   continue
c-----
c      read floating point number
c-----
        if(hase)then
            read(str,'(bn,e20.13)')fout
        else
            read(str,'(bn,f20.0)') fout
        endif
        return
        end

        subroutine tdomain(fnum,fden,f,g,p,r,nn,npts,dt,gwidth,theshift,
     1      verbose,maxbumps,lpositive,tol,sacbin,dotwice,rayp)
c-----
c      time domain deconvolution
c-----
c      fnum    Ch* - file name for numerator
c      fden    Ch* - file name for denominator
c      f   R*4 - numerator array
c      g   R*4 - denominator array
c      p   R*4 - array
c      r   R*4 - array
c      nn  I*4 - power of 2 for fft
c      npts    I*4 - number of points in initial time series
c      dt  R*4 - sampling interval
c      gwidth  R*4 - filter width
c      theshift R*4    - time delay for output of result
c      verbose L   - .true. for a lot of intermediate files
c      maxbumps I*4    - number of iterators
c      lpositive L - .true. only ouput positive pulses
c      tol R*4 - convergence tolerance
c      sacbin  L   - .true. sac files are binary else ASCII
c      dotwice L   - .true. use double length FFT to avoid wrap around
c      rayp    R   - ray parameter in sec/km
c-----
        parameter (LER=6, LIN=5, LOT=6)
        parameter (NSAMP=131072)
        character*80 fnum, fden
        real f(NSAMP),g(NSAMP), p(NSAMP), r(NSAMP)
        real gwidth, theshift, tol
        logical lpositive,  verbose, sacbin, dotwice
        integer maxbumps

        parameter (MAXG=1000)
        real amps(MAXG)
        integer shifts(MAXG)

        character*12 resfile, filename
c-----
c      check
c-----
        if(maxbumps .gt. MAXG)then
            maxbumps = MAXG
        endif
        if(maxbumps.lt.0)then
            maxbumps = MAXG
        endif
c-----
c      open the sac files and read in the data
c-----
        if(sacbin)then
            call brsac(1,NSAMP,fnum,f,nerr)
        else
            call arsac(1,NSAMP,fnum,f,nerr)
        endif
        if(nerr.ne.0)then
            write(LER,*)'Problem reading the numerator file'
            call usage()
        endif
c-----
c      get header values
c-----
        call getfhv('DELTA   ',cmpdt,nerr)
        call getnhv('NPTS    ',npts, nerr)
        call getfhv('B       ',btime_num,nerr)
        dt = cmpdt
        nn = npts
c-----
c      find power of 2 that >= number of observed points
c-----
        call npow2(nn)
c-----
c      if use double length FFT , nn -> 2*nn, but comparison is
c      for residuals and output is based on nn before doubling
c-----
        if(dotwice)then
            nwrout = nn
            nn = 2*nn
            nout = nn
        else
            nout = nn
            nwrout = nn
        endif
        call scmxmn(f,npts,depmax,depmin,depmen,indmax,indmin)
        write(LOT,*)'Max Min amplitudes of Numerator  : ',depmin,depmax
c-----
c      filter and save the numerator
c-----
        call gfilter(f,gwidth,npts,nn,dt)
        call scmxmn(f,nout,depmax,depmin,depmen,indmax,indmin)

        call setnhv('NPTS    ',nout, nerr)
        call setfhv('TIMMAX  ',btime_num+indmax*dt,nerr)
        call setfhv('TIMMIN  ',btime_num+indmin*dt,nerr)
        call setfhv('DEPMAX  ',depmax,nerr)
        call setfhv('DEPMIN  ',depmin,nerr)
        call setfhv('DEPMEN  ',depmen,nerr)
        call bwsac(1,NSAMP,'numerator',f)
        call bwsac(1,NSAMP,'observed' ,f)
c-----
c      work on denominator
c-----
        if(sacbin)then
            call brsac(1,NSAMP,fden,g,nerr)
        else
            call arsac(1,NSAMP,fden,g,nerr)
        endif
        if(nerr.ne.0)then
            write(LER,*)'Problem reading the denominator file'
            call usage()
        endif
        call getfhv('B       ',btime_num,nerr)
c-----
c      get header values
c-----
        call scmxmn(g,npts,depmax,depmin,depmen,indmax,indmin)
        write(LOT,*)'Max Min amplitudes of Denominator: ',depmin,depmax
c-----
c      filter and save the denominator
c-----
        call gfilter(g,gwidth,npts,nn,dt)
        call scmxmn(g,nout,depmax,depmin,depmen,indmax,indmin)
        call setnhv('NPTS    ',nout, nerr)
        call setfhv('TIMMAX  ',btime_den+indmax*dt,nerr)
        call setfhv('TIMMIN  ',btime_den+indmin*dt,nerr)
        call setfhv('DEPMAX  ',depmax,nerr)
        call setfhv('DEPMIN  ',depmin,nerr)
        call setfhv('DEPMEN  ',depmen,nerr)
        call bwsac(1,NSAMP,'denominator',g)
c-----
c      compute the power in the numerator for error scaling
c-----
        power = 0.0
        do 100 i=1, nn
            power = power + f(i)*f(i)
 100    continue
c-----
c    correlate the signals
c-----
        call fcorrelate(f,g,nn,nn,dt)
c-----
c    find the peak in the correlation
c-----
        maxlag = nout/2
        write(LOT,'(/,a27,f10.5)') 'The maximum spike delay is ', 
     1   real(maxlag) * dt
c
        if(lpositive) then
            call getmax(g,maxlag,amps(1),shifts(1))
        else
            call getabsmax(g,maxlag,amps(1),shifts(1))
        endif
c-----
c      divide by DELTA t since a unit area impulse consists of
c      a single point with amplitude 1/dt
c-----
        amps(1) = amps(1) / dt
        WRITE(6,*)'maxlag,amps,shifts:',maxlag,amps(1),shifts(1)
c
        nshifts = 1
c-----
c    read in the signals again
c-----
        call zero(f,nn)
        call zero(g,nn)
        call rsac1(fnum,f,ndummy,beg,delta,NSAMP,nerr)
        call rsac1(fden,g,ndummy,b,dt,NSAMP,nerr)
        npts = nn
c-----
c    compute the predicted deconvolution result
c-----
        call zero(p,nn)
        call build_decon(amps,shifts,nshifts,p,npts,nn,gwidth,dt)
        call setnhv('NPTS    ',nn, nerr)
        if(verbose) then
            call phs_shift(p,theshift,nn,nn,dt)      
          call wsac1('d001',p,npts,-theshift,dt,nerr)
            call phs_shift(p,-theshift,nn,nn,dt)      
        endif
c-----
c    convolve the prediction with the denominator signal
c----- 
        call convolve(p,g,npts,nn,dt)

c
        if(verbose) then
          call wsac1('p001',p,npts,beg,dt,nerr)
        endif
c-----
c    filter the signals
c-----     
        call gfilter(f,gwidth,npts,nn,dt)
        call gfilter(g,gwidth,npts,nn,dt)
c
        if(verbose)then
          write(resfile,'(a1,i3.3)') 'r',0
          call wsac1(resfile,f,nout,beg,dt,nerr)
        endif
c-----      
c    compute the residual (initial error is 1.0)
c-----
        call getres(f,p,nout,r,sumsq_ip1)
c
        sumsq_i = 1.0
        sumsq_ip1 = sumsq_ip1 / power
        d_error = 100*(sumsq_i - sumsq_ip1) 
c
          write(resfile,'(a1,i3.3)') 'r',1
        if(verbose)then
          call wsac1(resfile,r,nout,beg,dt,nerr)
        endif
c    
        write(LOT,1000)
        write(LOT,1001)
     &  resfile, dt*amps(1),(shifts(1)-1)*dt,100*sumsq_ip1,
     &  d_error
1000  format(/,1x,'File',9x,
     & 'Spike amplitude   Spike delay   Misfit   Improvement')
1001  format(1x,a10,2x,e16.9,2x,f10.3,3x,f7.2,'%',3x,f9.4,'%')
c
c   
        do while(d_error .gt. tol .and. nshifts .lt. (maxbumps))
c
          nshifts = nshifts + 1
        sumsq_i = sumsq_ip1
c
          call zero(g,nn)
          call rsac1(fden,g,ndummy,b,dt,NSAMP,nerr)
          call gfilter(g,gwidth,npts,nn,dt)
          call fcorrelate(r,g,nn,nn,dt)
        if(lpositive)then
         call getmax(g,maxlag,amps(nshifts),shifts(nshifts))
          else
           call getabsmax(g,maxlag,amps(nshifts),shifts(nshifts))
          endif
          amps(nshifts) = amps(nshifts) / dt
c
          call zero(p,nn)
          call build_decon(amps,shifts,nshifts,p,npts,nn,gwidth,dt)
        if(verbose)then
            write(filename,'(a1,i3.3)') 'd',nshifts
            call phs_shift(p,theshift,nn,nn,dt)      
          call wsac1(filename,p,nout,-theshift,dt,nerr)
            call phs_shift(p,-theshift,nn,nn,dt)      
          endif
c       
        call zero(g,nn)
          call rsac1(fden,g,ndummy,b,dt,NSAMP,nerr)
          call convolve(p,g,npts,nn,dt)
        if(verbose)then
            write(filename,'(a1,i3.3)') 'p',nshifts
          call wsac1(filename,p,nout,beg,dt,nerr)
          endif
c               
          call zero(f,nn)
          call rsac1(fnum,f,ndummy,beg,delta,NSAMP,nerr)
          call gfilter(f,gwidth,npts,nn,dt)
          call getres(f,p,nout,r,sumsq_ip1)
          
          sumsq_ip1 = sumsq_ip1/ power
          write(resfile,'(a1,i3.3)') 'r',nshifts
        if(verbose)then
            call wsac1(resfile,r,nout,beg,dt,nerr)
        endif
        d_error = 100*(sumsq_i - sumsq_ip1)
        
        write(LOT,1001)
     &   resfile,dt*amps(nshifts),(shifts(nshifts)-1)*dt,
     &   100*sumsq_ip1,d_error
c   
        enddo
c
c     
        write(LOT,1010) d_error
1010  format(/,1x,'Last Error Change = ',f9.4,'%',/)
c-----
c    if the last change made no difference, drop it
c-----      
        fit = 100 - 100*sumsq_ip1
c
        if(d_error .le. tol)then
           nshifts = nshifts - 1
           fit = 100 - 100*sumsq_i
           write(LOT,*)'Hit the min improvement tolerance - halting.'
        endif
c
        if(nshifts .ge. maxbumps)then
           write(LOT,*)'Hit the max number of bumps - halting.'
        endif
c
        write(LOT,*)'Number of bumps in final result: ', nshifts
        write(LOT,1011) fit
1011  format(1x,'The final deconvolution reproduces ',
     &    f6.1,'% of the signal.',/)
c-----
c    compute the final prediction
c-----
        call zero(p,nn)
        call build_decon(amps,shifts,nshifts,p,npts,nn,gwidth,dt)
        call zero(g,nn)
        call rsac1(fden,g,ndummy,b,dt,NSAMP,nerr)
        call convolve(p,g,npts,nn,dt)
        call wsac1('predicted',p,nout,beg,dt,nerr)
        call zero(g,nn)
c-----
c    write out the answer
c-----
        call zero(p,nn)
        call build_decon(amps,shifts,nshifts,p,npts,nn,gwidth,dt)
        call phs_shift(p,theshift,nn,nn,dt)      
c-----
c      we use the original numerator file to get header values,
c      related to event location and station, and then modify the fields
c      accordingly
c-----
        call newhdr()
        call rsac1(fnum,g,ndummy,b,dt,NSAMP,nerr)
        call scmxmn(p,nwrout,depmax,depmin,deomen,indmax,indmin)
        call setnhv('NPTS    ',nwrout,nerr)
        call setfhv('B       ',-theshift,nerr)
        call setfhv('TIMMAX  ',-theshift+indmax*dt,nerr)
        call setfhv('TIMMIN  ',-theshift+indmin*dt,nerr)
        call setfhv('DEPMAX  ',depmax,nerr)
        call setfhv('DEPMIN  ',depmin,nerr)
        call setfhv('DEPMEN  ',depmen,nerr)
        theend = -theshift + (nn-1)*dt
        call setfhv('E       ',theend,nerr)
        call setnhv('NZSEC   ',-12345,nerr)
        call setfhv('USER0   ',gwidth,nerr)
        call setkhv('KUSER0  ','Rftn    ',nerr)
        call setkhv('KUSER1  ','IT_DECON',nerr)
        call setkhv('KEVNM   ','Rftn    ',nerr)
        call setkhv('KEVNMC  ','        ',nerr)
        call setfhv('USER5   ',fit,nerr)
C       call setkhv('KUSER5  ','FIT     ',nerr)
        if(rayp.gt.0.0)then
            call setfhv('USER4   ', rayp, nerr)
C           call setkhv('KUSER4  ','p(km/s) ', nerr)
        endif
        call wsac0('decon.out',xdummy,p,nerr)
c-----
c    write out the gaussian filter
c-----
        if(verbose)then
          call newhdr()
          call zero(p,nn)
          p(1) = 1 / dt
          call phs_shift(p,theshift,nn,nn,dt)      
          call gfilter(p,gwidth,npts,nn,dt)
          call wsac1('thefilter',p,nn,beg,dt,nerr)
        endif
        
        return
        end

        subroutine build_decon(amps,shifts,nshifts,p,n,nn,gwidth,dt)
c-----
c      compute the predicted time series from a set of 
c      amplitudes and shifts
c-----
c      amps    R*4 - array of amplitudes
c      shifts  R*4 - array of time shifts
c      p   R*4 - array generated
c      n   I*4 - number of points in array
c      nn  I*4 - number of points rounded up to nearest power of 2
c      gwidth  R*4 - filter parameter
c      dt  R*4 - sampling interval 
c-----
        real p(n), amps(nshifts)
        integer shifts(nshifts)
        integer i, n, nshifts
        
        call zero(p,nn)
        do 1 i = 1, nshifts
            p(shifts(i)) = p(shifts(i)) + amps(i)
    1   continue
        call gfilter(p,gwidth,nn,nn,dt)

        return
        end

        subroutine phs_shift(x,theshift,n,nn,dt)
c-----
c      time shift a signal
c------
c      X   R*4 - signal to be shifted
c      theshift R*4    - time shift in seconds
c      n   I*4 - length of signal
c      nn  I*4 - length of signal to nearest power of 2
c      dt  R*4 - sampling interval
c-----
        real x(n), pi, two_pi, theshift, d_omega
        integer i, n
        integer forward, inverse

        parameter (NSAMP=131072)
        common/zval/z1, z2
        complex z1(NSAMP), z2(NSAMP)
c     
        forward = 1
        inverse = -1
        pi = 3.1415927
        two_pi = 2 * pi

        do 1000 i=1,nn
            if(i.le.n)then
                z1(i) = cmplx(x(i),0.0)
            else
                z1(i) = cmplx(0.0,0.0)
            endif
 1000   continue
c-----
c      get the Fourier transform
c-----
        call zfour(z1,nn,-1,dt,df)
        d_omega = two_pi * df
c-----
c      time shift
c-----
        n21 = nn / 2 + 1
        do 2000 i=1,n21
            ang = (i-1)*d_omega*theshift
            c = cos(ang)
            s = sin(ang)
            z1(i) = z1(i) * cmplx(c,-s)
            if(i.gt.1)then
                z1(nn + 2 - i) = conjg(z1(i))
            endif
 2000   continue
c-----
c      ensure Nyquist frequency element is real
c-----
        z1(n21) = cmplx(real(z1(n21)),0.0)
c-----
c      inverse Fourier transform to the time domain
c-----
        call zfour(z1,nn,+1,dt,df)
c-----
c      redonstitute the time series
c-----
        do 3000 i=1,nn
            x(i) = real(z1(i))
 3000   continue
        return
        end
                
        subroutine fcorrelate(f,g,n,nn,dt)
        real f(n), g(n)
        integer n, nn
        real dt
c-----
c    correlation routine - correlates f and g and replaces the 
c      g with the cross-correlation the value is normalized
c      by the zero-lag autocorrelation of g
c-----
c      f   R*4 - array 1
c      g   R*4 - array 2
c      n   I*4 - number of points in time series
c      nn  I*4 - length of signal to nearest power of 2
c      dt  R*4 - sampling interval
c-----  
        parameter (NSAMP=131072)
        common/zval/z1, z2
        complex z1(NSAMP), z2(NSAMP)

        real sum0
        integer n21, i
c-----
c
c      compute the zero-lag autocorrelation of g
c-----
        sum0 = 0.0
        do 1 i = 1, n
            sum0 = sum0 + g(i)*g(i)
    1   continue
        sum0 = sum0 * dt
c------
c      cross correlate by frequency domain multiplication
c      normalize result by zero lag autocorrelationc
c-----
        do 1000 i=1,nn
            if(i.le.n)then
                z1(i) = cmplx(f(i),0.0)
                z2(i) = cmplx(g(i),0.0)
            else
                z1(i) = cmplx(0.0,0.0)
                z2(i) = cmplx(0.0,0.0)
            endif
 1000   continue
c-----
c      get Fourier transforms
c-----
        call zfour(z1,nn,-1,dt,df)
        call zfour(z2,nn,-1,dt,df)
c-----
c      cross correlation = F G*
c-----
        n21 = nn / 2 + 1
        do 2000 i=1,n21
            z1(i) = z1(i) * conjg(z2(i))
            if(i.gt.1)then
                z1(nn + 2 - i ) = conjg(z1(i))
            endif
 2000   continue
c-----
c      ensure Nyquist frequency element is real
c-----
        z1(n21) = cmplx(real(z1(n21)),0.0)
c------
c      compute inverse Fourier transform
c-----
        call zfour(z1,nn,+1,dt,df)   
c-----
c      update the g array, normalizing in the process
c-----
        do 20 i = 1,nn
            g(i) = real(z1(i))/sum0
   20   continue
        return
        end   

        subroutine zero(x,n)
c-----
c      zero a real array
c-----
c      x   R*4 - array to be zeroed
c      n   I*4 - number of points in array
c-----
        real x(n)
        integer i,n
        
        do 1 i = 1,n
            x(i) = 0
    1   continue
        return
        end

        subroutine getmax(x,n,maxvalue,maxindex)
c-----
c      get maximum value in array and return index to array position
c-----
c      x   R*4 - array
c      n   I*4 - number of points in array
c      maxvalue R*4    - maximum value
c      maxindex I*4    - position of maximum value in array
c-----
        real x(n), maxvalue
        integer i,n,maxindex
        
        maxvalue = x(1)
        maxindex = 1
        do 20 i = 2, n
            if(x(i) .gt. maxvalue) then
                maxvalue = x(i)
                maxindex = i
            endif
   20   continue
        return
        end

        subroutine getabsmax(x,n,thevalue,maxindex)
c-----
c      get maximum value in array and return index to array position
c-----
c      x   R*4 - array
c      n   I*4 - number of points in array
c      thevalue R*4    - maximum value
c      maxindex I*4    - position of maximum value in array
c-----
        real x(n), maxvalue, thevalue
        integer i,n,maxindex
        
        maxvalue = abs(x(1))
        maxindex = 1
        thevalue = x(1)
        do 20 i = 2, n
            if(abs(x(i)) .gt. maxvalue) then
                maxvalue = abs(x(i))
                thevalue = x(i)
                maxindex = i
            endif
   20   continue
        return
        end

        subroutine gfilter(x,gwidth,n,nn,dt)
c-----
c      convolve a function with a unit-area Gaussian filter
c-----
c      x   R*4 - array to be filtered
c      Gwidth
c          R*4 - Filter is exp( - ( pi freq/gwidth)**2)
c      n   I*4 - number of points in time series
c      nn  I*4 - number of points at next largest power of 2
c      dt  R*4 - sampling interval
c-----
        real x(n)
        real gwidth

        parameter (NSAMP=131072)
        common/zval/z1, z2
        complex z1(NSAMP), z2(NSAMP)

        real fac
        integer i, n21

        do 1000 i = 1, nn
            if(i.le.n)then
                z1(i) = cmplx(x(i),0.0)
            else
                z1(i) = cmplx(0.0,0.0)
            endif
 1000   continue
c-----
c      get Fourier transform
c-----
        call zfour(z1,nn,-1,dt,df)
c-----
c      Gaussian filter
c-----
        n21 = nn / 2 + 1
        do 2000 i=1, n21
            freq = ( i - 1) * df
            fac = 3.1415927*freq/gwidth
            fac = fac * fac
            if(fac.gt.50.0)then
                fac = 0.0
            else
                fac = exp( -fac)
            endif
            z1(i) = z1(i) * fac
            if(i.gt.1)then
                z1(nn + 2 - i) = conjg(z1(i))
            endif
 2000   continue
c-----
c      ensure Nyquist frequency element is real
c-----
        z1(n21) = cmplx(real(z1(n21)),0.0)
c-----
c      get inverse Fourier transform
c-----
        call zfour(z1,nn,+1,dt,df)
c-----
c      reconstitute the series
c-----
        do 3000 i=1,nn
            x(i) = real(z1(i))
 3000   continue
        return
        end

        subroutine convolve(x,y,n,nn,dt)
c-----
c      convolve x and y, replacing the x array
c-----
c      x   R*4 - array
c      y   R*4 - array
c      n   I*4 - number of points in time series
c      nn  I*4 - number of points rounded up to next power of 2
c      dt  R*4 - sampling interval
c-----
        real x(n), y(n), dt
        integer n, nn

        parameter (NSAMP=131072)
        common/zval/z1, z2
        complex z1(NSAMP), z2(NSAMP)

c------
c      convolve  by frequency domain multiplication
c-----
        do 1000 i=1,nn
            if(i.le.n)then
                z1(i) = cmplx(x(i),0.0)
                z2(i) = cmplx(y(i),0.0)
            else
                z1(i) = cmplx(0.0,0.0)
                z2(i) = cmplx(0.0,0.0)
            endif
 1000   continue
c-----
c      get Fourier transforms
c-----
        call zfour(z1,nn,-1,dt,df)
        call zfour(z2,nn,-1,dt,df)
c-----
c      convolution = F  G
c-----
        n21 = nn / 2 + 1
        do 2000 i=1,n21
            z1(i) = z2(i) * z1(i)
            if(i.gt.1)then
                z1(nn + 2 - i ) = conjg(z1(i))
            endif
 2000   continue
c-----
c      ensure Nyquist frequency element is real
c-----
        z1(n21) = cmplx(real(z1(n21)),0.0)
c------
c      compute inverse Fourier transform
c-----
        call zfour(z1,nn,+1,dt,df)   
c-----
c      save the correlated function
c-----
        do 20 i = 1,nn
            x(i) = real(z1(i))
   20   continue
        return
        end
 
        subroutine getres(x,y,n,r,sumsq)
c-----
c      get sum square residuals between x and y arrays
c-----
c      x   R*4 - array
c      y   R*4 - array
c      n   I*4 - number of points
c      r   R*4 - array of residuals
c      sumsq   R*4 - sum of residuals squared
c-----
        real x(n), y(n), r(n), sumsq
        integer i,n
        
        sumsq = 0 
        do 20 i = 1, n
            r(i) = x(i) - y(i)
            sumsq = sumsq + r(i)*r(i) 
   20   continue
        return
        end

