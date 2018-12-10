        program sacdecon
c---------------------------------------------------------------------c
c                                                                    c
c     COMPUTER PROGRAMS IN SEISMOLOGY                                c
c     VOLUME V                                                       c
c                                                                    c
c     PROGRAM: SACDECON                                              c
c                                                                    c
c     COPYRIGHT 1998 R. B. Herrmann                                  c
c                                                                    c
c     Department of Earth and Atmospheric Sciences                   c
c     Saint Louis University                                         c
c     221 North Grand Boulevard                                      c
c     St. Louis, Missouri 63103                                      c
c     U. S. A.                                                       c
c                                                                    c
c---------------------------------------------------------------------c
c      CHANGES
c
c      15 JAN 2004 - added the Gaussian filter parameter ALPHA to the
c          code so that can comapre to results of Iterative Decon
c      09 JAN 2005 - thshift undefined in main routine - 
c         change dto the shift  baker@usgs.gov
c       23 JUL 2009 - ensured consistent call to FFT routine with
c               complex arguments. We now call zfour
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
        real alp, gwidth, theshift, wlevel
        logical dowin, sacbin
c-----
c      internal variables
c-----
        parameter (LER=6, LIN=5, LOT=6)
        parameter (NSAMP=131072)
        real y(NSAMP)
        complex zn(NSAMP), zd(NSAMP), zr(NSAMP)
c-----
c      call machine dependent initialization
c-----
        call mchdep()
c-----
c      parse command line arguments
c-----
        call gcmdln(fnum,fden,dowin,alp,sacbin,theshift,wlevel,gwidth)
        WRITE(LER,*)'Numerator            :', fnum
        WRITE(LER,*)'Denominator          :', fden
        WRITE(LER,*)'DO COSINE TAPER      :', dowin
        WRITE(LER,*)'SAC Binary           :',sacbin
        WRITE(LER,*)'alp (time domain)    :',alp
        WRITE(LER,*)'alpha (gauss param)  :',gwidth 
        WRITE(LER,*)'Time shift           :',theshift
        WRITE(LER,*)'Water level          :',wlevel
c-----
c      open the sac files and read in the data
c-----
        if(sacbin)then
            call brsac(1,NSAMP,fnum,y,nerr)
        else
            call arsac(1,NSAMP,fnum,y,nerr)
        endif
c-----
c      get header values
c-----
        call getfhv('DELTA   ',cmpdt,nerr)
        call getnhv('NPTS    ',npts, nerr)
        dt = cmpdt
        nn = npts
        call npow2(nn)
        nn2 = 2*nn
        afac = alp/nn
        fac = 1.0
        dfac = exp(-afac)
        call scmxmn(y,npts,depmax,depmin,depmen,indmax,indmin)
        write(0,*)'Numerator: ',depmin,depmax
c-----
c      zero fill to avoid problems of wrap around
c-----
        do 1101 i=1,nn2
            if(i.le.npts)then
                zn(i) = cmplx(fac*y(i),0.0)
                fac = fac * dfac
            else
                zn(i) = cmplx(0.0,0.0)
            endif
 1101   continue
        call zfour(zn,nn2,-1,dt,df)
c-----
c      work on denominator
c-----

        if(sacbin)then
            call brsac(1,NSAMP,fden,y,nerr)
        else
            call arsac(1,NSAMP,fden,y,nerr)
        endif
c-----
c      get header values
c-----
            afac = alp/nn2
            fac = 1.0
            dfac = exp(-afac)
        call scmxmn(y,npts,depmax,depmin,depmen,indmax,indmin)
        write(0,*)'Denominator: ',depmin,depmax
        do 1102 i=1,nn2
            if(i.le.npts)then
                zd(i) = cmplx(fac*y(i),0.0)
                fac = fac * dfac
            else
                zd(i) = cmplx(0.0,0.0)
            endif
 1102   continue
        call zfour(zd,nn2,-1,dt,df)
c-----
c      do the decon
c-----
        n21 = nn2/2 + 1
        zdmax = 0.0
        do 1103 i=1,n21
            if(cabs(zd(i)).gt.zdmax)then
                zdmax = cabs(zd(i))
            endif
 1103   continue
        if(zdmax.eq.0.0)zdmax = 1.0
        water = wlevel * zdmax
        write(0,*)'zdmax,water:',zdmax,water
c-----
c      decon, but also shift output by about 1/16 number of points
c-----
        write(0,*)'nn2,n21,theshift:',nn2,n21,theshift
        do 1104 i=1,n21
            omega = 6.2831853 * ( i-1)*df
            facd = cabs(zd(i))**2 + water**2
            zr(i) = zn(i)*conjg(zd(i))/facd
c-----
c      time shift
c-----
            omega=6.2831853*(i-1)*df
            theta=omega*theshift
            zr(i) = zr(i) * exp(-afac*theshift)
     1          *cmplx(cos(theta),-sin(theta))
c-----
c      cosine taper
c-----
            if(dowin)then
                tap = 0.5 * (1.0 + cos(3.1415927*(i-1)/nn2))
                zr(i) = zr(i) * tap
            endif
            if(i.gt.1)then
                zr(nn2+2-i)=conjg(zr(i))
            endif
 1104   continue
        zr(n21)=cmplx(real(zr(n21)),0.0)
        call zfour(zr,nn2,+1,dt,df)
        afac = alp/nn2
        fac = 1.0
        dfac = exp(+afac)
        do 1105 i=1,npts
            zr(i) = cmplx(fac*real(zr(i)),0.0)
            fac = fac * dfac
 1105   continue
c-----
c      Gaussian Filter
c-----
        if(gwidth.gt.0.0)then
            call zfour(zr,nn2,-1,dt,df)
            do 1106 i=1,n21
                freq =  ( i-1)*df
                fac = 3.1415927*freq/gwidth
                fac = fac * fac
                    if(fac.gt.50.0)then
                            fac = 0.0
                    else
                            fac = exp( -fac)
                    endif
                    zr(i) = zr(i) * fac
                if(i.gt.1)then
                    zr(nn2+2-i)=conjg(zr(i))
                endif
 1106       continue
            zr(n21)=cmplx(real(zr(n21)),0.0)
            call zfour(zr,nn2,+1,dt,df)
        endif
c-----
c      get the time series
c-----
        do 1107 i=1,npts
            y(i) = real(zr(i))
 1107   continue
        call scmxmn(y,npts,depmax,depmin,depmen,indmax,indmin)
        write(0,*)'Decon: ',depmin,depmax
c-----
c                              output the filtered trace
c-----
        call setfhv('B',-theshift,nerr)
        call setfhv('DEPMAX', depmax, ierr)
        call setfhv('DEPMIN', depmin, ierr)
        call setfhv('DEPMEN', depmen, ierr)
        call setfhv('TIMMAX  ',-theshift + indmax*dt  ,nerr)
        call setfhv('TIMMIN  ',-theshift + indmin*dt  ,nerr)
        theend = -theshift + (npts-1)*dt
        call setfhv('E',theend,nerr)
        call setnhv('NZSEC',-12345,nerr)
        Call setkhv('KUSER0','Decon',nerr)
        call setkhv('KUSER1  ','FD_DECON',nerr)
        if(sacbin)then
            call bwsac (2,LN,'sacdecon.bin',y)
        else
            call awsac (2,LN,'sacdecon.asc',y)
        endif
        end


        subroutine gcmdln(fnum,fden,dowin,alp,sacbin,theshift,
     1     wlevel,gwidth) 
c-----
c      parse command line arguments
c
c      requires subroutine targ() and funtion mnmarg()
c
c-----
c-----
c      fnum    - Ch*   file name for numerator
c      fden    - Ch*   file name for denominator
c      dowin   - L .true. use a cosine taper
c      alp - R time domain damping factor
c      sacbin  L       - .true. if SAC binary output, else ascii
c      theshift - R    time shift
c      wlevel  - R water level for deconvolution 0.001 default
c      gwidth  - R gaussian filter parameter
c-----
        
        character*80 name
        character fnum*(*), fden*(*)
        logical dowin
        logical sacbin
        real alp, theshift, wlevel
        integer mnmarg

        nmarg = mnmarg()
        dowin = .false.
        i = 0
        nfile = 0
        alp = 2.3
        gwidth = 1.0
        sacbin = .true.
        wlevel = 0.01
        theshift = 0.0
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
            else if(name(1:2).eq.'-T')then
                dowin = .true.
            else if(name(1:4).eq.'-ALP' .and. name(1:6).ne.'-ALPHA')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,f10.0)')alp
            else if(name(1:6).eq.'-ALPHA')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,f10.0)')gwidth
            else if(name(1:2).eq.'-W')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,f10.0)')wlevel
            else if(name(1:2).eq.'-D')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,f10.0)')theshift
            else if(name(1:2).eq.'-B')then
                sacbin = .true.
            else if(name(1:2).eq.'-A' .and. name(1:3).ne.'-ALP')then
                sacbin = .false.
            else if(name(1:2).eq.'-?')then
                call usage()
            else if(name(1:2).eq.'-h')then
                call usage()
            endif
        go to 11
   13   continue
        if(nfile.ne.2)call usage()
        return
        end

        subroutine usage()
        integer*4 LER, LIN, LOT
        parameter(LER=6, LIN=5, LOT=6)
        write(LER,*)'USAGE: sacdecon  ',
     1      ' -FN file_num -FD file_den [-W level]',
     2      ' [-T] [-ALP alp] [-D delay] [-B -A] [-?] [-h]',
     3      ' -ALPHA alpha ]'
        write(LER,*)
     1      ' -FN  file_num (default none) numerator'
        write(LER,*)
     1      ' -FD  file_den (default none) denominator'
        write(LER,*)
     1      ' -W   (default 0.001) water level'
        write(LER,*)
     1      ' -T   (default false) cosine taper in freq domain'
        write(LER,*)
     1      ' -ALP alp (default 2.3) complex frequency parameter'
        write(LER,*)
     1      ' -ALPHA alpha (default 1.0) Gaussian filter parameter'
        write(LER,*)
     1      ' -D   delay (0 sec) Begin output delay sec before t=0'
        write(LER,*)
     1      ' -A   (default false) data are SAC ascii'
        write(LER,*)
     1      ' -B   (default true) data are SAC binary'
        write(LER,*)
     1      ' -?   This help message'
        write(LER,*)
     1      ' -h   This help message'
        write(LER,*)'Output files:'
        write(LER,*)'  sacdecon.bin: SAC binary freq domain ',
     1          'deconvolution'
        write(LER,*)'  sacdecon.asc: SAC ASCII  freq domain ',
     1          'deconvolution'

        write(LER,*)'SAC header values'
        write(LER,*)' KUSER0:  Decon ',' KUSER1:  FD_DECON'
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
