c-----
c       23 JUL 2009 - ensured consistent call to FFT routine with
c               complex arguments. We now call zfour
c-----
c       General purpose FFT routines to
c       handle real time functions econominally.
c       These will use one complex N-point FFT to compute the 
c           FT of a 2N real 
c           time series (realft(data,eata,N,isign,dt,df)
c       and compute the FT/IFT of two N point time series 
c           with a single complex
c       N point FFT
c
c       These are similar to routines in the 
c           copyrighted Numerical Recipes Codes
c       but these are derived directly from Brigham.
c-----
        subroutine ftwofft(DATA1, DATA2, FFT1, FFT2, dt, df, N)
c-----
c       E. Oran Brigham
c       The Fast Fourier Transform and Its Applications
c       Prentice Hall, Englewood Cliffs, New Jersey
c       1988
c       ISBN 0-13-307505-2
c-----
c       Find the Discrete Fourier transform of two real time
c       series of length N with just one call to four()
c       Brigham (1988), Sec. 9.3
c-----
c       DATA1   R   - real array 1
c       DATA2   R   - real array 2
c       FFT1    C   - complex transform of DATA1 - output -
c       FFT2    C   - complex transform of DATA2 - output -
c       dt  R   - time sampling interval
c       df  R   - frequency damping interval = 1/(N*dt)
c       N   I   - number of points, power of 2
c-----
        real DATA1(N), DATA2(N)
        complex FFT1(N), FFT2(N)
        real rn, rnmn, in, inmn

        do 1000 i=1,N
            fft1(i) = cmplx( data1(i), data2(i))
 1000   continue
        call zfour(fft1,n,-1,dt,df)
c-----
c       now unwrap
c-----
        n21 = n / 2 + 1
        do 2000 i=1,n21
        if(i.eq.1)then
            rn =  real(fft1(1))
            in = aimag(fft1(1))
            fft1(1) = cmplx(rn, 0.0)
            fft2(1) = cmplx(in, 0.0)
        else
            rn   =  real(fft1(i    ))
            in   = aimag(fft1(i    ))
            rnmn =  real(fft1(n+2-i))
            inmn = aimag(fft1(n+2-i))
            fft1(i) = cmplx(0.5*(rn+rnmn), 0.5*(in-inmn))
            fft2(i) = cmplx(0.5*(in+inmn),-0.5*(rn-rnmn))
            fft1(n+2-i) = conjg(fft1(i))
            fft2(n+2-i) = conjg(fft2(i))
        endif
 2000   continue
        fft1(n21) = cmplx(real(fft1(n21)), 0.0)
        fft2(n21) = cmplx(real(fft2(n21)), 0.0)
        return
        end

        subroutine itwofft(DATA1, DATA2, FFT1, FFT2, dt, df, N)
c-----
c       E. Oran Brigham
c       The Fast Fourier Transform and Its Applications
c       Prentice Hall, Englewood Cliffs, New Jersey
c       1988
c       ISBN 0-13-307505-2
c-----
c       Find the inverse Discrete Fourier transform of two real time
c       series of length N with just one call to four()
c       Brigham (1988), Sec. 9.3
c-----
c       DATA1   R   - real array 1 - output -
c       DATA2   R   - real array 2 - output -
c       FFT1    C   - complex transform of DATA1
c       FFT2    C   - complex transform of DATA2
c       dt  R   - time sampling interval
c       df  R   - frequency damping interval = 1/(N*dt)
c       N   I   - number of points, power of 2
c-----
        real DATA1(N), DATA2(N)
        complex FFT1(N), FFT2(N)
        real rn, rnmn, in, inmn

c-----
c       now wrap
c-----
        n21 = n / 2 + 1
        do 2000 i=1,n
        if(i.eq.1)then
            rn =  real(fft1(1))
            in =  real(fft2(1))
            fft1(1) = cmplx( rn, in)
        else
            in   =  aimag(fft1(i    )) +  real(fft2(i     ))
            inmn = -aimag(fft1(i    )) +  real(fft2(i     ))
            rn   =   real(fft1(i    )) - aimag(fft2(i     ))
            rnmn =   real(fft1(i    )) + aimag(fft2(i     ))
            
            fft1(i)     = cmplx( rn  , in  )
C           fft1(n+2-i) = cmplx( rnmn, inmn)
        endif
 2000   continue

        call zfour(fft1,n,+1,dt,df)
        do 1000 i=1,N
            data1(i) =  real(fft1(i))
            data2(i) = aimag(fft1(i))
 1000   continue
        return
        end

        subroutine realft(data,eata,n,isign,dt)
c-----
c       compute inverse FFT for  2n real FFT with one 
c       N FFT operation isign = +1
c       derived following Brigham (1988), Sec. 9.3
c-----
c       data    R   - array of 2n values
c       n   I   - 2n / 2
c       isign   I   - -1 forward FFT, +1 inverse FFT
c
c       for ISIGN < 0
c       Input :  FORTRAN indexing
c           data(1), data(2), ..., data(2n-1)
c       Output: FORTRAN indexing
c           R(1), R(n21)   where n21 = 2n/2 +1
c           R(2), I(2)
c           ...   ...
c           R(2n),I(2n)
c       we can get others for real time series by Real Even, Imag odd
c
c       for ISIGN > 0
c       Output:  FORTRAN indexing
c           data(1), data(2), ..., data(2n-1)
c       Input : FORTRAN indexing
c           R(1), R(n21)   where n21 = 2n/2 +1
c           R(2), I(2)
c           ...   ...
c           R(2n),I(2n)
c       we can get others for real time series by Real Even, Imag odd
c-----
        real data(n+n)
        real eata(n+n)
        integer n, isign

        real re, ro, ie, io
        real*8 dtheta, s, c, ds, dc

c-----
c       The proper Delta T and Delta F will be imposed later
c       just force the FFT to use ddt = 1 for the FFT, and
c       ddf = 1 for the inverse FFT
c-----
        if(isign.lt.0)then
            ddt = 1.0
            ddf = 0.0
            df = 1.0/(2*n*dt)
            call rfour(data,n,isign,ddt,ddf)
        else
            df = 1.0/(2*n*dt)
            ddt = 0.0
            ddf = 1.0
        endif

c----
c       rearrange
c-----
        dtheta = 3.1415927/real(n)
        dc = cos(dtheta)
        ds = sin(dtheta)
        do 1000 i=1,N

        if(i.eq.1)then
            if(isign.lt.0)then
                eata(1) = (data(1) + data(2))
                eata(2) = (data(1) - data(2))
            else
                eata(1) = (data(1) + data(2)) /2.0
                eata(2) = (data(1) - data(2)) /2.0
            endif
            c = 1.0d+00
            s = 0.0d+00
        else
            jr = 2*i -1
            ji = 2*i
            kr = 2*n + 2*(2 - i) -1
            ki = 2*n + 2*(2 - i) 
            ts = s*dc + c*ds
            tc = c*dc - s*ds
            s = ts
            c = tc

            re =  ( data(jr) +  data(kr) )/2.0
            io =  ( data(ji) -  data(ki) )/2.0
            
            if(isign.lt.0)then
                ro =  ( data(jr) -  data(kr) )/2.0
                ie =  ( data(ji) +  data(ki) )/2.0
                eata(jr) =   re + c*ie - s*ro
                eata(ji) =   io - s*ie - c*ro
                eata(kr) =   re - c*ie + s*ro
                eata(ki) = - io - s*ie - c*ro
            else
                ie = 0.5 * 
     1              ( -s * (data(ji)+data(ki)) + c*(data(jr)-data(kr)))
                ro = 0.5 * 
     1              ( -c * (data(ji)+data(ki)) - s*(data(jr)-data(kr)))
                eata(jr) = re + ro
                eata(ji) = io + ie 
                eata(kr) = re - ro
                eata(ki) = ie - io 
            endif
        endif
 1000   continue
        do 2000 i=1,2*n
        if(isign.lt.0)then
            data(i) = eata(i) * dt
        else if(isign.gt.0)then
            data(i) = 2.0 * eata(i) * df
        endif
 2000   continue
        if(isign.gt.0)then
            call rfour(data,n,isign,ddt,ddf)
        endif
        return
        end

        subroutine rfour(rarr,nn,isign,dt,df) 
c-----
c     THE input is a real array
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
c-----
        real rarr(*) 
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
    1 tempr = rarr(j) 
        tempi = rarr(j+1) 
        rarr(j) = rarr(i) 
        rarr(j+1)=rarr(i+1) 
        rarr(i) = tempr 
        rarr(i+1) = tempi 
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
        tempr=wr*rarr(j)-wi*rarr(j+1) 
        tempi=wr*rarr(j+1)+wi*rarr(j) 
        rarr(j)=rarr(i)-tempr 
        rarr(j+1)=rarr(i+1)-tempi 
        rarr(i)=rarr(i)+tempr 
    8 rarr(i+1) = rarr(i+1)+tempi 
        tempr = wr 
        wr = wr*wstpr-wi*wstpi + wr 
    9 wi = wi*wstpr+tempr*wstpi + wi 
        mmax = istep 
        go to 6 
   10 continue 
        if(isign.lt.0) go to 1002 
c     frequency to time domain 
        do 1001 iiii = 1,n 
 1001 rarr(iiii) = rarr(iiii) * df 
        return 
 1002 continue 
c     time to frequency domain 
        do 1003 iiii = 1,n 
 1003 rarr(iiii) = rarr(iiii) * dt 
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
