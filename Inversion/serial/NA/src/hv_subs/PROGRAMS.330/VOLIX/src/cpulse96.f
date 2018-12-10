        program cpulse96
c-----
c
c     program cpulse96
c
c     cpulse96 -l ntau -t dt -n n -d ndelay -f fname 
c         -F pulsefile -M mult -ALL -EQEX -EQF
c
c
c     in terms of earthquake seismology, the output is ground velocity
c     for a steplike source. To get the same results with the input
c     pulse, we use the convention that the pulse is that of the
c     far-field displacement, e.g., one sided pulse. This is
c     differentiated in the program to make the output
c     effectively ground velocity
c
c     This program takes the output of seis81, and computes the
c     sixteen Green s functions. It also includes the source time
c     function of a parabolic pulse of total length 4 L dt
c-----
c     CHANGES
c     11 SEP 2000 - build in P, SV and SH first arrival times
c         Note this is a noop at present
c     03 OCT 2003 - error in output which would not permit the THF to be
c         written (Mario La Rocca <mlarocca@ov.ingv.it>)
c     05 FEB 2004 - modified to be slightly more tolerant about DT for
c         user supplied pulse - now issues WARNING and not termination
c         mlarocca@ov.ingv.it
c     07 FEB 2005 - add a -Z flag to indicate that the internal 
c         parabolic or triangular pulses are to be zero phase 
c     12 DEC 2007 - set header values evlat, evlon, stlat, stlon to -12345
c       for compatibility of resulting SAC traces
c     23 JUL 2009 - ensured consistent call to FFT routine with
c               complex arguments. We now call zfour
c-----
        parameter (MXDST=100, NSAMP=1024, NSAMP2=4096, NSAMP4=8192)
        parameter (LER=0, LIN=5, LOT=6)
        
        real*4 dst(MXDST), dt(MXDST), t0(MXDST),
     1      vred(MXDST)
        integer*4 npts(MXDST)
        character mname*80
        real*4 tmin(MXDST)
        complex strace(NSAMP2)
        common/grn/green
        real*4 green(NSAMP,16)
        character rfile*80
        real*4 mult
        integer*4 ieqex
        logical doq, dozero
        real duration
        integer*4 ntau

        character ostr*80

        integer jsrc(16)
        data tmin/100*1.0e+30/
        data jsrc/15*1,1*0/
c-----
c     call machine dependent initialization
c-----
        call mchdep()
c-----
        call gcmdln(ntau,ipt,alp,idva,
     1      iodva,ieqex,xmult,rfile,ostr,
     2      doq,dozero)
c-----
c     initialize
c-----
        do i=1,NSAMP2
            strace(i) = cmplx(0.0,0.0)
        enddo
c-----
c     define output Green s functions
c-----
        if(ieqex.eq.0)then
            jsrc(11) = 0
            jsrc(12) = 0
            jsrc(13) = 0
            jsrc(14) = 0
            jsrc(15) = 0
        else if(ieqex.eq.1)then
            jsrc( 1) = 0
            jsrc( 2) = 0
            jsrc( 3) = 0
            jsrc( 4) = 0
            jsrc( 5) = 0
            jsrc( 6) = 0
            jsrc( 7) = 0
            jsrc( 8) = 0
        endif
c-----
c     open waveform output file
c-----
        open(8,file='cseis96.amp',status='old',form='formatted',
     1            access='sequential')
        rewind 8
c-------
c     read through the data file to obtain the earliest
c     arrival time
c-----
        read(8,'(a)')mname
        read(8,100)ndst
        read(8,104)xsour,zsour
        do 1234 i=1,ndst
            read(8,103)dst(i),dt(i),npts(i),t0(i),vred(i)
            if(npts(i) .gt. NSAMP)then
                call usage('No more than -n 1024 points ')
            endif
 1234   continue
  200   continue
            read(8,116)ncode,ii,t,ax,az,phx,phz,pnew,mwave,ros,vel,
     1          sumtq,rsrf,vsrf
c-----
c     ncode  - ray number
c     ii     - ray type
c     t      - arrival time
c     ax     - horizontal amplitude
c     az     - vertical amplitude
c     phx    - horizontal phase
c     phz    - vertical phase
c     pnew   -
c     mwave  - 
c     ros    - density at source
c     vel    - velocity at source (can be P or S depending on ray)
c     sumtq  - t* for the ray
c     rsrf   - density at surface at receiver
c     vsrf   - velocity at surface at receiver
c------
            if(ncode.eq.0)goto 201
            if(t.lt.tmin(ii))tmin(ii)=t
        goto 200
  201   continue
c-----
c-----
        mdst=ndst
        pi = 3.1415927
c-----
c     make the synthetic seismograms as a function of distance
c-----
        do 500 jd=1,mdst
            n = npts(jd)
            do 501 i=1,16
                do 502 j=1,n
                    green(j,i)=0.0
  502           continue
  501       continue
            xr = dst(jd) - xsour
            yr = 0.0

            hs = zsour
            if(vred(jd).gt.0.0)then
                ti = t0(jd) + xr/vred(jd)
            else
                ti = t0(jd)
            endif
c-----
c     if the input is from a given pulse file
c     we read it in. We also fill it up in such a manner to
c     increase the time domain sampling
c-----
            if(rfile.ne.' ')then
                call rdtrc(strace,tau,rfile,
     1              npts(jd),ddt,dtin,mult)
                ndec = 1
                if(dtin.ne.dt(jd))then
                    write(LER,*)'Warning: RFILE dt',ddt,
     1              ' is not same as dt',dt(jd),
     1                  ' specified for distance',dst(jd)
                endif
            else
c-----
c         make analytic pulse
c
c         to save computations, we fill up the source
c         pulse array with the time function and its
c         hilbert transfort
c         but with a fine sampling
c-----
                tau = ntau * dt(jd)
                call maktrc(strace,tau,
     1              npts(jd),ddt,dt(jd),ndec,ipt,
     2              dozero,duration)
            endif
c-----
c     begin loop over distance - depth suites
c-----

            rewind 8
c-----
c     FORTRAN INPUT - NO CARRIAGE CONTROL
c-----
  100   format(26i3)
  103   format(2e12.4, i10, 2e12.4)
  104   format(8f10.5)
  116   format(2i5,6e11.4,i5,5e11.4)
c-----
            read(8,'(a)')mname
            read(8,100)ndst
            read(8,104)xsour,zsour
        do 1235 i=1,ndst
            read(8,103)dst(i),dt(i),npts(i),t0(i),vred(i)
 1235   continue
  600       continue
            tau = ntau*dt(i)
            n = npts(jd)
            read(8,116)ncode,ii,t,ax,az,phx,phz,pnew,mwave,ros,vel,
     1          sumtq,rsrf,vsrf
c-----
c     ncode   I*4 - sequential number of ray
c     ii  I*4 - distance index
c     ax  R*4 - amplitude of horizontal motion
c     az  R*4 - amplitude of vertical motion (P-SV only)
c     phx R*4 - phase of horizontal motion
c     phz R*4 - phase of vertical motion (P-SV only)
c     pnew    R*4 - source ray angle of incidence in radians 
c                 measured from horizontal 
c                 positive downward
c                 negative upward
c     mwave   I*4 - 1 ray leaves focal sphere as P
c               2 ray leaves focal sphere as SV
c               3 ray leaves focal sphere as SH
c     ros R*4 - density at source point
c     vel R*4 - wave velocity at source point
c               (P-vel for mwave =1, S-vel for mwave = 2 or 3)
c     sumtq   R*4 - t* for the ray at 1.0 Hz. 
c                 Frequency independent Q assumed
c     rsrf    R*4 - surface density
c     vsrf    R*4 - surface velocity of ray - 
c                 depends on the incident ray, e.g.,
c                 P or S on last leg
c     
c                 
c-----
            if(ncode.eq.0)goto 601
            if(ii.ne.jd)goto 600
c-----
c     to fix reversed radial direction to make radial positive away 
c     from source
c-----
            if(mwave.ne.3)phx = phx + pi
            tshift = t - ti
            angle = 3.1415927/2.0  + pnew
            if(mwave.eq.1 .or. mwave.eq.2)then
                call source(1,mwave,angle,amdd,phdd,ros,vel)
                call srccon(1,strace,az,phz,amdd,phdd,
     1              n,tshift,ddt,dt,ndec,jsrc(1),doq,sumtq)
                call source(2,mwave,angle,amdd,phdd,ros,vel)
                call srccon(2,strace,ax,phx,amdd,phdd,
     1              n,tshift,ddt,dt,ndec,jsrc(2),doq,sumtq)
                call source(3,mwave,angle,amds,phds,ros,vel)
                call srccon(3,strace,az,phz,amds,phds,
     1              n,tshift,ddt,dt,ndec,jsrc(3),doq,sumtq)
                call source(4,mwave,angle,amds,phds,ros,vel)
                call srccon(4,strace,ax,phx,amds,phds,
     1              n,tshift,ddt,dt,ndec,jsrc(4),doq,sumtq)
                call source(6,mwave,angle,amss,phss,ros,vel)
                call srccon(6,strace,az,phz,amss,phss,
     1              n,tshift,ddt,dt,ndec,jsrc(6),doq,sumtq)
                call source(7,mwave,angle,amss,phss,ros,vel)
                call srccon(7,strace,ax,phx,amss,phss,
     1              n,tshift,ddt,dt,ndec,jsrc(7),doq,sumtq)
                call source(9,mwave,angle,amex,phex,ros,vel)
                call srccon(9,strace,az,phz,amex,phex,
     1              n,tshift,ddt,dt,ndec,jsrc(9),doq,sumtq)
                call source(10,mwave,angle,amex,phex,ros,vel)
                call srccon(10,strace,ax,phx,amex,phex,
     1              n,tshift,ddt,dt,ndec,jsrc(10),doq,sumtq)
                call source(11,mwave,angle,amvf,phvf,ros,vel)
                call srccon(11,strace,az,phz,amvf,phvf,
     1              n,tshift,ddt,dt,ndec,jsrc(11),doq,sumtq)
                call source(12,mwave,angle,amvf,phvf,ros,vel)
                call srccon(12,strace,ax,phx,amvf,phvf,
     1              n,tshift,ddt,dt,ndec,jsrc(12),doq,sumtq)
                call source(13,mwave,angle,amhf,phhf,ros,vel)
                call srccon(13,strace,az,phz,amhf,phhf,
     1              n,tshift,ddt,dt,ndec,jsrc(13),doq,sumtq)
                call source(14,mwave,angle,amhf,phhf,ros,vel)
                call srccon(14,strace,ax,phx,amhf,phhf,
     1              n,tshift,ddt,dt,ndec,jsrc(14),doq,sumtq)
            else if(mwave.eq.3)then
                call source(5,mwave,angle,amdt,phdt,ros,vel)
                call srccon(5,strace,ax,phx,amdt,phdt,
     1              n,tshift,ddt,dt,ndec,jsrc(5),doq,sumtq)
                call source(8,mwave,angle,amst,phst,ros,vel)
                call srccon(8,strace,ax,phx,amst,phst,
     1              n,tshift,ddt,dt,ndec,jsrc(8),doq,sumtq)
                call source(15,mwave,angle,amht,phht,ros,vel)
                call srccon(15,strace,ax,phx,amht,phht,
     1              n,tshift,ddt,dt,ndec,jsrc(15),doq,sumtq)
            endif



            goto 600
  601       continue
c-----
c     now that computations are done
c     integrate point force solutions
c-----
                call ginteg(11,dt,n)
                call ginteg(12,dt,n)
                call ginteg(13,dt,n)
                call ginteg(14,dt,n)
                call ginteg(15,dt,n)
c-----
c     output
c-----

            xr = dst(jd) - xsour
            yr = 0.0

            hs = zsour
            nt = n

        call output(tau,idva,iodva,ieqex,
     1      ostr,
     3      xr, yr, hs, n, dt, ti -0.5*duration, mname)
  500   continue
        close(8)
        end

        subroutine rdtrc(strace,tau,rfile,n,ddt,dt,mult)
c-----
c     make the pulse and its Hilbert Transform
c     using a wavelet read in
c
c     green(NSAMP,10) - R*4   temporary storage
c     strace          - Complex
c      real strace(i) -    pulse
c      imag strace(i) -    Hilbert transform of pulse
c     tau     - R*4   (pulse has duration 4*tau)
c     n       - I*4   points in final time series
c     ddt     - I*4   sampling interval of pulse in seconds
c     dt      - R*4   sampling interval in seconds
c     mult        - R*4   scaling factor
c-----

        integer NSAMP, NSAMP2
        parameter(NSAMP = 1024,NSAMP2=4096)
        complex strace(NSAMP2)
        common/grn/green
        real*4 green(NSAMP,16)
        character rfile*80
        real*4 mult, tau
        complex z1, z2

c-----
c     read in the trace information
c-----
        open(1,file=rfile,status='unknown',form='formatted',
     1      access='sequential')
        rewind 1
        read(1,'(i5,f10.3)')npls,ddt
        read(1,10)(green(i,1),i=1,npls)
   10   format(4e15.7)
        close(1)
        tau = npls*ddt
c-----
c     now we play games with the difference between dt and ddt
c-----
        do 100 i=1,NSAMP2
            if(i.le.npls)then
                strace(i) = cmplx(mult*green(i,1),0.0)
            else
                strace(i) = cmplx(0.0,0.0)
            endif
  100   continue
        call zfour(strace,NSAMP2,-1,ddt,df)
        n21 = NSAMP2 / 2 + 1
c-----
c     form spectrum for the Hilbert and also shift it to the
c     center of the resultant time history. The flip causes
c     the time series to be shifted to the center.
c-----
        flip = 1.0
        do 200 i=1,n21
            om = 6.2831853*(i-1)*df
            z1 = strace(i)
            z2 = z1 * cmplx(0.0,om)
            strace(i) = 2.0*flip*z2
            if(i.gt.1)then
                j = NSAMP2 + 2 -i
                strace(j) = cmplx(0.0,0.0)
            endif
            flip = - flip
  200   continue
        call zfour(strace,NSAMP2,+1,ddt,df)
c-----
c     change the polarity of the quadrature component
c-----
        do 300 i=1,NSAMP2
            strace(i) = conjg(strace(i))
  300   continue
        return
        end

        subroutine maktrc(strace,tau,n,ddt,dt,ndec,ipt,
     1      dozero,duration)
c-----
c     make the pulse and its Hilbert Transform
c
c     strace          - Complex
c      real strace(i) -    pulse
c      imag strace(i) -    Hilbert transform of pulse
c     tau     - R*4   (pulse has duration 4*tau)
c     n       - I*4   points in final time series
c     ddt     - I*4   sampling interval of pulse in seconds
c     dt      - R*4   sampling interval in seconds
c     ndec        - I*4   for pulse resolution finer sampling
c     dozero      - L zero phase for tringular/parabolic pulse
c     duration    - R*4   duration of pulse for triangular/parabolic
c-----

        integer NSAMP, NSAMP2
        parameter(NSAMP = 1024,NSAMP2=4096)
        complex strace(NSAMP2)
        real*4 tau
        integer n, ndec, ipt
        real ddt, dt, duration
        logical dozero
        real tr, ti
c-----
c     first define ndec, so that we do not exceed array sizes
c     we trade off ultimate time resolution with fitting into
c     array size
c-----
        nmax = 8.0*tau /dt
        ndec = n/nmax
        if(ndec.gt.10)then
            ndec = 10
        else if(ndec.lt.1)then
            ndec = 1
        endif
        ddt = dt/real(ndec)
        
        do 150 i=1,NSAMP2
            td = (i-(NSAMP2/2))*ddt - 1.0*tau
            if(ipt.eq.0)then
                strace(i) = cmplx(tri(td,tau),trih(td,tau))/(2.*tau)
            else
c-----
c     the division by 2*tau is to create the proper
c     normalization for a parabolic pulse.
c-----
                tr=tri(td,tau) - tri(td-2.*tau,tau)
                ti=trih(td,tau)- trih(td-2.*tau,tau)
                strace(i) = cmplx(tr,ti)/(2.*tau)
            endif
  150   continue
        if(ipt.eq.0)then
            duration = 2.0*tau
        else if(ipt.eq.1)then
            duration = 4.0*tau
        endif
        if(.not. dozero)then
            duration = 0.0
        endif
        return
        end

c-----
c     build up seismograms
c
c     *****************************************************
c
        subroutine source(num,mwave,ang,am,ph,ros,vel)
c-----
c     num R*4 - Green's function desired
c     mwave   R*4 - 1 P-wave
c               2 SV-wave
c               3 SH-wave
c     angle   R*4 - ray angle measured from downward vertical in radians
c     am  R*4 - excitation amplitude
c     ph  R*4 - excitation phase
c     ros R*4 - density at source
c     vel R*4 - ray velocity at source
c-----
        r1 = sin(ang)
        h1 = cos(ang)
        rr = r1*r1
        hh = h1*h1
        hr = r1*h1
        den2 = 1.0/(12.5667*ros*vel**2)
        den3 = 1.0/(12.5667*ros*vel**3)
        if(num.eq. 1)then
            pamp  = (2.0*hh - rr)*den3
            svamp = 3.0*hr *den3
            shamp = 0.0
        else if(num.eq. 2)then
            pamp  = (2.0*hh - rr)*den3
            svamp = 3.0*hr *den3
            shamp = 0.0
        else if(num.eq. 3)then
            pamp  = -2.0*hr * den3
            svamp = (hh - rr) *den3
            shamp = 0.0
        else if(num.eq. 4)then
            pamp  = -2.0*hr * den3
            svamp = (hh - rr) *den3
            shamp = 0.0
        else if(num.eq. 5)then
            pamp  = 0.0
            svamp = 0.0
            shamp =  h1 * den3
        else if(num.eq. 6)then
            pamp  = rr * den3
            svamp = -hr * den3
            shamp = 0.0
        else if(num.eq. 7)then
            pamp  = rr * den3
            svamp = -hr * den3
            shamp = 0.0
        else if(num.eq. 8)then
            pamp = 0.0
            svamp = 0.0
            shamp = -r1 * den3
        else if(num.eq. 9)then
            pamp = den3
            svamp = 0.0
            shamp = 0.0
        else if(num.eq.10)then
            pamp  = den3
            svamp = 0.0
            shamp = 0.0
        else if(num.eq.11)then
            pamp  = - h1 * den2
            svamp = - r1 * den2
            shamp = 0.0
        else if(num.eq.12)then
            pamp  = - h1 * den2
            svamp = - r1 * den2
            shamp = 0.0
        else if(num.eq.13)then
            pamp  =   r1 * den2
            svamp = - h1 * den2
            shamp = 0.0
        else if(num.eq.14)then
            pamp  =   r1 * den2
            svamp = - h1 * den2
            shamp = 0.0
        else if(num.eq.15)then
            pamp  = 0.0
            svamp = 0.0
            shamp = - den2
        else if(num.eq.16)then
            pamp = 0.0
            svamp = 0.0
            shamp = 0.0
        endif
        am = 0.0
        ph = 0.0
        if(mwave.eq.1)then
            am = abs(pamp)
            if(pamp .lt. 0.0)ph = 3.1415927
        else if(mwave.eq.2)then
            am = abs(svamp)
            if(svamp .lt. 0.0)ph = 3.1415927
        else if(mwave.eq.3)then
            am = abs(shamp)
            if(shamp .lt. 0.0)ph = 3.1415927
        endif
        return
        end

        subroutine gcmdln(ntau,ipt,alp,idva,
     1      iodva,ieqex,xmult,rfile,ostr,
     2      doq,dozero)
c-----
c     parse command line arguments
c     requires subroutine mgtarg() and function mnmarg()
c-----
c     vred    R*4 - reduction velocity
c     t0  R*4 - reduced travel time tfirst = t - x/vred - t0
c     dt  R*4 - desired sample interval
c     n   I*4 - number of data points
c
c     ntau    I*4 - pulse duration factor for parabolic and triangular
c             pulses
c     ipt I*4 - pulse type
c             0 - triangular
c             1 - parabolic
c             2 - Ohnaka
c             3 - Dirac delta function
c             4 = user pulse in file rfile
c     alp R*4 - Ohnaka pulse shape factor, fc= alp / 2 pi
c     idva    I*4 - time history type
c             0 - displacement
c             1 - velocity
c             2 - acceleration
c     iodva   I*4 - force output time history type
c             Useful if internal pulse is not to reprsent a
c             step in fault displacement, explosion pressure,
c             or force
c             0 - displacement
c             1 - velocity
c             2 - acceleration
c     ieqex   I*4 -   Green's function choice
c             0 = earthquake and explosion Green's functions
c             1 = explosion and point force Green's functions
c             2 = earthquake, explosion and point force Green's func
c     xmult   R*4 - moment scaling factor
c     rfile   C*80- name of user provided pulse 
c             note dt must be the same as specified in dfile
c     doq L   - .true. do causal Q
c     dozero  L   .true. for triangular and parabolic pulses, center
c             at zero lag
c-----
        character*(*)  rfile
        integer*4 ntau, ipt, idva, iodva, ieqex 
        real*4 xmult, alp
        logical  doq, dozero
        character ostr*(*)

        integer*4 mnmarg
        character*80 name
c-----
c     initialization
c-----
        ntau = -1
        ipt = -1
        alp = -1.0
        idva = 1
        iodva = -1
        ieqex = 2
        rfile = ' '
        xmult = 1.0
        ostr = ' '
        doq = .false.
        dozero = .false.
        nmarg=mnmarg()
        i = 0
   11   i = i + 1
            if(i.gt.nmarg)go to 13
            call mgtarg(i,name)
            if(name(1:2).eq.'-l')then
                i=i+1
                call mgtarg(i,name)
                read(name,'(bn,i10)')ntau
            else if(name(1:2).eq.'-a')then
                i=i+1
                call mgtarg(i,name)
                read(name,'(bn,f10.0)')alp
            else if(name(1:2).eq.'-t')then
                ipt = 0
            else if(name(1:2).eq.'-p')then
                ipt = 1
            else if(name(1:2).eq.'-o')then
                ipt = 2
            else if(name(1:2).eq.'-i')then
                ipt = 3
            else if(name(1:2).eq.'-F')then
                ipt = 4
                        i=i+1
                        call mgtarg(i,rfile)
            else if(name(1:2).eq.'-D' )then
                idva = 0
            else if(name(1:2).eq.'-V' )then
                idva = 1
            else if(name(1:2).eq.'-A' .and. name(1:4).ne.'-ALL')then
                idva = 2
            else if(name(1:3).eq.'-OD')then
                iodva = 0
            else if(name(1:3).eq.'-OV')then
                iodva = 1
            else if(name(1:3).eq.'-OA')then
                iodva = 2
            else if(name(1:3).eq.'-EQ')then
                ieqex = 0
            else if(name(1:3).eq.'-EX')then
                ieqex = 1
            else if(name(1:4).eq.'-ALL')then
                ieqex = 2
            else if(name(1:2).eq.'-m')then
                i=i+1
                call mgtarg(i,name)
                read(name,'(bn,f20.0)')xmult
            else if(name(1:2).eq.'-Q')then
                doq = .true.
            else if(name(1:2).eq.'-Z')then
                dozero = .true.
            else if(name(1:2).eq.'-?')then
                call usage('Help')
            else if(name(1:2).eq.'-h')then
                call usage('Help')
            endif
        go to 11
   13   continue
c-----
c   do some error checking, e.g., we cannot generate velocity
c   for triangular pulse
c-----
        if(ipt.ge.0 .and. ipt .le.1 .and . ntau.le.0)ntau = 1
        if(ipt.eq.2 .and. alp .le.0.0)alp = 1.0
        if(ipt.eq.2 .and. alp.lt.0.0)
     1      call usage('No alpha for Ohnaka pulse')
        if(ipt.lt.0)
     1      call usage('No pulse shape defined')
        if(iodva.lt.0)iodva = idva
        lostr = len(ostr)
        ostr = 'cpulse96'
        do 14 i=1,nmarg
            call getarg(i,name)
            L = lgstr(name)
            lo = lgstr(ostr)
            lomax = lo + 1 + L
            if(lomax .lt. lostr)then
                ostr = ostr(1:lo)//' '//name(1:L)
            endif
   14   continue
        return
        end

        subroutine usage(str)
        parameter (LER=0, LIN=5, LOT=6)
        character str*(*)
        write(LER,*)'cpulse96:',str
        write(LER,*)'USAGE: ',
     1  'cpulse96 [ -v  ] [ -t -o -p -i ] -a alpha',
     2  ' -l L [ -D -V -A]  [-F rfile ] [ -m mult] ',
     3  ' [ -OD -OV -OA ] [-Z]  ',
     4  ' [-Q]',
     4  ' [-EQEX -EXF -ALL]  [-?] [-h]'

        write(LER,*)
     1  ' -v          Verbose output'
        write(LER,*)
     1  ' -t          Triangular pulse of base 2 L dt'
        write(LER,*)
     1  ' -p          Parabolic Pulse of base  4 L dt'
        write(LER,*)
     1  ' -l L        (default 1 )duration control parameter'
        write(LER,*)
     1  ' -D          Output is ground displacment'
        write(LER,*)
     1  ' -V          Output is ground velocity (default)'
        write(LER,*)
     1  ' -A          Output is ground acceleration'
        write(LER,*)
     1  ' -F rfile    User supplied pulse'
        write(LER,*)
     1  ' -m mult     Multiplier (default 1.0)'
        write(LER,*)
     1  ' -OD         Output is ground displacement'
        write(LER,*)
     1  ' -OV         Output is ground velocity'
        write(LER,*)
     1  ' -OA         Output is ground acceleration'
        write(LER,*)
     1  ' -Q          (default false) do causal Q'
        write(LER,*)
     1  ' -EXF        (default) Explosion/point force green functions'
        write(LER,*)
     1  ' -EQEX       Earthquake and double couple green functions'
        write(LER,*)
     1  ' -ALL        Earthquake, Explosion and Point Force '
        write(LER,*)
     1  ' -Z           (default false) zero phase ',
     2  'triangular/parabolic pulse'
        write(LER,*)
     1  ' -?          Write this help message'
        write(LER,*)
     1  ' -h          Write this help message'
        stop
        end

        subroutine ginteg(ngrn,dt,n)
        parameter (MXDST=100, NSAMP=1024, NSAMP2=4096, NSAMP4=8192)
        common/grn/green
        real*4 green(NSAMP,16)
        vfirst = green(1,ngrn)
        sum = 0.0
        do 1000 i=1,n
            sum = sum + (green(i,ngrn) - vfirst) * dt
            green(i,ngrn) = sum
 1000   continue
        return
        end

        subroutine srccon(k,strace,a,p,am,ph,n,
     1            tshift,ddt,dt,ndec,jsrc,doq,sumtq)
        parameter (MXDST=100, NSAMP=1024, NSAMP2=4096, NSAMP4=8192)
        common/grn/green
        real*4 green(NSAMP,16)
        complex strace(NSAMP2)
        logical doq
        real*4 sumtq
c-----
c     sumtq   R*4 - SUM T/Q for this ray, frequency independent Q
c-----

        complex zdata(NSAMP)
        complex ztmp
c-----
c     The pulse is centered on strace but has a sample interval
c     of ddt seconds, while the final trace has a sampling interval
c     of dt seconds. No interpolation is done. However, by making ddt
c     smaller than dt, we effectively have some interpolation
c-----
        if(jsrc.eq.0)return
        amp = a*am
        phz = p+ph
        c = amp*cos(phz)
        s = amp*sin(phz)
        if(ndec.lt.1)ndec = 1
        nshift= ((tshift/ddt ) + 0.5)
        if(nshift.lt.1)nshift = 1
        do 100 i=1,n
            j = ndec*i - nshift + (NSAMP2)/2
            if(j.ge.1.and.j.le.NSAMP2)then
                strac1=real(strace(j))
                strac2=aimag(strace(j))
                zdata(i) = cmplx( c*strac1 + s*strac2 , 0.0)
            else
                zdata(i) = cmplx(0.0, 0.0)
            endif
  100   continue
c-----
c     if we require causal Q, convolve with the Futterman operator here
c-----
        if(doq)then
            call zfour(zdata,n,-1,dt,df)
c-----
c         apply t*filter
c-----
            n21 = n/2 + 1
            do 300 i=1,n21
                freq = (i-1)*df
c-----
c             prevent bad things with logs of zero
c-----
                if(freq.lt.df)freq = df
                call causlq(freq,sumtq,0.0,ztmp)
                zdata(i) = zdata(i)*ztmp
                if(i.gt.1)then
                    zdata(n+2-i) = conjg(zdata(i))
                endif
  300       continue
c-----
c         Inverse FFT
c-----
            call zfour(zdata,n,+1,dt,df)
        endif
c-----
c     update the output Green s functions for this ray
c-----
        do 500 i = 1, n
            green(i,k) = green(i,k) + real(zdata(i))
  500   continue
        return
        end

        function tri(t,tau)
c-----
c     triangular source time function
c     a causal
c-----
        if(t.lt.-tau)then
                tri=0.0
        else if(t.ge. -tau .and. t .le. 0.0)then
                tri = (1. + (t/tau))/tau
        else if(t.ge.0.0 .and. t .le. tau)then
                tri = (1. - (t/tau))/tau
        else
                tri = 0.0
        endif
        return
        end

        function trih(t,tau)
c-----
c     hilbert transform of triangular pulse
c-----
        if(t.lt.0.0)then
                sign= +1.0
        else
                sign= -1.0
        endif
        tt = abs(t)
        if(tt.eq.0.0)then
                trih = 0.0
        else if(tt.gt.0.0.and.tt.lt.tau)then
                trih = (tt/tau)*alog( (tau*tau-tt*tt)/(tt*tt) ) -
     1                  alog ( (tau-tt)/(tau+tt) )
        else if(tt.eq.tau)then
                trih = 2.* alog(2.)
        else if(tt.gt.tau)then
                trih = (tt/tau)*alog( (tt*tt-tau*tau)/(tt*tt) ) -
     1                   alog( (tt-tau)/(tt+tau) )
        endif
        trih = sign * trih/(3.1415927*tau)
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
        subroutine npow2(nsamp,npts,npts21)
c-----
c     Given nsamp, find npts >= nsamp such that npts is a power of 2
c-----  
        integer*4 nsamp, npts, npts21
        npts = 1
 1000   continue
        npts = 2*npts
        if(npts.lt.nsamp)goto 1000
        npts21 = npts/2 + 1
        return
        end


        
        subroutine causlq(freq,ts,nu,temp)
            real freq, ts, nu
            complex temp
            fref = 1.0
            pi = 3.1415927
            facr = -pi*ts*freq**(1.0-nu)
            if(nu .eq. 0.0)then
                faci = 2.*pi*freq*ts*alog(freq/fref)/pi
            else
                cot = cos(0.5*pi*nu)/sin(0.5*pi*nu)
                faci = 2.*pi*freq*ts*cot*(1. - (fref/freq)**nu)/2.
            endif
            temp = cexp(cmplx(facr,faci))
        return
        end


        subroutine minteg(x,nt,dt,dodc)
        parameter (LER=0, LIN=5, LOT=6)
c-----
c     Form the integral x(t) dt numerically
c     x   R*4 - array to be integrated
c     nt  I*4 - number of data points in time series
c     dt  R*4 - sampling interval
c     dodc    L*4 - .true. remve x(1) before summing
c-----
        
        real x(nt)
        integer*4 nt
        real*4 dt
        logical dodc
        sum=0.0
        if(dodc)then
            do 100 k = 1,nt
                sum = sum + dt*(x(k) - x(1))
                x(k)=sum
  100       continue
        else
            do 200 k = 1,nt
                sum = sum + dt*x(k)
                x(k)=sum
  200       continue
        endif
        return
        end

        subroutine mderiv(y,n,dt)
c-----
c     Use centered difference to take derivative of time series
c-----
c     y   R*4 - time series to be differentiated and returned as y
c     n   I*4 - length of time series
c     dt  R*4 - sample interval
c-----
        real*4 y(1)
        b0=1.0/dt
        xm1=y(1)
        xm2=0.
        y1=0.
        do 11 i=1,n
            y2=b0*(y(i)-xm1) 
            xm2=xm1
            xm1=y(i)
            y(i)=y2
   11   continue
        return
        end
        

        subroutine output(tau,idva,iodva,ieqex,
     1      ostr,
     3      dist, yr, depths, npts, dt, tfirst, mname)
c-----
c     write the Green's functions for this distance
c-----
c     tau R*4 - pulse duration parameter
c     idva    I*4 - displacmeent, velocity or acceleration
c     dist    R*4 - distance
c     dt  R*4 - sample rate
c     npts    I*4 - number of points in time series
c     tfirst  R*4 - time of first sample
c     depthr  R*4 - receiver depth
c     depths  R*4 - source depth
c     iopen9  I*4 - indicate if temp file is used
c     kkl I*4 - = 1 indicates Love wave eigenfunctions
c     kkr I*4 - = 1 indicates Rayleigh wave eigenfunctions
c     ostr    C*(*)   - command line string
c     mname   C*80    - name of model file
c-----
c-----
c     variables in common blocks
c
c     iftype  I*4 File type
c             1 - single trace
c             3 - three component
c             16 - Green's function
c
c     iobsyn  I*4 1 - observed
c             2 - synthetic
c
c     itmfrq  I*4 1 - time series
c             2 - Fourier spectra (not implemented)
c
c     iunit   I*4 1   - counts 
c             2   - cm
c             3   - cm/sec
c             4   - cm/sec/sec
c             5   - m
c             6   - m/sec
c             7   - m/sec/sec
c             8   - microns
c             9   - microns/sec
c             10  - microns/sec/sec
c     junit   I*4
c             11  - Pa  (nt/m^2)
c             12  - MPa  (mega Pascal)
c
c     cfilt   C*80    comment on filtering operations
c
c     keyear  I*4 event year
c     kemon   I*4 event mon
c     keday   I*4 event day
c     kehour  I*4 event hour
c     kemin   I*4 event minute
c     esec    R*4 event second
c     evlat   R*4 event latitude
c     evlon   R*4 event longitude
c     evdep   R*4 event depth
c     
c     stname  C*8 station name
c             NOTE: AH(6), SEED(5), CSS(6), SAC(8)
c     
c     stlat   R*4 station latitude
c     stlon   R*4 station longitude
c     stelev  R*4 station elevation
c
c
c     distkm  R*4 epicentral distance in km
c     distdg  R*4 epicentral distance in degrees (111.195 km/degree)
c     evstaz  R*4 event -> station azimuth, degrees east of north
c     stevaz  R*4 station -> event azimuth, degrees east of north
c     
c     cpulse  C*80    pulse description
c
c     ccomnt  C*80    comment
c
c     jsrc    I*4 Array of indices to indicate if traces are 
c                 present
c             iftype =  1  jsrc(1) indicates if trace is 
c                 present
c             iftype =  3  jsrc(i),i=1,5 indicates if trace 
c                 is present
c                 i = Z, 2=N, 3=E, 4=R, 5=T
c                 (we could carry all five, but the real 
c                 purpose  is to permit no more than 
c                 3 traces, but 2 may all that are 
c                 present in a real data set)
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
     2      TP, TSV, TSH,
     3      sa, sc, sf, sl, sn, sr
        real*4 esec, distkm, distdg, evstaz, stevaz
        real*4 evlat, evlon, evdep 
        real*4 stlat, stlon, stelev
        REAL*4 TP, TSV, TSH
        REAL*4 SA, SC, SF, SL, SN, SR

        common/chdr96/ stname, cfilt, cpulse, ccomnt
        character*8 stname*8, cfilt*80, cpulse*80, ccomnt*80    

        character stcomp*8
        real*4 cmpaz, cmpinc, cmpdt
        real*4 ssec
c-----
c     internal variables
c-----
        parameter (LER=0, LIN=5, LOT=6)
        parameter (NGRN=16)
        integer*4  idva, iodva, ieqex, npts
        real*4 dt, tfirst, depths, depthr, tau
        character ostr*(*)
        character mname*80

        integer NSAMP
        parameter(NSAMP=1024)
        common/grn/green
        real*4 green(NSAMP,16)
        real*4 x(NSAMP)
        complex z
        integer iszrt(16)
        character*8 ost(16)
        data iszrt/1,4, 1,4,5, 1,4,5, 1,4, 1,4, 1,4,5,1/
        data ost/'ZDD     ', 'RDD     ', 'ZDS     ', 'RDS     ',
     1       'TDS     ', 'ZSS     ', 'RSS     ', 'TSS     ',
     2       'ZEX     ', 'REX     ', 'ZVF     ', 'RVF     ',
     3       'ZHF     ', 'RHF     ', 'THF     ', 'PEX     '/

c-----
c     establish the jsrc array for the file16(V) header
c-----
        do 100 i=1,21
            if(i.le.15)then
                jsrc(i) = iszrt(i)
                if(ieqex.eq.0 .and. i.ge.11)then
                    jsrc(i) = 0
                else if(ieqex.eq.1 .and. i.le.8)then
                    jsrc(i) = 0
                endif
            else
                jsrc(i) = 0
            endif
  100   continue
        los = lgstr(ostr)
        cpulse = ostr(1:los)
        if(idva.eq.0)then
            iunit = 2
        else if(idva.eq.1)then
            iunit = 3
        else if(idva.eq.2)then
            iunit = 4
        endif
        if(iodva.eq.0)then
            iunit = 2
        else if(iodva.eq.1)then
            iunit = 3
        else if(iodva.eq.2)then
            iunit = 4
        endif
        junit = 11
            DEPTHR = 0
            iftype = 16
            iobsyn = 2
            itmfrq = 1
            cfilt = 'None'
            keyear = 0
            kemon = 0
            keday = 0
            kehour = 0
            kemin = 0
            esec = 0.0
            evlat = -12345
            evlon = -12345
            evdep = depths

            stname = 'GRN16'
            stlat  = -12345
            stlon = -12345
            stelev = depthr
            distkm = dist
            distdg = dist/111.195
            evstaz = 0.0
            stevaz = 180.0
            lmnm = lgstr(mname)
            ccomnt = mname(1:lmnm)
c-----
c     TP TSV and TSH have no meaning in ordeinary sense
c-----
        TP = -12345.
        TSV = -12345.
        TSH = -12345.
        SA = 0.0
        SC = 0.0
        SF = 0.0
        SL = 0.0
        SN = 0.0
        SR = 0.0
c-----
c         Output header
c-----
            call wrhd96(LOT,nerr)
c-----
c     Now that all spectra for a given distance are in
c     make time history and output in desired order
c-----
        np2 = npts/2 + 1
        df = 1.0/(npts*dt)
        do 1300 jk=1,16
            if(jsrc(jk).ne.0)then
                do 1301 i=1,npts
                    x(i) = green(i,jk)
 1301           continue

                if(idva.eq.0)then
                    call minteg(x,npts,dt,.false.)
                else if(idva.eq.2)then
                    call mderiv(x,npts,dt)
                endif
c-----
c         vertical
c-----
                if(iszrt(jk).eq.1)then
                    cmpinc = -90.0
                    cmpaz  =   0.0
c-----
c         radial
c-----
                else if(iszrt(jk).eq.4)then
                    cmpinc = 0.0
                    cmpaz  =   0.0
c-----
c         transverse
c-----
                else if(iszrt(jk).eq.5)then
                    cmpinc = 0.0
                    cmpaz  =  90.0
                endif
                stcomp = ost(jk)
                cmpdt = dt
                ksyear = 0
                ksmon = 0
                ksday = 0
                kshour = 0
                ksmin = 0
                ssec = tfirst
                call wrtr96(LOT,stcomp,cmpinc,cmpaz,
     1              cmpdt, npts, ksyear, ksmon, 
     2              ksday, kshour, ksmin, ssec, 
     3              x,nerr,NSAMP)
            endif
 1300   continue
        return
        end
