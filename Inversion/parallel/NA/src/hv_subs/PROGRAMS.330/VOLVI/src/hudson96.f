        program hudson96
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME VI                                                      c
c                                                                     c
c      PROGRAM: HUDSON96                                              c
c                                                                     c
c      COPYRIGHT 2008                                                 c
c      R. B. Herrmann                                                 c
c      Department of Earth and Atmospheric Sciences                   c
c      Saint Louis University                                         c
c      221 North Grand Boulevard                                      c
c      St. Louis, Missouri 63013                                      c
c      U. S. A.                                                       c
c                                                                     c
c---------------------------------------------------------------------c
c       CHANGES
c       17 MAR 2008 - created this program for the purpose of 
c             synthesizing a teleseismic P-wave signal
c             for a high sample rate. To save time, the
c             approach of Hudson (1969) is used that
c             separates the problem into a source crust,
c             receiver crust and teleseism problems, with the source
c             and receiver crust response fcomputed using propagator
c             matrices.
c
c             This program may seem  complicated, primarily because
c             large sections of the programs hspec96p, hrftn96 and
c             time96 are used as the building blocks of this program.
c             However it does make synthetics fast
c          I recommend using the P wave synthetics from first arrival to
c          the fime of PP or PcP, which ever is earliest.
c
c          I do not recommend the use of Z and R (e.g., SV) S synthetics
c          since too much is missing.
c
c          The SH synthetics can be used up to the arrival of SS, or ScS,
c          whichever is earliest. However do not use part 80 degrees
c          this the timing and geometrical spreading are that of SKS and not S   
c
c       11 APR 2008 - corrected an error that led to incorrect synthetics
c          if the source was at the default 60 km of the source model
c       25 MAY 2008 - SH times are written into the header so that they 
c          appear in the SAC file
c       28 MAY 2008 - to properly set the arrival time, the start of the
c          teleseismic attenuation pulse is determined empirically using two
c          FFT's. This must be done using the same DT as for the synthetic
c          because of the frequency dispersion inherent in the causal T* 
c          operator
c       29 MAY 2008 - corrected indexing when merging source and receiver 
c          structure with teleseism structure. The halfspace was not
c          inserted correctly
c       18 FEB 2009 -  took care to add common/depref/refdep in frstar
c       17 MAR 2009 -  the default offset of P is 10 sec, and
c                      is NOW 20 sec not 10 sec for S
c                      corrected numerical underflow for S-wave incident
c                      at the receiver by using a compound matrix
c       19 MAR 2009 -  added -TSTAR tstar to override the computed T*
c                      to permit the use of the same code, if T* = 0
c                      we use 0.0001 internally
c       07 MAY 2009 -  Add -TTONLY flag to get travel times only
c       23 JUL 2009 - ensured consistent call to FFT routine with
c               complex arguments. We now call zfour
c       11 OCT 2010  - forced extension of source/receiver models to
c                    a predefined maximum depth
c        
c  TO DO - 
c  1. Redo the computation of the ray parameter which 
c     is simplistic. The computation should use the
c     hudsonsrc.mod down and hudsonrec.mod up. or do this
c     iteratively
c TODO     The ray parameter is computed using the source structure and the
c          receiver structure.  It is assumed that the ray turning point is
c          beneath the shallow source and receiver specific structures
c          receiver. This is done by modifying the travel time routine which
c          was common to other codes until this point by using the
c          models stored in the common blocks srcmod recmod
c       
c  2. For the S first arrival do not permit SKS
c  6. redo frstar so that the model read is in the fstarr which
c      then returns the medium properties of the original model
c      This wil make geometrical simpler
c  7. See effect of Earth flattening on p-tau.  Currently for
c     short distance and deep event, the ray parameters of P pP and sP
c     are different and hence the whole validity of things is questioned
c     the flattening is a kludge, but actually works for hspec96?
c     THIS IS REQUIRED
c  8. 
c----- 
        implicit none

c-----
c       input/output default LUN
c-----
        integer LER,LOT,LIN
        parameter(LER=0, LIN=5, LOT=6)
c-----
c       velocity model file
c-----
        integer NL
        parameter (NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL),
     2      frefp(NL), frefs(NL)
        real d,a,b, rho,qa,qb,etap,etas,frefp,frefs
        common/depref/refdep
        real refdep
        real refdepsrc, refdeprec

        integer mmax, iunit, iiso, iflsph, idimen, icnvel
        integer ierr
        character title*80
c-----
c       command line information
c-----
        character modeltel*100
        character modelsrc*100
        character modelrec*100

        real dt, gcarc, hs, offset, utstar
        integer npts

        logical dosrc, dorec, dotel, dop, dokjar, dottonly
c-----
c       teleseismic parameters
c-----
        real hstel, hrtel
        real zsrc, zrec
        real timep, timesv, pvelrec, svelrec, densrec
        real vsa, vsb, vsr
        real raypp, geomp, tstarp
        real rayps, geoms, tstars
        real rayp , geom , tstar 

        real tstaradj

        real TP, TSV, TSH

        real srcdelay, recdelay
c-----
c       Earth parameters
c-----
        common/earth/radius
        real radius
        real kmdeg
c-----
c       receiver region arrays
c-----
        integer MAXPTS
        parameter (MAXPTS=16384)
        complex recz(MAXPTS), recr(MAXPTS), rect(MAXPTS)
        complex srcp(21,MAXPTS)
c-----
c       teleseismic tstar operator
c-----
        complex zts(MAXPTS)
c-----
c       time series parameters
c-----
        real alpha
        real tshft
c-----
c       parameters for output of Green function
c-----
        real SR,  SL, SN, SA, SC, SF
        real VRA, VRB

c-----
c       routine parameters
c-----
        integer lmax
        real dph
        integer iprog
        integer n,n1,n2,n21,nyq2
        real fl, fu, df
        real rr, tt0, tshift, freq, omega, tmp
        complex ztmp, zomega, zpfac
        real dr, di

        real datar, datai
        integer jk, i

        integer lnblnk

c-----
c       ksrc = temporary array for jsrc for output
c       here jsrc != 0 reflects source information
c       if receiver is not in a fluid DO NOT output pressure field
c-----
        integer ksrc(21)
c-----
c       lsrc maps jsrc to output Green s functions. e.g., if
c       jsrc(8) = radial explosion term, but in final output it
c       occupies position 10, or jsrc(lsrc(10)) = computed
c-----
        integer*4 lsrc(21)
        data lsrc/1,2,3,4,13,5,6,14,7,8,9,10,11,12,15,16,17,18,19,20,21/

        

c-----
c       initialize
c-----
        radius = 6371.
        kmdeg=radius*3.1415927/180.0
        do i=1,21
            ksrc(i) = 0
        enddo
c-----
c       parse command line arguments
c-----
        call gcmdln(dt, gcarc, hs, npts,
     1    modeltel,modelsrc,modelrec,offset,
     2    dosrc,dorec,dotel,dop,dokjar,utstar,dottonly,
     3    zsrc,zrec)
       if(dop)then
c-----
c                force P for MT sources
c-----
             ksrc( 1) = 1
             ksrc( 2) = 1
             ksrc( 3) = 1
             ksrc( 4) = 1
             ksrc( 6) = 1
             ksrc( 7) = 1
             ksrc( 9) = 1
             ksrc(10) = 1
        else
             ksrc( 1) = 1
             ksrc( 2) = 1
             ksrc( 3) = 1
             ksrc( 4) = 1
             ksrc( 5) = 1
             ksrc( 6) = 1
             ksrc( 7) = 1
             ksrc( 8) = 1
             ksrc( 9) = 1
             ksrc(10) = 1
        endif
c-----
c       make npts a power of 2, and check the size
c-----
        call npow2(npts)
        if(npts.gt.MAXPTS)npts = MAXPTS
c-----
c       output command line 
c-----
        if(dop)then
             WRITE(LOT,*)'Synthetic seismogram parameters for P'
        else
             WRITE(LOT,*)'Synthetic seismogram parameters for S'
        endif
        WRITE(LOT,*)'Offset               : ',offset
        WRITE(LOT,*)'Source depth (km)    : ',hs
        WRITE(LOT,*)'Arc distance (deg)   : ',gcarc
        WRITE(LOT,*)'npts                 : ',npts
        WRITE(LOT,*)'dt                   : ',dt
        WRITE(LOT,*)'Teleseism model      : ',
     1                                 modeltel(1:lnblnk(modeltel))
        WRITE(LOT,*)'Source region model  : ',
     1                                 modelsrc(1:lnblnk(modelsrc))
        WRITE(LOT,*)'Receiver region model: ',
     1                                 modelrec(1:lnblnk(modelrec))
c-----
c       develop strategy
c       Eventually correct  the delay for very deep sources. Of course we could
c       We will propagate from the MAX(src depth, src crust, 60 km) position in
c       the teleseismic model to MAX( rec crust, 60 km) 
c-----
        call defineprop(modeltel,modelsrc,modelrec,hs,gcarc,
     1      hstel, hrtel, zsrc, zrec)
        WRITE(LOT,*)'Teleseism model used from ',
     1      hstel,' (km) at source'
        WRITE(LOT,*)'Teleseism model used to   ',
     1      hrtel,' (km) at receiver'
c-----
c       get the travel time etc for the teleseismic propagation
c-----
C ORDER IMPORTANT HERE FOR RAYP AND TSTAR
        call frstar(gcarc*111.195,hstel,hrtel,modeltel,2,
     1   TSV, pvelrec, svelrec, densrec, 
     2   vsa, vsb, vsr, rayps, geoms, tstars,
     1   .false.)
        TSH  = TSV
C        WRITE(6,*)'TSV,tstars,rayps,geoms:',TSV,tstars,rayps,geoms
        call frstar(gcarc*111.195,hstel,hrtel,modeltel,1,
     1   TP, pvelrec, svelrec, densrec, 
     2   vsa, vsb, vsr, raypp, geomp, tstarp,
     1   .false.)
C        WRITE(6,*)'TP,tstarp,raypp,geomp:',TP,tstarp,raypp,geomp
        WRITE(LOT,*)'Teleseismic propagation'
        if(dop)then
             WRITE(LOT,*)'P-wave travel time   : ',TP,'(s)'
             WRITE(LOT,*)'Ray parameter        : ',raypp,' (sec/km)'
             WRITE(LOT,*)'Geometrical spreading: ',geomp
             WRITE(LOT,*)'T*(P) from model     : ',tstarp
             rayp = raypp
             tstar = tstarp
             geom = geomp
        else 
             WRITE(LOT,*)'S-wave travel time   : ',TSV,'(s)'
             WRITE(LOT,*)'Ray parameter        : ',rayps,' (sec/km)'
             WRITE(LOT,*)'Geometrical spreading: ',geoms
             WRITE(LOT,*)'T*(S) from model     : ',tstars
             rayp = rayps
             tstar = tstars
             geom = geoms
        endif
        if(utstar.ge.0.0)then
           if(dop)then
             WRITE(LOT,*)'T*(P) used           : ',utstar
           else
             WRITE(LOT,*)'T*(S) used           : ',utstar
           endif
        endif
        WRITE(LOT,*)'P-velocity rec base  : ',pvelrec,' (km/s)'
        WRITE(LOT,*)'S-velocity rec base  : ',svelrec,' (km/s)'
        WRITE(LOT,*)'Density    rec base  : ',densrec,' (gm/cm^3)'
        WRITE(LOT,*)'P-velocity src base  : ',vsa,' (km/s)'
        WRITE(LOT,*)'S-velocity src base  : ',vsb,' (km/s)'
        WRITE(LOT,*)'Density    src base  : ',vsr,' (gm/cm^3)'

c-----
c       define alpha for time domain damping
c-----
        alpha = 2.5/(npts*dt)
c-----
c       get the t* attenuation
c-----
        if(dop)then
              if(utstar.eq.0.0)then     
              call telop(0.0001,TP ,geomp,dt,alpha,npts,zts,dotel,
     1             tstaradj,dokjar,offset)
              else if(utstar.gt.0.0)then     
              call telop(utstar,TP ,geomp,dt,alpha,npts,zts,dotel,
     1             tstaradj,dokjar,offset)
              else
              call telop(tstarp,TP ,geomp,dt,alpha,npts,zts,dotel,
     1             tstaradj,dokjar,offset)
              endif
        else
              if(utstar.eq.0.0)then     
              call telop(0.0001,TSV,geoms,dt,alpha,npts,zts,dotel,
     1             tstaradj,dokjar,offset)
              else if(utstar.gt.0.0)then     
              call telop(utstar,TSV,geoms,dt,alpha,npts,zts,dotel,
     1             tstaradj,dokjar,offset)
              else
              call telop(tstars,TSV,geoms,dt,alpha,npts,zts,dotel,
     1             tstaradj,dokjar,offset)
              endif
        endif
        WRITE(LOT,*)'T* pulse adjustment  : ',tstaradj,' (s)'
c-----
c       Get receiver region response
c-----
        WRITE(LOT,*)'Computing receiver region response'
        call getrec(dt,alpha,rayp,npts,recz,recr,rect,'hudsonrec.mod',
     1      hrtel,recdelay,dorec,pvelrec,svelrec,densrec,dop)
c-----
c       Get source region  input to teleseismic propagation
c-----
        WRITE(LOT,*)'Computing source region input'
        call getsrc(dt,alpha,rayp,npts,srcp,'hudsonsrc.mod',hs,hstel,
     1      srcdelay,dosrc,vsa,vsb,vsr,dop)
        if(dottonly)then
             WRITE(LOT,'(a,1x,f5.1,1x,f5.1,1x,f6.2)')
     1          'TT:',gcarc,HS,TP+srcdelay+recdelay
        else
c-----
c       output synthetics in Green function format
c-----
        open(unit=4,file='hspec96.grn',status='unknown',
     1      form='unformatted',access='sequential')

        rewind 4
c------
c       note I need the A B RHO at source depth and receiver depth
c       not at base of layers as above
c-----
        iprog = 4
        n =npts
        n1 = 1
        n2 = npts/2 + 1
        df = 1./(npts*dt)
        fl = (n1 -1 ) * df
        fu = (n2 -1 ) * df
        nyq2 = 2*(npts/2 + 1)
        
        
        write(4)iprog
        write(4) alpha,fl,fu,dt,n,n1,n2,df,nyq2
        write(4)modeltel
C        WRITE(6,*)'iprog:', iprog
C        WRITE(6,*)'alpha:', alpha
C        WRITE(6,*)'fl:', fl
C        WRITE(6,*)'fu:', fu
C        WRITE(6,*)'dt:', dt
C        WRITE(6,*)'n:', n
C        WRITE(6,*)'n1:', n1
C        WRITE(6,*)'n2:', n2
C        WRITE(6,*)'df:', df
C        WRITE(6,*)'nyq2:', nyq2
c-----
c       safety do an insert
c----
        call getmod(1,modelsrc,mmax,title,iunit,iiso,iflsph,
     1          idimen,icnvel,ierr,.false.)
        if(ierr .lt. 0)stop
        refdepsrc = refdep
            call srclay(hs,lmax,dph)

            VSB = b(lmax)
            VSA = a(lmax)
            VSR = rho(lmax)
        WRITE(LOT,*)'SRC:A,B,R:',VSA,VSB,VSR
c-----
c       define TI constants
c-----
            SR = VSR
            SL = VSR*VSB*VSB
            SN = VSR*VSB*VSB
            SA = VSR*VSA*VSA
            SC = VSR*VSA*VSA
            SF = SA - 2.0*SN
        call getmod(1,modelrec,mmax,title,iunit,iiso,iflsph,
     1          idimen,icnvel,ierr,.false.)
        if(ierr .lt. 0)stop
        refdeprec = refdep
            call srclay(0.0,lmax,dph)
            VRB = b(lmax)
            VRA = a(lmax)
        WRITE(LOT,*)'REC:A,B:',VRA,VRB

c-----
c       TSHFT is the absolute time of the first sample point of the
c       time series. We need this done carefully because of the
c       use of complex frequency (e.g., time damping) for control of
c       singularities. The time of the first sample is
c         -offset + timep + recdelay + srcdelay which accounts for the
c       initial offset, the propagation from the source down to the
c       teleseism model, through the deep earth model, and up from the
c       teleseism model at the receiver to the surface.
c       To do this correctly we must also synthetically delay the
c       teleseism part for propagation. Everything in unraveled in
c       hpulse96
c-----
        if(dop)then
             timep = TP + srcdelay + recdelay
             timesv = -12345.
             TSHFT = -offset + timep
        else
             timep = -12345.
             timesv = TSV + srcdelay + recdelay
             TSHFT = -offset + timesv
        endif
            write(4)gcarc*111.195,tshft,hs-refdepsrc,
     1          0.0-refdeprec,
     2          timep, timesv, timesv,
     3          SA, SC, SF, SL, SN, SR

        TSHFT = TSHFT + tstaradj
C        WRITE(6,*)'TSHFT:',TSHFT
        write(4)ksrc
        do i=n1,n2
                 freq = (i-1)*df
c-----
c                handle the first point as origin time
c-----
                 tmp = 6.2831853*freq*(TSHFT )
                 ztmp = cmplx(cos(tmp),sin(tmp))
       
                 do jk=1,21
                      if(ksrc(jk).eq.1)then
          
                      if(jk.eq.1 .or. jk.eq.3 .or. jk.eq.6.
     1   .or. jk.eq.9)then
                      datar = real(recz(i)*ztmp*zts(i)*srcp(jk,i))
                      datai = aimag(recz(i)*ztmp*zts(i)*srcp(jk,i))
C       WRITE(6,*)'i,freq,jk,r,t,T,s,R,I:',i,freq,jk,recz(i),
C     1     ztmp,zts(i),srcp(jk,i),datar,datai,TSHFT
                      else if(jk.eq.2 .or. jk.eq.4 .or. jk.eq.7.
     1   .or. jk.eq.10)then
                      datar = real(recr(i)*ztmp*zts(i)*srcp(jk,i))
                      datai = aimag(recr(i)*ztmp*zts(i)*srcp(jk,i))
C       WRITE(6,*)'i,freq,jk,r,t,T,s,R,I:',i,freq,jk,recr(i),
C     1     ztmp,zts(i),srcp(jk,i),datar,datai,TSHFT
                      else if(jk.eq.5 .or. jk.eq.8)then
                      datar = real(rect(i)*ztmp*zts(i)*srcp(jk,i))
                      datai = aimag(rect(i)*ztmp*zts(i)*srcp(jk,i))
C       WRITE(6,*)'i,freq,jk,r,t,T,s,R,I:',i,freq,jk,rect(i),
C     1     ztmp,zts(i),srcp(jk,i),datar,datai,TSHFT
                      endif
                      write(4)datar,datai
            
             endif
                  enddo
        enddo
        rr = -1.0
        tt0 = 0.0
        write(4)rr,tt0

        close(4)
c-----
c       
c-----
        WRITE(LOT,*)'Program Hudson completed'
        endif
        end

        subroutine telop(tstar,T,geom,dt,alpha,npts,zts,dotel,
     1     tstaradj,dokjar,offset)
c-----
c	compute zts, the Fourier transform of the teleseism
c       response to include the time shift due to propagation and
c       the anelastic attenuation operation.
c
c       We follow 
c       Kjartansson, E., Constant Q-wave propagation and
c             Attenuation, J. Geophys. Res., Vol 84, 4737-4748, 1979
c
c       tstar     R*4    - T*
c       T         R*4    - travel time, required for Kjartansson operator
c       geom      R*4    - teleseismic geometrical spreading factor
c       dt        R*4    - samling interval
c       alpha     R*4    - time domain attenuation operator
c       npts      I*4    - number of data points
c       zts       C*4    - array of complex attenuation operators
c       dotel     L      - .true. comp[ute complex T* operator
c       tstaradj  R*4    - adjustment to get initial pulse
c       dokjar    L      - .true. use Kjartasson operator, else use
c                          Futterman
c       offset    R      - source pulse offset - note this is
c                          use here only for the telop correction
c----- 
        implicit none

        integer MAXPTS
        parameter (MAXPTS=16384)
        real tstar, T, dt, alpha, geom, tstaradj
        integer npts
        complex zts(MAXPTS)
        logical dotel, dokjar
        real offset 
c-----
c       internal variables
c-----
        integer n21, i
        real freq, df
        complex ztmp
        real qi
        real fac, dfac, depmax
        integer ifound

        complex z(MAXPTS)
        real x(MAXPTS)
        real x0, y0, x1, y1, xi


        n21 = npts/2 + 1
        df = 1./(npts*dt)
        if(tstar.gt.0.0)then
        qi = tstar/T
        else
        qi = 0.0
        endif


        if(.not. dotel)then
c-----
c       turn off for testing
c-----
             do i=1,n21
                  zts(i) = cmplx(1.0,0.0)
             enddo
             T = 0.0
             tstaradj = 0.0
        else
              do i=1,n21
                 freq = (i-1)*df
                 if(freq.lt.df)freq = 0.01 * df
                 call causlq(freq,alpha,qi,T,ztmp,dokjar)
                 zts(i) = ztmp*geom
              enddo
c-----
c       compute the time delay
c           define a zero phase source pulse
c-----
            do i=1,npts
                z(i) = cmplx(0.0,0.0)
                if(i.eq.1)then
                    z(i) = cmplx(0.50/dt,0.0)
                else if(i.eq.2)then
                    z(i) = cmplx(0.25/dt,0.0) * exp(-alpha*dt)
                else if(i.eq.npts)then
                    z(i) = cmplx(0.25/dt,0.0) * exp(+alpha*dt)
                endif
            enddo
            call zfour(z,npts,-1,dt,df)

              do i=1,n21
                 freq = (i-1)*df
                 if(freq.lt.df)freq = 0.01 * df
                 call causlq(freq,alpha,qi,T,ztmp,dokjar)
                 z(i) = z(i) * ztmp
              enddo
c-----
c            apply the time offset
c-----
             do i=1,n21
                 freq = (i-1)*df
                 fac = 6.2831853*freq*(-offset+T)
                 ztmp = cmplx(cos(fac), sin(fac) )
                 z(i) = z(i) * ztmp
             enddo
c-----
c            force a real time series
c-----
             do i=1,n21
                  if(i.gt.1)then
                      z(npts -i + 2) = conjg(z(i))
                  endif
             enddo
             call zfour(z,npts,+1,dt,df)
c-----
c            undamp
c-----
             fac = exp(alpha*offset)
             dfac = exp(alpha*dt)
             depmax = -1.0e+37
             do i=1,npts
                   z(i) = z(i) * fac
                   fac = fac * dfac
                   x(i) = real(z(i))
                   if(x(i).gt.depmax)then
                       depmax = x(i)
                   endif
             enddo
c-----
c            attempt to get the initial time by finding the first point
c                   >= depmax/100.0
c-----
            tstaradj = 0
            ifound = 0
            tstaradj = 0.0
            do i=1,npts-1
               if(x(i+1).ge.0.01*depmax .and. ifound.eq.0)then
                     ifound = 1
c-----
c                    do a 2 point interpolation
c-----
                     x0 = i
                     y0 = x(i)
                     x1 = i+1
                     y1 = x(i+1)
                     if(y1 .ne. y0)then
                     xi = x0 + ( 0.01*depmax - y0)*( x1-x0)/(y1-y0)
                     tstaradj = - offset + (xi-1)*dt 
                     endif
C        WRITE(6,*)'xi,offset,tstaradj,dt:',xi,offset,tstaradj,dt
C        WRITE(6,*)'depmax:',depmax
C        WRITE(6,*)'x0    :',x0    
C        WRITE(6,*)'y0    :',y0    
C        WRITE(6,*)'x1    :',x1    
C        WRITE(6,*)'y1    :',y1    
               endif
            enddo
        endif

        
        return
        end

        subroutine getrec(dt,alp,rayp,npts,recz,recr,rect,
     1      modelrec,hrtel,recdelay,dorec,
     2      pvelrec,svelrec,densrec,pincident)
        implicit none
c-----
c       routine calling arguments
c-----
        real dt,alp,rayp,pvelrec,svelrec,densrec
        integer npts
        complex recz(*), recr(*), rect(*)
        character modelrec*(*)
        real hrtel, recdelay
        logical dorec
        logical pincident

c-----
c
c-----
        common/damp/alpha,ieqex
            real alpha
            integer ieqex
c-----
c       ray parameter values
c-----
        common/c/pmin,pmax,dp,pcntrl
        real pmin,pmax,dp,pcntrl
c-----
c       receiver model parameters
c-----
        integer LER, LIN, LOT
        parameter(LER=0, LIN=5, LOT=6)
c-----
c       earth model information
c-----
        integer NL
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL),
     2      frefp(NL), frefs(NL)
        real d, a, b, rho, qa, qb, etap, etas, frefp, frefs
        common/modlly/mmax
        integer mmax
        common/depref/refdep
        real refdep

        integer iunit, iiso, iflsph, idimen, icnvel, ierr
        character title*80
        common/lyrctl/lyrins
        logical lyrins
c-----
c       internal variables
c-----
        integer i, n21
        real freq, df
        complex Z, R, T
        integer lmaxr
        real dphr
        real eta
 
        alpha = alp
        pmin = rayp
        pmax = rayp
        dp = rayp
        pcntrl = -1.0
        lyrins = .true.

        n21 = npts/2 + 1
        df = 1.0/(npts*dt)
c-----
c	get the model at the receiver
c-----
        call getmod(1,modelrec,mmax,title,iunit,iiso,iflsph,
     1          idimen,icnvel,ierr,.false.)
        if(ierr .lt. 0)stop
c-----
c       make sure that we use 1/Q
c-----
        do i=1,mmax
            if(qa(i).lt.0.0)qa(i) = 0.0
            if(qb(i).lt.0.0)qb(i) = 0.0
            if(qa(i) .gt. 1.0)qa(i) = 1.0/qa(i)
            if(qb(i) .gt. 1.0)qb(i) = 1.0/qb(i)
            if(frefp(i) .le. 0.0)frefp(i) = 1.0
            if(frefs(i) .le. 0.0)frefs(i) = 1.0
        enddo
c-----
c       insert a boundary at hrtel and determine the index of the layer here
c       output for safety as a debug
c-----
        call insert(hrtel+refdep)
        call dezero()
        call srclyr(hrtel+refdep,lmaxr,dphr)
c-----
c       Change the maximum number of layers and
c       place teleseism model parameters into new base
c       layer - note in future preserve the Q etc
c-----
        mmax = lmaxr
        a(mmax)   = pvelrec
        b(mmax)   = svelrec
        rho(mmax) = densrec
c-----
c       should also carry the Q etc
c-----
        WRITE(LOT,'(a)')
     1     '  LAYER   H(km)      PVel      SVel      Dens'
        do i=1,lmaxr
              if(i.lt.lmaxr)then
              WRITE(LOT,'(i5,4f10.2)')i,d(i),a(i),b(i),rho(i)
              else
              WRITE(LOT,'(i5,4f10.2)')-i,d(i),a(i),b(i),rho(i)
              endif
        enddo
c-----
c	compute the vertical tau delay
c-----
        recdelay = 0.0
        do i=1,lmaxr-1
           if(pincident)then
                eta = sqrt(1.0 - rayp*rayp*a(i)*a(i))
                recdelay = recdelay + d(i)*eta/a(i)
           else
                eta = sqrt(1.0 - rayp*rayp*b(i)*b(i))
                recdelay = recdelay + d(i)*eta/b(i)
           endif
        enddo
        WRITE(LOT,*)'recdelay   :',recdelay
    

c-----
c       compute the response renaming the excit of hrftn96 to exitpw
c       for incident P
c-----
        do i=1,n21
             freq = (i-1)*df
             if(freq.lt.df) freq = 0.01*df
             call excitpw(freq,pincident,rayp,Z,R,T,pvelrec,svelrec)
C        IF(I.EQ.2648)THEN
C          WRITE(6,*)'Z,R,T:',Z,R,T
C        ENDIF
             recz(i) = Z
             recr(i) = R
             rect(i) = T
        enddo
c-----
c       turn off for testing
c-----
        if(.not. dorec)then
             recdelay = 0.0
             do i=1,n21
                   recz(i) = cmplx(1.0,0.0)
                   recr(i) = cmplx(1.0,0.0)
                   rect(i) = cmplx(1.0,0.0)
             enddo
        endif
        return
        end

        subroutine defineprop(modeltel,modelsrc,modelrec,hs,gcarc,
     1      hstel, hrtel, zsrc, zrec)
        implicit none
c-----
c       routine parameters
c-----
        character modeltel*(*)
        character modelsrc*(*)
        character modelrec*(*)
        real hs, gcarc
        real hstel, hrtel, zsrc, zrec
c-----
c       velocity model parameters
c------
        integer NL
        parameter (NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL),
     2      frefp(NL), frefs(NL)
        real d,a,b, rho,qa,qb,etap,etas,frefp,frefs
        common/depref/refdep
        real refdep

        integer mmax, iunit, iiso, iflsph, idimen, icnvel
        integer ierr
        character title*80

        real layerdepth(3*NL)

c-----
c       internal variables
c-----
        real depth
        integer i,j

c-----
c       determine the base of the receiver model if not modeltel
c       also it must be a minimum of zrec km deep which should be upper mantle
c-----
        if(modelrec.eq.modeltel)then
             hrtel = zrec
        else
             call getmod(1,modelrec,mmax,title,iunit,iiso,iflsph,idimen,
     1            icnvel, ierr,.false.)
c-----
c            compute the depth to the base of the model
c-----
             depth = 0.0
             do i=1,mmax
                 depth = depth + d(i)
             enddo
             if(depth.lt.zrec)then
                   hrtel = zrec
             else
                   hrtel = depth 
             endif
        endif
c-----
c       determine the base of the source model if not modeltel
c       also it must be a minimum of zsrc km deep which should be upper mantle
c-----
        if(modelsrc.eq.modeltel)then
             hstel = zsrc
        else
             call getmod(1,modelsrc,mmax,title,iunit,iiso,iflsph,idimen,
     1            icnvel, ierr,.false.)
c-----
c            compute the depth to the base of the model
c-----
             depth = 0.0
             do i=1,mmax
                 depth = depth + d(i)
             enddo
             if(depth.lt.zsrc)then
                   hstel = zsrc
             else
                   hstel = depth
             endif
        endif
        if(hstel.le.hs)then
             hstel = hs + zsrc
        endif
c-----
c       build the hudsonsrc.mod hudsonrec.mod model files
c       These must be spherical
c       These must have the same number of layers
c       read all models, compute depths, sort, make uniq
c
c   beware of refdep
c-----
        do i=1,3*NL
             layerdepth(i) = 0.0
        enddo
        j=1
c-----
c            get boundaries for teleseism model
c-----
             call getmod(1,modeltel,mmax,title,iunit,iiso,iflsph,idimen,
     1            icnvel, ierr,.false.)
             depth = 0
             do i=1,mmax 
                   layerdepth(j) = depth
                   j = j + 1
                   depth = depth + d(i)
             enddo
             layerdepth(j) = depth
c-----
c            get boundaries for receiver shallow model
c-----
             call getmod(1,modelrec,mmax,title,iunit,iiso,iflsph,idimen,
     1            icnvel, ierr,.false.)
             depth = 0
             j = j + 1
             do i=1,mmax 
                   layerdepth(j) = depth
                   j = j + 1
                   depth = depth + d(i)
             enddo
             layerdepth(j) = depth
             j = j + 1
             layerdepth(j) = zsrc
             j = j + 1
             layerdepth(j) = zrec
c-----
c            get boundaries for source shallow model
c-----
             call getmod(1,modelsrc,mmax,title,iunit,iiso,iflsph,idimen,
     1            icnvel, ierr,.false.)
             depth = 0
             j = j + 1
             do i=1,mmax 
                   layerdepth(j) = depth
                   j = j + 1
                   depth = depth + d(i)
             enddo
             layerdepth(j) = depth
c-----
c       now sort and get uniq
c-----
            call bsort(j,layerdepth,+1)

            call uniq(j,layerdepth)
c-----
c           output the layer thicknesses of the model
c-----
c-----
c       now do tel+rec
c       tel first and then top fill rec
c-----
        call modmerge(modeltel,modelsrc,modelrec, j, 
     1          layerdepth, zsrc, zrec)

        return
        end

        subroutine modmerge(modeltel,modelsrc,modelrec,
     1      numbdy,layerdepth, zsrc, zrec)
c-----
c       merge the models to create hudsonsrc.mod and hudsonrec.mod
c-----
c-----
c       velocity model parameters
c------
        integer NL
        parameter (NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL),
     2      frefp(NL), frefs(NL)
        real d,a,b, rho,qa,qb,etap,etas,frefp,frefs
        common/depref/refdep
        real refdep
c-----
c       temporary storage for teleseism model
c-----
        common/telmod/Td(NL),Ta(NL),Tb(NL),Trho(NL),
     1      Tqa(NL),Tqb(NL),Tetap(NL),Tetas(NL),
     2      Tfrefp(NL), Tfrefs(NL)
        real Td,Ta,Tb, Trho,Tqa,Tqb,Tetap,Tetas,Tfrefp,Tfrefs
        real Trefdep
        integer Tmmax
c-----
c       temporary storage for src model
c-----
        common/srcmod/Sd(NL),Sa(NL),Sb(NL),Srho(NL),
     1      Sqa(NL),Sqb(NL),Setap(NL),Setas(NL),
     2      Sfrefp(NL), Sfrefs(NL)
        real Sd,Sa,Sb, Srho,Sqa,Sqb,Setap,Setas,Sfrefp,Sfrefs
        real Srefdep
        integer Smmax
c-----
c       temporary storage for rec model
c-----
        common/recmod/Rd(NL),Ra(NL),Rb(NL),Rrho(NL),
     1      Rqa(NL),Rqb(NL),Retap(NL),Retas(NL),
     2      Rfrefp(NL), Rfrefs(NL)
        real Rd,Ra,Rb, Rrho,Rqa,Rqb,Retap,Retas,Rfrefp,Rfrefs
        real Rrefdep
        integer Rmmax
c-----
c       temporary storage for new model
c-----
        common/tmpmod/Nd(NL),Na(NL),Nb(NL),Nrho(NL),
     1      Nqa(NL),Nqb(NL),Netap(NL),Netas(NL),
     2      Nfrefp(NL), Nfrefs(NL)
        real Nd,Na,Nb, Nrho,Nqa,Nqb,Netap,Netas,Nfrefp,Nfrefs
        real Nrefdep
        integer Nmmax
c-----
c       command line parameters
c-----
        character modeltel*(*), modelsrc*(*), modelrec*(*)
        integer numbdy
        real layerdepth(3*NL), zsrc, zrec
c-----
c       parameters for reading the model
c-----
        integer mmax, iunit, iiso, iflsph, idimen, icnvel
        integer ierr
        character title*80
c-----
c       internal routine  variables
c-----
        integer i
        real dphlow, dphhgh
c-----
c	get the teleseism model and store separately
c-----
        call getmod(1,modeltel,mmax,title,iunit,iiso,iflsph,idimen,
     1            icnvel, ierr,.false.)
        call localmod(mmax,refdep,d,a,b,rho,qa,qb,etap,etas,frefp,frefs,
     1   Tmmax,Trefdep,Td,Ta,Tb,Trho,Tqa,Tqb,Tetap,Tetas,Tfrefp,Tfrefs,
     2   .true.)
c-----
c	get the src model and store separately
c-----
        call getmod(1,modelsrc,mmax,title,iunit,iiso,iflsph,idimen,
     1            icnvel, ierr,.false.)
        if(d(mmax).eq.0.0)d(mmax) = zsrc
        call localmod(mmax,refdep,d,a,b,rho,qa,qb,etap,etas,frefp,frefs,
     1   Smmax,Srefdep,Sd,Sa,Sb,Srho,Sqa,Sqb,Setap,Setas,Sfrefp,Sfrefs,
     2   .true.)
c-----
c	get the rec model and store separately
c-----
        call getmod(1,modelrec,mmax,title,iunit,iiso,iflsph,idimen,
     1            icnvel, ierr,.false.)
        if(d(mmax).eq.0.0)d(mmax) = zrec
        call localmod(mmax,refdep,d,a,b,rho,qa,qb,etap,etas,frefp,frefs,
     1   Rmmax,Rrefdep,Rd,Ra,Rb,Rrho,Rqa,Rqb,Retap,Retas,Rfrefp,Rfrefs,
     2   .true.)
c-----
c       now create the combined src/teleseismic model
c       using ugly linear searches
c
c       since the combined model will have a number of layers >=
c       to any original model, we focus on filling the combined
c       model
c-----
        Nrefdep = Trefdep
        Nmmax = numbdy -1
c-----
        
        do i=1,numbdy-1
          if(layerdepth(i+1).le. zsrc)then
              call getindex(layerdepth(i),layerdepth(i+1),Sd,Tmmax,k)
              if(k.gt.0)then
                Nd(i) = layerdepth(i+1) - layerdepth(i)
                Na(i) = Sa(k)
                Nb(i) = Sb(k)
                Nrho(i) = Srho(k)
                Nqa(i) = Sqa(k)
                Nqb(i) = Sqb(k)
                Netap(i) = Setap(k)
                Netas(i) = Setas(k)
                Nfrefp(i) = Sfrefp(k)
                Nfrefs(i) = Sfrefs(k)
              endif
          else
              call getindex(layerdepth(i),layerdepth(i+1),Td,Tmmax,k)
              if(k.gt.0)then
                Nd(i) = layerdepth(i+1) - layerdepth(i)
                Na(i) = Ta(k)
                Nb(i) = Tb(k)
                Nrho(i) = Trho(k)
                Nqa(i) = Tqa(k)
                Nqb(i) = Tqb(k)
                Netap(i) = Tetap(k)
                Netas(i) = Tetas(k)
                Nfrefp(i) = Tfrefp(k)
                Nfrefs(i) = Tfrefs(k)
              endif
          endif
        enddo
c-----
c   lets look at the result
c-----
C        WRITE(0,*)'MERGED TELESEISM/SRC MODEL TOP'
C        do i=1,20
C        WRITE(0,'(i5,4f10.4)')i,Nd(i),Na(i),Nb(i),NRho(i)
C        enddo
C        WRITE(0,*)'..... ......... ......... ......... .........'
C        do i=Nmmax -19, Nmmax
C        WRITE(0,'(i5,4f10.4)')i,Nd(i),Na(i),Nb(i),NRho(i)
C        enddo
c-----
c       put this into the array for output
c-----
        call localmod(mmax,refdep,d,a,b,rho,qa,qb,etap,etas,frefp,frefs,
     1   Nmmax,Nrefdep,Nd,Na,Nb,Nrho,Nqa,Nqb,Netap,Netas,Nfrefp,Nfrefs,
     2   .false.)
        iflsph = 1
        call putmod(1,'hudsonsrc.mod',Nmmax,
     1       'Merged teleseism/src/model',
     2       iunit,iiso,iflsph,idimen,icnvel,.false.)
c-----
c       now create the combined rec/teleseismic model
c       using ugly linear searches
c
c       since the combined model will have a number of layers >=
c       to any original model, we focus on filling the combined
c       model
c-----
        Nrefdep = Trefdep
        Nmmax = numbdy -1
c-----
        
        do i=1,numbdy-1
          if(layerdepth(i+1).le. zrec)then
              call getindex(layerdepth(i),layerdepth(i+1),Rd,Tmmax,k)
              if(k.gt.0)then
                Nd(i) = layerdepth(i+1) - layerdepth(i)
                Na(i) = Ra(k)
                Nb(i) = Rb(k)
                Nrho(i) = Rrho(k)
                Nqa(i) = Rqa(k)
                Nqb(i) = Rqb(k)
                Netap(i) = Retap(k)
                Netas(i) = Retas(k)
                Nfrefp(i) = Rfrefp(k)
                Nfrefs(i) = Sfrefs(k)
              endif
          else
              call getindex(layerdepth(i),layerdepth(i+1),Td,Tmmax,k)
              if(k.gt.0)then
                Nd(i) = layerdepth(i+1) - layerdepth(i)
                Na(i) = Ta(k)
                Nb(i) = Tb(k)
                Nrho(i) = Trho(k)
                Nqa(i) = Tqa(k)
                Nqb(i) = Tqb(k)
                Netap(i) = Tetap(k)
                Netas(i) = Tetas(k)
                Nfrefp(i) = Tfrefp(k)
                Nfrefs(i) = Tfrefs(k)
              endif
          endif
        enddo
c-----
c   lets look at the result
c-----
C        WRITE(0,*)'MERGED TELESEISM/REC MODEL TOP'
C        do i=1,20
C        WRITE(0,'(i5,4f10.4)')i,Nd(i),Na(i),Nb(i),NRho(i)
C        enddo
C        WRITE(0,*)'..... ......... ......... ......... .........'
C        do i=Nmmax -19, Nmmax
C        WRITE(0,'(i5,4f10.4)')i,Nd(i),Na(i),Nb(i),NRho(i)
C        enddo
c-----
c       put this into the array for output
c-----
        call localmod(mmax,refdep,d,a,b,rho,qa,qb,etap,etas,frefp,frefs,
     1   Nmmax,Nrefdep,Nd,Na,Nb,Nrho,Nqa,Nqb,Netap,Netas,Nfrefp,Nfrefs,
     2   .false.)
        iflsph = 1
        call putmod(1,'hudsonrec.mod',Nmmax,
     1       'Merged teleseism/rec/model',
     2       iunit,iiso,iflsph,idimen,icnvel,.false.)
        return
        end

        subroutine getindex(dlow,dhgh,d,mmax,k)
c-----
c       convert the thickness array d into depths and
c       then find the bounding index for depths between dlow and dhgh
c-----
        implicit none
        real dlow, dhgh, d(*)
        integer mmax, k

        integer i
        real dphl,dphh
        real dmid
        dmid = 0.5*(dlow+dhgh)
        dphl = 0.0
        do i=1,mmax 
          dphh = dphl + d(i)
          if(dmid.ge.dphl .and. dmid.le.dphh)then
              k = i
              return
          endif
          dphl = dphh
        enddo
        k = -1
        return
        end
       

        subroutine localmod(mmax,refdep,d,a,b,rho,qa,qb,etap,etas,
     1   frefp,frefs,
     1   Lmmax,Lrefdep,Ld,La,Lb,Lrho,Lqa,Lqb,Letap,Letas,Lfrefp,Lfrefs,
     2   tolocal)
c-----
c       put from getmod storage to local storage is tolocal == .true.
c       else get
c-----
        integer mmax
        real refdep, d(*), a(*), b(*), rho(*), qa(*), qb(*)
        real etap(*), etas(*), frefp(*), frefs(*)
        integer Lmmax
        real Lrefdep, Ld(*), La(*), Lb(*), Lrho(*), Lqa(*), Lqb(*)
        real Letap(*), Letas(*), Lfrefp(*), Lfrefs(*)

        logical tolocal

c-----
c     copy from getmod input to local storage
c-----
        if(tolocal)then
              Lmmax = mmax
              Lrefdep = refdep
              do i=1,Lmmax
                  Ld(i) = d(i)
                  La(i) = a(i)
                  Lb(i) = b(i)
                  Lrho(i) = rho(i)
                  Lqa(i) = qa(i)
                  Lqb(i) = qb(i)
                  Letap(i) = etap(i)
                  Letas(i) = etas(i)
                  Lfrefp(i) = frefp(i)
                  Lfrefs(i) = frefs(i)
              enddo
        else
c-----
c     copy form local storage to putmod arrays
c     also do QC so that VS != 0, and Qinv is output
c     assuming that Q is always greater than 1
c-----
              mmax = Lmmax
              refdep = Lrefdep
              do i=1,mmax
                  d(i) = Ld(i)
                  a(i) = La(i)
                  if(Lb(i).eq.0.0)then
                       b(i) = 0.0001
                  else
                       b(i) = Lb(i)
                  endif
                  rho(i) = Lrho(i)
                  if(Lqa(i).gt.1.0)then
                       qa(i) = 1.0/Lqa(i)
                  else
                       qa(i) = Lqa(i)
                  endif
                  if(Lqb(i).gt.1.0)then
                       qb(i) = 1.0/Lqb(i)
                  else
                       qb(i) = Lqb(i)
                  endif
                  etap(i) = Letap(i)
                  etas(i) = Letas(i)
                  frefp(i) = Lfrefp(i)
                  frefs(i) = Lfrefs(i)
              enddo
        endif
        return
        end


        subroutine uniq(n,x)
        integer n
        real x(n)
      
        m = 1
        do i=2,n
              if(x(i).ne.x(m))then
                  m = m + 1
                  x(m) = x(i)
              endif
        enddo
        n = m
        return
        end

        subroutine bsort(nn,x,isign)
c-----
c   http://linux.wku.edu/~lamonml/algor/sort/bubble.html
c   transliterated to FORTRAN
c-----
c       do bubble sort.
c       isign=+1  increase   =-1 decrease.
c
        real x(nn)
        real temp
        integer i,j
        do i= nn ,1, -1
             do j=1,i-1,1
             x0 = x(j+1)-x(j)
             if(isign.gt.0)x0 = - x0
             if( x0 .gt.0.0)then
                   temp = x(j)
                   x(j) = x(j+1)
                   x(j+1) = temp
             endif
             enddo
        enddo
      
        return
        end


        subroutine gcmdln(dt, gcarc, hs, npts,
     1    modeltel,modelsrc,modelrec,offset,
     2    dosrc,dorec,dotel,dop,dokjar,utstar,dottonly,zsrc,zrec)
c-----
c       parse command line arguments
c       requires subroutine mgtarg() and function mnmarg()
c-----
        implicit none

        integer LER,LOT,LIN
        parameter(LER=0, LIN=5, LOT=6)

        character modeltel*(*)
        character modelsrc*(*)
        character modelrec*(*)

        real dt, gcarc, hs, offset,utstar,zsrc,zrec
        integer npts
        logical dosrc,dorec,dotel,dop,dokjar,dottonly

        integer lnblnk

c-----
c       internal variables
c-----
        integer i, nmarg
        integer mnmarg
        character*100 name
        logical ext

        logical loffsetdefault

c-----
c       initialize
c-----
        modeltel = ' '
        modelsrc = ' '
        modelrec = ' '

        hs = 0.0
        dt = 1024
        npts = 1
        gcarc = 50
        offset = 10
        dosrc = .true.
        dorec = .true.
        dotel = .true.
        dop = .true.
        dokjar = .true.
        loffsetdefault = .true.
        utstar = -12345.
        dottonly = .false.
        zsrc = 100.0
        zrec =  60.0
   

        nmarg=mnmarg()
        i = 0
   11   i = i + 1
             if(i.gt.nmarg)goto 13
             call mgtarg(i,name)
             if(name(1:4).eq.'-TEL')then
                   i = i + 1
                   call mgtarg(i,modeltel)
             else if(name(1:2).eq.'-S' .and.
     1             name(1:4).ne.'-SRC')then
                   dop = .false.
             else if(name(1:4).eq.'-SRC')then
                   i = i + 1
                   call mgtarg(i,modelsrc)
             else if(name(1:2).eq.'-P')then
                   dop = .true.
             else if(name(1:4).eq.'-REC')then
                   i = i + 1
                   call mgtarg(i,modelrec)
             else if(name(1:3).eq.'-DT')then
                   i = i + 1
                   call mgtarg(i,name)
                   read(name,'(bn,f10.0)')dt
             else if(name(1:3).eq.'-HS')then
                   i = i + 1
                   call mgtarg(i,name)
                   read(name,'(bn,f10.0)')hs
             else if(name(1:2).eq.'-G')then
                   i = i + 1
                   call mgtarg(i,name)
                   read(name,'(bn,f10.0)')gcarc
             else if(name(1:6).eq.'-TSTAR')then
                   i = i + 1
                   call mgtarg(i,name)
                   read(name,'(bn,f10.0)')utstar
             else if(name(1:5).eq.'-ZSRC')then
                   i = i + 1
                   call mgtarg(i,name)
                   read(name,'(bn,f10.0)')zsrc
             else if(name(1:5).eq.'-ZREC')then
                   i = i + 1
                   call mgtarg(i,name)
                   read(name,'(bn,f10.0)')zrec
             else if(name(1:7).eq.'-TTONLY')then
                   dottonly = .true.
             else if(name(1:2).eq.'-O')then
                   i = i + 1
                   call mgtarg(i,name)
                   read(name,'(bn,f10.0)')offset
                   if(offset.le.0.0)then
                         offset = 10.0
                         loffsetdefault = .true.
                   else
                         loffsetdefault = .false.
                   endif
             else if(name(1:2).eq.'-N' .and.
     1      name(1:3).ne.'-NO')then
                   i = i + 1
                   call mgtarg(i,name)
                   read(name,'(bn,i10)')npts
             else if(name(1:6).eq.'-NOSRC' )then
                   dosrc = .false.
             else if(name(1:6).eq.'-NOREC' )then
                   dorec = .false.
             else if(name(1:6).eq.'-NOTEL' )then
                   dotel = .false.
             else if(name(1:2).eq.'-F')then
                   dokjar = .false.
             else if(name(1:2).eq.'-?')then
                call usage(' ')
             else if(name(1:2).eq.'-h')then
                call usage(' ')
             endif

             go to 11
   13   continue
c-----
c       safety checks
c-----
        if(modeltel .eq. ' ')call usage(' ')
        if(modelsrc .eq. ' ')modelsrc = modeltel
        if(modelrec .eq. ' ')modelrec = modeltel
        if (zsrc.le.0.0)zsrc = 100.0
        if (zrec.le.0.0)zrec =  60.0
        inquire(file = modeltel, exist=ext)
        if( .not. ext)then
            write(LER,*)' Model file is not valid: ',
     1      modelrec(1:lnblnk(modeltel))
            call usage(' ')
        endif
        inquire(file = modelsrc, exist=ext)
        if( .not. ext)then
            write(LER,*)' Source region file is not valid: ',
     1      modelrec(1:lnblnk(modeltel))
            call usage(' ')
        endif
        inquire(file = modelrec, exist=ext)
        if( .not. ext)then
            write(LER,*)' Receiver region file is not valid: ',
     1      modelrec(1:lnblnk(modeltel))
            call usage(' ')
        endif
         

        if(hs.lt.0.0 .or. hs.gt.750)then
            call usage(' Only depth between 0 - 750 km permitted')
        endif
        if( loffsetdefault .and. .not. dop )then
            offset = 20.0
        endif
        
        return
        end

        subroutine usage(ostr)
c------
c       write out program syntax
c-----
        character ostr*(*)
        integer LER,LOT,LIN
        parameter(LER=0, LIN=5, LOT=6)
        if(ostr.ne. ' ' )then
            lostr = lgstr(ostr)
            write(LER,*)ostr(1:lostr)
        endif
        write(LER,*)'USAGE: ',
     1  'hudson96 [-TEL modeltel] [-SRC modelsrc] [-REC modlrec] ',
     1      ' [ -HS hs ] [ -O offset ] [-P] [-S]',
     2      ' [-GCARC gcarc] [-DT dt] [-NPTS npts] ',
     3      ' [-NOSRC] [-NOREC] [-NOTEL] [-TSTAR tstar]',
     4      ' [-TTONLY] [-HSRC hsrc] [-HREC hrec] [-?] [-h]'
        write(LER,*)
     1  '-TEL modeltel (default none) Earth velocity model'
        write(LER,*)
     1  '-SRC modelsrc (default modeltel) Velocity model for source'
        write(LER,*)
     1  '-REC modelrec (default modeltel) Velocity model for receiver'
        write(LER,*)
     1  '-HS hs        (default 10 km) Source depth '
        write(LER,*)
     1  '-GCARC gcarc  (default 30 degrees) Arc distance '
        write(LER,*)
     1  '-NPTS npts    (default 1024) Number of points in time series'
        write(LER,*)
     1  '-DT dt        (default 0.05 sec) Sampling interval'
        write(LER,*)
     1  '-O  offset    (default 10.0 sec) Time offset ',
     1  'before signal (use more for S) '
        write(LER,*)
     1  '-P            (default true)  make P synthetic '
        write(LER,*)
     1  '-S            (default false)  make S synthetic '
        write(LER,*)
     1  '-NOSRC        (default false) Turn off source input'
        write(LER,*)
     1  '-NOREC        (default false) Turn off receiver part'
        write(LER,*)
     1  '-F            (default false) Use Futterman operator'
        write(LER,*)
     1  '                         else use Kjartansson (1979)'
        write(LER,*)
     1  '-NOTEL        (default false) Turn off teleseism ',
     2  ' geometrical'
        write(LER,*)
     1  '-TSTAR tstar  (default -12345.) If >=0 use this T* instead',
     2  ' of one computed from velocity model'
        write(LER,*)
     1  '-TTONLY       (default false) Output only travel times'
        write(LER,*)
     1  '                     spreading and attenuation'
        write(LER,*)
     1  '-ZSRC zsrc     (default see below',
     2     ' Maximum thickness of Src model merged with teleseismic'
        write(LER,*)
     1  '-ZREC zrec     (default see below',
     2     ' Maximum thickness of Rec model merged with teleseismic'
        write(LER,*) '    If halfspace of Src/Rec model has thickness',
     1     ' of 0 km, the halfspace is extended to the the depth of ',
     2     ' zsrc/zrec km. The default of 100 km for the source and ',
     3     ' 60 km for the receiver models.' 
        write(LER,*)
     1  '-?            Display this usage message'
        write(LER,*)
     1  '-h            Display this usage message'
        write(LER,*)
     1  '-----------------------------------------------'
        write(LER,*)
     1  'The program creates the file hspec96.grn for hpulse96.'
        write(LER,*)
     1  'The program creates the files hudsonsrc.mod and',
     2  ' hudsonrec.mod,'
        write(LER,*)
     1  'which are the complete models (model96 format) used ',
     2  'for the src and rec regions from the surface to interior.'
        write(LER,*)
     1  'these can be plotted using shwmod96.'
        write(LER,*)


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

        subroutine causlq(freq,alpha,qi,tp,ztmp,dokjar)
c-----
c       Kjartansson, E. (1979).
c          Constant Q-wave propagation and attenuation,
c       J. Geophys. Res. 84, 4737-4748.
c-----
            real freq, tp,  qi
            complex Ztmp
            logical dokjar

            real gam, fac, om1, fref
            real pi, pi2
            complex atn, omega, zfac
            real zr, zi
            complex at

            real ts

            ts = tp * qi
            omega = cmplx(6.2831853*freq, -alpha)
            fref = 1.0
            om1 = 6.2831853*fref
            pi = 3.1415927
            pi2 = pi / 2.

            

            if(dokjar)then

                   gam = atan(qi)/pi
                   if(gam.le.0.0)then
                      atn = cmplx(1.0,0.0)
                   else
                       fac = pi2*gam
                       rfac = sin(fac)/cos(fac)
                       atn =
     1                     (omega/om1)**dble(1.0-gam)  *
     2                     dcmplx(dble(rfac),1.0d+00)
                   endif
c-----
c           form the propagation term
c-----
                   
                   zfac = - om1*tp*atn 
            else
                 zfac =  omega*ts*clog(omega/om1)/pi + 
     1             cmplx(0.0,1.0)*omega*ts/2
                 zfac = zfac * cmplx(0.0,1.0) 
            endif
            zr = real(zfac)
            zi = aimag(zfac)
            if(zr.gt. -88)then
                 ztmp = exp(zr)*cmplx(cos(zi),sin(zi))
            else
                 ztmp = cmplx(0.0,0.0)
            endif
        return
        end

        subroutine getsrc(dt,alp,rayp,npts,srcp,modelsrc,hs,hstel,
     1      srcdelay,dosrc,vsa,vsb,vsr,pteleseismic)
        implicit none
c-----
c       routine calling arguments
c-----
        real dt,alp,rayp,vsa,vsb,vsr
        integer npts
        complex srcp(21,npts)
        character modelsrc*(*)
        real hs, hstel, srcdelay
        logical dosrc
        logical pteleseismic
c-----
c
c-----
        common/damp/alpha,ieqex
            real alpha
            integer ieqex
        common/jout/jsrc(21) , jbdrys, jbdryh
            integer jsrc, jbdrys, jbdryh
c-----
c       ray parameter values
c-----
        common/c/pmin,pmax,dp,pcntrl
        real pmin,pmax,dp,pcntrl
c-----
c       receiver model parameters
c-----
        integer LER, LIN, LOT
        parameter(LER=0, LIN=5, LOT=6)
c-----
c       earth model information
c-----
        integer NL
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL),
     2      frefp(NL), frefs(NL)
        real d, a, b, rho, qa, qb, etap, etas, frefp, frefs
        common/modlly/mmax
        integer mmax
        common/depref/refdep
        real refdep

        integer iunit, iiso, iflsph, idimen, icnvel, ierr
        character title*80

        common/lyrctl/lyrins
        logical lyrins
c-----
c       control of wavefield
c-----
        common/ctlrud/rpup,rpdn,rsup,rsdn,dorud
        logical rpup, rpdn, rsup, rsdn, dorud
        common/ctlsud/spup,spdn,ssup,ssdn,dosud
        logical spup, spdn, ssup, ssdn, dosud

        common/cntrl/ishank,hnkarg,dstcor,dokjar
        integer*4 dstcor
        real*4 hnkarg
        logical ishank,dokjar

        common/earth/radius
        real radius

c-----
c       internal variables
c-----
        integer i, n21, j
        real freq, df
        complex Z, R
        integer lmaxb, lmaxs
        real dphr
        real eta
        complex gg(21)
        real sp, cp
        complex zstapha
 
c-----
c       set the parameters normally used by hspec96p
c-----
        alpha = alp
        pmin = rayp
        pmax = rayp
        dp = rayp
        pcntrl = -1
        ieqex = 2
        jbdrys = 1
        jbdryh = 0
        lyrins = .true.

        dokjar = .true.
        if(pteleseismic)then
             rpdn = .true.
             rsdn = .false.
        else
             rsdn = .true.
             rpdn = .false.
        endif
        rsup = .false.
        rpup = .false.
        dorud = .true.
        dosud = .false.


        n21 = npts/2 + 1
        df = 1.0/(npts*dt)
c-----
c	get the model at the source
c-----
        call getmod(1,modelsrc,mmax,title,iunit,iiso,iflsph,
     1          idimen,icnvel,ierr,.false.)
        if(ierr .lt. 0)stop
c-----
c       make sure that we use 1/Q
c-----
        do i=1,mmax
            if(qa(i).lt.0.0)qa(i) = 0.0
            if(qb(i).lt.0.0)qb(i) = 0.0
            if(qa(i) .gt. 1.0)qa(i) = 1.0/qa(i)
            if(qb(i) .gt. 1.0)qb(i) = 1.0/qb(i)
            if(frefp(i) .le. 0.0)frefp(i) = 1.0
            if(frefs(i) .le. 0.0)frefs(i) = 1.0
        enddo
        call modcpy(.true.)
c-----
c       check model for inconsistencies
c-----
        call chkmod()
c-----
c       insert a boundary at hstel and determine the index of the layer here
c       output for safety as a debug
c-----
        call insert(hstel+refdep)
c-----
c       get the index for the source depth
c-----
        call insert(hs+refdep)
        call dezero()
c-----
c           check whether neighboring layers are identical
c           to avoid redundant evaluation
c----
        call equlck()
c-----
c       determine position of source and receive in new layers
c-----
        call srclyr(hstel+refdep,lmaxb,dphr)
        call srclyr(hs+refdep,lmaxs,dphr)
c-----
c       Change the maximum number of layers and
c       place teleseism model parameters into new base
c       layer - note in future preserve the Q etc
c-----
        mmax = lmaxb
        a(mmax)   = vsa
        b(mmax)   = vsb
        rho(mmax) = vsr
c-----
c       should also carry the Q etc
c-----
        WRITE(LOT,'(a)')
     1     '  LAYER   H(km)      PVel      SVel      Dens'
        do i=1,lmaxb
              if(i.ge.lmaxs .and. i.le.lmaxb)then
              WRITE(LOT,'(i5,4f10.2)')i,d(i),a(i),b(i),rho(i)
              else
              WRITE(LOT,'(i5,4f10.2)')-i,d(i),a(i),b(i),rho(i)
              endif
        enddo
c-----
c	compute the vertical tau delay
c-----
        srcdelay = 0.0
        do i=lmaxs,lmaxb-1
           if(pteleseismic)then
                eta = sqrt(1.0 - rayp*rayp*a(i)*a(i))
                srcdelay = srcdelay + d(i)*eta/a(i)
           else
                eta = sqrt(1.0 - rayp*rayp*b(i)*b(i))
                srcdelay = srcdelay + d(i)*eta/b(i)
           endif
        enddo
        WRITE(LOT,*)'srcdelay   :',srcdelay
        if(pteleseismic)then
        sp = rayp*a(lmaxb)
        cp = sqrt(1.0 -sp*sp)
        else
        sp = rayp*b(lmaxb)
        cp = sqrt(1.0 -sp*sp)
        endif

c-----
c       now that we have the material properties at the source and
c       receiver depths, apply the sphericity correction
cDEBUG
c-----
c-----
c       transform the spherical model to a flat model
c-----
        if(iflsph.ne.0)then
            call adosph()
        endif


    

c-----
c       compute the response renaming the excit of hrftn96 to exitpw
c       for incident P
c-----
        do i=1,n21
             freq = (i-1)*df
             if(freq.lt.df) freq = 0.01*df
c-----
c            get the stationary phase factor
c            i omega/V cos theta
c-----
             if(pteleseismic)then
                  zstapha = cmplx(0.0,1.0)*
     1               cmplx(6.2831853*freq,-alpha)*cp/a(lmaxb)
             else
                  zstapha = cmplx(0.0,1.0)*
     1               cmplx(6.2831853*freq,-alpha)*cp/b(lmaxb)
             endif
c-----
c            get the medium response
c-----
             call excitsrc(freq,rayp,gg,hs,hstel,lmaxs,lmaxb)
c-----
c	CONVERT TO P and SV FROM Z R
c-----
             do j=1,16
                srcp(j,i) = cmplx(0.0,0.0)
             enddo
             if(pteleseismic)then
c-----
c                 DD
c-----
                  srcp(1,i)  = (gg(1)*(-cp) + gg(2)*sp)*zstapha
                  srcp(2,i)  = (gg(1)*(-cp) + gg(2)*sp)*zstapha
c-----
c                 DS
c-----
                  srcp(3,i)  = (gg(3)*(-cp) + gg(4)*sp)*zstapha
                  srcp(4,i)  = (gg(3)*(-cp) + gg(4)*sp)*zstapha
c-----
c                 SS
c-----
                  srcp(6,i)  = (gg(5)*(-cp) + gg(6)*sp)*zstapha
                  srcp(7,i)  = (gg(5)*(-cp) + gg(6)*sp)*zstapha
c-----
c                 EX
c-----
                  srcp(9,i)  = (gg(7)*(-cp) + gg(8)*sp)*zstapha
                  srcp(10,i) = (gg(7)*(-cp) + gg(8)*sp)*zstapha
             else
c-----
c                 DD
c-----
                  srcp(1,i)  = (gg(1)*(sp) + gg(2)*cp)*zstapha
                  srcp(2,i)  = (gg(1)*(sp) + gg(2)*cp)*zstapha
c-----
c                 DS
c-----
                  srcp(3,i)  = (gg(3)*(sp) + gg(4)*cp)*zstapha
                  srcp(4,i)  = (gg(3)*(sp) + gg(4)*cp)*zstapha
                  srcp(5,i)  =                   gg(13)*zstapha
c-----
c                 SS
c-----
                  srcp(6,i)  = (gg(5)*(sp) + gg(6)*cp)*zstapha
                  srcp(7,i)  = (gg(5)*(sp) + gg(6)*cp)*zstapha
                  srcp(8,i)  =                   gg(14)*zstapha
c-----
c                 EX
c-----
                  srcp(9,i)  = (gg(7)*(sp) + gg(8)*cp)*zstapha
                  srcp(10,i) = (gg(7)*(sp) + gg(8)*cp)*zstapha
             endif
        enddo
c-----
c       turn off for testing
c-----
        if(.not. dosrc)then
             srcdelay = 0.0
             do j=1,21
                  do i=1,n21
                   srcp(j,i) = cmplx(1.0,0.0)
                  enddo
             enddo
        endif
        return
        end

        subroutine insert(dph)
        implicit none
        real dph
        integer LER, LIN, LOT
        parameter (LER=0, LIN=5, LOT=6)
        integer NL
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL),
     2      frefp(NL), frefs(NL)
        real d, a, b, rho, qa, qb, etap, etas, frefp, frefs
        common/modlly/mmax
        integer mmax
        common/shwave/bsh(NL), rhosh(NL), qbsh(NL)
        real bsh, rhosh, qbsh

        real dp, dphh, hsave, dep
        integer m, ls
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
                a(m+1) = a(m)
                b(m+1) = b(m)
                rho(m+1) = rho(m)
                qa(m+1) = qa(m)
                qb(m+1) = qb(m)
                etap(m+1) = etap(m)
                etas(m+1) = etas(m)
                frefp(m+1) = frefp(m)
                frefs(m+1) = frefs(m)
                bsh(m+1) = bsh(m)
                qbsh(m+1) = qbsh(m)
                rhosh(m+1) = rhosh(m)
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

        subroutine frstar(r,hs,hr,mname,ipsvsh,time,pvel,svel,den,
     1      vsa, vsb, vsr, rayp, geom, tstar, dolock)
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
c       pvel    R   Velocity of P wave at receiver
c       svel    R   Velocity of S wave at receiver
c       den     R   Density at receiver
c       vsa R   P-wave velocity at source
c       vsb R   S-wave velocity at source
c       vsr R   Density at source
c       rayp R   Ray parameter in sec/km
c       geom R   geometrical spreading factor
c       tstar R   geometrical spreading factor
c       dolock L .true. apply locked mode which means to ignore the
c                bottom layer
c-----
        real r, hs, hr, time, pvel, svel, vsa, vsb, vsr
        real rayp, geom, tstar
        logical dolock
        character mname*(*)
        integer ipsvsh
c-----
c-----
c       internal variables
c-----
        real depths, depthr
        real dphs, dphr, dphref
        integer lmaxs, lmaxr, lmaxref
        real tm, rp(2),ts
        real dpdx, deginrad, sinis, sinir, cosir, cosis, sindeg
        real vs, vr, rs, rr

        
        common/earth/radius
        real radius

        common/depref/refdep
        real refdep
c-----
c       compute the travel time
c-----
        call fstarr(r,time,lmaxs, lmaxr, lmaxref,
     1      hs, hr, ipsvsh,iflsph, rayp,
     2      tstar, dolock,
     3      mname, varec, vbrec, rhorec,
     4      vasrc, vbsrc, rhosrc)

        svel = vbrec
        pvel = varec
        den  = rhorec
        vsa  = vasrc
        vsb  = vbsrc
        vsr  = rhosrc
c-----
C       compute the geometrical spreading
C
C       since a flat earth is always used, we use the flat layered
C       geometrical spreading from Officer
C       The geometrical spreading is dimensionless and gives the decrease in amplitude from
C       a distance of 1 km from the source to the receiver
C                        2                            2
C            ( rhos vs  a sin Is Vs                  d T      )
C       sqrt |  -------------------------------     -----     |
C            |                                         2      |
C            ( rhor vr sin DEG  cos Ir Rs Cos Is     dx       )
C
C       where p (sec/km) = dT/dx
C       a = radius of sphere about source - we use 1.0 km
C       Is= incident angle at source
C       Ir= incident angle at receiver
C       Rs= distance from center of Earth to source
C       Rr= distance from center of Earth to receiver
C       DEG=epicental distance in degrees
C       rhos and vs are the density and wave velocity at the source depth
c       rhor and vr are the density and wave velocity at the receiver depth
c
c       To get the dp/dx, we determine p at different distances, and then form
c       an interpolating polynomial
c
              call fstarr(r-500,tm,lmaxs, lmaxr, lmaxref,
     1            hs+refdep, hr+refdep, ipsvsh,iflsph, rp(1),
     2            ts, dolock,
     3            mname, varec, vbrec, rhorec,
     4            vasrc, vbsrc, rhosrc)

              call fstarr(r+500,tm,lmaxs, lmaxr, lmaxref,
     1            hs+refdep, hr+refdep, ipsvsh,iflsph, rp(2),
     2            ts, dolock,
     3            mname, varec, vbrec, rhorec,
     4            vasrc, vbsrc, rhosrc)


              dpdx = abs(rp(1) - rp(2))/(500.0 - ( - 500.0))
              deginrad = r/radius
              sindeg = sin(deginrad)
  
              if(ipsvsh.eq.1)then
                  vs = vsa
                  vr = pvel
                  rs = vsr
                  rr = den
              else if(ipsvsh.eq.2)then
                  vs = vsb
                  vr = svel
                  rs = vsr
                  rr = den
              else if(ipsvsh.eq.3)then
                  vs = vsb
                  vr = svel
                  rs = vsr
                  rr = den
              else if(ipsvsh.eq.4)then
                  vs = vsa
                  vr = pvel
                  rs = vsr
                  rr = den
              else if(ipsvsh.eq.5)then
                  vs = vsb
                  vr = pvel
                  rs = vsr
                  rr = den
              endif
              sinis = rayp*vs
              cosis = sqrt(1.0-sinis*sinis)
              sinir = rayp*vr
              cosir = sqrt(1.0-sinir*sinir)

              fac = (rs*vs*sinis*vs*dpdx)/
     1           (rr*vr*sindeg*cosir*(radius-hs)*cosis)
              geom = sqrt(abs(fac))
              

c-----
        return
        end

        subroutine fstarr(dist,tfirst,lmaxs,lmaxr,lmaxref,
     1      hs,hr,ipsvsh,iflsph, rayp,
     2      tstar, dolock,
     3      mname, varec, vbrec, rhorec,
     4      vasrc, vbsrc, rhosrc)
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
c       hs      R   - depth of source
c       hs      R   - depth of receiver
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
c       mname   Ch* - name of the velocity model
c       varec   R - P velocity at receiver (untransformed)
c       vbrec   R - S velocity at receiver (untransformed)
c       rhorec  R - Density at receiver (untransformed)
c       vasrc   R - P velocity at source (untransformed)
c       vbsrc   R - S velocity at source (untransformed)
c       rhosrc  R - Density at source (untransformed)
c-----
c       since this routine is to be used for omega-k,
c       we will approximate the direct arrival
c
c       18 JAN 2008 - everything is straightforward. The addition of
c          the request for pP and sP changes the logic in that
c          the direct arrival is ignored, and that the upgoing refraction 
c          from the source is ignored. We handle this by just setting
c          a very large tfirst before trying to do the modified 
c          downward path refraction to avoid another level of
c          if/then/else/endif
c       24 MAR 2008 - modified to place the model read into this
c          routine instead of in frstar
c-----
        real dist, tfirst, dphs, dphr, hs, hr, depths, depthr
        real rayp
        real varec, vbrec, rhorec
        real vasrc, vbsrc, rhosrc
        logical dolock
        integer lmaxs, lmaxr, lmaxref
        character mname*(*)

        integer NL
        parameter (NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        common/modlly/mmmax
        integer mmmax
        common/depref/refdep
        real refdep
        common/shwave/bsh(NL), rhosh(NL), qbsh(NL)
        real bsh, rhosh, qbsh

        integer mmax

        integer iunit, iiso, iflsph, idimen, icnvel, ierr
        character title*80

        real v(NL), h(NL), qi(NL)

        real*8 c, s, t, x, p, tint, dxdp, vel, pnew, pupper
        real*8 ts
        real*8 sumx, sumt
        logical ext

        real tds, tdr
        common/earth/radius
        real radius

c-----
c       first read in the model and determine the medium parameters at the
c       source and receiver depths
c-----
c       get the earth model
c-----
        inquire(file=mname,exist=ext)
        if(.not. ext) call usage('Model file does not exist')
        l = lgstr(mname)

                call getmod(1,mname,mmax,title,iunit,iiso,iflsph,
     1          idimen,icnvel,ierr,.false.)
        mmmax = mmax
        if(ierr .lt. 0)return
                call adomod()
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
c-----
c       get the medium parameters at the source and reciever depths
c-----
        varec = a(lmaxr)
        vbrec = b(lmaxr)
        rhorec = rho(lmaxr)
        vasrc = a(lmaxs)
        vbsrc = b(lmaxs)
        rhosrc = rho(lmaxs)
c-----
c       prepare for the computations
c-----
        depths = hs + refdep
        depthr = hr + refdep

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
            call adosph()
            tds = radius*alog(radius/(radius-hs))
            tdr = radius*alog(radius/(radius-hr))
        else
            tds = depths
            tdr = depthr
        endif
c-----
c       now fill in velocity array according to desired first arrival
c       for SH there can be no water layer
c       for SV can be a water layer
c       Also define the Q for the T* analysis. Note we define
c        eventually q = 1/Q based on whether the given Q > or < 1
c-----
        do 100 i=1,mmax
            if(ipsvsh.eq.1)then
                v(i) = a(i)
                qi(i) = qa(i)
            else if(ipsvsh.eq.2)then
                v(i) = b(i)
                qi(i) = qb(i)
                if(b(i).le.0.001)then
                    v(i) = a(i)
                    qi(i) = qa(i)
                endif
            else if(ipsvsh.eq.3)then
                v(i) = bsh(i)
                qi(i) = qbsh(i)
            else if(ipsvsh.eq.4)then
                v(i) = a(i)
                qi(i) = qa(i)
            else if(ipsvsh.eq.5)then
                v(i) = a(i)
                qi(i) = qa(i)
            endif
            if(qi(i) .gt. 1.0)then
                qi(i) = 1.0 / qi(i)
            endif
            h(i) = d(i)
 100    continue
c-----
c       For the computations we look at four cases
c       1) direct path between source and receiver 
c          a) source and receiver in the same layer
c          b) source and receiver in different layers
c       2) refracted arrivals       
c          a) path is downward from source and then up to
c             receiver
c          b) path is upward from the source and then down to
c             receiver
c          This recognized the possibility that velocity does
c          not increase uniformly with depth
c-----
                    
c-----
c       direct arrival source/receiver at same layer
c-----
        if(v(lmaxs).eq.0.0)return
        if(v(lmaxr).eq.0.0)return
        if(lmaxs .eq. lmaxr)then
            tfirst = sqrt(dist**2 + abs(tds - tdr)**2)/
     1          v(lmaxs)
            rayp = (dist/sqrt(dist**2 + abs(tds - tdr)**2))/
     1          v(lmaxs)
            tstar = tfirst*qi(lmaxs)
        else
c-----
c       direct arrival source/receiver in different layers
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
            ps = 1.0/v(lmaxs)
            pr = 1.0/v(lmaxr)
            if(ps.lt.pr)then
                pupper = ps
            else
                pupper = pr
            endif
            do 1000 l=lmn,lmx
                if(v(l).eq.0.0)return
                p = 1.0/v(l)
                if(p.lt.pupper)pupper = p
 1000       continue
            p = 0.5d+00  * pupper
            do 1111 iter=1,10
                x = 0.0d+00
                t = 0.0d+00
                ts = 0.0d+00
                tint = 0.0d+00
                dxdp = 0.0d+00
                do 1500 l=lmn,lmx - 1
                    vel = dble(v(l))
                    s = p*vel
                    c = dsqrt(1.0d+00 - s*s)
                    t = t + dble(h(l)) /(vel*c)
                    x = x + dble(h(l)) * s / c
                    dxdp  = dxdp + dble(h(l)) *
     1                  vel/(c*c*c)
                    tint = tint + dble(h(l)) * c / vel
                    ts = ts + qi(l) * dble(h(l))/(c*vel)
                   

 1500           continue
                pnew = p - (x-dble(dist))/dxdp
c-----
c       safety - we must have a real ray, with upper bound
c       of  min[ 1/v(src), 1/v(rec)]
c-----  
                if(pnew .gt. pupper)then
                    if(iter.lt.10)then
                        pnew = 0.999d+00 * pupper
                    else
c-----
c       this is propably working like a refraction, so stop iterations
c-----  
                        t = tint + p * (dist)
                        go to 1112
                    endif
                endif
                p = pnew
 1111       continue
 1112       continue
            tfirst = t
            rayp = p
            tstar = ts
        endif
c-----
c       now proceed through the possible refracted arrivals
c       considering first upward rays from the source
c-----  
        if(lmn.gt.1)then
        do 3020 m=1,lmn-1
c-----
c       m is the refracting layer
c
c       get velocity of refractor
c-----
            vel = v(m)
            if(v(m).eq.0.0)return
            p = 1.0/vel
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
        if(v(lmn).ge.vel)go to 3020
        if(v(lmx).ge.vel)go to 3020
c-----
c       single leg
c-----
        sumt = 0.0
        sumx = 0.0
        ts = 0.0
            do 3021 l=1,lmx-1,lmn
                if(v(l).gt.vel)go to 3020
                cs = sqrt(abs(1.0 - p*p*v(l)*v(l)))
                sumt = sumt + h(l)*cs/v(l)
                sumx = sumx + h(l)*p/cs
                ts = ts + qi(l)*h(l)/(cs * v(l))
 3021       continue
            do 3022 l=m+1,lmn-1
                if(v(l).gt.vel)go to 3020
                cs = sqrt(abs(1.0 - p*p*v(l)*v(l)))
                sumt = sumt + 2.0*h(l)*cs/v(l)
                sumx = sumx + 2.0*h(l)*p/cs
                ts = ts + 2.0*qi(l)*h(l)/(cs * v(l))
 3022       continue
            tint = sumt
            tt = tint + dist / vel
            ts = ts + qi(m)*(dist-sumx)/v(m)
            if(tt .lt. tfirst .and. dist.ge.sumx)then
                  tfirst = tt
                  rayp = p
                 tstar = ts
            endif
 3020       continue
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
c       get velocity of refractor
c-----
            vel = v(m)
            if(v(m).eq.0.0)return
            p = 1.0/vel
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
c       refraction velocity
c-----
        if(v(lmn).ge.vel)go to 2020
        if(v(lmx).ge.vel)go to 2020
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
                      if(a(l).gt.vel)go to 2020
                      cs = sqrt(abs(1.0 - p*p*a(l)*a(l)))
                      sumt = sumt + 2.*h(l)*cs/a(l)
                      sumx = sumx + 2.*h(l)*p*a(l)/cs
                      if(qa(l).gt.1.0)qa(l) = 1.0/qa(l)
                      ts = ts + 2.*qa(l)*h(l)/(cs * a(l))
                  enddo
            else if(ipsvsh.eq.5)then
c-----
c               sP
c-----
                  do  l=lmaxref,lmaxs - 1
                      if(a(l).gt.vel)go to 2020
                      if(b(l).gt.vel)go to 2020
                      csa = sqrt(abs(1.0 - p*p*a(l)*a(l)))
                      csb = sqrt(abs(1.0 - p*p*b(l)*b(l)))
                      sumt = sumt + h(l)*csa/a(l)
     1                        +h(l)*csb/b(l)
                      sumx = sumx + 2.*h(l)*p*a(l)/csa
                      if(qa(l).gt.1.0)qa(l) = 1.0/qa(l)
                      if(qb(l).gt.1.0)qb(l) = 1.0/qb(l)
                      ts = ts + qa(l)*h(l)/(csa * a(l))
     1                        + qb(l)*h(l)/(csb * b(l))
                  enddo
            endif
c-----
c       continue
c-----
            do 2021 l=lmn,lmx - 1
                if(v(l).gt.vel)go to 2020
                if(v(l).eq.0.0)return
                cs = sqrt(abs(1.0 - p*p*v(l)*v(l)))
                sumt = sumt + h(l)*cs/v(l)
                sumx = sumx + h(l)*p*v(l)/cs
                ts = ts + qi(l)*h(l)/(cs * v(l))
 2021       continue
c-----
c       double leg
c-----
            do 2022 l=lmx,m-1
                if(v(l).gt.vel)go to 2020
                if(v(l).eq.0.0)return
                cs = sqrt(abs(1.0 - p*p*v(l)*v(l)))
                sumt = sumt + 2.0*h(l)*cs/v(l)
                sumx = sumx + 2.0*h(l)*p*v(l)/cs
                ts = ts + 2.*qi(l)*h(l)/(cs * v(l))
 2022       continue
            tint = sumt
            tt = tint + dist / vel
            ts = ts + qi(m)*(dist-sumx)/vel
            if(tt .lt. tfirst .and. dist.ge.sumx)then
                 tfirst = tt
                 rayp = p
                 tstar = ts
            endif
 2020       continue
             if(tfirst .eq. 1.0e+30)then
                tfirst = -12345.
                rayp   = -12345.
                tstar  = -12345.
             endif
        return
        end

        subroutine adosph()
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
c-----
c       mmax    I*4 number of layers
c       ipsvsh  I*4     1 - get P time
c                       2 - get SV time
c                       3 - get SH time
c-----
        integer NL
        parameter (NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        common/modlly/mmax
        integer mmax
        common/depref/refdep
        real refdep
        common/shwave/bsh(NL), rhosh(NL), qbsh(NL)
        real bsh, rhosh, qbsh

        double precision z0,z1,r0,r1,dr,ar,tmp

        common/earth/radius
        real radius

        ar=radius
        dr=0.0d0
        r0=ar + refdep
        d(mmax)=1.0
        do 10 i=1,mmax
            r1=r0-dble(d(i))
            z0=ar*dlog(ar/r0)
            z1=ar*dlog(ar/r1)
            d(i)=z1-z0
c-----
c        attempt 7 15 2007 - use standard rule but at mid layer depth as per DGH
c-----
            TMP=(ar+ar)/(r0+r1)

            a(i)=a(i)*tmp
            btmp = b(i)
            b(i)=btmp*tmp
            bsh(i)=btmp*tmp
            qbtmp = qb(i)
            qbsh(i)=qbtmp
            rhotmp=rho(i)
            rhosh(i) = rhotmp * tmp **(-5.0)
            rho(i) = rhotmp * tmp **(-2.275)
            r0 = r1
   10   continue
        d(mmax)=0.0
        return
        end

        subroutine adomod()
c-----
c       just fill the rhosh, bsh and qbsh arrays 
c-----
        integer NL
        parameter (NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        common/modlly/mmax
        integer mmax
        common/depref/refdep
        real refdep
        common/shwave/bsh(NL), rhosh(NL), qbsh(NL)
        real bsh, rhosh, qbsh

        do  i=1,mmax
            bsh(i)=b(i)
            qbsh(i)=qb(i)
            rhosh(i) = rho(i) 
        enddo
        return
        end

        subroutine srclyr(depth,lmax,dph)
        implicit none
        real depth, dph
        integer lmax
        integer LER, LIN, LOT
        parameter (LER=0, LIN=5, LOT=6)
        integer NL
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL),
     2      frefp(NL), frefs(NL)
        real d, a, b, rho, qa, qb, etap, etas, frefp, frefs
        common/modlly/mmax
        integer mmax
        common/shwave/bsh(NL), rhosh(NL), qbsh(NL)
        real bsh, rhosh, qbsh

        integer m
        real dep
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

        subroutine excitpw(freq,dop,rayp,Z,R,T,pvelrec,svelrec)
        real freq, omega, rayp,pvelrec,svelrec
        logical dop
        complex Z, R, T

        integer NL
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
            real d, a, b, rho, qa, qb, etap, etas, frefp, frefs
        common/modlly/mmax
            integer mmax
c-----
c       matrix components in layers and boundaries saved
c-----
        common/damp/alpha,ieqex
            real alpha
            integer ieqex
c-----
c       matrix components in layers and boundaries saved
c-----
        common/hamat/har
            complex*16 har(NL,4,4)
        common/damat/dar
            complex*16 dar(NL,5,5)
        common/hsrfr/hsr
            complex*16 hsr(2,5)
        common/gbrfr/gbr
            complex*16 gbr(2,5)
        common/hlmat/hal
            complex*16 hal(NL,2,2)
        common/hsrfl/hsl
            complex*16 hsl(2,2)
        common/gbrfl/gbl
            complex*16 gbl(2,2)
        common/hexex/hex
            real*8 hex(NL)
        common/hexexw/hexw
            real*8 hexw(NL)
        common/dexex/dex
            real*8 dex(NL)
        common/lexex/lex 
            real*8 lex(NL)
        common/water/iwater, iwats, iwatb
            integer iwater(NL), iwats(2), iwatb(2)
        common/updnsm/equalu, equald
            logical equalu(NL), equald(NL)
         logical compute

        complex*16 RA, RB
c-----
c       internal variables
c-----
        complex*16 om, om2, wvno, wvno2
c-----
c       P-SV quantities
c-----
        complex*16 aa(4,4), da(4,4), g(4,4)
        real*8 ex, exe, exa, exl
        complex*16   cd(5),e(5),fr, ca(5,5)
        complex*16 hl(2,2)
c-----
c       SH quantities
c-----
        complex *16 e1, e2, e21, e22,fl,  d11, d21
        real*8 exel, exb
        
        integer m, i, j
        integer jbdrys, jbdryh


        omega = 6.2831853*freq
        om =  dcmplx(dble(omega), dble(-alpha))
        om2 = om * om
        wvno  = dcmplx(dble(rayp),0.0d+00)*om
        wvno2 = wvno * wvno

        exe = 0.0d+00
        exa = 0.0d+00
        exl = 0.0d+00
        ex  = 0.0d+00
        exel= 0.0d+00
c-----
c       First evaluate all the Haskell matrices
c       Surface is free  , jbdrys = 1
c       Base is halfspace, jbdryh = 0
c-----
        jbdrys = 1
        jbdryh = 0
c-----
c       Check the model to avoid excessive computation
c-----
        call equlck()
c-----
c       evaluate the matrices
c-----  
        call evlmat(om,wvno,jbdrys,jbdryh,wvno2,om2)
        call evgmat(om, om2, wvno, wvno2, g, RA,RB)
c-----
c       initialize the halfspace
c-----
        do  i=1,5
            e(i) = gbr(1,i)
            cd(i) = e(i)
        enddo
        e1 = dcmplx(1.0d+00,0.0d+00)
        e2 = dcmplx(0.0d+00,0.0d+00)
c-----
c       now compute the Haskell Product from top down
c-----
        do 1000 j=1,4
            do 1100 i=1,4
                aa(i,j) = dcmplx(0.0d+00, 0.0d+00)
 1100       continue
            aa(j,j) = dcmplx(1.0d+00,0.0d+00)
 1000   continue
c-----
c       For the SH problem the equations are
c
c-----
                
        do 1340 m=1,mmax-1
c-----
c           get da(m) matrix
c-----
            call copy4(da,har,m,1,hex,ex)
            exl = exl + ex
            call copy2(hl,hal,m,1,lex,exb)
c-----
c           form da * aa
c-----
            call dmult(da,aa)
c-----
c           copy da to aa
c-----
            call copy(da,aa)
            d11=hl(1,1)*e1 + hl(1,2)*e2
            d21=hl(2,1)*e1 + hl(2,2)*e2
            e1 = d11
            e2 = d21
            exel = exel + exb

 1340   continue
        do 1341 m=mmax-1,1,-1
c-----
c       handle the compound matrices
c-----
            call copy5(ca,dar,m,1,dex,exa)
            call cmult(e,ca,exa,exe)
 1341   continue
c-----
c       put in the halfspace condition
c-----
        e(1) = e(1)*hsr(1,1) + e(2)*hsr(1,2) + e(3)*hsr(1,3)
     1      + e(4)*hsr(1,4) + e(5)*hsr(1,5)
            call dmult(g ,aa)
c-----
c       compute the receiver function
c-----
        if(dop)then
c-----
c           UZ
c-----
C            if(exl.lt.88)then
C            Z = -  g(2,1)*exp(-exl)/(g(1,1)*g(2,2) - g(2,1)*g(1,2))
            if((exe-exl).lt.88)then
            Z = -  g(2,1)*exp(exl-exe)/(e(1))
            else
                 Z = dcmplx(0.0d+00,0.0d+00)
            endif
c-----
c           add this to get displacement pulse
c-----
            Z = Z * pvelrec /(cmplx(0.0,1.0)*om)
c-----
c           Ur
c-----
C            if(exl.lt.88)then
C            R =   dcmplx(0.0d+00, 1.0d+00) * g(2,2)*exp(-exl)/
C     1          (g(1,1)*g(2,2) - g(2,1)*g(1,2))
            if((exe-exl).lt.88)then
            R =   dcmplx(0.0d+00, 1.0d+00) * g(2,2)*exp(exl-exe)/
     1          (e(1))
            else
                 R = dcmplx(0.0d+00,0.0d+00)
            endif
c-----
c           add this to get displacement pulse
c-----
            R = R * pvelrec /(cmplx(0.0,1.0)*om)
c-----
c           Ut
c-----
            T = cmplx(0.0,0.0)
        else
c-----
c           UZ
c-----
CRBH      IF(FREQ.GE.6.0 .and. FREQ.LE.7.0 .OR. FREQ.LT.01)then
CRBH          WRITE(6,*)'freq:',freq
CRBH          WRITE(6,*)'exl,g1212:',exl,g(1,1)*g(2,2) - g(2,1)*g(1,2)
CRBH          WRITE(6,*)'exe,e(1):',exe,e(1)*exp(exe-exl)
CRBH      ENDIF
C            if(exl.lt.88)then
C            Z =   dcmplx(0.0d+00, 1.0d+00) * g(1,1)*exp(-exl)/
C     1          (g(1,1)*g(2,2) - g(2,1)*g(1,2))
            if((exe-exl).lt.88)then
            Z =   dcmplx(0.0d+00, 1.0d+00) * g(1,1)*exp(exl-exe)/
     1          (e(1))
            else
                 Z = dcmplx(0.0d+00,0.0d+00)
            endif
c-----
c           add this to get displacement pulse
c-----
            Z = Z * svelrec/(cmplx(0.0,1.0)*om)
c-----
c           Ur
c-----
C            if(exl.lt.88)then
C            R = -                           g(1,2)*exp(-exl)/
C     1          (g(1,1)*g(2,2) - g(2,1)*g(1,2))
            if((exe-exl).lt.88)then
            R = -                           g(1,2)*exp(exl-exe)/
     1          (e(1))
            else
                 R = dcmplx(0.0d+00,0.0d+00)
            endif
c-----
c           add this to get displacement pulse
c-----
            R = - R * svelrec/(cmplx(0.0,1.0)*om)
c-----
c           Ut 
c FILL THIS IN - now the b must be causal
c-----
            if(exel.lt.88)then
                 T = 2.*exp(-exel)/(e1 +
     1         e2/(rho(mmax)*b(mmax)*b(mmax) *RB))
            else
                 T = dcmplx(0.0d+00,0.0d+00)
            endif
        endif
        return
        end

        subroutine evgmat(om, om2, wvno, wvno2, g, RA,RB)
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        common/modlly/mmax
        common/damp/alpha,ieqex
        complex*16 xka,xkb,ra,rb,gam,gamm1
        complex*16 atna,atnb
        real *8 ex,exa,exb
        complex*16 om,wvno,wvno2,om2
        complex*16 g(4,4)
        complex*16 zone
c-----
c       matrix components in layers and boundaries saved
c-----
        complex*16 CDSQRT
        zone  = dcmplx(1.0d+00,0.0d+00)
        m = mmax
            call aten(om,qa(m),qb(m),xka,xkb,
     1          alpha,a(m),b(m),atna,atnb,iwat,
     2          frefp(m),frefs(m))
            h =(dble(rho(m)*b(m)*b(m))*atnb*atnb)
            gam=dble(b(m))*(wvno/om)
            gam = gam * atnb
            gam = dcmplx(2.0d+00,0.0d+00)*gam*gam
            gamm1 = gam - zone
            ra=CDSQRT(wvno2-xka*xka)
            rb=CDSQRT(wvno2-xkb*xkb)
            g(1,1) =   gam/wvno
            g(2,1) = - gamm1/rb
            g(3,1) =   g(1,1)
            g(4,1) = - g(2,1)
            g(1,2) = - gamm1 / ra
            g(2,2) =   g(1,1)
            g(3,2) =   gamm1 / ra
            g(4,2) =   g(1,2)
            g(1,3) = - 1.0d+00/(rho(m)*om2)
            g(2,3) =   wvno / ( rb * rho(m)*om2)
            g(3,3) =   g(1,3)
            g(4,3) = - g(2,3)
            g(1,4) =   wvno / ( ra * rho(m)*om2)
            g(2,4) =   g(1,3)
            g(3,4) = - wvno / ( ra * rho(m)*om2)
            g(4,4) =   g(2,4)
            do i=1,4
                do j=1,4
                    g(i,j) = g(i,j)*0.5
                enddo
            enddo
        return
        end

        subroutine copy(da,aa)
c-----
c       copy da array to aa
c-----
        implicit none
        complex*16 da(4,4), aa(4,4)
        integer i, j
        
        do 1000 j=1,4
            do 1100 i=1,4
                aa(i,j) = da(i,j)
 1100       continue
 1000   continue
        return
        end

        subroutine excitsrc(freq,rayp,gg,hs,hr,lmaxs,lmaxr)

c-----
c     sample response for all wavenumbers at a given frequency
c     using Bouchon equal wavenumber sampling = dk
c     with offset of 0.218dk
c-----
        parameter(LER=0,LIN=5,LOT=6)
        parameter(NL=200)
        parameter(NSOURCE=100,NRECEIVER=100,NSR=100)
        integer lmaxs, lmaxr, mdpths, mdpthr
        common/bound/vmin,vamin,vamax,vbmin,vbmax
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        common/modlly/mmax
        common/jout/jsrc(21) , jbdrys, jbdryh
        common/damp/alpha,ieqex
c-----
c       control of wavefield
c-----
        common/ctlrud/rpup,rpdn,rsup,rsdn,dorud
        logical rpup, rpdn, rsup, rsdn, dorud
        common/ctlsud/spup,spdn,ssup,ssdn,dosud
        logical spup, spdn, ssup, ssdn, dosud

        common/cntrl/ishank,hnkarg,dstcor,dokjar
        integer*4 dstcor
        real*4 hnkarg
        logical ishank,dokjar

c-----
c       set up common blocks for wavenumber sampling at
c       suitable depths. This is necessary since the matrix
c       evaluation is done here for all source-receiver pairs
c       The source-receiver distance is important for the
c       wavenumber sampling at low frequencies
c-----
        common/kint1/gasymp
            logical gasymp(NSR)
        common/kint2/mkup
            integer mkup(NSR)
        common/kint3/wave
            real*4 wave(NSR,2)
        common/kint4/aa,bb,cc
            complex aa(NSR,21),bb(NSR,21),cc(NSR,21)

        complex*16 wvn,om, wvn2, om2
        complex gg(21)
        complex zeye
c-----
c       matrix components in layers and boundaries saved
c-----
        complex*16 har(NL,4,4), dar(NL,5,5), hsr(2,5), gbr(2,5), 
     1      hal(NL,2,2), hsl(2,2), gbl(2,2)
        real*8 hex(NL), lex(NL), dex(NL), hexw(NL)
        common/hamat/har
        common/damat/dar
        common/hsrfr/hsr
        common/gbrfr/gbr
        common/hlmat/hal
        common/hsrfl/hsl
        common/gbrfl/gbl
        common/hexex/hex
        common/hexexw/hexw
        common/dexex/dex
        common/lexex/lex 
        common/water/iwater(NL),iwats(2),iwatb(2)
        common/lyrctl/lyrins
        logical lyrins
        common/c/pmin,pmax,dp,pcntrl
        real pmin,pmax,dp,pcntrl


        zeye = cmplx(0.0,1.0)
        omega=6.2831853*freq
        om=dcmplx(dble(omega),-dble(alpha))
        om2 = om * om
c-----
c       output by ray parameter
c-----
            p  = rayp
c-----
c       safety
c-----
            if(p.eq.0.0)p = 0.01*dp
            wvn=dcmplx(dble(p),0.0d+00)*om
            wvn2 = wvn*wvn
c-----
c       evaluate matrices first
c-----
            call evlmat(om,wvn,jbdrys,jbdryh,wvn2,om2)
c-----
c       now evaluate for a specific source, receiver position
c-----
            call rshof(gg,om,wvn,lmaxs,
     1                  lmaxr,wvn2,om2)
c-----
c       To make radial look pulse like for small ray parameters
c       Also take time derivative of the point force and 
c           pressure fields
c       The multiplication by 'i' accomplishes this for 
c           responses with a linear term in wavenumber
c       Note that these are the integrands and not the 
c           output Green s functions. Thus
c       the TDS contribution (output=5) is integrand (13)
c-----
                    if(pcntrl .le.0.0)then
                        GG(2)  = - GG(2) * zeye
                        GG(3)  =   GG(3) * zeye
                        GG(5)  = - GG(5)
                        GG(6)  =   GG(6) * zeye
                        GG(14) =   GG(14) * zeye
                        GG(8)  = - GG(8) * zeye
                        GG(9)  =   GG(9) * zeye * om
                        GG(10) = - GG(10) * zeye * om * zeye
                        GG(11) =   GG(11) * zeye * om * zeye
                        GG(12) =   GG(12) * zeye * om 
                        GG(15) =   GG(15) * zeye * om 
                    endif
c-----
c       output
c-----
        return
        end

        subroutine aten(omega,qa,qb,xka,xkb,alpha,a,b,atna,atnb,iwat,
     2      frefp,frefs)
c-----
c       make velocities complex, using Futterman causality operator
c-----
        real*4 qa,qb,alpha,a,b
        complex*16 omega,at,atna,atnb,xka,xkb
        real*8 pi, om1p, om1s, oml, fac, pi2
        common/cntrl/ishank,hnkarg,dstcor,dokjar
        integer*4 dstcor
        real*4 hnkarg
        logical ishank,dokjar
        real*8 CDABS
        complex*16 CDLOG
c-----
c       reference frequency is fref hz
c-----
        om1p=6.2831853*frefp
        om1s=6.2831853*frefs
        pi2 = 1.5707963
        pi=3.1415927d+00
        if(dokjar)then
c-----
c       Kjartansson Constant Q, causal Q operator
c       Kjartansson, E. (1979). 
c          Constant Q-wave propagation and attenuation,
c       J. Geophys. Res. 84, 4737-4748.
c-----
            gama = atan(qa)/pi
            gamb = atan(qb)/pi
            if(gama.le.0.0)then
                atna = cmplx(1.0,0.0)
            else
                fac = pi2*gama
                rfac = sin(fac)/cos(fac)
                atna = dcmplx(1.0d+00,0.0d+00)/
     1              (( (omega/om1p)**dble(-gama) ) *
     2               dcmplx(1.0d+00,-dble(rfac)))
            endif
            if(b.gt.1.0e-04*a)then
                if(gamb.le.0.0)then
                    atnb = cmplx(1.0,0.0)
                else
                    fac = pi2*gamb
                    rfac = sin(fac)/cos(fac)
                    atnb = dcmplx(1.0d+00,0.0d+00)/
     1              (( (omega/om1s)**dble(-gamb) ) *
     2               dcmplx(1.0d+00,-dble(rfac)))
                endif
            endif
        else
c-----
c       Futterman Causal Q
c-----
c           low frequency cutoff is 0.01 hz
c-----
            oml=0.062831853d+00
            atna=dcmplx(1.0d+00,0.0d+00)
            atnb=dcmplx(1.0d+00,0.0d+00)
            if(qa.gt.0.0)then
              at=dcmplx(0.0d+00,0.0d+00)
              if(CDABS(omega).gt.oml) at=CDLOG(omega/om1p)/pi
              if(CDABS(omega).le.oml) then
                fac=dsqrt(oml*oml + dble(alpha*alpha))/oml
              at=CDLOG(dcmplx(dble(oml),-dble(alpha))/(om1p*fac))/pi
              endif
              atna=(1.+dble(qa)*at+dcmplx(0.0d+00,dble(qa/2.)))
            endif
            if(qb.gt.0.0)then
              at=dcmplx(0.0d+00,0.0d+00)
              if(CDABS(omega).gt.oml) at=CDLOG(omega/om1s)/pi
              if(CDABS(omega).le.oml) then
                fac=dsqrt(oml*oml + dble(alpha*alpha))/oml
              at=CDLOG(dcmplx(dble(oml),-dble(alpha))/(om1s*fac))/pi
              endif
               atnb=(1.+dble(qb)*at+dcmplx(0.0d+00,dble(qb/2.)))
            endif
        endif
        xka=omega/(dble(a)*atna)
        if(b.le.1.0e-04*a)then
            iwat = 1
            xkb = dcmplx(0.0d+00,0.0d+00)
        else
            iwat = 0
            xkb=omega/(dble(b)*atnb)
        endif
        return
        end

        subroutine cmult(e,ca,exa,exe)
        common/lwater/lfluid
        logical lfluid
c-----
c       FORM EC where e(1x5) c(5x5)
        complex*16 ca(5,5)
        real*8 exa,exe,eval
        real *8 xnorm
        complex*16 e(5)
        complex*16 c, ee(5)

        if(lfluid)then
            IUP = 2
        else
            IUP = 5
        endif
        do 1350 i=1,IUP
            c = dcmplx(0.0d+00,0.0d+00)
            do 1349 j=1,IUP
                c=c+ e(j) * ca(j,i)
 1349       continue
            ee(i)=c
 1350   continue
        exe = exe + exa
        call normc(ee,eval,xnorm)
        do 1351 i=1,IUP
            e(i) = ee(i)*xnorm
 1351   continue
        exe = exe + eval
        return
        end

        subroutine rcmult(y,c,exa,exe)
        common/lwater/lfluid
        logical lfluid
c-----
c       FORM YC where y(5x5) c(5x5) RETURN Y
c-----
        complex*16 c(5,5)
        real*8 exa,exe,eval
        real *8 xnorm
        complex*16 y(5,5)
        complex*16 ztmp, ee(5,5)
        if(lfluid)then
            IUP = 2
        else
            IUP = 5
        endif
        do 1350 i=1,IUP
            do 1351 j=1,IUP
                ztmp = dcmplx(0.0d+00,0.0d+00)
                do 1349 k=1,IUP
                    ztmp=ztmp+ y(i,k) * c(k,j)
 1349           continue
                ee(i,j)=ztmp
 1351       continue
 1350   continue
        exe = exe + exa
        call rnormc(ee,eval,xnorm)
        do 1353 j=1,IUP
            do 1352 i=1,IUP
                y(i,j) = ee(i,j)*xnorm
 1352       continue
 1353   continue
        exe = exe + eval
        return
        end

        subroutine dmult(da,aa)
c-----
c       propagator up
c       FORM D = DA
c-----
        complex*16 aa(4,4)
        complex*16 sumd,ea(4,4),da(4,4)
        do 1360 i=1,4
            do 1361 j=1,4
                sumd = dcmplx(0.0d+00,0.0d+00)
                do 1362 jj=1,4
                    sumd=sumd+da(i,jj) * aa(jj,j)
 1362           continue
                ea(i,j)=sumd
 1361       continue
 1360   continue
        do 1363 j=1,4
            do 1364 i=1,4
                da(i,j)=ea(i,j)
 1364       continue
 1363   continue
        return
        end

        subroutine dnka(ca,wvno,wvno2,om2,gam,rho,iwat,w,x,cosp,ex)
        complex*16 ca(5,5)
        complex*16 wvno, wvno2, om2
        real*4 rho
        complex*16 gam
        complex*16 w,x,cosp
        real*8 ex
        common/ ovrflw / a0,cpcq,cpy,cpz,cqw,cqx,xy,xz,wy,wz
        real *8 a0
        complex*16 cpcq,cpy,cpz,cqw,cqx,xy,xz,wy,wz

        complex*16 gam2,gamm1,gamm2,a0c,xz2,wy2,temp
        real*8 dfac
        complex*16 cqww2, cqxw2, g1wy2, gxz2, g2wy2, g2xz2
        complex*16 gg1, a0cgg1
        complex*16 zrho, zrho2
c-----
c        A11     A12     A13    -A13     A15     A16
c        A21     A22     A23    -A23     A25     A15
c        A31     A32     A33    1-A33   -A23    -A13
c       -A31    -A32    1-A33    A33     A23     A13
c        A51     A52    -A32     A32     A22     A12
c        A61     A51    -A31     A31     A21     A11
c-----
c       this will be multipled on the left by the G matrix
c
c       [ G11   G12 G13 -G13    G15 G16 ]
c
c-----
c       or on the right by
c
c       [ H11   H21 H31 -H31    H51 H61  ] ^T
c-----
c       the number of multiplications can be reduced from 36 to 25 
c       if we define a new matrices
c       related to the original matrices by
c-----
c         A11     A12     A13         A15     A16
c         A21     A22     A23         A25     A15
c        2 A31   2 A32   2 A33 -1   -2 A23  -2 A13
c         A51     A52    -A32         A22     A12
c         A61     A51    -A31         A21     A11
c-----
c
c       [ G11   G12  G13    G15 G16  ]
c       [ H11   H21 2 H31   H51 H61  ] ^T
c
c-----
c       this means that some of the original definitions of the 
c       Aij elements must be changed for the
c       definition of the modified 5x5 compount A matrix
c
c       old 6x6                 new 5x5
c       A11 = 1 - A33 + A22     1 - (1/2)(new A33 + 1) + new A2)
c       A53 = -A32              A43 = - (1/2) new A32
c       A63 = -A31              A53 = - (1/2) new A31
c-----
c       To recover the needed elements, we note that the 
c       old G14 = -old G14 = new G13
c-----
c-----
        zrho = dcmplx(dble(rho),0.0d+00)
        if(iwat.eq.1)then
c-----
c       fluid layer
c-----
            do 100 j=1,5
                do 101 i=1,5
                    ca(i,j) = dcmplx(0.0d+00, 0.0d+00)
  101           continue
  100       continue
            if(ex.gt.35.0d+00)then
                dfac = 0.0d+00
            else
                dfac = dexp(-ex)
            endif
            ca(3,3) = dfac
            ca(1,1) = cosp
            ca(5,5) = cosp
            ca(1,2) = - x/(zrho*om2)
            ca(2,1) = - zrho*w*om2
            ca(2,2) = cosp
            ca(4,4) = cosp
            ca(4,5) = ca(1,2)
            ca(5,4) = ca(2,1)
        else
c-----
c       elastic layer
c-----
            zrho2= dcmplx(dble(rho*rho),0.0d+00)
            gam2  = gam*gam
            gamm1 = gam-1.
            gamm2 = gamm1*gamm1
            cqww2 = cqw * wvno2
            cqxw2 = cqx / wvno2
            gg1 = gam*gamm1
            a0c  = dcmplx(2.0d+00,0.0d+00)*
     1          (dcmplx(a0,0.0d+00)-cpcq)
            xz2  = xz/wvno2
            gxz2 = gam*xz2
            g2xz2 = gam2 * xz2
            a0cgg1 = a0c*(gam+gamm1)
            wy2  = wy*wvno2
            g2wy2 = gamm2 * wy2
            g1wy2 = gamm1 * wy2

c-----
c       OK by symmetry
c----
            temp = a0c*gg1 + g2xz2 + g2wy2
            ca(3,3) = a0 + temp + temp
            ca(1,1) = cpcq-temp
            ca(1,2) = (-cqx + wvno2*cpy)/(zrho*om2)
            temp = dcmplx(0.5d+00,0.0d+00)*a0cgg1 + gxz2 + g1wy2
            ca(1,3) = wvno*temp/(zrho*om2)

            ca(1,4) = (-cqww2+cpz)/(zrho*om2)
            temp = wvno2*(a0c + wy2) + xz
            ca(1,5) = -temp/(zrho2*om2*om2)

            ca(2,1) = (-gamm2*cqw + gam2*cpz/wvno2)*zrho*om2
            ca(2,2) = cpcq
            ca(2,3) = (gamm1*cqww2 - gam*cpz)/wvno
            ca(2,4) = -wz
            ca(2,5)=ca(1,4)


            temp =dcmplx(0.5d+00,0.0d+00)*a0cgg1*gg1 
     1          + gam2*gxz2 + gamm2*g1wy2
            ca(3,1) = -dcmplx(2.0d+00,0.0d+00)*temp*zrho*om2/wvno
            ca(3,2) = -wvno*(gam*cqxw2 - gamm1*cpy)*
     1          dcmplx(2.0d+00,0.0d+00)

            ca(3,4)=-dcmplx(2.0d+00,00d+00)*ca(2,3)
            ca(3,5)=-dcmplx(2.0d+00,00d+00)*ca(1,3)

            ca(4,1) = (-gam2*cqxw2 + gamm2*cpy)*zrho*om2
            ca(4,2) = -xy
            ca(4,3)= -ca(3,2)/(dcmplx(2.0d+00,00d+00))
            ca(4,4)=ca(2,2)
            ca(4,5)=ca(1,2)

            temp = gamm2*(a0c*gam2 + g2wy2) + gam2*g2xz2
            ca(5,1) = -zrho2*om2*om2*temp/wvno2
            ca(5,2)=ca(4,1)
            ca(5,3)=-ca(3,1)/(dcmplx(2.0d+00,00d+00))
            ca(5,4)=ca(2,1)
            ca(5,5)=ca(1,1)
        endif
        return
        end

        subroutine hska(aa,w,x,y,z,cosp,cosq,wvno,wvno2,gam,gamm1,rho,
     1      iwat,ex,om2)
        complex*16 wvno,wvno2
        complex*16 aa(4,4),w,x,y,z,cosp,cosq,gam,gamm1
        complex*16 cpq, gcpq, zw2, gzw2, g1w, g1y, gx
        complex*16 zrho
        real*8 ex, dfac
        complex*16 om2
        zrho = dcmplx(dble(rho),0.0d+00)
        if(iwat.eq.1)then
c-----
c       fluid layer
c-----
            do 100 j=1,4
                do 101 i=1,4
                    aa(i,j) = dcmplx(0.0d+00,0.0d+00)
  101           continue
  100       continue
            if(ex.gt.35.0d+00)then
                dfac = 0.0d+00
            else
                dfac = dexp(-ex)
            endif
            aa(1,1) = dfac
            aa(4,4) = dfac
            aa(2,2) = cosp
            aa(3,3) = cosp
            aa(2,3) = -x/(zrho*om2)
            aa(3,2) = - zrho*w*om2
        else
c-----
c       elastic layer
c-----
            cpq = cosp-cosq
            gcpq = gam*cpq
            zw2 = z/wvno2
            gzw2 = gam*zw2
            g1w = gamm1*w
            g1y = gamm1*y
            gx = gam*x
            aa(1,1)=   gcpq + cosq
                aa(1,3)= - wvno * cpq/(zrho*om2)
            aa(1,2)=   wvno*(-g1w+gzw2)
            aa(1,4)=   (wvno2*w-z)/(zrho*om2)
            aa(2,1)=   (gx - wvno2*g1y)/wvno
            aa(2,2)= - gcpq + cosp
            aa(2,3)=   (-x+wvno2*y)/(zrho*om2)
            aa(2,4)= - aa(1,3)
            aa(3,1)=   zrho*om2*gamm1*gcpq/wvno
            aa(3,2)=   zrho*om2*((-gamm1*g1w)+(gam*gzw2))
            aa(3,3)=   aa(2,2)
            aa(3,4)= - aa(1,2)
            aa(4,1)=   zrho*om2*(((gam*gx)/wvno2) - (gamm1*g1y))
            aa(4,2)= - aa(3,1)
            aa(4,3)= - aa(2,1)
            aa(4,4)=   aa(1,1)
        endif
        return
        end

        subroutine hskl(hl,cosql,yl,zl,h,iwat)
        complex*16 hl(2,2)
        complex*16 cosql,yl,zl,h
        integer iwat
c-----
c       cosql = cosh ( nu d )
c       zl    = nu sinh ( nu d )
c       yl    = sinh ( nu d ) / nu
c       h     = rho beta ^ 2
c-----
        if(iwat.eq.0)then   
            hl(1,1) = cosql
            hl(2,1) = zl*h
            hl(1,2) = yl/h
            hl(2,2) = cosql
        else
            hl(1,1) = dcmplx(1.0d+00,0.0d+00)
            hl(1,2) = dcmplx(0.0d+00,0.0d+00)
            hl(2,1) = dcmplx(0.0d+00,0.0d+00)
            hl(2,2) = dcmplx(1.0d+00,0.0d+00)
        endif
        return
        end 

        subroutine lmult(d11,d12,d21,d22,hl,iwat,exel,exb,icomp)
c-----
c       multiply SH matrix by a row vector on left
c-----
        complex*16 d11,d12,d21,d22,hl(2,2),e1,e2
        real*8 exel, exb
        logical icomp
c-----
c       fluid layer do nothing, just return, 
c       equivalent to multiplying by
c       identity matrix
c-----
        if(iwat.eq.0)then
c-----
c       elastic layer
c-----
            e1=d11
            e2=d12
c-----
c           a11 = cosql
c           a12 = yl
c           a21 = zl
c           a22 = cosql
c-----
            d11=e1*hl(1,1) + e2*hl(2,1)
            d12=e1*hl(1,2) + e2*hl(2,2)
            exel = exel + exb
            if(icomp)then
                e1=d21
                e2=d22
                d21=e1*hl(1,1) + e2*hl(2,1)
                d22=e1*hl(1,2) + e2*hl(2,2)
            endif
        endif
        return
        end

        subroutine normc(e,ex,xnorm)
        common/lwater/lfluid
        logical lfluid
        real*8 ex
        real *8 test,testt,x,y,fac,xnorm
        complex*16 e(*)
        real*8 DREAL
        test = 0.0D+00
        testt = 0.0D+00
        if(lfluid)then
            IUP = 2
        else
            IUP = 5
        endif
        do 2 i = 1,IUP
            if(dabs(dreal(e(i))).gt.testt) testt =dabs(dreal(e(i)))
            if(dabs(dimag(e(i))).gt.testt) testt =dabs(dimag(e(i)))
    2   continue
        if(testt.lt.1.0e-30)testt=1.0
        do 1 i =1,IUP
            x=dreal(e(i))/testt
            y=dimag(e(i))/testt
            fac = x*x + y*y
            if(test.lt.fac) test = fac
    1   continue
        test = testt*dsqrt(test)
        if(test.lt.1.0d-30) test=1.0
        xnorm = 1./test
        ex =-dlog(xnorm)
        return
        end

        subroutine rnormc(e,ex,xnorm)
        common/lwater/lfluid
        logical lfluid
        real*8 ex
        real *8 test,testt,x,y,fac,xnorm
        complex*16 e(5,5)
        real*8 DREAL
        test = 0.0D+00
        testt = 0.0D+00
        if(lfluid)then
            IUP = 2
        else
            IUP = 5
        endif
        do 3 j=1,IUP
            do 2 i = 1,IUP
            if(dabs(dreal(e(i,j))).gt.testt) 
     1          testt =dabs(dreal(e(i,j)))
            if(dabs(dimag(e(i,j))).gt.testt) 
     1          testt =dabs(dimag(e(i,j)))
    2       continue
    3   continue
        if(testt.lt.1.0e-30)testt=1.0
        do 4 j=1,IUP
            do 1 i =1,IUP
                x=dreal(e(i,j))/testt
                y=dimag(e(i,j))/testt
                fac = x*x + y*y
                if(test.lt.fac) test = fac
    1       continue
    4   continue
        test = testt*dsqrt(test)
        if(test.lt.1.0d-30) test=1.0
        xnorm = 1./test
        ex =-dlog(xnorm)
        return
        end

        subroutine evlmat(om,wvno,jbdrys,jbdryh,wvno2,om2)
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        common/modlly/mmax
        parameter(NSOURCE=100,NRECEIVER=100,NSR=100)
        common/source/depths(NSOURCE),lmaxs(NSOURCE),mdpths
        common/receiv/depthr(NRECEIVER),lmaxr(NRECEIVER),mdpthr
        real*4 depths, depthr
        integer lmaxs, lmaxr, mdpths, mdpthr
        common/damp/alpha,ieqex
        complex*16 om,xka,xkb,ra,rb,gam,gamm1,p,q,w,x,y
        complex*16 cosq,z,cosp
        complex*16  ca(5,5), hl(2,2)
        complex*16 atna,atnb
        complex*16 yl,zl,cosql
        real *8 ex,exa,exb
        complex*16 wvno,wvno2,om2
        complex*16 h
        complex*16 aa(4,4)
        complex*16 zone
c-----
c       matrix components in layers and boundaries saved
c-----
        complex*16 har(NL,4,4), dar(NL,5,5), hsr(2,5), gbr(2,5), 
     1      hal(NL,2,2), hsl(2,2), gbl(2,2)
        real*8 hex(NL), lex(NL), dex(NL), hexw(NL)
        common/hamat/har
        common/damat/dar
        common/hsrfr/hsr
        common/gbrfr/gbr
        common/hlmat/hal
        common/hsrfl/hsl
        common/gbrfl/gbl
        common/hexex/hex
        common/hexexw/hexw
        common/dexex/dex
        common/lexex/lex 
        common/water/iwater(NL),iwats(2),iwatb(2)
        common/updnsm/equalu(NL), equald(NL)
        logical equalu, equald
        logical compute
        complex*16 CDSQRT
        zone  = dcmplx(1.0d+00,0.0d+00)
c-----
c       evaluate the G matrix 
c           gbr(1,x) is for normal stack
c           gbr(2,x) is for inverted stack
c-----
        call evalg(jbdryh,mmax,mmax-1,gbr,gbl,1,wvno,om,om2,wvno2)
        call evalg(jbdrys,1,   2,     gbr,gbl,2,wvno,om,om2,wvno2)
c-----
c       evaluate the H matrix
c           hsr(1,x) is for normal stack
c           hsr(2,x) is for inverted stack
c-----
        call evalh(jbdrys,1,   2,     hsr,hsl,1,wvno,om,om2,wvno2)
        call evalh(jbdryh,mmax,mmax-1,hsr,hsl,2,wvno,om,om2,wvno2)
c-----

c-----
c       matrix multiplication from bottom layer upward
c-----
        do 1340 m = 1,mmax,1
c-----
c       first check to see if computations already done
c-----
            if(equald(m))then
                compute = .false.
            else
                compute = .true.
            endif
            if(compute)then
                call aten(om,qa(m),qb(m),xka,xkb,
     1              alpha,a(m),b(m),atna,atnb,iwat,
     2              frefp(m),frefs(m))
                h =(dble(rho(m)*b(m)*b(m))*atnb*atnb)
                gam=dble(b(m))*(wvno/om)
                gam = gam * atnb
                gam = dcmplx(2.0d+00,0.0d+00)*gam*gam
                gamm1 = gam - zone
                ra=CDSQRT(wvno2-xka*xka)
                rb=CDSQRT(wvno2-xkb*xkb)
                p=ra*dble(d(m))
                q=rb*dble(d(m))
                call var(p,q,ra,rb,w,x,y,z,cosp,cosq,
     1              ex,exa,exb,yl,zl,cosql,iwat)
                call dnka(ca,wvno, wvno2, om2,gam,rho(m),
     1                  iwat,w,x,cosp,ex)
                call hska(aa,w,x,y,z,cosp,cosq,wvno,wvno2,
     1                  gam,gamm1,rho(m),
     1              iwat,ex,om2)
                call hskl(hl,cosql,yl,zl,h,iwat)
            endif
            iwater(m) = iwat
            call copy5(ca,dar,m,0,dex,exa)
            call copy4(aa,har,m,0,hex,ex)
            call copy2(hl,hal,m,0,lex,exb)
 1340   continue
        return
        end

        subroutine equlck()
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        common/modlly/mmax
        common/updnsm/equalu(NL), equald(NL)
        logical equalu, equald
c-----
c       To avoid repeated computation, check to see if 
c       neighboring layers
c       are identical, once for going up and another for going down
c-----
c       First check top down
c-----
        do 100 m=1,mmax
            if(m.eq.1)then
                equald(m) = .false.
            else if(m.gt.1
     1          .and. a(m).eq.a(m-1) 
     2          .and. b(m).eq.b(m-1)
     3          .and. d(m).eq.d(m-1) 
     4          .and. rho(m).eq.rho(m-1)
     5          .and. qa(m).eq.qa(m-1)
     6          .and. qb(m).eq.qb(m-1) )then
                equald(m) = .true.
            else
                equald(m) = .false.
            endif
  100   continue
c-----
c       check bottom up
c-----
        do 200 m=1,mmax
            if(m.eq.mmax)then
                equalu(m) = .false.
            else if(m.lt.mmax
     1          .and. a(m).eq.a(m+1) 
     2          .and. b(m).eq.b(m+1)
     3          .and. d(m).eq.d(m+1) 
     4          .and. rho(m).eq.rho(m+1)
     5          .and. qa(m).eq.qa(m+1)
     6          .and. qb(m).eq.qb(m+1) )then
                equalu(m) = .true.
            else
                equalu(m) = .false.
            endif
  200   continue
        return
        end

        subroutine srclay(depth,lmax,dph)
        parameter (LER=0, LIN=5, LOT=6)
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        common/modlly/mmax
        common/lyrctl/lyrins
        logical lyrins
        if(.not.lyrins)then
            call modcpy(.false.)
            call insert(depth)
        endif
        call srclyr(depth,lmax,dph)
        return
        end


        subroutine copy5(ca,dar,m,itofrm,dex,exa)
        parameter (NL=200)
        complex*16 dar(NL,5,5)
        complex*16 ca(5,5)
        integer itofrm
        real*8 dex(NL)
        real*8 exa
c-----
c       copy from ca to dar
c-----
        if(itofrm.eq.0)then
            do 100 j=1,5
                do 110 i=1,5
                    dar(m,i,j) = ca(i,j)
  110           continue
                dex(m) = exa
  100       continue
c-----
c       copy from dar to ca
c-----
        else
            do 200 j=1,5
                do 210 i=1,5
                    ca(i,j) = dar(m,i,j)
  210           continue
                exa = dex(m)
  200       continue
        endif
        return
        end

        subroutine copy4(aa,har,m,itofrm,hex,ex)
        parameter (NL=200)
        complex*16 har(NL,4,4)
        complex*16 aa(4,4)
        integer itofrm
        real*8 hex(NL)
        real*8 ex
c-----
c       copy from aa to har
c-----
        if(itofrm.eq.0)then
            do 100 j=1,4
                do 110 i=1,4
                    har(m,i,j) = aa(i,j)
  110           continue
                hex(m) = ex
  100       continue
c-----
c       copy from har to aa
c-----
        else
            do 200 j=1,4
                do 210 i=1,4
                    aa(i,j) = har(m,i,j)
  210           continue
                ex = hex(m)
  200       continue
        endif
        return
        end

        subroutine copy2(hl,hal,m,itofrm,lex,exb)
        parameter (NL=200)
        complex*16 hal(NL,2,2)
        complex*16 hl(2,2)
        integer itofrm
        real*8 lex(NL)
        real*8 exb
c-----
c       copy from hl to hal
c-----
        if(itofrm.eq.0)then
            do 100 j=1,2
                do 110 i=1,2
                    hal(m,i,j) = hl(i,j)
  110           continue
                lex(m) = exb
  100       continue
c-----
c       copy from hal to hl
c-----
        else
            do 200 j=1,2
                do 210 i=1,2
                    hl(i,j) = hal(m,i,j)
  210           continue
                exb = lex(m)
  200       continue
        endif
        return
        end

        subroutine evalg(jbdry,m,m1,gbr,gbl,in,
     1      wvno,om,om2,wvno2)
        complex*16 gbr(2,5), gbl(2,2)
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        common/modlly/mmax
        common/damp/alpha,ieqex
        complex*16 om,xka,xkb,ra,rb,gam,gamm1
        complex*16 atna,atnb
        complex*16 wvno,wvno2,om2
        complex*16 h
        common/lwater/lfluid
        logical lfluid
        complex*16 CDSQRT
c-----
c       set up halfspace conditions
c-----
        call aten(om,qa(m),qb(m),xka,xkb,alpha,a(m),
     1      b(m),atna,atnb,iwat,
     2      frefp(m),frefs(m))
        ra=CDSQRT(wvno2-xka*xka)
        rb=CDSQRT(wvno2-xkb*xkb)
        gam = dble(b(m))*wvno/om
        gam = gam * atnb
        gam = 2.0d+00 * (gam * gam)
        gamm1 = gam - dcmplx(1.0d+00,0.0d+00)
        h =(dble(rho(m)*b(m)*b(m))*atnb*atnb)
c-----
c       set up halfspace boundary conditions
c
c       jbdry   = -1  RIGID
c           =  0  ELASTIC
c           = +1  FREE SURFACE
c
c-----
        if(jbdry.lt.0)then
c-----
c       RIGID - check properties of layer above
c-----
            if(b(m) .gt. 0.0)then
c-----
c               ELASTIC ABOVE - RIGID
c-----
                gbr(in,1) = dcmplx(1.0d+00,0.0d+00)
                gbr(in,2) = dcmplx(0.0d+00,0.0d+00)
                gbr(in,3) = dcmplx(0.0d+00,0.0d+00)
                gbr(in,4) = dcmplx(0.0d+00,0.0d+00)
                gbr(in,5) = dcmplx(0.0d+00,0.0d+00)
                gbl(in,1) = dcmplx(1.0d+00,0.0d+00)
                gbl(in,2) = dcmplx(0.0d+00,0.0d+00)
            else
c-----
c               FLUID ABOVE - RIGID
c-----
                gbr(in,1) = dcmplx(0.0d+00,0.0d+00)
                gbr(in,2) = dcmplx(0.0d+00,0.0d+00)
                gbr(in,3) = dcmplx(0.0d+00,0.0d+00)
                gbr(in,4) = dcmplx(0.0d+00,0.0d+00)
                gbr(in,5) = dcmplx(0.0d+00,0.0d+00)
                if(lfluid)then
                    gbr(in,1) = dcmplx(1.0d+00,0.0d+00)
                else
                    gbr(in,4) = dcmplx(1.0d+00,0.0d+00)
                endif
c-----
c               (pseudo SH)
c-----
                gbl(in,1) = dcmplx(1.0d+00,0.0d+00)
                gbl(in,2) = dcmplx(0.0d+00,0.0d+00)
            endif
        else if(jbdry.eq.0)then
c-----
c       HALFSPACE
c-----
            if(iwat.eq.0)then
c-----
c               ELASTIC HALFSPACE
c-----
c       multiply G of Herrmann 2001 by - rho^2 om^4 k^2 ra rb
c       should have no effect since it is in both numerator and
c       denominator -- however will not give the correct potential
c       coefficient -- so rethink?
c-----



                gbr(in,1)=dble(rho(m)*rho(m))*om2*om2*
     1              (-gam*gam*ra*rb+wvno2*gamm1*gamm1)
                gbr(in,2)=-dble(rho(m))*(wvno2*ra)*om2
                gbr(in,3)=-dble(rho(m))*(-gam*ra*rb+wvno2*gamm1)
     1              *om2*wvno
                gbr(in,4)=dble(rho(m))*(wvno2*rb)*om2
                gbr(in,5)=wvno2*(wvno2-ra*rb)
        gbr(in,1) = 0.25*gbr(in,1)/(-rho(m)*rho(m)*om2*om2*wvno2*ra*rb)
        gbr(in,2) = 0.25*gbr(in,2)/(-rho(m)*rho(m)*om2*om2*wvno2*ra*rb)
        gbr(in,3) = 0.25*gbr(in,3)/(-rho(m)*rho(m)*om2*om2*wvno2*ra*rb)
        gbr(in,4) = 0.25*gbr(in,4)/(-rho(m)*rho(m)*om2*om2*wvno2*ra*rb)
        gbr(in,5) = 0.25*gbr(in,5)/(-rho(m)*rho(m)*om2*om2*wvno2*ra*rb)

C               gbl(in,1) = dble(rho(m))*rb
C               gbl(in,2) = dcmplx(1.0d+00,0.0d+00)
C     1             /(dble(b(m)*b(m))*atnb*atnb)
                gbl(in,1) =  dble(rho(m))*(dble(b(m)*b(m))*atnb*atnb)*rb
                gbl(in,2) =  dcmplx(1.0d+00,0.0d+00)
            else if(iwat.eq.1)then
c-----
c               FLUID HALFSPACE
c-----
                if(lfluid)then
                    gbr(in,1) = dble(0.5) / ra
                    gbr(in,2) = dcmplx(0.5d+00,0.0d+00)/
     1                  (-dble(rho(m))*om2)
                    gbr(in,3) = dcmplx(0.0d+00,0.0d+00)
                    gbr(in,4) = dcmplx(0.0d+00,0.0d+00)
                    gbr(in,5) = dcmplx(0.0d+00,0.0d+00)
                else
                    gbr(in,1) = dcmplx(0.0d+00,0.0d+00)
                    gbr(in,2) = dcmplx(0.0d+00,0.0d+00)
                    gbr(in,3) = dcmplx(0.0d+00,0.0d+00)
                    gbr(in,4) = dble(0.5*rho(m)*om2) / ra
                    gbr(in,5) = dcmplx(-0.5d+00,0.0d+00)
                endif
                gbl(in,1) = dcmplx(0.0d+00,0.0d+00)
                gbl(in,2) = dcmplx(1.0d+00,0.0d+00)
            endif
        else if(jbdry.eq.1)then
c-----
c       FREE - check properties of layer above
c-----
            if(b(m) .gt. 0.0)then
                gbr(in,1) = dcmplx(0.0d+00,0.0d+00)
                gbr(in,2) = dcmplx(0.0d+00,0.0d+00)
                gbr(in,3) = dcmplx(0.0d+00,0.0d+00)
                gbr(in,4) = dcmplx(0.0d+00,0.0d+00)
                gbr(in,5) = dcmplx(1.0d+00,0.0d+00)
                gbl(in,1) = dcmplx(0.0d+00,0.0d+00)
                gbl(in,2) = dcmplx(1.0d+00,0.0d+00)
                
            else
                gbr(in,1) = dcmplx(0.0d+00,0.0d+00)
                gbr(in,2) = dcmplx(0.0d+00,0.0d+00)
                gbr(in,3) = dcmplx(0.0d+00,0.0d+00)
                gbr(in,4) = dcmplx(0.0d+00,0.0d+00)
                gbr(in,5) = dcmplx(0.0d+00,0.0d+00)
                if(lfluid)then
                    gbr(in,2) = dcmplx(1.0d+00,0.0d+00)
                else
                    gbr(in,5) = dcmplx(1.0d+00,0.0d+00)
                endif
                gbl(in,1) = dcmplx(0.0d+00,0.0d+00)
                gbl(in,2) = dcmplx(1.0d+00,0.0d+00)
            endif
        endif
        return
        end

        subroutine evalh(jbdry,m,m1,hsr,hsl,in,wvno,om,om2,wvno2)
        complex*16 hsr(2,5), hsl(2,2)
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        common/modlly/mmax
        common/damp/alpha,ieqex
        complex*16 om,xka,xkb,ra,rb,gam,gamm1
        complex*16 atna,atnb
        complex*16 wvno,wvno2,om2
        complex*16 h
        complex*16 CDSQRT
c-----
c       set up surface conditions
c-----
        call aten(om,qa(m),qb(m),xka,xkb,alpha,a(m),
     1      b(m),atna,atnb,iwat,
     2      frefp(m),frefs(m))
        ra=CDSQRT(wvno2-xka*xka)
        rb=CDSQRT(wvno2-xkb*xkb)
        gam = dble(b(m))*wvno/om
        gam = gam * atnb
        gam = 2.0d+00 * (gam * gam)
        gamm1 = gam - dcmplx(1.0d+00,0.0d+00)
        h =(dble(rho(m)*b(m)*b(m))*atnb*atnb)
c-----
c       do top surface condition
c
c           = -1  RIGID
c           =  0  ELASTIC
c           = +1  FREE SURFACE
c-----
        if(iwat.eq.0)then
            if(jbdry.eq.0)then
c-----
c           ELASTIC ABOVE SOLID
c-----
                hsr(in,1) = (wvno2 - ra*rb)
                hsr(in,2) = dble(rho(m))*om2*rb
                hsr(in,3) = 2.0d+00*
     1              wvno*dble(rho(m))*om2*
     1              ( gamm1 - gam*ra*rb/wvno2)
                hsr(in,4) = -dble(rho(m))*om2*ra
                hsr(in,5) = dble(rho(m)*rho(m))*om2*om2
     1              *(gamm1*gamm1 - 
     1              gam*gam*ra*rb/wvno2)
                hsl(in,1) = dcmplx(0.5d+00,0.0d+00)
                hsl(in,2) = 0.5d+00*h*rb
            else if(jbdry.eq.-1)then
c-----
c           RIGID ABOVE SOLID
c-----
                hsr(in,1) = dcmplx(0.0d+00,0.0d+00)
                hsr(in,2) = dcmplx(0.0d+00,0.0d+00)
                hsr(in,3) = dcmplx(0.0d+00,0.0d+00)
                hsr(in,4) = dcmplx(0.0d+00,0.0d+00)
                hsr(in,5) = dcmplx(1.0d+00,0.0d+00)
                hsl(in,1) = dcmplx(0.0d+00,0.00d+00)
                hsl(in,2) = dcmplx(1.0d+00,0.00d+00)
            else if(jbdry.eq.1)then
c-----
c           FREE ABOVE SOLID
c-----
                hsr(in,1) = dcmplx(1.0d+00,0.0d+00)
                hsr(in,2) = dcmplx(0.0d+00,0.0d+00)
                hsr(in,3) = dcmplx(0.0d+00,0.0d+00)
                hsr(in,4) = dcmplx(0.0d+00,0.0d+00)
                hsr(in,5) = dcmplx(0.0d+00,0.0d+00)
                hsl(in,1) = dcmplx(1.0d+00,0.00d+00)
                hsl(in,2) = dcmplx(0.0d+00,0.00d+00)
            endif
        else if(iwat.gt.0)then
            if(jbdry.eq.0)then
c-----
c           HALFSPACE ABOVE FLUID
c-----
                hsr(in,1) = ra
                hsr(in,2) = -dble(rho(m))*om2
                hsr(in,3) = dcmplx(0.0d+00,0.0d+00)
                hsr(in,4) = dcmplx(0.0d+00,0.0d+00)
                hsr(in,5) = dcmplx(0.0d+00,0.0d+00)
                hsl(in,1) = dcmplx(1.0d+00,0.00d+00)
                hsl(in,2) = dcmplx(0.0d+00,0.00d+00)
            else if(jbdry.eq.-1)then
c-----
c           RIGID ABOVE FLUID
c-----
                hsr(in,1) = dcmplx(0.0d+00,0.0d+00)
                hsr(in,2) = dcmplx(1.0d+00,0.0d+00)
                hsr(in,3) = dcmplx(0.0d+00,0.0d+00)
                hsr(in,4) = dcmplx(0.0d+00,0.0d+00)
                hsr(in,5) = dcmplx(0.0d+00,0.0d+00)
                hsl(in,1) = dcmplx(1.0d+00,0.00d+00)
                hsl(in,2) = dcmplx(0.0d+00,0.00d+00)
            else if(jbdry.eq.1)then
c-----
c           FREE ABOVE FLUID
c-----
                hsr(in,1) = dcmplx(1.0d+00,0.0d+00)
                hsr(in,2) = dcmplx(0.0d+00,0.0d+00)
                hsr(in,3) = dcmplx(0.0d+00,0.0d+00)
                hsr(in,4) = dcmplx(0.0d+00,0.0d+00)
                hsr(in,5) = dcmplx(0.0d+00,0.0d+00)
                hsl(in,1) = dcmplx(1.0d+00,0.00d+00)
                hsl(in,2) = dcmplx(0.0d+00,0.00d+00)
            endif
        endif
        return
        end

        subroutine var(p,q,ra,rb,w,x,y,z,cosp,cosq,ex,
     1      exa,exl,yl,zl,cosql,iwat)
c     not modified for negative p,q
c     this assumes that real p and real q have same signs
        common/ovrflw/a0,cpcq,cpy,cpz,cqw,cqx,xy,xz,wy,wz
        complex*16 cpcq,cpy,cpz,cqw,cqx,xy,xz,wy,wz
        complex*16 p,q,ra,rb,w,x,y,z,cosp,cosq
        complex*16 yl,zl,cosql
        complex*16 eqp,eqm,epp,epm,sinp,sinq
        real *8 a0,pr,pi,qr,qi,fac,qmp,ex,exa,exl
c-----
c       form terms such as cos(p), sin(p), cos(q), sin(q)
c       and cos(p)*cos(q)
c
c       Introduce a factorization of exponentials to
c       make a pseudo floating point system
c
c       ex is the exponent in cosp
c       exl is the exponent in cosq for SH
c       exa is the exponent in cosp*cosq
c-----
        real*8 DREAL
      ex=0.0d+00
      exl = 0.0d+00
      a0=0.0d+00
      pr=dreal(p)
      pi=dimag(p)
      epp=dcmplx(dcos(pi),dsin(pi))/2.
      epm=dconjg(epp)
      ex=pr
      fac=0.0
      if(pr.lt.15.) fac=dexp(-2.*pr)
      cosp=epp + fac*epm
      sinp=epp - fac*epm
      w=sinp/ra
      x=ra*sinp
        if(iwat.eq.1)then
c-----
c       fluid layer
c-----
            a0 = 1.0d+00
            exa = ex
            cosq = 1.0d+00
            y = 0.0d+00
            z = 0.0d+00
            cosql = 1.0d+00
            yl = 0.0d+00
            zl = 0.0d+00
            exl = 0.0d+00
        else
c-----
c       elastic layer
c-----
            qr=dreal(q)
            qi=dimag(q)
            eqp=dcmplx(dcos(qi),dsin(qi))/2.
            eqm=dconjg(eqp)
            exl=qr
            fac=0.0d+00
            if(qr.lt.15.) fac=dexp(-2.*qr)
            cosql=eqp + fac*eqm
            sinq=eqp - fac*eqm
            yl=sinq/rb
            zl=rb*sinq
c-----
c       form factors for compound P-SV matrix
c-----
            exa=pr + qr
            cpcq=cosp*cosql
            cpy=cosp*yl
            cpz=cosp*zl
            cqw=cosql*w
            cqx=cosql*x
            xy=x*yl
            xz=x*zl
            wy=w*yl
            wz=w*zl
            fac=0.0d+00
            qmp=qr-pr
            if(qmp.gt.-40.) fac=dexp(qmp)
            cosq=cosql*fac
            y=fac*yl
            z=fac*zl
            fac=0.0d+00
            if(exa.lt.60.) a0=dexp(-exa)
        endif
        return
        end

        subroutine modcpy(totmp) 
        logical totmp
c-----
c       copy model to temporary array
c-----
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        common/modlly/mmax
        common/modelt/dt(NL),at(NL),bt(NL),rhot(NL),
     1      mmaxt,qat(NL),qbt(NL),etapt(NL),etast(NL),
     1      frefpt(NL),frefst(NL)
c-----
c       copy to temporary array
c-----
        if(totmp)then
            do 20 i = 1,mmax 
                dt(i) = d(i)
                at(i) = a(i)
                bt(i) = b(i)
                rhot(i) = rho(i)
                qat(i) = qa(i)
                qbt(i) = qb(i)
                etapt(i) = etap(i)
                etast(i) = etas(i)
                frefpt(i) = frefp(i)
                frefst(i) = frefs(i)
   20       continue 
            mmaxt = mmax
        else
            do 30 i = 1,mmaxt 
                d(i) = dt(i)
                a(i) = at(i)
                b(i) = bt(i)
                rho(i) = rhot(i)
                qa(i) = qat(i)
                qb(i) = qbt(i)
                etap(i) = etapt(i)
                etas(i) = etast(i)
                frefp(i) = frefpt(i)
                frefs(i) = frefst(i)
   30       continue 
            mmax = mmaxt
        endif
        return 
        end 

        subroutine dezero()
        parameter (LER=0, LIN=5, LOT=6)
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        common/modlly/mmax
c-----
c       ultimately get rid of zero thickness layers - this
c       will require readjusting the model from top down, and
c       also readjusting the source and receiver indices.
c----
c       Here just guarantee that the halfspace is not of zero thickness
c-----
        dmin = 1.0e+30
        do 100 i=1,mmax-1
            if(d(i) .lt. dmin .and. d(i).gt.0.0)dmin = d(i)
  100   continue
c       if(d(mmax).le.0.0)then
c           d(mmax) = 0.1*dmin
c       endif
        return
        end

        subroutine rshof(gg,om,wvno, lmaxs, lmaxr, wvno2, om2) 
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        common/modlly/mmax
        common/bound/vmin,vamin,vamax,vbmin,vbmax
        common/damp/alpha,ieqex
c-----
c       gus - surface displacements or potentials or top of layer
c-----
        complex*16 gus(21)
        complex*16 wvno,wvno2
        complex*16 cd(5),da(4,4),fr,y(4,4)
        complex gg(21) 
        complex*16 om,fourpo,ka2,kb2 , om2, ka, kb
        complex*16 d11,d12,fl 
        complex*16 s21,s32,s14,s34,s32e,s34e 
        complex*16 s24,s33
        complex*16 atna,atnb 
        complex*16 wv4pi
        complex*16 zero
        real *8 fact,exe,exl,exel,exll,elj
        real *8 exwu
        common/cntrl/ishank,hnkarg,dstcor,dokjar
        integer*4 dstcor
        real*4 hnkarg
        logical ishank,dokjar
        common/ctlsud/spup,spdn,ssup,ssdn,dosud
        logical spup, spdn, ssup, ssdn, dosud
        complex*16 ra,rb,gam,gamm1,h
        complex*16 sdd(4), sds(4), sss(4), sep(4), svf(4), shf(4)
        complex*16 zone
        common/jout/jsrc(21) , jbdrys, jbdryh
        complex*16 haa(4,4), saa(2,2)
        real*8 fourpi
        real*8 CDABS
        complex*16 CDSQRT
        real*8 DREAL

c-----
c       Initialization
c-----
        fourpi=12.5663706d+00
        fourpo=12.5663706d+00*om*om
        zero  = dcmplx(0.0d+00,0.0d+00)
        zone  = dcmplx(1.0d+00,0.0d+00)
c-----
c       do not evaluate for wvno = 0.0
c-----
        if(CDABS(wvno).eq.0.0d+00) then
            do 102 i=1,21
                gg(i) = cmplx(0.0,0.0)
  102       continue
        else
c-----
c       process for this wavenumber and frequency
c-----
            do 101 i = 1,21
                gus(i) = zero
  101       continue
            wv4pi = 2.0d+00 * wvno / fourpi
            call aten(om,qa(lmaxs),qb(lmaxs),ka ,kb ,alpha,
     1          a(lmaxs),b(lmaxs),atna,atnb,iwats,
     2          frefp(lmaxs),frefs(lmaxs)) 
            h =(dble(rho(lmaxs)*b(lmaxs)*b(lmaxs))*atnb*atnb)
            gam=dble(b(lmaxs))*(wvno/om)
            gam = gam * atnb
            gam = dcmplx(2.0d+00,0.0d+00)*gam*gam
            gamm1 = gam - zone
            ka2=ka *ka  
            kb2=kb *kb  
            ra=CDSQRT(wvno2-ka2)
            rb=CDSQRT(wvno2-kb2)
            if(b(lmaxr).eq.0.0)then
                iwatr = 1
            else
                iwatr = 0
            endif
            if(dosud)then
                call hsupdn(haa,saa,wvno2,om,om2,wvno,lmaxs,
     1              lmaxr,lmaxs,1)
            endif
c-----
c           this follows Wang and Herrmann (1980), 
c               "A numerical study of P-, SV-, and SH-wave 
c               generation in a plane layerd medium,
c               Bull. Seism. Soc. Am. 70, 1015-1036.
c           For p-SV we use equation (8). Note below in "do 60" 
c               that the Z component (odd index) is reversed 
c               from (8) so that +z is up
c
c       using the trick of reducing the 6x6 compound matrices to 5x5, 
c       the following correspondence between what is given here
c        and that in Wang and Herrmann is in effect
c-----
c       Compound Matrix (6x6)   Here (5x5)
c
c       X|12/12 ==      cd(1)
c       X|12/13 ==      cd(2)
c       X|12/14 ==      cd(3)
c       X|12/23 ==      -cd(3)
c       X|12/24 ==      cd(4)
c       X|12/34 ==      cd(5)
c       X|12/21 = -X|12/12
c       X|12/ij = -X|12/ji
c
c       where
c
c        |12
c       X|  = X|12/ij
c        |ij
c-----
            call scoef(cd,da,fr,om,exe,exl,exwu,wvno,
     1          fl,d11,d12,exel,exll,lmaxs,lmaxr,wvno2,om2) 

c-----
c       Form X|12/ij x Zj2 in Equation 8 of Wang and Herrmann
c----
c KLUDGE to CHANGE ORDER AND ALSO GET UR CORRECT FOR WATER LAYER
            do 50 k=1,2
c-----
c       k =1 UR     k=2 UZ
c       jj=1 = UZ, jj=2 = UR
c-----
                if(k.eq.1)then
                    jj = 2
                else if(k.eq.2)then
                    jj = 1
                endif
                if(iwatr.eq.1 .and. jj.eq.2)then
c-----
c       compute Tz and then get Ur from UZ
c-----
c-----
                    j = 3
                else
                    j = k
                endif
                y(1,jj)= cd(1)*da(2,j) + cd(2)*da(3,j) 
     1              + cd(3)*da(4,j)
                y(2,jj)= -cd(1)*da(1,j) - cd(3)*da(3,j) 
     1              + cd(4)*da(4,j)
                y(3,jj)=-cd(2)*da(1,j) + cd(3)*da(2,j) 
     1              + cd(5)*da(4,j)
                y(4,jj)= -cd(3)*da(1,j) - cd(4)*da(2,j) 
     1              - cd(5)*da(3,j)
   50       continue
            if(iwatr.eq.1)then
                y(1,2) = - wvno*y(1,2)/(rho(lmaxr)*om*om)
                y(2,2) = - wvno*y(2,2)/(rho(lmaxr)*om*om)
                y(3,2) = - wvno*y(3,2)/(rho(lmaxr)*om*om)
                y(4,2) = - wvno*y(4,2)/(rho(lmaxr)*om*om)
            endif
c-----
c           evaluate different Green's functions
c           apply source terms
c----- 
c-----
c           START OF P-SV
c-----
c           First compute the DELTA displacement-stress source terms
c           for inverted model, the UZ, TR elements change, These will
c           be only those required for dipoles and forces
c           
c           Stress-displacement discontinuities for Green's functions
c           Green   dUr dUz dTz dTr
c           DD      s32     s34
c           DS  s21
c           SS              s14
c           EX      s32e        s34e
c           VF          s33
c           HF              s24
c-----
            if(iwats.eq.1)then
                s14  = zero
                s21  = zero
                s24  = zero
                s32  = zero
                s33  = zero
                s34  = zero
                s34e = zero
            else
                s14  = -wv4pi
                s21  = 2.0d+00*kb2/(dble(rho(lmaxs))*fourpo)
                s24  = -2.0d+00/fourpi
                s32  = 4.*ka2/(dble(rho(lmaxs))*fourpo)
                s33  = dcmplx(-2.0d+00/fourpi, 0.0)
                s34  = wv4pi*( (4.*ka2/kb2) - 3.0d+00)
                s34e = 2.0d+00*wv4pi*(ka2/kb2)
            endif
            s32e=2.0d+00*ka2/(dble(rho(lmaxs))*fourpo) 
c-----
c           receiver beneath the source
c-----
            if(lmaxr .gt. lmaxs)then
                s14  = - s14
                s24  = - s24
                s32  = - s32
                s32e = - s32e
                s34  = - s34
                s34e = - s34e
            endif
c-----
c           For complete wavefield computation do not
c           waste cycles computing W matrix elements
c-----
            if(.not. dosud)then
                do 61 j=1,2
c       DD
                    gus(j   )=s32 *y(2,j)+ s34*y(4,j)
c       DS
                    gus(j+ 2)=s21 *y(1,j)             
c       SS
                    gus(j+ 4)=             s14*y(4,j)
c       EX
                    gus(j+ 6)=s32e*y(2,j)+s34e*y(4,j)
c       VF
                    gus(j+ 8)=s33 *y(3,j)
c       HF
                    gus(j+10)=s24 *y(4,j)
   61           continue
            ELSE
                sdd(1) = zero
                sdd(2) = s32
                sdd(3) = zero
                sdd(4) = s34
                sds(1) = s21
                sds(2) = zero
                sds(3) = zero
                sds(4) = zero
                sss(1) = zero
                sss(2) = zero
                sss(3) = zero
                sss(4) = s14
                sep(1) = zero
                sep(2) = s32e
                sep(3) = zero
                sep(4) = s34e
                svf(1) = zero
                svf(2) = zero
                svf(3) = s33
                svf(4) = zero
                shf(1) = zero
                shf(2) = zero
                shf(3) = zero
                shf(4) = s24
c-----
c               Change the source vector for up/down going wavefields
c-----

                call svupdn(sdd,haa)
                call svupdn(sds,haa)
                call svupdn(sss,haa)
                call svupdn(sep,haa)
                call svupdn(svf,haa)
                call svupdn(shf,haa)
c-----
c           compute the response
c-----
                do 64 j=1,2
                    gus(j   ) = sdd(1)*y(1,j) + sdd(2)*y(2,j) 
     1                  + sdd(3)*y(3,j) + sdd(4)*y(4,j)
                    gus(j+ 2) = sds(1)*y(1,j) + sds(2)*y(2,j) 
     1                  + sds(3)*y(3,j) + sds(4)*y(4,j)
                    gus(j+ 4) = sss(1)*y(1,j) + sss(2)*y(2,j) 
     1                  + sss(3)*y(3,j) + sss(4)*y(4,j)
                    gus(j+ 6) = sep(1)*y(1,j) + sep(2)*y(2,j) 
     1                  + sep(3)*y(3,j) + sep(4)*y(4,j)
                    gus(j+ 8) = svf(1)*y(1,j) + svf(2)*y(2,j) 
     1                  + svf(3)*y(3,j) + svf(4)*y(4,j)
                    gus(j+10) = shf(1)*y(1,j) + shf(2)*y(2,j) 
     1                  + shf(3)*y(3,j) + shf(4)*y(4,j)
   64           continue
            ENDIF
c-----
c       invert the vertical
c-----
            do 62 j=1,12,1
                gus(j) = -gus(j)
   62       continue
c-----
c           if receiver beneath the source unflip radial
c-----
            if(lmaxr .gt. lmaxs)then
                do 63 j=2,12,2
                    gus(j) = - gus(j)
   63           continue
            endif
c-----
c           If the receiver is in the water and the source 
c           is an explosion generate pressure time history at receiver
c           Also in water layer radial time history 
c           is generated differently
c-----
            if(iwatr.eq.1)then
                gus(16) =  dble(rho(lmaxr))*gus(8) /wvno
                gus(17) =  dble(rho(lmaxr))*gus(2) /wvno
                gus(18) =  dble(rho(lmaxr))*gus(4) /wvno
                gus(19) =  dble(rho(lmaxr))*gus(6) /wvno
                gus(20) =  dble(rho(lmaxr))*gus(10) /wvno
                gus(21) =  dble(rho(lmaxr))*gus(12) /wvno
            endif
c-----
c           END OF P-SV
c-----
c           START OF SH
c-----
            if(iwats.eq.0 .and. iwatr.eq.0)then
                sds(1) = - 2.0d+00/(rho(lmaxs)*12.5663706d+00*
     1              b(lmaxs)*b(lmaxs)*atnb*atnb)
                sds(2) = zero
                sss(1) = zero
                sss(2) =  2.0d+00*wvno/12.5663706d+00
                shf(1) = zero
                shf(2) =  2.0d+00/12.5663706d+00
                if(.not. dosud )then
                    gus(13) = - ( d11*sds(1)           )
                    gus(14) = - (            d12*sss(2))
                    gus(15) = - (            d12*shf(2))
                ELSE
                    call shupdn(sds,saa)
                    call shupdn(sss,saa)
                    call shupdn(shf,saa)
                    gus(13) = - ( d11*sds(1) + d12*sds(2))
                    gus(14) = - ( d11*sss(1) + d12*sss(2))
                    gus(15) = - ( d11*shf(1) + d12*shf(2))
                ENDIF
                if(lmaxr .gt. lmaxs)then
                    gus(13) = - gus(13)
                endif
            endif
c-----
c           END OF SH
c-----
c-----
c           do final scaling for exponential
c-----
c-----
c           SV
c-----
c           fix for radial derived from vertical for fluid
c-----
            do 71 k=1,2
                elj = -exe + exl 
                fact = 0.0D+00
                if(elj.gt.-55.) fact=dexp(elj)
                do 72 i=0,10,2
                    j = i + k
                    gg(j) = ( gus(j) * fact/fr)
c-----
c           flip UZ to make vertical positive up
c-----
                    if(k.eq.1)then
                        gg(j) = -gg(j)
                    endif
   72           continue
c----
c           do pressure field
c----
                if(k.eq.1)then
                    gg(16) = - (gus(16)*fact/fr)
                    gg(17) = - (gus(17)*fact/fr)
                    gg(18) = - (gus(18)*fact/fr)
                    gg(19) = - (gus(19)*fact/fr)
                    gg(20) = - (gus(20)*fact/fr)
                    gg(21) = - (gus(21)*fact/fr)
                endif
   71       continue
c-----
c           SH
c-----
            elj=-exel+exll
            if(iwats.eq.0)then
                if(elj.gt.-55.) then
                    fact = dexp(elj)
                    gg(13)=(gus(13)*fact)/(fl)
                    gg(14)=(gus(14)*fact)/(fl)
                    gg(15)=(gus(15)*fact)/(fl)
                else
                    gg(13) = cmplx(0.0,0.0)
                    gg(14) = cmplx(0.0,0.0)
                    gg(15) = cmplx(0.0,0.0)
                endif
            else
                gg(13) = cmplx(0.0,0.0)
                gg(14) = cmplx(0.0,0.0)
                gg(15) = cmplx(0.0,0.0)
            endif
        endif
        return
        end

        subroutine hsupdn(aa,saa,wvno2,om,om2,wvno,m,
     1      lmaxr,lmaxs,isr)
c-----
c       aa  C*16    P-SV matrix
c       saa C*16    SH matrix
c       wvno2   C*16    wavenumber squared
c       om  C*16    complex angular frequency
c       wvno    C*16    wavenumber 
c       m   I   layer index for matrix computation
c               (perform wavefield separation in this layer)
c       lmaxr   I   index of receiver layer
c       lmaxs   I   index of source layer
c               (used for definition of up/down)
c       isr I   0 - compete for receiver layer
c               1 - compute for source
c-----
        complex*16 aa(4,4), saa(2,2), wvno2, om, wvno,om2
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        common/modlly/mmax
        common/damp/alpha,ieqex
        common/ctlrud/rpup,rpdn,rsup,rsdn,dorud
        logical rpup, rpdn, rsup, rsdn, dorud
        common/ctlsud/spup,spdn,ssup,ssdn,dosud
        logical spup, spdn, ssup, ssdn, dosud
        complex*16 h, gam, gamm1, ra, rb, atna, atnb, xka, xkb
        complex*16 zero, zone 
        logical pup, pdn, sup, sdn, doud
        real*8 drho
        complex*16 CDSQRT

c-----
c       define necessary constants
c-----
        zero = dcmplx(0.0d+00,0.0d+00)
        zone = dcmplx(1.0d+00,0.0d+00)
        drho = dble(rho(m))
c-----
c       if receiver is beneath source, it is necessary to redefine
c       local concept of up/down
c-----
        if(isr.eq.0)then
            if(lmaxs .gt. lmaxr)then
                pup = rpup
                pdn = rpdn
                sup = rsup
                sdn = rsdn
            else
                pup = rpdn
                pdn = rpup
                sup = rsdn
                sdn = rsup
            endif
            doud = dorud
        else
            if(lmaxs .gt. lmaxr)then
                pup = spup
                pdn = spdn
                sup = ssup
                sdn = ssdn
            else
                pup = spdn
                pdn = spup
                sup = ssdn
                sdn = ssup
            endif
            doud = dosud
        endif
        call aten(om,qa(m),qb(m),xka,xkb,
     1      alpha,a(m),b(m),atna,atnb,iwat,
     2      frefp(m),frefs(m))
        h =(dble(rho(m)*b(m)*b(m))*atnb*atnb)
        gam=dble(b(m))*(wvno/om)
        gam = gam * atnb
        gam = dcmplx(2.0d+00,0.0d+00)*gam*gam
        gamm1 = gam - zone
        ra=CDSQRT(wvno2-xka*xka)
        rb=CDSQRT(wvno2-xkb*xkb)
        do 100 j=1,4
            do 101 i = 1,4
                aa(i,j) = zero
  101       continue
  100   continue
        saa(1,1) = zero
        saa(1,2) = zero
        saa(2,1) = zero
        saa(2,2) = zero
c-----
c       now test wave field conditions
c-----
        if(pup .and. pdn .and. sup .and. sdn)then
            aa(1,1)  = zone
            aa(2,2)  = zone
            aa(3,3)  = zone
            aa(4,4)  = zone
            saa(1,1) = zone
            saa(2,2) = zone
        else
c-----
c           water layer
c-----
            if(iwat.ne.0)then
c-----
c               coefficients of exp(nua*h)
c-----
                if(pup)then
                    aa(2,2) = aa(2,2) + 0.5d+00*zone
                    aa(2,3) = aa(2,3) - 
     1                  ra/(dble(2.0*rho(m))*om2)
                    aa(3,2) = aa(3,2) - 
     1                  drho*om2/(2.0d+00*ra)
                    aa(3,3) = aa(3,3) + 0.5d+00*zone
                endif
c-----
c               coefficients of exp(-nua*h)
c-----
                if(pdn)then
                    aa(2,2) = aa(2,2) + 0.5d+00*zone
                    aa(2,3) = aa(2,3) + 
     1                  ra/(dble(2.0*rho(m))*om2)
                    aa(3,2) = aa(3,2) + 
     1                  drho*om2/(2.0d+00*ra)
                    aa(3,3) = aa(3,3) + 0.5d+00*zone
                endif
c-----
c           elastic layer
c-----
            else
c-----
c               coefficients of exp(nua*h)
c-----
                if(pup)then
                    aa(1,1) = aa(1,1) + gam
                    aa(1,2) = aa(1,2) - wvno*gamm1/ra
                    aa(1,3) = aa(1,3) - wvno/(drho*om2)
                    aa(1,4) = aa(1,4) + 
     1                  wvno2/(drho*ra*om2)
                    aa(2,1) = aa(2,1) + gam*ra/wvno
                    aa(2,2) = aa(2,2) - gamm1
                    aa(2,3) = aa(2,3) - ra/(drho*om2)
                    aa(3,1) = aa(3,1) + 
     1                  drho*gam*gamm1*om2/wvno
                    aa(3,2) = aa(3,2) - drho*om2*
     1                  gamm1*gamm1/ra
                    aa(4,1) = aa(4,1) + drho*om2*
     1                  gam*gam*ra/wvno2
                endif
c-----
c               coefficients of exp(-nua*h)
c-----
                if(pdn)then
                    aa(1,1) = aa(1,1) + gam
                    aa(1,2) = aa(1,2) + wvno*gamm1/ra
                    aa(1,3) = aa(1,3) - wvno/(drho*om2)
                    aa(1,4) = aa(1,4) - 
     1                  wvno2/(drho*ra*om2)
                    aa(2,1) = aa(2,1) - gam*ra/wvno
                    aa(2,2) = aa(2,2) - gamm1
                    aa(2,3) = aa(2,3) + ra/(drho*om2)
                    aa(3,1) = aa(3,1) + 
     1                  drho*gam*gamm1*om2/wvno
                    aa(3,2) = aa(3,2) + drho*om2*
     1                  gamm1*gamm1/ra
                    aa(4,1) = aa(4,1) - drho*om2*
     1                  gam*gam*ra/wvno2
                endif
c-----
c               coefficients of exp(nub*h)
c-----
                if(sup)then
                    aa(1,1) = aa(1,1) - gamm1
                    aa(1,2) = aa(1,2) + gam*rb/wvno
                    aa(1,3) = aa(1,3) + wvno/(drho*om2)
                    aa(1,4) = aa(1,4) - rb/(drho*om2)
                    aa(2,1) = aa(2,1) - wvno*gamm1/rb
                    aa(2,2) = aa(2,2) + gam
                    aa(2,3) = aa(2,3) + 
     1                  wvno2/(drho*rb*om2)
                    aa(3,1) = aa(3,1) - 
     1                  drho*gam*gamm1*om2/wvno
                    aa(3,2) = aa(3,2) + 
     1                  drho*om2*gam*gam*rb/wvno2
                    aa(4,1) = aa(4,1) - drho*om2*
     1                  gamm1*gamm1/rb
                    saa(1,1) = saa(1,1) + zone
                    saa(1,2) = saa(1,2) + zone/(h*rb)
                    saa(2,1) = saa(2,1) + (h*rb)
                    saa(2,2) = saa(2,2) + zone
                endif
c-----
c               coefficients of exp(-nub*h)
c-----
                if(sdn)then
                    aa(1,1) = aa(1,1) - gamm1
                    aa(1,2) = aa(1,2) - gam*rb/wvno
                    aa(1,3) = aa(1,3) + wvno/(drho*om2)
                    aa(1,4) = aa(1,4) + rb/(drho*om2)
                    aa(2,1) = aa(2,1) + wvno*gamm1/rb
                    aa(2,2) = aa(2,2) + gam
                    aa(2,3) = aa(2,3) - 
     1                  wvno2/(drho*rb*om2)
                    aa(3,1) = aa(3,1) - 
     1                  drho*gam*gamm1*om2/wvno
                    aa(3,2) = aa(3,2) - 
     1                  drho*om2*gam*gam*rb/wvno2
                    aa(4,1) = aa(4,1) + drho*om2*
     1                  gamm1*gamm1/rb
                    saa(1,1) = saa(1,1) + zone
                    saa(1,2) = saa(1,2) - zone/(h*rb)
                    saa(2,1) = saa(2,1) - (h*rb)
                    saa(2,2) = saa(2,2) + zone
                endif
                aa(2,4) = -aa(1,3)
                aa(3,3) =  aa(2,2)
                aa(3,4) = -aa(1,2)
                aa(4,2) = -aa(3,1)
                aa(4,3) = -aa(2,1)
                aa(4,4) =  aa(1,1)
c-----
c       special case clean up
c-----
                if(pup.and.sup .or. pdn.and.sdn)then
                    aa(1,1) = zone
                    aa(2,2) = zone
                    aa(3,3) = zone
                    aa(4,4) = zone
                    aa(1,3) = zero
                    aa(2,4) = zero
                    aa(3,1) = zero
                    aa(4,2) = zero
                endif
                do 200 j=1,4
                    do 201 i = 1,4
                        aa(i,j) = aa(i,j)/dble(2.0)
  201               continue
  200           continue
                saa(1,1) = saa(1,1) / dble(2.0)
                saa(1,2) = saa(1,2) / dble(2.0)
                saa(2,1) = saa(2,1) / dble(2.0)
                saa(2,2) = saa(2,2) / dble(2.0)
            endif
        endif
        return
        end

        subroutine scoef(cd,da,fr,om,exe,exl,exwu,wvno,
     1      fl,d11,d12,exel,exll,llmaxs,llmaxr,wvno2, om2)
c-----
c       llmaxs  I*4 - source layer index in original model
c       llmaxr  I*4 - receiver layer index in original model
c-----
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        common/modlly/mmax
        common/damp/alpha,ieqex
        complex*16 om, om2
        complex*16  da(4,4),ca(5,5)
        complex*16   cd(5),e(5),fr
        complex*16 d11,d12,e1,e2,e21, e22, fl
        real *8 exe,exl,exel,exll,ex,exa,exb,exwu
        real*8 dzero
        complex*16 wvno2, wvno
        complex*16 zdum
        complex*16 aa(4,4)
        complex*16 cy(5,5)
        complex*16 zero, zone
        complex*16 y11, y12, y21, y22, sd11, sd21
        complex*16 hl(2,2)
c-----
c       matrix components in layers and boundaries saved
c-----
        complex*16 har(NL,4,4), dar(NL,5,5), hsr(2,5), gbr(2,5), 
     1      hal(NL,2,2), hsl(2,2), gbl(2,2)
        real*8 hex(NL), lex(NL), dex(NL), hexw(NL)
        common/hamat/har
        common/damat/dar
        common/hsrfr/hsr
        common/gbrfr/gbr
        common/hlmat/hal
        common/hsrfl/hsl
        common/gbrfl/gbl
        common/hexex/hex
        common/hexexw/hexw
        common/dexex/dex
        common/lexex/lex 
        common/water/iwater(NL),iwats(2),iwatb(2)
        common/updnsm/equalu(NL), equald(NL)
        logical equalu, equald
        logical retrieve
c-----
c       check for decomposition at wavefield at receiver
c-----
        common/ctlrud/rpup,rpdn,rsup,rsdn,dorud
        logical rpup, rpdn, rsup, rsdn, dorud
        common/ctlsud/spup,spdn,ssup,ssdn,dosud
        logical spup, spdn, ssup, ssdn, dosud
        complex*16 haa(4,4), saa(2,2)

C       COMMON/DEBUG/VERBY
C       LOGICAL VERBY
c-----
c-----
c       this routine computes the layer response. 
c          To simplify the mathematics
c       of the case of receiver above or beneath the source, the
c       layer is internally flipped.
c
c       llmaxs  I*4 - source layer index in original model
c       llmaxr  I*4 - receiver layer index in original model
c       mm  I*4 - pointer to layer in original model
c
c       lmaxs   I*4 - source layer index in current model
c       lmaxr   I*4 - receiver layer index in current model
c       mm  I*4 - pointer to layer in current model
c
c       in  I*4 - 1 use current model (source beneath receiver)
c               - 2 use inverted model (receiver beneath source)
c-----
c       initialize matrices
c-----
        zero = dcmplx(0.0d+00,0.0d+00)
        zone  = dcmplx(1.0d+00,0.0d+00)
        dzero = 0.0d+00
        exe=0.0
        exl=0.0
        exwu = 0.0
        do 2 j = 1,4
            do 3 i = 1,4
                da(i,j)=zero
    3       continue
            da(j,j) = zone
    2   continue
        do 12 j=1,5
            do 13 i=1,5
                cy(i,j) = zero
   13       continue
            cy(j,j) = zone
   12   continue
        y11 = zone
        y12 = zero
        y21 = zero
        y22 = zone
        exel = 0.0
        exll = 0.0
c-----
c     set up halfspace conditions
c-----
        if(llmaxs .ge. llmaxr)then
            in = 1
        else
            in = 2
        endif
        do 100 i=1,5
            e(i) = gbr(in,i)
  100   continue
        if(dreal(om).lt.0.12.and.dreal(om).gt.0.06)then
        endif
        e1 = gbl(in,1)
        e2 = gbl(in,2)
        do 11 i=1,5
            cd(i)=e(i)
   11   continue
        d11=e1
        d12=e2
c-----
c       set up limits on the layer stacking
c-----
        if(llmaxs .ge. llmaxr)then
            lmaxs = llmaxs
            lmaxr = llmaxr
        else
            lmaxs = mmax - llmaxs + 2
            lmaxr = mmax - llmaxr + 2
        endif
c-----
c       matrix multiplication from bottom layer upward
c-----
        do 1340 mm = mmax,1,-1
            if(llmaxs .ge. llmaxr)then
                m = mm
                if(equalu(m))then
                    retrieve = .false.
                else
                    retrieve = .true.
                endif
            else
                m = mmax + 1 - mm
                if(equald(m))then
                    retrieve = .false.
                else
                    retrieve = .true.
                endif
            endif
            iwat = iwater(m)
            if(retrieve)then
                call copy5(ca,dar,m,1,dex,exa)
                call copy2(hl,hal,m,1,lex,exb)
                call copy4(aa,har,m,1,hex,ex)
            endif

            call cmult(e,ca,exa,exe)
            call lmult(e1,e2,e21,e22,hl,iwat,exel,exb,.false.)
            if(mm.lt.lmaxr)then
                call rcmult(cy,ca,exa,exl)
                call lmult(y11,y12,y21,y22,hl,iwat,
     1              exll,exb,.true.)
            else if(mm.ge.lmaxr .and. mm.lt.lmaxs) then
                call dmult(da,aa)
                    exl = exl + ex
c-----
c       save values at top of source layer
c-----
            else if(mm.eq.lmaxs) then
                    do 1352 i=1,5
                    cd(i)=e(i)
 1352           continue
                exl=exe
                exll = exel
                d11=e1
                d12=e2
            endif
            if(mm.eq.1)then
                do 200 i=1,5
                    ca(i,1) = hsr(in,i)
  200           continue
                sd11 = hsl(in,1)
                sd21 = hsl(in,2)
                zdum = e1
                e1 = zdum*sd11 + e2*sd21
c               e2 = zdum*sd11 - e2*sd21
                zdum = y11
                y11 = zdum*sd11 + y12*sd21
                zdum = y21
                y21 = zdum*sd11 + y22*sd21
                zdum = dcmplx(0.0,0.0)
                do 1402 i=1,5
                    zdum = zdum + e(i)*ca(i,1)
 1402           continue
                e(1) = zdum
                call rcmult(cy,ca,dzero,exl)
            endif
 1340 continue
c-----
c       get final matrices
c-----
c-SH
        fl=e1
c-P-SV
c       form x(l,m)y(ij|12)
c-----take care of x(i,j) y(1j|12) and replace the da
            aa(1,1) =   zero
            aa(2,1) =   cy(1,1)
            aa(3,1) =   cy(2,1)
            aa(4,1) =   cy(3,1) /(2.0d+00 )
c change sign 0430 1200
            aa(1,2) = - aa(2,1)
            aa(2,2) =   zero
            aa(3,2) = - aa(4,1)
            aa(4,2) =   cy(4,1)
            aa(1,3) = - aa(3,1)
            aa(2,3) = - aa(3,2)
            aa(3,3) =   zero
            aa(4,3) =   cy(5,1)
            aa(1,4) = - aa(4,1)
            aa(2,4) = - aa(4,2)
            aa(3,4) = - aa(4,3)
            aa(4,4) =   zero
            call dmult(da,aa)
        fr=e(1)
c-----
c       if decomposion of wavefield at receiver modify the
c       returned matrices
c-----
        if(dorud)then
            
            call hsupdn(haa,saa,wvno2,om,om2,wvno,llmaxr, 
     1          llmaxr,llmaxs,0)
c-----
c           SH
c           B = SAA B
c-----
            d11 = d11*(saa(1,1)*y11 + saa(1,2)*y21)
            d12 = d12*(saa(1,1)*y11 + saa(1,2)*y21)
c-----
c           SV
c           BT = ST CD X HAAT
c           Form X HAAT
c-----
            call trans4(haa)
            call dmult(da,haa)
        else
            d11 = y11*d11
            d12 = y11*d12
        endif
        return
        end

        subroutine svupdn(s,haa)
        complex*16 s(4), haa(4,4)
        complex*16 t(4), tmp
        do 100 i=1,4
            tmp = dcmplx(0.0d+00,0.0d+00)
            do 101 j=1,4
                tmp = tmp + haa(i,j)*s(j)
  101       continue
            t(i) = tmp
  100   continue
        do 200 i=1,4
            s(i) = t(i)
  200   continue
        return
        end

        subroutine shupdn(s,saa)
        complex*16 s(2), saa(2,2)
        complex*16 tmp
            tmp = saa(1,1)*s(1) + saa(1,2)*s(2)
            s(2) = saa(2,1)*s(1) + saa(2,2)*s(2)
            s(1) = tmp
        return
        end

        subroutine trans4(a)
c-----
c       from a transpose 
c-----
        complex*16 a(4,4), zdum
        do 100 i=1,4
            do 101 j=i,4
                zdum = a(i,j)
                a(i,j) = a(j,i)
                a(j,i) = zdum
  101       continue
  100   continue
        return
        end

        subroutine chkmod()
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        common/modlly/mmax
        common/jout/jsrc(21) , jbdrys, jbdryh
        common/lwater/lfluid
        logical lfluid
c-----
c       check model for inconsistencies
c-----
c-----
c       Model cannot consist entirely of water layers
c       Also determine first solid layer from top
c-----
        iw = 0  
        isoldd = 0
        do 100 i=1,mmax
            if(b(i).eq.0.0)then
                iw = iw + 1
            else
                if(isoldd .eq.0)isoldd=i
            endif
  100   continue
        if(iw .eq. mmax)then
            lfluid = .true.
C           call werror('MODEL CONSISTS ONLY OF LIQUID LAYERS')
        else
            lfluid = .false.
        endif
c-----
c       Determine first solid layer from bottom
c-----
        iw = 0  
        isoldu = 0
        do 101 i=mmax,1,-1
            if(b(i).eq.0.0)then
                iw = iw + 1
            else
                if(isoldu .eq.0)isoldu=i
            endif
  101   continue
c-----
c       Check for interior water layer
c-----
        if(iw.gt.0 .and. .not. lfluid)then
            do 102 i=isoldd,isoldu
                if(b(i).eq.0.0)then
                call werror('MODEL HAS INTERIOR  FLUID LAYERS')
                endif
  102       continue
        endif
c-----
c       If boundary condition is rigid, and the adjacent layer is
c       fluid, terminate 
c-----
C       if(b(1).le.1.0e-04 .and. jbdrys.eq.-1 .and. lfluid)then
C           call werror('TOP LAYER IS FLUID AND RIGID BOUNDARY')
C       endif
C       if(b(mmax).le.1.0e-04 .and. jbdryh.eq.-1 .and. lfluid)then
C           call werror('BOTTOM LAYER IS FLUID AND RIGID BOUNDARY')
C       endif
        return
        end

        subroutine werror(ostr)
c-----
c       output error message and terminate program
c-----
        parameter(LER=0, LIN=5, LOT=6)
        character ostr*(*)
        write(LER,*)'PROGRAM TERMINATION'
        write(LER,*)ostr
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
