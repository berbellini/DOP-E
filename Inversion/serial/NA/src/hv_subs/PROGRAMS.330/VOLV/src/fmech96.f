        program fmech96
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME V                                                       c
c                                                                     c
c      PROGRAM: FMECH96                                               c
c                                                                     c
c      COPYRIGHT 1996 R. B. Herrmann                                  c
c                                                                     c
c      Department of Earth and Atmospheric Sciences                   c
c      Saint Louis University                                         c
c      221 North Grand Boulevard                                      c
c      St. Louis, Missouri 63103                                      c
c      U. S. A.                                                       c
c---------------------------------------------------------------------c
c-----
c       11 SEP 2000 add TP, TSV TSH to FILE96 format
c       14 JAN 2001 add A, C, F, L, N and density information
c       30 JAN 2002 ensure that AZ is written into f96 header
c       08 JUN 2002 changed e20.13 to e20.0 format in subroutine chtofp
c       09 JAN 2007 put in stevaz correctly based on the BAZ
c       23 MAY 2009 In an effort to both comply with and promulgate
c          international standards, the NEIC Policy and Procedures
c          Working Group has made the decision to adopt the
c          international standard for converting scalar moment
c          to moment magnitude proposed by the IASPEI Commission on
c          Seismological Observation and Interpretation (CoSOI)
c          Working Group on Magnitude --
c
c         old     Mw = 2/3 log Mo - 10.7    Mo in dyne-cm
c         new     Mw = 2/3(log Mo - 16.1)   Mo in dyne-cm
c                 Mw = 2/3(log Mo -  9.1)   Mo in N-m
c       23 MAY 2009 - correctly set cmpinc and cmpaz for the horizontals
c         for ZRT
c       12 JUN 2009 - NOTE if this is to be used for TI media be very
c         careful about the sue of SA AN SL SC for the dislocation
c-----
        parameter (LER=0, LIN=5, LOT=6)
c-----
c       variables in common blocks
c
c       iftype  I*4 File type
c               1 - single trace
c               3 - three component
c               16 - Green's function
c               21 - Green's function
c
c       iobsyn  I*4 1 - observed
c               2 - synthetic
c
c       itmfrq  I*4 1 - time series
c               2 - Fourier spectra (not implemented)
c
c       iunit   I*4 1   - counts 
c               2   - cm
c               3   - cm/sec
c               4   - cm/sec/sec
c               5   - m
c               6   - m/sec
c               7   - m/sec/sec
c               8   - microns
c               9   - microns/sec
c               10  - microns/sec/sec, junit
c       junit   I*4
c               11  - Pa  (nt/m^2)
c               12  - MPa  (mega Pascal)
c
c       cfilt   C*80    comment on filtering operations
c
c       keyear  I*4 event year
c       kemon   I*4 event mon
c       keday   I*4 event day
c       kehour  I*4 event hour
c       kemin   I*4 event minute
c       esec    R*4 event second
c       evlat   R*4 event latitude
c       evlon   R*4 event longitude
c       evdep   R*4 event depth
c       
c       stname  C*8 station name
c               NOTE: AH(6), SEED(5), CSS(6), SAC(8)
c       
c       stlat   R*4 station latitude
c       stlon   R*4 station longitude
c       stelev  R*4 station elevation
c
c
c       distkm  R*4 epicentral distance in km
c       distdg  R*4 epicentral distance in degrees (111.195 km/degree)
c       evstaz  R*4 event -> station azimuth, degrees east of north
c       stevaz  R*4 station -> event azimuth, degrees east of north
c       
c       cpulse  C*80    pulse description
c
c       ccomnt  C*80    comment
c
c       jsrc    I*4 Array of indices to indicate if traces are 
c                   present
c               iftype =  1  jsrc(1) indicates if trace is 
c                   present
c               iftype =  3  jsrc(i),i=1,5 indicates if trace 
c                   is present
c                   i = Z, 2=N, 3=E, 4=R, 5=T
c                   (we could carry all five, but the real 
c                   purpose  is to permit no more than 
c                   3 traces, but 2 may all that are 
c                   present in a real data set)
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
     2      TP, TSV, TSH, SA, SC, SF, SL, SN, SR
        real*4 esec, distkm, distdg, evstaz, stevaz
        real*4 evlat, evlon, evdep 
        real*4 stlat, stlon, stelev 
        REAL*4 TP, TSV, TSH
        REAL*4 SA, SC, SF, SL, SN, SR

        common/chdr96/ stname, cfilt, cpulse, ccomnt
        character*8 stname*8, cfilt*80, cpulse*80, ccomnt*80    

        character stcomp*8
        real*4 cmpaz, cmpinc, cmpdt
        integer*4 npts
        real*4 ssec
c-----
c       command line arguments
c-----
        real*4 dip, rake, strike, xmom, az ,baz
        real*4 forcex, forcey, forcez
        real*4 xmt(3,3)
        logical verby
        logical lrot
        integer isds
c-----
c       internal program variables
c-----
        integer*4 NSAMP
        parameter (NSAMP=16384)
        common/xx/x
        real*4 x(NSAMP)
        common/xxr/xr
        real*4 xr(NSAMP)
        common/xxt/xt
        real*4 xt(NSAMP)
        common/xxz/xz
        real*4 xz(NSAMP)
        common/xxp/xp
        real*4 xp(NSAMP)
        real*4 fz(4),fr(4),ft(4)
        integer jjsrc(21)
        logical lhavez, lhaver, lhavet
c-----
c       call machine dependent initialization
c-----
        call mchdep()
c-----
c       parse command line arguments
c-----
        call gcmdln(dip,rake,strike,isds,verby,xmom,az,baz,lrot,
     1      xmt,forcex,forcey,forcez)
c-----
c       INITIALIZE
c-----
        degrad = 0.01745329
        sint = sin(baz*degrad)
        cost = cos(baz*degrad)
        sina = sin(az*degrad)
        cosa = cos(az*degrad)
        sin2a = sin(2.*az*degrad)
        cos2a = cos(2.*az*degrad)
            do 1002 i=1,4
                fz(i) = 0.0
                fr(i) = 0.0
                ft(i) = 0.0
 1002       continue
        if(isds .eq.0)then
            if(verby)then
                write(LER,*)'Double couple:'
                write(LER,*)'   Dip =',dip,' Rake= ',rake,
     1              ' Strike= ',strike
                write(LER,*)'   Moment= ',xmom,' az=',az,
     1              ' baz=',baz
                write(LER,*)'   Rotation = ', lrot
            endif
c-----
c       STRIKE DIP SLIP MOMENT SPECIFIED
c-----
            str=strike
            slp=rake
            call trans(dip,str,slp,f1,f2,f3,v1,v2,v3)
            fz(1) = f3 * v3
            fr(1) = fz(1)
            fz(2) = (f1*v3+f3*v1)*cosa + (f2*v3+f3*v2)*sina
            fr(2) = fz(2)
            fz(3) = (f1*v1-f2*v2)*cos2a + (f1*v2+f2*v1)*sin2a
            fr(3) = fz(3)
            ft(1) = 0.0
            ft(2) = (f1*v3+f3*v1)*sina - (f2*v3+f3*v2)*cosa
            ft(3) = (f1*v1-f2*v2)*sin2a - (f1*v2+f2*v1)*cos2a
            ft(4) = 0.0
            do 2001 i=1,4
                fz(i) = fz(i) * xmom
                fr(i) = fr(i) * xmom
                ft(i) = ft(i) * xmom
 2001       continue
        else if(isds.eq.1)then
            if(verby)then
                write(LER,*)'Moment tensor:'
                write(LER,*)'   Mxx =',xmt(1,1),' Mxy = ',
     1              xmt(1,2),' Mxz = ',xmt(1,3)
                write(LER,*)'   Myy =',xmt(2,2),' Myz = ',
     1              xmt(2,3),' Mzz = ',xmt(3,3)
                write(LER,*)'   az=',az,' baz=',baz
                write(LER,*)'   Rotation = ', lrot
            endif
c-----
c       MOMENT TENSOR SPECIFIED
c-----
            fz(1) = -(xmt(1,1)+xmt(2,2))/6.0 + xmt(3,3)/3.0
            fr(1) = fz(1)
            fz(2) = xmt(1,3)*cosa + xmt(2,3)*sina
            fr(2) = fz(2)
            fz(3) = 0.5*(xmt(1,1)-xmt(2,2))*cos2a + xmt(1,2)*sin2a
            fr(3) = fz(3)
            fz(4) = (xmt(1,1)+xmt(2,2)+xmt(3,3))/3.0
            fr(4) = fz(4)
            ft(1) = 0.0
            ft(2) = -xmt(2,3)*cosa + xmt(1,3)*sina
            ft(3) = 0.5*(xmt(1,1)-xmt(2,2))*sin2a - xmt(1,2)*cos2a
            ft(4) = 0.0
        else if(isds.eq.2)then
            if(verby)then
                write(LER,*)'Explosion'
                write(LER,*)'   az=',az,' baz=',baz,'exp=',isds
                write(LER,*)'   Rotation = ', lrot
            endif
c-----
c       EXPLOSION SPECIFIED
c-----
            fz(4) = xmom
            fr(4) = xmom
        endif
        fzh = (forcex*cosa + forcey*sina)
        fth = (forcex*sina - forcey*cosa)

c-----
c following the notation of Wang and Herrmann (1980)
c-----
C       if(verby)then
C           write(LER,*)' fz=',fz
C           write(LER,*)' fr=',fr
C           write(LER,*)' ft=',ft
C       endif
c-----
c       process Green's Functions
c-----
 9995   continue
        lhavez = .false.
        lhaver = .false.
        lhavet = .false.
        call rdhd96(LIN,nerr)
        if(nerr .lt. 0)go to 9999
c-----
c       redefine the azimuth
c-----
        evstaz = az
        stevaz = baz
c-----
c       save the input field description, since the output will
c       be three component and the input must be 16 or 21
c-----
        iftyin = iftype
        if(iftyin .ne. 16 .and. iftyin.ne.21)go to 9999
c-----
c       copy the jsrc array, and reset
c-----
        do 1004 i=1,21
            if(i.le.iftyin)then
                jjsrc(i) =jsrc(i)
                jsrc(i)  = 0
            else
                jjsrc(i) = 0
                jsrc(i) = 0
            endif
 1004   continue
        if(jjsrc(1) .ne. 0)lhavez = .true.
        if(jjsrc(3) .ne. 0)lhavez = .true.
        if(jjsrc(6) .ne. 0)lhavez = .true.
        if(jjsrc(9) .ne. 0)lhavez = .true.
        if(jjsrc(11) .ne. 0)lhavez = .true.
        if(jjsrc(13) .ne. 0)lhavez = .true.

        if(jjsrc(2) .ne. 0)lhaver = .true.
        if(jjsrc(4) .ne. 0)lhaver = .true.
        if(jjsrc(7) .ne. 0)lhaver = .true.
        if(jjsrc(10) .ne. 0)lhaver = .true.
        if(jjsrc(12) .ne. 0)lhaver = .true.
        if(jjsrc(14) .ne. 0)lhaver = .true.

        if(jjsrc(5) .ne. 0)lhavet = .true.
        if(jjsrc(8) .ne. 0)lhavet = .true.
        if(jjsrc(15) .ne. 0)lhavet = .true.

        if(lhavez)jsrc(1) = 1
        if(lhaver)jsrc(2) = 4
        if(lhavet)jsrc(3) = 5
        if(lhavet .or. lhaver)then
            if(.not. lrot)then
                jsrc(2) = 2
                jsrc(3) = 3
            endif
        endif

        do 1005 i=1,NSAMP
            xz(i)=0.0
            xr(i)=0.0
            xt(i)=0.0
            xp(i)=0.0
 1005   continue
c-----
c       read in the 16 basic solutions
c-----
        do 1001 j=1,iftype
c-----
c       initialize input array
c-----
            do 300 i=1,NSAMP
                x(i)=0.0
  300       continue
c-----
c       PUT XMOM ABOVE IN COEFFICIENTS
c-----
            if(jjsrc(j).ne.0)then
                call rdtr96(LIN,stcomp,cmpinc,cmpaz,
     1              cmpdt, npts, ksyear, ksmon, 
     2              ksday, kshour, ksmin, ssec, 
     3              x,nerr,NSAMP)
                    if(nerr .lt. 0 )go to 9999
                do 9200 i=1,npts
                    if(j.eq.1)then
                        xz(i)=xz(i) + fz(1)*x(i)
                    elseif(j.eq.2)then
                        xr(i)=xr(i) + fr(1)*x(i)
                    elseif(j.eq.3)then
                        xz(i)=xz(i) + fz(2)*x(i)
                    elseif(j.eq.4)then
                        xr(i)=xr(i) + fr(2)*x(i)
                    elseif(j.eq.5)then
                        xt(i)=xt(i) + ft(2)*x(i)
                    elseif(j.eq.6)then
                        xz(i)=xz(i) + fz(3)*x(i)
                    elseif(j.eq.7)then
                        xr(i)=xr(i) + fr(3)*x(i)
                    elseif(j.eq.8)then
                        xt(i)=xt(i) + ft(3)*x(i)
                    else if(j.eq.9)then
                        xz(i)=xz(i) + fz(4)*x(i)
                    else if(j.eq.10)then
                        xr(i)=xr(i) + fr(4)*x(i)
                    else if(j.eq.11)then
                        xz(i)=xz(i) + forcez*x(i)
                    else if(j.eq.12)then
                        xr(i)=xr(i) + forcez*x(i)
                    else if(j.eq.13)then
                        xz(i)=xz(i) + fzh*x(i)
                    else if(j.eq.14)then
                        xr(i)=xr(i) + fzh*x(i)
                    else if(j.eq.15)then
                        xt(i)=xt(i) + fth*x(i)
                    else if(j.eq.16)then
                        xp(i)=xp(i) + fz(4)*x(i)
                    else if(j.eq.17)then
                        xp(i)=xp(i) + fz(1)*x(i)
                    else if(j.eq.18)then
                        xp(i)=xp(i) + fz(2)*x(i)
                    else if(j.eq.19)then
                        xp(i)=xp(i) + fz(3)*x(i)
                    else if(j.eq.20)then
                        xp(i)=xp(i) + forcez*x(i)
                    else if(j.eq.21)then
                        xp(i)=xp(i) + fzh*x(i)
                    endif
 9200           continue
            endif
 1001   continue
c-----
c       rotate horizontal components into receiver coordinates N and E
c-----
        if(.not. lrot)then
            do 9210 i = 1,npts
                yr = xr(i)
                yt = xt(i)
                xr(i) = - cost*yr + sint*yt
                xt(i) = -sint*yr - cost*yt
 9210       continue
        endif
c-----
c       output the header, and then the traces
c-----  
        iftype = 3
c-----
c       redefine the jsrc array
c-----
        call wrhd96(LOT,nerr)
c-----
c       output the Vertical Component
c-----
                    if(jsrc(1).ne.0)then
                    stcomp = 'Z       '
                    cmpinc = -90.0
                    cmpaz = 0.0
                    call wrtr96(LOT,stcomp,cmpinc,cmpaz,
     1                  cmpdt, npts, ksyear, ksmon, 
     2                  ksday, kshour, ksmin, ssec, 
     3                  xz,nerr,NSAMP)
                endif
c-----
c       output the Radial/North Component
c-----
                if(jsrc(2).ne.0)then
                    cmpinc = 0.0
                    if(lrot)then
                        stcomp = 'R       '
                        cmpaz = az 
                        cmpaz = amod(cmpaz,360.0)
                        if(cmpaz.lt.0.0)cmpaz=cmpaz+360.0
                    else
                        stcomp = 'N       '
                        cmpaz = 0.0
                    endif
                    call wrtr96(LOT,stcomp,cmpinc,cmpaz,
     1                  cmpdt, npts, ksyear, ksmon, 
     2                  ksday, kshour, ksmin, ssec, 
     3                  xr,nerr,NSAMP)
                endif
c-----
c       output the Transverse/East Component
c-----
                if(jsrc(3).ne.0)then
                    if(lrot)then
                        cmpaz = az  + 90.0
                        cmpaz = amod(cmpaz,360.0)
                        if(cmpaz.lt.0.0)cmpaz=cmpaz+360.0
                        stcomp = 'T       '
                    else
                        stcomp = 'E       '
                        cmpaz = 90.0
                    endif
                    call wrtr96(LOT,stcomp,cmpinc,cmpaz,
     1                  cmpdt, npts, ksyear, ksmon, 
     2                  ksday, kshour, ksmin, ssec, 
     3                  xt,nerr,NSAMP)
                endif
        go to 9995
 9999   continue
        end

        subroutine trans(dip,stk,rake,f1,f2,f3,v1,v2,v3)
            degrad=0.01745329
            sins=sin(stk*degrad)
            coss=cos(stk*degrad)
            sind=sin(dip*degrad)
            cosd=cos(dip*degrad)
            sinr=sin(rake*degrad)
            cosr=cos(rake*degrad)
            a11=cosr*coss + sinr*cosd*sins
            a12=cosr*sins - sinr*cosd*coss
            a13= - sinr*sind
            a21= -sins*sind
            a22= coss*sind
            a23= - cosd
            f1=a11
            f2=a12
            f3=a13
            v1=a21
            v2=a22
            v3=a23
        return
        end

        subroutine gcmdln(dip,rake,strike,isds,verby,xmom,
     1      az,baz,lrot,x,fx,fy,fz)
        parameter (LER=0, LIN=5, LOT=6)
c-----
c       parse command line arguments
c
c       requires subroutine mgtarg() and function mnmarg()
c
c-----
c       dip R*4 - dip of fault
c       rake    R*4 - rake angle of motion on fault
c       strike  R*4 - strike of fault
c       isds    I*4 - indicator of couple source description
c                 -1 none given
c                  0 strike, dip, rake
c                  1 moment tensor
c                  2 explosion
c       verby   L   - verbose output on standard error
c       xmom    R*4 - seismic moment of strike, dip, rake form
c                 1.0 == 1.0 dyne-cm for km,gm,km/s system
c       az  R*4 - source -> receiver azimuth
c       baz R*4 - receiver -> source azimuth
c       lrot    L   - false, form north and east components using 
c                   back azimuth
c                 true, output radial and transverse motion
c       x(3,3)  R*4 - moment tensor components (units of dyne-cm)
c       fx  R*4 - point force components
c       fy  R*4   1.0 == 1.0  dyne for km,gm,km/s system
c       fz  R*4
c       scale   R*4 - overall amplitude scaling to facilitate
c                 input of x(3,3)
c-----
        real*4 x(3,3)
        character*25 name
        logical lrot, verby
        integer*4 mnmarg
        dip = 0.0
        rake = 0.0
        strike = 0.0
        isds=-1
        verby=.false.
        xmom=1.0
        az=0.0
        baz=0
        fx = 0.0
        fy = 0.0
        fz = 0.0
        do 120 i=1,3
            do 121 j=1,3
                x(j,i) = 0.0
  121       continue
  120   continue
        lrot = .false.
        nmarg=mnmarg()
        if(nmarg.le.0)then
            call usage()
        endif
        i = 0
   11   i = i + 1
        if(i.gt.nmarg)go to 13
            call mgtarg(i,name)
            if(name(1:2).eq.'-D')then
                i=i+1
                call mgtarg(i,name)
                call chtofp(name,dip)
                isds = 0
            else if(name(1:2).eq.'-S')then
                i=i+1
                call mgtarg(i,name)
                call chtofp(name,strike)
                isds = 0
            else if(name(1:2).eq.'-R'.and.name(1:3).ne.'-RO')then
                i=i+1
                call mgtarg(i,name)
                call chtofp(name,rake)
                isds = 0
            else if(name(1:3).eq.'-M0' .or. name(1:3).eq.'-MO')then
                i=i+1
                call mgtarg(i,name)
                call chtofp(name,xmom)
            else if(name(1:3).eq.'-MW' .or. name(1:3).eq.'-Mw')then
                i=i+1
                call mgtarg(i,name)
                call chtofp(name,xmom)
                xmom = 10.**(1.5*xmom + 16.10)
            else if(name(1:2).eq.'-A')then
                i=i+1
                call mgtarg(i,name)
                call chtofp(name,az)
            else if(name(1:2).eq.'-B')then
                i=i+1
                call mgtarg(i,name)
                call chtofp(name,baz)
            else if(name(1:2).eq.'-E')then
                isds = 2
            else if(name(1:3).eq.'-RO')then
                lrot = .true.
            else if(name(1:3).eq.'-xx' .or. name(1:3).eq.'-XX')then
                i=i+1
                call mgtarg(i,name)
                call chtofp(name,x(1,1))
                isds = 1
            else if(name(1:3).eq.'-yy' .or. name(1:3).eq.'-YY')then
                i=i+1
                call mgtarg(i,name)
                call chtofp(name,x(2,2))
                isds = 1
            else if(name(1:3).eq.'-zz' .or. name(1:3).eq.'-ZZ')then
                i=i+1
                call mgtarg(i,name)
                call chtofp(name,x(3,3))
                isds = 1
            else if(name(1:3).eq.'-xy' .or. name(1:3).eq.'-XY')then
                i=i+1
                call mgtarg(i,name)
                call chtofp(name,x(1,2))
                x(2,1) = x(1,2)
                isds = 1
            else if(name(1:3).eq.'-xz' .or. name(1:3).eq.'-XZ')then
                i=i+1
                call mgtarg(i,name)
                call chtofp(name,x(1,3))
                x(3,1) = x(1,3)
                isds = 1
            else if(name(1:3).eq.'-yz' .or. name(1:3).eq.'-YZ')then
                i=i+1
                call mgtarg(i,name)
                call chtofp(name,x(2,3))
                x(3,2) = x(2,3)
                isds = 1
            else if(name(1:3).eq.'-fx' .or. name(1:3).eq.'-FX')then
                i=i+1
                call mgtarg(i,name)
                call chtofp(name,fx)
            else if(name(1:3).eq.'-fy' .or. name(1:3).eq.'-FY')then
                i=i+1
                call mgtarg(i,name)
                call chtofp(name,fy)
            else if(name(1:3).eq.'-fz' .or. name(1:3).eq.'-FZ')then
                i=i+1
                call mgtarg(i,name)
                call chtofp(name,fz)
C           else if(name(1:2).eq.'-V')then
C               verby = .true.
            else if(name(1:2).eq.'-?')then
                call usage()
            else if(name(1:2).eq.'-h')then
                call usage()
            endif
        goto 11
   13   continue
        if(lrot)then
               baz = amod(az+180.0,360.0)
        endif
c-----
c       convert everything to the units required for the synthetics.
c       The synthetics are generated using KM for distance and
c       layer thickness, KM/SEC for velocity, and GM/CC for density.
c-----
        fx = fx / 1.0E+15
        fy = fy / 1.0E+15
        fz = fz / 1.0E+15
        xmom = xmom / 1.0E+20
        x(1,1) = x(1,1) / 1.0E+20
        x(1,2) = x(1,2) / 1.0E+20
        x(1,3) = x(1,3) / 1.0E+20
        x(2,3) = x(2,3) / 1.0E+20
        x(3,3) = x(3,3) / 1.0E+20
        x(2,2) = x(2,2) / 1.0E+20
        return
        end
        
        subroutine usage()
        parameter (LER=0, LIN=5, LOT=6)
        write(LER,*)'fmech96 -D Dip -S Stk -R Rake -M0 Mom -E'
        write(LER,*)' -MW mw  -A Az -B Baz  -ROT'
        write(LER,*)'  -XX Mxx -YY Myy -ZZ Mzz -XY -Mxy -XZ Mxz'
        write(LER,*)'  -YZ Myz  -fx Fx -fy Fy -fz Fz'
        write(LER,*)
     1  ' '
        write(LER,*)
     1  ' -D dip               dip of fault plane'
        write(LER,*)
     1  ' -S Strike            strike of fault plane'
        write(LER,*)
     1  ' -R Rake              slip angle on fault plane'
        write(LER,*)
     1  ' -M0 Moment (def=1.0) Seismic moment in units of dyne-cm'
        write(LER,*)
     1  ' -MW mw               Moment Magnitude  '
        write(LER,*)
     1  ' -E                   Explosion'
        write(LER,*)
     1  ' -A Az                Source to Station Azimuth'
        write(LER,*)
     1  ' -B Baz               Station to Source azimuth'
C       write(LER,*)
C     1 ' -V                   Verbose output'
        write(LER,*)
     1  ' -ROT      (default ZNE) generate Z, R, T 3 component series'
        write(LER,*)
     1  ' -fx FX -fy Fy -fZ fz  Point force amplitudes ',
     2  ' (N,E,down) in  dynes'
        write(LER,*)
     1  ' -XX Mxx -YY Myy -ZZ Mzz  Moment tensor elements in units of'
        write(LER,*)
     1  ' -XY Mxy -XZ Mxz -YZ Myz    dyne-cm'
        write(LER,*)
     1  ' -?                   This online help'
        write(LER,*)
     1  ' -h                   This online help'
        stop
        end

        subroutine chtofp(str,fout)
c------
c       routine to convert string to floating point
c       The E format is accepted as well as free form
c       input
c
c       This assumes that the string is 20 characters long
c-----
        character*(*) str
        real*4 fout
        integer*4 lgstr
        logical hase
        l = lgstr(str)
c------
c       If the string str contains an E or e, then
c       we parse the result using an E format.
c------
        hase = .false.
        do 100 i=1,l
            if(str(i:i) .eq.'e' .or. str(i:i).eq.'E')then
                hase = .true.
            endif
  100   continue
c-----
c       read floating point number
c-----
        if(hase)then
            read(str,'(bn,e20.0)')fout
        else
            read(str,'(bn,f20.0)') fout
        endif
        return
        end
