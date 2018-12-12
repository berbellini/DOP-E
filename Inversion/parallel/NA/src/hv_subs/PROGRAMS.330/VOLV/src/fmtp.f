      program fmtp
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME V                                                       c
c                                                                     c
c      PROGRAM: FMTP                                                  c
c                                                                     c
c      COPYRIGHT 2002 R. B. Herrmann                                  c
c                                                                     c
c      Department of Earth and Atmospheric Sciences                   c
c      Saint Louis University                                         c
c      221 North Grand Boulevard                                      c
c      St. Louis, Missouri 63103                                      c
c      U. S. A.                                                       c
c                                                                     c
c---------------------------------------------------------------------c
c       Changes:
c       04 SEP 2002 - improper test and variable definition in tpdss
c-----
c     This program converts from a strike, dip, rake
c     description of an earthquake fault motion to provide
c     the orientations of the pressure and tension axes. As a
c     double-check, the two nodal planes corresponding to the
c     T- and P-axes are determined.
c-----
c     This program also will accept the trend and plunge
c     angles of the tension and pressure axes to determine
c     the nodal planes.
c
c     The output consists of the orientations of the X, Y, Z
c     T and P axes in cartesian and spherical coordinates,
c     the strike, dip and rake angles for each nodal plane.
c
c     30 AUG 2009 - added a temporary flag to output the nodal planes as
c     SDRPLANES  stk1 dip1 rake1 stk2 dip2 rake2 for flag with a grep
c
c-----
        implicit none
        integer LER, LIN, LOT
        parameter (LER=0, LIN=5, LOT=6)
        logical totp
        real dip, rake, strike
        real ptrend, pplunge, ttrend, tplunge
        logical doextra
c-----
c       machine dependent initialization
c-----
        call mchdep()
c----
        call gcmdln(dip, rake, strike, 
     1      ptrend, pplunge, ttrend, tplunge,
     2      totp,doextra)
        if(totp)then
            call trans(dip,strike,rake,ttrend,
     1          tplunge,ptrend,pplunge,0)
            call tpdss(ttrend,tplunge,ptrend,pplunge,doextra)
        else
            call tpdss(ttrend,tplunge,ptrend,pplunge,doextra)
        endif
        end

      subroutine trans(dip,stk,rake,stkt,plnt,stkp,plnp,iverby)
      parameter (LER=0, LIN=5, LOT=6)
c-----
c     This takes strike, dip, and rake angles and defines a coordinate
c     transformation in which the rake vector, fault normal vector
c     and the normal vector of the auxiliary plane in terms of the
c     local NS, EW, and vertical cartesian coordinate system.
c-----
c     dip   angle of dip measured from horizontal plane.
c           0 < dip < 90
c     stk   direction of nodal plane strike. When looking
c           along stirke, fault plane dips downward to the
c           right
c           0 < stk < 360
c     rake  direction of movement along nodal plane.
c           If -180 < rake < 180, then the P-wave
c           vertically downward will a first motion polarity
c           equal to the sign of rake
c-----
c
c     REFERENCE:
c
c     Herrmann, R. B. (1975). A student's guide to the use of
c     P and S wave data, Earthquake Notes 46, 29-39
c
c-----
c     The X, Y and Z axes form a right handed coordinate
c     system, such that the compressional quadrant is
c     occurs whenever the xy > 0. The Z axis is the
c     null axis. The Y-axis is normal to the fault
c     plane and the X-axis is normal to the auxiliary plane
c     Note that for the input convention on dip, the
c     Y-axis will always point in the negative z-direction
c
c-----
   10 format(' ',a1,':','(',f10.2,',',f10.2,',',f10.2,' )')
   11 format(' ',a1,':','(',f10.2,',',f10.2,         ' )')
   12 format(' '/' ',6x,'  X-DIR   ','  Y-DIR   ','  Z-DIR   ',/)
   13 format(' '/' ',6x,'  TREND   ','  PLUNGE  ',/)
      degrad=0.017452329
      sins=sin(stk*degrad)
      coss=cos(stk*degrad)
      sind=sin(dip*degrad)
      cosd=cos(dip*degrad)
      sinf=sin(rake*degrad)
      cosf=cos(rake*degrad)
c-----
c     X-axis
c-----
      a11=cosf*coss+sinf*cosd*sins
      a12=cosf*sins-sinf*cosd*coss
      a13=-sinf*sind
c-----
c     Y-axis
c-----
      a21=-sins*sind
      a22=coss*sind
      a23=-cosd
c-----
c     Z-axis
c-----
      a31=coss*sinf-cosd*cosf*sins
      a32=cosd*cosf*coss+sinf*sins
      a33=cosf*sind
      t1=a11+a21
      t2=a12+a22
      t3=a13+a23
      tnorm=sqrt(2.)
      if(t3.le.0.0) tnorm=-tnorm
      t1=t1/tnorm
      t2=t2/tnorm
      t3=t3/tnorm
      p1=a11-a21
      p2=a12-a22
      p3=a13-a23
      pnorm=sqrt(2.)
      if(p3.le.0.0) pnorm=-pnorm
      p1=p1/pnorm
      p2=p2/pnorm
      p3=p3/pnorm
      xnorm=1.0
      if(a13.lt.0.0) xnorm=-1.0
      ynorm=1.0
      if(a23.lt.0.0) ynorm=-1.0
      znorm=1.0
      if(a33.lt.0.0) znorm=-1.0
      a11=a11/xnorm
      a12=a12/xnorm
      a13=a13/xnorm
      a21=a21/ynorm
      a22=a22/ynorm
      a23=a23/ynorm
      a31=a31/znorm
      a32=a32/znorm
      a33=a33/znorm
c------
c     now all vectors point into the lower hemisphere.
c     To get the correct orientations, we require
c     that if the center of the focal sphere is a compression
c     that the X,Y, Z axes form a right handed coordinate system in the
c     lower hemisphere, otherwise it will form a left handed
c     coordinate system
c-----
c     obtain P-wave polarity at the center of the focal
c     sphere
c-----
      xy = a13*a23
c-----
c     determine if right handed or left handed coordinate system
c-----
      z3 = a11*a22 - a12*a21
c-----
c     make right handed coordinate system
c-----
      if(z3.lt.0.0)then
            tmp1=a11
            tmp2=a12
            tmp3=a13
            a11=a21
            a12=a22
            a13=a23
            a21=tmp1
            a22=tmp2
            a23=tmp3
      endif
      if(sign(1.0,xy).ne.sign(1.0,rake))then
            tmp1=a11
            tmp2=a12
            tmp3=a13
            a11=a21
            a12=a22
            a13=a23
            a21=tmp1
            a22=tmp2
            a23=tmp3
      endif
      call phase(a11,a12,a13,delx,betx)
      call phase(a21,a22,a23,dely,bety)
      call phase(a31,a32,a33,delz,betz)
      call phase(t1,t2,t3,delt,bett)
      call phase(p1,p2,p3,delp,betp)
      plx = 90-betx
      ply = 90-bety
      plz = 90-betz
      plt=90-bett
      plp=90-betp
      stkt=delt
      plnt=90.0-bett
      stkp=delp
      plnp=90.0-betp
      if(iverby.eq.1)then
            write(LOT,12)
            write(LOT,10)'X:',a11,a12,a13
            write(LOT,10)'Y:',a21,a22,a23
            write(LOT,10)'Z:',a31,a32,a33
            write(LOT,10)'T:',t1,t2,t3
            write(LOT,10)'P:',p1,p2,p3
            write(LOT,13)
            write(LOT,11)'X',delx,plx
            write(LOT,11)'Y',dely,ply
            write(LOT,11)'Z',delz,plz
            write(LOT,11)'T',stkt,plnt
            write(LOT,11)'P',stkp,plnp
      endif
      return
      end

      subroutine phase(x,y,z,del,bet)
c-----
c     This expresses a three component vector in terms
c     of two angles in a spherical coordinate system
c-----
c     x     x-coordinate of vector
c     y     y-coordinate of vector
c     z     z-coordinate of vector
c-----
c     del   direction of projection of the vector on a
c           horizontal plane, measured positively, clockwise
c           from north - the trend of the vector.
c     bet   angle that the 3-D vector makes with the
c           downward vertical. This is just 90 - plunge
c-----
      implicit real*4 (a-h,o-z)
      degrad=0.017452329
c     due to numerical inaccuracy we need:
      if(abs(z).le.1.0) then
            bet=acos(z)
      else if(z.ge.1.0) then
            bet = 0.0
c       write(6,*) ' ACOS(Z) = 0; Z = ',z
      else if(z.le.-1.0) then
            bet = 3.14159265
c       write(6,*) ' ACOS(Z) = 0; Z = ',z
      endif
      bet=bet/degrad
      if(x.eq. 0.0 .and. y.eq. 0.0)then
            del = 0.0
      else
            del = atan2(y,x)/degrad
            if(del.gt.360.) del=del-360.
            if(del.lt.0.0) del=del+360.
      endif
      return
      end

      subroutine tpdss(ttrend,tplunge,ptrend,pplunge,doextra)
      parameter (LER=0, LIN=5, LOT=6)
c-----
c     This takes the orientations of the tension and pressure
c     axes, and determines the two nodal planes
c-----
c
c     convert (P,T) to (dip,rake,stk) and (X,Y,Z).
c-----
        real ttrend, tplunge,ptrend,pplunge
        logical doextra
      double precision t1,t2,t3,p1,p2,p3,f1,f2,f3,xn1,xn2,xn3
      double precision deg,stkt,plnt,stkp,plnp
      double precision stk,rake,dip,stkk,rakep,dipp,xx,yy
      double precision stk0,rake0,dip0
      double precision stk0sv, rake0sv,dip0sv
   14 format(/' ',' NODAL PLANES '/)
   15 format(' ',' STK= ',f10.2/' ',' DIP= ',f10.2/' ','RAKE= ',f10.2)
   16 format('SDRPLANES:',3f10.2,5x,3f10.2)
      deg=3.141592653/180.0
        dstkt = ttrend
        dplnt = tplunge
        dstkp = ptrend
        dplnp = pplunge
        
c-----
c     strike angle is measured from the north clockwise.
c     plunge angle is measured from the horizon downward.
c-----
      if(dstkt.lt.-999.0d+00) go to 200
      stkp=dstkp*deg
      plnp=dplnp*deg
      stkt=dstkt*deg
      plnt=dplnt*deg
      t1=dcos(plnt)*dcos(stkt)
      t2=dcos(plnt)*dsin(stkt)
      t3=dsin(plnt)
      p1=dcos(plnp)*dcos(stkp)
      p2=dcos(plnp)*dsin(stkp)
      p3=dsin(plnp)
      f1=t1+p1
      f2=t2+p2
      f3=t3+p3
      yy=dsqrt(f1*f1+f2*f2+f3*f3)
      f1=f1/yy
      f2=f2/yy
      f3=f3/yy
      xn1=t1-p1
      xn2=t2-p2
      xn3=t3-p3
      yy=dsqrt(xn1*xn1+xn2*xn2+xn3*xn3)
      xn1=xn1/yy
      xn2=xn2/yy
      xn3=xn3/yy
      if(xn3.gt.0.0)then
            xn1=-xn1
            xn2=-xn2
            xn3=-xn3
            f1=-f1
            f2=-f2
            f3=-f3
      endif
      dip=dacos(-xn3)
      stk=datan2(-xn1,xn2)
      xx=f1*xn2-f2*xn1
      rake=datan2(-f3,xx)
      dip0=dip/deg
      rake0=rake/deg
      if(rake0.lt.-180.)rake0=rake0+180.
      if(rake0.gt.180.)rake0=rake0-360.
      stk0=stk/deg
      if(stk0.lt.0.0) stk0=360.0+stk0
      write(LOT,14)
      write(LOT,*) ' '
      write(LOT,15)stk0,dip0,rake0
      stk0sv = stk0
      dip0sv = dip0
      rake0sv = rake0
      if(f3.gt.0.0)then
            xn1=-xn1
            xn2=-xn2
            xn3=-xn3
            f1=-f1
            f2=-f2
            f3=-f3
      endif
      dipp=dacos(-f3)
      stkk=datan2(-f1,f2)
      xx=f2*xn1-f1*xn2
      rakep=datan2(-xn3,xx)
      dip0=dipp/deg
      rake0=rakep/deg
      stk0=stkk/deg
      if(rake0.lt.-180.)rake0=rake0+360.
      if(rake0.gt.180.)rake0=rake0-360.
      if(stk0.lt.0.0) stk0=360.0+stk0
      write(LOT,*)' '
      write(LOT,*) '            OR'
      write(LOT,*)' '
      write(LOT,15)stk0,dip0,rake0
      if(doextra)then
        write(LOT,*) ' '
        write(LOT,16)stk0sv,dip0sv,rake0sv,stk0,dip0,rake0
      endif
      call trans(sngl(dip0),sngl(stk0),sngl(rake0),
     1  sngl(stkt),sngl(plnt),sngl(stkp),sngl(plnp),1)
200   continue
      return
      end

        subroutine gcmdln(dip, rake, strike, 
     1      ptrend, pplunge, ttrend, tplunge,
     2      totp,doextra)
c-----
c       parse command line arguments and return control
c       parameters
c
c       totp    L   - convert from strike dip rake to t P axes
c       dip R   - fault plane dip
c       rake    R   - fault plane rake
c       strike  R   - fault plane strike
c       ptrend  R   - trend of P axis
c       pplunge R   - plunge of P axis - horiz - 0 vertical doen = +90
c       ttrend  R   - trend of T axis
c       tplunge R   - plunge of T axis - horiz - 0 vertical doen = +90
c       doextra L   - output the line
c              SDRPLANES  stk1 dip1 rake1 stk2 dip2 rake2 (default false)
c-----

        implicit none
        real dip, rake, strike
        real ptrend, pplunge, ttrend, tplunge
        logical totp
        logical doextra

        integer i
        integer nmarg
        logical doit

        character*10 names
        integer*4 mnmarg
        nmarg = mnmarg()
        totp = .false.
        doit = .false.
        dip = 0
        rake =0
        strike = 0
        ptrend = 0
        pplunge = 0
        ttrend = 0
        tplunge = 90
         doextra = .false.

        i = 0
   11   i = i + 1
            if(i.gt.nmarg)goto 13
            call mgtarg(i,names)
            if(names(1:2).eq.'-D')then
                i = i + 1
                call mgtarg(i,names)
                read(names,'(bn,f10.0)')dip
                totp = .true.   
                doit = .true.
            else if(names(1:2).eq.'-R')then
                i = i + 1
                call mgtarg(i,names)
                read(names,'(bn,f10.0)')rake
                totp = .true.   
                doit = .true.
            else if(names(1:2).eq.'-S')then
                i = i + 1
                call mgtarg(i,names)
                read(names,'(bn,f10.0)')strike
                totp = .true.   
                doit = .true.
            else if(names(1:4).eq.'-PTR')then
                i = i + 1
                call mgtarg(i,names)
                read(names,'(bn,f10.0)')ptrend
                totp = .false.  
                doit = .true.
            else if(names(1:4).eq.'-PPL')then
                i = i + 1
                call mgtarg(i,names)
                read(names,'(bn,f10.0)')pplunge
                totp = .false.  
                doit = .true.
            else if(names(1:4).eq.'-TTR')then
                i = i + 1
                call mgtarg(i,names)
                read(names,'(bn,f10.0)')ttrend
                totp = .false.  
                doit = .true.
            else if(names(1:4).eq.'-TPL')then
                i = i + 1
                call mgtarg(i,names)
                read(names,'(bn,f10.0)')tplunge
                totp = .false.  
                doit = .true.
            else if(names(1:2).eq.'-X')then
                doextra = .true.
            else if(names(1:2).eq.'-?')then
                call usage(' ')
            else if(names(1:2).eq.'-h')then
                call usage(' ')
            endif
        go to 11
   13   continue
        if(.not. doit)then
            call usage('No command line arguments specified')
        endif
        return
        end

        subroutine usage(str)
        implicit none   
        integer LER, LIN, LOT
        parameter (LER=0, LIN=5, LOT=6)
        character str*(*)
        if(str.ne.' ' )then
            write(LER,*)'fmtp:',str
        endif
        write(LER,*)'USAGE: ',
     1  'fmtp  -S stk -R rake -D dip', '  OR '
        write(LER,*)'       ',
     1  'fmtp  -PTR ptrend -PPP pplunge -TTR ttrend -TPP tplunge', 
     2  '  OR '
        write(LER,*)'       ',
     1  'fmtp  -? -h'
        write(LER,*)
     1  ' -D dip               dip of fault plane'
        write(LER,*)
     1  ' -S Strike            strike of fault plane'
        write(LER,*)
     1  ' -R Rake              rake angle on fault plane'
        write(LER,*)
     1  ' -PTR ptrend          P-axis trend, 0 = N'
        write(LER,*)
     1  ' -PPL pplunge         P-axis plunge, 90 = down'
        write(LER,*)
     1  ' -TTR ptrend          T-axis trend, 0 = N'
        write(LER,*)
     1  ' -TPL pplunge         T-axis plunge, 90 = down'
        write(LER,*)
     1  ' -X     (default false) output extra line SDRPLANES sdr sdr'
        write(LER,*)
     1  ' -?           Write this help message'
        write(LER,*)
     1  ' -h           Write this help message'
        stop
        end
