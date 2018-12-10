        

      subroutine  mtdec(mtensor,m0,ounit,stk0,dip0,rak0,
     1  stk1, dip1, rak1,xmw,pcvld) 
c---------------------------------------------------------------------c
c                                                                   c
c    COMPUTER PROGRAMS IN SEISMOLOGY                                c
c    VOLUME VIII                                                    c
c                                                                   c
c    PROGRAM: MTDEC                                                 c
c                                                                   c
c    COPYRIGHT 1988 M. L. Jost                                      c
c                                                                   c
c    Department of Earth and Atmospheric Sciences                   c
c    Saint Louis University                                         c
c    221 North Grand Boulevard                                      c
c    St. Louis, Missouri 63103                                      c
c    U. S. A.                                                       c
c                                                                   c
c---------------------------------------------------------------------c
c     Changes:
c     07 JAN 2003 - replaced Numerical Recipe routines
c     09 JAN 2005 - slip not defined in subroutine trans1
c         removed everything related to X Y and B axes baker@usgs.gov
c     21 JAN 2006 - for special case of a pure double couple, the
c         returned STK DIP RAKE MW were wrong - fixed if statement in
c         subroutine mjrcpl
c     04 OCT 2006 - In an effort to both comply with and promulgate 
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
c-----
      parameter (LER=0, LIN=5, LOT=6)
      implicit real*8 (a-h,o-z)
        integer index(3)
        real*4 pcvld
      real*8 z(3,3),ev(3),evd(3),evp(3)
      real*8 a1a1(3,3),a2a2(3,3),a3a3(3,3)
        real*8 zz(3,3)
        real*4 mtensor(6), m0,stk0, dip0, rak0, xmw
        real*4 stk1, dip1, rak1
        real*8 xmt(3,3)
        integer ounit
      tol = 1.0e-10
      np=3
c-----
c     mtensor(i) =>  i=1(m11);=2(m22);=3(m12);=4(m13);=5(m23);=6(m33)
c-----
        xmt(1,1) = mtensor(1)*m0
        xmt(2,2) = mtensor(2)*m0
        xmt(1,2) = mtensor(3)*m0
        xmt(2,1) = mtensor(3)*m0
        xmt(1,3) = mtensor(4)*m0
        xmt(3,1) = mtensor(4)*m0
        xmt(2,3) = mtensor(5)*m0
        xmt(3,2) = mtensor(5)*m0
        xmt(3,3) = mtensor(6)*m0

c-----
c
c   mtdec -mc -cc -dc -vd -clvd -a
c
c   This program performs a general moment tensor decomposition.
c   First, the isotropic part is taken out of the moment tensor.
c   The null-, T-, and P-axis are determined from
c   the eigenvalues and eigenvectors of a seismic moment tensor.
c   Consequently, the output of MTEIG (VIII) is used as input.
c   The eigenvector to the largest positive eigenvalue is the T
c   axis. The eigenvector to the largest negative eigenvalue is
c   the P-axis. 
c -mc The smallest eigenvalue in the
c   absolute sense is set to zero defining the major double 
c   couple. The corresponding eigenvector defines the null-axis.
c   From these axes, the focal mechanism i.e. strike, dip, and slip
c   as well as the trend and plunge of the X-, Y-, null-, T-, and P-
c   axes is determined following Herrmann (1975, see XYZTP (II)).
c   The minor couple is then constructed by letting the eigenvector
c   of the largest eigenvalue (in the absolute sense) be the null axis
c   (Jost and Herrmann, 1989).
c -cc The deviatoric moment tensor is decomposed into a 
c   double couple and a
c   CLVD (Knopoff and Randall, 1970).
c -dc The moment tensor is decomposed into three double couples 
c   following
c   Ben-Menahem and Singh (1981, eq.4.57). For each double couple,
c   the corresponding focal mechanism is determined (XYZTP(II)).
c -vd The moment tensor is decomposed into three vector 
c   dipoles following
c   Ben-Menahem and Singh (1981, eq.4.55).
c -clvd The moment tensor is decomposed into three 
c   compensated linear vector 
c   dipoles following Ben-Menahem and Singh (1981, eq.4.56).
c -a  all decompositions are considered.
c----- 
c-----
c     decompose the moment tensor
c-----
        ndata = 3
        write(ounit,'(a)') 'MT    INPUT MOMENT TENSOR'
        do 4 i=1,ndata
            write(ounit,55)(xmt(i,j),j=1,3)
    4   continue
   55   format('MT ',3e15.7)
c---- SEISMIC MOMENT OF THE GIVEN MOMENT TENSOR
      xmom = 0.0
        do 11 i=1,ndata
            do 12 j=1,ndata
                xmom=xmom+xmt(i,j)*xmt(i,j)
   12       continue
   11   continue
        xmom = sqrt(xmom/2.0)
    5   format(5e15.7)
        call tred2(np,np,xmt,ev,evd,zz)
        call tql2(np,np,ev,evd,zz,ierr)
c-----
c     use jordan and silver moment iffespdctive of major/minor DC
c-----
        xmom = ev(1)**2 + ev(2)**2 + ev(3)**2
        xmom = sqrt(xmom/2.0)
c-----
c     m0 is the moment scaling factor fro the Green's functions
c     For Computer programs in seismology with earth model in
c     km/sec and gm/cc, m0 == 1.0+20
c-----
c-----
c     ASSUME GREEN's FUNCTIONS ARE IN SAME PHYSICAL UNITS AS
c     OBSERVED - but that Greens are for a Moment of 1.0e+20 dyne cm
c-----
            xmw =  (dlog10(xmom) - 16.10)/1.5


      write(ounit,'(a)') 'MT  '
      write(ounit,'(a)') 'MT   EIGENVALUES                EIGENVECTORS'
      do 6 i=1,ndata
      write(ounit,7)ev(i),(zz(j,i),j=1,3)
    6 continue
    7 format('MT ',e15.7,2x,5f11.7)
c-----
c     output different representations
c-----
        call isotr(ndata,index,tol,thresh,xmom,zmom,ev,evp,
     1      evd,zz,ounit,eps,eps1,eps2)
            pcvld = eps1
        call dyad(a1a1,a2a2,a3a3,zz)
        call dcclvd(ndata,np,index,tol,thresh,xmom,zmom,ev,evp,
     1          evd,zz,ounit,
     1      amom,dcstk0,dcdip0,dcrak0,dcstk1,dcdip1,dcrak1)
            stk0 = dcstk0
            dip0 = dcdip0
            rak0 = dcrak0
            stk1 = dcstk1
            dip1 = dcdip1
            rak1 = dcrak1
      return
      end

      subroutine isotr(ndata,index,tol,thresh,xmom,zmom,ev,evp,evd,z,
     1      ounit,eps,eps1,eps2)
      parameter (LER=0, LIN=5, LOT=6)
      implicit real*8 (a-h,o-z)
        integer ounit
      dimension z(3,3),ev(3),evp(3),evd(3),index(3)
      real tevd(3)
      write(ounit,'(a)')'MT '
      write(ounit,'(a,e15.7)') 
     $'MT   SEISMIC MOMENT OF THE INPUT MOMENT TENSOR: ',xmom
      write(ounit,'(a)') 'MT  ' 
      expl=0.
      do 1 i=1,ndata
      expl=expl+ev(i)
    1 continue
      write(ounit,'(a,e15.7)') 
     $'MT   ISOTROPIC COMPONENT (TRACE OF MOMENT TENSOR): ',expl
      write(ounit,'(a)') 'MT  '
      zmom = expl/3.0
      wmax=0.
      do 2 i=1,ndata
      if(abs(ev(i)).gt.wmax) wmax=abs(ev(i))
    2 continue
      thresh=tol*wmax
      if(expl.gt.thresh) then
      write(ounit,'(a)')
     $'MT   INPUT MOMENT TENSOR IS PARTIALLY DUE TO AN EXPLOSION'
      elseif(expl.lt.-thresh) then
      write(ounit,'(a)')
     $'MT   INPUT MOMENT TENSOR IS PARTIALLY DUE TO AN IMPLOSION'
      else
      write(ounit,'(a)') 
     $'MT   INPUT MOMENT TENSOR IS PURELY DEVIATORIC'
      expl = 0.0
      endif
      do 3 i=1,ndata
      evd(i)=ev(i)-expl/float(ndata)
    3 continue
      write(ounit,'(a)') 
     $'MT   EIGENVALUES OF THE PURELY DEVIATORIC MOMENT TENSOR'
      write(ounit,56) (evd(i),i=1,ndata)
      write(ounit,'(a)') 'MT  '
c---- sorting the eigenvalues (principal stresses) into ascending order by absolute value
      do 4 i=1,ndata
      index(i)=i
      evp(i)=abs(evd(i))
      tevd(i) = evd(i)
    4 continue
      do 8 j=2,ndata
           a=evp(j)
           b=tevd(j)
           nin=index(j)
           do 9 i=j-1,1,-1
                if(evp(i).le.a) goto 10
                evp(i+1)=evp(i)
                tevd(i+1)=tevd(i)
                index(i+1)=index(i)
    9      continue
           i=0
   10      evp(i+1)=a
           tevd(i+1)=b
           index(i+1)=nin
    8 continue
      write(ounit,'(a)') 'MT   SORTED ABS(EIGENVALUES)'
      write(ounit,56) (evp(i),i=1,ndata)
      write(ounit,'(a)') 'MT   SORTED EIGENVALUES'
      write(ounit,56) (tevd(i),i=1,ndata)
c---- Dziewonski, Chou, Woodhouse,  JGR 1981 2825-2852
c---- eps=0 double couple
c---- eps=0.5 vector dipole
c-----
c     in the manner of Dreger and Ford
      eps=evp(1)/evp(ndata)
      eps1=eps*200.0
      eps2=100.0-eps1
      write(ounit,'(a)') 'MT  |EPSILON|    % OF CLVD     % OF DC'
      write(ounit,11) eps, eps1, eps2
c     in the manner of Dreger and Ford
      epsilon2 = - 2.* tevd(1)/abs(tevd(3))
      xk = zmom /(abs(zmom) + abs(tevd(3)) )
      write(ounit,'(a)') 'MT      2 Epsilon      k (Hudson)'
      write(ounit,'(a3,5x,2f10.3)')'MT ',epsilon2, xk
    6 format(5e15.7)
   56 format('MT ',5e15.7)
   11 format('MT ',f10.6,2f12.6)
      return 
      end

      subroutine dcclvd(ndata,np,index,tol,thresh,xmom,zmom,
     1  ev,evp,evd,z,ounit,
     1  amom,dcstk0,dcdip0,dcrak0,dcstk1,dcdip1,dcrak1)
      parameter (LER=0, LIN=5, LOT=6)
      implicit real*8 (a-h,o-z)
        integer ounit
      dimension z(3,3),ev(3),evd(3),evp(3)
      dimension index(3),iev(3)
      jj=0
      do 1 i=1,ndata
            if(abs(evd(i)).lt.thresh) jj=jj+1
    1 continue
      write(ounit,*) ' '
      write(ounit,'(a)') 
     $'  --------------------------------------------------------------'
      write(ounit,*) ' '
      write(ounit,'(a)') 
     $'  DECOMPOSITION INTO A DOUBLE COUPLE AND A CLVD'
      if(jj.eq.1) then
            write(ounit,*) ' '
            write(ounit,'(a)') 
     $      '  DEVIATORIC MOMENT TENSOR IS DUE TO A DOUBLE COUPLE'
            call xyztp2(evd,z,np,ndata,tol,dip0,stk0,rak0,
     1      dip1,stk1,rak1,ounit)
            dcdip0 = dip0
            dcstk0 = stk0
            dcrak0 = rak0
            dcdip1 = dip1
            dcstk1 = stk1
            dcrak1 = rak1
      elseif(jj.eq.2) then
            write(ounit,*) ' '
            write(ounit,'(a)') 
     $      '  DEVIATORIC MOMENT TENSOR IS DUE TO A VECTOR DIPOLE'
      elseif(jj.eq.0) then
            write(ounit,*) ' '
            write(ounit,'(a)') 
     $      '  SEISMIC MOMENTS OF ISOTROPIC PART, DC, CLVD'
            eps=evp(1)/evp(3)
            amj = evp(3) * (1.0 - 2.0 * eps)
            bmj = evp(3) * eps
            amom = amj
            write(ounit,5) zmom,amj,bmj
c---- MAJOR COUPLE
            ev(index(1)) = 0.0
            ev(index(2)) = -evd(index(3))
            ev(index(3)) = evd(index(3))
            write(ounit,*) ' '
            write(ounit,'(a)') ' EIGENVALUES AND EIGENVECTORS OF THE DC'
            do 3 i=1,ndata
                  write(ounit,4) ev(i),(z(j,i),j=1,ndata)
    3       continue
    4 format(e15.7,2x,5f11.7)
    5 format(5e15.7)
            call xyztp2(ev,z,np,ndata,tol,dip0,stk0,rak0,
     1      dip1,stk1,rak1,ounit)
            dcdip0 = dip0
            dcstk0 = stk0
            dcrak0 = rak0
            dcdip1 = dip1
            dcstk1 = stk1
            dcrak1 = rak1
            iev(1) = -1
            iev(2) = -1
            iev(3) = -1
            iev(index(3)) = 2
            write(ounit,*) ' '
            write(ounit,'(a)') '  EIGENVALUE          AXES OF CLVD'
            do 6 i=1,ndata
                  write(ounit,7) iev(i),(z(j,i),j=1,3)
    6       continue
    7 format(i9,2x,5f11.7)
      endif
      return
      end

      subroutine dyad(a1a1,a2a2,a3a3,z)
      implicit real*8 (a-h,o-z)
      dimension z(3,3),a1a1(3,3),a2a2(3,3),a3a3(3,3)
      a1a1(1,1) = z(1,1) * z(1,1)
      a1a1(2,2) = z(2,1) * z(2,1)
      a1a1(3,3) = z(3,1) * z(3,1)
      a1a1(1,2) = z(1,1) * z(2,1)
      a1a1(1,3) = z(1,1) * z(3,1)
      a1a1(2,3) = z(2,1) * z(3,1)
      a1a1(2,1) = a1a1(1,2)
      a1a1(3,1) = a1a1(1,3)
      a1a1(3,2) = a1a1(2,3)
      a2a2(1,1) = z(1,2) * z(1,2)
      a2a2(2,2) = z(2,2) * z(2,2)
      a2a2(3,3) = z(3,2) * z(3,2)
      a2a2(1,2) = z(1,2) * z(2,2)
      a2a2(1,3) = z(1,2) * z(3,2)
      a2a2(2,3) = z(2,2) * z(3,2)
      a2a2(2,1) = a2a2(1,2)
      a2a2(3,1) = a2a2(1,3)
      a2a2(3,2) = a2a2(2,3)
      a3a3(1,1) = z(1,3) * z(1,3)
      a3a3(2,2) = z(2,3) * z(2,3)
      a3a3(3,3) = z(3,3) * z(3,3)
      a3a3(1,2) = z(1,3) * z(2,3)
      a3a3(1,3) = z(1,3) * z(3,3)
      a3a3(2,3) = z(2,3) * z(3,3)
      a3a3(2,1) = a3a3(1,2)
      a3a3(3,1) = a3a3(1,3)
      a3a3(3,2) = a3a3(2,3)
      return
      end

      subroutine xyztp2(ev,xmt,np,nmt,tol,dip0,stk0,rak0,
     1      dip1,stk1,rak1,ounit)
c---- This program gives the nodal planes of the earthquake
c---- rupture together with the orientation of the pressure
c---- and tension axes, as well as the null vector.
c---- The eigenvectors of the moment tensor are converted
c---- to the null-vector and the tension and pressure axes.
      implicit real*8 (a-h,o-z)
        integer ounit
      dimension ev(np),xmt(np,np)
      call trans1(stkt,plnt,stkp,plnp,ev,xmt,np,tol,nmt)
      call tpdss(stkt,plnt,stkp,plnp,dip0,stk0,rak0,
     1      dip1,stk1,rak1,ounit)
      return
      end
 
      subroutine trans1(stkt,plnt,stkp,plnp,ev,xmt,np,tol,nmt)
      parameter (LER=0, LIN=5, LOT=6)
      implicit real*8 (a-h,o-z)
      dimension ev(np),xmt(np,np)
      wmax=0.
      do 1 i=1,nmt
      if(abs(ev(i)).gt.wmax) wmax=abs(ev(i))
    1 continue
      thresh=tol*wmax
      do 2 i=1,nmt
      if(abs(ev(i)).lt.thresh) then
      a31=xmt(1,i)
      a32=xmt(2,i)
      a33=xmt(3,i)
      endif
      if(ev(i).lt.-thresh) then
      p1=xmt(1,i)
      p2=xmt(2,i)
      p3=xmt(3,i)
      endif
      if(ev(i).gt.thresh) then
      t1=xmt(1,i)
      t2=xmt(2,i)
      t3=xmt(3,i)
      endif
    2 continue
      tnorm=1.0
      if(t3.le.0.0) tnorm=-tnorm
      t1=t1/tnorm
      t2=t2/tnorm
      t3=t3/tnorm
      pnorm=1.0
      if(p3.le.0.0) pnorm=-pnorm
      p1=p1/pnorm
      p2=p2/pnorm
      p3=p3/pnorm
      xnorm=sqrt(2.)
      call phase(t1,t2,t3,delt,bett)
      call phase(p1,p2,p3,delp,betp)
      plt=90.-bett
      plp=90.-betp
      stkt=delt
      plnt=90.0-bett
      stkp=delp
      plnp=90.0-betp
      return
      end

      subroutine phase(x,y,z,del,bet)
c-----
c   This expresses a three component vector in terms
c   of two angles in a spherical coordinate system
c-----
c   x     x-coordinate of vector
c   y     y-coordinate of vector
c   z     z-coordinate of vector
c-----
c   del   direction of projection of the vector on a
c         horizontal plane, measured positively, clockwise
c         from north - the trend of the vector.
c   bet   angle that the 3-D vector makes with the
c         downward vertical. This is just 90 - plunge
c-----
      implicit real*8 (a-h,o-z)
      degrad=0.017452329
c   due to numerical inaccuracy we need:
      if(abs(z).le.1.0) then
            bet=acos(z)
      else if(z.ge.1.0) then
            bet = 0.0
c     write(6,*) ' ACOS(Z) = 0; Z = ',z
      else if(z.le.-1.0) then
            bet = 3.14159265
c     write(6,*) ' ACOS(Z) = 0; Z = ',z
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

      subroutine tpdss(dstkt,dplnt,dstkp,dplnp,dip0,stk0,rak0,
     1      dip1,stk1,rak1,ounit)
      parameter (LER=0, LIN=5, LOT=6)
c-----
c   This takes the orientations of the tension and pressure
c   axes, and determines the two nodal planes
c-----
c
c   convert (P,T) to (dip,slip,stk) and (X,Y,Z).
c-----
      implicit real*8 (a-h,o-z)
        integer ounit
      deg=3.141592653/180.0
c-----
c   strike angle is measured from the north clockwise.
c   plunge angle is measured from the horizon downward.
c-----
      if(stkt.lt.-999.) go to 200
      stkp=dstkp*deg
      plnp=dplnp*deg
      stkt=dstkt*deg
      plnt=dplnt*deg
      t1=cos(plnt)*cos(stkt)
      t2=cos(plnt)*sin(stkt)
      t3=sin(plnt)
      p1=cos(plnp)*cos(stkp)
      p2=cos(plnp)*sin(stkp)
      p3=sin(plnp)
      f1=t1+p1
      f2=t2+p2
      f3=t3+p3
      yy=sqrt(f1*f1+f2*f2+f3*f3)
      f1=f1/yy
      f2=f2/yy
      f3=f3/yy
      xn1=t1-p1
      xn2=t2-p2
      xn3=t3-p3
      yy=sqrt(xn1*xn1+xn2*xn2+xn3*xn3)
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
c   due to numerical inaccuracy we need:
      if(abs(xn3).le.1.0) then
        dip=acos(-xn3)
        elseif(xn3.ge.1.0) then
        dip = 3.14159265
c     write(6,*) ' DIP(xn3) = 0; xn3 = ',xn3
        elseif(xn3.le.-1.0) then
        dip = 0.0
c     write(6,*) ' DIP(xn3) = 0; xn3 = ',xn3
        endif
      stk=atan2(-xn1,xn2)
      xx=f1*xn2-f2*xn1
      slip=atan2(-f3,xx)
      dip0=dip/deg
      rak0=slip/deg
      if(rak0.lt.-180.)rak0=rak0+180.
      if(rak0.gt.180.)rak0=rak0-360.
      stk0=stk/deg
      if(stk0.lt.0.0) stk0=360.0+stk0
      if(stk0.lt.0.001) stk0=0.
      if(abs(dip0).lt.0.01) dip0=0.
      if(abs(rak0).lt.0.01) rak0=0.
      write(ounit,14)
   14 format(/' ',' NODAL PLANES ')
      write(ounit,'(a,f10.5)') '   STK= ',stk0
      write(ounit,'(a,f10.5)') '   DIP= ',dip0
      write(ounit,'(a,f10.5)') '  SLIP= ',rak0
      if(f3.gt.0.0)then
            xn1=-xn1
            xn2=-xn2
            xn3=-xn3
            f1=-f1
            f2=-f2
            f3=-f3
      endif
c   due to numerical inaccuracy we need:
      if(abs(f3).le.1.0) then
        dipp=acos(-f3)
        elseif(f3.ge.1.0) then
        dipp = 3.14159265
c     write(6,*) ' DIPP(f3) = 0; f3 = ',f3
        elseif(f3.le.-1.0) then
        dipp = 0.0
c     write(6,*) ' DIPP(f3) = 0; f3 = ',f3
        endif
      stkk=atan2(-f1,f2)
      xx=f2*xn1-f1*xn2
      slipp=atan2(-xn3,xx)
      dip1=dipp/deg
      rak1=slipp/deg
      stk1=stkk/deg
      if(rak1.lt.-180.)rak1=rak1+360.
      if(rak1.gt.180.)rak1=rak1-360.
      if(stk1.lt.0.0) stk1=360.0+stk1
      if(stk1.lt.0.001) stk1=0.
      if(abs(dip1).lt.0.01) dip1=0.
      if(abs(rak1).lt.0.01) rak1=0.
      write(ounit,*) '  OR'
      write(ounit,'(a,f10.5)') '   STK= ',stk1
      write(ounit,'(a,f10.5)') '   DIP= ',dip1
      write(ounit,'(a,f10.5)') '  SLIP= ',rak1
      call trans(dip1,stk1,rak1,stkt,plnt,stkp,plnp,1,ounit)
200   continue
      return
      end

      subroutine trans(dip,stk,slip,stkt,plnt,stkp,plnp,iverby,ounit)
      parameter (LER=0, LIN=5, LOT=6)
      implicit real*8 (a-h,o-z)
        integer ounit
c-----
c   This takes strike, dip, and slip angles and defines a coordinate
c   transformation in which the slip vector, fault normal vector
c   and the normal vector of the auxiliary plane in terms of the
c   local NS, EW, and vertical cartesian coordinate system.
c-----
c   dip   angle of dip measured from horizontal plane.
c         0 < dip < 90
c   stk   direction of nodal plane strike. When looking
c         along stirke, fault plane dips downward to the
c         right
c         0 < stk < 360
c   slip  direction of movement along nodal plane.
c         If -180 < slip < 180, then the P-wave
c         vertically downward will a first motion polarity
c         equal to the sign of slip
c-----
c
c   REFERENCE:
c
c   Herrmann, R. B. (1975). A student's guide to the use of
c   P and S wave data, Earthquake Notes 46, 29-39
c
c-----
c   The X, Y and Z axes form a right handed coordinate
c   system, such that the compressional quadrant is
c   occurs whenever the xy > 0. The Z axis is the
c   null axis. The Y-axis is normal to the fault
c   plane and the X-axis is normal to the auxiliary plane
c   Note that for the input convention on dip, the
c   Y-axis will always point in the negative z-direction
c
c-----
      degrad=0.017452329
      sins=sin(stk*degrad)
      coss=cos(stk*degrad)
      sind=sin(dip*degrad)
      cosd=cos(dip*degrad)
      sinf=sin(slip*degrad)
      cosf=cos(slip*degrad)
c-----
c   X-axis
c-----
      a11=cosf*coss+sinf*cosd*sins
      a12=cosf*sins-sinf*cosd*coss
      a13=-sinf*sind
c-----
c   Y-axis
c-----
      a21=-sins*sind
      a22=coss*sind
      a23=-cosd
c-----
c   Z-axis
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
c   now all vectors point into the lower hemisphere.
c   To get the correct orientations, we require
c   that if the center of the focal sphere is a compression
c   that the X,Y, Z axes form a right handed coordinate system in the
c   lower hemisphere, otherwise it will form a left handed
c   coordinate system
c-----
c   obtain P-wave polarity at the center of the focal
c   sphere
c-----
      xy = a13*a23
c-----
c   determine if right handed or left handed coordinate system
c-----
      z3 = a11*a22 - a12*a21
c-----
c   make right handed coordinate system
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
C      x1234 = 1.0
C      if(sign(x1234,xy).ne.sign(x1234,slip))then
C            tmp1=a11
C            tmp2=a12
C            tmp3=a13
C            a11=a21
C            a12=a22
C            a13=a23
C            a21=tmp1
C            a22=tmp2
C            a23=tmp3
C      endif
      call phase(a11,a12,a13,delx,betx)
      call phase(a21,a22,a23,dely,bety)
      call phase(a31,a32,a33,delz,betz)
      call phase(t1,t2,t3,delt,bett)
      call phase(p1,p2,p3,delp,betp)
      plx = 90.-betx
      ply = 90.-bety
      plz = 90.-betz
      plt=90.-bett
      plp=90.-betp
      stkt=delt
      plnt=90.0-bett
      stkp=delp
      plnp=90.0-betp
      if(iverby.eq.1)then
      if(abs(a11).lt.0.01) a11=0.
      if(abs(a12).lt.0.01) a12=0.
      if(abs(a13).lt.0.01) a13=0.
      if(abs(a21).lt.0.01) a21=0.
      if(abs(a22).lt.0.01) a22=0.
      if(abs(a23).lt.0.01) a23=0.
      if(abs(a31).lt.0.01) a31=0.
      if(abs(a32).lt.0.01) a32=0.
      if(abs(a33).lt.0.01) a33=0.
      if(abs(t1).lt.0.01) t1=0.
      if(abs(t2).lt.0.01) t2=0.
      if(abs(t3).lt.0.01) t3=0.
      if(abs(p1).lt.0.01) p1=0.
      if(abs(p2).lt.0.01) p2=0.
      if(abs(p3).lt.0.01) p3=0.
            write(ounit,*) ' '
      write(ounit,*) '           X-DIR          Y-DIR          Z-DIR'
            write(ounit,'(a,3f15.7)') '  X: ',a11,a12,a13
            write(ounit,'(a,3f15.7)') '  Y: ',a21,a22,a23
            write(ounit,'(a,3f15.7)') '  Z: ',a31,a32,a33
            write(ounit,'(a,3f15.7)') '  T: ',t1,t2,t3
            write(ounit,'(a,3f15.7)') '  P: ',p1,p2,p3
            write(ounit,*) ' '
            write(ounit,*) '       TREND       PLUNGE'
      if(abs(delx).lt.0.01) delx=0.
      if(abs(plx).lt.0.01) plx=0.
      if(abs(dely).lt.0.01) dely=0.
      if(abs(ply).lt.0.01) ply=0.
      if(abs(delz).lt.0.01) delz=0.
      if(abs(plz).lt.0.01) plz=0.
      if(abs(stkt).lt.0.01) stkt=0.
      if(abs(plnt).lt.0.01) plnt=0.
      if(abs(stkp).lt.0.01) stkp=0.
      if(abs(plnp).lt.0.01) plnp=0.
            write(ounit,'(a,3f12.5)') '  X ',delx,plx
            write(ounit,'(a,3f12.5)') '  Y ',dely,ply
            write(ounit,'(a,3f12.5)') '  Z ',delz,plz
            write(ounit,'(a,3f12.5)') '  T ',stkt,plnt
            write(ounit,'(a,3f12.5)') '  P ',stkp,plnp
      endif
      return
      end

      subroutine mteig(ndata,np,ev,xmt,zz,ounit)
      parameter (LER=0, LIN=5, LOT=6)
      implicit real*8 (a-h,o-z)
        integer ounit
      real*8 xmt(np,np),ev(np),ev1(3)
        real*8 zz(3,3)
c
c   This program determines the eigenvalues and eigenvectors of
c   a seismic moment tensor. The program uses the Householder 
c   transformation with further QL decomposition (Press et.al. 1987).
c
      do 3 i=1,ndata
      write(ounit,4) (xmt(i,j),j=1,ndata)
    3 continue
    4 format(5e15.7)
      write(ounit,*) ' '
        call tred2(np,np,xmt,ev,ev1,zz)
        call tql2(np,np,ev,ev1,zz,ierr)
      call out3v(ev,zz,np,ndata,ounit)
      return
      end
  
      subroutine out3v(ev,z,np,ndata,ounit)
      parameter (LER=0, LIN=5, LOT=6)
      implicit real*8 (a-h,o-z)
        integer ounit
      dimension z(np,np),ev(np)
      write(ounit,*) ' '
      write(ounit,'(a)') '  EIGENVALUES                EIGENVECTORS'
      do 1 j=1,ndata
      write(ounit,2) ev(j),(z(i,j),i=1,ndata)
    1 continue
    2 format(e15.7,2x,5f11.7)
      return
      end
      
