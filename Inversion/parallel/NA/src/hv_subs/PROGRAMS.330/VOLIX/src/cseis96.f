c   program seis81
c
c   **************
c
c   Program seis81 is designed for a two point ray tracing and
c   ray synthetic seismogram computation in a 2-D laterally
c   inhomogeneous media with curved interfaces and block structures.
c
c   ****************************************************************
c
c
c-----
c     ARRAY PARAMETERS
c     NRYSEG Number of layer segments in a given ray description.
c         Used in manually input ray description
c     NDIST  Number of distances for which synthetic is generated
c     NBDLYR Number of boundaries in model
c     NNODE  Maximum number of nodes defining each boundary
c     KRAYSG Number of elements of aan individual ray
c         used in the ay(10,KRAYSG) dimension
c         for each layer there will be more than one such
c         because of the ray is computed each dt seconds or less
c-----
        integer NRYSEG
        parameter (NRYSEG=100)
        integer NDIST
        parameter (NDIST=100)
        integer NBDLYR
        parameter (NBDLYR=15)
        integer NNODE
        parameter (NNODE=30)
        integer KRAYSG
        parameter (KRAYSG=1200)
        
      integer code
      real vel(6),y(3)
      common/auxi/intr,int1,iout,irefr,lay,itype,nder,iprint,mprint,ntr
      common/auxx/mmx(20),mmy(20),mmxy(20),maux
      common/cod/code(NRYSEG),kref,kc
      common/den/rho1(NBDLYR),rho2(NBDLYR),nrho
      common/pq/qp1(NBDLYR),qp2(NBDLYR),etap1(NBDLYR),etap2(NBDLYR)
      common/sq/qs1(NBDLYR),qs2(NBDLYR),etas1(NBDLYR),etas2(NBDLYR)
      common/eq/q(KRAYSG), eta(KRAYSG)
      common/eqdoit/idoit
      common/dist/dst(NDIST),ndst,idist,astart,astep,afin,reps
      common/dout/dti(NDIST),npts(NDIST),t0(NDIST),vred(NDIST)
      common/intrf/a1(NNODE,NBDLYR),b1(NNODE,NBDLYR),
     1  c1(NNODE,NBDLYR),d1(NNODE,NBDLYR),x1(NNODE,NBDLYR),
     2  brd(2),iii(NNODE,NBDLYR),npnt(NBDLYR),nint
      common/medim/v(300),sx(350),sy(350),nx(20),ny(20),nvs(NBDLYR),
     1ptos(NBDLYR)
      common/ray/ay(10,KRAYSG),ds(9,NRYSEG),kint(NRYSEG),
     1  ros,mreg,n,iref,ind,ind1,iy(KRAYSG),ips(KRAYSG)
      common/src/azm,fdip,rake
      common/vcoef/a02(300),a20(300),a22(300)
c-----
c     ay(i,j) = i = value, j = index for the point on the ray
c         ay(1,i)...the time corresponding to the point.
c         ay(2,i)...x-coordinate of the point.
c         ay(3,i)...z-coordinate of the point.
c         ay(4,i)...the direction of the ray, specified by an angle
c                   fi. the angle is determined in the same way as
c                   fi0. See input data card no.12.
c         ay(5,i)...velocity v at the point.
c         ay(6,i)...dv/dx
c         ay(7,i)...dv/dz.
c         ay(8,1)...second derivative of v with respect to x.
c         ay(9,i)...second derivative of v with respect to x and z.
c         ay(10,i)..second derivative of v with respect to z.
c         iy(i).....layer corresponding to this spatial point
c-----
        parameter (LER=0, LIN=5, LOT=6)
        logical ext
        logical verbose, ldoray
        character mname*80
c-----
c     call machine dependent initialization
c-----
        call mchdep()
c-----
c     parse command line arguments
c-----
        call gcmdln(verbose,ldoray)
c-----
c     open cseis96.dat
c-----
        inquire(file='cseis96.dat',exist=ext)
        if(.not. ext)then
            write(LER,*)'Ray description file cseis96.dat not located'
            go to 9000
        endif
        open(4,file='cseis96.dat',access='sequential',form='formatted',
     1      status='unknown')
        rewind 4
      irun=0
        LU1 = 1
        read(4,'(a)')mname
        lm = lgstr(mname)
        write(LOT,*)'Model: ',mname(1:lm)
      read(4,*)mprint
      write(LOT,110)mprint
      open (LU1,file='cseis96.trc',access='sequential',
     1      form='unformatted',status='unknown')
c-----
c   NEW ADDITION RBH 26 AUG 84
c-----
      open(8,file='cseis96.amp',access='sequential',form='formatted',
     1            status='unknown')
      rewind 8
      if(irun.eq.0.and.LU1.ne.0)rewind LU1
      irun=1
c
c   specification of the model
c
      call model(xmin,xmax,zmin,zmax)
c
c   specification of the synthetic seismograms
c
      naux=0
        icont = 1
        read(4,*)mep,mout,mdim,method,
     1      mreg,itmax
      if(naux.eq.0)write(LOT,101)mep,mout,mdim,method,mreg,itmax
      if(nrho.eq.1)mreg=1
      naux=1
      if(itmax.eq.0)itmax=20
      irec=1
      iup=0
      ius=0
      mreg=0
      ndr=7
      nder=0
      if(mdim.ne.0)nder=1
      if(mdim.ne.0)ndr=10
c
      if(mep.gt.0)go to 3
      ndst=-mep
        do 1234 i=1,ndst
            read(4,*)dst(i),dti(i),npts(i),t0(i),vred(i)
            write(LOT,103)dst(i),dti(i),npts(i),t0(i),vred(i)
 1234   continue
        if(ndst.gt.1)then
            rstep=(dst(ndst)-dst(1))/float(ndst-1)
        else
            rstep = 0.0
        endif
      go to 4
    3 ndst=mep
      read(4,*)rmin,rstep
      write(LOT,104)rmin,rstep
      a=rmin
      DO 13 I=1,mep
      dst(i)=a
   13 a=a+rstep
c
    4 read(4,*)xsour,zsour,tsour,reps
      write(LOT,104)xsour,zsour,tsour,reps
      read(4,*)dt,amin1,astep1,amax1,amin2,astep2,amax2,ac
      write(LOT,104)dt,amin1,astep1,amax1,amin2,astep2,amax2,ac
      if(dt.lt.0.0000001)dt=1.
      if(ac.lt.0.0000001)ac=0.0001
      tmax=10000.
      if(reps.lt..00001)reps=.05
      ind=-1
        write(LOT,*)'xsour,zsour,tsour,amin1,ac',
     1      xsour, zsour, tsour, amin1, ac
      call ray2(xsour,zsour,tsour,amin1,x,z,t,fi,tmax,dt,ac)
      if(ind.eq.10)write(LOT,111)ind
      if(ind.eq.10)go to 99
      lay=ind
      isour=ind
c-----
c     get the density at the source
c-----
      y(1)=xsour
      y(2)=zsour
      y(3)=amin1
        itype = 1
      call veloc(y,vel)
      ros=rho1(ind)+rho2(ind)*vel(1)
c-----
c   ADDITION RBH
c-----
        lm = lgstr(mname)
        write(8,'(a)')mname(1:lm)
        write(8,100)ndst
        write(8,104)xsour,zsour,tsour,rstep
        do 1235 i=1,ndst
            write(8,103)dst(i),dti(i),npts(i),t0(i),vred(i)
 1235   continue
c-----
c   generate file LU1 for plotting of ray diagrams,travel times and
c   amplitudes
c-----
        if(LU1.gt.0)then
            write(LU1)lm
            write(LU1)mname(1:lm)
            write(LU1)icont,xmin,xmax,zmin,zmax,ldoray
            write(LU1)nint,(npnt(i),i=1,nint)
            do 22 iij=1,nint
                ni=npnt(iij)-1
                 write(LU1)(a1(j,iij),b1(j,iij),
     1              c1(j,iij),d1(j,iij),
     1              x1(j,iij),iii(j,iij),j=1,npnt(iij))
   22       continue
            write(LU1)xsour,zsour
            write(LU1)ndst
            write(LU1)(dst(i),i=1,ndst)
        endif
c
c    loop for elementary waves
c
      if(mout.eq.0)write(LOT,105)
        ncode = 0
   20   continue
        read(4,*)kc,kref,(code(i),i=1,kref)
        ncode = ncode + 1
        if(LU1.gt.0 )then
            write(LU1)ncode,kc,kref,(code(i),i=1,kref)
        endif
        if(kref.eq.0)go to 49
        if(mout.ne.0)write(LOT,107)
        write(LOT,106)ncode,kc,kref,(code(i),i=1,kref)
        if(mout.ne.0)write(LOT,108)
c
c
      astart=amin1
      astep=astep1
      afin=amax1
      if(kc.lt.0)astart=amin2
      if(kc.lt.0)astep=astep2
      if(kc.lt.0)afin=amax2
      call tramp(xsour,zsour,tsour,dt,ac,tmax,itmax,mout,
     1  mdim,ncode,LU1,method,ldoray)
      go to 20
c
c   end of loop for elementary waves
c
   49 jj=0
      aa=0.
c-----
c   ADDITION BY RBH
c-----
      write(8,116)jj,jj,aa,aa,aa,aa,aa,aa,jj,aa,aa,aa,aa,aa
c-----
c     terminate the run
c-----
        icont = 0
        if(LU1.ne.0)then
            write(LU1)icont,xmin,xmax,zmin,zmax,ldoray
        endif
        go to 99
c
  100 format(26i3)
  101 format(2x,26i3)
  103 format(2e12.4, i10, 2e12.4)
  104 format(8f10.5)
  105 format(//2x,'listing of wave codes'//2x,'int.code',5x,'e x t e r n
     1 a l   c o d e')
  106 format(1x,12i5/(11x,10i5))
  107 format(2x,'int.code',5x,'e x t e r n a l   c o d e'/
     1  '  NCODE   KC  KREF CODE(i=1,KREF)')
  108 format('  CONDITION   IND IN1 ITER       X         Z',
     2  '         T           pnew        dd II INDX'   )
  110 format(3i2,a68)
  111 format(/2x,'ind=',i5,/)
  116 format(2i5,6e11.4,i5,5e11.4)
c
   99 close (LU1)
        close (4)
 9000   continue
            end
c
c   **************************************************************
c
      subroutine model(xmin,xmax,zmin,zmax)
c
c   approximation of interfaces and velocity distribution in
c   individual layers.
c
c   ****************************************************************
c
c
c   A short description of the routine for the approximation of
c   interfaces and velocity distribution within individual layers
c   model().
c   *************
c
c   The routine model consists of two principal parts:
c   1) In the first part, the input data concerning interfaces
c      and velocity distribution inside individual layers are
c      read in. Interfaces are interpolated by cubic splines
c      (routine splin), the velocity distribution by bicubic
c      splines (routines biap and splin). The coefficients thus
c      obtained are stored in commons /intrf/, /medim/ and
c      /vcoef/. The input data concerning the density distri-
c      bution are stored in common/den/.
c   2) In the second part, certain tables describing the model are
c      optionally printed on the line printer,see the detailed
c      description in the section on output tables.
c   Routines used: biap, splin, veloc. The routine veloc (see its
c   detailed description following) is used only in the
c   second part of the program, for the construction of the printer
c   plot.
c   Note that instead of the system of routines model, biap, and
c   splin, another system of routines based on some other appro-
c   ximation of the data may be used. In this new system it would
c   only be necessary to generate the values of the quantities
c   stored in commons /intrf/, /medim/, /vcoef/, /den/.
c   The quantities stored in common/auxx/and the quantities lay,
c   itype,nder from common/auxi/ are used for internal
c   communication between the routines model, biap, splin, veloc,
c   see the description of commons /auxi/ and /auxx/.
c
c
c   *************************************************************
        integer NBDLYR
        parameter (NBDLYR=15)
        integer NNODE
        parameter (NNODE=30)
      integer ipr(61)
      real x(30),fx(30),aux(12),vel(6),y(3)
      common/medim/v(300),sx(350),sy(350),nx(20),ny(20),nvs(NBDLYR),
     1ptos(NBDLYR)
      common/intrf/a1(NNODE,NBDLYR),b1(NNODE,NBDLYR),
     1  c1(NNODE,NBDLYR),d1(NNODE,NBDLYR),x1(NNODE,NBDLYR),
     2  brd(2),iii(NNODE,NBDLYR),npnt(NBDLYR),nint
      common/den/rho1(NBDLYR),rho2(NBDLYR),nrho
      common/pq/qp1(NBDLYR),qp2(NBDLYR),etap1(NBDLYR),etap2(NBDLYR)
      common/sq/qs1(NBDLYR),qs2(NBDLYR),etas1(NBDLYR),etas2(NBDLYR)
      common/auxi/intr,int1,iout,irefr,lay,itype,nder,iprint,mprint,ntr
      common/auxx/mmx(20),mmy(20),mmxy(20),maux
c
c   input of data - interfaces
c
        read(4,*)nint,(npnt(i),i=1,nint)
        write(6,102)nint,(npnt(i),i=1,nint)
        do 11 i=1,nint
c-----
c         For each interface get the x and Z coordinates
c-----
            nc=npnt(i)
            read(4,*)(x(j),fx(j),iii(j,i),j=1,nc)
            write(6,103)(x(j),fx(j),iii(j,i),j=1,nc)
c
c   determination of coefficients of cubic parabolas
c   approximating interfaces
c
            do 1 j=1,nc
                x1(j,i)=x(j)
                a1(j,i)=fx(j)
    1       continue
            j1=1
            nmin=1
    2       do 3 j=j1,nc
                j2=j
C               if(iii(j,i))4,3,6
                if(iii(j,i) .lt. 0)then
                    go to 4
                else if(iii(j,i) .eq. 0)then
                    go to 3
                else
                    go to 6
                endif
    3       continue
    4       if(nmin.eq.j2)go to 5
            fx(nmin)=a1(nmin,i)
            call splin(x,fx,nmin,j2)
            key=0
            go to 8
    5       if(j2.eq.nc)go to 11
            j1=j2+1
            nmin=j2
            go to 2
    6       if(nmin.eq.j2)go to 7
            fx(nmin)=a1(nmin,i)
            call splin(x,fx,nmin,j2)
            key=1
            go to 8
    7       in=iii(j2,i)
            x1(j2,i)=x1(in,i-1)
            a1(j2,i)=a1(in,i-1)
            b1(j2,i)=b1(in,i-1)
            c1(j2,i)=c1(in,i-1)
            d1(j2,i)=d1(in,i-1)
            if(j2.eq.(nc-1))go to 11
            j1=j2+1
            nmin=j1
            go to 2
    8       if((j2-nmin).eq.1)go to 10
            j3=j2-1
            do 9 j=nmin,j3
                h=x(j+1)-x(j)
                d=(a1(j+1,i)-a1(j,i))/h
                d1(j,i)=(fx(j+1)-fx(j))/(3.*h)
                c1(j,i)=fx(j)
                b1(j,i)=d-.333333*h*(fx(j+1)+2.*fx(j))
    9       continue
C           if(key)5,5,7
            if(key.le.0)then
                go to 5
            else
                go to 7
            endif
   10       d1(nmin,i)=0.
            c1(nmin,i)=0.
            b1(nmin,i)=(a1(j2,i)-a1(nmin,i))/(x(j2)-x(nmin))
C           if(key)5,5,7
            if(key.le.0)then
                go to 5
            else
                go to 7
            endif
   11   continue
c-----
c     END PROCESSING INTERFACES
c-----
      nc=npnt(1)
      brd(1)=x1(1,1)
      brd(2)=x1(nc,1)
c
c   input data - velocity distribution in individual
c   layers
c
        mx2=0
        my2=0
        mxy2=0
        nc=nint-1
        do 13 i=1,nc
            WRITE(6,*)'Velocity Distribution:',i
            mx1=mx2+1
            my1=my2+1
            mxy1=mxy2+1
            read(4,*)mx,my
            write(6,102)mx,my
            nx(i)=mx
            ny(i)=my
            mx2=mx1+mx-1
            my2=my1+my-1
            mxy2=mxy1+mx*my-1
            read(4,*)(sx(j),j=mx1,mx2)
            read(4,*)(sy(j),j=my1,my2)
            write(6,104)(sx(j),j=mx1,mx2)
            write(6,104)(sy(j),j=my1,my2)
            m1=mxy1
            do 12 l=1,mx
                m2=m1+my-1
                read(4,*)(v(j),j=m1,m2)
                write(6,104)(v(j),j=m1,m2)
                m1=m2+1
   12       continue
c
c   determination of coefficients of bicubic parabolas
c   approximating velocity distribution
c
            call biap(mx1,mx,my1,my,mxy1)
        mmx(i)=mx1
        mmy(i)=my1
        mmxy(i)=mxy1
   13 continue
c
c   densities and s velocities
c
        read(4,*)(rho1(i),rho2(i),i=1,nc)
        write(6,104)(rho1(i),rho2(i),i=1,nc)
        read(4,*)(qp1(i),qp2(i),i=1,nc)
        write(6,104)(qp1(i),qp2(i),i=1,nc)
        read(4,*)(etap1(i),etap2(i),i=1,nc)
        write(6,104)(etap1(i),etap2(i),i=1,nc)
        read(4,*)(qs1(i),qs2(i),i=1,nc)
        write(6,104)(qs1(i),qs2(i),i=1,nc)
        read(4,*)(etas1(i),etas2(i),i=1,nc)
        write(6,104)(etas1(i),etas2(i),i=1,nc)
c-----
c     check for default density definition and also get correct Q
c-----
        do 32 irho = 1,nc
            if(rho1(irho).eq.0.0 .and. rho2(irho).eq.0.0)then
                rho1(irho) = 1.7
                rho2(irho) = 0.2
            endif
            nvs(irho) = 0
   32   continue
        read(4,*)(ptos(i),i=1,nc)
        write(6,104)(ptos(i),i=1,nc)
        do 33 irho=1,nc
   33 if(ptos(irho).lt.sqrt(2.))ptos(irho)=1.732
        nrho=0
        if(ptos(1).ge.100.0)then
            rho1(1)=1.0
            rho2(1)=0.0
            nrho = 1
        endif
c-----
c   printer plot of the velocity distribution
c-----
        read(4,*)vmin,vmax,bmin,bmax
        write(6,104)vmin,vmax,bmin,bmax
        if(abs(bmin).lt. 0.00001)bmin=a1(1,1)
        if(abs(bmax).lt. 0.00001)bmax=a1(1,nint)
c-----
c     output for ray plotting program
c-----
        xmin = brd(1)
        xmax = brd(2)
        zmin = bmin
        zmax = bmax
c-----
c   numerical form of interfaces
c-----
      aux1=(brd(2)-brd(1))/25.
        if(mprint.ge.2)write(6,105)nint
        if(mprint.ge.2)write(6,107)
        do 19 i=1,nint
        nc=npnt(i)-1
        if(mprint.ge.2)write(6,112)
        if(mprint.ge.2)write(6,108)i
        do 15 j=1,nc
        if(mprint.ge.2)
     1write(6,109)d1(j,i),c1(j,i),b1(j,i),a1(j,i),x1(j,i),x1(j+1,i),iii(
     2j,i)
   15 continue
        if(mprint.ge.2)write(6,110)
        aux2=brd(1)
        aux3=aux2
        aux(1)=aux2
        L=1
        j=1
   16 aux3=aux3+aux1
        if(aux3.gt.brd(2))go to 17
        L=L+1
        aux(L)=aux3
        if(L.lt.12)go to 16
   17 if(mprint.ge.2)write(6,111)(aux(m),m=1,L)
        k=1
   18 if(aux(k).gt.x1(j+1,i))j=j+1
        if(aux(k).gt.x1(j+1,i))go to 18
        aux4=aux(k)-x1(j,i)
        aux(k)=((d1(j,i)*aux4+c1(j,i))*aux4+b1(j,i))*aux4+a1(j,i)
        k=k+1
        if(k.le.L)go to 18
        if(mprint.ge.2)write(6,111)(aux(m),m=1,L)
        if(mprint.ge.2)write(6,112)
        if((aux3+aux1).gt.brd(2))go to 19
        l=0
        go to 16
   19 continue
c-----
c     numerical form of velocity distribution
c     for printer plot
c-----
c     MPNX = number of x positions
c     MPNY = number of y positions
c-----
        MPNX = 60
        MPNY = 50
        MPNX1 = MPNX + 1
        MPNY1 = MPNY + 1
        lay=1
        itype=1
        y(3)=1.
        vmm=vmax-vmin
        vmmm=vmm/10.
        if(mprint.ge.1)write(6,114)vmin,vmax,vmmm
        dy=(bmax-bmin)/real(MPNY)
        dx=(brd(2)-brd(1))/4.
        y(2)=bmin-dy
        k1=1
        aux(1)=brd(1)
        do 20 i=2,5
   20 aux(i)=aux(i-1)+dx
        if(mprint.ge.1)write(6,113)(aux(i),i=1,5)
        dx=(brd(2)-brd(1))/real(MPNX)
        maux=0
        nder=0
        do 29 k=1,MPNY1
            y(2)=y(2)+dy
            y(1)=brd(1)-dx
            do 28 l=1,MPNX1
                y(1)=y(1)+dx
                if(lay.ge.nint)goto 24
   21 lay1=lay+1
                nc=npnt(lay1)-1
                do 22 i=1,nc
                    j=i
                    if(y(1).lt.x1(i+1,lay1))go to 23
   22           continue
   23           a2=a1(j,lay1)
                b2=b1(j,lay1)
                c2=c1(j,lay1)
                d2=d1(j,lay1)
                x2=x1(j,lay1)
                aux1=y(1)-x2
                zint=((d2*aux1+c2)*aux1+b2)*aux1+a2
                if(y(2).ge.zint)lay=lay+1
                if(lay.ge.nint)go to 27
                if(y(2).gt.zint)go to 21
                if(lay.le.0)go to 27
   24           nc=npnt(lay)-1
                do 25 i=1,nc
                    j=i
                    if(y(1).lt.x1(i+1,lay))go to 26
   25           continue
   26           a2=a1(j,lay)
                b2=b1(j,lay)
                c2=c1(j,lay)
                d2=d1(j,lay)
                x2=x1(j,lay)
                aux1=y(1)-x2
                zint=((d2*aux1+c2)*aux1+b2)*aux1+a2
                if(y(2).lt.zint)lay=lay-1
                if(lay.le.0)go to 27
                if(y(2).lt.zint)go to 24
   27           if(lay.le.0.or.lay.ge.nint)ipr(l)=9
                if(lay.le.0.or.lay.ge.nint)go to 28
                call veloc(y,vel)
                aux1=10.*(vel(1)-vmin)/vmm
                ipr(l)=ifix(aux1)
                if(ipr(l).lt.0)ipr(l) = 0
                if(ipr(l).gt.9)ipr(l) = 9
        itype=1
   28 continue
        if(k1.eq.k.and.mprint.ge.1)write(6,115)y(2),(ipr(jj),jj=1,MPNX1)
        if(k1.ne.k.and.mprint.ge.1)write(6,116)(ipr(jj),jj=1,MPNX1)
        if(k1.eq.k)k1=k1+10
   29 continue
c
  102 format(16i5)
  103 format(3(2f10.4,i5))
  104 format(8f10.4)
  105 format(/1x,'m o d e l   o f   m e d i u m'/2x,'number of interf
     1aces - ',i3)
  107 format(/1x,'i n t e r f a c e s'/1x,'interfaces are approx
     1imated by cubic parabolas  z=d*(x-x1)**3+c*(x-x1)**2+b*(x-x1)+a  b
     2etween  x1  and  x2'/1x,'coefficients of parabolas are determined
     3 by cubic spline interpolation')
  108 format(/1x,'coefficients of parabolas approximating interface',i3/
     1/15x,'d',14x,'c',14x,'b',14x,'a',2x,'in interval',5x,'from x1',
     27x,'to x2',5x,'index')
  109 format(1x,4e15.5,15x,f10.5,f12.5,i10)
  110 format(1x,'numerical form of interface'/1x,'upper row - x-coor
     1dinates of points of interface, lower row - corresponding z-coordi
     2nates of points of interface'/)
  111 format(1x,13f9.4)
  112 format(/)
  113 format(8x,5(f11.3,1x))
  114 format(1x,'velocity distribution'/1x,'isolines constructed from',
     1  f10.4,'  to ',f10.4,'  with increment',f10.4/)
  115 format(f7.2,1x,61i1)
  116 format(8x,61i1)
c
      return
      end

      subroutine biap(mx1,mx,my1,my,mxy1)
c
c
c   this routine performs the determination of coefficients of
c   the bicubic spline interpolation
c
        integer NBDLYR
        parameter (NBDLYR=15)
      real x(30),fx(30)
      common/medim/v(300),sx(350),sy(350),nx(20),ny(20),nvs(NBDLYR),
     1ptos(NBDLYR)
      common/vcoef/a02(300),a20(300),a22(300)
c
        do 1 j=1,mx
        l=mx1+j-1
    1   x(j)=sx(l)
        do 3 i=1,my
            do 2 j=1,mx
            k=mxy1+(j-1)*my+i-1
    2       fx(j)=v(k)
            call splin(x,fx,1,mx)
            do 3 j=1,mx
                k=mxy1+(j-1)*my+i-1
    3   a20(k)=fx(j)
c
        do 4 i=1,my
            l=my1+i-1
    4   x(i)=sy(l)
        do 6 j=1,mx
            do 5 i=1,my
                k=mxy1+(j-1)*my+i-1
    5   fx(i)=v(k)
            call splin(x,fx,1,my)
            do 6 i=1,my
            k=mxy1+(j-1)*my+i-1
    6   a02(k)=fx(i)
c
        do 7 j=1,mx
        l=mx1+j-1
    7   x(j)=sx(l)
        do 9 i=1,my
        do 8 j=1,mx
        k=mxy1+(j-1)*my+i-1
    8   fx(j)=a02(k)
        call splin(x,fx,1,mx)
        do 9 j=1,mx
        k=mxy1+(j-1)*my+i-1
    9   a22(k)=fx(j)
c
        return
        end
c
c   *************************************************************
c
      subroutine splin(x,fx,nmin,nmax)
c
c   cubic spline interpolation with zero curvatures at
c   the end points
c
      real a(30),b(30),h(30),f(30),x(30),fx(30)
c
      if((nmax-nmin).eq.1)go to 4
      nmin1=nmin+1
      nmax1=nmax-1
      h(nmin1)=x(nmin1)-x(nmin)
      d2=(fx(nmin1)-fx(nmin))/h(nmin1)
      do 1 i=nmin1,nmax1
      h(i+1)=x(i+1)-x(i)
      d1=d2
      d2=(fx(i+1)-fx(i))/h(i+1)
      b(i)=h(i)+h(i+1)
      fx(i)=3.*(d2-d1)/b(i)
      a(i)=h(i)/b(i)
    1 b(i)=h(i+1)/b(i)
    4 fx(nmin)=0.
      fx(nmax)=0.
      if((nmax-nmin).eq.1)return
      h(nmin)=0.
      f(nmin)=0.
      do 2 i=nmin1,nmax1
      xpom=2.+a(i)*h(i-1)
      h(i)=-b(i)/xpom
    2 f(i)=(fx(i)-a(i)*f(i-1))/xpom
      do 3 i=nmin,nmax1
      j=nmax1-(i-nmin)
    3 fx(j)=h(j)*fx(j+1)+f(j)
      return
      end
c
c    s u b r o u t i n e   r a y 2
c
c    *****************************************************************
c
c   *****************************************************************
c
c
c   A short description of the universal 2-D ray tracing routine
c   ray2(x0,z0,t0,fi0,x,z,t,fi,tmax,dt,ac)
c   **************************************
c
c   The routine ray2 computes the ray specified by initial condi-
c   tions x0,z0,t0,fi0. It returns the quantities x,z,t,fi at the
c   point of termination of the ray. The computation is con-
c   trolled by the quantities tmax,dt,ac.
c   The meaning of these parameters is as follows:
c
c         x0,z0...  coordinates of the initial point of the ray.
c         t0...     initial time.
c         fi0...    initial angle.
c         x,z...    coordinates of the point at which the computation
c                   of the ray is terminated.
c         t...      corresponding time.
c         fi...     corresponding angle of the ray.
c         tmax,dt,ac... the same meaning as described in the input
c                   data card no.9.(ray81) card12 (seis81)
c       tmax      seis81 tmax=10000. at most
c
c   The history of the ray is stored in common/ray/, see above. The
c   quantity ind in this common specifies the reason for the termina-
c   tion of the ray.
c   Routines used: rkgs, fct, out, veloc,root.
c   To specify the model, the routine model, with routines biap and
c   splin, has to be used. The information about the model is trans-
c   mitted to the ray tracing package by common/intrf/, common/den/,
c   common/medim/ and common/vcoef/. Routine veloc is used to deter-
c   mine the velocity and its derivatives at any arbitrary point in
c   the medium.
c   The numerical code of the wave is transmitted to the package by
c   common/cod/.
c
c   ****************************************************************
      subroutine ray2(x0,z0,t0,fi0,x,z,t,fi,tmax,dt,ac)
c
c   ray tracing by the runge-kutta's method
c
        integer NRYSEG
        parameter (NRYSEG=100)
        integer NBDLYR
        parameter (NBDLYR=15)
        integer NNODE
        parameter (NNODE=30)
        integer KRAYSG
        parameter (KRAYSG=1200)
      integer code
      real y(3),dery(3),din(3),prm(5),aux(8,3)
      common/intrf/a1(NNODE,NBDLYR),b1(NNODE,NBDLYR),
     1  c1(NNODE,NBDLYR),d1(NNODE,NBDLYR),x1(NNODE,NBDLYR),
     2  brd(2),iii(NNODE,NBDLYR),npnt(NBDLYR),nint
      common/auxi/intr,int1,iout,irefr,lay,itype,nder,iprint,mprint,ntr
      common/cod/code(NRYSEG),kref,kc
      common/ray/ay(10,KRAYSG),ds(9,NRYSEG),kint(NRYSEG),
     1  ros,mreg,n,iref,ind,ind1,iy(KRAYSG),ips(KRAYSG)
      common/eqdoit/idoit
c
c
      y(1)=x0
      y(2)=z0
      y(3)=fi0
      irefr=0
      n=0
      iref=1
      iout=0
      if(ind.gt.0) go to 6
c
c   for ind=-1: determination of the number of the layer, in which
c   the initial point (x0,z0) is situated
c
      if(y(1).lt.brd(1).or.y(1).gt.brd(2)) go to 4
      int1=0
      intr=1
    1 nl=npnt(intr)-1
      do 2 i=1,nl
      j=i
      if(y(1).lt.x1(i+1,intr)) go to 3
    2 continue
    3 a2=a1(j,intr)
      b2=b1(j,intr)
      c2=c1(j,intr)
      d2=d1(j,intr)
      x2=x1(j,intr)
      au =y(1)-x2
      zint=((d2*au +c2)*au +b2)*au +a2
      if(y(2).lt.zint.and.intr.eq.1) go to 4
      if(abs(y(2)-zint).lt..0000001)int1=intr
      if(y(2).lt.zint) goto 5
      isour=intr
      if(abs(y(2)-zint).lt..0000001.and.intr.eq.nint)go to 5
      if(intr.eq.nint)go to 4
      intr=intr+1
      go to 1
    4 ind=10
      return
c
c   determination of the initial conditions for the Runge-Kutta
c   procedure
c
    5 iint1=int1
      if(ind.ge.0) go to 6
      ind=isour
      return
c
    6 lay=isour
      int1=iint1
      if(isour.ne.iabs(code(1)))ind=14
      if(isour.ne.iabs(code(1)))return
      ind=0
      itype=isign(1,code(1))
        idoit = itype
      prm(1)=t0
      prm(2)=tmax
      prm(3)=dt
      prm(4)=ac
      t=prm(1)
      din(1)=.33
      din(2)=.33
      din(3)=.34
    9 do 10 i=1,3
   10 dery(i)=din(i)
c
c   computation of the ray
c
      call rkgs(prm,y,dery,3,ihlf,aux)
      if(ihlf.eq.11)ind=5
      if(ihlf.eq.12)ind=6
      if(ihlf.eq.13)ind=7
      if(ind.ge.5.and.ind.le.7)return
      if(abs(prm(5)).gt..0001) go to 11
c
c   integration from the point of reflection/transmission with
c   a specified regular time step
c
      t=ay(1,n)
      deltt=prm(3)
      if(t.ge.tmax)go to 15
      go to 14
   11 if(abs(prm(5)-1.1).gt..0001) go to 13
c
c   loop for the iterative determination of the point of
c   incidence (not used in this version)
c
      prm(1)=ay(1,n)
      prm(3)=prm(3)/(2.**(ihlf+1))
      do 12 i=1,3
   12 y(i)=ay(i+1,n)
      t=ay(1,n)
      n=n-1
      go to 9
   13 x=y(1)
      z=y(2)
      t=ay(1,n)
      if(abs(prm(5)-2.).gt..0001.and.irefr.eq.1)ind1=-ind1
      if(abs(prm(5)-2.).gt..0001)go to 20
c
c   integration from the point of reflection/transmission to the
c   closest point of the ray corresponding to the regular time
c   step along the ray
c
      prm(1)=ay(1,n)
      i=int((prm(1)-t0)/dt)
      prm(2)=float(i+1)*dt+t0
      go to 9
   14 prm(1)=prm(2)
      prm(2)=tmax
      prm(3)=dt
      n=n-1
      go to 9
   15 ind=12
      x=y(1)
      z=y(2)
      if(irefr.eq.1)ind1=-ind1
c
c
   20 fi=ay(4,n)
c
      return
      end
c
c   ****************************************************************
c
      subroutine rkgs(prmt,y,dery,ndim,ihlf,aux)
      real y(3),dery(3), aux(8,3),a(4),b(4),c(4),prmt(5)
c
c   standard ibm ssp routine to solve a system of ordinary
c   differential equations of the first order by the runge-
c   kutta's method
c
        do 1 i=1,ndim
            aux(8,i)=.0666667*dery(i)
    1   continue
        x=prmt(1)
        xend=prmt(2)
        h=prmt(3)
        prmt(5)=0.
        call fct(x,y,dery)
        tval = h*(xend-x)
        if(tval.lt.0.0)then
            go to 38
        else if(tval.eq.0.0)then
            go to 37
        else
            go to 2
        endif
   2  a(1)=.5
      a(2)=.2928932
      a(3)=1.707107
      a(4)=.1666667
      b(1)=2.
      b(2)=1.
      b(3)=1.
      b(4)=2.
      c(1)=.5
      c(2)=.2928932
      c(3)=1.707107
      c(4)=.5
      do 3 i=1,ndim
      aux(1,i)=y(i)
      aux(2,i)=dery(i)
      aux(3,i)=0.
    3 aux(6,i)=0.0
      irec=0
      h=h+h
      ihlf=-1
      istep=0
      iend=0
 4      continue
        ttval = (x+h-xend)*h
        if(ttval.lt.0.0)then
            go to 7
        else if(ttval.eq.0.0)then
            go to 6
        else 
            go to 5
        endif
    5 h=xend-x
    6 iend=1
    7 call outp(x,y,dery,irec,ndim,prmt)
        if(prmt(5).eq.0.0)then
            go to 8
        else
            go to 40
        endif
    8 itest=0
    9 istep=istep+1
      j=1
   10 continue
            aj=a(j)
            bj=b(j)
            cj=c(j)
            do 11 i=1,ndim
                r1=h*dery(i)
                r2=aj*(r1-bj*aux(6,i))
                y(i)=y(i)+r2
                r2=r2+r2+r2
                aux(6,i)=aux(6,i)+r2-cj*r1
   11       continue
        if(j.ge.4)then
            go to 15
        endif
        j=j+1
        if(j.ne.3)then
            x = x + 0.5*h
        endif
        call fct(x,y,dery)
        go to 10
   15   continue
        if(itest .gt. 0)then
            go to 20
        endif
            do 17 i=1,ndim
                aux(4,i)=y(i)
   17       continue
            itest=1
            istep=istep+istep-2
   18 ihlf=ihlf+1
      x=x-h
      h=.5*h
        do 19 i=1,ndim
            y(i)=aux(1,i)
            dery(i)=aux(2,i)
            aux(6,i)=aux(3,i)
   19   continue
        go to 9
   20   imod=istep/2
        iiii = istep-imod-imod
        if(iiii.ne.0)then
            call fct(x,y,dery)
            do 22 i=1,ndim
                aux(5,i)=y(i)
                aux(7,i)=dery(i)
   22       continue
            go to 9
        endif
        delt=0.
        do 24 i=1,ndim
            delt=delt+aux(8,i)*abs(aux(4,i)-y(i))
   24   continue
        xxxx = delt-prmt(4)
        if(xxxx .le. 0.0)then
            go to 28
        endif
        if(ihlf.ge.10)then
            go to 36
        endif
        do 27 i=1,ndim
            aux(4,i)=aux(5,i)
   27   continue
        istep=istep+istep-4
        x=x-h
        iend=0
        go to 18
   28 call fct(x,y,dery)
        do 29 i=1,ndim
            aux(1,i)=y(i)
        aux(2,i)=dery(i)
        aux(3,i)=aux(6,i)
        y(i)=aux(5,i)
        dery(i)=aux(7,i)
   29   continue
        call outp(x-h,y,dery,ihlf,ndim,prmt)
CFIX      if(prmt(5))40,30,40
        if(prmt(5).ne.0.0)then
            return
        else
            do 31 i=1,ndim
                y(i)=aux(1,i)
                dery(i)=aux(2,i)
   31       continue
            irec=ihlf
            if(iend.le.0)then
                ihlf=ihlf-1
                istep=istep/2
                h=h+h
                if(ihlf .ge. 0)then
                    imod=istep/2
                    iiii=istep-imod-imod
                    if(iiii.eq.0)then
                        xxxx=delt-.02*prmt(4)
                        if(xxxx .le. 0.0)then
                            ihlf=ihlf-1
                            istep=istep/2
                            h=h+h
                        endif
                    endif
                endif
                go to 4
            else
                go to 39
            endif
        endif
CFIX      if(iend)32,32,39
CFIX   32 ihlf=ihlf-1
CFIX      if(ihlf)4,33,33
CFIX   33 imod=istep/2
CFIX      if(istep-imod-imod)4,34,4
CFIX   34 if(delt-.02*prmt(4))35,35,4
CFIX   35 ihlf=ihlf-1
CFIX      istep=istep/2
CFIX      h=h+h
CFIX      go to 4
   36 ihlf=11
      call fct(x,y,dery)
      go to 39
   37 ihlf=12
      go to 39
   38 ihlf=13
   39 call outp(x,y,dery,ihlf,ndim,prmt)
   40 return
      end
c
c   **************************************************************
c
      subroutine fct(x,y,dery)
c
c   computation of the right-hand sides in the ray tracing system
c
      real vel(6),y(3),dery(3)
c
      call veloc(y,vel)
      dery(1)=vel(1)*cos(y(3))
      dery(2)=vel(1)*sin(y(3))
      dery(3)=vel(2)*sin(y(3))-vel(3)*cos(y(3))
      return
      end
c
c
c   ****************************************************************
c
      subroutine outp(x,y,dery,ihlf,ndim,prmt)
c
c   external routine in the runge-kutta's ray tracing.
c   performs various computations at each point of the ray.
c
        integer NRYSEG
        parameter (NRYSEG=100)
        integer NBDLYR
        parameter (NBDLYR=15)
        integer NNODE
        parameter (NNODE=30)
        integer KRAYSG
        parameter (KRAYSG=1200)
      integer code
      real y(3),dery(3),prmt(5),yold(2),yint(2),vel(6),y1(2)
      common/cod/code(NRYSEG),kref,kc
      common/intrf/a1(NNODE,NBDLYR),b1(NNODE,NBDLYR),
     1  c1(NNODE,NBDLYR),d1(NNODE,NBDLYR),x1(NNODE,NBDLYR),
     2  brd(2),iii(NNODE,NBDLYR),npnt(NBDLYR),nint
      common/ray/ay(10,KRAYSG),ds(9,NRYSEG),kint(NRYSEG),
     1  ros,mreg,n,iref,ind,ind1,iy(KRAYSG),ips(KRAYSG)
      common/auxi/intr,int1,iout,irefr,lay,itype,nder,iprint,mprint,ntr
      common/den/rho1(NBDLYR),rho2(NBDLYR),nrho
      common/eqdoit/idoit
c
c
        inspl=0
        n=n+1
        nrr=n
        ntr=1
        if(n.gt.KRAYSG)then
            ind=13
            prmt(5)=3.
            yan=y(3)
            y(3)=101.
            call veloc(y,vel)
            y(3)=yan
            return
        endif
        intr=lay
c
c   check the position of the point of the ray
c
      nl=npnt(intr)-1
      do 2 i=1,nl
      j=i
      if(y(1).lt.x1(i+1,intr)) goto 3
    2 continue
    3 if(iref.eq.1)go to 27
      if(kint(iref-1).eq.(n-1))go to 21
   27 a2=a1(j,intr)
      b2=b1(j,intr)
      c2=c1(j,intr)
      d2=d1(j,intr)
      x2=x1(j,intr)
      aux=y(1)-x2
      zint=((d2*aux+c2)*aux+b2)*aux+a2
      if(y(2).lt.zint) go to 9
      intr=lay+1
      nl=npnt(intr)-1
      do 4 i=1,nl
      j=i
      if(y(1).lt.x1(i+1,intr)) go to 5
    4 continue
    5 a2=a1(j,intr)
      b2=b1(j,intr)
      c2=c1(j,intr)
      d2=d1(j,intr)
      x2=x1(j,intr)
      aux=y(1)-x2
      zint=((d2*aux+c2)*aux+b2)*aux+a2
      if(y(2).gt.zint) go to 9
   21 if(y(1).lt.brd(1).or.y(1).gt.brd(2))go to 7
c
c   the ray did not cross an interface
c
        ay(1,n)=x
        do 6 i=1,3
            ay(i+1,n)=y(i)
    6   continue
        yan=y(3)
        y(3)=101.
        call veloc(y,vel)
        y(3)=yan
      return
c
c   the ray crossed one of the vertical boundaries
c
    7 if(y(1).lt.brd(1)) ind=1
      if(y(1).gt.brd(2)) ind=2
   33 aux=(brd(ind)-ay(2,n-1))/(y(1)-ay(2,n-1))
      ay(1,n)=ay(1,n-1)+aux*(x-ay(1,n-1))
      do 8 i=2,3
    8 ay(i+1,n)=ay(i+1,n-1)+aux*(y(i)-ay(i+1,n-1))
      ay(2,n)=brd(ind)
      x=ay(1,n)
      y(1)=ay(2,n)
      y(2)=ay(3,n)
      y(3)=ay(4,n)
      prmt(5)=3.
      kint(iref)=n
      ntr=2
      go to 32
c
c   the ray crossed an interface.
c
c   determine the point of incidence
c
    9 yold(1)=ay(2,n-1)
      yold(2)=ay(3,n-1)
      irr=iref
      call root(yold,y,yint,anyx,anyz,rc,intr,j,ind,iroot)
      if(ind.eq.1.or.ind.eq.2)go to 33
      ind1=intr
      if(iroot.eq.100)go to 30
c
c   the loop for the iterative determination of the point
c   of incidence (not used in this version)
c
      if(iout.ne.0)go to 11
   10 n=n-1
      y1(1)=yint(1)
      y1(2)=yint(2)
      prmt(5)=1.1
      iout=1
      return
   11 iac=iac+1
      if(iac.eq.10)go to 30
      aux1=y1(1)-yint(1)
      aux2=y1(2)-yint(2)
      aux=aux1*aux1+aux2*aux2
      if(aux.gt..00001) go to 10
c
c   determination of parameters of incident wave at the point
c   of incidence
c
   30 iac=0
      iout=0
      ay(2,n)=yint(1)
      ay(3,n)=yint(2)
      xdif=y(1)-ay(2,n-1)
      ydif=y(2)-ay(3,n-1)
      if(abs(xdif).lt.0.001)go to 60
      aux=(ay(2,n)-ay(2,n-1))/xdif
      go to 62
   60 if(abs(ydif).lt.0.001)go to 61
      aux=(ay(3,n)-ay(3,n-1))/ydif
      go to 62
   61 aux=0.
   62 continue
      ay(1,n)=ay(1,n-1)+aux*(x-ay(1,n-1))
      ay(4,n)=ay(4,n-1)+aux*(y(3)-ay(4,n-1))
      x=ay(1,n)
      do 12 i=1,3
   12 y(i)=ay(i+1,n)
      ntr=3
      if(iref.eq.1.and.kc.eq.1.and.lay.eq.intr)go to 26
      ntr=4
      if(iref.eq.1.and.kc.eq.(-1).and.lay.ne.intr)go to 26
      if(iref.eq.1)go to 49
      nnn=kint(iref-1)
      if(nnn.le.0)go to 49
      if(abs(ay(2,nnn)-yint(1)).gt.0.0001.or.intr.ne.int1)go to 49
      irr=irr-1
        ntr=5
      if(aindex.gt.1.)go to 28
        ntr=6
      go to 26
c
c
   49 itype1=itype
        idoit = 0
      itype=1
      call veloc(y,vel)
      ds(4,iref)=vel(1)
      ds(8,iref)=rho1(lay)+rho2(lay)*vel(1)
      itype=-1
      call veloc(y,vel)
      ds(5,iref)=vel(1)
      itype=itype1
        idoit = itype
      ds(1,iref)=0.
      ds(3,iref)=0.
      ds(6,iref)=0.
      ds(7,iref)=0.
      ds(9,iref)=0.
      if(kref.gt.1)ic1=code(iref)
      lay1=lay
      px=cos(y(3))
      pz=sin(y(3))
      yan=y(3)
      aux=anyx*px+anyz*pz
      tang=-1.57079632
      if(-anyx*anyz.gt.0)tang=-tang
      if(abs(anyz).gt.0.00000001*abs(anyx))tang=atan(-anyx/anyz)
      tang=atan(-anyx/anyz)
      if(aux.le.0.) go to 13
      tang=tang+3.14159265
      anyx=-anyx
      anyz=-anyz
      aux=-aux
      rc=-rc
   13 aa=3.14159265-tang+y(3)
      ds(2,iref)=aa
      if(aa.gt.3.14159265)ds(2,iref)=aa-6.2831852
      if(kref.le.1)go to 31
c
c   multiply reflected wave
c
      ntr=7
      if((iref+1).gt.kref) go to 24
      ncd=code(iref+1)-code(iref)
      ntr=8
      if(intr.eq.int1.and.ncd.ne.0)go to 26
      if(intr.ne.int1.or.ncd.ne.0)go to 46
      irefr=1
      kint(iref)=0
      iref=iref+1
      irr=iref
      ds(2,iref)=ds(2,iref-1)
      ds(4,iref)=ds(4,iref-1)
      ds(5,iref)=ds(5,iref-1)
      ds(8,iref)=ds(8,iref-1)
      do 53 i=1,9
   53   ds(i,iref-1)=0.
      ncd=code(iref+1)-code(iref)
      ntr=9
      if((iref+1).gt.kref) go to 24
   46 int1=intr
      ntr=10
      if(iref.eq.1.and.kc.eq.(-1).and.intr.ne.lay) go to 26
      ntr=11
      if(iref.eq.1.and.kc.eq.1.and.intr.eq.lay) go to 26
      j11=j
      i11=intr
      l11=lay1
      if(ncd.ne.0) go to 22
c
c   reflection of unconverted wave
c
      iiii=iii(j,intr)
      ntr=12
      if(iiii.eq.(-2))go to 48
      y(3)=101.
      call veloc(y,vel)
      y(3)=yan
      aux2=pz-2.*aux*anyz
      aux1=px-2.*aux*anyx
      iref=iref+1
      v1=vel(1)
      v2=v1
      aindex=1.0
      go to 18
c
c   refracted wave
c
   31 ntr=13
      if(intr.eq.lay.and.lay.eq.1) go to 24
      ntr=14
      if(intr.gt.lay.and.intr.eq.nint) go to 24
      ncd=1
      int1=intr
      go to 14
c
c   refraction or reflection of converted wave
c
   22 ncd=iabs(code(iref+1))-iabs(code(iref))
      iiii=iii(j,intr)
      ntr=15
      if(ncd.eq.0.and.iiii.eq.(-2)) go to 48
      ic2=code(iref+1)
   14 y(3)=101.
      call veloc(y,vel)
      y(3)=yan
      v1=vel(1)
   45 iref=iref+1
c
c   check for transmission at an interface which coincides
c   with another interface
c
      if(kref.gt.1)then
            itype=isign(1,code(iref))
            idoit = itype
        endif
      if(ncd.eq.0)go to 16
      if(intr.eq.lay)go to 15
      ntr=16
      if(ncd.lt.0)go to 26
      ntr=17
      if(intr.eq.nint) go to 26
      lay=lay+1
      go to 29
   15 ntr=18
        if(ncd.gt.0.and.kref.gt.1)go to 26
      if(lay.eq.1.and.kref.eq.1)go to 24
      ntr=19
      if(lay.eq.1)go to 26
      lay=lay-1
   29  intra=intr
      ja=j
      if(intr.eq.lay1)go to 36
      nc=npnt(intr+1)
      do 34 i=1,nc
      jj1=i
      ii1=iii(i,intr+1)
      if(j.eq.ii1)go to 35
   34 continue
      go to 16
   35 intr=intr+1
      j=jj1
      go to 37
   36 ii1=iii(j,intr)
      if(ii1.le.0)go to 16
      intr=intr-1
      j=ii1
c
c   transmission at an interface which coincides with
c   another interface
c
   37 n=n+1
      n1=n+1
      ntr=20
      if(n1.gt.KRAYSG) go to 50
      ndr=7
      if(nder.eq.1)ndr=10
      do 38 i=n,n1
      do 38 l=1,ndr
   38 ay(l,i)=ay(l,n-1)
      n=n1
      lay1=lay
      ind1=intr
      kint(iref)=-1
      do 47 i=1,9
   47 ds(i,iref)=0.
      ntr=21
      if(kref.gt.1.and.(iref+1).gt.kref)go to 24
      int1=intr
      if(kref.gt.1)ncd=iabs(code(iref+1))-iabs(code(iref))
      ntr=22
      if(ncd.eq.0.and.iii(j,intr).gt.0)go to 23
      if(kref.gt.1)ic2=code(iref+1)
      go to 45
   16 lay2=lay
      call veloc(y,vel)
      v2=vel(1)
      if(kref.gt.1)ncd=iabs(ic2)-iabs(ic1)
      aindex=v2/v1
      aux1=1.-aindex*aindex*(1.-aux*aux)
      if(aux1.gt.0.) go to 17
      ntr=23
      go to 28
   17 if(ncd.eq.0)aux1=aindex*aux-sqrt(aux1)
      if(ncd.ne.0)aux1=aindex*aux+sqrt(aux1)
      aux2=aindex*pz-aux1*anyz
      aux1=aindex*px-aux1*anyx
      inspl=0
   18 if(abs(aux1).gt..000001) go to 19
      y(3)=1.57079632*sign(1.,pz)
      if(kref.le.1.or.ncd.ne.0) go to 20
      y(3)=-y(3)
      go to 20
   19 y(3)=atan2(aux2,aux1)
c
   20 prmt(5)=2.
      kint(irr)=nrr
      ds(1,irr)=rc
      aux=3.14159265-tang+y(3)
      ds(3,irr)=aux
      if(aux.gt.3.14159265)ds(3,irr)=aux-6.2831852
c
c   determination of parameters of generated waves at a
c   point of incidence
c
      if(ncd.ne.0)go to 52
c
c   check whether the reflecting interface coincides with
c   another interface
c
      ja=j11
      intra=i11
      lay2=l11
      if(intra.eq.l11)go to 43
   40 lay2=lay2+1
      ntr=24
      if(intra.eq.nint.and.nder.ne.0)go to 51
      if(intra.eq.nint)return
      nc=npnt(intra+1)
      do 41 i=1,nc
      jj1=i
      ii1=iii(i,intra+1)
      if(ja.eq.ii1)go to 42
   41 continue
      go to 44
   42 intra=intra+1
      ja=jj1
      go to 40
   43 lay2=lay2-1
      ii1=iii(ja,intra)
      if(ii1.le.0)go to 44
      intra=intra-1
      ja=ii1
      go to 43
   44 if(lay2.ne.0)go to 52
      ds(6,irr)=0.
      ds(7,irr)=0.
      ds(9,irr)=0.
      return
   52 itype1=itype
        idoit = 0
      itype=1
      l11=lay
      lay=lay2
      call veloc(y,vel)
      ds(6,irr)=vel(1)
      ds(9,irr)=rho1(lay2)+rho2(lay2)*vel(1)
      itype=-1
      call veloc(y,vel)
      ds(7,irr)=vel(1)
      itype=itype1
        idoit = itype1
      lay=l11
      return
c
c   termination of computation
   48 ind=16
      go to 25
   51 ind=15
      go to 25
   23 ind=17
      go to 25
   26 ind=8
      go to 25
   24 ind=intr+100
      if(intr.eq.1)ind=3
      if(intr.eq.nint)ind=4
      go to 25
   28 ind=9
   25 prmt(5)=3.
      kint(irr)=nrr
      iref=irr
      if(ind.eq.3.or.ind.eq.4)go to 32
      ds(1,irr)=rc
      aa=3.14159265-tang+y(3)
      ds(2,irr)=aa
      if(aa.gt.3.14159265)ds(2,irr)=aa-6.2831852
      ds(3,irr)=0.
      go to 32
   50 ind=13
      prmt(5)=3.
   32 if(ind.eq.9.or.ind.eq.8)return
      yan=y(3)
      y(3)=101.
      call veloc(y,vel)
      y(3)=yan
      return
      end
c   *************************************************************
c
      subroutine root(yold,ynew,yint,anyx,anyz,rc,intr,j,ind,iroot)
c
c   determination of the point of incidence at an interface
c   (intersection of the ray with the interface).
c
        integer NBDLYR
        parameter (NBDLYR=15)
        integer NNODE
        parameter (NNODE=30)
      real yold(2),ynew(2),yint(2)
      common/intrf/a1(NNODE,NBDLYR),b1(NNODE,NBDLYR),
     1  c1(NNODE,NBDLYR),d1(NNODE,NBDLYR),x1(NNODE,NBDLYR),
     2  brd(2),iii(NNODE,NBDLYR),npnt(NBDLYR),nint
c
      sh(x9)=(exp(x9)-exp(-x9))/2.
      ch(x9)=(exp(x9)+exp(-x9))/2.
      arccos(x9)=atan(sqrt(1.-x9*x9)/x9)
      arcch(x9)=alog(x9+sqrt(x9*x9-1.))
      arcsh(x9)=alog(x9+sqrt(x9*x9+1.))
c
c
      iroot=100
c
c   determination of the element of the interface with the point
c   of incidence
c
      nc=npnt(intr)
      xa=yold(1)
      xb=ynew(1)
      aux=xb-xa
      x2=x1(j,intr)
      a2=a1(j,intr)
      if(abs(aux).gt.0.000001)go to 1
      yint(1)=yold(1)
      go to 2
    1 p=(ynew(2)-yold(2))/aux
      q=yold(2)-p*xa
      if(aux.lt.0.)go to 21
      if(xa.ge.x2)go to 2
      x3=x1(j+1,intr)
      a3=a1(j+1,intr)
      y2=p*x2+q
      y3=p*x3+q
   20 aux1=(y3-a3)*(y2-a2)
      if(aux1.le.0.)go to 2
      j=j-1
      x3=x2
      x2=x1(j,intr)
      a3=a2
      a2=a1(j,intr)
      if(xa.ge.x2)go to 2
      y3=y2
      y2=p*x2+q
      go to 20
   21 x3=x1(j+1,intr)
      if(xa.le.x3)go to 2
      a3=a1(j+1,intr)
      y2=p*x2+q
      y3=p*x3+q
   22 aux1=(y3-a3)*(y2-a2)
      if(aux1.le.0.)go to 2
      j=j+1
      x2=x3
      x3=x1(j+1,intr)
      a2=a3
      if(xa.le.x3)go to 2
      a3=a1(j+1,intr)
      y2=y3
      y3=p*x3+q
      go to 22
c
    2 b2=b1(j,intr)
      c2=c1(j,intr)
      d2=d1(j,intr)
      if(abs(aux).le. 0.0001)go to 17
      xa=x2
      xb=x1(j+1,intr)
      aux=xb-xa
c
c   determination of coefficients of cubic equation d*x**3+c*x**2+b*x+
c   roots are looked for in interval (0,1)
c
      d=d2*aux*aux*aux
      c=c2*aux*aux
      b=(b2-p)*aux
      a=a2-p*xa-q
c
      if(abs(d).lt. 0.000001)go to 10
c
c   transformation of cubic equation into form y**3+3*p*y+q
c   substituting y=x+c/(3*d)
c
      aux1=c/(3.*d)
      q=aux1*aux1*aux1-.5*(b*aux1-a)/d
      p=b/(3.*d)-aux1*aux1
      diskr=q*q+p*p*p
      if(q.eq. 0.0) go to 8
      if(p.eq. 0.0) go to 7
      r=sign(1.,q)*sqrt(abs(p))
      ax=q/(r*r*r)
      if(p.gt. 0.0) go to 6
        if(diskr .le. 0.0)then
            go to 3
        else if (diskr.gt.0.0)then
            go to 5
        endif
CFIX      if(diskr)3,3,5
c
c   p.lt.0.and.diskr.le.0
c
    3 d=arccos(ax)/3.
      xr=-2.*r*cos(d)-aux1
      xr1=2.*r*cos(1.047198-d)-aux1
      xr2=2.*r*cos(1.047198+d)-aux1
   25 nrt=3
CFIX      if(xr.ge.0..and.xr.le.1.)go to 18
      if(xr.lt.0.0 .or. xr.gt.1.0)then
            nrt=nrt-1
            xr=xr2
        endif
CFIX   18 if(xr1.ge.0..and.xr1.le.1.)go to 19
        if(xr1.lt.0.0 .or. xr1.gt.1.0)then
            nrt=nrt-1
            xr1=xr2
        endif
CFIX   19 if(xr2.ge.0..and.xr2.le.1.)go to 4
        if(xr2.lt.0.0 .or. xr2.gt.1.0)then
            nrt=nrt-1
        endif
      if(nrt.eq.1)go to 15
      rr=xa+xr*aux-yold(1)
      if(abs(rr).lt..001)rr=1000000.
      rr1=xa+xr1*aux-yold(1)
      if(abs(rr1).lt..001)rr1=1000000.
      if(abs(rr1).lt.abs(rr))xr=xr1
      if(abs(rr1).lt.abs(rr))rr=rr1
      if(nrt.eq.2)go to 15
      rr1=xa+xr2*aux-yold(1)
      if(abs(rr1).lt..001)rr1=1000000.
      if(abs(rr1).lt.abs(rr))xr=xr2
      go to 15
c
c   p.lt.0..and.diskr.gt.0.
c
    5 xr=-2.*r*ch(arcch(ax)/3.)-aux1
      go to 15
c
c   p.gt.0
c
    6 xr=-2.*r*sh(arcsh(ax)/3.)-aux1
      go to 15
c
c   p.eq.0.
c
    7 xr=-sign(1.,q)*exp(alog(2.*abs(q))/3.)-aux1
      go to 15
c
c   q.eq.0
c
    8 xr=-aux1
        if(p.lt.0.0)then
            xr1=sqrt(-3.*p)-aux1
            xr2=-xr1-2.*aux1
            go to 25
        else
            go to 15
        endif
CFIX      if(p)9,15,15
c
c   q.eq.0..and.p.lt.0.
c
CFIX    9 xr1=sqrt(-3.*p)-aux1
CFIX      xr2=-xr1-2.*aux1
CFIX      go to 25
c
c   reduction of cubic equation to quadratic equation
c
   10 if(abs(c).lt. 0.000001) go to 14
      diskr=b*b-4.*c*a
      p=-b/(2.*c)
CFIX      if(diskr)11,11,12
        if(diskr.le.0.0)then
            xr=p
            go to 15
        else
            go to 12
        endif
c
CFIX   11 xr=p
CFIX      go to 15
c
   12 q=sqrt(diskr)/(2.*c)
      xr=p+q
      xr1=p-q
      nrt=2
CFIX      if(xr.ge.0..and.xr.le.1.)go to 23
        if(xr.lt.0.0 .or. xr.gt.1.0)then
            nrt=nrt-1
            xr=xr1
        endif
CFIX   23 if(xr1.ge.0..and.xr1.le.1.)go to 24
        if(xr1.lt.0.0 .or. xr1.gt.1.0)then
            nrt=nrt-1
        endif
CFIX   24 if(nrt.eq.1)go to 15
        if(nrt.ne.1)then
            rr=xa+xr*aux-yold(1)
            if(abs(rr).lt..001)rr=1000000.
            rr1=xa+xr1*aux-yold(1)
            if(abs(rr1).lt..001)rr1=1000000.
            if(abs(rr1).lt.abs(rr))xr=xr1
        endif
        go to 15
c
c   reduction of cubic equation to linear equation
c
   14 continue
      if(abs(b).eq.0.0)then
            ind = 1
      else
            xr=-a/b
      endif
   15 yint(1)=xa+xr*(xb-xa)
      if(yint(1).lt.brd(1))ind=1
      if(yint(1).gt.brd(2))ind=2
      if(ind.eq.1.or.ind.eq.2)return
   17 aux=yint(1)-x2
      yint(2)=((d2*aux+c2)*aux+b2)*aux+a2
      aux1=(3.*d2*aux+2.*c2)*aux+b2
      aux2=1./(aux1*aux1+1.)
      aux2=sqrt(aux2)
      anyx=-aux1*aux2
      anyz=aux2
      rc=6.*d2*aux+2.*c2
      rc=rc*anyz*anyz*anyz
      return
      end

      subroutine rsvp(iflag,azm,fdip,rake,toa,p,sv)
c**************************************************
c  p-wave and sv-wave radiation patterns for arbitrary
c  orientation of a double-couple source in a homogeneous,
c  elastic full-space.
c-----input parameters:
c   iflag=0  for p-wave radiation only.
c        =1  for sv-wave radiation only.
c        =2  for both.
c   azm  =   azimuth in degrees to station minus
c            azimuth in degrees of fault strike
c            (both positive clockwise from north).
c   fdip =   dip of fault in degrees from horizontal
c            (positive down to the right when
c             facing in direction of fault strike).
c   rake =   angle in degrees between slip vector and
c            line made by intersection of fault plane
c            and horizontal plane (measured in fault
c            plane, positive counterclockwise from
c            horizontal line pointing in the fault strike
c            direction).
c   toa  =   ray take-off angle measured in degrees from
c            downward vertical.
c   note on coordinate system: see figure 4.20 in text
c            by Aki and Richards, p.114.
c-----output:
c   p    =   dimensionless p-wave radiation amplitude.
c   sv   =   dimensionless sv-wave radiation amplitude.
c   note on formulas: computations use equations 4.84 and
c            4.85 from Aki and Richards, p.115.
c-----written by g.Zandt, June 1982.
c***************************************************
        integer NRYSEG
        parameter (NRYSEG=100)
        integer KRAYSG
        parameter (KRAYSG=1200)
      common/ray/ay(10,KRAYSG),ds(9,NRYSEG),kint(NRYSEG),
     1  ros,mreg,n,iref,ind,ind1,iy(KRAYSG),ips(KRAYSG)
c-----
c   ay(5,1) is velocity at the source
c   ros is the density at the source
c-----
      fact = 4.*3.1415927*ros*ay(5,1)**3
c-----
c   fact is the source normalization for dislocation sources
c   it is 4 pi rho alpha^3 for P
c         4 pi rho  beta^3 for S
c   The knowledge of P,S is not directly here since
c   the calling sequence makes this choice and
c   since the velocity at the source is wavetype specific
c-----
c- initialize p and sv.
      p=9.9
      sv=9.9
c- convert degrees to radians
      rad=0.01745329
      az=azm*rad
      fd=fdip*rad
      rk=rake*rad
      ta=toa*rad
c- compute common terms
c-----azimuth terms
      sinaz=sin(az)
      cosaz=cos(az)
      sinaz2=sinaz*sinaz
      sin2az=2.0*sinaz*cosaz
c-----dip terms
      sinfd=sin(fd)
      cosfd=cos(fd)
      sin2fd=2.0*sinfd*cosfd
      cos2fd=1.0-2.0*sinfd*sinfd
c-----rake terms
      sinrk=sin(rk)
      cosrk=cos(rk)
c-----take-off angle terms
      sinta=sin(ta)
      costa=cos(ta)
      sinta2=sinta*sinta
      costa2=costa*costa
      sin2ta=2.0*sinta*costa
      cos2ta=costa2-sinta2
c- compute radiation amplitudes
      if(iflag.eq.1) go to 20
c-----compute p-wave radiation if iflag=0 or 2
      term1=cosrk*sinfd*sinta2*sin2az
      term2=cosrk*cosfd*sin2ta*cosaz
      term3=sinrk*sin2fd*(costa2-sinta2*sinaz2)
      term4=sinrk*cos2fd*sin2ta*sinaz
      p=term1-term2+term3+term4
      p=p/fact
c
 20   if(iflag.eq.0) go to 30
c-----compute sv-wave radiation if iflag=1 or 2
      term1=sinrk*cos2fd*cos2ta*sinaz
      term2=cosrk*cosfd*cos2ta*cosaz
      term3=0.5*cosrk*sinfd*sin2ta*sin2az
      term4=0.5*sinrk*sin2fd*sin2ta*(1.0+sinaz2)
      sv=term1-term2+term3-term4
      sv=sv/fact
c- computations complete, return
 30   return
      end
c
c   *****************************************************
c
      subroutine tramp(xsour,zsour,tsour,dt,ac,tmax,itmax,mout,
     1mdim,ncode,LU1,method,ldoray)
c
c   two-point ray tracing
c
c   ****************************************************************
c
c
c   A short description of the universal two-point ray tracing
c   routine
c   tramp(xsour,zsour,dt,ac,tmax,tstop,itmax,mout,mdim,
c   ncode,LU1,method).
c   ******************************
c
c   The routine tramp is designed for a two-point ray tracing. It
c   computes the rays which arrive at a system of receivers distri-
c   buted regularly or irregularly along the earths surface. The
c   shooting method is used. The basic system of the initial angles
c   of rays at the source is specified in the input data card no.12.
c   To be sure that no receivers are outside the range of initial
c   angles, the whole range of angles (e.g. from 3.1416 to -3.1416)
c   can be used. To save the computing time, the range may be selected
c   narrower. From the computed termination points of these rays, the
c   rays arriving at specified receivers are iteratively determined.
c   The travel time curve may be arbitrary, with loops, caustic,
c   shadows, etc.
c   A special care is devoted to certain irregular situations in the
c   ray field. The routine is automatically looking for all boundary
c   rays(i.e. the rays corresponding to the boundaries of shadow
c   zones at the earths surface, to the left-hand side and right-hand
c   side upper corners of the model, etc.). A special care is also
c   devoted to the problem of slightly reflected waves (which have
c   rays similar to head waves). As the routine uses only a single
c   arithmetics, all the rays of this wave, especially those rays
c   which are very close to the critical ray, cannot be determined.
c   In a future modification of the program, it would be useful to
c   supplement the program by some simple bending procedure for this
c   situation.
c
c   The meaning of the input and output parameters is as follows:
c
c      xsour,zsour,tsour,dt,ac... see input card no. 11 and 12
c      tmax...  maximum travel time permitted. In seis81, the travel-
c               time selection is not performed at this stage.
c               Formally , tmax=10000.
c      itmax,mout,method... see input card no.9.
c      LU1,..See input card no.1.
c      ncode... Internal wave code. See details in the description
c               of routine gener.
c
c   The information on the system of receiver positions along the earths
c   surface, on the basic system of initial angles for ray shooting,
c   and on the required accuracy of computations (see reps, input card
c   no.11) is transmitted to tramp by the coomon /dist/.
c   Routines used : ray2, ampl,jpar.
c   The routine does not need any information on the velocity distri-
c   bution in the model and on the history of rays. It uses only the
c   information on the x-coordinate of the termination point, on the
c   initial angles of rays and on the parameter ind.In other words,
c   the routine can work with any other initial value ray tracing
c   routine which provides these quantities. In fact, it has been
c   successfully aplied in various other programs, even where the
c   ray parameter has a diferent meaning than the initial angle (see
c   e.g. the program syns81 in this package).
c     The routine generates two files, first is on the logical unit
c   LU1 (two-point ray diagrams, travel timers and amplitudes),
c   and the second on the logical unit lu2 (impulse synthetic seismo-
c   grams). The content of files is described elsewhere. Similarly,
c   the output on the line printer is described elsewhere.
c
c   ******************************************************************
        logical ldoray
        integer NBDLYR
        parameter (NBDLYR=15)
        integer NRYSEG
        parameter (NRYSEG=100)
        integer NDIST
        parameter (NDIST=100)
        integer KRAYSG
        parameter (KRAYSG=1200)
      integer code
      real time(200),xcoor(200),ampx(200),ampz(200),zcoor(200)
        integer indi(200)
      common/ray/ay(10,KRAYSG),ds(9,NRYSEG),kint(NRYSEG),
     1  ros,mreg,n,iref,ind,ind1,iy(KRAYSG),ips(KRAYSG)
      common/dist/dst(NDIST),ndst,idist,astart,astep,afin,reps
      common/cod/code(NRYSEG),kref,kc
      common/eq/q(KRAYSG), eta(KRAYSG)
      common/den/rho1(NBDLYR),rho2(NBDLYR),nrho
        common/medim/v(300),sx(350),sy(350),nx(20),ny(20),nvs(NBDLYR),
     1      ptos(NBDLYR)
c
c
c
       aa=astart-astep
      index=0
      inum=0
      inds=0
c
c   loop in ray parameters aa, from astart to afin, with the step
c   astep
c
    1 aa=aa+astep
      if(astep.gt.0..and.aa.gt.afin)go to 99
      if(astep.lt.0..and.aa.lt.afin)go to 99
      ind=inds
      call ray2(xsour,zsour,tsour,aa,x,z,t,fi,tmax,dt,ac)
      if(mout.eq.2)write(6,100)ind,ind1,x,z,t,aa
      if(inum.ne.0)go to 2
      aold=aa
      iold=ind
      xold=x
      inum=1
      go to 1
c
c   parameters for the preceding ray: aa=aold, x=xold, ind=iold
c   parameters for the new ray: aa=anew, x=xnew, ind=inew
c
    2 inew=ind
      anew=aa
      xnew=x
      if(inew.eq.3.and.iold.eq.3)go to 40
      if(inew.eq.3)go to 50
      if(iold.eq.3)go to 55
      if(inew.eq.9.and.iold.ne.9)go to 30
      if(inew.ne.9.and.iold.eq.9)go to 35
c
c   no iterations, take a new ray in the loop
c
    3 iold=inew
      xold=xnew
      aold=anew
      go to 1
c
c   regular case: iold=3, inew=3
c
   40 xxnew=xnew
      xxold=xold
      aanew=anew
      aaold=aold
   41 if(xxnew.gt.xxold)go to 46
c
c   regular case: xxnew.le.xxold, itrend=-1 (reverse branch)
c
      itrend=-1
      do 42 i=1,ndst
      ndd=ndst-i+1
      if(dst(ndd).gt.xxold)go to 42
      if(dst(ndd).lt.xxnew)go to 3
      dd=dst(ndd)
      ii=ndd
      go to 43
   42 continue
      go to 3
c
c   regular case: xxnew.gt.xxold, itrend=1 (regular branch)
c
   46 itrend=1
      do 44 i=1,ndst
      if(dst(i).lt.xxold)go to 44
      if(dst(i).gt.xxnew)go to 3
      dd=dst(i)
      ii=i
      go to 43
   44 continue
      go to 3
   43 p1=aaold
      p2=aanew
      x1=xxold
      x2=xxnew
c
c   i t e r a t i o n s
c
   60 iter=0
      isign=1
   61 iter=iter+1
      isign=-isign
      pnew=0.5*(p1+p2)
      if(isign.gt.0.and.method.eq.1)pnew=p1+(dd-x1)*(p2-p1)/
     1(x2-x1)
      if(method.eq.0)pnew=p1+(dd-x1)*(p2-p1)/(x2-x1)
      if(iter.gt.itmax)go to 80
      ind=inds
      call ray2(xsour,zsour,tsour,pnew,x,z,t,fi,tmax,dt,ac)
      if(mout.eq.2)write(6,101)ind,ind1,iter,x,z,t,pnew,dd
      if(ind.ne.3)p2=pnew
      if(ind.ne.3)go to 61
      if(abs(x-dd).lt.reps)go to 65
      if(x1.lt.x2.and.dd.gt.x)go to 63
      if(x1.gt.x2.and.dd.lt.x)go to 63
      p2=pnew
      x2=x
      go to 61
   63 p1=pnew
      x1=x
      go to 61
c
c   successfull iterations
c
   65 index=index+1
      if(mout.eq.2)write(6,102)ind,ind1,iter,x,z,t,pnew,dd,ii,index
      if(mout.eq.2.and.mdim.eq.0)write(6,114)
      if(mout.eq.1.and.mdim.eq.0)write(6,113)dd,x,z,t,pnew,ind,
     1ind1,iter,ii,inew
c-----
c     output the ray description
c     consisting of
c     
c     n = number of segments
c     ind =
c     ay(2,i) = X - coordinate
c     ay(3,i) = Z - coordinate
c     ay(1,i) = Time
c     iy(i)   = layer number
c     Q(i)    = Q at this point
c     eta = frequency dependence of Q for this ray element
c-----
        if(LU1.ne.0 )then
                write(LU1)n,ind,ii,ncode
            if(n.gt.0 .and. ldoray)then
                write(LU1)(ay(2,i),ay(3,i),ay(1,i),iy(i),
     1          ips(i),q(i),eta(i),i=1,n)
            endif
        endif
      time(index)=t
      xcoor(index)=x
      zcoor(index)=z
      if(mdim.eq.0)go to 80
      if(mdim.eq.1)spr=1.
      if(mdim.eq.1)go to 10
      n=0
      q0=0.
      p0=1./ay(5,1)
      call jpar(q0,p0,spr,p,0)
   10 call ampl(spr,ax,az,phx,phz,at,pht,ish)
        sumtq = 0.0
      if(mdim.lt.3)go to 12
        aux=0.0
c-----
c   solution to dynamic ray tracing in the case of gradient
c   or constant velocity. e.g., q = integral v ds
c
c     at the same time get the INT T/Q for Q
c-----
        do 11 i=2,n
            tdif=ay(1,i)-ay(1,i-1)
            vaver=0.5*(ay(5,i)+ay(5,i-1))
            vaver2=vaver*vaver
            aux=aux+vaver2*tdif

            sumtq = sumtq + 0.5*tdif*(q(i-1) + q(i))
   11   continue
        if(aux.lt.0.000001)aux=0.001
        aux=sqrt(ay(5,1)/aux)
        spr=spr/aux
        ax=ax*aux
        az=az*aux
        at=at*aux
   12 continue
        if(code(1).lt.0)then
                mwave = 2
        else
                mwave = 1
        endif
c-----
c   at this point we have everything except the source
c   to do the source we need to know the ray parameter,
c   and the velocity and density at the source
c   This will be output for use by a rayto10 program
c   after passing through a source pulse filter
c-----
c     first get velocity and density at surface
c-----
        vsrf = ay(5,n)
        ipors = ips(n)
        if(ipors.gt.0)then
            vpsrf = vsrf
        else
            vpsrf = vsrf * ptos(1)
        endif
        rsrf = rho1(1) + vpsrf*rho2(1)
        write(8,116)ncode,ii,t,ax,az,phx,phz,pnew,mwave,ros,
     1      ay(5,1),sumtq,rsrf,vsrf
c-----
c   if SH computed output it also
c-----
      if(ish.gt.0)then
           mshwv=3
           taz=0.0
           tphz=0.0
           write(8,116)ncode,ii,t,at,taz,pht,tphz,pnew,mshwv,ros,
     1      ay(5,1),sumtq,rsrf,vsrf
      endif
  116 format(2i5,6e11.4,i5,5e11.4)
c------
c     if the source is not a dislocation or
c     if the source is given by a table
c     then use the factor for amplitude s/amplitude p ratio
c     otherwise, this will be taken care of in the subroutine
c     source
c-----
      ampz(index)=az
      ampx(index)=ax
      indi(index)=ind
        if(mout.eq.1)then
                write(6,105)dd,x,z,t,pnew,spr,ax,phx,az,phz,ind,ind1
     1                  ,iter,ii,inew,n
        elseif(mout.eq.2)then
                write(6,110)ax,phx,az,phz,spr
        endif
c
c
c
   80 p1=pnew
      x1=x
      if(iter.gt.itmax)p1=aaold
      if(iter.gt.itmax)x1=xxold
      p2=aanew
      x2=xxnew
      if(itrend.eq.(-1))go to 85
      ii=ii+1
      if(ii.gt.ndst)go to 3
      dd=dst(ii)
      if(dd.gt.xxnew)go to 3
      go to 60
   85 ii=ii-1
      if(ii.lt.1)go to 3
      dd=dst(ii)
      if(dd.lt.xxnew)go to 3
      go to 60
c
c   e n d   o f    i t e r a t i o n s
c
c
c    boundary rays. case iold.ne.3, inew=3
c
   50 xxnew=xnew
      aanew=anew
      p1=aold
      p2=anew
   54 ires=0
      iter=0
   51 pnew=0.5*(p1+p2)
      iter=iter+1
      ind=inds
      call ray2(xsour,zsour,tsour,pnew,x,z,t,fi,tmax,dt,ac)
      if(mout.eq.2)write(6,103)ind,ind1,iter,x,z,t,pnew
      if(ind.eq.3)go to 52
      p1=pnew
      go to 53
   52 ires=1
      xxold=x
      aaold=pnew
      p2=pnew
   53 if(iter.ge.itmax.and.ires.eq.1)go to 41
      if(iter.ge.itmax)go to 3
      go to 51
c
c   boundary rays. case iold=3, inew.ne.3
c
   55 xxold=xold
      aaold=aold
      p1=aold
      p2=anew
      ires=0
      iter=0
   56 pnew=0.5*(p1+p2)
      iter=iter+1
      ind=inds
      call ray2(xsour,zsour,tsour,pnew,x,z,t,fi,tmax,dt,ac)
      if(mout.eq.2)write(6,103)ind,ind1,iter,x,z,t,pnew
      if(mout.eq.1.and.iter.eq.20)write(6,107)x,z,t,pnew,ind,ind1,iter
      if(ind.eq.3)go to 57
      p2=pnew
      go to 58
   57 ires=1
      xxnew=x
      aanew=pnew
      p1=pnew
   58 if(iter.ge.itmax.and.ires.eq.1)go to 41
      if(iter.ge.itmax)go to 3
      go to 56
c
c   critical rays. case iold.ne.9, iold.ne.9, inew=9
c
   30 iter=0
      p1=aold
      p2=anew
      ires=0
   31 pnew=0.5*(p1+p2)
      iter=iter+1
      ind=inds
      call ray2(xsour,zsour,tsour,pnew,x,z,t,fi,tmax,dt,ac)
c-----
c     ind 1   termination at left border
c     ind 2   termination at right border
c     ind 3   termination at surface
c     ind 4   termination at bottom
      if(mout.eq.2)write(6,104)ind,ind1,iter,x,z,t,pnew
      if(mout.eq.1.and.iter.eq.20)write(6,111)x,z,t,pnew,ind,ind1,iter
      if(ind.eq.9.or.ind.gt.100)go to 32
      if(ind.eq.3)go to 33
      p1=pnew
      go to 34
   32 p2=pnew
      go to 34
   33 ires=1
      p1=pnew
   34 if(iter.ge.itmax.and.ires.eq.0)go to 3
      if(iter.lt.itmax)go to 31
      p2=pnew
      xxnew=x
      aanew=p2
      p1=aold
      go to 54
c
c   critical rays. case iold=9, inew.ne.9, inew.ne.3.
c   (this version not yet  written)
c
   35 continue
      go to 3
  100 format(2x,'            ',1x,2i3,3x,3f10.4, f15.10)
  101 format(2x,'iterations  ',1x,3i3,   3f10.4,2f15.10)
  102 format(2x,'successful  ',1x,3i3,   3f10.4, f15.10,f10.4,2i3)
  103 format(2x,'boundary ray',1x,3i3,   3f10.4, f15.10)
  104 format(2x,'critical ray',1x,3i3,   3f10.4, f15.10)
  105 format(4f10.5,f15.10,e15.8,e10.3,f10.5,e10.3,f10.5,5i3,i5)
  107 format(10x,3f10.5,f15.10,3i3,5x,'boundary ray')
  110 format(25x,e10.3,f10.5,e10.3,f10.5,e15.8,/)
  111 format(10x,3f10.5,f15.10,3i3,5x,'critical ray')
  113 format(4f10.5,f15.10,5i3)
  114 format(/)
c
c
   99 n=0
      naux=0
c-----
c     indicate the end of the ray run
c-----
        write(LU1)n,naux,naux,ncode
c-----
      if(LU1.ne.0 )write(LU1)index
      if(index.eq.0)return
      if(LU1.ne.0 )then
            write(LU1)(indi(i),xcoor(i),zcoor(i),time(i),
     1          ampx(i),ampz(i),i=1,index)
        endif
      return
      end

      subroutine ampl(fj,ax,az,phx,phz,at,pht,ish)
c
c   computation of amplitudes
c
c   ****************************************************************
c
c   A short description of the universal routine for the evaluation
c   of ray amplitudes of seismic body waves along a known ray
c   ampl(fj,ax,az,phx,phz).
c   ***********************
c
c   The routine ampl computes the ray amplitudes and phase shifts
c   of arbitrary multiply-reflected (possibly converted) seismic
c   body waves in a 2-D laterally inhomogeneous media with curved
c   interfaces. The computation is performed along a known ray, and
c   all the information about the ray is transmitted to ampl by
c   common/ray/. It is also assumed that the geometrical spreading
c   was determined independently (e.g.by the procedure jpar, see
c   above). The geometrical spreading may correspond to an arbitrary
c   type of source, e.g. to a point source, line source.
c   The meaning of the input and output parameters is as follows:
c
c         fj...     geometrical spreading.
c         ax,az...  amplitudes of horizontal and vertical components.
c         phx,phz...phase shifts of horizontal and vert. components.
c------
c   addition RBH
c
c      at        amplitude tangential component
c        pht       phase tangential component
c
c       ish       = 0 no SH computed
c                 > 0 SH computed
c-----
c
c   The numerical code of the wave is transmitted to the routine
c   ampl by common/cod/.
c   Routines used: coef8.
c   Note that the model information used in ampl and coef8
c   is transmitted to these routines by common/ray/.
c   In other words, routine ampl can be used with any other
c   ray tracing routine, when the routine generates this
c   common, and when the geometrical spreading is also determined in
c   advance (it can be determined, e.g., by jpar).
c
c
c   *****************************************************************
        integer NRYSEG
        parameter (NRYSEG=100)
        integer KRAYSG
        parameter (KRAYSG=1200)
      integer code
      common/cod/code(NRYSEG),kref,kc
      common/ray/ay(10,KRAYSG),ds(9,NRYSEG),kint(NRYSEG),
     1  ros,mreg,n,iref,ind,ind1,iy(KRAYSG),ips(KRAYSG)
c
      ph=0.
      pht=0.0
      q=1.
      tq=1.
c-----
c   initially set ish = true
c   if any part of path is a P then ish = 0
c-----
      ish = 1
      v=1.
      tv=1.
      al=1.
      tal=1.
      if(kref.le.1)ic1=code(1)
      if(kref.le.1)ic2=code(1)
   13 i1=kint(iref)
      if(i1.le.0)iref=iref-1
      if(i1.le.0)go to 13
      dst=ay(2,i1)-ay(2,1)
      iref1=iref-1
      if(iref1.eq.0)go to 10
c
c   loop for interfaces
c
      do 5 i=1,iref1
      i1=kint(i)
      if(i1.le.0)go to 5
      sn1=ds(2,i)
      sn2=ds(3,i)
      p=cos(sn1)
      vp1=ds(4,i)
      vs1=ds(5,i)
      ro1=ds(8,i)
      vp2=ds(6,i)
      vs2=ds(7,i)
      ro2=ds(9,i)
      if(kref.le.1)ic1=code(1)
      if(kref.le.1)ic2=code(1)
      if(kref.le.1)go to 11
      ic1=code(i)
      ii=i
   12 ii=ii+1
      if(kint(ii).le.0)go to 12
      ic2=code(ii)
      if((iabs(ic2)-iabs(ic1)).eq.0)go to 2
c-----
c   transmitted waves
c-----
   11 if(ic1.lt.0)go to 1
      ish = 0
      p=abs(p/vp1)
      av=ro1*vp1
      if(ic2.gt.0)nc=3
      if(ic2.gt.0)av=(ro2*vp2)/av
      if(ic2.lt.0)nc=4
      if(ic2.lt.0)av=(ro2*vs2)/av
      nct = -1
      tav = 0.0
      go to 4
    1 p=abs(p/vs1)
      av=ro1*vs1
      if(ic2.gt.0)then
            nc=7
            av=(ro2*vp2)/av
            ish=0
      else
            nc=8
            av=(ro2*vs2)/av
            tav=av
            nct=10
      endif
      go to 4
c-----
c   reflected waves
c-----
    2 if(abs(vp2).lt..00001)ro2=0.
      if(ic1.lt.0)go to 3
      ish = 0
      p=abs(p/vp1)
      if(ic2.gt.0)nc=1
      if(ic2.gt.0)av=1.
      if(ic2.lt.0)nc=2
      if(ic2.lt.0)av=vs1/vp1
      nct = -1
        tav = 0.0
      go to 4
    3 p=abs(p/vs1)
      if(ic2.gt.0)then
            nc=5
            av=vp1/vs1
            ish=0
            nct=-1
      else
            nc=6
            av=1.
            tav=1.
            nct=9
      endif
      go to 4
    4 v=v*av
      tv=tv*tav
      nd=0
      if(sn1.gt.1.57079632)nd=1
      if(sn1.lt.-1.57079632)nd=1
c----
c   CHANGED
c    call coef8(p,vp1,vs1,ro1,vp2,vs2,ro2,nc,nd,r,phs)
c----
      if(nct.gt.0)then
            call coef8(p,vp1,vs1,ro1,vp2,vs2,ro2,nct,nd,tr,tph)
      endif
      call coef8(p,vp1,vs1,ro1,vp2,vs2,ro2,nc,nd,r,phs)
      q=q*r
      tq=tq*tr
      ph=ph+phs
      pht=pht + tph
      sn1=sin(sn1)
      sn2=sin(sn2)
      al=al*abs(sn1/sn2)
    5 continue
c
c   end of loop
c
   10 v0=ay(5,1)
      v=v*ros*v0
      tv=tv*ros*v0
      ro1=ds(8,iref)
      if(ro1.lt.0.5)ro1=1.7+0.2*v0
      i=kint(iref)
      v=v/(ro1*ay(5,i))
      tv=tv/(ro1*ay(5,i))
      v=sqrt(v)
      tv=sqrt(tv)
c-----
c   Cerveny and Psencik have - 1.5707 since
c   their definition of a Fourier transform
c   is different. Here we use
c   G(f) = int g(t) exp(- j 2 pi f t ) dt
c-----
      if(fj.lt.0.)ph=ph+1.57079632
      if(fj.lt.0.)pht=pht+1.57079632
      aux=abs(fj)
      al=al*aux
      fj=sqrt(al)
      u=(v*q)/fj
      tu=(tv*tq)/fj
c
c
      if(mreg.ge.0)go to 20
      ax=0.
      az=u
      at=tu
      phx=0.
      phz=ph
      return
   20 continue
c
      if(ind.ne.3)go to 8
      if(mreg.eq.1)go to 8
      vp1=ds(4,iref)
      vs1=ds(5,iref)
      ro2=0.
      sn1=ds(2,iref)
      p=cos(sn1)
c------
c sense of nd changed RBH
c-----
      nd=1
      if(sn1.lt.1.57079632.and.dst.ge.0.)nd=0
      if(sn1.lt.1.57079632.and.dst.lt.0.)nd=0
      if(ic2.lt.0)go to 6
      ish = 0
      p=abs(p/vp1)
c-----
c   RBH fix
c-----
c   CHANGED
c    call coef8(p,vp1,vs1,ro1,0.00001,0.00001,0.00001,5,nd,rx,phx)
c----
c   CHANGED
c    call coef8(p,vp1,vs1,ro1,0.00001,0.00001,0.00001,6,nd,rz,phz)
c----
      call coef8(p,vp1,vs1,ro1,0.00001,0.00001,0.00001,4,nd,rx,phx)
      call coef8(p,vp1,vs1,ro1,0.00001,0.00001,0.00001,3,nd,rz,phz)
      rt = 0.0
      tph=0.0
      go to 7
    6 p=abs(p/vs1)
c-----
c   RBH fix
c-----
      call coef8(p,vp1,vs1,ro1,0.00001,0.00001,0.00001,8,nd,rx,phx)
      call coef8(p,vp1,vs1,ro1,0.00001,0.00001,0.00001,7,nd,rz,phz)
c----
c   CHANGED
c    call coef8(p,vp1,vs1,ro1,0.00001,0.00001,0.00001,7,nd,rx,phx)
c    call coef8(p,vp1,vs1,ro1,0.00001,0.00001,0.00001,8,nd,rz,phz)
c----
       call coef8(p,vp1,vs1,ro1,0.00001,0.00001,0.00001,10,nd,rt,tph)
    7 ax=u*rx
      az=u*rz
      at=tu*rt
      phx=phx+ph
      phz=phz+ph
      pht=pht+tph
      go to 9
    8 ax=u*cos(ay(4,i))
      az=-u*sin(ay(4,i))
      at = tu
      phx=ph
      phz=ph
    9 return
      end
c
c   ****************************************************************
c
      subroutine coef8(p,vp1,vs1,ro1,vp2,vs2,ro2,ncode,nd,rmod,rph)
c   *****************************************************************
c
c
c   A short description of the universal routine for the deter-
c   mination of reflection/transmission coefficients
c   coef8(p,vp1,vs1,ro1,vp2,vs2,ro2,ncode,nd,rmod,rph).
c   ***************************************************
c
c   The routine coef8 is a standard routine for the computation of
c   reflection/transmission coefficients of plane waves at a plane
c   interface between two homogeneous media or at a free surface
c   of a homogeneous halfspace. in the second case, the con-
c   version coefficients can also be optionally computed.
c   The codes of individual coefficients are specified by the
c   following numbers
c   a/ Interface between two solid halfspaces
c   p1p1...1       p1s1...2       p1p2...3       p1s2...4
c   s1p1...5       s1s1...6       s1p2...7       s1s2...8
c   sh1sh1 9       sh1sh2 10
c   b/ Free surface (for ro2.lt.0.00001)
c   pp.....1       px.....5       px,pz...x- and z- components of the
c   ps.....2       pz.....6       coef.of conversion,incident P wave
c   sp.....3       sx.....7       sx,sz...x- and z- components of the
c   ss.....4       sz.....8       coef.of conversion,incident S wave
c
c   shsh   9       shshs  10
c-----------
c   NOTE TO FIX MAJOR PROBLEM IN AMPL
c   the correspondence is now
c
c   b/ Free surface (for ro2.lt.0.00001)
c   pp.....1     sp.....5
c     ps.....2     ss.....6
c     pz.....3     sz.....7
c     px.....4     sx.....8
c
c    shsh    9     shshs  10
c     c
c     R. B. Herrmann -- Saint Louis University 8/23/84
c
c    SH coeficients are due to R Nowack, MIT
c
c
c-----------
c   I n p u t   p a r a m e t e r s
c         p...ray parameter
c         vp1,vs1,ro1...parameters of the first halfspace
c         vp2,vs2,ro2...parameters of second halfspace. for the free
c                  surface take ro2.lt.0.00001,eg.ro2=0., and
c                  arbitrary values of vp2 and vs2
c         ncode...code of the computed coefficient
c         nd...=0  when the positive direction of the ray
c                  and the x-axis make an acute angle
c              =1  when the wave impinges on the interface
c                  against the positive direction of the x-axis
c
c   O u t p u t   p a r a m e t e r s
c         rmod,rph...modul and argument of the coefficient
c
c   N o t e s
c   1/ Positive p...in the direction of propagation
c   2/ Positive s...to the left from p
c   3/ Time factor of incident wave ... exp(-i*omega*t)
c   4/ Formulae are taken from Cerveny ,Molotkov, Psencik, Ray Method
c      in Seismology, pages 30-35. Due to the note 2, the signs at
c      certain coefficients are opposite
c
c      Written by V.C. 1976.
c
c
c   ***************************************************************
      complex b(4),rr,c1,c2,c3,c4,h1,h2,h3,h4,h5,h6,h,hb,hc
      real prmt(4),d(4),dd(4)
c-----
c   for sh components ncode .ge.9
c-----
      if(ncode.ge.9)goto 400
      if(ro2.lt.0.000001)go to 150
      prmt(1)=vp1
      prmt(2)=vs1
      prmt(3)=vp2
      prmt(4)=vs2
      a1=vp1*vs1
      a2=vp2*vs2
      a3=vp1*ro1
      a4=vp2*ro2
      a5=vs1*ro1
      a6=vs2*ro2
      q=2.*(a6*vs2-a5*vs1)
      pp=p*p
      qp=q*pp
      x=ro2-qp
      y=ro1+qp
      z=ro2-ro1-qp
      g1=a1*a2*pp*z*z
      g2=a2*x*x
      g3=a1*y*y
      g4=a4*a5
      g5=a3*a6
      g6=q*q*pp
      do 21 i=1,4
      dd(i)=p*prmt(i)
   21 d(i)=sqrt(abs(1.-dd(i)*dd(i)))
      if(dd(1).le.1..and.dd(2).le.1..and.dd(3).le.1..and.dd(4).le.1.)
     1go to 100
c
c   complex coefficients
      do 22 i=1,4
      if(dd(i).gt.1.)go to 23
      b(i)=cmplx(d(i),0.)
      go to 22
   23 b(i)= cmplx(0.,d(i))
   22 continue
      c1=b(1)*b(2)
      c2=b(3)*b(4)
      c3=b(1)*b(4)
      c4=b(2)*b(3)
      h1=g1
      h2=g2*c1
      h3=g3*c2
      h4=g4*c3
      h5=g5*c4
      h6=g6*c1*c2
      h=1./(h1+h2+h3+h4+h5+h6)
      hb=2.*h
      hc=hb*p
      go to (1,2,3,4,5,6,7,8),ncode
    1 rr=h*(h2+h4+h6-h1-h3-h5)
      go to 26
    2 rr=vp1*b(1)*hc*(q*y*c2+a2*x*z)
      if(nd.ne.0)rr=-rr
      go to 26
    3 rr=a3*b(1)*hb*(vs2*b(2)*x+vs1*b(4)*y)
      go to 26
    4 rr=-a3*b(1)*hc*(q*c4-vs1*vp2*z)
      if(nd.ne.0)rr=-rr
      go to 26
    5 rr=-vs1*b(2)*hc*(q*y*c2+a2*x*z)
      if(nd.ne.0)rr=-rr
      go to 26
    6 rr=h*(h2+h5+h6-h1-h3-h4)
      go to 26
    7 rr=a5*b(2)*hc*(q*c3-vp1*vs2*z)
      if(nd.ne.0)rr=-rr
      go to 26
    8 rr=a5*b(2)*hb*(vp1*b(3)*y+vp2*b(1)*x)
      go to 26
c   real coefficients
  100 e1=d(1)*d(2)
      e2=d(3)*d(4)
      e3=d(1)*d(4)
      e4=d(2)*d(3)
      s1=g1
      s2=g2*e1
      s3=g3*e2
      s4=g4*e3
      s5=g5*e4
      s6=g6*e1*e2
      s=1./(s1+s2+s3+s4+s5+s6)
      sb=2.*s
      sc=sb*p
      go to (101,102,103,104,105,106,107,108),ncode
  101 r=s*(s2+s4+s6-s1-s3-s5)
      go to 250
  102 r=vp1*d(1)*sc*(q*y*e2+a2*x*z)
      if(nd.ne.0)r=-r
      go to 250
  103 r=a3*d(1)*sb*(vs2*d(2)*x+vs1*d(4)*y)
      go to 250
  104 r=-a3*d(1)*sc*(q*e4-vs1*vp2*z)
      if(nd.ne.0)r=-r
      go to 250
  105 r=-vs1*d(2)*sc*(q*y*e2+a2*x*z)
      if(nd.ne.0)r=-r
      go to 250
  106 r=s*(s2+s5+s6-s1-s3-s4)
      go to 250
  107 r=a5*d(2)*sc*(q*e3-vp1*vs2*z)
      if(nd.ne.0)r=-r
      go to 250
  108 r=a5*d(2)*sb*(vp1*d(3)*y+vp2*d(1)*x)
      go to 250
c
c   earths surface,complex coefficients and coefficients of conversion
  150 a1=vs1*p
      a2=a1*a1
      a3=2.*a2
      a4=2.*a1
      a5=a4+a4
      a6=1.-a3
      a7=2.*a6
      a8=2.*a3*vs1/vp1
      a9=a6*a6
      dd(1)=p*vp1
      dd(2)=p*vs1
      do 151 i=1,2
  151 d(i)=sqrt(abs(1.-dd(i)*dd(i)))
      if(dd(1).le.1..and.dd(2).le.1.)go to 200
      do 154 i=1,2
      if(dd(i).gt.1.)go to 155
      b(i)=cmplx(d(i),0.)
      go to 154
  155 b(i)= cmplx(0.,d(i))
  154 continue
      h1=b(1)*b(2)
      h2=h1*a8
      h=1./(a9+h2)
c    go to (161,162,163,164,165,166,167,168),ncode
c-----
c   fix to make the ncodes correspond
c-----
      go to (161,162,166,165,163,164,168,167),ncode
  161 rr=(-a9+h2)*h
      go to 26
  162 rr=a5*b(1)*h*a6
      if(nd.ne.0)rr=-rr
c-----
c   fix by RBH
c-----
      rr = -rr
      go to 26
  163 rr=a5*b(2)*h*a6*vs1/vp1
      if(nd.ne.0)rr=-rr
      go to 26
  164 rr=-(a9-h2)*h
      go to 26
  165 rr=a5*h1*h
      if(nd.ne.0)rr=-rr
      go to 26
  166 rr=a7*b(1)*h
      go to 26
  167 rr=a7*b(2)*h
      go to 26
  168 rr=-a5*vs1*h1*h/vp1
      if(nd.ne.0)rr=-rr
   26 z2=real(rr)
      z3=aimag(rr)
      if(z2.eq.0..and.z3.eq.0.)go to 157
      rmod=sqrt(z2*z2+z3*z3)
      rph=atan2(z3,z2)
        goto 300
  157 rmod=0.
      rph=0.
        goto 300
c
c   earths surface,real coefficients and coefficients of conversion
  200 s1=d(1)*d(2)
      s2=a8*s1
      s=1./(a9+s2)
c    go to (201,202,203,204,205,206,207,208),ncode
c-----
c   fix to make ncodes correspond
c-----
      go to (201,202,206,205,203,204,208,207),ncode
  201 r=(-a9+s2)*s
      go to 250
  202 r=a5*d(1)*s*a6
      if(nd.ne.0)r=-r
c----
c   fix by RBH
c----
      r = -r
      go to 250
  203 r=a5*d(2)*s*a6*vs1/vp1
      if(nd.ne.0)r=-r
      go to 250
  204 r=(s2-a9)*s
      go to 250
  205 r=a5*s1*s
      if(nd.ne.0)r=-r
      go to 250
  206 r=a7*d(1)*s
      go to 250
  207 r=a7*d(2)*s
      go to 250
  208 r=-a5*vs1*s1*s/vp1
      if(nd.ne.0)r=-r
      goto 250
  250   continue
        if(r.lt.0.0)then
                rmod = -r
                rph = -3.1415927
        else
                rmod = r
                rph = 0.0
        endif
      go to 300
  400 continue
      if(ro2.lt.0.00001)then
            rph = 0.0
            if(ncode.eq.9)then
                  r = +1
            elseif(ncode.eq.10)then
                  r = +2
            endif
            rmod = r
            rph  = 0.0
      else
            a1=vs1*ro1
            a2=vs2*ro2
            dd(1)=p*vs1
            dd(2)=p*vs2
            d(1)=sqrt(abs(1.- (dd(1))**2))
            d(2)=sqrt(abs(1.- (dd(2))**2))
            if(dd(1).le.1.0.and.dd(2).le.1.0)then
                  s=a1*d(1) + a2*d(2)
                  if(ncode.eq.9)then
                        r = (a1*d(1)-a2*d(2))/s
                  elseif(ncode.eq.10)then
                        r = (2.*a1*d(1))/s
                  endif
                  rph=0.0
                  rmod=abs(r)
                  if(r.lt.0.0)then
                        rph=3.1415927
                  endif
            else
                  do 302 i=1,2
                        if(dd(i).gt.1.)then
                              b(i)=cmplx(0.0,-d(i))
                        else
                              b(i)=cmplx(d(i),0.0)
                        endif
  302             continue
                  h1 = a1*b(1) + a2*b(2)
                  if(ncode.eq.9)then
                        rr =(a1*b(1)-a2*b(2))/h1
                  elseif(ncode.eq.10)then
                        rr = 2.*a1*b(1)/h1
                  endif
                  z2=real(rr)
                  z3=aimag(rr)
                  if(z2.eq.0.0.and.z3.eq.0.0)then
                        rmod=0.0
                        rph=0.0
                  else
                        rmod=cabs(rr)
                        rph=atan2(z3,z2)
                  endif
            endif
      endif
  300   continue
c------
c   reverse phase for Herrmann convention of Fourier Transform
c-----
      if(ncode.lt.9)rph = -rph
c       write(6,*)'coef8:p,vp1,vs1,ro1,vp2,vs2,ro2,ncode,nd,rmod,rph'
c       write(6,*)p,vp1,vs1,ro1,vp2,vs2,ro2,ncode,nd,rmod,rph
        return
        end
c
c   ***************************************************************
c
      subroutine jpar(q0,p0,q,p,iprint)
c
c   dynamic ray tracing
c
c   ******************************************************************
c
c
c   A short description of the universal  routine for 2-D dynamic
c   ray tracing along a known ray jpar(q0,p0,q,p,iprint).
c   *****************************************************
c
c
c   When the ray is known and its history is stored in common/ray/,
c   2-D dynamic ray tracing can be performed along the ray by
c   the procedure jpar.
c   The meaning of the input and output parameters is as follows:
c
c      q0,p0... initial values for the dynamic ray tracing
c               system.
c      q,p...   the solutions of the dynamic ray tracing system
c               at the termination point of the ray. The quan-
c               tity q represents the spreading function for a
c               line source. The spreading function for the
c               point source can be easily obtained by well-known
c               methods
c      iprint.. controls the printing of dynamic ray tracing
c               along the ray on a line printer. In seis81, no
c               print on the line printer of dynamic ray tracing
c               results is performed, iprint=0.
c
c   The numerical code of the wave is transmitted to the routine
c   jpar by common/cod/.
c   Routines functions used: rk,fcta.
c   Note that all the information about the model used in jpar,rk
c   and fcta is transmitted to these routines only by common/ray/.
c   In other words, the dynamic ray tracing package can also be
c   used with other ray tracing routines, when they generate
c   common/ray/ as described above.
c
c
c   ****************************************************************
        integer NRYSEG
        parameter (NRYSEG=100)
        integer KRAYSG
        parameter (KRAYSG=1200)
      integer code
      real y(2)
      common/ray/ay(10,KRAYSG),ds(9,NRYSEG),kint(NRYSEG),
     1  ros,mreg,n,iref,ind,ind1,iy(KRAYSG),ips(KRAYSG)
      common/cod/code(NRYSEG),kref,kc
c
      irf=0
      x=ay(1,1)
      y(1)=q0
      y(2)=p0
    2   continue
        i1=kint(iref)
        if(i1.le.0)iref=iref-1
        if(i1.le.0)go to 2
c
    8   continue
        irf=irf+1
      i1=kint(irf)
      if(i1.le.0)go to 8
    1   continue
            i2=i1+1
            call rk(x,y,irf,iprint)
      if(n.eq.kint(iref))go to 4
      d11=0.5*ds(1,irf)
      if(kref.le.1)go to 6
      ic1=code(irf)
    6   continue
            alfa=ds(2,irf)
            cs1=cos(alfa)
            sn1=sin(alfa)
            beta=ds(3,irf)
            sn2=sin(beta)
            delta=1.57079632-ay(4,i1)
            cs=cos(delta)
            sn=sin(delta)
            v1=ay(5,i1)
            vx=ay(6,i1)
            vy=ay(7,i1)
            akapa1=(vx*cs-vy*sn)/v1
            vs1=vx*sn+vy*cs
            delta=1.57079632-ay(4,i2)
            cs=cos(delta)
            sn=sin(delta)
            v2=ay(5,i2)
            vx=ay(6,i2)
            vy=ay(7,i2)
            akapa2=(vx*cs-vy*sn)/v2
            vs2=vx*sn+vy*cs
            s1=2.*(akapa1*sn1-akapa2*sn2)*cs1/v1
            s1=s1+2.*d11*(sn1/v1-sn2/v2)
            s1=s1+(vs1-vs2)*(cs1*cs1/(v1*v1))
            y(2)=(y(2)*sn1-y(1)*s1/sn1)/sn2
            y(1)=y(1)*sn2/sn1
    7       continue
                irf=irf+1
                i1=kint(irf)
            if(i1.le.0)go to 7
            if(kref.le.1)go to 1
            ic2=code(irf)
      if(iabs(ic2).ne.iabs(ic1))go to 1
      y(1)= -y(1)
      y(2)= -y(2)
      go to 1
c
    4   continue
        q=y(1)
        p=y(2)
        return
        end
c
c   **************************************************************
c
      subroutine rk(x,y,irf,iprint)
c
c   modified euler's method to solve a system of of ordinary
c   differential equations of first order
c
        integer NRYSEG
        parameter (NRYSEG=100)
        integer KRAYSG
        parameter (KRAYSG=1200)
      real y(2),dery(2),y1(2),y2(2)
      common/ray/ay(10,KRAYSG),ds(9,NRYSEG),kint(NRYSEG),
     1  ros,mreg,n,iref,ind,ind1,iy(KRAYSG),ips(KRAYSG)
c
      if(iprint.eq.3)write(6,100)x,y(1),y(2)
      i1=kint(irf)
      n=n+1
    1 h=ay(1,n+1)-ay(1,n)
      x=x+h
      call fcta(y,dery,2)
      do 2 i=1,2
      y1(i)=dery(i)
    2 y2(i)=y(i)+h*y1(i)
      n=n+1
      call fcta(y2,dery,2)
      do 3 i=1,2
    3 y(i)=y(i)+.5*h*(y1(i)+dery(i))
      if(iprint.eq.3)write(6,100)x,y(1),y(2)
  100 format(2x,3f15.7)
      if(n.eq.i1) go to 4
      go to 1
    4 return
      end
c
c   ***************************************************************
c
      subroutine fcta(y,dery,ndim)
c
c   computation of the right-hand sides of the dynamic ray tracing
c   system
c
        integer NRYSEG
        parameter (NRYSEG=100)
        integer KRAYSG
        parameter (KRAYSG=1200)
      real y(ndim),dery(ndim)
      common/ray/ay(10,KRAYSG),ds(9,NRYSEG),kint(NRYSEG),
     1  ros,mreg,n,iref,ind,ind1,iy(KRAYSG),ips(KRAYSG)
c
      v=ay(5,n)
      delta=1.57079632-ay(4,n)
      cs=cos(delta)
      sn=sin(delta)
      vxx=ay(8,n)
      vxy=ay(9,n)
      vyy=ay(10,n)
      v22=vxx*cs*cs-2.*vxy*cs*sn+vyy*sn*sn
      dery(1)=v*v*y(2)
      dery(2)=-(v22*y(1))/v
      return
      end
c
c   *************************************************************
c
        subroutine veloc(y,vel)
c
c   determination of velocity and its derivatives
c   for bicubic polynomial approximation
c
c   *************************************************************
c
c
c   A short description of the routine for the determination of
c   velocity and its derivatives at an arbitrary point of the
c   medium  veloc(y,vel).
c   *********************
c
c   In this routine, the specified point at which the velocity
c   is to be determined is first localized in the given velo-
c   city network. Then the values of velocity and of its first
c   partial derivatives are determined. if isr.ne.0, the
c                     (mdim.gt.1 seis81)
c   second derivatives are also determined.
c   The meaning of the input and output parameters is as follows
c   follows:
c
c     y(1),y(2)...x and z coordinates of the point at which the
c                 velocity is to be determined.
c     y(3)...     value of the corresponding angle of the ray. it
c                 has only a formal meaning in the routine.
c     vel(1)...   velocity at the point y(1),y(2).
c     vel(2),vel(3)... first partial derivatives vx and vz
c     vel(4),vel(5),vel(6)... second partial derivatives vxx, vxz
c                 and vzz.
c
c   Note that for the proper operation of routine veloc, all the
c   parameters of the common/auxx/ must be specified a priori
c   (see the description of common/auxx/). Also the parameters
c   lay,itype and nder,(see the description of common/auxi/), must
c   be specified before the use of routine veloc.
c
c
c   *****************************************************************
        integer NBDLYR
        parameter (NBDLYR=15)
        integer NRYSEG
        parameter (NRYSEG=100)
        integer KRAYSG
        parameter (KRAYSG=1200)
        real y(3),vel(6)
        common/medim/v(300),sx(350),sy(350),nx(20),ny(20),nvs(NBDLYR),
     1      ptos(NBDLYR)
        common/vcoef/a02(300),a20(300),a22(300)
        common/ray/ay(10,KRAYSG),ds(9,NRYSEG),kint(NRYSEG),
     1      ros,mreg,n,iref,ind,ind1,iy(KRAYSG),ips(KRAYSG)
        common/auxx/mmx(20),mmy(20),mmxy(20),maux
      common/auxi/intr,int1,iout,irefr,lay,itype,nder,iprint,mprint,ntr
      common/pq/qp1(NBDLYR),qp2(NBDLYR),etap1(NBDLYR),etap2(NBDLYR)
      common/sq/qs1(NBDLYR),qs2(NBDLYR),etas1(NBDLYR),etas2(NBDLYR)
      common/eq/q(KRAYSG), eta(KRAYSG)
      common/eqdoit/idoit
c
        mx=nx(lay)
        my=ny(lay)
        i2 = 0
        j2 = 0
        do 1 i=2,mx
            k=mmx(lay)+i-1
            if(y(1).lt.sx(k))go to 2
    1   continue
    2   i1=k
        do 3 i=2,my
            k=mmy(lay)+i-1
            if(y(2).lt.sy(k))go to 4
    3   continue
    4   j1=k
        if(maux.eq.0 .or. i1.ne.i2 .or. j1.ne.j2)then
            i=i1-mmx(lay)
            j=j1-mmy(lay)
            i2=i1
            j2=j1
            xx=sx(i1-1)
            yy=sy(j1-1)
            mm=mmxy(lay)-1
            k=mm+(i-1)*my+j
            b20=a20(k)
            b02=a02(k)
            b22=a22(k)
            b00=v(k)
            k=mm+i*my+j
            c20=a20(k)
            c02=a02(k)
            c22=a22(k)
            c00=v(k)
            k=mm+(i-1)*my+j+1
            d20=a20(k)
            d02=a02(k)
            d22=a22(k)
            d00=v(k)
            k=mm+i*my+j+1
            e20=a20(k)
            e02=a02(k)
            e22=a22(k)
            e00=v(k)
            hx=sx(i1)-xx
            hy=sy(j1)-yy
            xa=3.*hx
            ya=3.*hy
            xb=hx/3.
            yb=hy/3.
            d32=(e22-d22)/xa
            d30=(e20-d20)/xa
            b30=(c20-b20)/xa
            b32=(c22-b22)/xa
            d12=(e02-d02)/hx-xb*(e22+2.*d22)
            d10=(e00-d00)/hx-xb*(e20+2.*d20)
            b10=(c00-b00)/hx-xb*(c20+2.*b20)
            b12=(c02-b02)/hx-xb*(c22+2.*b22)
            b03=(d02-b02)/ya
            b13=(d12-b12)/ya
            b23=(d22-b22)/ya
            b33=(d32-b32)/ya
            b01=(d00-b00)/hy-yb*(d02+2.*b02)
            b11=(d10-b10)/hy-yb*(d12+2.*b12)
            b21=(d20-b20)/hy-yb*(d22+2.*b22)
            b31=(d30-b30)/hy-yb*(d32+2.*b32)
            maux=1
        endif
        ax=y(1)-xx
        az=y(2)-yy
        aux1=((b33*az+b32)*az+b31)*az+b30
        aux2=((b23*az+b22)*az+b21)*az+b20
        aux3=((b13*az+b12)*az+b11)*az+b10
        aux4=((b03*az+b02)*az+b01)*az+b00
        vel(1)=((aux1*ax+aux2)*ax+aux3)*ax+aux4
        vel(2)=(3.*aux1*ax+2.*aux2)*ax+aux3
        if(nder.ne.0)then
            vel(4)=6.*aux1*ax+2.*aux2
        endif
        aux1=(3.*b33*az+2.*b32)*az+b31
        aux2=(3.*b23*az+2.*b22)*az+b21
        aux3=(3.*b13*az+2.*b12)*az+b11
        aux4=(3.*b03*az+2.*b02)*az+b01
        vel(3)=((aux1*ax+aux2)*ax+aux3)*ax+aux4
        if(nder.ne.0)then
            vel(5)=(3.*aux1*ax+2.*aux2)*ax+aux3
            aux1=3.*b33*az+b32
            aux2=3.*b23*az+b22
            aux3=3.*b13*az+b12
            aux4=3.*b03*az+b02
            vel(6)=2.*(((aux1*ax+aux2)*ax+aux3)*ax+aux4)
        endif
c-----
c     get the P velocity
c-----  
        if(idoit.ne.0)then
        if(nvs(lay).eq.0)then
            vp = vel(1)
        else
            vp = vel(1) / ptos(lay)
        endif
        endif

        if(itype.gt.0.and.nvs(lay).eq.0)go to 10
        if(itype.lt.0.and.nvs(lay).eq.1)go to 10
        if(itype.gt.0.and.nvs(lay).eq.1)go to 9
        if(ptos(lay).ge.100.)go to 12
c-----
c     compute S velocity from P-velocity, e.g., nvs = 0
c-----
        vel(1)=vel(1)/ptos(lay)
        vel(2)=vel(2)/ptos(lay)
        vel(3)=vel(3)/ptos(lay)
        if(nder.eq.0)go to 10
        vel(4)=vel(4)/ptos(lay)
        vel(5)=vel(5)/ptos(lay)
        vel(6)=vel(6)/ptos(lay)
        go to 10
   12   vel(1)=0.
        vel(2)=0.
        vel(3)=0.
        if(nder.eq.0)go to 10
        vel(4)=0.
        vel(5)=0.
        vel(6)=0.
        go to 10
    9   vel(1)=vel(1)*ptos(lay)
        vel(2)=vel(2)*ptos(lay)
        vel(3)=vel(3)*ptos(lay)
        if(nder.eq.0)go to 10
        vel(4)=vel(4)*ptos(lay)
        vel(5)=vel(5)*ptos(lay)
        vel(6)=vel(6)*ptos(lay)
   10   if(y(3).lt.100.)return
c-----
c     compute the Q and frequency dependence at this point
c     also whether P or S
c
c     itype > 0 = P
c     itype < 0 = S
c-----
            if(idoit.gt.0)then
                ips(n) = idoit
                q(n) = qp1(lay) + vp * qp2(lay)
                if(q(n).gt.1.0)q(n) = 1.0/q(n)
                eta(n) = etap1(lay) + vp * etap2(lay)
            else if(idoit.lt.0)then
                ips(n) = idoit
                q(n) = qs1(lay) + vp * qs2(lay)
                if(q(n).gt.1.0)q(n) = 1.0/q(n)
                eta(n) = etas1(lay) + vp * etas2(lay)
            endif
        ay(5,n)=vel(1)
        ay(6,n)=vel(2)
        ay(7,n)=vel(3)
        if(nder.eq.0)return
        ay(8,n)=vel(4)
        ay(9,n)=vel(5)
        ay(10,n)=vel(6)
        iy(n) = lay
        return
        end

      subroutine gcmdln(verbose,ldoray)
c-----
c   parse command line arguments and return control
c   parameters
c
c   requires subroutine mgtarg(i,name) to return
c         the i'th argument in the string name
c
c   and the function mnmarg() to return the number
c         of arguments excluding program name
c         The first argument is i = 1
c
c-----
c     verbose L   - verbose printer output
c     ldoray  L   - output ray tracing information
c-----
        logical verbose, ldoray
        character name*20
c-----
c   set up defaults for poor usage test
c-----
        verbose = .false.
        ldoray = .false.
        nmarg = mnmarg()
        i = 0
   11   i = i + 1
        if(i.gt.nmarg)goto 13
            call mgtarg(i,name)
            if(name(1:2).eq.'-R')then
                ldoray = .true.
            else if(name(1:2).eq.'-V'.or.name(1:2).eq.'-v')then
                verbose = .true.
            else if(name(1:2) .eq. '-?')then
                call usage(' ')
            else if(name(1:2) .eq. '-h')then
                call usage(' ')
            endif
        go to 11
   13   continue
      return
      end

        subroutine usage(str)
        character str*(*)
        integer LER
        parameter (LER=0)
        write(LER,*)'cseis96:',str
        write(LER,*)'Usage cseis96 [-v] [-P] [-?] [-h]'
        write(LER,*)
     1  ' -v    (default false) verbose output'
        write(LER,*)
     1  ' -R    (default false) generate file for CRAY96'
        write(LER,*)
     1  ' -?    (default false) this help screen'
        write(LER,*)
     1  ' -h    (default false) this help screen'
        stop
        end

