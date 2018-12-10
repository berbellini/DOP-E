c----------------------------------------------------------------------c
c                                                                    c
c    COMPUTER PROGRAMS IN SEISMOLOGY                                 c
c    VOLUME IV                                                       c
c                                                                    c
c    PROGRAM: TRFTNPV96                                               c
c                                                                    c
c    COPYRIGHT 2002                                                  c
c    R. B. Herrmann                                                  c
c    Department of Earth and Atmospheric Sciences                    c
c    Saint Louis University                                          c
c    221 North Grand Boulevard                                       c
c    St. Louis, Missouri 63103                                       c
c    U. S. A.                                                        c
c                                                                    c
c----------------------------------------------------------------------c
        program rftnpv96
c-----
c     CHANGES
c     21 MAY 2002 Corrected scaling for depth - previously
c             if layer thickness inversion used, only
c             the second model thicknesses were used
c     07 JUN 2002 Removed limitation on number of receivers functions
c     23 NOV 2002 - Carry through weights for observed 
c          RFTN's in tmpsrfi.18
c-----
c
c     This program plots the observed and predicted
c     receiver functions as well as the initial and current
c     current model
c-----
c     FIXES
c     04 APR 2002 - corrects trace plots when more than 26 traces
c-----
        parameter(LER=0,LIN=5,LOT=6)
        parameter(NL=200,NL2=NL+NL)
c-----
c     LIN - unit for FORTRAN read from terminal
c     LOT - unit for FORTRAN write to terminal
c     LER - unit for FORTRAN error output to terminal
c     NL  - number of layers in model
c     NL2 - number of columns in model (first NL2/2 are
c         - velocity parameters, second NL2/2 are Q values)
c-----
        integer nf10(NL)
c-----
c     machine dependent initialization
c-----
        call mchdep()
c-----
c     get the command line 
c-----
       call gcmdln(nid,zmax)
c-----
        nid = 0
c-----
c     iprog is a binary OR: 2= rftn, 1=surf
c-----
        call gttmp0(iprog,itot,nf1,nf2,nf34,nf5,nf67,nf10,nfilt,
     1      nup, dlam,qaqb,wref,invcsl,
     2      lstinv,twnmin,twnmax,iter,nurftn,pval,
     3      sigv,sigr,sigg,
     4      idtwo,idum2,idum3,idum4,idum5,
     5      rdum1,rdum2,rdum3,rdum4,rdum5)
c-----
c     iprog   = inversion type. Logical or of  rftn 2, surf 1 - if does
c             not match, terminate this run
c     itot    = total number of iterations [37]
c     nf1     = 1 estimated stdev computed from residuals
c               0 no scaling by residuals
c     nf2     = TOTAL number of Love Wave Gamma Modes
c               0 DO NOT PROCESS Love Wave Gamma Data for Q
c     nf34    = TOTAL number of Love Wave modes to process 
c          (from C and U)
c     nf5     = TOTAL number of Rayleigh Wave Gamma Modes
c               0 DO NOT PROCESS Rayleigh Wave Gamma Data for Q
c     nf67    = TOTAL number of Rayleigh Wave modes to process 
c          (from C and U)
c     nf10    = Input Format (from model file)
c               0 - Inversion a, rho fixed
c               1 - Inversion Poisson Ratio Fixed, Rho computed from Vp
c     nfilt   = smoothing parameter 
c               0  No model Weight  No smoothing
c               1  Model Weight     No smoothing
c               2  No model weight  Smoothing
c               3  Model weight     Smoothing
c     nup     = state
c             =1 partials computed 
c             =0 surfinv run =1 before
c             =2 on update of invcsl>=1 or 3
c     dlam    = damping [32]
c     qaqb    = 2.25 default  [34]
c     wref    = reference frequency 1.0 is default [33]
c     invcsl  = 0 acausal, 1 uncoupled causel, 2 coupled causal [35]
c     invdep  = 0 last inversion was for depth
c             = 1 last inversion was for velocity and Q inverse
c     lstinv  = 2,3,4,5 depending on the last inversion
c             invdep = 1 for 2,3,4 and 0 for 5
c     twnmin  These give the receiver function window
c     twnmax
c     iter    Current iteration
c     nurftn  Number of receiver functions to be read
c-----
        if((iprog/2).ne.1)then
            WRITE(LOT,*)'rftnpv96 requires a receiver ',
     1         'function inversion'
            STOP
        endif
            

c-----
c-----
c     get solution
c-----
            call pinitf('TRFTNPV96.PLT')
c-----
c     now plot the velocity/Q inverse model(s) for reference
c-----
        call pltmod(nid,zmax)
c-----
c     Now systematically plot the Receiver functions
c-----
        call pltrcv(nurftn,twnmin,twnmax,nid,zmax)
        call newpen(1)
        call pend()
        close(1,status='keep')
        end

        subroutine pltrcv(nurftn,twnmin,twnmax,nid,zmax)
c-----
c     plot the receiver functions
c     predicted is blue in background
c     observed is red on top
c-----
        integer NL, NL2
        parameter(NL=200,NL2=NL+NL)
        real dd(NL)
        integer NSAMP
        parameter(NSAMP=8192)
        real o(NSAMP), p(NSAMP)
        
        integer i, j, itrace, npage
        real vmax, dy
        real xx, yy, x0, y0
        real ylen, xlen, yht

        character ostr*17

        character kstnm*8
        integer nzyear, nzjday, nzhour, nzmin, nzmon, nzday

        integer nurftn
        real twnmin, twnmax
c-----
c     open file with all receiver functions and partial derivatives
c-----
            open(4,file='tmpsrfi.18',form='unformatted',
     1          access='sequential',status='unknown')
            rewind 4
c-----
c     INITIAL TEST VERSION 06 JAN 2002 plot up to 15 per column, 
c          2 columns
c-----
c-----
c     define plot limits for a nice plot
c     Basically we can place 28 on a page, but for a small
c     number, such as for documentation, rearrange vertically
c-----
        if(nurftn.le.8)then
            y0 = 5.5
            x0 = 2.5
            itmod = 14
            xlen = 3.0
            ylen = 0.5
        else if(nurftn.gt.8 .and. nurftn.le.26)then
            y0 = 7.0
            x0 = 2.5
            itmod = 14
            xlen = 3.0
            ylen = 0.5
        else if(nurftn.gt.26 .and. nurftn.le.39)then
            y0 = 7.0
            x0 = 2.500 
            itmod = 14 
            xlen = 1.5
            ylen = 0.5
        else if(nurftn.gt.39)then
            y0 = 7.0
            x0 = 2.500 
            itmod = 14 
            xlen = 1.5
            ylen = 0.5
        endif
            
        itrace = 0
        jtrace = 0
        npage = 1
 1000   continue
        if(jtrace.gt.0 .and. mod(jtrace,39).eq.0
     1      .and.nurftn.gt.npage*39)then
            call frame()
            call pltmod(nid,zmax)
            itrace = 0
            npage = npage + 1
        endif
        read(4,end=9000)kstnm
        read(4,end=9000)nzyear, nzjday, nzhour, nzmin, nzmon, nzday,rayp
        read(4,end=9000)nlow,nhgh,mm,mm2,iid,
     1          gaussalp,sum2o,sum2p,sum2r,redv,invwgt
        vmax = 0.0
        j = 0
        do 1100 iii=nlow,nhgh
                j = j + 1
                read(4)(dd(i),i=1,mm),o(j),p(j),dy
                if( abs(o(j)) .gt.vmax) vmax = abs( o(j) )
                if( abs(p(j)) .gt.vmax) vmax = abs( p(j) )
 1100   continue
        npt = j
        yht = ylen / 7.0
        yy = y0 - ylen*mod(itrace,itmod-1)
        xx = x0 + ((itrace)/(itmod-1))*(xlen + 1.0)
        call newpen(2)
c-----
c       plot predicted
c-----
        call pltit(xx,yy,xlen,ylen,vmax,npt,p)
C       WRITE(6,*)'XX,YY,XLEN,YLEN,VMAX,NPT:',xx,yy,xlen,ylen,vmax,npt
        call newpen(1)
        itrace = itrace + 1
        jtrace = jtrace + 1
c-----
c     put in the time scales
c-----
        call newpen(1)
        if(itrace.eq.1)then
        if(nurftn.lt.13)then
        call dolinx(x0,y0-ylen*nurftn,xlen,twnmax,twnmin,0.07,
     1      .false.,.false.,.true.,10,'Time (sec)')
        else

        call dolinx(x0,0.5,xlen,twnmax,twnmin,0.07,
     1      .false.,.false.,.true.,10,'Time (sec)')
        endif
        endif
        call dolinx(x0,y0+0.7,xlen,twnmax,twnmin,0.07,.true.,
     1      .false.,.false.,10,'Time (sec)')
        if(itrace.eq. 14)then
        call dolinx(x0+1.0+xlen,0.5,xlen,twnmax,twnmin,0.07,.false.,
     1      .false.,.true.,10,'Time (sec)')
        call dolinx(x0+1.0+xlen,y0+0.7,xlen,twnmax,twnmin,0.07,.true.,
     1      .false.,.false.,10,'Time (sec)')
        endif
        if(itrace.eq. 27)then
        call dolinx(x0+2.0+xlen+xlen,0.5,xlen,twnmax,twnmin,
     1      0.07,.false.,.false.,.true.,10,'Time (sec)')
        call dolinx(x0+2.0+xlen+xlen,y0+0.7,xlen,twnmax,twnmin,
     1      0.07,.true.,.false.,.false.,10,'Time (sec)')
        endif
c-----
c     annotate
c-----
        write(ostr,1)nzyear,nzmon,nzday,nzjday,nzhour,nzmin
    1   format(i4.4,i2.2,i2.2,'(',i3.3,')',i2.2,i2.2)
        call gright(xx+xlen+0.5,yy+ylen-3.0*yht,yht,ostr,0.0)
        call gleft (xx-0.30    ,yy+ylen-2.0*yht,yht,kstnm,0.0)
        write(ostr,'(f5.2)')gaussalp
        call gleft (xx-0.30    ,yy+ylen-3.5*yht, yht,ostr,0.0)
        write(ostr,'(f4.0)')redv
        call gleft (xx-0.30    ,yy+ylen-5.0*yht, yht,ostr,0.0)
        write(ostr,'(f6.3)')rayp
        call gleft (xx-0.30    ,yy+ylen-6.5*yht, yht,ostr,0.0)
c-----
c       plot observed
c-----
        call newpen(4)
        call pltit(xx,yy,xlen,ylen,vmax,npt,o)
        go to 1000
 9000   continue
            close(4)
        return
        end

        subroutine pltit(x0,y0,xlen,ylen,vmax,npt,x)
        real x0, y0, xlen, ylen, vmax
        integer npt
        real x(npt)
        real xx, yy
        do 1000 i=1,npt
            xx = x0 + (i-1)*xlen/npt
            yy = y0 + ylen*x(i)/vmax
            if(i.eq.1)then
                call plot(xx,yy,3)
            else
                call plot(xx,yy,2)
            endif
 1000   continue
        return
        end


      subroutine gttmp0(iprog,itot,nf1,nf2,nf34,nf5,nf67,nf10,nfilt,
     1      nup,dlam,qaqb,wref,invcsl,
     2      lstinv,twnmin,twnmax,iter,nurftn,pval,
     3      sigv,sigr,sigg,
     3      idtwo,idum2,idum3,idum4,idum5,
     4      rdum1,rdum2,rdum3,rdum4,rdum5)
        integer NL
        parameter (NL=200)
        integer nf10(NL)
c-----
c     read control file
c-----
      open(1,file='tmpsrfi.00',form='unformatted',access='sequential')
      rewind 1
      read(1) iprog,itot,nf1,nf2,nf34,nf5,nf67,nf10,nfilt,nup,dlam
     1      ,qaqb,wref,invcsl,lstinv,twnmin,twnmax,iter,nurftn,pval,
     3      sigv,sigr,sigg,
     3      idtwo,idum2,idum3,idum4,idum5,
     4      rdum1,rdum2,rdum3,rdum4,rdum5
      close(1,status='keep')
        return
        end

        subroutine pltmod(nid,dmax)
c-----
c     plot the velocity or the Q inverse model
c     
c     nid I*4 0 for velocity 1 for gamma/Q
c-----
c-----
c     common for igetmod
c-----
        integer NL, NL2, NLAY
        parameter(NL=200,NLAY=200,NL2=NL+NL)

        common/isomod/d(NLAY),va(NLAY),vb(NLAY),rho(NLAY),
     1      qa(NLAY),qb(NLAY),etap(NLAY),etas(NLAY), 
     2      frefp(NLAY), frefs(NLAY)
        common/depref/refdep
        integer mmax, iunit, iiso, iflsph, idimen, icnvel
        common/modtit/title
        character title*80
c-----
c     LIN - unit for FORTRAN read from terminal
c     LOT - unit for FORTRAN write to terminal
c     LER - unit for FORTRAN error output to terminal
c     NL  - number of layers in model
c-----
        real barr(4), qbarr(4), darr(4)
c-----
c     set the plot
c-----  
        x0 = 0.5
        y0 = 3.0
        call plot(x0,y0,-3)
        xlen = 1.5
        ylen = 4.0
c-----
c     get the velocity model twice, once for the
c     current model and then for the initial model
c     however we use both to get the plot extremes
c-----
c-----
c-----
c     get bounds for the plot
c-----
        bmin = 1.0e+38
        bmax = 0.0
        qbmin = 1.0e+38
        qbmax = 0.0
        zmin = 0.0
        zmax = 0.0
        do 2000 jmd=1,2
c-----
c         get velocity model and q(beta) inverse model
c-----
            if(jmd.eq.1)then
        call getmod(2,'tmpsrfi.17',mmax,title,iunit,iiso,iflsph,
     1      idimen,icnvel,ierr,.false.)
            else
        call getmod(2,'tmpmod96.000',mmax,title,iunit,iiso,iflsph,
     1      idimen,icnvel,ierr,.false.)
            endif
        mlyr = mmax
            depth = 0.0
        do 1000 i=1,mlyr
            if(nid.eq.0)then
                if(vb(i).gt.bmax)bmax = vb(i)
                if(vb(i).lt.bmin)bmin = vb(i)
            else if(nid.eq.1)then
                if(qb(i) .gt.1.0)qb(i) = 1.0/qb(i)
                if(qb(i).gt.qbmax)qbmax = qb(i)
                if(qb(i).lt.qbmin)qbmin = qb(i)
            endif
            if(i.lt.mlyr)depth = depth + d(i)
 1000   continue
        if(depth+d(mlyr-1).gt.zmax)zmax = depth+d(mlyr-1)
 2000   continue
        if(dmax .gt. 0.0)then
            zmax = dmax
        endif
c-----
c     safety test - force a minimum separation in km/sec units
            bmax = bmax + 0.125
            bmin = bmin - 0.125
        
        darr(1) = 0.0
        darr(2) = zmax
        WRITE(6,*)dmax,zmax
c-----
c     now scale the velocity/Q inverse axis
c-----
        barr(1) = bmin
        barr(2) = bmax
        call dvscale(bmin,bmax)
        barr(3) = bmin
        barr(4) = (bmax - bmin)/xlen
        qbarr(1) = 0.95*qbmin
        qbarr(2) = 1.05*qbmax
        call gscale(qbarr,xlen,2,1)
        qbmin = qbarr(3)
        qbmax = qbarr(3) + xlen * qbarr(4)

        yylen = 10.0
        call gscale(darr,yylen,2,1)
        dmin = 0.0
C       dmax = darr(3) + yylen*darr(4)
        dmax = darr(2)
c-----
c     create bounding box
c-----
        call box(0.0,0.0,xlen,ylen)
c-----
c     label the depth axis
c-----
            call labz(dmin,dmax,10,0.0,ylen,ylen,0,dd)
c-----
c     plot the velocity/q inverse axis
c-----
        if(nid.eq.0)then
            call dolinx(0.0,ylen,xlen,bmax,bmin,
     1          0.05,.true.,.true.,.true.,9,'VS (km/s)')
C           call myaxis(0.0,ylen,'VS (KM/S)', 9,xlen,0.0, 
C     1             barr(3),barr(4))
        else
            call myaxis(0.0,ylen,'QBINV'    , 5,xlen,0.0,
     1          qbarr(3),qbarr(4))
        endif
c-----
c     plot the velocity and Qb models
c-----
c-----
c     now read the model file two times, once for original
c     and then for current
c-----
	call gclip('ON',0.0,0.0,xlen,ylen)
        do 3000 jmd=1,2
            if(jmd.eq.1)then
        call getmod(2,'tmpmod96.000',mmax,title,iunit,iiso,iflsph,
     1      idimen,icnvel,ierr,.false.)
            call newpen(4)
            else
        call getmod(2,'tmpsrfi.17',mmax,title,iunit,iiso,iflsph,
     1      idimen,icnvel,ierr,.false.)
            call newpen(2)
            endif
        mlyr = mmax
        if(nid .eq. 0)then
            dx = 1.0/barr(4)
            call vplt(0.0,ylen,mlyr,vb,d,bmin,bmax,dx,
     1          dmin,dmax,dd,jmd)
        else if(nid.eq.1)then
            dx = 1.0/qbarr(4)
            call vplt(0.0,ylen,mlyr,qb,d,qbmin,qbmax,dx,
     1          dmin,dmax,dd,jmd)
        endif
 3000   continue
	call gclip('OFF',0.0,0.0,xlen,ylen)
c-----
c     annotate
c-----
        call newpen(2)
        call plot(0.0,-0.25,3)
        call plot(0.25*xlen, -0.25, 2)
        call symbol(0.30*xlen, -0.25, 0.10, ' ACI Model',0.0,11)
        call newpen(4)
        call plot(0.0,-0.50,3)
        call plotd(0.25*xlen,-0.50,21,0.05)
        call symbol(0.30*xlen, -0.50, 0.10, ' CIA Model',0.0,11)
c-----
c     reset the plot
c-----
        call plot(-x0,-y0,-3)
        return
        end

        subroutine box(xl,yl,xu,yu)
        call plot(xl,yl,3)
        call plot(xu,yl,2)
        call plot(xu,yu,2)
        call plot(xl,yu,2)
        call plot(xl,yl,2)
        return
        end
c---------------------------------------------------------------------c
c                                                                   c
c    COMPUTER PROGRAMS IN SEISMOLOGY                                c
c    VOLUME I                                                       c
c                                                                   c
c    PROGRAM: AXIS                                                  c
c                                                                   c
c    COPYRIGHT (C)  1986 R. B. Herrmann                             c
c                                                                   c
c    Department of Earth and Atmospheric Sciences                   c
c    Saint Louis University                                         c
c    221 North Grand Boulevard                                      c
c    St. Louis, Missouri 63103                                      c
c    U. S. A.                                                       c
c                                                                   c
c---------------------------------------------------------------------c
        subroutine myaxis(xpage,ypage,ititl,nchar,axlen,angle,
     1          firstv,deltav)
        real HTNUM,SYMHT,TICHT
       parameter (HTNUM=0.07,SYMHT=0.10,TICHT=0.07)
c-----
c     xpage,ypage  coordinates of starting point of axis in
c                  inches relative to current origin
c     ititl   axis title
c     nchar   number of characters in title
c             >0 for tic marks, numbers and title on
c                counterclockwise side
c             <0 for tic marks, numbers and title on
c                clockwise side
c     axlen   floating point axis length in inches
c     angle   angle (degrees) that x axis makes with
c             horizontal x-direction
c     firstv  scale value at (xpage,ypage)
c     deltav  change in slace between tic marks per unit
c             length
c             <0 value at xpage,ypage is a max, and values
c             decrease along the axis
c-----
c-----
c     we do not worry about the vagaries of storage of characters
c     since the address is passed on down to subroutine symbol
c-----
        character*(*) ititl
        a=1.0
        kn=nchar
        if(kn.lt.0)then
                a= -a
                kn= -kn
        endif
c-----
c     if deltav is too large invoke scientific notation
c-----
        ex = 0.0
        if(deltav.ne.0.0)then
                yex= alog10(abs(deltav))
                if(yex.lt.-2.0)then
                ex = aint(yex+0.01) - 1.0
            else if(yex.lt.-1.0)then
                ex = aint(yex+0.01) 
            endif
                if(yex.ge.2.0) ex = aint(yex + 0.01)
        endif
        xval = firstv*10.0**(-ex)
        xdel = deltav*10.0**(-ex)
        ct = cos(angle*0.01745329)
        st = sin(angle*0.01745329)
        ntic = axlen + 1.0
c-----
c     first put in numbers and title
c     adjust offset for numbers
c-----
        dx = - HTNUM
        dy = 1.5*a*HTNUM - 0.5*HTNUM
c-----
c     find initial position given rotation
c-----
        xn = xpage + dx*ct - dy*st
        yn = ypage + dy*ct + dx*st
        mtic = ntic/2
        do 100 i=1,ntic
            if(i.eq.1)then
                if(xval.lt.0.0)then
                    dll = 2.0*HTNUM
                else
                    dll = 1.0*HTNUM
                endif
                if(abs(xval).ge.10.0)dll = dll + HTNUM
                dxx = dll*ct
                dyy = dll*st
                    call number(xn+dxx,yn+dyy,HTNUM,xval,angle,2)
            else if(i.eq.ntic)then
                if(xval.lt.0.0)then
                    dll = 4.0*HTNUM
                else
                    dll = 3.0*HTNUM
                endif
                if(abs(xval).ge.10.0)dll = dll + HTNUM
                dxx = -dll*ct
                dyy = -dll*st
                    call number(xn+dxx,yn+dyy,HTNUM,xval,angle,2)
            else
                    call number(xn,yn,HTNUM,xval,angle,2)
            endif
                xval=xval+xdel
                xn=xn+ct
                yn=yn+st
c-----
c     halfway down axis put in title
c-----
                if(i.eq.mtic)then
                        z=kn
                        if(ex.ne.0.0)z=z+7.0
                        dx = -0.5*SYMHT*z + axlen*0.5
                        dy = (2.5*a-0.5)*SYMHT
                        xt=xpage +dx*ct-dy*st
                        yt=ypage +dy*ct+dx*st
                        call symbol(xt,yt,SYMHT,ititl,angle,kn)
                        if(ex.ne.0.0)then
                                z=kn+2
                                xt=xt+z*ct*SYMHT
                                yt=yt+z*st*SYMHT
                                call symbol(xt,yt,SYMHT,'*10',angle,3)
                                xt=xt+(3.0*ct-0.8*st)*SYMHT
                                yt=yt+(3.0*st+0.8*ct)*SYMHT
                                call number(xt,yt,0.7*SYMHT,ex,angle,
     1                                  -1)
                        endif
                endif
  100   continue
c-----
c     now put in tic marks
c-----
        call plot(xpage+axlen*ct,ypage+axlen*st,3)
        dx = - TICHT*st*a
        dy =   TICHT*ct*a
        a = ntic -1
        xn = xpage + a*ct
        yn = ypage + a*st
        do 200 i=1,ntic
                call plot(xn,yn,2)
                call plot(xn+dx,yn+dy,2)
                call plot(xn,yn,2)
                xn=xn-ct
                yn=yn-st
  200   continue
        return
        end

        subroutine labz(dmin,dmax,nd,x0,y0,ylen,id,ddd)
c-----
c     label the depth axis - however 
c-----
c     dmin    R*4 - minimum depth
c     dmax    R*4 - maximum depth
c     nd  I*4 - number of segments
c     x0  R*4
c     Y0  R*4 - upper left corner
c     ylen    R*4 - length of Y-axis
c     id  I*4 - 0 single
c               1 double
c-----
        character ostr*10, fmt*10
        dd = (dmax - dmin) / real (nd)
        if(id.eq.1)then 
            ndo = 2
        else
            ndo = 1
        endif
        if(ndo.eq.1)then
            dz = ylen / real(nd)
        else
            dz = 0.5 * ylen / real(nd)
        endif
        ddd = dz / dd
        if(dmax.gt.10 .and. dmax.lt.100)then
            ndec = 2
            fmt = '(f10.2)'
        else if(dmax.ge.100)then
            ndec = 0
            fmt = '(f10.0)'
        else
            ndec = 3
            fmt = '(f10.3)'
        endif
        do 100 n=1,ndo
            do 110 m=1,nd-1
                ypos = ylen - real(m)*dz - (n-1)*0.5*ylen
                yval = dmin + m*dd
                call plot(x0,ypos,3)
                call plot(x0+0.10,ypos,2)
                write(ostr,fmt)yval
                ht = 0.05   
                call symbol(x0-11.0*ht,ypos-0.5*ht,ht,ostr,0.0,10)
  110       continue
  100   continue
        hht = 0.07
            nchar = 10
            call symbol(x0-9.0*ht,ylen/2.0 - nchar*0.5*hht,
     1          hht,'Depth (km )',90.0,nchar)
        return
        end
                
        subroutine gcmdln(nid,zmax)
c-----
c     parse the command line parameters
c-----
c     nid I*4 0 for velocity dispersion and model (default)
c             1 for gamma dispersion and Qb inv model
c-----
        integer*4 nid
        character names*20
        real zmax
        nmarg = mnmarg()
c-----
c     initialize defaults
c-----
        nid = 0
        zmax = -1.0
c-----
c     loop through command line arguments
c-----

        i = 0
 1000   continue
            i = i + 1
            if(i.gt.nmarg)go to 2000
            call mgtarg(i,names)
            if(names(1:2).eq.'-V' .or. names(1:2).eq.'-v')then
                nid = 0
            else if (names(1:2).eq.'-G' .or. names(1:2).eq.'-g')then
                nid = 1
            else if (names(1:2).eq.'-Z' .or. names(1:2).eq.'-z')then
                i = i + 1
                call mgtarg(i,names)
                read(names,'(bn,f10.0)')zmax
            endif
        go to 1000
 2000   continue
        return
        end
                
        subroutine vplt(xpos,ypos,n, x, y,xmin,xmax,dx,
     1          ymin,ymax,dy,iskp)
c-----
c     y   R*4 - array of parameters to be plotted along y-axis
c     x   R*4 - array of parameters to be plotted along x-axis
c     iskp    I*4 - 2 solid line
c               1 dashed line
c-----
c     plotting velocity model is odd since have to force the steps
c-----
        parameter(NL=200,NL2=NL+NL)
        real*4 x(*), y(*)
        dphl = 0.0
        dphu = y(1)
        do 100 i=2,n
            xx1 = xpos + (x(i-1) - xmin)*dx
            xx2 = xpos + (x(i  ) - xmin)*dx
            yy1 = ypos - (dphl - ymin)*dy
            yy2 = ypos - (dphu - ymin)*dy
            call plot(xx1,yy1,3)
            if(iskp.eq.2)then
                call plot(xx1,yy2,2)
                call plot(xx2,yy2,2)
            else
                call plotd(xx1,yy2,21,0.05)
                call plotd(xx2,yy2,21,0.05)
            endif
            dphl = dphu
            dphu = dphl + y(i)
  100   continue
        dphu = 1.10*dphl
        yy2 = ypos - (dphu - ymin)*dy
        if(iskp.eq.2)then
            call plot(xx2,yy2,2)
        else
            call plotd(xx2,yy2,21,0.05)
        endif
        return
        end

        subroutine dvscale(vmin,vmax)
c-----
c     rescale xval so that a plot will look nice
c     this is currently for velocity which can never be negative
c-----
c     vmin    R   - value to be scaled
c     vmax    R   - value to be scaled
c-----
c-----
c     first get size
c-----
        size = max(vmin,vmax)
        if(size.gt.1)then
            lpow = alog10(size)
            tmp = 10.0**lpow
            if(size.lt.tmp)then
                tmp = 0.1 * tmp
            endif
        else
            lpow = alog10(size) + 1
            tmp = 10.0**lpow
            if(size.lt.tmp)then
                tmp = 0.1 * tmp
            endif
        endif
        vmn = vmin
        do 1000 i=0,20
            xx = real(i) * 0.5 *tmp
            if(xx .lt. vmin)then
                vmn = xx 
            endif
 1000   continue
        vmx = vmax
        do 2000 i=20,0,-1
            xx = real(i) * 0.5 *tmp
            if(xx .gt. vmax)then
                vmx = xx 
            endif
 2000   continue
C       write(6,*)vmin,vmax,vmn,vmx,size,lpow,tmp
        if(vmx.gt.vmn)then
            vmin = vmn
            vmax = vmx
        endif
C       write(6,*)vmin,vmax,vmn,vmx,size,lpow,tmp

        return
        end
