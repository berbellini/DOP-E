        program srfphr96
c-----
c     CHANGES
c     21 FEB 2011 - add option for P/S inversion for lstinv
c-----
c       general routine to plot current model
c       and resolution kernels
c-----
        parameter(NL=200,NL2=NL+NL,NLP2=NL+2)
        common/ctrl/ d(NLP2),a(NL),b(NLP2),r(NL),rat(NL),
     $      v(NL2,NL2),qbinv(NLP2),qainv(NLP2),
     2      wc(NL2)
        logical wc
c-----
c       temporary data array
c-----
        real*4 darr(4), aarr(4), barr(4), qbarr(4) 
        real*4 od(NLP2), oa(NLP2), ob(NLP2), oqbinv(NLP2)
c-----
c       open data file
c-----
        open(2,file='tmpsrfi.02', access='sequential',form='formatted',
     1      status='unknown')
        rewind 2
c-----
c       get data set
c-----
        read(2,'(7i5,2e16.8)',end=9999,err=9999)
     1      iunit,mlyr,itype,invdep,mlw,mup,lstinv,dlam
c-----
c       iunit   R*4 - 0 km, 1 ft, 2 m
c       mlyr    I*4 - number of layers in model
c       itype   I*4
c       dlam    R*4 - current damping value
c       invdep  I*4 - 1 - inversion was for beta, Q(s) inverse or both
c                 0 - inversion was for layer thickness
c       mlw i*4
C       mup I*4 - lower and upper limits of Resolution matrix output
c       lstinv  I*4 - Current resolution kernels are for
c                 2 - inversion for S velocity
c                 3 - inversion for inverse Q(s) 
c                 4 - inversion for both S velocity and  inverse Q(s) 
c                 5 - layer thickness only
c                 6 - P and S velocity inversion in shallow96
c-----

        amin =  1.0e+37
        amax = -1.0e+37
        bmin =  1.0e+37
        bmax = -1.0e+37
        qbmin =  1.0e+37
        qbmax = -1.0e+37
        hmin =  1.0e+37
        hmax = -1.0e+37
        do 145 ii=mlw, mup
            if(ii.gt.mlyr)then
                i = ii - mlyr
            else
                i = ii
            endif
            if(invdep.eq.1)then
                read(2,'(i5,4e16.8/5x,5e16.8)',end=9999,err=9999)i,
     1              dum,d(i),a(i),b(i),
     2              qbinv(i),oa(i),ob(i),oqbinv(i)
                    od(i) = d(i)
            WRITE(0,*)i,a(i),b(i),qbinv(i),oa(i),ob(i),oqbinv(i)
            else
                read(2,'(i5,4e16.8/5x,2e16.8)',end=9999,err=9999)i,
     1              dum,d(i),a(i),b(i),
     2              qbinv(i),od(i)
                ob(i) = b(i)
                oa(i) = a(i)
                oqbinv(i) = qbinv(i)
            endif
            read(2,'(5e16.8)',end=9999,err=9999)(v(ii,j),j=mlw,mup)
c-----
c       beware of layer thickness change
c-----
            if(invdep.eq.1)then
                call amxmn(ob(i),bmax,bmin)
                call amxmn(oa(i),amax,amin)
                call amxmn(oqbinv(i),qbmax,qbmin)
            endif
                call amxmn(b(i),bmax,bmin)
                call amxmn(a(i),amax,amin)
                call amxmn(qbinv(i),qbmax,qbmin)
  145   continue
        close (2)
c-----
c       get total model thickness
c-----  
        dph = 0.0
        odph = 0.0
        do 146 i=1,mlyr-1
            dph = dph + d(i) 
            if(invdep.eq.0)then
                odph = odph + od(i) 
            endif
  146   continue
        if(invdep.eq.0)then
            if(odph .gt. dph)dph = odph
        endif
        dph = dph * 1.1
        darr(1) = 0.0
        darr(2) = dph
c-----
c       set up plot area
c-----
        xorg = 0.5
        yorg = 0.5
c-----
c       plot it all
c-----
        call pinitf('SRFPHR96.PLT')
c-----
c       annotate the figure
c-----
        call plot( xorg,yorg,-3)
        ylen = 6.0
        ht = 0.14
        if(lstinv.eq.6)then
              call symbol(0.0,ylen+7.0*ht,ht,
     1            'Velocity Model and Resolution ',0.0,30)
        else
              call symbol(0.0,ylen+7.0*ht,ht,
     1            'Velocity/Q Inverse Model and Resolution ',0.0,40)
        endif
        call symbol(0.0,ylen+5.5*ht,ht,
     1      ' (Red/dashed=current,Blue/solid=next)',0.0,37)
        call symbol(0.0,ylen+4.0*ht,ht,'Damping=',0.0,8)
        call number(999.0,999.0,ht,dlam,0.0,1003)
c-----
c       put in the velocity model
c-----
        xlen = 3.0
        ylen = 6.0
c-----
c       get axis scaling
c-----
        barr(1) = bmin
        barr(2) = bmax
        call dvscale(bmin,bmax)
        barr(3) = bmin
        barr(4) = (bmax - bmin)/xlen

        aarr(1) = amin
        aarr(2) = amax
        call dvscale(amin,amax)
        aarr(3) = amin
        aarr(4) = (amax - amin)/xlen

        qbarr(1) = 0.95*qbmin
        qbarr(2) = 1.05*qbmax
        call gscale(qbarr,xlen,2,1)
        qbmin = qbarr(3)
        qbmax = qbarr(3) + xlen * qbarr(4)

        yylen = 10.0
        call gscale(darr,yylen,2,1)
        dmin = 0.0
        dmax = darr(3) + yylen*darr(4)

        if(lstinv.eq.4)then
            call labz(dmin,dmax,10,0.0,ylen,ylen,1,dd,iunit)
        else if(lstinv.eq.6)then
            call labz(dmin,dmax,10,0.0,ylen,ylen,1,dd,iunit)
        else
            call labz(dmin,dmax,10,0.0,ylen,ylen,0,dd,iunit)
        endif

        call box(0.0,0.0,xlen,ylen)
        apos = ylen
        bpos = ylen
        qbpos = ylen
        if(lstinv.eq.4)then
            call newpen(1)
            call plot(0.0 ,ylen/2.0,3)
            call plot(xlen,ylen/2.0,2)
            call newpen(1)
            qbpos = ylen / 2.0
        else if(lstinv.eq.6)then
            call newpen(1)
            call plot(0.0 ,ylen/2.0,3)
            call plot(xlen,ylen/2.0,2)
            call newpen(1)
            apos = ylen / 2.0
        endif
        WRITE(0,*)amin,amax,bmin,bmax,qbmin,qbmax
        if(lstinv.ne.3)then
            if(iunit.eq.0)then
            call dolinx(0.0,ylen,xlen,bmax,bmin,
     1          0.05,.true.,.true.,.true.,9,'VS (km/s)')
            else if(iunit.eq.1)then
            call dolinx(0.0,ylen,xlen,bmax,bmin,
     1          0.05,.true.,.true.,.true.,9,'VS (ft/s)')
            else if(iunit.eq.2)then
            call dolinx(0.0,ylen,xlen,bmax,bmin,
     1          0.05,.true.,.true.,.true.,9,'VS (m/s)')
            endif
        endif
        if(lstinv.eq.3 .or. lstinv.eq.4)then
            call myaxis(0.0,0.0,'QBINV',-5,xlen,0.0,qbarr(3),qbarr(4))
        else if(lstinv.eq.6)then
            call dolinx(0.0,0.0,xlen,amax,amin,
     1          0.05,.false.,.false.,.true.,9,'VP (km/s)')
        endif
c-----
c       plot the current velocity and Qb models
c-----
        call newpen(2)
        if(lstinv.eq.2 .or. lstinv.eq.4 .or. lstinv.eq.5)then
            dx = 1.0/barr(4)
            call vplt(0.0,bpos,mlyr,ob,od,bmin,bmax,dx,
     1          dmin,dmax,dd,1)
        else if(lstinv.eq.6)then
            dx = 1.0/barr(4)
            call vplt(0.0,bpos,mlyr,ob,od,bmin,bmax,dx,
     1          dmin,dmax,dd,1)
            dx = 1.0/aarr(4)
            call vplt(0.0,apos,mlyr,oa,od,amin,amax,dx,
     1          dmin,dmax,dd,1)
        else
            dx = 1.0/qbarr(4)
            call vplt(0.0,qbpos,mlyr,oqbinv,od,qbmin,qbmax,dx,
     1          dmin,dmax,dd,1)
        endif
c-----
c       plot the future model if the current damping is used
c-----
        call newpen(4)
        if(lstinv.eq.2 .or. lstinv.eq.4 .or. lstinv.eq.5)then
            dx = 1.0/barr(4)
            call vplt(0.0,bpos,mlyr, b, d,bmin,bmax,dx,
     1          dmin,dmax,dd,0)
        endif
        if(lstinv.eq.3 .or. lstinv.eq.4)then
            dx = 1.0/qbarr(4)
            call vplt(0.0,qbpos,mlyr, qbinv, d,qbmin,qbmax,dx,
     1          dmin,dmax,dd,0)
        else if(lstinv.eq.6)then
            dx = 1.0/barr(4)
            call vplt(0.0,bpos,mlyr, b, d,bmin,bmax,dx,
     1          dmin,dmax,dd,0)
            dx = 1.0/aarr(4)
            call vplt(0.0,apos,mlyr, a, d,amin,amax,dx,
     1          dmin,dmax,dd,0)
        endif
        call newpen(1)

        call plot(-xorg,-yorg,-3)
c-----
c       put in the resolution kernels
c-----
        call plot( xorg+3.0, yorg,-3)
        xlen = 6.0
        ylen = 6.0
        ylen = 6.0
        ht = 0.14
        call box(0.0,0.0,xlen,ylen)
        ht = 0.10
        nchar=28
        call symbol(0.5*xlen-0.5*nchar*ht,ylen+2.0*ht,ht,
     1      'NORMALIZED RESOLUTION MATRIX',0.0,nchar)
        if(lstinv.eq.4)then
            call newpen(4)
            call plot(0.0 ,ylen/2.0,3)
            call plot(xlen,ylen/2.0,2)
            call plot(xlen/2.0,ylen,3)
            call plot(xlen/2.0, 0.0,2)
            call newpen(1)
        endif
c-----
c       normalize them all
c-----
        ntr = mup - mlw + 1
        yamp = ylen/(ntr + 1)
        if(yamp.lt.0.05)yamp=0.05
c BEWARE of DIVIDE by ZERO
        dx = xlen/real(ntr+1)
        dy = ylen/real(ntr-1)
c-----
c       put in matrix columns
c-----
        jjj = 0
        ndo = mup - mlw +1
        do 157 jj=mlw,mup
            jjj = jjj + 1
            xx = x0 + (jj - mlw+1)*dx
            if(jj.ge.mlw .and. jj.le.mup .and. ndo.le.20)then
                call newpen(4)
                call plot(xx,0.0,3)
                call plotd(xx,ylen,21,0.05)
                call newpen(1)
            endif
c-----
c           hack since 100 is the maximum number of layers
c-----
            if(jj.lt.100)then
                ii = jj/10
                i  = mod(jj,10) 
                ht = 0.07
                if(jj.gt.0)then
                call number(xx-0.5*ht,-3.5*ht,ht,real(ii),0.0,-1)
                endif
                call number(xx-0.5*ht,-5.0*ht,ht,real(i ),0.0,-1)
            else
                call number(xx-0.5*ht,-2.0*ht,ht,1.0,0.0,-1)
                call number(xx-0.5*ht,-3.5*ht,ht,0.0,0.0,-1)
                call number(xx-0.5*ht,-5.0*ht,ht,0.0,0.0,-1)
            endif
  157   continue
c-----
c       plot the resolution kernels
c       by COLUMN
c-----
        JJJ = 0
        if(ndo.le.20)then
            ht = 0.07
            ht2 = 0.07
        else
            ht = 0.05
            ht2 = 0.05
        endif
        do 155 JJ=mlw, mup
            JJJ = JJJ + 1
            x0 = 0.0
            rmax = -1.0e+37
c-----
c           get max so that we can plot relative amplitudes
c-----
            vm = 0.0
            III = 0
            call newpen(1)
            do 159 II=mlw,mup
                if(abs(v(ii,jj)).gt.vm)vm=abs(v(ii,jj))
                III = III + 1
                y0 = (mup - II  )*dy
            call number(xlen+0.15,y0-0.035,0.07,real(III),0.0,-1)
  159       continue
            call newpen(2)
            do 156 II=mlw,mup
                yval = v(ii,jj)
                y0 = (mup - II  )*dy
                yy = y0
                xx = x0 + (JJ - mlw+1)*dx + yamp*yval/vm
                if(II.eq.mlw)then
                    call plot(xx,yy,3)
                else
                    call plot(xx,yy,2)
                endif
                if(yval.gt.rmax)then
                    rmax = yval
                    xxsv = xx
                    yysv = yy
                endif
  156       continue
            call newpen(1)
            
            call symbol(xxsv,yysv,ht,char(1),0.0,-1)
            call number(xxsv-6*ht2,yysv+0.0*ht2,ht2,rmax,0.0,3)
            call plot(xlen-0.10,y0,3)
            call plot(xlen     ,y0,2)
  155   continue
CRBHc-----
CRBHc   plot the resolution kernels
CRBHc   by ROW
CRBHc-----
CRBH    iii = 0
CRBH    if(ndo.le.20)then
CRBH        ht = 0.07
CRBH        ht2 = 0.07
CRBH    else
CRBH        ht = 0.05
CRBH        ht2 = 0.05
CRBH    endif
CRBH    do 155 ii=mlw, mup
CRBH        iii = iii + 1
CRBH        x0 = 0.0
CRBH        y0 = (mup - ii + 1 )*dy
CRBH        rmax = -1.0e+37
CRBH        call newpen(2)
CRBHc-----
CRBHc       get max so that we can plot relative amplitudes
CRBHc-----
CRBH        vm = 0.0
CRBH        do 159 jj=mlw,mup
CRBH            if(abs(v(ii,jj)).gt.vm)vm=abs(v(ii,jj))
CRBH  159       continue
CRBH        do 156 jj=mlw,mup
CRBH            yval = v(ii,jj)
CRBH            yy = y0 + yamp * yval/vm
CRBH            xx = x0 + (jj - mlw)*dx
CRBH            if(jj.eq.mlw)then
CRBH                call plot(xx,yy,3)
CRBH            else
CRBH                call plot(xx,yy,2)
CRBH            endif
CRBH            if(yval.gt.rmax)then
CRBH                rmax = yval
CRBH                xxsv = xx
CRBH                yysv = yy
CRBH            endif
CRBH  156       continue
CRBH        call newpen(1)
CRBH        
CRBH        call symbol(xxsv,yysv,ht,char(1),0.0,-1)
CRBH        call number(xxsv+ht2,yysv+ht2,ht2,rmax,0.0,4)
CRBH        call plot(xlen-0.10,y0,3)
CRBH        call plot(xlen     ,y0,2)
CRBH        call number(xlen+0.15,y0-0.035,0.07,real(iii),0.0,-1)
CRBH  155   continue
            
        call newpen(1)
        call plot(-xorg-3.0,-yorg,-3)
        call newpen(1)
        call pend()
 9999   continue
        end

        subroutine box(xl,yl,xu,yu)
        call plot(xl,yl,3)
        call plot(xu,yl,2)
        call plot(xu,yu,2)
        call plot(xl,yu,2)
        call plot(xl,yl,2)
        return
        end

        subroutine labz(dmin,dmax,nd,x0,y0,ylen,id,ddd,iunit)
c-----
c       label the depth axis - however 
c-----
c       dmin    R*4 - minimum depth
c       dmax    R*4 - maximum depth
c       nd  I*4 - number of segments
c       x0  R*4
c       Y0  R*4 - upper left corner
c       ylen    R*4 - length of Y-axis
c       id  I*4 - 0 single
c                 1 double
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
        if(iunit.eq.0)then
            nchar = 10
            call symbol(x0-7.0*ht,ylen/2.0 + nchar*0.5*hht,
     1          hht,'DEPTH (km)',90.0,nchar)
        else if(iunit.eq.1)then
            nchar = 10
            call symbol(x0-7.0*ht,ylen/2.0 + nchar*0.5*hht,
     1          hht,'DEPTH (ft)',90.0,nchar)
        else if(iunit.eq.2)then
            nchar = 9
            call symbol(x0-7.0*ht,ylen/2.0 + nchar*0.5*hht,
     1          hht,'DEPTH (m)',90.0,nchar)
        endif
        return
        end
                
        subroutine vplt(xpos,ypos,n, x, y,xmin,xmax,dx,
     1          ymin,ymax,dy,iskp)
c-----
c       y   R*4 - array of parameters to be plotted along y-axis
c       x   R*4 - array of parameters to be plotted along x-axis
c       iskp    I*4 - 0 solid line
c                 1 dashed line
c-----
c       plotting velocity model is odd since have to force the steps
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
            if(iskp.eq.0)then
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
        if(iskp.eq.0)then
            call plot(xx2,yy2,2)
        else
            call plotd(xx2,yy2,21,0.05)
        endif
        return
        end

c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: AXIS                                                  c
c                                                                     c
c      COPYRIGHT (C)  1986 R. B. Herrmann                             c
c                                                                     c
c      Department of Earth and Atmospheric Sciences                   c
c      Saint Louis University                                         c
c      221 North Grand Boulevard                                      c
c      St. Louis, Missouri 63103                                      c
c      U. S. A.                                                       c
c                                                                     c
c---------------------------------------------------------------------c
        subroutine myaxis(xpage,ypage,ititl,nchar,axlen,angle,
     1          firstv,deltav)
        real HTNUM,SYMHT,TICHT
       parameter (HTNUM=0.07,SYMHT=0.10,TICHT=0.07)
c-----
c       xpage,ypage  coordinates of starting point of axis in
c                    inches relative to current origin
c       ititl   axis title
c       nchar   number of characters in title
c               >0 for tic marks, numbers and title on
c                  counterclockwise side
c               <0 for tic marks, numbers and title on
c                  clockwise side
c       axlen   floating point axis length in inches
c       angle   angle (degrees) that x axis makes with
c               horizontal x-direction
c       firstv  scale value at (xpage,ypage)
c       deltav  change in slace between tic marks per unit
c               length
c               <0 value at xpage,ypage is a max, and values
c               decrease along the axis
c-----
c-----
c       we do not worry about the vagaries of storage of characters
c       since the address is passed on down to subroutine symbol
c-----
        character*(*) ititl
        a=1.0
        kn=nchar
        if(kn.lt.0)then
                a= -a
                kn= -kn
        endif
c-----
c       if deltav is too large invoke scientific notation
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
c       first put in numbers and title
c       adjust offset for numbers
c-----
        dx = - HTNUM
        dy = 1.5*a*HTNUM - 0.5*HTNUM
c-----
c       find initial position given rotation
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
c       halfway down axis put in title
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
c       now put in tic marks
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

        subroutine dvscale(vmin,vmax)
c-----
c       rescale xval so that a plot will look nice
c       this is currently for velocity which can never be negative
c-----
c       vmin    R   - value to be scaled
c       vmax    R   - value to be scaled
c-----
c-----
c       first get size
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

        subroutine amxmn(x,xmax,xmin)
c----
c       given x, update the xmax and xmin
c-----
            if(x.gt.xmax)xmax=x
            if(x.lt.xmin)xmin=x
            return
       end
