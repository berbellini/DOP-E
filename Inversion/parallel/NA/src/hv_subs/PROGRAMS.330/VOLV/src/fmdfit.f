        program fmdfit
c----------------------------------------------------------------------c
c                                                                      c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                 c
c      VOLUME V                                                        c
c                                                                      c
c      PROGRAM: fmdfit                                                 c
c                                                                      c
c      COPYRIGHT 2002                                                  c
c      R. B. Herrmann                                                  c
c      Department of Earth and Atmospheric Sciences                    c
c      Saint Louis University                                          c
c      221 North Grand Boulevard                                       c
c      St. Louis, Missouri 63103                                       c
c      U. S. A.                                                        c
c                                                                      c
c----------------------------------------------------------------------c
c Revision history:
c       31 JUL 2002 - programs created
c       21 APR 2004 - corrected typo in subroutine usage
c       11 OCT 2004 - Added optional focal mechanism plot
c               note radius chosen internally works with
c               fmdfit -HMN 0 -HMX 26 -MECH
c       09 JAN 2005 - extra long lines split
c
c       Plot focal mechanism grid search best fits results from 
c       srdisp96 and ..... as a function of given depth
c----------------------------------------------------------------------c
        implicit none
        integer LER, LIN, LOT
        parameter (LIN=5,LOT=6,LER=0)                          
        


        real x0, xlen, y0, ylen
        
        real hmn, hmx 
        real fmn, fmx
        character instr*170
        integer ls, lsep
        real h, stk, dip, rak, rr, rl, xmwr, xmwl, best
        logical doit
        integer lgstr
        real x6, x7
        logical domech

        call gcmdln(hmn,hmx,
     1      x0,y0,xlen,ylen,domech)

        fmn = 0.0
        fmx = 1.0
        
        call pinitf('FMDFIT.PLT')
c-----
c       draw a box
c-----
        call gbox(x0,y0,x0+xlen,y0+ylen)
c-----
c       put in the axes
c-----
        call dolinx(x0,y0,xlen,hmn,hmx,
     1      0.10,.true.,.false.,.true., 10, 'Depth (km)')
        call dolinx(x0,y0+ylen,xlen,hmn,hmx,
     1      0.10,.false.,.false.,.false., 10, 'Depth (km)')
        call doliny(x0,y0,ylen,fmn,fmx,
     1      0.10,.false.,.true.,.true.,3,'Fit')
        call doliny(x0+xlen,y0,ylen,fmn,fmx,
     1      0.10,.true.,.true.,.false.,3,'Fit')
c-----
c       read data and plot only within the search bounds
c-----
 1000   continue
        read(5,'(a)',end=9000)instr
c-----
c       check format
c-----
            ls = lgstr(instr)
            doit = .false.
            lsep = index(instr,'SRFGRD96')
            if(lsep.gt.0 .and.ls.gt.8)then
                read(instr(lsep+8:ls),*)h,stk,dip,rak,rr,rl,
     1          xmwr,xmwl,best
                doit = .true.
            endif
            lsep = index(instr,'WVFGRD96')
            if(lsep.gt.0 .and.ls.gt.8)then
                read(instr(lsep+8:ls),*)h,stk,dip,rak,xmwr,best
                doit = .true.
            endif
            lsep = index(instr,'WVFMTD96')
            if(lsep.gt.0 .and.ls.gt.8)then
                read(instr(lsep+8:ls),*)h,stk,dip,rak,xmwr,x6,x7,best
                doit = .true.
            endif
            lsep = index(instr,'WVFMT96')
            if(lsep.gt.0 .and.ls.gt.7)then
                read(instr(lsep+7:ls),*)h,stk,dip,rak,xmwr,x6,x7,best
                doit = .true.
            endif
c-----
c       if have valid entry, plot it
c-----
            if(doit .and. best.gt.0.0)then
                call plotit(h,best,fmn,fmx,hmn,hmx,
     1              x0,y0,xlen,ylen,stk,dip,rak,domech)
            endif
        go to 1000
 9000   continue
        call pend()
        end

        subroutine plotit(h,best,fmn,fmx,hmn,hmx,
     1     x0,y0,xlen,ylen,stk,dip,rak,domech)
c-----
c       h   R   - depth
c       best    R   - best value [-1,1]
c       fmn R   - minimum fit value for y-axis
c       fmx R   - maximum fit value for y-axis
c       hmn R   - minimum depth value for x-axis
c       hmx R   - maximum depth value for x-axis
c       x0  R   - x-coordinate of lower left corner
c       y0  R   - y-coordinate of lower left corner
c       xlen    R   - length of x-axis
c       ylen    R   - length of y-axis
c       stk R   - strike of focal mechanism
c       dip R   - dip of focal mechanism
c       rak R   - rake of focal mechanism
c       domech  L   - plot focal mechanism
c-----
        implicit none
        real h,best,hmn,hmx,fmn,fmx
        real x0,y0,xlen,ylen
        logical domech
        real xx, yy
        integer jpen
        real dip,stk,rak,xc,yc,radc
        integer kolor, ityp
c-----
c       map best [0,1] to colors [BLUE,RED] = [1100, 1000]
        jpen = 1100 - 100.0*best
        xx = x0 + xlen * (h - hmn)/(hmx - hmn)
        yy = y0 + ylen * (best - fmn)/(fmx - fmn)
        call newpen(jpen)
        call fillit('CI',0.05,xx,yy)
        call newpen(1)
        call curvit('CI',0.05,xx,yy)
        if (domech)then
           xc = xx
           yc = yy+0.30
            kolor = 1
            radc = 0.1
            ityp = 1
            call dofmplt(dip,stk,rak,xc,yc,kolor,radc,ityp)
        endif
        return
        end



        subroutine gcmdln(hmn,hmx,
     1      x0,y0,xlen,ylen,domech)
c-----
c       parse command line arguments
c
c       requires subroutine mgtarg() and function mnmarg()
c
c-----
c       hmn,hmn         R*4 - depth range search parameters
c-----
        implicit none
        integer LER, LIN, LOT
        parameter (LER=0, LIN=5, LOT=6)
c-----
c       parameters for grid search
c-----
        real hmn, hmx
        real x0, y0, xlen, ylen
        logical domech

        character*25 names
        integer*4 mnmarg
        integer nmarg
        integer i
c-----
c       initialize variables
c-----
        hmn = 0
        hmx = 40
        
        x0 = 1.0
        y0 = 1.0
        xlen = 8.0
        ylen = 6.0
        domech = .false.
c-----
c       process command line arguments
c-----
        nmarg=mnmarg()
        i = 0
   11   i = i + 1
        if(i.gt.nmarg)go to 13
            call mgtarg(i,names)
            if(names(1:4).eq.'-HMN')then
                i=i+1
                call mgtarg(i,names)
                call chtofp(names,hmn)
            else if(names(1:4).eq.'-HMX')then
                i=i+1
                call mgtarg(i,names)
                call chtofp(names,hmx)
            else if(names(1:2).eq.'-M')then
                domech = .true.
            else if(names(1:2).eq.'-?')then
                call usage(' ')
            else if(names(1:2).eq.'-h')then
                call usage(' ')
            endif
        goto 11
   13   continue
        return
        end
        
        subroutine usage(ostr)
        implicit none
        character ostr*(*)
        integer LER, LIN, LOT
        parameter (LER=0, LIN=5, LOT=6)
        integer ls
        integer lgstr
        ls = lgstr(ostr)
        if(ostr.ne.' ') write(LER,*)ostr(1:ls)
        write(LER,*)'fmdfit -HMN hmn -HMX hmx -M '
        write(LER,*)
     1  ' -HMN hmn    (default   0)  Minimum depth'
        write(LER,*)
     1  ' -HMX hmx    (default  40)  Maximum depth'
        write(LER,*)
     1  ' -M          (default no) plot mechanism '
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
        implicit none
        character*(*) str
        real*4 fout
        integer*4 lgstr
        integer i, l
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
            read(str,'(bn,e20.13)')fout
        else
            read(str,'(bn,f20.0)') fout
        endif
        return
        end

        subroutine dofmplt(dip,stk,rake,xc,yc,kolor,radc,ityp)
c-----
c       ityp    I*4 - 0 plot mechanism
c-----
        real dip,stk,rake,xc,yc,radc
        integer kolor, ityp
        integer NTR
        parameter (NTR=90 )
        real*4 val(NTR,2), tval(NTR)
        real*4 mom(3,3), gamm(3), phi(3), theta(3)
        common/pltcon/x0,y0,rad,eqarea
        logical eqarea
        common/pltmap/tr0,dtr,pl0,dpl

        parameter (mp=3,np=3)
        real*8 xmt(np,np),ev(np),ev1(np)
        real*8 z(3,3)

        character title*80
        character subtit*80

        logical ifm
        logical fmamp
        integer*4 fmfill
        logical fmtype
        logical verbos
        logical nullit
        logical hastit
        logical hassub
        logical pltmch
        logical hemis
        
        real*4 titsiz

        character sym(5)*3
        data sym/'  P',' SV',' SH','POL',' S '/
        data mom/9*0.0/

        x0 = xc
        y0 = yc
        rad = radc

        eqarea = .true.
        itype = 1
        fmfill = 1
        ifm = .false.
        fmamp = .false.
        fmtype = .false.
        verbos = .false.
        nullit = .true.
        title = ' '
        hastit = .false.
        subtit = ' '
        hassub = .false.
        pltmch = .true.
        titsiz = 0.05
        hemis = .true.

        call tensor(stk,dip,rake,1.0,mom)
c-----
c       force symmetry of moment tensor - this is valid for
c       double couples. 
c       If the field of a single couple is desired, do not force
c       this symmetry
c-----
        mom(3,1) = mom(1,3)
        mom(3,2) = mom(2,3)
        mom(2,1) = mom(1,2)

c-----
c       compute eigenvalues and eigenvectors of moment tensor matrix
c       Get index of largest and smallest eigenvalues
c-----
        do 100 i=1,3
            do 110 j=1,3
                xmt(i,j) = dble(mom(i,j))
  110       continue
  100   continue
        call tred2(np,np,xmt,ev,ev1,z)
        call tql2(np,np,ev,ev1,z,ierr)
c-----
c       get the largest and smallest eigenvectors
c-----
        elg = -1.0e+38
        esm =  1.0e+38
        ilg = 1
        ism = 1
C       write(6,*)' EIGENVALUE, AND EIGENVECTOR OF M(i,j)'
C    2  format(' ',e11.4,'  (',e11.4,',',e11.4,',',e11.4,')')
        do 120 j=1,3
C           write(6,2)ev(j),(xmt(i,j),i=1,3)
            if(ev(j).gt.elg)then
                elg = ev(j)
                ilg = j
            endif
            if(ev(j).lt.esm)then
                esm = ev(j)
                ism = j
            endif
  120   continue
c-----
c       if(ityp.ne.0)just output the mechanism as the P or T axis, but
c       temper everything with the apparent dip.
c-----

c-----
c       sample the radiation pattern on focal sphere, 
c       and save value of amplitude on the focal sphere
c-----
        wid = 0.02*rad
        if(wid .lt. 0.005)wid = 0.005
        if(fmfill .gt. 0)then
            npl = 3*(rad/wid )
        else
            npl = 20
        endif
        dtr = 360.0 / (NTR - 1)
        dpl = 90.0 / (npl - 1)
        tr0 = 0.0
        pl0 = 0.0
        open(4,
     1  status='scratch',form='unformatted',access='sequential')
        rewind 4
        vmax = -1.0e+38
        vmin =  1.0e+38
        if(pltmch)then
        do 1000 ip = 1, npl
            if(fmfill .gt. 0)then
                r = ip*wid/3.0
                if(eqarea)then
                    pl = 2.0*asin(0.707*r/rad)
                else
                    pl = 2.0*atan(r/rad)
                endif
            else
                pl = pl0 + (ip-1)*dpl
                pl = pl * 3.1415927/180.0
            endif
            do 1100 it = 1, NTR
                tr = tr0 + (it-1)*dtr
c-----
c           convert to spherical coordinates on lower unit hemisphere
c-----
                tr = tr * 3.1415927/180.0
                if(hemis)then
                call getmtn(tr,pl,gamm,phi,theta,mom,sump,sumsv,sumsh)
                else
                call getmtn(tr,-pl,gamm,phi,theta,mom,sump,sumsv,sumsh)
                endif
c-----
c       get  P-wave pattern for itype = 1 
c       get SV-wave pattern for itype = 2 
c       get SH-wave pattern for itype = 3 
c       get S  polarization for itype = 4
c       get S  wave pattern for itype = 5
c-----
            if(itype.eq.1)then
                sum = sump
            else if(itype.eq.2)then
                sum = sumsv
            else if(itype.eq.3)then
                sum = sumsh
            else if(itype.eq.4)then
                sum = sump
                if(mod(it,5).eq.1 .and. mod(ip,3).eq.0)then
                if(sumsh.ne.0.0 .or. sumsv.ne.0.0)then
                    pol = atan2(sumsh,sumsv)
                    call pltpol(ip,it,pol,tr)
                endif
                endif
            else if(itype.eq.5)then
                sum = sqrt(sumsh**2 + sumsv**2)
            endif

            tval (it) = sum
            if(sum.gt.vmax)vmax = sum
            if(sum.lt.vmin)vmin = sum
 1100       continue
            write(4)(tval(it),it=1,NTR)
 1000   continue
        endif
c-----
c       define symbol height
c-----
        ht = 0.14*rad
        hht = titsiz
        call gwidth(0.5*wid)
c-----
c       put on the focal mechanism type
c-----
        call gfont(3)
c-----
c       plot focal circle here, so that it overwrites the title
c       string
c-----
        call gwidth(wid)
        call circle(rad,x0,y0,nullit)
c-----
c       go through focal mechanism options - they are not 
c       necessarily all exclusive
c-----
c-----
c       plot in the regions of positive functional value
c-----
        cval = 0.0
        if(itype.ge.1 .and. (itype.le.3 .or. itype.eq.5)
     1       .and. .not. ifm)then
            if(abs(vmin).gt.vmax)vmax = abs(vmin)
            if(fmfill.gt.0)then
                call newpen(kolor)
                call conshd(val,tval,1,NTR,npl,cval,NTR)
            else if(fmfill.lt.0)then
                call consh1(val,tval,1,NTR,npl,vmax,NTR)
            endif
        endif
c-----
c       
c       Otherwise just plot the zero contour
c-----
            cval = 0.0
            call newpen(1)
            call contur(val,tval,2,NTR,npl,cval,NTR)
            call circle(rad,x0,y0,.false.)
        close(4)
        return
        end

        subroutine getmtn(tr,pl,gamm,phi,theta,mom,sump,sumsv,sumsh)
        real gamm(3), theta(3), phi(3), mom(3,3)
                cost = cos(tr)
                sint = sin(tr)
                cosp = cos(pl)
                sinp = sin(pl)
                gamm(1) = cost*cosp
                gamm(2) = sint*cosp
                gamm(3) = sinp
                theta(1) = sinp*cost
                theta(2) = sinp*sint
                theta(3) = -cosp
                phi(1) = -sint
                phi(2) = cost
                phi(3) = 0.0
                sump = 0.0
                sumsv = 0.0
                sumsh = 0.0
                do 1200 i=1,3
                    do 1201 j=1,3
                        sump=sump+gamm(j)*gamm(i)*
     1                      mom(i,j)
                        sumsv=sumsv+theta(i)*gamm(j)*
     1                      mom(i,j)
                        sumsh=sumsh+phi(i)*gamm(j)*
     1                      mom(i,j)
 1201               continue
 1200           continue
        return
        end

        subroutine contur(f,g,i1,imax,jmax,val,npts)
        integer LIN, LOT, LER
        parameter (LIN=5, LER=0, LOT=6)
        integer  NTR
        parameter (NTR=90)
        real*4 f(NTR,2),g(NTR)
c-----
c       contouring program
c           plot a contour of level val from array f(imax,jmax)
c-----
        do 1000 jj = 2,jmax
            call getf(f,g,jj,i1,imax,npts)
        do 1000 i = i1,imax
c-------
c
c       F(i-1,j-1)      F(i-1,j)
c
c       F(i,j-1)    F(i,j)
c
c-----
c-----
c       to save memory and to permit this to
c       work in a non-virtual machine, the  function f(time,freq)
c       are stored sequentially on UNIT 01. At any one time we
c       need only access 2 columns, so we use j = 2
c-----
            j = 2
c-----
c       get center value, assuming Laplacian = 0
c-----
            cval = (F(i,j)+F(i-1,j)+F(i-1,j-1)+F(i,j-1))/4.0
c-----
c       now sequentially process the four triangles
c-----
            call lookt(i,jj,i  ,j  ,i-1,j  ,val,cval,F)
            call lookt(i,jj,i-1,j  ,i-1,j-1,val,cval,F)
            call lookt(i,jj,i-1,j-1,i  ,j-1,val,cval,F)
            call lookt(i,jj,i  ,j-1,i  ,j  ,val,cval,F)
 1000   continue
        return
        end

        subroutine lookt(ii,jj,i1 ,j1 ,i2,j2 ,val,cval,F)
        integer  NTR
        parameter (NTR=90)
        real*4 f(NTR,2)
        real*4 xc(3), yc(3)
        logical inarea
        logical interp
        nc = 0
c-----
c       sequentially go through corners of triangle
c-----
        inarea = interp(F(i1,j1),F(i2,j2),val,pos)
        if(inarea)then
            nc = nc + 1
            xc(nc) =  (i1 + (i2-i1)*pos)
            yc(nc) =  (j1 + (j2-j1)*pos) + (jj -2)
        endif
        inarea = interp(cval,F(i2,j2),val,pos)
        if(inarea)then
            nc = nc + 1
            xc(nc) =  (ii-0.5 + (i2-(ii-0.5))*pos)
            yc(nc) =  (1.5 + (j2-1.5)*pos) + (jj -2)
        endif
        inarea = interp(F(i1,j1),cval,val,pos)
        if(inarea)then
            nc = nc + 1
            xc(nc) =  (i1 + ((ii-0.5)-i1)*pos)
            yc(nc) =  (j1 + (1.5-j1)*pos) + (jj -2)
        endif
        ipen = 3
        iret = 1
        if(nc.gt.1)then
            do 100 i=1,nc
            call pplot(xc(i),yc(i),x,y)
            if(iret.gt.0)then
                call plot(x,y,ipen)
                ipen = 2
            endif
  100       continue
        endif
        return
        end

        function interp(y0,y1,val,x)
        logical interp
        real y1,y0,val,x
        denom = y1-y0
        if(denom.eq.0.0 .and. val.ne.y0)then
            interp = .false.
        elseif(denom.eq.0.0 .and. val.eq.y0)then
            interp = .true.
            x = 0.0
        else
            x = (val-y0)/denom
            if(x.lt.0.0 .or. x.gt.1.0)then
                interp = .false.
            else
                interp = .true.
            endif
        endif
        return
        end

        subroutine getf(f,g,jj,n1,n,npts)
      integer LIN, LOT, LER
      parameter (LIN=5, LER=0, LOT=6)
      integer  NTR
      parameter (NTR=90)
        real*4 f(NTR,2),g(NTR)
c-----
c       get the traces envelopes from UNIT 01
c       then update the array f(4096,2) so that column 1 refers
c       always to the leftmost column, even after updating
c-----
c-----
c       use special care for the first entry, jj = 2
c-----
        if(jj.eq.2)then
            rewind 4
            read(4)(g(j),j=1,npts)
            DO 5009 KK = 1,npts
                f(KK,2) = g(KK)
 5009       continue
        endif
c-----
c       normal processing, no check for EOF on read
c
c       1. read in new column
c       2. determine maximum amplitude of new column
c       3. normalize new column, put in f(KK,2)
c          move old f(KK,2) to f(KK,1)
c-----
            read(4)(g(j),j=1,npts)
            DO 5011 KK = 1,npts
                temp = f(KK,2)
                f(KK,2) = g(KK)
                f(KK,1) = temp
 5011       continue
        return
        end

        subroutine pplot(xi,yi,xx,yy)
        real*4 xx,yy
        common/pltcon/x0,y0,rad,eqarea
        logical eqarea
        common/pltmap/tr0,dtr,pl0,dpl
c-----
c       This routine converts an xi, yi coordinate to equivalent value
c       of trend and plunge, which are then projected onto a plane
c       using an equal area or stereographic projection
c-----
c       xi  R*4 -
c       yi  R*4 -
c       tr0 R*4 - value of trend corresponding to index xi = 1.0
c       dtr R*4 - trend increment
c       pl0 R*4 - value of plunge corresponding to index yi = 1.0
c       dpl R*4 - plunge increment
c       xx  R*4 - returned plot coordinate
c       yy  R*4 - returned plot coordinate
c       rad R*4 - radius of sphere
c       x0  R*4 - plot origin offset in x
c       y0  R*4 - plot origin offset in y
c       eqarea  L   - .true. equal area projection
c                   - .false. stereographic projection
c-----
        tr = tr0 + (xi - 1.0)*dtr
        pl = pl0 + (yi - 1.0)*dpl
        tr = tr * 3.1415927/180.0
        pl = pl * 3.1415927/180.0
c-----
c       convert to spherical coordinates on lower unit hemisphere
c-----
        x = cos(tr)*cos(pl)
        y = sin(tr)*cos(pl)
        z = sin(pl)
c-----
c       obtain angle between downward Z axis and vector
c       This is really just pi/2 - pl, but we will compute here
c       because there may be times when just (x,y,z) are given
c-----
        theta = atan2(sqrt(x**2 + y**2), abs(z))
c-----
c       get projection of the point onto the z=0 plane, determining
c       its distance from the origin
c-----
        if(eqarea)then
            r = sqrt(2.0)* sin(theta/2.0)
        else
            r = sin(theta/2.0)/cos(theta/2.0)
        endif
        if(x .ne. 0.0 .and. y .ne. 0.0)then
            alpha=atan2(y,x)
        else
            alpha = 0.0
        endif
c-----
c       get plotting point, noting that plotting system (xx,yy) is with
c       positive xx to right and yy up
c       and the fact that the trend =0 is equal to yy direction
c-----
        yy = y0 + rad * r * cos(alpha)
        xx = x0 + rad * r * sin(alpha)
        return
        end

        subroutine consh1(f,g,i1,imax,jmax,vmax,npts)
c-----
c       plot polarity information as + or - signs, with the amplitude
c       keyed to the amplitude on the focal sphere
c-----
        integer LIN, LOT, LER
        parameter (LIN=5, LER=0, LOT=6)
        integer  NTR
        parameter (NTR=90)
        real*4 f(NTR,2),g(NTR)
        common/pltcon/x0,y0,rad,eqarea
        logical eqarea


        do 1000 jj = 2,jmax
            call getf(f,g,jj,i1,imax,npts)
c-----
c       search for regions of f greater than cval
c       If this condition is met, interpolate and plot
c       We must be careful since the f(,) is a 2-D array, and
c       the first and last j indices do not wrap around, thus
c       we must handle jj=2 case carefully
c-----
            if(mod(jj,2).eq.0)then
            do 1100 ii=1,imax,3
            xi = ii
            xj = jj
            call pplot(xi,xj,xx,yy)
            cal = f(ii,1)
            ssize = 0.05 * rad * abs(f(ii,1))/vmax
            dx = ssize/2.0
            dy = ssize/2.0
c-----
c       Make a plus sign
c_____
            if(cal.gt.0.0)then
                call plot(xx+dx,yy,3)
                call plot(xx-dx,yy,2)
                call plot(xx,yy+dy,3)
                call plot(xx,yy-dy,2)
                call plot(xx,yy,3)
c-----
c       Make a minus sign
c-----
            elseif(cal.lt.0.0)then
                call plot(xx+dx,yy,3)
                call plot(xx-dx,yy,2)
            endif

 1100   continue
        endif
 1000   continue
        return
        end

        subroutine conshd(f,g,i1,imax,jmax,cval,npts)
        integer LIN, LOT, LER
        parameter (LIN=5, LER=0, LOT=6)
        integer  NTR
        parameter (NTR=90)
        real*4 f(NTR,2),g(NTR)
        do 1000 jj = 2,jmax
            call getf(f,g,jj,i1,imax,npts)
c-----
c       search for regions of f greater than cval
c       If this condition is met, interpolate and plot
c       We must be careful since the f(,) is a 2-D array, and
c       the first and last j indices do not wrap around, thus
c       we must handle jj=2 case carefully
c-----
            if(jj.eq.2)then
                call pltpos(f,1,i1,jj-1,imax,cval)
            endif
            call pltpos(f,2,i1,jj,imax,cval)
 1000   continue
        return
        end

        subroutine pltpos(f,j,i1,jj,imax,cval)
        integer  NTR
        parameter (NTR=90)
        real*4 f(NTR,2)
        logical interp
        yi = jj
        xi = i1
        ipen = 3
        call pplot(xi,yi,xold,yold)
        call plot(xold,yold,ipen)
        vold = f(i1,j)
        do 1100 i = i1+1,imax
            vnew = f(i,j)
            xi = i
            if(vold.ge.cval .and. vnew.ge.cval)then
                call pplot(xi,yi,xx,yy)
                call plot(xold,yold,3)
                call plot(xx,yy,2)
            else if(vold.lt. cval .and. vnew.le.cval)then
                call pplot(xi,yi,xx,yy)
                call plot(xx,yy,3)
            else if(interp(vold,vnew,cval,x))then
                call pplot(xi+x-1.0,yi,xx,yy)
                if(vold.ge.cval)then
                    call plot(xold,yold,3)
                    call plot(xx,yy,2)
                    call plot(xx,yy,3)
                else
                    call plot(xx,yy,3)
                    call pplot(xi,yi,xx,yy)
                    call plot(xx,yy,2)
                endif
            endif
            xold = xx
            yold = yy
            vold = vnew
 1100   continue
        return
        end

        subroutine gttrpl(ev,xmt1,xmt2,xmt3,hemis,tr,pl)
        real*8 ev, xmt1, xmt2, xmt3
        logical hemis
        
        common/pltmap/tr0,dtr,pl0,dpl
        common/pltcon/x0,y0,rad,eqarea
        logical eqarea
        real*8 r
c-----
c       plot pole on focal sphere
c-----
c       convert to trend and plunge, lower hemisphere
c-----
        degrad = 3.1415927/180.0
        if(xmt3 .lt. 0.0d+00)then
            xmt1 = - xmt1
            xmt2 = - xmt2
            xmt3 = - xmt3
        endif
        r = xmt1**2 + xmt2**2 + xmt3**2
        r = dsqrt(r)
        xmt1 = xmt1/r
        xmt2 = xmt2/r
        xmt3 = xmt3/r
        pl = 90.0 - sngl(dacos(xmt3))/degrad
        r = xmt1**2 + xmt2**2 
        r = dsqrt(r)
        if(r.le.1.0d-06)then
            tr = 0.0
        else
            tr = sngl(datan2(xmt2,xmt1))/degrad
        endif
        if(.not.hemis)then
            pl =  - pl
            tr = 180 + tr
        endif
        return
        end


        subroutine plpole(ev,xmt1,xmt2,xmt3,sym,hemis)
        real*8 ev, xmt1, xmt2, xmt3
        character sym*1
        logical hemis
        
        common/pltmap/tr0,dtr,pl0,dpl
        common/pltcon/x0,y0,rad,eqarea
        logical eqarea
        real*8 r
c-----
c       plot pole on focal sphere
c-----
c       convert to trend and plunge, lower hemisphere
c-----
        degrad = 3.1415927/180.0
        if(xmt3 .lt. 0.0d+00)then
            xmt1 = - xmt1
            xmt2 = - xmt2
            xmt3 = - xmt3
        endif
        r = xmt1**2 + xmt2**2 + xmt3**2
        r = dsqrt(r)
        xmt1 = xmt1/r
        xmt2 = xmt2/r
        xmt3 = xmt3/r
        pl = 90.0 - sngl(dacos(xmt3))/degrad
        r = xmt1**2 + xmt2**2 
        r = dsqrt(r)
        if(r.le.1.0d-06)then
            tr = 0.0
        else
            tr = sngl(datan2(xmt2,xmt1))/degrad
        endif
        if(.not.hemis)then
            pl =  - pl
            tr = 180 + tr
        endif
        xi = (tr - tr0)/dtr + 1.0
        yi = (pl - pl0)/dpl + 1.0
        call pplot(xi,yi,xx,yy)
        call symfil(xx,yy,0.05*rad)
        ht = 0.12*rad
        call symbol(xx+0.6*ht,yy+ht,ht,sym,0.0,1)
        return
        end

        subroutine tensor(stk,dip,slip,xmom,m)
c-----
c       calculate moment tensor for a double couple mechanism
c-----
c       stk - R*4   - strike, measured clockwise from north
c                 when looking down, e.g., east = 90
c       dip - R*4   - dip, measured from horizontal when looking
c                 in direction of strike, fault dips down to right
c                 0 <= dip <= 90
c       slip    - R*4   - slip. Measured from horizontal with respect
c                 to strike. If -180 < slip <= 180 is taken as
c                 the convention, then a negative slip, implies 
c                 that the P-wave first motion in the center of the
c                 focal sphere is negative (e.g., dilatation)
c       xmom    - R*4   - moment value
c       m(3,3)  - R*4   - moment tensor
c-----
        integer LOT
        parameter (LOT=6)
        real stk, dip, slip, xmom, m(3,3)
        real degrad, tol, dp, st, sl
        real sinp, cosp, sint, cost, sinp2, cosp2, sinlp, coslp
        real sint2, cost2
        
            degrad=0.0174532925
            tol = 1.0e-7
            dp = degrad*dip
            st = degrad*stk
            sl = degrad*slip
            sinp=sin(dp)
            cosp=cos(dp)
            sinlp=sin(sl)
            coslp=cos(sl)
            sint=sin(st)
            cost=cos(st)
            sinp2=sin(2.*dp)
            cosp2=cos(2.*dp)
            sint2=sin(2.*st)
            cost2=cos(2.*st)
            m(1,1)=(-sinp*coslp*sint2-sinp2*sinlp*sint*sint)*xmom
            m(2,2)=(sinp*coslp*sint2-sinp2*sinlp*cost*cost)*xmom
            m(3,3)=(sinp2*sinlp)*xmom
            m(1,2)=(sinp*coslp*cost2+0.5*sinp2*sinlp*sint2)*xmom
            m(1,3)=(-cosp*coslp*cost-cosp2*sinlp*sint)*xmom
            m(2,3)=(-cosp*coslp*sint+cosp2*sinlp*cost)*xmom
            m(2,1) = m(1,2)
            m(3,1) = m(1,3)
            m(3,2) = m(2,3)
c-----
c       clean up small values
c-----
        xmax= -1.0e+38
        do 4 i=1,3
            do 5 j=1,3
                if(abs(m(i,j)).gt.xmax)xmax = abs(m(i,j))
    5       continue
    4   continue
        thresh = tol * xmax
        do 6 i=1,3
            do 7 j=1,3
                if(abs(m(i,j)).lt.thresh) m(i,j) = 0.0
    7       continue
    6   continue

c-----
c       write out the information
c-----
C       write(LOT,*)' DIP       =',dip
C       write(LOT,*)' SLIP      =',slip
C       write(LOT,*)' STRIKE    =',stk
C       write(LOT,*)' MOMENT    =',xmom
C       write(LOT,*)' MOMENT TENSOR:'
C       write(LOT,*)' M(1,1)    =',m(1,1)
C       write(LOT,*)' M(1,2)    =',m(1,2)
C       write(LOT,*)' M(1,3)    =',m(1,3)
C       write(LOT,*)' M(2,1)    =',m(2,1)
C       write(LOT,*)' M(2,2)    =',m(2,2)
C       write(LOT,*)' M(2,3)    =',m(2,3)
C       write(LOT,*)' M(3,1)    =',m(3,1)
C       write(LOT,*)' M(3,2)    =',m(3,2)
C       write(LOT,*)' M(3,3)    =',m(3,3)
        return
        end

        subroutine pltpol(ip,it,pol,tr)
c-----
c       plot polarization angle
c-----
        common/pltcon/x0,y0,rad,eqarea
        logical eqarea
            xi = it
            yi = ip
            call pplot(xi,yi,xx,yy)
            ang = tr + pol 
            c = cos(ang)
            s = sin(ang)
            xxx = 0.05*rad*s
            yyy = 0.05*rad*c
            call plot(xx+xxx,yy+yyy,3)
            call plot(xx-xxx,yy-yyy,2)
            call plot(xx-xxx,yy-yyy,3)
        return
        end
 
c-----
chttp://gams.nist.gov/serve.cgi/ModuleComponent/1070/Fullsource/ITL/RS.f
c
c       Modified RBHerrmann, Saint Louis University 07 JAN 2003 to
c       use DOUBLE PRECISION INSTEAD OF REAL. THIS ALSO
c       REQUIRED AMIN1 -> DMIN1 AMAX1 -> DMAX1
c-----
      SUBROUTINE TQL2(NM,N,D,E,Z,IERR)
C***BEGIN PROLOGUE  TQL2
C***DATE WRITTEN   760101   (YYMMDD)
C***REVISION DATE  830518   (YYMMDD)
C***CATEGORY NO.  D4A5,D4C2A
C***KEYWORDS  EIGENVALUES,EIGENVECTORS,EISPACK
C***AUTHOR  SMITH, B. T., ET AL.
C***PURPOSE  Compute eigenvalues and eigenvectors of symmetric
C            tridiagonal matrix.
C***DESCRIPTION
C
C     This subroutine is a translation of the ALGOL procedure TQL2,
C     Wilkinson.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 227-240(1971).
C
C     This subroutine finds the eigenvalues and eigenvectors
C     of a SYMMETRIC TRIDIAGONAL matrix by the QL method.
C     The eigenvectors of a FULL SYMMETRIC matrix can also
C     be found if  TRED2  has been used to reduce this
C     full matrix to tridiagonal form.
C
C     On Input
C
C        NM must be set to the row dimension of two-dimensional
C          array parameters as declared in the calling program
C          dimension statement.
C
C        N is the order of the matrix.
C
C        D contains the diagonal elements of the input matrix.
C
C        E contains the subdiagonal elements of the input matrix
C          in its last N-1 positions.  E(1) is arbitrary.
C
C        Z contains the transformation matrix produced in the
C          reduction by  TRED2, if performed.  If the eigenvectors
C          of the tridiagonal matrix are desired, Z must contain
C          the identity matrix.
C
C      On Output
C
C        D contains the eigenvalues in ascending order.  If an
C          error exit is made, the eigenvalues are correct but
C          unordered for indices 1,2,...,IERR-1.
C
C        E has been destroyed.
C
C        Z contains orthonormal eigenvectors of the symmetric
C          tridiagonal (or full) matrix.  If an error exit is made,
C          Z contains the eigenvectors associated with the stored
C          eigenvalues.
C
C        IERR is set to
C          Zero       for normal return,
C          J          if the J-th eigenvalue has not been
C                     determined after 30 iterations.
C
C     Calls PYTHAG(A,B) for sqrt(A**2 + B**2).
C
C     Questions and comments should be directed to B. S. Garbow,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C     ------------------------------------------------------------------
C***REFERENCES  B. T. SMITH, J. M. BOYLE, J. J. DONGARRA, B. S. GARBOW,
C                 Y. IKEBE, V. C. KLEMA, C. B. MOLER, *MATRIX EIGEN-
C                 SYSTEM ROUTINES - EISPACK GUIDE*, SPRINGER-VERLAG,
C                 1976.
C***ROUTINES CALLED  PYTHAG
C***END PROLOGUE  TQL2
C
        implicit none
      INTEGER I,J,K,L,M,N,II,L1,L2,NM,MML,IERR
      DOUBLE PRECISION D(N),E(N),Z(NM,N)
      DOUBLE PRECISION B,C,C2,C3,DL1,EL1,F,G,H,P,R,S,S2
      DOUBLE PRECISION PYTHAG
C
C***FIRST EXECUTABLE STATEMENT  TQL2
      IERR = 0
      IF (N .EQ. 1) GO TO 1001
C
      DO 100 I = 2, N
  100 E(I-1) = E(I)
C
      F = 0.0E0
      B = 0.0E0
      E(N) = 0.0E0
C
      DO 240 L = 1, N
         J = 0
         H = ABS(D(L)) + ABS(E(L))
         IF (B .LT. H) B = H
C     .......... LOOK FOR SMALL SUB-DIAGONAL ELEMENT ..........
         DO 110 M = L, N
            IF (B + ABS(E(M)) .EQ. B) GO TO 120
C     .......... E(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
C                THROUGH THE BOTTOM OF THE LOOP ..........
  110    CONTINUE
C
  120    IF (M .EQ. L) GO TO 220
  130    IF (J .EQ. 30) GO TO 1000
         J = J + 1
C     .......... FORM SHIFT ..........
         L1 = L + 1
         L2 = L1 + 1
         G = D(L)
         P = (D(L1) - G) / (2.0E0 * E(L))
         R = PYTHAG(P,1.0D+00)
         D(L) = E(L) / (P + SIGN(R,P))
         D(L1) = E(L) * (P + SIGN(R,P))
         DL1 = D(L1)
         H = G - D(L)
         IF (L2 .GT. N) GO TO 145
C
         DO 140 I = L2, N
  140    D(I) = D(I) - H
C
  145    F = F + H
C     .......... QL TRANSFORMATION ..........
         P = D(M)
         C = 1.0E0
         C2 = C
         EL1 = E(L1)
         S = 0.0E0
         MML = M - L
C     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
         DO 200 II = 1, MML
            C3 = C2
            C2 = C
            S2 = S
            I = M - II
            G = C * E(I)
            H = C * P
            IF (ABS(P) .LT. ABS(E(I))) GO TO 150
            C = E(I) / P
            R = SQRT(C*C+1.0E0)
            E(I+1) = S * P * R
            S = C / R
            C = 1.0E0 / R
            GO TO 160
  150       C = P / E(I)
            R = SQRT(C*C+1.0E0)
            E(I+1) = S * E(I) * R
            S = 1.0E0 / R
            C = C * S
  160       P = C * D(I) - S * G
            D(I+1) = H + S * (C * G + S * D(I))
C     .......... FORM VECTOR ..........
            DO 180 K = 1, N
               H = Z(K,I+1)
               Z(K,I+1) = S * Z(K,I) + C * H
               Z(K,I) = C * Z(K,I) - S * H
  180       CONTINUE
C
  200    CONTINUE
C
         P = -S * S2 * C3 * EL1 * E(L) / DL1
         E(L) = S * P
         D(L) = C * P
         IF (B + ABS(E(L)) .GT. B) GO TO 130
  220    D(L) = D(L) + F
  240 CONTINUE
C     .......... ORDER EIGENVALUES AND EIGENVECTORS ..........
      DO 300 II = 2, N
         I = II - 1
         K = I
         P = D(I)
C
         DO 260 J = II, N
            IF (D(J) .GE. P) GO TO 260
            K = J
            P = D(J)
  260    CONTINUE
C
         IF (K .EQ. I) GO TO 300
         D(K) = D(I)
         D(I) = P
C
         DO 280 J = 1, N
            P = Z(J,I)
            Z(J,I) = Z(J,K)
            Z(J,K) = P
  280    CONTINUE
C
  300 CONTINUE
C
      GO TO 1001
C     .......... SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 30 ITERATIONS ..........
 1000 IERR = L
 1001 RETURN
      END
      SUBROUTINE TRED2(NM,N,A,D,E,Z)
C***BEGIN PROLOGUE  TRED2
C***DATE WRITTEN   760101   (YYMMDD)
C***REVISION DATE  830518   (YYMMDD)
C***CATEGORY NO.  D4C1B1
C***KEYWORDS  EIGENVALUES,EIGENVECTORS,EISPACK
C***AUTHOR  SMITH, B. T., ET AL.
C***PURPOSE  Reduce real symmetric matrix to symmetric tridiagonal
C            matrix using and accumulating orthogonal transformation
C***DESCRIPTION
C
C     This subroutine is a translation of the ALGOL procedure TRED2,
C     NUM. MATH. 11, 181-195(1968) by Martin, Reinsch, and Wilkinson.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
C
C     This subroutine reduces a REAL SYMMETRIC matrix to a
C     symmetric tridiagonal matrix using and accumulating
C     orthogonal similarity transformations.
C
C     On Input
C
C        NM must be set to the row dimension of two-dimensional
C          array parameters as declared in the calling program
C          dimension statement.
C
C        N is the order of the matrix.
C
C        A contains the real symmetric input matrix.  Only the
C          lower triangle of the matrix need be supplied.
C
C     On Output
C
C        D contains the diagonal elements of the tridiagonal matrix.
C
C        E contains the subdiagonal elements of the tridiagonal
C          matrix in its last N-1 positions.  E(1) is set to zero.
C
C        Z contains the orthogonal transformation matrix
C          produced in the reduction.
C
C        A and Z may coincide.  If distinct, A is unaltered.
C
C     Questions and comments should be directed to B. S. Garbow,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C     ------------------------------------------------------------------
C***REFERENCES  B. T. SMITH, J. M. BOYLE, J. J. DONGARRA, B. S. GARBOW,
C                 Y. IKEBE, V. C. KLEMA, C. B. MOLER, *MATRIX EIGEN-
C                 SYSTEM ROUTINES - EISPACK GUIDE*, SPRINGER-VERLAG,
C                 1976.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  TRED2
C
        implicit none
      INTEGER I,J,K,L,N,II,NM,JP1
      DOUBLE PRECISION A(NM,N),D(N),E(N),Z(NM,N)
      DOUBLE PRECISION F,G,H,HH,SCALE
C
C***FIRST EXECUTABLE STATEMENT  TRED2
      DO 100 I = 1, N
C
         DO 100 J = 1, I
            Z(I,J) = A(I,J)
  100 CONTINUE
C
      IF (N .EQ. 1) GO TO 320
C     .......... FOR I=N STEP -1 UNTIL 2 DO -- ..........
      DO 300 II = 2, N
         I = N + 2 - II
         L = I - 1
         H = 0.0E0
         SCALE = 0.0E0
         IF (L .LT. 2) GO TO 130
C     .......... SCALE ROW (ALGOL TOL THEN NOT NEEDED) ..........
         DO 120 K = 1, L
  120    SCALE = SCALE + ABS(Z(I,K))
C
         IF (SCALE .NE. 0.0E0) GO TO 140
  130    E(I) = Z(I,L)
         GO TO 290
C
  140    DO 150 K = 1, L
            Z(I,K) = Z(I,K) / SCALE
            H = H + Z(I,K) * Z(I,K)
  150    CONTINUE
C
         F = Z(I,L)
         G = -SIGN(SQRT(H),F)
         E(I) = SCALE * G
         H = H - F * G
         Z(I,L) = F - G
         F = 0.0E0
C
         DO 240 J = 1, L
            Z(J,I) = Z(I,J) / H
            G = 0.0E0
C     .......... FORM ELEMENT OF A*U ..........
            DO 180 K = 1, J
  180       G = G + Z(J,K) * Z(I,K)
C
            JP1 = J + 1
            IF (L .LT. JP1) GO TO 220
C
            DO 200 K = JP1, L
  200       G = G + Z(K,J) * Z(I,K)
C     .......... FORM ELEMENT OF P ..........
  220       E(J) = G / H
            F = F + E(J) * Z(I,J)
  240    CONTINUE
C
         HH = F / (H + H)
C     .......... FORM REDUCED A ..........
         DO 260 J = 1, L
            F = Z(I,J)
            G = E(J) - HH * F
            E(J) = G
C
            DO 260 K = 1, J
               Z(J,K) = Z(J,K) - F * E(K) - G * Z(I,K)
  260    CONTINUE
C
  290    D(I) = H
  300 CONTINUE
C
  320 D(1) = 0.0E0
      E(1) = 0.0E0
C     .......... ACCUMULATION OF TRANSFORMATION MATRICES ..........
      DO 500 I = 1, N
         L = I - 1
         IF (D(I) .EQ. 0.0E0) GO TO 380
C
         DO 360 J = 1, L
            G = 0.0E0
C
            DO 340 K = 1, L
  340       G = G + Z(I,K) * Z(K,J)
C
            DO 360 K = 1, L
               Z(K,J) = Z(K,J) - G * Z(K,I)
  360    CONTINUE
C
  380    D(I) = Z(I,I)
         Z(I,I) = 1.0E0
         IF (L .LT. 1) GO TO 500
C
         DO 400 J = 1, L
            Z(I,J) = 0.0E0
            Z(J,I) = 0.0E0
  400    CONTINUE
C
  500 CONTINUE
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION PYTHAG(A,B)
        implicit none
C***BEGIN PROLOGUE  PYTHAG
C***REFER TO  EISDOC
C
C     Finds sqrt(A**2+B**2) without overflow or destructive underflow
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  PYTHAG
      DOUBLE PRECISION A,B
C
      DOUBLE PRECISION P,Q,R,S,T
C***FIRST EXECUTABLE STATEMENT  PYTHAG
      P = DMAX1(ABS(A),ABS(B))
      Q = DMIN1(ABS(A),ABS(B))
      IF (Q .EQ. 0.0E0) GO TO 20
   10 CONTINUE
         R = (Q/P)**2
         T = 4.0E0 + R
         IF (T .EQ. 4.0E0) GO TO 20
         S = R/T
         P = P + 2.0E0*P*S
         Q = Q*S
      GO TO 10
   20 PYTHAG = P
      RETURN
      END
        subroutine  boxrot(xc,yc,xwidt,xleng,trn,kolor)
c-----
c       put in a BOX
c-----
        sl = xleng*sin(trn*3.1415927/180.0)
        cl = xleng*cos(trn*3.1415927/180.0)
        sw = xwidt*sin(trn*3.1415927/180.0)
        cw = xwidt*cos(trn*3.1415927/180.0)
        x1 = xc + sl + cw
        x2 = xc + sl - cw
        x3 = xc - sl + cw
        x4 = xc - sl - cw
        y1 = yc + cl - sw
        y2 = yc + cl + sw
        y3 = yc - cl - sw
        y4 = yc - cl + sw
        call newpen(kolor)
        call shadet(x1,y1,x2,y2,x3,y3,0,0,0.01,0.01)
        call shadet(x2,y2,x3,y3,x4,y4,0,0,0.01,0.01)
        call newpen(1)
        call plot(x1,y1,3)
        call plot(x2,y2,2)
        call plot(x4,y4,2)
        call plot(x3,y3,2)
        call plot(x1,y1,2)
        call plot(x1,y1,3)
        return
        end

        subroutine boxzob(xc,yc,radc,trnp,plnp,trnt,plnt,
     1      trnb,plnb,kolor)
        character regime*2
        integer red, green, blue
c-----
c       first determine the type of mechanism according to Zoback 1992)
c       dir = direction of maximum compressive stress
c-----
        red = 1000
        green = 1050
        blue = 1100
        regime = 'XX'
        if(plnp.lt.40. .and. plnb.ge.45 .and. plnt.le.20)then
            regime = 'SS'
            dir = trnt + 90.
            kolor1 = green
            kolor2 = green
        else if(plnp.le.20. .and. plnb.ge.45 .and. plnt.le.40)then
            regime = 'SS'
            dir = trnp 
            kolor1 = green
            kolor2 = green
        else if(plnp .ge. 52 .and. plnb.le.45 .and. plnt .le. 35.0)then
            regime = 'NF'
            dir = trnb
            kolor1 = red
            kolor2 = red
        else if(plnp.ge.40. .and. plnp.lt.52 .and. plnb.le.45 .and.
     1      plnt.le.20.)then
            regime = 'NS'
            dir = trnp + 90.0
            kolor1 = red
            kolor2 = green
        else if(plnp.le.20. .and. plnb.le.45 .and. plnt.ge.40. .and.
     1      plnt.le.52.)then
            regime = 'TS'
            dir = trnp
            kolor1 = blue
            kolor2 = green
        else if(plnp.le.35. .and. plnb.le.45 .and. plnt.ge.52.)then
            regime = 'TF'
            dir = trnp
            kolor1 = blue
            kolor2 = blue
        endif
        write(6,'(a,1x,a2,2x,6f6.0)')
     1      'Regime:',regime,trnp,plnp,trnt,plnt,trnb,plnb
        
c-----
c       put in a BOX
c-----
        if(regime.ne.'XX')then
            sl = 1.2*radc*sin(dir*3.1415927/180.0)
            cl = 1.2*radc*cos(dir*3.1415927/180.0)
            sw = 0.15*radc*sin(dir*3.1415927/180.0)
            cw = 0.15*radc*cos(dir*3.1415927/180.0)
            x1 = xc + sl + cw
            x2 = xc + sl - cw
            x3 = xc - sl + cw
            x4 = xc - sl - cw
            y1 = yc + cl - sw
            y2 = yc + cl + sw
            y3 = yc - cl - sw
            y4 = yc - cl + sw
            call newpen(kolor1)
            call shadet(x1,y1,x2,y2,x3,y3,0,0,0.01,0.01)
            call newpen(kolor2)
            call shadet(x2,y2,x3,y3,x4,y4,0,0,0.01,0.01)
            call newpen(1)
            call plot(x1,y1,3)
            call plot(x2,y2,2)
            call plot(x4,y4,2)
            call plot(x3,y3,2)
            call plot(x1,y1,2)
            call plot(x1,y1,3)
        endif
        return
        end


        subroutine circle(rad,x0,y0,nullit)
        real*4 x0, y0, rad
c-----
c       draw a circle of radius rad centered at (x0,y0)
c-----
        logical nullit
        degrad = 3.1415927/180.0
c-----
c       if necessary, reset the region behind the focal circle to the
c       background color
c-----
        if(nullit)then
            call newpen(0)
            ipatx=0
            ipaty=0
            xlen=0.01
            ylen=0.01
            x1 = x0
            y1 = y0
            x2 = x0 + rad
            y2 = y0
            do 100 i=0,360,5
                arg = i*degrad
                x3 = x0 + 1.15*rad*cos(arg)
                y3 = y0 + 1.15*rad*sin(arg)
                call shadet(x1,y1,x2,y2,x3,y3,
     1              ipatx,ipaty,xlen,ylen)
                x2 = x3
                y2 = y3
  100       continue
            call newpen(1)
        endif
        ipen = 3
        do 200 i=0,360,5
            arg = i*degrad
            x = x0 + rad*cos(arg)
            y = y0 + rad*sin(arg)
            call plot(x,y,ipen)
            ipen = 2
  200   continue
        call plot(x,y,3)
        return
        end


        subroutine symfil(xx,yy,ht)
        real*4 x(9),y(9)
        data x/0.9238,0.3827,-0.3827,-0.9238,-0.9238,-0.3827,
     1      0.3827,0.9238,0.9238/
        data y/0.3827,0.9238,0.9238,0.3827,-0.3827,-0.9238,
     1      -0.9238,-0.3827,0.3827/
c-----
c       put in a filled octagon
c-----
        do 100 i=1,8
            xxx = xx + ht*x(i)
            yyy = yy + ht*y(i)
            xxxx= xx + ht*x(i+1)
            yyyy= yy + ht*y(i+1)
            call newpen(1)
            call shadet(xx,yy,xxx,yyy,xxxx,yyyy,0,0,0.01,0.01)
            call newpen(1000)
            call plot(xxx,yyy,3)
            call plot(xxxx,yyyy,2)
  100   continue
        call newpen(1)
        return
        end

