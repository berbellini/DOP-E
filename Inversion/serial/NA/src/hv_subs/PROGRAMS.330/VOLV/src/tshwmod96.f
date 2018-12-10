        program shwmod96
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME V                                                       c
c                                                                     c
c      PROGRAM: SHWMOD96                                              c
c                                                                     c
c      COPYRIGHT 2000 R. B. Herrmann                                  c
c                                                                     c
c      Department of Earth and Atmospheric Sciences                   c
c      Saint Louis University                                         c
c      221 North Grand Boulevard                                      c
c      St. Louis, Missouri 63103                                      c
c      U. S. A.                                                       c
c                                                                     c
c---------------------------------------------------------------------c
c       This program plot model96 models files
c
c       shwmod96 -XLEN xlen -YLEN ylen -X0 x0 -Y0 y0 
c           -VMIN vmin -VMAX vmax -ZMIN zmin -ZMAX zmax
c           -K kolor -W width
c           model96_file(s)
c-----
c       Changes
c
c       07 SEP 2001 permit choice of P-velocity, 
c           S-velocity or density for plot
c       NOTE TO add Q info, must plot 1/Q or Q consistently
c       27 JAN 2002 added line width option
c-----
c-----
c-----
c       LIN I*4 - logical unit for standard input
c       LOT I*4 - logical unit for standard output
c-----
        integer LIN, LOT, LER
        parameter (LIN=5,LOT=6,LER=0)
c-----
c       command line arguments
c-----
        integer NMODMX
        parameter (NMODMX=500)
        real x0, y0, xlen, ylen, vmin, vmax, zmin, zmax
        character*80  mname(NMODMX)
        integer kolor
        integer kpar
        real wid
        logical dolabx, dolaby, doleg
c-----
c       parse the command line
c-----
        call gcmdln(x0,y0,xlen,ylen,vmin,vmax,zmin,zmax,
     1      nmod,mname,NMODMX,kolor,kpar,wid,dolabx,dolaby,doleg)
        if(nmod.gt.0)then
c-----
c       open plot file
c-----
            call pinitf('SHWMOD96.PLT')
            call newpen(1)
            call doframe(x0,y0,xlen,ylen,zmin,zmax,vmin,vmax,kpar,
     1          dolabx,dolaby)
c-----
c           plot all the models
c-----
            call pltmdl(x0,y0,xlen,ylen,zmin,zmax,vmin,vmax,
     1          mname,nmod,NMODMX,kolor,kpar,wid,doleg)
            call pend()
        endif
        end

        subroutine doframe(x0,y0,xlen,ylen,zmin,zmax,vmin,vmax,kpar,
     1      dolabx,dolaby)
        real*4 x0,y0,xlen,ylen,zmin,zmax,vmin,vmax
        integer kpar
        logical dolabx, dolaby
        call gbox(x0,y0,x0+xlen,y0+ylen)
        if(dolabx)then
            call dolinx(x0,y0     ,xlen,vmin,vmax,0.10,
     1          .true.,.true.,.false.,
     2          1,' ')
            if(kpar.eq.1)then
                call dolinx(x0,y0+ylen,xlen,vmin,vmax,0.10,
     1          .false.,.true.,.true. ,
     2          17,'P-Velocity (km/s)')
            else if(kpar.eq.2)then
                call dolinx(x0,y0+ylen,xlen,vmin,vmax,0.10,
     1          .false.,.true.,.true. ,
     2          17,'S-Velocity (km/s)')
            else if(kpar.eq.3)then
                call dolinx(x0,y0+ylen,xlen,vmin,vmax,0.10,
     1          .false.,.true.,.true. ,
     2          17,'Density (gm/cm^3)')
            else if(kpar.eq.6)then
                call dolinx(x0,y0+ylen,xlen,vmin,vmax,0.10,
     1          .false.,.true.,.true. ,
     2           3,'Eta')
            endif
        else
            call dolinx(x0,y0     ,xlen,vmin,vmax,0.10,
     1          .true.,.false.,.false.,
     2          1,' ')
            call dolinx(x0,y0+ylen,xlen,vmin,vmax,0.10,
     1          .false.,.false.,.false. ,
     2          1 ,' ')
        endif
        if(dolaby)then
            call doylin(x0     ,y0,ylen,-zmax,-zmin,0.10,
     1          .false.,.true. ,.true. ,10,'Depth (km)',.true.)
            call doylin(x0+xlen,y0,ylen,-zmax,-zmin,0.10,
     1          .true. ,.true. ,.false.,1,' ',.true.)
        else
            call doliny(x0     ,y0,ylen,-zmax,-zmin,0.10,
     1          .false.,.true. ,.false.,1 ,' ')
            call doliny(x0+xlen,y0,ylen,-zmax,-zmin,0.10,
     1          .true. ,.true. ,.false.,1,' ')
        endif
        return
        end

        subroutine pltmdl(x0,y0,xlen,ylen,zmin,zmax,vmin,vmax,
     1          mname,nmod,NMODMX,kolor,kpar,wid,doleg)
c-----
c-----
        real x0, y0, xlen, ylen, zmin, zmax, vmin, vmax
        integer nmod, NMODMX
        character*80 mname(NMODMX)
        integer kolor
        integer kpar
        real wid
        logical doleg
c-----
c       kpar    I*4 1 plot P-velocity (km/sec)
c       kpar    I*4 2 plot S-velocity (km/sec)
c       kpar    I*4 3 plot density    (gm/cc)
c       kpar    I*4 6 plot Eta = F/(A-2L)
c-----

        logical ext
        integer ls, lt
c-----
c       LIN I*4 - logical unit for standard input
c       LOT I*4 - logical unit for standard output
c-----
        integer LIN, LOT, LER
        parameter (LIN=5,LOT=6,LER=0)
c-----
c       NLAY    I*4 - maximum number of layers in the model
c       NBDY    I*4 - maximum number of boundaries in the 
c                   ray description
c                   potentially this in NLAY+2
c-----
        integer NLAY, NBDY
        parameter (NLAY=100, NBDY=NLAY+2)
c-----
c       mname   C*80    - name of the model file
c       title   C*80    - title of the model file
c       ext L   - logical variable to see if the model file 
c                   exists
c       ierr    I*4 - 0 model file read in correctly
c                 -1 model file not read in
c-----
c       mmax    I*4 - number of layers in the model, last layer is
c                    halfspace
c       iunit   I*4 - 0 Kilometer, Gram, Sec
c       iiso    I*4 - 0 isotropic 
c                 1 transversely anisotropic 
c                 2 general anisotropic 
c       iflsph  I*4 - 0 flat earth model
c                 1 spherical earth model
c-----
        integer*4 mmax, iunit, iiso, iflsph
        character*80 title
        integer ierr
c-----
c       get the earth model(s)
c-----
c  note if the depth limits are not given search 
c           the model for everything
        yleg = y0 + ylen
        xleg = x0 + xlen
        do 1000 i=1,nmod
            inquire(file=mname(i),exist=ext)
            if(.not. ext)go to 9000
            ls = lgstr(mname(i))

            write(LOT,*)'Model name: ',mname(i)(1:ls)
            if(kolor.ge.0)then
                ipen = kolor
            else
                if(nmod.eq.1)then
                    ipen = 2
                else
                    ipen = 1000 + real(100 * 
     1                  real(i-1)/real(nmod-1)) 
                    if(ipen.gt.1100)ipen = 1100
                    if(ipen.lt.1000)ipen = 1000
                endif
            endif
            if(doleg)then
                call newpen(ipen)
                call plot(xleg+0.1,yleg,3)
                call gwidth(wid)
                call plot(xleg+0.3,yleg,2)
                call newpen(1)
                dy = 0.14
                ht = 0.71 * dy
                call symbol(xleg+0.4,yleg-0.5*ht,ht,
     1              mname(i),0.0,ls)
                yleg = yleg - dy
                call gwidth(0.001)
            endif
            call newpen(ipen)
            call getmod(1,mname(i),mmax,title,iunit,iiso,iflsph,
     1          idimen,icnvel,ierr,.true.)
            lt = lgstr(title)
            write(LOT,*)'title:',title(1:lt)
            write(LOT,*)'mmax:',mmax
            write(LOT,*)'iunit:',iunit
            write(LOT,*)'iiso:',iiso
            write(LOT,*)'iflsph:',iflsph
            write(LOT,*)'idimen:',idimen
            write(LOT,*)'icnvel:',icnvel
c-----
c       error checking
c-----
            if(idimen.ne.1)then
                write(LER,*)'1-D velocity model required'
                go to 9000
            endif
            if(icnvel.ne.0)then
                write(LER,*)'constant velocity model required'
                go to 9000
            endif
c-----
c           model is OK to plot, and everything necessary is in common
c-----
            if(kpar.eq.1)then
                if(iiso.eq.1)then
            call dopltmd(iiso,x0,y0,xlen,ylen,zmin,zmax,
     1          vmin,vmax,mmax,mname(i),kpar,wid,'TA')
                endif
            call dopltmd(iiso,x0,y0,xlen,ylen,zmin,zmax,
     1          vmin,vmax,mmax,mname(i),kpar,wid,'TC')
            else if(kpar.eq.2)then
                if(iiso.eq.1)then
            call dopltmd(iiso,x0,y0,xlen,ylen,zmin,zmax,
     1          vmin,vmax,mmax,mname(i),kpar,wid,'TN')
                endif
            call dopltmd(iiso,x0,y0,xlen,ylen,zmin,zmax,
     1          vmin,vmax,mmax,mname(i),kpar,wid,'TL')
            else if(kpar.eq.3)then
            call dopltmd(iiso,x0,y0,xlen,ylen,zmin,zmax,
     1          vmin,vmax,mmax,mname(i),kpar,wid,'TR')
            else if(kpar.eq.6)then
            call dopltmd(iiso,x0,y0,xlen,ylen,zmin,zmax,
     1          vmin,vmax,mmax,mname(i),kpar,wid,'ET')
            endif
            call newpen(1)
 9000       continue
 1000   continue
        return
        end

        subroutine dopltmd(iiso,x0,y0,xlen,ylen,zmin,zmax,
     1          vmin,vmax,mmax,mname,kpar,wid,cmd)
        real*4 x0, y0, xlen, ylen, vmin, vmax, zmin, zmax
        integer*4 mmax
        integer kpar
        real wid
        character cmd*(*)
c-----
c       iiso    I*4 - 0 isotropic
c                 1 transversely anisotropic
c                 2 general anisotropic
c       kpar    I*4 1 plot P-velocity (km/sec)
c       kpar    I*4 2 plot S-velocity (km/sec)
c       kpar    I*4 3 plot density    (gm/cc)
c       kpar    I*4 6 plot Eta = F/(A-2L)
c       wid R   line width
c       cmd C*  Indicates what to plot TA TC TN TL Density Eta
c-----
        character mname*(*)
c-----
c       model parameters from common block
c-----
        integer NL
        parameter(NL=200)
        common/timodel/D(NL),TA(NL),TC(NL),TL(NL),TN(NL),TF(NL),
     1      TRho(NL),
     2      qa(NL),qb(NL),etap(NL),etas(NL), 
     3      frefp(NL), frefs(NL)
        real D, TA, TC, TF, TL, TN, TRho
        real qa, qb, etap, etas, frefp, frefs
        common/depref/refdep
        real refdep

        real pval
        logical dodash
c-----
c       At this point the program igetmod has made all
c       layer thicknesses positive and indicated where the
c       0 reference depth (- = above) is located
c-----
c-----
c       lets plot the parameter
c       For isotropic media the A and N values are plotted with dashes
c       These corresponding to the VPH and VSH (quasi-SH propagating
c       horizontally)
c-----
        if(cmd.eq.'TA' .or. cmd.eq.'TN')then
            dodash =.true.
        else
            dodash = .false.
        endif
        dvdx =  xlen/(vmax-vmin)
        dzdz = -ylen/(zmax-zmin)
        xx0 = x0
        yy0 = y0 + ylen

c-----
c       to do a discontinuous model, we need
c-----
        z0 = 0.0
        zz = zmin
        if(cmd.eq.'TA')then
            pval = sqrt(TA(1)/TRho(1))
        else if(cmd.eq.'TC')then
            pval = sqrt(TC(1)/TRho(1))
        else if(cmd.eq.'TN')then
            pval = sqrt(TN(1)/TRho(1))
        else if(cmd.eq.'TL')then
            pval = sqrt(TL(1)/TRho(1))
        else if(cmd.eq.'TR')then
            pval = TRho(1)
        else if(cmd.eq.'ET')then
            pval = TF(1)/(TA(1)-2.0*TL(1))
        endif
        xx = xx0 + (pval - vmin)*dvdx
c-----
c       spatially clip
c-----
        call gclip('on',x0,y0,x0+xlen,y0+ylen)

        zz = -refdep
        if(zz .lt. zmin)then
            yy = yy0
        else
            yy = yy0 + (zz   - zmin)*dzdz
        endif
        call gwidth(wid)
        call plot(xx,yy,3)
        do 1000 i=2,mmax
            zz = zz + d(i-1) 

            if(zz.gt.zmax)go to 1000
            if(zz .lt. zmin)then
                yy = yy0
            else
                yy = yy0 + (zz   - zmin)*dzdz
            endif
            if(dodash)then
                call plotd(xx,yy,21,0.1)
            else
                call plot(xx,yy,2)
            endif
        if(cmd.eq.'TA')then
            pval = sqrt(TA(i)/TRho(i))
        else if(cmd.eq.'TC')then
            pval = sqrt(TC(i)/TRho(i))
        else if(cmd.eq.'TN')then
            pval = sqrt(TN(i)/TRho(i))
        else if(cmd.eq.'TL')then
            pval = sqrt(TL(i)/TRho(i))
        else if(cmd.eq.'TR')then
            pval = TRho(i)
        else if(cmd.eq.'ET')then
            pval = TF(i)/(TA(i)-2.0*TL(i))
        endif
            xx = xx0 + (pval - vmin)*dvdx
            if(dodash)then
                call plotd(xx,yy,21,0.1)
            else
                call plot(xx,yy,2)
            endif
 1000   continue
        if(dodash)then
            call plotd(xx,y0,21,0.1)
        else
            call plot(xx,y0,2)
        endif
        call plot(xx,y0,3)
        call gwidth(0.001)
        call gclip('off',x0,y0,x0+xlen,y0+ylen)
            
        return
        end
        

        subroutine gcmdln(x0,y0,xlen,ylen,vmin,vmax,zmin,zmax,
     1      nmod,mname,NMODMX,kolor,kpar,wid,dolabx,dolaby,doleg)
c-----
c       parse the command line arguments
c-----
c       x0  R*4 - lower left corner of plot frame
c       y0  R*4 - lower left corner of plot frame
c       xlen    R*4 - width  of plot frame on page
c       ylen    R*4 - height of plot frame on page
c       kolor   I   - pen color for plto, if < 0 use rainbow
c                   and annotate the plot
c       vmin    R*4 - minimum velocity for plot
c       vmax    R*4 - maximum velocity for plot
c       zmin    R*4 - minimum depth for plot
c       zmax    R*4 - mazimum depth for plot
c       nmod    I*4 - maximum number of model files to overlay
c       mname   Ch*80   - array of model files
c       NMODMX  I*4 - array dimension of mname
c       kpar    I*4 - Plot parameter control
c                   1 P-velocity
c                   2 S-velocity (default)
c                   3 density
c                   4 QP
c                   5 QS (note fix this for 1/Q Q stuff!!!)
c                   6 ETA = F/(A-2L)
c       wid R   - line width in inches
c       dolabx  L   - .true. label X-axis with velocity/density
c       dolaby  L   - .true. label Y-axis with depth
c       doleg   L   - .true. put in file legend to right of plot
c-----
        real*4 x0, y0, xlen, ylen, zmin, zmax, vmin, vmax
        real wid
        logical dolabx, dolaby, doleg
        integer nmod, NMODMX
        integer kolor
        character*80 mname(NMODMX)

        integer mnmarg
        character name*80

        x0 = 2.0
        y0 = 1.0
        xlen = 6.0
        ylen = 6.0
        kolor = 1
        vmin = 2.0
        vmax = 5.0
        zmin = 0.0
        zmax = 60.0
        nmod = 0
        kpar = 2
        wid = 0.001
        dolabx = .true.
        dolaby = .true.
        doleg = .false.
        
        nmarg = mnmarg()
        i = 0
 1000   continue
            i = i + 1
            if(i.gt.nmarg)go to 2000
            call mgtarg(i,name)
            if(name(1:2).eq.'-K')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,i10)')kolor
            else if(name(1:5).eq.'-VMIN')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,f10.0)')vmin
            else if(name(1:5).eq.'-VMAX')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,f10.0)')vmax
            else if(name(1:5).eq.'-ZMIN')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,f10.0)')zmin
            else if(name(1:5).eq.'-ZMAX')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,f10.0)')zmax
            else if(name(1:3).eq.'-X0')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,f10.0)')x0
            else if(name(1:3).eq.'-Y0')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,f10.0)')y0
            else if(name(1:5).eq.'-XLEN')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,f10.0)')xlen
            else if(name(1:5).eq.'-YLEN')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,f10.0)')ylen
            else if(name(1:2) .eq. '-P')then
                kpar = 1
            else if(name(1:2) .eq. '-S')then
                kpar = 2
            else if(name(1:2) .eq. '-D')then
                kpar = 3
            else if(name(1:2) .eq. '-E')then
                kpar = 6
            else if(name(1:2).eq.'-W')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,f10.0)')wid
                if(wid .le. 0.0)wid = 0.01
            else if(name(1:7).eq.'-NOLABX')then
                dolabx = .false.
            else if(name(1:7).eq.'-NOLABY')then
                dolaby = .false.
            else if(name(1:4).eq.'-LEG')then
                doleg = .true.
            else if(name(1:2) .eq. '-?')then
                call usage()
            else if(name(1:2) .eq. '-h')then
                call usage()
            else
                if(nmod.lt.NMODMX)then
                    nmod = nmod +1
                    mname(nmod) = name
                endif
            endif
        go to 1000
 2000   continue
c-----
c       safety
c-----
        if(kpar.lt.1.or.kpar.gt.7)then
            kpar = 2
        endif
        return
        end

        subroutine usage()
        integer LER
        parameter (LER=0)
        write(LER,*)'Usage: shwmod96 -XLEN xlen -YLEN ylen',
     1      ' -X0 x0 -Y0 y0 ',
     2      '-VMIN vmin -VMAX vmax -ZMIN zmin -ZMAX zmax',
     3      '-K kolor [-P -S -D ] [ -W width ] ',
     4      ' [-NOLABX -NOLABY] [-LEG] model96_file[s]'
        write(LER,*)
     1  '-XLEN xlen (default 6.0 )  Length of horizontal axis'
        write(LER,*)
     1  '-YLEN ylen (default 6.0 )  Length of depth axis'
        write(LER,*)
     1  '-X0 x0     (default 2.0 )  (x0,y0) are lower left corner'
        write(LER,*)
     1  '-Y0 y0     (default 1.0 )  '
        write(LER,*)
     1  '-K kolor   (default  1  )  Profile in color',
     2  ' if kolor < 0 use red->blue progression'
        write(LER,*)
     1  '-VMIN vmin (default 2.0 )  Minimum value of horizontal'
        write(LER,*)
     1  '-VMAX vmax (default 5.0 )  Minimum value of horizontal'
        write(LER,*)
     1  '-ZMIN zmin (default 0.0 )  Minimum value of depth axis'
        write(LER,*)
     1  '-ZMAX zmax (default 60. )  Minimum value of h axisorizontal'
        write(LER,*)
     1  '-W   width (default 0.001) Width of line (inch) for model plot'
        write(LER,*)
     1  '-NOLABX    (default label X) Do not label X-axis'
        write(LER,*)
     1  '-NOLABY    (default label Y) Do not label Y-axis'
        write(LER,*)
     1  '-LEG       (default none) Put in file legend'
        write(LER,*)
     1  '-P         (default S )  plot P-velocity '
        write(LER,*)
     1  '-S         (default S )  plot S-velocity '
        write(LER,*)
     1  '-D         (default S )  plot density '
        write(LER,*)
     1  '-E         (default S )  plot ETA=F/(A-2L) '
        write(LER,*)
     1  '-?         (default none )  this help message '
        write(LER,*)
     1  '-h         (default none )  this help message '
        stop
        end
        
        subroutine doylin(x0,y0,yleng,ymax,ymin,
     1      sizey,ticlft,lablft,dopow,lly,titley,dorev)
c-----
c       put in tic marks along the y-axis 
c
c       x0  R*4 - position of bottom side of axis
c       y0  R*4 - position of bottom side of axis
c       yleng   R*4 - length of Y axis
c       ymax    R*4 - maximum value of number corresponding to
c                   top
c       ymax    R*4 - maximum value of number corresponding to
c                   bottom
c       sizey   R*4 - size of numbers, e.g., 10
c       ticlft  L   - .true. tics go left
c                 .false. tics go right
c       lablft  L   - .true. numbers goes to left of axis
c                 .false. goes to right
c       dopow   L   - .true. put in numbers labels
c       ly  I*4 length of y-axis title string
c       titley  Ch  y-axis title string
c       dorev   L   .true. numbers increase downward
c               just change the sign of the plotted number
c-----
        real*4 x0, y0, yleng, ymax, ymin
        logical ticlft, lablft, dopow
        integer ly
        character titley*(*)
        character title*80
        logical dorev

        logical dosci

        character cmd*10

        scaly=1.0
        maxdig = 2
        ly = lly
c-----
c       put in the axis
c-----
        yl = y0 + yleng
        ycen = y0 + 0.5*yleng
        call plot(x0,y0,3) 
        call plot(x0,yl,2) 
c-----
c       set up positioning for tics
c-----
        tlen1 = 0.6*sizey
        tlen2 = 1.2*sizey
        if(lablft)then
            if(ticlft)then
                xpos   = x0 - 2*sizey
                naxdig = -(maxdig+3+2)
                xpos2  = x0 + naxdig* sizey -1.2*sizey
            else
                xpos   = x0 - sizey
                naxdig = -(maxdig+3+1)
                xpos2  = x0 + naxdig* sizey -1.2*sizey
            endif
            cmd = 'RIGHT'
        else
            if(ticlft)then
                xpos   = x0 + 1.00*sizey
                naxdig =  (maxdig+3+1)
                xpos2  = x0 + naxdig* sizey +1.0*sizey
            else
                xpos   = x0 + 2.0*sizey
                naxdig =  (maxdig+3+2)
                xpos2  = x0 + naxdig* sizey +1.0*sizey
            endif
            cmd = 'LEFT'
        endif
c-----
c       compute limits for scale
c-----
        ymn = min(ymin,ymax)
        ymx = max(ymin,ymax)
        yyv = ymn
c-----
c       safety
c-----
        if(ymn.eq.ymx)ymx = ymn + 1.0
c-----
c       get some idea of the size of the number to be plotted
c-----
        yval = max(abs(ymn),abs(ymx))
c-----
c       we set the maximum number of decimal digits as 3, 
c       if the number is greater than 100 or less than 0.01, use
c       scientific notation
c       get the maximum number of digits for the number
c-----
        if(yval.le.0.01 .and.yval.gt.0.0 )then
            dosci = .true.
            ynorm = alog10(1.0/yval)
            nynorm = ynorm+1
            ynorm = 10.0**( nynorm)
            nscaly = -nynorm
        else if( yval.ge.100.0)then
            dosci = .true.
            ynorm = alog10(yval)
            nynorm = ynorm
            ynorm = 10.0**(-nynorm)
            nscaly =  nynorm
        else
            dosci = .false.
            ynorm = 1.0
            nynorm = 0.0
            nscaly = 0
        endif
c-----
c       choose the maximum number of tics somehow on
c       yleng, and size of numbers and necessary decimals
c-----
        dyy = ynorm*(ymx - ymn)/5.0
        if(dyy.eq.0.0)dyy = 1.0
c-----
c       get start limits for search - we will place at most 10 tics
c       here
c-----
        n = alog10(dyy)
        if(dyy.lt.1.0)n = n -1
        ifac = (10.0**(-n)*dyy)
        
c-----
c       ifac should be a number between 1 and 9
c-----
        if(ifac.eq.3)then 
            ifac = 4
        else if(ifac.ge.6)then 
            ifac = 10
        endif
c----- 
c       now get a new minimum value that is some multiple 
c       of ifac 10**n 
c----- 
        dyy = ifac * 10.0**(n)
        ilow = ynorm*ymn/dyy -1
        iup  = ynorm*ymx/dyy +1
        if(ifac.eq.1)then
            ndyyy = 5
        else if(ifac.eq.2)then
            ndyyy = 4
        else
            ndyyy = ifac
        endif
        dyyy = dyy/ndyyy

c-----
c       loop for labels
c-----
c       ye is the end of the previous label
c-----
        yc = ymn
        call plot(x0,y0,3)
c-----
c       here
c       yy is the position of the plot corresponding to a value of yyv
c       yend is the position of the end of the last number plotted
c-----
        ye = y0 
        do 1000 i=ilow,iup
            do 1001 j=0,ndyyy-1
            
            yyv = i*dyy + j*dyyy
            if(i.eq.ilow .and. j.eq.0)then
                if(yyv.lt.0.0)then
                    yyy = 0.5*(maxdig+3)*sizey
                else
                    yyy = 0.5*(maxdig+2)*sizey
                endif
            endif
            yyv = yyv/ynorm
            if(j.eq.0)then
                tlen = tlen2
            else
                tlen = tlen1
            endif
            if(yyv.ge.ymn .and. yyv .le.ymx)then
                yy = y0 + (yyv - ymn)*yleng/(ymx - ymn)
                call plot(x0,yy,3)
                if(ticlft)then  
                    call plot(x0-tlen,yy,2)
                else
                    call plot(x0+tlen,yy,2)
                endif
                call plot(x0,yy,3)
                if(j.eq.0 .and. dopow)then
                    yp = yy - 0.5*sizey
                    
                    yyy = 0.5*sizey
                    yb = yy - yyy
                    yend = yy + yyy
                    if(yb.gt.amin1(ye,yl) .and. 
     1                 yend.le.amax1(ye,yl))then
                    if(dorev)then
                    call mynum(xpos,yp,sizey,-yyv*ynorm,0.0,maxdig,cmd)
                    else
                    call mynum(xpos,yp,sizey,yyv*ynorm,0.0,maxdig,cmd)
                    endif
                    endif
                endif
            endif
 1001       continue
 1000   continue
        if(dosci)then
            scaly = 1.0/ynorm
        endif
c-----
c       put in the title if we put in the numbers
c-----
        if(dopow )then
            sizeyy = 1.2*sizey
            title = ' '
            if(ly.gt.0)then
                title = titley(1:ly)
            else
                ly = 0
            endif
            if(dosci)then
c-----
c       note need space for this
c-----
                title(ly +1:ly+5)=' *10 '
                if(lablft)then
                    call gcent(xpos2,ycen,sizeyy,title, 90.0)
                    ye = ycen + 0.5*sizeyy*(ly+4)
                    call number(xpos2-0.7*sizey,ye,0.7*sizeyy,
     1                  real(nscaly), 90.0,-1)
                else
                    call gcent(xpos2,ycen,sizeyy,title,-90.0)
                    ye = ycen - 0.5*sizeyy*(ly+4)
                    call number(xpos2+0.7*sizeyy,ye,0.7*sizeyy,
     1                  real(nscaly),-90.0,-1)
                endif
            else
                if(lablft)then
                    call gcent(xpos2,ycen,sizeyy,title, 90.0)
                else
                    call gcent(xpos2,ycen,sizeyy,title,-90.0)
                endif
            endif
        endif
        return
        end

