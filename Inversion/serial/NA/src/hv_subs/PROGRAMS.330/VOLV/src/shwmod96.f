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
c       07 SEP 2001 permit choice of P-velocity, S-velocity 
c           or density for plot
c       NOTE TO add Q info, must plot 1/Q or Q consistently
c       27 JAN 2002 added line width option
c       08 APR 2003 slightly modified width option
c       10 NOV 2009 added option to place the legend inside the plot frame
c       05 MAY 2010 added dashed line option - but be careful with automatic
c            -K -1
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
        logical dolabx, dolaby, doleg, dolegin
        integer lintyp
c-----
c       parse the command line
c-----
        call gcmdln(x0,y0,xlen,ylen,vmin,vmax,zmin,zmax,
     1      nmod,mname,NMODMX,kolor,kpar,wid,dolabx,dolaby,
     2      doleg,dolegin,lintyp)
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
     1          mname,nmod,NMODMX,kolor,kpar,wid,doleg,dolegin,lintyp)
            call pend()
        endif
        end

        subroutine doframe(x0,y0,xlen,ylen,zmin,zmax,vmin,vmax,kpar,
     1          dolabx,dolaby)
        real*4 x0,y0,xlen,ylen,zmin,zmax,vmin,vmax
        integer kpar
        logical dolabx, dolaby

        if(xlen.le.2.5)then
              ht = xlen/25.0
        else
              ht = 0.10
        endif
        call gbox(x0,y0,x0+xlen,y0+ylen)
        if(dolabx)then
            call dolinx(x0,y0     ,xlen,vmin,vmax,ht,
     1          .true.,.true.,.false.,
     2          1,' ')
            if(kpar.eq.1)then
                call dolinx(x0,y0+ylen,xlen,vmin,vmax,ht,
     1          .false.,.true.,.true. ,
     2          17,'P-Velocity (km/s)')
            else if(kpar.eq.2)then
                call dolinx(x0,y0+ylen,xlen,vmin,vmax,ht,
     1          .false.,.true.,.true. ,
     2          17,'S-Velocity (km/s)')
            else if(kpar.eq.3)then
                call dolinx(x0,y0+ylen,xlen,vmin,vmax,ht,
     1          .false.,.true.,.true. ,
     2          17,'Density (gm/cm^3)')
            else if(kpar.eq.4)then
                call dolinx(x0,y0+ylen,xlen,vmin,vmax,ht,
     1          .false.,.true.,.true. ,
     2          4,'1/Qp')
            else if(kpar.eq.5)then
                call dolinx(x0,y0+ylen,xlen,vmin,vmax,ht,
     1          .false.,.true.,.true. ,
     2          4,'1/Qs')
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
            call dnliny(x0     ,y0,ylen,-zmax,-zmin,0.10,
     1          .false.,.true. ,.true. ,10,'Depth (km)')
            call dnliny(x0+xlen,y0,ylen,-zmax,-zmin,0.10,
     1          .true. ,.true. ,.false.,1,' ')
        else
            call dnliny(x0     ,y0,ylen,-zmax,-zmin,0.10,
     1          .false.,.true. ,.false.,1 ,' ')
            call dnliny(x0+xlen,y0,ylen,-zmax,-zmin,0.10,
     1          .true. ,.true. ,.false.,1,' ')
        endif
        return
        end

        subroutine pltmdl(x0,y0,xlen,ylen,zmin,zmax,vmin,vmax,
     1          mname,nmod,NMODMX,kolor,kpar,wid,doleg,dolegin,lintyp)
c-----
c-----
        real x0, y0, xlen, ylen, zmin, zmax, vmin, vmax
        integer nmod, NMODMX
        character*80 mname(NMODMX)
        integer kolor
        integer kpar
        real wid
        logical doleg,dolegin
        integer lintyp
c-----
c       kpar    I*4 1 plot P-velocity (km/sec)
c                   2 plot S-velocity (km/sec)
c                   3 plot density    (gm/cc)
c                   4 plot 1/Qp
c                   5 plot 1/Qs
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
c       model parameters from common block
c-----
        integer NL
        parameter (NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        common/depref/refdep
c-----
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
c  note if the depth limits are not given search the 
c           model for everything
        if(dolegin)then
              yleg = y0 + 0.6*ylen
              xleg = x0 + 0.05*xlen
        else
              yleg = y0 + ylen
              xleg = x0 + xlen
        endif
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
                call gwidth(0.001)
                call symbol(xleg+0.4,yleg-0.5*ht,ht,
     1              mname(i),0.0,ls)
                yleg = yleg - dy
            endif
            call newpen(ipen)
            call getmod(1,mname(i),mmax,title,iunit,iiso,iflsph,
     1          idimen,icnvel,ierr,.false.)
            lt = lgstr(title)
            write(LOT,*)'title:',title(1:lt)
            write(LOT,*)'mmax:',mmax
            write(LOT,*)'iunit:',iunit
            write(LOT,*)'iiso:',iiso
            write(LOT,*)'iflsph:',iflsph
            write(LOT,*)'idimen:',idimen
            write(LOT,*)'icnvel:',icnvel
c-----
c           convert Q to 1/Q assuming that Q >= 1
c-----
           do j=1,mmax
                if(qa(j).gt.1.0)qa(j) = 1./qa(j)    
                if(qb(j).gt.1.0)qb(j) = 1./qb(j)    
           enddo
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
            call dopltmd(x0,y0,xlen,ylen,zmin,zmax,
     1          vmin,vmax,mmax,mname(i),kpar,wid,lintyp)
            call newpen(1)
 9000       continue
 1000   continue
        return
        end

        subroutine dopltmd(x0,y0,xlen,ylen,zmin,zmax,
     1          vmin,vmax,mmax,mname,kpar,wid,lintyp)
        real*4 x0, y0, xlen, ylen, vmin, vmax, zmin, zmax
        integer*4 mmax
        integer kpar
        real wid
        integer lintyp
c-----
c       kpar    I*4 1 plot P-velocity (km/sec)
c                   2 plot S-velocity (km/sec)
c                   3 plot density    (gm/cc)
c                   4 plot 1/Qp
c                   5 plot 1/Qs
c       wid R   line width
c-----
        character mname*(*)
c-----
c       model parameters from common block
c-----
        integer NL
        parameter (NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        common/depref/refdep
c-----
c       At this point the program igetmod has made all
c       layer thicknesses positive and indicated where the
c       0 reference depth (- = above) is located
c-----
c-----
c       lets plot the parameter
c-----
        dvdx =  xlen/(vmax-vmin)
        dzdz = -ylen/(zmax-zmin)
        xx0 = x0
        yy0 = y0 + ylen

c-----
c       to do a discontinuous model, we need
c-----
        z0 = 0.0
        zz = zmin
C       write(0,*)'vmin,vmax,dvdx,dzdz:',vmin,vmax,dvdx,dzdz
        if(kpar.eq.1)then
            xx = xx0 + (a(1) - vmin)*dvdx
        else if(kpar.eq.2)then
            xx = xx0 + (b(1) - vmin)*dvdx
        else if(kpar.eq.3)then
            xx = xx0 + (rho(1) - vmin)*dvdx
        else if(kpar.eq.4)then
            xx = xx0 + (qa(1) - vmin)*dvdx
        else if(kpar.eq.5)then
            xx = xx0 + (qb(1) - vmin)*dvdx
        endif
C       write(0,*)'xx,xx0,a(1),b(1),rho(1),vmin,dvdx:',
c           xx,xx0,a(1),b(1),rho(1),vmin,dvdx
C       write(0,*)'kpar,refdep:',kpar,refdep
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
C       write(0,*)'xx,yy:',xx,yy
        do 1000 i=2,mmax
            zz = zz + d(i-1) 

            if(zz.gt.zmax)go to 1000
            if(zz .lt. zmin)then
                yy = yy0
            else
                yy = yy0 + (zz   - zmin)*dzdz
            endif
            if(lintyp.eq.1)then
                  ipat = 3
                  patlen = 0.05
                  call plotd(xx,yy,ipat,patlen)
            else if(lintyp.eq.2)then
                  ipat = 7
                  patlen = 0.05
                  call plotd(xx,yy,ipat,patlen)
            else 
                  call plot(xx,yy,2)
            endif

            if(kpar.eq.1)then
                xx = xx0 + (a(i) - vmin)*dvdx
            else if(kpar.eq.2)then
                xx = xx0 + (b(i) - vmin)*dvdx
            else if(kpar.eq.3)then
                xx = xx0 + (rho(i) - vmin)*dvdx
            else if(kpar.eq.4)then
                xx = xx0 + (qa(i) - vmin)*dvdx
            else if(kpar.eq.5)then
                xx = xx0 + (qb(i) - vmin)*dvdx
            endif
            if(lintyp.eq.1)then
                  ipat = 3
                  patlen = 0.05
                  call plotd(xx,yy,ipat,patlen)
            else if(lintyp.eq.2)then
                  ipat = 7
                  patlen = 0.05
                  call plotd(xx,yy,ipat,patlen)
            else 
                  call plot(xx,yy,2)
            endif
 1000   continue
            if(lintyp.eq.1)then
                  ipat = 3
                  patlen = 0.05
                  call plotd(xx,y0,ipat,patlen)
            else if(lintyp.eq.2)then
                  ipat = 7
                  patlen = 0.05
                  call plotd(xx,y0,ipat,patlen)
            else 
                  call plot(xx,y0,2)
            endif
        call plot(xx,y0,3)
        call gwidth(0.001)
        call gclip('off',x0,y0,x0+xlen,y0+ylen)
            
        return
        end
        

        subroutine gcmdln(x0,y0,xlen,ylen,vmin,vmax,zmin,zmax,
     1      nmod,mname,NMODMX,kolor,kpar,wid,dolabx,dolaby,
     2      doleg,dolegin,lintyp)
c-----
c       parse the command line arguments
c-----
c       x0  R*4 - lower left corner of plot frame
c       y0  R*4 - lower left corner of plot frame
c       xlen    R*4 - width  of plot frame on page
c       ylen    R*4 - height of plot frame on page
c       kolor   I   - pen color for plot, if < 0 use rainbow
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
c       wid R   - line width in inches
c       dolabx  L   - .true. label X-axis with velocity/density
c       dolaby  L   - .true. label Y-axis with depth
c       doleg   L   - .true. put in file legend to right of plot
c       delegin L   - .true. put legend insides requires -L flag
c       lintyp  I   - 0 solid, 1 short dashes, 2 long dashes
c-----
        real*4 x0, y0, xlen, ylen, zmin, zmax, vmin, vmax
        real wid
        logical dolabx, dolaby, doleg,dolegin
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
        dolegin = .false.
        lintyp = 0
        
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
            else if(name(1:3).eq.'-DT')then
                i = i + 1
                call mgtarg(i,name)
                if(name(1:6).eq."short")then
                      lintyp = 1
                else if(name(1:5).eq."long")then
                      lintyp = 2
                else if(name(1:6).eq."solid")then
                      lintyp = 0
                endif
            else if(name(1:2) .eq. '-P')then
                kpar = 1
            else if(name(1:2) .eq. '-S')then
                kpar = 2
            else if(name(1:2) .eq. '-D')then
                kpar = 3
            else if(name(1:3) .eq. '-QP')then
                kpar = 4
            else if(name(1:3) .eq. '-QS')then
                kpar = 5
            else if(name(1:2) .eq. '-W')then
                i = i + 1
                call mgtarg(i,name)
                read(name,'(bn,f10.0)')wid
                if(wid .le. 0.0)wid = 0.01
            else if(name(1:7).eq.'-NOLABX')then
                dolabx = .false.
            else if(name(1:7).eq.'-NOLABY')then
                dolaby = .false.
            else if(name(1:4).eq.'-LEG')then
                if(name(1:6).eq.'-LEGIN')then 
                   dolegin = .true.
                else
                   dolegin = .false.
                endif
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
        if(kpar.lt.1.or.kpar.gt.5)then
            kpar = 2
        endif
        if(kpar.eq.4 .or. kpar.eq.5)then
             if(vmin .eq. 2.0)then
                  vmin = 0.0
             endif
             if(vmax .eq. 5.0)then
                  vmax = 1.0
             endif
        endif
        return
        end

        subroutine usage()
        integer LER
        parameter (LER=0)
        write(LER,*)'Usage: shwmod96 -XLEN xlen -YLEN ylen',
     1      ' -X0 x0 -Y0 y0 ',
     2      '-VMIN vmin -VMAX vmax -ZMIN zmin -ZMAX zmax ',
     3      '-K kolor [-P -S -D -QP -QS ] [ -W width ] ',
     4      ' [-NOLABX -NOLABY] [-LEG] model96_file[s]',
     5      ' [ -DT lintyp ]'
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
     1  '-ZMAX zmax (default 60. )  Minimum value of horizontal axis'
        write(LER,*)
     1  '-W   width (default 0.001) Width of line (inch) for model plot'
        write(LER,*)
     1  '-NOLABX    (default label X) Do not label X-axis'
        write(LER,*)
     1  '-NOLABY    (default label Y) Do not label Y-axis'
        write(LER,*)
     1  '-LEG       (default none) Put in file legend'
        write(LER,*)
     1  '-LEGIN     (default none) Put in file legend inside frame'
        write(LER,*)
     1  '-P         (default S )  plot P-velocity '
        write(LER,*)
     1  '-S         (default S )  plot S-velocity '
        write(LER,*)
     1  '-D         (default S )  plot density '
        write(LER,*)
     1  '-QP         (default S )  plot 1/QP '
        write(LER,*)
     1  '-QS         (default S )  plot 1/QS '
        write(LER,*)
     1  '-DT linetype (default solid) linetype= solid short long '
        write(LER,*)
     1  '-?         (default none )  this help message '
        write(LER,*)
     1  '-h         (default none )  this help message '
        stop
        end
        
