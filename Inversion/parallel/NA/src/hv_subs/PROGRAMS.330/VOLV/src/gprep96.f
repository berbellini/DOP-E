        program gprep96
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME V                                                       c
c                                                                     c
c      PROGRAM: GPREP96                                               c
c                                                                     c
c      COPYRIGHT 1996 R. B. Herrmann                                  c
c                                                                     c
c      Department of Earth and Atmospheric Sciences                   c
c      Saint Louis University                                         c
c      221 North Grand Boulevard                                      c
c      St. Louis, Missouri 63103                                      c
c      U. S. A.                                                       c
c                                                                     c
c---------------------------------------------------------------------c
c-----
c       program to prepare input for genray91(V)
c       This program permits automatic ray generation
c
c       CHANGES
c          19 JUN 2009 - changed order ofat line 92 so that the data
c              statement follows the last type declaration
c-----
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
        parameter (NLAY=200, NBDY=NLAY+2)
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
        character mname*80, title*80, deny*80, bounce*80
        logical ext
        integer*4 mmax, iunit, iiso, iflsph
c-----
c       depths  R*4 - source depth
c       depthr  R*4 - receiver depth
c-----
        real*4 depths, depthr
        logical lsrcbd, lrecbd

c-----
c       ibdmap(NBDY)    I*4 - mapping array from original
c                   layers to boundaries
c                   The boundaries have been incremented
c                       by the source and receiver
c       permit(NBDY)    L  - if .false. deny reflections and conversions
c                   This is used to simulate gradients
c                   in internal layers
c       idmap(NBDY) I*4 - mapping into ibdmap used with layer denial
c-----
        integer*4 ibdmap(NBDY)
        logical permit(NBDY)
        integer*4 idmap(NBDY)

c-----
c       command line information
c-----
        integer*4 mxseg
        logical dop, dosv, dosh, doconv, dotran, dorefl
        common/debug/dodbg
        logical dodbg
c-----
c       model
c-----
        common/isomod/d(NLAY),a(NLAY),b(NLAY),rho(NLAY),
     1      qa(NLAY),qb(NLAY),
     1      etap(NLAY),etas(NLAY), frefp(NLAY), frefs(NLAY)
        common/depref/refdep
c-----
c       information on reverberations permitted
c-----
        integer*4 reverb(NLAY)
        data reverb/NLAY*1000/

c-----
c       parse the command line
c-----
        call gcmdln(mxseg,dop,dosv,dosh,doconv, dotran, dorefl,
     1      mname,depths,depthr,deny,bounce)
c-----
c       get the earth model
c-----
        inquire(file=mname,exist=ext)
        if(.not. ext)go to 9000
        l = lgstr(mname)

        write(LOT,*)'Model name: ',mname(1:l)
        call getmod(1,mname,mmax,title,iunit,iiso,iflsph,
     1      idimen,icnvel,ierr,.true.)
        write(LOT,*)'title:',title
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
c       open file for the ray description file
c-----
        open(2,file='genray96.ray',form='formatted',access='sequential',
     1      status='unknown')
        rewind 2
c-----
c       Enter the source and receiver depths
c-----
        write(LOT,*)'Source   Depth:',depths
        write(LOT,*)'Receiver Depth:',depthr
        depths = depths + refdep
        depthr = depthr + refdep
        write(2,'(a)')mname
        write(2,*)depths,depthr
c-----
c       get the layer index for the source and receiver in the model
c-----
        call getlyr(depths,idxsrc,mmax,lsrcbd,dsrcly)
        call getlyr(depthr,idxrec,mmax,lrecbd,drecly)
        write(LOT,10)'Source   :',depths,idxsrc,lsrcbd,dsrcly,
     1           'Receiver :',depthr,idxrec,lrecbd,drecly
   10   format(/' ',10x,'     Depth','     Layer',' Boundary '
     1      ,'PosInLayer'/
     2      ' ',a10,f10.3,i10,l10,f10.3/
     3      ' ',a10,f10.3,i10,l10,f10.3)
c-----
c       Enter the type of ray description desired
c-----
        do 1010 i=1,NBDY
            permit(i) = .true.
 1010   continue
        if(deny.ne.' ')then
            call denial(permit,mmax,deny)
        endif
        if(bounce.ne.' ')then
            call bnce(reverb,mmax,bounce)
        endif
c-----
c       depending on the answer construct the working boundary model
c-----
c-----
c       insert the boundaries into the model
c-----  
        call insert(mmax,mbdmax,idxsrc,idxrec,
     1      dsrcly,drecly,ibdmap,permit,idmap,mwork,msrc,mrec,
     2      mrsrc,mrrec)
        
c-----
c       at this point the source and receiver positions have
c       been placed remapped into artificial layer boundaries
c-----
        call pinitf('GPREP96.PLT')
        call pltmod(mbdmax,msrc,mrec,ibdmap,mname,title,
     1      dsrcly,drecly,permit,
     1      depths,depthr,reverb)
        call gwidth(0.05)
c-----
c       define the minimum number of segments, which is
c       related to the minimum required for a direct ray
c-----
        minseg = abs(mrsrc - mrrec) 
        dcolor = 100.0 /(mxseg - minseg + 1)
        do 4000 maxseg = minseg,mxseg
            kolor = 1000 + (maxseg-minseg + 1)*dcolor
            call newpen(kolor)
            call doray2(mwork,mrsrc,mrrec,maxseg,idmap,mbdmax,
     1          ibdmap,depths,depthr,
     2          dop,dosh,dosv,doconv,dotran, dorefl,
     3          reverb)
 4000   continue
        call gwidth(0.0)
        call newpen(1)
        call pend()
 9000   continue
        close (2)
        end

        subroutine pltmod(mmax,msrc,mrec,ibdmap,mname,title,
     1      dsrcly,drecly,permit,depths,depthr,reverb)
c-----
c       mmax    I*4 - maximum number of intefaces
c       msrc    I*4 - index of the source layer
c       mrec    I*4 - index of the receiver layer in the augmented model
c       ibdmap(NBDY)    I*4 - mapping array from original
c                   layers to boundaries
c                   The boundaries have been incremented
c                   by the source and receiver
c       mname   C*80    - name of the model file
c       title   C*80    - title of the model file
c       dsrcly  R*4 - depth of source within this layer
c       drecly  R*4 - depth of receiver in this layer
c       permit  L   - if .false. deny reflections and conversions
c                   This is used to simulate gradients
c                   in internal layers
c       depths  R*4 - source depth
c       depthr  R*4 - receiver depth
c-----
        integer NLAY, NBDY
        parameter (NLAY=200, NBDY=NLAY+2)
        integer*4 mmax, msrc, mrec
        integer*4 ibdmap(NBDY)
        character mname*80, title*80
        real*4 dsrcly, drecly
        logical permit(NBDY)
        real*4 depths, depthr

        common/isomod/d(NLAY),a(NLAY),b(NLAY),rho(NLAY),
     1      qa(NLAY),qb(NLAY),
     1      etap(NLAY),etas(NLAY), frefp(NLAY), frefs(NLAY)
        common/depref/refdep

c-----
c       information on reverberations permitted
c-----
        integer*4 reverb(NLAY)

        y0 = 7.0
        x0 = 0.8
        dy = (y0 - 1.0)/(mmax-1)
        xl = x0
        xh = xl + 6.0
        
        l = lgstr(title)
        xc = 0.5*(xh + xl - l*0.20)
        call symbol(xc,y0+0.40,0.20,title(1:l),0.0,l)
        call symbol(xh-0.20,y0+0.40,0.10,'DEPTH SRC=',0.0,10)
        call number(999.0  ,y0+0.40,0.10,depths-refdep,0.0,3)
        call symbol(999.0  ,y0+0.40,0.10,' REC=',0.0,5)
        call number(999.0  ,y0+0.40,0.10,depthr-refdep,0.0,3)

        call symbol(xh+0.10,y0+0.10,0.10,'H'  ,0.0,1)
        call symbol(xh+1.00,y0+0.10,0.10,'VP' ,0.0,2)
        call symbol(xh+1.90,y0+0.10,0.10,'VS' ,0.0,2)
        call symbol(xh+2.80,y0+0.10,0.10,'RHO',0.0,3)
        ibd = 0
        dep = 0.0
        do 100 i=1,mmax
            yy = y0 - (i-1)*dy
            yyy = yy - 0.5*dy
            call plot(xl,yy,3)
            if(ibdmap(i).gt.0)then
                ii = ibdmap(i)
c-----
c               This is a real layer boundary
c-----
                if(permit(i))then
                    call plotd(xh,yy,21,0.05)
                else
                    call plotd(xh,yy,16,0.05)
                endif
                call plotd(xh,yy,21,0.05)
                dm1 = dep
                dep = dep + d(ii)
                if(refdep.ge.dep)then
                call number(xh+0.10,yyy,0.10,-d(ii),0.0,3)
                else
                call number(xh+0.10,yyy,0.10, d(ii),0.0,3)
                endif
                call number(xh+1.00,yyy,0.10,a(ii)  ,0.0,3)
                call number(xh+1.90,yyy,0.10,b(ii)  ,0.0,3)
                call number(xh+2.80,yyy,0.10,rho(ii),0.0,3)
                ibd = ibd + 1
                call number(x0-0.5,yy-0.05,0.10,float(ibd),0.0,-1)
                if(i.ne.mmax)then
                    call number(xl-0.20,yyy,0.10,
     1                  float(reverb(ii)),0.0,-1)
                endif
            
            else
c-----
c               This is an artificial boundary for the source/receiver
c-----
                call plotd(xh,yy,7,0.05)
                if(ibdmap(i).lt.-1000)then
                    call symbol(x0-0.5,yy-0.05,0.10,'REC',0.0,3)
                    call symbol(xh+0.1,yyy,0.10,'(',0.0,1)
                    call number(999.0 ,yyy,0.10,drecly,0.0,3)
                    call symbol(999.0 ,yyy,0.10,')',0.0,1)
                else
                    call symbol(x0-0.5,yy-0.05,0.10,'SRC',0.0,3)
                    call symbol(xh+0.1,yyy,0.10,'(',0.0,1)
                    call number(999.0 ,yyy,0.10,dsrcly,0.0,3)
                    call symbol(999.0 ,yyy,0.10,')',0.0,1)
                endif
            endif
  100   continue
        yy = y0 - (msrc-1)*dy
        call symbol(xl,yy,0.10,char(1),0.0,-1)
        yy = y0 - (mrec-1)*dy
        call symbol(xh,yy,0.10,char(1),0.0,-1)
c-----
c       annotate the plot with description of line segments
c-----
        l = lgstr(mname)
        xc = x0 + 0.5*(xh - xl + l*0.20)
        call symbol(6.2,0.6,0.07,mname(1:l),0.0,l)

        call plot  (x0     ,0.60,3)
        call plotd (x0+1.00,0.60, 7,0.05)
        call symbol(x0+1.1,0.55,0.10,'ARTIFICIAL BOUNDARY ',0.0,20)
        call symbol(999.0 ,0.55,0.10,'FOR SOURCE-RECEIVER',0.0,19)
        call plot  (x0    ,0.40,3)
        call plotd (x0+1.0,0.40,21,0.05)
        call symbol(x0+1.1,0.35,0.10,'LAYER BOUNDARY'     ,0.0,14)
        call plot  (x0    ,0.20,3)
        call plotd (x0+1.0,0.20,16,0.05)
        call symbol(x0+1.1,0.15,0.10,
     1      'LAYER BOUNDARY (NO CONVERSIONS)',0.0,31)
        return
        end

        subroutine doray2(mmax,msrc,mrec,maxseg,idmap,mbdmax,
     1      ibdmap,depths,depthr,
     2      dop,dosh,dosv,doconv,dotran, dorefl,
     3          reverb)
c-----
c       mmax    I*4 - maximum number of boundaries in modified model
c       msrc    I*4 - boundary index for the source
c       mrec    I*4 - boundary index for the receiver
c       maxseg  I*4 - maximum number of components for the ray
c       idmap(NBDY) I*4 - mapping into ibdmap used with layer denial
c       mbdmax  I*4 - number of augmented boundaries, used for ploting
c       ibdmap(NBDY)    I*4 - mapping array from original
c                       layers to boundaries
c                       The boundaries have been incremented
c                       by the source and receiver
c-----
        integer NLAY, NBDY
        parameter (NLAY=100, NBDY=NLAY+2)
        integer*4 idmap(NBDY)
        integer*4 ibdmap(NBDY)

        logical dop, dosh, dosv,doconv, dotran, dorefl
        integer*4 reverb(NLAY)
        logical doplot
        integer*4 mmax

c-----
c       LIN I*4 - logical unit for standard input
c       LOT I*4 - logical unit for standard output
c-----
        integer LIN, LOT, LER
        parameter (LIN=5,LOT=6,LER=0)

        common/debug/dodbg
        logical dodbg

c-----
c       id(NBDY)    I*4 - Ray description
c-----
        integer id(NBDY)

        integer*4 gtbit2
        integer*4 ib(NBDY)

        write(LOT,*)'MAX_BDY=',mmax,' MSRC=',msrc,' MREC=',mrec, 
     1      ' MAXSEG=',maxseg
c-----
c       process over the maximum number of ray segments
c-----
        ilow = -1
        iup = 2**maxseg
        i = -1
 2000   continue
            i = i + 1
            if(i.lt.ilow .or. i.ge.iup)return
c-----
c           parse ray geometry
c-----
            id(1) = msrc
            do 2100 j=1,maxseg 
                iud = gtbit2(i,maxseg-j)
                ib(j) = iud
                if(iud.eq.0)then
                    id(j+1) = id(j) - 1
                else
                    id(j+1) = id(j) + 1
                endif
                if(id(j+1).lt.1 )then
                    go to 2000
                endif
                if(id(j+1).gt.mmax)then
                    go to 2000
                endif
 2100       continue
            if(id(maxseg+1).eq.mrec)then
c-----
c               The ray terminates at the receiver, now ignore
c               any direction changes at the source or receiver
c               layer that do not appear at the ends
c-----
                if(ichndr(id,maxseg,msrc,mrec,i).lt.0)go to 2000
c-----
c               In the special case of a receiver at the surface
c               do not permit the free surface reflection
c-----
                if(ircsrf(id,maxseg,depthr,i).lt.0)go to 2000
                if( depthr.eq.0.0.and.id(maxseg).eq.1)go to 2000
c-----
c               this is a successful ray
c-----
                if(dodbg)then
                write(LOT,'(a,i5,a,20i2)')
     1              ' ray             ',i,'*****',
     1              (id(j),j=1,maxseg+1)
                write(LOT,'(a,i5,a,20i2)')
     1              '  idmap(i)       ',i,'*****',
     1              (idmap(id(j)),j=1,maxseg+1)
                write(LOT,'(a,i5,a,20i6)')
     1              ' ibdmap(idmap(i))',i,'*****',
     1              (ibdmap(idmap(id(j))),j=1,maxseg+1)
                endif
c-----
c               now convert back to layers and process for
c               different wavetypes, not the number of types
c               is maxseg 
c-----
                mxseg1 = maxseg+1
                call remap2(id,mxseg1,idmap,ibdmap,
     1              msrc,mrec,
     2              depths,depthr,
     3              dop,dosv,dosh,doconv,dotran, dorefl,
     4              doplot,reverb,mmax)
                if(doplot)call pltray(id,maxseg+1,idmap,mbdmax)
                write(LOT,*)' '
            endif
        go to 2000
        end

        function gtbit2(i,j)
        integer gtbit2
c-----
c       In base 2, return the value of bit j of the number i
c       i   I*4 - number
c       j   I*4 - bit position 0,...,N-1
c-----
            k = 2**j
            gtbit2 = mod(i/k , 2)
        return
        end

        subroutine pltray(id,mseg,idmap,mbdmax)
        integer NLAY, NBDY
        parameter (NLAY=100, NBDY=NLAY+2)
        integer*4 idmap(NBDY)
        integer*4 id(1)
        x0 = 0.8
        y0 = 7.0
        dy = (y0 - 1.0)/(mbdmax-1)
        dx = 6.0/(mseg-1)
        ipen = 3
        do 200 i=1,mseg
            ii = idmap(id(i))
            yy = y0 - (ii-1)*dy
            xx = x0 + (i-1)*dx
            call plot(xx,yy,ipen)
            ipen = 2
  200   continue
        return
        end


        subroutine iyesno(iyes)
        logical iyes
c-----
c       LIN I*4 - logical unit for standard input
c       LOT I*4 - logical unit for standard output
c-----
        integer LIN, LOT, LER
        parameter (LIN=5,LOT=6,LER=0)
        character ans*2

 1000   continue
            read(LIN,'(a)')ans
            if(ans(1:1).eq.'y' .or. ans(1:1).eq.'Y')then
                iyes = .true.
                return
            else if(ans(1:1).eq.'n' .or. ans(1:1).eq.'N')then
                iyes = .false.
                return
            endif
            write(LOT,*)'Enter (yn) only'
        go to 1000
        end
        
        subroutine getlyr(depthp,idxpnt,mmax,lbdy,dlyr)
        real*4 depthp, dlyr
        integer*4 idxpnt, mmax
        logical lbdy

        integer NLAY
        parameter (NLAY=200)
        common/isomod/d(NLAY),a(NLAY),b(NLAY),rho(NLAY),
     1      qa(NLAY),qb(NLAY),
     1      etap(NLAY),etas(NLAY), frefp(NLAY), frefs(NLAY)
        common/depref/refdep

c-----
c       Determine the layer in which the depth places a point
c       This is used to find the position of the source and receiver
c
c       depthp  R*4 - Desired depth
c       idxpnt  I*4 - Layer index
c       mmax    I*4 - Maximum number of layers in a model
c       lbdy    L   - .true. depth is on a layer boundary
c       dlyr    R*4 - depth in layer
c
c       NOTE. In case the point is at a layer boundary, it will be 
c       placed at the top of the layer beneath
c-----
        lbdy = .false.
        if(depthp .eq. 0.0)then
            idxpnt = 1
            dlyr = 0.0
            lbdy = .true.
        else
            dtop = 0.0
            do 100 i=1,mmax-1
                dbot = dtop + d(i)
                if(depthp .ge.dtop .and. depthp.lt.dbot)then
                    idxpnt = i
                    if(depthp .eq.dtop)then
                        lbdy = .true.
                    else
                        lbdy = .false.
                    endif
                    dlyr = depthp - dtop
                    return
                endif
                dtop = dbot
  100       continue
c-----
c           If get to this point, the depth point is in the halfspace
c           so add a layer
c-----
            d(mmax) = depthp - dtop
c           d(mmax+1) = 0.0
c           a(mmax+1) = a(mmax)
c           b(mmax+1) = b(mmax)
c           rho(mmax+1) = rho(mmax)
c           qa(mmax+1) = qa(mmax)
c           qb(mmax+1) = qb(mmax)
c           etap(mmax+1) = etap(mmax)
c           etas(mmax+1) = etas(mmax)
c           mmax = mmax + 1
c           dlyr = 0.0
        dlyr = depthp - dtop
            lbdy = .true.
            idxpnt = mmax
        endif
        return
        end

        subroutine insert(mmax,mbdmax,idxsrc,idxrec,
     1      dsrcly,drecly,ibdmap,permit,idmap,mwork,msrc,mrec,
     2      mrsrc,mrrec)
c-----
c       mmax    I*4 - number of layers in the original model
c       mbdmax  I*4 - number of reflecting boundaries
c       idxsrc  I*4 - layer index for source 
c       idxrec  I*4 - layer index for receiver 
c       dsrcly  R*4 - depth of source within this layer
c       drecly  R*4 - depth of receiver in this layer
c       ibdmap  I*4 - mapping array from boundary to layer index
c       permit  L   - if .false. deny reflections and conversions
c                       This is used to simulate gradients
c                       in internal layers
c       idmap(NBDY) I*4 - mapping into ibdmap used with layer denial
c       mwork       I*4 - number of working boundaries in reduced model
c       msrc        I*4 - boundary index of source in model
c       mrec        I*4 - boundary index of receiver in model
c       mrsrc       I*4 - boundary index of source in reduced model
c       mrrec       I*4 - boundary index of receiver in reduced model
c-----
        integer NLAY, NBDY
        parameter (NLAY=200, NBDY=NLAY+2)
        common/isomod/d(NLAY),a(NLAY),b(NLAY),rho(NLAY),
     1      qa(NLAY),qb(NLAY),
     1      etap(NLAY),etas(NLAY), frefp(NLAY), frefs(NLAY)
        common/depref/refdep


        integer*4 ibdmap(NBDY)
        logical permit(NBDY)
        integer*4 idmap(NBDY)

c-----
c       LIN I*4 - logical unit for standard input
c       LOT I*4 - logical unit for standard output
c-----
        integer LIN, LOT, LER
        parameter (LIN=5,LOT=6,LER=0)

        common/debug/dodbg
        logical dodbg

        logical srcrec
c-----
c       first check the trivial case of the source or receiver
c-----
c-----
c       now do the mapping - use a negative number to 
c       indicate an artificial boundary
c-----

c-----
c       initialize to the original model
c-----
        do 1000 i=1,mmax
            ibdmap(i) = i
 1000   continue
        mbdmax = mmax
        if(dodbg)then
            write(LOT,*)'IBDMAP'
            write(LOT,*)(ibdmap(i),i=1,mbdmax)
        endif

c-----
c       now add the artificial boundary for the source and
c       receiver. The extra care is to ensure that the plot
c       looks nice as well as the model listing
c-----      
        if(idxsrc .eq. idxrec)then
            if(dsrcly .gt. drecly)then
                call putin(ibdmap,permit,idxsrc,mbdmax,.true.)
                call putin(ibdmap,permit,idxrec,mbdmax,.false.)
            else
                call putin(ibdmap,permit,idxrec,mbdmax,.false.)
                call putin(ibdmap,permit,idxsrc,mbdmax,.true.)
            endif
        else
            call putin(ibdmap,permit,idxsrc,mbdmax,.true.)
            call putin(ibdmap,permit,idxrec,mbdmax,.false.)
        endif
c-----
c       now verify the mapping
c-----
c-----
c       now output the model according to the boundaries
c       in addition compute the reduced model
c
c       The printed output at the end of the do 5000 is not
c       necessary for the ray tracing. It does verify the
c       layer insertion code. The special case of the source
c       and receiver in the same layer, requires some added
c       care
c-----
        if(idxsrc.eq.idxrec)then
            d1 = min(dsrcly,drecly)
            d2 = max(dsrcly,drecly) - d1
            d3 = d(idxsrc)  - max(dsrcly,drecly)
            srcrec = .true.
            kount = 0
        else
            srcrec = .false.
        endif
            

    1   format(' BDRY      LAYER             H      P-VEL     S-VEL   ',
     1      'DENSITY CONV  '/
     2  ' ____________________________________________________________',
     3  '_____')
    2   format(' ',5x, 5x,i5,5x,4f10.3,1x)
    3   format(' ',i5, 1x, a4, 50x,1x,l1)
    4   format(' ',
     1  ' _____________________________________________________________'
     3  ,'____'/
     2  ' Working Boundaries=',i5/
     3  ' Source Boundary   =',i5/
     4  ' Receiver Boundary =',i5)
        write(LOT,1)
c-----
c       mwork is the working boundary index for those boundaries that
c       permit reflections
c-----
        mwork = 0
        do 5000 ibdy=1,mbdmax
            if(ibdmap(ibdy).gt.0)then
                write(LOT,3)ibdy,'----', permit(ibdy)
            else if(ibdmap(ibdy).lt.-1000)then
                write(LOT,3) -ibdy,'REC ', permit(ibdy)
            else
                write(LOT,3) -ibdy,'SRC ', permit(ibdy)
            endif
            if(permit(ibdy) .or. ibdmap(ibdy).lt.0)then
                mwork = mwork + 1
                idmap(mwork) = ibdy
            endif
            if(ibdmap(ibdy).lt.0)then
                if(ibdmap(ibdy) .lt. -1000)then
                    mrrec = mwork
                    mrec = ibdy
                else
                    mrsrc = mwork
                    msrc = ibdy
                endif
            endif
            i = abs(ibdmap(ibdy))
            if(ibdmap(ibdy).lt.0.0)then
                if(i.gt.1000)then
                    j = i - 1000
                    dd = d(j) - drecly
                else
                    j = i
                    dd = d(j) - dsrcly
                endif
                if(srcrec)then
                    if( kount .eq.1)then
                        dd = d2
                    else if(kount .eq. 2)then
                        dd = d3
                    endif
                    kount = kount + 1
                endif
                write(LOT,2)j,dd,a(j),b(j),rho(j)
            else
                if(i.eq.idxsrc)then
                    dd = dsrcly
                else if(i.eq.idxrec)then
                    dd = drecly
                else
                    dd = d(i)
                endif
                if(srcrec .and. kount.eq.0)then
                    dd = d1
                    kount = kount + 1
                endif
                write(LOT,2)i,dd,a(i),b(i),rho(i)
            endif
 5000   continue
        write(LOT,4)mwork, mrsrc, mrrec
        write(LOT,*)(idmap(i),i=1,mwork)
        return
        end

        subroutine denial(permit,mmax,deny)
c-----
c       Interactive routine for user to reduce the number of
c       ray conversions in internal layers
c
c       permit(NBDY)    L   - 
c             if .false. deny reflections and conversions
c                       This is used to simulate gradients
c                       in internal layers
c       mmax        I*4 - maximum number of layers in model
c       deny        C*80    - data file with conversion denial
c-----
        integer LIN, LOT, LER
        parameter (LIN=5,LOT=6,LER=0)

        integer NLAY, NBDY
        parameter (NLAY=200, NBDY=NLAY+2)
        common/isomod/d(NLAY),a(NLAY),b(NLAY),rho(NLAY),
     1      qa(NLAY),qb(NLAY),
     1      etap(NLAY),etas(NLAY), frefp(NLAY), frefs(NLAY)
        common/depref/refdep

        logical permit(NBDY)
        character deny*80

        logical ext

c-----
c       'The purpose is to specify the layers in which no'
c       'conversions are permitted at layer boundaries. '
c       'This means that there will be no reflected multiples'
c       'solely within the layer. For example, if there are'
c       '3 layers over a halfspace, and multiples are denied'
c       'in layer 2, then a ray path 1 2 2 3 2 1 is '
c       'prohibited, but 1 2 3 2 1 is permitted.'
c       'In addition, there will be no conversion between'
c       'P and SV on transmission or reflection.'
c       ' '
c       'This permits the user to specify rays of interest'
c       ' '
c       'Note that reflections will always be permitted  '
c       'off the top  model boundaries.'
c
c       'This model has ',mmax+1, ' boundaries'
c-----

        inquire(file=deny,exist=ext)
        if(.not. ext)return
        open(1,file=deny,status='unknown',form='formatted',
     1          access='sequential')
        rewind 1

 1000   continue
        read(1,*,end=9999,err=1000)i
            if(i.le.0)return
            if(i.ge. 1 .and. i.le.mmax+1)then
                write(LOT,*)i
                permit(i) = .false.
            endif
        go to 1000
 9999   continue
        close(1)
        return
        end
            
        subroutine bnce(reverb,mmax,bounce)
c-----
c       Interactive routine for user to reduce the number of
c       ray conversions in internal layers
c
c       permit(NBDY)    L   - 
c            if .false. deny reflections and conversions
c                       This is used to simulate gradients
c                       in internal layers
c       mmax        I*4 - maximum number of layers in model
c       deny        C*80    - data file with conversion denial
c-----
        integer LIN, LOT, LER
        parameter (LIN=5,LOT=6,LER=0)

        integer NLAY, NBDY
        parameter (NLAY=100, NBDY=NLAY+2)

        integer*4 reverb(NLAY)
        character bounce*80

        logical ext


        inquire(file=bounce,exist=ext)
        if(.not. ext)return
        open(1,file=bounce,status='unknown',form='formatted',
     1          access='sequential')
        rewind 1

        write(LOT,*)'Layer Maximum_Reverberations'
 1000   continue
        read(1,*,end=9999,err=1000)i,j
c-----
c       layer i, number of bounces j
c-----
            if(i.le.0)return
            if(i.ge. 1 .and. i.le.mmax)then
                write(LOT,*)i,j
                reverb(i) = j
            endif
        go to 1000
 9999   continue
        close(1)
        return
        end
            

        function ichndr(id,maxseg,msrc,mrec,iray)
c-----
c       id  I*4 - array of connected boundaries
c       maxseg  I*4 - number of segments, the end is maxseg+1
c       msrc    I*4 - source boundary index
c       mrec    I*4 - receiver boundary index
c       iray    I*4 - ray index
c
c       ichndr      +1 keep this ray
c               -1 drop this ray
c-----
        integer NLAY, NBDY
        parameter (NLAY=100, NBDY=NLAY+2)
        integer id(NBDY)
        
c-----
c       LIN I*4 - logical unit for standard input
c       LOT I*4 - logical unit for standard output
c-----
        integer LIN, LOT, LER
        parameter (LIN=5,LOT=6,LER=0)

        common/debug/dodbg
        logical dodbg


        ichndr = 1
c-----
c       look at the Laplacian
c
c       if ray index is 2 3 4 then Laplacian is 0
c       a 2 3 2 gives non zero and represents a reflection
c       If the center element if a source or receiver layer, then
c       this is artificial
c-----
        do 1000 i=2,maxseg
            ilap = id(i-1) - 2*id(i) + id(i+1)
            if(ilap .ne. 0 .and. id(i).eq.msrc)then
                ichndr = -1
            endif
            if(ilap .ne. 0 .and. id(i).eq.mrec)then
                ichndr = -1
            endif
 1000   continue
        if(ichndr .lt. 0 .and. dodbg)then
            write(LOT,'(a,i10,a,20i2)')' ray',iray,' XXXXX',
     1          (id(j),j=1,maxseg+1)
        endif
        return
        end
            
        function ircsrf(id,maxseg,depthr,iray)
c-----
c       For the special case of the receiver at the
c       surface the generalized ray will use the free surface
c       Zoepritz coefficients rather than adding together the
c       direct and reflected rays
c
c       So if the receiver is at the surface and the previous
c       boundary is 1, then ignore the ray
c-----
c       id  I*4 - array of connected boundaries
c       maxseg  I*4 - number of segments, the end is maxseg+1
c       depthr  R*4 - receiver depth
c       iray    I*4 - ray index
c
c       ircsrf  I*4 - 1 use this ray
c                 -1 drop this ray
c-----
        integer NLAY, NBDY
        parameter (NLAY=100, NBDY=NLAY+2)
        integer id(NBDY)

c-----
c       LIN I*4 - logical unit for standard input
c       LOT I*4 - logical unit for standard output
c-----
        integer LIN, LOT, LER
        parameter (LIN=5,LOT=6,LER=0)


        common/debug/dodbg
        logical dodbg

        ircsrf = 1
        if( depthr .eq. 0.0 .and. id(maxseg).eq.1)ircsrf = -1
        if(ircsrf .lt. 0)then
            if(dodbg)then
            write(LOT,'(a,i10,a,20i2)')' ray',iray,' YYYYY',
     1          (id(j),j=1,maxseg+1)
            endif
        endif
        return
        end

        subroutine remap2(id,mseg,idmap,ibdmap,
     1              msrc,mrec,
     2              depths,depthr,
     3              dop,dosv,dosh,doconv,dotran, dorefl,
     4              doplot,reverb,mmax)
c-----
c       LIN I*4 - logical unit for standard input
c       LOT I*4 - logical unit for standard output
c-----
        integer LIN, LOT, LER
        parameter (LIN=5,LOT=6,LER=0)

        integer NLAY, NBDY
        parameter (NLAY=100, NBDY=NLAY+2)
        integer*4 id(NBDY)
        integer mseg
        integer*4 idmap(NBDY), ibdmap(NBDY)
        integer  msrc, mrec
        real*4 depths, depthr
        logical dop, dosv, dosh,doconv, dotran, dorefl
        logical doplot
        integer*4 reverb(NLAY)
        integer*4 mmax

        common/debug/dodbg
        logical dodbg

        integer i1(100), i2(100)
        integer j1(100), j2(100)
        integer*4 treverb(NLAY)

        integer gtbit2

        logical srcup, outref, outtrn
        logical isfld
c-----
c       first determine if both the source and the receiver are in the
c       same layer. if they are, then determine if the ray is up or down
c-----
        if(id(2) .gt. id(1))then
            srcup = .false.
        else
            srcup = .true.
        endif
c-----
c       go through the id list and drop out any source or receiver
c       index that is not the first or last
c-----
        j = 1
        i1(j) = idmap(msrc)
        do 1000 i=2,mseg-1
            k = idmap(id(i))
            if(k.ne. idmap(msrc) .and. k.ne.idmap(mrec))then
                j = j + 1
                i1(j) = k
            endif
 1000   continue
        j = j + 1
        i1(j) = idmap(mrec)
        jmax = j
c-----
c       at this stage we have the unique elements 
c       now put in the ray type
c-----
        do  987 j=1,jmax
            j1(j) = j
  987   continue
        if(dodbg)then
            write(LOT,*)'jmax=',jmax,' I1=',(i1(i),i=1,jmax)
            write(LOT,*)'     ',jmax,' J1=',(j1(i),i=1,jmax)
        endif
c-----
c       now fill in the missing elements
c-----
        j = 1
        k = 1
        i2(k) = i1(1)
        j2(k) = j1(1)
 1100   continue
            ilw = i1(j)
            iup = i1(j+1)
            jold = j1(j)
            if(ilw.gt.iup)then
                inc = -1
            else
                inc =  1
            endif
            do 1200 i=ilw+inc,iup,inc
                k = k + 1
                i2(k) = i
                j2(k) = jold 
 1200       continue
            j = j + 1
            if(j.eq.jmax)go to 1300
        go to 1100
 1300   continue
        kmax = k
        if(dodbg)then
            write(LOT,*)'kmax=',kmax,' I2=',(i2(i),i=1,kmax)
            write(LOT,*)'     ',kmax,' J2=',(j2(i),i=1,kmax)
        endif
c-----
c       now repurge the source and receiver boundaries on c
c           interior of ray
c-----
        j = 1
        i1(j) = ibdmap(i2(1))
        j1(j) = j2(1)
        do 1400 i=2,kmax-1
            k = ibdmap(i2(i))
            if(k .gt. 0)then
                j = j + 1
                i1(j) = k
                j1(j) = j2(i)
            endif
 1400   continue
        j = j + 1
        i1(j) = ibdmap(i2(kmax))
        j1(j) = j2(kmax)
        jmax = j
        if(dodbg)then
            write(LOT,*)'jmax=',jmax,' I1=',(i1(i),i=1,jmax)
            write(LOT,*)'     ',jmax,' J1=',(j1(i),i=1,jmax)
        endif
c-----
c       now convert from boundaries to layers
c-----
        if(jmax.eq.2 .and. depths.lt.depthr)then
            ndeg = -1
        else
            ndeg = 1
        endif
c-----
c       first and last are special because of the partial segment in
c       layer
c
c       first ray is always in the layer of the source
c-----
        i2(1) = abs(i1(1))
        j2(1) = j1(2)

            do 2000 j=2,jmax-2
                if(i1(j+1).lt.i1(j))then
                    i2(j) = i1(j+1)
                else
                    i2(j) = i1(j)
                endif
                j2(j) = j1(j+1)
 2000       continue
        if(jmax.gt.2)then
            i2(jmax -1) = abs(i1(jmax)+1000)
            j2(jmax -1)    = j1(jmax)
        endif
        nn = jmax - 1
c-----
c       now that the ray description is defined, determine whether the
c       source and receiver are in the same layer, and if the number
c       of segments is greater than 1, determine if the ray went
c       upward or downward
c-----
        if(nn.gt.1 .and. i2(1).eq.i2(nn))then
            if(srcup)then
                ndeg = -ndeg
            endif
        else if(nn.eq.1)then
                ndeg = 1
        endif
            
        if(dodbg)then
            write(LOT,113) nn, ndeg
            write(LOT,114) (i2(j),j=1,nn)
            write(LOT,114) (j2(j),j=1,nn) 
            write(LOT,*)' '
        endif
c-----
c       Now we have the layers through which the rays propagate, e.g.,
c       the i2 array, and the ray type variable, the j2 values
c
c       Now test for the maximum number of reverberations in each layer
c       If that test is passed, then output the rays, else return
c-----
        doplot = .false.
        do 2100 i=1,mmax
            treverb(i) = reverb(i)
 2100   continue
        do 2150 j=1,nn
            k = treverb(i2(j))
            treverb(i2(j)) = k - 1
 2150   continue
c-----
c       now check to see if the treverb < 0 in any layer, in which case
c       the layer has too many ray segments
c-----
        do 2160 i=1,mmax
            if(treverb(i) .lt. 0)return
 2160   continue
        doplot = .true.
c-----
c       Now output ray descriptions for all rays
c-----
c-----
c       now fill SH
c-----
        if(dosh)then
            call fillit(j1,nn,4)
            call chkfld(mmax,nn,i2,j1,isfld)
            if(.not. isfld)then
                write(2,103)nn, ndeg
                write(2,104) (i2(j),j=1,nn)
                write(2,104) (j1(j),j=1,nn)
                write(LOT,113)nn, ndeg
                write(LOT,114) (i2(j),j=1,nn)
                write(LOT,114) (j1(j),j=1,nn) 
            endif
        endif
c-----
c       now fill P-SV
c-----
        if(dop .and. dosv .and. doconv)then
            maxpow = j2(nn)
            iup = 2**maxpow
            do 3000 i=0,iup-1
                do 3100 j=1,nn
                    ibt = j2(j)-1
                    ips = gtbit2(i,ibt)
                    if(ips.gt.0)then
                        j1(j) = 3
                    else
                        j1(j) = 5
                    endif
 3100           continue
                outref = .true.
                outtrn = .true.
c-----
c               if conversions are desired only on
c               reflection then only output those
c-----
                if(dorefl)then
                    call chkref(nn,i2,j1,outref)
                endif
c-----
c               if conversions are desired only on
c               transmission then only output those
c-----
                if(dotran)then
                    call chktrn(nn,i2,j1,outtrn)
                endif
                if(outtrn .and. outref)then
                    call chkfld(mmax,nn,i2,j1,isfld)
                    if(.not. isfld)then
                        write(2,103)nn, ndeg
                        write(2,104) (i2(j),j=1,nn)
                        write(2,104) (j1(j),j=1,nn)
                        write(LOT,113)nn, ndeg
                        write(LOT,114) (i2(j),j=1,nn)
                        write(LOT,114) (j1(j),j=1,nn) 
                    endif
                endif
 3000       continue
        else
            if(dop)then
                call fillit(j1,nn,5)
                call chkfld(mmax,nn,i2,j1,isfld)
                if(.not. isfld)then
                    write(2,103)nn, ndeg
                    write(2,104) (i2(j),j=1,nn)
                    write(2,104) (j1(j),j=1,nn)
                    write(LOT,113)nn, ndeg
                    write(LOT,114) (i2(j),j=1,nn)
                    write(LOT,114) (j1(j),j=1,nn) 
                endif
            endif
            if(dosv)then
                call fillit(j1,nn,3)
                call chkfld(mmax,nn,i2,j1,isfld)
                if(.not. isfld)then
                    write(2,103)nn, ndeg
                    write(2,104) (i2(j),j=1,nn)
                    write(2,104) (j1(j),j=1,nn)
                    write(LOT,113)nn, ndeg
                    write(LOT,114) (i2(j),j=1,nn)
                    write(LOT,114) (j1(j),j=1,nn) 
                endif
            endif
        endif
c-----
c       output to file does not require carriage control on MSDOS/UNIX
c-----
  103   format(2i5)
  104   format(5x,14i5) 
c-----
c       UNIX FORTRAN OUTPUT DOES NOT REQUIRE CARRIAGE CONTROL
c-----
c  113  format(2i5)
c  114  format(5x,14i5) 
c-----
c       MSDOS FORTRAN OUPUT REQUIRES CARRIAGE CONTROL
c-----
  113   format(2i5)
  114   format(5x,14i5) 
        return
        end

        subroutine gcmdln(mxseg,dop,dosv,dosh,doconv,dotran,dorefl,
     1      mname,depths,depthr,deny,bounce)
c-----
c       parse the command line arguments
c-----
c       mxseg   R*4 - upper limit for the maximum number of segments
c                   considered
c       dop L   - include P rays
c       dosh    L   - include SH
c       dosv    L   - include SV
c       doconv  L   - permit P-SV conversions
c       dorefl  L   - permit P-SV only on reflection 
c       dotran  L   - permit P-SV only on transmission 
c       mname   C*80    - name of model file
c       depths  R*4 - source depth
c       depthr  R*4 - receiver depth
c       deny    C*80    - file name of boundary interactions to deny
c       bounce  C*80    - file name of maximum internal bounces
c-----
        logical dop, dosh, dosv, doconv, dotran, dorefl
        character mname*80, deny*80, bounce*80

        common/debug/dodbg
        logical dodbg

        character name*40
        integer mnmarg
        mxseg = 12
        dop = .false.
        dosh = .false.
        dosv = .false.
        doconv = .false.
        dorefl = .false.
        dotran = .false.
        dodbg = .false.
        depths = 0.0
        depthr = 0.0
        deny = ' '
        mname = ' '
        bounce = ' '
        nmarg = mnmarg()
        i = 0
 1000   continue
            i = i + 1
            if(i.gt.nmarg)go to 2000
                call mgtarg(i,name)
                if(name(1:4).eq.'-DOP')then
                    dop = .true.
                else if(name(1:5).eq.'-DOSV')then
                    dosv = .true.
                else if(name(1:5).eq.'-DOSH')then
                    dosh = .true.
                else if(name(1:4).eq.'-DOA')then
                    dop = .true.
                    dosv = .true.
                    dosh = .true.
                else if(name(1:5).eq.'-DOCO')then
                    doconv = .true.
                else if(name(1:5).eq.'-DORE')then
                    dorefl = .true.
                    doconv = .true.
                else if(name(1:5).eq.'-DOTR')then
                    dotran = .true.
                    doconv = .true.
                else if(name(1:6).eq.'-DEBUG')then
                    dodbg = .true.
                else if(name(1:2).eq.'-M')then
                    i = i + 1
                    call mgtarg(i,mname)
                else if(name(1:5).eq.'-DENY')then
                    i = i + 1
                    call mgtarg(i,deny)
                else if(name(1:2).eq.'-R')then
                    i = i + 1
                    call mgtarg(i,bounce)
                else if(name(1:2).eq.'-N')then
                    i = i + 1
                    call mgtarg(i,name)
                    read(name,'(bn,i10)')mxseg
                    if(mxseg.lt.1 .or. 
     1                  mxseg.gt.31)mxseg=12
                else if(name(1:3).eq.'-HS')then
                    i = i + 1
                    call mgtarg(i,name)
                    read(name,'(bn,f10.0)')depths
                else if(name(1:3).eq.'-HR')then
                    i = i + 1
                    call mgtarg(i,name)
                    read(name,'(bn,f10.0)')depthr
                else if(name(2:2) .eq. '?')then
                    call usage()
                else if(name(2:2) .eq. 'h')then
                    call usage()
                endif
        go to 1000
 2000   continue
        if(mname .eq. ' ')call usage()
        if( .not. dop .and. .not. dosh .and. .not. dosv)call usage()
        return
        end

        subroutine usage()
        integer LER
        parameter (LER=0)
        write(LER,*)'Usage: gprep96 -M model  [-DOP] [-DOSV]',
     1      ' [-DOSH] [-DOALL] [-DOCONV] [-DOREFL] ',
     2      ' [-DOTRAN] [-DEBUG] [-DENY deny ]',
     3      ' [-R reverb] [-N maxseg] [-HS depths] [-HR depthr]'
        write(LER,*)
     1  '-M model  (default none )  Earth model file name'
        write(LER,*)
     1  '-N nseg   (default 12   )  Maximum number of ray segments '
        write(LER,*)
     1  '-DOP      (default false)  Generate P ray description'
        write(LER,*)
     1  '-DOSV     (default false)  Generate SV ray description'
        write(LER,*)
     1  '-DOSH     (default false)  Generate SH ray description'
        write(LER,*)
     1  '-DOALL    (default false)  Generate P, SV, SH ray description'
        write(LER,*)
     1  '-DOCONV   (default false)  Permit P-SV conversions'
        write(LER,*)
     1  '-DOREFL   (default false)  P-SV conversions on  reflection'
        write(LER,*)
     1  '-DOTRAN   (default false)  P-SV conversions on  transmission'
        write(LER,*)
     1  '-DENY     (default none )  file with layer conversion denial'
        write(LER,*)
     1  '-R reverb (default none )  file with maximum number ',
     2      'of multiples in layer'
        write(LER,*)
     1  '-HS depths(default 0.0  )  source depth '
        write(LER,*)
     1  '-HR depthr(default 0.0  )  receiver depth '
        write(LER,*)
     1  '-?        (default none )  this help message '
        write(LER,*)
     1  '-h        (default none )  this help message '
        stop
        end
        
        function issrrc(ibdy,issrc)
c-----
c       parse the boundary index to determine the
c       boundary and also whether the boundary is that
c       of the source or the receiver
c
c       ibdy    I*4 - index of the boudary
c       issrc   L   - .true. this is the source boundary
c       issrrc  I*4 - boundary index
c-----
        integer*4 ibdy
        logical issrc
        if(ibdy.lt.0)then
            if(ibdy.lt.-1000)then
                issrrc =  (1000-ibdy)
                issrc = .false.
            else
                issrrc = -ibdy
                issrc = .true.
            endif
        else
            issrrc = -1
        endif
        return
        end

        subroutine putin(ibdmap,permit,idx,mbdmax,issrc)
c-----
c       actually put in the extra boundary into the model
c-----
c       ibdmap  I*4 - mapping array from boundary to layer index
c       permit  L   - if .false. deny reflections and conversions
c                       This is used to simulate gradients
c                       in internal layers
c       idx I*4 - layer index for the source or reciever
c       mbdmax  I*4 - number of reflecting boundaries
c       issrc   L   - .true. is the source
c                 .false. is the receiver
c-----
c-----
c       LIN I*4 - logical unit for standard input
c       LOT I*4 - logical unit for standard output
c-----
        integer LIN, LOT, LER
        parameter (LIN=5,LOT=6,LER=0)

        integer NLAY, NBDY
        parameter (NLAY=100, NBDY=NLAY+2)
        integer*4 ibdmap(NBDY)
        logical permit(NBDY)
        integer*4 idx
        integer*4 mbdmax
        logical issrc

        common/debug/dodbg
        logical dodbg

        logical done

            done = .false.
            jjsrrc = mbdmax+1
            do 1100 i=1,mbdmax
                if(abs(ibdmap(i)).eq.idx .and. .not. done)then
                    jjsrrc = i
                    done = .true.
                endif
 1100       continue
            if(dodbg)then
                if(issrc)then
                    write(LOT,*)'jjsrc=',jjsrrc
                else
                    write(LOT,*)'jjrec=',jjsrrc
                endif
            endif
            do 1150 i=mbdmax,jjsrrc,-1
                ibdmap(i+1) = ibdmap(i)
                permit(i+1) = permit(i)
 1150       continue
            jjsrrc = jjsrrc + 1
            if(issrc)then
                ibdmap(jjsrrc) = - abs(ibdmap(jjsrrc))
            else
                ibdmap(jjsrrc) = - abs(ibdmap(jjsrrc)) - 1000
            endif
            permit(jjsrrc) = .true.
            mbdmax = mbdmax + 1

            if(dodbg)then
                write(LOT,*)'IBDMAP'
                write(LOT,*)(ibdmap(i),i=1,mbdmax)
            endif

        return
        end

        subroutine fillit(j1,nn,intval)
        integer NLAY, NBDY
        parameter (NLAY=100, NBDY=NLAY+2)
        integer*4  j1(NBDY)
            do 100 i=1,nn
                j1(i) = intval
 100        continue
        return
        end

        subroutine chkref(nn,i2,j1,output)
c-----
c       parse the ray description to permit only
c       ray conversion on transmission
c
c-----
c       nn  I*4 - Number of ray segments
c       i2  I*4 - layer in which the ray occurs
c       j1  I*4 - wave type SV(3) or P(5)
c       output  L   - .true. output the ray information
c                 .false. do not output
c-----
        integer*4 nn, i2(nn), j1(nn)
        logical output


        do 1000 i=2,nn
            if(i2(i-1).ne.i2(i))then
c-----
c           first deny conversion on transmission
c-----
                if( j1(i-1).ne.j1(i))then
                    output = .false.
                endif
c-----
c         deny conversion on upward reflection, except for free surface
c-----
            endif
 1000   continue
        return
        end

        subroutine chktrn(nn,i2,j1,output)
c-----
c       parse the ray description to permit only
c       ray conversion on transmission
c
c-----
c       nn  I*4 - Number of ray segments
c       i2  I*4 - layer in which the ray occurs
c       j1  I*4 - wave type SV(3) or P(5)
c       output  L   - .true. output the ray information
c                 .false. do not output
c-----
        integer*4 nn, i2(nn), j1(nn)
        logical output


        do 1000 i=2,nn
            if(i2(i-1).eq.i2(i))then
c-----
c         first deny conversion on transmission
c-----
                if( j1(i-1).ne.j1(i))then
                    output = .false.
                endif
c-----
c         deny conversion on upward reflection, except for free surface
c-----
            endif
 1000   continue
        return
        end

        subroutine chkfld(mmax,nn,i2,j1,isfld)
c-----
c       If a layer is a fluid, e.g., a zero S velocity, and if the
c       ray segment is S, then the ray cannot be used
c-----  
        
        integer NLAY
        parameter (NLAY=200)
        common/isomod/d(NLAY),a(NLAY),b(NLAY),rho(NLAY),
     1      qa(NLAY),qb(NLAY),
     1      etap(NLAY),etas(NLAY), frefp(NLAY), frefs(NLAY)
        common/depref/refdep


        integer i2(100), j1(100)
        logical isfld
        isfld = .false.
        do 100 i=1,nn

            lay = i2(i)
            if(b(lay) .le. 0.000001 .and. j1(i).ne.5)then
                isfld = .true.
                return
            endif
  100   continue
        return
        end

