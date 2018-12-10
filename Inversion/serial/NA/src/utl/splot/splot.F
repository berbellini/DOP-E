c------------------------------------------------------------------------
c
c	Program splot - plots projections of multidimensional points 
c			onto a 2-D plane. Reads points (models) from
c			Neighbourhood algorithm direct access file (NAD). 
c
c       Comments
c
c		ndmax 		: Maximum number of dimensions 
c				  including special parameters
c		nplmax 		: Maximum number of plots 
c		nxmax 		: Maximum number of x panels in a plot 
c		nymax 		: Maximum number of y panels in a plot 
c		ngridmax 	: Maximum size of plotting grid array
c
c	Calling sequence:
c
c		All arguments in brackets are optional.
c
c					M. Sambridge, RSES, May 1999.
c
c------------------------------------------------------------------------
c
	Program splot
c
#if NA_F77_COMPATIBILITY
        include         '../../NA_src/na_param.inc'

        parameter       (ndmax=nd_max,nhmax=nh_max,nemax=nmod_max)
#endif

	parameter	(nxmax = 100)
	parameter	(nymax = 100)
	parameter	(nplmax = 50)
        parameter       (ndismax = 100)
        parameter       (nmodmax = 20)
        parameter       (ngridmax = 20000)
c
        character*128   title
#if NA_F77_COMPATIBILITY

c                                                      NAD file
        real*4          models(ndmax*nemax)
        real*4          data(nemax)
        integer         rank(nemax)
        character       header(nhmax)
c                                                      frame arrays
	real*4          df(0:ndmax)
	real*4          ticl(0:ndmax)
        real*4          refmodel(0:ndmax)
	integer		nsub(0:ndmax)
	integer		ndec(0:ndmax)
	character*72	str(0:ndmax)
#else
c                                                      NAD file
        real*4, allocatable :: df(:)        ! 0->ndmax
        real*4, allocatable :: ticl(:)      ! 0->ndmax
        real*4, allocatable :: refmodel(:)  ! 0->ndmax
        integer, allocatable :: nsub(:)     ! 0->ndmax
        integer, allocatable :: ndec(:)     ! 0->ndmax
        character*72, allocatable :: str(:) ! 0->ndmax
c                                                      frame arrays

        real*4, allocatable :: models(:)      ! ndmax*nemax
        real*4, allocatable :: data(:)        ! nemax
        integer, allocatable :: rank(:)       ! nemax
        character, allocatable :: header(:)   ! nhmax
#endif
	character*72	tstr
	real*4		dxinc(nplmax)
	real*4		dyinc(nplmax)
	real*4		dxinct(nplmax)
	real*4		dyinct(nplmax)

	real*4		range(2,0:nplmax)
	real*4		xb(4)
	real*4		yb(4)
	real*4		arr(ngridmax)

	integer		jx(nplmax),jy(nplmax)

	logical		model
	logical		pbest
	logical		postscript
	logical		finetune
	logical		grid
	logical		inbounds

        character*256   fnme
	character*256	string

	itick = 0
	do i=1,ngridmax
	   arr(i) = 0.0
	end do
	write(*,*)
	write(*,*)' S-plot: '
	write(*,*)

	grid = .false.
	nargs=iargc()
	if(nargs.eq.1)then
	   call getarg(1,string)
	   if(string(1:2).eq.'-g'.or.string(1:2).eq.'-G')then
	     grid = .true.
	   end if
	end if
c                                               read in ensemble of models
c                                               from direct access nad file
        lu_in = 15

c						take an advanced look
c						into command file to get
c						nad filename

	open(lu_in,file='splot.cmd',status='old')
	read(lu_in,*)
	read(lu_in,*)
	read(lu_in,*)
        read(lu_in,fmt='(a72)')string
        call fnamelen(string,72,n)
        if(n.eq.0)then
           write(*,*)
           write(*,*)' Error no input filename found in splot.cmd'
           write(*,*)
           stop
        end if
        read(string(1:n),fmt='(a72)')fnme

	close(lu_in)

	write(*,100)fnme
 100    format(1x,' Reading direct access file',24x,': ',a40)
#if NA_F77_COMPATIBILITY
        iform = 0
        call read_nad
     &               (lu_in,fnme,nh,nhu,iform,
     &                nd,ne,nhmax,ndmax,nemax,
     &                header,data,models)
#else
        iform = -1
c                                               read nad header 
        call read_nad
     &               (lu_in,fnme,nh,nhu,iform,
     &                nd,ne,nhmax,ndmax,nemax,
     &                header,data,models)
c                                               allocate memory to dynamic arrays
        allocate(df(0:nd))
        allocate(ticl(0:nd))
        allocate(refmodel(0:nd))
        allocate(nsub(0:nd))
        allocate(ndec(0:nd))
        allocate(str(0:nd))
        allocate(models(nd*ne))
        allocate(data(ne))
        allocate(rank(ne))
        allocate(header(nh))
        nhmax = nh
        nemax = ne
        ndmax = nd
c                                               read nad
        call read_nad
     &               (lu_in,fnme,nh,nhu,iform,
     &                nd,ne,nhmax,ndmax,nemax,
     &                header,data,models)
#endif

c                                               read command file 
c                                               to overide defaults
c
	call readcmds
     &               (nplmax,ndmax,nd,np,jx,jy,range,
     &		      ipenmin,ipenmax,fmin,fmax,ltype,
     &		      ipenref,ipenbest,pbest,title,
     &                frametyp,df,ticl,
     &                nsub,ndec,str,xtit,ytit,rmult,
     &                smultps,smultx,tmultps,tmultx,cmult)

	write(*,*)' Number of models in ensemble ',
     &            '                     = ',ne
	write(*,*)' Number of plots             ',
     &            '                      = ',np
	write(*,*)' Number of variables           ',
     &            '                    = ',nd
c	write(*,*)' Number of special variables   ',
c    &            '                    = ',nsp

	if(ne.gt.11000.and..not.grid)then
	   write(*,*)
	   write(*,*)' Note that for large numbers of models'
	   write(*,*)' more compact postscript files',
     &               ' and'
	   write(*,*)' quicker X graphics are produced'
	   write(*,*)' if the grid option is used.'
	   write(*,*)
	   write(*,*)' syntax: '
	   write(*,*)'     splot-x -g (or splot-p -g)'
	   write(*,*)
	end if
c						read in reference model
c						if one exists
        model = .false.
	finetune = .false.
        lu3 = 15
        open(lu3,file='model.ref',status='old',err=5)
	read(lu3,*)
	read(lu3,*)
	read(lu3,*)
        do i=1,nd
           read(lu3,*)refmodel(i)
        end do
        model  = .true.
	if(ipenref.lt.0)model = .false.

 5      continue
	
c	
	dmin = data(1)
	imin = 1
	dmax = dmin
	do i=1,ne
	   if(data(i).lt.dmin)imin = i
	   dmin = min(dmin,data(i))
	   dmax = max(dmax,data(i))
	end do
	write(*,*)' Minimum data value             ',
     &            '                   =',dmin
	write(*,*)' Index of model with min data   ',
     &            '                   =',imin
	write(*,*)' Maximum data value             ',
     &            '                   =',dmax
	if(grid)write(*,*)' Plotting discretized image'
c						rank models by misfit
	if(ltype.eq.0)call indexx(ne,data,rank)
c
c						choose layout of plot.
	xsize = 26.8
	ysize = 18.8
	xsize = 24.12
	ysize = 16.92
	ysize = 15.
c
	call arrangefig(np,xsize,ysize,nx,ny,xo,yo,xlen,ylen,dx,dy)

	xo = xo + 1.5
	yo = yo + 1.0
	if(frametyp.eq.2)dy = dy*1.2
        if(np.le.6)del = 0.3
c
        open(30,file='splot.ps',status='unknown')
c
c                                       A4
        call xname ('S-plot',6)
        call hplots(1,0,30,0)
c
c ---------------------------- Start plotting -------------------------
c
c       Read in colour table (if required)
c
	lu_col = 20
        open(lu_col,file='pal.in',status='unknown')
        call ldcolr(lu_col)
        call typset(0.0,0.0)
c       call zpick(3,0,0)
        call zpick(0,0,0)
c						determine whether
c						we are using postscript
	n = -1
	postscript = .false.
	call plottype(n)
	if(n.eq.0) postscript = .true.
	   smult=smultx
	   if(postscript)smult=smultps
	   tmult=tmultx
	   if(postscript)tmult=tmultps
c						
           ipen1 = 1
           ithk = 2
           ithk = 0
	   radmin = 0.07
	   csize = cmult*xlen*0.07
	   size = smult*xlen*0.07
	   rad = rmult*xlen/100.
	   if(nx.eq.2)rad = rmult*xlen/110.
	   if(nx.eq.1)rad = rmult*xlen/130.
	   if(rad.lt.radmin.and..not.grid)rad = radmin
	   ticlen = xlen/10.
	   xinc = 0.
	   xinc2 = 0.
	   yinc = 0.
	   yinc2 = 0.
	   if(ny.gt.1)xinc =0.5*ticlen 
	   if(ny.ge.3)xinc =0.7*ticlen 
	   if(ny.ge.4)xinc =1.5*ticlen 
	   if(ny.ge.5)xinc =2.0*ticlen 
	   if(nx.eq.3)xinc =0.3*ticlen 
	   if(ny.gt.1)yinc =0.4*ticlen 
	   if(ny.gt.3)yinc =0.9*ticlen 
	   if(ny.gt.3)yinc2 =0.8*ticlen 
	   if(ny.gt.3)size =0.4*ticlen 
           size = size*smult
	   tsize = tmult*0.05*xsize
           call pen(1,ithk)

	   if(postscript)then
             if(.not.grid)rad = rad*15.
             size = size*1.3
	     if(ny.ge.4)size = size*2.
	     if(ny.ge.3)size = size*1.4
             tsize = tsize *0.8
             yinc = 0.
             yinc2 = 0.
             xinc = 0
	     ytit = ytit+1
	     yo = yo + 1.
	     if(ny.eq.1)then
               size=size*1.2
	       yinc = -0.3*ticlen
	       yinc2 = -0.2*ticlen
	       xinc = -0.2*ticlen
	       xinc2 = -0.5*ticlen
	     end if
c						read in special 
c						fine tuning increments
c						for postscript figure

	     open(31,file='finetune',status='old',err=10)
	     finetune = .true.
	     do i=1,23
	        read(31,*)
	     end do
	     do i=1,np
	        read(31,*,err=10)dxinc(i),dxinct(i),
     &                           dyinc(i),dyinct(i)
	     end do
	     do i=1,np
	     end do
	     write(*,*)' Finetuning of axes read in' 
 10          continue
	     close(31)
           end if
	   nup = 200
	   if(ny.gt.1)nup = 100
c	   ne = 0
c	   fmax = log(fmax)
c	   fmin = log(fmin)
	   imap = 0
	   imap = 1
	   eps = 0.0001
c                                               plot panels
           ii = 0
           do iy=1,ny
              ys = yo - (iy-1)*(ylen+dy)
              do ix=1,nx
                 xs = xo + (ix-1)*(xlen+dx)
                 xb(1) = xs
                 xb(2) = xs+xlen
                 xb(3) = xs+xlen
                 xb(4) = xs
                 yb(1) = ys
                 yb(2) = ys
                 yb(3) = ys+ylen
                 yb(4) = ys+ylen
                 ii = ii + 1
                 i = ii

		 kx = jx(i)
		 ky = jy(i)
	         if(kx.eq.0)then
                    range(1,0) = 1. 
                    range(2,0) = ne
	            ndec(0) = 0
                    df(0) = ne
                    refmodel(0) = ne/2 
                    nsub(0) = 10
                    str(0) = 'No. of models'
	            d2x = xlen/(range(2,0)-range(1,0))
	         else if(kx.lt.0.or.kx.gt.nd)then
                    write(*,*)
                    write(*,*)' Error in splot.cmd' 
                    write(*,*)
                    write(*,*)' Input variable out of range'
                    write(*,*)' Trying to plot variables       :',kx,ky
                    write(*,*)' Range of variables in NAD file :',1,nd
                    write(*,*)
                    stop
	         else
	            d2x = xlen/(range(2,kx)-range(1,kx))
                 end if
	         if(ky.le.0.or.ky.gt.nd)then
                    write(*,*)
                    write(*,*)' Error in splot.cmd' 
                    write(*,*)
                    write(*,*)' Input variable out of range'
                    write(*,*)' Trying to plot variables       :',kx,ky
                    write(*,*)' Range of variables in NAD file :',1,nd
                    write(*,*)
                    stop
                 end if
	         d2y = ylen/(range(2,ky)-range(1,ky))
	         d2p = (ipenmax-ipenmin+1)/(fmax-fmin)

	         if(grid)then
c						fill grid values
c						and plot matrix
	            nsx = int(xlen/rad)
	            nsx = min(nsx,nup)
	            nsy = nsx
	            dsx = xlen/real(nsx)
	            dsy = ylen/real(nsy)
	            if(nsx*nsy.gt.ngridmax)then
		    write(*,*)
		    write(*,*)' Error array sizes not large enough'
		    write(*,*)' nsx      = ',nsx
		    write(*,*)' nsy      = ',nsy
		    write(*,*)' ngridmax = ',ngridmax
		    write(*,*)
		    write(*,*)' Remedy'
		    write(*,*)' increase size of ngridmax',
     &                        ' to at least ',nsx*nsy
		    stop
	            end if
	            s11 = ipenmin
	            s22 = ipenmax
		    do j=1,ne
                       if(kx.eq.0)then
                          x = (j-range(1,kx))*d2x
		       else
	                  k = kx + nd*(j-1)
                          x = (models(k)-range(1,kx))*d2x
	               end if
                       k = ky + nd*(j-1)
                       y = (models(k)-range(1,ky))*d2y
                       dat = data(j)
		       if(ltype.eq.0)then
                          dat = rank(j)*100./real(ne)
		       end if
	               dinc = (dat-fmin)*d2p
c		       if(mod(dinc,1.).eq.0.0)dinc = dinc-0.5
	               if(dinc.lt.1.+eps.and.
     &                 dinc.gt.1.0-eps)dinc = dinc-0.5
	               ip = ipenmin + dinc
                       if(dat.gt.fmax)ip = ipenmax
                       if(dat.lt.fmin)ip = ipenmin
		       iix = 1+x/dsx
		       iiy = 1+y/dsy
	               iv = iix + (iiy-1)*nsx
	               if(ii.ge.1.and.ii.le.ngridmax)then
		          arr(iv) = ip
		       end if
		    end do
c						plot grid
c
                    call picol
     &                   (xs,ys,xlen,ylen,nsx,nsy,
     &                    arr,ipenmin,ipenmax,s11,s22,imap)
		    do js=1,nsy
		    do is=1,nsx
	               k = is + (js-1)*nsx
	               arr(k) = 0.0
		    end do
		    end do

		 else

c						plot points
		 iplast = 0
 	         do j = 1,ne
	            inbounds = .true.
                    if(kx.eq.0)then
	               k1 = kx + nd*(j-1)
                       x = xs + (j-range(1,kx))*d2x
		    else
	               k1 = kx + nd*(j-1)
                       x = xs + (models(k1)-range(1,kx))*d2x
	            end if
	            k2 = ky + nd*(j-1)
                    y = ys + (models(k2)-range(1,ky))*d2y
c	            dat = log(data(j))
 	            dat = data(j)
		    if(ltype.eq.0)then
 	               dat = rank(j)*100./real(ne)
		    end if
	            dinc = (dat-fmin)*d2p
c		    if(mod(dinc,1.).eq.0.0)dinc = dinc-0.5
	            if(dinc.lt.1.+eps.and.
     &                 dinc.gt.1.0-eps)dinc = dinc-0.5
	            ip = ipenmin + dinc
	            if(dat.gt.fmax)ip = ipenmax
	            if(dat.lt.fmin)ip = ipenmin
c	            write(*,*)j,dat,dinc,ip
		    if(ip.ne.iplast)call pen(ip,ithk) 
	            iplast = ip
c						check if model is
c						within plot range
c
                    if(kx.eq.0)then
		    if(models(k2).lt.range(1,ky).or.
     &		       models(k2).gt.range(2,ky))inbounds = .false.
	            else
		    if(models(k1).lt.range(1,kx).or.
     &		       models(k2).lt.range(1,ky).or.
     &		       models(k1).gt.range(2,kx).or.
     &		       models(k2).gt.range(2,ky))inbounds = .false.
	            end if
	 	    if(inbounds)then
	               if(postscript)then
 	                  call pcirclef(x,y,rad)
	               else 
 	                  call csymbl(x,y,3,rad,21)
		       end if
	            end if
		 end do

	  	 end if

c						plot frame and labels

		 if(frametyp.le.1)then

		 call pen(ipen1,2) 
                 call plot_box(xb,yb,1)
                 call pen(ipen1,0)

                 if(kx.lt.10)then
                    if(ky.lt.10)then
                       write(string,fmt='(i1," and ",i1)')kx,ky
                    else if(ky.lt.100)then
                       write(string,fmt='(i1," and ",i2)')kx,ky
                    else if(ky.lt.1000)then
                       write(string,fmt='(i1," and ",i3)')kx,ky
                    end if
                 else if(kx.lt.100)then
                    if(ky.lt.10)then
                       write(string,fmt='(i2," and ",i1)')kx,ky
                    else if(ky.lt.100)then
                       write(string,fmt='(i2," and ",i2)')kx,ky
                    else if(ky.lt.1000)then
                       write(string,fmt='(i2," and ",i3)')kx,ky
                    end if
                 else if(kx.lt.1000)then
                    if(ky.lt.10)then
                       write(string,fmt='(i3," and ",i1)')kx,ky
                    else if(ky.lt.100)then
                       write(string,fmt='(i3," and ",i2)')kx,ky
                    else if(ky.lt.1000)then
                       write(string,fmt='(i3," and ",i3)')kx,ky
                    end if
                 end if
c
                 nstr = 11
                 call pen(ipen1,ithk)
	         xn = xs + xlen*0.5-0.22*nstr*size
                 call typstr(xn,ys+1.05*ylen,size,string,0.0,nstr)

		 if(frametyp.eq.1)then
c						add range labels
		    
	            sl = size*0.6
                    nxd = ndec(kx)
	            off = nxd*sl*0.3
	            x = xs -off
	            y = ys - 2.0*sl
                    call typnum
     &              (X,Y,sl,range(1,kx),0.,nxd)
	            x = xs + xlen-off 
                    call typnum
     &              (X,Y,sl,range(2,kx),0.,nxd)
c
	 	    nyd = ndec(ky)
	            off = 3*sl + nyd*sl
	            x = xs - off
	            y = ys 
                    call typnum
     &              (X,Y,sl,range(1,ky),0.,nyd)
	            y = ys + ylen 
                    call typnum
     &              (X,Y,sl,range(2,ky),0.,nyd)
	           
		 end if

		 else if(frametyp.eq.2)then
c						plot full 
c                                               annotated frame
	           tstr = ' '
	           sizl = size*0.8
	           sizt = size*0.8
	           ipenf1 = 1
	           ipenf2 = 1
	           ipenf3 = 1
	           ticlx = d2x*df(kx)
	           ticly = d2y*df(ky)
	           f0x = range(1,kx)
	           f0y = range(1,ky)
                   xincs = xinc
                   yincs = yinc
                   xincs2 = xinc2
                   yincs2 = yinc2
	           if(finetune)then
                     xincs = xinc+dxinc(ii)
                     yincs = yinc+dyinc(ii)
                     xincs2 = xinc2+dxinct(ii)
                     yincs2 = yinc2+dyinct(ii)
	           end if
                   call framer(xs,ys,xlen,ylen,
     &                         f0x,ticlx,
     &                         df(kx),nsub(kx),ndec(kx),
     &                         f0y,ticly,
     &                         df(ky),nsub(ky),ndec(ky),
     &                         str(kx),str(ky),tstr,sizl,sizt,
     &                         ipenf1,ipenf2,ipenf3,
     &                         xincs,yincs,yincs2,xincs2,
     &                         ticlen,postscript)
		   call pen(ipen1,2) 
                   call plot_box(xb,yb,1)
                   call pen(ipen1,1)
		 end if


	         if(model)then
c                                               plot ref model

	           x = xs + (refmodel(kx)-range(1,kx))*d2x
		   y = ys + (refmodel(ky)-range(1,ky))*d2y
	           call pen(0,ithk)
	           xb(1) = x - csize*0.5
	           yb(1) = y - csize*0.5
	           xb(2) = x + csize*0.5
	           yb(2) = y - csize*0.5
	           xb(3) = x + csize*0.5
	           yb(3) = y + csize*0.5
	           xb(4) = x - csize*0.5
	           yb(4) = y + csize*0.5
                   call fillpoly(xb,yb,4)
	           call pen(ipenref,0)
c		   call csymbl(x,y,3,csize,4) 
	           call plot(x+csize*0.45,y+csize*0.45,3)
	           call plot(x-csize*0.45,y-csize*0.45,2)
	           call plot(x+csize*0.45,y-csize*0.45,3)
	           call plot(x-csize*0.45,y+csize*0.45,2)
		 end if
c
	         if(pbest)then
c                                               plot best model

	           k = kx + nd*(imin-1)
	           x = xs + (models(k)-range(1,kx))*d2x
	           k = ky + nd*(imin-1)
		   y = ys + (models(k)-range(1,ky))*d2y
	           call pen(0,0)
	           xb(1) = x - csize*0.5
	           yb(1) = y - csize*0.5
	           xb(2) = x + csize*0.5
	           yb(2) = y - csize*0.5
	           xb(3) = x + csize*0.5
	           yb(3) = y + csize*0.5
	           xb(4) = x - csize*0.5
	           yb(4) = y + csize*0.5
                   call fillpoly(xb,yb,4)
	           call pen(ipenbest,0)
c		   call csymbl(x,y,3,csize,2) 
	           call plot(x+csize*0.45,y,3)
	           call plot(x-csize*0.45,y,2)
	           call plot(x,y-csize*0.45,3)
	           call plot(x,y+csize*0.45,2)
		 end if
              end do
           end do

c						write title
           call pen(1,2)
	   yt = yo + (ylen) + 1.5
	   if(ny.eq.1)yt = yo + ylen + 2.0
	   nchar=ilenr(title)
	   xt = xo + 0.5*(nx*xlen+(nx-1)*dx) - 0.27*nchar*0.5*tsize
           call typstr(xtit,ytit,tsize,title,0.0,nchar)
           call pen(1,0)
c						draw legend
	   wleg = 5.0
	   hleg = 0.75
	   scleg = 0.6*size
	   if(ny.gt.1)scleg = 0.8*size
           xl = xo + 0.5*(nx*xlen+(nx-1)*dx) 
           yl = yo - (ny-1)*(ylen+dy) - 1.5*hleg-0.3
	   if(ny.eq.1)yl = yl - 1.

           if(frametyp.eq.2)yl = yl - hleg
           call draw_legend_fill
     &           (xl,yl,fmax,fmin,ipenmin,ipenmax,
     &            wleg,hleg,scleg,ltype)

	write(*,*)
c
c ---------------------------- Stop plotting -------------------------
c
c       Call hplots to terminate plotting
c
        call hplots (0,0,0,0)

	if(postscript)then
           write(*,*)' Created postscript file : splot.ps'
           write(*,*)
	end if

	stop
	end
c
c ---------------------------------------------------------------------
c
c       arrangefig - decides how to arrange a multiple plot 
c		     on a landscape page
c
c       Input:
c		np		: number
c		xsize		: x-length of plot in cm
c		ysize		: y-length of plot in cm
c
c       Output:
c		xo,yo		: origin of plot
c		nx,ny		: Number of min-plots in x and y directions
c		xlen,ylen	: dimensions of mini plots
c		dx,dy		: distance between mini plots.
c
c	Comments
c		All dimensions are relative to an A3 page, although
c		A postscript plot will automatically scale to A4. 
c
c					M. Sambridge, RSES, March. 1999
c
c ---------------------------------------------------------------------
c
	subroutine arrangefig
     &             (np,xsize,ysize,nx,ny,xo,yo,xlen,ylen,dx,dy)

	xo = 0.2
	yo = 1.8

	if(np.eq.1)then
	   nx = 1
	   ny = 1
	else if(np.eq.2)then
	   nx = 2
	   ny = 1
	else if(np.eq.3)then
	   nx = 3
	   ny = 1
	else if(np.eq.4)then
	   nx = 2
	   ny = 2
	else if(np.gt.4.and.np.le.6)then
	   nx = 3
	   ny = 2
	else if(np.gt.6.and.np.le.8)then
	   nx = 4
	   ny = 2
	else if(np.gt.8.and.np.le.12)then
	   nx = 4
	   ny = 3
	else if(np.gt.12.and.np.le.16)then
	   nx = 4
	   ny = 4
	else if(np.gt.16.and.np.le.20)then
	   nx = 5
	   ny = 4
	else if(np.gt.20.and.np.le.30)then
	   nx = 6
	   ny = 5
	else if(np.gt.30.and.np.le.36)then
	   nx = 6
	   ny = 6
	else if(np.gt.36.and.np.le.42)then
	   nx = 7
	   ny = 6
	else if(np.gt.42.and.np.le.48)then
	   nx = 8
	   ny = 6
	else
	   write(*,*)' Cannot handle more than 49 plots'
	end if

	if(np.eq.0)then
	   nx = 1
	   ny = 1
	   xlen = xsize
	   ylen = ysize
	   dx = 0.
	   dy = 0.
	else if(np.eq.-1)then
	   nx = 1
	   ny = 1
	   xlen = min(xsize,ysize)*.8
	   ylen = xlen
	   xo = 0.2 + (xsize-xlen)/2. 
	   yo = 0.2 + 1
	   dx = 0.
	   dy = 0.
	else
	   xv = xsize/nx
	   yv = ysize/ny
	   xlen = min(xv,yv)*.8
	   ylen = xlen
	   dx = xv - xlen
	   dy = yv - ylen
	   xo = xo + 0.5*dx
	   yo = yo + ysize - yv +0.5*dy
	   if(np.eq.4)then
             dx = dx*0.5
           end if
	end if

	if(ny.eq.2)dy = 0.8*dy

	return
	end
c
c ---------------------------------------------------------------------
c
c       plot_box - plots a rectangle
c
c       Input:
c             xb                : x-coordinates of rectangle
c             yb                : y-coordinates of rectangle
c
c       Output:
c
c ---------------------------------------------------------------------
c
	Subroutine plot_box(xb,yb,ipen)
	
	real*4	xb(4),yb(4)

	call pen(ipen,0)
	call plot(xb(1),yb(1),3)
	call plot(xb(2),yb(2),2)
	call plot(xb(3),yb(3),2)
	call plot(xb(4),yb(4),2)
	call plot(xb(1),yb(1),2)


	return
	end
c
c ---------------------------------------------------------------------
c
c       readcmds - reads in extra commands controlling plots
c
c       Input:
c             nplmax               : maximum number of plots 
c             ndmax                : maximum number of variables 
c
c       Output:
c             np                   : number of plots 
c             jx                   : x-variables for ith plot 
c             jy                   : y-variables for ith plot 
c	      ...		   : and many other things
c
c       Comments:
c
c                                       M. Sambridge, RSES, August 1998.
c
c ---------------------------------------------------------------------
c

	Subroutine readcmds
     &             (nplmax,ndmax,nd,np,jx,jy,range,
     &		    ipenmin,ipenmax,fmin,fmax,ltype,
     &		    ipenref,ipenbest,pbest,title,
     &              frametyp,df,ticl,nsub,ndec,str,
     &		    xtit,ytit,rmult,smultps,smultx,
     &              tmultps,tmultx,cmult)

        integer         jx(*)
        integer         jy(*)
	real*4		range(2,0:*)
	logical		pbest
	character*256	string	
	character*50	title	
c	character*1	yesno	
c							frame arrays
	real*4		df(0:*)
	real*4		ticl(0:*)
	integer		nsub(0:*)
	integer		ndec(0:*)
	character*72	str(0:*)


        open(50,file='splot.cmd',status='old')

  5     read(50,fmt='(a72)',end=6)string
        if(string(1:1).eq.'#')go to 5
        read(50,*)np
	if(np.gt.nplmax)then
	   write(*,*)
	   write(*,*)' Too many plots requested'
	   write(*,*)
	   write(*,*)' Number of plots requested : ',np
	   write(*,*)' Maximum number of plots   : ',nplmax
	   write(*,*)
	   stop
	end if
           
	do i=1,np
	   read(50,*)jx(i),jy(i)
	end do
	read(50,*)ipenref
	read(50,*)ipenbest
	pbest = .true.
	if(ipenbest.lt.0)pbest = .false.
        ltype = 0
	read(50,*)fmin,fmax,ltype
	if(ltype.ne.0)ltype = 1
	read(50,*)ipenmin,ipenmax
	read(50,fmt='(a44)')title
	rmult = 1.0
	smult = 1.0
	cmult = 1.0
	tmult = 1.0
	read(50,*)xtit,ytit,tmultps,tmultx
	read(50,*)smultps,smultx
	read(50,*)cmult
	read(50,*)rmult
	read(50,*)frametyp

  6     continue

        close(50)

	write(*,*)' Options read in from command file:'

c						read range data
c						and frame parameters
c
        open(51,file='ranges.in',status='old')
	read(51,*)
	read(51,*)
	read(51,*)
	read(51,*)
	do i=1,nd
	   read(51,*)range(1,i),range(2,i),ndec(i)
	end do
	if(frametyp.eq.2)then
	   read(51,*)
	   read(51,*)
	   read(51,*)
	   read(51,*)
	   read(51,*)
	   do i=1,nd
	      read(51,fmt='(a70)')string
	      read(string,*)df(i),nsub(i)
              n = ilenr(string)
	      k1 = 0
	      k2 = 0
	      do j=1,n
                 if(string(j:j).eq.':')then
                    k1 = k2
                    k2 = j
		 end if
              end do
	      if(k1.eq.0.or.k2.eq.0)then
		 write(*,*)' Error reading file ranges.in'
		 write(*,*)' Frame label character string',
     &                     ' of variable ',i
                 write(*,*)' is not between a pair of colons'
	         stop
	      end if
	      read(string(k1+1:k2-1),fmt='(a30)')str(i)
c 	      write(*,*)df(i),nsub(i),str(i)
           end do
	   close (51)
	end if

	return
	end
c
c       Function ilenr returns the position of the first non blank character
c       in string `char' searching from the right.
c
        Function ilenr(char)
        character*(*) char
        ilenr = 0
c                               determine length of character string
        n = len(char)
c                               find first non blank character from the right.
        do 10 i = n,1,-1
           if(char(i:i).ne.' ')go to 11
 10     continue
        write(6,*)'Warning subroutine ilenr no non',
     &              ' blank characters found'
        ilenr = 0
        return
 11     ilenr = i
        return
        end

c
c-----------------------------------------------------------------------
c
c       Numerical recipes routine
c
c-----------------------------------------------------------------------
c
      SUBROUTINE indexx(n,arr,indx)
      INTEGER n,indx(n),M,NSTACK
      REAL arr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
      REAL a
      do 11 j=1,n
        indx(j)=j
11    continue
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 13 j=l+1,ir
          indxt=indx(j)
          a=arr(indxt)
          do 12 i=j-1,1,-1
            if(arr(indx(i)).le.a)goto 2
            indx(i+1)=indx(i)
12        continue
          i=0
2         indx(i+1)=indxt
13      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        itemp=indx(k)
        indx(k)=indx(l+1)
        indx(l+1)=itemp
        if(arr(indx(l+1)).gt.arr(indx(ir)))then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l)).gt.arr(indx(ir)))then
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l+1)).gt.arr(indx(l)))then
          itemp=indx(l+1)
          indx(l+1)=indx(l)
          indx(l)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx(l)
        a=arr(indxt)
3       continue
          i=i+1
        if(arr(indx(i)).lt.a)goto 3
4       continue
          j=j-1
        if(arr(indx(j)).gt.a)goto 4
        if(j.lt.i)goto 5
        itemp=indx(i)
        indx(i)=indx(j)
        indx(j)=itemp
        goto 3
5       indx(l)=indx(j)
        indx(j)=indxt
        jstack=jstack+2
        if(jstack.gt.NSTACK)pause 'NSTACK too small in indexx'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END
c
c --------------------------------------------------------------------
c
c     Subroutine draw_legend_fill - draws a filled legend 
c
c --------------------------------------------------------------------
c
      Subroutine draw_legend_fill
     &           (xl,yl,fmax,fmin,ipen1,ipen2,
     &            wleg,hleg,scleg,ltype)

      character*10      xlab1,xlab2
      real*4    xa(4),ya(4)

      npens = ipen2-ipen1+1
      call pen(1,0)
      xo = xl - 0.5*wleg
      yo = yl
      dec = 0.3*hleg
c
      dx = wleg/npens
      df = (fmax-fmin)/npens
c                                       fill legend
      do i=1,npens
         ipen = ipen1 + (i-1)
         x = xo + (i-1)*dx
         y = yo
         xa(1) = x
         ya(1) = y
         xa(2) = x+dx
         ya(2) = y
         xa(3) = x+dx
         ya(3) = y+hleg
         xa(4) = x
         ya(4) = y+hleg
         call pen(ipen,0)
         call fillpoly(xa,ya,4)
      end do
c                                       draw legend frame
      call pen(1,0)
      call plot(xo,yo-dec,3)
      call plot(xo,yo+hleg,2)
      call plot(xo+wleg,yo+hleg,2)
      call plot(xo+wleg,yo-dec,2)
      call plot(xo+wleg,yo,2)
      call plot(xo,yo,2)
c                                       write labels
 
      write(xlab1,fmt='(f3.0,"%")')fmin
      write(xlab2,fmt='(f3.0,"%")')fmax
      k = ilenr(xlab1)
      xt = xo - k*0.1*scleg
      yt = yo-dec-1.1*scleg
c     call typstr(xt,yt,scleg,xlab1,0.0,k)
      idec = 1
      if(fmin.ne.0.0)then
         if(abs(fmin).lt.0.1)idec = 2
         if(abs(fmin).lt.0.01)idec = 3
         if(abs(fmin).lt.0.001)idec = 4
         if(abs(fmin).lt.0.0001)idec = 5
      end if
      if(ltype.eq.0)idec = -1
      
      call typnum(xt,yt,scleg,fmin,0.0,idec)

      xt = xo + wleg - (k-1)*0.2*scleg
c     call typstr(xt,yt,scleg,xlab2,0.0,k)
      idec = 1
      if(fmax.ne.0.0)then
         if(abs(fmax).lt.0.1)idec = 2
         if(abs(fmax).lt.0.01)idec = 3
         if(abs(fmax).lt.0.001)idec = 4
         if(abs(fmax).lt.0.0001)idec = 5
      end if
      if(ltype.eq.0)idec = -1
      call typnum(xt,yt,scleg,fmax,0.0,idec)

      if(ltype.eq.0)then
        k = 10
        xt = xo + 0.5*wleg - k*0.5*0.5*scleg
        call typstr(xt,yt,scleg,'Percentile',0.0,k)
      else
        k = 6
        xt = xo + 0.5*wleg - k*0.5*0.5*scleg
        call typstr(xt,yt,scleg,'Misfit',0.0,k)
      end if

      return
      end

C===================================================================
C
      SUBROUTINE framer
     &           (X0,Y0,XASL,YASL,FX0,XTICL,DFX,NXSUB,NXDEC,FY0,
     &            YTICL,DFY,NYSUB,NYDEC,xstr,ystr,tstr,sizl,sizt,
     &            ipen1,ipen2,ipen3,xinc,yinc,
     &            yinc2,xinc2,ticlen,postscript)
C
c
C     This subroutine draws a frame
C     Parameters:
C     X0,Y0 = coordinates of lower left corner (cm)
C     XASL,YASL length of x- and y-axis (cm)
C     FX0,FY0 = functional values at lower left corner
C     XTICL,YTICL = distance (cm) between large tics
C     DFX,DFY = step in function values between two large tics
C     NXSUB,NYSUB = number of subdivisions with small tics between any two
C                   large tics
C     NXDEC,NYDEC = number of decimals in large tic labels (set -1 for integer)
C     ITX(18),ITY(18)= text for x and y axis respectively (centered at char 18
C     ITT(18)= title text (max 72 char, start at ITT(1), unlike ITX,ITY)
C     sizl = labelsize in cm
C     sizt = textletter size in cm
C     ipen1 frame-box pen     
C     ipen2 tics and symbols pen                       
C     ipen3 text pen
c
c
C
      character*72 xstr,ystr,tstr
      logical		postscript
c     DATA CM/2.54/
C
C
C  Formerly convert to inches
      st=sizt
      sl=sizl
      X0I=X0
      Y0I=Y0
      xlow=x0i
      ylow=y0i
      XASLI=XASL
      YASLI=YASL
      xup=x0i+xasli
      yup=y0i+yasli
C
      call Pen(ipen1,0)
C     Move to lower left corner and draw frame
      X=X0I
      Y=Y0I
      call Plot(X,Y,3)
      call Plot(x,yup,2)
      call Plot(xup,yup,2)
      call Plot(xup,y,2)
      call Plot(x,y,2)
c
C
C     Compute subdivision parameters
      NX=XASL/XTICL+1.01
      NY=YASL/YTICL+1.01
      NXSUB=MAX0(1,NXSUB)
      NYSUB=MAX0(1,NYSUB)
      DX=XTICL/real(NXSUB)
      DY=YTICL/real(NYSUB)
c
c     Reorigin in function values relative to (fx,fy)=(0,0)
      if(fx0.lt.0.and.dfx.gt.0) rox=(ifix(fx0/dfx)-1)*dfx
      if(fx0.ge.0.and.dfx.gt.0) rox=(ifix(fx0/dfx)  )*dfx
      if(fx0.lt.0.and.dfx.lt.0) rox=(ifix(fx0/dfx)  )*dfx
      if(fx0.ge.0.and.dfx.lt.0) rox=(ifix(fx0/dfx)-1)*dfx
*                                            
      if(fy0.lt.0.and.dfy.gt.0) roy=(ifix(fy0/dfy)-1)*dfy
      if(fy0.ge.0.and.dfy.gt.0) roy=(ifix(fy0/dfy)  )*dfy
      if(fy0.lt.0.and.dfy.lt.0) roy=(ifix(fy0/dfy)  )*dfy
      if(fy0.ge.0.and.dfy.lt.0) roy=(ifix(fy0/dfy)-1)*dfy
*                                            
      x0i=x0i-xticl/dfx *(fx0-rox)
      y0i=y0i-yticl/dfy *(fy0-roy)
c                  
      x=x0i
      y=ylow
      F=rox
      call Plot(xlow,ylow,3)
C
*
      call Pen(ipen2,0)
*
C     Plot tics and numbers lower x-axis
      DO 20 I=1,NX+1
      DO 10 J=1,NXSUB
      X=X+DX
      if(x.lt.xlow) goto 10
      if(x.gt.xup)  goto 21
      call Plot(X,Y,3)
      call Plot(X,Y+.25*ticlen,2)
      call Plot(X,Y,3)
   10 CONTINUE
      F=F+DFX
      call Plot(x,y,3)
      call Plot(X,Y+.40*ticlen,2)
      nxd=nxdec
      call nfig(f,nxd,nf)
c							x-axis
c     write(*,*)' x-axis nf = ',nf,' nxd = ',nxd
      call typnum(X-sl*0.2*nf,Y-1.1*sl-yinc,0.8*sl,F,0.,NXD)
      call Plot(X,Y,3)
   20 CONTINUE
c  21 call numc(xstr,nct)
   21 continue
      call Pen(ipen3,0)
      k = ilenr(xstr)
c     xback = 18*0.5*sl
      xback = 0.5*k*0.7*sl
      if(postscript)xback = 0.5*k*0.5*sl
      call typstr(Xlow+.5*XASLI-xback,
     &            Ylow-2.5*sl-yinc-yinc2,sl,xstr,0.,72)
      call Pen(ipen2,0)
C
      X=Xlow
      Y=Y0I
      call Plot(Xlow,ylow,3)
      F=roy
C
C     Plot left y-axis
      DO 40 I=1,NY+1
      DO 30 J=1,NYSUB
      Y=Y+DY
      if(y.lt.ylow) goto 30
      if(y.gt.yup)  goto 41
      call Plot(X,Y,3)
      call Plot(X+.25*ticlen,Y,2)
      call Plot(X,Y,3)
   30 CONTINUE
      F=F+DFY
      call Plot(x,y,3)
      call Plot(X+.40*ticlen,Y,2)
      nyd=nydec
c							y-axis
      call nfig(f,nyd,nf)
c     write(*,*)' y-axis nf = ',nf,' nyd = ',nyd
      call typnum(X-sl*0.55*nf-xinc,Y-sl*0.23,0.8*sl,F,0.,NYD)
c     call typnum(X-.35,Y-sl*nf/3.,0.8*sl,F,90.,NYD)
      call Plot(X,Y,3)
   40 CONTINUE
c  41 call numc(ystr,nct)
   41 continue
      call Pen(ipen3,0)
c     call typstr(Xlow-sl,Ylow+.5*YASLI-2.40,sl,ystr,90.,72)
      k = ilenr(ystr)
c     zback = 18*0.5*sl
      zback = 0.5*k*0.5*sl
      call typstr(Xlow-2.5*sl-xinc-xinc2,
     &            Ylow+.5*YASLI-zback,sl,ystr,90.,72)
      call Pen(ipen2,0)
C
C     Move to upper left corner
      X=X0I
      Y=Yup
      call Plot(Xlow,ylow,3)
C                                                 
C     Plot top x-axis
      DO 60 I=1,NX+1
      DO 50 J=1,NXSUB
      X=X+DX
      if(x.lt.xlow) goto 50
      if(x.gt.xup)  goto 61
      call Plot(X,Y,3)
      call Plot(X,Y-.25*ticlen,2)
      call Plot(X,Y,3)
   50 CONTINUE
      call Plot(X,Y-.40*ticlen,2)
      call Plot(X,Y,3)
   60 CONTINUE
   61 continue
C
C
C     Plot right y-axis
      X=Xup
      Y=Y0I
      call Plot(xlow,ylow,3)
      DO 80 I=1,NY+1
      DO 70 J=1,NYSUB
      Y=Y+DY
      if(y.lt.ylow) goto 70
      if(y.gt.yup)  goto 81
      call Plot(X,Y,3)
      call Plot(X-.25*ticlen,Y,2)
      call Plot(X,Y,3)
   70 CONTINUE
      call Plot(X-.40*ticlen,Y,2)
      call Plot(X,Y,3)
   80 CONTINUE
   81 continue
C
C     title text
      call Pen(ipen3,0)
      call Plot(xup,yup,3)
      X=Xlow
c     call numc(tstr,nct)
      call Plot(X,Yup+1.50,3)
      call typstr(X,Yup+1.50,st,tstr,0.,72)
      call Pen(1,0)
      RETURN
      END
      subroutine nfig(f,ndec,nf)
      if(abs(f).lt.0.00001) then
            x=0
      else
            x=alog10(abs(f))
            if(abs(f-1).lt.0.0001) x=0
      endif
      if(x.ge.0) then
            nf=x+2+ndec
      else
            ng=abs(x)+2
            nf=max0(ng,ndec+2)
      endif
      if(f.lt.0) nf=nf+1
c  plotpak use
c     ndec= nf*10+ndec
      return
      end
      subroutine numc(c,nc)
      character*(*) c
      m=len(c)
      nc=m+1
10      nc=nc-1
      if(c(nc:nc).eq.' '.and.nc.gt.1) goto 10
      return
      end
c
c ---------------------------------------------------------------------
c
c       read_nad - read a direct access file in NAD format. 
c
c	Input:
c	      lu	      : logical unit of file
c	      fnme	      : filename 
c	      iform	      : If iform=-1 on input then just read header and return 
c	      nhmax	      : maximum size of array header
c	      ndmax	      : maximum size of array data
c	      nemax	      : maximum size of array models
c
c	Output:
c	      nh	      : total length in bytes of file header 
c	      nhu	      : length in bytes of user written part of header 
c	      iform	      : format for nad file (1=multi-record,0=single)
c	      mul	      : number of records in the header (with padding)
c	      nd	      : dimension of parameter space
c	      ne	      : number of models in ensemble
c	      header          : header character string of length nh (char)
c	      data(nd)	      : array of data values for each model (real*4)
c	      models(nd,ne)   : array of model values  (real*4)
c
c	Comments:
c                The direct access NAD file format:
c
c                VARIABLE       TYPE            SIZE IN BYTES
c                nd             int             4
c                ne             int             4
c                nh             int             4
c                nhu            int             4
c                header         character       nh
c                models         real*4          4*nd*ne
c                data           real*4          4*nd
c                tail           real*4          4
c
c                File must contain a single record of length
c                [4x(4+nd*ne+ne+1) + nh] bytes
c
c		Calls are made to subroutine read_da. 
c
c                This routine assumes that direct access
c                files are opened with the record length specified
c                in bytes. This is the default for most machines
c                but not Dec machines. (A compiler option is
c                available on the Compaq Alpha to use bytes rather than 
c                4-byte words.)
c
c					M. Sambridge, RSES, Nov. 2001
c
c ---------------------------------------------------------------------
c
      	Subroutine read_nad
     &                     (lu,fnme,nh,nhu,iform,nd,ne,
     &                      nhmax,ndmax,nemax,
     &                      header,data,models)
c
      	real*4            models(*)
      	real*4            data(*)
      	character*(*)     header
        logical           headeronly

      	character*256 	  fnme

	mul = 0
        headeronly = .false.
        if(iform.eq.-1)then
          headeronly = .true.
          iform = 0
        end if
	iform = 0
c						read in size of header
        len = 4
        open(lu,file=fnme,status='unknown',
     &       form='unformatted',access='direct',recl=len,err=100)
        read(lu,rec=1)mul
        close(lu)
c						read in number of models,
c						number of dimensions and 
c						size of header
	if(mul.gt.0)then
c						we are in single record format
          len = 16
          open(lu,file=fnme,status='unknown',
     &         form='unformatted',access='direct',recl=len)
          read(lu,rec=1)nd,ne,nh,nhu
          mul = 0
          close(lu)
	else
c						we are in multi record format
          len = 20
          open(lu,file=fnme,status='unknown',
     &         form='unformatted',access='direct',recl=len)
          read(lu,rec=1)mul,nd,ne,nh,nhu
          mul = -mul
          close(lu)
	  iform = 1
	end if
        if(headeronly)return
c						check array sizes
	if(nd.gt.ndmax)then
          write(*,*)
          write(*,*)
     &    'Error - maxmimum number of dimensions is too small'
          write(*,*)
     &    '        current maximum number of dimensions =',ndmax
          write(*,*)
     &    '        number of dimensions in input file   =',nd
          write(*,*)'        '
          write(*,*)'Remedy - increase size of parameter ndmax'
          write(*,*)'         to at least this value and recompile'
          write(*,*)'        '
          stop
        end if

	if(ne.gt.nemax)then
          write(*,*)
          write(*,*)
     &    'Error - maxmimum number of models is too small'
          write(*,*)
     &    '        current maximum number of models =',nemax
          write(*,*)
     &    '        number of models in input file   =',ne
          write(*,*)'        '
          write(*,*)'Remedy - increase size of parameter nemax'
          write(*,*)'         to at least this value and recompile'
          write(*,*)'        '
          stop
        end if

	if(nh.gt.nhmax)then
          write(*,*)
          write(*,*)
     &    'Error - maxmimum size of header is too small'
          write(*,*)
     &    '        current maximum size of header =',nhmax
          write(*,*)
     &    '        size of header in input file   =',nh
          write(*,*)'        '
          write(*,*)'Remedy - increase size of parameter nhmax'
          write(*,*)'         to at least this value and recompile'
          write(*,*)'        '
          stop
        end if

 	call read_da(lu,fnme,nd,ne,nh,mul,header,models,data)

      	return
 100    write(*,*)
        write(*,*)' Error reading file : ',fnme
        write(*,*)' Does it exist ? '
        write(*,*)
	stop

      	end

c
c ---------------------------------------------------------------------
c
c       read_da - read a direct access file containing
c                 ne models with dimension nd, ne data 
c                 values and a header of size nh.
c
c	Input:
c	      lu		: logical unit of file
c	      fnme		: filename 
c	      nh		: length in bytes of file header (minus padding)
c	      mul		: number of records in the header (with padding)
c	      nd		: dimension of parameter space
c	      ne		: number of models in ensemble
c
c	Output:
c	      header		: header character string of length nh (char)
c	      data(nd)		: array of data values for each model (real*4)
c	      models(nd,ne)	: array of model values  (real*4)
c
c	Comments:
c		 The direct access NAD file format:
c
c		 VARIABLE	TYPE		SIZE IN BYTES
c		 nd		int     	4
c		 ne		int     	4
c		 nh		int     	4
c		 header		character     	nh
c		 models		real*4     	4*nd*ne
c		 data		real*4     	4*nd
c		 tail		real*4     	4
c	
c	         File must contain a single record of length 
c		 [4x(3+nd*ne+ne+1) + nh] bytes
c
c
c					M. Sambridge, RSES, Nov. 2001
c
c ---------------------------------------------------------------------
c
      	Subroutine read_da(lu,fnme,nd,ne,nh,mul,header,models,data)
c
      	real*4            models(nd,ne)
      	real*4            data(ne)
      	character         header(nh)
c     	real*4            tail

      	character*256 	  fnme
c					
c                                              calculate length of header
 	iform = 0
 	if(mul.ne.0)iform = 1

 	lenh  = 4*5+nh
        len   = 4*(nd+1)

 	if(iform.eq.1)then

           open(lu,file=fnme,status='unknown',
     &          form='unformatted',access='direct',recl=lenh)

       	   read(lu,rec=1)idum1,idum2,idum3,idum4,idum5,header
	
 	   close(lu)

c	   write(*,*)' record length               ',len
c	   write(*,*)' header length               ',lenh
c	   write(*,*)' Number of records in header ',mul
c
           open(lu,file=fnme,status='unknown',
     &          form='unformatted',access='direct',recl=len)

           do i=1,ne
              call rnad(lu,i+mul,nd,models(1,i),data(i))
           end do

           close(lu)

 	else
c					       read in file as 
c					       a single record
	   len = 4+nd*ne+ne
           len = 4*len+nh

           open(lu,file=fnme,status='unknown',
     &          form='unformatted',access='direct',recl=len)

           read(lu,rec=1)i,j,k,kk,header,models,data

           close(lu)

 	end if
c
      	return
      	end

      	Subroutine rnad(lu,i,nd,models,data)

      	real*4            models(nd)
      	real*4            data

      	read(lu,rec=i)models,data

      	return
      	end

        Subroutine fnamelen (string,m,n)

        character*1     string(m)

        do i=1,m
          if(string(i).eq.' ')then
             n = i-1
             return
          end if
        end do

        return
        end

