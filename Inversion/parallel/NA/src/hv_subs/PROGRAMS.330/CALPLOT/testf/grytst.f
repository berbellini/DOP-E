c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: GRYTST                                                c
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
c-----
c	this program tests an improvement to the
c	CALPLOT device drivers to permit gray scale plots under
c	certain circumstances. This will be useful for plotting
c	seismic traces together with a trace attribute modulating
c	either the shading or the color of the background
c
c	To do this we invoke solid rectangle shading with ipatx = ipaty = 0
c	but with the color in newpen(ipen) > 1000. The present gray
c	shading will fill solid with ipen=1000, and with an
c	areal density decreasing with increasing ipen
c-----
	call pinitf('GRYTST.PLT')
	call gunit('in')
c-----
c	program to test gray shading
c-----
	xwid = 6.0
	ywid = 6.0
	nx = 40
	ny = 40
	dx = xwid/nx
	dy = ywid/ny
	ampmin = -1.0
	ampmax = 1.0
	ishdrv = 1
c-----
c	jplot = 1 multiply sine functions
c	jplot = 2 multiply triangular functions
c-----
	do 200 jplot=1,1
		call plot(1.0,1.0,-3)
		do 100 ix=1,nx
			x = (ix-1)*dx
			if(jplot.eq.1)then
				xfac = sin(3.1415927*x/xwid)
			else
				xfac = tri(x,xwid)
			endif
			do 101 iy=1,ny
				y = (iy-1)*dy
				if(jplot.eq.1)then
					yfac=sin(6.2831853*y/ywid)
				else
					yfac = tri(y+y,ywid)
				endif
				amp = xfac*yfac
				call putshd(x,x+dx,y,y+dy,amp,ishdrv,ampmin,ampmax)
  101			continue
  100		continue
		call box(0.0,0.0,xwid,ywid)
c-----
c		put out a nice gray scale according to ipen
c-----
		call graysc(ampmin,ampmax,ishdrv)
		call plot(-1.0,-1.0,-3)
		if(jplot.eq.1)call frame()
  200	continue
	call pend()
	stop
	end

	function tri(x,xl)
c-----
c	function to produce triangular function with a period
c	of 2*xl
c-----
c-----
c	first reduce periodicity to variable xval between [0,1]
c-----
		xval = amod(x,xl+xl)/(xl+xl)
		if(xval.le.0.25)then
			tri = xval
		elseif(xval.gt.0.25 .and. xval.le.0.75)then
			tri = 0.50 - xval
		elseif(xval.gt.0.75 .and. xval.le.1.00)then
			tri = xval - 1.0
		endif
		tri = 4.0*tri
	return
	end

	
	subroutine box(xl,yl,xh,yh)
		call plot(xl,yl,3)
		call plot(xl,yh,2)
		call plot(xh,yh,2)
		call plot(xh,yl,2)
		call plot(xl,yl,2)
		call plot(xl,yl,3)
	return
	end

	subroutine putshd(xl,xh,yl,yh,amp,ishdrv,ampmin,ampmax)
		xcol = (amp - ampmin)/(ampmax-ampmin)
		xcol = 100.0*xcol  + 0.5
		ipen = 1000 + xcol
		if(ishdrv .lt.0)ipen = 1100 - xcol
		if(ipen .lt. 1000)ipen = 1000
		if(ipen .gt. 1100)ipen=1100
		ipatx = 0
		ipaty = 0
		xlen = 0.02
		ylen = 0.02
c-----
c	ipen = 1000 is red, 1100 = blue or 1000 = dark, 1100 = light halftone
c-----
		call newpen(ipen)
		call shader(xl,yl,xh,yh,ipatx,ipaty,xlen,ylen)
	return
	end

	subroutine graysc(ampmin,ampmax,ishdrv)
		nbox=21
		damp = (ampmax-ampmin)/(nbox-1)
		xl=7.0
		xh=7.5
		dy = 6.0/nbox
		do 100 i=1,nbox
			amp = ampmin + (i-1)*damp
			yl=(i-1)*dy
			yh=yl+dy
			call putshd(xl,xh,yl,yh,amp,ishdrv,ampmin,ampmax)
			call newpen(1)
			call number(xh+0.1,yh-dy/2.,0.07,amp,0.0,2)
  100		continue
		call box(7.0,0.0,7.5,6.0)
	return
	end
