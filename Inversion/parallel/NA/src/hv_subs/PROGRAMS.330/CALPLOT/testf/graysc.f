c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: GRAYSC                                                c
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
	real*4 x(41),y(41),amp(41)
	call pinitf('GRAYSC.PLT')
	call gunit('in')
c-----
c	draw some colored lines
c-----
	do 100 i=0,20
		call newpen(i)
		xx = 8.5
		yy = 7.0 - 0.3*i
		call gwidth(0.10)
		call plot(xx+0.5,yy,3)
		call plot(xx+1.00,yy,2)
		call newpen(1)
		call gwidth(0.01)
		call number(xx+1.1,yy,0.10,real(i),0.0,-1)
  100	continue
	call gwidth(0.01)
c-----
c	put out a nice gray scale according to ipen
c-----
	call graysc(1.5,1.0,2.0,7.0,1000,1020,1)
	call graysc(3.0,1.0,3.5,7.0,1020,1040,1)
	call graysc(4.5,1.0,5.0,7.0,1040,1060,1)
	call graysc(6.0,1.0,6.5,7.0,1060,1080,1)
	call graysc(7.5,1.0,8.0,7.0,1080,1100,1)
	call pend()
	stop
	end

	subroutine box(x0,y0,x1,y1)
		call plot(x0,y0,3)
		call plot(x0,y1,2)
		call plot(x1,y1,2)
		call plot(x1,y0,2)
		call plot(x0,y0,2)
	return
	end

	subroutine graysc(x0,y0,x1,y1,ipen0,ipen1,ipinc)
	call box(x0,y0,x1,y1)
	ipnmn=ipen0
	ipnmx=ipen1
	ipnd=ipnmx-ipnmn
	ylen = y1 - y0
	dy =  ylen/(ipnd+1)
	do 100 ipen=ipnmn,ipnmx,ipinc
		xl=x0
		xh=x1
		yh = y1 - (ipen-ipnmn)*dy
		yl = yh - dy*ipinc
		if(yl.lt.y0)yl=y0
		ipatx=0
		ipaty=0
		xlen=0.01
		ylen=0.01
		call newpen(ipen)
		call shader(xl,yl,xh,yh,ipatx,ipaty,xlen,ylen)
		call newpen(1)
		call number(xh+0.1,yh-dy/2.,0.10,real(ipen),0.0,-1)
  100	continue
	return
	end
