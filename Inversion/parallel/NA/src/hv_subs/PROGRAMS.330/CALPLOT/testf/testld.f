c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: TESTLD                                                c
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
c----
c	program to test the scale and lined routines
c----
	real xaxlen,yaxlen
c	parameter (xaxlen=5.0,yaxlen=1.0)
	dimension x(120),c(120),s(120)
	dimension cc(120),ss(120)
	xaxlen=5.0
	yaxlen=1.0
	xx = xaxlen
	yy = yaxlen
	do 100 i=1,100
		ang=0.01*6.2831853*(i-1)
		x(i)=ang
		c(i)=1234.0*cos(ang)
		cc(i)=0.015 * cos(ang)
		ss(i)=0.015 * cos(ang)
		s(i)=1234.0*sin(ang)
  100	continue
	call pinitf('TESTLD.PLT')
	call gunit('in')
	call plot(1.5,0.0,-3)
	do 200 j=1,3
		if(j.eq.1)then
			call plot(1.0,5.0,-3)
			call gscale(x,xaxlen,100,1)
			call gscale(c,yaxlen,100,1)
			call gscale(s,yaxlen,100,1)
			xlen = 0.05
			ipat = 5
			call axis(0.0,0.0,'X-AXIS',-6,xaxlen,0.0,
     1				x(101),x(102))
			call axis(0.0,0.0,'Y-AXIS',6,yaxlen,90.0,
     1				s(101),s(102))
			call lined(x,c,100,1,0,0,ipat,xlen)
			call lined(x,s,100,1,-10,1,ipat,xlen)
		elseif(j.eq.2)then
			call plot(0.0,-2.0,-3)
			call gscale(x,xaxlen,50,-2)
			call gscale(c,yaxlen,50,2)
			call gscale(s,yaxlen,50,2)
			xlen = 0.05
			ipat = 21
			call axis(0.0,0.0,'X-AXIS',-6,xaxlen,0.0,
     1				x(101),x(103))
			call axis(0.0,0.0,'Y-AXIS',6,yaxlen,90.0,
     1				s(101),s(103))
			call lined(x,c,50,2,-2,2,ipat,xlen)
			call lined(x,s,50,2,0,0,ipat,xlen)
		elseif(j.eq.3)then
			call plot(0.0,-2.0,-3)
			call gscale(x,xaxlen,100,1)
			call gscale(cc,yaxlen,100,1)
			call gscale(ss,yaxlen,100,1)
			xlen = 0.05
			ipat = 7
			call axis(0.0,0.0,'X-AXIS',-6,xaxlen,0.0,
     1				x(101),x(102))
			call axis(0.0,0.0,'Y-AXIS',6,yaxlen,90.0,
     1				ss(101),ss(102))
			call lined(x,cc,100,1,0,0,ipat,xlen)
			call lined(x,ss,100,1,-10,1,ipat,xlen)
		endif
  200	continue
	call pend()
	stop
	end
