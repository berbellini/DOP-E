c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: DTEST                                                 c
c                                                                     c
c      COPYRIGHT (C) 1987 R. B. Herrmann                              c
c                                                                     c
c      Department of Earth and Atmospheric Sciences                   c
c      Saint Louis University                                         c
c      221 North Grand Boulevard                                      c
c      St. Louis, Missouri 63103                                      c
c      U. S. A.                                                       c
c                                                                     c
c---------------------------------------------------------------------c
c-----
c	program to test dashed line routines
c-----
	call pinitf('DTEST.PLT')
	call gunit('in')
	call plot(0.0,0.0,-3)
c-----
c	test line algorithms
c-----
	do 2030 ipage = 1,3
		if(ipage.eq.1)then
			xlen = 0.05
		elseif(ipage.eq.2)then
			xlen = 0.1
		else
			xlen = 0.2
		endif
		call symbol(1.4,8.6,0.14,'TEST OF PLOTD: XLEN= ',0.0,21)
		call number(999.,999.,0.14,xlen,0.0,2)
		do 2010 ipat=0,31
			yy = 8.5 - 0.2*ipat
			xi = ipat
			call number(1.0,yy-0.05,0.10,xi,0.0,-1)
			call plot(1.4,yy,3)
			call plotd(7.5,yy,ipat,xlen)
 2010		continue
c-----
c		test ability to get correct pattern for a line
c		composed of many small segments
c-----
		do 2011 i=0,61
			yy = 2.0
			xx = 1.4 +  i*0.1
			if(i.eq.0)then
				call plot(xx,yy,3)
			else
				call plotd(xx,yy,31,xlen)
			endif
 2011		continue
		call axis(1.4,1.6,'LENGTH',-6,6.0,0.0,0.0,1.0)
		if(ipage.ne.3)call frame()
 2030	continue
c-----
c	draw a sequence of polygons to test algorithm
c	ability to correctly dash non-vertical and non-horizontal
c	line segments
c-----
	call frame()
	do 1000 k=0,31,3
		rad =  (k+1)*0.1
			x0 = 4.0
			y0 = 5.0
			xlen = 0.05
			inc = 15
			ipat = k
		do 100 i=0,360,inc
			ang = i*3.1415927/180.0
			xx = x0 + rad*cos(ang)
			yy = y0 + rad*sin(ang)
			if(i.eq.0)call plot(xx,yy,3)
			call plotd(xx,yy,ipat,xlen)
  100		continue
 1000	continue
	call pend()
	stop
	end
