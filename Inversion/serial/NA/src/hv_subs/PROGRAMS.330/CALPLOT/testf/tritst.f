c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: TRITST                                                c
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
	call pinitf('TRITST.PLT')
	call gunit('in')
	ipatx = 4
	ipaty = 4
	xlen = 0.05
	ylen = 0.05
	call shadet(0.0,0.0,0.0,1.0,1.0,1.0,ipatx,ipaty,xlen,ylen)
	call plot(1.0,1.0,-3)
	call factor(0.5)
	call shadet(0.0,0.0,0.0,1.0,1.0,1.0,ipatx,ipaty,xlen,ylen)
	call factor(1.0)
	call plot(1.0,0.0,-3)
	ipatx = 0
	ipaty = 7
	call shadet(0.0,0.0,0.0,1.0,1.0,1.0,ipatx,ipaty,xlen,ylen)
	call plot(-2.0,-1.0,-3)
	call plot(0.0,2.0,-3)
	call factor(2.0)
	ipatx = 0
	call shadet(0.0,0.0,0.0,1.0,1.0,1.0,ipatx,ipaty,xlen,ylen)
	call pend()
	stop
	end
