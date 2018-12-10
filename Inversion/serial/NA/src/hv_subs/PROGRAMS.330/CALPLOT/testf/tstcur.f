c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: TSTCUR                                                c
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
	real*4 x(22),y(22)
	character*2 ic
	character*10 ttlx,ttly
	character*40 outstr
	call pinitf('INTER')
	call gunit('in')
	n = 20
	do 100 i=1,n
		x(i) = i
		y(i) = i*i*i
  100	continue
	do 1000 i=1,4
		call plot(2.0,1.0,-3)
		if(i.eq.1)then
			nocx=2
			nocy=3
			fct = 1.0
			inteq = 0
			lintyp = 0
		elseif(i.eq.2)then
			nocx=0
			nocy=2
			fct = 0.8
			inteq = i
			lintyp = -1
		elseif(i.eq.3)then
			nocy=0
			nocx=2
			fct = 1.0
			inteq = 1
			lintyp = +1
		elseif(i.eq.4)then
			nocx=0
			nocy=0
			fct = 0.8
			inteq = i
			lintyp = 0
		endif
		call factor(fct)
		yaxlen=6.0
		xaxlen=5.0
		call pltscl(x,xaxlen,n,x1,deltax,nocx)
		call pltscl(y,yaxlen,n,y1,deltay,nocy)
		ttlx = 'X-AXIS'
		mtx = 6
		ttly = 'Y-AXIS'
		mty = 6
		call algaxe(xaxlen,yaxlen,nocx,nocy,ttlx,ttly,mtx,mty,
     1		x1,y1,deltax,deltay)
		if(i.le.2)then
		call pltlog(x,y,n,x1,y1,deltax,deltay,lintyp,inteq,ht,
     1		nocx,nocy)
		else
		call pltlgd(x,y,n,x1,y1,deltax,deltay,lintyp,inteq,ht,
     1		nocx,nocy,8*i -1 , 0.10)
		endif
		call curuxy(xx,yy,x1,y1,deltax,deltay,nocx,nocy,ic)
		write(outstr,'(2f10.3,1x,a)')xx,yy,ic
		call gwrtxt(0.0,5.0,outstr,0)
		call factor(1.0)
		call plot(-2.0,-1.0,-3)
		call frame()
 1000	continue
	call pend()
	stop
	end
