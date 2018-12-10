c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: TABL                                                  c
c                                                                     c
c      COPYRIGHT (C)  2004 R. B. Herrmann                             c
c                                                                     c
c      Department of Earth and Atmospheric Sciences                   c
c      Saint Louis University                                         c
c      221 North Grand Boulevard                                      c
c      St. Louis, Missouri 63103                                      c
c      U. S. A.                                                       c
c                                                                     c
c---------------------------------------------------------------------c
c
c     this routine produces a symbol table, which shows the characters
c     available in the symbol routine.
c-----
	real znum(4)
	real z, xs, ys, x, y
	character*2 cc
	integer m, ia, ib, ndec, mdec
	data znum/ 10293.84756,0.019375,-204.86,-12345.6789/
	call pinitf('TABL.PLT')
	call gunit('in')
	call rect(0.4,0.4,9.6,7.5)
	call rect(0.6,0.6,9.4,7.1)
	call gfont(1)
	call symbol(2.15,7.25,0.15, 
     1		'characters available in symbol routine',0.0,38)
	z=0.0
	m=0
	x=0.6
	y=6.57
	xs=0.75
	ys=0.25
	call gfont(0)
	do 1000	ia=0,7
		do 1100 ib=0,13
			call number(x+.1,y+.18,.10,z,0.0,-1)
				cc(1:1) = char(m)
			call symbol(x+xs,y+ys ,.30 ,cc,0.0,-1)
				z=z+1.0
				m=m+1
				y=y-0.45
 1100		continue
		x=x+1.1
		call plot(x,0.6,3)
		call plot(x,7.1,2)
		y=6.57
		xs=.55
		ys=.05
 1000	continue
	call gfont(1)
	call symbol(2.05,0.45,.10,
     1 	'integer for use in symbol call shown to left of each symbol'
     2				,0.0,59)
c-----
c	begin a new page
c-----
	call frame()
c-----
c	the following tests the number subroutine for precision
c-----
	call symbol(1.0,2.5,.14,
     1		'example of number subroutine',90.0,28 )
	call symbol(2.62,7.5,.14,
     1		'call number(xx,yy,ht,fpn,ang,ndec)',0.0,34)
	call symbol( 1.5,7.2,.10,'ndec',0.0,4)
	call symbol( 2.0,7.2,.10,'number',0.0,6)
	call symbol( 4.5,7.2,.10,'ndec',0.0,4)
	call symbol( 5.0,7.2,.10,'number',0.0,6)
	call symbol( 7.5,7.2,.10,'ndec',0.0,4)
	call symbol( 8.0,7.2,.10,'number',0.0,6)
	y=7.0
	do 2000 ia=1,4
		do 2100 ib=1,11
			ndec = ib-6
			call number(1.5,y,.07,float(ndec),0.0,-1)
			call number(2.0,y,.07,znum(ia),0.0,ndec)

			ndec = 999 + ib
			mdec = ndec + 1000
			call number(4.5,y,.07,float(ndec),0.0,-1)
			call number(5.0,y,.07,znum(ia),0.0,ndec)
			call number(7.5,y,.07,float(mdec),0.0,-1)
			call number(8.0,y,.07,znum(ia),0.0,mdec)
			y=y-0.14
 2100		continue
		y=y-0.2
 2000	continue
	call pend()
	end

	subroutine rect( llx,  lly,  urx,  ury)
	real llx, lly, urx, ury
	call plot(llx,lly,3)
	call plot(urx,lly,2)
	call plot(urx,ury,2)
	call plot(llx,ury,2)
	call plot(llx,lly,2)
	call plot(llx,lly,3)
	return
	end
