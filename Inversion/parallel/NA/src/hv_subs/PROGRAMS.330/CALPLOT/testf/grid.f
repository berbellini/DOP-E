
	integer NARR
	parameter(NARR=16)

	real xarr(NARR) 
	real yarr(NARR) 

	character strout*100

	
	integer i,j,ipen
	real  x1,y1,x2,y2,x3,y3,rad,ang
	call pinitf('GRID.PLT')
c-----
c	put out a sequence of tics 
c-----
	call newpen(1)
	do i=0,100,5
		if(mod(i,10) .eq.0)then
			call plot(i*0.1, 0.0 ,3)
			call plot(i*0.1, 8.0 ,2)
		endif
	enddo
	do i=0,80,5
		if(mod(i,10) .eq.0)then
			call plot( 0.0 ,i*.1,3)
			call plot(10.00,i*.1,2)
		endif
	enddo
c------
c	put in little squares 
c	draw a box for orientation of X and Y axes 
c-----

	ipen = 1000
	do j=0,7,2
		y1 = j 
		y2 = y1 + 1.0
		do i=0,9,2
			x1 = i 
			x2 = x1 + 1.0
			ipen = ipen + 5
			call newpen(ipen)
			call shader(x1,y1,x2,y2, 0, 0, 0.01, 0.01)
		enddo
	enddo
	call newpen(0)
	call plot( 9.0,7.0,3)
	call plot( 8.0,6.0,2)


	ipen = 1000
	do j=0,7,2
		y1 = j  + 1.0
		y2 = y1 + 1.0
		y3 = y1
		do i=0,9,2
			x1 = i  + 1.0
			x2 = x1 + 0.5
			x3 = x1 + 1.0
			ipen = ipen + 5
			call newpen(ipen)
			call shadet(x1,y1,x2,y2,x3,y3, 0, 0, 0.01,0.01)
		enddo
	enddo
	

	call newpen(2)
	do i=1,NARR
		ang = 6.2831853*(i-1)/(NARR-1)
		xarr(i) = 4.0 + 1.0*cos(ang)
		yarr(i) = 4.0 + 1.0*sin(ang)
	enddo
	call shadep(NARR,xarr,yarr)
		


	call pend()
	end

	subroutine box(xl, yl, xh, yh)
	real xl, yl, xh, yh
	call plot(xl,yl,3)
	call plot(xh,yl,2)
	call plot(xh,yh,2)
	call plot(xl,yh,2)
	call plot(xl,yl,2)
	call plot(xl,yl,3)
	return
	end

