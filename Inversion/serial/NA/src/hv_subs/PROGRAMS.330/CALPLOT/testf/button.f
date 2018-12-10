	subroutine mkmenu()
	common/menu/nm,mxl(5),myl(5),mxh(5),myh(5),mstr(5),mactio(5)
	integer*4 nm
	real*4 mxl, myl, mxh, myh
	character mstr*8
	integer*4 mactio
	nm = 5
	do 100 i=1,nm
		mxl(i) = -1
		myl(i) = -1
		mxh(i) = -1
		myh(i) = -1
		mactio(i) = i
  100	continue
	mstr(1) = '  Mode  '
	mstr(2) = ' Select '
	mstr(3) = ' Delete ' 
	mstr(4) = '  Copy  '
	mstr(5) = '  Exit  '
	return
	end


	program button
	common/ctrl/black
	integer*4 black
	common/menu/nm,mxl(5),myl(5),mxh(5),myh(5),mstr(5),mactio(5)
	integer*4 nm
	real*4 mxl, myl, mxh, myh
	character mstr*8
	integer*4 mactio

	real*4 xl,yl,xh,yh
	integer*4 i
	integer*4 cmd
	real*4 xv, yv
	character c*2
	character ostr*8
	call pinitf('INTER')
	call mkmenu()
	border = 8.0/60.0
	black = 1

	call drawit()
	
	do 1000 i=1,nm
		ostr = mstr(i)
		call draw_button(8.5,3.0 -i*0.4,ostr,xl,yl,xh,yh)
		mxl(i) = xl
		myl(i) = yl
		mxh(i) = xh
		myh(i) = yh
 1000	continue
c-----
c	event loop for values 
c-----
	cmd = -1
 2000	continue
		call gcursor('Arrow')
		call curaxy(xv, yv, c)
		cmd = -1
		do 2100 i=1, nm
			if(inside(xv,yv,mxl(i),myl(i),
     1				mxh(i),myh(i)).ne.0)then
     				cmd = mactio(i)
				call newpen(1)
				ostr = mstr(i)
				call gmesg(ostr)
			endif
 2100		continue
		if(cmd .eq. 5)then
			goto 9999
		else if(cmd .eq. 4)then
			call pinitf('button.plt')
			call drawit()
			call pend()
			call pinitf('INTER')
			call gcursor('Arrow')
		endif
	go to 2000
 9999	continue
	call gmesg('Closing Session - Any Key to End')
	call pend()
	end

	integer function inside( xv,  yv,  xlb, ylb,  xhb,  yhb)
	real*4 xv, yv, xlb, ylb, xhb, yhb
	if(xv .ge. xlb .and. xv .le. xhb .and. 
     1		yv .ge. ylb .and. yv .le. yhb)then
		inside = 1
	else
		inside = 0
	endif
	return
	end

	subroutine draw_button( xl,  yl, str,  xlw,  ylw, xup,  yup)
	real*4 xl, yl
	character str*(*)
	real*4 xlw, ylw, xup, yup
	common/ctrl/black
	integer*4 black
	integer*4 lstr
	integer*4 nc
	integer*4 IBcolr
	integer*4 IBLigh
	integer*4 IBDark
	integer*4 IBFore
	real*4 bdr, title, x(6), y(6)
	bdr = 0.04
	title = 20.0/60.0
	lstr = len(str)
c-----
c	add a space
c-----
	nc = lstr*8
	IBDark  = 1100
	IBcolr       = 1090
	IBLigh = 1070
c-----
c	color
c-----
	if(black.eq.1)then
		IBFore  =    0
	else
		IBFore  =    1
	endif
	call newpen(IBDark)
	x(1) = xl
	y(1) = yl
	x(2) = xl + lstr*0.125
	y(2) = yl
	x(3) = xl + lstr*0.125
	y(3) = yl + title 
	x(4) = xl + lstr*0.125 - bdr
	y(4) = yl + title - bdr 
	x(5) = xl + lstr*0.125 - bdr
	y(5) = yl + bdr 
	x(6) = xl + bdr
	y(6) = yl + bdr 
	call shadep(6,x,y)

	call newpen(IBLigh)
	x(1) = xl
	y(1) = yl
	x(2) = xl + bdr
	y(2) = yl + bdr 
	x(3) = xl + bdr
	y(3) = yl + title - bdr 
	x(4) = xl + lstr*0.125 - bdr
	y(4) = yl + title - bdr 
	x(5) = xl + lstr*0.125 
	y(5) = yl + title  
	x(6) = xl  
	y(6) = yl + title  
	call shadep(6,x,y)

	call newpen(IBcolr)
	x(1) = xl + bdr
	y(1) = yl + bdr
	x(2) = xl + lstr*0.125 - bdr
	y(2) = yl + bdr
	x(3) = xl + lstr*0.125 - bdr
	y(3) = yl + title - bdr
	x(4) = xl + bdr
	y(4) = yl + title - bdr
	call shadep(4,x,y)
	call newpen(IBFore)
	call gwrtxt(xl+0.0625,yl+0.015,str,0)
	xlw = xl
	ylw = yl
	xup = xl + lstr*0.125
	yup = yl + title
	return
	end

	subroutine drawit()
	integer*4 i, ipen
	real*4 xx, yy
	call box(3.0,4.0,6.0,6.0)
	call gclip('on',3.0,4.0,6.0,6.0)
	ipen = 3
	call newpen(4)
	do 1000 i=0, 100, 1
		xx = 3.0 + (i)*(6.0-3.0)/100.0
		yy = 5.0 + 2.0*sin(10.*xx)+cos(10.*xx)
		call plot(xx,yy,ipen)
		ipen = 2
 1000	continue
	call gclip('off',3.0,4.0,6.0,6.0)
	return
	end

	subroutine box(xl, yl, xh, yh)
	real*4 xl, yl, xh, yh
	call plot(xl,yl,3)
	call plot(xh,yl,2)
	call plot(xh,yh,2)
	call plot(xl,yh,2)
	call plot(xl,yl,2)
	call plot(xl,yl,3)
	return
	end

	function lgstr(str)
c-----
c	function to find the length of a string
c	this will only be used with file system path names
c	thus the first blank 
c	indicates the end of the string
c-----
	character*(*) str
	integer*4 lgstr
	n = len(str)
	lgstr = 1
	do 1000 i=n,1,-1
		lgstr = i
		if(str(i:i).ne.' ')go to 100
 1000	continue
  100   continue
        return
	end
