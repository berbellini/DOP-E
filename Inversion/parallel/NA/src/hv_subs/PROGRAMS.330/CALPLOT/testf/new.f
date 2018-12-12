

	character strout*100
	
	integer*4 narr	
	integer*4 i
	integer*4 ix, iy
	character txt*20
	character txt1*10
	character tetx*80

	character ch*2
	real*4  xv, yv
	integer HasMouse 
	real*4 XminDev, YminDev, 
     1		XmaxDev, YmaxDev, XminClip, 
     2		YminClip, XmaxClip, YmaxClip
	integer*4 Color
	real*4 xarr(7), yarr(7)
	data xarr/ 1.0, 2.0, 3.0, 2.8, 2.5, 2.0, 1.0/
	data yarr/ 1.0, 1.0, 2.0, 2.5, 2.0, 3.0, 2.0/

	narr = 7
	call ginitf('INTER','NEW')
c-----
c	/* put out a sequence of tics */
c-----
	call newpen(1)
	do 1000 i=0, 100, 5
		if(mod(i,10).ne.0)then
			call plot(i*0.1,0.0 ,3)
			call plot(i*0.1,0.15,2)
		else 
			call plot(i*0.1,0.0 ,3)
			call plot(i*0.1,0.3 ,2)
			call number(i*0.1-0.07,0.35,0.14,(0.1*i),0.0,-1)
		endif
 1000	continue
	do 2000 i=0 , 100, 5
		if(mod(i,10).ne.0)then
			call plot(0.0 ,i*0.1,3)
			call plot(0.15,i*0.1,2)
		else
			call plot(0.0 ,i*0.1,3)
			call plot(0.3 ,i*0.1,2)
			call number(0.35,i*0.1-0.07,0.14,(0.1*i),0.0,-1)
		endif
 2000	continue
c-----
c	/* draw a box for orientation of X and Y axes */
c-----
	call newpen(1)
	call plot(0.0,0.0,3) 
	call plot(1.0,0.0,2) 
	call plot(1.0,1.0,2) 
	call plot(0.0,1.0,2)
	call plot(0.0,0.0,2)
	call symbol(0.5,0.10,0.10,'X',0.0,1)
	call symbol(0.10,0.5,0.10,'Y',0.0,1)
	call newpen(2)
	call plot(0.0,0.0,3) 
	call plot(2.0,0.0,2) 
	call plot(2.0,1.5,2) 
	call plot(0.0,1.5,2)
	
c-----
c	/* draw clip box */
c-----
	call newpen(1000)
 
	call plot (2.0,2.0,3)
	call plotd(3.0,2.0,9,0.05)
	call plotd(3.0,4.0,9,0.05)
	call plotd(2.0,4.0,9,0.05)
	call plotd(2.0,2.0,9,0.05)

	call newpen(4)
c-----
c	/* plot a border around the polygon - note this is
c		not a polygon clip */
c-----
	call gwidth(0.05)
	do 3000 i=1, narr
		if(i .eq. 1 )then
			call plot(xarr(i),yarr(i),3)
		else
			call plot(xarr(i),yarr(i),2)
		endif
 3000	continue
	call plot(xarr(1),yarr(1),2)
c-----
c	/* turn on clip box */
c-----
	call gclip('on',2.0,2.0,3.0,4.0)

c-----
c	/* shade a triangle */
c-----
	call shadet(1.0,1.0,2.0,5.0,1.0,0.5,0,0,0.05,0.05)

c-----
c	/* now do a polygon fill in green */
c-----
	call newpen(1050)
	call shadep(narr,xarr,yarr)
	call newpen(4)
	call gwidth(0.001)
c-----
C/* crudely get clip region */
c-----
	do 4000 i=250,10000,250
		yv = real(i)/1000.0
		call plot(0.0,yv,3) 
		call plot(10.0,yv,2)
 4000	continue

	
c-----
c	/* draw a straight line to see how it is clipped */
c-----
	call newpen(5)
	call gwidth(0.10)
	call plot(0.5,1.5,3)
	call plot(4.0,4.0,2)
	call gwidth(0.001)
	call gclip('on',0.0,0.0,5.0,5.0)
	call plot(-10.0,-10.0,3)
	call plot(5.0,10.0,2)
	call plot(40.0,40.0,3)
	call gclip('off',-1.0,-1.0,-1.0,-1.0)
	call box(5.0,4.0,8.0,6.0)
	call gclip('on',5.0,4.0,8.0,6.0)
	call gmesg('Crosshair Curso in Box Arrow Outsider')
	call gcursor('Cross')
	call cross(ix,iy,ch)
    1	format('ix: ',i5,' iy: ',i5,' c:',i5)
	write(strout,1)ix, iy, ichar(ch(1:1))
	call gwrtxt(5.5,5.0,strout,1)

	call ginfo(HasMouse, XminDev, YminDev, 
     1		XmaxDev, YmaxDev, XminClip, 
     1		YminClip, XmaxClip, YmaxClip, Color)
    2	format('HasMouse = ',i5)
	write(strout,2)HasMouse
	call gwrtxt(0.5,7.5,strout,1)
    3	format('Device Limits: (',f10.3,',',f10.3,') - (',
     1	 f10.3,',',f10.3,')')
	write(strout,3)XminDev,YminDev,XmaxDev,YmaxDev
	write(6,'(a)')strout
	call gwrtxt(0.5,7.0,strout,1)
    4	format('Clip Limits: (',f10.3,',',f10.3,') - (',
     1	 f10.3,',',f10.3,')')
	write(strout,4)XminClip,YminClip,XmaxClip,YmaxClip
	call gwrtxt(0.5,6.5,strout,1)
	write(strout,5)Color
    5	format('Color = ',i5)
	call gwrtxt(0.5,6.0,strout,1)
c-----
c	/* turn off the clip region to see if all is recovered */
c-----

	call gclip('off',5.0,4.0,8.0,6.0)
	call box(5.0,6.9,8.0,7.9)
	call gmesg('Arrow Cursor')
	call gcursor('Arrow')
		call cross(ix,iy,ch)
	write(strout,1)ix, iy, ichar(ch(1:1))
	call gwrtxt(5.5,5.0,strout,1)
		call newpen(2)
		call gwrtxt(5.5,7.5,strout,1)
	call gmesg('Crosshair Cursor')
	call gcursor('Cross')
		call cross(ix,iy,ch)
		write(strout,1)ix, iy, ichar(ch(1:1))
		call newpen(3)
		call gwrtxt(5.5,7.0,strout,1)
	call gmesg('Plus Cursor')
	call gcursor('Plus')
		call cross(ix,iy,ch)
		write(strout,1)ix, iy, ichar(ch(1:1))
		call newpen(3)
		call newpen(4)
		call gwrtxt(5.5,6.5,strout,1)
	call gmesg('XOR Arrow Cursor')
	call gcursor('XORArrow')
		call cross(ix,iy,ch)
		write(strout,1)ix, iy, ichar(ch(1:1))
		call newpen(3)
		call gwrtxt(5.5,6.0,strout,1)
	call gwrtxt(3.0,2.0,'HELLO WORLD   ',0)
	call gwrtxt(3.0,1.5,'HELLO WORLD   ',1)
	call gwrtxt(3.0,1.0,'HELLO WORLD   ',2)
	call gwrtxt(3.0,0.5,'ENTER STRING  ',2)
	call gwrtxt(3.0,0.5,'              ',2)
	call plot(3.0,0.5,3)
	call gcursor('Plus')
	call grdtxt(txt1,10)
	call gwrtxt(3.0,0.2,txt1,2)

	call frame()
	call pend()
	end
	
	
	
	subroutine box( xl,  yl,  xh,  yh)
	call plot(xl,yl,3)
	call plot(xh,yl,2)
	call plot(xh,yh,2)
	call plot(xl,yh,2)
	call plot(xl,yl,2)
	call plot(xl,yl,3)
	return
	end
