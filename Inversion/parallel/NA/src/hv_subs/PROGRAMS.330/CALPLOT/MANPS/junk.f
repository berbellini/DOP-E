	real*4 x(5), y(5)
	data x / 7.0, 11.0, 11.0, 10.0,  9.0 /
	data y / 6.0,  7.0, 11.0, 12.0,  11.0/

	integer i
	real xx, yy
c	call pinitf('LINE.PLT')
	call pinit()
	call newpen(1000)
	call shader(8.0,1.0,12.0,2.0,0,0,0.01,0.01)
	call newpen(1050)
	call shadet(1.0,7.0,1.5,9.0,2.0,6.0,0,0,0.01,0.01)
	call newpen(1100)
	call shadep(5,x,y)
	call newpen(1)
	call plot(-10.0,-8.0,3)
	call plot( 20.0,16.0,2)
	do 100 i=-10,20,1
		xx = i
		yy = i
		call number(1.0,xx,0.20,real(i),0.0,1)
		call number(xx,1.0,0.20,real(i),0.0,1)
  100	continue
	call box(0.5,0.5,9.5,7.5)
	call box(-0.5,-0.5,10.5,8.5)
	call newpen(2)
	call box(0.0,0.0,10.0,8.0)
	call pend()
	end

	subroutine box(xl, yl, xh,yh)
	call plot(xl,yl,3)
	call plot(xh,yl,2)
	call plot(xh,yh,2)
	call plot(xl,yh,2)
	call plot(xl,yl,2)
	call plot(xl,yl,3)
	return
	end
