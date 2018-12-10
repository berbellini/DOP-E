

	call ginitf('INTER','TXTMOD')
	call newpen(1050)
	call shader(5.5,6.5,7.0,2.0, 0, 0, 0.01, 0.01)
	call dorow(6.0)
	call newpen(3000)
	call dorow(4.8)
	call newpen(3001)
	call dorow(3.6)
	call dorow(3.6)
	call newpen(3001)
	call dorow(2.4)
	call pend()
	end

	subroutine dorow ( yy)
	real yy
	integer i
	integer fillcolor(3) 
	data fillcolor/ 1000, 1050, 1100 /
	integer linecolor(5) 
	data linecolor / 0, 1, 2, 3, 4 /
	data nfill /3/
	data nline /5/
	real xx, yv
	do i=0, nfill -1
		xx = 1.0 + i
		call newpen(fillcolor(i))
		call shader(xx,yy,xx+1.1 ,yy-1.0, 0, 0, 0.01, 0.01)
	enddo
	call newpen(2)
	call gwrtxt(6.0,yy     ,"This is a string - 0",0)
	call gwrtxt(6.0,yy-0.25,"This is a string - 1",1)
	call gwrtxt(6.0,yy-0.50,"This is a string - 2",2)
	call gwidth(0.05)
	do i=0,nline-1
		yv = yy - 0.1 - i*1.0/real(nline)
		call newpen(linecolor(i))
		call plot(0.5,yv,3)
		call plot(1.0+nfill  + 0.6,yv,2)
		
	enddo
	call gwidth(0.001)
	return
	end

