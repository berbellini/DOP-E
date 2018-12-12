	call pinitf('WIDTH.PLT')
	call gunit('in')
	do 1000 ifont=1,4
		call gfont(ifont)
		do 100 i = 16,111
			ix = mod((i-1),16)
			iy = (i-1)/16
			xx = 1.0 + ix*0.5
			yy = 7.0 - (iy-1)*1.0
			ang = 0.0
			call gwidth(0.05)
			call symbol(xx,yy,0.5,char(i),ang,-1)
			call gwidth(0.0)
  100		continue
		call frame()
 1000	continue
	call pend()
	end
