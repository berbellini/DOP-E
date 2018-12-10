	character*10 txt
	character*10 txt1
	character*80 tetx
	character*1 ic
	call pinitf('INTER')
	call gunit('in')
	call newpen(1)
	call gwrtxt(1.0,2.0,'HELLO WORLD   ',0)
	call gwrtxt(1.0,1.5,'HELLO WORLD   ',1)
	call gwrtxt(1.0,1.0,'HELLO WORLD   ',2)
	txt='TODAY    :'
	call newpen(2)
	call gwrtxt(3.0,5.0,txt,0)
	call gwrtxt(3.0,4.5,txt,1)
	call gwrtxt(3.0,4.0,txt,2)
	call grdtxt(txt1,10)
	call newpen(4)
	call gwrtxt(0.0,7.0,txt1,2)
c-----
c	put a cursor on the screen, get coordinates when a key is
c	hit. Then plot this at the position. Note that the (0,0)
c	of the character is its lower left corner
c-----
        call gcursor('CROSS')
	call currxy(xx,yy,ic)
c-----
c	echo character input -- beware when doing this with curuxy
c-----
c------
c	safety on input character
c-----
	if(ichar(ic(1:1)).lt.32)then
		if(ichar(ic(1:1)).eq.1)then
			call gwrtxt(xx,yy,'Left mouse button',0)
		else if(ichar(ic(1:1)).eq.2)then
			call gwrtxt(xx,yy,'Right mouse button',0)
		else if(ichar(ic(1:1)).eq.3)then
			call gwrtxt(xx,yy,'Center mouse button',0)
		else
			ic = ' '
			call gwrtxt(xx,yy,ic,0)
		endif
	else if(ichar(ic(1:1)).gt.126)then
		ic=' '
		call gwrtxt(xx,yy,ic,0)
	else
		call gwrtxt(xx,yy,ic,0)
	endif
c-----
c	put legend at the bottom
c-----
	write(tetx,1)xx,yy,ic(1:1)
    1	format('For input xx=',f10.2,' yy=',f10.2,' character=',a,' ')
	call gwrtxt(0.0,0.1,tetx,1)
	call gwrtxt(8.0,0.1,'CR TO CONTINUE',2)
	call grdtxt(txt1,1)
	call frame()
	call pend()
	end
	
