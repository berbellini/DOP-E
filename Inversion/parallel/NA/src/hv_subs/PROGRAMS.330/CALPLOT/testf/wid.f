
	call pinitf('WID.PLT')
	call gunit('in')
	do 100 i = 1, 2
		if(i.eq.1)then
			call gunit('in')
			call gwidth(0.25)
			call plot(3.0,3.0,3)
			call plot(5.0,5.0,2)
			call plot(3.0,3.0,3)
			call plot(5.0,4.0,2)
			call plot(3.0,3.0,3)
			call plot(4.0,5.0,2)
			call plot(3.0,3.0,3)
			call plot(3.0,5.0,2)
			call plot(3.0,3.0,3)
			call plot(5.0,3.0,2)
		else
			call gunit('cm')
			call gwidth(1.0)
			call plot(1.0,1.0,3)
			call plot(6.0,6.0,2)
			call plot(1.0,1.0,3)
			call plot(6.0,3.0,2)
			call plot(1.0,1.0,3)
			call plot(3.0,6.0,2)
			call plot(1.0,1.0,3)
			call plot(6.0,1.0,2)
			call plot(1.0,1.0,3)
			call plot(1.0,6.0,2)
		endif
  100	continue
	call pend()
	end
		
