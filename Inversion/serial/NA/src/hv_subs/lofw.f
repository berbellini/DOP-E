	integer function lofw(word)
c
	character*(*) word
c
	lw=len(word)
	k=0
	do i=1,lw
	  if (word(i:i).eq.' ') go to 99
	  k=k+1
	end do
99	lofw=k
c
	return
	end
