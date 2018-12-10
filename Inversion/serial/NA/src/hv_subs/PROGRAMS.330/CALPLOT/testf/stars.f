	program star
	call pinitf('INTER')
	call dowheel()
	call pend()
	end

	subroutine dowheel()
	integer*4 i
	integer*4 color
	real*4 dang, ang
	real*4 rad
	real*4 x0, y0
c-----
c	define the center 
c-----
	call plot(5.0,4.0,-3)
	call gclip('on',-4.0,-2.0,4.0,2.0)
	x0 = 0.0
	y0 = 0.0
	dang = 0.0
	rad = 5.0
	do 1000 i=1, 1000
		x0 = 0.0 + 8.0*(1.0-0.02*random(100))
		y0 = 0.0 + 8.0*(1.0-0.02*random(100))
		rad = 0.5 * 0.01*random(100)
		dang = random(30)
		ipen = 1000 + random(100)
		call newpen(ipen)
		call dostar(x0,y0,rad,dang )
 1000	continue
	end

	subroutine  dostar(x0,y0,rad,dang )
	real*4 x0, y0, rad, dang
	parameter(NPOLY=10)
	real*4 xarr(NPOLY), yarr(NPOLY)

	integer*4 i
	real*4 degrad 
	real*4  ang
	real*4 radius
	degrad = 3.1415927/180.0
	do 1000 i=1, NPOLY
		ang = dang + (i-1)*36
		ang = ang * degrad
		if(mod(i,2).eq.1)then
			radius = rad
		else
			radius = 0.3819*rad
		endif
		xarr(i) = x0 + radius*cos(ang)
		yarr(i) = y0 + radius*sin(ang)
 1000	continue
	call shadep(NPOLY, xarr, yarr)
	return 
	end

	function random(i)
		random = i*ran2(idum)
	return
	end

      FUNCTION ran2(idum)
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     *IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
     *NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software 3,5.
