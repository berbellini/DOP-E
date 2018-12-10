c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: PLTTST                                                c
c                                                                     c
c      COPYRIGHT (C)  1986 R. B. Herrmann                             c
c                                                                     c
c      Department of Earth and Atmospheric Sciences                   c
c      Saint Louis University                                         c
c      221 North Grand Boulevard                                      c
c      St. Louis, Missouri 63103                                      c
c      U. S. A.                                                       c
c                                                                     c
c---------------------------------------------------------------------c
      program plttst
      call pinitf('PLTTST.PLT')
	call gunit('in')
      call symbol(5.0,0.5,0.25,'X',0.0,1)
      call symbol(0.5,5.0,0.25,'Y',0.0,1)
      do 100 i=1,40
      icolor = i
      call newpen(icolor)
      rmax= i*0.5
      call number(0.25,rmax,0.10,rmax,0.0,+1)
      call number(rmax,0.25,0.10,rmax,0.0,+1)
      do 100 j=0,90
      ang = 3.1415927*float(j)/180.
      xx = rmax*cos(ang)
      yy= rmax*sin(ang)
      if(j.gt.0)then
            ipen=2
      else
            ipen=3
      endif
      call plot(xx,yy,ipen)
  100 continue
      call plot(0.0,0.0,3)
      call plot(8.0,0.0,2)
      call plot(8.0,8.0,2)
      call plot(0.0,8.0,2)
      call plot(0.0,0.0,2)
      call frame()
      call newpen(1)
      do 400 i=1,40
      xmax= i*0.5
      ymax=xmax
      call number(0.25,ymax,0.10,ymax,0.0,+1)
      call number(xmax,0.25,0.10,xmax,0.0,+1)
      call plot(0.0,0.0,3)
      call plot(xmax,0.0,2)
      call plot(xmax,ymax,2)
      call plot(0.0,ymax,2)
      call plot(0.0,0.0,2)
  400 continue
      call symbol(5.0,0.5,0.25,'X',0.0,1)
      call symbol(0.5,5.0,0.25,'Y',0.0,1)
      call pend()
      stop
      end
