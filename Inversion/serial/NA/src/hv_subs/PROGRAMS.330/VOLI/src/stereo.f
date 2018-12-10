      program stereo
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME II                                                      c
c                                                                     c
c      PROGRAM: STEREO                                                c
c                                                                     c
c      COPYRIGHT 1986 R. B. Herrmann                                  c
c                                                                     c
c      Department of Earth and Atmospheric Sciences                   c
c      Saint Louis University                                         c
c      221 North Grand Boulevard                                      c
c      St. Louis, Missouri 63103                                      c
c      U. S. A.                                                       c
c                                                                     c
c---------------------------------------------------------------------c
c
c     stereo -d dip -s slip -a strike -r rad -E -S
c---------------------------------------------------------------------c
c
      parameter (LER=0, LIN=5, LOT=6)
      real*4 dip,slip,strike,rad,x0,y0
      integer*4 iarea
c-----
c     machine dependent initialization
c-----
      call mchdep()
c-----
c     initialize plot interface
c-----
      call pinitf('STEREO.PLT')
      call plot(0.0,0.0,-3)
c-----
c     get command line information
c-----
      call gcmdln(dip,slip,strike,rad,iarea,x0,y0)
c-----
c     invoke stereo plot
c-----
      if(iarea.eq.0 .or. iarea.eq.1)then
            call strplt(dip,slip,strike,rad,iarea,x0,y0)
      else if (iarea.eq.2)then
            call blnkplt(rad,x0,y0)
      endif
c-----
c     close plot interface
c-----
      call pend()
      stop
      end


      subroutine strplt(dip,slip,strike,rad,iarea,x0,y0)
      real*4 dip,slip,strike,rad
      integer*4 iarea
      real*4 x(361),y(361),c(361),s(361)
c-----
c     perform stereographic projection
c       plotting planes of equal dip, slip
c       annotate circumference with tics marks each strike degrees
c
c       dip          dip increment in degrees
c       slip         slip increment in degrees
c       strike       strike increment in degrees
c       rad          radius of circle in inches
c       iarea        = 0 (default)  equal area
c                    .ne. 0  stereograpic projection
c-----
c-----
c     initialize a table of sin and cos
c-----
      do 100 i=0,360
            ang = i*3.1415927/180.0
            s(i+1) = sin(ang)
            c(i+1) = cos(ang)
  100 continue
c-----
c       draw dip planes
c-----
c-----
c     convert dip increment (real*4) into integer
c-----
      idip = dip
      if(dip.le.0)idip = 5
        do 200 id = 0, 180, idip
            do 201 i = 1,181
            if(id.eq.0)then
                  k = 181
                  x(i) = x0 + rad*s(i)
                  y(i) = y0 + rad*c(i)
            elseif(id.eq.90)then
                  k = 2
                  x(1) = x0
                  y(1) = y0 + rad
                  x(2) = x0
                  y(2) = y0 - rad
            elseif(id .eq. 180)then
                  k = 181
                  x(i) = x0 - rad * s(i)
                  y(i) = y0 + rad * c(i)
            else
                  k = 181
                  xs = c(i)
                  ys = s(i)*c(id+1)
                  zs = s(i)*s(id+1)
                  xj = atan2(sqrt(xs**2 + ys**2),zs)
                  if(iarea .eq. 0)then
                        r = 1.4142135*sin(xj/2.)
                  else
                        r = tan(xj/2.)
                  endif
                  alpha = atan2(ys,xs)
                  y(i) = y0 + rad*r*cos(alpha)
                  x(i) = x0 + rad*r*sin(alpha)
            endif
  201       continue
            call liner(k,x,y)
  200 continue
c-----
c     plot curves of equal slip
c-----
      islip = slip
      if(islip.le.0)islip = 5
      do 300 is=0, 180, islip
            do 301 id=0, 181, 1
                  k = 181
                  xs = c(is+1)
                  ys = s(is+1)*c(id+1)
                  zs = s(is+1)*s(id+1)
                  xj = atan2(sqrt(xs**2 + ys**2),zs)
                  if(iarea.eq.0)then
                        r = 1.4142135*sin(xj/2.)
                  else
                        r = tan(xj/2.)
                  endif
                  alpha=atan2(ys,xs)
                  y(id+1) = y0 + rad*r*cos(alpha)
                  x(id+1) = x0 + rad*r*sin(alpha)
  301       continue
            call liner(k,x,y)
  300       continue
c-----
c     put in azimuth tics
c-----
      istr = strike
      if(istr.le.0)istr = 5
      do 400 i=0, 359, istr
            xs = x0 + (rad + 0.20)*s(i+1)
            ys = y0 + (rad + 0.20)*c(i+1)
            if(i.gt.180)then
                  xs = xs - 0.20
            endif
            strk = i
            call number(xs,ys,0.025*rad,strk,0.0,-1)
  400 continue
c-----
c       now place a little circle at the center of the plot
c-----
        call curvit('CI',0.025*rad,x0,y0)
      return
      end

      subroutine blnkplt(rad,x0,y0)
      real*4 rad
      integer*4 iarea
      real*4 c(361),s(361)
c-----
c     perform stereographic projection
c       plotting planes of equal dip, slip
c       annotate circumference with tics marks each strike degrees
c
c       dip          dip increment in degrees
c       slip         slip increment in degrees
c       strike       strike increment in degrees
c       rad          radius of circle in inches
c       iarea        = 0 (default)  equal area
c                    .ne. 0  stereograpic projection
c-----
c-----
c     initialize a table of sin and cos
c-----
      do 100 i=0,360
            ang = i*3.1415927/180.0
            s(i+1) = sin(ang)
            c(i+1) = cos(ang)
  100 continue
c-----
c       circle
c-----
      ipen = 3
      do  i=0,360
            xs = x0 + (rad )*s(i+1)
            ys = y0 + (rad )*c(i+1)
           call plot(xs,ys,ipen)
           ipen = 2
      enddo
c-----
c     put in azimuth tics
c-----
      do 400 i=0, 359, 5
            if(i.eq.0 .or. i.eq.90 .or. i.eq.180 .or. i.eq.270)then
                  xs = x0 + (rad - 0.05)*s(i+1)
                  ys = y0 + (rad - 0.05)*c(i+1)
                  call plot(xs,ys,3)
                  xs = x0 + (rad + 0.05)*s(i+1)
                  ys = y0 + (rad + 0.05)*c(i+1)
                  call plot(xs,ys,2)
            else
                  xs = x0 + (rad - 0.05)*s(i+1)
                  ys = y0 + (rad - 0.05)*c(i+1)
                  call plot(xs,ys,3)
                  xs = x0 + (rad + 0.00)*s(i+1)
                  ys = y0 + (rad + 0.00)*c(i+1)
                  call plot(xs,ys,2)
            endif
C            if(mod(i,10).eq.0)then
C                  call plot(xs,ys,3)
C                  xs = x0 + (rad + 0.20)*s(i+1)
C                  ys = y0 + (rad + 0.20)*c(i+1)
C                  if(i.gt.180)then
C                        xs = xs - 0.20
C                  endif
C                  strk = i
C                  call number(xs,ys,0.025*rad,strk,0.0,-1)
C            endif
  400 continue
c-----
c       now place a little plus at the center of the plot
c-----
        call plot(x0-0.05,y0,3)
        call plot(x0+0.05,y0,2)
        call plot(x0,y0-0.05,3)
        call plot(x0,y0+0.05,2)
        call plot(x0,y0,3)
      return
      end

      subroutine liner(k,x,y)
      real*4 x(*),y(*)
c-----
c     draw the curve given by the x,y pairs
c-----
      x0 = x(1)
      y0 = y(1)
      call plot(x0,y0,3)
      do 100 i=1,k
            x0 = x(i)
            y0 = y(i)
            call plot(x0,y0,2)
  100 continue
      return
      end
      subroutine gcmdln(dip,slip,strike,rad,iarea,x0,y0)
      real*4 dip,slip,strike,rad,x0,y0
      integer*4 iarea
c-----
c     parse command line arguments and return control
c     parameters
c
c     requires subroutine mgtarg(i,name) to return
c           the i'th argument in the string name
c
c     and the function numarg() to return the number
c           of arguments excluding program name
c           The first argument is i = 1
c
c-----
      character*20 name
      integer*4 mnmarg
      dip = 5.0
      slip = 5.0
      strike = 5.0
      iarea = 0
      rad = 2.0
        x0 = 4.0
        y0 = 4.0

      nmarg = mnmarg()
      i = 0
   11 i = i + 1
      if(i.gt.nmarg)goto 13
            call mgtarg(i,name)
            if(name(1:3).eq.'-DD')then
                  i=i+1
                  call mgtarg(i,name)
                  read(name,'(bn,f20.0)')dip
            else if(name(1:3).eq.'-DA')then
                  i=i+1
                  call mgtarg(i,name)
                  read(name,'(bn,f20.0)')strike
            else if(name(1:3).eq.'-DS')then
                  i=i+1
                  call mgtarg(i,name)
                  read(name,'(bn,f20.0)')slip
            else if(name(1:3).eq.'-EQ')then
                  iarea = 0
            else if(name(1:3).eq.'-ST')then
                  iarea = 1
            else if(name(1:2).eq.'-B')then
                  iarea = 2
            else if(name(1:4).eq.'-RAD')then
                  i = i + 1
                  call mgtarg(i,name)
                  read(name,'(bn,f20.0)')rad
            else if(name(1:3).eq.'-X0')then
                  i = i + 1
                  call mgtarg(i,name)
                  read(name,'(bn,f20.0)')x0
            else if(name(1:3).eq.'-Y0')then
                  i = i + 1
                  call mgtarg(i,name)
                  read(name,'(bn,f20.0)')y0
            else if(name(1:2).eq.'-h' .or. name(1:2).eq.'-?')then
            call usage()
            endif
            go to 11
   13 continue
      return
      end

        subroutine usage()
        integer LOT
        parameter (LOT=6)
        write(LOT,*)
     1  'stereo -DD dip_inc -DA az_inc -X0 x0 -Y0 y0',
     1  ' [-EQ | -ST ] [-B] [ -h | - ?]'
        WRITE(LOT,*)
     1  '[create a stereo net]'
        WRITE(LOT,*)
     1  '-DD  dip_inc  (default 5 )   Dip increment in degrees'
        WRITE(LOT,*)
     1  '-DA   az_inc  (default 5 )   Azimuth increment in degrees'
        WRITE(LOT,*)
     1  '-DS slip_inc  (default 5 )   Slip increment in degrees'
        WRITE(LOT,*)
     1  '-X0   x0      (default 4.0)  Center of circle on page'
        WRITE(LOT,*)
     1  '-Y0   y0      (default 4.0)  Center of circle on page'
        WRITE(LOT,*)
     1  '-RAD rad      (default 2.0)  Radius of circle on page'
        WRITE(LOT,*)
     1  '-ST                          Stereographic projection'
        WRITE(LOT,*)
     1  '-EQ           (default)      Equal Area projection'
        WRITE(LOT,*)
     1  '-B                           Plot circle only for data'
        WRITE(LOT,*)
     1  '-h                           Usage query, but no execution'
        WRITE(LOT,*)
     1  '-?                           Usage query, but no execution'
        stop
        end
