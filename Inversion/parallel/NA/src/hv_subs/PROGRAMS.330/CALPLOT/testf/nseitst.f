c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: NSEITST                                               c
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
c-----
c       this program tests an improvement to the
c       CALPLOT device drivers to permit gray scale plots under
c       certain circumstances. This will be useful for plotting
c       seismic traces together with a trace attribute modulating
c       either the shading or the color of the background
c
c       To do this we invoke solid rectangle shading with ipatx = ipaty = 0
c       but with the color in newpen(ipen) > 1000. The present gray
c       shading will fill solid with ipen=1000, and with an
c       areal density decreasing with increasing ipen
c-----
        parameter (NPTS=81)
        real*4 x(NPTS),y(NPTS),amp(NPTS)
        call pinitf('NSEITST.PLT')
        call gunit('in')
c-----
c       program to test seismic trace shading
c-----
c       fill up trace array and the attribute array
c-----
        ampmin = 1.0e+38
        ampmax =-1.0e+38
        tramin = 1.0e+38
        tramax =-1.0e+38
        do 100 i=0,NPTS-1
            x(i+1) = i*0.05
            amp(i+1) = abs(0.25 * cos(3.1415927*x(i+1)/2.))
            if(ampmin.gt.amp(i+1))ampmin=amp(i+1)
            if(ampmax.lt.amp(i+1))ampmax=amp(i+1)
            y(i+1) = amp(i+1)*sin(3.1415927*x(i+1)*5./2.)
            if(tramin.gt.abs(y(i+1)) )tramin=abs(y(i+1))
            if(tramax.lt.abs(y(i+1)) )tramax=abs(y(i+1))
  100   continue
        call puttrc(1.0,1.0,x,y,NPTS)
        call symbol(5.2,1.0,0.1,'TRACE',0.0,5)
        ishdrv = +1
        call putshd(1.0,1.5,x,amp,NPTS,ishdrv,ampmin,ampmax)
        call symbol(5.2,1.5,0.1,'SHADE',0.0,5)
        call plot(1.0,2.0,-3)
        do 400 i=2,NPTS
            ampav = (amp(i-1)+amp(i))/2
            xl = x(i-1)
            xh = x(i)
            xcol = (ampav - ampmin)/(ampmax -ampmin)
            if(i.eq.2)then
                ipen = 3
            else
                ipen = 2
            endif
            call plot(xl,xcol,ipen)
            call plot(xh,xcol,2)
  400   continue
        call plot(-1.0,-2.0,-3)
        call symbol(5.2,2.0,0.10,'ATTRIBUTE',0.0,9)
        ishdrv =+1
        call putshd(1.0,3.5,x,amp,NPTS,ishdrv,ampmin,ampmax)
        call puttrc(1.0,3.5,x,y,NPTS)
        call symbol(5.2,3.5,0.10,'ZERO ATTRIB=DARK',0.0,16)
        call symbol(5.2,3.3,0.10,'ZERO ATTRIB=RED ',0.0,16)
        ishdrv = -1
        call putshd(1.0,4.5,x,amp,NPTS,ishdrv,ampmin,ampmax)
        call puttrc(1.0,4.5,x,y,NPTS)
        call symbol(5.2,4.5,0.10,'MAX  ATTRIB=DARK',0.0,16)
        call symbol(5.2,4.3,0.10,'MAX  ATTRIB=RED ',0.0,16)
c-----
c       put out a nice gray scale according to ipen
c-----
        call graysc()
c-----
c       put up a trace, but have the shading keyed according to the
c       peak amplitude
c       This is tricky, since we must first search for local maxima
c-----
        call puttra(1.0,5.3,x,y,NPTS,tramax,tramin,0)
        call symbol(5.2,5.3,0.10,'POSITIVE  TRACE SHADE',0.0,21)
        call puttra(1.0,5.8,x,y,NPTS,tramax,tramin,1)
        call symbol(5.2,5.8,0.10,'NEGATIVE  TRACE SHADE',0.0,21)
        call puttra(1.0,6.3,x,y,NPTS,tramax,tramin,2)
        call symbol(5.2,6.3,0.10,'TWO-SIDED TRACE SHADE',0.0,21)
        call pend()
        stop
        end

        subroutine puttrc(x0,y0,x,y,n)
        real*4 x(1), y(1)
        call plot(x0,y0,-3)
        call shdsei(0.0,0.0,0,1,0)
        do 200 i=1,n
            xx = x(i)
            yy = y(i)
            if(i.eq.1)then
                call plot(xx,yy,3)
            else
                call plot(xx,yy,2)
            endif
  200   continue
        call shdsei(0.0,0.0,0,0,0)
        call plot(-x0,-y0,-3)
        return
        end

        subroutine puttra(x0,y0,x,y,n,ampmx,ampmn,iplmn)
        real*4 x(1), y(1)
        parameter (NPTS=81)
        integer kolor(NPTS)
c-----
c       The idea here is to follow the trace output, saving the
c       current maximum. At a zero crossing, we output the
c       color corresponding to this maximum, and reset
c       we must have the color attribute defined beforehand
c
c       note that the trace extrema, ampmx and ampmn,
c       are used to color the trace, not the extrema of the attribute
c-----
        ampmax = abs(y(1))
        iold = 1
        do 100 i=2,n
            if(abs(y(i)) .gt.ampmax)ampmax = abs(y(i))
            if(sign(1.0,y(i)) .ne. sign(1.0,y(i-1)))then
                color =  1100.0 - 100.0*(ampmax-ampmn)/
     1          (ampmx - ampmn)
                kolr = color
                if(kolr.lt.1000)kolr = 1000
                if(kolr.gt.1100)kolr = 1100
                do 101 j=iold,i-1
                    kolor(j) = kolr
  101           continue
                iold = i
                ampmax = 0.0
            endif
  100   continue
        kolor(n) = kolor(n-1)
        call plot(x0,y0,-3)
        call shdsei(0.0,0.0,0,1,iplmn)
        do 200 i=1,n
            xx = x(i)
            yy = y(i)
            if(i.eq.1)then
                call newpen(kolor(1))
                call plot(xx,yy,3)
            else
                if(sign(1.0,y(i)) .ne. sign(1.0,y(i-1)))then
                    fac = (y(i) - y(i-1))/(x(i) - x(i-1))
                    x2 = x(i-1) - y(i-1)/fac
                    call plot(x2,0.0,2)
                    call newpen(kolor(i-1))
c-----
c       force plotting of segment in this color
c-----
                    call shdsei(0.0,0.0,0,0,0)
c-----
c       initialize plotting of segment for next color
c-----
                    call shdsei(0.0,0.0,0,1,iplmn)
                    call plot(x2,0.0,3)
                    call newpen(kolor(i))
                    call plot(xx,yy,2)
                else
                    call plot(xx,yy,2)
                endif
            endif
  200   continue
        call shdsei(0.0,0.0,0,0,0)
c-----
c       now plot the trace in black so it can be seen
c-----
        call newpen(1)
        do 400 i=1,n
            xx = x(i)
            yy = y(i)
            if(i.eq.1)then
                call plot(xx,yy,3)
            else
                call plot(xx,yy,2)
            endif
  400   continue
        call plot(-x0,-y0,-3)
        return
        end

        subroutine putshd(x0,y0,x,amp,n,ishdrv,ampmin,ampmax)
        real*4 x(1),amp(1)
        call plot(x0,y0,-3)
        do 300 i=2,n
            xl = x(i-1)
            xh = x(i)
            ampav = (amp(i-1) + amp(i))/2
            xcol = (ampav - ampmin)/(ampmax-ampmin)
            xcol = 100.0*xcol  + 0.5
            ipen = 1000 + xcol
            if(ishdrv .lt.0)ipen = 1100 - xcol
            if(ipen .lt. 1000)ipen = 1000
            if(ipen .gt. 1100)ipen=1100
            yl = -0.25
            yh =  0.25
            ipatx = 0
            ipaty = 0
            xlen = 0.02
            ylen = 0.02
c-----
c       ipen = 1000 is red, 1100 = blue or 1000 = dark, 1100 = light halftone
c-----
            call newpen(ipen)
            call shader(xl,yl,xh,yh,ipatx,ipaty,xlen,ylen)
  300   continue
c-----
c       undo special shading
c-----
        call newpen(1)
        call plot(-x0,-y0,-3)
        return
        end

        subroutine graysc()
        call plot(7.0,5.0,-3)
        call plot(0.0,-4.0,2)
        call plot(0.5,-4.0,2)
        call plot(0.5,0.0,2)
        call plot(0.0,0.0,2)
        ipnmn=1000
        ipnmx=1100
        ipnd=ipnmx-ipnmn
        ipinc=5
        dy = -4.0/ipnd
        do 100 ipen=ipnmn,ipnmx-1,ipinc
            xl=0.0
            xh=0.5
            yl=(ipen-ipnmn)*dy
            yh=yl+dy*ipinc
            ipatx=0
            ipaty=0
            xlen=0.01
            ylen=0.01
            call newpen(ipen)
            call shader(xl,yl,xh,yh,ipatx,ipaty,xlen,ylen)
            call newpen(1)
            call number(xh+0.1,yh-dy/2.,0.07,real(ipen),0.0,-1)
  100   continue
        call plot(-7.0,-5.0,-3)
        return
        end
