c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: SCALE                                                 c
c                                                                     c
c      COPYRIGHT (C)  1997 R. B. Herrmann                             c
c                                                                     c
c      Department of Earth and Atmospheric Sciences                   c
c      Saint Louis University                                         c
c      3507 Laclede Avenue                                            c
c      St. Louis, Missouri 63103                                      c
c      U. S. A.                                                       c
c                                                                     c
c---------------------------------------------------------------------c
        subroutine gscale(x,xaxlen,n,inc)
c       purpose of scale is to find extreme values of array x
c       and set parameters for a nice plot by a combined use
c       of axis and line.
c-----
c       x       array to be scaled dimension n*inc+inc +1
c       xaxlen  axis length in inches into which data must fit
c       n       total number of points to be fit, spaced inc apart
c       inc     in the original array, each |inc|'th point is
c               used for scaling. Others are ignored. If array is
c               not well bounded, parameter will not be the best
c               but will suffice for the actual points used
c               The idea is to use the line routine to get a fast
c               plot
c               inc < 1, the first value is assumed to be a maximum
c-----
c       The object of scale is to define two parameters firstx and 
c       deltax
c       inc > 0 firstx is a min and deltax is positive
c       inc < 0 firstx is a max and deltax is negative
c-----
        dimension x(1),save(7)
c        data save/1.0,2.0,4.0,5.0,8.0,10.0,20.0/
        save(1) = 1.0
        save(2) = 2.0
        save(3) = 4.0
        save(4) = 5.0
        save(5) = 8.0
        save(6) = 10.0
        save(7) = 20.0
c-----
c       for a nice plot we want scale to start at some nice number
c       raised to some power
c-----
        k = iabs(inc)
        nn=n*k
c-----
c       get largest and smallest numbers in the array
c-----
        xmin = x(1)
        xmax=xmin
        do 100 i=1,nn,k
                xx=x(i)
                if(xx.gt.xmax)xmax=xx
                if(xx.lt.xmin)xmin=xx
  100   continue
        firstv = xmin
        deltav = (xmax-xmin)/xaxlen
c-----
c       make deltav a multiple of one of the save(i)
c-----
        if(deltav.le.0.0)then
                deltav=2.0*firstv
                if(deltav.le.0.0)deltav=1.0
            tmax = firstv + deltav
        else
                xlog=alog10(deltav)
c               if(xlog.lt.0.0)xlog=xlog-1.0
                i = xlog + 500
c-----
c assume computer never has floating point number outside the
c range > 10**500 or < 10**(-500)
c-----
                p=10.0**(i-500)
                deltav=deltav/p - 0.01
                is=0
                do 50 i=1,6
                        if(save(i).lt.deltav)is=i
   50           continue
c-----
c       construct final deltav
c-----
 1000   continue
            is = is + 1
            if(is.gt.7)go to 1001
                deltav = save(is)*p
                range = aint(xaxlen+0.49)*deltav
                xmid= 0.5*(xmin+xmax)
            if(xmid.lt.0.0)xmid = xmid - deltav/2.0
c-----
c       find scaled element nearest the center
                xscal = deltav*aint(xmid/(0.99*deltav))
            firstv= xscal - 0.5*deltav*aint(xaxlen+0.49)
c-----
c       test total range
c-----
                tmax = firstv + range
            if(tmax .lt. 1.01*xmax .or. firstv.gt.0.99*xmin)
     1          go to 1000
                if(xmin*firstv.le.0.0)firstv=0.0
        endif
 1001   continue
c-----
c       invoke sign of inc
c-----
        if(inc.gt.0)then
                x(nn+1)=firstv
                x(nn+k+1)=deltav
        else
                x(nn+1)=tmax
                x(nn+k+1)= -deltav
        endif
        return
        end
