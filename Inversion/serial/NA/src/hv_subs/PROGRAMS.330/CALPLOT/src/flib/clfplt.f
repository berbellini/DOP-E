c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: CLFPLT                                                c
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

        subroutine curixy(ix,iy,ic)
        character ic*(*)
        integer *4 jx, jy
c-----
c     convert screen coordinates back into plot coordinates
c-----
            call cross(jx,jy,ic)
            ix = jx
            iy = jy
        return
        end

        subroutine curaxy(xx,yy,ic)
        character ic*(*)
c-----
c       returns coordinates of cursor position in user units (inch/cm) 
c       with respect to the absolute origin
c-----
        common/Zcplot/iunit,ifont
            call curixy(ix,iy,ic)
            xx = ix/1000.0
            yy = iy/1000.0
            if(iunit.eq.1)then
                xx = xx * 2.54
                yy = yy * 2.54
            endif
        return
        end

        subroutine currxy(xx,yy,ic)
        character ic*(*)
      common/Xcplot/xold,yold,xcur,ycur
        common/Ycplot/xstp,ystp
c-----
c       returns cursor position relative to current origin,
c       taking into account any user scaling ,e.g., call factor
c-----
        common/Zcplot/iunit,ifont
            call curaxy(ax,ay,ic)
            if(iunit.eq.1)then
                ax = ax * 0.3937
                ay = ay * 0.3937
            endif
            xx = (ax-xold)/xstp
            yy = (ay-yold)/ystp
        return
        end

