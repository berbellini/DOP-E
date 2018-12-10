c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: GUNIT                                                 c
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
        subroutine gunit(str)
        common/Ycplot/xstp,ystp
        common/Zcplot/iunit,ifont
        character str*(*)
            if(str(1:2) .eq. 'cm' .or. str(1:2).eq.'CM')then
                junit = 1
            else
                junit = 0
            endif
            if(junit.ne.iunit)then
                if(junit.eq.1)then
                    xstp = xstp * 0.3937
                    ystp = ystp * 0.3937
                else
                    xstp = xstp / 0.3937
                    ystp = ystp / 0.3937
                endif
            endif
            iunit = junit
        return
        end
