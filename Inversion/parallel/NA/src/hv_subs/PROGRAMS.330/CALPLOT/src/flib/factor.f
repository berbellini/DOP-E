c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: FACTOR                                                c
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
      subroutine factor (fact) 
c-----
c       multiply all pen moves by fact
c       useful in reducing size of plots, just
c       do a call factor immediately after
c       the call plots or call pinit
c----
        common/Zcplot/iunit,ifont
        if(fact.le.0.0)then
            factnw=1.0
        else
            factnw = fact
        endif
        if(iunit.eq.1)then
            factnw = factnw * 0.3937
        endif
        call plot(factnw,factnw,1001)
        return 
        end 
