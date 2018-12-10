c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: PEND                                                  c
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
        subroutine pend()
c-----
c----- TERMINATE PLOT
c-----
c-----
c      First reset width and pen to the default
c-----
       call gwidth(0.0)
       call newpen(1)
        call dfclos(0)
        return
        end

        subroutine gend(mode)
        integer mode

c----- TERMINATE PLOT with mode = 0 for normal, mode = 1
c       for automatic close of window
c-----

c-----
c      First reset width and pen to the default
c-----
       call gwidth(0.0)
       call newpen(1)
        if(mode.lt.0 .or. mode.gt.1)then
            call dfclos(0)
        else
            call dfclos(mode)
        endif
        return
        end
