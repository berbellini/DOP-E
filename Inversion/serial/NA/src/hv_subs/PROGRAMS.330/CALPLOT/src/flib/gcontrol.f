c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: GRONTROL                                              c
c                                                                     c
c      COPYRIGHT (C)  2005 R. B. Herrmann                             c
c                                                                     c
c      Department of Earth and Atmospheric Sciences                   c
c      Saint Louis University                                         c
c      3507 Laclede Avenue                                            c
c      St. Louis, Missouri 63103                                      c
c      U. S. A.                                                       c
c                                                                     c
c---------------------------------------------------------------------c
c       send information to low level graphics
        subroutine gcontrol(type, p1, p2, p3, p4)
        integer type
        real p1, p2, p3, p4
        integer i1, i2, i3, i4
        i1 = p1 * 10000.0
        i2 = p2 * 10000.0
        i3 = p3 * 10000.0
        i4 = p4 * 10000.0
        call dfcontrol(type, i1, i2, i3, i4)
        return
        end
