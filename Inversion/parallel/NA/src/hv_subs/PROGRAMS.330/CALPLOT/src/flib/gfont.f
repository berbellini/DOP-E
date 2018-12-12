c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: GFONT                                                 c
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
cRevision History:
c        27 DEC 2000 -- saved the font in the Zcplot named common 
c---------------------------------------------------------------------c
        subroutine gfont(ifont)
        common/Zcplot/iunit,jfont
        integer*4 ifont
        integer*4 i4font
        integer jfont
        i4font = ifont
        jfont = ifont
        if(i4font.lt.0)i4font = 0
        call dffont(i4font)
        return
        end
