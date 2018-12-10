c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: GPHTXT                                                c
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

        

c-----
c       routines to implement string output to terminals
c       on some terminals this will use the terminal character 
c       generator, e.g.,
c       tektronix, on others it will use a bit map, e.g., PC, and others
c       will use a vector draw
c-----

            subroutine gwrtxt(xx,yy,text,flgbln)
            real*4 xx,yy
            character text*(*)
            integer*4 flgbln
            integer*4 ipen
            integer  lstr
            ipen = 2000 + flgbln
            lstr = len(text)

                call plot(xx,yy,3)
                call newpen(ipen)
                call dfgott(lstr,text)
            return
            end

            subroutine grdtxt(text,ltext)
            character text*(*)
            integer   ltext
            integer*4 lst
                text=' '
                lst = ltext
                call dfgint(lst,text)
            return
            end
