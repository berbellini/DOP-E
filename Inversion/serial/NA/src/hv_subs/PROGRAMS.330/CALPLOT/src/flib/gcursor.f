c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: GCURSOR                                               c
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

C     interface to subroutine dfcurs
C    + [c] (cmd[reference])
C     integer*4 cmd
C     end

        subroutine  gcursor(cmd)
        character cmd*(*)

        integer*4 icmd
        if(cmd(1:1) .eq. 'A' .or. cmd(1:1) .eq. 'a')then
c           Arrow
            icmd = 0
        else if(cmd(1:1) .eq. 'X' .or. cmd(1:1) .eq. 'x')then
c           Arrow
            icmd = 1
        else if(cmd(1:1) .eq. 'C' .or. cmd(1:1) .eq. 'c')then
c           Crosshairs
            icmd = 2
        else if(cmd(1:1) .eq. 'P' .or. cmd(1:1) .eq. 'p')then
c           Plus
            icmd = 3
        else if(cmd(1:1) .eq. 'B' .or. cmd(1:1) .eq. 'b')then
c           Box
            icmd = 4
        else if(cmd(1:1) .eq. 'R' .or. cmd(1:1) .eq. 'r')then
c           RubberBand
            icmd = 5
        else if(cmd(1:1) .eq. 'O' .or. cmd(1:1) .eq. 'o')then
c           Off
            icmd = 6
        else if(cmd(1:1) .eq. 'H' .or. cmd(1:1) .eq. 'h')then
c           Hyperbola
            icmd = 7
        else 
c           none
            icmd = 0
        endif
        call dfcurs(icmd)
        return
        end
