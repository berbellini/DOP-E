c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: GINFO                                                 c
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

        subroutine  ginfo(HasMouse, XminDev, YminDev, 
     1      XmaxDev, YmaxDev, XminClip, 
     2      YminClip, XmaxClip, YmaxClip, Color)
        integer*4 HasMouse
        real*4 XminDev, YminDev, XmaxDev, YmaxDev
        real*4 XminClip, YminClip, XmaxClip, YmaxClip
        integer*4 Color

        integer*4  gHasMouse, gXminDev, gYminDev, gXmaxDev, gYmaxDev
        integer*4  gXminClip, gYminClip, gXmaxClip, gYmaxClip
        integer*4 gColor
        common/Xcplot/xold,yold,xcur,ycur
        common/Ycplot/xstp,ystp

        call dfinfo(gHasMouse, gXminDev, gYminDev, 
     1      gXmaxDev, gYmaxDev, gXminClip, 
     2      gYminClip, gXmaxClip, gYmaxClip,gColor)
        HasMouse = gHasMouse 
C       For Device boundaries, use absolute coordinates 
        XminDev = gXminDev 
        YminDev = gYminDev 
        XmaxDev = gXmaxDev 
        YmaxDev = gYmaxDev 
C       return clip coordinates with respect to current origin 
        XminClip = gXminClip /1000.0 
        YminClip = gYminClip /1000.0 
        XmaxClip = gXmaxClip /1000.0 
        YmaxClip = gYmaxClip /1000.0 
        XminClip = (XminClip-xold)/xstp
        YminClip = (YminClip-yold)/ystp
        XmaxClip = (XmaxClip-xold)/xstp
        YmaxClip = (YmaxClip-yold)/ystp
        Color = gColor
        return
        end
