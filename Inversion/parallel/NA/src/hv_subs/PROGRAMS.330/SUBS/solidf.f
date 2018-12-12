
        subroutine gsolid(x, y, h, isymb)
C-----
C   solid fill with 
C       symb    symbol
C       0   square
C       1   circle
C       2   triangle
C       3   diamond
C       4   octagon
C       5   inverted triangle
c-----
          
        real x, y, h
        integer isymb
        real dang, ang
        integer nside, i
        real rad
        real xp(17)
        real yp(17)
        
        jsymb = mod(isymb,6)
        if(jsymb.eq.0)then
c-----
C square
c-----
                rad = 1.414 * h
                nside = 4
                dang = 45.0
        else if(jsymb.eq.1)then
c-----
C circle
c-----
                rad =  h
                nside = 16
                dang = 0.0
        else if(jsymb.eq.2)then
c-----
C triangle
c-----
                rad =  1.5*h
                nside = 3
                dang = 90.0
        else if(jsymb.eq.3)then
c-----
C diamond
c-----
                rad =  1.5*h
                nside = 4
                dang = 0.0
        else if(jsymb.eq.4)then
c-----
C octagon 
c-----
                rad = 1.0823 * h
                nside = 8
                dang = 22.5
        else if(jsymb.eq.2)then
c-----
C inverted triangle
c-----
                rad =  1.5*h
                nside = 3
                dang = 270.0
        endif
        do 100 i=1,nside
            ang = dang + (i-1)*360.0/real(nside)
            ang = ang * 0.017453293
            xp(i) = x + rad * cos(ang)
            yp(i) = y + rad * sin(ang)
 100    continue
        call shadep(nside,xp,yp)
        return
        end

        subroutine gpoly(x, y, h, isymb)
C-----
C   draw polygon fill with 
C       symb    symbol
C       0   square
C       1   circle
C       2   triangle
C       3   diamond
C       4   octagon
C       5   inverted triangle
c-----
          
        real x, y, h
        integer isymb
        real dang, ang
        integer nside, i
        real rad
        real xp(17)
        real yp(17)
        
        jsymb = mod(isymb,6)
        if(jsymb.eq.0)then
c-----
C square
c-----
                rad = 1.414 * h
                nside = 4
                dang = 45.0
        else if(jsymb.eq.1)then
c-----
C circle
c-----
                rad =  h
                nside = 16
                dang = 0.0
        else if(jsymb.eq.2)then
c-----
C triangle
c-----
                rad =  1.5*h
                nside = 3
                dang = 90.0
        else if(jsymb.eq.3)then
c-----
C diamond
c-----
                rad =  1.5*h
                nside = 4
                dang = 0.0
        else if(jsymb.eq.4)then
c-----
C octagon 
c-----
                rad = 1.0823 * h
                nside = 8
                dang = 22.5
        else if(jsymb.eq.2)then
c-----
C inverted triangle
c-----
                rad =  1.5*h
                nside = 3
                dang = 270.0
        endif
        do 100 i=1,nside
            ang = dang + (i-1)*360.0/real(nside)
            ang = ang * 0.017453293
            xp(i) = x + rad * cos(ang)
            yp(i) = y + rad * sin(ang)
 100    continue
        ipen = 3
        do 1000 i=1,nside
            call plot(xp(i),yp(i),ipen)
            ipen = 2
 1000   continue
        call plot(xp(1),yp(1),ipen)
        call plot(xp(1),yp(1),3)
            
        return
        end
