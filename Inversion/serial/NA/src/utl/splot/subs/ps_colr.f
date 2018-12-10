C-------------------------------------------------------------ms-
      SUBROUTINE FILLPOLY(x,y,n)
C
C     Draws a filled polygon in current pen colour.
C
c-------------------------------------------------------------------

      COMMON/P00000/LPLOT,IROT,IL34,A,B,C,D,ASP,THET

      real*4  x(n),y(n)

      if(n.gt.2)then
         call plot(x(1),y(1),3)
	 do 10 i = 2,n
	    call plot(x(i),y(i),2)
 10      continue
      end if
      write(lplot,*)' closepath fill'

      return
      end

C    
C-------------------------------------------------------------ms-
      SUBROUTINE PEN(IPEN,ITHK)
C
C     allows choice of pen colour and thickness  
C    
C	modified 10/3/93 so that pen only controls colour
C	and ithk only controls pen thickness.
C
      COMMON/P00000/LPLOT,IROT,IL34,A,B,C,D,ASP,THET
      COMMON/PC0000/colr(256),colg(256),colb(256)
      COMMON/L00000/PSCA,xo,yo
      WRITE(LPLOT,*) ' stroke'
      if(il34.eq.0) then
        thik = 1.5
        thin = 1.0
        thin1= 0.7
        thin2= 0.5
        thin3= 0.35
      elseif(il34.eq.1) then
        thik = 1.0
        thin = 0.7
        thin1= 0.5
        thin2= 0.35
        thin3= 0.25
      endif
c
      thw = thin
      if(ithk.eq.-3)write(lplot,*)' thin4'
      if(ithk.eq.-2)write(lplot,*)' thin3'
      if(ithk.eq.-1)write(lplot,*)' thin2'
      if(ithk.eq.0)write(lplot,*)' thin1'
      if(ithk.eq.1)write(lplot,*)' thik1'
      if(ithk.eq.2)write(lplot,*)' thik2'
      if(ITHK.gt.1.and.ithk.ne.2)then
        thw = thw*float(ITHK)
        write(lplot,*) thw,' setlinewidth'
      end if
      ip = IPEN+1
      if(ipen.ge.0)then
         rc = colr(ip)
         gc = colg(ip)
         bc = colb(ip)
c        write(LPLOT,73) thw,rc,gc,bc
c						If pen is undefined then
c						set to pen0 (white)
c
         if(rc.lt.0..or.gc.lt.0..or.bc.lt.0.)then
            write(LPLOT,74)0
         else
            if(ipen.ge.0.and.ipen.lt.10)then
               write(LPLOT,74) ipen
            else if(ipen.lt.100)then
               write(LPLOT,75) ipen
            else if(ipen.lt.1000)then
               write(LPLOT,76) ipen
            end if
         end if
      else
         write(LPLOT,*) ' 0.0 setgray'
      end if
      xv=psca*xo
      yv=psca*yo
c     write(LPLOT,fmt='(f9.3,1x,f9.3,a)') xv,yv,' pM'
 71   format(1x,f5.2,a)
 73   format(1x,f5.2,' setlinewidth ',3f8.5,' setrgbcolor')
 74   format(1x,' pen',i1)
 75   format(1x,' pen',i2)
 76   format(1x,' pen',i3)
C
      RETURN
      END
C--------------------------------------------------------------bk--
      SUBROUTINE COLINT
C
C     set up colour values for colours up to 31 
C     (white background)
C
      COMMON/P00000/LPLOT,IROT,IL34,A,B,C,D,ASP,THET
      COMMON/PC0000/colr(256),colg(256),colb(256)
      COMMON/L00000/PSCA,xo,yo
C 
      dimension ICOL0(96)
      data icol0/
     &           255,255,255,  0,  0,  0,  0,  0,  0,255, 31, 31,
     &            63, 63,255, 31,255, 31,225,225, 31,225,139, 31,
     &           191,127, 63,225,191, 97,159, 63,159,225,139,139,
     &           225,139,225,139,139,225,139,225,139,159,159,159,
     &            15, 15, 15, 31, 31, 31, 47, 47, 47, 63, 63, 63, 
     &            79, 79, 79, 95, 95, 95,111,111,111,127,127,127,
     &           143,143,143,159,159,159,175,175,175,191,191,191,
     &           207,207,207,223,223,223,239,239,239,255,255,255/
      dimension ICOL1(108)
      data icol1/
     &             0,250,250, 50,250,250,100,250,250,150,250,250,
     &           200,250,250,250,250,250,250,200,250,250,150,250,
     &           250,100,250,250, 50,250,250,  0,250,255,255,255, 
     &             0,250,  0, 50,250, 50,100,250,100,150,250,150,
     &           200,250,200,250,250,250,250,250,200,250,250,150,
     &           250,250,100,250,250, 50,250,250,  0,255,255,255, 
     &             0,  0,250, 50, 50,250,100,100,250,150,150,250,
     &           200,200,250,250,250,250,250,200,200,250,150,150,
     &           250,100,100,250, 50, 50,250,  0,  0,255,255,255/ 
      dimension ICOL2(108)
      data icol2/
     &             0,100,100, 50,150,150,100,200,200,150,250,250,
     &           200,250,250,250,250,250,250,200,250,250,150,250,
     &           200,100,200,150, 50,150,100,  0,100,255,255,255, 
     &             0,100,  0, 50,150, 50,100,200,100,150,250,150,
     &           200,250,200,250,250,250,250,250,200,250,250,150,
     &           200,200,100,150,150, 50,100,100,  0,255,255,255, 
     &             0,  0,100, 50, 50,150,100,100,200,150,150,250,
     &           200,200,250,250,250,250,250,200,200,250,150,150,
     &           200,100,100,150, 50, 50,100,  0,  0,155,155,155/ 
       do 10 ic=1,32           
         ik = (ic-1)*3
         colr(ic) = float(icol0(ik+1))/255.
         colg(ic) = float(icol0(ik+2))/255.  
         colb(ic) = float(icol0(ik+3))/255.
 10    continue
c      do 11 ic=33,68           
c        ik = (ic-33)*3
c        colr(ic) = float(icol1(ik+1))/255.
c        colg(ic) = float(icol1(ik+2))/255.  
c        colb(ic) = float(icol1(ik+3))/255.
c11    continue
c      do 12 ic=69,104           
c        ik = (ic-69)*3
c        colr(ic) = float(icol2(ik+1))/255.
c        colg(ic) = float(icol2(ik+2))/255.  
c        colb(ic) = float(icol2(ik+3))/255.
c12    continue
c      do 20 ic=105,256
c        colr(ic) = 1.0
c        colg(ic) = 1.0
c        colb(ic) = 1.0
c20    continue

       do 30 ic=1,32
         if(ic-1.ge.0.and.ic-1.lt.10)then 
            write(lplot,100)ic-1,colr(ic),colg(ic),colb(ic)
         else if(ic-1.lt.100)then 
            write(lplot,101)ic-1,colr(ic),colg(ic),colb(ic)
         else if(ic-1.lt.1000)then 
            write(lplot,102)ic-1,colr(ic),colg(ic),colb(ic)
         end if
 30    continue
       do ic=33,256
         colr(ic) = -1.
         colg(ic) = -1.
         colb(ic) = -1.
       end do
 100  format(1x,' /pen',i1,' { ',3(f6.4,1x),' setrgbcolor} def')
 101  format(1x,' /pen',i2,' { ',3(f6.4,1x),' setrgbcolor} def')
 102  format(1x,' /pen',i3,' { ',3(f6.4,1x),' setrgbcolor} def')
       RETURN
       END
c------------------------------------------------------------------
      subroutine LDCOLR(lunit)

	character*256	string
c
c     load a colour map for pens > 32 from lunit
c
      COMMON/PC0000/colr(256),colg(256),colb(256)
      COMMON/P00000/LPLOT,IROT,IL34,A,B,C,D,ASP,THET
  5   read(lunit,200) string
      if(string(1:1).eq.'#')go to 5
      read(string,*) ncol
      do 10 j=1,ncol
  15    read(lunit,200) string
        if(string(1:1).eq.'#')go to 15
        read(string,*) k,ir,ig,ib
        colr(k+1) = float(ir)/255.
        colg(k+1) = float(ig)/255.
        colb(k+1) = float(ib)/255.
        if(k.ge.0.and.k.lt.10)then 
           write(lplot,100)k,colr(k+1),colg(k+1),colb(k+1)
        else if(k.lt.100)then 
           write(lplot,101)k,colr(k+1),colg(k+1),colb(k+1)
        else if(k.lt.1000)then 
           write(lplot,102)k,colr(k+1),colg(k+1),colb(k+1)
        end if
 10   continue

 100  format(1x,' /pen',i1,' { ',3(f6.4,1x),' setrgbcolor} def')
 101  format(1x,' /pen',i2,' { ',3(f6.4,1x),' setrgbcolor} def')
 102  format(1x,' /pen',i3,' { ',3(f6.4,1x),' setrgbcolor} def')
 200  format(a72)

      RETURN
      END
C-----------------------------------------------------------gh--bk
      subroutine HPLOTS(ION,IRO,LPL,ILS)
C
C    Initialises plotter : if ILS .eq. 0  -  mapped from A4 paper
C                          if ILS .eq. 1  -  mapped from A3 paper
C      and establishes handshaking characteristics (ION=1).
C      Terminates plot file if  ION .ne. 1
C
      COMMON/P00000/LPLOT,IROT,IL34,A,B,C,D,ASP,THET
      COMMON/PC0000/colr(256),colg(256),colb(256)
      COMMON/L00000/PSCA,xo,yo
c     character*40 iword
      IF(ION.eq.0)GO TO 80
C
      A=1.0
      B=0.0
      C=1.0
      D=0.0
      asp = 0.6666
      IL34 = ils
      LPLOT=LPL
      IROT=IRO
      do 31 j=1,256
        colr(j) = 0.0
        colg(j) = 0.0
        colb(j) = 0.0
 31   continue
C
c                                   Postscript initialisation
      write(LPLOT,fmt='("%!")')
      write(LPLOT,*)'gsave'
      write(LPLOT,*)'/pM {stroke newpath moveto} def'
      write(LPLOT,*)'/pL {lineto} def'
      write(LPLOT,*)'/rM {rmoveto} def'
      write(LPLOT,*)'/rL {rlineto} def'
      write(LPLOT,*)'2 setlinejoin'
      write(LPLOT,*)'/hV {/Helvetica findfont '
      write(LPLOT,*)' exch scalefont setfont} def'
      write(LPLOT,*)'/thick {1.6 setlinewidth} def '
      write(LPLOT,*)' /norm  {1.0 setlinewidth} def'
      write(LPLOT,*)' /thin  {0.4 setlinewidth} def'
      if(il34.eq.0)then
         write(LPLOT,*)' /thik1 {1.5 setlinewidth} def'
         write(LPLOT,*)' /thik2 {2.1 setlinewidth} def'
         write(LPLOT,*)' /thin1 {1.0 setlinewidth} def'
         write(LPLOT,*)' /thin2 {0.7 setlinewidth} def'
         write(LPLOT,*)' /thin3 {0.5 setlinewidth} def'
         write(LPLOT,*)' /thin4 {0.35 setlinewidth} def'
      elseif(il34.eq.1) then
         write(LPLOT,*) ' /thik1 {1.0 setlinewidth} def'
         write(LPLOT,*) ' /thik2 {2.1 setlinewidth} def'
         write(LPLOT,*) ' /thin1 {0.7 setlinewidth} def'
         write(LPLOT,*) ' /thin2 {0.5 setlinewidth} def'
         write(LPLOT,*) ' /thin3 {0.35 setlinewidth} def'
         write(LPLOT,*) ' /thin4 {0.25 setlinewidth} def'
      end if
      write(LPLOT,*) ' /bblue  {0 0 0.7 setrgbcolor} def'
      write(LPLOT,*) ' /bgreen  {0 0.4 0 setrgbcolor} def'
      write(LPLOT,*) ' /black  { 0 0 0 setrgbcolor} def'
      write(LPLOT,*) ' /blue  {0 0.6 1.0 setrgbcolor} def'
      write(LPLOT,*) ' /orange {1 0.8 0 setrgbcolor} def'
      write(LPLOT,*) ' /yellow {1 1 0 setrgbcolor} def'
      write(LPLOT,*) ' /white  {1 1 0.8 setrgbcolor} def'
      write(LPLOT,*) ' /brown  {0.75 0.5 0.25 setrgbcolor} def'
      write(LPLOT,*) ' /red    {1 0.4 0.4 setrgbcolor} def'
      write(LPLOT,*) ' /pink   {1 0.0 0.6 setrgbcolor} def'
      write(LPLOT,*) ' /green  {0.4 1 0.4 setrgbcolor} def'
      write(LPLOT,*) ' /cyan   {0 1 0.8 setrgbcolor} def'
      write(LPLOT,*) '590 5 translate '
      write(LPLOT,*) '90 rotate '
      write(LPLOT,*) 'newpath '
c
      if ( IL34 .eq. 0 )  then
        PSCA = 72.0/2.54
C
      else if (IL34 .eq. 1) then
        PSCA = 72.0/(SQRT(2.0)*2.54)
C
      else if (IL34 .eq. 2) then
C    initialise plotter to allow convenient screen mapping
C    (1cm = 10 points)
C
        PSCA = 10.0
C
      end if
C
C                                   Load colours,Select pen #1
C
      call colint
      call pen(1,0)
C
      RETURN
C
C
 80    continue

C                                   Write NA logo 
        xlogo = 26
        ylogo = 0.7
        size = 0.4
        call logo_NA(xlogo,ylogo,size)

C                                   Postscript close
       write(LPLOT,*) 'stroke'
       write(LPLOT,*) 'showpage'
       write(LPLOT,*) 'grestore'
       close(LPLOT)
       RETURN
C
      END
C------------------------------------------------------------bk--
      SUBROUTINE PICOL(xori,yori,xxl,yyl,nsxx,nsyy,arr,
     &                 nth11,nth22,sth11,sth22,imap)
C----------
C        picol         colour image plotting
C                      xori - x origin for image
C                      yori - y origin for image
C                      xxl  - x dimension of block
C                      yyl  - y dimension of block
C                      arr(nsxx,nsyy) - data array
C                      nth11 - lowest colour
C                      nth22 - highest colour
C                      sth11 - array value for bottom of lowest colour
C                      sth22 - array value for top of highest colour
C
C       The array is now organised with the origin at the bottom left
C	The x-direction is to the right an dthe y-direction is up.
C
C       For uneven spacings PIMASK can be used with 
C       the pen set to the current colour   
C
C	Colour map can be reversed by swapping pen limits around.
C--------- 
C    
      COMMON/P00000/LPLOT,IROT,IL34,A,B,C,D,ASP,THET
      COMMON/PC0000/colr(256),colg(256),colb(256)
      COMMON/L00000/PSCA,xo,yo
      real*4 xl,yl,sth1,sth2
      integer*4 nsx,nsy,nth1,nth2
      dimension arr(nsxx,nsyy)
      integer nc
C      
      xl = xxl
      yl = yyl
      nsx = nsxx
      nsy = nsyy
      nth1 = nth11
      nth2 = nth22
      nmax = max(nth11,nth22)
      nmin = min(nth11,nth22)
      sth1 = sth11
      sth2 = sth22
c     write(6,*) 'picol:',xl,yl,nsx,nsy,nth1,nth2,sth1,sth2
c 
      xso = xori*psca
c
c     yso = (yori+yyl)*psca  
c				 move origin to bottom left
      yso = yori*psca
      dsx = xl*psca/float(nsx)
      dsy = yl*psca/float(nsy)
c                                write Postscript commands                      
      write(LPLOT,*) 'gsave'
      write(LPLOT,81) xso,yso
 81   format(1x,2f8.3,' pM')
c     ysa = yso-dsy
      ysa = yso
      do 10 j=1,nsy
        xsa = xso
        do 20 i=1,nsx
          nc = int(float(nth1)+
     &         float(nth2-nth1)*(arr(i,j)-sth1)/(sth2-sth1))
c
c					map array values outside 
c					of limits to end pens
c
	  if(imap.eq.0)then
 	     if(nc.lt.nmin)nc = nmin
 	     if(nc.gt.nmax)nc = nmax
	  else
 	     if(nc.lt.nmin)go to 19
 	     if(nc.gt.nmax)go to 19
	  end if
c
          call pen(nc,0)
          write(LPLOT,81) xsa,ysa
          write(LPLOT,82) xsa,ysa+dsy,xsa+dsx,ysa+dsy,xsa+dsx,ysa
 82       format(1x,2f8.3,' pL',2f8.3,' pL',2f8.3,' pL')
          write(LPLOT,83) xsa,ysa
 83       format(1x,2f8.3,' pL closepath fill ')
 19       xsa = xsa+dsx	
 20     continue
c       ysa = ysa-dsy
        ysa = ysa+dsy
 10   continue
c
      WRITE(LPLOT,*) 'grestore'
      RETURN      
      END
c------------------------------------------------------------bk--
      subroutine pimask (xori,yori,xxl,yyl,nsxx,nsyy,arr,
     &                   nk11,nk22,nk33,sk11,sk22,kpen1)
c------
c         pimask         bitmap imagemask plotting in postscript 
c                        (1 bit representation)
c                        xori - x origin for bit image
c                        yori - y origin for bit image
c                        xxl  - x dimension of block
c                        yyl  - y dimension of block
c                        arr(nsxx,nsyy) - data array
c                        nk11  - mask on
c                        nk22  - mask on ) alternate
c                        nk33  - mask off) rows
c                        sk11,sk22  - limits for band 
c                        kpen  - pen number
c     nk1,nk2,nk3 must all lie in 1,16 and 16-(nk.-1) will
c     be interpreted as 4 1bit operations in terms of the mask
c
c     if(kpen .le. 16 )then
c       pens translated to grey scale with
c       0 - white 16 - black
c     elseif(kpen .gt. 16) then
C       pencolour k used for mask
c     endif
c    
c------
      real*4 xl,yl,sk1,sk2
      integer*4 nsx,nsy,nk1,nk2,nk3
      character*1 car,chex(16)
      dimension car(256,256)
      dimension arr(nsxx,nsyy)
      COMMON/P00000/LPLOT,IROT,IL34,A,B,C,D,ASP,THET
      COMMON/PC0000/colr(256),colg(256),colb(256)
      COMMON/L00000/PSCA,xo,yo
      COMMON/T00000/zxb0,zyb0,fchar,ffac
      character*4 fchar
      data chex/'0','1','2','3','4','5','6','7','8','9',
     &          'A','B','C','D','E','F'/
c    
      do 100 j=1,256
        do 110 k=1,256
          car(k,j) = ' '
 110    continue
 100  continue
c                                     create Hex map
      xl=xxl
      yl=yyl
      nsx=nsxx
      nsy=nsyy
      nk1=nk11
      nk2=nk22
      nk3=nk33
      sk1=sk11
      sk2=sk22
      kpen = kpen1
      write(6,*) 'pimask:',xl,yl,nsx,nsy,nk1,nk2,nk3,sk1,sk2,kpen 
c
      nsx1 = nsx
      nr = MOD(nsx,2)
      if(nr.eq.1) nsx1=nsx1+1
      nsx2 = nsx1/2
      do 10 j=1,nsy
        do 20 i=1,nsx1
          if(mod(j,2).eq.0) then
            npix = nk2
          elseif(mod(j,2).eq.1) then
            npix = nk3
          endif
          if(arr(i,j).lt.sk1) npix = nk1                  
          if(arr(i,j).gt.sk2) npix = nk1
          car(i,j) = chex(17-npix)
c        write(6,*) i,j, arr(i,j),npix,car(i,j)
20      continue
10    continue     
c100   format(Z1)
c                                     write Postscript commands 
      if(kpen.le.16) then
        gray = float(kpen)/16.
        write(8,80) gray    
 80     format(1x,f8.4,' setgray')
      elseif(kpen.gt.16) then
c       write(6,*) kpen
        call pen(kpen,1)
      endif
      nsx4 = nsx*4
      xso = xori*psca
      yso = yori*psca
      xscl = xl*psca
      yscl = yl*psca
      write(8,*) 'gsave'
      write(8,81) xso,yso
 81    format(1x,2f8.3,' translate')
      write(8,82) xscl,yscl
 82    format(1x,2f8.3,' scale')
      write(8,83) nsx4,nsy,nsx4,nsy
 83    format(2i4,' true  [',i4,' 0 0',i4,' 0 0 ]')
      write(8,*) '{<'
      write(8,85) ((car(ip,jp),ip=1,nsx1),jp=1,nsy)
 85    format(1x,64a1)
      write(8,*) '>}  imagemask'
      write(8,*) 'grestore'
c
      return 
      end
c-------------ms-------------------------------------
	subroutine xname(string,n)
c
c	This is a dummy routine for postscript library
c       It is only meaningful when called using Xpak
c
	character*256	string
	return
	end
C-------------------------------------------------------------ms-
      SUBROUTINE logo_NA(x,y,size)
C
C     Draws NA-plot logo in postscript file
C
c-------------------------------------------------------------------

      real*4  x,y
      real*4  xb(4),yb(4)

      call zpick(8,0,0)
      xsize = 3.4*size
      ysize = 1.1*size
      xb(1) = x
      xb(2) = x+xsize
      xb(3) = x+xsize
      xb(4) = x
      yb(1) = y
      yb(2) = y
      yb(3) = y+ysize
      yb(4) = y+ysize

c     call plot_box(xb,yb,1)
      call pen(1,0)
      call fillpoly(xb,yb,4)
      call pen(0,1)
      call typstr
     &     (x+0.23,y+0.13,size*0.8,'S-plot',0.0,6)
c     call typstr
c    &     (x+0.23,y+0.13,size*0.8,'NA-plot',0.0,7)

      return
      end
c-------------ms-------------------------------------
	subroutine plottype(n)
c
c	This is a routine to tell the user whether
c	the postscript library has been compiled. 
c
	n = 0
	return
	end 
c-------------ms-------------------------------------
