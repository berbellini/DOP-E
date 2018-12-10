      subroutine cfopen(lunit,ncol,id,ir,ig,ib)
c 
      integer*4 id(255),ir(255),ig(255),ib(255),ncol
      integer*4 lunit
	character*256	string
c
  5   read(lunit,fmt='(a72)')string
      if(string(1:1).eq.'#')go to 5
      read(string,*)ncol
      do 10 i=1,ncol
  15    read(lunit,fmt='(a72)')string
        if(string(1:1).eq.'#')go to 15
        read(string,*)id(i),ir(i),ig(i),ib(i)
  10  continue 
c
      return
      end     
c
      subroutine logo_NA(x,y,s)
c 
c				A dummy routine for xpak
c
      return
      end     
