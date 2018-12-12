*
*     write a "bar chart" in ascii
*
*     n = number of values
*     labels -> labels for each value
*     values -> values to chart
*     nl = number of characters in label (to print)
*     nx = number of x's in maximum value
*
      subroutine bar_chart(n,labels,values,nl,nx,ounit)
      integer n, nl, nx,ounit
      character*256 labels(n)
      real values(n)
      
       write(ounit,*) ' '
       write(ounit,*) 'Misfit Summary:'
       write(ounit,*) ' '
       write(ounit,'(a18,1x,a1,49x,a1)')'Fractional Misfit:','0','1'
       write(ounit,100) '|','|','|','|','|'
100    format(19x,a1,9x,a1,9x,a1,9x,a1,9x,a1,9x,a1)      
       do 1 i = 1, n
         call do_bc_line(labels(i),values(i),12,50,ounit)
1      continue
       write(ounit,100) '|','|','|','|','|','|'
       write(ounit,*) ' '
       
       return
       end
*
*   label = string to plot on left
*   misfit = fractional length of x string
*   nl = length of string to print
*   nx = maximum number of x's (for a value of 1.0)
*

      subroutine do_bc_line(label,misfit,nl,nx,ounit)
        character*256 label 
        real      misfit 
        integer   nl, nx, ounit, i, n
        character*16 myformat
*     
*      build a format statement
*       
       if(nl .lt.10) then
           write(myformat,'(a2,i1,a11)') '(a', nl, ',1x,f5.2,$)'
       else
           write(myformat,'(a2,i2,a11)') '(a', nl, ',1x,f5.2,$)'
       endif
*
       write(ounit,myformat) label,  misfit
*
*      clips the value to 1.0
*      
       n = ifix(misfit * nx)
       if(n .gt. nx) n = nx
*      
       write(ounit,'(a2,$)')' |'
       do 100 i = 1, n
          write(ounit,'(a1,$)')'x'
100    continue
*
       do 105 i = n+1, nx-1
          write(ounit,'(a1,$)')' '
105    continue
       write(ounit,'(a1)')'|'

       return
       end
      
