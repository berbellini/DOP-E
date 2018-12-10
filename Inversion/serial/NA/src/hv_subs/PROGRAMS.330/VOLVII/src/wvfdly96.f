        program wvfdly96
c-----
c       The purpose of this program is to analyze the time shifts
c       arising from the waveform matching of wvfgrd96.  The time shift 
c       is placed in the USER9 header field of the predicted trace, which 
c       has the form SLMZ.pre
c
c       The time shifts for waveform matching arise for several reasons:
c       o The origin time and epicentral distance are incorrect
c       o The velocity model used for hte inversion is incorrect
c       o The velocity model used to define the P-arrival time is not the
c         same as the velocity model used for the waveform inversion
c         (assuming that the initial trace alignment is based on the
c          P arrival time)
c        
c       This program considers the azimuthal variation of the Z R and T 
c       components, which are contained in the files Z.dat, R.dat and T.dat
c       The files consist of one line per station with an entry
c         Azimuth Time_Shift
c
c       The convention for the time shift is that a positive value
c       indicates that the synthetic traces must be later for alignment,
c       due to a later origin time than assumed, a slower velocity
c       model than assumed, or that the actual epicentral distance is
c       greater than assumed
c
c       These are fit to a functional form
c         Time_shift = A + B cos Az + C Sin Az
c
c       For small changes in epicentral location, 
c         A is the change in origin time
c       The second two terms can be rewritten, for interpretation, as
c          R/V cos ( Az - Theta)
c       where R is the shift of the epicenter in km and Theta + 180 is the
c       direction of the epicenter shift
c
c       Consider the following example with time shifts indicated by -, 0 +
c       where I is the initial epicenter and T is the true epicenter
c               
c                       0          delay           * *
c                       |            |           *    *
c                       |            |          *       *
c               + ------I--T--- -    **---|----|----|----| Az
c                       |            | * 90  *180  270  360
c                       |            |  *   *
c                       |            |    *
c                       0            |
c
c       This pattern would indicate a functional fit of
c          time_shift = A + C sin Az  with C negative
c       The remove the time shift the epicenter should be moved
c       to the east, e.g., 270 + 180 degrees
c
c       Among the output of the program are the entries
c          Origin time shift:   1.2931364    
c          px               :  0.32548723    
c          py               :  0.21361648    
c          ang              :   213.27681    
c          Rayl delay (s)   :  0.38932496    
c          Rayl shift (km)  :   1.2458400    
c          Average delay    :   1.4797298    
c          T0SHIFT=   1.2931364    
c          RSHIFT=   1.2069074    
c          AZSHIFT=   213.27679 
c
c          T0SHIFT is the change in origin time, to be ADDED to the
c                  assumed origin time
c          The epicenter is to be shifted RSHIFT (km) in the 
c          AZSHIFT (degree) direction.
c
c       To estimate the RSHIFT from the B and C coefficients, 
c       it is assumed that the Rayleigh wave group velocity (vrayl) is
c       3.1 km/s. The Love wave group velocity is assumed to be greater
c       by a factor of 1/0.92.
c
c       The transverse component time delays, e.g., Love, are converted to
c       pseudo Rayleigh wave (Z and R component) time delays by using the
c       model
c        
c         Time_shift = A + 0.92 B cos Az + 0.92 C Sin Az
c-----
c       CHANGES
c       16 SEP 2010 created
c-----
        real a(3,3), x(3), y(3)
        integer i,j

        real azobs(1000), delayobs(1000), vrtovl(1000)
c-----
c       machine dependent initialization
c-----
        call mchdep()
c-----
c       define the Rayleigh wae group velocity
c-----
        vrayl=3.1
c-----
c       initialize
c-----
        do i=1,3
              y(i) = 0.0
              x(i) =0.0
           do j=1,3
              a(i,j) = 0.0
           enddo
        enddo
c-----
c       read in Z.dat if it exists
c-----
        n = 0
        sumd = 0.0
        call fill('Z.dat',1.0,a,y,n,sumd,azobs,delayobs,vrtovl)
        call fill('R.dat',1.0,a,y,n,sumd,azobs,delayobs,vrtovl)
        call fill('T.dat',0.92,a,y,n,sumd,azobs,delayobs,vrtovl)
c-----
c       put in some stability to the inverse
c-----
        do i=1,3
            a(i,i) = a (i,i) + 0.0001
        enddo
        do i=1,3
          write(6,'(i5, 3f20.10,5x,f20.10)')i,(a(i,j),j=1,3),y(i)
        enddo
c-----
c       invert the matrix
c-----
        call matinv(a,3)
        do i=1,3
          write(6,'(i5, 3f20.10)')i,(a(i,j),j=1,3)
        enddo
c-----
c       get the solution vector
c-----
        do i=1,3
           x(i) = 0.0
           do j=1,3
              x(i) = x(i) + a(i,j)*y(j)
           enddo
        enddo
        ang = atan2(-x(3),-x(2))*180./3.1415927
        if(ang.lt.0.0)ang = ang + 360.
c-----
c       interpreted the results
c-----
        write(6,*)'Origin time shift:',x(1)
        write(6,*)'px               :',x(2)
        write(6,*)'py               :',x(3)
        write(6,*)'ang              :',ang
        write(6,*)'Rayl delay (s)   :',sqrt(x(2)*x(2) + x(3)*x(3))
        write(6,*)'Rayl shift (km)  :',sqrt(x(2)*x(2) + x(3)*x(3))*vrayl
        write(6,*)'Average delay    :',sumd/float(n)
        write(6,*)'T0SHIFT=',x(1)
        write(6,*)'RSHIFT=',sqrt(x(2)*x(2) + x(3)*x(3))*vrayl
        write(6,*)'AZSHIFT=',amod(180.0+atan2(x(3),x(2))*180./3.1415927
     1      ,360.)

c------
c       output the residual
c-----
        rms = 0.0
        do i=1,n
            saz= sin(3.1415926*azobs(i)/180.0)*vrtovl(i)
            caz= cos(3.1415926*azobs(i)/180.0)*vrtovl(i)
            
            tpre = x(1) + x(2)*caz + x(3)*saz
            write(6,'(i5,3f10.1)'),i,azobs(i),delayobs(i),
     1        delayobs(i)-tpre
            rms = rms + (delayobs(i)-tpre)**2
        enddo
        write(6,*)'RMS              :',sqrt(rms/n)
        end

        subroutine fill(fname,vrvl,a,y,n,sumd,azobs,delayobs,vrtovl)
        character fname*(*)
        real vrvl
        real a(3,3), y(3)
        real azobs(*), delayobs(*), vrtovl(*)
        integer n
        logical ext
        inquire(file=fname,exist=ext)
        if(ext)then
            open(1,file=fname,access='sequential',form='formatted',
     1        status='old')
            rewind 1
 1001       continue
            read(1,*,end=1002)az,delay
            n = n + 1
            azobs(n) = az
            delayobs(n) = delay
            vrtovl(n) = vrvl
            saz= sin(3.1415926*az/180.0)
            caz= cos(3.1415926*az/180.0)
c-----
c           correction for Love which has smaller slowness than Rayleigh
c-----
            saz = saz * vrvl
            caz = caz * vrvl
            a(1,1) = a(1,1) + 1.0*1.0
            a(1,2) = a(1,2) + caz*1.0
            a(1,3) = a(1,3) + saz*1.0
            a(2,1) = a(2,1) + 1.0*caz
            a(2,2) = a(2,2) + caz*caz
            a(2,3) = a(2,3) + saz*caz
            a(3,1) = a(3,1) + 1.0*saz
            a(3,2) = a(3,2) + caz*saz
            a(3,3) = a(3,3) + saz*saz
            y(1) = y(1) + delay * 1.0
            y(2) = y(2) + delay * caz
            y(3) = y(3) + delay * saz
            sumd = sumd + delay
            go to 1001
 1002       continue
            close(1)
        endif
        return
        end


        subroutine matinv(a,n) 
        real*4 a(3,3),v(3) 
c-----
c    this algorithm determines the inverse of a symmetric n x n matrix 
c    the inverse is determined and returned in the same location as a 
c-----
        nm1 = n - 1 
        do 35 k = 1,n 
            pivot = 1.0/a(1,1) 
            do 20 i = 2,n 
                v(i-1) = a(1,i) 
   20       continue
            do 30 i = 1,nm1 
                yy = -v(i) * pivot 
                a(i,n) = yy 
                do 31 j = 1,nm1 
                    a(i,j) = a(i+1,j+1) + v(j) * yy 
   31           continue
   30       continue
            a(n,n) = - pivot 
   35   continue
        do 40 i = 1,n 
            do 41 j = 1,n 
                a(i,j) = - a(i,j) 
   41       continue
   40   continue
        do 50 i = 2,n 
            ii = i - 1 
            do 51 j = 1,ii 
                a(i,j) = a(j,i) 
   51       continue
   50   continue
        return 
        end 
