        subroutine shwrfw()
c-----
c       display receiver function weights
c
c       For simplicity just read it from the file
c       instead of passing all arguments through the
c       command line
c-----
        implicit none
        integer LER, LIN, LOT
        parameter(LER=0,LIN=5,LOT=6)

        character fname*80
        integer iprog, itot, nf1, nf2, nf34, nf5, nf67, nfilt
        integer nup, ivcs, lsin
        real dlam, qaqb, wref
        real twmn,twmx
        integer iter,nurftn,indp
        real pval, sigv, sigr, sigg
        integer id2 ,idum2,idum3,idum4,idum5
        real rdum1,rdum2,rdum3,rdum4,rdum5

        character kstnm*8
        integer nzyear, nzjday, nzhour, nzmin, nzmon, nzday
        integer n
        real dt, rayp, gaussalp, delay, invwgt

        integer i, ls
        integer lgstr

        integer NL
        parameter (NL=200)
        integer nf10(NL)
        data nf10/NL*1/

        call gttmp0(iprog,itot,nf1,nf2,nf34,nf5,nf67,nf10,nfilt,
     1      nup,dlam,qaqb,wref,ivcs,
     2      lsin,twmn,twmx,iter,nurftn,indp,pval,
     4              sigv,sigr,sigg,
     3      id2 ,idum2,idum3,idum4,idum5,
     4      rdum1,rdum2,rdum3,rdum4,rdum5)
c----
c       open the file with the receiver function information
c-----
        open(3,file='tmpsrfi.15',access='sequential',form='formatted',
     1      status='unknown')
        rewind 3
        write(LOT,1)
        write(LOT,2)
        do 1000 i=1,nurftn
            read(3,'(a)',end=9000)fname
            read(3,'(a)')kstnm
            read(3,*)nzyear, nzjday, nzhour, nzmin, nzmon, nzday
            read(3,*)n, dt, rayp, gaussalp, delay, invwgt
            ls = lgstr(fname)
        write(LOT,3)i, kstnm, gaussalp, invwgt, fname(1:ls)

 1000   continue
 9000   continue
        write(LOT,2)
        close(3)
    1   format('NUM ','Station ','     Gauss','   Rftn Wt',
     1      ' File Name')
    2   format('----------------------------------------'   ,
     1  '----------------------------------------'   )
    3   format(i3,1x,a8,f10.2,f10.2,1x,a)
        WRITE(LOT,*)'Use option 50 to change RFTN inversion weight'
        WRITE(LOT,*)'Use option 49 to redisplay this menu'

        return
        end

        subroutine chgrfw(numa,i)
c-----
c       change the receiver function weight
c-----
        implicit none
        integer LER, LIN, LOT
        parameter(LER=0,LIN=5,LOT=6)

c-----
c       numa    number of arguments on command line
c       i   current command line argument
c-----
        integer numa, i
c-----
c       name    command line argument
c       j   internal counter
c       lgstr   definition of function type
c       ls  length of string 
c       ival    index to receiver function to change
c       fval    new value of weight
c       murftn  actual number of lists entries read
c-----
        character*20 name
        integer j
        integer lgstr
        integer ls
        integer ival
        real fval
        integer murftn

c-----
c       variables for tmpsrfi.00
c-----
        integer iprog,itot,nf1,nf2,nf34,nf5,nf67,nup,nfilt
        real dlam,qaqb,wref
        integer ivcs,lsin,iter,nurftn,indp
        real twmn,twmx,pval
        real sigv,sigr,sigg
        integer id2 ,idum2,idum3,idum4,idum5
        real rdum1,rdum2,rdum3,rdum4,rdum5
        integer NL
        parameter (NL=200)
        integer nf10(NL)
        data nf10/NL*1/
c-----
c       temporary parameters for storing information about 
c       receiver function inversion
c       from tmpsrfi.15
c-----
        integer MXRF
        parameter (MXRF=1000)
        character*80 fname(MXRF)
        character*8 kstnm(MXRF)
        integer nzyear(MXRF), nzjday(MXRF), nzhour(MXRF)
        integer  nzmin(MXRF), nzmon(MXRF), nzday(MXRF)
        integer n(MXRF)
        real dt(MXRF), rayp(MXRF), gaussalp(MXRF), 
     1                  delay(MXRF), invwgt(MXRF)

c-----
c       get the number of receiver functions
c-----
        call gttmp0(iprog,itot,nf1,nf2,nf34,nf5,nf67,nf10,nfilt,
     1      nup,dlam,qaqb,wref,ivcs,
     2      lsin,twmn,twmx,iter,nurftn,indp,pval,
     3              sigv,sigr,sigg,
     4      id2 ,idum2,idum3,idum4,idum5,
     5      rdum1,rdum2,rdum3,rdum4,rdum5)
c----
c       open the file with the receiver function information
c-----
        open(3,file='tmpsrfi.15',access='sequential',form='formatted',
     1      status='unknown')
        rewind 3
        murftn = 0
        do 1000 j=1,nurftn
            read(3,'(a)',end=9000)fname(j)
            read(3,'(a)')kstnm(j)
            read(3,*)nzyear(j), nzjday(j), nzhour(j), 
     1          nzmin(j), nzmon(j), nzday(j)
            read(3,*)n(j), dt(j), rayp(j), gaussalp(j), 
     1                  delay(j), invwgt(j)
            murftn = murftn + 1
 1000   continue
 9000   continue
        close(3)
c-----
c       now process the command line or get the interactive input
c-----
        if(numa.ne.0)then
            i=i+1
            name = ' '
            call mgtarg(i,name)
            read(name,'(i10)')ival
            i=i+1
            name=' '
            call mgtarg(i,name)
            read(name,'(f10.0)')fval
            if(fval.lt.0.0001)fval=0.0001
            if(ival.lt.1.or.ival.gt.nurftn)then
                return
            else
                invwgt(ival) = fval
            endif
        else
            write(LOT,*)'Enter i'
            read(LIN,*)ival
            if(ival.lt.1.or.ival.gt.nurftn)then
                write(LOT,*)'ival out of range'
                return
            else
                ls = lgstr(fname(ival))
                write(LOT,*)'Weight for ',
     1              fname(ival)(1:ls),' is ',invwgt(ival)
                write(LOT,*)'Enter new weight'
                read(LIN,*)fval
                if(fval.lt.0.0001)fval=0.0001
            endif
            invwgt(ival) = fval
        endif
c-----
c       now update the tmpsrfi.15 file with the new weight
c-----
        open(3,file='tmpsrfi.15',access='sequential',form='formatted',
     1      status='unknown')
        rewind 3
        do 2000 j=1,murftn
            ls = lgstr(fname(j))
            write(3,'(a)')fname(j)(1:ls)
            write(3,'(a)')kstnm(j)
            write(3,'(i4,1x,i3,1x,i2,1x,i2,1x,i2,1x,i2)')
     1          nzyear(j), nzjday(j), nzhour(j), 
     2          nzmin(j), nzmon(j), nzday(j)
            write(3,'(i5,f10.3,f10.3,f10.3,f10.3,f10.3)')
     1          n(j), dt(j), rayp(j), gaussalp(j), 
     2          delay(j), invwgt(j)
 2000   continue
        close(3)
        return
        end
