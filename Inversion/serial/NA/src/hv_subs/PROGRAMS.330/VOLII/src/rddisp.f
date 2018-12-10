        subroutine rddisp(fildsp,lun,ierr,ndsp,NOBS,
     1      jlorr,jobs,jobsyn,jmode,
     2      fper,fobs,fobserr)
c-----
c       CHANGES
c       19 JAN 2002 - Make obsyn=2 for synthetic, 1 for observed
c                     to be compatible with f96subs.f 2=syn/1=obs
c       05 MAR 2002 - synthetic denoted by 'T' and observed by 'X'
c-----
c       read a SURF96 dispersion file
c-----
c       fildsp  Ch* name of dispersion file
c       lun I*4 logical unit for read
c       ierr    I*4 Error code: -1 file does not exist
c       ndsp    I*4 number of dispersion values
c
c       ftype   A   should be SURF96
c       jlorr() I*4 1 = Love, 2 = Rayleigh
c       jobs()  I*4 1 = phase velocity, 2=group velocity, 3=gamma
c       jobsyn() I*4    1 = observed, 2 synthetic
c       jmode() I*4 mode: 0=Fund, 1=1 st, etc
c       fper()  R*4 period
c       fobs()  R*4 Array of ground velocities
c       fobserr()   R*4 Error in u()
c       
c
c
CSURF96 R U T 0  350.000     2.694     0.635  
c-----
c-----
c       command line variables
c-----
        character fildsp*(*)
        integer lun, ierr, ndsp,  NOBS
        integer*4 jlorr(NOBS), jobs(NOBS), jobsyn(NOBS), jmode(NOBS)  
        real*4  fper(NOBS)   , fobs(NOBS)   , fobserr(NOBS)
c-----
c       internal variables
c-----
        logical ext
        character instr*132
        character ic*1
        integer ilorr, imode, iobs
        ierr = -1
        ext = .false.
        inquire(file=fildsp,exist=ext)
        if(.not. ext)then
            ierr = -1
        else
            open(lun,file=fildsp,access='sequential',
     1      form='formatted',status='old')
            rewind lun
c-----
c           read values and parse them correctly
c-----
            ndsp = 0
 1000       continue
            read(lun,'(a)',end=9999)instr
            ls = lgstr(instr)
c-----
c           do the parsing
c-----
            if(instr(1:6).eq.'surf96' .or.
     1          instr(1:6).eq.'SURF96')then
c-----
c               now get wave type
c-----      
                lsep = 6
                call getblnk(instr,lsep,ls,lnobl)
                ic = instr(lnobl:lnobl)
                if(ic(1:1).eq.'R'. or. ic(1:1).eq.'r')then
                    ilorr = 2
                else if(ic(1:1).eq.'L'. or. ic(1:1).eq.'l')then
                    ilorr = 1
                else if(ic(1:1).eq.'A'. or. ic(1:1).eq.'a')then
                    ilorr = 0
                else
                    go to 1000
                endif
c-----
c               now get observation type
c-----
                lsep = lnobl
                call getblnk(instr,lsep,ls,lnobl)
                ic = instr(lnobl:lnobl)
                if(ic(1:1).eq.'C'. or. ic(1:1).eq.'c')then
                    iobs = 1
                else if(ic(1:1).eq.'U'. or. ic(1:1).eq.'u')then
                    iobs = 2
                else if(ic(1:1).eq.'G'. or. ic(1:1).eq.'g')then
                    iobs = 3
                else
                    go to 1000
                endif
c-----
c               now get whether observed or synthetic
c-----
                lsep = lnobl
                call getblnk(instr,lsep,ls,lnobl)
                ic = instr(lnobl:lnobl)
                if(ic(1:1).eq.'T'. or. ic(1:1).eq.'t')then
                    iobsyn = 2
c-----
c       the F designation is deprecated
c-----
                else if(ic(1:1).eq.'F'. or. ic(1:1).eq.'f')then
                    iobsyn = 1
                else if(ic(1:1).eq.'X'. or. ic(1:1).eq.'x')then
                    iobsyn = 1
                else
                    go to 1000
                endif
c-----
c-----
c               now get the values using list directed IO
c-----
                lsep = lnobl
                call getblnk(instr,lsep,ls,lnobl)
                read(instr(lnobl:ls),*)
     1          imode,per,obs,obserr
c-----
c       test for array dimensions
c-----
                if(ndsp.eq.NOBS)go to 9999
                ndsp = ndsp + 1
                jlorr(ndsp) = ilorr
                jobs(ndsp)  = iobs
                jobsyn(ndsp)  = iobsyn
                jmode(ndsp) = imode
                fper(ndsp)  = per
                fobs(ndsp)  = obs
                fobserr(ndsp)  = obserr
            endif
            
            go to 1000
        endif
 9999       continue
            close (lun)
            ierr = 0
        return
        end


        subroutine getblnk(instr,lsep,ls,lnobl)
c-----
c       determine first non-blank character
c
c       instr   Ch* Character string to be parsed
c       lsep    I*4 index of last non blank character
c       ls  I*4 length of input string
c       lnobl   I*4 index of first non blank character
c-----
        character instr*(*)
        integer lsep,ls,lnobl
        character tab*1
        tab=char(9)
        lnobl = lsep+1
        igotit = 0
        do 1000 i=lsep+1,ls
            if(igotit.eq.0)then
            if(instr(i:i).ne.' ' .and. instr(i:i).ne.tab)then
                lnobl = i
                igotit = 1
            endif
            endif
 1000   continue
        return
        end
