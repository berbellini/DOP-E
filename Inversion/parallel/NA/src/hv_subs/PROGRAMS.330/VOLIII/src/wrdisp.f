        subroutine wrdisp(fildsp,lun,ndsp,
     1      jlorr,jobs,jobsyn,jmode,
     2      fper,fobs,fobserr,NOBS)
c-----
c       CHANGES
c       19 JAN 2002 - Make obsyn=2 for synthetic, 1 for observed
c                     to be compatible with f96subs.f 2=syn/1=obs
c       05 MAR 2002 - synthetic denoted by 'T' and observed by 'X'
c-----
c       read a SURF96 dispersion file
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
c       jobsyn()I*4 1 = observed, 2 theoretical
c       jmode() I*4 mode: 0=Fund, 1=1'st, etc
c       fper()  R*4 period
c       fobs()  R*4 Array of ground velocities
c       fobserr()   R*4 Error in u()
c       
c
c
CSURF96 R U T -1   350.000     2.694     0.635  
CSURF96 R U F  0  290.0000   4.18800   1.27200  
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
        character ostr*12
        integer ilorr, imode, iobs
        logical ext
        ext = .false.
        inquire(file=fildsp,exist=ext)
        if(.not. ext)then
            open(lun,file=fildsp,access='sequential',
     1      form='formatted',status='new')
        else
            open(lun,file=fildsp,access='sequential',
     1      form='formatted',status='old')
        endif
            rewind lun
c-----
c           write values and parse them correctly
c-----
            do 1000 i=1,ndsp
                ostr(1:6) = 'SURF96'
c-----
c               now set wave type
c-----      
                ilorr = jlorr(i) 
                iobs  = jobs(i) 
                iobsyn  = jobsyn(i) 
                mode  = jmode(i) 
                per    = fper(i)
                obs    = fobs(i)
                obserr = fobserr(i)
                if(ilorr.eq.0)then
                    ostr(7:8) = ' A'
                else if(ilorr.eq.1)then
                    ostr(7:8) = ' L'
                else if(ilorr.eq.1)then
                    ostr(7:8) = ' R'
                endif
c-----
c               now set observation type
c-----
                if(iobs.eq.1)then
                    ostr(9:10) = ' C'
                else if(iobs.eq.2)then
                    ostr(9:10) = ' U'
                else if(iobs.eq.3)then
                    ostr(9:10) = ' G'
                endif
c-----
c               now set observation type
c-----
                if(iobsyn.eq.2)then
                    ostr(11:12) = 'T'
                else if(iobsyn.eq.1)then
                    ostr(11:12) = 'X'
                endif
                write(lun,1)mode,per,obs,obserr
 1000       continue
    1   format(a12,1x,i3,1x,3g11.4)
        return
        end

