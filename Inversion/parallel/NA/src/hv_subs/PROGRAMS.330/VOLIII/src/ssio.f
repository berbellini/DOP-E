c-----
c       common routines for dispersion files
c-----

        subroutine gtsmdl(lun,mmax,d,a,b,rho,qa,qb,nper,
     1      mname,ipar,fpar)
        parameter (NL=1000)
        real *4 d(NL), a(NL), b(NL), rho(NL), qa(NL), qb(NL)
        character mname*80
        integer ipar(10)
        real*4 fpar(10)
            read(lun)mname
            read(lun)ipar
            read(lun)fpar
            read(lun)mmax,(d(i),a(i),b(i),rho(i),
     2          i=1,mmax)
CTMP     1          qa(i),qb(i),
            read(lun)nper
        return
        end

        subroutine ptsmdl(lun,mmax,d,a,b,rho,qa,qb,nper,
     1      mname,ipar,fpar)
        parameter (NL=1000)
        real *4 d(NL), a(NL), b(NL), rho(NL), qa(NL), qb(NL)
        character mname*80
        integer ipar(10)
        real*4 fpar(10)
            write(lun,*)mname
            write(lun,*)ipar
            write(lun,*)fpar
            write(lun,*)mmax,(d(i),a(i),b(i),rho(i),
     1          i=1,mmax)
CTMP     1          qa(i),qb(i),
            write(lun,*)nper
        return
        end

        subroutine gtshed(lun,ifunc,kmode,t0,ierr)
            real*8 t0
            ierr = 0
            read(lun,end=200,err=1001) ifunc,kmode,t0
            return
  200       continue
                ierr = 200
                return
 1001       continue
                ierr = 1001
        return
        end

        subroutine ptshed(lun,ifunc,kmode,t0)
            real*8 t0
            write(lun,*) ifunc,kmode,t0
        return
        end

        subroutine gtsval(lun,cp,kmode,ierr)
        parameter(MAXMOD=2000)
        real*8 cp(MAXMOD)
            read(lun,end=2001,err=1001) (cp(i),i=1,kmode)
            return
 2001       continue
                ierr = 200
                return
 1001       continue
                ierr = 1001
            return
        end

        subroutine ptsval(lun,cp,kmode)
        parameter(MAXMOD=2000)
        real*8 cp(MAXMOD)
            write(lun,*) (cp(i),i=1,kmode)
            return
        end


c-----
c       common routines for eigenfunction files
c-----

        subroutine getmdl(lun,mmax,d,a,b,rho,qa,qb,nper,dphs,dphr,
     1      mname,ipar,fpar)
        parameter (NL=1000)
        real *4 d(NL), a(NL), b(NL), rho(NL), qa(NL), qb(NL)
        character mname*80
        integer ipar(10)
        real*4 fpar(10)
            read(lun)mname
            read(lun)ipar
            read(lun)fpar
            read(lun)mmax,(d(i),a(i),b(i),rho(i),
     1          qa(i),qb(i),i=1,mmax)
            read(lun)nper,dphs,dphr
        return
        end

        subroutine putmdl(lun,mmax,d,a,b,rho,qa,qb,nper,dphs,dphr,
     1      mname,ipar,fpar)
        parameter (NL=1000)
        real *4 d(NL), a(NL), b(NL), rho(NL), qa(NL), qb(NL)
        character mname*80
        integer ipar(10)
        real*4 fpar(10)
            write(lun,*)mname
            write(lun,*)ipar
            write(lun,*)fpar
            write(lun,*)mmax,(d(i),a(i),b(i),rho(i),
     1          qa(i),qb(i),i=1,mmax)
            write(lun,*)nper,dphs,dphr
        return
        end

        subroutine gethed(lun,ifunc,kmode,t0,ierr)
            ierr = 0
            read(lun,end=200,err=1001) ifunc,kmode,t0
            return
  200       continue
                ierr = 200
                return
 1001       continue
                ierr = 1001
        return
        end

        subroutine puthed(lun,ifunc,kmode,t0)
            real*4 t0
            write(lun,*) ifunc,kmode,t0
        return
        end

        subroutine getegn(lun,lorr,intext,wvno,u,gamma,
     1      sur,sdur,suz,sduz,sare,wvnsrc,sur0,
     2      rur,rtr,ruz,rtz,rare,wvnrec,rur0,
     3      sumkr,sumgr,sumgv,ierr)
c-----
c       lun I*4 - Logical unit number
c       lorr    I*4 - 1 = Love, 2 = Rayleigh
c       intext  I*4 - 1 write (lr)eigen85 files
c               - 0 write internal files for slat2d96
c       wvno    R*4 - wavenumber for path
c       u   R*4 - group velocity for path
c       gamma   R*4 - anelastic atenuation coefficient for path
c
c       sur R*4 - UR for Rayleigh or UT for Love at source depth
c       sdur    R*4 - d UR/dz for Rayleigh or d UT/dz 
c                   for Love at source depth
c       suz R*4 - UZ for Rayleigh at source depth
c       sduz    R*4 - d UZ/dz for Rayleigh at source depth
c       sare    R*4 - energy integral at source
c       wvnsrc  R*4 - wavenumber at source
c       sur0    R*4 - free surface ellipticity ur/uz at source
c
c       wvnrec  R*4 - wavenumber at receiver
c       rur R*4 - UR for Rayleigh or UT for Love at receiver depth
c       rtr R*4 - d UR/dz for Rayleigh or d UT/dz 
c                   for Love at receiver depth
c       ruz R*4 - UZ for Rayleigh at receiver depth
c       rtz R*4 - d UZ/dz for Rayleigh at receiver depth
c       rare    R*4 - energy integral at receiver
c       rur0    R*4 - ellipticity at receiver
c
c       sumkr   R*4 - sum of kr for slat2d96
c       sumgv   R*4 - sum of r/u for slat2d96
c       sumgr   R*4 - sum of gamma*r for slat2d96
c
c       ierr    I*4 - 1001 EOR or ERR on read
c-----
        integer LER
        parameter (LER=0)
        ierr = 0
        if(intext .ne. 0)then
            if(lorr.eq.1) then
                read(lun,end=1001,err=1001)wvno,u,gamma,
     1              sur,sdur,
     1              rur,rtr,
     1              sare
            else if(lorr.eq.3) then
                read(lun,end=1001,err=1001)wvno,u,gamma,
     1              sur,sdur,
     1              rur,rtr,
     1              wvnsrc,wvnrec,sare,rare
            else if(lorr.eq.2) then
                read(lun,end=1001,err=1001)wvno,u,gamma,
     1              sur,sdur,suz,sduz,
     1              rur,rtr,ruz,rtz,
     2              sur0, rur0,sare
            else if(lorr.eq.4) then
                read(lun,end=1001,err=1001)wvno,u,gamma,
     1              sur,sdur,suz,sduz,
     1              rur,rtr,ruz,rtz,
     1              wvnsrc,wvnrec,sur0,rur0,sare,rare
            endif
        else
            if(lorr.eq.1) then
                read(lun,end=1001,err=1001)wvno,are,u,gamma,
     1              sur,sdur,
     1              rur,rtr,
     1              sare,rare,
     1              wvnsrc,wvnrec,
     2              sumkr, sumgr,sumgv
            else if(lorr.eq.2) then
                read(lun,end=1001,err=1001)wvno,are,u,gamma,
     1              sur,sdur,suz,sduz,
     1              rur,rtr,ruz,rtz,
     2              sur0, rur0,sare,rare,
     2              wvnsrc,wvnrec,
     2              sumkr, sumgr,sumgv
            endif
        endif
        return
 1001   continue
            ierr = 1001
            return
        end

        subroutine putegn(lun,lorr,intext,wvno,u,gamma,
     1      sur,sdur,suz,sduz,sare,wvnsrc,sur0,
     2      rur,rtr,ruz,rtz,rare,wvnrec,rur0,
     2      sumkr,sumgr,sumgv)
c-----
c       lun I*4 - Logical unit number
c       lorr    I*4 - 1 = Love, 2 = Rayleigh
c       intext  I*4 - 1 write (lr)eigen85 files
c               - 0 write internal files for slat2d96
c       wvno    R*4 - wavenumber for path
c       u   R*4 - group velocity for path
c       gamma   R*4 - anelastic atenuation coefficient for path
c
c       sur R*4 - UR for Rayleigh or UT for Love at source depth
c       sdur    R*4 - d UR/dz for Rayleigh or d UT/dz 
c                   for Love at source depth
c       suz R*4 - UZ for Rayleigh at source depth
c       sduz    R*4 - d UZ/dz for Rayleigh at source depth
c       sare    R*4 - energy integral at source
c       wvnsrc  R*4 - wavenumber at source
c       sur0    R*4 - free surface ellipticity ur/uz at source
c
c       wvnrec  R*4 - wavenumber at receiver
c       rur R*4 - UR for Rayleigh or UT for Love at receiver depth
c       rtr R*4 - d UR/dz for Rayleigh or d UT/dz 
c                   for Love at receiver depth
c       ruz R*4 - UZ for Rayleigh at receiver depth
c       rtz R*4 - d UZ/dz for Rayleigh at receiver depth
c       rare    R*4 - energy integral at receiver
c       rur0    R*4 - ellipticity at receiver
c
c       sumkr   R*4 - sum of kr for slat2d96
c       sumgv   R*4 - sum of r/u for slat2d96
c       sumgr   R*4 - sum of gamma*r for slat2d96
c-----
        if(intext .ne. 0)then
            if(lorr.eq.1) then
        write(LUN,*)'wvno,u,gamma/sut,sdut,rut,rdut/sale'
                write(lun,*)wvno,u,gamma
                write(lun,*)sur,sdur,rur,rtr
                write(lun,*)sare
            else if(lorr.eq.3) then
                write(lun,*)wvno,u,gamma
                write(lun,*)sur,sdur,rur,rtr
                write(lun,*)wvnsrc,wvnrec,sare,rare
            else if(lorr.eq.2) then
        write(LUN,*)'wvno,u,gamma/sur,sdur,suz,sduz/',
     1  'rur,rtr,ruz,rtz,sur0,rur0,sare'
                write(lun,*)wvno,u,gamma
                write(lun,*)sur,sdur,suz,sduz
                write(lun,*)rur,rtr,ruz,rtz
                write(lun,*)sur0, rur0,sare
            else if(lorr.eq.4) then
                write(lun,*)wvno,u,gamma
                write(lun,*)sur,sdur,suz,sduz
                write(lun,*)rur,rtr,ruz,rtz
                write(lun,*)sur0, rur0,sare,rare
                write(lun,*)wvnsrc,wvnrec
            endif
        else
            if(lorr.eq.1) then
                write(lun,*)wvno,are,u,gamma,ur,dur
CTMP     1          ,arer,wvnor,
CTMP     1              sumkr, sumgr,sumgv,
CTMP     2              wvnsrc,wvnrec
            else if(lorr.eq.2) then
                write(lun,*)wvno,are,u,gamma,ur,dur,uz,duz
CTMP     1          ,arer,wvnor,sumkr, sumgr,sumgv, ur0r,
CTMP     2          wvnsrc,wvnrec.sur0,rur0
            endif
        endif
        return
        end

