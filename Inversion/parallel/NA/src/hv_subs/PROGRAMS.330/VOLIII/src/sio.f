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
            write(lun)mname
            write(lun)ipar
            write(lun)fpar
            write(lun)mmax,(d(i),a(i),b(i),rho(i),
     1          i=1,mmax)
CTMP     1          qa(i),qb(i),
            write(lun)nper
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
            write(lun) ifunc,kmode,t0
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
            write(lun) (cp(i),i=1,kmode)
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
            write(lun)mname
            write(lun)ipar
            write(lun)fpar
            write(lun)mmax,(d(i),a(i),b(i),rho(i),
     1          qa(i),qb(i),i=1,mmax)
            write(lun)nper,dphs,dphr
C           write(6,*)mname
C           write(6,*)ipar
C           write(6,*)fpar
C           write(6,*)mmax
C           write(6,1)(d(i),a(i),b(i),rho(i),
C     1         qa(i),qb(i),i=1,mmax)
C           write(6,*)nper,dphs,dphr
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
            write(lun) ifunc,kmode,t0
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
c       sdu R*4 - d UZ/dz for Rayleigh at source depth
c       sare    R*4 - energy integral at source
c       wvnsrc  R*4 - wavenumber at source
c       sur0    R*4 - free surface ellipticity ur/uz at source
c
c       wvnrec  R*4 - wavenumber at receiver
c       rur R*4 - UR for Rayleigh or UT for Love at receiver depth
c       rtr R*4 - Tr for Rayleigh or Tt
c                   for Love at receiver depth
c       ruz R*4 - UZ for Rayleigh at receiver depth
c       rtz R*4 - Tz for Rayleigh at receiver depth
c       rare    R*4 - energy integral at receiver
c       rur0    R*4 - ellipticity at receiver
c
c       sumkr   R*4 - sum of kr for slat2d96
c       sumgv   R*4 - sum of r/u for slat2d96
c       sumgr   R*4 - sum of gamma*r for slat2d96
c
c       ierr    I*4 - 1001 EOR or ERR on read
c-----
c       tare    R*4 - place holder
c-----
        integer LER
        parameter (LER=0)
        ierr = 0
        if(intext .ne. 0)then
            sumkr = 0.0
            sumgr = 0.0
            sumgv = 0.0
            if(lorr.eq.1) then
                read(lun,end=1001,err=1001)wvno,u,gamma,
     1              sur,sdur,
     1              rur,rtr,
     1              sare
                wvnsrc = wvno
                wvnrec = wvno
                rare = sare
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
                wvnsrc = wvno
                wvnrec = wvno
                rare = sare
            else if(lorr.eq.4) then
                read(lun,end=1001,err=1001)wvno,u,gamma,
     1              sur,sdur,suz,sduz,
     1              rur,rtr,ruz,rtz,
     1              wvnsrc,wvnrec,sur0,rur0,sare,rare
            endif
        else
            if(lorr.eq.1) then
                read(lun,end=1001,err=1001)wvno,tare,u,gamma,
     1              sur,sdur,
     1              rur,rtr,
     1              sare,rare,
     1              wvnsrc,wvnrec,
     2              sumkr, sumgr,sumgv
            else if(lorr.eq.2) then
                read(lun,end=1001,err=1001)wvno,tare,u,gamma,
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
     3      sumkr,sumgr,sumgv)
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
c       rtr R*4 - Tr for Rayleigh or Tt
c                   for Love at receiver depth
c       ruz R*4 - UZ for Rayleigh at receiver depth
c       rtz R*4 - Tz for Rayleigh at receiver depth
c       rare    R*4 - energy integral at receiver
c       rur0    R*4 - ellipticity at receiver
c
c       sumkr   R*4 - sum of kr for slat2d96
c       sumgv   R*4 - sum of r/u for slat2d96
c       sumgr   R*4 - sum of gamma*r for slat2d96
c-----
c       tare    R*4 - place holder
c-----
        if(intext .ne. 0)then
            if(lorr.eq.1) then
                write(lun)wvno,u,gamma,
     1              sur,sdur,rur,rtr,sare
            else if(lorr.eq.3) then
                write(lun)wvno,u,gamma,
     1              sur,sdur,rur,rtr
     1              ,wvnsrc,wvnrec,sare,rare
            else if(lorr.eq.2) then
                write(lun)wvno,u,gamma,
     1              sur,sdur,suz,sduz,
     2              rur,rtr,ruz,rtz,
     3              sur0, rur0,sare
            else if(lorr.eq.4) then
                write(lun)wvno,u,gamma,
     1              sur,sdur,suz,sduz,rur,rtr,ruz,rtz
     1              ,wvnsrc,wvnrec,sur0,rur0,sare,rare
            endif
        else
            if(lorr.eq.1) then
                write(lun)wvno,tare,u,gamma,
     1              sur,sdur,
     1              rur,rtr,
     1              sare,rare,
     1              wvnsrc,wvnrec,
     2              sumkr, sumgr,sumgv
            else if(lorr.eq.2) then
                write(lun)wvno,tare,u,gamma,
     1              sur,sdur,suz,sduz,
     1              rur,rtr,ruz,rtz,
     2              sur0, rur0,sare,rare,
     2              wvnsrc,wvnrec,
     2              sumkr, sumgr,sumgv
            endif
        endif
        return
        end

c-----
c       routines for depth dependent eigenfunctions
c       note these will NEVER be used with SLAT96
c-----
        subroutine getder(lun,lorr,wvno,u,gamma,
     1      sur,sdur,sd2ur,suz,sduz,sd2uz,sare,wvnsrc,sur0,
     2      rur,rtr,ruz,rtz,rare,wvnrec,rur0,
     3      sumkr,sumgr,sumgv,ierr,mmax,dcdh,dcda,dcdb,dcdr,
     4      ur,tur,uz,tuz,ipar)
c-----
c       lun I*4 - Logical unit number
c       lorr    I*4 - 5 = Love, 6 = Rayleigh
c       wvno    R*4 - wavenumber for path
c       u   R*4 - group velocity for path
c       gamma   R*4 - anelastic atenuation coefficient for path
c
c       sur R*4 - UR for Rayleigh or UT for Love at source depth
c       sdur    R*4 - d UR/dz for Rayleigh or d UT/dz 
c                   for Love at source depth
c       suz R*4 - UZ for Rayleigh at source depth
c       sdu R*4 - d UZ/dz for Rayleigh at source depth
c       sare    R*4 - energy integral at source
c       wvnsrc  R*4 - wavenumber at source
c       sur0    R*4 - free surface ellipticity ur/uz at source
c
c       wvnrec  R*4 - wavenumber at receiver
c       rur R*4 - UR for Rayleigh or UT for Love at receiver depth
c       rtr R*4 - Tr for Rayleigh or Tt
c                   for Love at receiver depth
c       ruz R*4 - UZ for Rayleigh at receiver depth
c       rtz R*4 - Tz for Rayleigh at receiver depth
c       rare    R*4 - energy integral at receiver
c       rur0    R*4 - ellipticity at receiver
c
c       sumkr   R*4 - sum of kr for slat2d96
c       sumgv   R*4 - sum of r/u for slat2d96
c       sumgr   R*4 - sum of gamma*r for slat2d96
c
c       mmax    I*4 - number of layers in model
c
c       dcdh    R*4 - array of layer thickness partials
c       dcda    R*4 - array of layer p-velocity partials
c       dcdb    R*4 - array of layer s-velocity partials
c       dcdr    R*4 - array of layer density partials
c
c       ur  R*4 - array of radial eigenfunction
c       tur R*4 - array of radial stress eigenfunction
c       uz  R*4 - array of vertical eigenfunction
c       tuz R*4 - array of vertical stress eigenfunction
c
c       ierr    I*4 - 1001 EOR or ERR on read
c       ipar    I*4 - array of integer controls
c-----
        integer LER
        parameter (LER=0)
        parameter (NL=1000)
        real*4 dcdh(NL), dcda(NL), dcdb(NL), dcdr(NL)
        real*4 ur(NL), tur(NL), uz(NL), tuz(NL)
        integer*4 ipar(10)
        ierr = 0
c-----
c       initialize
c-----
            sumkr = 0.0
            sumgr = 0.0
            sumgv = 0.0
            do 185 i=1,mmax
                ur(i) = 0.0
                tur(i) = 0.0
                uz(i) = 0.0
                tuz(i) = 0.0
                dcdh(i) = 0.0
                dcda(i) = 0.0
                dcdb(i) = 0.0
                dcdr(i) = 0.0
  185       continue
c-----
c       read in the data stream
c-----

            if(lorr.eq.5) then
                read(lun,end=1001,err=1001)wvno,u,gamma,
     1              sur,sdur,sd2ur,rur,rtr,sare
            rare = sare
                if(ipar(4).eq.1)then
                read(lun,end=1001,err=1001)(ur(i),i=1,mmax)
                read(lun,end=1001,err=1001)(tur(i),i=1,mmax)
                endif
                if(ipar(5).eq.1)then
                read(lun,end=1001,err=1001)(dcdh(i),i=1,mmax)
                endif
                if(ipar(7).eq.1)then
                read(lun,end=1001,err=1001)(dcdb(i),i=1,mmax)
                endif
                if(ipar(8).eq.1)then
                read(lun,end=1001,err=1001)(dcdr(i),i=1,mmax)
                endif
            else if(lorr.eq.6) then
                read(lun,end=1001,err=1001)wvno,u,gamma,
     1              sur,sdur,sd2ur,suz,sduz,sd2uz,
     1              rur,rtr,ruz,rtz,
     2              sur0, rur0,sare
            rare = sare
                if(ipar(4).eq.1)then
                read(lun,end=1001,err=1001)(ur(i),i=1,mmax)
                read(lun,end=1001,err=1001)(tur(i),i=1,mmax)
                read(lun,end=1001,err=1001)(uz(i),i=1,mmax)
                read(lun,end=1001,err=1001)(tuz(i),i=1,mmax)
                endif
                if(ipar(5).eq.1)then
                read(lun,end=1001,err=1001)(dcdh(i),i=1,mmax)
                endif
                if(ipar(6).eq.1)then
                read(lun,end=1001,err=1001)(dcda(i),i=1,mmax)
                endif
                if(ipar(7).eq.1)then
                read(lun,end=1001,err=1001)(dcdb(i),i=1,mmax)
                endif
                if(ipar(8).eq.1)then
                read(lun,end=1001,err=1001)(dcdr(i),i=1,mmax)
                endif
            endif
        return
 1001   continue
            ierr = 1001
            return
        end


        subroutine putder(lun,lorr,wvno,u,gamma,
     1      sur,sdur,sd2ur,suz,sduz,sd2uz,sare,wvnsrc,sur0,
     2      rur,rtr,ruz,rtz,rare,wvnrec,rur0,
     3      sumkr,sumgr,sumgv,mmax,dcdh,dcda,dcdb,dcdr,
     4      ur,tur,uz,tuz,ipar)
c-----
c       lun I*4 - Logical unit number
c       lorr    I*4 - 5 = Love, 6 = Rayleigh
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
c       rtr R*4 - Tr for Rayleigh or Tt
c                   for Love at receiver depth
c       ruz R*4 - UZ for Rayleigh at receiver depth
c       rtz R*4 - Tz for Rayleigh at receiver depth
c       rare    R*4 - energy integral at receiver
c       rur0    R*4 - ellipticity at receiver
c
c       sumkr   R*4 - sum of kr for slat2d96
c       sumgv   R*4 - sum of r/u for slat2d96
c       sumgr   R*4 - sum of gamma*r for slat2d96
c
c       mmax    I*4 - number of layers in model
c
c       dcdh    R*4 - array of layer thickness partials
c       dcda    R*4 - array of layer p-velocity partials
c       dcdb    R*4 - array of layer s-velocity partials
c       dcdr    R*4 - array of layer density partials
c
c       ur  R*4 - array of radial eigenfunction
c       tur R*4 - array of radial stress eigenfunction
c       uz  R*4 - array of vertical eigenfunction
c       tuz R*4 - array of vertical stress eigenfunction
c       ipar    I*4 - array of integer controls
c
c-----
        parameter (NL=1000)
        real*4 dcdh(NL), dcda(NL), dcdb(NL), dcdr(NL)
        real*4 ur(NL), tur(NL), uz(NL), tuz(NL)
        integer*4 ipar(10)
            if(lorr.eq.5) then
                write(lun)wvno,u,gamma,
     1              sur,sdur,sd2ur,rur,rtr,sare
                if(ipar(4).eq.1)then
                write(lun)(ur(i),i=1,mmax)
                write(lun)(tur(i),i=1,mmax)
                endif
                if(ipar(5).eq.1)then
                write(lun)(dcdh(i),i=1,mmax)
                endif
                if(ipar(7).eq.1)then
                write(lun)(dcdb(i),i=1,mmax)
                endif
                if(ipar(8).eq.1)then
                write(lun)(dcdr(i),i=1,mmax)
                endif
            else if(lorr.eq.6) then
                write(lun)wvno,u,gamma,
     1              sur,sdur,sd2ur,suz,sduz,sd2uz,
     2              rur,rtr,ruz,rtz,
     3              sur0, rur0,sare
                if(ipar(4).eq.1)then
                write(lun)(ur(i),i=1,mmax)
                write(lun)(tur(i),i=1,mmax)
                write(lun)(uz(i),i=1,mmax)
                write(lun)(tuz(i),i=1,mmax)
                endif
                if(ipar(5).eq.1)then
                write(lun)(dcdh(i),i=1,mmax)
                endif
                if(ipar(6).eq.1)then
                write(lun)(dcda(i),i=1,mmax)
                endif
                if(ipar(7).eq.1)then
                write(lun)(dcdb(i),i=1,mmax)
                endif
                if(ipar(8).eq.1)then
                write(lun)(dcdr(i),i=1,mmax)
                endif
            endif
        return
        end
