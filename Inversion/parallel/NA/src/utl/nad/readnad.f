c----------------------------------------------------------------------
c
c	readnad - Program to read Neighbourhood Algorithm Direct 
c	          access file and write out information to standard out.
c
c	Comments:
c		  This program also acts as an example of how to call the 
c		  routine read_nad which reads a direct access nad
c		  file and puts output into arrays models, data, header
c		  and variables nd,ne,nh.
c
c	Calling sequence:
c
c			 readnad file.in 
c
c	where,
c	     `file.in'  = input direct access `*.nad' file (e.g. na.nad)
c	
c----------------------------------------------------------------------
c
 	Program readnad
c
c	parameter	(nemax=201000,ndmax=24,nhmax=10000)
 	parameter	(nemax=41000,ndmax=450,nhmax=10000)
	
	real*4		models(ndmax*nemax)
	real*4		data(nemax)
	real*4		ranges(2,ndmax)
	real*4		scales(ndmax+1)
	character	header(nhmax)

      	character*256 	fnme
      	character*256 	answer

        write(*,100)

        m = iargc()
        if(m.ne.1)then
           write(*,*)' Incorrect number of arguments'
           write(*,*)
           write(*,*)' Syntax:'
           write(*,*)'        readnad filein'
           write(*,*)
           stop
        end if

c						read in filenames
        call getarg(1,answer)
        read(answer,fmt='(a20)')fnme
c						read in ensemble of models
c						from direct access nad file 
        lu_in = 15
c	write(*,*)' Reading direct access file...'

	call read_nad
     &               (lu_in,fnme,nh,nhu,iform,nd,ne,
     &                nhmax,ndmax,nemax,
     &                header,data,models)

c	write(*,*)' Finished reading direct access file'
c
c						write out ensemble of models
c						in rfi ascii format
	nha = nh - nhu
	write(*,*)' Summary of nad file read in'
	write(*,*)
	write(*,*)' Number of variables                  = ',nd
	write(*,*)' Number of models                     = ',ne
	write(*,*)' Total length of header               = ',nh
	write(*,*)' Length of reserved header            = ',nha
	write(*,*)' Length of user header                = ',nhu
	a = data(1)
	am = data(1)
	avg = 0.
	mopt = 1
	do i=1,ne
	   if(data(i).lt.a)then
              a = data(i)
	      mopt = i
	   end if
	   if(data(i).gt.am)then
              am = data(i)
	      moptm = i
	   end if
	   avg = avg + data(i)
	end do
	avg = avg/real(ne)
	write(*,*)
	write(*,*)' Minimum data value                   = ',a
	write(*,*)' Maximum data value                   = ',am
	write(*,*)' Average misfit over all models       = ',avg
	write(*,*)' Model with minimum data value        = ',mopt
	write(*,*)' Model with maximum data value        = ',moptm

	open(50,file='model.opt',status='unknown')
	do i=1,nd
	   jj = i+(mopt-1)*nd
	   write(50,*)models(jj)
	end do
c
c						read nad header

	if(nhu.eq.nh)then
	   write(*,*)' No NA header in NAD file'
	else
	   lu = 20
 	   call read_NAheader
     &          (lu_in,fnme,nd,nh,nhu,iform,
     &           ranges,scales,ns1,ns,nr)

	   write(*,*)' Number of samples in first iteration = ',ns1
	   write(*,*)' Number of samples per iteration      = ',ns
	   write(*,*)' Number of cells per iteration        = ',nr

	   write(*,*)
	   write(*,*)' scale parameter                      = ',
     &               scales(1)
	   write(*,*)
	   write(*,*)' Boundaries of parameter space:'
	   write(*,*)
	   write(*,300)
	   do i=1,nd
	      write(*,200)i,ranges(1,i),ranges(2,i),scales(i+1)
	   end do

	end if

 100    format(/"  Program readnad - reads",
     &          " direct access nad file "/
     &          20x,"and writes out information to screen."/)

 200    format(1x,i3,2x,':',2x,f9.5,2x,f9.5,2x,f9.5)
 300    format(13x,'Min',8x,'Max',4x,'Scale factor')

	stop
	end
c
c ---------------------------------------------------------------------
c
c       read_nad - read a direct access file in NAD format. 
c
c	Input:
c	      lu	      : logical unit of file
c	      fnme	      : filename 
c	      iform	      : If iform=-1 on input then just read header and return 
c	      nhmax	      : maximum size of array header
c	      ndmax	      : maximum size of array data
c	      nemax	      : maximum size of array models
c
c	Output:
c	      nh	      : total length in bytes of file header 
c	      nhu	      : length in bytes of user written part of header 
c	      iform	      : format for nad file (1=multi-record,0=single)
c	      mul	      : number of records in the header (with padding)
c	      nd	      : dimension of parameter space
c	      ne	      : number of models in ensemble
c	      header          : header character string of length nh (char)
c	      data(nd)	      : array of data values for each model (real*4)
c	      models(nd,ne)   : array of model values  (real*4)
c
c	Comments:
c                The direct access NAD file format:
c
c                VARIABLE       TYPE            SIZE IN BYTES
c                nd             int             4
c                ne             int             4
c                nh             int             4
c                nhu            int             4
c                header         character       nh
c                models         real*4          4*nd*ne
c                data           real*4          4*nd
c                tail           real*4          4
c
c                File must contain a single record of length
c                [4x(4+nd*ne+ne+1) + nh] bytes
c
c		Calls are made to subroutine read_da. 
c
c                This routine assumes that direct access
c                files are opened with the record length specified
c                in bytes. This is the default for most machines
c                but not Dec machines. (A compiler option is
c                available on the Compaq Alpha to use bytes rather than 
c                4-byte words.)
c
c					M. Sambridge, RSES, Nov. 2001
c
c ---------------------------------------------------------------------
c
      	Subroutine read_nad
     &                     (lu,fnme,nh,nhu,iform,nd,ne,
     &                      nhmax,ndmax,nemax,
     &                      header,data,models)
c
      	real*4            models(*)
      	real*4            data(*)
      	character*(*)     header
        logical           headeronly

      	character*256 	  fnme

	mul = 0
        headeronly = .false.
        if(iform.eq.-1)then
          headeronly = .true.
          iform = 0
        end if
	iform = 0
c						read in size of header
        len = 4
        open(lu,file=fnme,status='unknown',
     &       form='unformatted',access='direct',recl=len,err=100)
        read(lu,rec=1)mul
        close(lu)
c						read in number of models,
c						number of dimensions and 
c						size of header
	if(mul.gt.0)then
c						we are in single record format
          len = 16
          open(lu,file=fnme,status='unknown',
     &         form='unformatted',access='direct',recl=len)
          read(lu,rec=1)nd,ne,nh,nhu
          mul = 0
          close(lu)
	else
c						we are in multi record format
          len = 20
          open(lu,file=fnme,status='unknown',
     &         form='unformatted',access='direct',recl=len)
          read(lu,rec=1)mul,nd,ne,nh,nhu
          mul = -mul
          close(lu)
	  iform = 1
	end if
        if(headeronly)return
c						check array sizes
	if(nd.gt.ndmax)then
          write(*,*)
          write(*,*)
     &    'Error - maxmimum number of dimensions is too small'
          write(*,*)
     &    '        current maximum number of dimensions =',ndmax
          write(*,*)
     &    '        number of dimensions in input file   =',nd
          write(*,*)'        '
          write(*,*)'Remedy - increase size of parameter ndmax'
          write(*,*)'         to at least this value and recompile'
          write(*,*)'        '
          stop
        end if

	if(ne.gt.nemax)then
          write(*,*)
          write(*,*)
     &    'Error - maxmimum number of models is too small'
          write(*,*)
     &    '        current maximum number of models =',nemax
          write(*,*)
     &    '        number of models in input file   =',ne
          write(*,*)'        '
          write(*,*)'Remedy - increase size of parameter nemax'
          write(*,*)'         to at least this value and recompile'
          write(*,*)'        '
          stop
        end if

	if(nh.gt.nhmax)then
          write(*,*)
          write(*,*)
     &    'Error - maxmimum size of header is too small'
          write(*,*)
     &    '        current maximum size of header =',nhmax
          write(*,*)
     &    '        size of header in input file   =',nh
          write(*,*)'        '
          write(*,*)'Remedy - increase size of parameter nhmax'
          write(*,*)'         to at least this value and recompile'
          write(*,*)'        '
          stop
        end if

 	call read_da(lu,fnme,nd,ne,nh,mul,header,models,data)

      	return
 100    write(*,*)
        write(*,*)' Error reading file : ',fnme
        write(*,*)' Does it exist ? '
        write(*,*)
	stop

      	end
c
c ---------------------------------------------------------------------
c
c       read_da - read a direct access file containing
c                 ne models with dimension nd, ne data 
c                 values and a header of size nh.
c
c	Input:
c	      lu		: logical unit of file
c	      fnme		: filename 
c	      nh		: length in bytes of file header (minus padding)
c	      mul		: number of records in the header (with padding)
c	      nd		: dimension of parameter space
c	      ne		: number of models in ensemble
c
c	Output:
c	      header		: header character string of length nh (char)
c	      data(nd)		: array of data values for each model (real*4)
c	      models(nd,ne)	: array of model values  (real*4)
c
c	Comments:
c		 The direct access NAD file format:
c
c		 VARIABLE	TYPE		SIZE IN BYTES
c		 nd		int     	4
c		 ne		int     	4
c		 nh		int     	4
c		 header		character     	nh
c		 models		real*4     	4*nd*ne
c		 data		real*4     	4*nd
c		 tail		real*4     	4
c	
c	         File must contain a single record of length 
c		 [4x(3+nd*ne+ne+1) + nh] bytes
c
c
c					M. Sambridge, RSES, Nov. 2001
c
c ---------------------------------------------------------------------
c
      	Subroutine read_da(lu,fnme,nd,ne,nh,mul,header,models,data)
c
      	real*4            models(nd,ne)
      	real*4            data(ne)
      	character         header(nh)
c     	real*4            tail

      	character*256 	  fnme
c					
c                                              calculate length of header
 	iform = 0
 	if(mul.ne.0)iform = 1

 	lenh  = 4*5+nh
        len   = 4*(nd+1)

 	if(iform.eq.1)then

           open(lu,file=fnme,status='unknown',
     &          form='unformatted',access='direct',recl=lenh)

       	   read(lu,rec=1)idum1,idum2,idum3,idum4,idum5,header
	
 	   close(lu)

c	   write(*,*)' record length               ',len
c	   write(*,*)' header length               ',lenh
c	   write(*,*)' Number of records in header ',mul
c
           open(lu,file=fnme,status='unknown',
     &          form='unformatted',access='direct',recl=len)

           do i=1,ne
              call rnad(lu,i+mul,nd,models(1,i),data(i))
           end do

           close(lu)

 	else
c					       read in file as 
c					       a single record
	   len = 4+nd*ne+ne
           len = 4*len+nh

           open(lu,file=fnme,status='unknown',
     &          form='unformatted',access='direct',recl=len)

           read(lu,rec=1)i,j,k,kk,header,models,data

           close(lu)

 	end if
c
      	return
      	end
c
c ---------------------------------------------------------------------
c
c       read_NAheader - extracts the NA program information placed
c			in the header block of a NAD file 
c
c ---------------------------------------------------------------------
c
 	Subroutine read_NAheader
     &             (lu,fnme,nd,nh,nhu,iform,
     &              ranges,scales,ns1,ns,nr)

	real*4		ranges(2,nd)
	real*4		scales(nd+1)
	character*256   fnme	
	
	if(nhu.ge.nh)return

	if(iform.eq.0)then
	   len = 7+nd*3 + 1
           len = 4*len
           open(lu,file=fnme,status='unknown',
     &          form='unformatted',access='direct',recl=len)
           read(lu,rec=1)i,j,k,kk,ns1,ns,nr,ranges,scales
           close(lu)
	else
           len = 4*5+nh
           open(lu,file=fnme,status='unknown',
     &          form='unformatted',access='direct',recl=len)
       	   read(lu,rec=1)idum1,idum2,idum3,idum4,idum5,
     &                   ns1,ns,nr,ranges,scales
           close(lu)
	end if

	return
	end

      	Subroutine rnad(lu,i,nd,models,data)

      	real*4            models(nd)
      	real*4            data

      	read(lu,rec=i)models,data

      	return
      	end

