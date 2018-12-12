c----------------------------------------------------------------------
c
c	nad2asc - Program to read Neighbourhodd Algorithm Direct 
c		  access file and write out an ascii `models' file. 
c
c	Comments:
c		  This program also acts as an example of how to call the 
c		  routine read_nad which reads a direct access nad
c		  file and puts output into arrays models, data, header
c		  and variables nd,ne,nh.
c
c	Calling sequence:
c
c			 nad2asc file.in file.out 
c
c	where,
c	     `file.in'  = input direct access `*.nad' file (e.g. `NAV22.nad')
c	     `file.out' = output ascii file 
c	
c----------------------------------------------------------------------
c
 	Program nad2asc
c
	parameter	(nemax=11000,ndmax=24,nhmax=10000)
	
	real*4		models(ndmax*nemax)
	real*4		data(nemax)
	character	header(nhmax)
        real*4          ranges(2,ndmax)
        real*4          scales(ndmax+1)

      	character*256 	fnme
      	character*256 	fnme2
      	character*256 	answer

        write(*,100)

        m = iargc()
        if(m.ne.2)then
           write(*,*)' Incorrect number of arguments'
           write(*,*)
           write(*,*)' Syntax:'
           write(*,*)'        nad2asc filein fileout'
           write(*,*)'        (e.g.   NAV22.nad NAV22.mdl)'
           write(*,*)
           stop
        end if

c						read in filenames
        call getarg(1,answer)
        read(answer,fmt='(a20)')fnme
        call getarg(2,answer)
        read(answer,fmt='(a20)')fnme2
c						read in ensemble of models
c						from direct access nad file 
        lu_in = 15
	write(*,*)' Reading direct access file...'

	call read_nad
     &               (lu_in,fnme,nh,nhu,iform,nd,ne,
     &                nhmax,ndmax,nemax,
     &                header,data,models)

	write(*,*)' Finished reading direct access file'
c
c                                               read nad header

        if(nhu.eq.nh)then
           write(*,*)' No NA header in NAD file'
        else
           call read_NAheader
     &          (lu_in,fnme,nd,nh,nhu,iform,
     &           ranges,scales,ns1,ns,nr)
	end if
c
c						write out ensemble of models
c						in rfi ascii format
	lu = 10

	write(*,*)
	write(*,*)' Writing ascii file...'
	write(*,*)

        call write_asc
     &             (lu,fnme,fnme2,nd,ne,nh,nhu,iform,
     &              ranges,scales,ns1,ns,nr,models,data)

 100    format(/" Program nad2asc - converts",
     &          " direct access nad file to"/
     &          18x," ascii receiver function model file."/)

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

      	Subroutine rnad(lu,i,nd,models,data)

      	real*4            models(nd)
      	real*4            data

      	read(lu,rec=i)models,data

      	return
      	end

c
c ---------------------------------------------------------------------
c
c       write_rfi - write ensemble of models produced by
c                   narfi and garfi programs.
c
c ---------------------------------------------------------------------
c
        Subroutine write_asc
     &             (lu,fnme,fnme2,nd,ne,nh,nhu,iform,
     &              ranges,scales,ns1,ns,nr,models,data)

        real*4          models(nd,ne)
        real*4          data(nd)
        real*4          ranges(2,nd)
        real*4          scales(nd+1)

      	character*256 	fnme2
      	character*256 	fnme

	open(lu,file=fnme2,status='unknown')

	write(lu,*)
	write(lu,fmt='("Ascii version of NAD file : ",a30)')fnme
	write(lu,*)
	write(lu,*)nd,'    : number of variables '
	write(lu,*)ne,'    : number of models in ensemble '
	write(lu,*)nh,'    : Length of header in bytes '
c	write(lu,*)nhu,'    : Length of user header in bytes '
	write(lu,*)ns1,'    : Number of samples at iteration 1'
	write(lu,*)ns,'    : Number of samples at other iterations'
	write(lu,*)nr,'    : Number of cells to resample'
	write(lu,*)scales(1),' : value of scales (1) '
	write(lu,*)
	write(lu,*)' Ranges and scale factors '
	do k=1,nd
           write(lu,*)k,ranges(1,k),ranges(2,k),scales(k+1)
	end do
	write(lu,*)
	write(lu,*)' Model and data values'
	write(lu,*)
	do i=1,ne
           write(lu,*)i,(models(k,i),k=1,nd),data(i)
	   write(lu,*)
	end do

	close(lu)

	return
	end
c ---------------------------------------------------------------------
c
c       read_NAheader - extracts the NA program information placed
c                       in the header block of a NAD file
c
c ---------------------------------------------------------------------
c
        Subroutine read_NAheader
     &             (lu,fnme,nd,nh,nhu,iform,
     &              ranges,scales,ns1,ns,nr)

        real*4          ranges(2,nd)
        real*4          scales(nd+1)
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

