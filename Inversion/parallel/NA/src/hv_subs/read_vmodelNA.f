c------------------------------------------------------------------------
c
c	Subroutine read_vmodelNA - reads in parameters associated with
c				 velocity model and 
c				 reads in all definitions to set up 
c				 the parameter space
c
c	Note_1: a priori model co-variance (scale factors) for i-th
c		parameter given by scales(i+1). scales(1) determines
c		the type of relative axis scaling in parameter space 
c		If scales(1) = 0 then; no scaling all other scale values ignored
c			         (This is only sensible when all parameters
c				  have the same physical units)
c		If scales(1) = -1 then; all parameter ranges are normalized 
c		                  and scales(j) j>1 are ignored 
c				  (Effectively the parameter range becomes
c				   the a priori co-variance for each parameter)
c		If scales(1) = anything else; then scales(j+1) is taken as the
c			       a priori co-variance for parameter j.
c
c	Note_2: Calls no other routines
c
c       Note_3:
c               range(1,1:nlayer): lower bound of thickness
c               range(2,1:nlayer): upper bound of thickness
c               range(1,nlayer+1:2*nlayer): lower bound of velocity_1
c               range(2,nlayer+1:2*nlayer): upper bound of velocity_1
c               range(1,2*nlayer+1:3*nlayer): lower bound of velocity_2
c               range(2,2*nlayer+1:3*nlayer): upper bound of velocity_2
c                 velocity_1 : S-wave velocity at upper interface
c                 velocity_2 : S-wave velocity at lower interface
c               range(1,3*nlayer+1:4*nlayer): lower bound of Vp/Vs ratio
c               range(2,3*nlayer+1:4*nlayer): upper bound of Vp/Vs ratio
c
c-----------------------------------------------------------------------
c
	subroutine read_vmodelNA(
     &          lu_vel, range, scales, moddim, 
     &          q_alpha, q_beta )
c
        include 'hv_param.inc'
c
	real*4		range(2,*)
c
	real*4		q_alpha(*),
     &			q_beta(*)
c
	real		scales(*)
c
	read(lu_vel,*) nlayer,scales(1)
c
	do i=1,nlayer
	  i2=nlayer+i
	  read(lu_vel,*) range(1,i), range(2,i), scales(i+1),
     &                   range(1,i2), range(2,i2), scales(i2+1)
	end do
	moddim=nlayer*2
c
	return
	end

