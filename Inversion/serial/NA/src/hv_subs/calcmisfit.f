c-------------------------------------------------------------------------
c
c	Subroutine calcmisfit - calculates a misfit value between 
c				observed data and predicted data
c				plus model roughness
c
c	Note: Calls no other routines
c
c------------------------------------------------------------------------
c
	subroutine calcmisfit(
     &		predicted_data, observed_data, sd, ndata, 
     &		weight, nwave, misfitval )
c
c
	include 'hv_param.inc'
c
c
	real*4		predicted_data(maxdata,maxwave),
     &			observed_data(maxdata,maxwave),
     &			sd(maxdata,maxwave),
     &			weight(maxwave),
     &			misfitval, aval
c
	integer         ndata(maxwave)
c
c						misfit between observed and
c						predicted
	misfitval=0.0


	do iw=1,nwave
c
	  do i=1,ndata(iw)
c	    write(*,*) observed_data(i,iw), predicted_data(i,iw)		
	    aval=((observed_data(i,iw)-predicted_data(i,iw))**2/(sd(i,iw))**2 )* weight(iw)
	    misfitval=misfitval+aval

	  end do

	end do

	

	return
	end
