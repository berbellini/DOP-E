c---------------------------------------------------------------------------
c
c	Subroutine read_seed - reads in a seed for the random number generator
c			       RAN0.
c			       To re-initialize call ran3 with a negative
c			       seed.
c
c	Note: Calls are made to : ran3
c---------------------------------------------------------------------------
c
	Subroutine read_seed(lu_seed, lu_out, iseed)
c
	common /GAswitches/ iwriteGA, idivGA
c
        read(lu_seed,*)iseed
	iseedo = iseed
	write(lu_out,*)' iseed = ',iseed
	write(lu_out,*)
	write(lu_out,*) '--------------------------------------------'
	write(lu_out,*)
	adummy = ran3(iseedo)
c
	return
	end
