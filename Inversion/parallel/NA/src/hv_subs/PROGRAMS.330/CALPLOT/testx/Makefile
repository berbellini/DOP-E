
CCMP	= cc
CFLAGS	= -O -I/usr/X11R6/include
LIB	= -L/usr/X11R6/lib -lX11 -lm

.c.o:
	$(CCMP) $(CFLAGS) -c $<

testx:	testx.o
	$(CCMP) $(LIB) testx.o  -o testx

clean:
	rm -f *.o
