
#include "calplot.h"
main()
{
	int i;
	pinitf("WID.PLT");
	gunit("cm");
	for(i=1;i<=2;i++){
		if(i == 1){
			gunit("in");
			gwidth(0.25);
			plot(3.0,3.0,3);
			plot(5.0,5.0,2);
			plot(3.0,3.0,3);
			plot(5.0,4.0,2);
			plot(3.0,3.0,3);
			plot(4.0,5.0,2);
			plot(3.0,3.0,3);
			plot(3.0,5.0,2);
			plot(3.0,3.0,3);
			plot(5.0,3.0,2);
		} else {
			gunit("cm");
			gwidth(1.0);
			plot(1.0,1.0,3);
			plot(6.0,6.0,2);
			plot(1.0,1.0,3);
			plot(6.0,3.0,2);
			plot(1.0,1.0,3);
			plot(3.0,6.0,2);
			plot(1.0,1.0,3);
			plot(6.0,1.0,2);
			plot(1.0,1.0,3);
			plot(1.0,6.0,2);
		}
	}
	pend();
}
