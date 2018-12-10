#include <stdio.h>

main()
{
	int i, j ;
	int c;
	for (i=0;i < 32; i++){
		printf("%2d ",i);
		for (j=5; j>=0 ; j--){
			c = ( ((char)i & 0x3f ) >> j) & 01 ;
			if(c == 0 )
				printf(".");
			else
				printf("X");
		}
		for (j=5; j>=0 ; j--){
			c = ( ((char)i & 0x3f ) >> j) & 01 ;
			if(c == 0 )
				printf(".");
			else
				printf("X");
		}
		printf("\n");

	}
}
