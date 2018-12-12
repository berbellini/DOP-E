#include        <stdio.h>
#include        <stdlib.h>
#include        <string.h>
#include	<ctype.h>

#include        "gsac.h"
#include	"gsac_docommand.h"

/* CHANGES
	04 SEP 2006 - command line aprsing ignores everything
		including and after a # - this is to permit
		comments
*/

int gsac_docommand(int ncmd, char  **cmdstr, char *input_lineptr);

/* get a line from the input 
 * separate it into different tokens (e.g., words) in the manner of
 * argc, argv;
 */
int gsac_tokenize_line(char *line, int *nstr, char **str);
int nstr;
char *gsac_strupr(char *s);
char *mystrtok (char *str, char *delims);
int  fgetline(FILE *fp, char s[], int lim, char *prompt);
int  sfgetline(FILE *fp, char s[], int lim);

#define NTOKEN 100
char *strcmd[NTOKEN];
#define NCHAR 2000
char input_lineptr[NCHAR] ;
char save_lineptr[NCHAR] ;

int gsac_parse_command(void)
{
	int lstr;
	int quit = 0 ;
	/* note the !quit test must be evaluated  first  */
	while(!quit && (fgetline(stdin,input_lineptr,NCHAR, "GSAC> ") != EOF) ){
		lstr = strlen(input_lineptr);
		if(lstr > 0){
			strcpy(save_lineptr,input_lineptr);
			gsac_tokenize_line(input_lineptr, &nstr, strcmd);
			quit =gsac_docommand(nstr, strcmd,save_lineptr);
			return (!quit);
		}
	}
	/* LATER WE WILL CALLOC AND REALLOC THIS 
		free(input_lineptr);
	*/
	return (quit);
}


int gsac_tokenize_line(char *line, int *nstr, char **str)
{
	/* separate the line by tokens
	 * determine the number of tokens
	 * increase the memory allocation of **str and *str
	 */
	char *break_set = "\t\n\'\" ";
	char *tokp, *sp;
	int n = 0 ;

	sp = line;
	/* get the number of tokens */
	/* this works since the fgetline filles an existing array.
	 * the mystrtok will actually replace the space-tab-quote
	 * by a NULL - and then sets pointers to the beginning of
	 * each substring. So no alloc or realloc is required */
	while ((tokp = mystrtok(sp,break_set)) != NULL )
	{
		str[n] = tokp;
		n++;
		sp = NULL ;
	}
	*nstr = n;
	return(n);
}

/*  http://web.cs.mun.ca/courses/cs3710-w98/c/overheads/examples/strtok.c */

/* modifies 04 SEP 2006 so that anything after and including a # is
	ignored unless qithin quotes or doublequotes
*/

char *mystrtok (char *str, char *delims)
{
  static char *pos = (char *) 0;
  char *start = (char *) 0;
  int foundquote = 0 ;
  int founddoublequote = 0 ;

  if (str)			/* Start a new string? */
    pos = str;

  if (pos) {
      /* Skip delimiters */
      while (*pos && strchr (delims, *pos)){
	/* this get a leading # for the situations
		one # two
		one #two
	   but not for a trailing
		one#two
	*/
	if(*pos == '#'){
		start = (char *) 0;
		return start;
	}
	if(*pos == '\''){
		foundquote = 1;
/*
	fprintf("found quote\n");
*/
	}
	if(*pos == '\"'){
		founddoublequote = 1;
/*
	fprintf("found doublequote\n");
*/
	}
	pos++;
      }
      if (*pos) {
	  start = pos;
	  /* if quote or doublequote skip until the next is found */
	  /* Skip non-delimiters */
	  if(foundquote == 1){
	  		while (*pos && !strchr ("\'", *pos))
	    		pos++;
	  		if (*pos)
	    			*pos++ = '\0';
			foundquote = 0;
          } else if(founddoublequote == 1){
	  		while (*pos && !strchr ("\"", *pos))
	    		pos++;
	  		if (*pos)
	    			*pos++ = '\0';
			founddoublequote = 0;
          } else {
	       while (*pos && !strchr (delims, *pos) )
	             pos++;
	          if ( *pos == '#' )
	             *pos = '\0';
	          else if ( *pos )
	             *pos++ = '\0';
	  
          }
	}
    }
    return start;
}

char *gsac_strupr(char *s)
{
	if (s != NULL )
	{
		char *p;

		for( p=s ; *p; p++)
			*p = toupper(*p);
	}
	return s;
}

#ifdef READLINE_LIBRARY
#  include "readline.h"
#  include "history.h"
int fgetline(FILE *fp, char s[], int len, char *prompt) {
	char *temp;
	s[0] = '\0';
	temp = readline (prompt);
	if (!temp)
		return EOF;
	if (*temp){
		if(strlen(temp) > len-1){
			printf("Line too long (%d characters versus %d in input buffer). Split line\n",
				(int)strlen(temp),len-1);
		} else {
			add_history (temp);
			strcpy(s,temp);
		}
	}
	free(temp);
	return 1;
}

#else

/* http://www.math.ncu.edu.tw/~shann/Teach/comp/c-samp/getline.htm */

int fgetline(FILE *fp, char s[], int len, char *prompt) {
    char *t;
    int c, lim;

	printf("%s", prompt);
    lim=len;

    t = s;
    while (--lim>1 && (c=getc(fp)) != EOF && c != '\n'){
        *s++ = c;
    }
    if (c == '\n')
        *s++ = '\0';
    else if (lim == 1) {
	*s++ = '\0';
	fprintf(stderr, "WARNING. fgetline: Line too long, split.\n");
    } else if(c == EOF) {
	    return EOF;
    }
    *s = '\0';
    return s - t;
}

#endif

/* I need to read a line and parse tokens in the gsac_trans.c However I
	do not require the Prompts 
*/
/* http://www.math.ncu.edu.tw/~shann/Teach/comp/c-samp/getline.htm */
int sfgetline(FILE *fp, char s[], int len) {
    char *t;
    int c, lim;

    lim=len;

    t = s;
    while (--lim>1 && (c=getc(fp)) != EOF && c != '\n'){
        *s++ = c;
    }
    if (c == '\n')
        *s++ = '\0';
    else if (lim == 1) {
	*s++ = '\0';
	fprintf(stderr, "WARNING. fgetline: Line too long, split.\n");
    } else if(c == EOF) {
	    return EOF;
    }
    *s = '\0';
    return s - t;
}

