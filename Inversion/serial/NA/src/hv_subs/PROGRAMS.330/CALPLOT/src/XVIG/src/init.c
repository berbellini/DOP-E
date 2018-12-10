/* File>>> init.c
--
-- %M% -- version %I% (IMEC)            last updated: %E%
--
-- Copyright (c) 1993
-- IMEC vzw
-- Kapeldreef 75
-- B-3001 LEUVEN
-- BELGIUM
--
-- Author   : A. Demaree
--
-- Date     : February 1, 1993
--
-- Function :
--
-- Comment  :
--
-- Review   :
--
*/


/*------------------------------------------------------------------------------
-- Global include files
------------------------------------------------------------------------------*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/file.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#ifdef SOLARIS
	/* beg RBH */
#include <libgen.h>
#include <X11/Xw32defs.h>
#endif
#include <unistd.h>
	/* end RBH */

/*------------------------------------------------------------------------------
-- Local include files
------------------------------------------------------------------------------*/

#define MAIN
#include "xviglocal.h"
#undef MAIN

/*------------------------------------------------------------------------------
-- Static variable declarations
------------------------------------------------------------------------------*/


/*------------------------------------------------------------------------------
-- Local function declarations
------------------------------------------------------------------------------*/


/*------------------------------------------------------------------------------
--
--
--
-------- ---------------------------------------------------------------------*/
static char *my_pathfind( char *path, const char *name, const char *mode);

char *my_pathfind( char *path, const char *name, const char *mode)
{
	char *fpath;
	char *token;
	token = (char *)strtok(path,":");
	while(token != NULL){
		fpath = calloc(strlen(token)+strlen(name)+2,sizeof(char));
		strcpy(fpath,token);
		strcat(fpath,"/");
		strcat(fpath,name);
		if(access(fpath, R_OK)==0)return(fpath);
		free(fpath);
		token = strtok(NULL,":");
	}
return(NULL);

}

int XviG_Init(char *classname,
              int color_array[][3],
              int nr_of_colors)
{
  XGCValues xgcv;
	/* beg RBH */
	int depth;
	int colorset;
	char *PATH;
	char *tpath;
	char *path;
	tpath = getenv("PATH");
	path=calloc(strlen(tpath)+1,sizeof(char));
	strcpy(path,tpath);
 
        XviG_pathname = (char *)my_pathfind(path, "calxvig", "rx");
fprintf(stderr," \n");
fprintf(stderr,"%s\n",XviG_pathname);
fprintf(stderr," \n");
	if(!XviG_pathname ){
		fprintf(stderr,
			"ERROR : Cannot find executable 'calxvig' in PATH\n");
		return 0;
	}

	/* end RBH */

  printf("\n>>> XviG Graphics Package   Copyright Imec (c) 1993 <<<\n\n");
  printf("Executing: %s\n",XviG_pathname);

  /*
  -- Get the environment variable 'XVIG' and check it

  if (!(XviG_pathname = getenv("XVIG")))
  {
    fprintf(stderr,
            "ERROR : Environment variable 'XVIG' has not been defined.\n");
    return 0;
  }

  if (access(XviG_pathname, F_OK) == -1)
  {
    fprintf(stderr, "ERROR : File '%s' not found.\n", XviG_pathname);
    return 0;
  }

  if (access(XviG_pathname, X_OK) == -1)
  {
    fprintf(stderr, "ERROR : File '%s' is not executable.\n", XviG_pathname);
    return 0;
  }
*/

  /*
  -- Open the display
  */

  if (!(XviG_display = XOpenDisplay(NULL)))
  {
    char *dname;

    if (!(dname = getenv("DISPLAY")))
      fprintf(stderr,
              "ERROR : Environment variable 'DISPLAY' has not been defined.\n");
    else
      fprintf(stderr,
              "ERROR : Cannot open display '%s'.\n", dname);

    return 0;
  }

  XviG_screen_nr = DefaultScreen(XviG_display);
  depth = DisplayPlanes(XviG_display,XviG_screen_nr);


  /*
  -- Create a dummy window
  */

  if (!(XviG_dummy_window = XCreateSimpleWindow(XviG_display,
                                   RootWindow(XviG_display, XviG_screen_nr),
                                   0, 0, 10, 10, 0,
                                   WhitePixel(XviG_display, XviG_screen_nr),
                                   BlackPixel(XviG_display, XviG_screen_nr))))
  {
    fprintf(stderr, "ERROR : Cannot create dummy window.\n");
    XCloseDisplay(XviG_display);

    return 0;
  }

  XSelectInput(XviG_display, XviG_dummy_window, KeyPressMask | ButtonPressMask);

  /*
  -- About the Graphical Context, create a separate one for drawing simple
  -- objects and for drawing filled objects to avoid interference
  */

  xgcv.graphics_exposures = False;
  XviG_gc = XCreateGC(XviG_display, XviG_dummy_window,
                      GCGraphicsExposures, &xgcv);
  XviG_gc_fill = XCreateGC(XviG_display, XviG_dummy_window,
                           GCGraphicsExposures, &xgcv);

  /*
  -- Create the colors and the fillpatterns
  */

  colorset = XviG_CreateColors(color_array, nr_of_colors,depth);
  if(colorset < 5)
	depth = 1;
  XviG_CreateFillpatterns();

  /*
  -- Create the Atoms (unique numbers) for the ClientMessage events
  */

  XviG_MESSAGE_INIT =
                  XInternAtom(XviG_display, MESSAGE_INIT, False);
  XviG_MESSAGE_KEY =
                  XInternAtom(XviG_display, MESSAGE_KEY, False);
  XviG_MESSAGE_BUTTON =
                  XInternAtom(XviG_display, MESSAGE_BUTTON, False);
  XviG_MESSAGE_KEY_BUTTON =
                  XInternAtom(XviG_display, MESSAGE_KEY_BUTTON, False);
  XviG_MESSAGE_SIZE =
                  XInternAtom(XviG_display, MESSAGE_SIZE, False);
  XviG_MESSAGE_SENSEKBD_ON =
                  XInternAtom(XviG_display, MESSAGE_SENSEKBD_ON, False);
  XviG_MESSAGE_SENSEKBD_OFF =
                  XInternAtom(XviG_display, MESSAGE_SENSEKBD_OFF, False);
  XviG_MESSAGE_SENSEKBD =
                  XInternAtom(XviG_display, MESSAGE_SENSEKBD, False);
  XviG_MESSAGE_CURSOR =
                  XInternAtom(XviG_display, MESSAGE_CURSOR, False);
  XviG_MESSAGE_QUIT =
                  XInternAtom(XviG_display, MESSAGE_QUIT, False);

  /* extensions rbh */
  XviG_MESSAGE_CLIP =
                  XInternAtom(XviG_display, MESSAGE_CLIP, False);
  XviG_MESSAGE_REVERSE =
                  XInternAtom(XviG_display, MESSAGE_REVERSE, False);
  XviG_MESSAGE_TOGRAY =
                  XInternAtom(XviG_display, MESSAGE_TOGRAY, False);
  XviG_MESSAGE_LCMAP_SIZE =
                  XInternAtom(XviG_display, MESSAGE_LCMAP_SIZE, False);
  XviG_MESSAGE_LCMAP_ENTRY =
                  XInternAtom(XviG_display, MESSAGE_LCMAP_ENTRY, False);
  XviG_MESSAGE_LCMAP_RESET =
                  XInternAtom(XviG_display, MESSAGE_LCMAP_RESET, False);
  XviG_MESSAGE_BOUNDS =
                  XInternAtom(XviG_display, MESSAGE_BOUNDS, False);
  XviG_MESSAGE_ORIGIN_HYPERBOLA =
                  XInternAtom(XviG_display, MESSAGE_ORIGIN_HYPERBOLA, False);

  /*
  -- About the Resource Manager
  */

  if (classname)
  {
    XviG_class = STRNEW(classname);
    XviG_GetRMDatabase();
  }
  else
    XviG_class = (char *) NULL;

  XviG_cursor_mode = False;

  return depth;
}

/*------------------------------------------------------------------------------
--
--
--
------------------------------------------------------------------------------*/

void XviG_Exit(void)
{
  XviG_DeleteRMDatabase();

  XviG_CleanupWindows();
  XviG_CleanupCursors();
  XviG_CleanupFillpatterns();
  XviG_CleanupFonts();

  XFreeGC(XviG_display, XviG_gc);
  XFreeGC(XviG_display, XviG_gc_fill);

  XDestroyWindow(XviG_display, XviG_dummy_window);
  XCloseDisplay(XviG_display);
}

