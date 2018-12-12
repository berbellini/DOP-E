/* File>>> xviglocal.h
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


#ifndef __XVIGLOCAL_H
#define __XVIGLOCAL_H


#include "xvig.h"
#include "commondata.h"


/*
-- Global variables
*/

#ifdef MAIN
#define EXTERN
#else
#define EXTERN extern
#endif  /* MAIN */

EXTERN Display *XviG_display;
EXTERN int XviG_screen_nr;
EXTERN Window XviG_dummy_window, XviG_window;
EXTERN Pixmap XviG_pixmap;
EXTERN GC XviG_gc, XviG_gc_fill;
EXTERN XEvent XviG_event;
EXTERN unsigned int XviG_window_width, XviG_window_height;
EXTERN char *XviG_pathname, *XviG_class;
EXTERN Bool XviG_cursor_mode;
EXTERN unsigned long XviG_xhair_color, XviG_save_color;
EXTERN unsigned int XviG_save_linewidth;
EXTERN int XviG_save_linestyle;

EXTERN Atom XviG_MESSAGE_INIT,
            XviG_MESSAGE_KEY,
            XviG_MESSAGE_BUTTON,
            XviG_MESSAGE_KEY_BUTTON,
            XviG_MESSAGE_SIZE,
            XviG_MESSAGE_SENSEKBD_ON,
            XviG_MESSAGE_SENSEKBD_OFF,
            XviG_MESSAGE_SENSEKBD,
            XviG_MESSAGE_CURSOR,
            XviG_MESSAGE_QUIT;

/* rbh extension */

EXTERN Atom XviG_MESSAGE_CLIP,
            XviG_MESSAGE_REVERSE,
            XviG_MESSAGE_TOGRAY,
            XviG_MESSAGE_LCMAP_SIZE,
            XviG_MESSAGE_LCMAP_ENTRY,
            XviG_MESSAGE_LCMAP_RESET,
	    XviG_MESSAGE_ORIGIN_HYPERBOLA,
            XviG_MESSAGE_BOUNDS;

EXTERN	int RBH_xoring_state;	/* to get around color map problems for xoring */

/*
-- General use macros
*/

#define Abs(n)    ((n) < 0 ? -(n) : (n))
#define Max(a,b)  ((a) > (b) ? (a) : (b))
#define Min(a,b)  ((a) < (b) ? (a) : (b))
#define STRNEW(str)      (strcpy((char *)malloc(strlen(str)+1),(str)))
#define STRNNEW(str,len) (strncpy((char *)malloc((len)+1),(str),(len)))


/*
-- For the initialisation
*/

extern int XviG_CreateColors(int colmap[][3], int nr, int depth);
extern void XviG_CreateFillpatterns(void);


/*
-- For cleanup when exiting
*/

extern void XviG_CleanupWindows(void);
extern void XviG_CleanupCursors(void);
extern void XviG_CleanupFillpatterns(void);
extern void XviG_CleanupFonts(void);


/*
-- For accessing the Resource Manager database
*/

extern void XviG_GetRMDatabase(void);
extern void XviG_DeleteRMDatabase(void);
extern char *XviG_GetRMGeometry(char *name);
extern int XviG_GetRMMaxFonts(void);
extern char *XviG_GetRMFontName(int nr);


/*
-- For maintaining the window list
*/

typedef struct WinItem aWinItem, *WinItem;

struct WinItem {
  char *name;
  int pid;
  Window window;
  Pixmap pixmap;
  unsigned int width, height;
  WinItem prev;
  WinItem next;
};

EXTERN WinItem selected_winitem;

extern void New_WinItem(char *name, int pid);
extern void Delete_WinItem(WinItem curwin);
extern WinItem Get_WinItem(char *name);


/*
-- For maintaining the list of user cursors
*/

typedef struct CurItem aCurItem, *CurItem;

struct CurItem {
  Pixmap pixmap;
  unsigned int width, height;
  int hot_x, hot_y;
  CurItem prev;
  CurItem next;
};

extern CurItem New_CurItem(Pixmap pixmap,
                           unsigned int width, unsigned int height,
                           int hot_x, int hot_y);
extern void Delete_CurItem(CurItem curcur);


#undef EXTERN

#endif  /* __XVIGLOCAL_H */
