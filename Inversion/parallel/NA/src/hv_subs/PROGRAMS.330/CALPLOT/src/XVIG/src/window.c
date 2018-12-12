/* File>>> window.c
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
#include <unistd.h>
#include <string.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>

/*------------------------------------------------------------------------------
-- Local include files
------------------------------------------------------------------------------*/

#include "xviglocal.h"

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
------------------------------------------------------------------------------*/

int XviG_OpenWindow(char *name,
                    int x,
                    int y,
                    unsigned int *width,
                    unsigned int *height)
{
  char *geometry;
  int bitmask;
  int rm_x, rm_y;
  unsigned int rm_width, rm_height;
  char window_id[32], init_x[32], init_y[32], init_width[32], init_height[32];
  int pid;

  if (XviG_cursor_mode)
  {
    fprintf(stderr, "ERROR : Window not opened, close the cursor first.\n");
    return 0;
  }

  if (Get_WinItem(name))
  {
    fprintf(stderr, "ERROR : Window '%s' is already open.\n", name);
    return 0;
  }

/* rbh 09/14/2009 
  if (geometry = XviG_GetRMGeometry(name))
  gcc -I ../include -Wall
	window.c:82: warning: suggest parentheses around assignment used as truth value
*/
  if ( (geometry = XviG_GetRMGeometry(name)) )
  {
    bitmask = XParseGeometry(geometry, &rm_x, &rm_y, &rm_width, &rm_height);

    if (bitmask & WidthValue)
      *width = rm_width;

    if (bitmask & HeightValue)
      *height = rm_height;

/* 09/14/2009 
window.c:97: warning: suggest explicit braces to avoid ambiguous ‘else’
window.c:104: warning: suggest explicit braces to avoid ambiguous ‘else’
*/
    if (bitmask & XValue){
      if (bitmask & XNegative){
        x = DisplayWidth(XviG_display, XviG_screen_nr)
            + rm_x - *width - 2*BORDER_WIDTH;
      } else {
        x = rm_x;
      }
    }

    if (bitmask & YValue) {
      if (bitmask & YNegative){
        y = DisplayHeight(XviG_display, XviG_screen_nr)
            + rm_y - *height - 2*BORDER_WIDTH;
      } else {
        y = rm_y;
      }
    }
  }

  sprintf(window_id, "%d", (int)XviG_dummy_window);
  sprintf(init_x, "%d", x);
  sprintf(init_y, "%d", y);
  sprintf(init_width, "%d", *width);
  sprintf(init_height, "%d", *height);

  /*
  -- Start the subprocess to create the window
  */

  pid = fork();

  switch (pid)
  {
    case 0 : /* Child process */
             if (execl(XviG_pathname, "calxvig", XviG_VERSION, name, window_id,
                       init_x, init_y, init_width, init_height, NULL) == -1)
             {
               fprintf(stderr, "ERROR : Cannot run execl.\n");
               exit(1);
             }
    case -1 : /* Error */
              fprintf(stderr, "ERROR : Cannot fork.\n");
              return 0;
  }

  /*
  -- Ask for the ClientEvent to get the real window ID, the pixmap ID
  -- and the window size
  */

  while (1)
  {
    XNextEvent(XviG_display, &XviG_event);

    if (XviG_event.type == ClientMessage)
    {
      if (XviG_event.xclient.message_type == XviG_MESSAGE_INIT)
      {
        XviG_window = (Window) XviG_event.xclient.window;
        XviG_pixmap = (Pixmap) XviG_event.xclient.data.l[0];
        XviG_window_width = *width
                          = (unsigned int) XviG_event.xclient.data.l[1];
        XviG_window_height = *height
                           = (unsigned int) XviG_event.xclient.data.l[2];

        break;
      }
      else
        printf("WARNING : Wrong ClientMessage received .....\n");
    }
    /*
    else
      printf("WARNING : Other event received .....\n");
    */
  }

  /*
  -- Store the window info in the window list
  */

  New_WinItem(name, pid);

  return 1;
}

/*------------------------------------------------------------------------------
--
--
--
------------------------------------------------------------------------------*/

int XviG_CloseWindow(char *name)
{
  WinItem curwin;

  if (XviG_cursor_mode)
  {
    fprintf(stderr, "ERROR : Window not closed, close the cursor first.\n");
    return 0;
  }

  if (!(curwin = Get_WinItem(name)))
  {
    fprintf(stderr, "ERROR : Window '%s' does not exist.\n", name);
    return 0;
  }

  XviG_event.xclient.message_type = XviG_MESSAGE_QUIT;
  XviG_event.xclient.format = 8;
  XviG_event.type = ClientMessage;

  if (!XSendEvent(XviG_display, curwin->window,
                  False, NoEventMask, &XviG_event))
  {
    fprintf(stderr, "ERROR : Cannot send ClientMessage 'quit'.\n");
    return 0;
  }

  XFlush(XviG_display);

  Delete_WinItem(curwin);

  return 1;
}

/*------------------------------------------------------------------------------
--
--
--
------------------------------------------------------------------------------*/

int XviG_SelectWindow(char *name)
{
  WinItem curwin;

  if (XviG_cursor_mode)
  {
    fprintf(stderr, "ERROR : Window not selected, close the cursor first.\n");
    return 0;
  }

  if (!(curwin = Get_WinItem(name)))
  {
    fprintf(stderr, "ERROR : Window '%s' does not exist.\n", name);
    return 0;
  }

  XviG_window = curwin->window;
  XviG_pixmap = curwin->pixmap;
  XviG_window_width = curwin->width;
  XviG_window_height = curwin->height;

  selected_winitem = curwin;

  return 1;
}

/*------------------------------------------------------------------------------
--
--
--
------------------------------------------------------------------------------*/

void XviG_WindowSize(unsigned int *width,
                     unsigned int *height)
{
  XviG_event.xclient.message_type = XviG_MESSAGE_SIZE;
  XviG_event.xclient.format = 8;
  XviG_event.type = ClientMessage;

  if (!XSendEvent(XviG_display, XviG_window, False, NoEventMask, &XviG_event))
  {
    fprintf(stderr, "ERROR : Cannot send ClientMessage 'size'.\n");

    *width = XviG_window_width;
    *height = XviG_window_height;

    return;
  }

  while (1)
  {
    XNextEvent(XviG_display, &XviG_event);

    if (XviG_event.type == ClientMessage)
    {
      if (XviG_event.xclient.message_type == XviG_MESSAGE_SIZE)
      {
        XviG_window_width = *width
                          = (unsigned int) XviG_event.xclient.data.l[0];
        XviG_window_height = *height
                           = (unsigned int) XviG_event.xclient.data.l[1];
        XviG_pixmap = (Pixmap) XviG_event.xclient.data.l[2];

        /*
        -- Update the currently selected window item in the window list
        */

        selected_winitem->width = XviG_window_width;
        selected_winitem->height = XviG_window_height;
        selected_winitem->pixmap = XviG_pixmap;

        break;
      }
      else
        printf("WARNING : Wrong ClientMessage received .....\n");
    }
    /*
    else
      printf("WARNING : Other event received .....\n");
    */
  }
}

/*------------------------------------------------------------------------------
--
--
--
------------------------------------------------------------------------------*/

void XviG_WindowPosition(int *x,
                         int *y)
{
  Window root_rtn;
  unsigned int width_rtn, height_rtn, bwidth_rtn, depth_rtn;

  if (!XGetGeometry(XviG_display, XviG_window, &root_rtn, x, y,
                    &width_rtn, &height_rtn, &bwidth_rtn, &depth_rtn))
  {
    printf("WARNING : Cannot get position of window.\n");
    *x = 0;
    *y = 0;
  }
}

/*------------------------------------------------------------------------------
--
--
--
------------------------------------------------------------------------------*/

void XviG_ClearWindow(void)
{
  XClearWindow(XviG_display, XviG_window);

  XSetForeground(XviG_display, XviG_gc,
                 BlackPixel(XviG_display, XviG_screen_nr));
  /* XSetFillStyle(XviG_display, XviG_gc, FillSolid); This is the default */
  XFillRectangle(XviG_display, XviG_pixmap, XviG_gc,
                 0, 0, XviG_window_width, XviG_window_height);

  /*
  -- Set the original color back in the graphical context
  */

  XSetForeground(XviG_display, XviG_gc, XviG_save_color);
}
