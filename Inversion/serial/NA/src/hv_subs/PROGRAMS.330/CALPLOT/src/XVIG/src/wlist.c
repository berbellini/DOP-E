/* File>>> wlist.c
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
#include <X11/Xlib.h>

/*------------------------------------------------------------------------------
-- Local include files
------------------------------------------------------------------------------*/

#include "xviglocal.h"

/*------------------------------------------------------------------------------
-- Static variable declarations
------------------------------------------------------------------------------*/

static WinItem startWindow = (WinItem) NULL,
               stopWindow = (WinItem) NULL;

/*------------------------------------------------------------------------------
-- Local function declarations
------------------------------------------------------------------------------*/


/*------------------------------------------------------------------------------
--
--
--
------------------------------------------------------------------------------*/

void New_WinItem(char *name,
                 int pid)
{
  WinItem curwin;

  curwin = (WinItem) malloc(sizeof(aWinItem));

  curwin->name = STRNEW(name);
  curwin->pid = pid;
  curwin->window = XviG_window;
  curwin->pixmap = XviG_pixmap;
  curwin->width = XviG_window_width;
  curwin->height = XviG_window_height;
  curwin->next = (WinItem) NULL;

  if (!startWindow)
  {
    curwin->prev = (WinItem) NULL;
    startWindow = stopWindow = curwin;
  }
  else
  {
    curwin->prev = stopWindow;
    stopWindow->next = curwin;
    stopWindow = curwin;
  }

  selected_winitem = curwin;
}

/*------------------------------------------------------------------------------
--
--
--
------------------------------------------------------------------------------*/

void Delete_WinItem(WinItem curwin)
{
  if (curwin == startWindow)
  {
    if (curwin == stopWindow)
      startWindow = (WinItem) NULL;
    else
    {
      curwin->next->prev = (WinItem) NULL;
      startWindow = curwin->next;
    }
  }
  else
  if (curwin == stopWindow)
  {
    curwin->prev->next = (WinItem) NULL;
    stopWindow = curwin->prev;
  }
  else
  {
    curwin->prev->next = curwin->next;
    curwin->next->prev = curwin->prev;
  }

  /*
  -- If it is the currently selected window, change all the global IDS
  */

  if ((curwin->window == XviG_window) && (curwin->prev || curwin->next))
  {
    WinItem newwin;

    newwin = curwin->prev ? curwin->prev : curwin->next;

    XviG_window = newwin->window;
    XviG_pixmap = newwin->pixmap;
    XviG_window_width = newwin->width;
    XviG_window_height = newwin->height;

    selected_winitem = newwin;
  }

  /*
  -- Free it
  */

  free(curwin->name);
  free(curwin);
}

/*------------------------------------------------------------------------------
--
--
--
------------------------------------------------------------------------------*/

WinItem Get_WinItem(char *name)
{
  WinItem loopwin;

  for (loopwin = startWindow; loopwin; loopwin = loopwin->next)
  {
    if (!strcmp(loopwin->name, name))
      return loopwin;
  }

  return (WinItem) NULL;
}

/*------------------------------------------------------------------------------
--
--
--
------------------------------------------------------------------------------*/

void XviG_CleanupWindows(void)
{
  WinItem loopwin, nextwin;

  for (loopwin = startWindow; loopwin; loopwin = nextwin)
  {
    nextwin = loopwin->next;
    (void) XviG_CloseWindow(loopwin->name);
  }
}
