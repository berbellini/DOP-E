/* File>>> snskbd.c
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


/*------------------------------------------------------------------------------
-- Local function declarations
------------------------------------------------------------------------------*/


/*------------------------------------------------------------------------------
--
--
--
------------------------------------------------------------------------------*/

void XviG_SetSenseKbd(int sense_char)
{
  /*
  -- Send the ClientMessage to set the `sense keyboard' functionality on/off
  */

  if (sense_char)
  {
    XviG_event.xclient.message_type = XviG_MESSAGE_SENSEKBD_ON;
    XviG_event.xclient.data.b[0] = (char) sense_char;
  }
  else
    XviG_event.xclient.message_type = XviG_MESSAGE_SENSEKBD_OFF;

  XviG_event.xclient.format = 8;
  XviG_event.type = ClientMessage;

  if (!XSendEvent(XviG_display, XviG_window, False, NoEventMask, &XviG_event))
    fprintf(stderr, "ERROR : Cannot send ClientMessage 'sense_kbd_on/off'.\n");
}

/*------------------------------------------------------------------------------
--
--
--
------------------------------------------------------------------------------*/

int XviG_SenseKbd(void)
{
  int sense_kbd;

  /*
  -- First send the ClientMessage to ask for a `sense keyboard'
  */

  XviG_event.xclient.message_type = XviG_MESSAGE_SENSEKBD;
  XviG_event.xclient.format = 8;
  XviG_event.type = ClientMessage;

  if (!XSendEvent(XviG_display, XviG_window, False, NoEventMask, &XviG_event))
  {
    fprintf(stderr, "ERROR : Cannot send ClientMessage 'sense_kbd'.\n");
    return 0;
  }

  /*
  -- Now ask for the `sense keyboard' event
  */

  while (1)
  {
    XNextEvent(XviG_display, &XviG_event);

    if (XviG_event.type == ClientMessage)
    {
      if (XviG_event.xclient.message_type == XviG_MESSAGE_SENSEKBD)
      {
        sense_kbd = XviG_event.xclient.data.b[0];
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

  return sense_kbd;
}
