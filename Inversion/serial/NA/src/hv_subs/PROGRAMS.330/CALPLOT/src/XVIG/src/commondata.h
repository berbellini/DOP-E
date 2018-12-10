/* File>>> commondata.h
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
-- Function : Definition of the common data of the user program (XviG library)
--            and the XviG child program.
--
-- Comment  :
--
-- Review   :
--
*/


#ifndef __COMMONDATA_H
#define __COMMONDATA_H


/*
-- The version string
*/

#define XviG_VERSION  "Version_1.1.0"


/*
-- The border width for the windows
*/

#define BORDER_WIDTH  4


/*
-- The Atom strings for the ClientMessage events
*/

#define MESSAGE_INIT          "_XVIG_MESSAGE_INIT"
#define MESSAGE_KEY           "_XVIG_MESSAGE_KEY"
#define MESSAGE_BUTTON        "_XVIG_MESSAGE_BUTTON"
#define MESSAGE_KEY_BUTTON    "_XVIG_MESSAGE_KEY_BUTTON"
#define MESSAGE_SIZE          "_XVIG_MESSAGE_SIZE"
#define MESSAGE_SENSEKBD_ON   "_XVIG_MESSAGE_SENSEKBD_ON"
#define MESSAGE_SENSEKBD_OFF  "_XVIG_MESSAGE_SENSEKBD_OFF"
#define MESSAGE_SENSEKBD      "_XVIG_MESSAGE_SENSEKBD"
#define MESSAGE_CURSOR        "_XVIG_MESSAGE_CURSOR"
#define MESSAGE_QUIT          "_XVIG_MESSAGE_QUIT"

/* extension rbh */
#define MESSAGE_CLIP          "_XVIG_MESSAGE_CLIP"
#define MESSAGE_REVERSE       "_XVIG_MESSAGE_REVERSE"
#define MESSAGE_TOGRAY        "_XVIG_MESSAGE_TOGRAY"
#define MESSAGE_BOUNDS        "_XVIG_MESSAGE_BOUNDS"
#define MESSAGE_LCMAP_SIZE    "_XviG_MESSAGE_LCMAP_SIZE"
#define MESSAGE_LCMAP_ENTRY   "_XviG_MESSAGE_LCMAP_ENTRY"
#define MESSAGE_LCMAP_RESET   "_XviG_MESSAGE_LCMAP_RESET"
#define MESSAGE_ORIGIN_HYPERBOLA "_XviG_MESSAGE_ORIGIN_HYPERBOLA"



/*
-- The data for the ClientMessage events
*/

#define DATA_SENSEKBD_YES  ((char) 1)
#define DATA_SENSEKBD_NO   ((char) 0)
#define DATA_CURSOR_ARROW  	0L
#define DATA_CURSOR_XORARROW  	1L
#define DATA_CURSOR_XHAIR  	2L
#define DATA_CURSOR_PLUS   	3L
#define DATA_CURSOR_BOX   	4L
#define DATA_CURSOR_RUBBER   	5L
#define DATA_CURSOR_OFF   	6L
#define DATA_CURSOR_HYPERBOLA  	8L


#endif  /* __COMMONDATA_H */
