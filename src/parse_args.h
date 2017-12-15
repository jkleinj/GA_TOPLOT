/*==============================================================================
$Id: parse_args.h,v 1.1 2008/02/05 17:17:32 jkleinj Exp $
parse_args.h : parse command line arguments
(C) 2006-2007 Alsseandro Pandini and Jens Kleinjung

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.
     This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
==============================================================================*/

#if !defined(PARSE_ARGS_H)
#define PARSE_ARGS_H

#include "ga.h"
#include "minset.h"

void license( void );
void usage( void );
void parse_args(int argc, char **argv, Gapar *gapar, Minset *ms, FILE *outfile);

#endif
