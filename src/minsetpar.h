/*==============================================================================
$Id: minsetpar.h,v 1.1 2008/02/05 17:17:32 jkleinj Exp $
minsetpar.h : MinSet default parameters
(C) 2006-2007 Jens Kleinjung and Alessandro Pandini

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

#if !defined(MINSETPAR_H)
#define MINSETPAR_H

/* stuctural code */
#define SA_ALPHABET
/* define generation of background distribution for most fitted */
/*#define GENERATE_BACKGROUND*/

/* number of randomly generated pools for background distribution */
#define NUMRAN_POOL 2500 

/* subset parameters */
#define BASESET "masterfilelist" /* base set of proteins: list of sequence file names */
#define SEQDIR "fastas/" /* relative path to directory containing sequence files */
#define ALPHABET "TOP2006" /* coding alphabet */
#define SUBSETSIZE 20. /* target size of subset relative to base set size (in % units) */
#define KWORDLENGTH 2 /* selection k-word (fragment) length */

#endif

