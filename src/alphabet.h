/*==============================================================================
$Id: alphabet.h,v 1.1 2008/02/05 17:17:32 jkleinj Exp $
alphabet.c : alphabet assignment routine
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

#if !defined(ALPHABET_H)
#define ALPHABET_H

/*___________________________________________________________________________*/
/* structure for alphabet data */
typedef struct {
    char name[200]; /* alphabet name */
    char *codeOrder; /* string defining order of code characters */
    float *freq; /* frequency of code characters */
    int size; /* code array size (spannning all character ASCII values) */
    int codeLength; /* number of code characters */
} Alphabet;

/*____________________________________________________________________________*/
/* prototypes */
void print_alphabet(Alphabet *alphabet);
void set_alphabet(Alphabet *alphabet);

#endif
