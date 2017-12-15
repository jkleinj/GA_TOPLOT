/*==============================================================================
minset.h : program for database subsetting
(C) 2006-2017 Jens Kleinjung
(C) 2006-2007 Alesseandro Pandini

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

#if !defined(MINSET_H)
#define MINSET_H

/*____________________________________________________________________________*/
/* includes */
#include "alphabet.h"
#include "ga.h"

/*____________________________________________________________________________*/
/* defines */

/* max macro */
#define max(a,b)  (((a) > (b)) ? (a) : (b))

/*___________________________________________________________________________*/
typedef struct
{
    char *name; /* protein (file)name */
    char *description; /* description in header of fastafile */
    char *seq; /* seq from fastafile */
    float entropy; /* entropy */
    float score; /* score */
} ProteinEntry;

/*___________________________________________________________________________*/
typedef struct
{
    ProteinEntry *protein; /* data for each protein */
    int n_prot; /*number of proteins */
    float entropy_sum; /* sum of single protein entropy in aa code */
} Prots;

/*____________________________________________________________________________*/
typedef struct 
{
	char basesetFileName[200]; /* input file: list of protein FASTA filenames */
	char seqdir[200]; /* relative path to directory holding sequence files */
	Alphabet alphabet; /* code alphabets */
	float subsetsize; /* target percentage of minset/baseset size */
	int kword_len; /* selection k-word (substring) length */

    /*____________________________________________________________________________*/
	Alphabet bg_freq;

    /*____________________________________________________________________________*/
	Prots prots; /* list of proteins */

    /*____________________________________________________________________________*/
	char *setfasta; /* string of all (concatenated) sequences of base set */
	int *setfasta_charCount; /* array of counts of single-character code symbols */
	char *polyfasta; /* string of subset of (concatenated) sequences */
	int *polyfasta_charCount; /* array of counts of single-character code symbols */

	float kl_distance; /* Kullback-Leibler distance */
	int total_len; /* total string length of concatenated sequences */
	int n_selected; /* number of selected proteins */

	char *subsetfasta; /* string of all (concatenated) sequences of subset */
    FILE *subsetOutFile; char *subsetOutFileName; /* file containing subset information */
} Minset;

/*____________________________________________________________________________*/
void initialise_minset(Pool *pool, Gapar *gapar, Minset *ms);
float run_minset(Pool *pool, Gapar *gapar, Minset *ms, int j, int k, int l, int ix);
void finalise_minset(Minset *ms);
int read_sequence(FILE *aafile, Prots *prots, int k);
void parametrise_minset(Minset *ms);
void print_subset(Pool *pool, Gapar *gaPar, Minset *ms);

#endif
