/*===============================================================================
$Id: ga.h,v 1.1 2008/02/05 17:17:32 jkleinj Exp $
ga : Genetic Algorithm
(C) 2003-2007 Jens Kleinjung
(C) 2006-2007 Alessandro Pandini

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
*==============================================================================*/

#if !defined(GA_H)
#define GA_H

/*____________________________________________________________________________*/
/* includes */
#include <alloca.h>
#include <assert.h>
#include <ctype.h>
#include <getopt.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/*____________________________________________________________________________*/
/* define debug mode and debug 'dump' function */
/*#define DEBUG*/
#define dump(x, fmt) fprintf(stderr, "dump@%s:%u: %s=" fmt"\n", __FILE__, __LINE__, #x, x);
#define dump2(x1, fmt1, x2, fmt2) fprintf(stderr, "dump2@%s:%u: %s=" fmt1"\t%s=" fmt2"\n", __FILE__, __LINE__, #x1, x1, #x2, x2);

/*____________________________________________________________________________*/
/* structures */

/* bit field for binary gene */
#ifdef BIT
typedef struct
{
	unsigned int state: 1; /* 1 bit integer */
} Bitgene;
#endif

/* genome pool */
typedef struct
{
#ifdef BIT
	Bitgene *bitgenome; /* bitgenome is an array of binary parameters=bitgenes */
#else
    int *genome; /* genome is an array of parameters=genes */
#endif
    float fitness; /* fitness of genome */
} Pool;

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/
/* basic GA parameters */
typedef struct
{
	/* GA scopes */
	int popsize; /* population size = genome number */
	int fitmate; /* number of fittest genomes to mate */
	int genenum; /* number of genes to optimize */
	int generation; /* number of optimisation cycles */
	int lowlim; /* lower limit of gene value */
	int uplim; /* upper limit of gene value */

	/* fitness sorting */
	int minimize; /* lowest fitness is top of list */
	int maximize; /* highest fitness is top of list */

	/* initial pool */
	int random; /* random values for genes of initial pool */
	int seeded; /* seeded values for genes of initial pool */
	int maxvar; /* maximal variation of seeded value */

	/* breeding modes: single choice mandatory */
	int crossover; /* breeding mode: crossover */
	int equilibrium;  /* breeding mode: equilibrium */

	/* repeats */
	int jackknife; /* splitting of database in 'JACKKNIFE' parts (1 = no jackknife) */
	int repeat; /* repeat entire GA 'REPEAT' times (1 = no repeat) */
} Gapar;

/*____________________________________________________________________________*/
/* prototypes */
FILE *safe_open(const char *name, const char *mode);
extern void *safe_malloc(size_t), *safe_realloc(void *, size_t);
int mem_fitness(Pool *pool, Gapar *gapar, int ix);
void constrain_genome(Pool *pool, Gapar *gapar, int ix, int maxgenes);
void set_outfilename(char *outfilename, char c0, char c2, char c4, char c6);

#endif
