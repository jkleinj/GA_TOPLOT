/*===============================================================================
ga : Genetic Algorithm
(C) 2003-2017 Jens Kleinjung
(C) 2006-2007 Alessandro Pandini

Usage: The MANUAL gives a brief description of the program.
	Sections for user modification are indicated by >>>>> delimiters.
	Consult the MANUAL how to add an application to the GA program.

References:
	Users publishing results obtained with the GA program and its applications 
	should acknowledge its use by one of the following citations:

	- J. Kleinjung, J. Romein, K. Lin and J. Heringa
	Contact-based sequence alignment.
	Nucleic Acids Res. 32(8) (2004) 2464-2473. 

    - A. Pandini, L. Bonati, F. Fraternali and J. Kleinjung
    MinSet: A general approach to derive maximally representative database
    subsets by using fragment dictionaries and its application to the SCOP
    database.
	Bioinformatics 23 (2007) 515-516.

    - A. Pandini and J. Kleinjung
    MinSet: automatic derivation of maximally representative subsets of
    structural domains.
	Nucleic Acids Res. (2007) submitted.

License:
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

/*____________________________________________________________________________*/
/* includes */
#include "ga.h"
#include "gapar.h"
#include "minset.h"
#include "parse_args.h"

/*____________________________________________________________________________*/
/* safe allocation */
static void *check_non_null(void *ptr)
{
    assert (ptr != 0);
    return ptr;
}

void *safe_malloc(size_t size)
{
    return check_non_null(malloc(size));
}

void *safe_realloc(void *ptr, size_t size)
{
    return check_non_null(realloc(ptr, size));
}

/*____________________________________________________________________________*/
/* open file */
FILE *safe_open(const char *name, const char *mode)
{
    FILE *file = fopen(name, mode);

    assert (file != 0);
    return file;
}

/*____________________________________________________________________________*/
/* generate random integer number within range  0 to 'gaPar.uplim' */
int get_rand(int limit)
{
    double dranum; /* random number in double format */
    int iranum; /* random number in integer format */
    dranum = rand() / (double)INT_MAX; /* generate rand and scale to range 0-1 */
    iranum = (int)((limit + 1) * dranum); /* scale to integer of range 0-limit */
    if (iranum == (limit + 1)) -- iranum; /* exclude upper limit */
    return iranum;
}

/*____________________________________________________________________________*/
/* generate random integer number within range 'low' to 'up' */
int get_rand_lu(unsigned int low, unsigned int up)
{
    assert(low < up);

    double dranum; /* random number in double format */
    int iranum; /* random number in integer format */
    dranum = rand() / (double)INT_MAX; /* generate rand and scale to range 0-1 */
    iranum = (int)(((up - low + 1) * dranum) + low); /* scale to range [low,up] */
    if (iranum == (up + 1)) -- iranum; /* exclude upper limit */
    return iranum;
}

/*____________________________________________________________________________*/
/* fitness comparison routine used by 'qsort'  */
int cmp_fitness(const void *pp1, const void *pp2)
{
    Pool *p1 = (Pool*)pp1;
    Pool *p2 = (Pool*)pp2;

#ifdef MINIMIZE
    return (p1->fitness > p2->fitness ? -1 : 1); /* MINIMIZE */
#else
    return (p1->fitness < p2->fitness ? -1 : 1);  /* MAXIMIZE */
#endif
}

/*____________________________________________________________________________*/
/* sort pool by fitness */
void sort_fitness(Pool *pool, Gapar *gaPar)
{
    qsort(&pool[0], gaPar->popsize, sizeof(Pool), cmp_fitness);
}

/*____________________________________________________________________________*/
/* bitfield genome: set gene state to either 0 or 1, depending on passed integer 'i' */
#ifdef BIT
void set_bitgene(Bitgene *bitgen, int i) 
{
	bitgen->state = i > 0 ? 1 : 0;
}
#endif

/*____________________________________________________________________________*/
/* bitfield genome: return integer 0 or 1, depending on gene bit 'state' */
#ifdef BIT
int get_bitgene(Bitgene *bitgen) 
{
	return (int)bitgen->state > 0 ? 1 : 0;
}
#endif

/*____________________________________________________________________________*/
/* vary value */
int vary_value(int val, Gapar *gaPar)
{
    int newval; /* new value */
    int mp; /* minus/plus */

    mp = get_rand(2);

    if (mp == 1)
        while ((newval = val - get_rand(gaPar->maxvar)) < gaPar->lowlim);
    else
        while ((newval = val + get_rand(gaPar->maxvar)) > gaPar->uplim);
    return newval;
}

/*____________________________________________________________________________*/
/* initialise random pool and fitness*/
static void init_random_pool(Pool *pool, Gapar *gaPar)
{
    unsigned int i, j;

    for (i = 0; i < gaPar->popsize; ++ i) /* genomes */
    {
        for (j = 0; j < gaPar->genenum; ++ j) /* genes */
		{
#ifdef BIT
            set_bitgene(&pool[i].bitgenome[j], get_rand(gaPar->uplim)); /* random value */
#else
            pool[i].genome[j] = get_rand(gaPar->uplim); /* random value */
#endif
		}

        pool[i].fitness = 0; /* initialise fitness */
    }
}

/*____________________________________________________________________________*/
/* initialise seeded pool and fitness */
void init_seeded_pool(Pool *pool, Gapar *gaPar)
{
    unsigned int i, j;

    for (i = 0; i < gaPar->popsize; ++ i) /* genomes */
    {
		/* assign defined value to each gene */
		/*
		pool[i].genome[0] = 4;
		*/

		pool[i].fitness = 0; /* initialise fitness */
    }
    
    /* vary values, not first genome */
    for (i = 1; i < gaPar->popsize; ++ i)
    {
		for (j = 0; j < gaPar->genenum; ++ j)
		{
#ifdef BIT
			set_bitgene(&pool[i].bitgenome[j], vary_value(get_bitgene(&pool[i].bitgenome[j]), gaPar));
#else
			pool[i].genome[j] = vary_value(pool[i].genome[j], gaPar);
#endif
		}
    }
}

/*____________________________________________________________________________*/
/* equilibrium: calc average parameter values in fittest genomes */
void equilibrium(Pool *pool, Gapar *gaPar, int *average)
{
    unsigned int i, j;

    for (j = 0; j < gaPar->genenum; ++ j)
    {
		average[j] = 0;
		for (i = 0; i < gaPar->fitmate; ++ i) /* sum up over the fittest genomes */
		{
#ifdef BIT
			average[j] += get_bitgene(&pool[i].bitgenome[j]);
#else
			average[j] += pool[i].genome[j];
#endif
		}
		average[j]  = (int)(average[j]/gaPar->fitmate); /* normalize */
    }
}

/*____________________________________________________________________________*/
/* breed_equilibrium: generate equilibrium genomes */ 
void breed_equilibrium(Pool *pool, Gapar *gaPar, int *average, int ix)
{
    unsigned int j;

    /*maxvar = gaPar.uplim / gaPar.generation;*/ /* decay to zero variation */

    for (j = 0; j < gaPar->genenum; ++ j)
	{
#ifdef BIT
		set_bitgene(&pool[ix].bitgenome[j], vary_value(get_bitgene(&pool[ix].bitgenome[j]), gaPar->maxvar));
#else
pool[ix].genome[j] = vary_value(pool[ix].genome[j], gaPar);
#endif
	}
    pool[ix].fitness = 0; /* initialise fitness */
}

/*____________________________________________________________________________*/
/* calculate sum of genome: n of selected items */
int gene_sum(Pool *pool, Gapar *gaPar, int ix)
{
	unsigned int gsum = 0;
    unsigned int i;

    for(i = 0; i < gaPar->genenum; ++ i)
	{
#ifdef BIT
		gsum += get_bitgene(&pool[ix].bitgenome[i]);
#else
        gsum += pool[ix].genome[i];
#endif
	}

	return gsum;
}

/*____________________________________________________________________________*/
/* constrain individual genome to specified maximal gene number */
void constrain_genome(Pool *pool, Gapar *gaPar, int ix, int maxgenes)
{
    unsigned int j;

	 /* check that total number of selected genes does not exceed 'maxgenes' */
    while(gene_sum(pool, gaPar, ix) > maxgenes)
	{
        /* deselect a random protein*/
#ifdef BIT
        do {
			j = (int) (rand()/(double)INT_MAX * gaPar->genenum);
        } while(get_bitgene(&pool[ix].bitgenome[j]) == 0);
        set_bitgene(&pool[ix].bitgenome[j], 0); 
#else
        do {
			j = (int) (rand()/(double)INT_MAX * gaPar->genenum);
        } while(pool[ix].genome[j] == 0);
        pool[ix].genome[j] = 0; 
#endif
    }
}

/*____________________________________________________________________________*/
/* breed_crossover: crossover between gaPar.fitmate pairs and create genome 'ix' */ 
void breed_crossover(Pool *pool, Gapar *gaPar, int iy)
{
    int ia, ib, j, yn;

    ia = -1;
    ib = -1;
    for (j = 0; j < gaPar->genenum; ++ j)
    {
		ia = get_rand(gaPar->fitmate); /* principal parent (1)*/

		do
			ib = get_rand(gaPar->fitmate); /* crossover parent (2) */
		while (ib == ia);

		yn = get_rand(2); /* choose crossover yes/no */
#ifdef BIT
		if (yn == 1)
			set_bitgene(&pool[iy].bitgenome[j], get_bitgene(&pool[ib].bitgenome[j]));
		else
			set_bitgene(&pool[iy].bitgenome[j], get_bitgene(&pool[ia].bitgenome[j]));
#else
		if (yn == 1)
			pool[iy].genome[j] = pool[ib].genome[j];
		else
			pool[iy].genome[j] = pool[ia].genome[j];
#endif
    }

    pool[iy].fitness = 0; /* initialise fitness */
}

/*____________________________________________________________________________*/
/* mem_fitness: get fitness from identical genome that has already been calculated */
int mem_fitness(Pool *pool, Gapar *gaPar, int ix)
{
    unsigned int i, j;
    unsigned int id;

    for (i = 0; i < ix; ++ i)
    {
		id = 1;

		/* test whether any gene of 'i' and 'ix' differs */
		for (j = 0; j < gaPar->genenum; ++ j)
		{
#ifdef BIT
			if (get_bitgene(&pool[i].bitgenome[j]) != get_bitgene(&pool[ix].bitgenome[j]))
				id = 0;
#else
			if (pool[i].genome[j] != pool[ix].genome[j])
				id = 0;
#endif
		}
		/* if identical and already fitnessed: copy fitness */
		if ((id == 1) && (pool[i].fitness > 0))
		{
			pool[ix].fitness = pool[i].fitness;
			return 1;
		} 
    }

    return 0;
}

/*____________________________________________________________________________*/
/* check_convergence: all genomes identical */
int check_convergence(Pool *pool, Gapar *gaPar)
{
    unsigned int i, j;

    for (i = 1; i < gaPar->popsize; ++ i)
		for (j = 0; j < gaPar->genenum; ++ j)
			if (pool[i].genome[j] != pool[i-1].genome[j])
				return 1;

    return 0;
}

/*____________________________________________________________________________*/
/* print pool in binary format */
void print_pool_bin(Pool *pool, Gapar *gaPar, FILE *outfile)
{
    unsigned int i, j;
    char ci;

    for (i = 0; i < gaPar->popsize; ++ i)
    {
        fwrite(&i, sizeof(int), 1, outfile);

        for (j = 0; j < gaPar->genenum; ++ j)
        {
#ifdef BIT
            ci = get_bitgene(&pool[i].bitgenome[j]) ? '1':'0';
#else
            ci = pool[i].genome[j] ? '1':'0';
#endif
            fwrite(&ci, sizeof(char), 1, outfile);
        }

        fwrite(&pool[i].fitness, sizeof(float), 1, outfile);

    }
    fflush(outfile);
}

/*____________________________________________________________________________*/
/* print pool in ASCII format */
/* argument 'fit' is a switch to print only the fittest genomes */
void print_pool_ascii(Pool *pool, Gapar *gaPar, FILE *outfile, int l, int fit)
{
    unsigned int i, j;
	unsigned int end;

    fprintf(outfile, "#%d\n", l);

	if (fit != 0)
		end = gaPar->fitmate;
	else
		end = 0;

    for (i = 0; i < end; ++ i)
    {
		fprintf(outfile, "%3d: ", i);

		for (j = 0; j < gaPar->genenum; ++ j)
		{
#ifdef BIT
			fprintf(outfile, "%1d ", get_bitgene(&pool[i].bitgenome[j]));
#else
			fprintf(outfile, "%1d ", pool[i].genome[j]);
#endif
		}

		fprintf(outfile, "%6.4f\n", pool[i].fitness);

    }
    fprintf(outfile, "\n");
    fflush(outfile);
}

/*____________________________________________________________________________*/
/* set output filename */
void set_outfilename(char *outfilename, char c0, char c2, char c4, char c6)
{
	outfilename[0] = c0;
	outfilename[2] = c2;
	outfilename[4] = c4;
	outfilename[6] = c6;
}

/*____________________________________________________________________________*/
/* parametrise GA */
void parametrise_ga(Gapar *gaPar)
{
	/* scopes */
	gaPar->popsize = (int)POPSIZE; assert(gaPar->popsize > 1);
	gaPar->fitmate = (int)FITMATE; assert(gaPar->fitmate < gaPar->popsize);
	gaPar->genenum = (int)GENENUM; assert(gaPar->genenum > 0);
	gaPar->generation = (int)GENERATION; assert(gaPar->generation > 0);
	gaPar->lowlim = (int)LOWLIM;
	gaPar->uplim = (int)UPLIM;
	assert(gaPar->uplim > gaPar->lowlim);
#ifdef BIT
	assert(gaPar->lowlim == 0 && gaPar->uplim == 1);
#endif

	/* fitness sorting */
	gaPar->minimize = (int)MINIMIZE;
	gaPar->maximize = (int)MAXIMIZE;
	assert((gaPar->minimize == 0 && gaPar->maximize == 1) || \
			(gaPar->minimize == 1 && gaPar->maximize == 0));

	/* pool initialisation */
	gaPar->random = (int)RANDOM;
	gaPar->seeded = (int)SEEDED;
	assert((gaPar->random == 0 && gaPar->seeded == 1) || \
			(gaPar->random == 1 && gaPar->seeded == 0));

	/* breedmode */
	gaPar->crossover = (int)CROSSOVER;
	gaPar->equilibrium = (int)EQUILIBRIUM;
	gaPar->maxvar = (int)MAXVAR;
	assert((gaPar->equilibrium == 0 && gaPar->crossover == 1) || \
			(gaPar->equilibrium == 1 && gaPar->crossover == 0));

	/* split database and/or repeat GA */
	gaPar->jackknife = (int)JACKKNIFE; assert(gaPar->jackknife > 0);
	gaPar->repeat = (int)REPEAT; assert(gaPar->repeat > 0);
}

/*____________________________________________________________________________*/
/* parse GA command line long_options */
void parse_args_ga(int argc, char **argv, Gapar *gaPar)
{
	int c;

	static struct option long_options[] =
	{
		{"popsize", required_argument, 0, 1},
		{"fitmate", required_argument, 0, 2},
		{"genenum", required_argument, 0, 3},
		{"generation", required_argument, 0, 4},
		{"lowlim", required_argument, 0, 5},
		{"uplim", required_argument, 0, 6},
		{"minimize", required_argument, 0, 7},
		{"maximize", required_argument, 0, 8},
		{"random", required_argument, 0, 9},
		{"seeded", required_argument, 0, 10},
		{"maxvar", required_argument, 0, 11},
		{"crossover", required_argument, 0, 12},
		{"equilibrium", required_argument, 0, 13},
		{"jackknife", required_argument, 0, 14},
		{"repeat", required_argument, 0, 15},
		{0, 0, 0, 0}
	};

	while ((c = getopt_long (argc, argv, "1:2:3:4:5:6:7:8:9:10:11:12:13:14:15:", long_options, NULL)) != -1)
	{
		switch(c)
		{
			case 1:
				gaPar->popsize = atoi(optarg);
				fprintf(stdout, "POPSIZE set to value %d\n", gaPar->popsize);
				break;
			case 2:
				gaPar->fitmate = atoi(optarg);
				fprintf(stdout, "FITMATE set to value %d\n", gaPar->fitmate);
				break;
			case 3:
				gaPar->genenum = atoi(optarg);
				fprintf(stdout, "GENENUM set to value %d\n", gaPar->genenum);
				break;
			case 4:
				gaPar->generation = atoi(optarg);
				fprintf(stdout, "GENERATION set to value %d\n", gaPar->generation);
				break;
			case 5:
				gaPar->lowlim = atoi(optarg);
				fprintf(stdout, "LOWLIM set to value %d\n", gaPar->lowlim);
				break;
			case 6:
				gaPar->uplim = atoi(optarg);
				fprintf(stdout, "UPLIM set to value %d\n", gaPar->uplim);
				break;
			case 7:
				gaPar->minimize = atoi(optarg);
				fprintf(stdout, "MINIMIZE set to value %d\n", gaPar->minimize);
				break;
			case 8:
				gaPar->maximize = atoi(optarg);
				fprintf(stdout, "MAXIMIZE set to value %d\n", gaPar->maximize);
				break;
			case 9:
				gaPar->random = atoi(optarg);
				fprintf(stdout, "RANDOM set to value %d\n", gaPar->random);
				break;
			case 10:
				gaPar->seeded = atoi(optarg);
				fprintf(stdout, "SEEDED set to value %d\n", gaPar->seeded);
				break;
			case 11:
				gaPar->maxvar = atoi(optarg);
				fprintf(stdout, "MAXVAR set to value %d\n", gaPar->maxvar);
				break;
			case 12:
				gaPar->crossover = atoi(optarg);
				fprintf(stdout, "CROSSOVER set to value %d\n", gaPar->crossover);
				break;
			case 13:
				gaPar->equilibrium = atoi(optarg);
				fprintf(stdout, "EQUILIBRIUM set to value %d\n", gaPar->equilibrium);
				break;
			case 14:
				gaPar->jackknife = atoi(optarg);
				fprintf(stdout, "JACKKNIFE set to value %d\n", gaPar->jackknife);
				break;
			case 15:
				gaPar->repeat = atoi(optarg);
				fprintf(stdout, "REPEAT set to value %d\n", gaPar->repeat);
				break;
			default:
				/*usage();*/
				break;	
		}
	}
}

/*____________________________________________________________________________*/
/* main function */
int main(int argc, char **argv)
{
    int i = 0, j = 0, k = 0, l = 0, ix = 0, iy = 0; /* counters */
	FILE *outfile = 0;
    /* digit encoding of 'pool' output files: */
    /* 1.: this repat, 2. gaPar.repga, 3. this fraction, 4. gaPar.jackknife */
    char outfilename[13] = "0_0.0_0.ga";
	int *average = 0;

    /*____________________________________________________________________________*/
	/* print program license */
	license();

    /*____________________________________________________________________________*/
	/* parametrise GA */
	Gapar gaPar; /* GA parameters */
	parametrise_ga(&gaPar); /* assign default GA parameters */

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/
	/* parametrise application */
	/* command line arguments of GA and application are processed together */
	Minset ms;
	parametrise_minset(&ms); /* assign default MINSET parameters */
	parse_args(argc, &(argv[0]), &gaPar, &ms, outfile); /* parse command line arguments */
    /*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*____________________________________________________________________________*/
	/* initialise GA */
    Pool pool[gaPar.popsize]; /* gene pool */
	/* average values for equilibrium */
	if (gaPar.equilibrium)
		average = safe_malloc(gaPar.popsize * sizeof(int)); 

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/
	/* initialise application */
	initialise_minset(&pool[0], &gaPar, &ms);
    /*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*____________________________________________________________________________*/
	/* allocate genomes */
    for(i = 0; i < gaPar.popsize; ++ i) /* allocate memory to genomes */
	{
#ifdef BIT
		pool[i].bitgenome = safe_malloc(sizeof(Bitgene) * gaPar.genenum); /* bitgene pool */
#else
        pool[i].genome = safe_malloc(sizeof(int) * gaPar.genenum); /* gene pool */
#endif
	}

    /*____________________________________________________________________________*/
	/* run GA */
    /* repeat entire GA */
    for (j = 0;  j < gaPar.repeat; ++ j)
    {
        /* set the starting point for pseudorandom integers */
        srand((unsigned int)(time(NULL)%10000));

        /* repeat for each gaPar.jackknife fraction */
        for (k = 0;  k < gaPar.jackknife; ++ k)
        {
            /*____________________________________________________________________________*/
			/* initialise pool */
			if (gaPar.random)
				init_random_pool(&pool[0], &gaPar); /* random pool */
			if (gaPar.seeded)
				init_seeded_pool(&pool[0], &gaPar); /* seeded pool */

			/* if required, constrain the size the newly generated pool */
			if (ms.subsetsize < 100.)
				for (iy = 0; iy < gaPar.popsize; ++ iy)
					constrain_genome(&pool[0], &gaPar, iy, ms.n_selected);

			/* set output file */
			outfile = safe_open(&outfilename[0], "w");

            /*____________________________________________________________________________*/
			/* for 'l' gaPar.generations */
			fprintf(stdout, "jackknife (max %d) / generation (max %d) :\n",
				gaPar.jackknife - 1, gaPar.generation - 1);
			/* first generation starts at ix=0, following generations start at ix=gaPar.fitmate */
            for (l = 0, ix = 0; l < gaPar.generation; ++ l, ix = gaPar.fitmate)
            {
				fprintf(stdout, "%d/%d ", k, l);
				fflush(stdout);

				/* for the population size (minus gaPar.fitmate) */
                for ( ; ix < gaPar.popsize; ++ ix)
                {
					/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/
					/* run application */
					/* application returns genome fitness unless identical genome already in memory */
					if (mem_fitness(&pool[0], &gaPar, ix) == 0)
						pool[ix].fitness = run_minset(&pool[0], &gaPar, &ms, j, k, l, ix);
					/*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
                }

				/*____________________________________________________________________________*/
				/* sort pool */
                sort_fitness(&pool[0], &gaPar); /* sort pool by fitness */
                
				if (check_convergence(&pool[0], &gaPar) == 0)
				{
					fprintf(stdout, "converged\n");
					/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/
					/* print results */
					print_pool_ascii(&pool[0], &gaPar, outfile, l, 1); /* print to file */
					fclose(outfile);
					print_subset(&pool[0], &gaPar, &ms);
					/*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
					exit(0);
				}

				/*____________________________________________________________________________*/
				/* breed new generation */
                for (iy = gaPar.fitmate; iy < gaPar.popsize; ++ iy)
                {
					if (gaPar.crossover)
						breed_crossover(&pool[0], &gaPar, iy); 
					if (gaPar.equilibrium)
					{
						equilibrium(&pool[0], &gaPar, &average[0]); /* equilibrium values */
						breed_equilibrium(&pool[0], &gaPar, &average[0], iy);
					}
                }

				/* if required, constrain the size the newly generated pool */
				if (ms.subsetsize < 100.)
					for (iy = gaPar.fitmate; iy < gaPar.popsize; ++ iy)
						constrain_genome(&pool[0], &gaPar, iy, ms.n_selected);
            }
        }
    }

	/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/
	/* print results */
	print_pool_ascii(&pool[0], &gaPar, outfile, l, 1); /* print to file */
    fclose(outfile);
	print_subset(&pool[0], &gaPar, &ms);
	/*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*____________________________________________________________________________*/
	/* finalise GA */
    /* free memory */
    for (i = 0; i < gaPar.popsize; ++ i)
	{
#ifdef BIT
        free(pool[i].bitgenome);
#else
        free(pool[i].genome);
#endif
	}

	/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/
	/* finalise application */
	finalise_minset(&ms);
    /*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    return 0;
}
