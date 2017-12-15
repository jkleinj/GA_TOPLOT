/*==============================================================================
parse_args.c : parse command line arguments
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

#include <assert.h>
#include "ga.h"
#include "minset.h"
#include "suffix_tree.h"

/*____________________________________________________________________________*/
/* license */
void license()
{
	fprintf(stdout,
		"\nMinSet, Copyright (C) 2006-2007 Jens Kleinjung and Alessandro Pandini\n"
		"MinSet comes with ABSOLUTELY NO WARRANTY.\n"
		"This is free software, and you are welcome to redistribute it under certain conditions.\n"
		"Read the MANUAL for details.\n\n");
}

/*____________________________________________________________________________*/
/* usage */
void usage(Gapar *gapar, Minset *ms)
{
	fprintf(stderr, "Usage: minset --baseset <baseset [FILE]>\n");
	fprintf(stderr, "\nOptions:\n"
		"\t--tag         \t [type] \t default \t description\n"
		"\t__________________________________________________________________________\n"
		"  GA\n"
		"\t--popsize     \t [INT]   \t %3d \t\t population size (number of genomes)\n"
		"\t--fitmate     \t [INT]   \t %3d \t\t fittest genomes : parents of next generation\n"
		"\t--genenum     \t [INT]   \t %3d \t\t number of genes (optimisation parameters)\n"
		"\t--generation  \t [INT]   \t %3d \t\t number of generations (optimisation loops)\n"
		"\t--lowlim      \t [INT]   \t %3d \t\t lower limit of parameter value range\n"
		"\t--uplim       \t [INT]   \t %3d \t\t upper limit of parameter value range\n"
		"\t--minimize    \t [BOOL]  \t %3d \t\t lowest score values are fittest\n"
		"\t--maximize    \t [BOOL]  \t %3d \t\t highest score values are fittest\n"
		"\t--random      \t [BOOL]  \t %3d \t\t random starting population\n"
		"\t--seeded      \t [BOOL]  \t %3d \t\t seeded starting population\n"
		"\t--maxvar      \t [INT]   \t %3d \t\t maximal parameter variation (mutation) in seeded population\n"
		"\t--crossover   \t [BOOL]  \t %3d \t\t crossover breeding\n"
		"\t--equilibrium \t [BOOL]  \t %3d \t\t equilibrium breeding\n"
		"\t--jackknife   \t [INT]   \t %3d \t\t split into 'JACKKNIFE' parts\n"
		"\t--repeat      \t [INT]   \t %3d \t\t repeat GA 'REPEAT' times\n"
		"  MINSET\n"
        "\t--baseset    \t [CHAR]  \t %s \t base set of proteins\n"
        "\t--seqdir      \t [CHAR]  \t %s \t relative path to directory holding sequence files\n"
        "\t--alphabet    \t [CHAR]  \t %s \t coding alphabet\n"
        "\t--subsetsize  \t [FLOAT] \t %3.0f \t\t target size of subset in percent units relative to base set\n"
        "\t--kwordlength \t [INT]   \t %3d \t\t selection k-word (fragment) length\n"
		"\n\talphabet choices:\n"
		"\tMV2000: amino acids (frequencies taken from Mueller and Vingron (2000) J.Comp.Biol.)\n"
		"\tCGT2004: structural fragments (character frequencies taken from Camproux et al. (2004) J.Mol.Biol.)\n"
		"\tTOP2006: topology code (character frequencies derived from ASTRAL40 1.69, Kleinjung (2006) unpublished)\n\n"
		,
		gapar->popsize, gapar->fitmate, gapar->genenum, gapar->generation,
		gapar->lowlim, gapar->uplim, gapar->minimize, gapar->maximize,
		gapar->random, gapar->seeded, gapar->maxvar,
		gapar->crossover, gapar->equilibrium, 
		gapar->jackknife, gapar->repeat, ms->basesetFileName, ms->seqdir,
		ms->alphabet.name, ms->subsetsize, ms->kword_len); 

	exit(1);
}

/*____________________________________________________________________________*/
/* print parameters */
void print_pars(Gapar *gapar, Minset *ms, FILE *outfile)
{
	FILE *parFile; char parFileName[] = "parameters";
	parFile = safe_open(parFileName, "w");

	fprintf(parFile,
		"popsize %3d\n"
		"fitmate %3d\n"
		"genenum %3d\n"
		"generation %3d\n"
		"lowlim %3d\n"
		"uplim %3d\n"
		"minimize %3d\n"
		"maximize %3d\n"
		"random %3d\n"
		"seeded %3d\n"
		"maxvar %3d\n"
		"crossover %3d\n"
		"equilibrium %3d\n"
		"jackknife %3d\n"
		"repeat %3d\n"
        "baseset %s\n"
        "seqdir %s\n"
        "alphabet %s\n"
        "subsetsize %3.0f\n"
        "kwordlength %3d\n",
		gapar->popsize, gapar->fitmate, gapar->genenum, gapar->generation,
		gapar->lowlim, gapar->uplim, gapar->minimize, gapar->maximize,
		gapar->random, gapar->seeded, gapar->maxvar,
		gapar->crossover, gapar->equilibrium, 
		gapar->jackknife, gapar->repeat, ms->basesetFileName, ms->seqdir,
		ms->alphabet.name, ms->subsetsize, ms->kword_len); 

	fclose(parFile);
}

/*____________________________________________________________________________*/
/* parse GA command line long_options */
void parse_args(int argc, char **argv, Gapar *gaPar, Minset *ms, FILE *outfile)
{
	/*extern char *optarg;
	extern int optind;*/
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
        {"baseset", required_argument, 0, 101},
        {"seqdir", required_argument, 0, 102},
        {"alphabet", required_argument, 0, 103},
        {"subsetsize", required_argument, 0, 104},
        {"kwordlength", required_argument, 0, 105},
        {"help", no_argument, 0, 1001},
		{0, 0, 0, 0}
	};

	while ((c = getopt_long (argc, argv, "1:2:3:4:5:6:7:8:9:10:11:12:13:14:15:101:102:103:104:105:1001", long_options, NULL)) != -1)
	{
		switch(c)
		{
			case 1:
				gaPar->popsize = atoi(optarg); assert(gaPar->popsize > 1);
				fprintf(stdout, "POPSIZE set to value %d\n", gaPar->popsize);
				break;
			case 2:
				gaPar->fitmate = atoi(optarg); assert(gaPar->fitmate < gaPar->popsize);
				fprintf(stdout, "FITMATE set to value %d\n", gaPar->fitmate);
				break;
			case 3:
				gaPar->genenum = atoi(optarg); assert(gaPar->genenum > 0);
				fprintf(stdout, "GENENUM set to value %d\n", gaPar->genenum);
				break;
			case 4:
				gaPar->generation = atoi(optarg); assert(gaPar->generation > 0);
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
				gaPar->jackknife = atoi(optarg); assert(gaPar->jackknife > 0);
				fprintf(stdout, "JACKKNIFE set to value %d\n", gaPar->jackknife);
				break;
			case 15:
				gaPar->repeat = atoi(optarg); assert(gaPar->repeat > 0);
				fprintf(stdout, "REPEAT set to value %d\n", gaPar->repeat);
				break;
            case 101:
                strcpy(ms->basesetFileName, optarg); assert (strlen(ms->basesetFileName) > 1);
                fprintf(stdout, "PROTLIST set to name %s\n", ms->basesetFileName);
                break;
            case 102:
                strcpy(ms->basesetFileName, optarg); assert (strlen(ms->seqdir) > 1);
                fprintf(stdout, "SEQDIR set to name %s\n", ms->seqdir);
                break;
            case 103:
                strcpy(ms->alphabet.name, optarg); assert (strlen(ms->alphabet.name) > 1);
                fprintf(stdout, "ALPHABET set to name %s\n", ms->alphabet.name);
                break;
            case 104:
                ms->subsetsize = atof(optarg); assert (ms->subsetsize > 0 && ms->subsetsize <= 100);
                fprintf(stdout, "PERCENTAGE set to value %f\n", ms->subsetsize);
                break;
            case 105:
                ms->kword_len = atoi(optarg); assert (ms->kword_len > 0);
                fprintf(stdout, "KWORDLENGTH set to value %d\n", ms->kword_len);
                break;
			default:
				usage(gaPar, ms);
				break;	
		}
	}

	assert(gaPar->uplim > gaPar->lowlim);
	assert((gaPar->minimize == 0 && gaPar->maximize == 1) || \
		(gaPar->minimize == 1 && gaPar->maximize == 0));
	assert((gaPar->random == 0 && gaPar->seeded == 1) || \
		(gaPar->random == 1 && gaPar->seeded == 0));
	assert((gaPar->equilibrium == 0 && gaPar->crossover == 1) || \
		(gaPar->equilibrium == 1 && gaPar->crossover == 0));


	print_pars(gaPar, ms, outfile);
	fprintf(stdout, "\n");
}

