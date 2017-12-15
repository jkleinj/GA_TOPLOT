/*=============================================================================
minset.c : a method for database subsetting
(C) 2006-2017 Jens Kleinjung
(C) 2006-2007 Alessandro Pandini

References:
    - A. Pandini, L. Bonati, F. Fraternali and J. Kleinjung
    MinSet: A general approach to derive maximally representative database
    subsets by using fragment dictionaries and its application to the SCOP
    database. Bioinformatics 23 (2007) 515-516.
    - A. Pandini and J. Kleinjung
    MinSet: automatic derivation of maximally representative subsets of
    structural domains. Nucleic Acids Res., submitted.

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
=============================================================================*/

#include "alphabet.h"
#include "getseqs.h"
#include "parse_args.h"
#include "suffix_tree.h"
#include "minset.h"
#include "minsetpar.h"

/*____________________________________________________________________________-*/  
/* get single-character code counts*/
void get_counts(char *seq, int *count, int alphabet_array_size)
{
    int i;

	/* initialise */
    for (i = 0; i < alphabet_array_size; ++ i)
        count[i] = 0;

    while ((*seq) != 0)
	{
		if (isalpha(*seq))
			++ count[(*seq)-'A'];
		++ seq;
	}
}

/*____________________________________________________________________________-*/  
/* Shannon formula */
float shannon(float p)
{
    return (-1 * (p * log(p) / log(2)));
}

/*____________________________________________________________________________-*/  
/* Shannon entropy of code */
float shannon_entropy(float *p, int length)
{
    unsigned int i;
    float H;

    for (i = 0, H = 0; i < length; ++ i)
    {
        if (p[i] == 0)
            continue;
        H += shannon(p[i]);
    }

    return H;
}

/*____________________________________________________________________________-*/  
/* relative entropy of code (Kullback-Leibler distance) */
float relative_entropy(float *p, float *q, int length)
{
    unsigned int i;
    float D;

    for (i = 0, D = 0; i < length; ++ i)
    {
        if (p[i] == 0 || q[i] == 0)
            continue;
        D += p[i] * log(p[i] / q[i]) / log(2);
    }

    return D;
}

/*____________________________________________________________________________-*/  
/* word entropy from suffix tree: */
/* the entropy of the subset based on a constant-length word alphabet */
float word_entropy(char *ref_string, char *search_string, int *p_nsymbols, int *p_kword_len)
{
    SUFFIX_TREE* tree; /* the suffix tree */
    DBL_WORD ref_len = 0; /* length of reference string */
    DBL_WORD search_len = 0; /* length of search (sub)string */
    DBL_WORD sub_len = 0; /* selection k-word length */
    unsigned int i = 0;
    char *sub_pc; /* pointer to char for search substring */
    unsigned int n_all = 0; /* total number of substrings */
    //DBL_WORD position; /* position of search (sub)string in reference string */
    float entropy;

    ref_len = (DBL_WORD)strlen(ref_string); /* the string to become the ST */
    search_len = (DBL_WORD)strlen(search_string); /* the string to use for searching */
	sub_len = (DBL_WORD)*p_kword_len; /* the length of the search words in the search string */
	/*fprintf(stderr, "ref_len %lu, search_len %lu, sub_len %lu\n", 
							ref_len, search_len, sub_len);*/
	assert(ref_len > 0);
	assert(search_len > 0);
	assert(sub_len > 0);
	assert(search_len <= ref_len);
	assert(sub_len <= search_len);

    tree = ST_CreateTree(ref_string, ref_len); /* generate the suffix tree */
    ST_InitTreeHits(tree);

    /* for the length of the search string */
    while (i < search_len - (sub_len - 1))
    {
        /* move pointer along search string to create substrings */
        sub_pc = search_string + i;

        /* skip all substrings containing the inter-sequence delimiter '-' */
        if (sub_pc[sub_len - 1] == '-')
        {
            i += sub_len; /* increment pointer position by substring length */
            continue;
        }

        ++ n_all; /* count all subSetSequences */

        /* search substring in suffix tree */
        //position = ST_FindSubstring(tree, sub_pc, sub_len);

        ++ i; /* increment pointer position */
    }

    assert(n_all == (tree->allhit + tree->allmiss));

    entropy = ST_TreeEntropy(tree);
    *p_nsymbols = tree->allsymbol; 

    ST_DeleteTree(tree);

    return entropy;
}

/*____________________________________________________________________________-*/  
/* score the seq by compression ratio */
float score_compress(char *instring, int insize, int total_len)
{
    unsigned char *out;
    unsigned int *work;
    unsigned int outsize = 0; 
    int maxoutsize;
    float compression_ratio;
    float zipscore;
    
    /* maximum output size w/ lz : (257/256)*insize + 1. */
    maxoutsize = ceil(257.0/256.0*insize) + 1;

/*  work - Pointer to internal working buffer, which must be able to hold (insize+65536) unsigned integers.*/
    work = safe_malloc(sizeof(unsigned int) * (insize + 65536));
    out = safe_malloc(sizeof(unsigned char) * maxoutsize);

/*  intersize = Huffman_Compress(instring, out, insize); */
/*    outsize   = LZ_CompressFast((unsigned char *) instring, out, insize, work); */

    if (insize == 0)
    {
        free(work);
        free(out);
        return 0.0; /* exclude null chromosome cases */
    }

    compression_ratio = (float) (insize - outsize) / insize; 
    compression_ratio = max(compression_ratio, 0.1);

    zipscore = (log(1 + insize / total_len) / log(2)) / compression_ratio; 

#ifdef DEBUG
        fprintf(stderr, "--> %d\t%d\t %d\t%5.1f\t%5.1f\t", 
				insize, outsize, total_len, compression_ratio, zipscore*100);
#endif

    free(work);
    free(out);
    return zipscore;
}

/*____________________________________________________________________________-*/  
/* parse the protein list file 
	the format is simply a list of filenames */
void parse_proteinlist(Minset *ms)
{
    FILE *listfile;
    char line[256];
    unsigned int i;
    unsigned int k = 0;
    int allocated = 64;

	/*____________________________________________________________________________-*/  
    /* record the protein names in 'prots' array */
    ms->prots.protein = safe_malloc(allocated * sizeof(ProteinEntry));

    listfile = safe_open(ms->basesetFileName, "r");

    while (fgets(line, 256, listfile))
    { 
		ms->prots.protein[k].name = safe_malloc(256 * sizeof(char));
		for (i = 0; i < 256 && isprint(line[i]) > 0; ++ i)
			ms->prots.protein[k].name[i] = line[i];
		ms->prots.protein[k].name[i] = '\0';

		if (i > 0)
			++ k;
		else
			fprintf(stderr, "Warning: empty protein name at line %d\n", k);

		if (k == allocated)
		{
			allocated += 64;
			ms->prots.protein = safe_realloc(ms->prots.protein, allocated * sizeof(ProteinEntry));
		}
    }

    fclose(listfile);

    ms->prots.n_prot = k;
}

/*____________________________________________________________________________*/
/* score the concatenated sequence */
float score_seq(Minset *ms, char *subSetSeq)
{
    int i;
    int strlengthSub;
	int n_symbol; /* number of k-word alphabet symbols (k-words) */
    float H = 0.; /* entropy */
    float E = 0.; /* expected entropy */
    float D = 0.; /* relative entropy */
    float score = 0.; /* fitness score */
    int charCount[ms->alphabet.size];
    float *p_count; /* relative frequency of chode character */
    float *p_bg; /* same for background distribution */

    p_count = safe_malloc(sizeof(float) * ms->alphabet.size);
    p_bg = safe_malloc(sizeof(float) * ms->alphabet.size);

    strlengthSub = strlen(subSetSeq); /* count before moving the seq pointer on*/

    get_counts(subSetSeq, charCount, ms->alphabet.size);

    for (i = 0; i < ms->alphabet.size; ++ i)
	{
        p_count[i] = charCount[i] / (float) strlengthSub; /* foregrund frequencies */
        p_bg[i] = ms->alphabet.freq[i]; /* background frequencies */
	}

    /*H = shannon_entropy(p_count, alphabet_array_len);*/ /* single-character alphabet */

    H = word_entropy(subSetSeq, subSetSeq, &n_symbol, &ms->kword_len); /* contant-length word alphabet */
    E = shannon_entropy(p_bg, ms->alphabet.size);
    D = relative_entropy(p_count, p_bg, ms->alphabet.size);

    score = H/E * (1 - D);

#ifdef DEBUG
    fprintf(stderr, "l: %d, H: %f, E: %f, D: %f, score: %f\n",
            strlengthSub, H, E, D, score);
#endif

    free(p_count);
    free(p_bg);

    return score;
}

/*____________________________________________________________________________*/
/* fill protein table with sequences and their entropies */
void fill_protein_table(Minset *ms)
{
	unsigned int k;
	FILE *fastaFile; char *fastaFileName;
	FILE *setFile;

	/* initialise counters */
	ms->total_len = 0;

	ms->setfasta_charCount = safe_malloc(ms->alphabet.size * sizeof(int));

	/*____________________________________________________________________________*/
	/* read all FASTA sequences */
    for (k = 0; k < ms->prots.n_prot; ++ k)
    {
		/* read sequence 'k' */
        fastaFileName = safe_malloc(strlen(ms->seqdir) + strlen(ms->prots.protein[k].name) + 6);
		strcpy(fastaFileName, "");
        sprintf(fastaFileName, "%s%s.tseq", ms->seqdir, ms->prots.protein[k].name);
	/*fprintf(stdout, "Reading sequence: %s\n", fastaFileName);*/
	fastaFile = safe_open(fastaFileName, "r");
        read_sequence(fastaFile, &(ms->prots), k);
        fclose(fastaFile);
		free(fastaFileName);

		/* add code character frequencies of this sequence to overall count */
        get_counts(ms->prots.protein[k].seq, &ms->setfasta_charCount[0], ms->alphabet.size);

		/* add length of this sequence to overall length */
        ms->total_len += strlen(ms->prots.protein[k].seq);

		/* compute the entropy of this sequence */
#ifndef COMPRESS_SCORE
        ms->prots.protein[k].entropy = score_seq(ms, ms->prots.protein[k].seq);
#endif
#ifdef COMPRESS_SCORE
        ms->prots.protein[k].entropy = score_compress(ms->prots.protein[k].seq, strlen(ms->prots.protein[k].seq));
#endif
    }

	/*____________________________________________________________________________*/
	/* concatenate all base set sequences to 'setfasta' with '-' delimiter */
    ms->setfasta = safe_malloc((ms->total_len + ms->prots.n_prot + 1) * sizeof(char));
    strcpy(ms->setfasta, ms->prots.protein[0].seq);
    for (k = 1; k < ms->prots.n_prot; ++ k)
    {
        strcat(ms->setfasta, "-");
        strcat(ms->setfasta, ms->prots.protein[k].seq);
    }

	/*____________________________________________________________________________*/
	/* print concatenated 'setfasta' sequence */
    setFile = safe_open("baseset.seq", "w");
	fprintf(setFile, "%s", ms->setfasta);
	fclose(setFile);

	free(ms->setfasta_charCount);
	free(ms->setfasta);
}

/*____________________________________________________________________________*/
/* calculate fitness of (concatenated) selected protein sequences */
float calculate_fitness(Pool *pool, Gapar *gaPar, int ix, Minset *ms)
{
    int i;
    int allocated = 1;

    ms->polyfasta = safe_malloc(allocated);
    strcpy(ms->polyfasta, "");

    for (i = 0; i < gaPar->genenum; ++ i)
    {
        if (pool[ix].genome[i] == 1)
        {
			allocated += strlen(ms->prots.protein[i].seq) + 1;
			ms->polyfasta = safe_realloc(ms->polyfasta, allocated * sizeof(char));
			strcat(ms->polyfasta, ms->prots.protein[i].seq);
			strcat(ms->polyfasta, "-");
        }
    }

#ifndef COMPRESS_SCORE
    pool[ix].fitness = score_seq(ms, ms->polyfasta);
#endif
#ifdef COMPRESS_SCORE
    pool[ix].fitness = score_compress(ms->polyfasta, strlen(ms->polyfasta), ms->total_len);
#endif
#ifdef DEBUG
	dump2(ms->polyfasta, "%s", pool[ix].fitness, "%f");
#endif

	ms->polyfasta_charCount = safe_malloc(ms->alphabet.size * sizeof(int));
    get_counts(ms->polyfasta, ms->polyfasta_charCount, ms->alphabet.size);

	free(ms->polyfasta);
	free(ms->polyfasta_charCount);

	return pool[ix].fitness;
}

/*____________________________________________________________________________*/
/* parse MinSet command line long_options */
void parse_args_ms(int argc, char **argv, Minset *ms)
{
    int c;

    /* MinSet command line long_options */
    static struct option long_options[] =
    {
        {"basesetFileName", required_argument, 0, 1},
        {"alphabet", required_argument, 0, 2},
        {"subsetsize", required_argument, 0, 3},
		{"kwordlength", required_argument, 0, 4},
        {0, 0, 0, 0}
    };

    while ((c = getopt_long (argc, argv, "1:2:3:4:", long_options, NULL)) != -1)
    {
        switch(c)
        {
            case 1:
                strcpy(ms->basesetFileName, optarg);
                fprintf(stdout, "PROTLIST set to name %s\n", ms->basesetFileName);
                break;
            case 2:
                /*strcpy(ms->alphabet.name, optarg);*/
                fprintf(stdout, "ALPHABET set to name %s\n", ms->alphabet.name);
                break;
            case 3:
                ms->subsetsize = atof(optarg);
                fprintf(stdout, "PERCENTAGE set to value %f\n", ms->subsetsize);
                break;
            case 4:
                ms->kword_len = atoi(optarg);
                fprintf(stdout, "KWORDLENGTH set to value %d\n", ms->kword_len);
                break;
            default:
                /*usage();*/
                break;
		}
	}
}

/*____________________________________________________________________________*/
/* parametrise MINSET */
void parametrise_minset(Minset *ms)
{
	strcpy(ms->basesetFileName, BASESET); assert (strlen(ms->basesetFileName) > 1);
	strcpy(ms->seqdir, SEQDIR); assert (strlen(ms->seqdir) > 1);
	strcpy(ms->alphabet.name, ALPHABET); assert (strlen(ms->alphabet.name) > 1);
	ms->subsetsize = (float)SUBSETSIZE; assert (ms->subsetsize > 0 && ms->subsetsize <= 100);
	ms->kword_len = (int)KWORDLENGTH; assert (ms->kword_len > 0);
}

/*____________________________________________________________________________*/
/* print subset composition and string */
void print_subset(Pool *pool, Gapar *gaPar, Minset *ms)
{
	unsigned int i, j;
	/*unsigned int k;*/
	FILE *subsetSeqFile;

    subsetSeqFile = safe_open("subset.seq", "w");

	/* print table: protein names and gene values of fittest genomes */
	for (j = 0; j < gaPar->genenum; ++ j)
	{
		/* sequence of protein names is identical to gene sequence */
		/* and print concatenated subset sequences delimited with '+' */
		fprintf(ms->subsetOutFile, "%7s ", ms->prots.protein[j].name);

		for (i = 0; i < gaPar->fitmate; ++ i)
		{
			fprintf(ms->subsetOutFile, " %1d", pool[i].genome[j]);

			if (i == 0 && pool[i].genome[j] == 1)
				fprintf(subsetSeqFile, "%s+", ms->prots.protein[j].seq);
		}

		fprintf(ms->subsetOutFile, "\n");
	}
	
	/*____________________________________________________________________________*/
	/* concatenate all subset sequences to 'subsetfasta' with '+' delimiter */
	/*
    ms->subsetfasta = safe_malloc((ms->total_len + ms->prots.n_prot + 1) * sizeof(char));
    strcpy(ms->subsetfasta, "");
    for (k = 0; k < ms->prots.n_prot; ++ k)
    {
		if (pool[0].genome[k] != 0)
		{
			strcat(ms->subsetfasta, ms->prots.protein[k].seq);
			strcat(ms->subsetfasta, "+");
		}
	}
	*/
	
	fclose(subsetSeqFile);
	free(ms->subsetfasta);
}

/*____________________________________________________________________________*/
void initialise_minset(Pool *pool, Gapar *gaPar, Minset *ms)
{
    /* set character frequencies of selected alphabet */
	set_alphabet(&(ms->alphabet));

    /*____________________________________________________________________________*/
    /* read and process sequence input */
	/* file containing the list of sequence filenames */
	parse_proteinlist(ms);

	/* read sequences, calculate entropy score per sequence */
	/* compute overall code character counts and overall sequence length */
	/* concatenate sequences to total 'setfasta' sequence */
	fill_protein_table(ms);

	gaPar->genenum = ms->prots.n_prot; /* overwrite number of genes in GA */
	/* set number of genes that are switched ON: to achieve subset target size */
	ms->n_selected = (int)floorf(gaPar->genenum * ms->subsetsize / 100);
	assert(ms->n_selected > 1);

    /*____________________________________________________________________________*/
	/* basic output file (name might be modified in 'run_minset') */
	ms->subsetOutFileName = safe_malloc(13 * sizeof(char));
	strcpy(ms->subsetOutFileName, "subset.list");
	ms->subsetOutFile = safe_open(ms->subsetOutFileName, "w");
}

/*____________________________________________________________________________*/
/* run minset */
float run_minset(Pool *pool, Gapar *gaPar, Minset *ms, int j, int k, int l, int ix)
{
	/* compute fitness of genome 'ix' */ 
	return calculate_fitness(&pool[0], gaPar, ix, ms);
}

/*____________________________________________________________________________*/
/* finalise minset */
void finalise_minset(Minset *ms)
{
	unsigned int i;

	/* close last output file */
	fclose(ms->subsetOutFile);

	/* alphabet */
	free(ms->alphabet.codeOrder);
	free(ms->alphabet.freq);

	/* filenames */
	free(ms->subsetOutFileName);

	/* protein set */
    for (i = 0; i < ms->prots.n_prot; ++ i)
    {
        free(ms->prots.protein[i].name);
        free(ms->prots.protein[i].description);
        free(ms->prots.protein[i].seq);
	}
    free(ms->prots.protein);
}

