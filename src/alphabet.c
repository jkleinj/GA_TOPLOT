/*==============================================================================
alphabet.c : alphabet assignment routine
(C) 2006-2017 Jens Kleinjung
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
==============================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "alphabet.h"

/*___________________________________________________________________________*/
/* AminoAcid Frequency (ARNDCQEGHILKMFPSTWYV):
	Mueller, T. & Vingron, M., J Comp Biol (2000) 7:761-776. [Table 3]  */
static const float MV2000[] =
{
	0.0771 /* A */,
	0.0501 /* R */,
	0.0462 /* N */,
	0.0538 /* D */,
	0.0146 /* C */,
	0.0409 /* Q */,
	0.0634 /* E */,
	0.0656 /* G */,
	0.0219 /* H */,
	0.0592 /* I */,
	0.0976 /* L */,
	0.0592 /* K */,
	0.0221 /* M */,
	0.0414 /* F */,
	0.0477 /* P */,
	0.0707 /* S */,
	0.0568 /* T */,
	0.0127 /* W */,
	0.0324 /* Y */,
	0.0669 /* V */
};

/*___________________________________________________________________________*/
/* Camproux Structural Alphabet Frequency (aAVWZBCDEOSRQIFUPHGYJKLNMTX):
	Camproux, A.C., Gautier, R. and Tuffery, P.,  J Mol Bol (2004) 339:591-605. [Table 1]  */
static const float CGT2004[] =
{
	0.026 /* a */,
	0.126 /* A */,
	0.056 /* V */,
	0.053 /* W */,
	0.045 /* Z */,
	0.047 /* B */,
	0.018 /* C */,
	0.020 /* D */,
	0.020 /* E */,
	0.015 /* O */,
	0.032 /* S */,
	0.017 /* R */,
	0.041 /* Q */,
	0.029 /* I */,
	0.019 /* F */,
	0.020 /* U */,
	0.044 /* P */,
	0.027 /* H */,
	0.034 /* G */,
	0.020 /* Y */,
	0.020 /* J */,
	0.041 /* K */,
	0.051 /* L */,
	0.049 /* N */,
	0.053 /* M */,
	0.030 /* T */,
	0.047 /* X */
};

/*___________________________________________________________________________*/
/* topology alphabet of TOPLOT program (ABCDEFGHIabcdefghi):
	frequencies derived from 'aha40.ls' list of scop1.67 using the
	'char_freq' program on the concatenated topology strings */
static const float TOP2006[] =
{
	0.113757 /* A */, 
	0.015368 /* B */,
	0.004562 /* C */,
	0.087518 /* D */,
	0.007862 /* E */,
	0.003041 /* F */,
	0.090818 /* G */,
	0.006147 /* H */,
	0.002653 /* I */,
	0.104665 /* J */,
	0.025689 /* K */,
	0.007150 /* L */,
	0.156950 /* a */,
	0.041737 /* b */,
	0.021904 /* c */,
	0.107221 /* d */,
	0.020674 /* e */,
	0.011065 /* f */,
	0.093536 /* g */,
	0.020480 /* h */,
	0.009318 /* i */,
	0.027113 /* j */,
	0.013459 /* k */,
	0.007312 /* l */
};

/*___________________________________________________________________________*/
void print_alphabet(Alphabet *alphabet)
{
	int i;

	fprintf(stdout, "Alphabet %s:\n", alphabet->name);

	for (i = 0; i < alphabet->size; ++ i)
		if (alphabet->freq[i] != 0)
			fprintf(stdout, "%d\t%c\t%f\n", i, (i + 'A'), alphabet->freq[i]);
}

/*___________________________________________________________________________*/
void set_alphabet(Alphabet *alphabet)
{
	int i;

	/* data structure of coding alphabet */
	static const struct myAlphabet
	{
		const char *name; /* name srting of alphabet */
		const char *codeOrder; /* sequence of coding characters */
		const float *freq; /* frequency of coding characters */
		const int   size; /* range of ASCII values of coding characters */
		const int   codeLength; /* number of coding characters */
	}
	/* array of currently implemented alphabets */
	alphabet_array[] =
	{
		/* AminoAcid Frequency:
			Mueller, T. & Vingron, M., J Comp Biol (2000) 7:761-776. [Table 3]  */
		{"MV2000", "ARNDCQEGHILKMFPSTWYV", MV2000, 26, 20},
		/* Camproux Structural Alphabet Frequency:
			Camproux, A.C., Gautier, R. and Tuffery, P.,  J Mol Bol (2004) 339:591-605. [Table 1]  */
		{"CGT2004", "aAVWZBCDEOSRQIFUPHGYJKLNMTX", CGT2004, 33, 27},
		/* topology alphabet of TOPLOT program
			frequencies derived from 'aha40.ls' list of scop1.67 using the
			'char_freq' program on the concatenated topology strings */
		{"TOP2006", "ABCDEFGHIJKLabcdefghijkl", TOP2006, 44, 24}
	};

	const struct myAlphabet *myAlphabet;

	/* match alphabet with specified alphabet name */
    for (myAlphabet = alphabet_array; strcmp(alphabet->name, myAlphabet->name) != 0; myAlphabet ++)
		if (myAlphabet == alphabet_array + sizeof alphabet_array / sizeof *alphabet_array - 1)
		{
			fprintf(stderr, "\nExiting: Alphabet '%s' not implemented!\nAvailable alphabets:\n"
							"\tMV2000 : amino acid frequencies\n"
							"\t\tMueller and Vingron (2000) J. Comp. Biol. 7, 761-776. [Table 3]\n"
							"\tCGT2004 : structural alphabet residue frequencies\n"
							"\t\tCamproux, Gautier and Tuffery (2004) J. Mol. Bol. 339, 591-605. [Table 1]\n\n"
							"\tTOP2006 : topology alphabet sec.str. pair frequencies\n"
							"\t\tKleinjung 2006\n\n",
							alphabet->name);
			exit(1);
		}

	/* set alphabet-specific values */
	alphabet->codeOrder = malloc(((int)strlen(myAlphabet->codeOrder) + 1) * sizeof(char));
	strcpy(alphabet->codeOrder, myAlphabet->codeOrder);
	alphabet->size = myAlphabet->size;
	alphabet->codeLength = strlen(myAlphabet->codeOrder);
	alphabet->freq = malloc(myAlphabet->size * sizeof(float));

    /* assign code character frequencies */
    for (i = 0; i < alphabet->codeLength; ++ i)
		alphabet->freq[alphabet->codeOrder[i] - 'A'] = myAlphabet->freq[i];

    /* set unused characters frequencies to zero */
    for (i = 0; i < alphabet->size; ++ i)
		if (strchr(alphabet->codeOrder, (i + 'A')) == 0)
			alphabet->freq[i] = 0;

	/*print_alphabet(alphabet);*/
}
