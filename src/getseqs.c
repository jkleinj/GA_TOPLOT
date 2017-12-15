/*==============================================================================
getseqs.c : Read alignment of FASTA sequences
(C) 2004 John Romein
==============================================================================*/

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>

#include "minset.h"

/*____________________________________________________________________________*/
static int read_sequence_name(FILE *file, ProteinEntry *protein)
{
    /* Reads a line starting with '>' and followed by the sequence name */
    /* Leading white space is skipped.  Reading proceeds until the end  */
    /* of the line is reached.  The name is read into sequence->name.   */

    int ch, length = 0, allocated = 64;

    while ((ch = getc(file)) != EOF && isspace(ch))
	;

    if (ch != '>')
		return 0;

    protein->description = safe_malloc(allocated);

    do
	{
		protein->description[length ++] = ch;

		if (length == allocated)
			protein->description = safe_realloc(protein->description, allocated += 64);
    } while ((ch = getc(file)) != EOF && ch != '\n' && isprint(ch));

    protein->description[length] = '\0';

	return 1;
}

/*____________________________________________________________________________*/
static void read_sequence_residues(FILE *file, ProteinEntry *protein)
{
    /* Reads the residues in a sequence, up to (but not including) the      */
    /* next sequence header (starting with '>'), or up to end of file.      */
    /* Residues may span multiple lines.  White space and gaps are skipped. */
    /* Nonalpha characters are rejected, resulting in an error message.     */
    /* Alpha characters are NOT converted to upper case.  The string is read*/
    /* into sequence->residues and zero-terminated; the length is stored    */
    /* into sequence->length.                                               */

    int ch, length = 0, allocated = 64;

    protein->seq = safe_malloc(allocated);

    while ((ch = getc(file)) != EOF && ch != '>')
		if (isalpha(ch))
		{
			protein->seq[length ++] = toupper(ch);

			if (length == allocated)
				protein->seq = safe_realloc(protein->seq, allocated += 64);
		}
		else if (ch != '.' && !isspace(ch))
		{
			fprintf(stderr, "illegal character '%c' in protein sequence\n", ch);
			exit(1);
		}

    if (ch == '>')
		ungetc('>', file);

    if (length == 0)
	{
		fprintf(stderr, "zero-sized sequence\n");
		exit(1);
    }

    protein->seq[length] = '\0';
}

/*____________________________________________________________________________*/
int read_sequence(FILE *file, Prots *prots, int k)
{
    if (read_sequence_name(file, &prots->protein[k]))
	{
		read_sequence_residues(file, &prots->protein[k]);
		return 1;
    }
	else
		return 0;
}

