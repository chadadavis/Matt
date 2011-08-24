/* Protein.h -- Prototypes and global constants for dealing with proteins and
 * amino acid/nucleotide sequences.
 *
 * Author: Matt Menke (2007)
 *
 * Matt is licensed under the GNU public license version 2.0. 
 *
 * If you would like to license Matt in an enviroment where the GNU public license is 
 * unacceptable (such as inclusion in a non-GPL software package) comercial Matt
 * licensing is available through the MIT and Tufts offices of Technology Transfer. 
 * Contact betawrap@csail.mit.edu or cowen@cs.tufts.edu for more information.
 */

#ifndef PROTEIN_H
#define PROTEIN_H


#include "util.h"

#define NAME_RESIDUE 1
#define NAME_DNA 2
#define NAME_BOTH (NAME_RESIDUE | NAME_DNA)

char LookupName(char *name, int type);

const extern int lookupTable[26];

const extern char ShortNames[42];
const extern char LongNames[42][4];

/* Gets the residue symbol from a letter.
 * Works with both upper and lower case.
 * Returns unknown code for everything that's not
 * a letter.
 */
int GetSymbol(char c);

struct SEQUENCE {
	char *name;
	int length;
	char *seq;
};
typedef struct SEQUENCE Sequence;

void AddCharacter(Sequence *s, char residue);
Sequence *CreateSequence(char *name);
void FreeSequence(Sequence *s);

/* Outputs the sequence in fasta format to f.  Can output
 * multiple sequences to the same file.  Returns 1 on succes,
 * 0 on failure.  If first or last are given, only displays
 * sequences from first to last.  If name is given, uses
 * that as the FASTA entry's name instead of the sequence's name.
 * A negative value for last means output to the end.
 * Assumes that the sequence uses the internal codes from 0 to 21.
 */

struct SEQUENCE_ALIGNMENT {
	Sequence ** seqs;
	int numSeqs;
};
typedef struct SEQUENCE_ALIGNMENT SequenceAlignment;

void FreeSequenceAlignment(SequenceAlignment *s);

int WriteSequenceFASTA(FILE *file, Sequence *s);
int WriteSequenceAlignmentFASTA(FILE *file, SequenceAlignment *s);

int WriteSequencePIR(FILE *file, Sequence *s);
int WriteSequenceAlignmentPIR(FILE *file, SequenceAlignment *s);

void Reorder(SequenceAlignment *s);

#endif
