/* Protein.c -- Functions and global constants for dealing with proteins and
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

#include "Protein.h"
#include "util.h"

#include <string.h>
#include <malloc.h>
#include <stdlib.h>

const int lookupTable[26] =
{  5,21,10,16,17,2,9,8,1,20,18,7,6,13,20,14,15,19,11,12,20,0,4,20,3,22};

const char ShortNames[42]   = {  'V',   'I',   'F',   'Y',
						   'W',   'A',   'M',   'L',
						   'H',   'G',   'C',   'S',
						   'T',   'N',   'P',   'Q',
						   'D',   'E',   'K',   'R',   'X',
						   'B',   'Z',
						   'A',   'T',   'G',   'C',   'U',  'I',
						   'A',   'T',   'G',   'C',   'U',  'I',
						   'A',   'T',   'G',   'C',   'U',  'I'};

const char LongNames[42][4] = {"VAL", "ILE", "PHE", "TYR",
						 "TRP", "ALA", "MET", "LEU",
						 "HIS", "GLY", "CYS", "SER",
						 "THR", "ASN", "PRO", "GLN",
						 "ASP", "GLU", "LYS", "ARG", "UNK",
						 "ASX", "GLZ",
						 "  A", "  T", "  G", "  C", "  U", "  I",
						 " +A", " +T", " +G", " +C", " +U", " +I",
						 " DA", " DT", " DG", " DC", " DU", " DI"};

int GetSymbol(char c) {
	char i = UCASE(c)-'A';
	if (i<0 || i>=36) return lookupTable['X'-'A'];
	return lookupTable[i];
}

char LookupName(char *name, int type) {
	int i;
	if (type & NAME_RESIDUE) {
		for (i=0; i<20; i++) {
			if (!strcmp(name, LongNames[i])) return i;
		}
		for (i=21; i<23; i++) {
			if (!strcmp(name, LongNames[i])) return i;
		}
	}
	if (type & NAME_DNA) {
		for (i=23; i<35; i++) {
			if (!strcmp(name, LongNames[i])) return i;
		}
	}
	return 20;
}


void FreeSequence(Sequence *s) {
	if (s) {
		if (s->name) free(s->name);
		if (s->seq) free(s->seq);
		free(s);
	}
}

void AddCharacter(Sequence *s, char residue) {
	if (s->length % 16 == 0) {
		s->seq = (char*) realloc(s->seq, (s->length+16) * sizeof(char));
	}
	s->seq[s->length++] = residue;
}

Sequence *CreateSequence(char *name) {
	Sequence *out = (Sequence*) calloc(sizeof(Sequence), 1);
	if (name)
		out->name = strdup(name);
	return out;
}

int WriteSequenceFASTA(FILE * file, Sequence *s) {
	char *name = "Unnamed";
	int i;
	if (!file) return 0;
	if (s->name) name = s->name;
	fprintf(file, ">%s\n", name);
	for (i=0; i < s->length; i += 60) {
		if (s->length >= i+60) {
			fwrite(s->seq+i, sizeof(char), 60, file);
		}
		else {
			fwrite(s->seq+i, sizeof(char), s->length-i, file);
		}
		fputc('\n', file);
	}
	return 1;
}

int WriteSequencePIR(FILE * file, Sequence *s) {
	char *name = "Unnamed";
	int i, j;
	if (!file) return 0;
	if (s->name) name = s->name;
	fprintf(file, ">P1;%s\n%s\n", name, name);
	for (i=0; i < s->length; i += 60) {
		if (i) fputc('\n', file);
		for (j=i; j < i + 60 && j < s->length; j+= 10) {
			if (j == i)
				fputc(' ', file);
			fputc(' ', file);
			if (s->length >= j+10) {
				fwrite(s->seq+j, sizeof(char), 10, file);
			}
			else {
				fwrite(s->seq+j, sizeof(char), s->length-j, file);
			}
		}
	}
	fprintf(file, "*\n");
	return 1;
}

void FreeSequenceAlignment(SequenceAlignment *s) {
	int i;
	for (i=0; i<s->numSeqs; i++) {
		FreeSequence(s->seqs[i]);
	}
	free(s->seqs);
	free(s);
}

int WriteSequenceAlignmentFASTA(FILE *file, SequenceAlignment *s) {
	int i;
	if (!file) return 0;
	for (i=0; i<s->numSeqs; i++) {
		if (!WriteSequenceFASTA(file, s->seqs[i])) return 0;
	}
	return 1;
}

int WriteSequenceAlignmentPIR(FILE *file, SequenceAlignment *s) {
	int i;
	if (!file) return 0;
	for (i=0; i<s->numSeqs; i++) {
		if (!WriteSequencePIR(file, s->seqs[i])) return 0;
		fputc('\n', file);
	}
	return 1;
}

#define IsBlank(x) ((x) == '-' || (x) == '.' || (x) == ' ')

void Reorder(SequenceAlignment *s) {
	if (s->numSeqs) {
		int i;
		for (i=0; i<s->seqs[0]->length; i++) {
			int j=i;
			while (j) {
				int c;
				for (c=0; c<s->numSeqs; c++) {
					if (!IsBlank(s->seqs[c]->seq[j])) break;
					if (!IsBlank(s->seqs[c]->seq[j-1])) break;
				}
				if (!IsBlank(s->seqs[c]->seq[j-1])) break;
				for (; c<s->numSeqs; c++) {
					if (IsBlank(s->seqs[c]->seq[j])) continue;
					if (IsBlank(s->seqs[c]->seq[j-1])) continue;
					break;
				}
				if (c != s->numSeqs) {
					break;
				}
				for (c=0; c<s->numSeqs; c++) {
					char temp = s->seqs[c]->seq[j];
					s->seqs[c]->seq[j] = s->seqs[c]->seq[j-1];
					s->seqs[c]->seq[j-1] = temp;
				}
				j--;
			}
		}
	}
}
