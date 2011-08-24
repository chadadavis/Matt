#include "MultipleAlignmentOutput.h"
#include <malloc.h>
#include <time.h>

void ReorderAlignment(MultipleAlignment *ma, int *order) {
	int i, j;
	for (i=0; i<ma->numChains; i++) {
		for (j=i; 1; j++) {
			if (order[i] == ma->chains[j]->id) {
				ResiduePositions *chain = ma->chains[i];
				WeightedResiduePositions residues = ma->residues[i];

				ma->chains[i] = ma->chains[j];
				ma->residues[i] = ma->residues[j];

				ma->chains[j] = chain;
				ma->residues[j] = residues;
				break;
			}
		}
	}
}

/* Converts a multiple alignment (possibly with gaps) to a bunch of
 * aligned single-letter seqeunces.  Also reorders aligned residues.
 */
SequenceAlignment *GetSequenceAlignment(MultipleAlignment *ma) {
	SequenceAlignment *out = (SequenceAlignment *) malloc(sizeof(SequenceAlignment));
	int chain, block, chain2;
	int *indices = (int*) calloc(ma->numChains, sizeof(int));
	out->numSeqs = ma->numChains;
	out->seqs = (Sequence**) malloc(sizeof(Sequence*) * ma->numChains);

	for (chain=0; chain<ma->numChains; chain++) {
		out->seqs[chain] = CreateSequence(ma->chains[chain]->pdb->idString);
	}
	for (block=0; block < ma->numBlocks; block++) {
		int first = ma->blocks[block].first;
		int last = ma->blocks[block].last;
		for (chain=0; chain<ma->numChains; chain++) {
			int find = first;
			int index;
			while (find <= last && !ma->residues[chain].res[find].exists) find++;
			if (find > last) break;
			index = ma->residues[chain].res[find].index;
			while (indices[chain] < index) {
				for (chain2=0; chain2<ma->numChains; chain2++) {
					if (chain2 != chain) {
						AddCharacter(out->seqs[chain2], '-');
					}
					else {
						AddCharacter(out->seqs[chain2], ma->chains[chain2]->res[indices[chain2]++].residue);
					}
				}
			}
		}
		while (first <= last) {
			for (chain=0; chain<ma->numChains; chain++) {
				if (ma->residues[chain].res[first].exists) {
					AddCharacter(out->seqs[chain], ma->chains[chain]->res[indices[chain]++].residue);
				}
				else {
					AddCharacter(out->seqs[chain], '-');
				}
			}
			first++;
		}
	}
	for (chain=0; chain<ma->numChains; chain++) {
		while (indices[chain] < ma->chains[chain]->length) {
			for (chain2=0; chain2<ma->numChains; chain2++) {
				if (chain2 != chain) {
					AddCharacter(out->seqs[chain2], '-');
				}
				else {
					AddCharacter(out->seqs[chain2], ma->chains[chain2]->res[indices[chain2]++].residue);
				}
			}
		}
	}
	free(indices);
	Reorder(out);
	return out;
}

int OutputAlignmentBentPDB(MultipleAlignment *ma, const char *fileName, int *order) {
	Atom **atoms = (Atom**) malloc(sizeof(Atom*)*ma->numChains);
	int i;
	int result = 0;
	for (i=0; i<ma->numChains; i++) {
		atoms[i] = (Atom*) memdup(ma->chains[i]->pdb->atoms, sizeof(Atom) * ma->chains[i]->pdb->numAtoms);
	}
	TransformPDB(ma);
	result = OutputAlignmentPDB(ma, fileName, order);
	for (i=0; i<ma->numChains; i++) {
		free(ma->chains[i]->pdb->atoms);
		ma->chains[i]->pdb->atoms = atoms[i];
	}
	free(atoms);
	return result;
}

int OutputAlignmentPDB(MultipleAlignment *ma, const char *filename, int *order) {
	int i, j;
	PDBChain **tempChains;
	FILE *out;
	if (!order) return 0;
	tempChains = (PDBChain**) malloc(sizeof(PDBChain*) * ma->numChains);
	out = fopen(filename, "wb");
	if (!out) return 0;

	for (i=0; i < ma->numChains; i++) {
		for (j=0; j<ma->numChains; j++) {
			if (order[i] == ma->chains[j]->id) {
				tempChains[i] = ma->chains[j]->pdb;
				break;
			}
		}
	}

	OutputChains(tempChains, ma->numChains, out, 1, 1, 0, 0);

	free(tempChains);

	fclose(out);
	return 1;
}


int OutputAlignmentRasmol(MultipleAlignment *ma, const char *filename) {
	static const char *colors[] = {"red", "blue", "green", "purple", "cyan", "orange", "yellow", "greenblue", "magenta", "[255,171,187]"};
	FILE *out = fopen(filename, "wb");
	int i, j, k;

	if (!out) return 0;

	fprintf(out, "select all\ncartoon\n");

	for (j=0; j<ma->numChains; j++) {
		int count = 0;
		fprintf(out, "select ");
		for (i = 0; i < ma->numBlocks; i++) {
			int res1 = ma->blocks[i].first;
			int res2 = ma->blocks[i].last;
			char char1[2] = {' ', '\0'};
			char char2[2] = {' ', '\0'};

			for (k=0; k<ma->numChains; k++) {
				int l;
				for (l=res1; l<=res2; l++) {
					if (!ma->residues[k].res[l].exists) break;
				}
				if (l <= res2) break;
			}
			if (k < ma->numChains) {
				continue;
			}
			if (count) fputc(',', out);

			res1 = ma->residues[j].res[res1].index;
			res2 = ma->residues[j].res[res2].index;

			char1[0] = ma->chains[j]->pdb->residues[res1].resLetter;
			res1 = ma->chains[j]->pdb->residues[res1].resNum;
			char2[0] = ma->chains[j]->pdb->residues[res2].resLetter;
			res2 = ma->chains[j]->pdb->residues[res2].resNum;

			if (char1[0] == ' ') char1[0] = 0;
			if (char2[0] == ' ') char2[0] = 0;
			if (ma->chains[j]->pdb->chainName != ' ')
				fprintf(out, "%i%s-%i%s:%c", res1, char1, res2, char2, ma->chains[j]->pdb->chainName);
			else
				fprintf(out, "%i%s-%i%s", res1, char1, res2, char2);
			count ++;
			if (count == 10) {
				fprintf(out, "\ncolor %s\n", colors[j%(sizeof(colors)/sizeof(char*))]);
				fprintf(out, "select ");
				count = 0;
			}
		}
		fprintf(out, "\ncolor %s\n", colors[j%(sizeof(colors)/sizeof(char*))]);
	}
	fprintf(out, "select all\n");
	fclose(out);
	return 1;
}

int OutputAlignmentSpdbv(MultipleAlignment *ma, const char *filename) {
	static const char *colors[] = {"red", "blue", "green", "purple", "cyan", "orange", "yellow", "<0.18,0.54,0.34>", "<1.0,0.0,1.0>", "<1.0,0.67,0.73>"};
	FILE *out = fopen(filename, "wb");
	int i, j, k;

	if (!out) return 0;

	fprintf(out, "please do\ncolor ribbon,res,side of all in grey;\n");

	for (j=0; j<ma->numChains; j++) {
		int output = 0;
		int count = 0;
		const char *color = colors[j%(sizeof(colors)/sizeof(char*))];
		fprintf(out, "color ribbon,res,side of chain \"%c\" and num ", ma->chains[j]->pdb->chainName);
		for (i = 0; i < ma->numBlocks; i++) {
			int res1 = ma->blocks[i].first;
			int res2 = ma->blocks[i].last;
			for (k=0; k<ma->numChains; k++) {
				int l;
				for (l=res1; l<=res2; l++) {
					if (!ma->residues[k].res[l].exists) break;
				}
				if (l <= res2) break;
			}
			if (k < ma->numChains) {
				continue;
			}
			res1 = ma->residues[j].res[res1].index;
			res2 = ma->residues[j].res[res2].index;
			for (k = res1; k <= res2; k++) {
				if (output) fputc(',', out);
				output = 1;
				fprintf(out, "%i", ma->chains[j]->pdb->residues[k].resNum);
			}
			count ++;
			if (count == 30) {
				if (color[0] == '<') {
					fprintf(out, " by");
				}
				else {
					fprintf(out, " in");
				}
				fprintf(out, " %s;\n", color);
				fprintf(out, "color ribbon,res,side of chain \"%c\" and num ", ma->chains[j]->pdb->chainName);
				count = 0;
				output = 0;
			}
		}
		if (color[0] == '<') {
			fprintf(out, " by");
		}
		else {
			fprintf(out, " in");
		}
		fprintf(out, " %s;\n", color);
	}
	fprintf(out, "thank you\n");
	fclose(out);
	return 1;
}

int OutputAlignmentFasta(MultipleAlignment *ma, const char *filename) {
	FILE *out = fopen(filename, "wb");
	SequenceAlignment *sa;
	int result = 0;
	if (out) {
		sa = GetSequenceAlignment(ma);
		result = WriteSequenceAlignmentFASTA(out, sa);
		FreeSequenceAlignment(sa);
		fclose(out);
	}
	return result;
}
int OutputAlignmentPir(MultipleAlignment *ma, const char *filename) {
	FILE *out = fopen(filename, "wb");
	SequenceAlignment *sa;
	int result = 0;
	if (out) {
		sa = GetSequenceAlignment(ma);
		result = WriteSequenceAlignmentPIR(out, sa);
		FreeSequenceAlignment(sa);
		fclose(out);
	}
	return result;
}

int OutputAlignmentMsf(MultipleAlignment *ma, char *filename) {
	FILE *out = fopen(filename, "wb");
	SequenceAlignment *sa;
	int i, j, k;
	if (out) {
		char t[100];
		time_t now = time(0);
		char ** names = (char**) malloc(sizeof(char*) * ma->numChains);
		char *period = strrchr(filename, '.');
		sa = GetSequenceAlignment(ma);
		*period = 0;
		fprintf(out, "MSF of: %s from:    1 to: %4i\n", filename, sa->seqs[0]->length);
		*period = '.';

		strftime(t, 100, "%d-%b-%y %X", localtime(&now));
		/* asctime is not thread-safe, and Windows has no thread-safe version of it. */
		fprintf(out, "%s MSF: %4i  Type: P %s  Check: 9999  ..\n\n", filename, sa->seqs[0]->length, t);

		for (i=0; i<ma->numChains; i++) {
			char *id = ma->chains[i]->pdb->idString;
			int len = strlen(id);
			if (len > 13) len = 12;
			names[i] = (char*) calloc(16, sizeof(char));
			memset(names[i], ' ', 13 * sizeof(char));
			memcpy(names[i], id, len);
		}
		for (i=ma->numChains-1; i>=0; i--) {
			char *id = ma->chains[i]->pdb->idString;
			int chainCount = 0;
			for (j=i-1; j>=0; j--) {
				if (!strcmp(names[i], names[j])) {
					names[j][14] = 1;
					names[i][14] = 1;
					if (ma->chains[i]->pdb->chainName == ma->chains[j]->pdb->chainName) {
						chainCount++;
						names[j][15] = 1;
						names[i][15] = 1;
					}
				}
			}
			if (names[i][15] || names[i][14]) {
				int len = 13;
				char *end = strchr(id, ':');
				char chain = ma->chains[i]->pdb->chainName;
				if (!chain || chain == 0) chain = '_';
				if (!end) end = strchr(id, ' ');
				if (end) len = (int)(end-id);
				if (names[i][15]) {
					if (len > 8) len = 8;
					if (chainCount > 100) len --;
					if (chainCount > 1000) len --;
					if (chainCount > 10000) len --;
					chainCount %= 100000;
					sprintf(names[i] + len, ":%c%i", ma->chains[i]->pdb->chainName, chainCount);
				}
				else {
					if (len > 10) len = 10;
					sprintf(names[i] + len, ":%c", ma->chains[i]->pdb->chainName);
				}
			}
		}

		for (i=0; i<ma->numChains; i++) {
			int count = 0;
			int check = 0;

			for (j=0; j<sa->seqs[i]->length; j++) {
				count++;
				check += count * UCASE(sa->seqs[i]->seq[j]);
				count = count%57;
			}
			check %= 10000;
			fprintf(out, "Name:  %-13s  Len:  %4i  Check:  %4i  Weight:  1.00\n", names[i], sa->seqs[i]->length, check);
		}
		fprintf(out, "\n//\n\n");

		for (i=0; i < sa->seqs[0]->length; i+=50) {
			for (j=0; j < ma->numChains; j++) {
				fprintf(out, "%-13s", names[j]);
				for (k=i; k<i+50 && k < sa->seqs[j]->length; k++) {
					char c = sa->seqs[j]->seq[k];
					if (c == '-') c = '.';
					fputc(c, out);
					if (k % 10 == 9 && k<i+49)
						fputc(' ', out);
				}
				fputc('\n', out);
			}
			fputc('\n', out);
		}
		for (i=0; i<ma->numChains; i++) {
			free(names[i]);
		}
		free(names);
		FreeSequenceAlignment(sa);
		fclose(out);
		return 1;
	}
	return 0;
}

void OutputAlignmentOrder(MultipleAlignment *ma, FILE *file, const AssemblyNode *node) {
	int i;
	if (node->id >= 0) {
		for (i=0; i<ma->numChains; i++) {
			if (ma->chains[i]->id == node->id) {
				fprintf(file, "%s", ma->chains[i]->pdb->idString);
				break;
			}
		}
	}
	else {
		fprintf(file, "(");
		OutputAlignmentOrder(ma, file, node->left);
		fprintf(file, ", ");
		OutputAlignmentOrder(ma, file, node->right);
		fprintf(file, ")");
	}
}

int OutputAlignmentText(MultipleAlignment *ma, const char *filename, int ref) {
	FILE *out = fopen(filename, "wb");
	int i;
	int nameLen = 0;
	int coreCount = 0;
	char format[30];
	SequenceAlignment *sa;
	if (!out) return 0;
	OutputAlignmentOrder(ma, out, ma->order->root);
	fprintf(out, "\n");
	for (i=0; i<ma->numChains; i++) {
		int len = (int)strlen(ma->chains[i]->pdb->idString);
		if (len > nameLen) nameLen = len;
	}

	for (i=0; i<ma->numResidues;i++) {
		int chain;
		for (chain=0; chain<ma->numChains; chain++) {
			if (!ma->residues[chain].res[i].exists) break;
		}
		coreCount += (chain == ma->numChains);
	}

	fprintf(out, "Core Residues: %i\n", coreCount);
	fprintf(out, "Core RMSD: %0.3lf\n", ma->rmsd);
	fprintf(out, "Raw Score: %0.3lf\n", ma->score);
	if (ma->numChains == 2 && ma->pvalue >= 0) {
		if (ma->pvalue >= 0.01)
			fprintf(out, "P-value: %0.4lf\n", ma->pvalue);
		else
			fprintf(out, "P-value: %0.4e\n", ma->pvalue);
	}
	fprintf(out, "Reference structure: %s\n\n", ma->chains[ref]->pdb->idString);
	sprintf(format, "%%-%is   ", nameLen);

	sa = GetSequenceAlignment(ma);
	for (i = 0; i<sa->seqs[0]->length; i+=60) {
		int end = i + 60;
		int j;

		if (end > sa->seqs[0]->length) end = sa->seqs[0]->length;
		for (j=0; j<ma->numChains; j++) {
			char *seq = sa->seqs[j]->seq;
			int p = 0;
			int k;
			fprintf(out, format, ma->chains[j]->pdb->idString);
			for (k = 0; k < i; k++) p += (seq[k] != '-');

			for (k = i; k<end; k++) {
				fprintf(out, "%c", seq[k]);
				p += (seq[k] != '-');
			}
			if (p) p--;
			fprintf(out, " %5i%c (%c)\n", ma->chains[j]->pdb->residues[p].resNum, ma->chains[j]->pdb->residues[p].resLetter, ma->chains[j]->pdb->chainName);
		}
		fprintf(out, "\n");
	}

	FreeSequenceAlignment(sa);
	fclose(out);
	return 1;
}
