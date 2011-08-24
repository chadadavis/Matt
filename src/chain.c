/* chain.c -- Contains functions for creating and manipulating pdb chains.
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

#include "chain.h"
#include <math.h>

#define Qconst (-27888.0f)
#define MaxHDist2 (9.0f*9.0f)
#define MinHDist2 (0.5f*0.5f)

void RelabelChains(PDBChain **chains, int numChains) {
	static const char names[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890~!@#$%^&()-+=_<>/'\\?\"[]{}abcdefghijklmnopqrstuvwxyz";
	int i;
	for (i=0; i<numChains; i++) {
		chains[i]->chainName = names[i%(sizeof(names)/sizeof(char))];
	}
}

void InitResidue(Residue *res, char *name) {
	/* Special value for not aligned to anything. */
	res->resNum = -10000000;
	res->resLetter = 0;
	strcpy(res->resName, name);
	res->firstAtom = res->lastAtom = -1;
}

void AddResidue(PDBChain *chain, char *resName) {
	if (chain->length % 16 == 0)
		chain->residues = (Residue*) realloc(chain->residues, sizeof(Residue)*(chain->length+16));
	InitResidue(&chain->residues[chain->length], resName);
	chain->length++;
}

int GetAtomIndex(PDBChain *chain, int pos, char*desc) {
	int i;
	if (pos >= chain->length || pos < 0 || chain->residues[pos].firstAtom==-1) return -1;
	for (i=chain->residues[pos].firstAtom; i <= chain->residues[pos].lastAtom; i++) {
		if (strcmp(desc, chain->atoms[i].name) == 0)
			return i;
	}
	return -1;
}

Atom * GetAtom(PDBChain *chain, int pos, char *desc) {
	int i = GetAtomIndex(chain, pos, desc);
	if (i==-1) return 0;
	return &chain->atoms[i];
}

Residue * GetResidue(PDBChain *chain, int pos) {
	if (pos >= chain->length || pos<0) return 0;
	return &chain->residues[pos];
}

void RotateChain(PDBChain *chain, const Matrix *m) {
	int i;
	for (i=chain->numAtoms-1; i>=0; i--) {
		Vector v = chain->atoms[i].pos;
		transformVect(&chain->atoms[i].pos, m, &v);
	}
}


struct FixAlignmentInfo {
	char prevGapped;
	char prevSkipped;
	int numGaps;
	int numSkips;
	int skipCount;
};

struct RES_TO_SEQRES_ALIGN_INFO {
	int alignedTo;
	int prevResIndex;
	int prevIndex;
	int prevSkip;
	int gapCount;
};
typedef struct RES_TO_SEQRES_ALIGN_INFO ResToSeqresAlignInfo;

struct RES_TO_SEQRES_ALIGN_LIST {
	int count;
	ResToSeqresAlignInfo *align;
};
typedef struct RES_TO_SEQRES_ALIGN_LIST ResToSeqresAlignList;

struct RES_ALIGN_INFO {
	ResToSeqresAlignList *list;
};
typedef struct RES_ALIGN_INFO ResAlignInfo;


void RemoveHetsWithoutCAs(PDBChain *chain) {
	int res, atomsToRemove = 0;
	for (res = 0; res<chain->numTempResidues; res++) {
		if (chain->tempResidues[res].het == 0 || chain->tempResidues[res].alphaCarbon) break;
		while (atomsToRemove < chain->numAtoms &&
			chain->atoms[atomsToRemove].resLetter == chain->tempResidues[res].resLetter &&
			chain->atoms[atomsToRemove].res == chain->tempResidues[res].res) {
				atomsToRemove++;
		}
	}
	if (res) {
		memmove(chain->tempResidues, chain->tempResidues+res, (chain->numTempResidues-res) * sizeof(TempResidue));
		chain->numTempResidues-= res;
		memmove(chain->atoms, chain->atoms+atomsToRemove, (chain->numAtoms-atomsToRemove) * sizeof(Atom));
		chain->numAtoms-=atomsToRemove;
	}

	for (res = chain->numTempResidues-1; res >= 0; res--) {
		if (chain->tempResidues[res].het == 0 || chain->tempResidues[res].alphaCarbon) break;
		while (chain->numAtoms &&
			chain->atoms[chain->numAtoms-1].resLetter == chain->tempResidues[res].resLetter &&
			chain->atoms[chain->numAtoms-1].res == chain->tempResidues[res].res) {
				chain->numAtoms--;
		}
	}
	chain->numTempResidues = res+1;
}

void SeqresResAlign(PDBChain *chain, int maxSkip) {
	int res, seqres, skip, i;
	int finalSkip;
	int pos = 0;
	int removedAtoms = 0;
	ResAlignInfo* info = (ResAlignInfo*) malloc(sizeof(ResAlignInfo) * chain->numTempResidues);
	info[0].list = (ResToSeqresAlignList *) calloc(chain->numTempResidues * (maxSkip+1), sizeof(ResToSeqresAlignInfo));

	RemoveHetsWithoutCAs(chain);
	if (maxSkip > chain->numTempResidues) maxSkip = chain->numTempResidues;
	if (!chain->numTempResidues) return;

	/* As many or more gaps as residues is bad. In that case, code to ignore SEQRES entries should work fine. */
	if (chain->numTempResidues-1 < maxSkip) maxSkip = chain->numTempResidues-1;

	for (res=1; res < chain->numTempResidues; res++) {
		info[res].list = info[res-1].list + maxSkip+1;
	}
	for (res=0; res < chain->numTempResidues; res++) {
		for (seqres = 0; seqres<chain->length; seqres++) {
			ResToSeqresAlignInfo best;
			if (strcmp(chain->tempResidues[res].resName, chain->residues[seqres].resName)) continue;
			best.alignedTo = seqres;
			best.gapCount = chain->length;

			for (skip=0; skip<maxSkip; skip++) {
				int delta;
				if (skip == res) {
					best.gapCount = 0;
					best.prevResIndex = -1;
					best.prevSkip = -1;
					best.prevIndex = -1;
					break;
				}
				for (delta = 1; delta <= 1+skip; delta++) {
					int res2 = res-delta;
					int i = skip-delta+1;
					int index;
					for (index = 0; index < info[res2].list[i].count; index++) {
						int seqres2 = info[res2].list[i].align[index].alignedTo;
						int gapCount = info[res2].list[i].align[index].gapCount;
						if (seqres2 >= seqres) break;
						if (seqres2+1 != seqres) {
							if (seqres - seqres2 != chain->tempResidues[res].res - chain->tempResidues[res2].res)
								gapCount++;
						}
						if (gapCount < best.gapCount) {
							best.gapCount = gapCount;
							best.prevResIndex = res2;
							best.prevSkip = i;
							best.prevIndex = index;
						}
					}
				}
				if (best.gapCount < chain->length) break;
			}
			if (best.gapCount == chain->length) {
				continue;
			}
			if ((info[res].list[skip].count & 0xFF) == 0) {
				info[res].list[skip].align = (ResToSeqresAlignInfo*)realloc(info[res].list[skip].align, sizeof(ResToSeqresAlignInfo) * (info[res].list[skip].count + 0x100));
			}
			info[res].list[skip].align[info[res].list[skip].count++] = best;
		}
		for (skip = 0; skip <= maxSkip; skip++) {
			if (info[res].list[skip].align) {
				info[res].list[skip].align = (ResToSeqresAlignInfo*)realloc(info[res].list[skip].align, sizeof(ResToSeqresAlignInfo) * info[res].list[skip].count);
			}
		}
	}

	for (skip=0; skip<=maxSkip; skip++) {
		finalSkip = 0;
		res = chain->numTempResidues-1;
		while (!info[res].list[skip].count && res && chain->tempResidues[res].het == 1 && !chain->tempResidues[res].alphaCarbon) res--;
		if (info[res].list[skip].count) break;
		for (finalSkip = 1; finalSkip <= skip; finalSkip++) {
			if (info[res-finalSkip].list[skip-finalSkip].count) break;
		}
		if (finalSkip <= skip) break;
	}

	if (skip > maxSkip) {
		res = chain->numTempResidues-1;
		while (!info[res].list[skip].count && res && chain->tempResidues[res].het == 1 && !chain->tempResidues[res].alphaCarbon) res--;
	}
	for (i=0; i<=res; i++) {
		while (pos < chain->numAtoms) {
			if (chain->atoms[pos].res != chain->tempResidues[i].res ||
				chain->atoms[pos].resLetter != chain->tempResidues[i].resLetter) break;
			pos++;
		}
	}
	/* Remove extra residues and their atoms. */
	chain->numTempResidues = res+1;
	chain->numAtoms = pos;
	chain->atoms = (Atom*) realloc(chain->atoms, sizeof(Atom)*chain->numAtoms);

	if (skip <= maxSkip) {
		int pick = 0;
		res -= finalSkip;
		skip-=finalSkip;
		for (i=1; i<info[res].list[skip].count; i++) {
			if (info[res].list[skip].align[pick].gapCount > info[res].list[skip].align[i].gapCount)
				pick = i;
		}
		chain->tempResidues[res].finalIndex = info[res].list[skip].align[pick].alignedTo;
		while (res > 0) {
			ResToSeqresAlignInfo *inf = &info[res].list[skip].align[pick];
			res = inf->prevResIndex;
			skip = inf->prevSkip;
			pick = inf->prevIndex;
			if (res < 0) break;
			chain->tempResidues[res].finalIndex = info[res].list[skip].align[pick].alignedTo;
		}
		for (res = 0; res < chain->numTempResidues; res++) {
			if (chain->tempResidues[res].finalIndex == -1) {
				int start, start2, end, end2, pos;
				int range;
				int range2;
				int endLen;
				start = res-1;
				pos = res-1;
				while (pos >= start-1 && pos > 0) {
					int delta = chain->tempResidues[pos].finalIndex - chain->tempResidues[pos-1].finalIndex;
					if (delta != 1 && delta != chain->tempResidues[pos].res - chain->tempResidues[pos-1].res) {
						start = pos-1;
						/* Run off end.  If only single mismatch near ends, don't delete leading/trailing
						 * seqres entries.  Otherwise, delete them.
						 */
						if (!start) start = -1;
					}
					pos--;
				}
				end = res+1;
				pos = res+1;
				while (pos <= end+1 && pos < chain->numTempResidues-1) {
					int delta = chain->tempResidues[pos+1].finalIndex - chain->tempResidues[pos].finalIndex;
					if (delta != 1 && delta != chain->tempResidues[pos+1].res - chain->tempResidues[pos].res) {
						end = pos+1;
						if (end == chain->numTempResidues-1) end++;
					}
					pos++;
				}
				/* Want to replace from start2 to end2 in SEQRES with start to end from ATOM/HETATM entries.
				 */
				if (start < 0) {
					start2 = 0;
					start = 0;
				}
				else start2 = chain->tempResidues[start].finalIndex;
				if (end >= chain->numTempResidues) {
					end = chain->numTempResidues-1;
					end2 = chain->length-1;
				}
				else end2 = chain->tempResidues[end].finalIndex;
				range = end-start;
				range2 = end2-start2;
				for (res = end+1; res < chain->numTempResidues; res++) {
					if (chain->tempResidues[res].finalIndex != -1) {
						chain->tempResidues[res].finalIndex += range-range2;
					}
				}
				endLen = (chain->length-end2-1)*sizeof(Residue);
				if (range2 < range) {
					chain->residues = (Residue*) realloc(chain->residues, sizeof(Residue) * (chain->length += range - range2));
					memmove(chain->residues+start2+range+1, chain->residues+end2+1, endLen);
				}
				else if (range2 > range) {
					memmove(chain->residues+start2+range+1, chain->residues+end2+1, endLen);
					chain->residues = (Residue*) realloc(chain->residues, sizeof(Residue) * (chain->length += range - range2));
				}
				for (res = start; res <= end; res++) {
					chain->tempResidues[res].finalIndex = start2 + res-start;
					InitResidue(&chain->residues[start2+res-start], chain->tempResidues[res].resName);
				}
			}
		}
		for (res=0; res < chain->numTempResidues; res++) {
			int res2 = chain->tempResidues[res].finalIndex;
			chain->residues[res2].resNum = chain->tempResidues[res].res;
			chain->residues[res2].resLetter = chain->tempResidues[res].resLetter;
		}

		/* Number everything before first ATOM record. */
		for (res = chain->tempResidues[0].finalIndex-1; res >= 0; res--) {
			chain->residues[res].resNum = chain->residues[res+1].resNum-1;
			chain->residues[res].resLetter = ' ';
		}
		/* Number everything after last ATOM record. */
		for (res = chain->tempResidues[chain->numTempResidues-1].finalIndex+1; res < chain->length; res++) {
			chain->residues[res].resNum = chain->residues[res-1].resNum+1;
			chain->residues[res].resLetter = ' ';
		}
		/* Number everything in gaps.  Number from lower number.
		 * If atom numbering doesn't fit with SEQRES (That is, more residues in a gap
		 * than allowed numbers), then leave the other residues basically random low numbers.
		 * Used to add insertion codes, but it didn't do much good.
		 */
		for (res = 1; res<chain->length; res++) {
			int res2 = res;
			int first = 1;
			if (chain->residues[res].resNum != -10000000) continue;
			while (res2 < chain->length && chain->residues[res2].resNum == -10000000) {
				chain->residues[res2].resNum = 1 + chain->residues[res2-1].resNum;
				chain->residues[res2].resLetter = ' ';
				res2++;
			}
			if (res2 == chain->length) break;
			/* Make sure that if residues have to have duplicate numbers, they won't have greater numbers
			 * than residues later in the sequence.  Note that this won't be true if numbering isn't
			 * monotinically increasing.
			 */
			while (res2>res && chain->residues[res2-1].resNum >= chain->residues[res2].resNum) {
				char letter = chain->residues[res2].resLetter;
				int num = chain->residues[res2].resNum;
				if (letter <= 'A' || letter > 'Z' || first) {
					letter = 'Z';
					num = chain->residues[res2].resNum-1;
					first = 0;
				}
				else letter --;
				chain->residues[res2-1].resLetter = letter;
				chain->residues[res2-1].resNum = num;
				res2--;
			}
			res = res2;
		}
	}
	else {
		free(chain->residues);
		chain->residues = 0;
		chain->length = 0;
	}

	for (res=0; res<chain->numTempResidues; res++) {
		for (skip=0; skip<=maxSkip; skip++) {
			free(info[res].list[skip].align);
		}
	}
	free(info[0].list);
	free(info);
}


void Terminate(PDBChain *chain, int maxGaps, int ter) {
	int i;
	if (chain->terminated) return;
 	chain->terminated = 1;

	/* If no ATOM records, do nothing. */
	if (chain->numAtoms) {
		int pos = 0;
		if (chain->chainName == 'H') {
			chain=chain;
		}
		if (chain->length)
			SeqresResAlign(chain, maxGaps);

		/* If no SEQRES records, populate them from tempResidues.
		 * This state may or may not be caused by SeqresAlign.
		 */
		if (chain->length == 0) {
			for (i=0; i<chain->numTempResidues; i++) {
				AddResidue(chain, chain->tempResidues[i].resName);
				chain->residues[i].resNum = chain->tempResidues[i].res;
				chain->residues[i].resLetter = chain->tempResidues[i].resLetter;
			}
		}
		free(chain->tempResidues);
		chain->numTempResidues = 0;
		chain->tempResidues = 0;

		/* Match up atom records with residues. */
		for (i=0; i<chain->numAtoms; i++) {
			while (!(chain->residues[pos].resNum == chain->atoms[i].res && chain->residues[pos].resLetter == chain->atoms[i].resLetter)) {
				pos++;
			}

			chain->atoms[i].resIndex = pos;
			if (chain->residues[pos].firstAtom == -1) chain->residues[pos].firstAtom = i;
			chain->residues[pos].lastAtom = i;
		}
	}
	for (i=0; i<chain->length; i++) {
		AddCharacter(chain->seq, LookupName(chain->residues[i].resName, NAME_BOTH));
	}
}

void End(PDBChain *chain, int maxGaps) {
	/* Require explicit terminate (checked elsewhere) or seqres
	 * or at least one atom (Not HETATM) record to consider a chain an actual chain
	 * otherwise, assumed to be just a bunch of molelcules.
	 */
	if (!chain->length && !chain->tempAtoms) chain->terminated = 1;

	if (!chain->terminated) {
		Terminate(chain, maxGaps, 0);
	}

	free(chain->tempResidues);
	chain->tempResidues = 0;
	chain->numTempResidues = 0;
}

void TerminateChainOutput(PDBChain *chain, FILE *out, int *index, int res, char c) {
	int num = chain->residues[res].resNum;
	fprintf(out, "TER   %5i      %-3s %c%4i%c%53s\n", index[0]++, chain->residues[res].resName, c, num, chain->residues[res].resLetter, "");
}

int OutputChain(PDBChain *chain, FILE *out, int seqres, int het, int *index) {
	int i;
	int terminated = 0;
	int lastAtom = -1;
	if (!out) return 0;
	if (!chain->numAtoms) return 1;
	for (i=0; i<chain->numAtoms; i++) {
		char *type;
		char charge[3] = "  ";
		Residue * res = &chain->residues[chain->atoms[i].resIndex];
		if (!chain->length) {
			if (!terminated) {
				if (!i) return 0;
				TerminateChainOutput(chain, out, index, chain->atoms[i-1].resIndex, chain->chainName);
				terminated = 1;
			}
		}

		if (chain->seq->seq[chain->atoms[i].resIndex] != lookupTable['X'-'A'] || !strcmp(chain->residues[chain->atoms[i].resIndex].resName, "UNK")) {
			type = "ATOM  ";
		}
		else {
			type = "HETATM";
		}
		if (chain->atoms[i].charge>0) {
			charge[0] = '0' + charge[1];
			charge[1] = '+';
		}
		else if (chain->atoms[i].charge<0) {
			charge[0] = '0' - charge[1];
			charge[1] = '-';
		}

		fprintf(out, "%s%5i %-4s %-3s %c%4i%c   %8.3lf%8.3lf%8.3lf%6.2f%6.2f          %c%c%s\n", type, index[0]++, chain->atoms[i].name, res->resName, chain->chainName, res->resNum, res->resLetter, chain->atoms[i].pos.x, chain->atoms[i].pos.y, chain->atoms[i].pos.z, chain->atoms[i].occupancy, chain->atoms[i].tempFactor, chain->atoms[i].element[0], chain->atoms[i].element[1], charge);
	}
	if (!terminated) TerminateChainOutput(chain, out, index, chain->atoms[i-1].resIndex, chain->chainName);

	return 1;
}

void OutputChainAlphas(PDBChain *chain, FILE *out, int het, int *numHelices, char **helixIDs) {
	int i;
	for (i=0; i < chain->numAlphaHelices; i++) {
		int len = 0;
		AlphaHelix *helix = &chain->alphaHelices[i];
		Residue *start = &chain->residues[helix->start];
		Residue *end = &chain->residues[helix->start + chain->alphaHelices[i].length - 1];
		char id[4];
		if (helix->id[0] && (!*helixIDs || !strstr(*helixIDs, helix->id))) {
			strcpy(id, helix->id);
		}
		else {
			int i = *numHelices;
			if (!*helixIDs || !strstr(*helixIDs, id)) {
				while (i) {
					sprintf(id, "%03i", i-1);
					if (strstr(*helixIDs, id)) break;
					i--;
				}
				sprintf(id, "%03i", i);
			}
			else {
				i=0;
				sprintf(id, "%03i", i);
				while (strstr(*helixIDs, id)) {
					i++;
					sprintf(id, "%03i", i);
				}
			}
		}
		fprintf(out, "HELIX  %3i %3s %s %c %4i%c %s %c %4i%c%2i%30s %5i    \n",
			++numHelices[0], id, start->resName, chain->chainName, start->resNum, start->resLetter, end->resName, chain->chainName, end->resNum, end->resLetter, helix->type, helix->comment, helix->length);
		if (*helixIDs) {
			len = (int)strlen(*helixIDs);
			helixIDs[0][len++] = '\n';
		}
		*helixIDs = (char*) realloc(*helixIDs, len+4);
		strcpy(helixIDs[0]+len, id);
	}
}

void OutputChainSheets(PDBChain *chain, FILE *out, int het, int *numSheets, char **sheetIDs, int *sheetLines) {
	int i, j;
	int s, e;
	Residue *start, *end, *current, *prev;
	for (i=0; i < chain->numBetaSheets; i++) {
		int len = 0;
		BetaSheet *b = &chain->betaSheets[i];
		char id[4];
		if (b->id[0] && (!*sheetIDs || !strstr(*sheetIDs, b->id))) {
			strcpy(id, b->id);
		}
		else {
			int i = *numSheets;
			if (!*sheetIDs || !strstr(*sheetIDs, id)) {
				while (i) {
					sprintf(id, "%03i", i-1);
					if (strstr(*sheetIDs, id)) break;
					i--;
				}
				sprintf(id, "%03i", i);
			}
			else {
				i=0;
				sprintf(id, "%03i", i);
				while (strstr(*sheetIDs, id)) {
					i++;
					sprintf(id, "%03i", i);
				}
			}
		}
		if (*sheetIDs) {
			len = (int)strlen(*sheetIDs);
			sheetIDs[0][len++] = '\n';
		}
		*sheetIDs = (char*) realloc(*sheetIDs, len+4);
		strcpy(sheetIDs[0]+len, id);
		s = b->strands[0].res1;
		e = b->strands[0].res1 + b->strands[0].length - 1;
		start = &chain->residues[s];
		end = &chain->residues[e];
		sheetLines[0]++;
		fprintf(out, "SHEET  %3i %3s%2i %s %c%4i%c %s %c%4i%c 0%30s          \n", 1, id, b->numStrands+1, start->resName, chain->chainName, start->resNum, start->resLetter, end->resName, chain->chainName, end->resNum, end->resLetter, "");
		for (j = 0; j < b->numStrands; j++) {
			char c1 = 'N', c2 = 'O';
			HydrogenBond *hbond;
			if (b->strands[j].direction == 1) {
				s = b->strands[j].res2;
				e = b->strands[j].res2 + b->strands[j].length - 1;
			}
			else {
				e = b->strands[j].res2;
				s = b->strands[j].res2 - b->strands[j].length + 1;
			}
			if (j < b->numStrands-1) {
				if (s > b->strands[j+1].res1) s = b->strands[j+1].res1;
				if (e < b->strands[j+1].res1 + b->strands[j+1].length-1) e = b->strands[j+1].res1 + b->strands[j+1].length-1;
			}
			start = &chain->residues[s];
			end = &chain->residues[e];
			current = &chain->residues[b->strands[j].res2];
			prev = &chain->residues[b->strands[j].res1];
			if (b->strands[j].direction == -1) {
				if (!AreHbondedAsymmetric(chain, b->strands[j].res1, b->strands[j].res2, &hbond)) {
					if (AreHbondedAsymmetric(chain, b->strands[j].res2, b->strands[j].res1, &hbond)) {
						c1 = 'O';
						c2 = 'N';
					}
					else if (AreHbondedAsymmetric(chain, 1+b->strands[j].res2, b->strands[j].res1-1, &hbond)) {
						current++;
						prev--;
					}
				}
				fprintf(out, "SHEET  %3i %3s%2i %s %c%4i%c %s %c%4i%c%2i  %c  %s %c%4i   %c  %s %c%4i%8s   \n", j+2, id, b->numStrands+1, start->resName, chain->chainName, start->resNum, start->resLetter, end->resName, chain->chainName, end->resNum, end->resLetter, b->strands[j].direction, c1, current->resName, chain->chainName, current->resNum, c2, prev->resName, chain->chainName, prev->resNum, "");
				sheetLines[0]++;
			}
			else {
				if (AreHbondedAsymmetric(chain, b->strands[j].res2+1, b->strands[j].res1, &hbond)) {
					current++;
				}
				else if (AreHbondedAsymmetric(chain, b->strands[j].res2+2, b->strands[j].res1+1, &hbond)) {
					current+=2;
					prev++;
				}
				fprintf(out, "SHEET  %3i %3s%2i %s %c%4i%c %s %c%4i%c%2i  %c  %s %c%4i%c  %c  %s %c%4i%c%8s  \n", j+2, id, b->numStrands+1, start->resName, chain->chainName, start->resNum, start->resLetter, end->resName, chain->chainName, end->resNum, end->resLetter, b->strands[j].direction, c1, current->resName, chain->chainName, current->resNum, current->resLetter, c2, prev->resName, chain->chainName, prev->resNum, prev->resLetter, "");
				sheetLines[0]++;
			}
		}
	}
}

void OutputChainSSbonds(PDBChain *chain, FILE *out, int het, int *numSSbonds) {
	int i;
	for (i=0; i<chain->numSSbonds; i++) {
		SSbond *ss = chain->ssbonds + i;
		Residue *res1 = &chain->residues[ss->res1];
		Residue *res2 = &chain->residues[ss->res2];
		if (het || (chain->seq->seq[ss->res1] < 20 && chain->seq->seq[ss->res2] < 20)) {
			fprintf(out, "SSBOND %3i %3s %c %4i%c   %3s %c %4i%c                       ", ++numSSbonds[0], res1->resName, chain->chainName, res1->resNum, res1->resLetter, res2->resName, chain->chainName, res2->resNum, res2->resLetter);
			if (ss->sym1) {
				fprintf(out, "%6i ", ss->sym1);
			}
			else fprintf(out, "       ");
			if (ss->sym2) {
				fprintf(out, "%6i", ss->sym2);
			}
			else fprintf(out, "      ");
			fprintf(out, "        \n");
		}
	}
}

void OutputChainCisPeps(PDBChain *chain, FILE *out, int het, int *numCisPeps, int model) {
	int i;
	for (i=0; i<chain->numCisPeps; i++) {
		CisPep *cispep = chain->cisPeps + i;
		Residue *res1 = &chain->residues[cispep->res1];
		Residue *res2 = &chain->residues[cispep->res2];
		if (het || (chain->seq->seq[cispep->res1] < 20 && chain->seq->seq[cispep->res2] < 20)) {
			fprintf(out, "CISPEP %3i %3s %c %4i%c   %3s %c %4i%c       %3i       %6.2lf   %18s\n", ++numCisPeps[0], res1->resName, chain->chainName, res1->resNum, res1->resLetter, res2->resName, chain->chainName, res2->resNum, res2->resLetter, model, cispep->angle, "");
		}
	}
}

int OutputChainFasta(PDBChain *chain, FILE *out, int het) {
	int count = 0;
	int lastAtom = -1;
	int i;
	if (!out) return 0;
	fprintf(out, ">%s\n", chain->idString);
	for (i=0; i<chain->numAtoms; i++) {
		Residue * res = &chain->residues[chain->atoms[i].resIndex];
		if (!het && (chain->seq->seq[chain->atoms[i].resIndex] == lookupTable['X'-'A'] && strcmp(chain->residues[chain->atoms[i].resIndex].resName, "UNK"))) {
			continue;
		}
		fputc(ShortNames[chain->seq->seq[chain->atoms[i].resIndex]], out);
		count ++;
		if (count && count %60 == 0) fputc('\n', out);
		while (i+1 < chain->numAtoms && chain->atoms[i].resIndex == chain->atoms[i+1].resIndex) i++;
	}
	if (!count || count%60) fputc('\n', out);
	return 1;
}

int OutputChainSeqres(PDBChain *chain, FILE *out, int het) {
	int i;
	int terminated = 0;
	int lastAtom = -1;
	int count = chain->length;
	int pos = 0;

	if (!out) return 0;
	if (!het) {
		count = 0;
		for (i=0; i<chain->length; i++) {
			if (chain->seq->seq[i] < 20) count++;
		}
	}
	for (i=0; i<chain->length; i++) {
		if (!het && chain->seq->seq[i] >= 20) continue;
		if (pos%13 == 0) {
			if (pos) fprintf(out, "         \n");
			fprintf(out, "SEQRES  %2i %c %4i  ", pos/13+1, chain->chainName, count);
		}
		fprintf(out, "%s ", chain->residues[i].resName);
		pos ++;
	}
	while (pos % 13) {
		fprintf(out, "    ");
		pos++;
	}
	if (pos) fprintf(out, "         \n");
	return 1;
}

void CleanupChain(PDBChain *chain) {
	int i;
	for (i=0; i<chain->numBetaSheets; i++) {
		free(chain->betaSheets[i].strands);
	}
	free(chain->betaSheets);
	FreeSequence(chain->seq);
	free(chain->cisPeps);
	free(chain->seqAdvs);
	free(chain->dbrefs);
	free(chain->ssbonds);
	free(chain->idString);
	free(chain->alphaHelices);
	free(chain->hbonds);
	free(chain->betaPairs);
	free(chain->residues);
	free(chain->atoms);
	free(chain->tempResidues);
	free(chain);
}

PDBChain *CreateChain(char *name, char chain) {
	PDBChain *out = (PDBChain*) calloc(1, sizeof(PDBChain));

	if (!name) {
		name = "Unnamed";
	}
	if (chain == ' ' || !chain) {
		out->idString = strdup(name);
	}
	else {
		out->idString = (char*) malloc(strlen(name)+3);
		sprintf(out->idString, "%s:%c", name, chain);
	}

	out->seq = CreateSequence(out->idString);
	out->chainName = chain;
	return out;
}

int GetResidueIndex(PDBChain *chain, int *position, int pdbPosition, char insertionCode) {
	int i;
	int min = 0, max = chain->length-1;
	/* normally use (min+max)/2, but guessing that the pdb numbering
	 * actually makes sense for the first guess seems a good idea.
	 */
	int current = pdbPosition - 1;
	if (current < 0 || current >= chain->length) current = (min+max/2);
	while (min <= max)
	{
		/* compare numbers, and adjust range if no match. */
		if (chain->residues[current].resNum > pdbPosition)
			max = current - 1;
		else if (chain->residues[current].resNum < pdbPosition)
			min = current + 1;
		/* Numbers match, now decide what to do. */

		/* Want first occurance of code */
		else if (insertionCode == -1) {
			if (current && chain->residues[current-1].resNum == pdbPosition) {
				max = current - 1;
			}
			else {
				*position = current;
				return 1;
			}
		}
		/* Otherwise compare insertion codes, just like we did the numbers, and adjust
		 * range if necessary.
		 */
		else if (chain->residues[current].resLetter > insertionCode)
			max = current - 1;
		else if (chain->residues[current].resLetter < insertionCode)
			min = current + 1;
		/* Otherwise, we're where we want to be. */
		else {
			if (current && chain->residues[current-1].resNum == pdbPosition && chain->residues[current-1].resLetter == insertionCode) {
				max = current - 1;
			}
			else {
				*position = current;
				return 1;
			}
		}
		current = (min+max)/2;
	}
	/* PDB numbering might not actually make sense...In which case, logical searches
	 * aren't guaranteed to work.  Also, insertion codes might not be in alphabetical
	 * order (Or ' ' codes might occur after lettered ones), in which case,
	 * the above search fails.
	 */
	for (i=0; i<chain->length; i++) {
		if (chain->residues[i].resNum == pdbPosition &&
			(chain->residues[i].resLetter == insertionCode ||
			 insertionCode == -1)) {
			*position = current;
			return 1;
		}
	}
	return GetResidueFirstIndex(chain, position, pdbPosition, insertionCode);
}

int GetResidueFirstIndexString(PDBChain *chain, int *position, char * rawString) {
	char code = rawString[4];
	int out;
	int index = atoi(rawString);
	rawString[4] = 0;
	out = GetResidueFirstIndex(chain, position, index, code);
	rawString[4] = code;
	return out;
}

int GetResidueFirstIndex(PDBChain *chain, int *position, int pdbPosition, char insertionCode) {
	/* PDB numbering might not actually make sense...In which case, logical searches
	 * aren't guaranteed to work.  Also, insertion codes might not be in alphabetical
	 * order (Or ' ' codes might occur after lettered ones), in which case,
	 * the above search fails.
	 */
	int i;
	if (insertionCode == 0) insertionCode = ' ';
	for (i=0; i<chain->length; i++) {
		if (chain->residues[i].resNum == pdbPosition &&
			(chain->residues[i].resLetter == insertionCode ||
			 insertionCode == -1)) {
			*position = i;
			return 1;
		}
	}
	return 0;
}

void CalcChainSecondary(PDBChain *chain) {
	if (!chain->secondaryCalculated) {
		CalcHbonds(chain);
		CalcBetaSheets(chain);
		CalcHelices(chain);
		chain->secondaryCalculated = 1;
		AssembleSheets(chain);
	}
}

void CalcHbonds(PDBChain *chain) {
	Atom *CA, *N, *O, *C, *C2, *O2, *CA2, *O3, *C3;

	double H_O2_dist, H_C2_dist, N_C2_dist, N_O2_dist;
	int i, j;

	free(chain->hbonds);
	chain->hbonds = 0;
	chain->numHbonds = 0;

	for (i=1; i<chain->length; i++) {
		double score1 = -500, score2 = -500;
		int hbond1 = -1, hbond2 = -1;
		Vector H;
		if (ShortNames[chain->seq->seq[i]] == 'P' ||
			!(C3 = GetAtom(chain, i-1, " C  ")) ||
			!(O3 = GetAtom(chain, i-1, " O  ")) ||
			!(CA = GetAtom(chain, i, " CA ")) ||
			!( N = GetAtom(chain, i,  " N  ")) ||
			!( C = GetAtom(chain, i,  " C  ")) ||
			!( O = GetAtom(chain, i,  " O  "))) continue;

		subVect(&H, &C3->pos, &O3->pos);
		normalizeVect(&H);
		addVect(&H, &H, &N->pos);
		/* Alternative method of predicting hydrogen position.
		 * Results aren't that much different, no clue if they're
		 * any better.  First two normalize operations might well not
		 * be needed, at a slight cost to accuracy.
		 */
		for (j=0; j<chain->length; j++) {
			double result = -9900;
			if (j==i-1 || j==i) continue;
			if (!( CA2 = GetAtom(chain, j,  " CA ")) ||
				!( C2 = GetAtom(chain, j,  " C  ")) ||
				!( O2 = GetAtom(chain, j,  " O  "))) continue;

			if (distSquaredVect(&CA->pos, &CA2->pos) > MaxHDist2) continue;

			if ((H_O2_dist = distSquaredVect(&H, &O2->pos)) >= MinHDist2 &&
				(H_C2_dist = distSquaredVect(&H, &C2->pos)) >= MinHDist2 &&
				(N_O2_dist = distSquaredVect(&N->pos, &O2->pos)) >= MinHDist2 &&
				(N_C2_dist = distSquaredVect(&N->pos, &C2->pos)) >= MinHDist2) {

				H_O2_dist = sqrt(H_O2_dist);
				H_C2_dist = sqrt(H_C2_dist);
				N_O2_dist = sqrt(N_O2_dist);
				N_C2_dist = sqrt(N_C2_dist);

				result = Qconst*(1/H_O2_dist - 1/H_C2_dist + 1/N_C2_dist - 1/N_O2_dist);
			}
			if (result < score2) {
				if (result < score1) {
					score2 = score1;
					hbond2 = hbond1;
					score1 = result;
					hbond1 = j;
				}
				else {
					score2 = result;
					hbond2 = j;
				}
			}
		}
		if (score1 < -500) {
			AddHbond(chain, i, hbond1);
			if (score2 < -500)
				AddHbond(chain, i, hbond2);
		}
	}
}

void AddHbond(PDBChain *chain, int res1, int res2) {
	int i = chain->numHbonds;
	chain->hbonds = (HydrogenBond*)realloc(chain->hbonds, sizeof(HydrogenBond) * (chain->numHbonds+1));

	while (i) {
		if (res1 > chain->hbonds[i-1].res1 ||
			(res1 == chain->hbonds[i-1].res1 && res2 >= chain->hbonds[i-1].res2)) break;
		chain->hbonds[i] = chain->hbonds[i-1];
		i--;
	}
	chain->hbonds[i].res1 = res1;
	chain->hbonds[i].res2 = res2;
	chain->hbonds[i].atom1 = GetAtomIndex(chain, res1, " N  ");
	chain->hbonds[i].atom2 = GetAtomIndex(chain, res2, " O  ");
	chain->numHbonds++;
}

int AreHbonded(PDBChain *chain, int res1, int res2, HydrogenBond **hbond) {
	if (AreHbondedAsymmetric(chain, res1, res2, hbond)) return 1;
	return AreHbondedAsymmetric(chain, res2, res1, hbond);
}

int AreHbondedAsymmetric(PDBChain *chain, int res1, int res2, HydrogenBond **hbond) {
	int min = 0;
	int max = chain->numHbonds-1;
	int mid;
	while (min < max) {
		mid = (min+max)/2;
		if (chain->hbonds[mid].res1 > res1) {
			max = mid-1;
		}
		else if (chain->hbonds[mid].res1 < res1) {
			min = mid+1;
		}
		else if (chain->hbonds[mid].res2 > res2) {
			max = mid-1;
		}
		else if (chain->hbonds[mid].res2 < res2) {
			min = mid+1;
		}
		else {
			min = max = mid;
		}
	}
	if (min == max && chain->hbonds[min].res1 == res1 && chain->hbonds[min].res2 == res2) {
		if (hbond) *hbond = &chain->hbonds[min];
		return 1;
	}
	return 0;
}

void CalcBetaSheets(PDBChain *chain) {
	int hbond;
	free(chain->betaPairs);
	chain->betaPairs = 0;
	chain->numBetaPairs = 0;

	for (hbond = 0; hbond < chain->numHbonds; hbond++) {
		int s1 = chain->hbonds[hbond].res1;
		int s2 = chain->hbonds[hbond].res2;
		int e1 = s1, e2 = s2;
		/* Check if hbond is already in argeement with a
		 * parallel (s2+1) or antiparallel (s2) strand alignment
		 */
		if (AreAlignedPairs(chain, s1, s2) ||
			AreAlignedPairs(chain, s1, s2+1)) continue;

		/* Check for an antiparallel pair of strands. */
		if (AreHbondedAsymmetric(chain, s2, s1, 0)) {
			while (AreHbondedAsymmetric(chain, e1+2, e2-2, 0)) {
				e1++;
				e2--;
				if (!AreHbondedAsymmetric(chain, e2-1, e1+1, 0)) {
					break;
				}
				e1++;
				e2--;
			}
			while (AreHbondedAsymmetric(chain, s2+2, s1-2, 0)) {
				s1--;
				s2++;
				if (!AreHbondedAsymmetric(chain, s1-1, s2+1, 0)) {
					break;
				}
				s1--;
				s2++;
			}
			if (s1 != e1) {
				AddStrandPair(chain, s1, s2, -1, e1-s1+1);
				continue;
			}
		}

		/* Check for a parallel pair of strands. */
		s2++;
		e2 = s2;
		if (AreHbondedAsymmetric(chain, s2+1, s1, 0) &&
			AreHbondedAsymmetric(chain, s1, s2-1, 0)) {
			while (AreHbondedAsymmetric(chain, e1+2, e2+1, 0)) {
				e1++;
				e2++;
				if (!AreHbondedAsymmetric(chain, e2+2, e1+1, 0)) break;
				e1++;
				e2++;
			}
			while (AreHbondedAsymmetric(chain, s2-1, s1-2, 0)) {
				s1--;
				s2--;
				if (!AreHbondedAsymmetric(chain, s1-1, s2-2, 0)) break;
				s1--;
				s2--;
			}
			if (s1 != e1) {
				AddStrandPair(chain, s1, s2, 1, e1-s1+1);
			}
		}
	}
}

void AddAlphaHelix(PDBChain *chain, AlphaHelix *helix) {
	int i = chain->numAlphaHelices;
	chain->alphaHelices = (AlphaHelix*) realloc(chain->alphaHelices, sizeof(AlphaHelix)*(chain->numAlphaHelices+1));
	while (i && chain->alphaHelices[i-1].start > helix->start) {
		chain->alphaHelices[i] = chain->alphaHelices[i-1];
		i--;
	}
	chain->alphaHelices[i] = *helix;
	chain->numAlphaHelices++;
}

void CalcStandardHelixType(PDBChain *chain, int skip, int type) {
	int i;
	for (i=0; i<chain->length-skip; i++) {
		if (AreHbondedAsymmetric(chain, i, i+skip, 0)) {
			int strandlen = 1;
			while (i+strandlen < chain->length-skip && AreHbondedAsymmetric(chain, i+strandlen, i+strandlen+skip, 0)) {
				strandlen++;
			}
			if (strandlen >= 2) {
				AlphaHelix helix;
				helix.start = i;
				helix.length = strandlen+skip-1;
				i+= strandlen-1;
				helix.type = type;
				helix.id[0] = helix.comment[0] = 0;
				AddAlphaHelix(chain, &helix);
				continue;
			}
		}
		if (AreHbondedAsymmetric(chain, i+skip, i, 0)) {
			int strandlen = 1;
			while (i + strandlen < chain->length - skip && AreHbondedAsymmetric(chain, i+strandlen+skip, i+strandlen, 0)) {
				strandlen++;
			}
			if (strandlen >= 2) {
				AlphaHelix helix;
				helix.start = i;
				helix.length = strandlen+skip;
				i+= strandlen-1;
				helix.type = type;
				helix.id[0] = helix.comment[0] = 0;
				AddAlphaHelix(chain, &helix);
				continue;
			}
		}
	}
}

void CalcHelices(PDBChain *chain) {
	int i, j;
	free(chain->alphaHelices);
	chain->alphaHelices = 0;
	chain->numAlphaHelices = 0;
	CalcStandardHelixType(chain, 4, RIGHT_ALPHA);
	CalcStandardHelixType(chain, 3, RIGHT_310);
	CalcStandardHelixType(chain, 5, RIGHT_PI);
	for (i=1; i<chain->numAlphaHelices; i++) {
		AlphaHelix helix = chain->alphaHelices[i];
		j = i;
		while (j && chain->alphaHelices[j-1].start > helix.start) {
			chain->alphaHelices[j] = chain->alphaHelices[j-1];
			j--;
		}
		chain->alphaHelices[j] = helix;
	}
}

int CompareDirections(PDBChain *chain, int p1, int p2, int p3, double *out) {
	Atom *a1 = GetAtom(chain, p1, " CA ");
	Atom *a2 = GetAtom(chain, p2, " CA ");
	Atom *a3 = GetAtom(chain, p3, " CA ");
	if (a1 && a2 && a3) {
		Vector v1, v2;
		subVect(&v1, &a2->pos, &a1->pos);
		subVect(&v2, &a3->pos, &a2->pos);
		normalizeVect(&v1);
		normalizeVect(&v2);
		*out = dotVect(&v1, &v2);
		return 1;
	}
	return 0;
}

void AssembleSheets(PDBChain *chain) {
	int i, j;
	int added;
	char *temp = (char*) calloc(chain->numBetaPairs, sizeof(char));
	BetaPair pairs[4];

	while (chain->numBetaSheets) {
		free(chain->betaSheets[--chain->numBetaSheets].strands);
	}
	free(chain->betaSheets);
	chain->betaSheets = 0;


	added = 0;
	for (i = 0; i<chain->numBetaPairs; i++) {
		int left = 0, right = 0;
		BetaSheet *sheet;
		if (temp[i]) continue;
		temp[i] = i+1;
		chain->betaSheets = (BetaSheet*) realloc(chain->betaSheets, sizeof(BetaSheet) * (chain->numBetaSheets+1));
		sheet = chain->betaSheets + chain->numBetaSheets;
		chain->numBetaSheets++;
		sheet->id[0] = 0;
		sheet->numStrands = 1;
		sheet->strands = (BetaSheetPair*) malloc(sizeof(BetaSheetPair));
		sheet->strands[0].direction = chain->betaPairs[i].direction;
		sheet->strands[0].length = chain->betaPairs[i].length;
		sheet->strands[0].res1 = chain->betaPairs[i].start1;
		sheet->strands[0].res2 = chain->betaPairs[i].start2;
		pairs[0] = pairs[1] = chain->betaPairs[i];
		FlipPair(pairs+1);
		while (left >= 0 || right >= 0) {
			int p;
			left = right = -1;
			for (p = 0; p < chain->numBetaPairs; p++) {
				BetaPair *pair = &chain->betaPairs[p];
				BetaResiduePair res;
				if (temp[p] == i+1) continue;
				for (j = 0; j<pairs[0].length; j++) {
					if (FindFirstMatch(pair, pairs[0].start1 + j, &res)) {
						if (left<0 || (temp[left] && !temp[j])) {
							BetaPair tempPair = *pair;
							double d;
							FlipPair(&tempPair);
							if (CompareDirections(chain,
													pairs[0].start2 + j*pairs[0].direction,
													pairs[0].start1 + j,
													tempPair.start2 + j*tempPair.direction, &d) && d > 0) {
								left = p;
								pairs[2] = tempPair;
							}
						}
					}
					else if (FindSecondMatch(pair, pairs[0].start1 + j, &res)) {
						if (left<0 || (temp[left] && !temp[j])) {
							double d;
							if (CompareDirections(chain,
													pairs[0].start2 + j*pairs[0].direction,
													pairs[0].start1 + j,
													pair->start2 + j*pair->direction, &d) && d > 0) {
								left = p;
								pairs[2] = *pair;
							}
						}
					}
				}
				for (j = 0; j<pairs[1].length; j++) {
					if (FindFirstMatch(pair, pairs[1].start1 + j, &res)) {
						if (left<0 || (temp[left] && !temp[j])) {
							double d;
							if (CompareDirections(chain,
													pairs[1].start2 + j*pairs[1].direction,
													pairs[1].start1 + j,
													pair->start2 + j*pair->direction, &d) && d > 0) {
								right = p;
								pairs[3] = *pair;
							}
						}
					}
					else if (FindSecondMatch(pair, pairs[1].start1 + j, &res)) {
						if (left<0 || (temp[left] && !temp[j])) {
							BetaPair tempPair = *pair;
							double d;
							FlipPair(&tempPair);
							if (CompareDirections(chain,
													pairs[1].start2 + j*pairs[1].direction,
													pairs[1].start1 + j,
													tempPair.start1 + j*tempPair.direction, &d) && d > 0) {
								right = p;
								pairs[3] = tempPair;
							}
						}
					}
				}
			}
			if (left >= 0) {
				sheet->strands = (BetaSheetPair*) realloc(sheet->strands, sizeof(BetaSheetPair) * (sheet->numStrands+1));
				memmove(sheet->strands + 1, sheet->strands, sheet->numStrands * sizeof(BetaSheetPair));
				sheet->strands[0].direction = pairs[2].direction;
				sheet->strands[0].length = pairs[2].length;
				sheet->strands[0].res1 = pairs[2].start1;
				sheet->strands[0].res2 = pairs[2].start2;
				pairs[0] = pairs[2];
				temp[left] = i+1;
				sheet->numStrands++;
			}
			if (right >= 0) {
				sheet->strands = (BetaSheetPair*) realloc(sheet->strands, sizeof(BetaSheetPair) * (sheet->numStrands+1));
				sheet->strands[sheet->numStrands].direction = pairs[3].direction;
				sheet->strands[sheet->numStrands].length = pairs[3].length;
				sheet->strands[sheet->numStrands].res1 = pairs[3].start1;
				sheet->strands[sheet->numStrands].res2 = pairs[3].start2;
				pairs[1] = pairs[3];
				FlipPair(pairs+1);
				temp[right] = i+1;
				sheet->numStrands++;
			}
		}
	}
	free(temp);
}

void AddStrandPair(PDBChain *chain, int start1, int start2, int direction, int length) {
	chain->betaPairs = (BetaPair*) realloc(chain->betaPairs, sizeof(BetaPair) * (chain->numBetaPairs+1));
	/* want first strand in sequence to appear first in the pair. */
	chain->betaPairs[chain->numBetaPairs].direction = direction;
	chain->betaPairs[chain->numBetaPairs].length = length;
	if (start1 < start2) {
		chain->betaPairs[chain->numBetaPairs].start1 = start1;
		chain->betaPairs[chain->numBetaPairs].start2 = start2;
	}
	else if (direction == 1) {
		chain->betaPairs[chain->numBetaPairs].start1 = start2;
		chain->betaPairs[chain->numBetaPairs].start2 = start1;
	}
	/* Antiparallel strands work backwards from start2 and forward from start,
	 * so have to modify them.
	 */
	else {
		chain->betaPairs[chain->numBetaPairs].start1 = start2-length+1;
		chain->betaPairs[chain->numBetaPairs].start2 = start1+length-1;
	}
	chain->numBetaPairs++;
}

int AreAlignedPairs(PDBChain *chain, int res1, int res2) {
	int pair;
	if (res2 < res1) {
		int temp = res1;
		res1 = res2;
		res2 = temp;
	}

	for (pair = 0; pair < chain->numBetaPairs; pair++) {
		if (chain->betaPairs[pair].start1 <= res1 && res1 < chain->betaPairs[pair].start1 + chain->betaPairs[pair].length) {
			int dist = chain->betaPairs[pair].direction*(res1 - chain->betaPairs[pair].start1);
			if (chain->betaPairs[pair].start2 + dist == res2) return chain->betaPairs[pair].direction;
		}
	}
	return 0;
}

void Renumber(PDBChain *chain) {
	int i;
	for (i=0; i<chain->length; i++) {
		chain->residues[i].resLetter = ' ';
		chain->residues[i].resNum = i+1;
	}
}

PDBChain *DuplicateChain(PDBChain *chain) {
	PDBChain *out = (PDBChain*)calloc(1, sizeof(PDBChain));
	int i;
	out->numAtoms = chain->numAtoms;
	out->atoms = (Atom*)memdup(chain->atoms, sizeof(Atom) * chain->numAtoms);
	out->seq = CreateSequence(chain->seq->name);
	for (i=0; i<chain->seq->length; i++) {
		AddCharacter(out->seq, chain->seq->seq[i]);
	}
	out->length = chain->length;
	out->residues = (Residue*) memdup(chain->residues, sizeof(Residue) * out->length);
	out->idString = strdup(chain->idString);
	out->chainName = chain->chainName;
	return out;
}

