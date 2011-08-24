/* Extend.c -- Contains functions to extend an alignment.
 * Includes both the extension pase run every round and the final pass.
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

#include "Extend.h"
#include <malloc.h>

#include "Score.h"

void Extend(MultipleAlignment *ma, Matrix **matrices, int useSkips) {
	int block, block2, chain, chain2, i;
	if (ma->numChains == 1) {
		return;
	}
	for (block=0; block<=ma->numBlocks; block++) {
		int skips = 0;

		while (1) {
				double bestDist = MAX_EXTEND_DIST;
				int best = -1;

				/* 0 is backwards, 1 is forwards.
				 * 0 starts with element before high, etc.
				 */
				int minDir = 0;
				int maxDir = 1;

				int overlap = 0;

				int low, high;
				if (block < ma->numBlocks) {
					high = ma->blocks[block].first;
					for (chain = 0; chain<ma->numChains; chain++) {
						if (ma->residues[chain].res[high].exists == 0 ||
							ma->residues[chain].res[high].index == 0 ||
							ma->chains[chain]->res[ma->residues[chain].res[high].index-1].exists == 0) {

							minDir = 1;
							break;
						}
					}
				}
				else minDir = 1;

				if (block) {
					low = ma->blocks[block-1].last;
					for (chain = 0; chain<ma->numChains; chain++) {
						if (ma->residues[chain].res[low].exists == 0 ||
							ma->residues[chain].res[low].index+1 == ma->chains[chain]->length ||
							ma->chains[chain]->res[ma->residues[chain].res[low].index+1].exists == 0) maxDir = 0;
					}

				}
				else maxDir = 0;

				if (block && block < ma->numBlocks) {
					if (skips >= 5) break;
					for (chain = 0; chain<ma->numChains; chain++) {
						if (ma->residues[chain].res[low].index + 1 >= ma->residues[chain].res[high].index) {
							if (!useSkips) {
								maxDir = 0;
								minDir = 1;
							}
							else
								overlap = 1;
							break;
						}
					}
				}

				for (i=minDir; i<=maxDir; i++) {
					int delta = -1+2*i;
					int pos;
					double dst = 0;
					double weight = 0;
					if (i == 0) pos = high;
					else pos = low;
					for (chain = 0; chain<ma->numChains; chain++) {
						Vector v1, v2;
						int r1 = ma->residues[chain].res[pos].index+delta;
						transformVect(&v1, &matrices[chain][block-i], &ma->chains[chain]->res[r1].coords);
						for (chain2 = chain+1; chain2<ma->numChains; chain2++) {
							int r2 = ma->residues[chain2].res[pos].index+delta;
							double w;
							transformVect(&v2, &matrices[chain2][block-i], &ma->chains[chain2]->res[r2].coords);
							w = ma->residues[chain].weight * ma->residues[chain2].weight;
							weight += w;
							dst += w * distVect(&v1, &v2);
						}

					}
					dst /= weight;
					if (dst < bestDist) {
						bestDist = dst;
						best = i;
					}
				}

				if (best<0) break;
				skips += overlap;
				ma->blocks[block-best].last++;
				for (block2 = block-best+1; block2<ma->numBlocks; block2++) {
					ma->blocks[block2].first++;
					ma->blocks[block2].last++;
				}

				for (chain=0; chain<ma->numChains; chain++) {
					int pos;
					int delta = -1+2*best;
					int r1;
					Vector v;
					if (best == 0) pos = high;
					else pos = low;

					r1 = ma->residues[chain].res[pos].index+delta;
					ma->residues[chain].res = (ResiduePosition*) realloc(ma->residues[chain].res, (ma->numResidues+1) * sizeof(ResiduePosition));
					memmove(ma->residues[chain].res+(pos+best+1), ma->residues[chain].res+(pos+best), (ma->numResidues-pos-best)*sizeof(ResiduePosition));
					ma->residues[chain].res[pos+best] = ma->chains[chain]->res[r1];
					v = ma->residues[chain].res[pos+best].coords;
					transformVect(&ma->residues[chain].res[pos+best].coords, &matrices[chain][block-best], &v);
				}
				ma->numResidues++;

		}

	}
}

static void CalcBlockExtents(MultipleAlignment *ma, int block, int *starts, int *ends) {
	int  chain;
	for (chain = 0; chain<ma->numChains; chain++) {
		if (block)
			starts[chain] = ma->residues[chain].res[ma->blocks[block-1].last].index+1;
		else
			starts[chain] = 0;
		if (block<ma->numBlocks)
			ends[chain] = ma->residues[chain].res[ma->blocks[block].first].index-1;
		else
			ends[chain] = ma->chains[chain]->pdb->length-1;
	}
}

static double ScoreSet(int *out, MultipleAlignment *ma, int length) {
	double sum = 0;
	int chain, chain2, i, div=0;
	for (chain = 0; chain<ma->numChains; chain++) {
		if (out[chain] < 0) continue;
		for (chain2 = chain+1; chain2<ma->numChains; chain2++) {
			double d = 0;
			if (out[chain2] < 0) continue;
			for (i=0; i<length; i++) {
				ResiduePosition *p1 = &ma->chains[chain]->res[out[chain]+i];
				ResiduePosition *p2 = &ma->chains[chain2]->res[out[chain2]+i];
				Vector v;
				d += lengthSquaredVect(subVect(&v, &p1->coords, &p2->coords));
			}
			d = sqrt(d/length);
			sum += d;
			div++;
		}
	}
	return sum / div;
}
/* Greedily finds a set of adjacent residues in all other chains closest to the specified set of residues.
 * Could find a tighter group, in theory, with another pass or two, but this allows a bit more variety,
 * possibly useful in the pseudo-dynamic partial alignment pass.
 *
 * full = 2 forces residues from each chain.
 * full = 0 forces at least one chain to have no residues.
 * full = 1 isn't used, but allows both full and partial sets of residues.
 */
static double AssembleSet(int *out, MultipleAlignment *ma, int *starts, int *ends, int chain, int residue, double rmsd, int length, int full) {
	double sum = 0;
	double worstScore = -1;
	int worstChain = -1;
	/* Doesn't include initial chain. */
	int numChainsWithPositions = 0;
	int chain2, residue2, i;
	int div = 0;
	double test = 10e10;

	for (residue2=residue; residue2<length+residue;residue2++) {
		if (!ma->chains[chain]->res[residue2].exists) return 10e10;
	}

	if (full != 2) {
		test = rmsd * rmsd * length;
	}

	out[chain] = residue;
	for (chain2 = 0; chain2<ma->numChains; chain2++) {
		double best;
		if (chain2 == chain) continue;
		best = test;
		out[chain2] = -1;
		for (residue2 = starts[chain2]; residue2<=ends[chain2]-length+1;residue2++) {
			double d = 0;
			for (i=0; i<length; i++) {
				Vector v;
				ResiduePosition *p1 = &ma->chains[chain]->res[residue+i];
				ResiduePosition *p2 = &ma->chains[chain2]->res[residue2+i];

				if (!p2->exists) {
					d = 1000*length*ma->numChains;
					break;
				}
				d += lengthSquaredVect(subVect(&v, &p1->coords, &p2->coords));
			}
			if (d<best) {
				out[chain2] = residue2;
				best = d;
			}
		}
		if (out[chain2] >= 0) {
			numChainsWithPositions++;
			sum += sqrt(best/length);
			if (worstScore < best) {
				worstScore = best;
				worstChain = chain2;
			}
		}

	}

	/* Simple sanity check.  Note that sum is only based on one structure, but if it's too much higher than the
	 * real RMSD, can use a different base residue to get a better sum.
	 * Only fail second test when full is 2.
	 */
	if (!numChainsWithPositions || sum > 2 * numChainsWithPositions * rmsd) return 10e100;

	if (full == 2 && numChainsWithPositions != ma->numChains - 1) return 10e100;

	if (numChainsWithPositions == ma->numChains - 1 &&
		full == 0) {
			/* Shouldn't happen, the way the code is structured.
			 */
			if (ma->numChains == 2) {
				return 10e10;
			}
			out[worstChain] = -1;
			numChainsWithPositions --;
	}

	return ScoreSet(out, ma, length);
}

void FinalExtend(MultipleAlignment *ma, double rmsd) {
	int *ends = (int*) malloc(sizeof(int)*ma->numChains);
	int *starts = (int*) malloc(sizeof(int)*ma->numChains);
	int *bestStarts = (int*) malloc(sizeof(int)*ma->numChains);
	int *current = (int*) malloc(sizeof(int)*ma->numChains);
	int i;
	int chain, chain2, chain3, pos, pos2;
	int length, block;
	for (length = 4; length>=1; length--) {
		for (block = 0; block<=ma->numBlocks; block++) {
			double best = length*ma->numChains*100000;
			ResiduePosition *fuck = ma->residues[0].res;
			for (chain = 0; chain<ma->numChains; chain++) {
				if (block)
					starts[chain] = ma->residues[chain].res[ma->blocks[block-1].last].index+1;
				else
					starts[chain] = 0;
				if (block<ma->numBlocks)
					ends[chain] = ma->residues[chain].res[ma->blocks[block].first].index-1;
				else
					ends[chain] = ma->chains[chain]->pdb->length-1;
			}
			for (chain = 0; chain<ma->numChains; chain++) {
				for (pos = starts[chain]; pos<=ends[chain]-length+1;pos++) {
					double sum = 0;
					current[chain] = pos;
					for (chain2 = 0; chain2<ma->numChains; chain2++) {
						double best2;
						if (chain2 == chain) continue;
						best2 = 100000*length*ma->numChains;
						for (pos2 = starts[chain2]; pos2<=ends[chain2]-length+1;pos2++) {
							double d = 0;
							for (i=0; i<length; i++) {
								Vector v;
								ResiduePosition *p1 = &ma->chains[chain]->res[pos+i];
								ResiduePosition *p2 = &ma->chains[chain2]->res[pos2+i];

								if (!p1->exists || !p2->exists) {
									d = 1000*length*ma->numChains;
									break;
								}
								d += lengthVect(subVect(&v, &p1->coords, &p2->coords));
							}
							if (d<best2) {
								current[chain2] = pos2;
								best2 = d;
							}
						}
						sum += best2;
					}
					if (sum/(ma->numChains-1) < rmsd * 2) {
						int div = 0;
						sum = 0;
						for (chain2 = 0; chain2<ma->numChains; chain2++) {
							for (chain3 = chain2+1; chain3<ma->numChains; chain3++) {
								double d = 0;
								for (i=0; i<length; i++) {
									ResiduePosition *p1 = &ma->chains[chain2]->res[current[chain2]+i];
									ResiduePosition *p2 = &ma->chains[chain3]->res[current[chain3]+i];
									Vector v;
									d += lengthVect(subVect(&v, &p1->coords, &p2->coords));
								}
								d = d/length;
								sum += d;
								div++;
							}
						}
						sum /= div;
						if (sum < best) {
							best = sum;
							memcpy(bestStarts, current, sizeof(int)*ma->numChains);
						}
					}
				}
			}
			if (best > rmsd) continue;
			else {
				AlignedBlock *g;
				for (i=block; i<ma->numBlocks;i++) {
					ma->blocks[i].first+=length;
					ma->blocks[i].last+=length;
				}
				ma->blocks = (AlignedBlock *) realloc(ma->blocks, (ma->numBlocks+1) * sizeof(AlignedBlock));
				memmove(ma->blocks+block+1, ma->blocks+block, (ma->numBlocks-block)*sizeof(AlignedBlock));
				ma->numBlocks++;

				g = ma->blocks+block;;
				if (!block)
					g->first = 0;
				else
					g->first = ma->blocks[block-1].last+1;
				g->last = g->first+length-1;

				for (chain=0; chain<ma->numChains; chain++) {
					ma->residues[chain].res = (ResiduePosition*) realloc(ma->residues[chain].res, sizeof(ResiduePosition)*(ma->numResidues+length));
					memmove(ma->residues[chain].res + g->first+length, ma->residues[chain].res+g->first, sizeof(ResiduePosition)*(ma->numResidues-g->first));
					for (i=0; i<length; i++) {
						ma->residues[chain].res[g->first+i] = ma->chains[chain]->res[bestStarts[chain]+i];
					}
				}
				ma->numResidues += length;

				block--;
			}
		}
	}
	free(ends);
	free(starts);
	free(bestStarts);
	free(current);
}

struct RESIDE_SET {
	double score;
	double combinedScore;
	double sortBy;
	int *indices;
	/* Counts of all unused sets. */
	int numSetsBefore;
	int numSetsAfter;

	/* Counts of selected set of sets. */
	int numPrevSets;
	int bestPrevSet;
};
typedef struct RESIDE_SET ResidueSet;

double CalcSetRMSD(MultipleAlignment *ma, ResidueSet *sets1, ResidueSet *sets2, int *c) {
	int i, j;
	double sum = 0;
	int count = 0;
	for (i=0; i<ma->numChains; i++) {
		int index1 = sets1->indices[i];
		Vector *v1;
		if (index1 < 0) continue;
		v1 = &ma->chains[i]->res[index1].coords;
		for (j=0; j<ma->numChains; j++) {
			int index2 = sets2->indices[j];
			Vector v;
			if (index2 < 0) continue;
			subVect(&v, v1, &ma->chains[j]->res[index2].coords);
			sum += lengthVect(&v);
			count ++;
		}
	}
	*c = count;
	return sum / count;
}

struct SET_ALIGNMENT_INFO {
	double score;
	int alignedTo;
	int pairCount;
};

typedef struct SET_ALIGNMENT_INFO SetAlignmentInfo;


void AlignSets(MultipleAlignment *ma, int *starts, int *ends, ResidueSet **sets, int *numSets, ResidueSet *sets1, int numSets1, ResidueSet *sets2, int numSets2, double rmsd) {
	int s1, s2, i, j;
	int max = 0;
	SetAlignmentInfo * info;

	int mergeCount = 0;
	info = (SetAlignmentInfo*) malloc(numSets1 * numSets2 * sizeof(SetAlignmentInfo));
	if (numSets1 && numSets2) {
		SetAlignmentInfo *s1Info, *lastS1Info;
		for (s1=0; s1<numSets1; s1++) {
			s1Info = info + s1*numSets2;
			lastS1Info = s1Info - numSets2;
			for (s2=0; s2<numSets2; s2++) {
				int count;
				double setRmsd = CalcSetRMSD(ma, sets1+s1, sets2+s2, &count);
				SetAlignmentInfo *inf = s1Info + s2;
				inf->alignedTo = 0;
				inf->score = 0;
				inf->pairCount = 0;
				/* Aligned to nothing does best. */
				if (s1) {
					*inf = lastS1Info[s2];
					inf->alignedTo = -1;
				}
				/* Aligned to something before s2 does best. */
				if (s2 && s1Info[s2-1].score > inf->score) {
					*inf = s1Info[s2-1];
				}
				/* Check for new aligned set. */
				if (setRmsd < rmsd) {
					double score = count*(rmsd - setRmsd);
					/* Can add in previous residues. */
					int pc = 0;
					if (s1 && s2) {
						pc = lastS1Info[s2-1].pairCount;
						score += lastS1Info[s2-1].score;
					}
					if (score > inf->score) {
						inf->pairCount = pc+1;
						inf->score = score;
						inf->alignedTo = s2;
					}
				}
			}
		}
		mergeCount = s1Info[numSets2-1].pairCount;
	}
	s1 = numSets1 - 1;
	s2 = numSets2 - 1;
	numSets[0] = numSets1 + numSets2 - mergeCount;
	*sets = (ResidueSet*)malloc((sizeof(ResidueSet) + sizeof(int) * ma->numChains) * numSets[0]);
	for (i = numSets[0]-1; i>=0; i--) {
		SetAlignmentInfo *inf = &info[s1*numSets2 + s2];
		sets[0][i].indices = ((int*)(sets[0] + numSets[0])) + i * ma->numChains;
		if (s1 == -1 || (s2 == 0 && inf[0].score == 0) || (s2 > 0 && inf[-1].score == inf[0].score)) {
			memcpy(sets[0][i].indices, sets2[s2].indices, ma->numChains * sizeof(int));
			s2--;
			continue;
		}
		if (s2 == -1 || inf->alignedTo == -1) {
			memcpy(sets[0][i].indices, sets1[s1].indices, ma->numChains * sizeof(int));
			s1--;
			continue;
		}
		memcpy(sets[0][i].indices, sets1[s1].indices, ma->numChains * sizeof(int));
		for (j=0; j<ma->numChains; j++) {
			if (sets2[s2].indices[j] != -1) {
				sets[0][i].indices[j] = sets2[s2].indices[j];
			}
		}
		s2--;
		s1--;
	}

	free(sets1);
	free(sets2);
	free(info);
}

void SplitAndJoinSets(MultipleAlignment *ma, int *starts, int *ends, int *splitChains, ResidueSet **sets, int *numSets, double rmsd) {
	ResidueSet *sets1, *sets2;
	int numSets1, numSets2;
	int s1, j, k;
	s1 = numSets[0];
	sets1 = sets[0];
	numSets2 = 0;
	numSets1 = 0;

	sets2 = (ResidueSet*)malloc((sizeof(ResidueSet) + sizeof(int) * ma->numChains) * s1);
	for (j = 0; j<s1; j++) {
		int addS1 = 0, addS2 = 0;

		sets2[numSets2].indices = ((int*)(sets2 + s1)) + numSets2 * ma->numChains;

		for (k=ma->numChains-1; k >= 0; k--) {
			sets2[numSets2].indices[k] = -1;
			if (sets1[j].indices[k] != -1) {
				if (splitChains[k]) {
					addS2 = 1;
					sets2[numSets2].indices[k] = sets1[j].indices[k];
					sets1[j].indices[k] = -1;
				}
				else {
					addS1 = 1;
				}
			}
		}
		numSets2 += addS2;
		sets1[numSets1] = sets1[j];
		numSets1 += addS1;
	}
	AlignSets(ma, starts, ends, sets, numSets, sets1, numSets1, sets2, numSets2, rmsd);
}

void Merge(MultipleAlignment *ma, int *starts, int *ends, AssemblyNode *o, ResidueSet **sets, int *numSets, double rmsd) {
	if (o->id >= 0) {
		int index = o->id;
		int i;
		int max = 1+ends[index]-starts[index];
		if (max) {
			max = max;
		}
		*numSets = 0;
		*sets = (ResidueSet*)malloc((sizeof(ResidueSet) + sizeof(int) * ma->numChains) * max);
		memset(*sets, -1, (sizeof(ResidueSet) + sizeof(int) * ma->numChains) * max);
		for (i=0; i<max; i++) {
			if (!ma->chains[index]->res[i + starts[index]].exists) continue;
			sets[0][numSets[0]].indices = ((int*)(sets[0] + max)) + i * ma->numChains;
			sets[0][numSets[0]].indices[index] = starts[index]+i;
			numSets[0]++;
		}
	}
	else {
		ResidueSet *sets1, *sets2;
		int numSets1, numSets2;
		int s1, s2, i, c1, c2;
		int max = 0;
		int *splitChains;
		int *chains;

		Merge(ma, starts, ends, o->left, &sets1, &numSets1, rmsd);
		Merge(ma, starts, ends, o->right, &sets2, &numSets2, rmsd);

		AlignSets(ma, starts, ends, sets, numSets, sets1, numSets1, sets2, numSets2, rmsd);

		splitChains = (int*) malloc(sizeof(int) * ma->numChains*2);
		chains  = splitChains + ma->numChains;
		memset(chains, 0, sizeof(int) * ma->numChains);

		/* Used for minor optimization. */
		for (s1=0; s1<numSets[0]; s1++) {
			for (c1 = 0; c1 < ma->numChains; c1++) {
				chains[c1] |= (sets[0][s1].indices[c1] >= 0);
			}
		}

		c2 = 0;
		for (c1 = 0; c1 < ma->numChains; c1++) {
			c2 += chains[c1];
		}
		if (c2 > 2) {
			for (i=0; i<10; i++) {
				int count = numSets[0];
				/* Slowest part of partial alignment algorithm. Try to realign all single chains and all pairs jointly. */
				memset(splitChains,0,sizeof(int) * ma->numChains);
				for (c1 = 0; c1 < ma->numChains; c1++) {
					if (!chains[c1]) continue;
					splitChains[c1] = 1;
					SplitAndJoinSets(ma, starts, ends, splitChains, sets, numSets, rmsd);
					for (c2 = c1+1; c2<ma->numChains; c2++) {
						if (!chains[c2]) continue;
						splitChains[c2] = 1;

						SplitAndJoinSets(ma, starts, ends, splitChains, sets, numSets, rmsd);
						splitChains[c2] = 0;
					}
					splitChains[c1] = 0;
				}
				for (s1=0; s1<numSets[0]; s1++) {
					int count = 0;
					int different = 0;
					/* Check if anything to do. */
					for (c1 = 0; c1 < ma->numChains; c1++) {
						if (sets[0][s1].indices[c1] == -1) {
							different += splitChains[c1];
						}
						else {
							different += 1-splitChains[c1];
							count ++;
						}
					}
					if (count <= 2 || !different) continue;
					for (c1 = 0; c1 < ma->numChains; c1++) {
						splitChains[c1] = (sets[0][s1].indices[c1] != -1);
					}

					s2 = numSets[0];
					SplitAndJoinSets(ma, starts, ends, splitChains, sets, numSets, rmsd);
					if (numSets[0] < s2) {
						s1 -= s2 - numSets[0];
						if (s1 < -1) s1 = -1;
					}
				}

				/* If things haven't changed much, stop. */
				if (numSets[0] == count) break;
			}
		}
		free (splitChains);
	}
}

void CalcPartial(MultipleAlignment *ma, double rmsd) {
	int block;


	int *starts = (int*) calloc(ma->numChains*4, sizeof(int));
	int *ends = starts+ma->numChains;
	int *indices = ends+ma->numChains;
	int *indices2 = indices+ma->numChains;
	int resIndex = 0;
	int s1, chain;

	int div = (ma->numChains)*(ma->numChains-1);

	ResidueSet *sets = 0;

	/* Account for roundoff and make sure nothing's added to the final alignment. */
	rmsd -= 0.00001;

	for (block=0; block <= ma->numBlocks; block++) {
		int numSets;

		/* Calculate residue ranges. */
		CalcBlockExtents(ma, block, starts, ends);

		sets = 0;
		numSets = 0;
		Merge(ma, starts, ends, ma->order->root, &sets, &numSets, rmsd);


		if (numSets) {
			int s2;
			s2 = 0;
			for (s1=0; s1<numSets; s1++) {
				int count = 0;
				for (chain=0; chain<ma->numChains; chain++) {
					count += (sets[s1].indices[chain] >= 0);
				}
				if (count >= 2) {
					sets[s2++] = sets[s1];
				}
			}

			numSets = s2;
			ma->blocks = realloc(ma->blocks, (ma->numBlocks + numSets) * sizeof(AlignedBlock));
			if (numSets) {
				for (s1 = ma->numBlocks-1; s1 >= block; s1 --) {
					ma->blocks[s1+numSets] = ma->blocks[s1];
					ma->blocks[s1+numSets].first += numSets;
					ma->blocks[s1+numSets].last += numSets;
				}
				if (block == ma->numBlocks) {
					ma->blocks[block].first = ma->blocks[block-1].last+1;
				}
				ma->blocks[block].last = ma->blocks[block].first;
				for (s1 = block+1; s1 < block + numSets; s1++) {
					ma->blocks[s1].first = ma->blocks[s1].last = ma->blocks[s1-1].first+1;
				}

				for (chain = 0; chain < ma->numChains; chain++) {
					ma->residues[chain].res = (ResiduePosition*)realloc(ma->residues[chain].res, (ma->numResidues+numSets)*sizeof(ResiduePosition));
					memmove(ma->residues[chain].res + ma->blocks[block].first + numSets,
							ma->residues[chain].res + ma->blocks[block].first,
							sizeof(ResiduePosition) * (ma->numResidues-ma->blocks[block].first));
					for (s1 = 0; s1<numSets; s1++) {
						ResiduePosition *pos = &ma->residues[chain].res[ma->blocks[block].first + s1];
						if (sets[s1].indices[chain] == -1) {
							pos->exists = 0;
						}
						else {
							*pos = ma->chains[chain]->res[sets[s1].indices[chain]];
						}
					}
				}
				ma->numBlocks += numSets;
				ma->numResidues += numSets;
				block += numSets;
			}
		}

		free(sets);
	}
	free(starts);
}
