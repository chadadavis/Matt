/* MultipleAlignment.h -- Contains prototypes for creating and manipulating multiple
 * alighments.
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

#ifndef MULTIPLE_ALIGNMENT_H
#define MULTIPLE_ALIGNMENT_H

#include "pdb.h"
#include "AssemblyOrder.h"


/* Use these around instead of the PDB files.  Makes things simpler and faster,
 * as my PDB structure is rather ugly.
 */
struct RESIDUE_POSITION {
	Vector coords;
	union {
		int hasCoords;
		int exists;
	};
	int index;
	char residue;
};
typedef struct RESIDUE_POSITION ResiduePosition;

/* Used to keep track of source PDB files and original atomic coordinates */
struct RESIDUE_POSITIONS {
	ResiduePosition *res;
	PDBChain *pdb;
	int length;
	int id;
};
typedef struct RESIDUE_POSITIONS ResiduePositions;

/* Contains aligned residue positions, after 3D transformation */
struct WEIGHTED_RESIDUE_POSITIONS {
	ResiduePosition *res;

	/* NOTE: Weighting is currently used inconcistently and everything's
	 * always given the same weight, except when excluding a single
	 * structurein the realignment phase.
	 */
	double weight;
};
typedef struct WEIGHTED_RESIDUE_POSITIONS WeightedResiduePositions;

/* First and last indices in WeightedResiduePositions of a block */
struct ALIGNED_BLOCK {
	int first;
	int last;
};
typedef struct ALIGNED_BLOCK AlignedBlock;

struct MULTIPLE_ALIGNMENT {
	/* Order of the structure assembly, as a tree allocated in a single block of memory. */
	AssemblyOrder *order;

	/* Raw score. */
	double score;

	double rmsd;
	double pvalue;

	ResiduePositions **chains;
	int numChains;
	int numResidues;

	/* This actually has the alignment */
	WeightedResiduePositions *residues;

	/* Which blocks were aligned.  Residues from the intermediate extension phase
	 * are added to end of blocks, while residues from the final extension phase
	 * are given their own blocks.
	 */
	AlignedBlock *blocks;
	int numBlocks;

	/* List of conflicts between adjacent blocks. Only used by primary dynamic
	 * programming algorithm (Not by any extension routine).
	 */
	int** conflictMap;

	/* Weighted averaged positions.  Used to estimate displacements. */
	Vector* averages;
};
typedef struct MULTIPLE_ALIGNMENT MultipleAlignment;





/* Calculates RMSD using current structure alignment. */
double CalcAlignmentRMSD(MultipleAlignment *ma);

/* Only works correctly when only two structures are aligned and CalcNonWeightedRMSD
 * was called, but no final extension phase has been run.
 */
double CalcAlignmentPvalue(MultipleAlignment *ma);

/* Transforms the PDB files using the current (bent) alignment. Only modifies
 * atom records.
 */
void TransformPDB(MultipleAlignment *ma);

void CleanupAlignment(MultipleAlignment *ma);

/* Calculates the unbent alignment of the structures. Does not affect original PDBs. */
void UnbentRMSDAlign(MultipleAlignment *ma);

/* Main function */
MultipleAlignment* Align(PDBChain **chains, int numChains, int displayStatus);



#endif


